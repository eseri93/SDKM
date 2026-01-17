library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(quanteda)
library(R.matlab)
library(stringr)
library(tidytext)
library(ggrepel)

wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

data<-data_corpus_inaugural
data<-tolower(data)
data_tokens <- tokens(data, remove_punct = TRUE, remove_symbols = T, remove_separators = T, remove_numbers = TRUE) 
data_tokens <-tokens_replace(data_tokens, pattern = lexicon::hash_lemmas$token, replacement = lexicon::hash_lemmas$lemma) #lemmatizzazione
tokens_to_remove <- unique(unlist(tokens(data_tokens))) # Convert tokens to a vector of unique tokens
tokens_to_remove <- tokens_to_remove[str_length(tokens_to_remove) < 3 | str_length(tokens_to_remove) > 14]
data_tokens<- tokens_remove(data_tokens, pattern = tokens_to_remove)
data_dfm <- dfm(data_tokens)
data_dfm_nostops <- dfm_remove(data_dfm, pattern = stopwords("en")) 
topfeatures(data_dfm_nostops, 20)
data_dfm_nostops_trim<-dfm_trim(data_dfm_nostops, min_termfreq = 11)
dim(data_dfm_nostops_trim)
data_dfm_nostops_trim<-dfm_tfidf(data_dfm_nostops_trim, force = TRUE)

data_text <- convert(data_dfm_nostops_trim, to = "data.frame")
rownames(data_text)<-data_text[,1]
data_text<-data_text[,-1]
data_text<-as.data.frame(t(data_text))

#write.csv(data_text, "data_text.csv")
#save(data_text, file="data_text.Rdata")

# Load cluster assignments from MATLAB
wordClusterData <- readMat("wordClusterIdx.mat")
docClusterData <- readMat("docClusterIdx.mat")

# Extract cluster assignments
wordClusterIdx <- as.vector(unlist(wordClusterData$wordClusterIdx))
docClusterIdx <- as.vector(unlist(docClusterData$docClusterIdx))


#Distribution of Words Across Clusters----------------------------------------------------------

# Create a data frame for word clusters
word_clusters_df <- data.frame(Word = rownames(data_text), Cluster = factor(wordClusterIdx))

# Count words per cluster
cluster_counts <- word_clusters_df %>%
  group_by(Cluster) %>%
  summarise(Number_of_Words = n())

# Order clusters by number of words
cluster_counts <- cluster_counts %>%
  arrange(desc(Number_of_Words)) %>%
  mutate(Cluster = factor(Cluster, levels = Cluster))

# Calculate percentages
cluster_counts <- cluster_counts %>%
  mutate(Percentage = (Number_of_Words / sum(Number_of_Words)) * 100)

# Plot the distribution
ggplot(cluster_counts, aes(x = Cluster, y = Number_of_Words, fill = Cluster)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label = paste0(Number_of_Words, " (", round(Percentage, 1), "%)")),
            vjust = -0.5, size = 3) +
  scale_fill_brewer(palette = "Pastel1") +
  labs(title = "Distribution of Words Across Clusters",
       x = "Word Cluster",
       y = "Number of Words") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(hjust = 0.5))


#Top Words in Each Cluster--------------------------------------------------------------------

# Calculate word frequencies (sum across documents)
word_frequencies <- rowSums(data_text)

# Create a data frame with word frequencies and clusters
word_frequencies_df <- data.frame(
  Word = rownames(data_text),
  Frequency = word_frequencies,
  Cluster = factor(wordClusterIdx)
)

# Get top N words per cluster
top_n <- 20

top_words_per_cluster <- word_frequencies_df %>%
  group_by(Cluster) %>%
  top_n(top_n, wt = Frequency) %>%
  arrange(Cluster, Frequency)

# Order words by frequency within clusters
top_words_per_cluster <- top_words_per_cluster %>%
  group_by(Cluster) %>%
  arrange(desc(Frequency)) %>%
  ungroup()

# Plot top words per cluster
ggplot(top_words_per_cluster, aes(x = reorder_within(Word, Frequency, Cluster), y = Frequency, fill = Cluster)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  scale_x_reordered() +
  facet_wrap(~ Cluster, scales = "free_y") +
  labs(title = "Top Words per Word Cluster",
       x = NULL,
       y = "Frequency") +
  theme_minimal(base_size = 12) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 8),
        strip.background = element_rect(fill = "lightblue"),
        strip.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))


#Distribution of Documents Across Clusters--------------------------------------------------------------------

# Create a data frame for document clusters
doc_clusters_df <- data.frame(Document = colnames(data_text), Cluster = factor(docClusterIdx))

# Count documents per cluster
doc_cluster_counts <- doc_clusters_df %>%
  group_by(Cluster) %>%
  summarise(Number_of_Documents = n()) %>%
  arrange(desc(Number_of_Documents)) %>%
  mutate(Cluster = factor(Cluster, levels = Cluster))

# Plot the distribution
ggplot(doc_cluster_counts, aes(x = Cluster, y = Number_of_Documents, fill = Cluster)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label = Number_of_Documents), vjust = -0.5, size = 3) +
  scale_fill_brewer(palette = "Pastel2") +
  labs(title = "Distribution of Documents Across Clusters",
       x = "Document Cluster",
       y = "Number of Documents") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.title = element_text(hjust = 0.5))


#Documents Colored by Cluster--------------------------------------------------------------------

# Simplify document names if necessary
doc_clusters_df <- doc_clusters_df %>%
  mutate(Document_ID = 1:n())

# Add a count of 1 to each document (optional)
doc_clusters_df$Count <- 1

# Plot each document as a bar colored by cluster with document names
ggplot(doc_clusters_df, aes(x = reorder(Document, Cluster), y = Count, fill = Cluster)) +
  geom_bar(stat = "identity", width = 0.8, show.legend = TRUE) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Documents Colored by Cluster",
       x = "Document",
       y = NULL,
       fill = "Cluster") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.title.x = element_text(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  ) +
  scale_y_continuous(breaks = NULL)

#library(openxlsx)
#write.xlsx(doc_clusters_df, 'doc_clusters_df.xlsx')

#Top Words per Document Cluster--------------------------------------------------------------------

# Reshape data to long format
data_text_long <- data_text %>%
  rownames_to_column('Word') %>%
  gather(key = "Document", value = "Frequency", -Word)

# Add word cluster assignments
data_text_long <- data_text_long %>%
  left_join(word_clusters_df, by = "Word")

# Add document cluster assignments
data_text_long <- data_text_long %>%
  left_join(doc_clusters_df, by = "Document")

# Aggregate frequencies and ensure uniqueness
top_words_per_doc_cluster <- data_text_long %>%
  group_by(Cluster.y, Word) %>%  # Group by document cluster and word
  summarise(TotalFrequency = sum(Frequency),
            WordCluster = first(Cluster.x)) %>%  # Get the word's cluster
  ungroup() %>%
  group_by(Cluster.y) %>%
  top_n(30, wt = TotalFrequency) %>%
  arrange(Cluster.y, desc(TotalFrequency))

# Plot
ggplot(top_words_per_doc_cluster, aes(x = TotalFrequency, y = reorder_within(Word, TotalFrequency, Cluster.y), fill = WordCluster)) +
  geom_col(show.legend = T) +
  facet_wrap(~ Cluster.y, scales = "free_y", ncol = 1, strip.position = "left") +
  scale_y_reordered() +
  labs(x = "Total Frequency", y = "", title = "Top Words per Document Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 1),
        strip.background = element_rect(fill = "lightblue"),
        strip.placement = "outside")



#DATA FOR TIME LINE

# Reshape data_text to a long format to handle words per document
data_text_long <- data_text %>%
  rownames_to_column('Word') %>%
  gather(key = "Document", value = "Frequency", -Word)

# Combine the document cluster assignments
data_text_long <- data_text_long %>%
  left_join(doc_clusters_df, by = "Document")

data_text_long<-subset(data_text_long, select=-c(Document_ID, Count))
head(data_text_long)

library(dplyr)

data_text_long <- data_text_long %>%
  arrange(Document, desc(Frequency))

# Now, filter to get the top 5 words per Document
top_words_per_document <- data_text_long %>%
  group_by(Document) %>%
  slice_head(n = 5) %>%
  ungroup()

library(stringr)

top_words_per_document <- top_words_per_document %>%
  mutate(
    Year = str_extract(Document, "^\\d{4}"),  # Extracts the year
    Document = str_remove(Document, "^\\d{4}-")  # Removes the year and the following hyphen
  )

#library(openxlsx)
#write.xlsx(top_words_per_document, 'doc_clusters_df.xlsx')



































