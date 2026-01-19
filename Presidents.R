############################################################
# Presidents.R
# U.S. Presidential Inaugural Addresses:
# 1) Build a cleaned term-document matrix (TF-IDF) using quanteda
# 2) Export the matrix for MATLAB (SDKM/DKM) if needed
# 3) Import MATLAB cluster assignments and produce diagnostic plots
############################################################

## ----------------------------
## 0) Packages
## ----------------------------
suppressPackageStartupMessages({
  library(quanteda)
  library(lexicon)
  library(R.matlab)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(tidytext)   # for reorder_within / scale_x_reordered / scale_y_reordered
})

## ----------------------------
## 1) Working directory (robust)
## ----------------------------
if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(wd)
} else {
  message("rstudioapi not available: using current working directory = ", getwd())
}

## ----------------------------
## 2) Load corpus
## ----------------------------
# data_corpus_inaugural is provided by quanteda (as a corpus object)
data_corpus <- data_corpus_inaugural

## ----------------------------
## 3) Tokenization + preprocessing
## ----------------------------
data_tokens <- tokens(
  data_corpus,
  remove_punct = TRUE,
  remove_symbols = TRUE,
  remove_separators = TRUE,
  remove_numbers = TRUE
)

# Lowercase tokens (safer than trying to tolower() a corpus object)
data_tokens <- tokens_tolower(data_tokens)

# Lemmatize using lexicon::hash_lemmas
data_tokens <- tokens_replace(
  data_tokens,
  pattern = lexicon::hash_lemmas$token,
  replacement = lexicon::hash_lemmas$lemma
)

# Remove tokens with length < 3 or > 14
all_tokens <- unique(unlist(tokens(data_tokens)))
tokens_to_remove <- all_tokens[str_length(all_tokens) < 3 | str_length(all_tokens) > 14]
data_tokens <- tokens_remove(data_tokens, pattern = tokens_to_remove)

## ----------------------------
## 4) DFM construction and trimming
## ----------------------------
data_dfm <- dfm(data_tokens)

# Remove English stopwords
data_dfm_nostops <- dfm_remove(data_dfm, pattern = stopwords("en"))

# Inspect top features (optional)
topfeatures(data_dfm_nostops, 20)

# Trim rare terms (min term frequency threshold)
data_dfm_trim <- dfm_trim(data_dfm_nostops, min_termfreq = 11)
dim(data_dfm_trim)

# TF-IDF weighting (force = TRUE keeps zero-count features consistent)
data_dfm_tfidf <- dfm_tfidf(data_dfm_trim, force = TRUE)

## ----------------------------
## 5) Export matrix for MATLAB (optional)
## ----------------------------
# Convert to a plain data.frame:
# rows = documents, columns = terms (quanteda default)
data_text_doc_term <- convert(data_dfm_tfidf, to = "data.frame")

# Move docnames to rownames, then drop docname column
rownames(data_text_doc_term) <- data_text_doc_term[, 1]
data_text_doc_term <- data_text_doc_term[, -1, drop = FALSE]

# Transpose to match your MATLAB convention: rows = terms, cols = documents
data_text <- as.data.frame(t(data_text_doc_term))

# Uncomment if you want to write to CSV for MATLAB
# write.csv(data_text, "data_text.csv", row.names = TRUE)

## ----------------------------
## 6) Import cluster assignments from MATLAB
## ----------------------------
# Expected: wordClusterIdx length = nrow(data_text) (terms)
#           docClusterIdx  length = ncol(data_text) (documents)

word_mat <- readMat("wordClusterIdx.mat")
doc_mat  <- readMat("docClusterIdx.mat")

wordClusterIdx <- as.vector(unlist(word_mat$wordClusterIdx))
docClusterIdx  <- as.vector(unlist(doc_mat$docClusterIdx))

# Sanity checks
if (length(wordClusterIdx) != nrow(data_text)) {
  stop("Length mismatch: wordClusterIdx has length ", length(wordClusterIdx),
       " but nrow(data_text) = ", nrow(data_text))
}
if (length(docClusterIdx) != ncol(data_text)) {
  stop("Length mismatch: docClusterIdx has length ", length(docClusterIdx),
       " but ncol(data_text) = ", ncol(data_text))
}

## ----------------------------
## 7) Word cluster distribution
## ----------------------------
word_clusters_df <- data.frame(
  Word = rownames(data_text),
  WordCluster = factor(wordClusterIdx)
)

cluster_counts <- word_clusters_df %>%
  count(WordCluster, name = "Number_of_Words") %>%
  arrange(desc(Number_of_Words)) %>%
  mutate(
    WordCluster = factor(WordCluster, levels = WordCluster),
    Percentage = 100 * Number_of_Words / sum(Number_of_Words)
  )

ggplot(cluster_counts, aes(x = WordCluster, y = Number_of_Words, fill = WordCluster)) +
  geom_col(show.legend = FALSE) +
  geom_text(
    aes(label = paste0(Number_of_Words, " (", round(Percentage, 1), "%)")),
    vjust = -0.5, size = 3
  ) +
  scale_fill_brewer(palette = "Pastel1") +
  labs(
    title = "Distribution of Words Across Clusters",
    x = "Word Cluster",
    y = "Number of Words"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

## ----------------------------
## 8) Top words per word cluster (by total TF-IDF)
## ----------------------------
# Note: Since data_text is TF-IDF, "Frequency" here means aggregate TF-IDF weight.
word_weights <- rowSums(data_text)

word_weights_df <- data.frame(
  Word = rownames(data_text),
  Weight = word_weights,
  WordCluster = factor(wordClusterIdx)
)

top_n_words <- 20

top_words_per_cluster <- word_weights_df %>%
  group_by(WordCluster) %>%
  slice_max(order_by = Weight, n = top_n_words, with_ties = FALSE) %>%
  arrange(WordCluster, desc(Weight)) %>%
  ungroup()

ggplot(
  top_words_per_cluster,
  aes(x = reorder_within(Word, Weight, WordCluster), y = Weight, fill = WordCluster)
) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  scale_x_reordered() +
  facet_wrap(~ WordCluster, scales = "free_y") +
  labs(
    title = "Top Words per Word Cluster (aggregate TF-IDF weight)",
    x = NULL,
    y = "Aggregate TF-IDF"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 8),
    strip.background = element_rect(fill = "lightblue"),
    strip.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )

## ----------------------------
## 9) Document cluster distribution
## ----------------------------
doc_clusters_df <- data.frame(
  Document = colnames(data_text),
  DocCluster = factor(docClusterIdx)
)

doc_cluster_counts <- doc_clusters_df %>%
  count(DocCluster, name = "Number_of_Documents") %>%
  arrange(desc(Number_of_Documents)) %>%
  mutate(DocCluster = factor(DocCluster, levels = DocCluster))

ggplot(doc_cluster_counts, aes(x = DocCluster, y = Number_of_Documents, fill = DocCluster)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = Number_of_Documents), vjust = -0.5, size = 3) +
  scale_fill_brewer(palette = "Pastel2") +
  labs(
    title = "Distribution of Documents Across Clusters",
    x = "Document Cluster",
    y = "Number of Documents"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(hjust = 0.5)
  )

## ----------------------------
## 10) Documents colored by cluster
## ----------------------------
doc_clusters_df <- doc_clusters_df %>%
  mutate(Document_ID = row_number(), Count = 1)

ggplot(doc_clusters_df, aes(x = reorder(Document, DocCluster), y = Count, fill = DocCluster)) +
  geom_col(width = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Documents Colored by Cluster",
    x = "Document",
    y = NULL,
    fill = "Cluster"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  ) +
  scale_y_continuous(breaks = NULL)

## ----------------------------
## 11) Top words per document cluster
## ----------------------------
data_text_long <- data_text %>%
  rownames_to_column("Word") %>%
  pivot_longer(cols = -Word, names_to = "Document", values_to = "Weight")

data_text_long <- data_text_long %>%
  left_join(word_clusters_df, by = "Word") %>%
  left_join(doc_clusters_df,  by = "Document")

top_words_per_doc_cluster <- data_text_long %>%
  group_by(DocCluster, Word) %>%
  summarise(
    TotalWeight = sum(Weight),
    WordCluster = first(WordCluster),
    .groups = "drop"
  ) %>%
  group_by(DocCluster) %>%
  slice_max(order_by = TotalWeight, n = 30, with_ties = FALSE) %>%
  arrange(DocCluster, desc(TotalWeight)) %>%
  ungroup()

ggplot(
  top_words_per_doc_cluster,
  aes(x = TotalWeight, y = reorder_within(Word, TotalWeight, DocCluster), fill = WordCluster)
) +
  geom_col(show.legend = TRUE) +
  facet_wrap(~ DocCluster, scales = "free_y", ncol = 1, strip.position = "left") +
  scale_y_reordered() +
  labs(
    x = "Total TF-IDF (sum across documents in cluster)",
    y = "",
    title = "Top Words per Document Cluster"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 1),
    strip.background = element_rect(fill = "lightblue"),
    strip.placement = "outside"
  )

## ----------------------------
## 12) Timeline prep: top words per document
## ----------------------------
# Sort words by weight within each document and keep top 5 per document
top_words_per_document <- data_text_long %>%
  arrange(Document, desc(Weight)) %>%
  group_by(Document) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    Year = str_extract(Document, "^\\d{4}"),
    Document = str_remove(Document, "^\\d{4}-")
  )

# Optional export
# write.csv(top_words_per_document, "top_words_per_document.csv", row.names = FALSE)






