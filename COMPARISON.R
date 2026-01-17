library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(quanteda)
library(R.matlab)
library(stringr)
library(tidytext)
library(ggrepel)
library(clue)

# --- Set working directory (adjust as needed) ---
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

# --- Load the data_text matrix ---
load("data_text.Rdata") 
terms <- rownames(data_text)

############################################################################
############################################################################

# Load SDKM cluster assignments from MATLAB
wordClusterData <- readMat("wordClusterIdx.mat")
docClusterData  <- readMat("docClusterIdx.mat")
wordClusterIdx  <- as.vector(unlist(wordClusterData$wordClusterIdx))
docClusterIdx   <- as.vector(unlist(docClusterData$docClusterIdx))

table(wordClusterIdx)
table(docClusterIdx)

# Load DKM outputs from MATLAB
wordDKM <- readMat("Udkm.mat")[[1]]
varDKM  <- readMat("Vdkm.mat")[[1]]
wordDKM <- data.frame(wordDKM)
varDKM  <- data.frame(varDKM)

colSums(wordDKM[,c(1,2,3)])
colSums(varDKM)


############################################################################
#      CONTINGENCY TABLE OF ALL WORDS: DKM CLUSTERS VS. SDKM CLUSTERS
############################################################################

wordClusterIdx_DKM <- apply(wordDKM, 1, function(x) which(x == 1))
table(wordClusterIdx_DKM)

cont_table <- table(
  DKM  = wordClusterIdx_DKM,
  SDKM = wordClusterIdx
)

print(cont_table)


############################################################################
# PLOT TOP 10 WORDS (DKM_i vs. SDKM_i) IN ONE BAR, COLOR IF OVERLAP
############################################################################

library(dplyr)
library(ggplot2)
library(tidytext)

# ----- DKM: Top Words per Cluster (limit to 30 per cluster) -----
word_frequencies <- rowSums(data_text)
word_frequencies_df_DKM <- data.frame(
  Word = rownames(data_text),
  Frequency = word_frequencies,
  Cluster = factor(wordClusterIdx_DKM)
)

top_n <- 30 

top_words_per_cluster_DKM <- word_frequencies_df_DKM %>%
  group_by(Cluster) %>%
  top_n(top_n, wt = Frequency) %>%
  arrange(Cluster, desc(Frequency)) %>%
  ungroup()

# ----- SDKM: Top Words per Cluster -----
word_frequencies_df_SDKM <- data.frame(
  Word = rownames(data_text),
  Frequency = word_frequencies,
  Cluster = factor(wordClusterIdx)
)

top_words_per_cluster_SDKM <- word_frequencies_df_SDKM %>%
  group_by(Cluster) %>%
  top_n(top_n, wt = Frequency) %>%
  arrange(Cluster, desc(Frequency)) %>%
  ungroup()

plot_data <- data.frame()

for(i in 1:3) {
  # 1) Subset top 10 for DKM cluster i (AFTER swapped)
  dkm_sub <- top_words_per_cluster_DKM %>%
    filter(Cluster == i) %>%
    arrange(desc(Frequency)) %>%
    slice_head(n = top_n) %>%
    rename(freq_dkm = Frequency) %>%
    select(Word, freq_dkm)
  
  # 2) Subset top 10 for SDKM cluster i (unchanged)
  sdkm_sub <- top_words_per_cluster_SDKM %>%
    filter(Cluster == i) %>%
    arrange(desc(Frequency)) %>%
    slice_head(n = top_n) %>%
    rename(freq_sdkm = Frequency) %>%
    select(Word, freq_sdkm)
  
  # 3) Full join => each word appears once
  pair_df <- full_join(dkm_sub, sdkm_sub, by="Word") %>%
    mutate(
      freq_dkm  = if_else(is.na(freq_dkm),  0, freq_dkm),
      freq_sdkm = if_else(is.na(freq_sdkm), 0, freq_sdkm)
    )
  
  # 4) Mark status
  pair_df <- pair_df %>%
    mutate(
      Status = case_when(
        freq_dkm > 0 & freq_sdkm > 0 ~ "Overlap",
        freq_dkm > 0                ~ "DKM only",
        freq_sdkm > 0               ~ "SDKM only"
      )
    )
  
  # 5) Single bar height = freq_dkm + freq_sdkm
  pair_df <- pair_df %>%
    mutate(freqPlot = freq_dkm + freq_sdkm)
  
  # 6) Facet panel label
  pair_df$ClusterID <- i
  
  # 7) Append
  plot_data <- bind_rows(plot_data, pair_df)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~
# B) Plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Factor the cluster IDs in ascending order
plot_data$ClusterID <- factor(plot_data$ClusterID, levels=c(1,2,3))

p <- ggplot(plot_data, 
            aes(x = reorder_within(Word, freqPlot, ClusterID), 
                y = freqPlot,
                fill = Status)) +
  geom_col(show.legend = TRUE) +
  scale_x_reordered() +
  coord_flip() +
  facet_wrap(~ ClusterID, scales="free_y") +
  scale_fill_manual(values=c("DKM only"="#1B9E77",
                             "SDKM only"="#D95F02", 
                             "Overlap"="#8E44AD")) +
  labs(
    title = "",
    x = NULL, y = "Frequency",
    fill = "Word Source"
  ) +
  theme_minimal(base_size=12) +
  theme(
    strip.background = element_rect(fill="lightblue", color=NA),
    strip.text = element_text(size=12)
  )

print(p)



############################################################################
#      CONTINGENCY TABLE OF ALL WORDS: DKM CLUSTERS VS. SDKM CLUSTERS
#      WITH DKM 1 ↔ 3 LABELS SWAPPED
############################################################################

# We define a small helper function to swap cluster labels: 1->3, 3->1, 2->2
swap_dkm_label <- function(x) {
  ifelse(x == 1, 3,
         ifelse(x == 3, 1, x)
  )
}

# 1) Suppose you have:
#       wordClusterIdx_DKM : integer vector of length (#words), each in {1,2,3}
#       wordClusterIdx     : integer vector of length (#words), each in {1,2,3} for SDKM
#    We will recode DKM's cluster labels using the swap function:
swapped_dkm <- swap_dkm_label(wordClusterIdx_DKM)

# 2) Build the 3×3 table (DKM vs. SDKM) AFTER swapping
cont_table_swap <- table(
  DKM  = swapped_dkm,
  SDKM = wordClusterIdx
)

print(cont_table_swap)


############################################################################
# PLOT TOP 10 WORDS (DKM_i vs. SDKM_i) IN ONE BAR, COLOR IF OVERLAP
#          with DKM 1 ↔ 3 labels swapped
############################################################################

library(dplyr)
library(ggplot2)
library(tidytext)  # for reorder_within, scale_x_reordered()

# Suppose you have data frames:
#   top_words_per_cluster_DKM  with columns (Word, Frequency, Cluster)  # original
#   top_words_per_cluster_SDKM with columns (Word, Frequency, Cluster)
# each with cluster ∈ {1,2,3}, but we want to swap 1↔3 in the DKM data.

# ~~~~~~~~~~~~~~~~~~~~~~~~~
#  A) Create swapped DKM data
# ~~~~~~~~~~~~~~~~~~~~~~~~~
swapped_dkm_df <- top_words_per_cluster_DKM %>%
  mutate(Cluster = swap_dkm_label(Cluster))  # now old 1→3, old 3→1, 2→2

# If your `Cluster` column is a factor, you might want to convert it to numeric
# first, or you can do a factor recoding. For simplicity, we assume it's numeric.

top_n <- 30

plot_data_swap <- data.frame()

for(i in 1:3) {
  # 1) Subset top 10 for DKM cluster i (AFTER swapped)
  dkm_sub <- swapped_dkm_df %>%
    filter(Cluster == i) %>%
    arrange(desc(Frequency)) %>%
    slice_head(n = top_n) %>%
    rename(freq_dkm = Frequency) %>%
    select(Word, freq_dkm)
  
  # 2) Subset top 10 for SDKM cluster i (unchanged)
  sdkm_sub <- top_words_per_cluster_SDKM %>%
    filter(Cluster == i) %>%
    arrange(desc(Frequency)) %>%
    slice_head(n = top_n) %>%
    rename(freq_sdkm = Frequency) %>%
    select(Word, freq_sdkm)
  
  # 3) Full join => each word appears once
  pair_df <- full_join(dkm_sub, sdkm_sub, by="Word") %>%
    mutate(
      freq_dkm  = if_else(is.na(freq_dkm),  0, freq_dkm),
      freq_sdkm = if_else(is.na(freq_sdkm), 0, freq_sdkm)
    )
  
  # 4) Mark status
  pair_df <- pair_df %>%
    mutate(
      Status = case_when(
        freq_dkm > 0 & freq_sdkm > 0 ~ "Overlap",
        freq_dkm > 0                ~ "DKM only",
        freq_sdkm > 0               ~ "SDKM only"
      )
    )
  
  # 5) Single bar height = freq_dkm + freq_sdkm
  pair_df <- pair_df %>%
    mutate(freqPlot = freq_dkm + freq_sdkm)
  
  # 6) Facet panel label
  pair_df$ClusterID <- i
  
  # 7) Append
  plot_data_swap <- bind_rows(plot_data_swap, pair_df)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~
# B) Plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Factor the cluster IDs in ascending order
plot_data_swap$ClusterID <- factor(plot_data_swap$ClusterID, levels=c(1,2,3))

s <- ggplot(plot_data_swap, 
            aes(x = reorder_within(Word, freqPlot, ClusterID), 
                y = freqPlot,
                fill = Status)) +
  geom_col(show.legend = TRUE) +
  scale_x_reordered() +
  coord_flip() +
  facet_wrap(~ ClusterID, scales="free_y") +
  scale_fill_manual(values=c("DKM only"="#1B9E77",
                             "SDKM only"="#D95F02", 
                             "Overlap"="#8E44AD")) +
  labs(
    title = "",
    x = NULL, y = "Frequency",
    fill = "Word Source"
  ) +
  theme_minimal(base_size=12) +
  theme(
    strip.background = element_rect(fill="lightblue", color=NA),
    strip.text = element_text(size=12)
  )

print(s)





