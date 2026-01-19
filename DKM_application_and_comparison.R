############################################################
# DKM_application_and_comparison.R
#
# Purpose:
# - Load TF-IDF matrix (terms x documents)
# - Load DKM outputs (Udkm, Vdkm) and derive hard labels
# - Load SDKM labels (wordClusterIdx, docClusterIdx)
# - Produce comparable plots: top words per cluster + docs-by-cluster bars
# - Align DKM cluster labels to SDKM using assignment (Hungarian/LSAP)
# - Export a labeled word list and contingency tables
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(R.matlab)
  library(stringr)
  library(tidytext)  # reorder_within(), scale_x_reordered()
  library(clue)      # solve_LSAP (Hungarian)
})

## ----------------------------
## 0) Working directory (robust)
## ----------------------------
if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(wd)
} else {
  message("rstudioapi not available: using current working directory = ", getwd())
}

## ----------------------------
## 1) Load data_text (terms x docs)
## ----------------------------
load("data_text.Rdata")  # must contain object `data_text`
X <- as.matrix(data_text)
storage.mode(X) <- "double"

terms <- rownames(X)
docs  <- colnames(X)

if (is.null(terms) || is.null(docs)) {
  stop("data_text must have rownames (terms) and colnames (documents).")
}

cat("Loaded X:", nrow(X), "terms x", ncol(X), "docs\n")

# Word frequencies (sum across documents)
word_frequencies <- rowSums(X)

## ----------------------------
## 2) Helpers
## ----------------------------

# robust extraction from readMat()
get_mat_var <- function(mat_list, candidates) {
  nm <- names(mat_list)
  hit <- candidates[candidates %in% nm]
  if (length(hit) > 0) return(mat_list[[hit[1]]])
  # fallback: first element (older files sometimes store unnamed)
  if (length(mat_list) >= 1) return(mat_list[[1]])
  NULL
}

# convert membership matrix -> hard labels
mat_to_labels <- function(M) {
  if (!is.matrix(M)) M <- as.matrix(M)
  max.col(M, ties.method = "first")
}

# coerce possible factor labels to integer safely
as_int_label <- function(x) {
  if (is.factor(x)) as.integer(as.character(x)) else as.integer(x)
}

## ----------------------------
## 3) Load SDKM labels (wordClusterIdx/docClusterIdx)
## ----------------------------
wc_mat <- readMat("wordClusterIdx.mat")
dc_mat <- readMat("docClusterIdx.mat")

wordClusterIdx <- as.vector(get_mat_var(wc_mat, c("wordClusterIdx", "wordCluster", "wc")))
docClusterIdx  <- as.vector(get_mat_var(dc_mat, c("docClusterIdx",  "docCluster",  "dc")))

if (length(wordClusterIdx) != nrow(X)) {
  stop("wordClusterIdx length does not match #terms.")
}
if (length(docClusterIdx) != ncol(X)) {
  stop("docClusterIdx length does not match #docs.")
}

cat("SDKM: K =", length(unique(wordClusterIdx)),
    "| Q =", length(unique(docClusterIdx)), "\n")

## ----------------------------
## 4) 3D “SDKM document coordinates” (optional object)
## ----------------------------
# For each document, compute mean TF-IDF within each word-cluster
K_sdkm <- length(unique(wordClusterIdx))

# Precompute term indices per word cluster (avoid doing inside the doc loop)
terms_by_cluster <- lapply(sort(unique(wordClusterIdx)), function(k) which(wordClusterIdx == k))

doc_coords <- sapply(seq_len(ncol(X)), function(i) {
  v <- X[, i]
  vapply(terms_by_cluster, function(idx) mean(v[idx]), numeric(1))
})
doc_coords <- t(doc_coords)  # docs x K

doc_3d <- data.frame(
  Document   = docs,
  docCluster = factor(docClusterIdx)
)

# If K==3, keep x/y/z names; otherwise keep generic columns
if (ncol(doc_coords) == 3) {
  doc_3d$x <- doc_coords[,1]
  doc_3d$y <- doc_coords[,2]
  doc_3d$z <- doc_coords[,3]
} else {
  doc_3d <- bind_cols(doc_3d, as.data.frame(doc_coords))
}

## ----------------------------
## 5) Load DKM outputs (Udkm, Vdkm)
## ----------------------------
Udkm_list <- readMat("Udkm.mat")
Vdkm_list <- readMat("Vdkm.mat")

Udkm <- get_mat_var(Udkm_list, c("Udkm", "U", "UdkmO"))
Vdkm <- get_mat_var(Vdkm_list, c("Vdkm", "V", "VdkmO"))

Udkm <- as.matrix(Udkm)
Vdkm <- as.matrix(Vdkm)

# Sanity checks: Udkm should be (#terms x K_dkm), Vdkm should be (#docs x Q_dkm)
if (nrow(Udkm) != nrow(X)) {
  stop("Udkm rows must match #terms.")
}
if (nrow(Vdkm) != ncol(X)) {
  stop("Vdkm rows must match #docs.")
}

K_dkm <- ncol(Udkm)
Q_dkm <- ncol(Vdkm)

cat("DKM: K =", K_dkm, "| Q =", Q_dkm, "\n")

wordClusterIdx_DKM <- mat_to_labels(Udkm)  # length = #terms
docClusterIdx_DKM  <- mat_to_labels(Vdkm)  # length = #docs

## ----------------------------
## 6) PART 1 — DKM plots
## ----------------------------
top_n <- 20

word_freq_df_DKM <- data.frame(
  Word      = terms,
  Frequency = word_frequencies,
  Cluster   = factor(wordClusterIdx_DKM)
)

top_words_per_cluster_DKM <- word_freq_df_DKM %>%
  group_by(Cluster) %>%
  slice_max(order_by = Frequency, n = top_n, with_ties = FALSE) %>%
  arrange(Cluster, desc(Frequency)) %>%
  ungroup()

p1 <- ggplot(top_words_per_cluster_DKM,
             aes(x = reorder_within(Word, Frequency, Cluster),
                 y = Frequency, fill = Cluster)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  scale_x_reordered() +
  facet_wrap(~ Cluster, scales = "free_y") +
  labs(title = "Top 20 Words per Cluster (DKM)",
       x = NULL, y = "Frequency") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 8),
        strip.background = element_rect(fill = "lightblue"),
        strip.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
print(p1)

doc_clusters_df_DKM <- data.frame(
  Document = docs,
  Cluster  = factor(docClusterIdx_DKM),
  Count    = 1
)

p2 <- ggplot(doc_clusters_df_DKM,
             aes(x = reorder(Document, Cluster), y = Count, fill = Cluster)) +
  geom_bar(stat = "identity", width = 0.8, show.legend = TRUE) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Documents by Cluster (DKM)",
       x = "Document", y = NULL, fill = "Cluster") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_y_continuous(breaks = NULL)
print(p2)

## ----------------------------
## 7) PART 2 — SDKM plots
## ----------------------------
word_freq_df_SDKM <- data.frame(
  Word      = terms,
  Frequency = word_frequencies,
  Cluster   = factor(wordClusterIdx)
)

top_words_per_cluster_SDKM <- word_freq_df_SDKM %>%
  group_by(Cluster) %>%
  slice_max(order_by = Frequency, n = top_n, with_ties = FALSE) %>%
  arrange(Cluster, desc(Frequency)) %>%
  ungroup()

p3 <- ggplot(top_words_per_cluster_SDKM,
             aes(x = reorder_within(Word, Frequency, Cluster),
                 y = Frequency, fill = Cluster)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  scale_x_reordered() +
  facet_wrap(~ Cluster, scales = "free_y") +
  labs(title = "Top 20 Words per Cluster (SDKM)",
       x = NULL, y = "Frequency") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 8),
        strip.background = element_rect(fill = "lightblue"),
        strip.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
print(p3)

doc_clusters_df_SDKM <- data.frame(
  Document = docs,
  Cluster  = factor(docClusterIdx),
  Count    = 1
)

p4 <- ggplot(doc_clusters_df_SDKM,
             aes(x = reorder(Document, Cluster), y = Count, fill = Cluster)) +
  geom_bar(stat = "identity", width = 0.8, show.legend = TRUE) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Documents by Cluster (SDKM)",
       x = "Document", y = NULL, fill = "Cluster") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_y_continuous(breaks = NULL)
print(p4)

## ----------------------------
## 8) PART 3 — Align DKM labels to SDKM
## ----------------------------

## (A) Terms: use Jaccard overlap of top-words lists and Hungarian assignment
dkm_clusters_terms  <- sort(unique(top_words_per_cluster_DKM$Cluster))
sdkm_clusters_terms <- sort(unique(top_words_per_cluster_SDKM$Cluster))

top_words_list_DKM <- lapply(dkm_clusters_terms, function(cl) {
  unique(top_words_per_cluster_DKM %>% filter(Cluster == cl) %>% pull(Word))
})
names(top_words_list_DKM) <- as.character(dkm_clusters_terms)

top_words_list_SDKM <- lapply(sdkm_clusters_terms, function(cl) {
  unique(top_words_per_cluster_SDKM %>% filter(Cluster == cl) %>% pull(Word))
})
names(top_words_list_SDKM) <- as.character(sdkm_clusters_terms)

sim_mat <- matrix(0, nrow = length(dkm_clusters_terms), ncol = length(sdkm_clusters_terms),
                  dimnames = list(as.character(dkm_clusters_terms), as.character(sdkm_clusters_terms)))

for (i in seq_along(dkm_clusters_terms)) {
  for (j in seq_along(sdkm_clusters_terms)) {
    a <- top_words_list_DKM[[i]]
    b <- top_words_list_SDKM[[j]]
    inter <- length(intersect(a, b))
    uni   <- length(union(a, b))
    sim_mat[i, j] <- if (uni > 0) inter / uni else 0
  }
}

cost_mat <- 1 - sim_mat
assignment_terms <- solve_LSAP(cost_mat)  # returns column indices
mapping_terms <- sdkm_clusters_terms[as.vector(assignment_terms)]
names(mapping_terms) <- as.character(dkm_clusters_terms)

cat("\nMapping (terms) DKM -> SDKM:\n")
print(data.frame(DKM_cluster = dkm_clusters_terms,
                 SDKM_cluster = mapping_terms))

# Add aligned labels for plotting/comparisons
top_words_per_cluster_DKM$Cluster_reordered <- factor(
  mapping_terms[as.character(top_words_per_cluster_DKM$Cluster)],
  levels = as.character(sdkm_clusters_terms)
)
top_words_per_cluster_SDKM$Cluster_reordered <- factor(
  as.character(top_words_per_cluster_SDKM$Cluster),
  levels = as.character(sdkm_clusters_terms)
)

## (B) Documents: use contingency table + Hungarian assignment (more stable than majority)
doc_cont <- table(DKM = doc_clusters_df_DKM$Cluster, SDKM = doc_clusters_df_SDKM$Cluster)
# Convert to cost: maximize overlap => minimize negative overlap
cost_docs <- max(doc_cont) - doc_cont
assignment_docs <- solve_LSAP(cost_docs)

mapping_docs <- colnames(doc_cont)[as.vector(assignment_docs)]
names(mapping_docs) <- rownames(doc_cont)

cat("\nMapping (docs) DKM -> SDKM (LSAP on contingency):\n")
print(mapping_docs)

doc_clusters_df_DKM$Cluster_reordered <- factor(
  mapping_docs[as.character(doc_clusters_df_DKM$Cluster)],
  levels = colnames(doc_cont)
)
doc_clusters_df_SDKM$Cluster_reordered <- factor(
  as.character(doc_clusters_df_SDKM$Cluster),
  levels = colnames(doc_cont)
)

# Combine for comparison plot
doc_clusters_df_DKM$Method  <- "DKM"
doc_clusters_df_SDKM$Method <- "SDKM"
doc_clusters_combined <- bind_rows(doc_clusters_df_DKM, doc_clusters_df_SDKM)

p6 <- ggplot(doc_clusters_combined,
             aes(x = reorder(Document, Cluster_reordered),
                 y = Count, fill = Cluster_reordered)) +
  geom_bar(stat = "identity", width = 0.8, show.legend = TRUE) +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~ Method, ncol = 1) +
  labs(title = "Documents by Cluster (Aligned Labels)",
       x = "Document", y = NULL, fill = "Cluster") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_y_continuous(breaks = NULL)
print(p6)

## ----------------------------
## 9) Export: all words labeled (DKM + SDKM)
## ----------------------------
dkm_output <- word_freq_df_DKM %>%
  transmute(Word, Cluster = as.character(Cluster), Method = "DKM")

sdkm_output <- word_freq_df_SDKM %>%
  transmute(Word, Cluster = as.character(Cluster), Method = "SDKM")

all_words <- bind_rows(dkm_output, sdkm_output)

write.table(all_words,
            file = "All_words.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

cat("\nWrote All_words.txt with word labels for DKM and SDKM.\n")

## ----------------------------
## 10) Contingency table on ALL words (optional manual swap)
## ----------------------------
# If you already KNOW DKM labels 1<->3 should be swapped, you can do it here.
swap_dkm_label <- function(x) {
  x <- as_int_label(x)
  ifelse(x == 1, 3,
         ifelse(x == 3, 1, x))
}

swapped_dkm <- swap_dkm_label(wordClusterIdx_DKM)

cont_table_swapped <- table(DKM = swapped_dkm, SDKM = wordClusterIdx)
cat("\nContingency (DKM vs SDKM) with DKM swap 1<->3:\n")
print(cont_table_swapped)

## ----------------------------
## 11) Plot overlap of top words (DKM_i vs SDKM_i) after swap 1<->3
## ----------------------------
# This reproduces your “single bar per word” idea, coloring by overlap status.
top_n_overlap <- 30

swapped_top_dkm <- top_words_per_cluster_DKM %>%
  mutate(Cluster_num = swap_dkm_label(Cluster)) %>%
  mutate(Cluster_num = as_int_label(Cluster_num))

plot_data <- tibble()

for (i in sort(unique(wordClusterIdx))) {
  dkm_sub <- swapped_top_dkm %>%
    filter(Cluster_num == i) %>%
    arrange(desc(Frequency)) %>%
    slice_head(n = top_n_overlap) %>%
    transmute(Word, freq_dkm = Frequency)
  
  sdkm_sub <- top_words_per_cluster_SDKM %>%
    filter(as_int_label(Cluster) == i) %>%
    arrange(desc(Frequency)) %>%
    slice_head(n = top_n_overlap) %>%
    transmute(Word, freq_sdkm = Frequency)
  
  pair_df <- full_join(dkm_sub, sdkm_sub, by = "Word") %>%
    mutate(freq_dkm  = if_else(is.na(freq_dkm),  0, freq_dkm),
           freq_sdkm = if_else(is.na(freq_sdkm), 0, freq_sdkm),
           Status = case_when(
             freq_dkm > 0 & freq_sdkm > 0 ~ "Overlap",
             freq_dkm > 0                 ~ "DKM only",
             TRUE                         ~ "SDKM only"
           ),
           freqPlot  = freq_dkm + freq_sdkm,
           ClusterID = factor(i, levels = sort(unique(wordClusterIdx)))
    )
  
  plot_data <- bind_rows(plot_data, pair_df)
}

p_overlap <- ggplot(plot_data,
                    aes(x = reorder_within(Word, freqPlot, ClusterID),
                        y = freqPlot, fill = Status)) +
  geom_col(show.legend = TRUE) +
  scale_x_reordered() +
  coord_flip() +
  facet_wrap(~ ClusterID, scales = "free_y") +
  labs(title = "Top words overlap per cluster (DKM swapped 1↔3 vs SDKM)",
       x = NULL, y = "Frequency", fill = "Word Source") +
  theme_minimal(base_size = 12) +
  theme(strip.background = element_rect(fill = "lightblue", color = NA),
        strip.text = element_text(size = 12))
print(p_overlap)







