library(readtext)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(quanteda)
library(R.matlab)
library(stringr)
library(tidytext)
library(ggrepel)
library(lexicon)
library(Matrix)
library(quanteda.textstats)
library(stringi)


wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

# If your working directory is the SDKM folder, this relative path will work:
rt <- readtext("20news-18828/*/*", encoding = "UTF-8")

# Add the newsgroup label from the directory name as a document variable
#rt$group <- basename(dirname(rt$doc_id))

quanteda_options(threads = max(1, parallel::detectCores(logical = TRUE) - 1))

# ---- 1) Header & boilerplate cleanup (keep Subject content, drop "From:")
clean_20ng <- function(x) {
  # keep Subject content but remove the label
  x <- gsub("(?mi)^Subject:\\s*", "", x, perl = TRUE)
  # drop the whole From line
  x <- gsub("(?mi)^From:.*$", " ", x, perl = TRUE)
  
  # drop a bunch of headers entirely
  hdrs <- c(
    "Archive-name","Last-modified","Version","Lines","Organization",
    "NNTP-Posting-Host","Distribution","Reply-To","Keywords","Xref",
    "Path","Message-ID","References","Sender","In-Reply-To","Newsgroups",
    "Date","Content-Type","Content-Transfer-Encoding","Mime-Version"
  )
  x <- gsub(paste0("(?mi)^(", paste(hdrs, collapse = "|"), "):.*$"),
            " ", x, perl = TRUE)
  
  # quoted replies, signatures, PGP blocks
  x <- gsub("(?m)^>+.*$", " ", x, perl = TRUE)
  x <- gsub("(?m)^-- ?.*$", " ", x, perl = TRUE)
  x <- gsub("(?s)-----BEGIN PGP.*?END PGP-----", " ", x, perl = TRUE)
  
  # emails, urls, file paths, ids/numbers
  x <- gsub("<[^>\\s]+@[^>\\s]+>", " ", x, perl = TRUE)
  x <- gsub("[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\\.[A-Za-z]{2,}", " ", x, perl = TRUE)
  x <- gsub("https?://\\S+|www\\.[^\\s]+", " ", x, perl = TRUE)
  x <- gsub("(?:\\b[A-Za-z]:)?[\\\\/][\\w./-]+", " ", x, perl = TRUE)
  x <- gsub("\\b\\d+[A-Za-z0-9_-]*\\b", " ", x, perl = TRUE)
  
  x
}


rt$text <- vapply(rt$text, clean_20ng, "", USE.NAMES = FALSE)

# ---- 2) Corpus -> tokens (email-aware)
corp <- corpus(rt, text_field = "text")

toks <- tokens(
  corp, what = "word",
  remove_punct = TRUE, remove_symbols = TRUE,
  remove_numbers = TRUE, remove_separators = TRUE,
  split_hyphens = TRUE
)
toks <- tokens_tolower(toks)

# keep only alphabetic words and constrain length 3–14
toks <- tokens_keep(toks, pattern = "^[a-z]{3,14}$", valuetype = "regex", padding = FALSE)

# lemmatize (fast dictionary replace)
toks <- tokens_replace(
  toks,
  pattern     = lexicon::hash_lemmas$token,
  replacement = lexicon::hash_lemmas$lemma
)

# standard + domain stopwords (metadata words from emails)
meta_stop <- c(
  "subject","from","organization","lines","article","writes","wrote",
  "posting","host","nntp","path","message","id","references","sender",
  "reply","keywords","distribution","newsgroups","xref","gmt","re",
  "edu","com","uk","ca","apr","may","jun","jul","aug","sep","oct","nov","dec"
)
toks <- tokens_remove(toks, c(stopwords("en"), meta_stop, c("don","ll","ve","re","t","d","s")))

# ---- 3) DFM -> trim -> TF-IDF -> (optional) L2 normalize
D <- dfm(toks)
dim(D)
topfeatures(D, 20)

# Remove features that are only numbers (e.g., "2", "3", "1993", "3.1")
D <- dfm_remove(D, pattern = "^\\d+(?:\\.\\d+)?$", valuetype = "regex", case_insensitive = FALSE)

# keep terms in ≥1% docs, and in ≤50% of docs (drop ultra-common)
D <- dfm_trim(D, min_docfreq = 0.01, docfreq_type = "prop")
D <- dfm_trim(D, max_docfreq = 0.5, docfreq_type = "prop")

topfeatures(D, 20)
summary(textstat_summary(D))

# features that are PURE numbers (e.g., "486", "1993", "3.1")
num_feats   <- featnames(D)[stri_detect_regex(featnames(D), "^(\\p{N}|\\d|\\d[\\d.]+)$")]
length(num_feats); head(num_feats)

# features that CONTAIN any digit (e.g., "x11", "3com", "win3")
alnum_feats <- featnames(D)[stri_detect_regex(featnames(D), "\\d")]
length(alnum_feats); head(alnum_feats)

# ---- Remove empty / very short / very long documents
nz <- ntoken(D) > 0
cat("Empty docs after trimming:", sum(!nz), "\n")
D <- D[nz, ]

cat("Docs left:", ndoc(D), " | Terms left:", nfeat(D), "\n")

# tf-idf + L2 normalization (good for clustering)
D_tfidf <- dfm_tfidf(D, scheme_tf = "prop", scheme_df = "inverse", force = TRUE)

dn <- docnames(D_tfidf)
parts <- strsplit(dn, "/", fixed = TRUE)
macro_clean <- gsub("[^A-Za-z0-9]+", "_", vapply(parts, `[[`, "", 1))
doc_id_num  <- vapply(parts, function(p) p[[length(p)]], "")
docnames(D_tfidf) <- make.unique(paste0(macro_clean, "_", doc_id_num), sep = "_d")

# ---- 4) Export formats for SDKM

data_news <- convert(D_tfidf, to = "data.frame")
rownames(data_news)<-data_news[,1]
data_news<-data_news[,-1]
data_news<-as.data.frame(t(data_news))

## From your current objects
doc_ids <- colnames(data_news)                 # same as docnames(D_tfidf) after your rename
# strip the trailing "_<number>" and any uniqueness suffix like "_d1"
macro <- sub("_[0-9]+(?:_d\\d+)?$", "", doc_ids)

# How many & which topics?
length(unique(macro))
sort(unique(macro))
table((macro))

#write.csv(data_news, "data_news.csv")
#save(data_news, file="data_news.Rdata")


#-------------------------------------------------------------------------------
#SAMPLING-----------------------------------------------------------------------
#-------------------------------------------------------------------------------

load("data_news.Rdata")

set.seed(123)

p <- 0.15                      # 20% sample
min_per_topic <- 100           # ensure at least this many per topic

idx_by_topic <- split(seq_along(doc_ids), macro)

sampled_idx <- unlist(lapply(idx_by_topic, function(ix) {
  n <- length(ix)
  size <- ceiling(n * p)
  size <- max(size, min_per_topic)     # honor the floor
  size <- min(size, n)                 # can't exceed available
  sample(ix, size, replace = FALSE)
}), use.names = FALSE)

# subset columns (docs) and keep order stable
sampled_idx <- sort(sampled_idx)
data_news_samp <- data_news[, sampled_idx, drop = FALSE]

# check distribution fidelity
macro_samp <- macro[sampled_idx]
orig_prop  <- prop.table(table(macro))
samp_prop  <- prop.table(table(macro_samp))
print(round(orig_prop, 3))
print(round(samp_prop, 3))

# write MATLAB-ready CSV (first column = terms)
out <- cbind(term = rownames(data_news_samp), data_news_samp)
#write.csv(out, "data_news_sample.csv", row.names = FALSE)

# (optional) also save R object
#save(data_news_samp, file="data_news_samp.Rdata")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#20NEWS SDKM--------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Load cluster assignments from MATLAB

U<-readMat("Usdkm_news.mat")
Y<-readMat("Ymsdkm_news.mat")
V<-readMat("Vsdkm_news.mat")
f<-readMat("fsdkm_news.mat")

Umat <- U$Usdkm.news   # terms x 3
Vmat <- V$Vsdkm.news   # docs  x 2

# cluster labels
wordClusterIdx <- max.col(Umat, ties.method = "first")  # length = #terms
docClusterIdx  <- max.col(Vmat, ties.method = "first")  # length = #docs

# quick counts
table(wordClusterIdx)
table(docClusterIdx)

# Create a data frame for word clusters
word_clusters_df <- data.frame(Word = rownames(data_news_samp), Cluster = factor(wordClusterIdx))
# Create a data frame for document clusters
doc_clusters_df <- data.frame(Document = colnames(data_news_samp), Cluster = factor(docClusterIdx))

library(readr)
pseudoF_20news <- read_csv("pseudoF_20news.csv")
library(tidyverse)
pseudoF_20news<-pseudoF_20news %>% remove_rownames %>% column_to_rownames(var="Row")

pseudoF_20news <- round(pseudoF_20news, 2)

library(dplyr); library(tidyr); library(ggplot2)
library(aricode)      # NMI
library(proxy)        # cosine distances
library(cluster)      # silhouette

X  <- as.matrix(data_news_samp)                  # terms x docs (same order as MATLAB CSV)
terms <- rownames(X); docs <- colnames(X)

K <- ncol(Umat); Q <- ncol(Vmat)

# labels (1..K for words, 1..Q for docs)
wordClusterIdx <- max.col(Umat, ties.method = "first")
docClusterIdx  <- max.col(Vmat, ties.method = "first")

# Quick sanity
stopifnot(nrow(Umat)==nrow(X), nrow(Vmat)==ncol(X))
table_words <- table(wordClusterIdx); table_docs <- table(docClusterIdx)
print(table_words); print(table_docs)

# --------- 1) What each WORD-CLUSTER is about ----------
#   Top-N terms per word cluster by average TF-IDF across docs
top_terms_word_cluster <- function(X, terms, wordClusterIdx, top_n = 30){
  res <- lapply(1:max(wordClusterIdx), function(k){
    idx <- which(wordClusterIdx==k)
    sc  <- rowMeans(X[idx, , drop=FALSE])      # mean tf-idf over all docs
    ord <- order(sc, decreasing=TRUE)[1:min(top_n, length(sc))]
    tibble(cluster=k, rank=1:length(ord), term=terms[idx][ord], mean_tfidf=sc[ord])
  })
  bind_rows(res)
}
top_words <- top_terms_word_cluster(X, terms, wordClusterIdx, top_n = 50)

# Also: distinctiveness of each term for doc clusters (log-odds, smoothed)
top_terms_doc_side <- function(X, terms, docClusterIdx, top_n=50, alpha=0.5){
  Tn <- nrow(X)
  out <- lapply(1:max(docClusterIdx), function(q){
    tf_q <- rowSums(X[, docClusterIdx==q, drop=FALSE])
    tf_o <- rowSums(X[, docClusterIdx!=q, drop=FALSE])
    p_q  <- (tf_q+alpha)/(sum(tf_q)+alpha*Tn)
    p_o  <- (tf_o+alpha)/(sum(tf_o)+alpha*Tn)
    lo   <- log(p_q/p_o)
    ord  <- order(lo, decreasing=TRUE)[1:top_n]
    tibble(doc_cluster=q, rank=1:top_n, term=terms[ord], log_odds=lo[ord])
  })
  bind_rows(out)
}
top_doc_markers <- top_terms_doc_side(X, terms, docClusterIdx, top_n=50)

# --------- 2) How word clusters relate to doc clusters (heatmap) ----------
#   Average TF-IDF of *each word-cluster* inside each *doc-cluster*
WQ <- sapply(1:Q, function(q){
  col_sel <- which(docClusterIdx==q)
  rowMeans(X[, col_sel, drop=FALSE])           # avg tf-idf per term in doc cluster q
})                                             # terms x Q
WQ_by_wordcluster <- aggregate(WQ, by=list(wordClusterIdx), FUN=mean)
colnames(WQ_by_wordcluster) <- c("word_cluster", paste0("doc_cluster_", 1:Q))
WQ_long <- pivot_longer(WQ_by_wordcluster, -word_cluster, names_to="doc_cluster", values_to="avg_tfidf")

ggplot(WQ_long, aes(x=doc_cluster, y=factor(word_cluster), fill=avg_tfidf))+
  geom_tile() + scale_fill_viridis_c() +
  labs(x="Document cluster", y="Word cluster", fill="Avg TF-IDF") +
  theme_minimal()

# --------- 3) Representative documents for each doc cluster ----------
# centroid (mean) of docs in cluster q; then cosine similarity
doc_centroids <- sapply(1:Q, function(q){
  rowMeans(X[, docClusterIdx==q, drop=FALSE])
}) # terms x Q

# cosine to each centroid
cos_to_centroid <- sapply(1:Q, function(q){
  cc <- doc_centroids[,q]
  # cosine(a,b) = (a·b)/(|a||b|)
  num <- crossprod(cc, X)                       # 1 x docs
  den <- sqrt(sum(cc^2)) * sqrt(colSums(X^2))   # 1 x docs
  (num/den)[1,]
}) # docs x Q

top_docs_per_cluster <- lapply(1:Q, function(q){
  ord <- order(cos_to_centroid[,q], decreasing=TRUE)[1:50]
  tibble(doc_cluster=q,
         rank=1:length(ord),
         doc_id = docs[ord],
         cosine = cos_to_centroid[ord,q])
})
top_docs <- bind_rows(top_docs_per_cluster)

# --------- 4) Internal quality: silhouette (cosine) ----------
# DOCS (may be heavy if many docs)
d_docs <- proxy::dist(t(X), method="cosine")
sil_docs <- silhouette(docClusterIdx, d_docs)
avg_sil_docs <- summary(sil_docs)$avg.width
avg_sil_docs

# WORDS (use a cap to keep cost reasonable)
max_words <- 3000
sel_words <- if(nrow(X) > max_words){
  order(rowSums(X), decreasing=TRUE)[1:max_words]
} else seq_len(nrow(X))

d_words <- proxy::dist(X[sel_words, , drop=FALSE], method="cosine")
sil_words <- silhouette(wordClusterIdx[sel_words], d_words)
avg_sil_words <- summary(sil_words)$avg.width
avg_sil_words


cat("\nSummary:\n",
    "K (word clusters)   :", K, "\n",
    "Q (doc clusters)    :", Q, "\n",
    "Avg sil (docs)      :", round(avg_sil_docs,3), "\n",
    "Avg sil (words)     :", round(avg_sil_words,3), "\n")





