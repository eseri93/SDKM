library(Matrix)
library(movMF)      # mixtures of von Mises–Fisher
library(mclust)     # adjustedRandIndex
library(R.matlab)   # readMat
library(tibble)
library(purrr)

wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)


## --- load your data (terms x documents, TF-IDF)
load("data_text.Rdata")        # assumes object 'data_text' exists
X <- as.matrix(data_text)      # terms = rows (N), docs = cols (J)

## --- load SDKM cluster assignments (from MATLAB)
wc <- readMat("wordClusterIdx.mat")$wordClusterIdx
dc <- readMat("docClusterIdx.mat")$docClusterIdx
wordClusterIdx <- as.vector(wc)
docClusterIdx  <- as.vector(dc)

## --- helper: L2-normalize rows, keep mask of nonzero rows
row_l2_normalize <- function(M) {
  nrm <- sqrt(rowSums(M^2))
  keep <- is.finite(nrm) & nrm > 0
  M_unit <- M[keep, , drop = FALSE]
  M_unit <- M_unit / nrm[keep]
  list(X = M_unit, keep = keep, nrm = nrm)
}

## --- (A) TERMS: vMF with K = 3 (match SDKM K)
set.seed(123)
K <- 3
rnormed <- row_l2_normalize(X)             # observations = term vectors in doc-space
X_terms_unit <- rnormed$X                  # rows on unit sphere, dim n_keep x J
keep_terms <- rnormed$keep

fit_terms <- movMF(X_terms_unit, k = K, control = list(nruns = 20))  # multiple starts

## hard assignments via log-densities (includes kappa constants)
logdens_terms <- sapply(
  1:K,
  function(k) dmovMF(X_terms_unit, fit_terms$theta[k, ], log = TRUE) + log(fit_terms$alpha[k])
)
term_clusters_vMF <- max.col(logdens_terms, ties.method = "first")

## map back to all terms (NA for zero rows if any)
term_clusters_all <- rep(NA_integer_, nrow(X))
term_clusters_all[keep_terms] <- term_clusters_vMF

## --- (B) DOCUMENTS: vMF with Q = 2 (match SDKM Q)
Q <- 2
# documents are columns; transpose so observations are rows
X_docs <- t(X)                                # now rows = documents, cols = terms
cnormed <- row_l2_normalize(X_docs)
X_docs_unit <- cnormed$X
keep_docs <- cnormed$keep

fit_docs <- movMF(X_docs_unit, k = Q, control = list(nruns = 20))

logdens_docs <- sapply(
  1:Q,
  function(q) dmovMF(X_docs_unit, fit_docs$theta[q, ], log = TRUE) + log(fit_docs$alpha[q])
)
doc_clusters_vMF <- max.col(logdens_docs, ties.method = "first")

## map back (should be all TRUE unless a zero-column exists)
doc_clusters_all <- rep(NA_integer_, ncol(X))
doc_clusters_all[keep_docs] <- doc_clusters_vMF

## --- (C) Comparisons to SDKM
# Adjusted Rand Index
ari_terms <- adjustedRandIndex(wordClusterIdx[keep_terms], term_clusters_vMF)
ari_docs  <- adjustedRandIndex(docClusterIdx[keep_docs],  doc_clusters_vMF)

cat(sprintf("vMF vs SDKM — ARI (terms): %.3f\n", ari_terms))
cat(sprintf("vMF vs SDKM — ARI (docs) : %.3f\n", ari_docs))

# Contingency tables (ignore NA if any)
tab_terms <- with(
  data.frame(SDKM = wordClusterIdx[keep_terms], vMF = term_clusters_vMF),
  table(SDKM, vMF)
)
tab_docs <- with(
  data.frame(SDKM = docClusterIdx[keep_docs], vMF = doc_clusters_vMF),
  table(SDKM, vMF)
)

tab_terms
tab_docs

## (optional) inspect concentration ||theta_k|| and mean directions
kappa_terms <- sqrt(rowSums(fit_terms$theta^2))
kappa_docs  <- sqrt(rowSums(fit_docs$theta^2))
kappa_terms; kappa_docs



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


X <- data_text
term_names <- rownames(X)
doc_names  <- colnames(X)

# If you filtered for vMF, keep the same indexing
# X_terms <- X[keep_terms, , drop = FALSE]
# X_docs  <- X[, keep_docs, drop = FALSE]

cos_sim <- function(a, b) {
  den <- sqrt(sum(a^2)) * sqrt(sum(b^2))
  if (den == 0) return(0)
  sum(a * b) / den
}

# ---------- SDKM: top terms for each TERM cluster (K=3) ----------
get_top_terms_sdkm <- function(X, wordClusterIdx, top_n = 30) {
  map_dfr(sort(unique(wordClusterIdx)), function(k) {
    rows_k   <- which(wordClusterIdx == k)
    centroid <- colMeans(X[rows_k, , drop = FALSE])
    centroid <- centroid / sqrt(sum(centroid^2) + 1e-12)
    scores   <- vapply(rows_k, function(i) cos_sim(X[i, ], centroid), numeric(1))
    ord      <- order(scores, decreasing = TRUE)
    tibble(
      method  = "SDKM",
      cluster = k,
      term    = term_names[rows_k][ord][seq_len(min(top_n, length(ord)))],
      score   = scores[ord][seq_len(min(top_n, length(ord)))]
    )
  })
}

# ---------- vMF (terms fit): top terms for each component ----------
# Fit on term-vectors (rows of X); "mu_k" is in document space; rank terms by cosine to mu_k.
get_top_terms_vmf_terms <- function(X, term_clusters_vMF, fit_terms, top_n = 30) {
  theta <- fit_terms$theta
  kappa <- sqrt(rowSums(theta^2))
  mu    <- sweep(theta, 1, kappa, "/")  # unit mean directions
  K     <- nrow(mu)
  map_dfr(seq_len(K), function(k) {
    rows_k <- which(term_clusters_vMF == k)
    scores <- vapply(rows_k, function(i) cos_sim(X[i, ], mu[k, ]), numeric(1))
    ord    <- order(scores, decreasing = TRUE)
    tibble(
      method  = "vMF-terms",
      cluster = k,
      term    = term_names[rows_k][ord][seq_len(min(top_n, length(ord)))],
      score   = scores[ord][seq_len(min(top_n, length(ord)))]
    )
  })
}

# ---------- SDKM: top terms characterizing each DOC cluster (Q=2) ----------
# Average TF-IDF across documents in each cluster; rank terms by that mean.
get_top_terms_for_doc_clusters_sdkm <- function(X, docClusterIdx, top_n = 30) {
  map_dfr(sort(unique(docClusterIdx)), function(q) {
    cols_q <- which(docClusterIdx == q)
    avg    <- rowMeans(X[, cols_q, drop = FALSE])
    ord    <- order(avg, decreasing = TRUE)
    tibble(
      method      = "SDKM",
      doc_cluster = q,
      term        = term_names[ord][seq_len(min(top_n, length(ord)))],
      weight      = avg[ord][seq_len(min(top_n, length(ord)))]
    )
  })
}

# ---------- vMF (docs fit): top terms per DOC component ----------
# Fit on document-vectors (columns of X); mu_q lies in term space directly.
get_top_terms_for_doc_clusters_vmf <- function(fit_docs, term_names, top_n = 30) {
  theta <- fit_docs$theta
  kappa <- sqrt(rowSums(theta^2))
  mu    <- sweep(theta, 1, kappa, "/")  # unit mean directions (in term space)
  Q     <- nrow(mu)
  map_dfr(seq_len(Q), function(q) {
    ord <- order(mu[q, ], decreasing = TRUE)
    tibble(
      method      = "vMF-docs",
      doc_cluster = q,
      term        = term_names[ord][seq_len(min(top_n, length(ord)))],
      loading     = mu[q, ord][seq_len(min(top_n, length(ord)))]
    )
  })
}

# ---------- Documents per cluster (lists) ----------
docs_by_cluster_sdkm <- tibble(doc = doc_names, cluster = docClusterIdx) %>%
  arrange(cluster, doc)

docs_by_cluster_vmf <- tibble(doc = doc_names, cluster = doc_clusters_vMF) %>%
  arrange(cluster, doc)

# If you used keep_docs, align names like:
# docs_by_cluster_vmf <- tibble(doc = doc_names[keep_docs], cluster = doc_clusters_vMF) %>%
#   arrange(cluster, doc)

# Identify the mismatching documents (if any)
mismatch_docs <- tibble(
  doc   = doc_names,             # or doc_names[keep_docs]
  SDKM  = docClusterIdx,         # or docClusterIdx[keep_docs]
  vMF   = doc_clusters_vMF
) %>% 
  filter(SDKM != vMF)

# ---------- Run and get tidy tables ----------
top_terms_sdkm_terms   <- get_top_terms_sdkm(X, wordClusterIdx, top_n = 30)
top_terms_vmf_terms    <- get_top_terms_vmf_terms(X, term_clusters_vMF, fit_terms, top_n = 30)

top_terms_sdkm_docs    <- get_top_terms_for_doc_clusters_sdkm(X, docClusterIdx, top_n = 30)
top_terms_vmf_docs     <- get_top_terms_for_doc_clusters_vmf(fit_docs, term_names, top_n = 30)

terms_comparison<-cbind(top_terms_sdkm_terms,top_terms_vmf_terms)
docs_comparison<-cbind(top_terms_sdkm_docs, top_terms_vmf_docs)

# Example: side-by-side comparison for interpretation
# (merge by cluster index and term, or just print them per cluster)



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#SIMULATION COMPARISON----------------------------------------------------------
#-------------------------------------------------------------------------------

library(Matrix)
library(movMF)      # mixtures of von Mises–Fisher
library(mclust)     # adjustedRandIndex
library(R.matlab)   # readMat
library(dplyr)
library(purrr)
library(stringr)
library(tibble)

wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

dat_e01 <- readMat("OutputSimData/simdata_e0p10.mat")
dat_e035 <- readMat("OutputSimData/simdata_e0p35.mat")
dat_e050 <- readMat("OutputSimData/simdata_e0p50.mat")
dat_e075 <- readMat("OutputSimData/simdata_e0p75.mat")
dat_e090 <- readMat("OutputSimData/simdata_e0p90.mat")
dat_e0110 <- readMat("OutputSimData/simdata_e1p10.mat")
dat_e0135 <- readMat("OutputSimData/simdata_e1p35.mat")
dat_e0150 <- readMat("OutputSimData/simdata_e1p50.mat")
dat_e0175 <- readMat("OutputSimData/simdata_e1p75.mat")
dat_e0200 <- readMat("OutputSimData/simdata_e2p00.mat")


library(R.matlab)
library(movMF)
library(mclust)
library(tibble); library(dplyr); library(purrr); library(stringr)

safe_l2_normalize_rows <- function(M) M / pmax(sqrt(rowSums(M^2)), 1e-12)
onehot_to_labels <- function(M) max.col(M, ties.method = "first")

# Try to find a variable whose name starts with a given prefix and is either
# a 3-D array (runs on 3rd dim) or a list of 2-D matrices (one per run).
.pick_container <- function(mm, prefix) {
  cands <- names(mm)[grepl(paste0("^", prefix), names(mm), ignore.case = TRUE)]
  for (nm in cands) {
    obj <- mm[[nm]]
    if (is.array(obj) && length(dim(obj)) == 3) {
      return(list(kind = "array3d", name = nm, obj = obj))
    }
    if (is.list(obj) && length(obj) > 0 && is.matrix(obj[[1]])) {
      return(list(kind = "list", name = nm, obj = obj))
    }
  }
  NULL
}

# Unified accessor that returns the r-th run matrices X (n×J), U_true (n×K), V_true (J×Q)
.get_run <- function(mm, Xc, Uc, Vc, r) {
  if (Xc$kind == "array3d") X <- Xc$obj[,,r, drop = FALSE][,,1] else X <- Xc$obj[[r]]
  if (Uc$kind == "array3d") U <- Uc$obj[,,r, drop = FALSE][,,1] else U <- Uc$obj[[r]]
  if (Vc$kind == "array3d") V <- Vc$obj[,,r, drop = FALSE][,,1] else V <- Vc$obj[[r]]
  list(X = X, U = U, V = V)
}

evaluate_vmf_matfile <- function(matfile, vmf_nstarts = 5, measure_time = FALSE) {
  mm <- readMat(matfile)
  
  # Find containers for X, U_true, V_true (names can be X_stack / X_list / etc.)
  Xc <- .pick_container(mm, "X")
  Uc <- .pick_container(mm, "U")
  Vc <- .pick_container(mm, "V")
  if (is.null(Xc) || is.null(Uc) || is.null(Vc)) {
    stop(sprintf("Could not find X/U/V containers in %s. Available vars: %s",
                 basename(matfile), paste(names(mm), collapse=", ")))
  }
  
  # Determine number of runs R
  if (Xc$kind == "array3d") {
    dims <- dim(Xc$obj)  # n × J × R
    n <- dims[1]; J <- dims[2]; R <- dims[3]
  } else {
    # list of matrices
    R <- length(Xc$obj)
    n <- nrow(Xc$obj[[1]]); J <- ncol(Xc$obj[[1]])
  }
  
  ari_rows <- numeric(R); ari_cols <- numeric(R)
  t_rows <- if (measure_time) numeric(R) else NULL
  t_cols <- if (measure_time) numeric(R) else NULL
  
  for (r in seq_len(R)) {
    run <- .get_run(mm, Xc, Uc, Vc, r)
    X <- run$X; U_true <- run$U; V_true <- run$V
    
    # Rows (units)
    Xr <- safe_l2_normalize_rows(X)
    if (measure_time) { t0 <- proc.time() }
    fit_r <- movMF(Xr, k = ncol(U_true), control = list(nruns = vmf_nstarts))
    if (measure_time) { t_rows[r] <- (proc.time() - t0)[["elapsed"]] }
    zhat_r  <- predict(fit_r, Xr, type = "class")
    ztrue_r <- onehot_to_labels(U_true)
    ari_rows[r] <- adjustedRandIndex(zhat_r, ztrue_r)
    
    # Columns (variables)
    Xc2 <- safe_l2_normalize_rows(t(X))  # J × n
    if (measure_time) { t0 <- proc.time() }
    fit_c <- movMF(Xc2, k = ncol(V_true), control = list(nruns = vmf_nstarts))
    if (measure_time) { t_cols[r] <- (proc.time() - t0)[["elapsed"]] }
    zhat_c  <- predict(fit_c, Xc2, type = "class")
    ztrue_c <- onehot_to_labels(V_true)
    ari_cols[r] <- adjustedRandIndex(zhat_c, ztrue_c)
  }
  
  eps <- str_match(basename(matfile), "e(\\d+p\\d+)")[,2] |> str_replace("p", ".") |> as.numeric()
  
  tibble(
    file = basename(matfile), eps = eps,
    n = n, J = J,
    runs = R,
    ARI_rows_mean = mean(ari_rows), ARI_rows_median = median(ari_rows),
    ARI_rows_exact_n = sum(ari_rows > 1 - 1e-12),
    ARI_rows_exact_prop = mean(ari_rows > 1 - 1e-12),
    ARI_cols_mean = mean(ari_cols), ARI_cols_median = median(ari_cols),
    ARI_cols_exact_n = sum(ari_cols > 1 - 1e-12),
    ARI_cols_exact_prop = mean(ari_cols > 1 - 1e-12),
    time_rows_median_s = if (measure_time) median(t_rows) else NA_real_,
    time_cols_median_s = if (measure_time) median(t_cols) else NA_real_
  )
}

# --- Run across files (same as before) ---
mat_files <- c(
  "OutputSimData/simdata_e0p10.mat",
  "OutputSimData/simdata_e0p35.mat",
  "OutputSimData/simdata_e0p50.mat",
  "OutputSimData/simdata_e0p75.mat",
  "OutputSimData/simdata_e0p90.mat",
  "OutputSimData/simdata_e1p10.mat",
  "OutputSimData/simdata_e1p35.mat",
  "OutputSimData/simdata_e1p50.mat",
  "OutputSimData/simdata_e1p75.mat",
  "OutputSimData/simdata_e2p00.mat"
)

set.seed(123)
vmf_summary <- map_dfr(mat_files, ~ evaluate_vmf_matfile(.x, vmf_nstarts = 5, measure_time = FALSE)) %>%
  arrange(eps)

print(vmf_summary, n = Inf)
