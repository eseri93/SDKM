**SDKM**

Code and data repository for “Spherical Double K-Means: a co-clustering approach for textual data analysis” by Ilaria Bombelli, Domenica Fioredistella Iezzi, Emiliano Seri, and Maurizio Vichi. The repository includes MATLAB, R scripts and datasets to reproduce the analyses in the paper.

*MATLAB files*

- SDKM.m: implementation of the Spherical Double K-Means (SDKM) algorithm.
- pseudo_F_index.m: pseudo-F index used to select the number of clusters (𝐾,𝑄).
- GenDataSDKM.m: synthetic data generator for SDKM (used in simulations).
- scenarios_simulation.m: simulation study under different noise/error levels.
- bghungar.m: robust MATLAB implementation of the Hungarian algorithm for square assignment problems (used here to handle label switching by optimally matching cluster labels).
- application_SDKM.m: runs SDKM on data_text.csv (U.S. presidential inaugural addresses term–document matrix).
- applicationNEWS_SDKM.m: runs SDKM on data_news_sample.csv (sampled 20 Newsgroups term–document matrix).
- DKM.m: implementation of the Double K-Means (DKM) baseline algorithm.
- application_DKM.m: runs DKM on data_text.csv (U.S. presidential inaugural addresses term–document matrix).
- choose_RndStarts.m evaluates the impact of the number of random starts (Rndst) in SDKM on avoiding local maxima, by repeatedly running SDKM on synthetic data and comparing the achieved objective values to the ground-truth objective. This is a diagnostic/tuning script; not required to reproduce the paper’s main results.

*Data files*

- data_text.csv: cleaned term–document matrix for the U.S. presidential inaugural addresses dataset (from the quanteda R package).
- data_news_sample.csv: cleaned term–document matrix for a stratified sample of the 20 Newsgroups dataset. The full dataset is freely available online (e.g., Kaggle: https://www.kaggle.com/datasets/crawford/20-newsgroups).

*R files*

- Presidents.R: loads the U.S. presidential inaugural addresses corpus (quanteda), performs preprocessing, and exports the term–document matrix. After running SDKM in MATLAB, it imports cluster assignments and produces plots/interpretation outputs.
- 20news.R: loads and cleans the 20 Newsgroups email corpus, builds a trimmed TF-IDF term–document matrix (exported for MATLAB SDKM/DKM), optionally creates a stratified sample, and imports MATLAB SDKM results to produce basic diagnostics (cluster summaries, heatmaps, representative documents, and cosine-silhouette scores).
- movMFcomparison.R: fits mixtures of von Mises–Fisher distributions to term and document vectors (TF-IDF, L2-normalized) and compares the resulting hard partitions with SDKM assignments via ARI, contingency tables, and interpretable top-term summaries; it also optionally evaluates vMF cluster recovery on simulated datasets saved as .mat files.
- DKM_application_and_comparison.R: loads TF-IDF data and MATLAB outputs from DKM and SDKM, produces comparable plots for top words and document cluster assignments, aligns DKM labels to SDKM via optimal assignment (Hungarian/LSAP) using overlap information, and exports a tab-separated file listing each word with its DKM and SDKM cluster labels.

**How to run**
1) *Requirements*

- MATLAB (tested with recent versions)
  - Statistics and Machine Learning Toolbox (required; used for functions such as rmse, etc.)
- R (≥ 4.0 recommended) with packages used in the scripts (e.g., quanteda, readtext, R.matlab, tidyverse, movMF, mclust, etc.)

Make sure your working directory is the repository root (or the folder containing the scripts and data files).

2) *Presidential inaugural addresses (data_text)*

R preprocessing

1. Run Presidents.R to:
   - load the quanteda inaugural corpus,
   - preprocess text,
   - build the TF-IDF term–document matrix,
   - export data_text.csv (and optionally data_text.Rdata).

MATLAB clustering

2. Run application_SDKM.m to:
   - compute pseudo-F for (𝐾,𝑄)∈{2,…,6}
   - run SDKM with K=3, Q=2,
   - save Usdkm.mat, Vsdkm.mat, Ymsdkm.mat, fsdkm.mat.

3. (Optional baseline) Run applicationDKM.m to run DKM and save Udkm.mat, Vdkm.mat, Ymdkm.mat, fdkm.mat.

R post-processing / plots

4. Run Presidents.R again (plot/interpretation section) to load MATLAB outputs (e.g., wordClusterIdx.mat, docClusterIdx.mat, or U/V matrices depending on your export) and generate plots and summaries.

3) *20 Newsgroups (data_news_sample)*

R preprocessing + sampling

1. Run 20news.R to:
    - load/clean the 20 Newsgroups corpus,
    - build the TF-IDF term–document matrix,
    - optionally sample a stratified subset,
    - export data_news_sample.csv (and optionally data_news_samp.Rdata).

MATLAB clustering

2. Run applicationNEWS_SDKM.m to run SDKM on data_news_sample.csv and save:
   - Usdkm_news.mat, Vsdkm_news.mat, Ymsdkm_news.mat, fsdkm_news.mat.

R diagnostics / interpretation

3. Run the analysis section of 20news.R to import the MATLAB outputs and compute diagnostics (cluster summaries, heatmaps, representative documents, silhouette scores).

4) *Simulation study*

   1. Run scenarios_simulation.m to generate synthetic datasets and evaluate SDKM under different noise levels.
   2. Run the simulation section in movMFcomparison.R to evaluate vMF recovery on the simulated .mat datasets.

5) *DKM vs SDKM comparison (Presidents)*

Run DKM_application_and_comparison.R after producing both SDKM and DKM MATLAB outputs on data_text.csv to generate side-by-side plots and label-aligned comparisons.
