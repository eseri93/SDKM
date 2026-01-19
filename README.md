SDKM

Repository associated with the paper “Spherical Double K-Means: a co-clustering approach for textual data analysis”. It contains the MATLAB and R code for the proposed methodology, simulations, and applications.

MATLAB files

- SDKM.m: implementation of the Spherical Double K-Means (SDKM) algorithm.
- pseudo_F_index.m: pseudo-F index used to select the number of clusters (𝐾,𝑄).
- GenDataSDKM.m: synthetic data generator for SDKM (used in simulations).
- scenarios_simulation.m: simulation study under different noise/error levels.
- bghungar.m: robust MATLAB implementation of the Hungarian algorithm for square assignment problems (used here to handle label switching by optimally matching cluster labels).
- application_SDKM.m: runs SDKM on data_text.csv (U.S. presidential inaugural addresses term–document matrix).
- applicationNEWS_SDKM.m: runs SDKM on data_news_sample.csv (sampled 20 Newsgroups term–document matrix).
- DKM.m: implementation of the Double K-Means (DKM) baseline algorithm.
- applicationDKM.m: runs DKM on data_text.csv (U.S. presidential inaugural addresses term–document matrix).

Data files

- data_text.csv: cleaned term–document matrix for the U.S. presidential inaugural addresses dataset (from the quanteda R package).
- data_news_sample.csv: cleaned term–document matrix for a stratified sample of the 20 Newsgroups dataset. The full dataset is freely available online (e.g., Kaggle: https://www.kaggle.com/datasets/crawford/20-newsgroups).

R files

- Presidents.R: loads the U.S. presidential inaugural addresses corpus (quanteda), performs preprocessing, and exports the term–document matrix. After running SDKM in MATLAB, it imports cluster assignments and produces plots/interpretation outputs.
- 20news.R: loads and cleans the 20 Newsgroups email corpus, builds a trimmed TF-IDF term–document matrix (exported for MATLAB SDKM/DKM), optionally creates a stratified sample, and imports MATLAB SDKM results to produce basic diagnostics (cluster summaries, heatmaps, representative documents, and cosine-silhouette scores).
- movMFcomparison.R: fits mixtures of von Mises–Fisher distributions to term and document vectors (TF-IDF, L2-normalized) and compares the resulting hard partitions with SDKM assignments via ARI, contingency tables, and interpretable top-term summaries; it also optionally evaluates vMF cluster recovery on simulated datasets saved as .mat files.
- DKM_application_and_comparison.R: loads TF-IDF data and MATLAB outputs from DKM and SDKM, produces comparable plots for top words and document cluster assignments, aligns DKM labels to SDKM via optimal assignment (Hungarian/LSAP) using overlap information, and exports a tab-separated file listing each word with its DKM and SDKM cluster labels.
