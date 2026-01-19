# SDKM
Repository associated with the paper "Spherical Double K-Means: a co-clustering approach for textual data analysis". It contains all the matlab and R codes for the methodology and applications.

- SDKM.m contains the SDKM code.
- pseudo_F_index.m contains the pseudo-F code for choosing the number of clusters.
- GenDataSDKM.m contains the code to generate SDKM data (for simulations).
- scenarios_simulation.m run simulations with different level of errors.
- bghungar.m implements a robust MATLAB version of the Hungarian algorithm to solve square assignment problems (used here to handle label switching by optimally matching cluster labels).
- application_SDKM.m perform the SDKM algorithm on "data_text.csv", i.e. the cleaned terms-documents frequency matrix of the presidential inaugural addresses dataset.
- applicationNEWS_SDKM.m perform the SDKM algorithm on "data_news_sample.csv", i.e. the cleaned terms-documents frequency matrix of the 20 newsgroup dataset.
- DKM.m contains the DKM code
- application_DKM.m perform the DKM algorithm on "data_text.csv", i.e. the cleaned terms-documents frequency matrix of the presidential inaugural addresses dataset.

- data_text.csv contains the cleaned terms-documents frequency matrix of the presidential inaugural addresses dataset (from the quanteda R package).
- data_news_sample.csv contains the cleaned terms-documents frequency matrix of the 20 newsgroup dataset sampled. The full dataset is freely available online https://www.kaggle.com/datasets/crawford/20-newsgroups.

- Presidents.R contains the data loading from the quanteda package of the U.S. presidential inauguaral adresses, the data preprocessing and saving. Then, after the SDKM application in matlab, the cluster assigments are loaded for the plotting and interpretation of the results.
- 20news.R loads and cleans the 20 Newsgroups email corpus, builds a trimmed TF-IDF term–document matrix (exported for MATLAB SDKM/DKM), optionally creates a stratified sample, and imports MATLAB SDKM results to produce basic diagnostics (cluster summaries, heatmaps, representative documents, and cosine-silhouette scores).
- movMFcomparison.R contains the Comparison with mixtures of von Mises–Fisher (vMF) components on both U.S. presidential data and simulated data.
- COMPARISON.R contains the comparison between the cluster assignments of SDKM and DKM on the U.S. presidential data.
