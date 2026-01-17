# SDKM
Repository associated with the paper "Spherical Double K-Means: a co-clustering approach for textual data analysis". It contains all the matlab and R codes for the methodology and applications.

- SDKM_tesi.m contains the SDKM code,
- randPU.m it is called in SDKM_tesi.m to generate random partitions,
- pseudo_F_index.m contains the pseudo-F code for choosing the number of clusters.
- GenDataSDKM.m contains the code to generate SDKM data (for simulations).
- scenarios_simulation.m run simulations with different level of errors.
- application_SDKM.m perform the SDKM algorithm on "data_text.csv", i.e. the cleaned terms-documents frequency matrix of the presidential inaugural addresses dataset.
- applicationNEWS_SDKM.m perform the SDKM algorithm on "data_news_sample.csv", i.e. the cleaned terms-documents frequency matrix of the 20 newsgroup dataset.

- data_text.csv contains the cleaned terms-documents frequency matrix of the presidential inaugural addresses dataset (from the quanteda R package).
- data_news_sample.csv contains the cleaned terms-documents frequency matrix of the 20 newsgroup dataset sampled. The full dataset is freely available online https://www.kaggle.com/datasets/crawford/20-newsgroups.

- Presidents.R contains the data loading from the quanteda package of the U.S. presidential inauguaral adresses, the data preprocessing and saving. Then, after the SDKM application in matlab, the cluster assigments are loaded for the plotting and interpretation of the results.
- 20news.R contains the data loading and preprocessing of the 20newsgoup corpus. Such as for the presidents application, the SDKM is done with matlab, and then the cluster assigments are loaded for the plotting and interpretation of the results.
- movMFcomparison.R contains the Comparison with mixtures of von Mises–Fisher (vMF) components on both U.S. presidential data and simulated data.
