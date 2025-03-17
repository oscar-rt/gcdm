# GCDM
Genetic Co-expression Diffusion Maps (GCDM), an unsupervised learning method for dimensionality reduction based on Diffusion Maps, designed to incorporate gene co-expression information from local genetic networks by employing Gromov-Wasserstein distance for graphs.

Dimensionality reduction is one of the workhorses in the single-cell RNA sequencing (scRNA-seq) data analysis pipeline. It constitutes the first step of many scRNA-seq data analysis algorithms, particularly algorithms for clustering and trajectory inference. Although many methods have been proposed recently, these methods concentrate on the similarities between individual cells given by Euclidean distances. Due to a large number of drop-out observations, the similarity may not reflect the actual closeness of cells. Besides, the existing methods do not incorporate the genetic co-expression information into the process. We present Genetic Co-expression Diffusion Maps (GCDM), a dimensionality reduction method based on Diffusion Maps, designed to incorporate gene co-expression information from local genetic networks by employing Gromov-Wasserstein distance for graphs. The similarities are calculated by the distance between genetic networks in single cells. Numerical results show that our proposed GCDM algorithm is an effective method for dimensionality reduction with very good accuracy and robustness properties.

Refer to the included Jupyter notebook for an example.

## Authors

- [@oscar-rt](https://www.github.com/oscar-rt)
