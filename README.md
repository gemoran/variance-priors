
This repository contains R code to reproduce the results in Moran et al (2018).

This code uses the Spike-and-Slab Lasso (with both the fixed and unknown variance cases) which is implemented in the `SSLASSO` R package, available on [CRAN](https://CRAN.R-project.org/package=SSLASSO). 

### Index
- "section_3-2.R"

This R script runs the simulation study in Section 3.2 and generates Figure 1. 

- "section_6-6.R"

This R script runs the simulation study in Section 6.6 and generates Table 1 and Figure 2.

- "section_7.R" 

This R script runs the analysis of the protein activity data of Clyde and Parmigiani (1998) which is detailed in Section 7.

- `scaledSSL`
This package implements the "scaled Spike-and-Slab Lasso" as described in Section 6.5 of Moran et al (2018). It is provided for the sole purpose of reproducing the results in "section_6.6.R". If you would like to use the Spike-and-Slab Lasso with unknown variance, please use the package `SSLASSO` available on CRAN [here](https://CRAN.R-project.org/package=SSLASSO). To install `scaledSSL` in order to run the code in "section_6.6.R", use the following command:
```
install_github("gemma-e-moran/variance-priors/scaledSSL")
```

### Reference
- Moran, G. E., Rockova, V. and George, E. I. (2018) "Variance prior forms for high-dimensional Bayesian variable selection'' [[arXiv:1801.03019]](https://arxiv.org/abs/1801.03019)
- Clyde, M. A. and Parmigiani, G. (1998) "Protein construct storage: Bayesian variable selection and prediction with mixtures" *Journal of Biopharmaceutical Statistics*, 8, 431-443 
