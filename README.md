# SBFA
Integrative Analysis of Multi-omics and Imaging Data with Incorporation of Biological Information via Structural Bayesian Factor Analysis

## R session infomation
```console
R version 4.2.0 (2022-04-22 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19044)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] lamW_2.1.1      pracma_2.3.8    dplyr_1.0.9     stringr_1.4.0   plyr_1.8.7      reticulate_1.25

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.8.3       rstudioapi_0.13    magrittr_2.0.3     rappdirs_0.3.3     tidyselect_1.1.2  
 [6] here_1.0.1         lattice_0.20-45    R6_2.5.1           rlang_1.0.2        fansi_1.0.3       
[11] tools_4.2.0        grid_4.2.0         png_0.1-7          utf8_1.2.2         DBI_1.1.2         
[16] cli_3.3.0          ellipsis_0.3.2     RcppParallel_5.1.5 assertthat_0.2.1   rprojroot_2.0.3   
[21] tibble_3.1.7       lifecycle_1.0.1    crayon_1.5.1       Matrix_1.4-1       purrr_0.3.4       
[26] vctrs_0.4.1        glue_1.6.2         stringi_1.7.6      compiler_4.2.0     pillar_1.7.0      
[31] generics_0.1.2     jsonlite_1.8.0     pkgconfig_2.0.3 
```
## SBFA usage
Please refer to our [EXAMPLE](./R_main/Example.R).
Input arguments: 
- X: Multiomics dataset with (number of rows = number of features) and (number of columns = number of subjects). This is a vertical concatenation of multiomics datasets. 

```For example, if you want to extract the common information from presence of effect allele recoded genotyping data (0 = No effect allele, 1 = 1 or 2 effect alleles; assume we have 10 candidate SNVs), normalized gene expression data (5 genes), and normalized regional imaging volumetric changes measurements (20 region of interests), then the X matrix should contains (5+10+20=35) rows/features. Say if we can get the the data of all three modalities from 50 common subjects (no missing values), then the X matrix should contains 50 columns/subjects.```

- type: A vector of length equal to number of features, indicating the underlying distribution for each feature. Each element of the feature is 0, 1, or 2 (0=Normal distribution; 1=Binomial distribution; 2=Negative binomial distribution). 

```In our example, since X matrix has 35 columns (10 presence of effect allele recoded genotyping data following Bernoulli distribution; 5 gene expression data following Gaussian distribution; 20 regional imaging measurements following Gaussian distribution), the input for type parameter will be'''
'''coffee
c(rep(1,10),rep(0,5),rep(0,20))
```.
