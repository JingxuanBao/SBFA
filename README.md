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

```coffee
SBFA_EM <- function(X,type,param,E,L,v1,v2,a_omega,b_omega,m.init=1,scale=T,W.init=NULL,eps=1e-3,maxIter=500)
```

Input arguments: 
- X: Multiomics dataset with (number of rows = number of features) and (number of columns = number of subjects). This is a vertical concatenation of multiomics datasets. 

*For example, if you want to extract the common information from presence of effect allele recoded genotyping data (0 = No effect allele, 1 = 1 or 2 effect alleles; assume we have 10 candidate SNVs), normalized gene expression data (5 genes), and normalized regional imaging volumetric changes measurements (20 region of interests), then the X matrix should contains (5+10+20=35) rows/features. Say if we can get the the data of all three modalities from 50 common subjects (no missing values), then the X matrix should contains 50 columns/subjects.*

- type: A vector of length equal to number of features, indicating the underlying distribution for each feature. Each element of the feature is 0, 1, or 2 (0=Gaussian distribution; 1=Binomial distribution; 2=Negative binomial distribution). 

*In our example, since X matrix has 35 columns (10 presence of effect allele recoded genotyping data following Bernoulli distribution; 5 gene expression data following Gaussian distribution; 20 regional imaging measurements following Gaussian distribution), the input for type parameter will be*

```coffee
c(rep(1,10),rep(0,5),rep(0,20))
```

- param: A vector of length equal to number of features, indicating the necessary parameters of underlying distribution for each feature (Gaussian: the standard deviation of Gaussian distribution; Binomial: number of trails parameter; Negative binomial: number of successes parameter).

*In our example, since X matrix has 35 columns (10 presence of effect allele recoded genotyping data following Bernoulli distribution, which is also binomial distribution with number of trails equal to 1; 5 gene expression data following Gaussian distribution - since the data has been normalized, the standard deviation is 1; 20 regional imaging measurements following Gaussian distribution - since the data has been normalized, the standard deviation is 1), the input for type parameter will be*

```coffee
c(rep(1,10),rep(1,5),rep(1,20))
```
- E: A two-column matrix with number of rows equal to total number of edges. 

*In our example, we assume 1) SNV1 (1st row in X), SNV3 (3rd row in X), SNV4 (4th row in X) are in the same LD block; and all the other SNVs are in different LD block (also different from SNV1, SNV3, and SNV4). 2) gene1 (11th row in X) and gene3 (13th row in X) are in the same pathway or enriched by the same module; and all the other genes are in different pathways (also different from gene1 and gene3 pathway). 3) region3 (18th row in X) and region5 (20th row in X) are in the same brain tissue; region2 (17th row in X), region7 (22th row in X) and region10 (25th row in X) are in the another brain tissue; all the other regions are in different brain tissue and different from either regioin3/region5 tissue or region2/region7/region10 tissue. Then the E matrix would be:*

```console
       [,1] [,2]
  [1,]    1    3
  [2,]    1    4
  [3,]    3    4
  [4,]   11   13
  [5,]   18   20 
  [6,]   17   22
  [7,]   17   25
  [8,]   22   25
  [9,]    3    1
 [10,]    4    1
 [11,]    4    2
 [12,]   13   11
 [13,]   20   18
 [14,]   22   17
 [15,]   25   17
 [16,]   25   22
```
- L: A scalar, denotes the number of latent dimensions.

*You may want to try different L and select the one gives you the lowest BIC. In our example, we can try 2,3,4 separately. For the tuning precess, please refer to SBFA_EM_AUTOT function*

- v1: A scalar, controls the sparsity. 

*You may want to try different v1 and select the one gives you the lowest BIC. In our example, we can try 3,4,5. For the tuning precess, please refer to SBFA_EM_AUTOT function. For the specific meaning of v1, please refer to our SBFA paper "Prior for λ: Incorporating Biological Knowledge" section.*

- v2: A scalar, set as log2. 

*Typically, you don't need to tune this parameter (but you can tune this parameter in our SBFA_EM_AUTOT function). For the specific meaning of v2, please refer to our SBFA paper "Prior for λ: Incorporating Biological Knowledge" section.*

- a_omega: A sccalar, denotes how much weights you want to give to your graph information. Larger a_omega corresponds to larger weights. 

*You can tune this parameter. In our example, we can fix it to be 3. For the tuning precess, please refer to SBFA_EM_AUTOT function. For the specific meaning of a_omega, please refer to our SBFA paper "Variational EM algorithm" section.*

- b_omega: A sccalar, set as 1.

*Typically, you don't need to tune this parameter (but you can tune this parameter in our SBFA_EM_AUTOT function). For the specific meaning of b_omega, please refer to our SBFA paper "Variational EM algorithm" section.*

- m.init: A sccalar or a vector of length equal to number of features. It is the straitforward estimation of the location parameter *m* used as initialization. If it is a scalar, it can be 0 (all features have initialization location parameter 0), 1 (the location parameters for all features are estimated using trimed mean across all subjects), 2 (the location parameters for all features are estimated using median across all subjects), or 3 (the location parameters for all features are estimated using mean across all subjects). If it is a vector of length equal to number of features, each element corresponds to the initialization of location parameter for the corresponding feature. Default is 1.

- scale: A boolean variable indicating whether you want to scale the data or not. Default is True.

*When reading data, different variables may in different scale. We may want to put all the variables in the same scale. For example, when we want to add penalty term to our model, the smaller scale variable tends to have larger penalty. So, it is often the case that we want to scale the variable. For the continuous variable like Gaussian distribution, it is easy to change the scale by just standardizing alll the variable. But for the discrete variable, it is no longer trivial to normalize, because when you multiple something, it won't be discrete any more. So, in our case, when we scale the variables, we didn't scale the matrix X. Instead, we scale the matrix mu. Our approach is to calculate the standard deviation of transformed data of X.*

- W.init: Initialization of factor loading matrix W. Defalt is NULL, which means the algorithm will use SVD decomposition to initialize W.

*The factor analysis problems suffer from an identification problem: As long as WxZ is the same, there is no difference. For example, we can multiply some matrix to W and then multiple the inverse of the matrix to Z. Then the resulting matrix mu is the same. There are multiple solutions to the identification problem such as put some penalty (it may give unique solution but it cannot guarantee to find the unique solution because of lots of local optimal solution). Thus, because of the identification problem and the unique solution problem, it is important for us to start with a good initial value. Random initialization implies random solution.*

*In our setting, we choose W and Z matrices close to orthogonal matrices. We do svd to the original matrix X and assign the LHS to be W and RHS to be Z. In this way we can easily find the good initial value for W and Z. This is not always the case. But if we choose the W and Z to be close to orthogonal matrix, then the SVD gives a reasonably good initial values.*

- eps: Convergence criteria. Defalt is 1e-3.

- maxIter: Maximum iterations if not convergent. Defalt is 500.

