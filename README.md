
clggm
=====

`clggm` contatins methods to estimate the coefficients of the stochastic process (Ornstein–Uhlenbeck) *d**X*<sub>*t*</sub> = *B*(*X*<sub>*t*</sub> − *a*)*d**t* + *D**d**W*<sub>*t*</sub> from observations of the invariant distribution.

Installation
------------

``` r
devtools::install_github("gherardovarando/clggm")
```

### Dependencies

-   `MASS`
-   `stats`
-   `genlasso` (only for the `genglasso` function), it can be obtained with:

``` r
devtools::install_github("glmgen/genlasso")
```

-   `testthat` (only in dev)

Usage
-----

``` r
library(clggm)

### define coefficient matrices
B <- matrix(nrow = 4, c(-4, 1,   0,  0, 
                         0, -3,  1,  0,
                         0,  0, -2,  1,
                         0,  0,  0, -1), byrow = TRUE)
D <- diag(c(1,1,1,1))
C <- D %*% t(D)

### solve continuous Lyapunov equation 
### to obtain covariance of invariant 
### distribution
Sigma <- clyap(B, C) 

### obtain observations from the invariant distribution
sample <- MASS::mvrnorm(n = 1000, mu = rep(0,4),Sigma =  Sigma)
```

### Estimate *B* knowing *C*

``` r
B0 <- - 0.5 * C %*% solve(cov(sample))
### penalized maximum-likelihood solved with proximal gradient
results <- proxgradllB(cov(sample), B= B0, C = C, eps = 1e-12, 
                       maxIter = 1000, lambda = 0.05, job = 0)
results 
```

    ## $N
    ## [1] 4
    ## 
    ## $Sigma
    ##             [,1]        [,2]        [,3]        [,4]
    ## [1,] 0.137000528 0.009879739 0.003578185 0.001539638
    ## [2,] 0.009879739 0.188136977 0.055622289 0.021826746
    ## [3,] 0.003578185 0.055622289 0.331898223 0.116011024
    ## [4,] 0.001539638 0.021826746 0.116011024 0.498167118
    ## 
    ## $B
    ##           [,1]       [,2]       [,3]       [,4]
    ## [1,] -3.673413  0.3299177  0.0000000  0.0000000
    ## [2,]  0.000000 -2.8732919  0.7294278  0.0000000
    ## [3,]  0.000000  0.0000000 -1.7289162  0.6363551
    ## [4,]  0.000000  0.0000000  0.0000000 -1.0036793
    ## 
    ## $C
    ##      [,1] [,2] [,3] [,4]
    ## [1,]    1    0    0    0
    ## [2,]    0    1    0    0
    ## [3,]    0    0    1    0
    ## [4,]    0    0    0    1
    ## 
    ## $lambda
    ## [1] 0.05
    ## 
    ## $diff
    ## [1] 1.526381e-15
    ## 
    ## $logLikl1
    ## [1] -1.600184
    ## 
    ## $iter
    ## [1] 33
    ## 
    ## $job
    ## [1] 0

``` r
### penalized least square 0.5*||S(B) - Sigma||_2^2 + lambda * ||B||_1,off
results2 <- proxgradlsB(cov(sample), B= B0, C = C, eps = 1e-12, 
                       maxIter = 1000, lambda = 0.01, job = 0)
results2$B
```

    ##           [,1]      [,2]       [,3]       [,4]
    ## [1,] -3.735025  0.000000  0.0000000  0.0000000
    ## [2,]  0.000000 -2.775588  0.3512418  0.0000000
    ## [3,]  0.000000  0.000000 -1.5716932  0.5593472
    ## [4,]  0.000000  0.000000  0.0000000 -0.9702106

#### Solutions along a regularization path

``` r
results <- llBpath(Sigma, eps = 1e-12, maxIter = 1000, job = 10)
t(sapply(results, function(res) c(lambda = res$lambda, 
                                npar = sum(res$B!=0),
                                fp = sum(res$B!=0 & B==0),
                                tp = sum(res$B!=0 & B!=0) ,
                                fn = sum(res$B==0 & B!=0),
                                tn = sum(res$B==0 & B==0),
                                errs = sum(res$B!=0 & B==0) + 
                                  sum(res$B==0 & B!=0))))
```

    ##           lambda npar fp tp fn tn errs
    ##  [1,] 0.02500000    8  1  7  0  8    1
    ##  [2,] 0.07777778    6  0  6  1  9    1
    ##  [3,] 0.13055556    6  0  6  1  9    1
    ##  [4,] 0.18333333    5  0  5  2  9    2
    ##  [5,] 0.23611111    5  0  5  2  9    2
    ##  [6,] 0.28888889    5  0  5  2  9    2
    ##  [7,] 0.34166667    5  0  5  2  9    2
    ##  [8,] 0.39444444    4  0  4  3  9    3
    ##  [9,] 0.44722222    4  0  4  3  9    3
    ## [10,] 0.50000000    4  0  4  3  9    3
