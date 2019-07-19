
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
    ##             [,1]       [,2]        [,3]       [,4]
    ## [1,] 0.132282030 0.01075479 0.003694099 0.00287548
    ## [2,] 0.010754792 0.18691262 0.050635856 0.03563093
    ## [3,] 0.003694099 0.05063586 0.334360144 0.14786543
    ## [4,] 0.002875480 0.03563093 0.147865432 0.51882916
    ## 
    ## $B
    ##           [,1]       [,2]       [,3]        [,4]
    ## [1,] -3.809909  0.3703037  0.0000000  0.00000000
    ## [2,]  0.000000 -2.8515474  0.5970447  0.07741383
    ## [3,]  0.000000  0.0000000 -1.7972632  0.68260157
    ## [4,]  0.000000  0.0000000  0.1851307 -1.01647030
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
    ## [1] 1.357952e-16
    ## 
    ## $logLikl1
    ## [1] -1.635144
    ## 
    ## $iter
    ## [1] 29
    ## 
    ## $job
    ## [1] 0

``` r
### penalized least square 0.5*||S(B) - Sigma||_2^2 + lambda * ||B||_1,off
results2 <- proxgradlsB(cov(sample), B= B0, C = C, eps = 1e-12, 
                       maxIter = 1000, lambda = 0.01, job = 0)
results2$B
```

    ##           [,1]     [,2]       [,3]       [,4]
    ## [1,] -3.877387  0.00000  0.0000000  0.0000000
    ## [2,]  0.000000 -2.72376  0.2127283  0.1684361
    ## [3,]  0.000000  0.00000 -1.6548317  0.6623955
    ## [4,]  0.000000  0.00000  0.0000000 -0.9327754

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

    ##       lambda npar fp tp fn tn errs
    ##  [1,]  0.025    8  1  7  0  8    1
    ##  [2,]  0.050    7  0  7  0  9    0
    ##  [3,]  0.075    6  0  6  1  9    1
    ##  [4,]  0.100    6  0  6  1  9    1
    ##  [5,]  0.125    6  0  6  1  9    1
    ##  [6,]  0.150    6  0  6  1  9    1
    ##  [7,]  0.175    5  0  5  2  9    2
    ##  [8,]  0.200    5  0  5  2  9    2
    ##  [9,]  0.225    5  0  5  2  9    2
    ## [10,]  0.250    5  0  5  2  9    2

Related code
------------

-   Some code is from the `lyapunov` package (<https://github.com/gragusa/lyapunov>).

-   Using FORTRAN code of Algorithm 705 from the Collected Algorithms from ACM, TRANSACTIONS ON MATHEMATICAL SOFTWARE, VOL. 18, NO. 2, PP. 232-238.75.
