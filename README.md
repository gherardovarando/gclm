
# clggm

`clggm` contatins methods to estimate the coefficients of the stochastic
process (Ornstein–Uhlenbeck) \(dX_t = B(X_t - a)dt + DdW_t\) from
observations of the invariant distribution.

## Installation

``` r
devtools::install_github("gherardovarando/clggm")
```

### Dependencies

  - `MASS`
  - `stats`
  - `genlasso` (only for the `genglasso` function), it can be obtained
    with:

<!-- end list -->

``` r
devtools::install_github("glmgen/genlasso")
```

  - `testthat` (only in dev)

## Usage

``` r
library(clggm)

### define coefficient matrices
B <- matrix(nrow = 4, c(-4, 1,   0,  0, 
                         0, -3,  1,  0,
                         0,  0, -2,  1,
                         0,  0,  0, -1), byrow = TRUE)
D <- diag(c(1,2,1,2))
C <- D %*% t(D)

### solve continuous Lyapunov equation 
### to obtain covariance of invariant 
### distribution
Sigma <- clyap(B, C) 

### obtain observations from the invariant distribution
sample <- MASS::mvrnorm(n = 10000, mu = rep(0,4),Sigma =  Sigma)
```

### Estimate \(B\) knowing \(C\)

``` r
B0 <- - 0.5 * C %*% solve(cov(sample))
### penalized maximum-likelihood solved with proximal gradient
results <- proxgradllB(cov(sample), B= B0, C = C, eps = 1e-12, 
                       maxIter = 1000, lambda = 0.05, job = 1)
results 
```

    ## $N
    ## [1] 4
    ## 
    ## $Sigma
    ##             [,1]       [,2]       [,3]        [,4]
    ## [1,] 0.148555747 0.07503349 0.01129353 0.009925968
    ## [2,] 0.075033489 0.70185934 0.06816678 0.062713744
    ## [3,] 0.011293526 0.06816678 0.56431835 0.612677504
    ## [4,] 0.009925968 0.06271374 0.61267750 1.981907470
    ## 
    ## $B
    ##          [,1]        [,2]       [,3]       [,4]
    ## [1,] -3.71938  0.70015713  0.0000000  0.0000000
    ## [2,]  0.00000 -2.88544745  0.3693624  0.0000000
    ## [3,]  0.00000  0.08280547 -1.7404734  0.7777934
    ## [4,]  0.00000  0.00000000  0.3678453 -1.1228428
    ## 
    ## $C
    ##      [,1] [,2] [,3] [,4]
    ## [1,]    1    0    0    0
    ## [2,]    0    4    0    0
    ## [3,]    0    0    1    0
    ## [4,]    0    0    0    4
    ## 
    ## $lambda
    ## [1] 0.05
    ## 
    ## $diff
    ## [1] 1.112156e-15
    ## 
    ## $objective
    ## [1] 1.397566
    ## 
    ## $iter
    ## [1] 90
    ## 
    ## $job
    ## [1] 1

``` r
### penalized least square 0.5*||S(B) - Sigma||_2^2 + lambda * ||B||_1,off
results2 <- proxgradlsB(cov(sample), B= B0, C = C, eps = 1e-12, 
                       maxIter = 1000, lambda = 0.01, job = 0)
results2$B
```

    ##           [,1]       [,2]        [,3]       [,4]
    ## [1,] -3.354636  0.5166269  0.02006317  0.0000000
    ## [2,]  0.000000 -2.8779855  0.38798324  0.1583115
    ## [3,]  0.000000  0.0000000 -1.64592888  0.8053895
    ## [4,]  0.000000  0.0000000  0.31607524 -1.1036053

#### Solutions along a regularization path

``` r
results <- llBpath(Sigma, eps = 1e-12, maxIter = 1000, job = 11)
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
    ##  [1,]    0.1    8  1  7  0  8    1
    ##  [2,]    0.2    7  0  7  0  9    0
    ##  [3,]    0.3    5  0  5  2  9    2
    ##  [4,]    0.4    5  0  5  2  9    2
    ##  [5,]    0.5    5  0  5  2  9    2
    ##  [6,]    0.6    5  0  5  2  9    2
    ##  [7,]    0.7    5  0  5  2  9    2
    ##  [8,]    0.8    5  0  5  2  9    2
    ##  [9,]    0.9    5  0  5  2  9    2
    ## [10,]    1.0    5  0  5  2  9    2

### Estimate \(C\) knowing \(B\)

``` r
results <- graddsllc(cov(sample), B, C = diag(4), C0 = diag(4), 
                      eps = 1e-12, maxIter = 1000, lambda = 0.0005)
cbind(C = diag(C), Cest = diag(results$C))
```

    ##      C      Cest
    ## [1,] 1 0.9857538
    ## [2,] 4 3.8935243
    ## [3,] 1 1.0098422
    ## [4,] 4 3.9197485

### Estimate \(B\) and \(C\)

``` r
B0 <- - 0.5 * diag(4) %*% solve(cov(sample))
results <- pnllbc(cov(sample), B0, C = diag(4), C0 = diag(4), 
                      eps = 1e-14, maxIter = 1000, intitr = 1,
                      lambda = 0.05, lambdac = 0.005, job = 0)
## B estimated
results$B
```

    ##           [,1]       [,2]        [,3]        [,4]
    ## [1,] -3.667638  0.5338953  0.00000000  0.00000000
    ## [2,]  0.000000 -0.6985771  0.04917071  0.03277296
    ## [3,]  0.000000  0.1219746 -1.19039342  0.48089501
    ## [4,]  0.000000  0.0000000  0.00000000 -0.28209690

``` r
## C estimated
cbind(C = diag(C), Cest = diag(results$C))
```

    ##      C      Cest
    ## [1,] 1 0.9735743
    ## [2,] 4 0.9493660
    ## [3,] 1 0.7175995
    ## [4,] 4 1.1249885

## Related code

  - Some code is from the `lyapunov` package
    (<https://github.com/gragusa/lyapunov>).

  - Using FORTRAN code of Algorithm 705 from the Collected Algorithms
    from ACM, TRANSACTIONS ON MATHEMATICAL SOFTWARE, VOL. 18, NO. 2, PP.
    232-238.75.
