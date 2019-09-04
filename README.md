
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
    ##            [,1]       [,2]       [,3]       [,4]
    ## [1,] 0.15601606 0.08321924 0.01558787 0.01211664
    ## [2,] 0.08321924 0.73264757 0.08683278 0.07235061
    ## [3,] 0.01558787 0.08683278 0.55492104 0.60110874
    ## [4,] 0.01211664 0.07235061 0.60110874 1.97583949
    ## 
    ## $B
    ##          [,1]       [,2]       [,3]       [,4]
    ## [1,] -3.58569  0.7140809  0.0000000  0.0000000
    ## [2,]  0.00000 -2.7795702  0.4197185  0.0000000
    ## [3,]  0.00000  0.1439541 -1.7630976  0.7750345
    ## [4,]  0.00000  0.0000000  0.3400480 -1.1156806
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
    ## [1] 3.342472e-15
    ## 
    ## $objective
    ## [1] 1.461488
    ## 
    ## $iter
    ## [1] 95
    ## 
    ## $job
    ## [1] 1

``` r
### penalized least square 0.5*||S(B) - Sigma||_2^2 + lambda * ||B||_1,off
results2 <- proxgradlsB(cov(sample), B= B0, C = C, eps = 1e-12, 
                       maxIter = 1000, lambda = 0.01, job = 0)
results2$B
```

    ##           [,1]        [,2]       [,3]          [,4]
    ## [1,] -3.190091  0.56017938  0.0000000  2.298023e-15
    ## [2,]  0.000000 -2.78017258  0.4427526  1.578973e-01
    ## [3,]  0.000000  0.07089838 -1.6362058  7.768397e-01
    ## [4,]  0.000000  0.00000000  0.3521275 -1.118297e+00

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
    ## [1,] 1 1.0290756
    ## [2,] 4 4.0162866
    ## [3,] 1 0.9962571
    ## [4,] 4 3.9106229

### Estimate \(B\) and \(C\)

``` r
B0 <- - 0.5 * diag(4) %*% solve(cov(sample))
results <- pnllbc(cov(sample), B0, C = diag(4), C0 = diag(4), 
                      eps = 1e-14, maxIter = 1000, intitr = 1,
                      lambda = 0.05, lambdac = 0.005, job = 0)
## B estimated
results$B
```

    ##           [,1]       [,2]        [,3]       [,4]
    ## [1,] -3.508704  0.5414103  0.00000000  0.0000000
    ## [2,]  0.000000 -0.6687548  0.07730917  0.0404080
    ## [3,]  0.000000  0.1426657 -1.20064232  0.4779548
    ## [4,]  0.000000  0.0000000  0.00000000 -0.2974951

``` r
## C estimated
cbind(C = diag(C), Cest = diag(results$C))
```

    ##      C      Cest
    ## [1,] 1 0.9690831
    ## [2,] 4 0.9396312
    ## [3,] 1 0.6999000
    ## [4,] 4 1.1927048

## Related code

  - Some code is from the `lyapunov` package
    (<https://github.com/gragusa/lyapunov>).

  - Using FORTRAN code of Algorithm 705 from the Collected Algorithms
    from ACM, TRANSACTIONS ON MATHEMATICAL SOFTWARE, VOL. 18, NO. 2, PP.
    232-238.75.
