
# gclm

`gclm` contatins methods to estimate a sparse parametrization of
covariance matrix as solution of a continuous time Lyapunov equation
(CLE):

\[ B\Sigma + \Sigma B^t + C = 0 \]

Solving the following \(\ell_1\) penalized loss minimization
problem:

\[ \arg\min L(\Sigma(B,C)) + \lambda \rho_1(B) + \lambda_C ||C - C_0||^2_F   \]

subject to \(B\) stable and \(C\) diagonal, where \(\rho_1(B)\) is the
\(\ell_1\) norm of the off-diagonal elements of \(B\) and
\(||C - C_0||^2_F\) is the squared frobenius norm of the difference
between \(C\) and a fixed diagonal matrix \(C_0\) (usually the
identity).

## Installation

``` r
#devtools::install_github("gherardovarando/gclm")
```

### Dependencies

  - `MASS`
  - `stats`
  - `testthat` (only in dev)

## Usage

``` r
library(gclm)

### define coefficient matrices
B <- matrix(nrow = 4, c(-4, 2,   0,  0, 
                         0, -3,  1,  0,
                         0,  0, -2,  0.5,
                         0,  0,  0, -1), byrow = TRUE)
D <- diag(c(1,2,1,2))
C <- D %*% t(D)

### solve continuous Lyapunov equation 
### to obtain real covariance matrix
Sigma <- clyap(B, C) 

### obtain observations 
sample <- MASS::mvrnorm(n = 100, mu = rep(0,4),Sigma =  Sigma)


### Solve minimization

res <- gclm(cov(sample), lambda = 0.5, lambdac = 0.01)

res$B
```

    ##           [,1]       [,2]       [,3]        [,4]
    ## [1,] -1.016968  0.4153772  0.0000000  0.00000000
    ## [2,]  0.000000 -0.4907475  0.0000000  0.00000000
    ## [3,]  0.000000  0.0000000 -0.9178208  0.08170168
    ## [4,]  0.000000  0.0000000  0.0000000 -0.40498108

``` r
res$C
```

    ## [1] 0.1945109 0.8876580 0.4806690 1.2627573

#### Solutions along a regularization path

``` r
path <- gclm.path(cov(sample), lambdac = 0.01)
t(sapply(path, function(res) c(lambda = res$lambda, 
                                npar = sum(res$B!=0),
                                fp = sum(res$B!=0 & B==0),
                                tp = sum(res$B!=0 & B!=0) ,
                                fn = sum(res$B==0 & B!=0),
                                tn = sum(res$B==0 & B==0),
                                errs = sum(res$B!=0 & B==0) + 
                                  sum(res$B==0 & B!=0))))
```

    ##           lambda npar fp tp fn tn errs
    ##  [1,] 0.77048220    4  0  4  3  9    3
    ##  [2,] 0.68487307    5  0  5  2  9    2
    ##  [3,] 0.59926393    5  0  5  2  9    2
    ##  [4,] 0.51365480    6  0  6  1  9    1
    ##  [5,] 0.42804567    6  0  6  1  9    1
    ##  [6,] 0.34243653    7  1  6  1  8    2
    ##  [7,] 0.25682740    8  2  6  1  7    3
    ##  [8,] 0.17121827    8  2  6  1  7    3
    ##  [9,] 0.08560913    9  3  6  1  6    4
    ## [10,] 0.00000000   16  9  7  0  0    9

## Related code

  - SOme inspiration from the `lyapunov` package
    (<https://github.com/gragusa/lyapunov>).

  - Using FORTRAN LAPACK.
