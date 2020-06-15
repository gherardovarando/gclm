
# gclm

This package contains the implementation of the algorithm in

Varando G, Hansen NR (2020) [Graphical continuous Lyapunov
models](https://arxiv.org/abs/2005.10483)

`gclm` contains methods to estimate a sparse parametrization of
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
## version on CRAN
install.packages("gclm")
## development version from github
devtools::install_github("gherardovarando/gclm")
```

## Usage

``` r
library(gclm)

### define coefficient matrices
B <- matrix(nrow = 4, c(-4,  2,  0,   0, 
                         0, -3,  1,   0,
                         0,  0, -2, 0.5,
                         0,  0,  0,  -1), byrow = TRUE)
C <- diag(c(1,4,1,4))

### solve continuous Lyapunov equation 
### to obtain real covariance matrix
Sigma <- clyap(B, C) 

### obtain observations 
sample <- MASS::mvrnorm(n = 1000, mu = rep(0,4), Sigma =  Sigma)


### Solve minimization

res <- gclm(cov(sample), lambda = 0.4, lambdac = 0.01)

res$B
```

    ##            [,1]       [,2]       [,3]       [,4]
    ## [1,] -0.9652156  0.2237539  0.0000000  0.0000000
    ## [2,]  0.0000000 -0.5769872  0.0000000  0.0000000
    ## [3,]  0.0000000  0.0000000 -0.9093145  0.1653510
    ## [4,]  0.0000000  0.0000000  0.0000000 -0.3042669

``` r
res$C
```

    ## [1] 0.3846540 0.8765735 0.5021123 1.3070860

The CLE can be freely multiplied by a scalar and thus the \(B,C\)
parametrization can be rescaled. As an example we can impose
\(C_{11} = 1\) as in the true \(C\) matrix, obtaining the estimators:

``` r
C1 <- res$C / res$C[1]
B1 <- res$B / res$C[1]

B1 
```

    ##           [,1]       [,2]     [,3]       [,4]
    ## [1,] -2.509308  0.5817017  0.00000  0.0000000
    ## [2,]  0.000000 -1.5000159  0.00000  0.0000000
    ## [3,]  0.000000  0.0000000 -2.36398  0.4298693
    ## [4,]  0.000000  0.0000000  0.00000 -0.7910144

``` r
C1
```

    ## [1] 1.000000 2.278862 1.305361 3.398082

#### Solutions along a regularization path

``` r
path <- gclm.path(cov(sample), lambdac = 0.01, 
                  lambdas = 10^seq(0, -3, length = 10))
t(sapply(path, function(res) c(lambda = res$lambda, 
                                npar = sum(res$B!=0),
                                fp = sum(res$B!=0 & B==0),
                                tp = sum(res$B!=0 & B!=0) ,
                                fn = sum(res$B==0 & B!=0),
                                tn = sum(res$B==0 & B==0),
                                errs = sum(res$B!=0 & B==0) + 
                                  sum(res$B==0 & B!=0))))
```

    ##            lambda npar fp tp fn tn errs
    ##  [1,] 1.000000000    4  0  4  3  9    3
    ##  [2,] 0.464158883    6  0  6  1  9    1
    ##  [3,] 0.215443469    6  0  6  1  9    1
    ##  [4,] 0.100000000    9  3  6  1  6    4
    ##  [5,] 0.046415888   10  3  7  0  6    3
    ##  [6,] 0.021544347   12  5  7  0  4    5
    ##  [7,] 0.010000000   12  5  7  0  4    5
    ##  [8,] 0.004641589   14  7  7  0  2    7
    ##  [9,] 0.002154435   15  8  7  0  1    8
    ## [10,] 0.001000000   16  9  7  0  0    9

## Related code

  - Some inspiration is from the `lyapunov` package
    (<https://github.com/gragusa/lyapunov>).
