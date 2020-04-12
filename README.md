
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
C <- diag(c(1,4,1,4))

### solve continuous Lyapunov equation 
### to obtain real covariance matrix
Sigma <- clyap(B, C) 

### obtain observations 
sample <- MASS::mvrnorm(n = 100, mu = rep(0,4),Sigma =  Sigma)


### Solve minimization

res <- gclm(cov(sample), lambda = 0.3, lambdac = 0.01)

res$B
```

    ##           [,1]       [,2]       [,3]        [,4]
    ## [1,] -1.120452  0.2870877  0.0000000  0.00000000
    ## [2,]  0.000000 -0.5589620  0.0000000  0.08421313
    ## [3,]  0.000000  0.0000000 -0.9806411  0.14099902
    ## [4,]  0.000000  0.0000000  0.0000000 -0.33050276

``` r
res$C
```

    ## [1] 0.2902704 0.9119084 0.5234533 1.2340809

The CLE can be freely multiplied by a scalar and thus the \(B,C\)
parametrization can be rescaled. As an example we can impose
\(C_{11} = 1\) as in the true \(C\) matrix, obtaining the estimators:

``` r
C1 <- res$C / res$C[1]
B1 <- res$B / res$C[1]

B1 
```

    ##          [,1]       [,2]      [,3]       [,4]
    ## [1,] -3.86003  0.9890355  0.000000  0.0000000
    ## [2,]  0.00000 -1.9256599  0.000000  0.2901196
    ## [3,]  0.00000  0.0000000 -3.378371  0.4857506
    ## [4,]  0.00000  0.0000000  0.000000 -1.1386031

``` r
C1
```

    ## [1] 1.000000 3.141583 1.803330 4.251488

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

    ##          lambda npar fp tp fn tn errs
    ##  [1,] 0.9131007    4  0  4  3  9    3
    ##  [2,] 0.8116450    4  0  4  3  9    3
    ##  [3,] 0.7101894    4  0  4  3  9    3
    ##  [4,] 0.6087338    6  0  6  1  9    1
    ##  [5,] 0.5072782    6  0  6  1  9    1
    ##  [6,] 0.4058225    6  0  6  1  9    1
    ##  [7,] 0.3043669    7  1  6  1  8    2
    ##  [8,] 0.2029113    7  1  6  1  8    2
    ##  [9,] 0.1014556    8  2  6  1  7    3
    ## [10,] 0.0000000   16  9  7  0  0    9

## Related code

  - Some inspiration is from the `lyapunov` package
    (<https://github.com/gragusa/lyapunov>).
