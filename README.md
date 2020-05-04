
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
devtools::install_github("gherardovarando/gclm")
```

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

    ##           [,1]       [,2]       [,3]       [,4]
    ## [1,] -1.102875  0.2596904  0.0000000  0.0000000
    ## [2,]  0.000000 -0.5585923  0.0000000  0.0000000
    ## [3,]  0.000000  0.0000000 -0.9903715  0.2037577
    ## [4,]  0.000000  0.0000000  0.0000000 -0.3436066

``` r
res$C
```

    ## [1] 0.3306905 0.9122027 0.5053309 1.2051588

The CLE can be freely multiplied by a scalar and thus the \(B,C\)
parametrization can be rescaled. As an example we can impose
\(C_{11} = 1\) as in the true \(C\) matrix, obtaining the estimators:

``` r
C1 <- res$C / res$C[1]
B1 <- res$B / res$C[1]

B1 
```

    ##           [,1]       [,2]      [,3]       [,4]
    ## [1,] -3.335068  0.7852974  0.000000  0.0000000
    ## [2,]  0.000000 -1.6891695  0.000000  0.0000000
    ## [3,]  0.000000  0.0000000 -2.994859  0.6161582
    ## [4,]  0.000000  0.0000000  0.000000 -1.0390579

``` r
C1
```

    ## [1] 1.000000 2.758479 1.528108 3.644371

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
    ##  [1,] 0.86275198    4  0  4  3  9    3
    ##  [2,] 0.76689065    5  0  5  2  9    2
    ##  [3,] 0.67102932    5  0  5  2  9    2
    ##  [4,] 0.57516799    6  0  6  1  9    1
    ##  [5,] 0.47930666    6  0  6  1  9    1
    ##  [6,] 0.38344533    6  0  6  1  9    1
    ##  [7,] 0.28758399    6  0  6  1  9    1
    ##  [8,] 0.19172266    8  2  6  1  7    3
    ##  [9,] 0.09586133    8  2  6  1  7    3
    ## [10,] 0.00000000   16  9  7  0  0    9

## Related code

  - Some inspiration is from the `lyapunov` package
    (<https://github.com/gragusa/lyapunov>).
