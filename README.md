
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

    ##            [,1]       [,2]       [,3]        [,4]
    ## [1,] -0.9558366  0.2183980  0.0000000  0.07412284
    ## [2,]  0.0000000 -0.6297386  0.0000000  0.18236536
    ## [3,]  0.0000000  0.0000000 -0.9202061  0.22274327
    ## [4,]  0.0000000  0.0000000  0.0000000 -0.31707165

``` r
res$C
```

    ## [1] 0.3972847 0.8313878 0.4728965 1.2315232

The CLE can be freely multiplied by a scalar and thus the \(B,C\)
parametrization can be rescaled. As an example we can impose
\(C_{11} = 1\) as in the true \(C\) matrix, obtaining the estimators:

``` r
C1 <- res$C / res$C[1]
B1 <- res$B / res$C[1]

B1 
```

    ##           [,1]       [,2]      [,3]       [,4]
    ## [1,] -2.405923  0.5497266  0.000000  0.1865736
    ## [2,]  0.000000 -1.5851065  0.000000  0.4590294
    ## [3,]  0.000000  0.0000000 -2.316238  0.5606641
    ## [4,]  0.000000  0.0000000  0.000000 -0.7980968

``` r
C1
```

    ## [1] 1.000000 2.092675 1.190321 3.099850

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
    ##  [1,] 0.9249956    4  0  4  3  9    3
    ##  [2,] 0.8222183    5  0  5  2  9    2
    ##  [3,] 0.7194410    5  0  5  2  9    2
    ##  [4,] 0.6166637    7  2  5  2  7    4
    ##  [5,] 0.5138865    8  2  6  1  7    3
    ##  [6,] 0.4111092    8  2  6  1  7    3
    ##  [7,] 0.3083319    8  2  6  1  7    3
    ##  [8,] 0.2055546    8  2  6  1  7    3
    ##  [9,] 0.1027773    9  3  6  1  6    4
    ## [10,] 0.0000000   16  9  7  0  0    9

## Related code

  - Some inspiration is from the `lyapunov` package
    (<https://github.com/gragusa/lyapunov>).
