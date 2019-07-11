#' Difference graph
#'
#' @param B1 matrix
#' @param B2 matrix
#' @param legend logical 
#' @param ... parameters to be passed to the igraph plotting
#'
#' @return a graph with colored edges where
#'         in red are plotted the edges in \code{B1} and not in
#'         \code{B2} and in blue the edges in \code{B2} and not in
#'         \code{B1}.
#' @export
deltaGraph <- function(B1, B2, legend = FALSE,  ...){
  gr1 <- igraph::graph_from_adjacency_matrix(abs(sign(t(B1))))
  gr2 <- igraph::graph_from_adjacency_matrix(abs(sign(t(B2))))
  d1 <-  igraph::set_edge_attr(gr1 - gr2, name = "color", index = , value = "red" )
  d1 <- igraph::add_edges(d1, edges = c(t(igraph::get.edgelist(gr2 - gr1))),
                  attr = list(color = "blue"))
  ed <- c(t(igraph::get.edgelist(igraph::intersection(gr1, gr2))))
  d1 <- igraph::add_edges(d1, edges = ed,
                  attr = list(color = "black"))
  igraph::plot.igraph(d1, ...)
  if (legend){
    legend("left", legend = c("correct", "missing", "wrong"),
           col = c("black", "red", "blue"), lty = 1, lwd = 2 ,cex = 0.5)
  }
  return(invisible(d1))
}


#' Hamming distance between adjacency matrices
#'
#' @param A matrix
#' @param B matrix
#'
#' @return The hamming distance between the directed graphs induced by
#'  \code{A} and \code{B}.
#' @export
hamming <- function(A, B){
  sum( (A != 0) != (B != 0) )
}


#' Minus Log-Likelihood for Gaussian model
#'
#' @param P the inverse of the covariance matrix
#' @param S the empirical covariance
#'
#' @importFrom stats var
#' @export
mll <- function(P, S){
 -determinant(P, logarithm = TRUE)$modulus + sum(S * P)
}

#' Minus log-likelihood for the invariant distribution of OU
#'
#' @param B coefficent matrix
#' @param S covariance matrix
#' @param C noise matrix
#'
#' @export
mllB <- function(B, S, C = diag(nrow(B))){
  P <- solve(lyapunov::clyap(A = B, Q = C))
  mll(P, S)
}


#' Generate data from invariant distribution of O-U process
#'
#' Sample observations from the invariant distribution of an
#' Ornstein-Uhlenbeck process:  \deqn{dX_t = B(X_t - \mu)dt + DdW_t}
#'
#' @param n number of sample
#' @param B square matrix, coefficients of the O-U equation
#' @param D noise matrix
#' @param mean invariant mean
#'
#' @return  A list with components:
#'   -  \code{sample} the produced sample
#'   -  \code{Sigma} the covariance matrix of the invariant distribution
#'   -  \code{C} the  C matrix obtained as DD'
#'   -  \code{B} the coefficient square matrix
#'   -  \code{mean} the mean vector of the invariant distribution
#'
#'
#' @importFrom MASS mvrnorm
#'
#' @export
rOUinv <- function(n = 1, B,
                D = diag(nrow = nrow(B), ncol = ncol(B)),
                mean = rep(0, nrow(B)) ){
  C <- D %*% t(D)
  Sigma <- clyap(A = B, Q = C)
  S <- mvrnorm(n = n, Sigma = Sigma, mu = mean)
  return(list(data = S, Sigma = Sigma, C = C, B = B, mean = mean))
}




#' Generate a random skew-symmetric matrix
#'
#' @param p integer the dimension of the desired matrix
#' @param rfun a function with signature \code{function(n)} that
#' generates \code{n} numbers (e.g. \code{rnorm})
#'
#' @return a skew-symmetric matrix
#'
#' @importFrom stats rnorm
#' @export
rSkewSymm <- function(p, rfun = rnorm){
  M <- matrix(nrow = p, ncol = p, data = rfun(p * p))
  return(M - t(M) )
}

#' Generate stable Metzler matrix
#'
#' @param p integere the dimension of the matrix
#' @param d proportion of non zero entries
#' @param lower logical, should the matrix be lower triangular
#' @param rfun the random number generator for the matrix entries
#' @param rdiag the random number generator for the diagonal added error
#'
#' @return a square matrix with positive off-diagonal entries (Metzler)
#' and eigenvalues with strict negative real part (stable)
#'
#' @importFrom stats rnorm
#' @export
rStableMetzler <- function(p = 1, d = 1, lower = FALSE,
                           rfun = rnorm, rdiag = rnorm){
  if (lower) d <- d / 2
  D <- abs(rfun(p * p)) * sample(c(0,1), replace = TRUE, size = p * p,
                             prob = c(1 - d, d))
  A <- matrix(nrow = p, ncol = p, data = D)
  diag(A) <- 0
  if (lower){
    A[upper.tri(A)] <- 0
  }
  diag(A) <- -rowSums(A) - abs(rdiag(p))
  return(A)
}


#' Recover lower triangular B
#'
#' Recover the only lower triangular stable matrix B such that the
#' associate OU process admit invariant distribution with covariance equal
#' to the given matrix.
#'
#'
#' @param Sigma PSD matrix, the covariance matrix
#' @param P PSD matrix, the inverse of the covariance
#' @param C PSD matrix
#'
#' @return A stable lower triangular matrix
#'
#' @export
lowertriangB <- function(Sigma,
                         P = solve(Sigma),
                         C = diag(nrow = nrow(Sigma))){

  p <- nrow(Sigma)
  l <- listInverseBlocks(Sigma)
  blong <- anti_t(0.5 * C %*% P)[ upper.tri(Sigma)]
  blong <- blong[length(blong):1]
  b <- blong[1:(p-1)]
  w <- l[[1]] %*% b
  t <- p
  for (i in 2:length(l)){
    b <- blong[t:(t + p - i - 1)]
    t <- t + p - i
    AA <- matrix(nrow = length(b), ncol = length(w), 0)
    for (j in (i+1):p){
      for (k in 1:(i-1)){
        ix <- (k - 1) * p - (k - 1) * k / 2 + (i - k)
        AA[j - i, ix] <-  P[k, j]
      }
    }
    b <- b + tcrossprod(AA, t(w) )
    w <- c(w, l[[i]] %*% b)
  }
  W <- matrix(nrow = p, ncol = p, 0)
  W[upper.tri(W)] <- w[length(w) : 1]
  W <- anti_t(W)
  W <- W - t(W)
  B <- (W - 0.5 * C) %*% P
  B[upper.tri(B)] <- 0
  return(B)
}


anti_t <- function(m)function (m){
  p <- nrow(m)
  j <- matrix(ncol = p, nrow = p, data = 0)
  for (i in 1:p) {
    j[i, p - i + 1] <- 1
  }
  return(j %*% t(m) %*% j)
}


#' List all the inverse matrices of the leading sub-matrix of P=Sigma^{-1}
#' (INTERNAL)
#'
#' @param Sigma an invertible  matrix
#'
#' @return a list with the inverse of the leading sub-matrices of Sigma^(-1)
#' @keywords internal
listInverseBlocks <- function(Sigma){
   l <- list()
   l[[1]] <- Sigma
   for (i in 2:(nrow(Sigma) )){
     l[[i]] <-l[[i - 1]][ - 1, - 1] - (l[[i - 1]][ - 1, 1] %*%
                                              t(l[[i - 1]][1, - 1])) /
       l[[i - 1]][1, 1]
   }
   return(l[-1])
}




#' Hessian of log-likelihood for Lyapunov parametrization
#'
#' @param B the B matrix
#' @param C the C matrix
#' @param Sigma the empirical covariance matrix
#'
#' @return the Hessian matrix of the log-likelihood
#' @export
hessianll <- function(B, C = diag(nrow(B)), Sigma){
  p <- nrow(B)
  pp <- p^2
  allres <- clyap(A = B, Q = C, all = TRUE)
  S <- matrix(nrow = p, data = allres[[6]])
  AA <- matrix(nrow = p, data = allres[[4]])
  EE <- matrix(nrow = p, data = allres[[5]])
  WKV <- allres[[7]]
  P <- solve(S)
  E <- matrix (nrow = p, ncol = p, 0)
  Ds <- list()
  H <- matrix(nrow = pp, ncol = pp, data = 0)

  Delta <- P %*% (S - Sigma) %*% P
  for (i in 1:pp){
    E[i] <- 1
    Cp <- E %*% S + S %*% t(E)
    Ds[[i]] <- clyap2(A = AA, Q = Cp, E = EE, WKV = WKV)
    E[i] <- 0
  }

  for (i in 1:pp){
    E[i] <- 1
    R <- E
    E[i] <- 0
    for (j in 1:i){
      E[i] <- 1
      R <- E
      E[i] <- 0
      E[j] <- 1
      R <- R + E
      E[j] <- 0
      DD <- clyap2(A = AA, Q = R, E = EE, WKV = WKV)
      H[i, j] <- sum((Ds[[i]] %*% P) * (Ds[[j]] %*% P)) +
                sum(Delta * DD)  - 2 * sum(Delta * (Ds[[i]] %*% P %*% Ds[[j]]))
      H[j, i] <- H[i, j]
    }
  }
return(H)
}

#' @rdname hessianll
#' @export
diagHessianll <- function(B, C = diag(nrow(B)), Sigma){
  p <- nrow(B)
  pp <- p^2
  allres <- clyap(A = B, Q = C, all = TRUE)
  S <- matrix(nrow = p, data = allres[[6]])
  AA <- matrix(nrow = p, data = allres[[4]])
  EE <- matrix(nrow = p, data = allres[[5]])
  WKV <- allres[[7]]
  P <- solve(S)
  E <- matrix (nrow = p, ncol = p, 0)
  Ds <- list()
  DH <- matrix(nrow = 1, ncol = pp, data = 0)

  Delta <- P %*% (S - Sigma) %*% P
  for (i in 1:pp){
    E[i] <- 1
    Cp <- E %*% S + S %*% t(E)
    Ds[[i]] <- clyap2(A = AA, Q = Cp, E = EE, WKV = WKV)
    E[i] <- 0
  }

  for (i in 1:pp){
    E[i] <- 1
    R <- E
    E[i] <- 0
    j = i
    E[i] <- 1
    R <- E
    E[i] <- 0
    E[j] <- 1
    R <- R + E
    E[j] <- 0
    DD <- clyap2(A = AA, Q = R, E = EE, WKV = WKV)
    DH[i] <- sum((Ds[[i]] %*% P) * (Ds[[j]] %*% P)) +
        sum(Delta * DD)  - 2 * sum(Delta * (Ds[[i]] %*% P %*% Ds[[j]]))
  }
  return(DH)
}

#' Compute the ROC curve 
#' 
#' @param x an array with confidence levels between 0 and 1
#' @param y the true target values (0 and >0)
#' @param n number of points 
#' @export
ROC <- function(x, y, n = 100){
 pp <- seq(from =1, to = 0, length.out = n)
 D <- sapply(pp, FUN = function(p){
   x[x <= p] <- 0
   c(FPR = FPR(x, y), TPR = TPR(x, y))
 })
 return(t(D))
}

#'@rdname ROC
#'@export
AUROC <- function(x, y, n = 100){
 D <- ROC(x, y, n)
 temp <- 0
 for (i in 1:(n - 1)){
   ## trapezoidal quadrature rule
   temp <- temp + (D[i + 1,2] + D[i, 2]) * ( D[i + 1,1] - D[i,1]) / 2
 }
 return(temp)
}

#'@rdname ROC
#'@export
FPR <- function(x, y){
  (sum(x != 0 & y == 0)) / (sum(y == 0))
}


#'@rdname ROC
#'@export
TPR <- function(x, y){
   sum(x != 0 & y != 0) / sum(y !=0) 
}


#' Generate a naive stable matrix 
#' 
#' @param p dimension of the matrix
#' @return a stable matrix with off-diagonal entries equal to 1 and 
#' diagonal entries equal to \code{-p}
#' @export
B0 <- function(p){
  M <- matrix(nrow = p, ncol = p, 1)
  diag(M) <- -p
  return(M)
}
