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
  P <- solve(clyap(A = B, Q = C))
  mll(P, S)
}







anti_t <- function (m){
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
