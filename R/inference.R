#' Recover lower triangular B
#'
#' Recover the only lower triangular stable matrix B such that 
#' \code{Sigma} (\eqn{\Sigma})
#' is the solution of the associated continuous Lyapunov equation:
#' \deqn{B\Sigma + \Sigma B' + C = 0}
#'
#'
#' @param Sigma acovariance matrix
#' @param P the inverse of the  covariance matrix
#' @param C symmetric positive definite matrix 
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





