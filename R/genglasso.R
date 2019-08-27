#' generalized graphical lasso
#'
#' Solve the generalized graphical lasso problem with proximal gradient
#' \deqn{\hat{P} = \arg \min _ P -LL(P) + \lambda||GB||_1}
#' 
#' @param Sigma the covariance matriz
#' @param P initial precision matrix
#' @param G Matrix in the generalized penalization
#' @param lambda penalization coefficient
#' @param maxIter maximum number of iterations
#' @param eps threshold for stopping criteria
#' @param alpha param for line search
#' @param beta param for line search
#' @param trace integer, if >0 print info
#' @return The estimated precision matrix
#' @export
#' @importFrom genlasso genlasso
genGlasso <- function(Sigma, P = diag(nrow(Sigma)),
                      G = -diag(nrow(Sigma)),
                      lambda = 0.1, maxIter = 1000,
                      eps = 1e-5, alpha = 0.5, beta = 0.5, trace = 0){
  
  p <- nrow(Sigma)
  a <- Inf
  ixd <- 0:(p - 1) * (p) + 1:p ##index of diagonal elements
  ixl <- matrix(nrow = p, ncol = p, 1:p^2)[lower.tri(Sigma)]
  ixu <- matrix(nrow = p, ncol = p, 1:p^2)[upper.tri(Sigma)]
  GG <- diag(p) %x% G ## kronecker product
  GG <- GG[-ixd, ] ## not penalty on diagonal of B
  GG[, ixl] <- GG[, ixu] + GG[ , ixl]
  GG <- GG[, -ixu] ## eliminate strict upper part
  n <- 0
  while (abs(a) > eps && n < maxIter){
    n <- n + 1
    Pold <- P
    S <- solve(P)
    u <- S - Sigma
    ### line search
    f <- mll(P, Sigma) + sum(abs(G %*% P))
    fnew <- Inf
    alph <- alpha
    while ( (fnew   > f - sum(u * (P - Pold) ) +
             sum((P - Pold) ^ 2) / (2* alph) ) ){
      
      P <- Pold + alph * u
      
      ### proximal step (genLasso)
      b <- genlasso::genlasso(P[-ixu], 
                              D = GG, 
                              minlam = alph * lambda)$beta
      P[-ixu] <- b[, ncol(b)]
      P[ixu] <- P[ixl]
      if (all(eigen(P, symmetric = TRUE, only.values = TRUE)$values > 0)){
        fnew <- mll(P, Sigma) + sum(abs(G %*% P))
      }else{
        fnew <- Inf
      }
      
      alph <- alph * beta
    }
    
    a <- sqrt(sum((Pold - P)^2))
    if (trace > 1){
      message("Iteration: ", n, " ||diff||:", signif(a))
    }
  }
  if (trace > 0){
    message("Stop after ", n, " iterations, with ||diff||=", signif(a))
  }
  
  return(P)
}