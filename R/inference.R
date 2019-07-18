#' Estimate Continuous Lyapunov Gaussian GM
#'
#' Estimate from cross sectional data 
#' the coefficients and the noise term of an Ornstein-Uhlenbeck
#' process (\eqn{dX(t) = BX(t) + \sqrt{C} dW(t)}) 
#' using penalized maximum-likelihood.
#'
#' @param Sigma the empirical covariance matrix
#' @param B an initial guess for the coefficient matrix B
#' @param C the initial guess for the noise matrix C
#' @param eps stopping criteria
#' @param alpha parameter for line search
#' @param beta parameter for line search
#' @param maxIter maximum number of iterations
#' @param trace if >0 print info (1:termination message only, 
#'               2:messages at each step)
#' @param lambda penalization coefficient
#' @param r logical, if TRUE the set of non-zero parameters will be updated 
#'          at every iteration
#' @param h logical if TRUE only the non-zero entries of B will be updated
#'
#' @return the estimated B matrix (\code{estimateBLL}) or
#' the estiamted C matrix (\code{estiamteCLL}).
#' @export
proxgradB <- function(Sigma, B, C = diag(ncol(Sigma)), eps =  1e-2,
                        alpha = 0.2, beta = 0.5,
                        maxIter = 1000, trace = 0,
                        lambda = 0, r = FALSE, h = FALSE){
  p <- ncol(Sigma)
  ixd <- 0:(p - 1) * (p) + 1:p ##index of diagonal elements
  a <- Inf
  n <- 0
  E <- matrix(nrow = p, ncol = p, 0)
  if (h) ix <- (1: p^2 )[B != 0]
  else ix <- 1:(p^2)
  ixnd <- ix[! ix %in% ixd ] ## not diagonal elements
  IX <- rep(1, p * p) ###index for fortran
  IX[-ix] <- 0
  allres <- clyap(A = B, Q = C, all = TRUE)
  S <- matrix(nrow = p, data = allres[[6]])
  AA <- matrix(nrow = p, data = allres[[4]])
  EE <- matrix(nrow = p, data = allres[[5]])
  WKV <- allres[[7]]
  P <- solve(S)
  tk <- 1
  u <- rep(0, length(ix))
  f <- mll(P = P, S = Sigma)
  Cp <- matrix(nrow = p, ncol = p, 0)
  while (a > eps && n < maxIter) {
    if (r){
      ix <- (1: p^2 )[B != 0]
      IX[-ix] <- 0 
      ixnd <- ix[! ix %in% ixd ] ## not diagonal elements
    }
    n <- n + 1

    tmp <- P %*% Sigma %*% P - P
    
    u <- gradllB(AA, EE, tmp, S, WKV, IX)[ix]   
    
    if (trace > 3){
      col = rep(1, p^2)
      col[B != 0] <- 2
      plot(u, col =  col[ix])
    }


    Bold <- B

    #### Beck and Teboulle line search
    f <- mll(P, Sigma) + lambda * sum(abs(Bold[ixnd]))
    fnew <- Inf

    alph <- alpha
    while ( fnew   > f - sum(u * (B[ix] - Bold[ix]) ) +
            sum((B[ix] - Bold[ix]) ^ 2) / (2* alph)  || fnew > f) {

      B[ix] <- Bold[ix] + alph * u

      ### soft thres
      B[ixnd] <- sign(B[ixnd]) * (abs(B[ixnd]) - alph * lambda)
      B[ixnd][abs(B[ixnd]) < (alph * lambda)] <- 0

      ### Lyapunv solution
      allres <- clyap(A = B, Q = C, all = TRUE)
      S <- matrix(nrow = p, data = allres[[6]])
      AA <- matrix(nrow = p, data = allres[[4]])
      EE <- matrix(nrow = p, data = allres[[5]])
      WKV <- allres[[7]]
      if (all( (diag(AA) * diag(EE))  < 0 )){ ###checking if B is stable
        P <- solve(S)
        fnew <- mll(P, Sigma) + lambda * sum(abs(B[ixnd]))
      }else{
        fnew <- Inf
      }
      alph <- alph * beta
    }
    #ix <- ix[sort(u, decreasing = TRUE, index.return = TRUE)$ix[1:p]]
    #IX[-ix] <- 0 
    #ixnd <- ix[! ix %in% ixd ] 
    a <- (f - fnew ) / (abs(f))      
    if (trace > 1){
      message("Iteration: ", n, " ||diff||:", signif(a),
              " alpha:", alph / beta, "||B||_0:", sum(B!=0))
    }
  }
  if (trace > 0){
    message("Stop after ", n, " iterations, with ||diff||=", signif(a))
  }
  attr(B, ".proxgradB") <- list(
    param = c( eps = eps, alpha = alpha, 
               beta = beta, maxIter = maxIter, 
               lambda = lambda, r = r, h = h),
    optimres = c(finalDiff = a, iterations = n)
    )
  return(B)
}




#' @rdname proxgradB
#' @param C0 penalization matrix
#' @param t0 initial step length
#' @export
gradC <- function(Sigma, B, C = diag(ncol(Sigma)), C0 = diag(ncol(Sigma)),
                        eps =  1e-2,
                        alpha = 0.2,
                        beta = 0.5,
                        maxIter = 1000, trace = 0,
                        lambda = 0,
                        t0 = 1){
  p <- ncol(Sigma)

  a <- Inf
  n <- 0
  E <- matrix(nrow = p, ncol = p, 0)
  ix <- (1: p^2 )[B != 0]
  ixc <- (1: p^2 )[C != 0]
  allres <- clyap(A = B, Q = C, all = TRUE)

  AA <- matrix(nrow = p, data = allres[[4]])
  EE <- matrix(nrow = p, data = allres[[5]])
  WKV <- allres[[7]]
  S <- clyap2(A = AA, Q = C, E = EE,  WKV = WKV)
  P <- solve(S)
  while (a > eps && n < maxIter){
    u <- rep(0, length(ix))
    n <- n + 1

    tmp <- P %*% Sigma %*% P - P

    v <- vapply(1:length(ixc), function(i){
      E[ixc[i]] <- 1
      Cp <- E + t(E)
      D <- clyap2(A = AA, Q = Cp, E = EE,  WKV = WKV)
      E[ixc[i]] <- 0
      sum(tmp * D )
    }, FUN.VALUE = 1) - 2 * lambda * ( (C - C0)[ixc])

    Cold <- C

    ### backtracking
    f <- mll(P, Sigma) + lambda * sum((C -C0)^2)
    fnew <- Inf
    t <- t0
    while (fnew > f - t * alpha * sum(v^2)){
      C[ixc] <- Cold[ixc] + t * v
      if (all(diag(C) > 0)){
        S <- clyap2(A = AA, Q = C, E = EE,  WKV = WKV)
        P <- solve(S)
        fnew <- mll(P, Sigma) + lambda * sum((C -C0)^2)
      }else{
        fnew <- Inf
      }
      t <- t * beta
    }

    a <- (f - fnew) / abs(f)
    if (trace > 1){
      message("Iteration: ", n, " ||diff||:", signif(a))
    }
  }
  if (trace > 0){
    message("Stop after ", n, " iterations, with ||diff||=", signif(a))
  }

  return(C)
}



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





