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
    #tmp <- Sigma - S
    # delta <- Sigma - S
    # Cp <- delta %*% (B + C %*% P)
    # Cp <- 0.5 * (Cp + t(Cp))
    # AA <- S %x% diag(p)
    # ixl <- matrix(nrow = p, ncol = p, 1:p^2)[lower.tri(Sigma)]
    # ixu <- matrix(nrow = p, ncol = p, 1:p^2)[upper.tri(Sigma)]
    # AA[, ixl] <-   AA[ , ixl] - AA[, ixu]
    # AA <- AA[, -c(ixu, ixd)]
    # u <- lm.fit(AA, c(Cp %*% S), offset = FALSE)$residual
    # 
    # u <- u[ix]
    # u <- u / sqrt(sum(u ^ 2))
    
    #Binv <- solve(B)
    #u1 <- vapply(ix, function(i){
    #  E[ix[i]] <- 1
    #  Cp <- E %*% S + S %*% t(E)
    #  return(- 0.5 * sum(Binv * (Cp %*% tmp)))
    #}, FUN.VALUE = 1)
    
    u <- gradllB(AA, EE, tmp, S, WKV, IX)[ix]   

    #Jac <- jacllB(AA, EE, S, WKV)
    #u1 <- Jac %*% c(tmp)
    #plot(u1, u)
    #u <- u / sqrt(sum(u ^ 2))
    
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
#' @export
aproxgradB <- function(Sigma, B, C = diag(ncol(Sigma)), eps =  1e-2,
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
  f <- mll(P = P, S = Sigma) + lambda * sum(abs(B[ixnd]))
  Cp <- matrix(nrow = p, ncol = p, 0)
  while (a > eps && n < maxIter) {
    if (r){
      ix <- (1: p^2 )[B != 0]
      IX[-ix] <- 0 
      ixnd <- ix[! ix %in% ixd ] ## not diagonal elements
    }
    n <- n + 1
    jac <- jacllB(AA, EE, S, WKV) 
    fold <- f
    alph <- 1
    k <- 0
    while (alph > 1e-8 && k < 100){
      tmp <- P %*% Sigma %*% P - P
      u <- jac %*% c(tmp)
      k <- k + 1
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
    }
    
    
    a <- (fold - fnew ) / (abs(fold))
    if (trace > 1){
      message("Iteration: ", n, " ||diff||:", signif(a),
              " alpha:", alph / beta)
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




#' generalized graphical lasso
#'
#' Solve the generalized graphical lasso problem with proximal gradient
#' 
#' @param Sigma the covariance matriz
#' @param P initial precision matrix
#' @param G Matrix in the generlaized penalization
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
      b <- genlasso(P[-ixu], D = GG, minlam = alph * lambda)$beta
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


