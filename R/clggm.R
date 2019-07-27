#' Solve continuous-time Lyapunov/Sylvester equations
#' \code{clyap} solve the continuous-time Lyapunov equations
#' \deqn{AXE' + EXA'+ Q=0.}
#' 
#' @param A Square matrix
#' @param Q Square matrix
#' @param E Square matrix or \code{NULL}
#' @param WKV the working vector obtained by a first call to lyapunov solver
#' @param all logical
#' 
#' @return The solution matrix \eqn{X}.
#' 
#' @details 
#' `clyap` use `BKDIS`  
#' from ALGORITHM 705, COLLECTED ALGORITHMS FROM ACM
#' 
#' @examples 
#' 
#' @useDynLib clggm
#' @export
clyap <- function(A, C, all = FALSE) {
  N <- ncol(A)
  Q <- matrix(nrow = N, ncol = N, 0)
  out <- .Fortran('DGELYP', as.integer(N),
                  as.double(A), as.double(-C),
                  as.double(Q), as.integer(0), as.integer(0),
                  PACKAGE = "clggm")
  if (all){
    out
  }else{
    matrix(out[[3]], N, N) 
  }
}




#' @rdname clyap
#' @export
clyap2 <- function(A, Q, WKV, E) {
  N <- ncol(A)
  out <- .Fortran('SYLGCQ', as.integer(N), as.integer(N), as.integer(N),
                  as.double(A), as.double(E), as.double(Q),
                  as.double(WKV),
                  PACKAGE = "clggm")
  matrix(out[[6]], N, N)
}


#' @rdname clyap
#' @export
gradllB <- function(A, E, D, S, WKV, IX = NULL){
  N <- ncol(A)
  if (is.null(IX)){
    IX <- rep(1, N * N)
  }
  grad <- matrix(nrow = N, ncol = N, 0)
  out <- .Fortran('GRADB', as.integer(N),
                  as.double(A), as.double(E), as.double(D), as.double(S),
                  as.double(WKV), as.double(grad), as.integer(IX),
                  PACKAGE = "clggm")
  return(out[[7]])
}


#' @rdname clyap
#' @export
jacllB <- function(A, E, S, WKV){
  N <- ncol(A)
  jac <- matrix(nrow = N * N, ncol = N * N, 0)
  out <- .Fortran('JACLLB', as.integer(N),
                  as.double(A), as.double(E), as.double(S),
                  as.double(WKV), as.double(jac),
                  PACKAGE = "clggm")
  return(matrix(nrow = N *N, ncol = N * N, data = out[[6]]))
}

#' Penalized likelihood estimation of CLGGM
#' 
#' Optimize the B matrix of a continuous Lyapunov 
#' Gaussian graphical model (CLGGM) using proximal gradient. 
#' \deqn{\hat{B} = \arg \min_B LL(B,C) + \lambda||B||_1}
#' @export
proxgradllB <- function(Sigma, B, C = diag(ncol(Sigma)), eps =  1e-2,
                       alpha = 0.5, 
                       maxIter = 1000, 
                       lambda = 0, job = 0){
  
 out <- .Fortran("PRXGRDLLB",as.integer(ncol(Sigma)), as.double(Sigma), as.double(B), 
          as.double(C), as.double(lambda), as.double(eps),
          as.double(alpha), as.integer(maxIter),as.integer(job),
          PACKAGE = "clggm")
 names(out) <- c("N", "Sigma", "B", "C", "lambda", "diff", 
                 "logLikl1", "iter", "job")
 out$Sigma <- matrix(nrow = out$N, out$Sigma)
 out$B <- matrix(nrow = out$N, out$B)
 out$C <- matrix(nrow = out$N, out$C)
 return(out)
}

#' Penalized likelihood estimation of CLGGM
#' 
#' Optimize the B matrix of a continuous Lyapunov 
#' Gaussian graphical model (CLGGM) using proximal gradient. 
#' \deqn{\hat{B} = \arg \min_B LL(B,C) + \lambda||B||_1}
#' @export
proxgradlsB <- function(Sigma, B, C = diag(ncol(Sigma)), eps =  1e-2,
                      alpha = 0.5, 
                      maxIter = 1000, 
                      lambda = 0, job = 0){
  
  out <- .Fortran("PRXGRDLSB",as.integer(ncol(Sigma)), as.double(Sigma), 
                  as.double(B), 
                  as.double(C), as.double(lambda), as.double(eps),
                  as.double(alpha), as.integer(maxIter),as.integer(job),
                  PACKAGE = "clggm")
  names(out) <- c("N", "Sigma", "B", "C", "lambda", "diff", 
                  "logLikl1", "iter", "job")
  out$Sigma <- matrix(nrow = out$N, out$Sigma)
  out$B <- matrix(nrow = out$N, out$B)
  out$C <- matrix(nrow = out$N, out$C)
  return(out)
}


#' Penalized likelihood estimation of CLGGM
#' 
#' Optimize the B matrix of a continuous Lyapunov 
#' Gaussian graphical model (CLGGM) using proximal coordinate descent. 
#' \deqn{\hat{B} = \arg \min_B LL(B,C) + \lambda||B||_1}
#' @export
proxcdllB <- function(Sigma, B, C = diag(ncol(Sigma)), eps =  1e-2,
                        alpha = 0.5, 
                        maxIter = 1000, 
                        lambda = 0, job = 0){
  
  out <- .Fortran("PRXCDLLB",as.integer(ncol(Sigma)), as.double(Sigma), as.double(B), 
                  as.double(C), as.double(lambda), as.double(eps),
                  as.double(alpha), as.integer(maxIter),as.integer(job),
                  PACKAGE = "clggm")
  names(out) <- c("N", "Sigma", "B", "C", "lambda", "diff", 
                  "logLikl1", "iter", "job")
  out$Sigma <- matrix(nrow = out$N, out$Sigma)
  out$B <- matrix(nrow = out$N, out$B)
  out$C <- matrix(nrow = out$N, out$C)
  return(out)
  
}