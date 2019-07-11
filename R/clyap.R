#' Solve continuous-time Lyapunov/Sylvester equations
#' \code{clyap} solve the continuous-time Lyapunov equations
#' \deqn{AXE' + EXA'+Q=0.}
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
#' `clyap` use `BKDIS` from ALGORITHM 705, COLLECTED ALGORITHMS FROM ACM
#' 
#' @examples 
#' 
#' @useDynLib clggm
#' @export
clyap <- function(A, Q, E=NULL, all = FALSE) {
  N <- ncol(A)
  if (is.null(E)) {
    E <- diag(N)
  }
  WKV   <- rep(0, 2*N^2+3*N)
  out <- .Fortran('SYLGC', as.integer(N), as.integer(N), as.integer(N),
                  as.double(A), as.double(E), as.double(Q),
                  as.double(WKV), as.integer(0), as.integer(0),
                  PACKAGE = "clggm")
  if (all){
    out
  }else{
    matrix(out[[6]], N, N) 
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
  out <- .Fortran('GRADLLB', as.integer(N),
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

#' @rdname clyap
#' @export
fproxgradB <- function(Sigma, B, C = diag(ncol(Sigma)), eps =  1e-2,
                       alpha = 0.5, 
                       maxIter = 1000, 
                       lambda = 0, all = FALSE){
  
 out <- .Fortran("PRXGRDLLB",as.integer(ncol(Sigma)), as.double(Sigma), as.double(B), 
          as.double(C), as.double(lambda), as.double(eps),
          as.double(alpha), as.integer(maxIter), PACKAGE = "clggm")
 if  (all) return(out)
 matrix(nrow = ncol(Sigma), out[[3]])
}

