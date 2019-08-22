#' Solve continuous-time Lyapunov equations
#' 
#' \code{clyap} solve the continuous-time Lyapunov equations
#' \deqn{BX + XB' + C=0}
#' Using the Bartels-Stewart algorithm with Hessenberg–Schur 
#' decomposition.
#' 
#' @param B Square matrix
#' @param C Square matrix
#' @param Q Square matrix, the orthogonal matrix used 
#' to transfrom the original equation
#' @param all logical
#' 
#' @return The solution matrix \eqn{X}.
#' 
#' @details 
#' \code{clyap} uses lapack subroutines: 
#' \code{DGEES} and \code{DTRSYL}
#' 
#' @examples 
#' 
#' @useDynLib clggm
#' @export
clyap <- function(B, C, Q = NULL, all = FALSE) {
  N <- ncol(B)
  JOB <-  0
  if (is.null(Q)) Q <- matrix(nrow = N, ncol = N, 0)
  else JOB <- 1
  WK <- rep(0, 5 * N)
  out <- .Fortran('DGELYP', as.integer(N),
                  as.double(B), as.double(-C),
                  as.double(Q),
                  as.double(WK),
                  as.integer(JOB), as.integer(0),
                  PACKAGE = "clggm")
  if (all){
    list(B = matrix(out[[2]], N, N),
         X = matrix(out[[3]], N, N), 
         Q = matrix(out[[4]], N, N), 
         INFO = out[[6]] )
  }else{
    matrix(out[[3]], N, N) 
  }
}


#' Penalized likelihood estimation of CLGGM
#' 
#' Optimize the B matrix of a continuous Lyapunov 
#' Gaussian graphical model (CLGGM) using proximal gradient. 
#' \deqn{\hat{B} = \arg \min_B LL(B,C) + \lambda||B||_1}
#' @param Sigma the observed covariance matrix
#' @param B an initial B matrix
#' @param C the C matrix 
#' @param eps convergence threshold for the proximal gradient
#' @param alpha Beck and Tabulle line search rate
#' @param maxIter the maximum number of iterations
#' @param lambda penalization coefficient 
#' @param job integer, rules to select the entries of B to be updated:
#'            0: all entries, 1: non-zero entries at each iteration,
#'            10: non-zero entries of initial B, 
#'            11: starting from non-zero entries of initial B update at 
#'            each iteration.
#' @return a list with the output of the optimization:
#'         * \code{N}
#'         * \code{Sigma} the covariance of the estimated CLGGM
#'         * \code{B} the estimated B matrix
#'         * \code{C}
#'         * \code{lambda} 
#'         * \code{diff} the value of the last relative decrease
#'         * \code{objective} the value of the objective function
#'         * \code{iter} number of iterations
#'         * \code{job} 
#'         
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
                 "objective", "iter", "job")
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
#' @param Sigma the observed covariance matrix
#' @param B an initial B matrix
#' @param C the C matrix 
#' @param eps convergence threshold for the proximal gradient
#' @param alpha Beck and Tabulle line search rate
#' @param maxIter the maximum number of iterations
#' @param lambda penalization coefficient 
#' @param job integer, rules to select the entries of B to be updated:
#'            0: all entries, 1: non-zero entries at each iteration,
#'            10: non-zero entries of initial B, 
#'            11: starting from non-zero entries of initial B update at 
#'            each iteration.
#' @return a list with the output of the optimization:
#'         * \code{N}
#'         * \code{Sigma} the covariance of the estimated CLGGM
#'         * \code{B} the estimated B matrix
#'         * \code{C}
#'         * \code{lambda} 
#'         * \code{diff} the value of the last relative decrease
#'         * \code{objective} the value of the objective function
#'         * \code{iter} number of iterations
#'         * \code{job} 
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
                  "objective", "iter", "job")
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
#' @param Sigma the observed covariance matrix
#' @param B an initial B matrix
#' @param C the C matrix 
#' @param eps convergence threshold for the proximal gradient
#' @param alpha Beck and Tabulle line search rate
#' @param maxIter the maximum number of iterations
#' @param lambda penalization coefficient 
#' @param job integer, rules to select the entries of B to be updated:
#'            0: all entries, 1: non-zero entries at each iteration,
#'            10: non-zero entries of initial B, 
#'            11: starting from non-zero entries of initial B update at 
#'            each iteration.
#' @return a list with the output of the optimization:
#'         * \code{N}
#'         * \code{Sigma} the covariance of the estimated CLGGM
#'         * \code{B} the estimated B matrix
#'         * \code{C}
#'         * \code{lambda} 
#'         * \code{diff} the value of the last relative decrease
#'         * \code{objective} the value of the objective function
#'         * \code{iter} number of iterations
#'         * \code{job} 
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
                  "objective", "iter", "job")
  out$Sigma <- matrix(nrow = out$N, out$Sigma)
  out$B <- matrix(nrow = out$N, out$B)
  out$C <- matrix(nrow = out$N, out$C)
  return(out)
}


#' Penalized likelihood estimation of CLGGM
#' 
#' Optimize the diagonal C matrix of a continuous Lyapunov 
#' Gaussian graphical model (CLGGM) usign gradient descent. 
#' \deqn{\hat{C} = \arg \min_C LL(B,C) + \lambda||C - C_0||_2^2}
#' Subject to C diagonal positive definite. 
#' @param Sigma the observed covariance matrix
#' @param B an initial B matrix
#' @param C the C matrix 
#' @param C0 the C_0 matrix
#' @param eps convergence threshold for the proximal gradient
#' @param alpha backtracking parameter
#' @param beta backtracking parameter
#' @param maxIter the maximum number of iterations
#' @param lambda penalization coefficient 
#' @param job integer
#' @return a list with the output of the optimization:
#'         * \code{N}
#'         * \code{Sigma} the covariance of the estimated CLGGM
#'         * \code{B} 
#'         * \code{C} the estimated C matrix
#'         * \code{C0}
#'         * \code{lambda} 
#'         * \code{diff} the value of the last relative decrease
#'         * \code{objective} the value of the objective function
#'         * \code{iter} number of iterations
#'         * \code{job} 
#' @export
graddsllc <- function(Sigma, B, C = diag(ncol(Sigma)),
                      C0 = diag(ncol(Sigma)),
                      eps =  1e-2,
                      alpha = 0.5, 
                      beta = 0.2,
                      maxIter = 1000, 
                      lambda = 0, job = 0){
  out <- .Fortran("GRDDSLLC",as.integer(ncol(Sigma)), as.double(Sigma), as.double(B), 
                  as.double(diag(C)), as.double(diag(C0)), as.double(lambda), 
                  as.double(eps),
                  as.double(alpha), as.double(beta),
                  as.integer(maxIter),as.integer(job),
                  PACKAGE = "clggm")
  names(out) <- c("N", "Sigma", "B", "C", "C0", "lambda", "diff", "beta", 
                  "objective", "iter", "job")
  out$Sigma <- matrix(nrow = out$N, out$Sigma)
  out$B <- matrix(nrow = out$N, out$B)
  out$C <- diag(out$C)
  out$C0 <- diag(out$C0)
  return(out)
}