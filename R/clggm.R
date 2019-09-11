#' Solve continuous-time Lyapunov equations
#' 
#' \code{clyap} solve the continuous-time Lyapunov equations
#' \deqn{BX + XB' + C=0}
#' Using the Bartels-Stewart algorithm with Hessenberg–Schur 
#' decomposition. Optionally the Hessenberg-Schur 
#' decomposition can be returned. 
#' 
#' @param B Square matrix
#' @param C Square matrix
#' @param Q Square matrix, the orthogonal matrix used 
#' to transfrom the original equation
#' @param all logical
#' 
#' @return The solution matrix \code{X} if \code{all = FALSE}. If 
#'         \code{all = TRUE} a list with components \code{X}, \code{B}
#'         and \code{Q}. Where \code{B} and \code{Q} are the 
#'         Hessenberg-Schur form of the original matrix \code{B} 
#'         and the orthogonal matrix that performed the transformation.  
#' 
#' @details 
#' 
#' If the matrix \code{Q} is set then the matrix \code{B} 
#' is assumed to be in upper quasi-triangular form 
#' (Hessenberg-Schur canonical form),
#' as required by LAPACK subroutine \code{DTRSYL} and \code{Q} is 
#' the orthogonal matrix associated with the Hessenberg-Schur form 
#' of \code{B}. 
#' Usually the matrix \code{Q} and the appropriate form of \code{B}
#' are obtained by a first call to \code{clyap(B, C, all = TRUE)}
#' 
#' 
#' \code{clyap} uses lapack subroutines: 
#' 
#' * \code{DGEES}
#' * \code{DTRSYL} 
#' 
#' Moreover, some additional subroutines are used from:
#' 
#' * Algorithm 705; a FORTRAN-77 software package for 
#'   solving the Sylvester matrix equation AXBT + CXDT = E
#' 
#' @examples 
#' B <- matrix(data = rnorm(9), nrow = 3)
#' ## make B negative diagonally dominant, thus stable:
#' diag(B) <- - 3 * max(B) 
#' C <- diag(runif(3))
#' X <- clyap(B, C)
#' ## check it is a solution:
#' max(abs(B %*% X + X %*% t(B) + C)) 
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
#' 
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
#' 
#' * \code{N}
#' * \code{Sigma} the covariance of the estimated CLGGM
#' * \code{B} the estimated B matrix
#' * \code{C}
#' * \code{lambda} 
#' * \code{diff} the value of the last relative decrease
#' * \code{objective} the value of the objective function
#' * \code{iter} number of iterations
#' * \code{job} 
#'         
#' @export
proxgradllB <- function(Sigma, B, C = diag(ncol(Sigma)), eps =  1e-2,
                       alpha = 0.5, 
                       maxIter = 1000, 
                       lambda = 0, job = 0){
  
 out <- .Fortran("PRXGRDLLB",as.integer(ncol(Sigma)), as.double(Sigma), 
                 as.double(B), 
          as.double(C), as.double(lambda), as.double(eps),
          as.double(alpha), as.integer(maxIter),
          as.integer(job), 
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
#' \deqn{\hat{B} = \arg \min_B ||Sigma - S(B,C)||_2^2 + \lambda||B||_1}
#' 
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
#' 
#' * \code{N}
#' * \code{Sigma} the covariance of the estimated CLGGM
#' * \code{B} the estimated B matrix
#' * \code{C}
#' * \code{lambda} 
#' * \code{diff} the value of the last relative decrease
#' * \code{objective} the value of the objective function
#' * \code{iter} number of iterations
#' * \code{job} 
#' 
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
#' 
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
#' 
#' * \code{N}
#' * \code{Sigma} the covariance of the estimated CLGGM
#' * \code{B} the estimated B matrix
#' * \code{C}
#' * \code{lambda} 
#' * \code{diff} the value of the last relative decrease
#' * \code{objective} the value of the objective function
#' * \code{iter} number of iterations
#' * \code{job} 
#' 
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
#' 
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
#' 
#' * \code{N}
#' * \code{Sigma} the covariance of the estimated CLGGM
#' * \code{B} 
#' * \code{C} the estimated C matrix
#' * \code{C0}
#' * \code{lambda} 
#' * \code{diff} the value of the last relative decrease
#' * \code{objective} the value of the objective function
#' * \code{iter} number of iterations
#' * \code{job} 
#' 
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
  names(out) <- c("N", "Sigma", "B", "C", "C0", "lambda", "diff", "objective", 
                  "beta", "iter", "job")
  out$Sigma <- matrix(nrow = out$N, out$Sigma)
  out$B <- matrix(nrow = out$N, out$B)
  out$C <- diag(out$C)
  out$C0 <- diag(out$C0)
  return(out)
}


#' Penalized-likelihood estimation of CLGGM
#' 
#' @param Sigma empirical covariance matrix
#' @param B initial B matrix
#' @param C intial C matrix
#' @param C0 penalization matrix
#' @param eps convergence tolerance 
#' @param alpha parameter line search 
#' @param beta parameter line search
#' @param maxIter maximum number of iterations
#' @param intitr number of internal iterations
#' @param lambda penalization coefficient for B
#' @param lambdac penalization coefficient for C
#' @param job integer 0,1,10 or 11 
#' @return a list with the result of the optimization
#' @export
pnllbc <- function(Sigma, B, C = diag(ncol(Sigma)),
                      C0 = diag(ncol(Sigma)),
                      eps =  1e-2,
                      alpha = 0.5, 
                      beta = 0.2,
                      maxIter = 100,
                      intitr = 100,
                      lambda = 0,
                      lambdac = 0,  job = 0){
  out <- .Fortran("PNLLBC",as.integer(ncol(Sigma)), as.double(Sigma), 
                  as.double(B), 
                  as.double(diag(C)), as.double(diag(C0)), 
                  as.double(lambda), as.double(lambdac), 
                  as.double(eps),
                  as.double(alpha), as.double(beta),
                  as.integer(maxIter),as.integer(intitr), 
                  as.integer(job), 
                  PACKAGE = "clggm")
  names(out) <- c("N", "Sigma", "B", "C", "C0", "lambda", "lambdac",
                  "diff", "objective", 
                  "beta", "iter", "intiter", "job")
  out$Sigma <- matrix(nrow = out$N, out$Sigma)
  out$B <- matrix(nrow = out$N, out$B)
  out$C <- diag(out$C)
  out$C0 <- diag(out$C0)
  return(out)
}
