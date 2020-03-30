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
#' @useDynLib gclm
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
                  as.integer(0),
                  as.integer(JOB), as.integer(0),
                  PACKAGE = "gclm")
  if (all){
    list(B = matrix(out[[2]], N, N),
         X = matrix(out[[3]], N, N), 
         Q = matrix(out[[4]], N, N), 
         INFO = out[[6]] )
  }else{
    matrix(out[[3]], N, N) 
  }
}


#' Penalized-likelihood estimation of CLGGM
#' 
#' @param Sigma empirical covariance matrix
#' @param B initial B matrix
#' @param C diagonal of intial C matrix
#' @param C0 diagonal of penalization matrix
#' @param loss one of "loglik" (default) or "frobenius"
#' @param eps convergence tolerance 
#' @param alpha parameter line search 
#' @param maxIter maximum number of iterations
#' @param lambda penalization coefficient for B
#' @param lambdac penalization coefficient for C
#' @param job integer 0,1,10 or 11 
#' @return a list with the result of the optimization
#' @useDynLib gclm
#' @export
gclm <- function(Sigma, B = - 0.5 * diag(ncol(Sigma)), 
                      C = rep(1, ncol(Sigma)),
                      C0 = rep(1, ncol(Sigma)),
                      loss = "loglik",
                      eps =  1e-2,
                      alpha = 0.5, 
                      maxIter = 100,
                      lambda = 0,
                      lambdac = 0,  
                      job = 0){
  if (loss == "loglik"){
    out <- .Fortran("GCLMLL",as.integer(ncol(Sigma)), as.double(Sigma), 
                    as.double(B), 
                    as.double(C), as.double(C0), 
                    as.double(lambda), as.double(lambdac), 
                    as.double(eps),
                    as.double(alpha),
                    as.integer(maxIter),
                    as.integer(job), 
                    PACKAGE = "gclm")    
  }
  if (loss == "frobenius"){
    out <- .Fortran("GCLMLS",as.integer(ncol(Sigma)), as.double(Sigma), 
                    as.double(B), 
                    as.double(C), as.double(C0), 
                    as.double(lambda), as.double(lambdac), 
                    as.double(eps),
                    as.double(alpha),
                    as.integer(maxIter),
                    as.integer(job), 
                    PACKAGE = "gclm")
  }

  names(out) <- c("N", "Sigma", "B", "C", "C0", "lambda", "lambdac",
                  "diff", "loss", 
                   "iter", "job")
  out$Sigma <- matrix(nrow = out$N, out$Sigma)
  out$B <- matrix(nrow = out$N, out$B)
  return(out)
}