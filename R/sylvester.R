#' Link to fortran solver for Lyapunov equation
#' 
#' @useDynLib clggm
#' @keywords internal
sylvester_c <- function(A, E, QQ, cond.number = FALSE, all = FALSE) {
  dA  <- dim(A)
  dE  <- dim(E)
  dQQ <- dim(QQ)
  dd  <- c(dA,dE,dQQ)
  if((max(dd) - min(dd))!=0)
    stop('A, E, QQ must be conformable square matrices')
  if(cond.number)
    ierr = 1
  else
    ierr = 0
  N     <- NAE <- NC <- dd[1]
  WKV   <- rep(0, 2*N^2+3*N)
  RCOND <- 0
  out <- .Fortran('SYLGC', as.integer(NAE), as.integer(NC), as.integer(N),
                  as.double(A), as.double(E), as.double(QQ),
                  as.double(WKV), as.integer(RCOND), as.integer(ierr),
                  PACKAGE = "clggm")
  if (all){
    out
  }else{
    matrix(out[[6]], N, N) 
  }
}

#' @rdname sylvester_c
sylvester_cq <- function(A, E, QQ, WKV = rep(0, 2*N^2+3*N),
                         cond.number = FALSE) {
  dA  <- dim(A)
  dE  <- dim(E)
  dQQ <- dim(QQ)
  dd  <- c(dA,dE,dQQ)
  if((max(dd) - min(dd))!=0)
    stop('A, E, QQ must be conformable square matrices')
  if(cond.number)
    ierr = 1
  else
    ierr = 0
  N     <- NAE <- NC <- dd[1]
  #WKV   <- rep(0, 2*N^2+3*N)
  RCOND <- 0
  out <- .Fortran('SYLGCQ', as.integer(NAE), as.integer(NC), as.integer(N),
                  as.double(A), as.double(E), as.double(QQ),
                  as.double(WKV), 
                  PACKAGE = "clggm")
  matrix(out[[6]], N, N)
}

