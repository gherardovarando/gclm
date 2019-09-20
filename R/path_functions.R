
#' Path of B estimates (log-lik)
#' 
#' @param Sigma empirical covariance matrix
#' @param lambdas increasing sequence of lambdas
#' @param C the known C = DD' 
#' @param B0 initial B matrix
#' @param eps convergence criteria
#' @param maxIter maximum iterations for each proximal gradient
#' @param job integer flag (0,1,10,11) as in \code{\link{proxgradllB}}
#' @export
llBpath <- function(Sigma, lambdas = NULL, 
                    C = diag(nrow(Sigma)),
                    B0 = NULL,
                    eps = 1e-8, maxIter = 1000, 
                    job = 11){
  if (is.null(lambdas)) {
    lambdas = seq(0, max(diag(Sigma)), length = 10)
  }
  results <- list()
  if (is.null(B0)){
    B0 <- - 0.5 * C %*% solve(Sigma)    
  }
  for (i in 1:length(lambdas)){
    results[[i]] <- proxgradllB(Sigma, B0, C, eps, alpha = 0.5, 
                                maxIter = maxIter, lambda = lambdas[i], 
                                job = job)
    B0 <- results[[i]]$B
    
  }
  return(results)
} 




#' Path of B estimates (frobenius)
#' 
#' @param Sigma empirical covariance matrix
#' @param lambdas increasing sequence of lambdas
#' @param C the known C = DD' 
#' @param B0 the initial B matrix
#' @param eps convergence criteria
#' @param maxIter maximum iterations for each proximal gradient
#' @param job integer flag (0,1,10,11)
#' @export
lsBpath <- function(Sigma, lambdas = NULL, 
                    C = diag(nrow(Sigma)), B0 = NULL,
                    eps = 1e-8, maxIter = 1000, 
                    job = 10){
  if (is.null(lambdas)) {
    lambdas = seq(max(abs(Sigma))/20, max(abs(Sigma)) / 2, length = 10)
  }
  results <- list()
  if (is.null(B0)){
    B0 <- - 0.5 * C %*% solve(Sigma)    
  }
  for (i in 1:length(lambdas)){
    results[[i]] <-proxgradlsB(Sigma, B0, C, eps, alpha = 0.5, 
                                maxIter = maxIter, lambda = lambdas[i], 
                                job = job)
    B0 <- results[[i]]$B
    
  }
  return(results)
} 

