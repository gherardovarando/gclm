
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
                    eps = 1e-6, maxIter = 1000, 
                    job = 0){
  if (is.null(lambdas)) {
    lambdas = seq(0, max(diag(Sigma)), length = 10)
  }
  results <- list()
  if (is.null(B0)){
    B0 <- - diag(p)  
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

#' @rdname llBpath
#' @param replicates integer 
#' @param data the dataset 
#' @export
llBstabilitypath <- function(data, replicates = 100, lambdas = NULL, 
                             C = diag(ncol(data)),
                             B0 = NULL,
                             eps = 1e-8, maxIter = 1000, 
                             job = 11){
  if (is.null(lambdas)) {
    lambdas = seq(0, 1, length = 10)
  }
  N <- nrow(data)
  p <- ncol(data)
  results <- replicate(replicates, {
    idx <- sample(N, size = floor(N / 2), replace = FALSE)
    Strain <- cov(data[idx,])
    Stest <- cov(data[-idx,])
    path <- llBpath(Strain, lambdas = lambdas, 
            C = C,
            B0 = B0,
            eps = eps, maxIter = maxIter, 
            job = job)
    loss <- sapply(path, function(res){
      mll(solve(res$Sigma), Stest)  
    })
    bidx <- which.min(loss)
    return(sign(abs(path[[bidx]]$B)))
  })
  apply(results, MARGIN = c(1,2), mean)
} 

#' @rdname lsBpath
#' @param replicates integer 
#' @param data the dataset 
#' @export
lsBstabilitypath <- function(data, replicates = 100, lambdas = NULL, 
                             C = diag(ncol(data)),
                             B0 = NULL,
                             eps = 1e-8, maxIter = 1000, 
                             job = 11){
  if (is.null(lambdas)) {
    lambdas = seq(0, 1, length = 10)
  }
  N <- nrow(data)
  p <- ncol(data)
  results <- replicate(replicates, {
    idx <- sample(N, size = floor(N / 2), replace = FALSE)
    Strain <- cov(data[idx,])
    Stest <- cov(data[-idx,])
    path <- lsBpath(Strain, lambdas = lambdas, 
                    C = C,
                    B0 = B0,
                    eps = eps, maxIter = maxIter, 
                    job = job)
    loss <- sapply(path, function(res){
      sum((res$Sigma - Stest) ^ 2)  
    })
    bidx <- which.min(loss)
    return(sign(abs(path[[bidx]]$B)))
  })
  apply(results, MARGIN = c(1,2), mean)
} 