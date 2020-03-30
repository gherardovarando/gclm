#' @rdname gclm
#' @param lambdas sequence of lambda
#' @param ... additional arguments passed to \code{gclm}
#' @export
gclm.path <- function(Sigma, lambdas = NULL, 
                      B = - 0.5 * diag(ncol(Sigma)),
                      C = rep(1, ncol(Sigma)),
                      ...){
  if (is.null(lambdas)) {
    lambdas = seq(0, max(diag(Sigma)), length = 10)
  }
  results <- list()
  for (i in 1:length(lambdas)){
    results[[i]] <- gclm(Sigma, B = B, C = C,
                         lambda = lambdas[i],  ...)
    B <- results[[i]]$B
    C <- results[[i]]$C
  }
  return(results)
}
