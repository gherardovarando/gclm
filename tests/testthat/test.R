context("gclm optimization")
library(gclm)

test_that("Optimization fast if initialized in optimum (loglik)", {
  ## generate random stable matrix
  B <- matrix(nrow = 10, ncol = 10, rnorm(100))
  diag(B) <- apply(B, MARGIN = 1, function(x) -sum(abs(x)))
  Sigma <- clyap(B, diag(10))
  res <- gclm(Sigma, B, eps = 0, lambda = 0, lambdac = -1)
  expect_equal(res$B, B)
  expect_equal(res$diff, 0)
  expect_lt(res$iter, 10)
})


test_that("Optimization fast if initialized in optimum (frobenius)", {
  ## generate random stable matrix
  B <- matrix(nrow = 10, ncol = 10, rnorm(100))
  diag(B) <- apply(B, MARGIN = 1, function(x) -sum(abs(x)))
  Sigma <- clyap(B, diag(10))
  res <- gclm(Sigma, B, eps = 1e-16, lambda = 0, lambdac = -1, 
              loss = "frobenius")
  expect_equal(res$B, B)
  expect_equal(res$diff, 0)
  expect_lt(res$iter, 10)
})
