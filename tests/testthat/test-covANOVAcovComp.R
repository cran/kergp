context("covANOVA, comparison with covComp")

PRINT <- FALSE

## ============================================================================
## Checks relative to the 'covMat' method for the class "TP". These
## concern: the gradient, the derivative (w.r.t. 'x').
## The code can be used to get more tests by changing the seed, 'd', 'n', ...
## and the 'k1Fun1' function
## ============================================================================

# library(kergp)
## require (numDeriv)  ## now in 'Depends'
precision <- 1e-6

set.seed(12345)
n <- 6
d <- 2

k1Fun1 <- k1Fun1Matern5_2

k1 <- covTP(d = 1, k1Fun1 = k1Fun1, cov = "homo")
k2 <- covTP(d = 1, k1Fun1 = k1Fun1, cov = "homo")
coef(k1) <- as.vector(simulPar(k1, nsim = 1))
coef(k2) <- as.vector(simulPar(k2, nsim = 1))

inputNames(k1) <- "x1"
inputNames(k2) <- "x2"

isoVec <- c(TRUE, FALSE)

for (i in 1:2){
  
  if (isoVec[i]) {
    coef(k2)["theta_1"] <- coef(k1)["theta_1"]
  }
  theta_1 <- coef(k1)["theta_1"]
  tau2_1 <- coef(k1)["sigma2"]
  theta_2 <- coef(k2)["theta_1"]
  tau2_2 <- coef(k2)["sigma2"]
  
  myEnv <- new.env()
  assign("k1", k1, envir = myEnv)
  assign("k2", k2, envir = myEnv)
  myCovComp <- covComp(~ (1.0 + k1()) * (1.0 + k2()), where = myEnv)
                     
  myCovANOVA <- covANOVA(d = 2, k1Fun1 = k1Fun1, 
                         cov = "homo", iso = isoVec[i])
  if (isoVec[i]){
    coef(myCovANOVA) <- c(theta_1, tau2_1, tau2_2, 1)
  } else {
    coef(myCovANOVA) <- c(theta_1, theta_2, tau2_1, tau2_2, 1)
  }
                      

  inputNames(myCovANOVA) <- c("x1", "x2")

  X <- matrix(runif(n * d), ncol = d,
              dimnames = list(NULL, c("x1", "x2")))

  KC <- covMat(myCovComp, X, compGrad = FALSE)
  KA <- covMat(myCovANOVA, X, compGrad = FALSE)
  KA2 <- covMat(myCovANOVA, X, Xnew = X)

  test_that(desc = sprintf("covComp and covANOVA match, iso = %s", isoVec[i]),
            code = expect_true(max(abs(KC - KA)) < precision))
  test_that(desc = sprintf("covANOVA, with / without Xnew, iso = %s", isoVec[i]),
            code = expect_true(max(abs(KA2 - KA)) < precision))

}