context("covTP, case of a shape parameter")

PRINT <- FALSE

## ============================================================================
## Checks relative to the 'covMat' method for the class "TP". These
## concern: the gradient, the derivative (w.r.t. 'x').
## The code can be used to get more tests by changing the seed, 'd', 'n', ...
## and the 'k1Fun1' function
## ============================================================================

## require (numDeriv)  ## now in 'Depends'
precision <- 1e-6

set.seed(12345)
n <- 10
d <- 4


## ============================================================================
## Define a covTP with the chosen one-dimensional kernel
## ============================================================================

## choose a design and parameter values
X <- array(runif(n * d), dim = c(n, d),
           dimnames = list(NULL, paste("x", 1:d, sep = "")))

kerns <- c("PowExp")
k1Fun1s <- paste("k1Fun1", kerns, sep = "")


## ============================================================================
## Replicate the test for each target one-dimensional kernel if wanted
## ============================================================================

for (iso1 in 0:1) {

    for (iFun in 1:length(k1Fun1s)) {
        
        k1Fun1 <- match.fun(k1Fun1s[iFun])
        
        if (PRINT) {
            cat(sprintf("Kernel function : %s\n", k1Fun1s[iFun]))
        }
        
        myCov <- covTP(k1Fun1 = k1Fun1, d = d, cov = "homo", iso1 = iso1)
        coefLower(myCov)[1:d] <- 0.0

        theta <- as.vector(simulPar(object = myCov, n = 1L))
        
        ## use small ranges to allow small kernel values
        ## theta[1:d] <- theta[1:d] / 20

        coef(myCov) <- theta
        
        ## ====================================================================
        ## check the gradient w.r.t. parameters
        ## ====================================================================
    
        res <- covAsVec(theta, myCov, X)
        grad.check <- jacobian(covAsVec, theta, object = myCov, X = X)
        
        errGrad <- grad.check - attr(res, "gradient")
        
        if (PRINT) {
            cat(sprintf("    gradient:   %e\n", max(abs(errGrad))))
        } else {
            test_that(desc = sprintf("gradient, symmetric case %s iso1 = %d",
                          kerns[iFun], iso1),
                      code = expect_true(max(abs(errGrad)) < precision))
        }
        
        ## ====================================================================
        ## check the derivative w.r.t. 'x'
        ## ====================================================================
    
        nNew <- 1
        XNew <- array(runif(nNew * d), dim = c(nNew, d),
                      dimnames = list(NULL, inputNames(myCov)))
        
        resDer <- covAsVec2(XNew, myCov, X)
        der.check <- jacobian(covAsVec2, XNew, object = myCov, X = X)
        errDer <- der.check - attr(resDer, "der")
        
        if (PRINT) {
            cat(sprintf("    derivative: %e\n", max(abs(errDer))))
        } else {
            test_that(desc = sprintf(paste("derivative, non-symmetric case",
                          " %s, iso1 = %d"), kerns[iFun], iso1),
                      code = expect_true(max(abs(errDer)) < precision))
        }
    }
    
}

## ============================================================================
## Check that we have the same results as DiceKriging For now only the
## non-iso case is coped with.
## ============================================================================

if (PRINT) {
    cat("Comparison of 'covMat' with 'DiceKriging'\n")
}

for (iso1 in 0) {
    
    for (iFun in 1:length(k1Fun1s)) {
        
        k1Fun1 <- match.fun(k1Fun1s[iFun])
        covType <- tolower(kerns[iFun])
        
        if (PRINT) {
            cat(sprintf("Cov. type: %s\n", covType))
        }

        if (iso1 == 0) {
            myCov <- covTP(k1Fun1 = k1Fun1, d = d, cov = "homo", iso1 = iso1)
            coefLower(myCov)[1:d] <- 0.0
            coefUpper(myCov)[1:d] <- 2.0
            theta <- as.vector(simulPar(object = myCov, n = 1L))
            coef(myCov) <- theta
            thetaMod <- c(theta[(d + 1):(2 * d)], theta[1:d], theta[2 * d + 1])
            DKcov <-
                DiceKriging::covStruct.create(covtype = covType, d = d,
                                              known.covparam = "All",
                                              var.names = inputNames(myCov),
                                              coef.cov = thetaMod[1:(2 * d)],
                                              coef.var = thetaMod[2 * d + 1])
        } else {
            ## XXX  to be implemented later
        }
            
        K1kgp <- covMat(myCov, X = X)
        attr(K1kgp, "gradient") <- NULL
        K1DK <- DiceKriging::covMatrix(DKcov, X = X)$C
        errSym <- max(abs(K1kgp - K1DK))
        
        if (PRINT) {
            cat(sprintf("   Symmetric case:    %e\n", max(abs(errSym))))
        } else {
            test_that(desc = sprintf(paste("covMat comparison with",
                          "DiceKriging, ",
                          "symmetric case. Shape iso = %d"), iso1),
                      code = expect_true(errSym < precision))
        }
        
        K2kgp <- covMat(myCov, X = X, Xnew = XNew)
        K2DK <- DiceKriging::covMat1Mat2(DKcov, X1 = X, X2 = XNew) 
        errNonSym <- max(abs(K2kgp - K2DK))
        
        if (PRINT) {
            cat(sprintf("   Asymmetric case:   %e\n", max(abs(errNonSym))))
        } else {
            test_that(desc = sprintf(paste("covMat comparison with",
                          "DiceKriging, ",
                          "asymmetric case. Shape iso = %d"), iso1),
                      code = expect_true(errNonSym < precision))
        }
        
    }
}

