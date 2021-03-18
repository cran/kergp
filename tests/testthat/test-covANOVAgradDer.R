context("covANOVA, gradients and derivatives")

PRINT <- FALSE

## ============================================================================
## Checks relative to the 'covMat' method for the class "ANOVA". These
## concern: the gradient, the derivative (w.r.t. 'x').
## The code can be used to get more tests by changing the seed, 'd', 'n', ...
## and the 'k1Fun1' function
## ============================================================================

# require (numDeriv)  ## now in 'Depends'
precision <- 1e-6

set.seed(12345)
n <- 10
d <- 4


## ============================================================================
## Define a covANOVA with the chosen one-dimensional kernel
## ============================================================================

## choose a design and parameter values
X <- array(runif(n * d), dim = c(n, d),
           dimnames = list(NULL, paste("x", 1:d, sep = "")))

kerns <- c("Cos", "Exp", "Gauss", "Matern3_2", "Matern5_2")
k1Fun1s <- paste("k1Fun1", kerns, sep = "")

## ============================================================================
## Replicate the test for each target one-dimensional kernel if wanted
## ============================================================================

isoVec <- c(TRUE, FALSE)

for (i in 1:1) {

    for (iFun in 1:length(k1Fun1s)) {
    for (iIso in 1:2){    
        k1Fun1 <- match.fun(k1Fun1s[iFun])
        
        if (PRINT) {
            cat(sprintf("Kernel function : %s\n", k1Fun1s[iFun]))
        }
        
        myCov <- covANOVA( k1Fun1 = k1Fun1, d = d, 
                           cov = "homo", iso = isoVec[iIso])
        
        theta <- as.vector(simulPar(object = myCov, n = 1L))
        
        ## use small ranges to allow small kernel values
        ## theta[1:d] <- theta[1:d] / 20

        coef(myCov) <- theta
        
        ## =====================================================================
        ## check the gradient w.r.t. parameters
        ## =====================================================================
    
        res <- kergp:::covAsVec(theta, myCov, X)
        grad.check <- jacobian(covAsVec, theta, object = myCov, X = X)
        
        errGrad <- grad.check - attr(res, "gradient")
        
        if (PRINT) {
            cat(sprintf("    gradient:   %e\n", max(abs(errGrad))))
        } else {
            test_that(desc = sprintf("gradient, symmetric case %s, iso = %s",
                          kerns[iFun], isoVec[iIso]),
                      code = expect_true(max(abs(errGrad)) < precision))
        }
        
        ## ====================================================================
        ## check the derivative w.r.t. 'x'
        ## ====================================================================
    
        nNew <- 1
        XNew <- array(runif(nNew * d), dim = c(nNew, d),
                      dimnames = list(NULL, inputNames(myCov)))
        
        res <- covAsVec2(XNew, myCov, X)
        der.check <- jacobian(covAsVec2, XNew, object = myCov, X = X)
        errDer <- der.check - attr(res, "der")
        
        if (PRINT) {
            cat(sprintf("    derivative: %e\n", max(abs(errDer))))
        } else {
            test_that(desc = sprintf("derivative, non-symmetric case %s, 
                                     iso = %s",
                          kerns[iFun], iso = isoVec[iIso]),
                      code = expect_true(max(abs(errDer)) < precision))
        }
    }
    }
}

