context("covRadial")

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
n <- 20
d <- 4

## ============================================================================
## utility functions to use 'numDeriv'
## ============================================================================

covAsVec <- function(par, object) {
    coef(object) <- par
    C <- covMat(object = object, X = X, compGrad = TRUE)
    grad <- attr(C, "gradient")
    C <- as.vector(C)
    dim(grad) <- c(length(C), length(par))
    attr(C, "gradient") <- grad
    C
}
covAsVec2 <- function(xNew, object) {
    C <- covMat(object = object, X = X, Xnew = xNew, 
                deriv = TRUE)
    der <- attr(C, "der")
    C <- as.vector(C)
    dim(der) <- c(length(C), object@d)
    attr(C, "der") <- der
    C
}

## ============================================================================
## Define a covTP with the chosen one-dimensional kernel
## ============================================================================

## choose a design and parameter values
X <- array(runif(n * d), dim = c(n, d),
           dimnames = list(NULL, paste("x", 1:d, sep = "")))

kerns <- c("Cos", "Exp", "Gauss", "Matern3_2", "Matern5_2")
k1Fun1s <- paste("k1Fun1", kerns, sep = "")

## ============================================================================
## Replicate the test for each target one-dimensional kernel
## ============================================================================

for (i in 1:1) {

    for (iFun in 1:length(k1Fun1s)) {
        
        k1Fun1 <- match.fun(k1Fun1s[iFun])
        
        if (PRINT) {
            cat(sprintf("Kernel function : %s\n", k1Fun1s[iFun]))
        }
        
        myCov <- covRadial( k1Fun1 = k1Fun1, d = d, cov = "homo")
        
        theta <- as.vector(simulPar(object = myCov, n = 1L))
        
        ## use small ranges to allow small kernel values
        ## theta[1:d] <- theta[1:d] / 20

        coef(myCov) <- theta
        
        ## =====================================================================
        ## check the gradient w.r.t. parameters
        ## =====================================================================
    
        res <- covAsVec(theta, myCov)
        grad.check <- jacobian(covAsVec, theta, object = myCov)
        
        errGrad <- grad.check - attr(res, "gradient")
        
        if (PRINT) {
            cat(sprintf("    gradient:   %e\n", max(abs(errGrad))))
        } else {
            test_that(desc = "gradient, symmetric case",
                      code = expect_true(max(abs(errGrad)) < precision))
        }
        
        ## ====================================================================
        ## check the derivative w.r.t. 'x'
        ## ====================================================================
    
        nNew <- 1
        XNew <- array(runif(nNew * d), dim = c(nNew, d),
                      dimnames = list(NULL, inputNames(myCov)))
        
        res <- covAsVec2(XNew, myCov)
        der.check <- jacobian(covAsVec2, XNew, object = myCov)
        errDer <- der.check - attr(res, "der")
        
        if (PRINT) {
            cat(sprintf("    derivative: %e\n", max(abs(errDer))))
        } else {
            test_that(desc = "derivative, non-symmetric case",
                      code = expect_true(max(abs(errDer)) < precision))
        }
    }
    
}

## ============================================================================
## Check that we have the same results as DiceKriging for the Gauss
## (Square-Exponential) kernel since the 'covTP' and 'covRadial' are
## then identical
## ============================================================================

if (PRINT) {
    cat("Comparison of 'covMat' with 'DiceKriging'\n")
}

for (iFun in 3:3) {
    
    k1Fun1 <- match.fun(k1Fun1s[iFun])
    covType <- tolower(kerns[iFun])

    if (PRINT) {
        cat(sprintf("Cov. type: %s\n", covType))
    }

    myCov <- covRadial( k1Fun1 = k1Fun1, d = d, cov = "homo")
    theta <- as.vector(simulPar(object = myCov, n = 1L))
    coef(myCov) <- theta
        
    DKcov <- DiceKriging::covStruct.create(covtype = covType, d = d,
                                           known.covparam = "All",
                                           var.names = inputNames(myCov),
                                           coef.cov = coef(myCov)[1:d],
                                           coef.var = coef(myCov)[d + 1])

    K1kgp <- covMat(myCov, X = X)
    attr(K1kgp, "gradient") <- NULL
    K1DK <- DiceKriging::covMatrix(DKcov, X = X)$C
    errSym <- max(abs(K1kgp - K1DK))
    
    if (PRINT) {
        cat(sprintf("   Symmetric case:    %e\n", max(abs(errSym))))
    } else {
        test_that(desc = "covMat comparison with DiceKriging, symmetric case",
                  code = expect_true(errSym < precision))
    }

    K2kgp <- covMat(myCov, X = X, Xnew = XNew)
    K2DK <- DiceKriging::covMat1Mat2(DKcov, X1 = X, X2 = XNew) 
    errNonSym <- max(abs(K2kgp - K2DK))
        
    if (PRINT) {
        cat(sprintf("   Asymmetric case:   %e\n", max(abs(errNonSym))))
    } else {
        test_that(desc = "covMat comparison with DiceKriging, asymmetric case",
                  code = expect_true(errNonSym < precision))
    }



    
}

