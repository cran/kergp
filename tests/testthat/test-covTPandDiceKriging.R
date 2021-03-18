context("covTP, comparison with DiceKriging")

PRINT <- FALSE

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

kerns <- c("Exp", "Gauss", "Matern3_2", "Matern5_2")
k1Fun1s <- paste("k1Fun1", kerns, sep = "")


## ============================================================================
## Check that we have the same results as DiceKriging
## ============================================================================

if (PRINT) {
    cat("Comparison of 'covMat' with 'DiceKriging'\n")
}

for (iFun in 1:length(k1Fun1s)) {
    
    k1Fun1 <- match.fun(k1Fun1s[iFun])
    covType <- tolower(kerns[iFun])

    if (PRINT) {
        cat(sprintf("Cov. type: %s\n", covType))
    }

    myCov <- covTP( k1Fun1 = k1Fun1, d = d, cov = "homo")
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

    nNew <- 1
    XNew <- array(runif(nNew * d), dim = c(nNew, d),
                  dimnames = list(NULL, inputNames(myCov)))
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

