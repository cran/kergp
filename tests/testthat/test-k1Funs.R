context("k1Fun1*")

## require (numDeriv)  ## now in 'Imports'

## ============================================================================
## Note that using a rather limited precision is needed because of the
## 2nd order derivative for the Matern 5/2 kernel. Remind that the
## 'k1Fun1*' functions must not be used with negative arguments!
## ============================================================================

precision <- 1e-5

k1Fun1s <- c("Cos", "Exp", "Gauss", "Matern3_2", "Matern5_2")

k1Fun1s <- paste("k1Fun1", k1Fun1s, sep = "")

## Caution x close to zero will cause an error for the Exp and
## Matern3_2 kernels which are not C^2
x <- seq(from = -3, to = 3, length.out = 100)

for (fn in k1Fun1s) {
    f <- match.fun(fn)
    kappa <- f(x)
    J <- H <- errJ <- errH <- rep(NA, length(x))
    for (i in seq_along(x)) {
        J[i] <- jacobian(f, x[i])
        H[i] <- hessian(f, x[i])
    }
    errJ <- J - attr(kappa, "der")[, 1]
    errH <- H - attr(kappa, "der")[, 2]

    if (FALSE) {
        cat("o function  ", fn, "\n")
        cat("  gradient: ", max(abs(errJ)), "\n")
        cat("  Hessian:  ", max(abs(errH)), "\n\n")
    }
    
    test_that(desc = sprintf("Function %s, Jacobian", fn),
              code = expect_true(max(abs(errJ)) < precision))

    test_that(desc = sprintf("Function %s, Hessian", fn),
              code = expect_true(max(abs(errH)) < precision))
              
}
    
