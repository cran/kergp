% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/q1CompSymm.R
\name{corLevCompSymm}
\alias{corLevCompSymm}
\title{Correlation Matrix for the Compound Symmetry Structure}
\usage{
corLevCompSymm(par, nlevels, levels, lowerSQRT = FALSE, compGrad = TRUE,
  cov = FALSE, impl = c("C", "R"))
}
\arguments{
  
  \item{par}{

    Numeric vector of length \code{1} if \code{cov} is
    \code{TRUE} or with length \code{2} else. The first element is the
    correlation coefficient and the second one (when it exists) is the
    variance.
    
  }
  \item{nlevels}{

    Number of levels.

  }
  \item{levels}{

    Character representing the levels.
    
  }
  \item{lowerSQRT}{

    Logical. When \code{TRUE} the (lower) Cholesky
    root \eqn{\mathbf{L}}{L} of the correlation matrix
    \eqn{\mathbf{C}}{C} is returned instead of the correlation matrix.
    
  }
  \item{compGrad}{

    Logical. Should the gradient be computed?

  }
  \item{cov}{Logical.

    If \code{TRUE} the matrix is a covariance
    matrix (or its Cholesky root) rather than a correlation matrix and
    the last element in \code{par} is the variance.

  }
  \item{impl}{

    A character telling which of the C and R implementations
    should be chosen.
  }
}
\value{

  A correlation matrix (or its Cholesky root) with the
  optional \code{gradient} attribute.

}
\description{

  Compute the correlation matrix for a the compound symmetry structure.

}
\note{

  When \code{lowerSQRT} is \code{FALSE}, the implementation
  used is always in R because no gain would then result from an
  implementation in C.

}
\examples{
checkGrad <- TRUE
lowerSQRT <- FALSE
nlevels <- 12
set.seed(1234)
par <- runif(1L, min = 0, max = pi)

##============================================================================
## Compare R and C implementations for 'lowerSQRT = TRUE'
##============================================================================
tR <- system.time(TR <- corLevCompSymm(nlevels = nlevels, par = par,
                                       lowerSQRT = lowerSQRT, impl = "R"))
tC <- system.time(T <- corLevCompSymm(nlevels = nlevels, par = par,
                                      lowerSQRT = lowerSQRT))
tC2 <- system.time(T2 <- corLevCompSymm(nlevels = nlevels, par = par,
                                        lowerSQRT = lowerSQRT, compGrad = FALSE))
## time
rbind(R = tR, C = tC, C2 = tC2)

## results
max(abs(T - TR))
max(abs(T2 - TR))

##===========================================================================
## Compare the gradients
##===========================================================================

if (checkGrad) {

    library(numDeriv)

    ##=======================
    ## lower SQRT case only
    ##========================
    JR <- jacobian(fun = corLevCompSymm, x = par, nlevels = nlevels,
                   lowerSQRT = lowerSQRT, impl = "R", method = "complex")
    J <- attr(T, "gradient")

    ## redim and compare.
    dim(JR) <- dim(J)
    max(abs(J - JR))
    nG <- length(JR)
    plot(1:nG, as.vector(JR), type = "p", pch = 21, col = "SpringGreen3",
         cex = 0.8, ylim = range(J, JR),
         main = paste("gradient check, lowerSQRT =", lowerSQRT))
    points(x = 1:nG, y = as.vector(J), pch = 16, cex = 0.6, col = "orangered")
}

}
\author{
Yves Deville
}
