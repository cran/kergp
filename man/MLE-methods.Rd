\name{mle-methods}
\docType{methods}
\alias{mle,covAll-method}
\alias{mle-methods}

\title{
  Maximum Likelihood estimation of Gaussian Process model parameters
}

\description{

  Maximum Likelihood estimation of Gaussian Process models which
  covariance structure is described in a covariance kernel object.
  
}
\usage{

\S4method{mle}{covAll}(object, 
    y, X, F = NULL, beta = NULL,
    parCovIni = coef(object),
    parCovLower = coefLower(object), 
    parCovUpper = coefUpper(object),
    noise = TRUE, varNoiseIni = var(y)/10,
    varNoiseLower = 0, varNoiseUpper = Inf,
    parFixed = NULL,
    compGrad = TRUE,
    doOptim = TRUE,
    method = "L-BFGS-B",
    control = list(fnscale = -1, trace = 3, REPORT = 1),
    parTrack = FALSE, trace  = 0, checkNames = TRUE)
}

\arguments{
  \item{object}{

    An object representing a covariance kernel.

  }
  \item{y}{

    Response vector.

  }
  \item{X}{

    Spatial (or input) design matrix.

  } 
  \item{F}{

    Trend matrix.

  }
  \item{beta}{

    Vector of trend coefficients if known/fixed.

  }
  \item{parCovIni}{

    Vector with named elements giving the initial values for the
    parameters. See \bold{Details}.

  }
  \item{parCovLower}{

    Lower bounds for the parameters.

  }  
  \item{parCovUpper}{

    Upper bounds for the parameters.

  }
  \item{noise}{

    Logical. Should a noise be added to the error term?

  }  
  \item{varNoiseIni}{

    Initial value for the noise variance.

  }  
  \item{varNoiseLower}{

    Lower bound for the noise variance.

  }  
  \item{varNoiseUpper}{

    Upper bound for the noise variance.

  }  
  \item{parFixed}{

    Vector with named elements giving the fixed parameters (if any).
    Not implemented yet.

  }
  \item{compGrad}{

    Logical: compute and use the analytical gradient in optim? This
    is only possible when \code{object} provides the analytical gradient.
    %% Must be TRUE at this stage.

  }
  \item{doOptim}{

    Logical. If \code{FALSE} no optimization is done.

  }
  \item{method}{

    To be passed to \code{\link{optim}}.

  }
  \item{control}{

    To be passed to \code{\link{optim}}. The element \code{fnscale = -1}
    states that a maximization is required and must not be changed.

  }
  \item{parTrack}{

    If \code{TRUE}, the parameter vectors used during the optimization
    are tracked and returned as a matrix.

  }
  \item{trace}{

    If \code{TRUE}, some elements are printed.

  }
  \item{checkNames}{

    if \code{TRUE} (default), check the compatibility of \code{X} with
    \code{object}, see \code{\link{checkX}}.

  }
}

\details{
  
  The vector with element with kernel parameters; \code{parIni}
  \code{parFixed}, \code{parLower}, \code{parUpper} must have named
  elements, and have name corresponding to those in the covariance given
  in \code{object}. The control tasks are not yet completed.
  
}
\value{

  A list with elements hopefully having understandable names. 
  \item{opt}{

    Optimization results if it was successful, or an error
    object otherwise.

  }
  \item{coef.kernel}{

    The vector of 'kernel' coefficients. This will include one or
    several variance parameters.

  }
  \item{coef.trend}{

    Estimate of the vector \eqn{\boldsymbol{\beta}}{\beta} of the trend
    coefficients.

  }
  \item{parTracked}{

    A matrix with rows giving the successive iterates during
    optimization if the \code{parTrack} argument was set to \code{TRUE}.

  }

}

\author{

  Y. Deville, O. Roustant

}
\note{
  
  Most checks concerning the parameter names, dimensions of provided
  objects, \dots are not implemented yet.

}

\seealso{
  
  \code{\link{gp}} for various examples.

}

\examples{
##=======================================================================
## Example. A 4-dimensional "covMan" kernel
##=======================================================================
d <- 4
myCovMan <- 
      covMan(
         kernel = function(x1, x2, par) { 
         htilde <- (x1 - x2) / par[1]
         SS2 <- sum(htilde^2)
         d2 <- exp(-SS2)
         kern <- par[2] * d2
         d1 <- 2 * kern * SS2 / par[1]            
         attr(kern, "gradient") <- c(theta = d1,  sigma2 = d2)
         return(kern)
      },
      label = "myGauss",
      hasGrad = TRUE,
      d = 4,    
      parLower = c(theta = 0.0, sigma2 = 0.0),
      parUpper = c(theta = +Inf, sigma2 = +Inf),
      parNames = c("theta", "sigma2"),
      par = c(NA, NA)
      )
kernCode <- "
       SEXP kern, dkern;
       int nprotect = 0, d;
       double SS2 = 0.0, d2, z, *rkern, *rdkern;

       d = LENGTH(x1);
       PROTECT(kern = allocVector(REALSXP, 1)); nprotect++;
       PROTECT(dkern = allocVector(REALSXP, 2)); nprotect++;
       rkern = REAL(kern);
       rdkern = REAL(dkern);

       for (int i = 0; i < d; i++) {
         z = ( REAL(x1)[i] - REAL(x2)[i] ) / REAL(par)[0];
         SS2 += z * z; 
       }

       d2 = exp(-SS2);
       rkern[0] = REAL(par)[1] * d2;
       rdkern[1] =  d2; 
       rdkern[0] =  2 * rkern[0] * SS2 / REAL(par)[0];

       SET_ATTR(kern, install(\"gradient\"), dkern);
       UNPROTECT(nprotect);
       return kern;
     "

## inline the C function into an R function: MUCH MORE EFFICIENT!!!
\dontrun{
require(inline)
kernC <- cfunction(sig = signature(x1 = "numeric", x2 = "numeric",
                                   par = "numeric"),
                    body = kernCode)
myCovMan <- covMan(kernel = kernC, hasGrad = TRUE, label = "myGauss", d = 4,
                   parNames = c("theta", "sigma2"))
}

##=======================================================================
## Example (continued). Simulate data for covMan and trend
##=======================================================================
n <- 100; 
X <- matrix(runif(n * d), nrow = n)
colnames(X) <- inputNames(myCovMan)

coef(myCovMan) <- myPar <- c(theta = 0.5, sigma2 = 2)
C <- covMat(object = myCovMan, X = X,
            compGrad = FALSE,  index = 1L)

library(MASS)
## set.seed(456)
y <- mvrnorm(mu = rep(0, n), Sigma = C)
p <- rpois(1, lambda = 4)
if (p > 0) {
  F <- matrix(runif(n * p), nrow = n, ncol = p)
  beta <- rnorm(p)
  y <- F \%*\% beta + y
} else F <- NULL
par <- parCovIni <- c(0.6, 4)

##=======================================================================
## Example (continued). ML estimation. Note the 'partrack' argument
##=======================================================================           
est <- mle(object = myCovMan,
           parCovIni = parCovIni,
           y = y, X = X, F = F,
           parCovLower = c(0.05, 0.05), parCovUpper = c(10, 100),
           parTrack = TRUE, noise = FALSE, checkNames = FALSE)


##=======================================================================
## Example (continued). Grid for graphical analysis
##=======================================================================
\dontrun{
theta.grid <- seq(from = 0.1, to = 0.7, by = 0.2)
sigma2.grid <- seq(from = 0.3, to = 6, by = 0.4)
par.grid <- expand.grid(theta = theta.grid, sigma2 = sigma2.grid)
ll <- apply(as.matrix(par.grid), 1, est$logLikFun)
llmat <- matrix(ll, nrow = length(theta.grid),
                ncol = length(sigma2.grid))
}                

##=======================================================================
## Example (continued). Explore the surface ?
##=======================================================================
\dontrun{
require(rgl)
persp3d(x = theta.grid, y = sigma2.grid, z = ll,
        xlab = "theta", ylab = "sigma2", zlab = "logLik",
        col = "SpringGreen3", alpha = 0.6)
}
##=======================================================================
## Example (continued). Draw a contour plot for the log-lik 
##                        and show iterates
##=======================================================================
\dontrun{
contour(x = theta.grid, y = sigma2.grid, z = llmat,
        col = "SpringGreen3", xlab = "theta", ylab = "sigma2",
        main = "log-likelihood contours and iterates",
        xlim = range(theta.grid, est$parTracked[ , 1], na.rm = TRUE),
        ylim = range(sigma2.grid, est$parTracked[ , 2], na.rm = TRUE))
abline(v = est$coef.kernel[1], h = est$coef.kernel[2], lty = "dotted")
niter <- est$opt$counts[1]
points(est$parTracked[1:niter-1, ],
       col = "orangered", bg = "yellow", pch = 21, lwd = 2, type = "o")
points(est$parTracked[niter, , drop = FALSE],
       col = "blue", bg = "blue", pch = 21, lwd = 2, type = "o", cex = 1.5)

text(x = est$parTracked[ , 1], y = est$parTracked[ , 2],
     labels = 0L:(NROW(est$parTracked)-1L),
     pos = 4, cex = 0.8, col = "orangered")
points(x = myPar["theta"], y = myPar["sigma2"],
       bg = "Chartreuse3", col = "ForestGreen",
       pch = 22, lwd = 2, cex = 1.4)

legend("topright", legend = c("optim", "optim (last)", "true"),
       pch = c(21, 21, 22), lwd = c(2, 2, 2), lty = c(1, 1, NA),
       col = c("orangered", "blue", "ForestGreen"),
       pt.bg = c("yellow", "blue", "Chartreuse3"))

}
}


