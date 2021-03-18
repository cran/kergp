\name{as.list, covTP-method}
\alias{as.list,covTP-method}

\title{Coerce a \code{covTP} Object into a List}
\usage{
\S4method{as.list}{covTP}(x)
}
\arguments{
  
  \item{x}{

    A \code{covTP} object representing a Tensor-Product covariance
    kernel.

  }

}
\value{
  
  A list with length \code{d} or \code{d + 1} where \code{d} is the
  "dimension" slot \code{x@d} of the object \code{x}. The first \code{d}
  elements of the list are one-dimensional \emph{correlation} kernel
  objects with class \code{"covTP"}. When \code{x} is a
  \emph{covariance} kernel (as opposed to a \emph{correlation} kernel),
  the list contains one more element which gives the variance.

}
\description{
  
  Coerce a \code{covTP} object representing a Tensor-Product
  covariance kernel on the \eqn{d}-dimensional Euclidean space
  into a list containing \eqn{d} one-dimensional kernels.
  
}
\section{Caution}{
  
  When \code{x} is not a correlation kernel the
  \code{(d + 1)}-th element of the returned list may be different in
  future versions: it may be a constant covariance kernel.
  
}

\examples{
set.seed(123)
d <- 6
myCov1 <- covTP(d = d, cov = "corr")
coef(myCov1) <- as.vector(simulPar(myCov1, nsim = 1))
as.list(myCov1)

## more examples and check the value of a 'covMat'
L <- list()
myCov <- list()

myCov[[1]] <- covTP(d = d, cov = "corr")
coef(myCov[[1]]) <- as.vector(simulPar(myCov[[1]], nsim = 1))
L[[1]] <- as.list(myCov[[1]])

myCov[[2]] <- covTP(k1Fun1 = k1Fun1PowExp, d = d, cov = "corr")
coef(myCov[[2]]) <- as.vector(simulPar(myCov[[2]], nsim = 1))
L[[2]] <- as.list(myCov[[2]])

myCov[[3]] <- covTP(k1Fun1 = k1Fun1PowExp, d = d, iso1 = 0L, cov = "corr")
coef(myCov[[3]]) <- as.vector(simulPar(myCov[[3]], nsim = 1))
L[[3]] <- as.list(myCov[[3]])

n <- 10
X <- matrix(runif(n * d), nrow = n,
            dimnames = list(NULL, paste("x", 1:d, sep = "")))
for (iTest in 1:3) {
   C <- covMat(L[[iTest]][[1]], X[ , 1, drop = FALSE])
   for (j in 2:d) {
      C <- C * covMat(L[[iTest]][[j]], X[ , j, drop = FALSE])
   }
   CTest <- covMat(myCov[[iTest]], X)
   print(max(abs(abs(C - CTest))))
}

}
\seealso{
  \code{\link{covTP}} and \code{\link{covTP-class}}.
}
