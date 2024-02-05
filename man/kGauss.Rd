\name{kGauss}
\alias{kGauss}
\alias{kSE}
\title{
  
  Gauss (Squared-Exponential) Kernel

}
\description{

  Gauss (or squared exponential) covariance function.

}
\usage{
kGauss(d)
}
\arguments{
  \item{d}{
    
    Dimension.
  }
}
\value{

  An object of class \code{"covMan"} with default parameters: 1 for
  ranges and variance values.

}

\references{

  C.E. Rasmussen and C.K.I. Williams (2006), \emph{Gaussian Processes
  for Machine Learning}, the MIT Press,
  \doi{10.7551/mitpress/3206.001.0001}

}


\examples{
kGauss()  # default: d = 1, nu = 5/2
myGauss <- kGauss(d = 2)
coef(myGauss) <- c(range = c(2, 5), sigma2 = 0.1)
myGauss
}
