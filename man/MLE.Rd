\name{mle}
\alias{mle}

\title{
  Generic function for Maximum Likelihood estimation of a Gaussian Process
  model
}

\description{
  
  Generic function for the Maximum Likelihood estimation of a Gaussian Process
  model.

}
\usage{
mle(object, ...)
}

\arguments{
  \item{object}{

    An object representing a covariance kernel.

  }
  \item{\dots}{

    Other arguments for methods, typically including a response, a
    design, ...
  }
  
}

\value{

  An estimated model, typically a list. 
  
}

\seealso{
  
  \code{\link{mle-methods}} for examples.

}

