\name{inputNames}
\alias{inputNames}

\title{

  Generic function returning the names of the inputs of a covariance
  kernel

}

\description{
  
  Generic function returning the names of the inputs attached with a
  covariance kernel.
  
}

\usage{
inputNames(object, ...)
}

\arguments{
  \item{object}{

    A covariance kernel object.

  }
  \item{\dots}{

    Other arguments for methods.

  }
}

\value{

  A character vector with \emph{distinct} input names that are used e.g. in
  prediction.

}

%%\author{}

\note{

  The input names are usually checked to control that they match the
  colnames of a spatial design matrix.

}

\seealso{
  \code{\link{checkX}}
}
