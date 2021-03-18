\name{optimMethods}
\alias{optimMethods}
\title{Optimization Methods (or Algorithms) for the \code{mle}
  Method}
\usage{
optimMethods(optimMethod = NULL,
             optimFun = c("both", "nloptr::nloptr", "stats::optim"))
}
\arguments{
  \item{optimMethod}{
    
    A character string used to find a method in a possible approximated
    fashion, see \bold{Examples}.

  }
  \item{optimFun}{

    Value of the corresponding formal argument of the \code{mle} method,
    or \code{"both"}. In the later case the full list of algorithms will
    be obtained.

  }  
}
\value{
  
  A data frame with four character columns: \code{optimFun},
  \code{optimMethod}, \code{globLoc} and \code{derNo}. The column
  \code{globLoc} indicate whether the method is global (\code{"G"}) or
  local (\code{"L"}). The column \code{derNo} indicates whether the
  method uses derivatives (\code{D}) or not (\code{"N"}) or
  \emph{possibly} uses it (\code{"P"}). Only methods corresponding the
  \code{optimFun = "stats::optim"} can have the value \code{"P"} for
  \code{derNo}. The data frame can be zero-row if \code{optimMethod} is
  given and no method match.

}
\description{
  
  Optimization methods (or algorithms) for the \code{mle} method.

}
\section{Caution}{
  The optimization method given in the argument
  \code{optimMethod} of the \code{mle} method should be compliant
  with the \code{compGrad} argument. Only a small number of
  possibilities have been tested, including the default values.
}
\references{
  
  See \href{https://nlopt.readthedocs.io/en/latest/}{The NLopt website}.

}
\seealso{
  
  \code{\link{mle-methods}}, \code{\link[stats]{optim}},
  \code{\link[nloptr]{nloptr}}.
  
}
\examples{
optimMethods()
optimMethods(optimMethod = "cobyla")
optimMethods(optimMethod = "nelder")
optimMethods(optimMethod = "BFGS")
optimMethods("CMAES")
}
