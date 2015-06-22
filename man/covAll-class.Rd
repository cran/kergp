\name{covAll-class}
\Rdversion{1.1}
\docType{class}
\alias{covAll-class}
\alias{inputNames,covAll-method}
\alias{checkX,covAll,matrix-method}
%% \alias{npar,covAll-method}
\alias{simulPar,covAll-method}


\title{Virtual class "covAll" }
\description{

  Virtual class \code{"covAll"}, union of classes including
  \code{"covTS"}, \code{"covMan"}.

}

\section{Methods}{
  \describe{
    \item{checkX}{
      
      \code{signature(object = "covAll", X = "matrix")}: checks the compatibility of a design with a given covariance object.

    }
    \item{inputNames}{

      \code{signature(object = "covAll")}: returns the character vector of input names.

    }
    \item{simulPar}{

      \code{signature(object = "covTS")}: simulates random values for the parameters.

    }
  }
}

%% \author{Yves Deville, David Ginsbourger, Olivier Roustant}

\examples{
showClass("covAll")
}
\keyword{classes}
