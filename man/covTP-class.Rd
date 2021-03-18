\name{covTP-class}
\Rdversion{1.1}
\docType{class}
\alias{covTP-class}
\alias{coef,covTP-method}
\alias{coef<-,covTP,numeric-method}
\alias{coefLower,covTP-method}
\alias{coefLower<-,covTP-method}
\alias{coefUpper,covTP-method}
\alias{coefUpper<-,covTP-method}
\alias{covMat,covTP-method}
\alias{npar,covTP-method}
\alias{scores,covTP-method}
\alias{show,covTP-method}
\alias{varVec,covTP-method}

\title{Class \code{"covTP"}}
\description{

   S4 class representing a Tensor Product (TP) covariance kernel.
  
}
\section{Objects from the Class}{
  
  Objects can be created by calls of the form \code{new("covTP", ...)}
  or by using the \code{\link{covTP}} function.
  
}
\section{Slots}{
  \describe{
    \item{\code{k1Fun1}:}{

      Object of class \code{"function"} A function of a scalar numeric
      variable.

    }
    \item{\code{k1Fun1Char}:}{

      Object of class \code{"character"} describing the function in the
      slot \code{k1Fun1}.

    }
    \item{\code{hasGrad}:}{

      Object of class \code{"logical"}. Tells if the value returned by
      the function \code{kern1Fun} has an attribute named \code{"der"}
      giving the derivative(s).

    }
    \item{\code{cov}:}{

      Object of class \code{"integer"}. The value \code{0L} corresponds
      to a correlation kernel while \code{1L} is for a covariance
      kernel.
      
    }
    \item{\code{iso}:}{

      Object of class \code{"integer"}. The value \code{1L} corresponds
      to an isotropic covariance, with all the inputs sharing the same
      range value.

    }
    \item{\code{iso1}:}{

      Object of class \code{"integer"} used only when the function in
      the slot \code{k1Fun1} depends on parameters i.e. has more than
      one formal argument. NOT IMPLEMENTED YET.  

    }
    \item{\code{label}:}{

      Object of class \code{"character"}. Short description of the
      object.

    }
    \item{\code{d}:}{

      Object of class \code{"integer"}. Dimension, i.e. number of
      inputs.

    }
    \item{\code{inputNames}:}{

      Object of class \code{"optCharacter"}. Names of the inputs.

    }
    \item{\code{parLower}:}{

      Object of class \code{"numeric"}. Numeric values for the lower
      bounds on the parameters. Can be \code{-Inf}.

    }
    \item{\code{parUpper}:}{

      Object of class \code{"numeric"}. Numeric values for the upper
      bounds on the parameters. Can be \code{Inf}.

    }
    \item{\code{par}:}{

      Object of class \code{"numeric"}. Numeric values for the
      parameters. Can be \code{NA}.


    }
    \item{\code{kern1ParN1}:}{

      Object of class \code{"integer"}. The number of parameters in
      \code{k1Fun1} (such as a shape).

    }
    \item{\code{parN1}:}{
      
      Object of class \code{"integer"}. Number of parameters of the
      function \code{kern1Fun} (such as a shape).

    }
    \item{\code{parN}:}{

      Object of class \code{"integer"}. Number of parameters for the
      object. The include: \emph{direct} parameters in the function
      \code{kern1Fun}, ranges, and variance.

    }
    \item{\code{kern1ParNames}:}{

      Object of class \code{"character"}. Names of the \emph{direct}
      parameters.

    }
    \item{\code{kernParNames}:}{

      Object of class \code{"character"}. Names of the parameters.
      
    }
  
  }
}
\section{Extends}{

  Class \code{"\linkS4class{covAll}"}, directly.

}
\section{Methods}{
  \describe{
    \item{coef}{

      \code{signature(object = "covTP")}: Get the vector of values for
      the parameters.

    }
    \item{coef<-}{

      \code{signature(object = "covTP", value = "numeric")}: Set the
      vector of values for the parameters.
 
    }
    \item{coefLower}{

      \code{signature(object = "covTP")}: Get the vector of
      lower bounds on the parameters.

    }
    \item{coefLower<-}{

      \code{signature(object = "covTP")}: Set the vector of lower
      bounds on the parameters.

    }
    \item{coefUpper}{

      \code{signature(object = "covTP")}: Get the vector of upper
      bounds on the parameters.
      
    }
    \item{coefUpper<-}{

      \code{signature(object = "covTP")}: Set the vector of upper
      bounds on the parameters. 

    }
    \item{covMat}{

      \code{signature(object = "covTP")}: Compute the covariance
      matrix for given sites.

    }
    \item{npar}{

      \code{signature(object = "covTP")}: Get the number of
      parameters.

    }
    \item{scores}{

      \code{signature(object = "covTP")}: Compute the scores
      i.e. the derivatives w.r.t. the parameters of the contribution of
      the covariance in the log-likelihood of a \code{gp}.

    }
    \item{show}{

      \code{signature(object = "covTP")}: Print or show the object.

    }
    \item{varVec}{

      \code{signature(object = "covTP")}: Compute the variance
      vector for given sites.

    }
  }
}
%% \references{}
%% \author{ }
%% \note{ }
\seealso{
  
  \code{\linkS4class{covRadial}} which is a similar covariance class and
  \code{\link{covTP}} which is intended to be the standard creator
  function for this class.

}
\examples{
showClass("covTP")
}
\keyword{classes}
