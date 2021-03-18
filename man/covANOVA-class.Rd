\name{covANOVA-class}
\Rdversion{1.1}
\docType{class}
\alias{covANOVA-class}
\alias{coef,covANOVA-method}
\alias{coef<-,covANOVA,numeric-method}
\alias{coefLower,covANOVA-method}
\alias{coefLower<-,covANOVA-method}
\alias{coefUpper,covANOVA-method}
\alias{coefUpper<-,covANOVA-method}
\alias{covMat,covANOVA-method}
\alias{npar,covANOVA-method}
\alias{scores,covANOVA-method}
\alias{show,covANOVA-method}
\alias{varVec,covANOVA-method}

\title{Class \code{"covANOVA"}}
\description{

   S4 class representing a Tensor Product (ANOVA) covariance kernel.
  
}
\section{Objects from the Class}{
  
  Objects can be created by calls of the form \code{new("covANOVA", ...)}
  or by using the \code{\link{covANOVA}} function.
  
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

      Object of class \code{"integer"}. The value \code{1L} corresponds
      to a general covariance kernel. The value of \code{0L} sets the variance parameter to \code{1}, which does \emph{not} correspond to a correlation kernel. See Section 'details' of \code{\link{covANOVA}}.}
      
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

      \code{signature(object = "covANOVA")}: Get the vector of values for
      the parameters.

    }
    \item{coef<-}{

      \code{signature(object = "covANOVA", value = "numeric")}: Set the
      vector of values for the parameters.
 
    }
    \item{coefLower}{

      \code{signature(object = "covANOVA")}: Get the vector of
      lower bounds on the parameters.

    }
    \item{coefLower<-}{

      \code{signature(object = "covANOVA")}: Set the vector of lower
      bounds on the parameters.

    }
    \item{coefUpper}{

      \code{signature(object = "covANOVA")}: Get the vector of upper
      bounds on the parameters.
      
    }
    \item{coefUpper<-}{

      \code{signature(object = "covANOVA")}: Set the vector of upper
      bounds on the parameters. 

    }
    \item{covMat}{

      \code{signature(object = "covANOVA")}: Compute the covariance
      matrix for given sites.

    }
    \item{npar}{

      \code{signature(object = "covANOVA")}: Get the number of
      parameters.

    }
    \item{scores}{

      \code{signature(object = "covANOVA")}: Compute the scores
      i.e. the derivatives w.r.t. the parameters of the contribution of
      the covariance in the log-likelihood of a \code{gp}.

    }
    \item{show}{

      \code{signature(object = "covANOVA")}: Print or show the object.

    }
    \item{varVec}{

      \code{signature(object = "covANOVA")}: Compute the variance
      vector for given sites.

    }
  }
}
%% \references{}
%% \author{ }
%% \note{ }
\seealso{
  
  \code{\link{covTP}}.

}
\examples{
showClass("covANOVA")
}
\keyword{classes}
