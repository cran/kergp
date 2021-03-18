\name{covTP}
\alias{covTP}
\title{Creator for the Class \code{"covTP"}}
\usage{
covTP(k1Fun1 = k1Fun1Gauss,
      cov = c("corr", "homo"),
      iso = 0, iso1 = 1L,
      hasGrad = TRUE,
      inputs = NULL,
      d = NULL,
      parNames,
      par = NULL, parLower = NULL, parUpper = NULL,
      label = "Tensor product kernel",
      ...)
}
\arguments{
  
  \item{k1Fun1}{

    A kernel function of a \emph{scalar} numeric variable, and possibly
    of an extra "shape" parameter. This function can also return the
    first-order derivative or the two-first order derivatives as an
    attribute with name \code{"der"} and with a matrix content. When an
    extra shape parameter exists, the gradient can also be returned
    as an attribute with name \code{"gradient"}, see \bold{Examples}
    later. The name of the function can be given as a character string.
    
  }
  \item{cov}{
    
    A character string specifying the kind of covariance kernel:
    correlation kernel (\code{"corr"}) or kernel of a homoscedastic GP
    (\code{"homo"}). Partial matching is allowed.
    
  }
  \item{iso}{
    
    Integer. The value \code{1L} corresponds to an isotropic covariance,
    with all the inputs sharing the same range value.
    
  }
  \item{iso1}{
    
    Integer. This applies only when \code{k1Fun1} contains one or more
    parameters that can be called 'shape' parameters. At now, only one
    such parameter can be found in \code{k1Fun1} and consequently
    \code{iso1} must be of length one. With \code{iso1 = 0} the shape
    parameter in \code{k1Fun1} will generate \code{d} parameters in the
    \code{covTP} object with their name suffixed by the dimension. When
    \code{iso1} is \code{1} only one shape parameter will be created in
    the \code{covTP} object.
    
  }
  \item{hasGrad}{
    
    Integer or logical. Tells if the value returned by the function
    \code{k1Fun1} has an attribute named \code{"der"} giving the
    derivative(s).
    
  }
  \item{inputs}{
    
    Character. Names of the inputs.
    
  }
  \item{d}{
    
    Integer. Number of inputs.
    
  }
  \item{parNames}{
    
    Names of the parameters. By default, ranges are prefixed
    \code{"theta_"} in the non-iso case and the range is named
    \code{"theta"} in the iso case.
    
  }
  \item{par}{
    
    Numeric values for the parameters. Can be \code{NA}.

  }
  \item{parLower}{
    
    Numeric values for the lower bounds on the parameters. Can be
    \code{-Inf}.
    
  }
  \item{parUpper}{
    
    Numeric values for the upper bounds on the parameters. Can be
    \code{Inf}.
    
  }
  \item{label}{
    
    A short description of the kernel object.
    
  }
  \item{...}{
    
    Other arguments passed to the method \code{new}.
    
  }

}
\value{
An object with class \code{"covTP"}.

%% @note When \code{k1Fun1} has more than one formal argument, the
%% arguments with position \code{> 1} are assumed to be parameters of
%% the model. Examples are functions with formals \code{function(x,
%% shape = 1.0)} or \code{function(x, alpha = 2.0, beta = 3.0)},
%% corresponding to vector of parameter names \code{c("shape")} and
%% \code{c("alpha", "beta")}.
}
\description{
  
  Creator for the class \code{"covTP"}.

}
\details{

  A tensor-product kernel on the \eqn{d}-dimensional Euclidean space
  takes the form \deqn{K(\mathbf{x},\,\mathbf{x}') = \sigma^2
    \prod_{\ell = 1}^d \kappa(r_\ell)}{K(x1, x2) = \sigma^2 * prod_l
    k1(r[l])} where \eqn{\kappa(r)}{k1(r)} is a suitable correlation
  kernel for a one-dimensional input, and \eqn{r_\ell}{r[l]} is given by
  \eqn{r_\ell := [x_\ell - x'_\ell] / \theta_\ell}{r[l] = (x2[l] -
    x2[l]) / \theta[l]} for \eqn{\ell = 1}{l =1} to \eqn{d}.
  
  In this default form, the tensor-product kernel depends on \eqn{d + 1}
  parameters: the \emph{ranges} \eqn{\theta_\ell >0}{\theta[l] > 0} and
  the variance \eqn{\sigma^2}{\sigma^2}.
  
  An \emph{isotropic} form uses the same range \eqn{\theta}{\theta}
  for all inputs, i.e. sets \eqn{\theta_\ell =
    \theta}{\theta[l]=\theta} for all \eqn{\ell}{l}. This is obtained by
  using \code{iso = TRUE}.
  
  A \emph{correlation} version uses \eqn{\sigma^2 = 1}{\sigma^2 =
    1}. This is obtained by using \code{cov = "corr"}.
  
  Finally, the correlation kernel \eqn{\kappa(r)}{k1(r)} can depend on
  a "shape" parameter, e.g. have the form
  \eqn{\kappa(r;\,\alpha)}{k1(r, \alpha)}. The extra shape parameter
  \eqn{\alpha}{\alpha} will be considered then as a parameter of the
  resulting tensor-product kernel, making it possible to estimate it
  by ML along with the range(s) and the variance.

}
\examples{
\dontrun{
if (require(DiceKriging)) {
    ## a 16-points factorial design and the corresponding response
    d <- 2; n <- 16; x <- seq(from = 0.0, to = 1.0, length.out = 4)
    X <- expand.grid(x1 = x, x2 = x)
    y <- apply(X, 1, DiceKriging::branin)

    ## kriging model with matern5_2 covariance structure, constant
    ## trend. A crucial point is to set the upper bounds!
    mycov <- covTP(k1Fun1 = k1Fun1Matern5_2, d = 2, cov = "homo")
    coefUpper(mycov) <- c(2.0, 2.0, 1e10)
    mygp <- gp(y ~ 1, data = data.frame(X, y),
               cov = mycov, multistart = 100, noise = FALSE)

    nGrid <- 50; xGrid <- seq(from = 0, to = 1, length.out = nGrid)
    XGrid <- expand.grid(x1 = xGrid, x2 = xGrid)
    yGrid <- apply(XGrid, 1, DiceKriging::branin)
    pgp <- predict(mygp, XGrid)$mean

    mykm <- km(design = X, response = y)
    pkm <- predict(mykm, XGrid, "UK")$mean
    c("km" = sqrt(mean((yGrid - pkm)^2)),
      "gp" = sqrt(mean((yGrid - pgp)^2)))
    
}
}
}