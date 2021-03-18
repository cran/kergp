\name{warpNorm}
\alias{warpNorm}
\alias{warpPower}
\alias{warpUnorm}
\alias{warpSpline1}
\alias{warpSpline2}
\alias{warpFun}
%\docType{data}
\title{
Warpings for Ordinal Inputs
}
\description{
Given warpings for ordinal inputs.
}
\usage{
warpNorm
warpUnorm
warpPower
warpSpline1
warpSpline2
}
\format{
  The format is a list of 6:

 $ fun : the warping function. The second argument is the vector of
 parameters. The function returns a numeric vector with an attribute
 \code{"gradient"} giving the derivative w.r.t. the parameters.
 
 $ parNames  : names of warping parameters (character vector). 
 
 $ parDefault: default values of warping parameters (numeric vector).
 
 $ parLower  : lower bounds of warping parameters (numeric vector).
 
 $ parUpper  : upper bounds of warping parameters (numeric vector). 
 
 $ hasGrad : a boolean equal to \code{TRUE} if \code{gradient} is
 supplied as an attribute of \code{fun}.

} 
\details{
  
  See \code{\link{covOrd}} for the definition of a warping in this
  context.  At this stage, two warpings corresponding to cumulative
  density functions (cdf) are implemented:
  
  \itemize{
    
    \item Normal distribution, truncated to \eqn{[0,1]}: \deqn{F(x) = [N(x) -
      N(0)] / [N(1) - N(0)]} where \eqn{N(x) = \Phi([x - \mu] /
      \sigma)}{N(x) = Phi([x - mu] / sigma)} is the cdf of the normal
    distribution with mean \eqn{\mu}{mu} and standard deviation
    \eqn{\sigma}{sigma}.
    
    \item Power distribution on \eqn{[0, 1]}: \eqn{F(x) = x^{pow}}{F(x) = x^pow}.

  }
  
  Furthermore, a warping corresponding to unnormalized Normal cdf is implemented, 
  as well as spline warpings of degree 1 and 2. 
  Splines are defined by a sequence of \code{k} knots between 0 and 1. The first knot is 0, and the last is 1.
  A spline warping of degree 1 is a continuous piecewise linear function. 
  It is parameterized by a positive vector of length \code{k-1}, representing the increments at knots.
  A spline warping of degree 2 is a non-decreasing quadratic spline. It is obtained by integrating a spline of degree 1.
  Its parameters form a positive vector of length \code{k}, representing the derivatives at knots. The implementation relies on the function \code{scalingFun1d} of \code{DiceKriging} package.
}
