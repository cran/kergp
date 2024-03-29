\name{kergp-package}
\encoding{UTF8}
\alias{kergp-package}
\alias{kergp}
\docType{package}
\title{
  
  Gaussian Process Laboratory 

}
\description{
  
  Laboratory Package for Gaussian Process interpolation, regression and
  simulation, with an emphasis on user-defined covariance kernels.

}

\note{
  
  This package was developed within the frame of the ReDice Consortium,
  gathering industrial partners (CEA, EDF, IFPEN, IRSN, Renault) and
  academic partners (Mines Saint-Étienne, INRIA, and the University of
  Bern) around advanced methods for Computer Experiments.

}

\details{
  \packageDESCRIPTION{kergp}
}

\section{Warning}{
  
  As a lab, \pkg{kergp} may strongly evolve in its future life. Users
  interested in stable software for the Analysis of Computer Experiments
  are encouraged to use other packages such as \pkg{DiceKriging}
  instead.
  
}

\author{

  Yves Deville (Alpestat), David Ginsbourger (University of Bern),
  Olivier Roustant (INSA Toulouse), with contributions from
  Nicolas Durrande (Mines Saint-Étienne).

  Maintainer: Olivier Roustant, <roustant@insa-toulouse.fr>

}
\references{

  Nicolas Durrande, David Ginsbourger, Olivier Roustant (2012).
  "Additive covariance kernels for high-dimensional gaussian process modeling".
  \emph{Annales de la Faculté des Sciences de Toulouse}, 21 (3):
  481-499.
  \href{https://afst.centre-mersenne.org/item/AFST_2012_6_21_3_481_0}{link}
  
  Nicolas Durrande, David Ginsbourger, Olivier Roustant, Laurent Carraro
  (2013).
  "ANOVA kernels and RKHS of zero mean functions for model-based sensitivity analysis".
  \emph{Journal of Multivariate Analysis}, 115, 57-67.
  \href{https://www.sciencedirect.com/science/article/pii/S0047259X1200214X}{link}
  
  David Ginsbourger, Xavier Bay, Olivier Roustant, Laurent Carraro
  (2012).
  "Argumentwise invariant kernels for the approximation of invariant functions".
  \emph{Annales de la Faculté des Sciences de Toulouse}, 21 (3):
  501-527.
  \href{https://afst.centre-mersenne.org/item/AFST_2012_6_21_3_501_0/}{link}
  
  David Ginsbourger, Nicolas Durrande, Olivier Roustant (2013).
  "Kernels and designs for modelling invariant functions: From group invariance to additivity".
  \emph{mODa 10 - Advances in Model-Oriented Design and
  Analysis. Contributions to Statistics}, 107-115.
  \href{https://link.springer.com/book/10.1007/978-3-319-00218-7}{link}
  
  Olivier Roustant, David Ginsbourger, Yves Deville (2012).
  "DiceKriging, DiceOptim: Two R Packages for the Analysis of
  Computer Experiments by Kriging-Based Metamodeling and Optimization".
  \emph{Journal of Statistical Software}, 51(1), 1-55.
  \doi{10.18637/jss.v051.i01}
}
%%\keyword{ package }
%%\seealso{}

\examples{
## ------------------------------------------------------------------
## Gaussian process modelling of function with invariance properties, 
## by using an argumentwise invariant kernel
## ------------------------------------------------------------------

## -- define manually an argumentwise invariant kernel --

kernFun <- function(x1, x2, par) {
  h <- (abs(x1) - abs(x2)) / par[1]
  S <- sum(h^2)
  d2 <- exp(-S)
  K <- par[2] * d2
  d1 <- 2 * K * S / par[1]   
  attr(K, "gradient") <- c(theta = d1,  sigma2 = d2)
  return(K)
}

## ---------------------------------------------------------------
## quicker: with Rcpp; see also an example  with package inline
## in "gp" doc. file. Note that the Rcpp "sugar" fucntions are
## vectorized, so no for loops is required.
## ---------------------------------------------------------------

\dontrun{

    cppFunction('
        NumericVector cppKernFun(NumericVector x1, NumericVector x2, 
                                 NumericVector par){
        int n1 = x1.size();
        double S, d1, d2; 
        NumericVector K(1), h(n1);
        h = (abs(x1) - abs(x2)) / par[0];  // sugar function "abs"
        S = sum(h * h);                    // sugar "*" and "sum" 
        d2 = exp(-S);
        K[0] = par[1] * d2;
        d1 = 2 * K[0] * S / par[0];   
        K.attr("gradient") = NumericVector::create(Named("theta", d1),
                                                   Named("sigma2", d2));
        return K;
     }')

}

## ---------------------------------------------------------------
## Below: with the R-based code for the kernel namely 'kernFun'.
## You can also replace 'kernFun' by 'cppKernFun' for speed.
## ---------------------------------------------------------------

covSymGauss <- covMan(kernel = kernFun,
                      hasGrad = TRUE,
                      label = "argumentwise invariant",
                      d = 2,
                      parLower = c(theta = 0.0, sigma2 = 0.0),
                      parUpper = c(theta = Inf, sigma2 = Inf),
                      parNames = c("theta", "sigma2"),
                      par = c(theta = 0.5, sigma2 = 2))

covSymGauss

## -- simulate a path from the corresponding GP --

nGrid <- 24; n <- nGrid^2; d <- 2
xGrid <- seq(from = -1, to = 1, length.out = nGrid)
Xgrid <- expand.grid(x1 = xGrid, x2 = xGrid)

Kmat <- covMat(object = covSymGauss, X = Xgrid,
               compGrad = FALSE, index = 1L)

library(MASS)
set.seed(1)
ygrid <- mvrnorm(mu = rep(0, n), Sigma = Kmat)

## -- extract a design and the corr. response from the grid --

nDesign <- 25
tab <- subset(cbind(Xgrid, ygrid), x1 > 0 & x2 > 0)
rowIndex <- seq(1, nrow(tab), length = nDesign)
X <- tab[rowIndex, 1:2]
y <- tab[rowIndex, 3]

opar <- par(mfrow = c(1, 3))
contour(x = xGrid, y = xGrid,
        z = matrix(ygrid, nrow = nGrid, ncol = nGrid), 
        nlevels = 15)
abline(h = 0, v = 0, col = "SpringGreen3")
points(x2 ~ x1, data = X, type = "p", pch = 21,
       col = "orangered", bg = "yellow", cex = 0.8)
title("GRF Simulation")


## -- Fit the Gaussian process model (trend + covariance parameters) -- 
covSymGauss
symgp <- gp(formula = y ~ 1, data = data.frame(y, X),
            inputs = names(X),
            cov = covSymGauss,
            parCovIni = c(0.1, 2),
            varNoiseIni = 1.0e-8,
            varNoiseLower = 0.9e-8, varNoiseUpper = 1.1e-8)

# mind that the noise is not a symmetric kernel
# so varNoiseUpper should be chosen as small as possible.

summary(symgp)

## -- predict and compare --

predSymgp <- predict(object = symgp, newdata = Xgrid, type = "UK")

contour(x = xGrid, y = xGrid,
        z = matrix(predSymgp$mean, nrow = nGrid, ncol = nGrid),
        nlevels = 15)
abline(h = 0, v = 0, col = "SpringGreen3")
points(x2 ~ x1, data = X, type = "p", pch = 21,
       col = "orangered", bg = "yellow", cex = 0.8)
title("Kriging mean")

contour(x = xGrid, y = xGrid,
        z = matrix(predSymgp$sd, nrow = nGrid, ncol = nGrid),
        nlevels = 15)
abline(h = 0, v = 0, col = "SpringGreen3")
points(x2 ~ x1, data = X, type = "p", pch = 21,
       col = "orangered", bg = "yellow", cex = 0.8)
title("Kriging s.d.")

par(opar)
}
