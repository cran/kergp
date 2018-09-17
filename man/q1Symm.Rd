\name{q1Symm}
\alias{q1Symm}
\title{
  Qualitative Correlation or Covariance Kernel with one Input
  and General Symmetric Correlation
}
\usage{

  q1Symm(factor, input = "x", cov = c("corr", "homo", "hete"))

}
\arguments{
  \item{factor}{

    A factor with the wanted levels for the covariance kernel object.

  }
  \item{input}{

    Name of (qualitative) input for the kernel.

  }
  \item{cov}{
    
    Character telling if the result is a correlation kernel, an
    homoscedastic covariance kernel or an heteroscedastic covariance
    kernel with an arbitrary variance vector.
    
  }

}
\value{

  An object with class \code{"covQual"} with \code{d = 1} qualitative
  input.

}
\description{

  Qualitative correlation or covariance kernel with one input and
  general symmetric correlation.

}
\note{

  Correlation kernels are needed in tensor products because the tensor
  product of two covariance kernels each with unknown variance would not
  be identifiable.

}
\examples{
School <- factor(1L:3L, labels = c("Bad", "Mean" , "Good"))
myCor <- q1Symm(School, input = "School")
coef(myCor) <- c(theta_2_1 = pi / 3, theta_3_1 = pi / 4, theta_3_2 = pi / 8)

## Use a data.frame with a factor
set.seed(246)
newSchool <- factor(sample(1L:3L, size = 20, replace = TRUE),
                    labels = c("Bad", "Mean" , "Good"))
C1 <- covMat(myCor, X = data.frame(School = newSchool),
             compGrad = FALSE, lowerSQRT = FALSE)

## or use a matrix with integer content
C2 <- covMat(myCor, X = cbind(School = as.integer(newSchool)),
             compGrad = FALSE, lowerSQRT = FALSE)
max(abs(C1 - C2))
}
\seealso{
The \code{\link{corLevSymm}} function used to compute the
correlation matrix and its gradients w.r.t. the correlation
parameters.
}
