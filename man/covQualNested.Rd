\name{covQualNested}
\alias{covQualNested}
\title{Nested Qualitative Covariance}
\usage{

covQualNested(input = "x",
              group = NULL,
              nestedLevels = NULL,
              between = "Symm",
              within = "Diag",
              covBet = c("corr", "homo", "hete"), 
              covWith = c("corr", "homo", "hete"),
              compGrad = TRUE,
              contrasts = contr.helmod)

}
\arguments{

  \item{input}{
    
    Name of the input, i.e. name of the column in the
    data frame when the covariance kernel is evaluated with the
    \code{\link{covMat,covQual-method}} method.
    
  }
  \item{group}{

    A factor of vector giving the groups, see \bold{Examples}.

  }
  \item{nestedLevels}{
    
    A factor or a vector giving the (nested) levels within the group for
    each level of \code{group}. If this is missing, each element of
    \code{group} is assumed to correspond to one nested level within the
    group and the levels within the group are taken as integers in the order
    of \code{group} elements.
    
  }
  \item{between}{
    
    Character giving the type of structure to use for the \emph{between}
    part. For now this can be one of the three choices \code{"Diag"}, 
    the diagonal structure of \code{\link{q1Diag}}, \code{"Symm"} for
    the general covariance of \code{\link{q1Symm}}, or \code{"CompSymm"}
    for the Compound Symmetry covariance of \code{\link{q1CompSymm}}.
    Default is \code{Symm}, corresponding to a specific correlation value for each pair of        groups. On the other hand, \code{Diag} corresponds to a common correlation value for all pairs of groups.

  }
  \item{within}{
    
    Character vector giving the type of structure to use for the
    \emph{within} part. The choices are the same as for
    \code{between}. The character vector is recycled to have length
    \eqn{G} so the \emph{within} covariances can differ across groups.
    Default is \code{"Diag"}, corresponding to a compound symmetry matrix.
    
  }
  \item{covBet}{
  }
  \item{covWith}{
    
    Character vector indicating the type of covariance matrix to be used
    for the generator between- or within- matrices, 
    as in \code{\link{q1Diag}}, \code{\link{q1Symm}} or \code{\link{q1CompSymm}}:
    correlation ("corr"), homoscedastic ("homo") or heteroscedastic ("hete").
    Partial matching is allowed. 
    This is different from the form of the resulting covariance matrix, see section
    \bold{Caution}.
    
  }
  \item{compGrad}{
    
    Logical.
    
  }
  \item{contrasts}{
    
    Object of class \code{"function"}. This function is similar to the
    \code{\link[stats]{contr.treatment}} or
    \code{\link[stats]{contr.treatment}} functions, but it must return
    an \emph{orthogonal} matrix. For a given integer \code{n}, it
    returns a matrix with \code{n} rows and \code{n - 1} columns forming
    a basis for the supplementary of a vector of ones in the
    \eqn{n}-dimensional Euclidean space. The \code{\link{contr.helmod}}
    can be used to obtain an orthogonal matrix hence defining an
    orthonormal basis.

  }

}
\value{
  
  An object with class \code{"covQualNested"}.

}
\description{

  Nested Qualitative Covariance

}
\note{

  For now the replacement method such as \code{'coef<-'} are inherited
  from the class \code{covQuall}. Consequently when these methods are
  used they do not update the covariance structure in the \code{between}
  slot nor those in the \code{within} (list) slot.
  
  This covariance kernel involves \code{two} categorical (i.e. factor)
  inputs, but these are nested. It could be aliased in the future as
  \code{q1Nested} or \code{q2Nested}.
  
}
\section{Caution}{

  When \code{covBet} and \code{covWith} are zero, the resulting matrix
  \emph{is not a correlation matrix}, due to the mode of
  construction. The "between" covariance matrix is a correlation but
  diagonal blocks are added to the extended matrix obtained by re-sizing
  the "between" covariance into a \eqn{n \times n}{ n * n} matrix.
  
}

\examples{

### Ex 1. See the vignette "groupKernel" for an example 
###       inspired from computer experiments.

### Ex 2. Below an example in data analysis.

country <- c("B", "B", "B", "F", "F" ,"F", "D", "D", "D")
cities <- c("AntWerp", "Ghent" , "Charleroi", "Paris", "Marseille", "Lyon",
            "Berlin", "Hambourg", "Munchen")

## create a nested covariance. Note that it will be used with ONE factor
## input obtained by pasing the group and subgroup.

nest <- covQualNested(input = "ccities",
                      group = country, nestedLevels = cities,
                      between = "Symm", within = "Diag",
                      compGrad = TRUE,
                      covBet = "corr", covWith = "corr")

## 'show' method as automatically invocated
nest

## 'covMat' method: if provided, 'X' must be  a data frame with a
## suitable factor or integer column. The input name was given
## at the creation of the covariance structure

Z <- sample(levels(nest)[[1]], size = 30, replace = TRUE)
df <- data.frame(ccities = factor(Z))
covMat(nest, X = df)
covMat(nest)
\dontrun{
    if (require(corrplot)) corrplot(cov2cor(covMat(nest)))
}

## 'simulate' method. One can give a 'X' formal, but with no
## duplicated values!

Sim <- simulate(nest, nsim = 30)
levs <- levels(nest)[[1]]
matplot(Sim, type = "l", main = "Simulated paths", xaxt = "n")
axis(side = 1, at = seq_along(levs), labels = levs)

## another 'covMat'

cc <- as.factor(sample(paste(country, cities, sep = "/"),
                size = 100, replace = TRUE))
covMat(nest, X = data.frame(ccities = cc))

}
