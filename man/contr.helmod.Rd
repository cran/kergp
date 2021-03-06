% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/miscUtils.R
\name{contr.helmod}
\alias{contr.helmod}
\title{Modified Helmert Contrast Matrix}
\usage{
contr.helmod(n)
}
\arguments{
\item{n}{Integer.}
}
\value{
An orthogonal matrix with \code{n} rows and \code{n - 1}
columns. The columns form a basis of the subspace orthogonal to
a vector of \code{n} ones.
}
\description{
Modified Helmert contrast (or coding) matrix.
}
\details{
The returned matrix is a scaled version of
\code{contr.helemert(A)}.
}
\examples{
A <- contr.helmod(6)
crossprod(A)

}
