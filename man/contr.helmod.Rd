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
  
  An orthogonal matrix with \code{n} rows and \code{n - 1} columns. The
  columns form a basis of the subspace orthogonal to a vector of
  \code{n} ones.
  
}
\description{
  
  Modified Helmert contrast (or coding) matrix.

}
\details{
  
The returned matrix is a scaled version of \code{contr.helmert(A)}.

}
\examples{
A <- contr.helmod(6)
crossprod(A)

}
