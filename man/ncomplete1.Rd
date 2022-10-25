\name{ncomplete1.for}
\alias{ncomplete1.for}
\alias{ncomplete1}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Function called by ncomplete.for}
\description{
Function called by ncomplete.for}
\usage{
ncomplete1.for(Z,NDIR=10000)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{Z}{The data set Z has to be a matrix with ncol(Z)-1 <= nrow(Z) <= 10000.
           The first ncol(Z)-1 columns of Z are the design matrix X.
           The last column of Z is the binary response vector y (0/1).}
  \item{NDIR}{Maximal number of directions (integer)}
}
\value{
A list with components
  \item{Z}{data matrix}
  \item{NDIR}{10000 directions}
}
\author{Andreas Christmann, Peter J. Rousseeuw
\email{Christmann@statistik.uni-dortmund.de}}
\references{}
\seealso{\code{\link{ncomplete.for}}}
\keyword{regression}
\keyword{robust}
\examples{}
%
%
%
