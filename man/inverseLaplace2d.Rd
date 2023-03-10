% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inverseLaplace2d.R
\name{inverseLaplace2d}
\alias{inverseLaplace2d}
\title{Inverse Laplace transform in 2D}
\usage{
inverseLaplace2d(t1, t2, h.star, ..., n = 15, n.cores = 1)
}
\arguments{
\item{t1}{First argument to Laplace transform}

\item{t2}{Second argument to Laplace transform}

\item{h.star}{FUN Laplace transform to be inverted}

\item{...}{Additional arguments to h.star}

\item{n}{Accuracy parameter}

\item{n.cores}{Number of cores to use defaults to 1}
}
\value{
Inverse Laplace transform
}
\description{
Inverse Laplace transform in 2D
}
\examples{
inverseLaplace2d(c(1,2),c(1,2),function(s1,s2){1/((s1+1)*(s2+1))})
}
