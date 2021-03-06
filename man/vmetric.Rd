% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/padicmetric.r
\name{vmetric}
\alias{vmetric}
\title{Computes v metric to create a matrix of p-adic transformed data.}
\usage{
vmetric(data, reg_eval = TRUE, ncores = 1)
}
\arguments{
\item{data}{(numeric) data matrix with features as rows and items as columns.}

\item{ncores}{(Integer) Number of cores. Default is 1.}

\item{reg.eval}{(logical) Compute R^2 to evaluate fitting of logarithmic approximation. Default is TRUE.}
}
\value{
matrix of v values (p-adic metric):
}
\description{
Computes v metric to create a matrix of p-adic transformed data.
}
\examples{

test ...

}
\author{
Ivan Imaz \email{ii236@cam.ac.uk}
}
