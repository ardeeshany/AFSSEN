% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/norm_matrix_K.R
\name{norm_matrix_K}
\alias{norm_matrix_K}
\title{Computation of the K-norm of a set of functions}
\usage{
norm_matrix_K(M, eigenval)
}
\arguments{
\item{M}{matrix. \code{J} \eqn{\times} \code{N} matrix containing in the column
\eqn{n} the coefficients of the projection of a function \eqn{y_n} on the
\code{J} eigenfunctions of the kernel.}

\item{eigenval}{vector. \code{J}-length vector of the eigenvalues of the
kernel.}
}
\value{
vector of length \code{N} containing the K-norm of the functions
}
\description{
It computes the norm associated to the kernel K of a set of functions.
The functions are represented as their projection on a
kernel basis.
}
\examples{
data(SobolevKernel)
data(simulation)
norm_K = sum(norm_matrix_K(Y_matrix, eigenval))

}
