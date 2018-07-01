\name{legendre_polynomial}
\alias{legendre_polynomial}
\title{legendre_polynomial}
\usage{legendre_polynomial(x, sm)}
\arguments{
   \item{x}{input node points}
   \item{sm}{List of plde fit}
}

\description{\code{legendre_polynomial} gives the Legendre polynomial design matrix
over the input node points.}

\examples{
# clean up
rm(list = ls())
library(plde)
x = seq(-1, 1, length = 200)
L = legendre_polynomial(x, list(dimension = 10))
# Legendre polynomial basis for dimension 1 to 10
matplot(x, L, type = "l")
}

\author{
JungJun Lee, Jae-Hwan Jhong, Young-Rae Cho, SungHwan Kim, Ja-yong Koo
}

\references{
JungJun Lee, Jae-Hwan Jhong, Young-Rae Cho, SungHwan Kim and Ja-Yong Koo. "Penalized Log-density Estimation Using Legendre Polynomials."
Submitted to Communications in Statistics - Simulation and Computation (2017), in revision.
}
