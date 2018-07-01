\name{soft_thresholding}
\alias{soft_thresholding}
\title{Soft thresholding operator}
\usage{soft_thresholding(y, threshold)}
\arguments{
   \item{y}{input real value}
   \item{threshold}{threshold value}
}

\description{\code{soft_thresholding} gives the soft threshold value of \eqn{y} given
the \code{threshold}. When \code{threshold} increasing, \eqn{y} shrinks to zero.}

\examples{
# clean up
rm(list = ls())
library(plde)
# soft thresholding operater
soft_thresholding(3, 1)
soft_thresholding(-3, 1)
# if the threshold value is large enough, it shrinks to zero
soft_thresholding(-3, 4)
soft_thresholding(3, 4)
# Plot of the soft thresholding operater
y = seq(-3, 3, length = 100)
st = NULL
for (i in 1 : length(y))
   st[i] = soft_thresholding(y[i], 1)
plot(y, y, col = "gray", type = "l", ylab = "ST")
lines(y, st, col = "blue")
}

\author{
JungJun Lee, Jae-Hwan Jhong, Young-Rae Cho, SungHwan Kim, Ja-yong Koo
}

\references{
JungJun Lee, Jae-Hwan Jhong, Young-Rae Cho, SungHwan Kim and Ja-Yong Koo. "Penalized Log-density Estimation Using Legendre Polynomials."
Submitted to Communications in Statistics - Simulation and Computation (2017), in revision.
}
