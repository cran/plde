\name{basic_values}
\alias{basic_values}
\title{Compute basic values}
\usage{basic_values(sm)}
\arguments{
   \item{sm}{List of plde fit}
}

\description{Compute basic values}

\details{
\code{basic_values} function computes
transformed variable (\code{sm$X_transform}),
rectangular node points (\code{sm$nodes}) and weights (\code{sm$weights})
for numerical integrations,
coefficient vector (\code{sm$coefficients}),
basis matrix at node and data points (\code{sm$B_mat}, \code{sm$X_mat}),
and basis mean (\code{sm$B_mean}).
}

\seealso{\code{\link{legendre_polynomial}}}

\author{
JungJun Lee, Jae-Hwan Jhong, Young-Rae Cho, SungHwan Kim, Ja-yong Koo
}

\references{
JungJun Lee, Jae-Hwan Jhong, Young-Rae Cho, SungHwan Kim and Ja-Yong Koo. "Penalized Log-density Estimation Using Legendre Polynomials."
Submitted to Communications in Statistics - Simulation and Computation (2017), in revision.
}
