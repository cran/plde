\name{compute_fitted}
\alias{compute_fitted}
\title{compute_fitted}
\usage{compute_fitted(x, sm)}
\arguments{
   \item{x}{grid points}
   \item{sm}{List of plde fit}
}

\description{\code{compute_fitted} function gives the fitted values
over the input grid points for the fixed tuning parameter \eqn{\lambda}.}

\details{
\code{compute_fitted} function computes fitted values of estimates
having support for the given data by scaling back
and change of variable technique. For more details, see Section 3.2
of the reference.
}

\seealso{\code{\link{legendre_polynomial}}}

\author{
JungJun Lee, Jae-Hwan Jhong, Young-Rae Cho, SungHwan Kim, Ja-yong Koo
}

\references{
JungJun Lee, Jae-Hwan Jhong, Young-Rae Cho, SungHwan Kim and Ja-Yong Koo. "Penalized Log-density Estimation Using Legendre Polynomials."
Submitted to Communications in Statistics - Simulation and Computation (2017), in revision.
}
