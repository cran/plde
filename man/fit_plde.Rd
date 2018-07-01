\name{fit_plde}
\alias{fit_plde}
\title{Fit plde for a fixed tuning parameter}
\usage{fit_plde(sm)}
\arguments{
   \item{sm}{List of plde fit}
}

\description{\code{fit_plde} gives the plde fit for a fixed tuning parameter}

\details{
This is the coordinate descent algorithm for computing \eqn{\hat\theta^\lambda}
when the penalty parameter \eqn{\lambda} is fixed.
See Algorithm 1 in the reference for more details.
}

\seealso{
\code{\link{fit_plde_sub}}, \code{\link{min_q_lambda}}
}

\author{
JungJun Lee, Jae-Hwan Jhong, Young-Rae Cho, SungHwan Kim, Ja-yong Koo
}

\references{
JungJun Lee, Jae-Hwan Jhong, Young-Rae Cho, SungHwan Kim and Ja-Yong Koo. "Penalized Log-density Estimation Using Legendre Polynomials."
Submitted to Communications in Statistics - Simulation and Computation (2017), in revision.
}
