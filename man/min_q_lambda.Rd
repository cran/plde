\name{min_q_lambda}
\alias{min_q_lambda}
\title{Minimization of the quadratic approximation to objective function}
\usage{min_q_lambda(sm)}
\arguments{
   \item{sm}{List of plde fit}
}

\description{\code{min_q_lambda} function gives the coefficient vector
(\code{sm$coefficients}) updated by the coordinate descent algorithm iteratively
until the quadratic approximation to the objective function convergences.}

\seealso{\code{\link{q_lambda}}, \code{\link{update}}}

\author{
JungJun Lee, Jae-Hwan Jhong, Young-Rae Cho, SungHwan Kim, Ja-yong Koo
}

\references{
JungJun Lee, Jae-Hwan Jhong, Young-Rae Cho, SungHwan Kim and Ja-Yong Koo. "Penalized Log-density Estimation Using Legendre Polynomials."
Submitted to Communications in Statistics - Simulation and Computation (2017), in revision.
}
