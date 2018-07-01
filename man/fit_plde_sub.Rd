\name{fit_plde_sub}
\alias{fit_plde_sub}
\title{Fit plde for a fixed tuning parameter}
\usage{fit_plde_sub(sm)}
\arguments{
   \item{sm}{List of plde fit}
}

\description{\code{fit_plde_sub} function computes
the updated normalizng constant (\code{sm$c_coefficients}),
Legendre density function estimator (\code{sm$f})
and the negative of penalized log-likelihood function (\code{sm$pen_loglik})
for each iteration.}

\author{
JungJun Lee, Jae-Hwan Jhong, Young-Rae Cho, SungHwan Kim, Ja-yong Koo
}

\references{
JungJun Lee, Jae-Hwan Jhong, Young-Rae Cho, SungHwan Kim and Ja-Yong Koo. "Penalized Log-density Estimation Using Legendre Polynomials."
Submitted to Communications in Statistics - Simulation and Computation (2017), in revision.
}
