\name{model_selection}
\alias{model_selection}
\title{Optimal model selection}
\usage{model_selection(fit, method = "AIC")}
\arguments{
   \item{fit}{Entire list of plde fit by all tuning parameters}
   \item{method}{model selection criteria. 'AIC' or 'BIC' is used. Default is 'AIC'.}
}

\description{\code{model_selection} function gives the optimal model over the whole
plde fits based on information criterian (AIC, BIC).
The optimal model is saved at \code{fit$optimal}.}

\author{
JungJun Lee, Jae-Hwan Jhong, Young-Rae Cho, SungHwan Kim, Ja-yong Koo
}

\references{
JungJun Lee, Jae-Hwan Jhong, Young-Rae Cho, SungHwan Kim and Ja-Yong Koo. "Penalized Log-density Estimation Using Legendre Polynomials."
Submitted to Communications in Statistics - Simulation and Computation (2017), in revision.
}
