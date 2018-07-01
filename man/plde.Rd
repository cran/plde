\name{plde}
\alias{plde}
\title{Penalized Log-density Estimation Using Legendre Polynomials}
\usage{plde(X, initial_dimension = 100, number_lambdas = 200,
     L = -0.9, U = 0.9, ic = 'AIC', epsilon = 1e-5, max_iterations = 1000,
     number_rectangular = 1000, verbose = FALSE)}
\arguments{
   \item{X}{Input vector, of dimension \eqn{n}.}
   \item{initial_dimension}{Positive interger that decides
                            initial dimension of Legendre polynomials.
                            Default is 100.}
   \item{number_lambdas}{The number of tuning parameter \eqn{\lambda} values.
                         Default is 200.}
   \item{L}{Lower bound of transformed data. Default is -0.9.}
   \item{U}{Upper bound of transformed data. Default is +0.9.}
   \item{ic}{Model selection criteria. 'AIC' or 'BIC' is used. Default is 'AIC'.}
   \item{epsilon}{Positive real value that controls the
                  iteration stopping criteria. In general, the smaller the value,
                  convergence needs more iterations.
                  Default is 1e-5.}
   \item{max_iterations}{Positive integer value that decides the maximum number of
                         iterations. Default is 1000.}
   \item{number_rectangular}{Number of node points for numerical integration}
   \item{verbose}{verbose}
}
\value{A list contains the whole fits of all tuning parameter \eqn{\lambda} sequence.
       For example, {\code{fit$sm[[k]]}} indicates the fit of
       \eqn{k} th lambda.}

\description{
This function gives the penalized log-density estimation
using Legendre polynomials.
}

\details{
The basic idea of implementation is to approximate the negative
log-likelihood function by a quadratic function and then
to solve penalized quadratic optimization problem
using a coordinate descent algorithm.
For a clear exposition of coordinate-wise updating scheme,
we briefly explain a penalized univariate quadratic problem
and its solution expressed as soft-thresholding operator
\code{soft_thresholding}.
We use this univariate case algorithm to update parameter vector
coordinate-wisely to find a minimizer.
}

\seealso{\code{\link{basic_values}}, \code{\link{compute_lambdas}}, \code{\link{fit_plde}},
         \code{\link{model_selection}}}

\author{
JungJun Lee, Jae-Hwan Jhong, Young-Rae Cho, SungHwan Kim, Ja-yong Koo
}

\examples{
# clean up
rm(list = ls())
library(plde)
Eruption = faithful$eruptions
Waiting = faithful$waiting
n = length(Eruption)
# fit PLDE
fit_Eruption = plde(Eruption, initial_dimension = 30, number_lambdas = 50)
fit_Waiting = plde(Waiting, initial_dimension = 30, number_lambdas = 50)
x_Eruption = seq(min(Eruption), max(Eruption), length = 100)
x_Waiting = seq(min(Waiting), max(Waiting), length = 100)
fhat_Eruption = compute_fitted(x_Eruption, fit_Eruption$sm[[fit_Eruption$number_lambdas]])
fhat_Waiting = compute_fitted(x_Waiting, fit_Waiting$sm[[fit_Waiting$number_lambdas]])
# display layout
par(mfrow = c(2, 2), oma=c(0,0,2,0), mar = c(4.5, 2.5, 2, 2))
#=======================================
# Eruption
#=======================================
col_index = rainbow(fit_Eruption$number_lambdas)
plot(x_Eruption, fhat_Eruption, type = "n", xlab = "Eruption", ylab = "", main = "")
# all fit plot
for(i in 1 : fit_Eruption$number_lambdas)
{
   fhat = compute_fitted(x_Eruption, fit_Eruption$sm[[i]])
   lines(x_Eruption, fhat, lwd = 0.5, col = col_index[i])
}
k_Eruption = density(Eruption, bw = 0.03)
lines(k_Eruption$x, k_Eruption$y / 2, lty = 2)

# optimal model
hist_col = rgb(0.8,0.8,0.8, alpha = 0.6)
hist(Eruption, nclass = 20, freq = FALSE, xlim = c(1.1, 5.9),
     col = hist_col, ylab = "", main = "", ylim = c(0, 1.2))
fhat_optimal_Eruption = compute_fitted(x_Eruption, fit_Eruption$optimal)
lines(x_Eruption, fhat_optimal_Eruption, col = "black", lwd = 2)
#========================================
# Waiting
#========================================
col_index = rainbow(fit_Waiting$number_lambdas)
plot(x_Waiting, fhat_Waiting, type = "n", xlab = "Waiting", ylab = "", main = "")
# all fit plot
for(i in 1 : fit_Waiting$number_lambdas)
{
   fhat = compute_fitted(x_Waiting, fit_Waiting$sm[[i]])
   lines(x_Waiting, fhat, lwd = 0.5, col = col_index[i])
}
k_Waiting = density(Waiting, bw = 1)
lines(k_Waiting$x, k_Waiting$y / 2, lty = 2)

# optimal model
hist_col = rgb(0.8,0.8,0.8, alpha = 0.6)
hist(Waiting, nclass = 20, freq = FALSE, xlim = c(40, 100),
     col = hist_col, ylab = "", main = "", ylim = c(0, 0.055))
fhat_optimal_Waiting = compute_fitted(x_Waiting, fit_Waiting$optimal)
lines(x_Waiting, fhat_optimal_Waiting, col = "black", lwd = 2)
}

\source{
This package is built on R version 3.4.2.
}

\references{
JungJun Lee, Jae-Hwan Jhong, Young-Rae Cho, SungHwan Kim and Ja-Yong Koo. "Penalized Log-density Estimation Using Legendre Polynomials."
Submitted to Communications in Statistics - Simulation and Computation (2017), in revision.

Friedman, Jerome, Trevor Hastie, and Rob Tibshirani. "Regularization paths for generalized linear models via coordinate descent." Journal of statistical software 33.1 (2010): 1.
}

