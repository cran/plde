#===============================================================================
basic_values = function(sm)
{
   # find an affine map for transformation of the given data
   X_min = min(sm$X)
   X_max = max(sm$X)
   sm$a = (sm$U - sm$L) / (X_max - X_min)
   sm$b = sm$U - sm$a * X_max
   # transform the data
   sm$X_transform = sm$a * sm$X + sm$b
   # compute rectangular points and weights
   sm$nodes = seq(-1, 1, length = sm$number_rectangular)
   sm$weights = sm$nodes[2] - sm$nodes[1]
   # coefficient vectors
   sm$coefficients = c(rep(0, sm$dimension))
   # compute B_mat and X_mat
   sm$B_mat = legendre_polynomial(sm$nodes, sm)
   sm$X_mat = legendre_polynomial(sm$X_transform, sm)
   # compute B_mean
   sm$B_mean = colMeans(sm$X_mat)
   return(sm)
}
#===============================================================================
compute_fitted = function(x, sm)
{
   x_transform = sm$a * x + sm$b
   B_mat = legendre_polynomial(x_transform, sm)
   f_X_transform = exp(B_mat %*% sm$coefficients - sm$c_coefficients)
   f_X = f_X_transform * sm$a
   return(f_X)
}
#===============================================================================
compute_lambdas = function(sm)
{
   lambda_max = max(abs(sm$B_mean)) + sm$epsilon
   lambda_min = sm$epsilon * lambda_max
   sm$lambda_all = seq(log(lambda_max), log(lambda_min), length = sm$number_lambdas)
   sm$lambda_all = exp(sm$lambda_all)
   return(sm)
}
#===============================================================================
fit_plde = function(sm)
{
   sm = fit_plde_sub(sm)
   for(iter in 1 : sm$max_iterations)
   {
      old_coefficients = sm$coefficients
      old_pen_loglik = sm$pen_loglik
      # compute score and information
      sm$B_int = crossprod(sm$B_mat, sm$f)
      sm$score = sm$B_mean - sm$B_int
      sm$information = crossprod(sm$B_mat, sm$f * sm$B_mat) -
         sm$B_int %*% t(sm$B_int)
      # compute base model
      sm$base_coefficients = sm$coefficients
      # minimize Quadratic approximation
      sm = min_q_lambda(sm)
      sm = fit_plde_sub(sm)
      # check convergence
      if (abs(sm$pen_loglik - old_pen_loglik) < sm$epsilon)
         break
   }
   return(sm)
}
#===============================================================================
fit_plde_sub = function(sm)
{
   sm$f = exp(sm$B_mat %*% sm$coefficients)
   # calculate c_coefficients
   sm$c_coefficients = log(sum(sm$f * sm$weights))
   # calculate f
   sm$f = sm$f / exp(sm$c_coefficients)
   # calculate fw = f * weight
   sm$f = as.vector(sm$f * sm$weights)
   # calculate pen_loglik
   sm$loglik = sm$B_mean %*% sm$coefficients - sm$c_coefficients
   sm$pen_loglik = - sm$loglik + sm$lambda * sum(abs(sm$coefficients))
   return(sm)
}
#===============================================================================
legendre_polynomial = function(x, sm)
{
   x = as.vector(x)
   dimension_p1 = sm$dimension + 1
   P = matrix(0, length(x), dimension_p1)
   P[, 1 : 2] = cbind(1, x)
   for (n in 2 : sm$dimension)
      P[, n + 1] = ((2 * n - 1) * x * P[, n] - (n - 1) * P[, n - 1]) / n
   P = P[, 2 : dimension_p1]
   return(P)
}
#===============================================================================
min_q_lambda = function(sm)
{
   # fit the initial model
   sm = q_lambda(sm)
   # coordinate descent algorithm begins
   for (iter in 1 : sm$max_iterations)
   {
      if (sm$verbose)
         cat("iter = ", iter, "Q_lambda = ", sm$Q_lambda, "\n")
      # we update sequentially
      # an updated coefficient will be current coefficient for next update
      old_Q_lambda = sm$Q_lambda
      sm = update(sm)
      sm = q_lambda(sm)
      # check convergence
      if (abs(sm$Q_lambda - old_Q_lambda) < sm$epsilon) break
   }
   return(sm)
}
#===============================================================================
model_selection = function(fit, method = 'AIC')
{
   ic = rep(0, fit$number_lambdas)
   penalty = log(fit$sample_size)
   if (method == 'AIC')
      penalty = 2
   for (k in 1 : fit$number_lambdas)
   {
      p = sum(fit$sm[[k]]$coefficients != 0)
      ic[k] = -2 * fit$sample_size * fit$sm[[k]]$loglik + penalty * p
   }
   fit$optimal_idx = which.min(ic)
   fit$optimal = fit$sm[[fit$optimal_idx]]
   return(fit)
}
#===============================================================================
plde = function(X, initial_dimension = 100, number_lambdas = 200,
                L = -0.9, U = 0.9, ic = 'AIC',
                epsilon = 1e-5,
                max_iterations = 1000,
                number_rectangular = 1000,
                verbose = FALSE)
{
   cat("=================================================\n")
   cat("Penalized Log-density Estimation\n")
   cat("Using Legendre Polynomials\n")
   cat("Version 0.1 by SDMLAB (May 28, 2018)\n")
   cat("Department of Statistics, Korea University, Korea\n")
   cat("=================================================\n")
   cat("initial dimension = ", initial_dimension,
       "number of lambdas = ", number_lambdas, "\n")
   # setup the process
   sm = list(X = X,
             sample_size = length(X),
             dimension = initial_dimension,
             number_lambdas = number_lambdas,
             L = L,
             U = U,
             epsilon = epsilon,
             max_iterations = max_iterations,
             number_rectangular = number_rectangular,
             ic = ic,
             verbose = verbose)
   sm = basic_values(sm)
   # calculate lambda
   sm = compute_lambdas(sm)
   fit = list(number_lambdas = sm$number_lambdas, sample_size = sm$sample_size)
   for (k in 1 : fit$number_lambdas)
   {
      # save tuning parameter
      sm$lambda = sm$lambda_all[k]
      # fit plde
      sm = fit_plde(sm)
      fit$sm[[k]] = sm
   }
   fit = model_selection(fit, ic)
   return(fit)
}
#===============================================================================
q_lambda = function(sm)
{
   delta = sm$coefficients - sm$base_coefficients
   Q = 0.5 * sum(delta %*% sm$information %*% delta) - sum(sm$score * delta)
   sm$Q_lambda = Q + sm$lambda * sum(abs(sm$coefficients))
   return(sm)
}
#===============================================================================
soft_thresholding = function(y, threshold)
{
   if (abs(y) <= threshold)
      return(0)
   if (y > threshold)
      return(y - threshold)
   if (y < -threshold)
      return(y + threshold)
}
#===============================================================================
update = function(sm)
{
   for (j in 1 : sm$dimension)
   {
      delta = sm$coefficients - sm$base_coefficients
      Ijj = sm$information[j, j]
      dQ_j = sum(sm$information[j, ] * delta) - sm$score[j]
      y = sm$coefficients[j] - dQ_j / Ijj
      sm$coefficients[j] = soft_thresholding(y, sm$lambda / Ijj)
   }
   return(sm)
}
#===============================================================================
