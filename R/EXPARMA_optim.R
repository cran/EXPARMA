
# Takes 2 compulsory inputs, i.e., the data and order of model to be fitted.
# Optimizes the parameters of the model of the given order (p,q) by minimizing
# RSS and returns a fit with the optimized parameters. The fit is returned
# using the previous function EXPARMAfit(), with the parameters inputted being
# the optimized ones. All other inputs to this function are arguments for the
# optim() part.

EXPARMA_optim <- function(ts_data, order, init, opt_method = "BFGS") {

  p <- order[1]
  q <- order[2]
  N <- length(ts_data)
  k <- 2*(p+q+1)
  if(missing(init)) init <- init_val(ts_data = ts_data, order = order)

  # Prameter estimation using optim()
  RSS <- function(par) {
    return(EXPARMAfit(ts_data = ts_data, order = order, par)$RSS)
  }
  opt_result <- stats::optim(par = init, fn = RSS, method = opt_method)

  # Save parameter estimates
  ar_par_phi <- c(opt_result$par[1:p],
                  rep(NA, max((max(p,q)-p),0)))
  ar_par_pi <- c(opt_result$par[(p+1):(2*p)],
                 rep(NA, max((max(p,q)-p),0)))
  ma_par_theta <- c(opt_result$par[((2*p)+1):((2*p)+q)],
                    rep(NA, max((max(p,q)-q),0)))
  ma_par_delta <- c(opt_result$par[(((2*p)+q)+1):(2*(p+q))],
                    rep(NA, max((max(p,q)-q),0)))
  ar_scale_par <- c(opt_result$par[(2*(p+q))+1],
                    rep(NA, max(p,q)-1))
  ma_scale_par <- c(opt_result$par[(2*(p+q+1))],
                    rep(NA, max(p,q)-1))
  df_pars <- data.frame("Phi" = ar_par_phi, "Pi" = ar_par_pi,
                        "Theta" = ma_par_theta, "Delta" = ma_par_delta,
                        "Gamma" = ar_scale_par, "Omega" = ma_scale_par)

  # Re-fit the model with optimised parameters
  model.fit <- EXPARMAfit(ts_data = ts_data, order = order,
                          par = opt_result$par)

  # Print parameters, RSS, AIC, fitted values and residuals
  # Also save other outputs of optim
  out <- list(series = deparse(substitute(ts_data)), order = c(p,q), n=N, k=k,
              par = df_pars, RSS = opt_result$value, AIC = model.fit$AIC,
              Fitted = model.fit$Fitted, Residuals = model.fit$Residuals,
              counts = opt_result$counts, convergence = opt_result$convergence,
              message = opt_result$message)
  return(out)
}
