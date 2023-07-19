

# Takes 3 inputs, all compulsory, i.e., data, order and parameters of the model
# to be fitted and returns the fitted values, residuals, RSS and AIC for the
# fitted model. No optimisation is done. Useful for simulation of data with
# given order and parameters.

EXPARMAfit <- function(ts_data, order, par) {

  p <- order[1]
  q <- order[2]
  N <- length(ts_data)
  k <- 2*(p+q+1)

  ar_par_phi <- par[1:p]
  ar_par_pi <- par[(p+1):(2*p)]
  ma_par_theta <- par[((2*p)+1):((2*p)+q)]
  ma_par_delta <- par[(((2*p)+q)+1):(2*(p+q))]
  ar_scale_par <- par[(2*(p+q))+1]
  ma_scale_par <- par[(2*(p+q+1))]
  Predicted <- rep(NA, times = max(p,q))
  Residuals <- rep(0, times = max(p,q))
  ar_part <- ma_part <- NULL
  for (t in ((max(p,q)):(N-1))) {
    for (i in 1:p){
      ar_part[i] <- (ar_par_phi[i] + (ar_par_pi[i]*exp(-1*(ar_scale_par)*((ts_data[t])^2))))*ts_data[t-i+1]
    }
    for (j in 1:q) {
      ma_part[j] <- (ma_par_theta[j] + (ma_par_delta[j]*exp(-1*(ma_scale_par)*((Residuals[t])^2))))*Residuals[t-j+1]
    }
    Predicted[t+1] <- sum(c(ar_part, ma_part))
    Residuals[t+1] <- (ts_data[t+1] - Predicted[t+1])
  }
  RSS <- sum(Residuals^2)
  AIC <- (N*log(RSS/(N-k))) + (2*k)
  model.fit <- list(Fitted = Predicted, Residuals = Residuals,
                    RSS = RSS, AIC = AIC)
  return(model.fit)
}





