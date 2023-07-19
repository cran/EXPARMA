# Takes two inputs, i.e., the data and order of the model to be fitted and
# returns the initial values over which optimization is to be done for fitting
# the model. The initial values are taken from an ARIMA model with the same
# order and appropriate differencing. These will act as default initial values
# in the EXPARMA function, unless specified by user. Without these values, the
# optimization process may fail.

init_val <- function(ts_data, order) {
  p <- order[1]
  q <- order[2]
  N <- length(ts_data)
  k <- 2*(p+q+1)
  d <- forecast::auto.arima(ts_data)$arma[6]
  arima_param <- unname((forecast::Arima(ts_data, order = c(p,d,q)))$coef)
  arima_init <- c(arima_param[1:p], rep(0.5,p), arima_param[(p+1):(p+q)],
                  rep(0.5,(q+2)))
  return(arima_init)
}

