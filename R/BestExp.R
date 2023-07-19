# Takes only one input, the data. Checks the AIC values of all EXPARMA models
# from order (1,1) up to (max.p, max.q) by finding the optimal parameters of
# those models and returns a fit with the model having the least AIC. The
# maximum order up to which checks are set to 5 by default.

BestExp <- function(ts_data, opt_method = "BFGS", max.p = 5, max.q = 5) {
  AIC <- NULL
  P <- NULL
  Q <- NULL
  for (i in 1:max.p) {
    for (j in 1:max.q) {
      AIC <- append(AIC, (EXPARMA_optim(ts_data, order = c(i,j)))$AIC)
      P <- append(P, i)
      Q <- append(Q, j)
    }
  }
  df <- data.frame("p" = P, "q" = Q, "AIC" = AIC)
  order_best <- unname(unlist(df[which(df$AIC == min(df$AIC)), 1:2]))
  best_model <- EXPARMA_optim(ts_data, order = order_best,
                        opt_method = opt_method)
  message(paste0("Best model is of order p=", order_best[1],
                 " and q=", order_best[2], " on the basis of minimum AIC"))
  return(best_model)
}

