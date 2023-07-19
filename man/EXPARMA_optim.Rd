
\name{EXPARMA_optim}
\alias{EXPARMA_optim}

\title{Fitting of Exponential Autoregressive Moving Average (EXPARMA) Model}
\description{This function will fit EXPARMA model for given parameters.}
\usage{EXPARMA_optim(ts_data, order, init, opt_method)}
%-
\arguments{
  \item{ts_data}{
a univarite time series data

}
\item{order}{
order represents the values of autoregresive (p) and moving average (q) parameters of EXPARMA model. p and q will take integer values of 1 or more. By default these values are set as 1.

}
\item{init}{
init represents the initial values of the model EXPARMA. These values can be user specified or the output from init_val() can be used.
}

\item{opt_method}{
opt_method represents the optimiation method used for estimating the model parameters.All methods available under optim() can be used here. By default we have made use of "BFGS".
}
}
\details{Takes 2 compulsory inputs, i.e., the data and order of model to be fitted.
Optimizes the parameters of the model of the given order (p,q) by minimizing. RSS and returns a fit with the optimized parameters. The fit is returned using the previous function EXPARMAfit(), with the parameters inputted being the optimized ones. All other inputs to this function are arguments for the optim() part.

}
\value{It returns the model parameters, RSS, AIC, fitted values and residuals.}

\examples{
datats=c(17597,14074,11425,11691,11298,12351,14311,
12349,10537,11755,13154,11989,13022,12107,11172,10667,
10091,12204,12274,22343)
EXPARMA_optim(datats,order=c(1,1))
}

