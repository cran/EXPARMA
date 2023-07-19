\name{init_val}
\alias{init_val}

\title{Fitting of Exponential Autoregressive Moving Average (EXPARMA) Model}
\description{This function will fit EXPARMA model for given parameters.}

\usage{init_val(ts_data, order)}
%-
\arguments{
  \item{ts_data}{
a univarite time series data

}
\item{order}{
order represents the values of autoregresive (p) and moving average (q) parameters of EXPARMA model. p and q will take integer values of 1 or more. By default these values are set as 1.

}
}
\details{This function takes two inputs, i.e., the data and order of the model to be fitted and returns the initial values over which optimization is to be done for fitting the model. The initial values are taken from an ARIMA model with the same order and appropriate differencing. These will act as default initial values in the EXPARMA function, unless specified by user. Without these values, the optimization process may fail.
}
\value{It returns the initial values of the EXAPRMA modelver which optimization is to be done for fitting the model.}
\examples{
datats=c(17597,14074,11425,11691,11298,12351,14311,
12349,10537,11755,13154,11989,13022,12107,11172,10667,
10091,12204,12274,22343)
init_val(datats,order=c(1,1))
}
