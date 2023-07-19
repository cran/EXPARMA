\name{EXPARMAfit}
\alias{EXPARMAfit}

\title{Fitting of Exponential Autoregressive Moving Average (EXPARMA) Model}
\description{This function will fit EXPARMA model for given parameters.}

\usage{EXPARMAfit(ts_data, order, par)}
%-
\arguments{
  \item{ts_data}{
a univarite time series data

}
\item{order}{
order represents the values of autoregresive (p) and moving average (q) parameters of EXPARMA model. p and q will take integer values of 1 or more. By default these values are set as 1.

}
\item{par}{
par is the parametric space of the EXPARMA model. The total number of parametrs are 2*(p+q+1).

}

}
\details{This function takes  3 inputs, all compulsory, i.e., data, order and parameters of the model to be fitted and returns the fitted values, residuals, RSS and AIC for the fitted model. No optimisation is done. This function is useful for simulation of data with given order and parameters.
}
\value{It returns the fitted EXPARMA model.}
\examples{
datats=c(17597,14074,11425,11691,11298,12351,14311,
12349,10537,11755,13154,11989,13022,12107,11172,10667,
10091,12204,12274,22343)
EXPARMAfit(datats,order=c(1,1), par=c(1,0.5,0.1,0.5,0.5,0.5))
}
