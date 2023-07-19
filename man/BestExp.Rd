\name{BestExp}
\alias{BestExp}

\title{Fitting of Best Exponential Autoregressive Moving Average (EXPARMA) Model}
\description{This function will fit best EXPARMA model.}

\usage{BestExp(ts_data, opt_method, max.p, max.q)}
%-
\arguments{
  \item{ts_data}{
a univarite time series data

}
\item{opt_method}{
opt_method represents the optimization methods that can be used. All the methods available at optim() function can be used here. By default the method used is "BFGS".}
\item{max.p}{
max.p is the order of autoregressive parameter of the model. By default the value is 5.}
\item{max.q}{
max.q is the order of moving average parameter of the model. By default the value is 5.}
}
\details{Takes only one input, the data. Checks the AIC values of all EXPARMA models from order (1,1) up to (max.p, max.q) by finding the optimal parameters of those models and returns a fit with the model having the least AIC. The maximum order up to which checks are set to 5 by default.
}
\value{It returns the best fitted EXPARMA model.}
\author{
Saikat Das, Bishal Gurung, Achal Lama and Kn Singh
}
\references{
Haggan and Ozaki (1981). Modelling nonlinear random vibrations using an amplitude-dependent autoregressive time series model. Biometrika,68(1):189-199.<doi:10.1093/biomet/68.1.189>.
Gurung (2015). An exponential autoregressive (EXPAR) (EXPAR) model for the forecasting of all India annual rainfall. Mausam,66(4):847-849. <doi:10.54302/mausam.v66i4.594>.
}
\examples{
datats=c(17597,14074,11425,11691,11298,12351,14311,
12349,10537,11755,13154,11989,13022,12107,11172,10667,
10091,12204,12274,22343)
BestExp(datats,opt_method="BFGS",max.p=1, max.q=1)
}

