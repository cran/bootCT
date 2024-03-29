% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boot_ardl.R
\name{boot_ardl}
\alias{boot_ardl}
\title{Bootstrap ARDL}
\usage{
boot_ardl(
  data,
  yvar = NULL,
  xvar = NULL,
  fix.ardl = NULL,
  info.ardl = "AIC",
  fix.vecm = NULL,
  info.vecm = "AIC",
  maxlag = 5,
  a.ardl = 0.05,
  a.vecm = 0.05,
  nboot = 2000,
  case = 3,
  a.boot.H0 = c(0.05, 0.025, 0.01),
  print = T
)
}
\arguments{
\item{data}{Input dataset. Must contain a dependent and a set of independent variables.}

\item{yvar}{Name of the dependent variable, enclosed in quotation marks. If NULL, the first variable will be used.}

\item{xvar}{Vector of names of the independent variables, each enclosed in quotation marks. If NULL, all variables except the first will be used.}

\item{fix.ardl}{Fixed lagged differences for the short term part of the ARDL equation.}

\item{info.ardl}{Selection criterion for the auto_ardl function. Options are "AIC", "AICc", BIC, "R2", "adjR2", if fix.ardl is null. Defaults to AIC.}

\item{fix.vecm}{Fixed lagged differences for the short term part of the VECM equation.}

\item{info.vecm}{Selection criterion for the VARselect function. Options are "AIC", "HQ", "SC", "FPE", if fix.vecm is null. Defaults to AIC.}

\item{maxlag}{Max number of lags for the auto_ardl  and VARselect procedures, if fix.ardl or fix.vecm are null, respectively.}

\item{a.ardl}{Threshold significance for the short-term ARDL coefficients significance.}

\item{a.vecm}{Threshold significance for the short-term VECM coefficients significance.}

\item{nboot}{Number of bootstrap replications.}

\item{case}{Model case, pertaining to the treatment of intercept and trend. Must be integer from 1 to 5. Defaults to 3.}

\item{a.boot.H0}{Probability/ies by which the critical quantiles of the bootstrap distribution(s) must be calculated.}

\item{print}{Show the progress bar.}
}
\value{
List of several elements including \itemize{
\item \code{data}: the data used to perform estimation and testing
\item \code{ARDL}: the estimated ARDL conditional model
\item \code{VECM}: the estimated VECM unconditional model
\item \code{jo.testX}: Johansen cointegration test on the independent variables
\item \code{pssbounds}: the PSS bound test output
\item \code{smgbounds}: the SMG bound test critical values
\item \code{fov.stat}: the test statistics on the conditional Fov tests
\item \code{t.stat}: the test statistics on the conditional t test
\item \code{find.stat}: the test statistics on the conditional Find tests
\item \code{quantfov}: the bootstrap conditional F Overall test critical value(s)
\item \code{quantt}: the bootstrap conditional t-test critical value(s)
\item \code{quantfind}: the bootstrap conditional F Independent test critical value(s)
\item \code{fakecoint}: indication of the situation in which \eqn{a_{y.x}\neq 0} but \eqn{a_{y.x}^{UC}=0}, signaling absence of cointegration.}
}
\description{
This is the main function of the package. It performs the bootstrap version of the ARDL bound test for cointegration.
}
\examples{
\dontrun{
#LOAD DATA
data("ger_macro")
# DATA PREPARATION
LNDATA = apply(ger_macro[,-1], 2, log)
col_ln = paste0("LN", colnames(ger_macro)[-1])
LNDATA = as.data.frame(LNDATA)
colnames(LNDATA) = col_ln
LNDATA$DATE = ger_macro$DATE
#ARDL BOOT 
BCT_res = boot_ardl(data = LNDATA,
                   yvar = "LNCONS",
                   xvar = c("LNINCOME","LNINVEST"),
                   maxlag = 5,
                   a.ardl = 0.1,
                   a.vecm = 0.1,
                   nboot = 2000,
                   case = 3,
                   a.boot.H0 = c(0.05),
                   print = T)
summary(boot_res)
}
}
