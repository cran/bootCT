% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_vecm_ardl.R
\name{sim_vecm_ardl}
\alias{sim_vecm_ardl}
\title{Generate data from a VECM/ARDL equation}
\usage{
sim_vecm_ardl(
  nobs,
  case = 1,
  sigma.in = diag(3),
  gamma.in,
  Axx.in,
  ayxUC.in,
  ayy.in,
  mu.in,
  eta.in,
  azeroy.in = 0,
  aoney.in = 0,
  burn.in,
  seed.in = NULL
)
}
\arguments{
\item{nobs}{number of observations.}

\item{case}{case related to intercept and trend}

\item{sigma.in}{error covariance matrix.}

\item{gamma.in}{list of short-run parameter matrices}

\item{Axx.in}{long-run relationships between the independent variables}

\item{ayxUC.in}{long-run unconditional relationship between dependent and independent variables, \eqn{\mathbf a_{yx}^{(UC)}} .
The second component ayxC, derived from conditioning, is calculated as\eqn{\mathbf a_{yx}^{(C)}= - \boldsymbol\omega'\mathbf A_{xx}}}

\item{ayy.in}{long-run relationship for the dependent variable \eqn{a_{yy}}}

\item{mu.in}{VAR intercept vector}

\item{eta.in}{VAR trend vector}

\item{azeroy.in}{Conditional ARDL intercept. Overridden if CASE I or CASE II}

\item{aoney.in}{Conditional ARDL trend. Overridden if CASE IV}

\item{burn.in}{burn-in number of observations}

\item{seed.in}{optional seed number for random error generation.}
}
\value{
A list that includes \itemize{
\item \code{dims}: a vector with the dataset dimension
\item \code{case}: the case given as input
\item \code{data}: the generated data
\item \code{diffdata}: the data first difference
\item \code{ut}: the generated random error matrix.
\item \code{sigma}: the error covariance matrix \eqn{\boldsymbol\Sigma}.
\item \code{omega}: the \eqn{\boldsymbol\omega} vector of parameters generated via conditioning
\item \code{At}: the conditional long-run parameter matrix \eqn{\tilde{\mathbf A}}
\item \code{ayx1}: the unconditional subvector of the ARDL equation \eqn{\mathbf a_{y.x}^{UC}}
\item \code{ayx}: the conditional subvector of the ARDL equation \eqn{a_{y.x}=a_{y.x}^{UC}-\omega'A_{xx}}
\item \code{gammalist}: the list of unconditional \eqn{\boldsymbol\Gamma_j} parameter matrices
\item \code{psilist}: the list of conditional \eqn{\boldsymbol\psi_{y.x,j}} parameter matrices
\item \code{azero}: the unconditional VECM intercept
\item \code{azero.c}: the conditional VECM intercept
\item \code{interc.ardl}: the conditional ARDL intercept
\item \code{aone}: the unconditional VECM trend
\item \code{aone.c}: the conditional VECM trend
\item \code{interc.ardl}: the conditional ARDL trend
\item \code{vmu}: the VAR intercept
\item \code{veta}: the VAR trend}
}
\description{
This function generates a dataset following a VECM/ARDL model equation.
}
\examples{
#PARAMETERS

#Sigma
corrm = matrix(0, ncol = 3, nrow = 3)
corrm[2,1] = 0.25
corrm[3,1] = 0.4
corrm[3,2] = -0.25
Corrm = (corrm + t(corrm)) + diag(3)
sds = diag(c(1.3, 1.2, 1))
Sigma = (sds \%*\% Corrm \%*\% t(sds))

#Gamma
gammax=list()
gammax[[1]] = matrix(c(0.6, 0, 0.2, 0.1, -0.3, 0, 0, -0.3, 0.2), nrow = 3, ncol = 3, byrow = TRUE)
gammax[[2]] = matrix(c(0.2, 0, 0.1, 0.05, -0.15, 0, 0, 0, 0.1), nrow = 3, ncol = 3, byrow = TRUE)

#DATA GENERATION
data_sim = sim_vecm_ardl(nobs = 200,
                         case = 3,
                         sigma.in = Sigma,
                         gamma.in = gammax,
                         Axx.in = matrix(c(0.3, 0.5, 0.4, 0.3), nrow = 2, ncol = 2),
                         ayxUC.in = c(0.5,0.6),
                         ayy.in = 0.7,
                         mu.in = rep(0.3, 3),
                         eta.in = rep(0, 3),
                         azeroy.in = 0.4,
                         aoney.in = 0,
                         burn.in = 50,
                         seed.in = 10)

}
