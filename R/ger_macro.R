#' Investment, Income and Consumption dataset.
#'
#' The data set contains quarterly, seasonally adjusted time series for 
#' West German fixed investment, disposable income, and consumption expenditures 
#' in billions of DM from 1960Q1 to 1982Q4. It was produced from file E1 of
#'  the data sets associated with Lutkepohl (2007). Originally obtained from Deutsche Bundesbank.
#'
#' @format A data frame with 92 rows and 4 variables:
#' \describe{
#'   \item{DATE}{Quarter}
#'   \item{INVEST}{Fixed investment (DM Billions)}
#'   \item{INCOME}{Disposable income (DM Billions)}
#'   \item{CONS}{Consumption expenditures (DM Billions)}
#' }
#' @source \url{http://www.jmulti.de/download/datasets/e1.dat}
"ger_macro"