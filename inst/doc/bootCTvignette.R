## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_chunk$set(echo = TRUE, fig.width = 7, fig.height = 3)

## ----setup, warning=FALSE-----------------------------------------------------
library(bootCT)

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
library(reshape2)
library(ggplot2)
library(dplyr)

## ----echo=TRUE----------------------------------------------------------------
sigma.in = matrix(c(1.69, 0.39, 0.52,
                    0.39, 1.44, -0.3,
                    0.52, -0.3, 1.00),3,3)

## ----echo=TRUE----------------------------------------------------------------
gamma1 = matrix(c(0.60, 0.00, 0.20,
                  0.10, -0.3, 0.00,
                  0.00, -0.3, 0.20),3,3,T)
gamma2 = gamma1 * 0.3
gamma.in = list(gamma1,gamma2)

## ----echo=TRUE----------------------------------------------------------------
ayy.in = 0.6
ayx.uc.in = c(0.4, 0.4)
axx.in = matrix(c(0.30, 0.50,
                  -0.4, 0.30),2,2,T)

## ----echo=TRUE----------------------------------------------------------------
case = 2
mu.in = c(2,2,2)
eta.in = c(0,0,0)
azero.in = c(0,0,0)
aone.in = c(0,0,0)

## ----echo=TRUE----------------------------------------------------------------
data.vecm.ardl_2 =
  sim_vecm_ardl(nobs=200,
                case = 2,
                sigma.in = sigma.in,
                gamma.in = gamma.in,
                axx.in = axx.in,
                ayx.uc.in = ayx.uc.in,
                ayy.in = ayy.in,
                mu.in = mu.in,
                eta.in = eta.in,
                azero.in = azero.in,
                aone.in = aone.in,
                burn.in = 100,
                seed.in = 999)

## ----echo=TRUE----------------------------------------------------------------
head(data.vecm.ardl_2$data)
head(data.vecm.ardl_2$diffdata)

## ----echo=TRUE----------------------------------------------------------------
df = data.vecm.ardl_2$data
meltdf <- melt(df,id="time")
ggplot(meltdf,
       aes(x = time, y = value, colour = variable, group = variable)) +
       geom_line() + ggtitle("CASE II - Level") + theme_bw()
df = data.vecm.ardl_2$diffdata
meltdf <- melt(df,id="time")
ggplot(meltdf,
       aes(x = time, y = value, colour = variable, group = variable)) +
       geom_line() + ggtitle("CASE II - Diff") + theme_bw()

## ----echo=TRUE----------------------------------------------------------------
# case III
case = 3
mu.in = c(0,0,0)
eta.in = c(0,0,0)
azero.in = c(2,2,2)
aone.in = c(0,0,0)

# case IV
case = 4
mu.in = c(0,0,0)
eta.in = c(0.6,0.6,0.6)
azero.in = c(2,2,2)
aone.in = c(0,0,0)

# case V
case = 5
mu.in = c(0,0,0)
eta.in = c(0,0,0)
azero.in = c(2,2,2)
aone.in = c(0.6,0.6,0.6)

## ----echo=TRUE----------------------------------------------------------------
data("ger_macro")

# Data preparation (log-transform)
LNDATA = apply(ger_macro[,-1], 2, log)
col_ln = paste0("LN", colnames(ger_macro)[-1])
LNDATA = as.data.frame(LNDATA)
colnames(LNDATA) = col_ln
LNDATA = LNDATA %>% select(LNCONS, LNINCOME, LNINVEST)
LNDATA$DATE = ger_macro$DATE

# First difference
lagdf1 = lag_mts(as.matrix(LNDATA[,-4]), k = c(1,1,1))
DIFF.LNDATA = na.omit(LNDATA[,-4] - lagdf1)
colnames(DIFF.LNDATA) = paste0("D_", colnames(LNDATA)[-4])
DIFF.LNDATA$DATE = ger_macro$DATE[-1]

# Graphs
dfmelt = melt(LNDATA, id = "DATE")
dfmelt = dfmelt%>%arrange(variable,DATE)
diff.dfmelt = melt(DIFF.LNDATA, id = "DATE")
diff.dfmelt = diff.dfmelt%>%arrange(variable,DATE)

ggplot(dfmelt,
       aes(x = DATE, y = value, colour = variable, group = variable)) +
  geom_line() + ggtitle("Level Variables (log-scale)") + theme_bw()

ggplot(diff.dfmelt,
       aes(x = DATE, y = value, colour = variable, group = variable)) +
  geom_line() + ggtitle("Diff. Variables (log-scale)") + theme_bw()

## ----results='hide'-----------------------------------------------------------
BCT_res_CONS = boot_ardl(data = LNDATA,
                         yvar = "LNCONS",
                         xvar = c("LNINCOME","LNINVEST"),
                         case = 3,
                         fix.ardl = NULL,
                         info.ardl = "AIC",
                         fix.vecm = NULL,
                         info.vecm = "AIC",
                         maxlag = 5,
                         a.ardl = 0.1,
                         a.vecm = 0.1,
                         nboot = 2000,
                         a.boot.H0 = c(0.01,0.05,0.1),
                         print = F)

## ----echo=TRUE----------------------------------------------------------------
summary(BCT_res_CONS, out="ARDL")

## ----echo=TRUE----------------------------------------------------------------
summary(BCT_res_CONS, out="VECM")

## ----echo=TRUE----------------------------------------------------------------
summary(BCT_res_CONS, out="cointVECM")

## ----echo=TRUE----------------------------------------------------------------
summary(BCT_res_CONS, out="cointARDL")

