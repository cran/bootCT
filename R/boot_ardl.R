#' Bootstrap ARDL
#'
#' This is the main function of the package. It performs the bootstrap version of the ARDL bound test for cointegration.
#'
#' @param data Input dataset. Must contain a dependent and a set of independent variables.
#' @param yvar Name of the dependent variable, enclosed in quotation marks. If NULL, the first variable will be used.
#' @param Xvar Vector of names of the independent variables, each enclosed in quotation marks. If NULL, all variables except the first will be used.
#' @param difflags Fixed lagged differences for the short term part of the ARDL equation.
#' @param maxlag Max number of lags for the auto_ardl procedure.
#' @param p.ardl Threshold p-value for the short-term ARDL coefficients significance in the bootstrap procedure.
#' @param p.vecm Threshold p-value for the short-term VECM coefficients significance in the bootstrap procedure.
#' @param B Number of bootstrap replications.
#' @param case Model case, pertaining to the treatment of intercept and trend. Must be integer from 1 to 5. Defaults to 3.
#' @param crit.H0 Probability/ies by which the critical quantiles of the bootstrap distribution(s) must be calculated.
#' @param print Show the progress bar.

#' @return List of several elements including \itemize{
#'\item \code{ARDL}: the conditional and unconditional ARDL models applied on the data
#'\item \code{pssbounds}: the PSS bound test output
#'\item \code{smgbounds}: the SMG bound test critical values
#'\item \code{fov.st}: the test statistics on the conditional and unconditional Fov tests
#'\item \code{t.st}: the test statistics on the conditional and unconditional t tests
#'\item \code{find.st}: the test statistics on the conditional and unconditional Find tests
#'\item \code{quantFOV}: the bootstrap conditional and unconditional F Overall test critical value(s)
#'\item \code{quantt}: the bootstrap conditional and unconditional t test critical value(s)
#'\item \code{quantFIND}: the bootstrap conditional and unconditional F Independent test critical value(s)}
#'
#' @examples
#'\dontrun{
#' data(ger_macro)
#' LNDATA = as.data.frame(log(ger_macro[,-1]))
#' colnames(LNDATA) = c("LNINVEST","LNINCOME","LNCONS")
#'
#' boot_res = boot_ardl(LNDATA, yvar = "LNINCOME", Xvar = c("LNCONS","LNINVEST"), maxlag = 5, B = 2000)
#' summary(boot_res)
#'}
#' @export
boot_ardl =
  function(data,
           yvar = NULL,
           Xvar = NULL,
           difflags = NULL,
           maxlag = 5,
           p.ardl = 0.05,
           p.vecm = 0.05,
           B = 2000,
           case = 3,
           crit.H0 = c(0.05, 0.025, 0.01),
           print = T) {

    if(!any(inherits(data,"data.frame"))){
      warning("Not a data frame object, converting")
      data=data.frame(data)
    }

    #SELECT ONLY NUMERIC COLUMNS AND ARRANGE y AND X
    dfnames = colnames(data)
    if (is.null(yvar)) {
      yvar = dfnames[1]
    }
    if (is.null(Xvar)) {
      Xvar = dfnames[-1]
    }

    if (!(yvar %in% dfnames)) {
      stop("y variable not in dataset")
    }

    if (!all(Xvar %in% dfnames)) {
      stop("at least one of the X variables not in dataset.")
    }

    if (yvar %in% Xvar) {
      stop("y variable and X variables must be separated.")
    }

    if (p.ardl < 0 | p.ardl > 1) {
      warning("Invalid ARDL threshold p.value. Defaulting to 0.05.")
      p.ardl = 0.05
    }

    if (p.vecm < 0 | p.vecm > 1) {
      warning("Invalid VECM threshold p.value. Defaulting to 0.05.")
      p.vecm = 0.05
    }

    if (length(p.ardl) > 1) {
      warning("Expecting a scalar value. Defaulting to 0.05.")
      p.ardl = 0.05
    }

    if (length(p.vecm) > 1) {
      warning("Expecting a scalar value. Defaulting to 0.05.")
      p.vecm = 0.05
    }

    if (any(difflags %% 1 != 0)) {
      warning("Expecting integer values. Rounding to the nearest integer.")
      difflags = round(difflags)
    }

    if (any(maxlag %% 1 != 0)) {
      warning("Expecting integer value. Rounding to the nearest integer.")
      maxlag = round(maxlag)
    }

    if (B %% 1 != 0) {
      warning("Expecting integer value. Rounding to the nearest integer.")
      B = round(B)
    }

    if (case < 1 | case > 5 | case %% 1 != 0) {
      warning("Inadmissible case. Defaulting to case III.")
      case = 3
    }

    if (any(crit.H0 < 0) | any(crit.H0 > 1)) {
      warning("Invalid probability for critical values. Defaulting to 1%, 2.5% and 5%")
      crit.H0 = c(0.01, 0.025, 0.05)
    }

    if(B<200){
      warning("suppressing print: B should be at least 200")
    }

    numdata = data %>% dplyr::select(yvar, Xvar)
    d = ncol(numdata)

    #LAGGED AND FIRST DIFFERENCE
    lagdata = lag_mts(as.matrix(numdata), k = rep(1, d))
    dlag0 = numdata - lagdata
    colnames(dlag0) <- paste("D_", colnames(numdata), sep = "")

    if (is.null(difflags)) {
      #INITIAL ARDL MODEL DATAFRAME
      df.ardl = cbind(lagdata, dlag0)

      #ID FOR TREND
      df.ardl$timeid = 1:nrow(df.ardl)

      #FORMULA BASED ON CASE
      text.ardl = (paste(colnames(df.ardl)[(1 + d)], "~", paste(colnames(df.ardl)[-(1 +
                                                                                      d)], collapse = '+')))

      if (case == 1) {
        formula.ardl = as.formula(paste(substr(text.ardl, 1, nchar(text.ardl) -
                                                 7), "-1"))
        fixed = c(-1, rep(0, d), rep(-1, (d - 1)))
      }

      if (case == 2) {
        formula.ardl = as.formula(paste(substr(text.ardl, 1, nchar(text.ardl) -
                                                 7)))
        fixed = c(-1, rep(0, d), rep(-1, (d - 1)))
      }

      if (case == 3) {
        formula.ardl = as.formula(paste(substr(text.ardl, 1, nchar(text.ardl) - 7)))
        fixed = c(-1, rep(0, d), rep(-1, (d - 1)))
      }

      if (case >= 4) {
        formula.ardl = as.formula(text.ardl)
        fixed = c(-1, rep(0, d), rep(-1, (d - 1)), 0)
      }

      #SELECT BEST ARDL ORDER
      ardl_select = ARDL::auto_ardl(
        formula.ardl,
        data = na.omit(df.ardl),
        max_order = maxlag,
        fixed_order = fixed
      )
      best_ardl = ardl_select$best_order[c(1, (d + 2):(2 * d))]

      #CONSTRUCT FINAL ARDL DATAFRAME
      final_dlag = lag_mts(dlag0, k = best_ardl)
    } else{
      if (length(difflags) != d) {
        warning("number of variable is not equal to number of difflags, taking only the first value")
        difflags0 = rep(difflags[1], d)
        difflags = difflags0
      }
      final_dlag = lag_mts(dlag0, k = difflags)
    }

    df.ardl.l = list(df.ardl.c = na.omit(cbind(dlag0[, 1], lagdata, final_dlag, dlag0[, -1])),
                     df.ardl.uc = na.omit(cbind(dlag0[, 1], lagdata, final_dlag)))

    #EMPTY LISTS FOR OUTPUT
    mod.ardl.l = beta.ardl.l = omega.ardl.l = e.ardl.l = sigmae.ardl.l =
      Sigmabeta.ardl.l =
      F.overall.os = t.dep.os = F.ind.os =
      pss.b = pss.f01 = pss.f05 = pss.f10 =
      pss.t01 = pss.t05 = pss.t10 =
      pss.fi01 = pss.fi02 = pss.fi05 = pss.fi10 = list()

    #EMPTY OBJECTS FOR BOOT PROCEDURE
    F.overallc = t.depc = F.indc = matrix(0, nrow = 2, ncol = B)
    quantFOV = quantFIND = quantt = matrix(0, nrow = length(crit.H0), ncol = 2)

    for (cond in 1:2) {
      formula.ardl = ardl_select$best_model$full_formula
      conditional = ifelse(cond == 1, "Conditional", "Unconditional")
      colnames(df.ardl.l[[cond]])[1] = colnames(dlag0)[1]

      #FINAL ARDL ESTIMATE
      if (case >= 4) {
        df.ardl.l[[cond]]$timeid = df.ardl$timeid
      }

      formula.ardl = ardl_to_lm(formula.ardl, case = case, d = d)

      if (cond == 2) {
        formula.ardl = strip_formula_UC(formula.ardl,case,d)
      }

      #ESTIMATE UNRESTRICTED ARDL MODEL
      mod.ardl.l[[cond]] = lm(formula.ardl, data = na.omit(df.ardl.l[[cond]]))
      beta.ardl.l[[cond]] = coef(mod.ardl.l[[cond]])

      #OMEGA DEPENDS ON CONDITIONING
      pos.omega = ((length(beta.ardl.l[[cond]]) - d + 2):(length(beta.ardl.l[[cond]]))) -
        1 * (case >= 4)
      omega.ardl.l[[cond]] = (beta.ardl.l[[cond]][pos.omega]) * (cond == 1)

      #RESIDUALS
      e.ardl.l[[cond]] = residuals(mod.ardl.l[[cond]])
      sigmae.ardl.l[[cond]] = summary(mod.ardl.l[[cond]])$sigma ^ 2

      Sigmabeta.ardl.l[[cond]] = summary(mod.ardl.l[[cond]])$cov.unscaled *
        pracma::repmat((sigmae.ardl.l[[cond]]), length(coef(mod.ardl.l[[cond]])))

      #TEST STATISTICS RESTRICTIONS BASED ON CASE
      if (case == 1) {
        terms.ur = 1:d
      }
      if (case == 2) {
        terms.ur = 1:(d + 1)
      }
      if (case == 3) {
        terms.ur = 2:(d + 1)
      }
      if (case == 4) {
        terms.ur = c(2:(d + 1), length(beta.ardl.l[[cond]]))
      }
      if (case == 5) {
        terms.ur = 2:(d + 1)
      }

      wtestfo = aod::wald.test(Sigmabeta.ardl.l[[cond]], beta.ardl.l[[cond]], Terms = terms.ur)
      F.overall.os[[cond]] = wtestfo$result$chi2[1] / wtestfo$result$chi2[2]
      t.dep.os[[cond]] = beta.ardl.l[[cond]][(1 + 1 * (case > 1))] / sqrt(Sigmabeta.ardl.l[[cond]][(1 +
                                                                                                      1 * (case > 1)), (1 + 1 * (case > 1))])
      wtestfi = aod::wald.test(Sigmabeta.ardl.l[[cond]], beta.ardl.l[[cond]], Terms = ((2:d) +
                                                                                         ((case > 1) * 1)))
      F.ind.os[[cond]] = wtestfi$result$chi2[1] / wtestfi$result$chi2[2]

      num = nrow(df.ardl.l[[cond]])

      #PSS BOUND TEST
      pss.b[[cond]] = dynamac::pssbounds(
        obs = num,
        fstat = F.overall.os[[cond]],
        tstat = t.dep.os[[cond]],
        case = case,
        k = d - 1,
        object.out = T
      )

      smk_crit = bootCT::smk_crit
      smk_test = NA

      nsmk = 81
      #SMK BOUND TEST ON INDEPENDENT
      if (case %in% c(1, 3, 5)) {
        case0 = case
        if (num <= 80) {
          rem = (num - smk_crit$num[1:12])
          check.nsmk = smk_crit$num[which.min(rem[rem > 0])]
          if (length(check.nsmk) == 0) {
            nsmk = 30
          }else{
            nsmk = check.nsmk
          }
        }

        cols_sel = which(as.numeric(substr(
          colnames(smk_crit), nchar(colnames(smk_crit)), nchar(colnames(smk_crit))
        )[-c(1, 2, 3)]) == d - 1) + 3
        smk_test0 = smk_crit %>% dplyr::filter(case == case0, num == nsmk)
        smk_test0 = smk_test0[, c(1, 2, 3, cols_sel)]
        smk_test=smk_test0[order(-smk_test0$prob),]
      }

      pss.f01[[cond]] = c(pss.b[[cond]]$ftest.I0.p01, pss.b[[cond]]$ftest.I1.p01)
      pss.f05[[cond]] = c(pss.b[[cond]]$ftest.I0.p05, pss.b[[cond]]$ftest.I1.p05)
      pss.f10[[cond]] = c(pss.b[[cond]]$ftest.I0.p10, pss.b[[cond]]$ftest.I1.p10)

      if (case %in% c(1, 3, 5)) {
        #T AND FIND TESTS ON CASES I, III, V
        pss.t01[[cond]] = c(pss.b[[cond]]$ttest.I0.p01, pss.b[[cond]]$ttest.I1.p01)
        pss.t05[[cond]] = c(pss.b[[cond]]$ttest.I0.p05, pss.b[[cond]]$ttest.I1.p05)
        pss.t10[[cond]] = c(pss.b[[cond]]$ttest.I0.p10, pss.b[[cond]]$ttest.I1.p10)
        pss.fi01[[cond]] = c(smk_test[1, 3], smk_test[1, 4])
        pss.fi02[[cond]] = c(smk_test[2, 3], smk_test[2, 4])
        pss.fi05[[cond]] = c(smk_test[3, 3], smk_test[3, 4])
        pss.fi10[[cond]] = c(smk_test[4, 3], smk_test[4, 4])
      }

      for (sh in c("fov", "t", "find")) {
        H0=sh
        #F OVERALL IMPOSING THE NULL
        if (H0 == "fov") {
          formula.ardl.H0 = strip_formula_H0(formula.ardl, H0, case, d)
          mod.ardl.H0 = lm(formula.ardl.H0, data = df.ardl.l[[cond]])

          beta.ardl.H0 = coef(mod.ardl.H0) #RESTRICTED ESTIMATES

          #OMEGA IF CONDITIONAL
          if (cond == 1) {
            omega.ardl.H0 = beta.ardl.H0[((length(beta.ardl.H0) - d + 2):length(beta.ardl.H0)) -
                                           1 * (case == 5)]
          } else{
            omega.ardl.H0 = rep(0, d - 1)
          }

          #RESTRICTED RESIDUALS
          e.ardl.H0 = residuals(mod.ardl.H0)

          typevecm = ifelse(case == 1, "none", ifelse(case < 4, "const", "both"))

          #SELECT UNRESTRICTED VECM
          vecm = vars::VARselect(
            na.omit(dlag0),
            lag.max = maxlag,
            type = typevecm,
            exogen = na.omit(lagdata)
          )

          vecmsel = min(vecm$selection)

          #CONFIRM UNRESTRICTED VECM
          vecm = vars::VAR(
            na.omit(dlag0),
            p = vecmsel,
            type = typevecm,
            exogen = na.omit(lagdata)
          )

          nx = nrow(na.omit(dlag0))
          dx = d * vecmsel + d * (d - 1)

          #IMPOSE WEAK EXOGENEITY AND REMOVE
          #NON-SIGNIFICANT SHORT-RUN ESTIMATES
          vecm.ser = vars::restrict(vecm,
                                    method = "ser",
                                    thresh = qnorm(1 - p.vecm))
          mat.ser = vecm.ser$restrictions
          mat.ser[, ncol(mat.ser) - d + 1] = c(1, rep(0, d - 1))
          mat.ser[, ncol(mat.ser) - d] = 1
          mat.ser[, (ncol(mat.ser) - d + 2):ncol(mat.ser)] = 1
          vecm.ser = vars::restrict(vecm, method = "manual", resmat = mat.ser)

          #VECM RESIDUALS
          e.vecm = residuals(vecm.ser)
          if (length(e.ardl.H0) > nrow(e.vecm)) {
            diffrow = length(e.ardl.H0) - nrow(e.vecm)
            e.z = cbind(e.ardl.H0[-c(1:diffrow)], e.vecm[,-1])
          } else{
            if (length(e.ardl.H0) < nrow(e.vecm)) {
              diffrow = nrow(e.vecm) - length(e.ardl.H0)
              e.z = cbind(e.ardl.H0, e.vecm[-c(1:diffrow),-1])
            } else{
              e.z = cbind(e.ardl.H0, e.vecm[,-1])
            }
          }

          #STORE PARAMETERS
          GAMMAX = do.call(cbind, vars::Acoef(vecm.ser))
          A = -vars::Bcoef(vecm.ser)[, (ncol(mat.ser) - d + 1):ncol(mat.ser)]

          azero = auno = rep(0, d)

          if (case > 1) {
            azero = vars::Bcoef(vecm.ser)[, (ncol(mat.ser) - d - (1 * (case > 3)))] #VECM INTERCEPT
          }

          if (case > 3) {
            auno = vars::Bcoef(vecm.ser)[, (ncol(mat.ser) - d)] #VECM TREND
          }

          #RESTRICTED BETA DEPENDING ON CASE
          if (case == 1) {
            hbeta.ardl.H0 = beta.ardl.H0[1:length(beta.ardl.H0)] *
              (summary(mod.ardl.H0)[[4]][, 4] < p.ardl)[1:length(beta.ardl.H0)]
            join_nameG = which(names(GAMMAX[1,]) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            join_nameB = which(names(hbeta.ardl.H0) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            GAMMAX[1,-join_nameG] = 0
            GAMMAX[1, join_nameG] = hbeta.ardl.H0[join_nameB]
          }

          if (case == 2) {
            azero[1] = 0
            #INTERCEPT CONSTRAINED TO VECM
            hbeta.ardl.H0 = beta.ardl.H0[1:length(beta.ardl.H0)] *
              (summary(mod.ardl.H0)[[4]][, 4] < p.ardl)[1:length(beta.ardl.H0)]
            join_nameG = which(names(GAMMAX[1,]) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            join_nameB = which(names(hbeta.ardl.H0) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            GAMMAX[1,-join_nameG] = 0
            GAMMAX[1, join_nameG] = hbeta.ardl.H0[join_nameB]
          }

          if (case == 3) {
            azero[1] = beta.ardl.H0[1] #ARDL UNCONSTRAINED INTERCEPT
            hbeta.ardl.H0 = beta.ardl.H0[2:length(beta.ardl.H0)] *
              (summary(mod.ardl.H0)[[4]][, 4] < p.ardl)[2:length(beta.ardl.H0)]
            join_nameG = which(names(GAMMAX[1,]) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            join_nameB = which(names(hbeta.ardl.H0) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            GAMMAX[1,-join_nameG] = 0
            GAMMAX[1, join_nameG] = hbeta.ardl.H0[join_nameB]
          }

          if (case == 4) {
            #TREND CONSTRAINED TO VECM
            azero[1] = beta.ardl.H0[1] #ARDL UNCONSTRAINED INTERCEPT
            auno[1] = 0
            hbeta.ardl.H0 = beta.ardl.H0[2:length(beta.ardl.H0)] *
              (summary(mod.ardl.H0)[[4]][, 4] < p.ardl)[2:length(beta.ardl.H0)]
            join_nameG = which(names(GAMMAX[1,]) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            join_nameB = which(names(hbeta.ardl.H0) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            GAMMAX[1,-join_nameG] = 0
            GAMMAX[1, join_nameG] = hbeta.ardl.H0[join_nameB]
          }

          if (case == 5) {
            azero[1] = beta.ardl.H0[1] #ARDL UNCONSTRAINED INTERCEPT
            auno[1] = beta.ardl.H0[length(beta.ardl.H0)] #ARDL UNCONSTRAINED TREND
            hbeta.ardl.H0 = beta.ardl.H0[2:(length(beta.ardl.H0) - 1)] *
              (summary(mod.ardl.H0)[[4]][, 4] < p.ardl)[2:(length(beta.ardl.H0) -
                                                             1)]
            join_nameG = which(names(GAMMAX[1,]) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            join_nameB = which(names(hbeta.ardl.H0) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            GAMMAX[1,-join_nameG] = 0
            GAMMAX[1, join_nameG] = hbeta.ardl.H0[join_nameB]
          }

          #F OVERALL
          A[1, 1] = 0
          ay.x = rep(0, d - 1)
          A[1,-1] = ay.x
          A[-1, 1] = rep(0, d - 1)
        }

        if (H0 == "t") {
          formula.ardl.H0 = strip_formula_H0(formula.ardl, H0, case,d)

          mod.ardl.H0 = lm(formula.ardl.H0, data = df.ardl.l[[cond]]) #case III
          summary(mod.ardl.H0)

          beta.ardl.H0 = coef(mod.ardl.H0)

          if (cond == 1) {
            omega.ardl.H0 = beta.ardl.H0[((length(beta.ardl.H0) - d + 2):length(beta.ardl.H0)) -
                                           1 * (case == 5)]
          } else{
            omega.ardl.H0 = rep(0, d - 1)
          }

          e.ardl.H0 = residuals(mod.ardl.H0)

          typevecm = ifelse(case == 1, "none", ifelse(case < 4, "const", "both"))

          vecm = vars::VARselect(
            na.omit(dlag0),
            lag.max = maxlag,
            type = typevecm,
            exogen = na.omit(lagdata)
          )

          vecmsel = min(vecm$selection)

          vecm = vars::VAR(
            na.omit(dlag0),
            p = vecmsel,
            type = typevecm,
            exogen = na.omit(lagdata)
          )

          nx = nrow(na.omit(dlag0))
          dx = d * vecmsel + d * (d - 1)

          vecm.ser = vars::restrict(vecm,
                                    method = "ser",
                                    thresh = qnorm(1 - p.vecm))
          mat.ser = vecm.ser$restrictions
          mat.ser[, ncol(mat.ser) - d + 1] = c(1, rep(0, d - 1))
          mat.ser[, ncol(mat.ser) - d] = 1
          mat.ser[, (ncol(mat.ser) - d + 2):ncol(mat.ser)] = 1
          vecm.ser = vars::restrict(vecm, method = "manual", resmat = mat.ser)

          summary(vecm.ser)

          e.vecm = residuals(vecm.ser)
          if (length(e.ardl.H0) > nrow(e.vecm)) {
            diffrow = length(e.ardl.H0) - nrow(e.vecm)
            e.z = cbind(e.ardl.H0[-c(1:diffrow)], e.vecm[,-1])
          } else{
            if (length(e.ardl.H0) < nrow(e.vecm)) {
              diffrow = nrow(e.vecm) - length(e.ardl.H0)
              e.z = cbind(e.ardl.H0, e.vecm[-c(1:diffrow),-1])
            } else{
              e.z = cbind(e.ardl.H0, e.vecm[,-1])
            }
          }
          GAMMAX = do.call(cbind, vars::Acoef(vecm.ser))
          A = -vars::Bcoef(vecm.ser)[, (ncol(mat.ser) - d + 1):ncol(mat.ser)]

          azero = auno = rep(0, d)

          if (case > 1) {
            azero = vars::Bcoef(vecm.ser)[, (ncol(mat.ser) - d - (1 * (case > 3)))] #interc del vecm
          }

          if (case > 3) {
            auno = vars::Bcoef(vecm.ser)[, (ncol(mat.ser) - d)] #trend del vecm
          }

          if (case == 1) {
            ay.x = -beta.ardl.H0[1:(d - 1)]
            hbeta.ardl.H0 = beta.ardl.H0[d:length(beta.ardl.H0)] *
              (summary(mod.ardl.H0)[[4]][, 4] < p.ardl)[d:length(beta.ardl.H0)]
            join_nameG = which(names(GAMMAX[1,]) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            join_nameB = which(names(hbeta.ardl.H0) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            GAMMAX[1,-join_nameG] = 0
            GAMMAX[1, join_nameG] = hbeta.ardl.H0[join_nameB]
          }
          if (case == 2) {
            azero[1] = 0
            ay.x = -beta.ardl.H0[1:(d - 1)]
            hbeta.ardl.H0 = beta.ardl.H0[d:length(beta.ardl.H0)] *
              (summary(mod.ardl.H0)[[4]][, 4] < p.ardl)[d:length(beta.ardl.H0)]
            join_nameG = which(names(GAMMAX[1,]) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            join_nameB = which(names(hbeta.ardl.H0) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            GAMMAX[1,-join_nameG] = 0
            GAMMAX[1, join_nameG] = hbeta.ardl.H0[join_nameB]
          }
          if (case == 3) {
            azero[1] = beta.ardl.H0[1]
            ay.x = -beta.ardl.H0[2:d]
            hbeta.ardl.H0 = beta.ardl.H0[(d+1):length(beta.ardl.H0)] *
              (summary(mod.ardl.H0)[[4]][, 4] < p.ardl)[(d+1):length(beta.ardl.H0)]
            join_nameG = which(names(GAMMAX[1,]) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            join_nameB = which(names(hbeta.ardl.H0) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            GAMMAX[1,-join_nameG] = 0
            GAMMAX[1, join_nameG] = hbeta.ardl.H0[join_nameB]
          }
          if (case == 4) {
            azero[1] = beta.ardl.H0[1]
            auno[1] = 0
            ay.x = -beta.ardl.H0[2:d]
            hbeta.ardl.H0 = beta.ardl.H0[(d+1):length(beta.ardl.H0)] *
              (summary(mod.ardl.H0)[[4]][, 4] < p.ardl)[(d+1):length(beta.ardl.H0)]
            join_nameG = which(names(GAMMAX[1,]) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            join_nameB = which(names(hbeta.ardl.H0) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            GAMMAX[1,-join_nameG] = 0
            GAMMAX[1, join_nameG] = hbeta.ardl.H0[join_nameB]
          }
          if (case == 5) {
            azero[1] = beta.ardl.H0[1]
            auno[1] = beta.ardl.H0[length(beta.ardl.H0)]
            ay.x = -beta.ardl.H0[2:d]
            hbeta.ardl.H0 = beta.ardl.H0[(d+1):(length(beta.ardl.H0) - 1)] *
              (summary(mod.ardl.H0)[[4]][, 4] < p.ardl)[(d+1):(length(beta.ardl.H0) - 1)]
            join_nameG = which(names(GAMMAX[1,]) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            join_nameB = which(names(hbeta.ardl.H0) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            GAMMAX[1,-join_nameG] = 0
            GAMMAX[1, join_nameG] = hbeta.ardl.H0[join_nameB]
          }

          #T TEST
          A[1, 1] = 0
          A[1,-1] = ay.x
          A[-1, 1] = rep(0, d - 1)
        }

        if (H0 == "find") {

          formula.ardl.H0 = strip_formula_H0(formula.ardl, H0, case,d)

          mod.ardl.H0 = lm(formula.ardl.H0, data = df.ardl.l[[cond]]) #case III
          summary(mod.ardl.H0)

          beta.ardl.H0 = coef(mod.ardl.H0)

          if (cond == 1) {
            omega.ardl.H0 = beta.ardl.H0[((length(beta.ardl.H0) - d + 2):length(beta.ardl.H0)) -
                                           (1 * case == 5)]
          } else{
            omega.ardl.H0 = rep(0, d - 1)
          }

          e.ardl.H0 = residuals(mod.ardl.H0)

          typevecm = ifelse(case == 1, "none", ifelse(case < 4, "const", "both"))

          vecm = vars::VARselect(
            na.omit(dlag0),
            lag.max = maxlag,
            type = typevecm,
            exogen = na.omit(lagdata)
          )

          vecmsel = min(vecm$selection)

          vecm = vars::VAR(
            na.omit(dlag0),
            p = vecmsel,
            type = typevecm,
            exogen = na.omit(lagdata)
          )

          summary(vecm)
          nx = nrow(na.omit(dlag0))
          dx = d * vecmsel + d * (d - 1)

          vecm.ser = vars::restrict(vecm,
                              method = "ser",
                              thresh = qnorm(1 - p.vecm))
          mat.ser = vecm.ser$restrictions
          mat.ser[, ncol(mat.ser) - d + 1] = c(1, rep(0, d - 1))
          mat.ser[, ncol(mat.ser) - d] = 1
          mat.ser[, (ncol(mat.ser) - d + 2):ncol(mat.ser)] = 1
          vecm.ser = vars::restrict(vecm, method = "manual", resmat = mat.ser)

          summary(vecm.ser)

          e.vecm = residuals(vecm.ser)
          if (length(e.ardl.H0) > nrow(e.vecm)) {
            diffrow = length(e.ardl.H0) - nrow(e.vecm)
            e.z = cbind(e.ardl.H0[-c(1:diffrow)], e.vecm[,-1])
          } else{
            if (length(e.ardl.H0) < nrow(e.vecm)) {
              diffrow = nrow(e.vecm) - length(e.ardl.H0)
              e.z = cbind(e.ardl.H0, e.vecm[-c(1:diffrow),-1])
            } else{
              e.z = cbind(e.ardl.H0, e.vecm[,-1])
            }
          }
          GAMMAX = do.call(cbind, vars::Acoef(vecm.ser))
          A = -vars::Bcoef(vecm.ser)[, (ncol(mat.ser) - d + 1):ncol(mat.ser)]

          azero = auno = rep(0, d)

          if (case > 1) {
            azero = vars::Bcoef(vecm.ser)[, (ncol(mat.ser) - d - (1 * (case > 3)))]
          }

          if (case > 3) {
            auno = vars::Bcoef(vecm.ser)[, (ncol(mat.ser) - d)]
          }

          if (case == 1) {
            ayy = -beta.ardl.H0[1]
            hbeta.ardl.H0 = beta.ardl.H0[2:length(beta.ardl.H0)] *
              (summary(mod.ardl.H0)[[4]][, 4] < p.ardl)[2:length(beta.ardl.H0)]
            join_nameG = which(names(GAMMAX[1,]) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            join_nameB = which(names(hbeta.ardl.H0) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            GAMMAX[1,-join_nameG] = 0
            GAMMAX[1, join_nameG] = hbeta.ardl.H0[join_nameB]
          }
          if (case == 2) {
            azero[1] = 0
            ayy = -beta.ardl.H0[1]
            hbeta.ardl.H0 = beta.ardl.H0[2:length(beta.ardl.H0)] *
              (summary(mod.ardl.H0)[[4]][, 4] < p.ardl)[2:length(beta.ardl.H0)]
            join_nameG = which(names(GAMMAX[1,]) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            join_nameB = which(names(hbeta.ardl.H0) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            GAMMAX[1,-join_nameG] = 0
            GAMMAX[1, join_nameG] = hbeta.ardl.H0[join_nameB]
          }
          if (case == 3) {
            azero[1] = beta.ardl.H0[1]
            ayy = -beta.ardl.H0[2]
            hbeta.ardl.H0 = beta.ardl.H0[3:length(beta.ardl.H0)] *
              (summary(mod.ardl.H0)[[4]][, 4] < p.ardl)[3:length(beta.ardl.H0)]
            join_nameG = which(names(GAMMAX[1,]) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            join_nameB = which(names(hbeta.ardl.H0) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            GAMMAX[1,-join_nameG] = 0
            GAMMAX[1, join_nameG] = hbeta.ardl.H0[join_nameB]
          }
          if (case == 4) {
            azero[1] = beta.ardl.H0[1]
            auno[1] = 0
            ayy = -beta.ardl.H0[2]
            hbeta.ardl.H0 = beta.ardl.H0[3:length(beta.ardl.H0)] *
              (summary(mod.ardl.H0)[[4]][, 4] < p.ardl)[3:length(beta.ardl.H0)]
            join_nameG = which(names(GAMMAX[1,]) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            join_nameB = which(names(hbeta.ardl.H0) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            GAMMAX[1,-join_nameG] = 0
            GAMMAX[1, join_nameG] = hbeta.ardl.H0[join_nameB]
          }
          if (case == 5) {
            azero[1] = beta.ardl.H0[1]
            auno[1] = beta.ardl.H0[length(beta.ardl.H0)]
            ayy = -beta.ardl.H0[2]
            hbeta.ardl.H0 = beta.ardl.H0[3:(length(beta.ardl.H0) - 1)] *
              (summary(mod.ardl.H0)[[4]][, 4] < p.ardl)[3:(length(beta.ardl.H0) -
                                                             1)]
            join_nameG = which(names(GAMMAX[1,]) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            join_nameB = which(names(hbeta.ardl.H0) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            GAMMAX[1,-join_nameG] = 0
            GAMMAX[1, join_nameG] = hbeta.ardl.H0[join_nameB]
          }

          #F INDEP
          ay.x = rep(0, d - 1)
          A[1, 1] = ayy
          A[1,-1] = ay.x
          A[-1, 1] = rep(0, d - 1)
        }

        #BOOT PROCEDURE
        for (b in 1:B) {
          #RESIDUAL RESAMPLING
          resid.b = (dplyr::sample_n(
            as.data.frame(e.z),
            size = nrow(numdata) + vecmsel + 1,
            replace = T
          ))

          #CENTER RESIDUALS
          meanmat = matrix(rep(apply(resid.b, 2, mean), nrow(resid.b)), ncol = d, byrow = T)
          resid.b = resid.b - meanmat

          #TAKE FIRST OBSERVATIONS
          mat.start = as.matrix(numdata)
          lagdfm = dlagdfm = list()
          for (lsel in 1:(vecmsel + 1)) {
            lagdfm[[lsel]] = lag_mts(as.matrix(numdata),
                                     k = rep(lsel, d),
                                     last.only = T)
            if (lsel == 1) {
              dlagdfm[[lsel]] = numdata - lagdfm[[(lsel)]]
            } else{
              dlagdfm[[lsel]] = lagdfm[[(lsel - 1)]] - lagdfm[[(lsel)]]
            }
          }

          dlag_all = do.call(cbind, dlagdfm)
          colnames(dlag_all) = paste0("D_", colnames(dlag_all))
          mat.start = cbind(numdata, lagdfm[[1]], dlag_all, 1)
          start.z = gtools::na.replace(mat.start[1:(vecmsel + 1),], 0)

          #GENERATE DATA UNDER THE NULL
          boot.data = boot_ardl_c(
            r_in = as.matrix(resid.b),
            GAMMAX = GAMMAX,
            A = A,
            interc = azero,
            trend = auno,
            start_z = as.matrix(start.z),
            omegat = omega.ardl.H0
          )$df[, 1:d]

          #DATASET NAMES
          typev = c("dep", "ind")
          nind = 1:(ncol(resid.b) - 1)
          nlag = 0:vecmsel
          lev.df = data.frame(names = apply(
            expand.grid(typev, nind, nlag), 1, paste, collapse = "_"
          )) %>%
            dplyr::filter(!stringr::str_detect(names, "^dep_[2-9]"))
          lev.col = as.character(lev.df[1:(d * 2),])
          diff.df = data.frame(names = apply(
            expand.grid(typev, nind, nlag), 1, paste, collapse = "_"
          )) %>%
            dplyr::filter(!stringr::str_detect(names, "^dep_[2-9]"))
          d.col = paste0("d_",
                         as.character(
                            diff.df[1:(d * (vecmsel + 1)),]
                         ))
          cnames = c(lev.col, d.col, "const")
          colnames(boot.data) = c("y", paste0("x", 1:(d - 1)))

          #LAGS AND DIFFERENCES
          lagboot1 = lag_mts(as.matrix(boot.data), k = rep(1, d))
          dlagboot0 = boot.data - lagboot1
          colnames(dlagboot0) <-
            paste("D_", colnames(boot.data), sep = "")

          #BOOT MODEL
          boot.ardl = na.omit(cbind(lagboot1, dlagboot0))
          final_dlag = lag_mts(dlagboot0, k = best_ardl)

          if (cond == 1) {
            boot.ardl = na.omit(cbind(
              dlagboot0[, 1],
              lagboot1,
              final_dlag,
              dlagboot0[, -1],
              timeid = 1:nrow(lagboot1)
            ))
          } else{
            boot.ardl = na.omit(cbind(
              dlagboot0[, 1],
              lagboot1,
              final_dlag,
              timeid = 1:nrow(lagboot1)
            ))
          }

          colnames(boot.ardl)[1] = "D_y"

          #INTERCEPT AND TREND BASED ON CASE
          if (case == 1) {
            formula.ardl.b = as.formula("D_y ~ . -1-timeid")
            terms.ur = 1:d
          }

          if (case == 2) {
            formula.ardl.b = as.formula("D_y ~ . -timeid")
            terms.ur = 1:(d + 1)
          }

          if (case == 3) {
            formula.ardl.b = as.formula("D_y ~ . -timeid")
            terms.ur = 2:(d + 1)
          }

          if (case == 4) {
            formula.ardl.b = as.formula("D_y ~ .")
            terms.ur = c(2:(d + 1), length(beta.ardl.l[[cond]]))
          }

          if (case == 5) {
            formula.ardl.b = as.formula("D_y ~ .")
            terms.ur = 2:(d + 1)
          }

          #MODEL UNDER THE NULL
          mod.ardl.b = lm(formula.ardl.b, data = boot.ardl)
          summary(mod.ardl.b)
          beta.ardl.b = coef(mod.ardl.b)
          e.ardl.b = residuals(mod.ardl.b)
          sigmae.ardl.b = summary(mod.ardl.b)$sigma ^ 2
          Sigmabeta.ardl.b = summary(mod.ardl.b)$cov.unscaled * pracma::repmat((sigmae.ardl.b), length(coef(mod.ardl.b)))

          #BOOT STATISTICS
          if (H0 == "fov") {
            wtestfo = aod::wald.test(Sigmabeta.ardl.b, beta.ardl.b, Terms = terms.ur)
            F.overallc[cond, b] = wtestfo$result$chi2[1] / wtestfo$result$chi2[2]
          }
          if (H0 == "t") {
            t.depc[cond, b] = (beta.ardl.b[(1 + 1 * (case > 1))] / sqrt(Sigmabeta.ardl.b[(1 +
                                                                                            1 * (case > 1)), (1 + 1 * (case > 1))]))
          }
          if (H0 == "find") {
            wtestfi = aod::wald.test(Sigmabeta.ardl.b, beta.ardl.b, Terms = ((2:d) + ((case > 1) * 1)))
            F.indc[cond, b] = wtestfi$result$chi2[1] / wtestfi$result$chi2[2]
          }

          #PRINT MESSAGE
          if (print & B>=200) {
            signal = floor(B / 200)
            if (!(b %% signal)) {
              prog = floor(b / B * 50)
              cat(
                '\r',
                conditional,
                " ",
                H0,
                ": [",
                strrep("#", prog),
                strrep(" ", 50 - prog),
                "] ",
                2 * prog,
                "%",
                sep = ""
              )
              flush.console()
            }
            ## Pretend that some work is being done.
            Sys.sleep(0.001)
            if (b == B) {
              cat(" \n\r")
            }
          }
        }
      }

      if (sum(F.overallc[cond, ]) != 0) {
        quantFOV[, cond] = quantile(F.overallc[cond, ], probs = 1 - crit.H0)
      }
      if (sum(t.depc[, cond]) != 0) {
        quantt[, cond] = quantile(t.depc[cond, ], probs = crit.H0)
      }
      if (sum(F.indc[cond, ]) != 0) {
        quantFIND[, cond] = quantile(F.indc[cond, ], probs = 1 - crit.H0)
      }

    }

    results = list(
      ARDL = mod.ardl.l,
      pssbounds = pss.b,
      smgbounds = smk_test,
      fov.st = unlist(F.overall.os),
      t.st = unlist(t.dep.os),
      find.st = unlist(F.ind.os),
      quantfov = (cbind(crit.H0, quantFOV)),
      quantt = data.frame(cbind(crit.H0, quantt)),
      quantfind = data.frame(cbind(crit.H0, quantFIND))
    )

    class(results) = "bootCT"

    return(results)

  }
