#' Bootstrap ARDL
#'
#' This is the main function of the package. It performs the bootstrap version of the ARDL bound test for cointegration.
#'
#' @param data Input dataset. Must contain a dependent and a set of independent variables.
#' @param yvar Name of the dependent variable, enclosed in quotation marks. If NULL, the first variable will be used.
#' @param xvar Vector of names of the independent variables, each enclosed in quotation marks. If NULL, all variables except the first will be used.
#' @param fix.ardl Fixed lagged differences for the short term part of the ARDL equation.
#' @param info.ardl Selection criterion for the auto_ardl function. Options are "AIC", "AICc", BIC, "R2", "adjR2", if fix.ardl is null. Defaults to AIC.
#' @param fix.vecm Fixed lagged differences for the short term part of the VECM equation.
#' @param info.vecm Selection criterion for the VARselect function. Options are "AIC", "HQ", "SC", "FPE", if fix.vecm is null. Defaults to AIC.
#' @param maxlag Max number of lags for the auto_ardl  and VARselect procedures, if fix.ardl or fix.vecm are null, respectively.
#' @param a.ardl Threshold significance for the short-term ARDL coefficients significance.
#' @param a.vecm Threshold significance for the short-term VECM coefficients significance.
#' @param nboot Number of bootstrap replications.
#' @param case Model case, pertaining to the treatment of intercept and trend. Must be integer from 1 to 5. Defaults to 3.
#' @param a.boot.H0 Probability/ies by which the critical quantiles of the bootstrap distribution(s) must be calculated.
#' @param print Show the progress bar.
#' @return List of several elements including \itemize{
#'\item \code{data}: the data used to perform estimation and testing
#'\item \code{ARDL}: the estimated ARDL conditional model
#'\item \code{VECM}: the estimated VECM unconditional model
#'\item \code{jo.testX}: Johansen cointegration test on the independent variables
#'\item \code{pssbounds}: the PSS bound test output
#'\item \code{smgbounds}: the SMG bound test critical values
#'\item \code{fov.stat}: the test statistics on the conditional Fov tests
#'\item \code{t.stat}: the test statistics on the conditional t test
#'\item \code{find.stat}: the test statistics on the conditional Find tests
#'\item \code{quantfov}: the bootstrap conditional F Overall test critical value(s)
#'\item \code{quantt}: the bootstrap conditional t-test critical value(s)
#'\item \code{quantfind}: the bootstrap conditional F Independent test critical value(s)
#'\item \code{fakecoint}: indication of the situation in which \eqn{a_{y.x}\neq 0} but \eqn{a_{y.x}^{UC}=0}, signaling absence of cointegration.}
#'
#' @examples
#'\dontrun{
#' data(ger_macro)
#' LNDATA = as.data.frame(log(ger_macro[,-1]))
#' colnames(LNDATA) = c("LNINVEST","LNINCOME","LNCONS")
#'
#' boot_res = boot_ardl(data=LNDATA,
#'                      yvar = "LNINCOME", xvar = c("LNCONS","LNINVEST"),
#'                      maxlag = 5, nboot = 2000)
#' summary(boot_res)
#'}
#' @export
boot_ardl =
  function(data,
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
    if (is.null(xvar)) {
      xvar = dfnames[-1]
    }

    if (!(yvar %in% dfnames)) {
      stop("y variable not in dataset.")
    }

    if (!all(xvar %in% dfnames)) {
      stop("at least one of the X variables not in dataset.")
    }

    if (yvar %in% xvar) {
      stop("y variable and X variables must be separated.")
    }

    if (a.ardl < 0 | a.ardl > 1) {
      warning("Invalid ARDL threshold p.value. Defaulting to 0.05.")
      a.ardl = 0.05
    }

    if (a.vecm < 0 | a.vecm > 1) {
      warning("Invalid VECM threshold p.value. Defaulting to 0.05.")
      a.vecm = 0.05
    }

    if (length(a.ardl) > 1) {
      warning("Expecting a scalar a.ardl value. Defaulting to 0.05.")
      a.ardl = 0.05
    }

    if(length(a.vecm) > 1) {
      warning("Expecting a scalar a.vecm value. Defaulting to 0.05.")
      a.vecm = 0.05
    }

    if(!is.null(fix.ardl)) {
      if(any(fix.ardl %% 1 != 0)){
      warning("Expecting integer fix.ardl values. Rounding to the nearest integer.")
      fix.ardl = round(fix.ardl)
      }
    }

    if (!is.null(fix.vecm)) {
      if(fix.vecm %% 1 != 0){
      warning("Expecting integer fix.vecm value. Rounding to the nearest integer.")
      fix.vecm = round(fix.vecm)
      }
    }

    if (!is.null(maxlag)) {
    if (any(maxlag %% 1 != 0)) {
      warning("Expecting integer maxlag value. Rounding to the nearest integer.")
      maxlag = round(maxlag)
    }
    }

    if (nboot %% 1 != 0) {
      warning("Expecting integer nboot value. Rounding to the nearest integer.")
      nboot = round(nboot)
    }

    if (case < 1 | case > 5 | case %% 1 != 0) {
      warning("Inadmissible case. Defaulting to case III.")
      case = 3
    }

    if (any(a.boot.H0 < 0) | any(a.boot.H0 > 1)) {
      warning("Invalid probability for critical values. Defaulting to 1%, 2.5% and 5%")
      a.boot.H0 = c(0.01, 0.025, 0.05)
    }

    if(nboot<200){
      warning("Suppressing print: nboot should be at least 200")
    }

    #SETTING FLAGS FOR SPURIOUS COINTEGRATION
    flagspur = rep(0,length(a.boot.H0))
    computeuc = 0
    
    numdata = data %>% dplyr::select(all_of(yvar), all_of(xvar))
    d = ncol(numdata)

    #LAGGED AND FIRST DIFFERENCE
    lagdata = lag_mts(as.matrix(numdata), k = rep(1, d))
    dlag0 = numdata - lagdata
    colnames(dlag0) <- paste("D_", colnames(numdata), sep = "")

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
    }

    if (case == 2) {
      formula.ardl = as.formula(paste(substr(text.ardl, 1, nchar(text.ardl) -
                                               7)))
    }

    if (case == 3) {
      formula.ardl = as.formula(paste(substr(text.ardl, 1, nchar(text.ardl) -
                                               7)))
    }

    if (case >= 4) {
      formula.ardl = as.formula(text.ardl)
    }

    if(is.null(fix.ardl)){
      if (case == 1) {
        fixed = c(-1, rep(0, d), rep(-1, (d - 1)))
      }

      if (case == 2) {
        fixed = c(-1, rep(0, d), rep(-1, (d - 1)))
      }

      if (case == 3) {
        fixed = c(-1, rep(0, d), rep(-1, (d - 1)))
      }

      if (case >= 4) {
        fixed = c(-1, rep(0, d), rep(-1, (d - 1)), 0)
      }

      if(is.null(info.ardl) | !(info.ardl %in% c("AIC","AICc","BIC","R2","adjR2"))){
        info.ardl="AIC"
        warning("Default ARDL lag choice criterion: AIC")
      }

      #SELECT BEST ARDL ORDER
      ardl_select = ARDL::auto_ardl(
        formula.ardl,
        data = na.omit(df.ardl),
        max_order = maxlag,
        fixed_order = fixed,
        selection=info.ardl
      )
      best_ardl = ardl_select$best_order[c(1, (d + 2):(2 * d))]

      #CONSTRUCT FINAL ARDL DATAFRAME
      final_dlag = lag_mts(dlag0, k = best_ardl)
      formula.ardl2=ardl_select$best_model$full_formula
    } else{
      if (length(fix.ardl) != d) {
        warning("number of variable is not equal to number of fix.ardl, taking only the first value")
        difflag0 = rep(fix.ardl[1], d)
        fix.ardl = difflag0
      }

      final_dlag = lag_mts(dlag0, k = fix.ardl)
      best_ardl = fix.ardl
      text.ardl2 = paste0(colnames(final_dlag),collapse = '+')
      formula.ardl2=paste(as.character(formula.ardl)[2],
                          as.character(formula.ardl)[1],
                          as.character(formula.ardl)[3],"+",
                          text.ardl2)
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
    F.overallc = t.depc = F.indc = matrix(0, nrow = 2, ncol = nboot)
    quantFOV = quantFIND = quantt = matrix(0, nrow = length(a.boot.H0), ncol = 2)


    typevecm = ifelse(case < 2, "none", ifelse(case < 4, "const", "both"))

    #SELECT UNRESTRICTED VECM
    if(is.null(fix.vecm)){

      if(is.null(info.vecm) | !(info.vecm %in% c("AIC","HQ","SC","FPE"))){
        info.vecm="AIC"
        idv = 1
        warning("Default VECM lag choice criterion: AIC")
      }else{
        idv = 1 * (info.vecm == "AIC") +
          2 * (info.vecm == "HQ") +
          3 * (info.vecm == "SC") +
          4 * (info.vecm == "FPE")
      }

      vecm = vars::VARselect(
        na.omit(dlag0),
        lag.max = maxlag,
        type = typevecm,
        exogen = na.omit(lagdata))

      vecmsel = vecm$selection[idv]
    }else{
      vecmsel=fix.vecm
    }

    #CONFIRM UNRESTRICTED VECM
    vecm = vars::VAR(
      na.omit(dlag0),
      p = vecmsel,
      type = typevecm,
      exogen = na.omit(lagdata)
    )

    summary(vecm)

    nx = nrow(na.omit(dlag0))
    dx = d * vecmsel + d * (d - 1)

    #IMPOSE WEAK EXOGENEITY AND REMOVE
    #NON-SIGNIFICANT SHORT-RUN ESTIMATES
    vecm.ser = vars::restrict(vecm,
                              method = "ser",
                              thresh = qnorm(1 - a.vecm))
    mat.ser = vecm.ser$restrictions
    mat.ser[,(ncol(mat.ser) - d + 1)] = c(1, rep(0, d - 1))
    if(typevecm=="const") mat.ser[, (ncol(mat.ser) - d)] = 1
    if(typevecm=="trend") mat.ser[, (ncol(mat.ser) - d)] = 1
    if(typevecm=="both"){
      mat.ser[, (ncol(mat.ser) - d)] = 1
      mat.ser[, (ncol(mat.ser) - d - 1)] = 1
    }
    mat.ser[, (ncol(mat.ser) - d + 2):ncol(mat.ser)] = 1
    vecm.ser = vars::restrict(vecm, method = "manual", resmat = mat.ser)

    vecm.out=vecm.ser
    
    H0=list(H0.c = c("fov","t","find"), H0.uc = "find")

    for (cond in 1:2) {
      formula.ardlx = as.formula(formula.ardl2)
      conditional = ifelse(cond == 1, "Conditional", "Unconditional")
      colnames(df.ardl.l[[cond]])[1] = colnames(dlag0)[1]

      #FINAL ARDL ESTIMATE
      if (case >= 4) {
        df.ardl.l[[cond]]$timeid = df.ardl$timeid
      }

      formula.ardlx = ardl_to_lm(formula.ardlx, case = case, d = d)

      if (cond == 2) {
        formula.ardlx = strip_formula_UC(formula.ardlx,case,d)
      }

      #ESTIMATE UNRESTRICTED ARDL MODEL
      mod.ardl.l[[cond]] = lm(formula.ardlx, data = na.omit(df.ardl.l[[cond]]))
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

      #smk_crit = bootCT::smk_crit
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

      h=H0[[cond]]

      for (hid in 1:length(h)) {
        #F OVERALL IMPOSING THE NULL
        if (h[hid] == "fov") {
          formula.ardl.H0 = strip_formula_H0(formula.ardlx, h[hid], case, d)
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
              (summary(mod.ardl.H0)[[4]][, 4] < a.ardl)[1:length(beta.ardl.H0)]
            join_nameG = which(names(GAMMAX[1,]) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            join_nameB = which(names(hbeta.ardl.H0) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            GAMMAX[1,-join_nameG] = 0
            GAMMAX[1, join_nameG] = hbeta.ardl.H0[join_nameB]
          }

          if (case == 2) {
            azero[1] = 0
            #INTERCEPT CONSTRAINED TO VECM
            hbeta.ardl.H0 = beta.ardl.H0[1:length(beta.ardl.H0)] *
              (summary(mod.ardl.H0)[[4]][, 4] < a.ardl)[1:length(beta.ardl.H0)]
            join_nameG = which(names(GAMMAX[1,]) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            join_nameB = which(names(hbeta.ardl.H0) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            GAMMAX[1,-join_nameG] = 0
            GAMMAX[1, join_nameG] = hbeta.ardl.H0[join_nameB]
          }

          if (case == 3) {
            azero[1] = beta.ardl.H0[1] #ARDL UNCONSTRAINED INTERCEPT
            hbeta.ardl.H0 = beta.ardl.H0[2:length(beta.ardl.H0)] *
              (summary(mod.ardl.H0)[[4]][, 4] < a.ardl)[2:length(beta.ardl.H0)]
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
              (summary(mod.ardl.H0)[[4]][, 4] < a.ardl)[2:length(beta.ardl.H0)]
            join_nameG = which(names(GAMMAX[1,]) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            join_nameB = which(names(hbeta.ardl.H0) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            GAMMAX[1,-join_nameG] = 0
            GAMMAX[1, join_nameG] = hbeta.ardl.H0[join_nameB]
          }

          if (case == 5) {
            azero[1] = beta.ardl.H0[1] #ARDL UNCONSTRAINED INTERCEPT
            auno[1] = beta.ardl.H0[length(beta.ardl.H0)] #ARDL UNCONSTRAINED TREND
            hbeta.ardl.H0 = beta.ardl.H0[2:(length(beta.ardl.H0) - 1)] *
              (summary(mod.ardl.H0)[[4]][, 4] < a.ardl)[2:(length(beta.ardl.H0) -
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

        if (h[hid] == "t") {
          formula.ardl.H0 = strip_formula_H0(formula.ardlx, h[hid], case,d)

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
              (summary(mod.ardl.H0)[[4]][, 4] < a.ardl)[d:length(beta.ardl.H0)]
            join_nameG = which(names(GAMMAX[1,]) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            join_nameB = which(names(hbeta.ardl.H0) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            GAMMAX[1,-join_nameG] = 0
            GAMMAX[1, join_nameG] = hbeta.ardl.H0[join_nameB]
          }
          if (case == 2) {
            azero[1] = 0
            ay.x = -beta.ardl.H0[1:(d - 1)]
            hbeta.ardl.H0 = beta.ardl.H0[d:length(beta.ardl.H0)] *
              (summary(mod.ardl.H0)[[4]][, 4] < a.ardl)[d:length(beta.ardl.H0)]
            join_nameG = which(names(GAMMAX[1,]) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            join_nameB = which(names(hbeta.ardl.H0) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            GAMMAX[1,-join_nameG] = 0
            GAMMAX[1, join_nameG] = hbeta.ardl.H0[join_nameB]
          }
          if (case == 3) {
            azero[1] = beta.ardl.H0[1]
            ay.x = -beta.ardl.H0[2:d]
            hbeta.ardl.H0 = beta.ardl.H0[(d+1):length(beta.ardl.H0)] *
              (summary(mod.ardl.H0)[[4]][, 4] < a.ardl)[(d+1):length(beta.ardl.H0)]
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
              (summary(mod.ardl.H0)[[4]][, 4] < a.ardl)[(d+1):length(beta.ardl.H0)]
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
              (summary(mod.ardl.H0)[[4]][, 4] < a.ardl)[(d+1):(length(beta.ardl.H0) - 1)]
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

        if (h[hid] == "find") {

          formula.ardl.H0 = strip_formula_H0(formula.ardlx, h[hid], case,d)

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
              (summary(mod.ardl.H0)[[4]][, 4] < a.ardl)[2:length(beta.ardl.H0)]
            join_nameG = which(names(GAMMAX[1,]) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            join_nameB = which(names(hbeta.ardl.H0) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            GAMMAX[1,-join_nameG] = 0
            GAMMAX[1, join_nameG] = hbeta.ardl.H0[join_nameB]
          }
          if (case == 2) {
            azero[1] = 0
            ayy = -beta.ardl.H0[1]
            hbeta.ardl.H0 = beta.ardl.H0[2:length(beta.ardl.H0)] *
              (summary(mod.ardl.H0)[[4]][, 4] < a.ardl)[2:length(beta.ardl.H0)]
            join_nameG = which(names(GAMMAX[1,]) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            join_nameB = which(names(hbeta.ardl.H0) %in% intersect(names(GAMMAX[1,]), names(hbeta.ardl.H0)))
            GAMMAX[1,-join_nameG] = 0
            GAMMAX[1, join_nameG] = hbeta.ardl.H0[join_nameB]
          }
          if (case == 3) {
            azero[1] = beta.ardl.H0[1]
            ayy = -beta.ardl.H0[2]
            hbeta.ardl.H0 = beta.ardl.H0[3:length(beta.ardl.H0)] *
              (summary(mod.ardl.H0)[[4]][, 4] < a.ardl)[3:length(beta.ardl.H0)]
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
              (summary(mod.ardl.H0)[[4]][, 4] < a.ardl)[3:length(beta.ardl.H0)]
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
              (summary(mod.ardl.H0)[[4]][, 4] < a.ardl)[3:(length(beta.ardl.H0) -
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
        for (b in 1:nboot) {
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
          if (h[hid] == "fov") {
            wtestfo = aod::wald.test(Sigmabeta.ardl.b, beta.ardl.b, Terms = terms.ur)
            F.overallc[cond, b] = wtestfo$result$chi2[1] / wtestfo$result$chi2[2]
          }
          if (h[hid] == "t") {
            t.depc[cond, b] = (beta.ardl.b[(1 + 1 * (case > 1))] / sqrt(Sigmabeta.ardl.b[(1 +
                                                                                            1 * (case > 1)), (1 + 1 * (case > 1))]))
          }
          if (h[hid] == "find") {
            wtestfi = aod::wald.test(Sigmabeta.ardl.b, beta.ardl.b, Terms = ((2:d) + ((case > 1) * 1)))
            F.indc[cond, b] = wtestfi$result$chi2[1] / wtestfi$result$chi2[2]
          }

          #PRINT MESSAGE
          if (print & nboot>=200) {
            signal = floor(nboot / 200)
            if (!(b %% signal)) {
              prog = floor(b / nboot * 50)
              cat(
                '\r',
                conditional,
                " ",
                h[hid],
                ": [",
                strrep("#", prog),
                strrep(" ", 50 - prog),
                "] ",
                2 * prog,
                "%",
                sep = ""
              )
              utils::flush.console()
            }
            Sys.sleep(0.001)
            if (b == nboot) {
              cat(" \n\r")
            }
          }
        }
      }

      if (sum(F.overallc[cond, ]) != 0) {
        quantFOV[, cond] = quantile(F.overallc[cond, ], probs = 1 - a.boot.H0)
      }
      if (sum(t.depc[, cond]) != 0) {
        quantt[, cond] = quantile(t.depc[cond, ], probs = a.boot.H0)
      }
      if (sum(F.indc[cond, ]) != 0) {
        quantFIND[, cond] = quantile(F.indc[cond, ], probs = 1 - a.boot.H0)
      }
      if(cond==1){
        findflag.c =rep(unlist(F.ind.os)[1],length(a.boot.H0))>quantFIND[,1]
        fakecoint = cbind(a.boot.H0,0+(findflag.c))

        if(all(fakecoint[,2]==0)){
          break
        }else{
          print("Computing Find test on the UC model to detect spurious cointegration")
          computeuc = 1
        }
      }

    }
    
    #ENRICH OUTPUT IF UNCONDITIONAL MODEL IS CONSIDERED
    if(computeuc == 1){
      findflag.uc = rep(unlist(F.ind.os)[2],length(a.boot.H0))<quantFIND[,2]
      flagspur = (findflag.uc & findflag.c) + 0
      find.stat = unlist(F.ind.os)
      quant.find = data.frame(cbind(a.boot.H0, quantFIND[,1], quantFIND[,2]))
    }else{
      find.stat = F.ind.os[[1]]
      quant.find = data.frame(cbind(a.boot.H0, quantFIND[,1]))
    }

    data.out=data.frame(numdata)

    #extra: Johansen tests on X 
    
    if(case==1) ecdet="none"
    if(case==2) ecdet="const"
    if(case==3) ecdet="none"
    if(case==4) ecdet="trend"
    if(case==5) ecdet="none"

    VECM.ctrace = urca::ca.jo(numdata[,-1],
                            ecdet=ecdet,
                            K=vecmsel+1,
                            spec="transitory",
                            type="trace")

    VECM.ceigen = urca::ca.jo(numdata[,-1],
                            ecdet=ecdet,
                            K=vecmsel+1,
                            spec="transitory",
                            type="eigen")

    results = list(
      data = data.out,
      ARDL = mod.ardl.l[[1]],
      VECM = vecm.out,
      jo.test = list(trace=VECM.ctrace,eigen=VECM.ceigen),
      pssbounds = pss.b,
      smgbounds = smk_test,
      fov.stat = F.overall.os[[1]],
      t.stat = t.dep.os[[1]],
      find.stat = find.stat,
      quantfov = data.frame(cbind(a.boot.H0, quantFOV[,1])),
      quantt = data.frame(cbind(a.boot.H0, quantt[,1])),
      quantfind = quant.find,
      fakecoint = flagspur
    )

    class(results) = "bootCT"

    return(results)

  }
