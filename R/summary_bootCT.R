  #' Summary method
#'
#' This function summarizes the ARDL bootstrap test and all the other asymptotic procedures all together.
#'
#' @param object an object of class "\code{bootCT}"
#' @param ... additional arguments, e.g. \code{out}: subset of output to print. Options (can be multiple) are: "all", "ARDL", "VECM", "cointVECM", "cointARDL". Defaults to "all".
#' @return the function returns a list of summary statistics, already present in the function \code{boot_ardl}, and displays them in an appropriate manner. Depending on the "out" argument, ARDL/VECM estimation outputs and/or ARDL/VECM cointegration tests can be displayed.
#' @export
summary.bootCT = function(object,
                          ...) {
  arguments = list(...)
  names = names(arguments)
  
  if(length(arguments)==0){
    out="all"
  }
  if(length(arguments)>1){
      nargs=unlist(names(arguments))
      if(any(nargs=="out")){
        out=arguments[[which(nargs=="out")]]
        warning("Extra arguments ignored, only \"out\" appropriate with bootCT class")
      }else{
        warning("Extra arguments ignored, only \"out\" appropriate with bootCT class")
      }
  }
  if(length(arguments)==1){
    nargs=unlist(names(arguments))
    if(any(nargs=="out")){
      out=arguments[[1]]
    }else{
      warning("Extra arguments ignored, only \"out\" appropriate with bootCT class")
    }
  }
  
  if(!all(out%in%c("all", "ARDL", "VECM", "cointVECM", "cointARDL"))){
    warning("Invalid 'out' options: use 'all','ARDL', 'VECM', 'cointVECM','cointARDL'")
    out=out[which(out%in%c("all", "ARDL", "VECM", "cointVECM", "cointARDL"))]
    if(length(out)==0){
      out="all"
    }
  }
  
  if(any(out=="all")){
    out="all"
  }
  
  SUMMARY_ARDL = summary(object$ARDL)
  SUMMARY_VECM = summary(object$VECM)
  SUMMARY_PSS = object$pssbounds
  SUMMARY_SMG = object$smgbounds
  FOV = object$quantfov
  FIND = object$quantfind
  TDEP = object$quantt
  TSTAT = object$t.st
  FINDSTAT = object$find.st
  FAKECOINT = object$fakecoint
  
  addtext = ""
  
  if (any(FAKECOINT == 1)) {
    selfake = paste0(FOV[which(FAKECOINT == 1),1], sep = ",")
    addtext = cat(
      "\n Conditional model rejects Fov null.
  \n However, Unconditional model does not reject Fov null for significance levels
  ",
      selfake,
      "\n system NOT COINTEGRATED."
    )
  }
  
  if (out == "all" | any(out == "ARDL")) {
    cat("CONDITIONAL ARDL MODEL\n")
    print(SUMMARY_ARDL)
  }
  
  if (out == "all" | any(out == "VECM")) {
    cat("UNCONDITIONAL VECM MODEL\n")
    print(SUMMARY_VECM)
  }
  
  if (out == "all" | any(out == "cointVECM")) {
    UCLtrace=unclass(object$jo.test$trace)
    UCLeigen=unclass(object$jo.test$eigen)
    
    sval1 = as.numeric(attr(UCLtrace,"teststat"))
    testats1=cbind(round(sval1,2),attr(UCLtrace,"cval"))
    dimnames(testats1)[[2]][1]="test"
    cat("VECM MODEL JOHANSEN COINTEGRATION TESTS\n\n")
    cat("Test type:",attr(UCLtrace,"type"),",",attr(UCLtrace,"model"),"\n\n")
    cat("Eigenvalues (lambda):\n")
    print(attr(UCLtrace,"lambda"))
    cat("\n\n")
    cat("Values of test statistic and critical values of test:\n\n")
    print(testats1)
    cat("\n\n")
    cat("Eigenvectors, normalised to first column:\n")
    cat("(These are the cointegration relations)\n\n")
    print(attr(UCLtrace,"V"))
    cat("\n\n")
    cat("Weights W:\n")
    cat("(This is the loading matrix)\n\n")
    print(attr(UCLtrace,"W"))
    
    cat("\n\n=======================================\n\n")
    
    sval2 = as.numeric(attr(UCLeigen,"teststat"))
    testats2=cbind(round(sval2,2),attr(UCLeigen,"cval"))
    dimnames(testats2)[[2]][1]="test"
    cat("Test type:",attr(UCLeigen,"type"),",",attr(UCLeigen,"model"),"\n\n")
    cat("Eigenvalues (lambda):\n")
    print(attr(UCLeigen,"lambda"))
    cat("\n\n")
    cat("Values of test statistic and critical values of test:\n\n")
    print(testats2)
    cat("\n\n")
    cat("Eigenvectors, normalised to first column:\n")
    cat("(These are the cointegration relations)\n\n")
    print(attr(UCLeigen,"V"))
    cat("\n\n")
    cat("Weights W:\n")
    cat("(This is the loading matrix)\n\n")
    print(attr(UCLeigen,"W"))
  }
  
  if (out == "all" | any(out == "cointARDL")) {
    cat(
      "\n\n\nObservations:",
      SUMMARY_PSS[[1]]$obs,
      "\nNumber of Lagged Regressors (not including LDV) (k): ",
      SUMMARY_PSS[[1]]$k,
      "\nCase: ",
      SUMMARY_PSS[[1]]$case,
      "\n--------------------------------------------------------------------------
-                         PSS    Fov-test                                -
--------------------------------------------------------------------------
                 <------- I(0) ----- I(1) ----->
10% critical value       ",
      sprintf("%.3f", SUMMARY_PSS[[1]]$ftest.I0.p10),
      "     ",
      sprintf("%.3f", SUMMARY_PSS[[1]]$ftest.I1.p10),
      "\n5% critical value        ",
      sprintf("%.3f", SUMMARY_PSS[[1]]$ftest.I0.p05),
      "     ",
      sprintf("%.3f", SUMMARY_PSS[[1]]$ftest.I1.p05),
      "\n1% critical value        ",
      sprintf("%.3f", SUMMARY_PSS[[1]]$ftest.I0.p01),
      "     ",
      sprintf("%.3f", SUMMARY_PSS[[1]]$ftest.I1.p01),
      "\n\nF-statistic = ",
      sprintf("%.3f", SUMMARY_PSS[[1]]$fstat),
      "\n\n     Bootstrap critical values\n    ",
      paste0(paste(FOV[, 1] * 100, "%"), collapse = "      "),
      "\n  ",
      paste0(sprintf("%.3f", round(FOV[, 2], 2)), collapse = "     "),
      addtext
    )
    
    if (SUMMARY_PSS[[1]]$case %in% c(1, 3, 5)) {
      cat(
        "\n\n\nObservations:",
        SUMMARY_PSS[[1]]$obs,
        "\nNumber of Lagged Regressors (not including LDV) (k): ",
        SUMMARY_PSS[[1]]$k,
        "\nCase: ",
        SUMMARY_PSS[[1]]$case,
        "\n--------------------------------------------------------------------------
-                          PSS    t-test                                 -
--------------------------------------------------------------------------
                  <------- I(0) ----- I(1) ----->
10% critical value       ",
        sprintf("%.3f", SUMMARY_PSS[[1]]$ttest.I0.p10),
        "    ",
        sprintf("%.3f", SUMMARY_PSS[[1]]$ttest.I1.p10),
        "\n5% critical value        ",
        sprintf("%.3f", SUMMARY_PSS[[1]]$ttest.I0.p05),
        "    ",
        sprintf("%.3f", SUMMARY_PSS[[1]]$ttest.I1.p05),
        "\n1% critical value        ",
        sprintf("%.3f", SUMMARY_PSS[[1]]$ttest.I0.p01),
        "    ",
        sprintf("%.3f", SUMMARY_PSS[[1]]$ttest.I1.p01),
        "\n\nt-statistic = ",
        sprintf("%.3f", SUMMARY_PSS[[1]]$tstat),
        "\n\n     Bootstrap critical values\n    ",
        paste0(paste(TDEP[, 1] * 100, "%"), collapse = "       "),
        "\n  ",
        paste0(sprintf("%.3f", round(TDEP[, 2], 2)), collapse = "     ")
      )
    } else{
      cat(
        "\n\n\nObservations:",
        SUMMARY_PSS[[1]]$obs,
        "\nNumber of Lagged Regressors (not including LDV) (k): ",
        SUMMARY_PSS[[1]]$k,
        "\nCase: ",
        SUMMARY_PSS[[1]]$case,
        "
       \nPSS bound t-test not applicable for case ",
        SUMMARY_PSS[[1]]$case,
        "\n\nt-statistic = ",
        sprintf("%.3f", round(TSTAT, 3)),
        "\n\n     Bootstrap critical values\n    ",
        paste0(paste(TDEP[, 1] * 100, "%"), collapse = "       "),
        "\n  ",
        paste0(sprintf("%.3f", round(TDEP[, 2], 2)), collapse = "     ")
      )
    }
    
    if (SUMMARY_PSS[[1]]$case %in% c(1, 3, 5)) {
      cat(
        "\n\n\nObservations:",
        SUMMARY_PSS[[1]]$obs,
        "\nNumber of Lagged Regressors (not including LDV) (k): ",
        SUMMARY_PSS[[1]]$k,
        "\nCase: ",
        SUMMARY_PSS[[1]]$case,
        "\n-----------------------------------------------------------------------------
-                         SMG    FInd-test                                  -
-----------------------------------------------------------------------------
                 <------- I(0) ----- I(1) ----->
10% critical value       ",
        SUMMARY_SMG[1, 4],
        "     ",
        SUMMARY_SMG[1, 5],
        "\n5% critical value        ",
        SUMMARY_SMG[2, 4],
        "     ",
        SUMMARY_SMG[2, 5],
        "\n2.5% critical value      ",
        SUMMARY_SMG[3, 4],
        "     ",
        SUMMARY_SMG[3, 5],
        "\n1% critical value        ",
        SUMMARY_SMG[4, 4],
        "     ",
        SUMMARY_SMG[4, 5],
        "\n\nF-statistic = ",
        sprintf("%.3f", round(FINDSTAT[1], 3)),
        "\n\n     Bootstrap critical values\n    ",
        paste0(paste(FIND[, 1] * 100, "%"), collapse = "      "),
        "\n  ",
        paste0(sprintf("%.3f", round(FIND[, 2], 2)), collapse = "     ")
      )
    } else{
      cat(
        "\n\n\nObservations:",
        SUMMARY_PSS[[1]]$obs,
        "\nNumber of Lagged Regressors (not including LDV) (k): ",
        SUMMARY_PSS[[1]]$k,
        "\nCase: ",
        SUMMARY_PSS[[1]]$case,
        "\nSMG bound F-test on independent variables not applicable for case ",
        SUMMARY_PSS[[1]]$case,
        "\n\nF-statistic = ",
        sprintf("%.3f", round(FINDSTAT, 3)),
        "\n\n     Bootstrap critical values\n    ",
        paste0(paste(FIND[, 1] * 100, "%"), collapse = "      "),
        "\n  ",
        paste0(sprintf("%.3f", round(FIND[, 2], 2)), collapse = "     ")
      )
    }
  }
}

  #' @examples
