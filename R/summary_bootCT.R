#' Summary method
#'
#' This function summarizes the ARDL bootstrap test and all the other asymptotic procedures all together.
#'
#' @param object an object of class "\code{bootCT}"
#' @param ... not parsed, added for compatibility
#' @return the function returns a list of summary statistics, already present in the function \code{boot_ardl}, and displays them in an appropriate manner.
#' @export
summary.bootCT = function(object,...){
SUMMARY_ARDL =lapply(object$ARDL,function(x) summary(x))
SUMMARY_PSS = object$pssbounds
SUMMARY_SMG = object$smgbounds
FOV = object$quantfov
FIND = object$quantfind
TDEP = object$quantt
FOVSTAT = object$fov.st
TSTAT = object$t.st
FINDSTAT = object$find.st

cat("CONDITIONAL MODEL\n")
print(SUMMARY_ARDL[[1]])
cat("\n\n\nUNCONDITIONAL MODEL\n")
print(SUMMARY_ARDL[[2]])

cat(
  "\n\n\nObservations:",SUMMARY_PSS[[1]]$obs,
  "\nNumber of Lagged Regressors (not including LDV) (k): ",SUMMARY_PSS[[1]]$k,
  "\nCase: ",SUMMARY_PSS[[1]]$case,
  "\n--------------------------------------------------------------------------
-                         PSS    Fov-test                                -
--------------------------------------------------------------------------
                 <------- I(0) ----- I(1) ----->
10% critical value       ",sprintf("%.3f",SUMMARY_PSS[[1]]$ftest.I0.p10),"     ",sprintf("%.3f",SUMMARY_PSS[[1]]$ftest.I1.p10),
  "\n5% critical value        ",sprintf("%.3f",SUMMARY_PSS[[1]]$ftest.I0.p05),"     ",sprintf("%.3f",SUMMARY_PSS[[1]]$ftest.I1.p05),
  "\n1% critical value        ",sprintf("%.3f",SUMMARY_PSS[[1]]$ftest.I0.p01),"     ",sprintf("%.3f",SUMMARY_PSS[[1]]$ftest.I1.p01),
  "\n\nF-statistic (C) = ",sprintf("%.3f",SUMMARY_PSS[[1]]$fstat)," F-statistic (UC) = ",sprintf("%.3f",SUMMARY_PSS[[2]]$fstat),
  "\n\n     Bootstrap critical values\n    ",
  paste0(paste(FOV[,1]*100,"%"),collapse="      "),
  "\nC  ",paste0(sprintf("%.3f",round(FOV[,2],2)),collapse="     "),
  "\nUC ",paste0(sprintf("%.3f",round(FOV[,3],2)),collapse="     "))

if(SUMMARY_PSS[[1]]$case%in%c(1,3,5)){
  cat(
    "\n\n\nObservations:",SUMMARY_PSS[[1]]$obs,
    "\nNumber of Lagged Regressors (not including LDV) (k): ",SUMMARY_PSS[[1]]$k,
    "\nCase: ",SUMMARY_PSS[[1]]$case,
    "\n--------------------------------------------------------------------------
-                          PSS    t-test                                 -
--------------------------------------------------------------------------
                  <------- I(0) ----- I(1) ----->
10% critical value       ",sprintf("%.3f",SUMMARY_PSS[[1]]$ttest.I0.p10),"    ",sprintf("%.3f",SUMMARY_PSS[[1]]$ttest.I1.p10),
    "\n5% critical value        ",sprintf("%.3f",SUMMARY_PSS[[1]]$ttest.I0.p05),"    ",sprintf("%.3f",SUMMARY_PSS[[1]]$ttest.I1.p05),
    "\n1% critical value        ",sprintf("%.3f",SUMMARY_PSS[[1]]$ttest.I0.p01),"    ",sprintf("%.3f",SUMMARY_PSS[[1]]$ttest.I1.p01),
    "\n\nt-statistic (C)= ",sprintf("%.3f",SUMMARY_PSS[[1]]$tstat)," t-statistic (UC)= ",sprintf("%.3f",SUMMARY_PSS[[2]]$tstat),
    "\n\n     Bootstrap critical values\n    ",
    paste0(paste(TDEP[,1]*100,"%"),collapse="       "),
    "\nC  ",paste0(sprintf("%.3f",round(TDEP[,2],2)),collapse="     "),
    "\nUC ",paste0(sprintf("%.3f",round(TDEP[,3],2)),collapse="     "))}else{
      cat(
        "\n\n\nObservations:",SUMMARY_PSS[[1]]$obs,
        "\nNumber of Lagged Regressors (not including LDV) (k): ",SUMMARY_PSS[[1]]$k,
        "\nCase: ",SUMMARY_PSS[[1]]$case,"
       \nPSS bound t-test not applicable for case ",SUMMARY_PSS[[1]]$case,
        "\n\nt-statistic (C)= ",sprintf("%.3f",round(TSTAT[1],3))," t-statistic (UC)= ",sprintf("%.3f",round(TSTAT[2],3)),
        "\n\n     Bootstrap critical values\n    ",
        paste0(paste(TDEP[,1]*100,"%"),collapse="       "),
        "\nC  ",paste0(sprintf("%.3f",round(TDEP[,2],2)),collapse="     "),
        "\nUC ",paste0(sprintf("%.3f",round(TDEP[,3],2)),collapse="     "))
    }

if(SUMMARY_PSS[[1]]$case%in%c(1,3,5)){
  cat(
    "\n\n\nObservations:",SUMMARY_PSS[[1]]$obs,
    "\nNumber of Lagged Regressors (not including LDV) (k): ",SUMMARY_PSS[[1]]$k,
    "\nCase: ",SUMMARY_PSS[[1]]$case,
    "\n-----------------------------------------------------------------------------
-                         SMG    FInd-test                                  -
-----------------------------------------------------------------------------
                 <------- I(0) ----- I(1) ----->
10% critical value       ",SUMMARY_SMG[1,4],"     ",SUMMARY_SMG[1,5],
    "\n5% critical value        ",SUMMARY_SMG[2,4],"     ",SUMMARY_SMG[2,5],
    "\n2.5% critical value      ",SUMMARY_SMG[3,4],"     ",SUMMARY_SMG[3,5],
    "\n1% critical value        ",SUMMARY_SMG[4,4],"     ",SUMMARY_SMG[4,5],
    "\n\nF-statistic (C)= ", sprintf("%.3f",round(FINDSTAT[1],3))," F-statistic (UC)= ",sprintf("%.3f",round(FINDSTAT[2],3)),
    "\n\n     Bootstrap critical values\n    ",
    paste0(paste(FIND[,1]*100,"%"),collapse="      "),
    "\nC  ",paste0(sprintf("%.3f",round(FIND[,2],2)),collapse="     "),
    "\nUC ",paste0(sprintf("%.3f",round(FIND[,3],2)),collapse="     "))}else{
      cat(
        "\n\n\nObservations:",SUMMARY_PSS[[1]]$obs,
        "\nNumber of Lagged Regressors (not including LDV) (k): ",SUMMARY_PSS[[1]]$k,
        "\nCase: ",SUMMARY_PSS[[1]]$case,
        "\nSMG bound F-test on independent variables not applicable for case ",SUMMARY_PSS[[1]]$case,
        "\n\nF-statistic (C)= ", sprintf("%.3f",round(FINDSTAT[1],3))," F-statistic (UC)= ",sprintf("%.3f",round(FINDSTAT[2],3)),
        "\n\n     Bootstrap critical values\n    ",
        paste0(paste(FIND[,1]*100,"%"),collapse="      "),
        "\nC  ",paste0(sprintf("%.3f",round(FIND[,2],2)),collapse="     "),
        "\nUC ",paste0(sprintf("%.3f",round(FIND[,3],2)),collapse="     "))
    }
}
#' @examples
