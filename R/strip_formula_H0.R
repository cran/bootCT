#' Removes the variables imposed under the null from an ARDL formula object
#'
#' @param formula.ardl the ARDL formula object
#' @param H0 the null hypothesis to be tested
#' @param case ARDL model case for the treatment of intercept and trend
#' @param d number of variables
#' @keywords internal
#' @export
#'
strip_formula_H0 = function(formula.ardl,H0,case,d){
  char.eqn = as.character(formula.ardl)
  str.pieces = gsub(" ", "",strsplit(char.eqn[3],split = "[+]")[[1]], fixed = TRUE)
  str.pieces.H0 = str.pieces
  if(H0=="fov"){k=1:d}
  if(H0=="t"){k=1}
  if(H0=="find"){k=2:d}
  str.pieces.H0=str.pieces[-k]
  str.interc=NA
  if(case == 4){
    str.pieces.H0=str.pieces.H0[-length(str.pieces.H0)]
  }
  if(case <= 2){
    str.interc="-1"
  }
  ans0 = na.omit(c(char.eqn[2],char.eqn[1],paste0(str.pieces.H0,collapse="+"),str.interc))
  ans = as.formula(paste(ans0,collapse = ""))
  return(ans)
}
