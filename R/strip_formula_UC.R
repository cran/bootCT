#' Removes the unlagged differences imposed under the conditional ARDL model
#'
#' @param formula.ardl the ARDL formula object
#' @param case ARDL model case for the treatment of intercept and trend
#' @param d number of variables
#' @keywords internal
#' @export
#' 

strip_formula_UC = function(formula.ardl,case,d){
  char.eqn=as.character(formula.ardl)
  str.pieces=strsplit(char.eqn[3],split="[+]")[[1]]
  str.pieces.uc=str.pieces
  
  subvec_D=(str.pieces[which(substr(str.pieces,1,2)=="D_" | substr(str.pieces,2,3)=="D_")])
  subvec_L=(str.pieces[which(substr(str.pieces,1,2)!="D_" & substr(str.pieces,2,3)!="D_")])[1:d]
  subvec_DL=(subvec_D[which(substr(subvec_D,nchar(subvec_D)-3,nchar(subvec_D)-2)==".l")])
  subvec_D0 = (subvec_D[which(substr(subvec_D,nchar(subvec_D)-3,nchar(subvec_D)-2)!=".l")])
  timechar = ifelse(case>=4,"timeid",NA)
  ordered_formula = na.omit(c(subvec_L,subvec_DL,timechar))
  
  ans = as.formula(paste(char.eqn[2],char.eqn[1],paste0(ordered_formula,collapse="+")))
  
  return(ans)
}

