#' Convert ARDL to lm
#'
#' Converts an ARDL formula to a lm formula
#'
#' @param formula.ardl the ARDL formula object
#' @param case ARDL model case for the treatment of intercept and trend
#' @param d number of variables
#' @keywords internal
#' @export
ardl_to_lm= function(formula.ardl,case,d){
  char.eqn=as.character(formula.ardl)
  str.pieces=gsub(" ", "",strsplit(char.eqn[3],split="[+]")[[1]], fixed = TRUE)
  str.pieces.lm=str.pieces
  for(i in 1:length(str.pieces)){
    if((substr(str.pieces[i],1,2)!="L(")){
      str.pieces.lm[i]=str.pieces[i]
    }else{
      str.pieces.lm[i]=paste0(substr(str.pieces[i],3,nchar(str.pieces[i])-3),".l",substr(str.pieces[i],nchar(str.pieces[i])-1,nchar(str.pieces[i])-1),collapse="")
    }
  }

  subvec_D=(str.pieces.lm[which(substr(str.pieces.lm,1,2)=="D_")])
  subvec_L=(str.pieces.lm[which(substr(str.pieces.lm,1,2)!="D_")])[1:d]
  subvec_DL=(subvec_D[which(substr(subvec_D,nchar(subvec_D)-2,nchar(subvec_D)-1)==".l")])
  subvec_D0 = (subvec_D[which(substr(subvec_D,nchar(subvec_D)-2,nchar(subvec_D)-1)!=".l")])
  timechar = ifelse(case>=4,"timeid",NA)
  ordered_formula = na.omit(c(subvec_L,subvec_DL,subvec_D0,timechar))

  ans = as.formula(paste(char.eqn[2],char.eqn[1],paste0(ordered_formula,collapse="+")))

  return(ans)
}
