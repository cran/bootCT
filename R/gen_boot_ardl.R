#' Call to cpp function to generate data following a conditional ARDL/VECM model
#'
#' @param r_in matrix of residuals
#' @param GAMMAX short term parameters
#' @param A long term parameters
#' @param start_z initial observations
#' @param omegat parameter vector implied in conditioning
#' @param interc intercept vector
#' @param trend trend vector
#' @return No return value, calls cpp function \code{boot_ardl_c}
#' @keywords internal
#' @export
gen_boot_ardl = function(r_in,
                         GAMMAX,
                         A,
                         start_z,
                         omegat,
                         interc,
                         trend){

  return(.Call(`_bootCT_boot_ardl_c`, r_in, GAMMAX,A,start_z,omegat,interc,trend))
}
