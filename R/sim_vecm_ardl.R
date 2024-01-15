#' Generate data from a VECM/ARDL equation
#'
#' @param nobs number of observations.
#' @param case case related to intercept and trend
#' @param sigma.in error covariance matrix \eqn{\boldsymbol\Sigma}
#' @param gamma.in list of VECM short-run parameter matrices \eqn{\boldsymbol\Gamma_j}
#' @param axx.in long-run relationships between the independent variables \eqn{\mathbf A_{xx}}
#' @param ayx.uc.in long-run unconditional relationship between dependent and independent variables, \eqn{\mathbf a_{yx}}.
#' The second component ayxC, derived from conditioning, is calculated as \eqn{\mathbf a_{yx}^{(C)}= - \boldsymbol\omega'\mathbf A_{xx}}
#' @param ayy.in long-run relationship for the dependent variable \eqn{a_{yy}}
#' @param mu.in VAR intercept vector \eqn{\boldsymbol\mu} (CASE II)
#' @param eta.in VAR trend vector \eqn{\boldsymbol\eta} (CASE IV)
#' @param azero.in VECM intercept \eqn{\boldsymbol{\alpha}_{0}} (CASE III-IV-V)
#' @param aone.in VECM trend \eqn{\boldsymbol{\alpha}_{1}} (CASE V)
#' @param burn.in burn-in number of observations
#' @param seed.in optional seed number for random error generation.
#' @return A list that includes \itemize{
#'\item \code{dims}: a vector with the dataset dimension
#'\item \code{case}: the case given as input
#'\item \code{data}: the generated data
#'\item \code{diffdata}: the data first difference
#'\item \code{ut}: the generated random error matrix.
#'\item \code{sigma}: the error covariance matrix \eqn{\boldsymbol\Sigma}.
#'\item \code{omega}: the \eqn{\boldsymbol\omega} vector of parameters generated via conditioning
#'\item \code{at}: the conditional long-run parameter matrix \eqn{\widetilde{\mathbf A}}
#'\item \code{ayy}: the coefficient weighting the EC term, \eqn{a_{yy}}
#'\item \code{ayx.uc}: the unconditional subvector of the ARDL equation \eqn{\mathbf a_{yx}}
#'\item \code{ayx2}: the conditioning effect \eqn{\omega'A_{xx}}
#'\item \code{ayx.c}: the conditional subvector of the ARDL equation \eqn{\widetilde{a}_{y.x}=a_{yx}-\omega'A_{xx}}
#'\item \code{gammalist}: the list of unconditional \eqn{\boldsymbol\Gamma_j} parameter matrices
#'\item \code{psilist}: the list of conditional \eqn{\boldsymbol\gamma_{y.x,j}} parameter matrices
#'\item \code{vmu}: the VAR intercept \eqn{\boldsymbol\mu}
#'\item \code{veta}: the VAR trend \eqn{\boldsymbol\eta}
#'\item \code{azero}: the unconditional VECM intercept \eqn{\boldsymbol\alpha_0}
#'\item \code{aone}: the unconditional VECM trend \eqn{\boldsymbol\alpha_1}
#'\item \code{azero.c}: the conditional VECM intercept \eqn{\boldsymbol\alpha_0^c}
#'\item \code{aone.c}: the conditional VECM trend \eqn{\boldsymbol\alpha_1^c}
#'\item \code{interc.ardl}: the conditional ARDL intercept \eqn{\alpha_{0.y}} (case > 2)
#'\item \code{trend.ardl}: the conditional ARDL trend \eqn{\alpha_{1.y}} (case = 5)
#'\item \code{theta0}: the \eqn{\theta_0} coefficient in the EC term (case = 2)
#'\item \code{theta1}: the \eqn{\theta_1} coefficient in the EC term (case = 4)
#'\item \code{interc.ec}: the conditional ARDL intercept derived from the EC tem \eqn{\alpha^{EC}_{0.y}} (case = 2)
#'\item \code{trend.ec}: the conditional ARDL trend derived from the EC tem \eqn{\alpha^{EC}_{1.y}} (case = 4)}
#' @examples
#'#PARAMETERS
#'
#'#Sigma
#'corrm = matrix(0, ncol = 3, nrow = 3)
#'corrm[2,1] = 0.25
#'corrm[3,1] = 0.4
#'corrm[3,2] = -0.25
#'corrs = (corrm + t(corrm)) + diag(3)
#'sds = diag(c(1.3, 1.2, 1))
#'sigma = (sds %*% corrs %*% t(sds))
#'
#'#Gamma
#'gammax = list()
#'gammax[[1]] = matrix(c(0.6, 0, 0.2, 0.1, -0.3, 0, 0, -0.3, 0.2), nrow = 3, ncol = 3, byrow = TRUE)
#'gammax[[2]] = matrix(c(0.2, 0, 0.1, 0.05, -0.15, 0, 0, 0, 0.1), nrow = 3, ncol = 3, byrow = TRUE)
#'
#'#DATA GENERATION
#'data_sim = sim_vecm_ardl(nobs = 200,
#'                          case = 3,
#'                          sigma.in = sigma,
#'                          gamma.in = gammax,
#'                          axx.in = matrix(c(0.3, 0.5, 0.4, 0.3), nrow = 2, ncol = 2),
#'                          ayx.uc.in = c(0.5,0.6),
#'                          ayy.in = 0.7,
#'                          mu.in = rep(0, 3),
#'                          eta.in = rep(0, 3),
#'                          azero.in = rep(0.4, 3),
#'                          aone.in = rep(0, 3),
#'                          burn.in = 50,
#'                          seed.in = 10)
#'
#' @export
sim_vecm_ardl =
  function(nobs,
           case = 1,
           sigma.in,
           gamma.in,
           axx.in,
           ayx.uc.in,
           ayy.in,
           mu.in = NULL,
           eta.in = NULL,
           azero.in = NULL,
           aone.in = NULL,
           burn.in = nobs * 0.5,
           seed.in = NULL){

    if (!all(eigen(sigma.in)$values >= 0) |
        !isSymmetric.matrix(sigma.in)) {
      stop("Invalid covariance matrix.")
    }
    
    if (!inherits(gamma.in,"list")) {
      stop("gamma.in must be a list.")
    }
    
    if (case < 1 | case > 5 | case %% 1 != 0) {
      warning("Inadmissible case. Defaulting to case III.")
      case = 3
    }

    if (burn.in <= 0) {
      warning("Inadmissible burn-in value. Defaulting to n * 0.5.")
      burn.in = n * 0.5
    }
    
    if(is.null(mu.in)){
      mu.in = rep(0,ncol(sigma.in))
    }
    
    if(is.null(eta.in)){
      eta.in = rep(0,ncol(sigma.in))
    }
    
    if(is.null(azero.in)){
      azero.in = rep(0,ncol(sigma.in))
    }
    
    if(is.null(aone.in)){
      aone.in = rep(0,ncol(sigma.in))
    }

    if (case == 1 &&
        (any(mu.in != 0) || any(azero.in != 0) || any(eta.in != 0) || any(aone.in != 0))) {
      stop("CASE I implies no intercept and no trend, please check parameter values")
    }

    if (case == 2 && (any(azero.in != 0) || any(eta.in != 0) || any(aone.in != 0))) {
      stop("CASE II implies restricted intercept and no trend, please check parameter values")
    }

    if (case == 3 && (any(mu.in != 0) || any(eta.in != 0) || any(aone.in != 0))) {
      stop("CASE III implies unrestricted intercept and no trend, please check parameter values.")
    }

    if (case == 4 && (any(mu.in != 0) || any(aone.in != 0))) {
      stop(
        "CASE IV implies unrestricted intercept and restricted trend, please check parameter values."
      )
    }
    
    if (case == 5 && (any(mu.in != 0) || any(eta.in != 0))) {
      stop(
        "CASE V implies unrestricted intercept and unrestricted trend, please check parameter values."
      )
    }

    cs = ncol(sigma.in)
    rs = nrow(sigma.in)
    cg = unlist(lapply(gamma.in, ncol))
    rg = unlist(lapply(gamma.in, nrow))
    dmu = length(mu.in)
    deta = length(eta.in)
    dazero = length(azero.in)
    daone = length(aone.in)
    da1 = length(ayy.in)+length(ayx.uc.in)
    da2 = length(ayy.in)+ncol(axx.in)
    da3 = length(ayy.in)+nrow(axx.in)
    
    if(length(ayy.in)!=1){
      stop("Invalid ayy entry, must be a scalar")
    }

    dims = c(cs, rs, cg, rg, dmu, deta, dazero, daone, da1, da2, da3)
    
    dimsval=dims[dims>0]

    if (length(unique(dimsval)) != 1) {
      stop("dimensions of input elements do not match.")
    }

    if(!is.null(seed.in)){
      set.seed(seed.in)
      }

    d = unique(dimsval)

    sigma = sigma.in
    gammax = gamma.in
    Mlag = length(gammax)
    omegat = sigma[1,-1] %*% solve(sigma[-1,-1])
    psi = list()

    for (j in 1:Mlag) {
      psi[[j]] = matrix(gammax[[j]][1,], nrow = 1, ncol = d) -
        omegat %*% matrix(gammax[[j]][-1,], nrow = (d - 1), ncol = d)
    }
    
    ayy = ayy.in
    Axx = axx.in
    ayx1 = ayx.uc.in
    ayx2 = omegat %*% Axx
    ayx = ayx1 - ayx2
    Ax = rbind(ayx1, Axx)
    At= cbind(c(ayy.in, rep(0, d - 1)), Ax)
    F1 = diag(d) - Reduce('+', gammax)
    
    #VAR intercept
    vmu = mu.in * (case > 1)
    veta = eta.in * (case > 3)

    #unconditional VECM intercept
    azero = (At %*% vmu) * (case == 2) + azero.in * (case > 2)
    aone = (At %*% veta) * (case == 4) + aone.in * (case == 5)

    # conditional A first row
    A.c = At[1,] - matrix(c(0, omegat %*% At[-1,-1]), nrow = 1)
    F1.c = F1[1,] - matrix(c(0, omegat %*% F1[-1,-1]), nrow = 1)

    #conditional VECM intercept
    azero.c0 = azero 
    aone.c0 = aone
    azero.c0[1] = azero[1] - omegat %*% azero[-1]
    aone.c0[1] = aone[1] - omegat %*% aone[-1]
    
    theta0 = theta1 = 0
    thetav = rep(0, d - 1)
    
    if(case == 2){
      thetav = -ayx/(ayy + 0.000001)
      theta0 = mu.in[1] - thetav %*% mu.in[-1] 
    }
    
    if(case == 4){
      thetav = -ayx/(ayy + 0.000001)
      theta1 = eta.in[1] - thetav %*% eta.in[-1] 
    }
    
    nss = nobs + burn.in + Mlag + 1

    ut = matrix(rnorm(n = nss * d), nss, d) %*% chol(sigma)

    azero.c = azero.c0[1]
    aone.c = aone.c0[1]
    ut.c = ut[, 1] - ut[,-1] %*% t(omegat)

    #Example 3 variables and 2 max lags

    #dep_0 = yt
    #ind.1_0 = x1t
    #ind.2_0 = x2t
    #dep_1 = yt-1
    #ind.1_1 = x1t-1
    #ind.2_1 = x2t-1

    #d_dep = dyt
    #d_ind.1_0 = dx1t
    #d_ind.2_0 = dx2t
    #d_dep_1 = dyt-1
    #d_ind.1_1 = dx1t-1
    #d_ind.2_1 = dx2t-1
    #d_dep_2 = dyt-2
    #d_ind.1_2 = dx1t-2
    #d_ind.2_2 = dx2t-2

    df.oss = data.frame(matrix(0, nrow = nss, ncol = (2 * d + #levels + past
                                                        d * (Mlag + 1) + #differences
                                                        1))) #trend index

    typev = c("dep", "ind")
    nind = 1:(d - 1)
    nlag = 0:Mlag

    #names of level variables
    df.levcol = data.frame(names = apply(
      expand.grid(typev, nind, nlag), 1, paste, collapse = "_"
    )) %>%
      dplyr::filter(!stringr::str_detect(names, "^dep_[2-9]"))
    lev.col = as.character(df.levcol[1:(d * 2),])

    #names of differenced variables
    df.diffcol = data.frame(names = apply(
      expand.grid(typev, nind, nlag), 1, paste, collapse = "_"
    )) %>%
      dplyr::filter(!stringr::str_detect(names, "^dep_[2-9]"))
    d.col = paste0("d_", as.character(df.diffcol[1:(d * (Mlag + 1)),]
    ))

    cnames = c(lev.col, d.col, "time")
    colnames(df.oss) = cnames
    for (w in 2:nss) {
      if (w > 2) {
        for (j in 0:Mlag - 1) {
          for (n in 1:Mlag) {
            
            df.oss[w, stringr::str_detect(cnames, paste0("^d_(.*)[", j + n, "]$"))] =
              df.oss[w - n, stringr::str_detect(cnames, paste0("^d_(.*)[", j, "]$"))]
            df.oss[w, stringr::str_detect(cnames, paste0("^(?!d_)(.*)[", j + n, "]$"))] =
              df.oss[w - n, stringr::str_detect(cnames, paste0("^(?!d_)(.*)[", j, "]$"))]
          }
        }
      }

      #filling X level variables
      #select Z lagged diffs
      v.dlag = df.oss[w, stringr::str_detect(cnames, "^d_(.*)[1-9]$")]

      #select Z lagged levels
      v.lag = df.oss[w, stringr::str_detect(cnames, "^(?!d_)(.*)[1-9]$")]

      #create X instant diffs
      df.oss[w, stringr::str_detect(cnames, "^d_ind_\\d*_0$")] =
        azero[-1] + aone[-1] * (w - burn.in - Mlag - 1) + ut[w,-1] +   #errore + trend + intercetta per x
        do.call(cbind, gammax)[-1,] %*% matrix(as.numeric(v.dlag),
                                               nrow = ncol(v.dlag),
                                               ncol = 1) - #diff in lag per z
        At[-1,] %*% matrix(as.numeric(v.lag),
                           nrow = ncol(v.lag),
                           ncol = 1) #lag levels per z lungo periodo

      #create X levels
      df.oss[w, stringr::str_detect(cnames, "^(?!d)(.*)[0]$")] =
        df.oss[(w - 1), stringr::str_detect(cnames, "^(?!d)(.*)[0]$")] + #riga precedente
        df.oss[w, stringr::str_detect(cnames, "^d_ind_\\d*_0$")] #riga appena creata

      #filling y level variables

      #select X instant diffs
      x.diff = df.oss[w, stringr::str_detect(cnames, "^(d_ind)(.*)[_0]$")]

      #create y instant diffs
      df.oss$d_dep_1_0[w] =
        azero.c + #intercept
        aone.c * (w - burn.in - Mlag - 1) + #trend
        ut.c[w] - #error term
        A.c %*% (t(as.matrix(v.lag))) + #Z lagged levels
        omegat %*% t(as.matrix(x.diff)) + #X instant diffs
        c(unlist(psi)) %*% t(as.matrix(v.dlag)) #Z lagged diffs

      #create y levels
      df.oss$dep_1_0[w] =
        df.oss$dep_1_0[w - 1] + #riga precedente
        df.oss$d_dep_1_0[w] #riga appena creata

      #id for trend
      df.oss$time[w] = (w - burn.in - Mlag - 1)
    }

    #level variables matrix
    z = df.oss[(burn.in + Mlag + 1 + 1):nss, stringr::str_detect(cnames, "^(?!d_)(.*)[_0]$")]

    #diff variables matrix
    dz = df.oss[(burn.in + Mlag + 1 + 1):nss, stringr::str_detect(cnames, "^(d_)(.*)[_0]$")]

    #add time id
    z=cbind(z,time=df.oss[(burn.in + Mlag + 1 + 1):nss,ncol(df.oss)])
    dz=cbind(dz,time=df.oss[(burn.in + Mlag + 1 + 1):nss,ncol(df.oss)])

    return(list(
        dims = c(nrow(z), ncol(z)),
        case = case,
        data = z,
        diffdata = dz,
        ut = ut,
        sigma = sigma,
        omega = omegat,
        at = At,
        ayy = ayy.in,
        ayx.uc = ayx1,
        ayx.2 = ayx2,
        ayx.c = ayx,
        gammalist = gammax,
        psilist = psi,
        vmu = vmu,
        veta = veta,
        azero = azero,
        aone = aone,
        azero.c = azero.c0,
        aone.c = aone.c0,
        interc.ardl = azero.c * (case > 2),
        trend.ardl = aone.c * (case > 4),
        theta0 = theta0 * (case == 2),
        theta1 = theta1 * (case == 4),
        interc.ec = ayy * theta0 * (case == 2),
        trend.ec = ayy * theta1 * (case == 4)
      )
    )
  }

