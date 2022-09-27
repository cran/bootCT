#' Generate data from a VECM/ARDL equation
#'
#' @param nobs number of observations.
#' @param case case related to intercept and trend
#' @param sigma.in error covariance matrix.
#' @param gamma.in list of short-run parameter matrices
#' @param Axx.in long-run relationships between the independent variables
#' @param ayxUC.in long-run unconditional relationship between dependent and independent variables, \eqn{\mathbf a_{yx}^{(UC)}} .
#' The second component ayxC, derived from conditioning, is calculated as\eqn{\mathbf a_{yx}^{(C)}= - \boldsymbol\omega'\mathbf A_{xx}}
#' @param ayy.in long-run relationship for the dependent variable \eqn{a_{yy}}
#' @param mu.in VAR intercept vector
#' @param eta.in VAR trend vector
#' @param azeroy.in Conditional ARDL intercept. Overridden if CASE I or CASE II
#' @param aoney.in Conditional ARDL trend. Overridden if CASE IV
#' @param burn.in burn-in number of observations
#' @param seed.in optional seed number for random error generation.
#'
#' @return A list that includes \itemize{
#'\item \code{dims}: a vector with the dataset dimension
#'\item \code{case}: the case given as input
#'\item \code{data}: the generated data
#'\item \code{diffdata}: the data first difference
#'\item \code{ut}: the generated random error matrix.
#'\item \code{sigma}: the error covariance matrix \eqn{\boldsymbol\Sigma}.
#'\item \code{omega}: the \eqn{\boldsymbol\omega} vector of parameters generated via conditioning
#'\item \code{At}: the conditional long-run parameter matrix \eqn{\tilde{\mathbf A}}
#'\item \code{ayx1}: the unconditional subvector of the ARDL equation \eqn{\mathbf a_{y.x}^{UC}}
#'\item \code{ayx}: the conditional subvector of the ARDL equation \eqn{a_{y.x}=a_{y.x}^{UC}-\omega'A_{xx}}
#'\item \code{gammalist}: the list of unconditional \eqn{\boldsymbol\Gamma_j} parameter matrices
#'\item \code{psilist}: the list of conditional \eqn{\boldsymbol\psi_{y.x,j}} parameter matrices
#'\item \code{azero}: the unconditional VECM intercept
#'\item \code{azero.c}: the conditional VECM intercept
#'\item \code{interc.ardl}: the conditional ARDL intercept
#'\item \code{aone}: the unconditional VECM trend
#'\item \code{aone.c}: the conditional VECM trend
#'\item \code{interc.ardl}: the conditional ARDL trend
#'\item \code{vmu}: the VAR intercept
#'\item \code{veta}: the VAR trend}
#' @examples
#'#PARAMETERS
#'
#'#Sigma
#'corrm = matrix(0, ncol = 3, nrow = 3)
#'corrm[2,1] = 0.25
#'corrm[3,1] = 0.4
#'corrm[3,2] = -0.25
#'Corrm = (corrm + t(corrm)) + diag(3)
#'sds = diag(c(1.3, 1.2, 1))
#'Sigma = (sds %*% Corrm %*% t(sds))
#'
#'#Gamma
#'gammax=list()
#'gammax[[1]] = matrix(c(0.6, 0, 0.2, 0.1, -0.3, 0, 0, -0.3, 0.2), nrow = 3, ncol = 3, byrow = TRUE)
#'gammax[[2]] = matrix(c(0.2, 0, 0.1, 0.05, -0.15, 0, 0, 0, 0.1), nrow = 3, ncol = 3, byrow = TRUE)
#'
#'#DATA GENERATION
#'data_sim = sim_vecm_ardl(nobs = 200,
#'                          case = 3,
#'                          sigma.in = Sigma,
#'                          gamma.in = gammax,
#'                          Axx.in = matrix(c(0.3, 0.5, 0.4, 0.3), nrow = 2, ncol = 2),
#'                          ayxUC.in = c(0.5,0.6),
#'                          ayy.in = 0.7,
#'                          mu.in = rep(0.3, 3),
#'                          eta.in = rep(0, 3),
#'                          azeroy.in = 0.4,
#'                          aoney.in = 0,
#'                          burn.in = 50,
#'                          seed.in = 10)
#'
#' @export
sim_vecm_ardl =
  function(nobs,
           case = 1,
           sigma.in = diag(3),
           gamma.in,
           Axx.in,
           ayxUC.in,
           ayy.in,
           mu.in,
           eta.in,
           azeroy.in = 0,
           aoney.in = 0,
           burn.in,
           seed.in = NULL){

    if (case < 1 | case > 5 | case %% 1 != 0) {
      warning("Inadmissible case. Defaulting to case III.")
      case = 3
    }

    if (burn.in <= 0) {
      warning("Inadmissible burn-in value. Defaulting to n * 0.5.")
      burn.in = n * 0.5
    }

    if (case == 1 &&
        (any(mu.in != 0) || any(eta.in != 0) || azeroy.in != 0 || aoney.in != 0)) {
      warning("CASE I implies no intercept and no trend, parameter values overridden")
    }

    if (case == 2 &&
        (any(eta.in != 0) || azeroy.in != 0 || aoney.in != 0)) {
      warning("CASE II implies restricted intercept and no trend, parameter values overridden")
    }

    if (case == 3 && (any(eta.in != 0) || aoney.in != 0)) {
      warning("CASE III implies unrestricted intercept and no trend, parameter values overridden.")
    }

    if (case == 4 && (aoney.in != 0)) {
      warning(
        "CASE IV implies unrestricted intercept and restricted trend, parameter values overridden."
      )
    }

    if (!all(eigen(sigma.in)$values >= 0) |
        !isSymmetric.matrix(sigma.in)) {
      stop("Invalid covariance matrix.")
    }

    if (!inherits(gamma.in,"list")) {
      stop("gamma.in must be a list.")
    }

    if (length(azeroy.in) > 1 || length(aoney.in) > 1) {
      warning(
        "Intercept / trend values of conditional ARDL model must be scalar. Using their first element."
      )
      azeroy.in = azeroy.in[1]
      aoney.in = aoney.in[1]
    }

    cs = ncol(sigma.in)
    rs = nrow(sigma.in)

    cg = unlist(lapply(gamma.in, ncol))
    rg = unlist(lapply(gamma.in, nrow))
    dmu = length(mu.in)
    deta = length(eta.in)

    dims = c(cs, rs, cg, rg, dmu, deta)

    if (length(unique(dims)) != 1) {
      stop("dimensions of input elements do not match.")
    }

    if(!is.null(seed.in)){
      set.seed(seed.in)
      }

    d = unique(dims)

    sigma = sigma.in
    gammax = gamma.in
    Mlag = length(gammax)
    omegat = sigma[1,-1] %*% solve(sigma[-1,-1])
    psi = list()

    for (j in 1:Mlag) {
      psi[[j]] = matrix(gammax[[j]][1,], nrow = 1, ncol = d) -
        omegat %*% matrix(gammax[[j]][-1,], nrow = (d - 1), ncol = d)
    }

    Axx = Axx.in
    ayx1 = ayxUC.in
    ayx2 = omegat %*% Axx
    ayx = ayx1 - ayx2
    Ax = rbind(ayx1, Axx)
    At = cbind(c(ayy.in, rep(0, d - 1)), Ax)

    PI = -At

    Id = diag(ncol(PI))
    F1 = Id - Reduce('+', gammax)

    #VAR intercept
    vmu = mu.in * (case > 1)
    veta = eta.in * (case > 3)

    #unconditional VECM intercept
    azero = -PI %*% vmu + (PI + F1) %*% veta
    aone = -PI %*% veta

    azero.c0 = azero
    aone.c0 = aone

    # conditional A first row
    A.c = At[1,] - matrix(c(0, omegat %*% At[-1,-1]), nrow = 1)
    F1.c = F1[1,] - matrix(c(0, omegat %*% F1[-1,-1]), nrow = 1)

    #conditional ARDL intercept
    if (case == 2) {
      azero.c0[1] =  A.c[1, ] %*% vmu + (A.c[1, ] + F1.c[1, ]) %*% veta
    } else{
      azero.c0[1] = azeroy.in
    }

    if (case == 4) {
      aone.c0[1] = A.c[1,] %*% veta
    } else{
      aone.c0[1] = aoney.in
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
                                               ncol = 1) + #diff in lag per z
        PI[-1,] %*% matrix(as.numeric(v.lag),
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
        Sigma = sigma,
        omega = omegat,
        AMat = At,
        PIMat = PI,
        ayx1 = ayx1,
        ayx = ayx,
        Glist = gammax,
        Psilist = psi,
        azero = azero,
        azero.c = azero.c0,
        interc.ardl = azero.c[1],
        aone = aone,
        aone.c = aone.c0,
        trend.ardl = aone.c[1],
        vmu = vmu,
        veta = veta
      )
    )
  }

