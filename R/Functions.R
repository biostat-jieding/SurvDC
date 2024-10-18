
###############################################################################
# The Main Function:
#   - for both parametric and semiparametric sceniaros
#   - allow for the consideration of a cured fraction
#   - one of the margins can be modeled nonparametric
#   - references:
#       + Czado and Van Keilegom (2023): https://doi.org/10.1093/biomet/asac067
#       + Delhelle and Van Keilegom (in arXiv): https://doi.org/10.48550/arXiv.2403.07963
#       + Ding and Van Keilegom (Manuscript in Preparation)
###############################################################################

# --- the main function for model fitting --- ##
#' @title Survival Analysis under Dependent Censoring using Copula with Unknown Association
#'
#' @description Provide approaches that can be used to model right-censored survival data under dependent censoring (without covariates).
#'   The copula-based approach is adopted and there is no need to explicitly specify the association parameter.
#'   Both marginal distributions of survival and censoring times can be considered as fully parametric, or one of the margins can be modeled nonparametrically.
#'   The existence of a cured fraction concerning survival time can also be taken into consideration.
#'
#' @param yobs a numeric vector that indicated the observed survival times.
#' @param delta a numeric vector that stores the right-censoring indicators.
#' @param tm a numeric vector that contains interested non-negative time points at which the survival probabilities will be evluated.
#'   Note that if we omit the definition of this argument (the default value becomes \code{NULL}), our function will automatically output survival probabilities at all oberserved time points, that is, \code{yobs}.
#' @param copfam a character string that specifies the copula family.
#'   Currently, it supports Archimedean copula families, including \code{"frank"} (the default value), \code{"clayton"}, \code{"gumbel"}, and \code{"joe"}.
#'   The degenerated independent censoring case can be considered as well by setting "indep".
#'   (other options will be added in the near future!)
#' @param margins a list used to define the distribution structures of both the survival and censoring margins.
#'   Specifically, it contains the following elements:
#'   \describe{
#'     \item{\code{survfam}}{ a character string that defines the assumed distribution for the survival time random variable,
#'       including \code{"lnorm"} for log-normal distribution, \code{"weibull"} for weibull distribution (other options will be added in the near future).}
#'     \item{\code{censfam}}{ a character string that defines the assumed distribution for the censoring time random variable, and the details are the same as those shown in \code{survfam}.}
#'     \item{\code{survtrunc}}{ a positive numeric value thats denotes the value of truncation for the assumed distribution, that is, \code{survfam}.}
#'     \item{\code{censtrunc}}{ a positive numeric value thats denotes the value of truncation for the assumed distribution, that is, \code{censfam}.}
#'   }
#'   Note if one of the marginal distributions should be modeled nonparametrically, one can let the corresponding argument to be \code{NULL} directly.
#'   For example if a semiparametric framework that defines the survival margin to be nonparametric and the censoring margin to be parametric, say log-normal, is desired,
#'   we can let \code{survfam = NULL} and \code{censfam = "lnorm"}, which is indeed the default value.
#'   Furthermore, if no truncation is imposed in \code{survfam} (or \code{censfam}), one can directly omit the specification of \code{survtrunc} (or \code{censtrunc}), which is the default specification.
#'   We also remark here that when a cured fraction is included (\code{cure = TRUE}), if \code{survfam} is not \code{NULL} and \code{survtrunc = NULL}, we will automatically let \code{survtrunc} to be \code{max(yobs)}.
#'   If we wants to model the data with a non-truncated survival distribution when there is a cured fraction, we can set \code{survtrunc = Inf}.
#' @param cure a logical value that indicates whether the existence of a cured fraction should be considered.
#' @param Var a list that controls the execution of the bootstrap for variance estimation,
#'     and it contains two elements:
#'     \code{do} is a logical value with default \code{FALSE} to tell the function whether the boostrap-based variances should be calculated;
#'     \code{nboot} is a numeric integer that specifies the number of bootstrap samples.
#' @param control indicates more detailed control of the underlying model fitting procedures.
#'   It is a list of the following three arguments:
#'   \describe{
#'     \item{\code{maxit}}{ a positive integer that denotes the maximum iteration number in optimization. The default value is \code{300}.}
#'     \item{\code{eps}}{ a positive small numeric value that denotes the tolerance for convergence. The default value is \code{1e-6}.}
#'     \item{\code{trace}}{ a logical value that judges whereh the tracing information on the progress of the model fitting should be produced. The default value if \code{TRUE}.}
#'     \item{\code{ktau.inits}}{ a numeric vector that contains initial values of the Kendall's tau.
#'       The default value is \code{NULL}, meaning that a grids of initial values will be automatically generated within our function.}
#'   }
#'
#' @details
#'   Various specifications of marginal distributions can be considered by choosing different combinations of the provided arguments.
#'   Generally speaking, the following scenarios can be realized:
#'   \describe{
#'     \item{\code{parametric survival and censoring margins (without cure)}}{
#'       both \code{survfam} and \code{censfam} are not \code{NULL} and \code{cure = FALSE}.}
#'     \item{\code{parametric survival and censoring margins (with cure)}}{
#'       both \code{survfam} and \code{censfam} are not \code{NULL} and \code{cure = TRUE}.}
#'     \item{\code{nonparametric survival margin and parametric censoring margin (without cure)}}{
#'       \code{survfam = NULL}, \code{censfam} is not \code{NULL} and \code{cure = FALSE}.}
#'     \item{\code{nonparametric survival margin and parametric censoring margin (with cure)}}{
#'       \code{survfam = NULL}, \code{censfam} is not \code{NULL} and \code{cure = TRUE}.}
#'     \item{\code{parametric survival margin and nonparametric censoring margin (without cure)}}{
#'       \code{survfam} is not \code{NULL}, \code{censfam = NULL} and \code{cure = FALSE}.}
#'     \item{\code{parametric survival margin and nonparametric censoring margin (with cure)}}{
#'       \code{survfam} is not \code{NULL}, \code{censfam = NULL} and \code{cure = TRUE}.}
#'   }
#'   References of the underlying methods include: Czado and Van Keilegom (2023) <doi:10.1093/biomet/asac067> and Delhelle and Van Keilegom (2024) <doi:10.48550/arXiv.2403.07963>
#'   Semiparametric modeling framework is proposed by Ding and Van Keilegom with manuscript in preparation.
#'
#' @return
#'   A list of fitted results is returned.
#'   Within this outputted list, the following elements can be found:
#'   \describe{
#'     \item{\code{probs}}{survival probabilities of the survial margin at \code{tm}.}
#'     \item{\code{ktau}}{Kendall's tau.}
#'     \item{\code{parapar}}{estimation of all parameters (except Kendall's tau) contained in the parametric part.}
#'     \item{\code{GoF}}{goodness-of-test results.}
#'     \item{\code{curerate}}{cure rate. If \code{cure = FALSE}, it is \code{NULL}.}
#'   }
#'
#' @importFrom methods is
#' @importFrom stats dlnorm dweibull integrate plnorm pnorm pweibull sd uniroot na.omit
#'
#' @examples
#' \donttest{
#' #----------------------------------------------------------#
#' # Basic preparations before running subsequent examples ####
#' #----------------------------------------------------------#
#'
#' ## clear the environment
#' rm(list=ls(all=TRUE))
#'
#' ## library necessary packages
#' library(SurvDC)
#'
#' #------------------------------------------------------------------------#
#' # simulated data from Frank copula log-Normal margins (without cure)
#' #------------------------------------------------------------------------#
#'
#' ## generate the simulated data
#'
#' # - the sample size of the generated data
#' n <- 2000
#'
#' # - information on the used copula
#' copfam.true <- "frank"
#' ktau.true <- 0.5
#' coppar.true <- ktau.to.coppar(ktau=ktau.true,copfam=copfam.true)
#'
#' # - parameters of the underlying log-normal marginal distributions
#' survpar.true <- c(2.20,1.00)
#' censpar.true <- c(2.20,0.25)
#'
#' # - true underlying survival and censoring times
#' set.seed(1)
#' u.TC <- copula::rCopula(
#'   n        = n,
#'   copula   = copula::archmCopula(
#'     family = copfam.true,
#'     param  = coppar.true,
#'     dim    = 2
#'   )
#' )
#' yobs.T <- qlnorm(1-u.TC[,1],survpar.true[1],survpar.true[2])
#' yobs.C <- qlnorm(1-u.TC[,2],censpar.true[1],censpar.true[2])
#'
#' # - observations
#' yobs  <- pmin(yobs.T,yobs.C)
#' delta <- as.numeric(yobs.T<=yobs.C)
#' cat("censoring rate is",mean(1-delta))
#'
#' ## model the data under different scenarios
#'
#' # - scenario 1: parametric survival and censoring margins
#' set.seed(1)
#' sol.scenario1 <- SurvDC(
#'   yobs    = yobs,
#'   delta   = delta,
#'   tm      = quantile(yobs, c(0.25,0.50,0.75)),
#'   copfam  = copfam.true,
#'   margins = list(survfam = "lnorm", censfam = "lnorm"),
#'   Var     = list(do = TRUE, nboot = 50)
#' )
#' sol.scenario1$probs
#' sol.scenario1$ktau
#' sol.scenario1$parapar
#'
#' # - scenario 2: nonparametric survival margin and parametric censoring margin
#' set.seed(1)
#' sol.scenario2 <- SurvDC(
#'   yobs    = yobs,
#'   delta   = delta,
#'   tm      = quantile(yobs, c(0.25,0.50,0.75)),
#'   copfam  = copfam.true,
#'   margins = list(survfam = NULL, censfam = "lnorm"),
#'   Var     = list(do = TRUE, nboot = 50)
#' )
#' sol.scenario2$probs
#' sol.scenario2$ktau
#' sol.scenario2$parapar
#'
#' # - scenario 3: parametric survival margin and nonparametric censoring margin
#' set.seed(1)
#' sol.scenario3 <- SurvDC(
#'   yobs    = yobs,
#'   delta   = delta,
#'   tm      = quantile(yobs, c(0.25,0.50,0.75)),
#'   copfam  = copfam.true,
#'   margins = list(survfam = "lnorm", censfam = NULL),
#'   Var     = list(do = TRUE, nboot = 50)
#' )
#' sol.scenario3$probs
#' sol.scenario3$ktau
#' sol.scenario3$parapar
#'
#' #------------------------------------------------------------------------#
#' # simulated data from Frank copula log-Normal margins (with cure)
#' #------------------------------------------------------------------------#
#'
#' ## generate the simulated data
#'
#' # - true underlying curerate
#' curerate.true <- 0.2
#'
#' # - true underlying survival and censoring times
#' set.seed(1)
#' u.TC <- copula::rCopula(
#'   n        = n,
#'   copula   = copula::archmCopula(
#'     family = copfam.true,
#'     param  = coppar.true,
#'     dim    = 2
#'   )
#' )
#' yobs.T <- sapply(u.TC[,1],function(uT){
#'   if(uT<=curerate.true){ val <- Inf }else{
#'     val <- EnvStats::qlnormTrunc((1-uT)/(1-curerate.true),survpar.true[1],survpar.true[2],0,15)
#'   }
#'   return(val)
#' })
#' yobs.C <- qlnorm(1-u.TC[,2],censpar.true[1],censpar.true[2])
#' cat("cure rate is",mean(yobs.T==Inf))
#'
#' # - observations
#' yobs  <- pmin(yobs.T,yobs.C)
#' delta <- as.numeric(yobs.T<=yobs.C)
#' cat("censoring rate is",mean(1-delta))
#'
#' ## model the data under different scenarios (with cure)
#'
#' # - scenario 4: parametric survival and censoring margins
#' set.seed(1)
#' sol.scenario4 <- SurvDC(
#'   yobs    = yobs,
#'   delta   = delta,
#'   tm      = quantile(yobs, c(0.25,0.50,0.75)),
#'   copfam  = copfam.true,
#'   margins = list(survfam = "lnorm", censfam = "lnorm"),
#'   Var     = list(do = TRUE, nboot = 50),
#'   cure    = TRUE
#' )
#' sol.scenario4$probs
#' sol.scenario4$ktau
#' sol.scenario4$parapar
#' sol.scenario4$curerate
#'
#' # - scenario 5: nonparametric survival margin and parametric censoring margin
#' set.seed(1)
#' sol.scenario5 <- SurvDC(
#'   yobs    = yobs,
#'   delta   = delta,
#'   tm      = quantile(yobs, c(0.25,0.50,0.75)),
#'   copfam  = copfam.true,
#'   margins = list(survfam = NULL, censfam = "lnorm"),
#'   Var     = list(do = TRUE, nboot = 50),
#'   cure    = TRUE
#' )
#' sol.scenario5$probs
#' sol.scenario5$ktau
#' sol.scenario5$parapar
#' sol.scenario5$curerate
#'
#' # - scenario 6: parametric survival margin and nonparametric censoring margin
#' set.seed(1)
#' sol.scenario6 <- SurvDC(
#'   yobs    = yobs,
#'   delta   = delta,
#'   tm      = quantile(yobs, c(0.25,0.50,0.75)),
#'   copfam  = copfam.true,
#'   margins = list(survfam = "lnorm", censfam = NULL),
#'   Var     = list(do = TRUE, nboot = 50),
#'   cure    = TRUE
#' )
#' sol.scenario6$probs
#' sol.scenario6$ktau
#' sol.scenario6$parapar
#' sol.scenario6$curerate
#' }
#'
#' @export SurvDC
SurvDC <- function(
    yobs,
    delta,
    tm      = NULL,
    copfam  = "frank",
    margins = list(survfam = NULL,    survtrunc = NULL,
                   censfam = "lnorm", censtrunc = NULL),
    cure    = FALSE,
    Var     = list(do = TRUE, nboot = 50),
    control = list(
      maxit      = 300,
      eps        = 1e-6,
      trace      = TRUE,
      ktau.inits = NULL
    )
){

  ## preparations
  n <- length(yobs)
  margins$type <- ifelse(any(c(is.null(margins$survfam),is.null(margins$censfam))),
                         "semiparametric","parametric")
  control.pre <- control.arguments()
  control <- lapply(names(control.pre),function(argi){
    if(is.null(control[[argi]])){
      val <- control.pre[[argi]]
    }else{ val <- control[[argi]] }
    return(val)
  })
  names(control) <- names(control.pre)

  ## set the value of truncation point (if NULL and cure)
  if(is.null(margins$survfam)){
    margins$survtrunc <- NULL
  }else{
    if(is.null(margins$survtrunc) & cure==TRUE){
      margins$survtrunc <- max(yobs)
    }
  }
  if(is.null(margins$censfam)){margins$censtrunc <- NULL}

  ## prepare the constraint matrix and vector
  constraints <- Parameters.Constraints(copfam=copfam,margins=margins,cure=cure)

  ## prepare initial values of the Kendall's tau
  if(is.null(control$ktau.inits)==TRUE | copfam=="indep"){
    if(copfam == "indep"){ control$ktau.inits <- 0 }else{
      ktau.left  <- ifelse(copfam=="frank",-1+1e-5,1e-5)
      control$ktau.inits <- seq(ktau.left,1-1e-5,ifelse(copfam=="frank",0.4,0.2))[-1]
    }
  }

  ## fit the semi-parametric model across many initial values
  if(control$trace==TRUE){cat("- Fitting the Model ...\n")}
  fit.model  <- list(value = -Inf)
  for(iktau in 1:length(control$ktau.inits)){

    # fit the semi-parametric model within try() to avoid error
    cytry <- try({
      sol.c <- Likelihood.Profile.Solve(
        yobs=yobs,delta=delta,copfam=copfam,margins=margins,
        ktau.init=control$ktau.inits[iktau],parapar.init=NULL,
        cure=cure,curerate.init=NULL,constraints=constraints,
        maxit=control$maxit,eps=control$eps
      )
      convergence <- sol.c$convergence
      param       <- sol.c$param
    }, silent = TRUE) # end try

    # - judge the current fit
    if(is(cytry, "try-error") == FALSE){
      if(convergence == TRUE){
        # __ calculate the value of likelihood
        if(length(control$ktau.inits)>1){
          if(margins$type=="parametric"){
            value <- Likelihood.Parametric(param=param,yobs=yobs,delta=delta,copfam=copfam,margins=margins,cure=cure)
          }else{
            value <- Likelihood.Profile.Kernel(param=param,yobs=yobs,delta=delta,copfam=copfam,margins=margins,cure=cure)
          }
        }else{ value <- Inf }
        # __ judge the current fit formally
        if(value > fit.model$value){
          fit.model$param <- param
          fit.model$value <- value
        }
      }
    }

  }

  ## extract estimates and prepare other needed ones

  # parametric part
  if(is.null(fit.model$param)==TRUE){
    stop("Fitting Failed for Current Setting!")
  }else{
    (ktau   <- fit.model$param[1])
    if(cure==TRUE & is.null(margins$survfam)==FALSE){
      parapar <- fit.model$param[-c(1,2)]
      curerate <- fit.model$param[2]
    }else{
      parapar <- fit.model$param[-1]
      curerate <- NULL
    }

  }

  # survival probabilities of survival time at interested points (and cure rate)
  if(is.null(tm)==TRUE){tm <- yobs}
  if(is.null(margins$survfam)==TRUE){
    probs.all <- SurvFunc.CG(tm=tm,yobs=yobs,delta=delta,copfam=copfam,ktau=ktau)$Shat
    if(cure == TRUE){
      curerate <- SurvFunc.CG(tm=max(yobs[delta==1]),yobs=yobs,delta=delta,copfam=copfam,ktau=ktau)$Shat
      probs    <- (probs.all-curerate)/(1-curerate)
    }else{ curerate <- NULL; probs <- probs.all }
  }else{
    if(margins$survfam %in% c("lnorm","weibull")){ parapar.surv <- parapar[1:2] }
    probs <- 1-parafam.p(x=tm,parameter=parapar.surv,distribution=margins$survfam,truncation=margins$survtrunc)
  }

  # calculate goodness-of-fit statistic
  GoF <- SurvDC.GoF(yobs=yobs,delta=delta,copfam=copfam,margins=margins,ktau=ktau,parapar=parapar,cure=cure,curerate=curerate)

  ## bootstrap procedure for variance estimation
  if(Var$do==TRUE){

    if(control$trace==TRUE){cat("- Doing Bootstrap to Quantify Variation ...\n")}

    # preparations

    # start the bootstrap
    results.boots <- rep(list(list()),Var$nboot)
    iboot <- 1; nwrong <- 0
    while(iboot <= Var$nboot){

      if(nwrong > ceiling(Var$nboot*0.1)){stop("Too Many Failed Fittings in Bootstrap!")}
      if(control$trace==TRUE){if(iboot%%10==0){cat("  + Bootstrap",iboot,"/",Var$nboot,"\n")}}

      # - generate the current bootstrap sample
      idx.iboot   <- sort(sample(1:n,n,replace=TRUE))
      yobs.iboot  <- yobs[idx.iboot]
      delta.iboot <- delta[idx.iboot]

      # - fit the current model and store the results
      try.iboot <- try({

        # __ fit the model and extract estimates
        sol.iboot <- Likelihood.Profile.Solve(
          yobs=yobs.iboot,delta=delta.iboot,
          copfam=copfam,margins=margins,
          ktau.init=ktau,parapar.init=parapar,
          cure=cure,curerate.init=curerate,
          constraints=constraints,
          maxit=control$maxit,eps=control$eps
        )
        convergence.iboot <- sol.iboot$convergence
        ktau.iboot        <- sol.iboot$param[1]
        if(cure==TRUE & is.null(margins$survfam)==FALSE){
          parapar.iboot  <- sol.iboot$param[-c(1,2)]
          curerate.iboot <- sol.iboot$param[2]
        }else{
          parapar.iboot  <- sol.iboot$param[-1]
          curerate.iboot <- NULL
        }

        # __ survival probabilities at interested points (and cure rate)
        if(is.null(margins$survfam)==TRUE){
          probs.all.iboot <- SurvFunc.CG(tm=tm,yobs=yobs.iboot,delta=delta.iboot,copfam=copfam,ktau=ktau.iboot)$Shat
          if(cure == TRUE){
            curerate.iboot <- SurvFunc.CG(tm=max(yobs.iboot[delta.iboot==1]),
                                          yobs=yobs.iboot,delta=delta.iboot,copfam=copfam,ktau=ktau.iboot)$Shat
            probs.iboot    <- (probs.all.iboot-curerate.iboot)/(1-curerate.iboot)
          }else{ curerate.iboot <- NULL; probs.iboot <- probs.all.iboot }
        }else{
          if(margins$survfam %in% c("lnorm","weibull")){ parapar.surv.iboot <- parapar.iboot[1:2] }
          probs.iboot <- 1-parafam.p(x=tm,parameter=parapar.surv.iboot,distribution=margins$survfam,truncation=margins$survtrunc)
        }

        # __ calculate goodness-of-fit statistic
        GoF.iboot <- SurvDC.GoF(
          yobs=yobs.iboot,delta=delta.iboot,copfam=copfam,
          margins=margins,ktau=ktau.iboot,parapar=parapar.iboot,
          cure=cure,curerate=curerate.iboot)

      }, silent = T)

      # - judge the current fit
      if(is(try.iboot,"try-error") == FALSE){
        if(convergence.iboot == TRUE){
          results.boots[[iboot]]$ktau    <- ktau.iboot
          results.boots[[iboot]]$parapar <- parapar.iboot
          results.boots[[iboot]]$probs   <- probs.iboot
          results.boots[[iboot]]$curerate <- curerate.iboot
          results.boots[[iboot]]$GoF <- GoF.iboot
        }else{ nwrong <- nwrong + 1; next }
      }else{ nwrong <- nwrong + 1; next }

      # - prepare for next iteration
      iboot <- iboot + 1

    }

    # get estimtes of standard errors
    ktau.se       <- sd(do.call(c,lapply(1:Var$nboot,function(iboot){results.boots[[iboot]]$ktau})))
    parapar.se    <- apply(do.call(cbind,lapply(1:Var$nboot,function(iboot){results.boots[[iboot]]$parapar})),1,sd)
    probs.se      <- apply(do.call(cbind,lapply(1:Var$nboot,function(iboot){results.boots[[iboot]]$probs})),1,sd)
    if(cure==TRUE){
      curerate.se <- sd(do.call(c,lapply(1:Var$nboot,function(iboot){results.boots[[iboot]]$curerate})))
    }

    # calculate p-value for goodness-of-fit statistic
    GoF.boots  <- do.call(c,lapply(1:Var$nboot,function(iboot){results.boots[[iboot]]$GoF}))
    GoF.pvalue <- mean(GoF.boots>=GoF)

  }else{
    ktau.se <- parapar.se <- probs.se <- curerate.se <- GoF.pvalue <- NULL
  }

  ## summary all fitted results

  # - Kendall's tau
  ktau <- c(Est=ktau,SE=ktau.se,zvalue=ktau/ktau.se,pvalue=2*(1-pnorm(abs(ktau/ktau.se))))

  # - parametric parameters
  parapar <- as.data.frame(cbind(
    Est    = parapar,
    SE     = parapar.se,
    zvalue = parapar/parapar.se,
    pvalue = 2*(1-pnorm(abs(parapar/parapar.se)))
  ))

  # - survival probabilities
  probs <- as.data.frame(cbind(
    time   = tm,
    Est    = probs,
    SE     = probs.se,
    zvalue = probs/probs.se,
    pvalue = 2*(1-pnorm(abs(probs/probs.se)))
  ))[order(tm),,drop=FALSE]

  # - cure rate
  if(cure==TRUE){
    curerate <- c(Est=curerate,SE=curerate.se,zvalue=curerate/curerate.se,pvalue=2*(1-pnorm(abs(curerate/curerate.se))))
  }else{
    curerate <- NULL
  }

  # - goodness-of-fit
  GoF <- c(statistic=GoF,pvalue=GoF.pvalue)

  ## output
  out <- list(
    probs     = probs,
    ktau      = ktau,
    parapar   = parapar,
    GoF       = GoF,
    curerate  = curerate,
    inputs    = list(
      yobs    = yobs,
      delta   = delta,
      copfam  = copfam,
      margins = margins
    )
  )
  return(out)

}

###############################################################################
# Other Auxiliary Functions
###############################################################################

## --- an auxiliary function used to fit a (semi-) parametric model --- ##
Likelihood.Profile.Solve <- function(
    yobs,delta,copfam,margins,ktau.init,parapar.init,cure,curerate.init,constraints,maxit,eps
){

  # prepare initial values for parametric parameters (if it is not nonparametric)
  if(is.null(parapar.init)==TRUE | copfam=="indep"){

    # calculate initial estimates for survival or censoring part
    for(fam in c("survfam","censfam")){

      # - if this part is NULL, omit the following calculations
      if(is.null(margins[[fam]])==TRUE){next}

      # - preparations
      if(fam=="survfam"){
        truncation <- margins$survtrunc
      }else{
        truncation <- margins$censtrunc
      }
      delta.parametric <- ifelse(rep(fam=="survfam",length(yobs)),delta,1-delta)

      # - obtain initial values of parametric parameters by assuming independence
      sol.fam.init <- SurvMLE(
        yobs=yobs,delta=delta.parametric,distribution=margins[[fam]],
        truncation=truncation,cure=(cure & fam=="survfam"),maxit=maxit)
      parapar.init <- c(parapar.init,sol.fam.init$parapar)
      if(cure==TRUE & fam=="survfam"){
        curerate.init <- sol.fam.init$curerate
      }

    }

  }

  # fit the model formally
  if(copfam=="indep"){
    if(cure==TRUE & is.null(margins$survfam)==FALSE){
      param <- c(0,curerate.init,parapar.init)
    }else{
      param <- c(0,parapar.init)
    }
    convergence <- TRUE
  }else{

    # - initial values
    if(cure==TRUE & is.null(margins$survfam)==FALSE){
      param.init <- c(ktau.init,curerate.init,parapar.init)
    }else{
      param.init <- c(ktau.init,parapar.init)
    }

    # - fit the model formally
    if(margins$type=="parametric"){

      sol.model <- stats::constrOptim(
        theta   = param.init,
        f       = Likelihood.Parametric,
        grad    = NULL,
        ui      = constraints$A,
        ci      = constraints$b,
        control = list(trace=0,fnscale=-1,maxit=maxit),
        yobs    = yobs, delta = delta, copfam = copfam, margins = margins, cure = cure
      )
      param <- sol.model$par
      convergence <- (sol.model$convergence == 0)

    }else if(margins$type=="semiparametric"){

      param.old <- param.init
      delta.nonparametric <- ifelse(rep(is.null(margins$survfam),length(yobs)),delta,1-delta)
      Syobs.old <- SurvFunc.CG(tm=NULL,yobs=yobs,delta=delta.nonparametric,
                               copfam=copfam,ktau=param.old[1],coppar=NULL)$Shat
      numit     <- 1
      repeat{
        sol.model <- stats::constrOptim(
          theta   = param.old,
          f       = Likelihood.Semiparametric,
          grad    = NULL,
          ui      = constraints$A,
          ci      = constraints$b,
          control = list(trace=0,fnscale=-1,maxit=maxit),
          Syobs   = Syobs.old, yobs = yobs, delta = delta, copfam = copfam,
          margins = margins, cure = cure
        )
        param  <- sol.model$par
        Syobs  <- SurvFunc.CG(tm=NULL,yobs=yobs,delta=delta.nonparametric,copfam=copfam,ktau=param[1],coppar=NULL)$Shat
        epsmax <- max(c(abs(param-param.old),abs(Syobs-Syobs.old)))
        if( epsmax>eps & numit<maxit ){
          param.old <- param
          Syobs.old <- Syobs
          numit     <- numit + 1
        }else{
          break
        }
      }
      convergence <- (numit<maxit)

    }

  }

  # output
  out <- list(
    param       = param,
    convergence = convergence
  )
  return(out)

}

## --- calculate the goodness-of-fit test statistic --- ##
SurvDC.GoF <- function(yobs,delta,copfam,margins,ktau,parapar,cure=FALSE,curerate=NULL){

  coppar  <- ktau.to.coppar(ktau=ktau,copfam=copfam)

  # distribution-related functions
  if(margins$type=="parametric"){

    # - length of parameters for two parts
    if(margins$survfam %in% c("lnorm","weibull")){ parapar.surv.num <- 2 }

    # - probabilities for both survival and censoring parts
    STC <- NULL
    for(fam in c("survfam","censfam")){
      if(fam=="survfam"){
        parapar.temp <- parapar[1:parapar.surv.num]
        truncation.temp <- margins$survtrunc
      }else{
        parapar.temp <- parapar[-c(1:parapar.surv.num)]
        truncation.temp <- margins$censtrunc
      }
      S.temp <- 1-parafam.p(x=yobs,parameter=parapar.temp,distribution=margins[[fam]],truncation=truncation.temp)
      STC <- cbind(STC,S.temp)
    }

    # - define probabilities formally
    if(cure==TRUE){
      ST <- curerate + (1-curerate)*STC[,1]
    }else{
      ST <- STC[,1]
    }
    SC <- STC[,2]

  }else if(margins$type=="semiparametric"){

    # - probabilities for parametric part time (judge later)
    S.temp <- 1-parafam.p(
      x=yobs,parameter=parapar,
      distribution=c(margins$survfam,margins$censfam),
      truncation=c(margins$survtrunc,margins$censtrunc)
    )

    # - probabilities for both survival and censoring parts
    if(is.null(margins$survfam) == TRUE){
      ST <- SurvFunc.CG(tm=NULL,yobs=yobs,delta=delta,copfam=copfam,ktau=NULL,coppar=coppar)$Shat
      SC <- S.temp
    }else{
      if(cure==TRUE){
        ST <- curerate + (1-curerate)*S.temp
      }else{
        ST <- S.temp
      }
      SC <- SurvFunc.CG(tm=NULL,yobs=yobs,delta=1-delta,copfam=copfam,ktau=NULL,coppar=coppar)$Shat
    }

  }

  # goodness-of-fit statistic
  # if(margins$type=="parametric"){
  #   SY.Model <- ST + SC + copdist.Archimedean(x=cbind(1-ST,1-SC),copfam=copfam,ktau=NULL,coppar=coppar) - 1
  # }else if(margins$type=="semiparametric"){
  #   SY.Model <- copdist.Archimedean(x=cbind(ST,SC),copfam=copfam,ktau=NULL,coppar=coppar)
  # }
  SY.Model <- copdist.Archimedean(x=cbind(ST,SC),copfam=copfam,ktau=NULL,coppar=coppar)
  SY.Empir <- sapply(yobs,function(yobsi){mean(yobs>yobsi)})
  GoF      <- sum((SY.Model-SY.Empir)^2)

  # output
  return(GoF)

}


## --- likelihood for a given parametric distribution
SurvMLE.Likelihood <- function(param,yobs,delta,distribution,truncation=NULL,cure=FALSE){

  # extract parameters
  if(cure==TRUE){
    curerate <- param[1]
    parameter <- param[-1]
  }else{
    parameter <- param
  }

  # basic distribution function
  FT.0 <- parafam.p(x=yobs[delta==0],parameter=parameter,distribution=distribution,truncation=truncation)
  fT.1 <- parafam.d(x=yobs[delta==1],parameter=parameter,distribution=distribution,truncation=truncation)

  # revision based on the cure status
  if(cure==TRUE){
    FT.0 <- FT.0*(1-curerate)
    fT.1 <- fT.1*(1-curerate)
  }

  # calculate likelihood
  Lik1  <- sum(log(pmax(1e-20,fT.1)),na.rm=TRUE)
  Lik0  <- sum(log(pmax(1e-20,1-FT.0)),na.rm=TRUE)
  Lik   <- Lik1 + Lik0

  # return
  return(Lik)

}

## --- maximum likelihood estimator for a given parametric distribution
SurvMLE <- function(yobs,delta,distribution,truncation=NULL,cure=FALSE,maxit=300){

  # setup initial values and constraints
  if(distribution %in% c("lnorm")){
    A <- rbind(c(0,1))
    b <- 1e-5
    param.init <- c(0,1)
  }else if(distribution %in% c("weibull")){
    A <- rbind(c(1,0),c(0,1))
    b <- c(1e-5,1e-5)
    param.init <- c(1,1)
  }

  # add information about cure
  if(cure==TRUE){
    A = as.matrix(Matrix::bdiag(rbind(1,-1),A))
    b = c(1e-5,-1+1e-5,b)
    param.init <- c(0.5,param.init)
  }

  # maximizing the likelihood function
  sol.model <- stats::constrOptim(
    theta   = param.init,
    f       = SurvMLE.Likelihood,
    grad    = NULL,
    ui      = A,
    ci      = b,
    control = list(trace=0,fnscale=-1,maxit=maxit),
    yobs    = yobs, delta = delta, distribution = distribution,
    truncation = truncation, cure = cure
  )
  if(cure==TRUE){
    curerate <- sol.model$par[1]
    parapar  <- sol.model$par[-1]
  }else{
    curerate <- NULL
    parapar  <- sol.model$par
  }
  convergence <- (sol.model$convergence == 0)

  # output
  out <- list(
    parapar     = parapar,
    curerate    = curerate,
    convergence = convergence
  )

}

## --- adjustment value of truncation --- ##
parafam.trunc <- function(truncation,parameter,distribution){

  if(distribution=="lnorm"){
    val.adjust <- plnorm(truncation,parameter[1],parameter[2])
  }else if(distribution=="weibull"){
    val.adjust <- pweibull(truncation,parameter[1],parameter[2])
  }
  return(val.adjust)

}

## --- value of distribution function --- ##
parafam.p <- function(x,parameter,distribution,truncation=NULL){

  # values with no truncation
  if(distribution=="lnorm"){
    val <- plnorm(x,parameter[1],parameter[2])
  }else if(distribution=="weibull"){
    val <- pweibull(x,parameter[1],parameter[2])
  }

  # truncation
  if(is.null(truncation)==FALSE){
    val.adjust <- parafam.trunc(truncation=truncation,parameter=parameter,distribution=distribution)
    val <- pmin(val/val.adjust,1)
  }

  # output
  return(val)

}

## --- value of density function --- ##
parafam.d <- function(x,parameter,distribution,truncation=NULL){

  # values with no truncation
  if(distribution=="lnorm"){
    val <- dlnorm(x,parameter[1],parameter[2])
  }else if(distribution=="weibull"){
    val <- dweibull(x,parameter[1],parameter[2])
  }

  # truncation
  if(is.null(truncation)==FALSE){
    val.adjust <- parafam.trunc(truncation=truncation,parameter=parameter,distribution=distribution)
    val <- val/val.adjust
    val[x>truncation] <- 0
  }

  # output
  return(val)

  # EnvStats::plnormTrunc(x,parameter[1],parameter[2],0,truncation)
  # EnvStats::dlnormTrunc(x,parameter[1],parameter[2],0,truncation)

}

## --- simplified version of profiled likelihood function: for semiparametric modeling --- ##
Likelihood.Semiparametric <- function(param,Syobs,yobs,delta,copfam,margins,cure=FALSE){

  # extract current values
  if(cure==TRUE & is.null(margins$survfam)==FALSE){
    curerate <- param[2]
    parapar  <- param[-c(1,2)]
  }else{
    parapar <- param[-1]
  }
  ktau    <- param[1]
  coppar  <- ktau.to.coppar(ktau=ktau,copfam=copfam)

  # - probabilities for parametric part time (judge later)
  delta.parametric <- ifelse(rep(is.null(margins$survfam),length(yobs)),1-delta,delta)
  S.temp <- 1-parafam.p(x=yobs,parameter=parapar,distribution=c(margins$survfam,margins$censfam),truncation=c(margins$survtrunc,margins$censtrunc))
  f.temp <- parafam.d(x=yobs[delta.parametric==1],parameter=parapar,distribution=c(margins$survfam,margins$censfam),truncation=c(margins$survtrunc,margins$censtrunc))

  # distribution-related functions
  if(is.null(margins$survfam)==TRUE){

    # basic quantities
    ST   <- Syobs
    SC   <- S.temp
    fC.0 <- f.temp

    # calculate likelihood
    hCT.1 <- cophfunc(x=cbind(ST,SC)[delta==1,,drop=FALSE],coppar=coppar,copfam=copfam,condvar=1)
    hTC.0 <- cophfunc(x=cbind(ST,SC)[delta==0,,drop=FALSE],coppar=coppar,copfam=copfam,condvar=2)
    Lik1  <- sum(log(pmax(1e-20,hCT.1)),na.rm=TRUE)
    Lik0  <- sum(log(pmax(1e-20,fC.0*hTC.0)),na.rm=TRUE)
    Lik   <- Lik1 + Lik0

  }else{

    # basic quantities
    SC <- Syobs
    if(cure==TRUE){
      ST <- curerate + (1-curerate)*S.temp
      fT.1 <- (1-curerate)*f.temp
    }else{
      ST <- S.temp
      fT.1 <- f.temp
    }

    # calculate likelihood
    hCT.1 <- cophfunc(x=cbind(ST,SC)[delta==1,,drop=FALSE],coppar=coppar,copfam=copfam,condvar=1)
    hTC.0 <- cophfunc(x=cbind(ST,SC)[delta==0,,drop=FALSE],coppar=coppar,copfam=copfam,condvar=2)
    Lik1  <- sum(log(pmax(1e-20,fT.1*hCT.1)),na.rm=TRUE)
    Lik0  <- sum(log(pmax(1e-20,hTC.0)),na.rm=TRUE)
    Lik   <- Lik1 + Lik0

  }

  # return
  return(Lik)

}

## --- parametric likelihood function: for fully parametric modeling --- ##
Likelihood.Parametric <- function(param,yobs,delta,copfam,margins,cure=FALSE){

  # extract current values
  ktau <- param[1]
  if(cure==TRUE){
    curerate <- param[2]
    parapar  <- param[-c(1,2)]
  }else{
    parapar <- param[-1]
  }
  coppar <- as.vector(ktau.to.coppar(ktau=ktau,copfam=copfam))

  # - survival and censoring time
  if(margins$survfam %in% c("lnorm","weibull")){
    parapar.surv <- parapar[1:2]
    parapar.cens <- parapar[-c(1:2)]
  }

  # distribution-related functions

  # - for survival time

  # __ basic distribution function
  FT   <- parafam.p(x=yobs,parameter=parapar.surv,distribution=margins$survfam,truncation=margins$survtrunc)
  fT.1 <- parafam.d(x=yobs[delta==1],parameter=parapar.surv,distribution=margins$survfam,truncation=margins$survtrunc)

  # __ revision based on the cure status
  if(cure==TRUE){
    FT   <- FT*(1-curerate)
    fT.1 <- fT.1*(1-curerate)
  }

  # - for censoring time
  FC   <- parafam.p(x=yobs,parameter=parapar.cens,distribution=margins$censfam,truncation=margins$censtrunc)
  fC.0 <- parafam.d(x=yobs[delta==0],parameter=parapar.cens,distribution=margins$censfam,truncation=margins$censtrunc)

  # calculate likelihood
  hCT.1 <- cophfunc(x=cbind(1-FT,1-FC)[delta==1,,drop=FALSE],coppar=coppar,copfam=copfam,condvar=1)
  hTC.0 <- cophfunc(x=cbind(1-FT,1-FC)[delta==0,,drop=FALSE],coppar=coppar,copfam=copfam,condvar=2)
  Lik1  <- sum(log(pmax(1e-20,fT.1*hCT.1)),na.rm=TRUE)
  Lik0  <- sum(log(pmax(1e-20,fC.0*hTC.0)),na.rm=TRUE)
  Lik   <- Lik1 + Lik0

  # return
  return(Lik)

}

## --- generate constraints --- ##
Parameters.Constraints <- function(copfam,margins,cure){

  # produce constraints
  if(copfam == "indep"){
    constraints <- NULL
  }else{

    # - bounds for the Kendall's tau
    ktau.left  <- ifelse(copfam=="frank",-1+1e-5,1e-5)
    ktau.right <- 1-1e-5
    constraints <- list(
      A = rbind(1,-1),
      b = c(ktau.left,-ktau.right)
    )

    # - bounds for the cure rate
    if(cure == TRUE & is.null(margins$survfam) == FALSE){
      constraints <- list(
        A = as.matrix(Matrix::bdiag(constraints$A,rbind(1,-1))),
        b = c(constraints$b,1e-5,-1+1e-5)
      )
    }

    # - define the constraints for other parametric parts

    for(fam in c("survfam","censfam")){

      # - if this part is NULL, omit the following calculations
      if(is.null(margins[[fam]])==TRUE){next}

      # - add constraints based on different parametric families
      if(margins[[fam]] %in% c("lnorm")){
        constraints <- list(
          A = as.matrix(Matrix::bdiag(constraints$A,rbind(c(0,1)))),
          b = c(constraints$b,1e-5)
        )
      }else if(margins[[fam]] %in% c("weibull")){
        constraints <- list(
          A = as.matrix(Matrix::bdiag(constraints$A,rbind(c(1,0),c(0,1)))),
          b = c(constraints$b,1e-5,1e-5)
        )
      }

    }

  }

  # output
  return(constraints)

}

## --- convert kendall's tau into copula parameter --- ##
#' @title Convert the Kendall's tau into the copula parameter
#'
#' @param ktau a numeric value that denotes the Kendall's tau.
#' @param copfam a character string that denotes the copula family.
#'
#' @export ktau.to.coppar
ktau.to.coppar <- function(ktau,copfam){

  # ktau's allowed range:
  # frank   : (-1,1)
  # gumbel  : [0,1)
  # clayton : (0,1)
  # joe     : [0,1)

  if(copfam=="frank"){

    coppar <- uniroot(
      f=function(coppar,ktau){
        if(abs(coppar)>=1e-5){
          int  <- integrate(f=function(t){t/(exp(t)-1)},0,coppar)$value
          ktauc <- 1-4*(1-int/coppar)/coppar
        }else{ ktauc <- 0 }
        int  <- integrate(f=function(t){t/(exp(t)-1)},0,coppar)$value
        ktauc - ktau
      },
      interval=c(-1e3,1e3),ktau=ktau
    )$root

  }else if(copfam=="gumbel"){

    coppar <- 1/(1-ktau)

  }else if(copfam=="joe"){

    coppar <- uniroot(
      f=function(coppar,ktau){
        int <- integrate(f=function(t,coppar){exp(-2*t)*t*(1-exp(-t))^(2/coppar-2)},0,Inf,coppar)$value
        1 - 4*int/(coppar^2) - ktau
      },
      interval=c(1,1e3),ktau=max(ktau,1e-5)
    )$root

  }else if(copfam=="clayton"){

    coppar <- 2*ktau/(1-ktau)

  }else if(copfam=="indep"){

    coppar <- NULL

  }

  return(coppar)

}

## --- convert copula parameter into kendall's tau --- ##
#' @title Convert the copula parameter the Kendall's tau
#'
#' @param coppar a numeric value that denotes the copula parameter.
#' @param copfam a character string that denotes the copula family.
#'
#' @export ktau.to.coppar
coppar.to.ktau <- function(coppar,copfam){

  # coppar's allowed range:
  # frank   : [-infty,infty]
  # gumbel  : [1,infty)
  # clayton : (0,infty)
  # joe     : [1,infty)

  if(copfam=="frank"){

    if(abs(coppar)>=1e-5){
      int  <- integrate(f=function(t){t/(exp(t)-1)},0,coppar)$value
      ktau <- 1-4*(1-int/coppar)/coppar
    }else{
      ktau <- 0
    }

  }else if(copfam=="gumbel"){

    ktau <- (coppar-1)/coppar

  }else if(copfam=="joe"){

    int  <- integrate(f=function(t,coppar){exp(-2*t)*t*(1-exp(-t))^(2/coppar-2)},0,Inf,coppar)$value
    ktau <- 1-4*int/(coppar^2)

  }else if(copfam=="clayton"){

    ktau <- coppar/(coppar+2)

  }else if(copfam=="indep"){

    ktau <- 0

  }

  return(ktau)

}

## --- Archimedean copula: generator function --- ##
generator.Archimedean <- function(x,coppar,copfam,inverse=FALSE){

  if(copfam=="frank"){
    if(inverse==FALSE){
      val <- -log((exp(-coppar*x)-1)/(exp(-coppar)-1))
    }else{
      val <- -log(1+exp(-x)*(exp(-coppar)-1))/coppar
    }
  }else if(copfam=="gumbel"){
    if(inverse==FALSE){
      val <- (-log(x))^coppar
    }else{
      val <- exp(-x^(1/coppar))
    }
  }else if(copfam=="joe"){
    if(inverse==FALSE){
      val <- -log(1-(1-x)^coppar)
    }else{
      val <- 1-(1-exp(-x))^(1/coppar)
    }
  }else if(copfam=="clayton"){
    if(inverse==FALSE){
      val <- (x^(-coppar)-1)/coppar
    }else{
      val <- (1+x*coppar)^(-1/coppar)
    }
  }else if(copfam=="indep"){
    if(inverse==FALSE){
      val <- -log(x)
    }else{
      val <- exp(-x)
    }
  }

  return(val)

}

## --- copula: h-function --- ##
cophfunc <- function(x,coppar,copfam,condvar=1){

  # preparations
  x <- pmin(pmax(x,1e-20),1-1e-20)
  numx <- nrow(x)
  if(condvar==1){
    u <- x[,1]; v <- x[,2]
  }else{
    u <- x[,2]; v <- x[,1]
  }

  # get the value of h function under different copulas
  if(copfam=="frank"){
    val <- exp(-coppar*u)*(exp(-coppar*v)-1)/(exp(-coppar)-1+(exp(-coppar*u)-1)*(exp(-coppar*v)-1))
  }else if(copfam=="gumbel"){
    uvpar <- (-log(u))^coppar + (-log(v))^coppar
    val <- exp(-uvpar^(1/coppar))*(uvpar)^(1/coppar-1)*(-log(u))^(coppar-1)/u
  }else if(copfam=="joe"){
    val <- ((1-u)^coppar+(1-v)^coppar-(1-u)^coppar*(1-v)^coppar)^(1/coppar-1)*(1-u)^(coppar-1)*(1-(1-v)^coppar)
  }else if(copfam=="clayton"){
    val <- u^(-coppar-1)*(u^(-coppar)+v^(-coppar)-1)^(-1/coppar-1)
  }else if(copfam=="indep"){
    val <- v
  }

  # output
  return(val)

}

## --- Archimedean copula: distribution function --- ##
copdist.Archimedean <- function(x,copfam,ktau,coppar=NULL){

  # preparations
  u <- x[,1]; v <- x[,2]
  if(is.null(coppar)==TRUE){ coppar <- ktau.to.coppar(ktau=ktau,copfam=copfam) }

  # distribution function
  val <- generator.Archimedean(
    x=generator.Archimedean(x=u,coppar,copfam)+generator.Archimedean(x=v,coppar,copfam),
    coppar,copfam,inverse=TRUE)

  # output
  return(val)

}

## --- estimated survival function based on copula-graphic estimator (Archimedean copula only) --- ##
SurvFunc.CG <- function(tm=NULL,yobs,delta,copfam,ktau,coppar=NULL){

  # preparation
  n <- length(yobs)
  yobs.order <- order(yobs)
  if(is.null(coppar)==TRUE){ coppar <- ktau.to.coppar(ktau=ktau,copfam=copfam) }

  # calculate values of survival function
  pai.yobs      <- (n-rank(yobs,ties.method="first")+1)/n
  pai.yobs.left <- pmax((n-rank(yobs,ties.method="first"))/n,0)
  diff.yobs     <- ifelse(delta==0,0,generator.Archimedean(x=pai.yobs,coppar=coppar,copfam=copfam)-generator.Archimedean(x=pai.yobs.left,coppar=coppar,copfam=copfam))
  if(is.null(tm)==FALSE){
    sum.tm <- sapply(tm,function(tmi){sum(diff.yobs[yobs<=tmi])})
  }else{
    sum.tm <- cumsum(diff.yobs[yobs.order])
    sum.tm[yobs.order] <- sum.tm
  }
  Shat <- generator.Archimedean(x=-sum.tm,coppar=coppar,copfam=copfam,inverse=TRUE)

  # return
  out <- list(Shat = Shat)
  return(out)

}

## --- estimated survival function based on Kaplan-Meier estimator --- ##
SurvFunc.KM <- function(tm=NULL,yobs,delta,type="right"){

  ### preparation ###
  if(is.null(tm)==TRUE){tm <- yobs}
  N <- length(yobs)
  yobs.order <- order(yobs)
  yobs.sort <- yobs[yobs.order]
  delta.sort <- delta[yobs.order]
  yobs.1max <- max(yobs[delta==1])

  # calculate the values of KM at pre-specified time points tm
  prods <- 1-delta.sort/(N-(1:N)+1)
  if(type=="right"){
    KMt <- sapply(tm,function(tmi){prod(prods[yobs.sort<=tmi])}) # right-continuous
  }else{
    KMt <- sapply(tm,function(tmi){prod(prods[yobs.sort<tmi])}) # left-continuous
  }

  # output
  return(KMt)

}

## --- profiled likelihood function with kernel smoothing: for semiparametric modeling --- ##
Likelihood.Profile.Kernel <- function(param,yobs,delta,copfam,margins,cure=FALSE){

  # extract current values
  if(cure==TRUE & is.null(margins$survfam)==FALSE){
    curerate <- param[2]
    parapar <- param[-c(1,2)]
  }else{
    parapar <- param[-1]
  }
  ktau    <- param[1]
  coppar  <- ktau.to.coppar(ktau=ktau,copfam=copfam)

  # probabilities for parametric part time (judge later)
  delta.parametric <- ifelse(rep(is.null(margins$survfam),length(yobs)),1-delta,delta)
  S.temp <- 1-parafam.p(x=yobs,parameter=parapar,distribution=c(margins$survfam,margins$censfam),truncation=c(margins$survtrunc,margins$censtrunc))
  f.temp <- parafam.d(x=yobs[delta.parametric==1],parameter=parapar,distribution=c(margins$survfam,margins$censfam),truncation=c(margins$survtrunc,margins$censtrunc))

  # distribution-related functions
  if(is.null(margins$survfam)){

    # - for survival time
    ST           <- SurvFunc.CG(tm=NULL,yobs=yobs,delta=delta,copfam=copfam,ktau=NULL,coppar=coppar)$Shat
    yobs.1.order <- order(yobs[delta==1])
    yobs.1.sort  <- yobs[delta==1][yobs.1.order]
    ST.1.sort    <- ST[delta==1][yobs.1.order]
    jumpT.1      <- c(1-ST.1.sort[1],-diff(ST.1.sort))
    h            <- (8*sqrt(2)/3)^(1/5)*sd(yobs.1.sort)*length(yobs)^(-1/5)
    fT.1         <- sapply(yobs[delta==1],function(x){sum(Kernel((x-yobs.1.sort)/h)*jumpT.1)/h})

    # - for censoring time
    SC   <- S.temp
    fC.0 <- f.temp

  }else{

    # - for survival time
    if(cure==TRUE){
      ST <- curerate + (1-curerate)*S.temp
      fT.1 <- (1-curerate)*f.temp
    }else{
      ST   <- S.temp
      fT.1 <- f.temp
    }

    # - for censoring time
    SC           <- SurvFunc.CG(tm=NULL,yobs=yobs,delta=1-delta,copfam=copfam,ktau=NULL,coppar=coppar)$Shat
    yobs.0.order <- order(yobs[delta==0])
    yobs.0.sort  <- yobs[delta==0][yobs.0.order]
    SC.0.sort    <- SC[delta==0][yobs.0.order]
    jumpC.0      <- c(1-SC.0.sort[1],-diff(SC.0.sort))
    h            <- (8*sqrt(2)/3)^(1/5)*sd(yobs.0.sort)*length(yobs)^(-1/5)
    fC.0         <- sapply(yobs[delta==0],function(x){sum(Kernel((x-yobs.0.sort)/h)*jumpC.0)/h})

  }

  # calculate likelihood
  hCT.1 <- cophfunc(x=cbind(ST,SC)[delta==1,,drop=FALSE],coppar=coppar,copfam=copfam,condvar=1)
  hTC.0 <- cophfunc(x=cbind(ST,SC)[delta==0,,drop=FALSE],coppar=coppar,copfam=copfam,condvar=2)
  Lik1  <- sum(log(pmax(1e-20,fT.1*hCT.1)),na.rm=TRUE)
  Lik0  <- sum(log(pmax(1e-20,fC.0*hTC.0)),na.rm=TRUE)
  Lik   <- Lik1 + Lik0

  # return
  return(Lik)

}

## --- define different kernel functions
Kernel <- function(u, name="Gaussian"){
  # the kernel function (not re-scaled) with only one continuous variable

  if(name=="Uniform"){
    ret <- 0.5*(abs(u)<=1)
  }else if(name=="Gaussian"){
    ret <- exp(-u^2/2)/sqrt(2*pi)
  }else if(name=="Epanechnikov"){
    ret <- 0.75*(1-u^2)*(abs(u)<=1)
  }else if(name=="Quartic"){
    ret <- (15/16)*(1-u^2)^2*(abs(u)<=1)
  }
  return(ret)
}

## --- initial values within the control --- ##
control.arguments <- function(
    maxit      = 300,
    eps        = 1e-6,
    trace      = TRUE,
    ktau.inits = NULL
){
  return(list(
    maxit      = maxit,
    eps        = eps,
    trace      = trace,
    ktau.inits = ktau.inits
  ))
}
