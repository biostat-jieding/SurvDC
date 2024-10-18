# *SurvDC* *R* *Package*: Survival Analysis under Dependent Censoring using Copula with Unknown Association

This is an *R* package for providing approaches that can be used to model **right-censored survival data** under **dependent censoring** (without covariates). 
- Key Technique: **COPULA**. 
- There is ***no need** to explicitly specify the **association parameter***.
- Allow for *flexible modeling frameworks* (with more details shown later):
  - Both marginal distributions are parametric (**fully parametric scenario**).
  - One of the marginal distributions is nonparametric (**semiparametric scenario**).
- The existence of **a cured fraction** concerning survival time can be considered.
- References of the underlying methods include: Czado and Van Keilegom (2023) <doi:10.1093/biomet/asac067> and Delhelle and Van Keilegom (2024) <doi:10.48550/arXiv.2403.07963>. Semiparametric modeling framework is proposed by Ding and Van Keilegom with manuscript in preparation.

## Package description and included main functions

### Installation
Installation of this package can be done locally after downloading the package manually from this github website. 
We will also upload this package to the Comprehensive R Archive Network (CRAN) so that it can be downloaded as a standard R package. 
Currently, it can be loaded using *R* command
```R
devtools::install_github("biostat_jieding/SurvDC")
```

### Main function and its arguments
The main function included in our *R* package is *SurvDC()*  and it can be called via the following *R* command:
```R
SurvDC(
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
)
```
We refer to its help page for more detailed explanations of the corresponding arguments (typing *?SurvDC()*). 
Here, we provide a brief introduction of them:
- **yobs** a numeric vector that indicated the observed survival times.
- **delta** a numeric vector that stores the right-censoring indicators.
- **tm** a numeric vector that contains interested non-negative time points at which the survival probabilities will be evluated. Note that if we omit the definition of this argument (the default value becomes *NULL*), our function will automatically output survival probabilities at all oberserved time points, that is, *yobs*.
- **copfam** a character string that specifies the copula family. Currently, it supports Archimedean copula families, including *"frank"* (the default value), *"clayton"*, *"gumbel"*, and *"joe"*. The degenerated independent censoring case can be considered as well by setting *"indep"*. (other options will be added in the near future!)
- **margins** a list used to define the distribution structures of both the survival and censoring margins. Specifically, it contains the following elements:
  - *survfam* a character string that defines the assumed distribution for the survival time random variable, including *"lnorm"* for log-normal distribution, *"weibull"* for weibull distribution (other options will be added in the near future).
  - *censfam* a character string that defines the assumed distribution for the censoring time random variable, and the details are the same as those shown in *survfam*.
  - *survtrunc* a positive numeric value thats denotes the value of truncation for the assumed distribution, that is, *survfam*.
  - *censtrunc* a positive numeric value thats denotes the value of truncation for the assumed distribution, that is, *censfam*.
Note if one of the marginal distributions should be modeled nonparametrically, one can let the corresponding argument to be *NULL* directly. For example if a semiparametric framework that defines the survival margin to be nonparametric and the censoring margin to be parametric, say log-normal, is desired, we can let *survfam = NULL* and *censfam = "lnorm"*, which is indeed the default value. Furthermore, if no truncation is imposed in *survfam* (or *censfam*), one can directly omit the specification of *survtrunc* (or *censtrunc*), which is the default specification. We also remark here that when a cured fraction is included (*cure = TRUE*), if *survfam* is not *NULL* and *survtrunc = NULL*, we will automatically let *survtrunc* to be *max(yobs)*.
#'   If we wants to model the data with a non-truncated survival distribution when there is a cured fraction, we can set *survtrunc = Inf*.
- **cure** a logical value that indicates whether the existence of a cured fraction should be considered.
- **Var** a list that controls the execution of the bootstrap for variance estimation, and it contains two elements: *do* is a logical value with default \code{FALSE} to tell the function whether the boostrap-based variances should be calculated; *nboot* is a numeric integer that specifies the number of bootstrap samples.
- **control** indicates more detailed control of the underlying model fitting procedures. It is a list of the following three arguments:
  - *maxit* a positive integer that denotes the maximum iteration number in optimization. The default value is *300*.
  - *eps* a positive small numeric value that denotes the tolerance for convergence. The default value is *1e-6*.
  - *trace* a logical value that judges whereh the tracing information on the progress of the model fitting should be produced. The default value if *TRUE*.
  - *ktau.inits* a numeric vector that contains initial values of the Kendall's tau. The default value is *NULL*, meaning that a grids of initial values will be automatically generated within our function.

### Tractable modeling frameworks
We remark here that two essential arguments are *margins* and *cure* which can help us realize the follow different modeling frameworks:
- *parametric survival and censoring margins (without cure)*
  - both *survfam* and *censfam* are not *NULL and *cure = FALSE*.
- *parametric survival and censoring margins (with cure)*
  - both *survfam* and *censfam* are not *NULL* and *cure = TRUE*.
- *nonparametric survival margin and parametric censoring margin (without cure)*
  - *survfam = NULL*, *censfam* is not *NULL* and *cure = FALSE*.
- *nonparametric survival margin and parametric censoring margin (with cure)*
  - *survfam = NULL*, *censfam* is not *NULL* and *cure = TRUE*.
- *parametric survival margin and nonparametric censoring margin (without cure)*
  - *survfam* is not *NULL*, *censfam = NULL* and *cure = FALSE*.
- *parametric survival margin and nonparametric censoring margin (with cure)*
  - *survfam* is not *NULL*, *censfam = NULL* and *cure = TRUE*.


## Numerical illustrations

We are going to use simulated datasets to illustrate the useages of our *R* package from different perspectives.
Before executing the following examples, we need to library the package:
```R
library(SurvDC)
```

### simulated data from Frank copula log-Normal margins (without cure)

#### generate the simulated data

the sample size of the generated data
```R
n <- 2000
```

information on the used copula
```R
copfam.true <- "frank"
ktau.true <- 0.5
coppar.true <- ktau.to.coppar(ktau=ktau.true,copfam=copfam.true)
```

parameters of the underlying log-normal marginal distributions
```R
survpar.true <- c(2.20,1.00)
censpar.true <- c(2.20,0.25)
```

true underlying survival and censoring times
```R
set.seed(1)
u.TC <- copula::rCopula(
  n        = n,
  copula   = copula::archmCopula(
    family = copfam.true,
    param  = coppar.true,
    dim    = 2
  )
)
yobs.T <- qlnorm(1-u.TC[,1],survpar.true[1],survpar.true[2])
yobs.C <- qlnorm(1-u.TC[,2],censpar.true[1],censpar.true[2])
```

observations
```R
yobs  <- pmin(yobs.T,yobs.C)
delta <- as.numeric(yobs.T<=yobs.C)
cat("censoring rate is",mean(1-delta))
```

#### model the data under different scenarios

scenario 1: parametric survival and censoring margins
```R
set.seed(1)
sol.scenario1 <- SurvDC(
  yobs    = yobs, 
  delta   = delta, 
  tm      = quantile(yobs, c(0.25,0.50,0.75)), 
  copfam  = copfam.true,
  margins = list(survfam = "lnorm", censfam = "lnorm"),
  Var     = list(do = TRUE, nboot = 50)
)
sol.scenario1$probs
sol.scenario1$ktau
sol.scenario1$parapar
```

scenario 2: nonparametric survival margin and parametric censoring margin
```R
set.seed(1)
sol.scenario2 <- SurvDC(
  yobs    = yobs, 
  delta   = delta, 
  tm      = quantile(yobs, c(0.25,0.50,0.75)), 
  copfam  = copfam.true,
  margins = list(survfam = NULL, censfam = "lnorm"),
  Var     = list(do = TRUE, nboot = 50)
)
sol.scenario2$probs
sol.scenario2$ktau
sol.scenario2$parapar
```

scenario 3: parametric survival margin and nonparametric censoring margin
```R
set.seed(1)
sol.scenario3 <- SurvDC(
  yobs    = yobs, 
  delta   = delta, 
  tm      = quantile(yobs, c(0.25,0.50,0.75)), 
  copfam  = copfam.true,
  margins = list(survfam = "lnorm", censfam = NULL),
  Var     = list(do = TRUE, nboot = 50)
)
sol.scenario3$probs
sol.scenario3$ktau
sol.scenario3$parapar
```



### simulated data from Frank copula log-Normal margins (with cure)

#### generate the simulated data

true underlying curerate
```R
curerate.true <- 0.2
```

true underlying survival and censoring times
```R
set.seed(1)
u.TC <- copula::rCopula(
  n        = n,
  copula   = copula::archmCopula(
    family = copfam.true,
    param  = coppar.true,
    dim    = 2
  )
)
yobs.T <- sapply(u.TC[,1],function(uT){
  if(uT<=curerate.true){ val <- Inf }else{
    val <- EnvStats::qlnormTrunc((1-uT)/(1-curerate.true),survpar.true[1],survpar.true[2],0,15)
  }
  return(val)
})
yobs.C <- qlnorm(1-u.TC[,2],censpar.true[1],censpar.true[2])
cat("cure rate is",mean(yobs.T==Inf))
```

observations
```R
yobs  <- pmin(yobs.T,yobs.C)
delta <- as.numeric(yobs.T<=yobs.C)
cat("censoring rate is",mean(1-delta))
```

#### model the data under different scenarios (with cure)

scenario 4: parametric survival and censoring margins
```R
set.seed(1)
sol.scenario4 <- SurvDC(
  yobs    = yobs, 
  delta   = delta, 
  tm      = quantile(yobs, c(0.25,0.50,0.75)), 
  copfam  = copfam.true,
  margins = list(survfam = "lnorm", censfam = "lnorm"),
  Var     = list(do = TRUE, nboot = 50),
  cure    = TRUE
)
sol.scenario4$probs
sol.scenario4$ktau
sol.scenario4$parapar
sol.scenario4$curerate
```

scenario 5: nonparametric survival margin and parametric censoring margin
```R
set.seed(1)
sol.scenario5 <- SurvDC(
  yobs    = yobs, 
  delta   = delta, 
  tm      = quantile(yobs, c(0.25,0.50,0.75)), 
  copfam  = copfam.true,
  margins = list(survfam = NULL, censfam = "lnorm"),
  Var     = list(do = TRUE, nboot = 50),
  cure    = TRUE
)
sol.scenario5$probs
sol.scenario5$ktau
sol.scenario5$parapar
sol.scenario5$curerate
```

scenario 6: parametric survival margin and nonparametric censoring margin
```R
set.seed(1)
sol.scenario6 <- SurvDC(
  yobs    = yobs, 
  delta   = delta, 
  tm      = quantile(yobs, c(0.25,0.50,0.75)), 
  copfam  = copfam.true,
  margins = list(survfam = "lnorm", censfam = NULL),
  Var     = list(do = TRUE, nboot = 50),
  cure    = TRUE
)
sol.scenario6$probs
sol.scenario6$ktau
sol.scenario6$parapar
sol.scenario6$curerate
```
