# *SurvDC* R Package
## Survival Analysis under Dependent Censoring using Copula with Unknown Association


An R package for estimating survival function under dependent censoring.
- adopt copula-based approaches
- without the consideration of covariates
- impose fully parametric margins or let one of the margins to be nonparametric
- with or without the existence of a cured fraction
- In summary, we consider the following cases:
  - fully-parametric margins (without cure)
  - fully-parametric margins (with cure in survival time)
  - nonparametric survival margin and parametric censoring margin (without cure)
  - nonparametric survival margin and parametric censoring margin (with cure in survival time)
  - parametric survival margin and nonparametric censoring margin (without cure)
  - parametric survival margin and nonparametric censoring margin (with cure in survival time)

For details are on the way ... Thank you for your attention!


## Examples

```R


#------------------------------------------------------------------------#
# Description (Example File for Replication): 
#   Survival Analysis under Dependent Censoring
#     - based on copulas
#     - fully parametric or partially nonparametric margins
#     - without or with the inclusion of a cured fraction
#------------------------------------------------------------------------#
# Author  : Jie Ding (DLUT) and Ingrid Van Keilegom (KUL)
# E-Mail  : biostat_jieding@outlook.com (Jie Ding)
#------------------------------------------------------------------------#

#----------------------------------------------------------#
# Basic preparations before running subsequent examples ####
#----------------------------------------------------------#

## clear the environment 
rm(list=ls(all=TRUE))

## library package
library(SurvDC)

#------------------------------------------------------------------------#
# simulated data from Frank copula log-Normal margins (without cure)
#------------------------------------------------------------------------#

## generate the simulated data

# - the sample size of the generated data
n <- 2000

# - information on the used copula
copfam.true <- "frank"
ktau.true <- 0.5
coppar.true <- ktau.to.coppar(ktau=ktau.true,copfam=copfam.true)

# - parameters of the underlying log-normal marginal distributions
survpar.true <- c(2.20,1.00)
censpar.true <- c(2.20,0.25)

# - true underlying survival and censoring times
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

# - observations
yobs  <- pmin(yobs.T,yobs.C)
delta <- as.numeric(yobs.T<=yobs.C)
cat("censoring rate is",mean(1-delta))

## model the data under different scenarios

# - scenario 1: parametric survival and censoring margins
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

# - scenario 2: nonparametric survival margin and parametric censoring margin
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

# - scenario 3: parametric survival margin and nonparametric censoring margin
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

#------------------------------------------------------------------------#
# simulated data from Frank copula log-Normal margins (with cure)
#------------------------------------------------------------------------#

## generate the simulated data

# - true underlying curerate
curerate.true <- 0.2

# - true underlying survival and censoring times
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

# - observations
yobs  <- pmin(yobs.T,yobs.C)
delta <- as.numeric(yobs.T<=yobs.C)
cat("censoring rate is",mean(1-delta))

## model the data under different scenarios (with cure)

# - scenario 4: parametric survival and censoring margins
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

# - scenario 5: nonparametric survival margin and parametric censoring margin
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

# - scenario 6: parametric survival margin and nonparametric censoring margin
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
