library(rstudioapi); currentfile_path=rstudioapi::getActiveDocumentContext()$path
currentfile_path=paste0(unlist(strsplit(currentfile_path,"\\/"))[1:(length(unlist(strsplit(currentfile_path,"\\/")))-1)],collapse="/")
setwd(currentfile_path)

######
library(BayesianTools) # The mcmc solver, needs a few packages such as coda
library(Rcpp)       # For c++ integration
library(RcppEigen)  # For c++ integration
library(socialmixr) # For POLYMOD contact data

age_lim <- c(0,2)
#  MUST make sure the contact matrix and the observations data is stratified by the right age groups

# mu -> Number of birth per day in England + Wales (2017 ONS)
# pop->  Population of England (2017 ONS)
dem = list( mu = 1860.564,   pop = 58744600 )

#  POLYMOD data for contact matrices
data(polymod) # Get the polymod
poly <- contact_matrix(polymod, countries = "United Kingdom", age.limits = age_lim, symmetric = TRUE)

# Get the observation data (cleaned). 
# Rows are the number of positive RSV samples in a week, columns are age group. 
# Number or rows is 364 (7 years of data from 2010 -> 2017)
obsdata <- read.table("../rdms/RSV_RSDM_pos_trim.txt")

# List of data
data = list(  obsdata=obsdata,   poly=poly,  dem=dem,  agegroup=age_lim )

# call cpp code

sourceCpp("./src/RcppModule.cpp")
classEvaluateLogLikelihood <- new(EvaluateLogLikelihood, dem$mu, dem$pop, age_lim)
classEvaluateLogLikelihood$contactMatrix <- poly$matrix
classEvaluateLogLikelihood$observedData <- as.matrix(obsdata)

# list of parameters to be fitted
parnames <- c("ga0", "om", "a", "b", "phi", "psi", "I1", "I2", "I3", "d1", "ep1", "ep2", "alpha")

# Define the log likelihood
llikelihood <- function(params){
  ll <-classEvaluateLogLikelihood$evaluateLogLikelihoodCpp(params)
  if (ll == -Inf) # This is to stop infinite log likelihoods gettings stuck
    ll = -1000000
  
  return(ll)
}

# Define the log priors (density and samplers)
density = function(params){
  p1 = 0
  for (i in 1:13)
  {
    p1 <- p1 + dlogis(params[i], location = 0, scale = 1, log = TRUE)
  }
  return(p1)
}

# Define the sampler for the first time step
sampler = function(n=1){
  s1 <- rlogis(13, location = 0, scale = 1)
  return(s1)
}

# Create the prior
prior <- createPrior(density = density, sampler = sampler, lower = rep(-1000, 13), upper = rep(1000, 13)) 
# lower are the lower bounds of the prior and upper are the upper bounds of the prior

# Create Bayesian Setup
bayesianSetup <- createBayesianSetup(llikelihood, prior, names = parnames)

# Differential Evolution with Snooker update, can choose different samplers, see
settingsDezs <- list(iterations = 1e5, nrChains = 3, thin = 10)
outDezs <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", 
                   settings = settingsDezs) # Standard function in Bayesian tools

summary(outDezs)
plot(outDezs)
correlationPlot(outDezs)
marginalPlot(outDezs, prior = TRUE)

##### binomial distribution

# eg: x=8, z=10, epsilon=0.5
x=12; z=10; epsilon=0.5
choose(z,x)*(epsilon^x)*(1-epsilon)^(z-x)
