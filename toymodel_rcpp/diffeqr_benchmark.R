# Benchmarking of deSolve/Rcpp vs. diffeqr
# Fabienne Krauer
# 26.09.2020
# SIRS toy model with seasonal forcing and birth/death process

# 0. Housekeeping ----------------------------------------------------
library(ggplot2)
library(Rcpp); library(deSolve); library(diffeqr); library(BayesianTools); library(coda)
####
# go to folder; library(rstudioapi)
currentfile_path=rstudioapi::getActiveDocumentContext()$path
currentfile_path=paste0(unlist(strsplit(currentfile_path,"\\/"))[1:(length(unlist(strsplit(currentfile_path,"\\/")))-1)],collapse="/")
setwd(currentfile_path)

sourceCpp("modelfunc.cpp") # must source after deSolve

maxtime <- 3650.0; beta <- 0.25; gamma <- 1/5; mu <- 1/(70*365); sigma <- 1/365; k <- 5
eta <- 0.2 # scaling of sinusoidally forced FOI
phi <- 10 # phase shift of sinusoidally forced FOI

theta <- c(beta=beta, gamma=gamma, mu=mu, sigma=sigma, eta=eta, phi=phi, k=k)
init <- log(c(100000, 1, 1)) # must be log states because integration of log coordinates in model

#####
# create ode with cpp inline
cppFunction('List toymodel_inline(double time, NumericVector logstate, NumericVector par) {
  NumericVector state = exp(logstate); // exponentiate the log states
  List rparam(par);
  NumericVector dxdt(state.length());
  // params 
  double beta = as<double>(rparam["beta"]); 
  double gamma = as<double>(rparam["gamma"]);
  double mu = as<double>(rparam["mu"]);
  double sigma = as<double>(rparam["sigma"]);
  double eta = as<double>(rparam["eta"]);
  double phi = as<double>(rparam["phi"]);
  // states
  double S = state[0];   double I = state[1];  double R = state[2];  
  double N = S + I + R;
  // Seasonally forced FOI
  double beta_eff = beta * (1+eta*pow(sin(M_PI*(time-phi)/365.0),2)); // can also be another sinusoidal function
  double FOI = beta_eff*I/N;
  // change in states
  dxdt[0] = -FOI*S + mu*N - mu*S + sigma*R; 
  dxdt[1] = FOI*S - gamma*I - mu*I; 
  dxdt[2] = gamma*I - mu*R - sigma*R;
  return List::create(dxdt/state); 
}')
####

# 1. SIMULATE ----------------------------------------------------

# deSolve ----------------------------------------------------
func_desolve <- function(theta, init, times) {
  traj <- data.frame(ode(y=init, times=times, func=toymodel_inline, parms=theta, method="lsoda", verbose=F))
  return(traj)
}

traj_desolve <- func_desolve(theta, init, seq(0.0, maxtime, by=1))
traj_desolve[,2:4] <- exp(traj_desolve[,2:4])
colnames(traj_desolve) <- c("time", "S", "I", "R")
ggplot(traj_desolve) + geom_line(aes(x=time, y=I)) + 
  scale_x_continuous(breaks=seq(from=0, to=maxtime, by=365)) +
  theme(panel.grid.minor.x = element_blank())

# Julia  ----------------------------------------------------
enviro <- diffeqr::diffeq_setup() # takes a couple of mins

func_julia <- function(init, param, times) {
  state <- exp(init)
  # Named vectors not possible at the moment
  beta <- param[1];   gamma <- param[2];   mu <- param[3];   sigma <- param[4];  eta <- param[5];  phi <- param[6]
  
  S = state[1]
  I = state[2]
  R = state[3] 
  N = S + I + R
  
  beta_eff = beta * (1+eta*sin(pi*(times-phi)/365.0)^2)
  FOI = beta_eff*I/N
  
  dSdt = -FOI*S + mu*N - mu*S + sigma*R
  dIdt = FOI*S - gamma*I - mu*I 
  dRdt = gamma*I - mu*R - sigma*R

  return(c(dSdt/S, dIdt/I, dRdt/R))
}

problem <- enviro$ODEProblem(func_julia, init, c(0.0, maxtime), theta)
problemacc <- diffeqr::jitoptimize_ode(enviro,problem) # accelerate, takes a couple of mins
solution <- enviro$solve(problemacc, enviro$AutoVern7(enviro$KenCarp4()), 
                         saveat=1.0, abstol=1e-8, reltol=1e-6) # default tolerance leads to weird behaviour with this model (generation of individuals)
traj_julia <- as.data.frame(t(sapply(solution$u,identity)))
traj_julia <- exp(traj_julia)
traj_julia <- cbind(time=as.data.frame(sapply(solution$t,identity)), traj_julia)
colnames(traj_julia) <- c("time", "S", "I", "R")


# compare --------------------------------------------------
foo <- merge(traj_desolve, traj_julia, by="time", all=T)
foo$N.x <- foo$S.x+foo$I.x+foo$R.x
foo$N.y <- foo$S.y+foo$I.y+foo$R.y

summary(foo$N.x)
summary(foo$N.y)

ggplot(foo) + geom_line(aes(x=time, y=I.x), color="red") + 
  geom_line(aes(x=time, y=I.y), linetype="dashed")
ggplot(foo) + geom_line(aes(x=time, y=R.x), color="red") + 
  geom_line(aes(x=time, y=R.y), linetype="dashed")
ggplot(foo) + geom_line(aes(x=time, y=N.x), color="red") + 
  geom_line(aes(x=time, y=N.y), linetype="dashed") # some numerical imprecision, not a major problem


# 2. FIT -----------------------------------------------------------------

# generate some fake data
data <- traj_desolve[,c("time", "I")]
data$obs <- NA
for (i in 1:nrow(data)) {
  data$obs[i] <- rnbinom(1, 5, mu=round(data$I[i]))
}
ggplot(data) +  geom_point(aes(x=time, y=obs)) +
  geom_line(aes(x=time, y=I), color="red")+
  scale_x_continuous(breaks=seq(from=0, to=maxtime, by=365)) +
  theme(panel.grid.minor.x = element_blank())

lower <- c(beta=1e-2, gamma=1e-3, mu=3e-5, sigma=1e-3, eta=0, phi=0, k=1e-3)
upper <- c(beta=1, gamma=0.5, mu=9e-5, sigma=0.5, eta=1, phi=364, k=1e6)
all(lower < upper)

# ll function
loglik_desolve <- function(theta, init, times, data) { # BT= BayesianTools
  
  traj <- func_desolve(theta, init, times)
  traj <- exp(traj)
  
  datapoint <- data$obs
  modelpoint <- traj$X2
  
  # Minimize the (negative) log likelihood
  ll <- sum(dnbinom(x=datapoint,
                    mu=modelpoint,
                    size=1/theta[["k"]],
                    log=TRUE), na.rm=TRUE)
  
  return(ll)
}

# test
loglik_desolve(theta, init, seq(0.0,maxtime), data)


loglik_julia <- function(theta, data) {
  
  # Update the problem with new draw of theta
  problem_new <- enviro$remake(problemacc, p=theta)
  solution <- enviro$solve(problem_new, 
                           enviro$AutoVern7(enviro$KenCarp4()), saveat=1.0,
                           abstol=1e-8, reltol=1e-6)
  
  traj <- as.data.frame(t(sapply(solution$u,identity)))
  traj <- exp(traj)
  
  datapoint <- data$obs
  modelpoint <- traj$V2
  
  # Minimize the  (negative) log likelihood
  ll <- sum(dnbinom(x=datapoint,
                    mu=modelpoint,
                    size=1/theta[["k"]],
                    log=TRUE), na.rm=TRUE)
  
  return(ll)
}

# test
loglik_julia(theta, data)


# Fit with BayesianTools --------------------------------------------------

theta_BT <- theta
theta_BT[1:length(theta_BT)] <- NA
# theta_est <- theta_est[order(names(theta_est))] # sort to alphabetical order
# est.pars <- names(which(is.na(theta_est))) # sort to NA first
# theta <- sort(theta, na.last=FALSE) 
est.index <- c(1:length(which(is.na(theta_BT))))
niterations <- 5000 # tested both for speed with 5000. 1 million is needed for convergence (tested only with Julia: 3.5 hours)
nchains <- 1

loglikwrap_desolve <- function(par) {
  
  parX = theta_BT
  parX[est.index] = par
  
  return(loglik_desolve(data = data,
                        theta = parX,
                        times = seq(0,maxtime),
                        init = init,
                        BT = TRUE))
}  

loglikwrap_julia <- function(par) {
  
  parX = theta_BT
  parX[est.index] = par
  
  return(loglik_julia(data = data,
                      theta = parX))
}  


prior <- createUniformPrior(lower=lower, upper=upper)
mcmc.settings <- list(iterations = niterations, nrChains = nchains)

# desolve
bayesiansetup_desolve <- createBayesianSetup(prior = prior,
                                             likelihood = loglikwrap_desolve,
                                             names = names(theta_BT[est.index]),
                                             parallel = FALSE)

system.time({trace_desolve <- runMCMC(bayesianSetup = bayesiansetup_desolve, 
                                      sampler = "DEzs", settings = mcmc.settings)})

saveRDS(trace_desolve, "toy/trace_BT_desolve.rds")


# Julia
bayesiansetup_julia <- createBayesianSetup(prior = prior,
                                           likelihood = loglikwrap_julia,
                                           names = names(theta_BT[est.index]),
                                           parallel = FALSE)

system.time({trace_julia <- runMCMC(bayesianSetup = bayesiansetup_julia, 
                                    sampler = "DEzs", settings = mcmc.settings)})

saveRDS(trace_julia, "toy/trace_BT_julia.rds")

# check raw traces
plot(trace_desolve, parametersOnly = TRUE)
plot(trace_julia, parametersOnly = TRUE)

# burn-in
nburn <- 100000 # For 1 mill iterations
trace_desolve_burn <- getSample(trace_desolve, parametersOnly = TRUE, coda=TRUE, start=nburn)
plot(trace_desolve_burn, parametersOnly = TRUE, start =nburn)

trace_julia_burn <- getSample(trace_julia, parametersOnly = TRUE, coda=TRUE, start=nburn)
plot(trace_julia_burn, parametersOnly = TRUE, start =nburn)

# check convergence
gelmanDiagnostics(trace_desolve_burn, plot=TRUE, start=nburn)
gelmanDiagnostics(trace_julia_burn, plot=TRUE, start=nburn)

# param summary
summary(trace_desolve_burn)
summary(trace_julia_burn)

###############################################################
### c++ model
# compile: sudo g++ -o hello hello.cpp
# run: ./hello
# read in output of c++ script
library(readr); library(reshape2)
lorenz_output=read_csv('../odeint-v2/examples/lorenz.csv',col_names=c('t','x1','x2','x3'))
lorenz_output_tidy=melt(lorenz_output,id.vars='t')
  
# time implicit
# list.dirs('../odeint-v2',recursive=FALSE)
ggplot(data=lorenz_output,aes(x=x1,y=x2,color=t)) + geom_path() + scale_colour_gradient() + 
  theme_bw() + theme(axis.title=element_text(size=9),panel.grid=element_line(linetype="dashed",colour="black",size=0.1))
# time course
ggplot(lorenz_output_tidy,aes(x=t,y=value,color=variable)) + geom_line()
