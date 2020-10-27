# install.packages("devtools"); library("devtools")
# install_github("SineadMorris/shinySIR")

# install.packages("deSolve"); install.packages('odin')
library(tidyverse); library(reshape2) # library(dplyr); library(readr); library(stringr); library(ggplot2); 
# ode solving, maximum likelihood, rcpp
library(deSolve); library(bbmle); library(Rcpp) # library(GillespieSSA)
library(rstudioapi); # library(fitdistrplus)
currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
# functions
source('RSV_model_functions.R')
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# set plotting theme
standard_theme=theme(panel.grid=element_line(linetype="dashed",colour="black",size=0.1),
                     plot.title=element_text(hjust=0.5,size=16),axis.text.x=element_text(size=9,angle=90),axis.text.y=element_text(size=9),
                     axis.title=element_text(size=14), text=element_text(family="Calibri"))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Set up SIRS model struct --------------------------------------------------------

# without reinfections the model would be: dx/dt = [-id_matr;id_matr;0]*diag(I_vect)*C_m*S_vect + K_m*x(t)
# with reinfections more complicated, abstractly denoted: 
# dx/dt = F(I_vect,C_m,S_vect) + K_m*x(t)
# I_vect: vector of infectious compartments, S_vect: vector of susceptible compartments, C_m: contact matrix
# K_m: kinetic matrix for linear terms (aging, recovery, waning of immunity)
#
# steps of constructing infection terms (lambda):
# linear list of state variables
# X_vars=matrix(abs(rnorm(n_inf*n_age*n_compartment)),n_inf*n_age*n_compartment,1) # matrix(0,n_inf*n_age*n_compartment,1)
# inf_vars_stacked=do.call(cbind,lapply(inf_vars_inds, function(x){X_vars[x]})) # do.call(cbind,inf_vars_inds)
# inf_vars_stacked_fullsize=t(matrix(1,1,n_inf)%*%inf_vars_stacked) 
# # full lambda column vector
# lambda_vect=diag(array(delta_susc))%*%contmatr_rowvector%*%inf_vars_stacked_fullsize
# # infection vector
# infection_vect=diag(X_vars[unlist(susc_vars_inds)])%*%lambda_vect
# 
# # put together RHS of ODEs (this is to test, we need to do it within ODE function)
# # dX/dt = F(delta,S,lambda) + K_m*X
# F_vect=matrix(0,dim_sys,1)
# F_vect[c(unlist(susc_vars_inds),unlist(inf_vars_inds))]=rbind(-infection_vect,infection_vect)
# rhs_odes=birth_term + F_vect + K_m%*%X_vars

# construct ODE fcn --------------------------------------------------------
sirs_template <- function(t,X,parms){
  birth_term=parms[[1]]; K_m=parms[[2]]; contmatr_rowvector=parms[[3]]; inf_vars_inds=parms[[4]]; susc_vars_inds=parms[[5]]
  # stack I vars
  inf_vars_stacked=do.call(cbind,lapply(inf_vars_inds, function(x){X[x]}))
  inf_vars_stacked_fullsize=t(matrix(1,1,n_inf)%*%inf_vars_stacked)
  lambda_vect=diag(array(delta_susc))%*%contmatr_rowvector%*%inf_vars_stacked_fullsize
  infection_vect=diag(X[unlist(susc_vars_inds)])%*%lambda_vect
  F_vect=matrix(0,dim_sys,1); F_vect[c(unlist(susc_vars_inds),unlist(inf_vars_inds))]=rbind(-infection_vect,infection_vect)
  dXdt=birth_term + F_vect + K_m%*%X; list(dXdt) }
## with seasonal forcing
sirs_seasonal_forc <- function(t,X,parms){
  birth_term=parms[[1]]; K_m=parms[[2]]; contmatr_rowvector=parms[[3]]; inf_vars_inds=parms[[4]]; susc_vars_inds=parms[[5]]
  # stack I vars
  inf_vars_stacked=do.call(cbind,lapply(inf_vars_inds, function(x){X[x]}))
  inf_vars_stacked_fullsize=t(matrix(1,1,n_inf)%*%inf_vars_stacked)
  lambda_vect=diag(array(delta_susc))%*%contmatr_rowvector%*%inf_vars_stacked_fullsize
  infection_vect=diag(X[unlist(susc_vars_inds)])%*%lambda_vect
  F_vect=matrix(0,dim_sys,1); F_vect[c(unlist(susc_vars_inds),unlist(inf_vars_inds))]=rbind(-infection_vect,infection_vect)
  dXdt=birth_term + F_vect + K_m%*%X; list(dXdt) }

# SET PARAMETERS --------------------------------------------------------
# system dimension
# number of age groups, reinfections and variables (S,I,R)
n_age=3; n_compartment=3; n_inf=3; dim_sys=n_age*n_compartment*n_inf
# force of infection terms
# linear indices of the I variables
inf_vars_inds=lapply(1:n_age, function(x_age){ sapply(1:n_inf, function(x_inf){
  fun_sub2ind(x_inf,x_age,varname='I',varname_list,n_compartment,n_age,n_inf) }) })
# linear indices of S variables
susc_vars_inds=lapply(1:n_age, function(x_age){ sapply(1:n_inf, function(x_inf){
  fun_sub2ind(x_inf,x_age,varname='S',varname_list,n_compartment,n_age,n_inf) }) })
# contact matrix
C_m_vals=rnorm(n_age^2,1,0.2); if (sum(C_m_vals<0)>0){C_m_vals[C_m_vals<0]=abs(C_m_vals[C_m_vals<0])}
C_m=matrix(C_m_vals,n_age,n_age)
contmatr_rowvector=t(do.call(cbind, lapply(1:nrow(C_m), function(x){diag(C_m[x,])%*%matrix(1,n_age,n_inf)})))
# build kinetic matrix --------------------------------------------------------
# waning (immunity) terms: R_i_j -> S_min(i+1,n_inf)_j
omega=1/1e2; # 1/runif(1,60,200)
# recovery terms
rho=1/6; # 1/rho=rweibull(1, shape=4.1,scale=8.3)
K_m=fun_K_m_sirs_multiage(dim_sys,n_age,n_inf,n_compartment,rho,omega,varname_list)
# susceptibility
# rbeta(35.583,11.417)~0.75; B(22.829,3.171)~0.9; B(6.117,12.882)~0.32
delta_primary=c(1,0.7,0.5); delta_susc=cbind(delta_primary,delta_primary/2,delta_primary/4)/500
# birth term into S_1_1
birth_rate=1; birth_term=matrix(c(birth_rate,rep(0,dim_sys-1)),dim_sys,1) # 0.01
# parameter inputs
params=list(birth_term,K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds)
# INITIAL CONDITIONS
initvals_sirs_model=matrix(0,dim_sys,1); 
# INITIAL SUSCEPTIBLES
ind_init_susc=susc_vars_inds[[1]][1]; initvals_sirs_model[ind_init_susc]=1e3
# INITIAL INFECTION. All first infection groups: sapply(inf_vars_inds, '[[',1)
initvals_sirs_model[inf_vars_inds[[1]][1]]=1
# duration of simul
n_years=12; max_time<-n_years*n_days_year; timesteps <- seq(0,max_time,by=0.1)
###
# integrate ODE --------------------------------------------------------
ptm <- proc.time(); ode_solution <- lsoda(initvals_sirs_model,timesteps,func=sirs_template,parms=params); proc.time() - ptm
df_ode_solution=ode_solution %>% as.data.frame() %>% setNames(c("t",fun_sirs_varnames(varname_list,n_age,n_inf)))
# process simul output
df_ode_solution_nonzero=df_ode_solution[,colSums(df_ode_solution)>0]
df_ode_solution_tidy=reshape2::melt(df_ode_solution[,colSums(df_ode_solution)>0],id.vars='t')
df_ode_solution_tidy[c('compartment','infection','agegroup')]=
  sapply(1:3, function(x) {sapply(strsplit(as.character(df_ode_solution_tidy$variable),'_'),'[[',x)})
df_ode_solution_tidy$compartment=factor(df_ode_solution_tidy$compartment,levels=varname_list)
ymaxval=max(df_ode_solution_nonzero[,2:ncol(df_ode_solution_nonzero)]); t_maxval=max(df_ode_solution_nonzero$t)
# inspect output
t_comp=5e2; View(round(df_ode_solution_nonzero[seq(1,t_comp,10),colSums(df_ode_solution_nonzero[1:t_comp,]>1)>1],1))

# Plot time course --------------------------------------------------------
# set themes
# standard_theme$strip.text=element_text(size=14); standard_theme$strip.text.y=element_text(size=16,color="red")
# standard_theme$legend.text=element_text(size=14); standard_theme$legend.title=element_text(size=14)
# standard_theme$legend.box.background = element_rect(colour="black")
#
# time logarithmic?
time_opt=c('_xlog') # ''
ggplot(df_ode_solution_tidy, aes(x=t,y=value,group=variable,color=compartment)) + 
  geom_line(size=1.05) + theme_bw() + standard_theme + labs(fill = "# infection") +
  facet_grid(infection~agegroup,scales='free',labeller=label_both) + # 1 panel: 1 agegroup, 1 # infection, 3 vartype
  # facet_grid(compartment~infection,scales='free',labeller=label_both) + # 1 panel: 1 vartype, 1 # infection, 2 age groups
  # facet_grid(compartment~agegroup+infection,scales='free',labeller=label_both) + 
  # facet_wrap(~compartment+infection,ncol=3,scales='free',labeller=label_both) + # 1 panel: 1 vartype, 1 # inf, 2 age groups  
  # facet_wrap(~compartment+agegroup,ncol=2,scales='free') + # 1 panel: 1 vartype, 1 age groups, 3 #s infection, 
  scale_x_log10(limits=c(1,t_maxval),breaks=scales::trans_breaks("log10", function(x) 10^x),
  labels=scales::trans_format("log10", scales::math_format(10^.x))) + annotation_logticks(sides='b') + 
  xlab('days') + ylab('') + ggtitle('SIRS simulation')
# scale_y_log10(limits=c(0.1,ymaxval)) + 
# scale_size_manual(values=rev(as.numeric(unique(df_ode_solution_tidy$infection)))*0.3) +
## SAVE
if (birth_rate==0) {birth_tag='_nobirth'} else {birth_tag=paste0('_birth_rate',birth_rate)}; group_type='vartype_grouped'
n_age_tag=paste0('_',n_age,'agegroups') # init_tag="_init_S22" # check ind_init_susc
timecourse_filename=paste0("simul_output/toymodel_timecourse_",group_type,time_opt,birth_tag,n_age_tag,".png")
ggsave(timecourse_filename,width=24,height=18,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# simple ODE fcn examples --------------------------------------------------------
#
# define ODE as a fcn

# linear chain model: -->x1-->x2-->x3 flow, x2 accumulates bc its outflow parameter is small, its value exceeds x1
parvect=c(1e1,1e1,1,1e1); 
linear_chain_ode <- function(t,x,parms) {
  p0=parms[1];p1=parms[2];p2=parms[3];p3=parms[4]; birth_vect=c(p0,0,0)
  kin_matr=matrix(c(-p1,0,0, p1,-p2,0, 0,p2,-p3),3,3,byrow=T); dxdt=birth_vect+kin_matr%*%x; list(dxdt) }
ode_solution_toymodel <- lsoda(c(0,0,0),seq(0,5e1,0.1),func=linear_chain_ode,parms=parvect)
df_ode_solution_toymodel=ode_solution_toymodel %>% as.data.frame() %>% setNames(c("t",paste0('x',1:3)))
ggplot(melt(df_ode_solution_toymodel,id.vars='t'),aes(x=t,y=value,color=variable)) + geom_line() + standard_theme + theme_bw()

# sir model
sir <- function(t,y,parms) {
  beta <- parms[1]; gamma <- parms[2]; S <- y[1]; I <- y[2]
  return(list(c(S = -beta*S*I, I = beta*S*I - gamma*I))) }

sir_birth_death <- function(t,y,parms){
  mu=parms[1]; beta=parms[2]; gamma=parms[3]; S=y[1]; I=y[2]
  dxdt = c(-beta*S*I-mu*S+mu,beta*S*I - (gamma+mu)*I)
  list(dxdt) }
# t - time; y - current state vector of the ODE at time t; parms - Parameter vector used by the ODE system
# Returns list with 1 component being a vector of length two containing: dS(t)/dt and dI(t)/dt
####

# simulate simple ODE examples (fractional) -------------------------------
gamma=1/3; mu=1/60; beta=1.05; sigma=beta/(gamma+mu)
initvals_S_I=c(0.99,0.01); max_time<-150; timesteps <- seq(0, max_time, by=0.1)
# ode_solution <- rk4(y=initvals_S_I,times=timesteps,func=sir_birth_death,parms=c(mu,beta,gamma))
params=c(mu,beta,gamma)
ode_solution <- lsoda(initvals_S_I,timesteps,func=sir_birth_death,parms=params)
df_ode_solution=ode_solution %>% as.data.frame() %>% setNames(c("t", "S", "I"))

# time vs state vars or statevars-statevars
ggplot(df_ode_solution, aes(x=S,y=I,color=t)) + # df_ode_solution
  # geom_line(aes(y=S),color="darkred") + geom_line(aes(y=I), color="steelblue", linetype="twodash") +
  geom_path() + scale_colour_gradient() + # scale_colour_gradientn(colours=terrain.colors(length(timesteps))) + 
  ggtitle('SIR+birth model') + theme_bw() + theme(axis.title=element_text(size=9),
        panel.grid=element_line(linetype="dashed",colour="black",size=0.1)) + 
  # xlab('time') + ylab('S,I') + xlim(0,max(timesteps)) + ylim(0,1)
  xlab('S') + ylab(('I'))

####
# calculate sum of square error ----------------------------------------------------
data=df_ode_solution # 
data[,c('S','I')] = data[,c('S','I')]*(1+matrix(rnorm(2*nrow(df_ode_solution),mean=0,sd=0.03),ncol=2))
# cbind(df_ode_solution$t,matrix(rnorm(2*nrow(df_ode_solution),mean=0,sd=0.01),ncol=2))
# sum of squared error
sir_sse = function(logparams,y){
  df_simul <- lsoda(initvals_S_I,timesteps,func=sir_birth_death,parms=exp(logparams)) %>% 
    as.data.frame() %>% setNames(c("t", "S", "I"))
  I_sqerr = (data$I - df_simul$I)^2; S_sqerr = (data$S - df_simul$S)^2
  sse=sum(I_sqerr) + sum(S_sqerr); sse }

sir_sse(log(params),data)

# Negative log likelihood ----------------------------------------------------
# assuming poisson distrib
N=1e6 # popul size
# calculate negative loglikelihood
sir_nll = function(logparams,y,popul_size){
  df_simul <- lsoda(initvals_S_I,timesteps,func=sir_birth_death,parms=exp(logparams)) %>% 
    as.data.frame() %>% setNames(c("t", "S", "I"))
  logdensities_I=dpois(x=round(popul_size*data$I),lambda=round(popul_size*df_simul$I),log=TRUE)
  logdensities_S=dpois(x=round(popul_size*data$S),lambda=round(popul_size*df_simul$S),log=TRUE)
  nll=-sum(logdensities_I) + (-sum(logdensities_S))
  nll
}

# R squared (least squares)
SS_tot=sum((df_ode_solution[,c('S','I')]-colMeans(df_ode_solution[,c('S','I')]))^2)
SS_res=sir_sse(log(params),data) # sum((data[,c('S','I')]-df_ode_solution[,c('S','I')])^2)
R_sq = 1-SS_res/SS_tot

# NegLogLikelihood
par_initguess=c(0.01,1,0.1)
sir_nll(logparams=log(par_initguess),data,N)/N
# with true value
sir_nll(logparams=log(params),data,N)/N

####
# least-square fitting with Nelder-Mead ----------------------------------------------------
# true values: params=c(1/60, 1.05,1/3)
par_initguess=c(0.01,1,0.1); parnames=c('mu','beta','gamma')
# Initialize parameters object
optim_proc=capture.output(optim(log(par_initguess),fn=sir_sse,method='Nelder-Mead',y=data,control=c(trace=2)))
####
# extract convergence process
fcn_extract_optimresults <- function(optim_proc,parnames){
  optim_outputs=as.data.frame(
  t(trimws(gsub("\\s+", " ",str_replace_all(optim_proc[which(grepl('\\$',optim_proc))+1],'\\[1\\]','')))),stringsAsFactors=FALSE)
  colnames(optim_outputs)=str_replace_all(optim_proc[grepl('\\$',optim_proc)],'\\$','')
  optim_outputs=optim_outputs %>% separate(par,into=paste0('log(',parnames,')'),sep=' ')
  optim_outputs[,1:length(parnames)]=as.numeric(optim_outputs[,1:length(parnames)])
  optim_outputs
}
####
optim_outputs=fcn_extract_optimresults(optim_proc,parnames)
# convergence process
conv_proc=gsub("\\s+", " ",optim_proc[(which(grepl('Stepsize|Exit',optim_proc))[1]+1):(which(grepl('Stepsize|Exit',optim_proc))[2]-1)])
conv_proc=as.data.frame(conv_proc) %>% separate('conv_proc',c('step_type','n_step','sse1','sse2'),sep=' ')
# post-optim error
sir_sse(as.numeric(optim_outputs[,1:length(parnames)]),data)
# initguess error
sir_sse(log(par_initguess),data)

# with loglikelihood
optim_proc_max_llh=optim(log(par_initguess),fn=sir_nll,method='Nelder-Mead',y=data,popul_size=N,control=c(trace=2))
# yes!
# fractional error
abs(exp(optim_proc_max_llh$par)-params)/params

###
# Amadillo sample code ----------------------------------------------------
# (code from Stefan: https://github.com/StefanFlasche/simpleCorePneumoModel)
library("RcppArmadillo"); library("coda"); library("rootSolve") # install.packages("rootSolve")
sourceCpp("toymodel_rcpp/functions.cpp")

df_kilifi=data.frame(Setting="kilifi",
                     Age_groups=c("<1y","1-5y","6-14y","15-20y","21-49y","50+"),
                     Age_group_upper = c(1,6,15,21,50,62),
                     Population=c(9617,45170,68547,33289,72143,24214),
                     VT.prev=c(0.41,0.338,.146,.142,.072,.038),#PCV10
                     NVT.prev=c(.46,.44,.39,.25,.22,.20),
                     N.carr=c(39,127,82,56,97,104),
                     VT.clear=c(.062,.12,.34,.34,.34,.34)/7, #weekly conversion to daily
                     NVT.clear=c(.086,.15,.34,.34,.34,.34)/7)

# set parameters
df=df_kilifi
competition=0.1
no.agegps=dim(df)[1]
agegp.l=df$Age_group_upper - c(0,df$Age_group_upper[-no.agegps])
state_prop=c(agegp.l*100-2,rep(1,no.agegps),rep(1,no.agegps),rep(0,no.agegps))
res = runsteady(y=state_prop, fun=SIS_ode_cpp, parms=param.list, times=c(0,1e5))
