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
                     plot.title=element_text(hjust=0.5,size=16),axis.text.x=element_text(size=9,angle=90),
                     axis.text.y=element_text(size=9),
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
# model without seas forcing in fcns file
## with seasonal forcing
sirs_seasonal_forc <- function(t,X,parms){ 
  birth_term=parms[[1]]; K_m=parms[[2]]; contmatr_rowvector=parms[[3]]; inf_vars_inds=parms[[4]]; susc_vars_inds=parms[[5]]
  forcing_vector=parms[[6]]
  # stack I vars
  inf_vars_stacked=do.call(cbind,lapply(inf_vars_inds, function(x){X[x]}))
  inf_vars_stacked_fullsize=t(matrix(1,1,n_inf)%*%inf_vars_stacked)
  lambda_vect=diag(forcing_vector[(t+0.1)*10]*array(delta_susc))%*%contmatr_rowvector%*%inf_vars_stacked_fullsize
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
# duration of simul
n_years=10; max_time=n_years*n_days_year; timesteps <- seq(0,max_time,by=0.1)
# seasonal forcing
peak_day=60; st_dev_season=27 # 27
dist_from_peak=apply(data.frame(abs(timesteps %% 365 - peak_day),n_days_year-(timesteps %% 365)+peak_day),1,min)
basal_rate=0.1; forcing_vector=basal_rate + exp(-0.5*(dist_from_peak/st_dev_season)^2)
# ggplot(data.frame(timesteps,forcing_vector),aes(x=timesteps,y=forcing_vector)) + geom_line() + standard_theme + theme_bw()
# parameter inputs
params=list(birth_term,K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds,forcing_vector)
# INITIAL CONDITIONS
initvals_sirs_model=matrix(0,dim_sys,1); 
# INITIAL SUSCEPTIBLES
ind_init_susc=susc_vars_inds[[1]][1]; initvals_sirs_model[ind_init_susc]=1e3
# INITIAL INFECTION. All first infection groups: sapply(inf_vars_inds, '[[',1)
initvals_sirs_model[inf_vars_inds[[1]][1]]=1
###
# integrate ODE --------------------------------------------------------
# sirs_template, sirs_seasonal_forc
ptm<-proc.time(); ode_solution<-lsoda(initvals_sirs_model,timesteps,func=sirs_seasonal_forc,parms=params); proc.time() - ptm
df_ode_solution=ode_solution %>% as.data.frame() %>% setNames(c("t",fun_sirs_varnames(varname_list,n_age,n_inf)))
# process simul output
df_ode_solution_nonzero=df_ode_solution[,colSums(df_ode_solution)>0]
df_ode_solution_tidy=reshape2::melt(df_ode_solution[,colSums(df_ode_solution)>0],id.vars='t')
df_ode_solution_tidy[c('compartment','infection','agegroup')]=
  sapply(1:3, function(x) {sapply(strsplit(as.character(df_ode_solution_tidy$variable),'_'),'[[',x)})
df_ode_solution_tidy$compartment=factor(df_ode_solution_tidy$compartment,levels=varname_list)
ymaxval=max(df_ode_solution_nonzero[,2:ncol(df_ode_solution_nonzero)]); t_maxval=max(df_ode_solution_nonzero$t)
# inspect output
# t_comp=5e2; View(round(df_ode_solution_nonzero[seq(1,t_comp,10),colSums(df_ode_solution_nonzero[1:t_comp,]>1)>1],1))

# Plot time course --------------------------------------------------------
# set themes
# standard_theme$strip.text=element_text(size=14); standard_theme$strip.text.y=element_text(size=16,color="red")
# standard_theme$legend.text=element_text(size=14); standard_theme$legend.title=element_text(size=14)
# standard_theme$legend.box.background = element_rect(colour="black")
#
# time logarithmic?
time_opt=c('_xlog'); vars_to_show=!grepl('S',df_ode_solution_tidy$variable) # 1:nrow(df_ode_solution_tidy) 
ggplot(df_ode_solution_tidy[vars_to_show,], aes(x=t,y=value,group=variable,color=compartment,linetype=compartment)) + 
  geom_line(size=1.05) + theme_bw() + standard_theme + labs(fill = "# infection") +
  facet_grid(infection~agegroup,scales='free',labeller=label_both) + # 1 panel: 1 agegroup, 1 # infection, 3 vartype
  # facet_grid(compartment~infection,scales='free',labeller=label_both) + # 1 panel: 1 vartype, 1 # infection, 2 age groups
  # facet_grid(compartment~agegroup+infection,scales='free',labeller=label_both) + 
  # facet_wrap(~compartment+infection,ncol=3,scales='free',labeller=label_both) + # 1 panel: 1 vartype, 1 # inf, 2 age groups  
  # facet_wrap(~compartment+agegroup,ncol=2,scales='free') + # 1 panel: 1 vartype, 1 age groups, 3 #s infection, 
  # scale_x_log10(limits=c(1,t_maxval),breaks=scales::trans_breaks("log10", function(x) 10^x),
  # labels=scales::trans_format("log10", scales::math_format(10^.x))) + annotation_logticks(sides='b') + 
  scale_x_continuous(breaks=seq(timesteps[1],timesteps[length(timesteps)],n_days_year),
      minor_breaks=seq(timesteps[1],timesteps[length(timesteps)],n_days_year/12),
      limits=c(7*n_days_year,n_years*n_days_year)) + scale_y_continuous(limits=c(0,2e2)) +
  xlab('days') + ylab('') + ggtitle('SIRS simulation')
# scale_y_log10(limits=c(0.1,ymaxval)) + scale_size_manual(values=rev(as.numeric(unique(df_ode_solution_tidy$infection)))*0.3)
## SAVE
if (birth_rate==0) {birth_tag='_nobirth'} else {birth_tag=paste0('_birthrate',birth_rate)}; group_type='vartype_grouped'
n_age_tag=paste0('_','agegroups',n_age) # init_tag="_init_S22" # check ind_init_susc
timecourse_filename=paste0("simul_output/toymodel_timecourse_seasforcing_",group_type,time_opt,birth_tag,n_age_tag,".png")
ggsave(timecourse_filename,width=24,height=18,units="cm")


