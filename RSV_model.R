# install.packages("devtools"); library("devtools")
# install_github("SineadMorris/shinySIR")

# install.packages("deSolve"); install.packages('odin')
library(tidyverse) # library(reshape2) # library(dplyr); library(readr); library(stringr); library(ggplot2); 
# ode solving, maximum likelihood, rcpp
library(deSolve); library(bbmle); library(Rcpp) # library(GillespieSSA)
library(rstudioapi); # library(fitdistrplus)
# contact data from https://bisaloo.github.io/contactdata/index.html (Prem 2017)
library(contactdata); library(wpp2019)
# covidm
currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
# functions
source('RSV_model_functions.R')
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# set plotting theme
standard_theme=theme(panel.grid=element_line(linetype="dashed",colour="black",size=0.1),
                     plot.title=element_text(hjust=0.5,size=16),axis.text.x=element_text(size=9,angle=90),
                     axis.tlext.y=element_text(size=9),
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
  birth_term=parms[[1]];K_m=parms[[2]];contmatr_rowvector=parms[[3]];inf_vars_inds=parms[[4]];susc_vars_inds=parms[[5]]
  forcing_vector=parms[[6]]; elem_time_step=parms[[7]]
  # stack I vars
  inf_vars_stacked=do.call(cbind,lapply(inf_vars_inds, function(x){X[x]}))
  inf_vars_stacked_fullsize=t(matrix(1,1,n_inf)%*%inf_vars_stacked)
  lambda_vect=diag(forcing_vector[(t/elem_time_step)+1]*array(delta_susc))%*%contmatr_rowvector%*%inf_vars_stacked_fullsize
  infection_vect=diag(X[unlist(susc_vars_inds)])%*%lambda_vect
  F_vect=matrix(0,dim_sys,1); F_vect[c(unlist(susc_vars_inds),unlist(inf_vars_inds))]=rbind(-infection_vect,infection_vect)
  dXdt=birth_term + F_vect + K_m%*%X; list(dXdt) }

# SET PARAMETERS --------------------------------------------------------
# selected country
country_sel="Germany"
# time resolution (in days)
elem_time_step=0.25
# population data
standard_age_groups <- fun_cntr_agestr("Germany","2020",seq(0,75,5),c(seq(4,74,5),99))
# RSV age groups (population data from wpp2019)
rsv_age_groups_low=c(0,0.5,1,1.5, 2,3,4, 5,10,15, 20); rsv_age_group_sizes=c(rep(0.4,4),rep(0.9,3),rep(4,3),79)
rsv_age_groups=fun_rsv_agegroups(standard_age_groups,rsv_age_groups_low,rsv_age_group_sizes)
# population by age group
N_tot=sum(rsv_age_groups$value) # S_by_age=rsv_age_groups$value
# number of age groups, reinfections and variables (S,I,R)
n_age=nrow(rsv_age_groups); varname_list=c('S','I','R'); n_compartment=length(varname_list); n_inf=3
dim_sys=n_age*n_compartment*n_inf; n_days_year=365
# query variables: fun_sub2ind(1:3,11,"R",varname_list,n_age,n_inf)
# force of infection terms
# linear indices of the I & S variables
l_inf_susc=fun_inf_susc_index_lists(n_age,n_inf,varname_list); inf_vars_inds=l_inf_susc[[1]]; susc_vars_inds=l_inf_susc[[2]]
# CONTACT MATRIX
# random values: C_m_vals=rnorm(n_age^2,1,0.2); C_m=fun_symm_contmatr(C_m_vals,n_age)
# contact matrix from covidm ("home","work","school","other")
C_m_full=Reduce('+',fun_covidm_contactmatrix(country_sel))
# create for our age groups
C_m=fun_create_red_C_m(C_m_full,rsv_age_groups)
# bc of reinfections we need to input contact matrix repeatedly
contmatr_rowvector=t(do.call(cbind, lapply(1:nrow(C_m), function(x){diag(C_m[x,])%*%matrix(1,n_age,n_inf)})))

# build kinetic matrix --------------------------------------------------------
# waning (immunity) terms: R_i_j -> S_min(i+1,n_inf)_j
omega=1/1e2; # 1/runif(1,60,200)
# recovery terms
rho=1/6; # 1/rho=rweibull(1, shape=4.1,scale=8.3)
# KINETIC MATRIX (aging terms need to be scaled by duration of age groups!)
K_m=fun_K_m_sirs_multiage(dim_sys,n_age,n_inf,n_compartment,rho,omega,varname_list,rsv_age_groups)
# R_final=fun_sub2ind(i_inf = 1:n_inf ,j_age = n_age,"R",varname_list,n_age,n_inf); diag(K_m[R_final,R_final])=(-1/365)*10
# susceptibility
# rbeta(35.583,11.417)~0.75; B(22.829,3.171)~0.9; B(6.117,12.882)~0.32
delta_primary=c(0.5,0.25,0.125); delta_susc=sapply(1:n_age, function(x) {delta_primary/(1.05^x * rsv_age_groups$value[x])})
# delta_susc*matrix(rep(rsv_age_groups$value,times=3),3,byrow = T)
# birth term into S_1_1 (Germany 2019: 778e3 births)
birth_rate=2130; birth_term=matrix(c(birth_rate,rep(0,dim_sys-1)),dim_sys,1) # 0.01
# duration of simul
n_years=10; max_time=n_years*n_days_year; timesteps <- seq(0,max_time,by=elem_time_step)
# seasonal forcing
forcing_vector=fun_seas_forc(timesteps,peak_day,st_dev_season,basal_rate=0.02)
# ggplot(data.frame(timesteps,forcing_vector),aes(x=timesteps,y=forcing_vector)) + geom_line() + standard_theme + theme_bw()
# INITIAL CONDITIONS: SUSCEPTIBLES
df_ode_solution
initvals_sirs_model=matrix(0,dim_sys,1); initvals_sirs_model[sapply(susc_vars_inds,"[[",1)]=rsv_age_groups$value
# INITIAL INFECTION. All first infection groups: sapply(inf_vars_inds, '[[',1)
initvals_sirs_model[inf_vars_inds[[1]][1]]=1
# or introduce as init state a stationary state
initvals_sirs_model[,1]=readRDS("simul_output/statsol_40years.RDS")

###
# integrate ODE --------------------------------------------------------
# parameter inputs
params=list(birth_term,K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds,forcing_vector,elem_time_step)
# sirs_template, sirs_seasonal_forc
ptm<-proc.time(); ode_solution<-lsoda(initvals_sirs_model,timesteps,func=sirs_seasonal_forc,parms=params); proc.time() - ptm
list_simul_output=fun_process_simul_output(ode_solution,varname_list,n_age,n_inf,rsv_age_groups$value)
df_ode_solution=list_simul_output[[1]]; df_ode_solution_tidy=list_simul_output[[2]]
# check how age group totals change
ggplot(df_ode_solution_tidy %>% group_by(t,agegroup) %>% summarise(agegroup_total=sum(value)),
  aes(x=t/365,y=agegroup_total/1e6,group=agegroup)) + geom_line() + facet_wrap(~agegroup,scales="free") + theme_bw() + standard_theme
# save
# ggsave("simul_output/agegroup_totals.png",width=24,height=18,units="cm")
# Plot time course --------------------------------------------------------
# time logarithmic?
# ymaxval=max(df_ode_solution_nonzero[,2:ncol(df_ode_solution_nonzero)]); t_maxval=max(df_ode_solution_nonzero$t)
xval_lims=c(7*n_days_year,9.5*n_days_year); xval_breaks=seq(timesteps[1],timesteps[length(timesteps)],n_days_year/12);time_opt=c('_xlog')
vars_to_show=grepl('I_',df_ode_solution_tidy$name)&(df_ode_solution_tidy$t>xval_lims[1]&df_ode_solution_tidy$t<xval_lims[2])
#&  df_ode_solution_tidy$agegroup<11
ggplot(df_ode_solution_tidy[vars_to_show,], aes(x=t/365,y=value_fract,group=name,color=infection,linetype=infection)) + 
  geom_line(size=1.05) + theme_bw() + standard_theme + # guides(color=FALSE,linetype=FALSE) +
  # 1 panel: 1 agegroup, 1 # infection, 3 vartype
  facet_wrap(~agegroup,nrow=3) + # ,scales='free',labeller=label_both
  scale_x_continuous(breaks=seq(timesteps[1],timesteps[length(timesteps)],n_days_year)/365,minor_breaks=xval_breaks/365) + 
  # scale_y_log10(limits=c(0.1,1e6)) + 
  xlab('years') + ylab('') + ggtitle('# infections by age group (RSV SIRS simulation)')
  # scale_y_continuous(limits=c(0,2e2)) + # ,limits=c(0,n_years*n_days_year)
## SAVE
if (birth_rate==0) {birth_tag='_nobirth'} else {birth_tag=paste0('_birthrate',birth_rate)}; group_type='vartype_grouped'
n_age_tag=paste0(n_age,'agegroups'); fract_abs=c("_fractional","_absnum")[1]
timecourse_filename=paste0("simul_output/RSV_DE_seasforcing_",n_age_tag,fract_abs,"_overlaid.png")
ggsave(timecourse_filename,width=30,height=16,units="cm")
# scale_y_log10(limits=c(0.1,ymaxval)) + scale_size_manual(values=rev(as.numeric(unique(df_ode_solution_tidy$infection)))*0.3)
# facet_grid(compartment~infection,scales='free',labeller=label_both) + # 1 panel: 1 vartype, 1 # infection, 2 age groups
# facet_grid(compartment~agegroup+infection,scales='free',labeller=label_both) + 
# facet_wrap(~compartment+infection,ncol=3,scales='free',labeller=label_both) + # 1 panel: 1 vartype, 1 # inf, 2 age groups  
# facet_wrap(~compartment+agegroup,ncol=2,scales='free') + # 1 panel: 1 vartype, 1 age groups, 3 #s infection, 
# scale_x_log10(limits=c(1,t_maxval),breaks=scales::trans_breaks("log10", function(x) 10^x),
# labels=scales::trans_format("log10", scales::math_format(10^.x))) + annotation_logticks(sides='b') + 


