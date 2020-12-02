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
rm(list=ls()); source('RSV_model_functions.R')
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# set plotting theme
standard_theme=theme(panel.grid=element_line(linetype="dashed",colour="black",size=0.1),
                     plot.title=element_text(hjust=0.5,size=16),axis.text.x=element_text(size=9,angle=90),
                     axis.tlext.y=element_text(size=9),
                     axis.title=element_text(size=14), text=element_text(family="Calibri"))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Set up SIRS ODE model --------------------------------------------------------
# see description in "RSV_model_functions.R"
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
# WANING (immunity) terms: R_i_j -> S_min(i+1,n_inf)_j
omega=1/1.25e2; # 1/runif(1,60,200)
# RECOVERY
rho=1/7; # 1/rho=rweibull(1, shape=4.1,scale=8.3)
# KINETIC MATRIX (aging terms need to be scaled by duration of age groups!)
K_m=fun_K_m_sirs_multiage(dim_sys,n_age,n_inf,n_compartment,rho,omega,varname_list,rsv_age_groups)
# SUSCEPTIBILITY # rbeta(35.583,11.417)~0.75; B(22.829,3.171)~0.9; B(6.117,12.882)~0.32
delta_primary=c(0.5,0.25,0.125); delta_susc=sapply(1:n_age, function(x) {delta_primary/(1.05^x * rsv_age_groups$value[x])})
# delta_susc*matrix(rep(rsv_age_groups$value,times=3),3,byrow = T)
# BIRTH RATE into S_1_1 (Germany 2019: 778e3 births)
birth_rate=2130; birth_term=matrix(c(birth_rate,rep(0,dim_sys-1)),dim_sys,1) # 0.01
####
# DURATION of SIMULATION
n_years=10.5; max_time=n_years*n_days_year; timesteps <- seq(0,max_time,by=elem_time_step)
# seasonal forcing
start_week=49; forcing_vector=fun_seas_forc(timesteps,peak_day=start_week*7,st_dev_season=27,basal_rate=0.02) # matplot(timesteps,forcing_vector,type="l")
# shutdown season
shutdown_start_week=8*52+(start_week-4); shutdown_stop_week=shutdown_start_week+12; shutdown_scale=0.1; n_prec=0.001
shutdown_list=fun_shutdown_seasforc(shutdown_start_week,shutdown_stop_week,shutdown_scale,forcing_vector,
                                      elem_time_step,basal_rate,n_prec)
forcing_vector=shutdown_list[[1]]; shutdown_limits=shutdown_list[[2]]
# INITIAL CONDITIONS
# introduce stationary state as init state?
initvals_sirs_model=matrix(0,dim_sys,1); stationary_init=TRUE
if (stationary_init){# statsol=readRDS("simul_output/statsol_100years.RDS") # as.numeric(statsol[2:length(statsol)])
  initvals_sirs_model[,1]=round(as.numeric(df_ode_solution[nrow(df_ode_solution),2:ncol(df_ode_solution)]))} else {
# entire popul into susceptibles
  initvals_sirs_model[sapply(susc_vars_inds,"[[",1)]=rsv_age_groups$value }
# INITIAL INFECTION 
initvals_sirs_model[inf_vars_inds[[1]][1]]=1 # all first infection groups: sapply(inf_vars_inds, '[[',1)

### integrate ODE --------------------------------------------------------
# deSolve input
params=list(birth_term,K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds,forcing_vector,elem_time_step)
# sirs_template, sirs_seasonal_forc
ptm<-proc.time(); ode_solution<-lsoda(initvals_sirs_model,timesteps,func=sirs_seasonal_forc,parms=params); proc.time()-ptm
list_simul_output=fun_process_simul_output(ode_solution,varname_list,n_age,n_inf,rsv_age_groups)
df_ode_solution=list_simul_output[[1]]; df_ode_solution_tidy=list_simul_output[[2]]
####
# initial vs final popul sizes: fun_agegroup_init_final_pop(df_ode_solution_tidy)
# PLOT how age group totals change
scale_val=c('free_y','fixed')[2]; ggplot(df_ode_solution_tidy %>% group_by(t_years,agegroup_name) %>% summarise(agegroup_total=sum(value)),
  aes(x=t_years,y=agegroup_total,group=agegroup_name)) + geom_line() + facet_wrap(~agegroup_name,scales=scale_val) +
  scale_y_log10() + theme_bw() + standard_theme + xlab("year") + ylab("million popul")
# save
ggsave(paste("simul_output/agegroup_totals_",scale_val,"yscale.png",sep=""),width=28,height=16,units="cm")
####
# Plot time course --------------------------------------------------------
xval_lims=c(max_time/365-3.75,max_time/365-0.25); xval_breaks=seq(0,max_time/365,by=1/4)
vars_to_show=grepl('I_',df_ode_solution_tidy$name)&(df_ode_solution_tidy$t_years>xval_lims[1]&df_ode_solution_tidy$t_years<xval_lims[2])
# y axis fixed? | facet by AGE only or age+INFECTION? | abs values or fraction?
# all permutations
all_perms=permutations(n=2,r=3,repeats.allowed=T)
for (k in 1:nrow(all_perms)) {
scale_val=c('fixed','free_y')[all_perms[k,1]];facet2tag=c('','infection')[all_perms[k,2]];value_type=c("value","value_fract")[all_perms[k,3]]
if (grepl("fract",value_type)){y_axis_tag='fraction'} else {y_axis_tag="# cases"}; if (nchar(facet2tag)){nrow_val=6} else{nrow_val=3}
# PLOT
ggplot(df_ode_solution_tidy[vars_to_show,],aes_string(x="t_years",y=value_type,group="name",color="infection")) + #,linetype="infection"
  geom_line(size=1.05) + theme_bw() + standard_theme + theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=6)) + 
  facet_wrap(as.formula(paste('~',gsub("^\\+","",paste(facet2tag,'+agegroup_name',sep='')),sep='')),nrow=nrow_val,scales=scale_val) + 
  scale_x_continuous(breaks=xval_breaks,minor_breaks=seq(0,max_time/365,by=1/12)) + # scale_y_log10(limits=c(0.1,1e6)) + 
  # shutdown
  geom_rect(aes(xmin=shutdown_limits[1]/365,xmax=shutdown_limits[2]/365,ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
  xlab('years') + ylab(y_axis_tag) + ggtitle('RSV infections by age group (SIRS simulation)')
## SAVE
timecourse_filename=fun_create_filename("simul_output",facet2tag,value_type,n_age,gsub("_y","",scale_val),"png")
ggsave(timecourse_filename,width=32,height=22,units="cm")}

### RSV data  --------------------------------------------------------
resp_virus_data_uk=read_csv("data/Respiratory viral defections by any method UK Ages.csv")
resp_virus_data_uk_tidy=resp_virus_data_uk[,!colnames(resp_virus_data_uk) %in% "X1"] %>% pivot_longer(!c("Year","startweek","Age"))
resp_virus_data_uk_tidy$Age=factor(resp_virus_data_uk_tidy$Age,levels=unique(resp_virus_data_uk_tidy$Age))
truthvals_rsv=resp_virus_data_uk_tidy$name %in% "RSV" & (!resp_virus_data_uk_tidy$Year %in% c(2014,2020))
# means across years
averages_years=data.frame(resp_virus_data_uk_tidy[truthvals_rsv,] %>% group_by(startweek,Age) %>% 
                            summarise(value=mean(value)),Year='mean',name='RSV')
if (!any(resp_virus_data_uk_tidy$Year %in% "mean")){
resp_virus_data_uk_tidy=rbind(resp_virus_data_uk_tidy,averages_years[,colnames(resp_virus_data_uk_tidy)])}
truthvals_rsv=resp_virus_data_uk_tidy$name %in% "RSV" & (!resp_virus_data_uk_tidy$Year %in% c(2014,2020))
resp_virus_data_uk_tidy[,"type"]="indiv year"; resp_virus_data_uk_tidy[resp_virus_data_uk_tidy$Year %in% "mean","type"]="5-year average"
resp_virus_data_uk_tidy[,"width"]=1.01; resp_virus_data_uk_tidy[resp_virus_data_uk_tidy$Year %in% "mean","width"]=1.015
# plot
ggplot(resp_virus_data_uk_tidy[truthvals_rsv,],aes(x=startweek,y=value,group=Year,color=Year,linetype=factor(type),size=width)) + 
  geom_line() + geom_point(aes(shape=factor(type))) + facet_wrap(~Age,scales="free") + scale_linetype_manual(values=c("solid","dashed")) + 
  scale_size(range=c(1,1.5), guide=FALSE) + labs(shape="data type",linetype="data type") + theme_bw() + standard_theme + ylab("")
# SAVE
ggsave("simul_output/uk_rsv_data.png",width=32,height=16,units="cm")
# resp_virus_data_uk_tidy[truthvals_rsv,] %>% group_by(Year,Age) %>% filter(value==max(value))