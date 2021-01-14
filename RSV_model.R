# install.packages("devtools"); library("devtools")
# install_github("SineadMorris/shinySIR")

# ode solving, maximum likelihood, rcpp
# contact data from https://bisaloo.github.io/contactdata/index.html (Prem 2017)
# library(contactdata); library(fitdistrplus);  library(bbmle); library(Rcpp); library(GillespieSSA)
lapply(c("tidyverse","deSolve","gtools","rstudioapi","wpp2019","plotly","Rcpp"), library, character.only=TRUE)
# functions
rm(list=ls()); currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
source('RSV_model_functions.R')
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
  dXdt=birth_term + F_vect + K_m %*% X; list(dXdt) }

# SET PARAMETERS --------------------------------------------------------
# selected country
country_sel="Germany"
# time resolution (in days)
elem_time_step=0.25
# population data
standard_age_groups <- fun_cntr_agestr("Germany","2015",seq(0,75,5),c(seq(4,74,5),99))
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
cm_path="~/Desktop/research/models/epid_models/covid_model/lmic_model/covidm/"
list_contmatrs=fun_covidm_contactmatrix(country_sel="Germany",currentdir_path=currentdir_path,cm_path=cm_path)
C_m_full=Reduce('+',list_contmatrs)
# make matrix symmetric
C_m_fullrecipr=fun_recipr_contmatr(C_m_full,age_group_sizes=standard_age_groups$values)
# create for our age groups
C_m=fun_create_red_C_m(C_m_fullrecipr,rsv_age_groups,
                orig_age_groups_duration=standard_age_groups$duration,orig_age_groups_sizes=standard_age_groups$values)
# make it reciprocal for the larger group
C_m=fun_recipr_contmatr(C_m,age_group_sizes=rsv_age_groups$value)
# bc of reinfections we need to input contact matrix repeatedly
contmatr_rowvector=t(do.call(cbind, lapply(1:nrow(C_m), function(x){diag(C_m[x,]) %*% matrix(1,n_age,n_inf)})))

# build kinetic matrix --------------------------------------------------------
# WANING (immunity) terms: R_i_j -> S_min(i+1,n_inf)_j
omega=1/1.25e2; # 1/runif(1,60,200)
# RECOVERY
rho=1/7; # 1/rho=rweibull(1, shape=4.1,scale=8.3)
# KINETIC MATRIX (aging terms need to be scaled by duration of age groups!)
K_m=fun_K_m_sirs_multiage(dim_sys,n_age,n_inf,n_compartment,rho,omega,varname_list,rsv_age_groups)
# SUSCEPTIBILITY # rbeta(35.583,11.417)~0.75; B(22.829,3.171)~0.9; B(6.117,12.882)~0.32
# NORMALIZE by age group sizes (this is for the infection terms ~ delta*(I1+I2+...+In)*S_i/N_i)
# agedep_fact determines strength of age dependence, if agedep_fact>1, decreasing susceptibility with age
agedep_fact=1; delta_primary=c(0.5,0.25,0.125)
delta_susc=sapply(1:n_age, function(x) {delta_primary/((agedep_fact^x)*rsv_age_groups$value[x])})
delta_susc_prop=delta_susc*matrix(rep(rsv_age_groups$value,3),nrow=3,byrow=T)
### PLOT susceptibility ~ f(age,exposure)
suscept_agedep=fcn_suscept_agedeptable(rsv_age_groups,delta_susc,n_inf)
ggplot(suscept_agedep,aes(x=agegroup,y=value,group=name,color=name)) + geom_line(size=2) + theme_bw() + standard_theme +
  xlab("age group (years)") + ylab("susceptibility") + ggtitle('susceptibility ~ f(age,#infection')
# if (length(unique(round(suscept_agedep$value,6)))==3){dep_tag="expos_dep"}else{dep_tag="age_expos_dep"}
# ggsave(paste0("simul_output/suscept_",dep_tag,".png"),width=32,height=22,units="cm")
#
# BIRTH RATE into S_1_1 (Germany 2019: 778e3 births)
birth_rate=2130; birth_term=matrix(c(birth_rate,rep(0,dim_sys-1)),dim_sys,1) # 0.01
####
# DURATION of SIMULATION
n_years=10.5; max_time=n_years*n_days_year; timesteps <- seq(0,max_time,by=elem_time_step)
# seasonal forcing
start_week=49; forcing_vector=fun_seas_forc(timesteps,peak_day=start_week*7,st_dev_season=27,basal_rate=0.02) 
# shutdown season
shutdown_start_week=7*52+(start_week-4); shutdown_stop_week=shutdown_start_week+12; shutdown_scale=0.1; n_prec=0.001
shutdown_list=fun_shutdown_seasforc(shutdown_start_week,shutdown_stop_week,shutdown_scale,forcing_vector,
                                      elem_time_step,basal_rate=0.02,n_prec)
forcing_vector=shutdown_list[[1]]; shutdown_limits=shutdown_list[[2]] #matplot(timesteps/365,forcing_vector,type="l"); axis(1,at=0:10)
# INITIAL CONDITIONS
# introduce stationary state as init state?
initvals_sirs_model=matrix(0,dim_sys,1); stationary_init=FALSE
if (stationary_init){
  # statsol=readRDS("simul_output/statsol_100years.RDS") # as.numeric(statsol[2:length(statsol)])
  initvals_sirs_model[,1]=round(as.numeric(df_ode_solution[nrow(df_ode_solution),2:ncol(df_ode_solution)])) } else {
# at t=0 entire popul into susceptibles
  initvals_sirs_model[sapply(susc_vars_inds,"[[",1)]=rsv_age_groups$value }
# INITIAL INFECTION 
initvals_sirs_model[inf_vars_inds[[1]][1]]=10 # all first infection groups: sapply(inf_vars_inds, '[[',1)

### integrate ODE --------------------------------------------------------
# deSolve input
params=list(birth_term,K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds,forcing_vector,elem_time_step)
# sirs_template, sirs_seasonal_forc
ptm<-proc.time(); ode_solution<-lsoda(initvals_sirs_model,timesteps,func=sirs_seasonal_forc,parms=params); proc.time()-ptm
list_simul_output=fun_process_simul_output(ode_solution,varname_list,n_age,n_inf,rsv_age_groups)
df_ode_solution=list_simul_output[[1]]; df_ode_solution_tidy=list_simul_output[[2]]; rm(list_simul_output)
####
# PLOT how age group totals change (initial vs final popul sizes: fun_agegroup_init_final_pop(df_ode_solution_tidy))
scale_val=c('free_y','fixed')[2]; fcn_plotagegroup_totals(df_ode_solution_tidy,scale_val)
# ggsave(paste("simul_output/agegroup_totals_",scale_val,"yscale.png",sep=""),width=28,height=16,units="cm")
####
# Plot time course --------------------------------------------------------
xval_lims=c(floor(shutdown_limits/365)[1]-0.15,ceiling(shutdown_limits/365)[2]+0.15); xval_breaks=seq(0,max_time/365,by=1/4)
agegr_lim=7
plot_df_ode_sol=subset(df_ode_solution_tidy,grepl('I',name) & agegroup<=agegr_lim & (t_years>xval_lims[1] & t_years<xval_lims[2]))
for (k in 1:(2^3)) {
# tags: y-axis fixed/free | facet by AGE/(age&INFECTION) | abs values/fraction (all_perms=permutations(n=2,r=3,repeats.allowed=T))
  g(scale_val,facet2tag,value_type,y_axis_tag,nrow_val,height_div,facet_formula) %=% fun_tcourse_plottags(k,nval=2,rval=3)
# PLOT
ggplot(plot_df_ode_sol,aes_string(x="t_years",y=value_type,group="name",color="infection")) +
  geom_line(size=1.05) + theme_bw() + standard_theme + theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=6)) + 
  facet_wrap(as.formula(facet_formula),ncol=round(agegr_lim/as.numeric(height_div)),scales=scale_val) + 
  scale_x_continuous(breaks=xval_breaks,minor_breaks=seq(0,max_time/365,by=1/12)) + ## shutdown:
  geom_rect(aes(xmin=shutdown_limits[1]/365,xmax=shutdown_limits[2]/365,ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
  xlab('years') + ylab(y_axis_tag) + ggtitle('RSV infections by age group (SIRS simulation)') # + xlim(c(0,5))
## SAVE
if (length(unique(round(array(delta_susc_prop),5)))==length(delta_primary)){foldername="suscept_noagedep"}else{foldername="suscept_agedep"}
timecourse_filename=fun_create_filename(paste0("simul_output/",foldername),facet2tag,value_type,n_age,gsub("_y","",scale_val),"png")
ggsave(timecourse_filename,width=30,height=18,units="cm")
}

### map infections to symptomatic cases -----------------------------------
prop_symptom=1-c(mean(rbeta(1e3,shape1=3,shape2=30)),mean(rbeta(1e3,shape1=9,shape2=43)),mean(rbeta(1e3,shape1=38,shape2=35)),
                mean(rbeta(1e3,shape1=36,shape2=11))) # prop_symptom_age=cbind(c(0,1,5,15),c(0.9,4.9,14.9,99),prop_symptom_age); 
list_symptom_agegroups=list(1:2,3:7,8:9,10:11); 
# EXPOSURE DEPENDENT or not?
expos_agedep=1 # expos_agedep=0: no dependence on exposure. if 0<expos_agedep<1:sublinear dependence, expos_agedep=1: linear depend
df_symptom_prop=fun_propsymptom_table(list_symptom_agegroups,expos_dep_val=expos_agedep,rsv_age_groups$agegroup_name,n_inf=3)
# plot clinical fraction as fcn of age and #infection
# ggplot(df_symptom_prop,aes(x=agegroup_name,y=sympt_value,color=factor(n_inf),linetype=factor(n_inf),group=n_inf)) + geom_line(size=1.25) + 
#   theme_bw() + standard_theme + labs(color="# infection") + ylab("% symptomatic")
# ggsave("simul_output/severity_age_noexp_depend.png",width=22,height=16,units="cm")
# calculate symptom cases
df_ode_solution_tidy_cases=fun_symptom_table(df_ode_solution_tidy,df_symptom_prop,bindcolnames=c("infection","agegroup"))
# sum of 1,2,3rd infections
df_ode_solution_tidy_cases_sum=df_ode_solution_tidy_cases %>% group_by(t_years,compartment,agegroup,agegroup_name) %>% 
  summarise(symptom_cases=sum(symptom_cases),symptom_cases_fract=sum(symptom_cases_fract))
# time window of plot
vars_to_show=df_ode_solution_tidy_cases_sum$t_years>6.75 & df_ode_solution_tidy_cases_sum$t_years<10.15 & 
             df_ode_solution_tidy_cases_sum$agegroup<7; plot_perms=permutations(n=2,r=2,repeats.allowed=T)
for (k_plot in 1:nrow(permutations(n=2,r=2,repeats.allowed=T))){
scale_val=c("free","fixed")[plot_perms[k_plot,1]]; vartype=c("symptom_cases_fract","symptom_cases")[plot_perms[k_plot,2]]
if (grepl("fract",vartype)){y_axis_tag='fraction'} else {y_axis_tag="# cases"}
ggplot(df_ode_solution_tidy_cases_sum[vars_to_show,],aes(x=t_years,y=symptom_cases_fract,group=compartment,color=compartment)) +
  geom_line(size=1.05) + theme_bw() + standard_theme + theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=8)) + 
  facet_wrap(~agegroup_name,scales=scale_val) + # df_ode_solution_tidy_cases$agegroup
  scale_x_continuous(breaks=xval_breaks,minor_breaks=seq(0,max_time/365,by=1/12)) + scale_y_continuous(breaks=(0:10)/10) + 
  geom_rect(aes(xmin=shutdown_limits[1]/365,xmax=shutdown_limits[2]/365,ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
  xlab('years') + ylab(y_axis_tag) + ggtitle('symptomatic cases by age group (age-structured SIRS model)') + 
  theme(legend.position='none')
# SAVE
symptvals=length(unique(df_symptom_prop[,c("sympt_value")])); symptagegrval=length(unique(df_symptom_prop[,c("sympt_group")]))
if (symptvals>symptagegrval){dep_tag="_ageexpdep"}else{dep_tag="_agedep_only"}
if (length(unique(round(array(delta_susc_prop),5)))==length(delta_primary)){foldername="suscept_noagedep"}else{foldername="suscept_agedep"}
full_filename=paste0("simul_output/",foldername,paste0("/symptomcases_sever",dep_tag),"/RSV_DE_symptomcases_sever",
                     dep_tag,"_y",scale_val,'_',gsub("# ","",y_axis_tag),".png")
ggsave(full_filename,width=32,height=22,units="cm") }

### Age distrib before and after shutdown ----------------------------
df_ode_solution_tidy_cases_sum$season=findInterval(df_ode_solution_tidy_cases_sum$t_years,c(0,(1:10)+0.51))
season_peaks=df_ode_solution_tidy_cases_sum %>% group_by(agegroup,agegroup_name,season) %>% summarise(max_case=max(symptom_cases)) %>%
  group_by(season) %>% mutate(max_case_share_season=max_case/sum(max_case))
season_peaks[,"pre_post_shtdn"]="pre"; season_peaks$pre_post_shtdn[season_peaks$season>=ceiling(max(shutdown_limits/365))]="post"
season_truthvals=season_peaks$agegroup<=floor(shutdown_limits/365)[1] & !season_peaks$season==floor(shutdown_limits/365)[2] & season_peaks$season>3
ylimvals=c(0.05,0.16)
# PLOT
ggplot(season_peaks[season_truthvals,],aes(x=agegroup_name,y=max_case_share_season,group=season,color=factor(season))) + 
  geom_line(aes(size=pre_post_shtdn)) + theme_bw() + standard_theme + 
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),legend.text=element_text(size=12),legend.title=element_text(size=14)) + 
  xlab("") + ylab("fraction of all cases in season") + scale_size_manual(values=c(3,1)) + # range=c(1, 2), guide=FALSE 
  scale_y_continuous(breaks=seq(ylimvals[1],ylimvals[2],by=0.01),limits=ylimvals) +
  labs(color="season",shape="before/after shutdown",size="before/after shutdown") # + ylim(ylimvals)
# SAVE
if (length(unique(round(array(delta_susc_prop),5)))==length(delta_primary)){foldername="suscept_noagedep"}else{foldername="suscept_agedep"}
symptvals=length(unique(df_symptom_prop[,c("sympt_value")])); symptagegrval=length(unique(df_symptom_prop[,c("sympt_group")]))
if (symptvals>symptagegrval){dep_tag="_ageexpdep"}else{dep_tag="_agedep_only"}
ggsave(paste0("simul_output/",foldername,"/age_distrib_under5yr_suscept",dep_tag,".png"),width=32,height=22,units="cm")

### UK RSV data  --------------------------------------------------------
resp_virus_data_uk=read_csv("data/Respiratory viral defections by any method UK Ages.csv")
resp_virus_data_uk_tidy=resp_virus_data_uk %>% pivot_longer(!c("Year","startweek","Age"))
resp_virus_data_uk_tidy$Age=factor(resp_virus_data_uk_tidy$Age,levels=unique(resp_virus_data_uk_tidy$Age))
leaveout_year=c(2014); truthvals_rsv=resp_virus_data_uk_tidy$name %in% "RSV" & (!resp_virus_data_uk_tidy$Year %in% leaveout_year)
# means across years
averages_years=data.frame(resp_virus_data_uk_tidy[truthvals_rsv,] %>% group_by(startweek,Age) %>% 
                            summarise(value=mean(value)),Year='mean',name='RSV')
if (!any(resp_virus_data_uk_tidy$Year %in% "mean")){
resp_virus_data_uk_tidy=rbind(resp_virus_data_uk_tidy,averages_years[,colnames(resp_virus_data_uk_tidy)])}
truthvals_rsv=(resp_virus_data_uk_tidy$name %in% "RSV") & resp_virus_data_uk_tidy$Year %in% c("mean",2020)
  # (!resp_virus_data_uk_tidy$Year %in% c(leaveout_year,))
resp_virus_data_uk_tidy[,"type"]="indiv year"; resp_virus_data_uk_tidy[resp_virus_data_uk_tidy$Year %in% "mean","type"]="5-year average"
resp_virus_data_uk_tidy[,"width"]=1.01; resp_virus_data_uk_tidy[resp_virus_data_uk_tidy$Year %in% "mean","width"]=1.015
# plot
ggplot(resp_virus_data_uk_tidy[truthvals_rsv,],aes(x=startweek,y=value,group=Year,color=Year)) + # ,linetype=factor(type)
  geom_line(size=1.25) + geom_point(aes(shape=factor(type),size=factor(width))) + facet_wrap(~Age,scales="free") + 
  scale_x_continuous(breaks=(0:10)*5) + scale_linetype_manual(values=c("solid","dashed")) + theme_bw() + standard_theme + ylab("") +
  scale_size_manual(values = c(1.5,2),guide=FALSE) + labs(shape="data type") # ,linetype="data type"
# scale_linetype_manual(values=c("solid","dashed")) + 
# SAVE
ggsave("simul_output/uk_rsv_data2020.png",width=32,height=16,units="cm")
# resp_virus_data_uk_tidy[truthvals_rsv,] %>% group_by(Year,Age) %>% filter(value==max(value))