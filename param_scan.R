# this script is a 2D parscan to explore the spectrum of dynamical behaviors ~ f(basal_rate,onset/offset_NPI)
rm(list=ls())
lapply(c("tidyverse","deSolve","gtools","rstudioapi","wpp2019","Rcpp","zoo","wesanderson","qs","lubridate","tsibble"),library,character.only=TRUE)
currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
source('fcns/RSV_model_functions.R')
### define constant parameters ----
# time resolution (in days)
elem_time_step=0.5
# population data
country_sel="United Kingdom"; standard_age_groups <- fun_cntr_agestr(country_sel,"2015",seq(0,75,5),c(seq(4,74,5),99))
# RSV age groups (population data from wpp2019)
rsv_age_groups_low=c(0,0.5,1,1.5, 2,3,4, 5,10,15, 20); rsv_age_group_sizes=c(rep(0.4,4),rep(0.9,3),rep(4,3),79)
rsv_age_groups=fun_rsv_agegroups(standard_age_groups,rsv_age_groups_low,rsv_age_group_sizes)
# population by age group # number of age groups, reinfections and variables (S,I,R)
N_tot=sum(rsv_age_groups$value); n_age=nrow(rsv_age_groups); varname_list=c('S','I','R'); n_compartment=length(varname_list); n_inf=3
dim_sys=n_age*n_compartment*n_inf; n_days_year=365
# linear indices of the I & S variables
l_inf_susc=fun_inf_susc_index_lists(n_age,n_inf,varname_list); inf_vars_inds=l_inf_susc[[1]]; susc_vars_inds=l_inf_susc[[2]]
# CONTACT MATRIX
# contact matrix from covidm ("home","work","school","other")
cm_path="~/Desktop/research/models/epid_models/covid_model/lmic_model/covidm/"
list_contmatrs=fun_covidm_contactmatrix(country_sel="Germany",currentdir_path=currentdir_path,cm_path=cm_path); C_m_full=Reduce('+',list_contmatrs)
# make matrix symmetric
C_m_fullrecipr=fun_recipr_contmatr(C_m_full,age_group_sizes=standard_age_groups$values)
# create for our age groups
C_m=fun_create_red_C_m(C_m_fullrecipr,rsv_age_groups,
                       orig_age_groups_duration=standard_age_groups$duration,orig_age_groups_sizes=standard_age_groups$values)
# make it reciprocal for the larger group
C_m=fun_recipr_contmatr(C_m,age_group_sizes=rsv_age_groups$value)
# bc of reinfections we need to input contact matrix repeatedly
contmatr_rowvector=t(do.call(cbind, lapply(1:nrow(C_m), function(x){diag(C_m[x,]) %*% matrix(1,n_age,n_inf)}))); rm(cm_matrices)
# build kinetic matrix --------------------------------------------------------
# WANING (immunity) terms: R_i_j -> S_min(i+1,n_inf)_j | RECOVERY
omega=1/150; rho=1/7 # 1/rho=rweibull(1, shape=4.1,scale=8.3)
# KINETIC MATRIX (aging terms need to be scaled by duration of age groups!)
K_m=fun_K_m_sirs_multiage(dim_sys,n_age,n_inf,n_compartment,rho,omega,varname_list,rsv_age_groups)
### ### ### ### ### ### ### ###
# BIRTH RATE into S_1_1 (Germany 2019: 778e3 births)
birth_rate=713e3/365; birth_term=matrix(c(birth_rate,rep(0,dim_sys-1)),dim_sys,1) # 0.01
# length of simul
n_years=15.25; max_time=n_years*n_days_year; timesteps <- seq(0,max_time,by=elem_time_step)
# season limits (empirical, UK)
g(forcing_vector_npi,shutdwn_lims,seas_force,seas_lims) %=% fun_shutdown_seasforc(timesteps,elem_time_step=0.5,forcing_above_baseline=2,
    npi_strength=0.5,npi_year=round(n_years-4),peak_week=46,season_width=4,npi_on=2,npi_off=2,n_prec=0.01,n_sd=2)
seas_lims$on=floor(seas_lims$on)+41*7/365; seas_lims$off=round(seas_lims$off)+11*7/365
# params to scan in
scan_params=list(susc=seq(0.08,0.11,length.out=4),expdep=c(1.5,2,3),agedep=1,t_npi=c(1,3,5),seasf=c(1,3,5))
param_table=expand.grid(scan_params); colnames(param_table)=names(scan_params)
if (length(unique(param_table$agedep))==1) {sel_col=which(grepl("age",colnames(param_table))); foldertag="expdep_"} else {
  sel_col=which(grepl("exp",colnames(param_table))); foldertag="agedep_"}
# CREATE FOLDER
foldername <- paste0("simul_output/parscan/",foldertag,paste0(names(scan_params),sapply(scan_params, length),collapse = "_"))
dir.create(foldername); sink(paste0(foldername,"/timestamp.txt")); cat(format(Sys.time(), "%F %H-%M")); sink()
write_csv(param_table,paste0(foldername,"/param_table.csv"))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### PARSCAN ----
# cores=detectCores(); cl<-makeCluster(cores[1]-1); registerDoParallel(cl) # clusterExport(cl, c("fcn_somal_sir_singlesimul")); 
# output<-foreach(k=1:k_lim,.combine=rbind,.packages=c("tidyr","deSolve","dplyr","RcppRoll")) %dopar% {
# temp<-...}
for (k_paramscan in 1:nrow(param_table)) {
  if (k_paramscan==1){df_ode_sol_sympt_cases=list();
  initvals_sirs_model=fcn_set_initconds(init_set="previous",init_cond_src="file",ode_solution,init_seed=10,
                                        seed_vars="all","simul_output/df_ode_solution_UK_long.RDS")
      } # else {
    # initvals_sirs_model=fcn_set_initconds(init_set="previous",init_cond_src="output",ode_solution,init_seed=10,seed_vars="all",filename) }
# assign params
  g(suscept_abs_level,suscept_expdep,suscept_agedep,npi_lims,forcing_strength) %=% param_table[k_paramscan,]
# NPI onset/end
npi_reduc_strength=0.5; npi_year=round(n_years-4); preseas_npi_on=npi_lims; postseas_npi_off=npi_lims-3
g(forcing_vector_npi,shutdwn_lims,seas_force,seas_lims) %=% fun_shutdown_seasforc(timesteps,elem_time_step=0.5,
  forcing_strength,npi_reduc_strength,npi_year,peak_week=46,season_width=4,preseas_npi_on,postseas_npi_off,n_prec=0.01,n_sd=2)
# susceptibility
delta_primary=suscept_abs_level/suscept_expdep^(0:2)
delta_susc=sapply(1:n_age, function(x) {delta_primary/((suscept_agedep^(x-1))*rsv_age_groups$value[x])})
# SIMULATE
approx_seas_forc <- approxfun(data.frame(t=timesteps,seas_force=forcing_vector_npi))
approx_introd <- approxfun(data.frame(t=timesteps,as.numeric(timesteps %% 30==0)*100))
params=list(birth_term,K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds,elem_time_step,delta_susc)
tm<-proc.time(); ode_solution<-lsoda(initvals_sirs_model,timesteps,func=sirs_seasonal_forc,parms=params); 
# g(ode_solution,df_ode_solution_tidy) %=% fun_process_simul_output(ode_solution,varname_list,n_age,n_inf,rsv_age_groups,neg_thresh=-1e-3)
df_ode_solution_tidy=fun_process_simul_output(ode_solution,varname_list,n_age,n_inf,rsv_age_groups,neg_thresh=-1e-3)[[2]]
round(proc.time()-tm,2)
print(k_paramscan)
# calc symptomatic cases
agedep_fact=1; clinfract_expos=rep(1,n_inf) # c(0.8,0.4,0.2) # rep(0.75,3) # c(0.5,0.125,0.05) # rep(0.5,3)
clin_fract_age_exp=fcn_clin_fract_age_exp(agedep_fact,clinfract_expos,rsv_age_groups,n_age,n_inf)
### CALCULATE SYMPTOMATIC CASES & sum of 1,2,3rd infections -------------------------------
df_ode_sol_sympt_cases[[k_paramscan]]=left_join(
  data.frame(subset(fun_symptomcases_table(df_ode_solution_tidy,clin_fract_age_exp,c("infection","agegroup")),
                    t/365 > npi_year-2),parset_id=k_paramscan), 
  data.frame(round(param_table[k_paramscan,],2),parset_id=k_paramscan))
if (k_paramscan==nrow(param_table)) {df_ode_sol_sympt_cases=bind_rows(df_ode_sol_sympt_cases)}
}
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### end of parscan
# SAVE/LOAD
qsave(df_ode_sol_sympt_cases,paste0(foldername,"/df_ode_sol_sympt_cases_EXP_DEP_highersusc",nrow(param_table),"_parset.qs"))
# filter out when no signif outbreak
df_ode_sol_sympt_cases=df_ode_sol_sympt_cases %>% # select(-symptom_cases_fract) %>% 
  # group_by(parset_id) %>% filter(max(symptom_cases)>1e5) %>% 
  group_by(t,parset_id) %>% mutate(sumval=sum(symptom_cases)) %>% group_by(parset_id) %>% mutate(normval=symptom_cases/max(sumval))
####
# plot symptom cases time courses ------------------
xval_lims=c(npi_year-2.5,npi_year+3.25); seas_lims_plot=subset(seas_lims,on>xval_lims[1] & off<xval_lims[2]) %>% pivot_longer(cols=!season)
xval_breaks=seq(0,max_time/365,by=1/4) # npi_year=round(n_years-4)
# PLOT
selval="symptom_cases" # symptom_cases | normval
ggplot(subset(df_ode_sol_sympt_cases,t/365>xval_lims[1] & t/365<xval_lims[2] & t%%7==0),aes(x=t/365,y=get(selval)/1e6)) + 
  geom_area(aes(fill=agegroup_name),position=position_stack(reverse=T),color="black",size=0.1) + # symptom_cases
  facet_grid(expdep~t_npi~susc+seasf,scales="free",labeller=label_both) + theme_bw() + standard_theme + # parset_id ,ncol=6
  theme(axis.text.x=element_text(size=5,vjust=0.5),axis.text.y=element_text(size=6),strip.text=element_text(size=8),legend.position="bottom") + 
  scale_x_continuous(breaks=xval_breaks,minor_breaks=seq(0,max_time/365,by=1/12),expand=expansion(0,0)) + #
  scale_y_continuous(expand=expansion(0,0)) + xlab('years') + ylab(paste0("cases",gsub("symptom_cases"," (1e6)",selval))) + labs(fill="") + 
  geom_rect(aes(xmin=shutdwn_lims[1]/365,xmax=shutdwn_lims[2]/365,ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
  geom_vline(data=seas_lims_plot,aes(xintercept=value),color="blue",linetype="dashed",size=0.15)
# save
ggsave(paste0(foldername,"/timecourse_symptcases_",gsub("symptom_cases","absval",selval),nrow(param_table),"parsets",".png"),
       units="cm",height=20,width=40)

# PLOT single parset
ggplot(subset(df_ode_sol_sympt_cases,parset_id==min(parset_id) & t>3800 & t<4900),aes(x=t,y=sumval/1e6)) + 
  geom_area(aes(fill=agegroup_name),position=position_stack(reverse=T),color="black",size=0.1) + ylab("cases (1e6)") +
  scale_x_continuous(breaks=(37:50)*100,expand=expansion(0,0)) + scale_y_continuous(expand=expansion(0,0)) + 
  theme_minimal() + # theme_bw() + standard_theme +
  theme(axis.text.x=element_text(vjust=0.75)) + geom_vline(xintercept=4700,linetype="dashed",size=0.5)
# symptom cases at given timepoint
sympt_cases_t=subset(df_ode_sol_sympt_cases,parset_id==min(parset_id) & t==4700)$symptom_cases

uk_popul=left_join(subset(popF,name %in% "United Kingdom")[,c("age","2020")],
      subset(popM,name %in% "United Kingdom")[,c("age","2020")],by="age",suffix=c("F","M")) %>% 
  mutate(totalpop=`2020F`+`2020M`,lower=as.numeric(gsub("-\\d+","",age)),upper=as.numeric(gsub("\\d+-","",age))+0.9)
if (any(is.na(uk_popul$lower))){
uk_popul[grepl("95",uk_popul$age),c("2020F","2020M","totalpop")]=(uk_popul[grepl("95",uk_popul$age),c("2020F","2020M","totalpop")]+
  uk_popul[uk_popul$age=="100+",c("2020F","2020M","totalpop")]); uk_popul=uk_popul[-which(uk_popul$age=="100+"),] 
  uk_popul=uk_popul %>% mutate(mean_age=(lower+upper+0.1)/2,fraction_pop=totalpop/sum(totalpop)) }

agegroup_match=data.frame(model_agegroup=1:nrow(rsv_age_groups),age_low=rsv_age_groups$age_low,age_high=rsv_age_groups$age_high,
  wpp_agegroup_low=unlist(lapply(rsv_age_groups$age_low,function(x){which(x>=uk_popul$lower & x<=uk_popul$upper)})),
  wpp_agegroup_high=unlist(lapply(rsv_age_groups$age_high,function(x){which(x>=uk_popul$lower & x<=uk_popul$upper)}))) %>%
  mutate(age_high=ifelse(model_agegroup<max(model_agegroup),age_low[model_agegroup+1],age_high),
         mean_age_arithm=(age_low+age_high)/2, mean_age_weighted=sapply(1:max(model_agegroup),function(x) {
           sum((uk_popul$totalpop[wpp_agegroup_low[x]:wpp_agegroup_high[x]]*
              uk_popul$mean_age[wpp_agegroup_low[x]:wpp_agegroup_high[x]])/
                 sum(uk_popul$totalpop[wpp_agegroup_low[x]:wpp_agegroup_high[x]]))})) %>%
  mutate(mean_age_weighted=ifelse(wpp_agegroup_low==wpp_agegroup_high,mean_age_arithm,mean_age_weighted)) %>% select(-mean_age_arithm)
  

### ### ### ### ### ### ### ### ### ### ### ### ### ###  
# calculate age distributions ----------------------------
season_peaks_AUC=lapply(unique(df_ode_sol_sympt_cases$parset_id),
      function(x){data.frame(fcn_seas_agedistrib(subset(df_ode_sol_sympt_cases,parset_id==x),max_time,timesteps,seas_case_threshold=1e4) %>%
        pivot_longer(cols=!c(agegroup,agegroup_name,season,agegr_size,mean_age,pre_post_shtdn)),parset_id=x)})
plot_season_agedist=left_join(
  bind_rows(lapply(season_peaks_AUC,function(x){data.frame(
    fcn_agedistrib_calc_season(x,seaslims=c(npi_year-1,npi_year+4),selvar="max_case", # max_case auc_case
                              agegr_merge_min=which(rsv_age_groups$age_low==5)-1,rsv_age_groups),parset_id=unique(x$parset_id))})),
                              data.frame(round(param_table[,-sel_col],2),parset_id=1:nrow(param_table)),by="parset_id")
season_peaks_AUC=left_join(bind_rows(season_peaks_AUC),
  data.frame(round(param_table[,-sel_col],2),parset_id=1:nrow(param_table)),by="parset_id")
if (length(unique(plot_season_agedist$pre_post_shtdn))==3) {colorvals=c(2,1,3)} else {colorvals=c(2,3)}
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# plot age distribs  ----------------------------
xval_lims=c(npi_year-0.6,npi_year+4.25); seas_lims_plot=subset(seas_lims,on>xval_lims[1] & off<xval_lims[2]) %>% pivot_longer(cols=!season)
selvar="number" # number | share
ggplot(subset(plot_season_agedist,grepl(selvar,name)),aes(x=season,y=value,group=rev(agegroup_merged))) +
  geom_bar(aes(fill=factor(pre_post_shtdn),alpha=agegroup_merged),color="black",stat="identity") + 
  scale_fill_manual(values=gg_color_hue(3)[colorvals]) + 
  facet_grid(expdep~t_npi~susc+seasf,scales="free",labeller=label_both) + theme_bw() + standard_theme + 
  theme(axis.text.x=element_text(size=6,angle=90,vjust=0.75),legend.position="bottom",legend.text=element_text(size=10),
        axis.text.y=element_text(size=8),strip.text=element_text(size=8)) + # legend.title=element_text(size=14),
  scale_x_continuous(breaks=seq(round(xval_lims)[1],round(xval_lims)[2]+1,1)) + scale_y_continuous(expand=expansion(c(0,0.02))) + 
  labs(fill="",alpha="") + guides(linetype=guide_legend(override.aes=list(fill=c(NA,NA,NA)))) + xlab("epi-year") + ylab("") + coord_flip()
# geom_text(aes(x=season,y=value,label=paste0(gsub("age=","",agegroup_merged),"\n",round(value,2))),size=2,
#          position=position_stack(vjust=0.5),check_overlap=T) +
# SAVE
ggsave(paste0(foldername,"/age_distrib_exp_dep_",selvar,".png"),units="cm",height=20,width=40)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# plot shift in mean age
mean_age_perseason=left_join(bind_rows(lapply(unique(season_peaks_AUC$parset_id),
  function(x){data.frame(fcn_calc_mean_age(subset(season_peaks_AUC,parset_id==x),season_min=npi_year+1,dep_name="previous exposure"),parset_id=x)})),
                             data.frame(round(param_table[,-sel_col],2),parset_id=1:nrow(param_table)),by="parset_id")
# segment plot
var_sel="value"; # value norm_val
if(grepl("norm",var_sel)) {ylab_tag=" (normalised)"} else {ylab_tag=" (year)"}
xval_lims=c(npi_year-0.8,npi_year+4.25); 
seas_lims_plot=subset(seas_lims %>% pivot_longer(cols=!season),value>npi_year-1 & value<npi_year+4 & name %in% "on")
ggplot(subset(left_join(mean_age_perseason,seas_lims %>% mutate(season=season+1)),grepl("all_age",mean_age_type) & 
    season>npi_year & season<npi_year+5 & grepl("auc",name) & as.numeric(mean_age_type)>2) %>% mutate(year=season-1)) + 
  facet_grid(expdep~t_npi~susc+seasf,scales="free",labeller=label_both) + # facet_wrap(expdep~t_npi~susc+seasf,scales="free",switch="y") + 
  geom_segment(aes(x=on,xend=on+1-0.05,y=get(var_sel),yend=get(var_sel),color=dep),size=1.3) + 
  geom_rect(aes(xmin=shutdwn_lims[1]/365,xmax=shutdwn_lims[2]/365,ymin=-Inf,ymax=Inf),fill="pink",color=NA,alpha=0.2,show.legend=TRUE) +
  geom_vline(data=seas_lims_plot,aes(xintercept=value),linetype="dashed",size=0.3) + 
  geom_hline(yintercept=1,color="blue",size=0.3,linetype="dashed") +
 scale_x_continuous(breaks=0:ceiling(xval_lims[2]),expand=expansion(0,0),limits=summary(seas_lims_plot$value)[c("Min.","Max.")]) +
  scale_y_continuous(expand=expansion(0.2,0)) + theme_bw() + standard_theme + 
  theme(axis.text.x=element_text(size=10,vjust=0.75),axis.text.y=element_text(size=8),
        legend.title=element_text(size=12),legend.text=element_text(size=12)) +
  xlab("epi-year") + ylab(paste0("<age symptomatic cases>",ylab_tag)) + labs(color="dependence",fill="NPI")
# save
ggsave(paste0(foldername,"/mean_age_exp_dep_",var_sel,".png"),units="cm",height=25,width=40)
