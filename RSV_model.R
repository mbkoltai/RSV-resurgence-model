# install.packages("devtools"); library("devtools")
# install_github("SineadMorris/shinySIR")

# ode solving, maximum likelihood, rcpp
# contact data from https://bisaloo.github.io/contactdata/index.html (Prem 2017)
# functions
rm(list=ls()); currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
# library(contactdata); library(fitdistrplus);  library(bbmle); library(Rcpp); library(GillespieSSA)
lapply(c("tidyverse","deSolve","gtools","rstudioapi","wpp2019","plotly","Rcpp","zoo","lubridate","tsibble"), library, character.only=TRUE)
source('fcns/RSV_model_functions.R')
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# SET PARAMETERS --------------------------------------------------------
# selected country
country_sel="United Kingdom"
# time resolution (in days)
elem_time_step=0.5
# population data
standard_age_groups <- fun_cntr_agestr(country_sel,"2015",seq(0,75,5),c(seq(4,74,5),99))
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
omega=1/150 # 1/runif(1,60,200)
# RECOVERY
rho=1/7 # 1/rho=rweibull(1, shape=4.1,scale=8.3)
# KINETIC MATRIX (aging terms need to be scaled by duration of age groups!)
K_m=fun_K_m_sirs_multiage(dim_sys,n_age,n_inf,n_compartment,rho,omega,varname_list,rsv_age_groups)
### ### ### ### ### ### ### ###
# BIRTH RATE into S_1_1 (Germany 2019: 778e3 births)
birth_rate=713e3/365; birth_term=matrix(c(birth_rate,rep(0,dim_sys-1)),dim_sys,1) # 0.01
# SUSCEPTIBILITY 
# normalised by age group sizes, for infection terms ~ delta*(I1+I2+...+In)*S_i/N_i
# agedep_fact determines strength of age dependence, if agedep_fact>1, decreasing susceptibility with age
agedep_fact=1.4; delta_primary=rep(0.25,3) # c(0.1,0.03,0.01)
delta_susc=sapply(1:n_age, function(x) {delta_primary/((agedep_fact^(x-1))*rsv_age_groups$value[x])})
delta_susc_prop=delta_susc*matrix(rep(rsv_age_groups$value,3),nrow=3,byrow=T)
### PLOT susceptibility: 
fcn_plot_suscept_table(fcn_suscept_agedeptable(rsv_age_groups,delta_susc,n_inf)) # ggsave("simul_output/suscept.png",width=32,height=22,units="cm")
# calculate R0 (at max seasonal forcing=1)
R0_calc_SIRS(C_m,delta_susc_prop,rho,n_inf) # scale_f=3; delta_susc_prop=delta_susc_prop*scale_f; delta_susc=delta_susc*scale_f
####
# DURATION of SIMULATION
n_years=20.25; max_time=n_years*n_days_year; timesteps <- seq(0,max_time,by=elem_time_step)
# seasonal forcing (above baseline level of 1) | npi_reduc_strength: reduction from baseline 
# shutdown season (if x --> (x+1)th season shut down) | preseas_npi_on/postseas_npi_off: on/off NPI week before/after season onset
forcing_strength=1; npi_reduc_strength=0; npi_year=round(n_years-5); preseas_npi_on=2; postseas_npi_off=2
g(forcing_vector_npi,shutdwn_lims,seas_force,seas_lims) %=% fun_shutdown_seasforc(timesteps,elem_time_step=0.5,
      forcing_strength,npi_reduc_strength,npi_year,peak_week=46,season_width=4,preseas_npi_on,postseas_npi_off,n_prec=0.01,n_sd=2) 
# set seas lims from UK data: peak is weeks 49/50, on/off is 41,11
seas_lims$on=floor(seas_lims$on)+41*7/365; seas_lims$off=round(seas_lims$off)+11*7/365
# PLOT seasonal forcing with NPI
fcn_plot_seas_forc(timesteps,seas_force,forcing_vector_npi,shutdwn_lims,seas_lims) 
# SAVE: ggsave(paste0("simul_output/NPI_y",npi_year,"_on",preseas_npi_on,"w_off",postseas_npi_off,"w.png"),units="cm",height=10,width=20)

# INITIAL CONDITIONS. # introduce stationary state as init state? | init_set: "previous" or anything else
prev_or_new=c("previous","fromscratch")[1]; filename="simul_output/df_ode_solution_UK_long.RDS"
# init_cond_src: "file" or "output"? 
initvals_sirs_model=fcn_set_initconds(init_set=prev_or_new,init_cond_src="output",ode_solution,init_seed=10,seed_vars="all",filename)
# manually choose a timepoint from prev simul
# initvals_sirs_model=matrix(as.numeric(round(ode_solution[min(which(round(timesteps/365,3)==19.5)),2:ncol(ode_solution)])))

### integrate ODE --------------------------------------------------------
params=list(birth_term,K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds,elem_time_step,delta_susc)
# interpolation fcns for seas forcing & extern introds
approx_seas_forc<-approxfun(data.frame(t=timesteps,seas_force=forcing_vector_npi))
approx_introd<-approxfun(data.frame(t=timesteps,as.numeric(timesteps %% 30==0)*10))
tm<-proc.time(); ode_solution<-lsoda(initvals_sirs_model,timesteps,func=sirs_seasonal_forc,parms=params); round(proc.time()-tm,2)
# reshape data | # check size of objs: fcn_objs_mem_use(1)
g(ode_solution,df_ode_solution_tidy) %=% fun_process_simul_output(ode_solution,varname_list,n_age,n_inf,rsv_age_groups,neg_thresh=-1e-3)
####
# PLOT how age group totals change (initial vs final popul sizes: fun_agegroup_init_final_pop(df_ode_solution_tidy))
fcn_plotagegroup_totals(df_ode_solution_tidy,scale_val=c('free_y','fixed')[1])
# ggsave(paste("simul_output/agegroup_totals_",scale_val,"yscale.png",sep=""),width=28,height=16,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Plot time course --------------------------------------------------------
xval_lims=c(npi_year-3.27,npi_year+4.5); seas_lims_plot=subset(seas_lims,on>xval_lims[1] & off<xval_lims[2]); agegr_lim=8
# PLOT
# tags: y-axis fixed/free | facet by AGE/(age&INFECTION) | abs values/fraction
# k=7 --> (free scale, infects separate, absval). k=8 --> (free, infects separate, fractional)
g(scale_val,facet2tag,value_type,y_axis_tag,ncol_val,facet_formula,foldername,caption_txt,subtitle_str,timecourse_filename) %=% 
    fun_tcourse_plottags(k=7,nval=2,rval=3,n_inf,n_age,colvar,agegr_lim,delta_susc_prop,delta_primary,
    preseas_npi_on,postseas_npi_off,npi_reduc_strength,forcing_strength)
# PLOT time course absolute values  --------------------------------------------------------
ggplot(subset(df_ode_solution_tidy,grepl('I',name) & agegroup<=agegr_lim & t/365>xval_lims[1] & t/365<xval_lims[2]),
    aes(x=t/365,y=get(value_type),group=name)) + geom_area(aes(fill=infection),position=position_stack(reverse=T),color="black",size=0.25) +
  facet_wrap(~agegroup_name,ncol=2,scales=scale_val) + theme_bw() + standard_theme + theme(axis.text.x=element_text(size=11,vjust=0.5),
    axis.text.y=element_text(size=12),legend.position="top",legend.title=element_blank(),strip.text=element_text(size=12)) +
  scale_x_continuous(breaks=seq(0,max_time/365,by=1/4),expand=expansion(0.01,0)) + scale_y_continuous(expand=expansion(0.01,0)) +
  geom_rect(aes(xmin=shutdwn_lims[1]/365,xmax=shutdwn_lims[2]/365,ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
  geom_vline(data=seas_lims_plot %>% pivot_longer(cols=!season),aes(xintercept=value),color="blue",linetype="dashed",size=0.3) +
  xlab('years') + ylab(y_axis_tag) + labs(subtitle=subtitle_str,caption=caption_txt)
## SAVE
ggsave(timecourse_filename, width=32,height=28,units="cm")

# PLOT time course normalised by max per age group  --------------------------------------------------------
xval_lims=c(npi_year-1.27,npi_year+3.25); seas_lims_plot=subset(seas_lims,on>xval_lims[1] & off<xval_lims[2])
ggplot(subset(df_ode_solution_tidy,grepl('I',name) & agegroup<=agegr_lim & t/365>xval_lims[1] & t/365<xval_lims[2]) %>%
         group_by(t,agegroup) %>% mutate(sum_val=sum(value)) %>% group_by(agegroup) %>% mutate(value_max_norm=value/max(sum_val)),
       aes(x=t/365,y=value_max_norm,group=name)) + geom_area(aes(fill=infection),color="black",size=0.25,position=position_stack(reverse=T)) + 
  facet_wrap(~agegroup_name,ncol=2,scales=scale_val) + theme_bw() + standard_theme + theme(axis.text.x=element_text(size=11,vjust=0.5),
  axis.text.y=element_text(size=12),legend.title=element_blank(),strip.text=element_text(size=12)) + #legend.position="top",
  scale_x_continuous(breaks=seq(0,max_time/365,by=1/4),expand=expansion(0.01,0)) + scale_y_continuous(expand=expansion(0.01,0)) +
  geom_rect(aes(xmin=shutdwn_lims[1]/365,xmax=shutdwn_lims[2]/365,ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
  geom_vline(data=seas_lims_plot %>% pivot_longer(cols=!season),aes(xintercept=value),color="blue",linetype="dashed",size=0.3) + 
  xlab('years') + ylab("% of maximum incidence") + labs(subtitle=subtitle_str,caption=caption_txt)
# SAVE
ggsave(gsub(".png","_byage_maxnorm.png",timecourse_filename), width=30,height=20,units="cm")

### PLOT time course faceted by infection --------------------------------------------------------
xval_lims=c(npi_year-1.27,npi_year+3.25); seas_lims_plot=subset(seas_lims,on>xval_lims[1] & off<xval_lims[2])
ggplot(subset(df_ode_solution_tidy,grepl('I',name) & t/365>xval_lims[1] & t/365<xval_lims[2]) %>% 
         group_by(t,infection) %>% mutate(sum_val=sum(value)) %>% group_by(infection) %>% mutate(value_max_norm=value/max(sum_val)),
  aes(x=t/365,y=value_max_norm,group=name)) + geom_area(aes(fill=agegroup_name),color="black",size=0.25,position=position_stack(reverse=T)) + 
  facet_wrap(~infection,ncol=1,scales=scale_val) + theme_bw() + standard_theme + theme(axis.text.x=element_text(size=11,vjust=0.5),
      axis.text.y=element_text(size=12),legend.position="top",legend.title=element_blank(),strip.text=element_text(size=12)) +
  scale_x_continuous(breaks=xval_breaks,expand=expansion(0,0)) + scale_y_continuous(expand=expansion(0.01,0)) +
  geom_rect(aes(xmin=shutdwn_lims[1]/365,xmax=shutdwn_lims[2]/365,ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
  geom_vline(data=seas_lims_plot %>% pivot_longer(cols=!season),aes(xintercept=value),color="blue",linetype="dashed",size=0.3) + 
  xlab('years') + ylab("% of maximum incidence") + labs(subtitle=subtitle_str,caption=caption_txt)
# SAVE
ggsave(gsub(".png","_byinf_maxnorm.png",timecourse_filename), width=30,height=20,units="cm")

### share of infections --------------------------------------------------------
xval_lims=c(npi_year-2.25,npi_year+3.25)
ggplot(subset(df_ode_solution_tidy,grepl('I',name) & t/365>xval_lims[1] & t/365<xval_lims[2] & agegroup<9) %>% group_by(t,agegroup) %>% 
mutate(value_fract_inf=value/sum(value)) %>% group_by(name) %>% mutate(value_fract_inf_smooth=rollmean(value_fract_inf,k=21,align="center",fill=NA)), 
  aes(x=t/365,y=value_fract_inf_smooth,group=name)) + geom_area(aes(fill=infection),color="black",size=0.25,position=position_stack(reverse=T)) + 
  facet_wrap(~agegroup_name,ncol=4,scales=scale_val) + theme_bw() + standard_theme + theme(axis.text.x=element_text(size=11,vjust=0.5),
      axis.text.y=element_text(size=12),legend.position="top",legend.title=element_blank(),strip.text=element_text(size=12)) +
  scale_x_continuous(breaks=xval_breaks,expand=expansion(0,0)) + scale_y_continuous(expand=expansion(0.01,0)) +
  geom_rect(aes(xmin=shutdwn_lims[1]/365,xmax=shutdwn_lims[2]/365,ymin=0,ymax=Inf),fill="grey",color=NA,alpha=0.01) +
  geom_vline(data=subset(seas_lims %>% pivot_longer(cols=!season),value>xval_lims[1]&value<xval_lims[2]),aes(xintercept=value),color="black",linetype="dashed",size=0.3) + 
  xlab('years') + ylab("% of cases at t") + labs(subtitle=subtitle_str,caption=caption_txt)
# SAVE
ggsave(gsub(".png","_byinf_share.png",timecourse_filename), width=32,height=18,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# dependence on age (>1) | dependence on prev exposure
agedep_fact=1; clinfract_expos=rep(0.4,n_inf) # c(0.8,0.4,0.2) # rep(0.75,3) # c(0.5,0.125,0.05) # rep(0.5,3)
clin_fract_age_exp=fcn_clin_fract_age_exp(agedep_fact,clinfract_expos,rsv_age_groups,n_age,n_inf)
# plot: fcn_plot_clinfract_table(clin_fract_age_exp)
### CALCULATE SYMPTOMATIC CASES & sum of 1,2,3rd infections -------------------------------
df_ode_sol_cases_sum=fun_symptomcases_table(df_ode_solution_tidy,clin_fract_age_exp,c("infection","agegroup"),seas_lims)

# PLOT SYMPTOMATIC CASES -----------------------------------
# k_plot: free/fixed | fraction/cases
g(scale_val,y_axis_tag,plotvar,foldername,full_filename,caption_txt,subtitle_str) %=% fun_sumcase_plot_tags(n_val=2,r_val=2,k_plot=2,
        clin_fract_age_exp,delta_susc_prop,preseas_npi_on,postseas_npi_off,npi_reduc_strength,forcing_strength)
xval_lims=c(npi_year-1.5,npi_year+3.5); seas_lims_plot=subset(seas_lims,on>xval_lims[1]&off<xval_lims[2]) %>% pivot_longer(col=!season)
# plot
ggplot(subset(df_ode_sol_cases_sum,t/365>xval_lims[1] & t/365<xval_lims[2] & agegroup<=9),aes(x=t/365,y=get(plotvar))) + geom_line() +
        facet_wrap(~agegroup_name,scales=scale_val,ncol=3) + theme_bw() + standard_theme +
  theme(axis.text.x=element_text(size=12,vjust=0.5),axis.text.y=element_text(size=13),strip.text=element_text(size=15),legend.position='none') +
  scale_x_continuous(breaks=xval_breaks[xval_breaks>npi_year-3],minor_breaks=seq(0,max_time/365,by=1/12),expand=expansion(0,0)) + 
  scale_y_continuous(expand=expansion(0.01,0)) + # 
  geom_rect(aes(xmin=shutdwn_lims[1]/365,xmax=shutdwn_lims[2]/365,ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
  geom_vline(data=seas_lims_plot,aes(xintercept=value),color="blue",linetype="dashed",size=0.3) +
  xlab('years') + ylab(y_axis_tag) + labs(caption=caption_txt,subtitle=subtitle_str)
# SAVE
# full_filename=gsub('.png','_reducedresurg.png',full_filename)
ggsave(full_filename,width=31,height=22,units="cm")

### symptom cases stacked as area plot
ggplot(subset(df_ode_sol_cases_sum,t/365>npi_year-0.35 & t/365<npi_year+3.3),aes(x=t/365,y=get(plotvar))) + #  & agegroup<=9
  geom_area(aes(fill=agegroup_name),position=position_stack(reverse=T),color="black",size=0.25) + theme_bw() + standard_theme +
  theme(axis.text.x=element_text(size=12,vjust=0.5),axis.text.y=element_text(size=13),strip.text=element_text(size=15)) +
  scale_x_continuous(breaks=xval_breaks[xval_breaks>npi_year-3],minor_breaks=seq(0,max_time/365,by=1/12),expand=expansion(0,0)) + 
  scale_y_continuous(expand=expansion(0,0)) +
  geom_rect(aes(xmin=shutdwn_lims[1]/365,xmax=shutdwn_lims[2]/365,ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
  geom_vline(data=subset(seas_lims_plot,value>npi_year-0.275),aes(xintercept=value),color="blue",linetype="dashed",size=0.3) +
  xlab('years') + ylab(y_axis_tag) + labs(fill="",caption=caption_txt,subtitle=subtitle_str)
# SAVE
ggsave(gsub("symptomcases","symptomcases_areaplot",full_filename),width=31,height=22,units="cm")

### Age distrib before and after shutdown ----------------------------
season_peaks_AUC=fcn_seas_agedistrib(df_ode_sol_cases_sum,max_time,timesteps,seas_case_threshold=1e3) %>%
  pivot_longer(cols=!c(agegroup,agegroup_name,season,agegr_size,mean_age,pre_post_shtdn))
# PLOT
plot_season_agedist=subset(season_peaks_AUC,season>npi_year-1 & season<=npi_year+4 & grepl("max_case",name)) %>% #  & agegroup<8
  mutate(agegroup_merged=as.character(agegroup_name)) %>% 
  mutate(agegroup_merged=factor(replace(agegroup_merged,as.numeric(agegroup_name)>7,"age>5yr"))) %>%
  group_by(season,name,agegroup_merged,pre_post_shtdn) %>% summarise(value=sum(value)) %>% 
  group_by(season,name) %>% mutate(value=replace(value,grepl("share",name),value/sum(value))) %>%
  mutate(value=replace(value,value>1,value/1e6)) %>% mutate(name=replace(name,name %in% "max_case_share_season","share of cases at peak")) %>%
  mutate(name=replace(name,name %in% "max_case","number of cases at peak (million)"))
if (length(unique(plot_season_agedist$pre_post_shtdn))==3) {colorvals=c(2,1,3)} else {colorvals=c(2,3)}
# bar plot
ggplot(plot_season_agedist,aes(x=season,y=value,group=rev(agegroup_merged))) +
 geom_bar(aes(fill=factor(pre_post_shtdn),alpha=agegroup_merged),color="black",stat="identity") + 
 scale_fill_manual(values=gg_color_hue(3)[colorvals]) + facet_wrap(~name,scales="free",nrow=2) + theme_bw() + standard_theme + 
 theme(axis.text.x=element_text(size=12,angle=0),legend.position="bottom",legend.text=element_text(size=12),
       legend.title=element_text(size=14),axis.text.y=element_text(size=12),strip.text=element_text(size=13)) + 
 scale_x_continuous(breaks=seq(round(xval_lims)[1],round(xval_lims)[2]+1,1)) + scale_y_continuous(breaks=,expand=expansion(c(0, 0.02))) + 
 geom_text(aes(x=season,y=value,label=paste0(gsub("age=","",agegroup_merged),"\n",round(value,2))),size=3.5,
           position=position_stack(vjust=0.5),check_overlap=T) +labs(fill="",alpha="",subtitle=subtitle_str,caption=caption_txt) + 
  guides(linetype=guide_legend(override.aes=list(fill=c(NA,NA,NA)))) + xlab("epi-year") + ylab("") + coord_flip()
# SAVE
agedistr_filename=fcn_agedistrib_plot_tags(delta_susc_prop,delta_primary,plot_season_agedist,clin_fract_age_exp,
                                           preseas_npi_on,postseas_npi_off,npi_reduc_strength,seas_case_threshold=1e2)
ggsave(agedistr_filename,width=32,height=16,units="cm")

### ### ### ### ### 
### mean age of infections ----------------------------
mean_age_perseason=subset(season_peaks_AUC,grepl("share",name)) %>% group_by(season,name) %>%
  summarise(under2yr=sum(value[agegroup<5]*mean_age[agegroup<5]/sum(value[agegroup<5])),
            under3yr=sum(value[agegroup<6]*mean_age[agegroup<6]/sum(value[agegroup<6])),
            under5yr=sum(value[agegroup<7]*mean_age[agegroup<7]/sum(value[agegroup<7])),all_agegroups=sum(value*mean_age/sum(value)),
  pre_post_shtdn=unique(pre_post_shtdn)) %>% pivot_longer(cols=!c(season,name,pre_post_shtdn),names_to="mean_age_type") %>%
  mutate(mean_age_type=factor(mean_age_type,unique(mean_age_type))) %>% group_by(name,mean_age_type) %>% 
  mutate(mean_val=mean(value[pre_post_shtdn=="pre-NPI" & season>3])) %>% 
  mutate(norm_val=value/mean_val,dep=gsub("_dep","",gsub("suscept_","",foldername)))
# PLOT
# ggplot(subset(mean_age_perseason,season>3 & grepl("max",name)), #  & mean_age_type %in% "mean_age_under5"
#   aes(x=season,y=value,fill=factor(pre_post_shtdn),group=mean_age_type,alpha=factor(mean_age_type))) +
#   geom_bar(stat="identity",color="black",position="dodge") + scale_alpha_manual(values=1-(0:(length(unique(mean_age_perseason$mean_age_type))-1))*0.2) +
#   scale_x_continuous(breaks=0:max(mean_age_perseason$season)) + scale_y_continuous(breaks=0:max(mean_age_perseason$value),expand=expansion(c(0,0.1))) +
#   scale_fill_manual(values=gg_color_hue(3)[colorvals],guide='none') + theme_bw() + standard_theme +
#   theme(axis.text.x=element_text(size=9,vjust=0.5),axis.text.y=element_text(size=8),legend.position="bottom",legend.text=element_text(size=12)) +
#         labs(alpha="",subtitle=subtitle_str,caption=caption_txt) + ylab("mean age of symptomatic cases") +
#   geom_text(aes(label=round(value,2)),size=4,color="black",alpha=1,vjust=-0.5,position=position_dodge(width=1)) # + coord_flip()
# # SAVE
# ggsave(gsub("age_distrib","mean_age",agedistr_filename),width=22,height=14,units="cm")

# segment plot
var_sel="value"; if(grepl("norm",var_sel)) {ylab_tag=" (normalised)"} else {ylab_tag=" (year)"}
ggplot(subset(left_join(mean_age_perseason,seas_lims %>% mutate(season=season+1)),season>npi_year-1 & season<npi_year+6 & grepl("max",name) & 
                as.numeric(mean_age_type)>2) %>% mutate(year=season-1)) + 
  geom_segment(aes(x=on,xend=on+1-0.05,y=get(var_sel),yend=get(var_sel),color=dep),size=2) + facet_wrap(~mean_age_type,scales="free") + 
  geom_rect(aes(xmin=shutdwn_lims[1]/365,xmax=shutdwn_lims[2]/365,ymin=-Inf,ymax=Inf),fill="pink",color=NA,alpha=0.1,show.legend=TRUE) +
  geom_vline(data=subset(seas_lims %>% pivot_longer(cols=!season),value>npi_year-2 & value<npi_year+4 & name %in% "on"),
             aes(xintercept=value),linetype="dashed",size=0.3) +
  scale_x_continuous(breaks=0:ceiling(xval_lims[2]),expand=expansion(0,0)) + scale_y_continuous(expand=expansion(0.15,0)) + 
  theme_bw() + standard_theme + theme(axis.text.x=element_text(size=12,vjust=0.75),axis.text.y=element_text(size=12), 
  legend.title=element_text(size=12),legend.text=element_text(size=12)) + xlab("epi-year") +ylab(paste0("<age symptomatic cases>",ylab_tag)) +
  labs(color="dependence",fill="NPI",caption=caption_txt) # guide_legend(override.aes=list(colour="red"))
# SAVE
ggsave(gsub("pct","pct_geomsegment",gsub("age_distrib","mean_age",agedistr_filename)),width=22,height=14,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### UK RSV data  --------------------------------------------------------
resp_virus_data_uk=read_csv("data/Respiratory viral detections by any method UK Ages.csv")
resp_virus_data_uk_tidy = resp_virus_data_uk %>% pivot_longer(!c("Year","startweek","Age")) %>% 
  mutate(Age=factor(gsub(" Y","Y",resp_virus_data_uk_tidyAge),levels=unique(gsub(" Y","Y",resp_virus_data_uk_tidyAge))))
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
# convert weel number to date
resp_virus_data_uk_tidy[,"date"]=as.Date(paste(resp_virus_data_uk_tidy$Year,resp_virus_data_uk_tidy$startweek,1,sep="-"),"%Y-%U-%u")
resp_virus_data_uk_tidy[,"year_week"]=gsub("mean-","",paste0(resp_virus_data_uk_tidy$Year,"-",resp_virus_data_uk_tidy$startweek))
resp_virus_data_uk_tidy$year_week=factor(resp_virus_data_uk_tidy$year_week,unique(resp_virus_data_uk_tidy$year_week))
# PLOT all years
ggplot(subset(resp_virus_data_uk_tidy,name %in% "RSV" & grepl("indiv",type)),aes(x=year_week,y=value,group=Age)) + 
  geom_area(aes(fill=Age),position=position_stack(reverse=T)) + 
  # geom_line(aes(color=Age)) + geom_point(size=0.4) + facet_wrap(~Age,ncol=2,scales="free") + # 
  xlab("Year, week") + ylab("number of reported RSV cases (4-week period)") + theme_bw() + standard_theme + 
  theme(axis.text.x=element_text(size=12,vjust=0.5),axis.text.y=element_text(size=13),legend.position="bottom")
#  labs(caption="source: gov.uk/health-and-social-care/health-protection-infectious-diseases")
# ggsave("simul_output/uk_rsv_data2014_2020.png",width=32,height=16,units="cm")
# ggsave("simul_output/uk_rsv_data2014_2020_age_facet.png",width=32,height=16,units="cm")
# plot MEAN vs 2020
ggplot(resp_virus_data_uk_tidy[truthvals_rsv,],aes(x=startweek,y=value,group=Year,color=Year)) + # ,linetype=factor(type)
  geom_line(size=1.25) + geom_point(aes(shape=factor(type),size=factor(width))) + facet_wrap(~Age,scales="free") + 
  scale_x_continuous(breaks=(0:10)*5) + scale_linetype_manual(values=c("solid","dashed")) + theme_bw() + standard_theme + ylab("") +
  scale_size_manual(values = c(1.5,2),guide=FALSE) + labs(shape="data type") # ,linetype="data type"
# SAVE
ggsave("simul_output/uk_rsv_data2020_multiyearaverage.png",width=32,height=16,units="cm")

### data with weekly resolution (no age resol) ----
x=read_csv("data/Respiratory viral detections by any method UK.csv")
resp_detects_weekly_all_age= x %>% mutate(year_week=factor(paste0(Year,"-",Week),unique(paste0(Year,"-",Week))), 
    RSV_rolling_av=rollmean(RSV,k=7,align="center",fill=NA),section=ceiling((as.numeric(year_week)*0.92)/100) ) %>% 
  mutate(section=ifelse(section>3,3,section)) # year_week=factor(year_week,unique(year_week)),
# plot
ggplot(resp_detects_weekly_all_age,aes(x=year_week)) + geom_point(aes(y=RSV,group=1,color=factor(Year))) + 
  geom_line(aes(y=RSV_rolling_av,group=1)) + labs(color="Year") + scale_y_continuous(expand = expansion(0,0.05)) + 
  geom_vline(data=subset(resp_detects_weekly_all_age, Week %in% c(40,11)),aes(xintercept=year_week),color="red",linetype="dashed",size=0.5) + 
  facet_wrap(~section,nrow=3,scale="free_x") + theme_bw() + standard_theme + theme(axis.text.x = element_text(vjust=0.5,size=6))
# save
ggsave("simul_output/uk_rsv_data2020_allagegroups_weekly.png",width=32,height=20,units="cm")