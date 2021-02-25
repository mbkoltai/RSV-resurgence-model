# install.packages("devtools"); library("devtools")
# install_github("SineadMorris/shinySIR")

# ode solving, maximum likelihood, rcpp
# contact data from https://bisaloo.github.io/contactdata/index.html (Prem 2017)
# functions
rm(list=ls()); currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
# library(contactdata); library(fitdistrplus);  library(bbmle); library(Rcpp); library(GillespieSSA)
lapply(c("tidyverse","deSolve","gtools","rstudioapi","wpp2019","plotly","Rcpp","zoo"), library, character.only=TRUE)
source('RSV_model_functions.R')
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# SET PARAMETERS --------------------------------------------------------
# selected country
country_sel="Germany"
# time resolution (in days)
elem_time_step=0.5
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
omega=1/1.25e2 # 1/runif(1,60,200)
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
  xlab("age group (years)") + ylab("susceptibility") + ggtitle('susceptibility ~ f(age, # infection)')
# if (length(unique(round(suscept_agedep$value,6)))==3){dep_tag="expos_dep"}else{dep_tag="age_expos_dep"}
# ggsave(paste0("simul_output/suscept_",dep_tag,".png"),width=32,height=22,units="cm")
#
# calculate R0 (at max seasonal forcing=1)
R0_calc_SIRS(C_m,delta_susc_prop,rho,n_inf)
# BIRTH RATE into S_1_1 (Germany 2019: 778e3 births)
birth_rate=2130; birth_term=matrix(c(birth_rate,rep(0,dim_sys-1)),dim_sys,1) # 0.01
####
# DURATION of SIMULATION
n_years=7.25; max_time=n_years*n_days_year; timesteps <- seq(0,max_time,by=elem_time_step)
# seasonal forcing
basal_rate=0.15; peak_week=49; season_width=3.8 # (in weeks)
# shutdown season (if x --> (x+1)th season shut down) | preseas_npi_on/postseas_npi_off: on/off NPI week before/after season onset
npi_year=4; preseas_npi_on=2; postseas_npi_off=2 # shutdown_scale=0.2
shutdown_list=fun_shutdown_seasforc(timesteps,elem_time_step,basal_rate,npi_year,peak_week,season_width,
                                    preseas_npi_on,postseas_npi_off,n_prec=0.01,st_devs=2:3,n_sd=2)
forcing_vector_npi=shutdown_list[[1]];shutdwn_lims=shutdown_list[[2]];seas_force=shutdown_list[[3]]; seas_lims=shutdown_list[[4]]
shutdown_scale=forcing_vector_npi[round(mean(shutdwn_lims))]
# PLOT seasonal forcing with NPI
fcn_plot_seas_forc(timesteps,seas_force,forcing_vector_npi,shutdwn_lims,subset(seas_lims,sd==2))
# SAVE
# seasfor_filename=paste0("simul_output/seas_forcing_NPI_y",npi_year,"_on",preseas_npi_on,"w_off",postseas_npi_off,"w.png")
# ggsave(seasfor_filename,units="cm",height=10,width=20)

# INITIAL CONDITIONS
# introduce stationary state as init state? # initvals_sirs_model[unlist(lapply((0:10)*9,function(x ){(4:6) + x}))]
if (exists("df_ode_solution")){initsrc="output"} else {initsrc="file"; df_ode_solution=c()}
initvals_sirs_model=fcn_init_susc_vals(stationary_init=TRUE,from_file_or_output=initsrc,simul_output=df_ode_solution,
                          susc_vars_inds,agegr_sizes=rsv_age_groups$value,sim_filepath="simul_output/statsol_10years.RDS")
# INITIAL INFECTION (taking stationary sol normally contains some Inf.s so no need to re-seed it)
initvals_sirs_model[inf_vars_inds[[1]][1]]=10 # all first infection groups: sapply(inf_vars_inds, '[[',1)

### integrate ODE --------------------------------------------------------
# deSolve input
params=list(birth_term,K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds,forcing_vector_npi,elem_time_step,delta_susc)
# solve with interpolation
approx_seas_forc <- approxfun(data.frame(t=timesteps,seas_force=forcing_vector_npi))
ptm<-proc.time(); ode_solution_interpol<-lsoda(initvals_sirs_model,timesteps,func=sirs_seasonal_forc_interpol,parms=params); proc.time()-ptm  

# reshape data
list_simul_output=fun_process_simul_output(ode_solution,varname_list,n_age,n_inf,rsv_age_groups,neg_thresh=-0.001)
df_ode_solution=list_simul_output[[1]]; df_ode_solution_tidy=list_simul_output[[2]]; rm(list_simul_output)
# check size of objs: fcn_objs_mem_use()
####
# PLOT how age group totals change (initial vs final popul sizes: fun_agegroup_init_final_pop(df_ode_solution_tidy))
fcn_plotagegroup_totals(df_ode_solution_tidy,scale_val=c('free_y','fixed')[2])
# ggsave(paste("simul_output/agegroup_totals_",scale_val,"yscale.png",sep=""),width=28,height=16,units="cm")
####
# Plot time course --------------------------------------------------------
# xval_lims=c(floor(shutdwn_lims/365)[1]-0.15,ceiling(shutdwn_lims/365)[2]+1.25)
xval_lims=c(3.77,max_time/365); xval_breaks=seq(0,max_time/365,by=1/4); agegr_lim=6; colvar="inf"; n_sd=2
seas_lims_plot=subset(seas_lims,sd==n_sd & on>xval_lims[1] & off<xval_lims[2])
# PLOT
for (k in 7:8) {# 1:(2^3)
# tags: y-axis fixed/free | facet by AGE/(age&INFECTION) | abs values/fraction
# k=7 --> (free scale, infects separate, absval). k=8 --> (free, infects separate, fractional)
  g(scale_val,facet2tag,value_type,y_axis_tag,ncol_val,facet_formula,foldername,caption_txt,timecourse_filename) %=% 
    fun_tcourse_plottags(k,nval=2,rval=3,n_inf,n_age,colvar,agegr_lim,delta_susc_prop,delta_primary,
                         preseas_npi_on,postseas_npi_off,shutdown_scale,basal_rate)
# PLOT
ggplot(subset(df_ode_solution_tidy,grepl('I',name)&agegroup<=agegr_lim&(t_years>xval_lims[1]&t_years<xval_lims[2])),
 aes_string(x="t_years",y=value_type,group="name",color="infection")) + geom_line(size=1.025) + theme_bw() + standard_theme +
 theme(axis.text.x=element_text(size=8,vjust=0.5),axis.text.y=element_text(size=8),legend.position="top",legend.title=element_blank()) +
 facet_wrap(as.formula(facet_formula),ncol=as.numeric(ncol_val),scales=scale_val) + scale_x_continuous(breaks=xval_breaks) +
 geom_rect(aes(xmin=shutdwn_lims[1]/365,xmax=shutdwn_lims[2]/365,ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
 geom_vline(data=seas_lims_plot,aes(xintercept=on),color="blue",linetype="dashed",size=0.3) +
 geom_vline(data=seas_lims_plot,aes(xintercept=off),color="turquoise4",linetype="dashed",size=0.3) +
 ggtitle('RSV infections by age group') + xlab('years') + ylab(y_axis_tag) + labs(caption=caption_txt)
## SAVE
plot_x_y=c(40,18); if (grepl("inf",colvar)){plot_x_y=plot_x_y[c(2,1)]*c(1.4,1)}
ggsave(timecourse_filename, width=plot_x_y[1],height=plot_x_y[2],units="cm")
}

### map infections to symptomatic cases -----------------------------------
# prop_symptom_age=cbind(c(0,1,5,15),c(0.9,4.9,14.9,99),prop_symptom_age)
list_symptom_agegroups=list(1:2,3:7,8:9,10:11); prop_symptom=1-c(mean(rbeta(1e3,shape1=3,shape2=30)),mean(rbeta(1e3,shape1=9,shape2=43)),
                 mean(rbeta(1e3,shape1=38,shape2=35)),mean(rbeta(1e3,shape1=36,shape2=11))) 
# EXPOSURE DEPENDENT or not? expos_dep_val=0 -> no dependence on exposure. if expos_dep_val>0 -> dependence on expos)
df_symptom_prop=fun_propsymptom_table(list_symptom_agegroups,expos_dep_val=1,rsv_age_groups$agegroup_name,n_inf=3)
# plot (clinical fraction ~ f(age and #infection))
ggplot(df_symptom_prop,aes(x=agegroup_name,y=sympt_value,color=factor(n_inf),linetype=factor(n_inf),group=n_inf)) + 
  geom_line(size=1.25) + theme_bw() + standard_theme + labs(color="# infection",linetype="# infection") + ylab("% symptomatic")
# ggsave("simul_output/severity_age_exp_depend.png",width=22,height=16,units="cm")

# CALCULATE SYMPTOMATIC CASES & sum of 1,2,3rd infections
df_ode_sol_cases_sum=fun_symptomcases_table(df_ode_solution_tidy,df_symptom_prop, bindcolnames=c("infection","agegroup"),n_sd=2)

# PLOT SYMPTOMATIC CASES
for (k_plot in 1:2){ # nrow(permutations(n=2,r=2,repeats.allowed=T))
# free/fixed | fraction/cases
g(scale_val,y_axis_tag,plotvar,foldername,full_filename,caption_txt) %=% fun_sumcase_plot_tags(n_val=2,r_val=2,k_plot,
                                                df_symptom_prop,preseas_npi_on,postseas_npi_off,shutdown_scale,basal_rate)
seas_lims_plot=subset(seas_lims,sd==n_sd & on>xval_lims[1] & off<xval_lims[2])
# plot
ggplot(subset(df_ode_sol_cases_sum,t_years>xval_lims[1] & t_years<(max_time/365-0.05) & agegroup<=9),
       aes(x=t_years,y=get(plotvar))) + geom_line(size=1.02) + theme_bw() + standard_theme + 
  theme(axis.text.x=element_text(size=8,vjust=0.5),axis.text.y=element_text(size=9)) + facet_wrap(~agegroup_name,scales=scale_val,ncol=3) + 
  scale_x_continuous(breaks=xval_breaks,minor_breaks=seq(0,max_time/365,by=1/12)) + # scale_y_continuous(breaks=(0:10)/10) + # 
  geom_rect(aes(xmin=shutdwn_lims[1]/365,xmax=shutdwn_lims[2]/365,ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
  geom_vline(data=seas_lims_plot,aes(xintercept=on),color="blue",linetype="dashed",size=0.3) +
  geom_vline(data=seas_lims_plot,aes(xintercept=off),color="turquoise4",linetype="dashed",size=0.3) +
  xlab('years') + ylab(y_axis_tag) + ggtitle('Symptomatic cases by age group') + theme(legend.position='none') +
  labs(caption=caption_txt) # subtitle=paste0("NPI during year ",npi_year)
# SAVE
# full_filename=gsub('.png','_reducedresurg.png',full_filename)
ggsave(full_filename,width=31,height=22,units="cm")
}

### Age distrib before and after shutdown ----------------------------
season_peaks_AUC=fcn_seas_agedistrib(df_ode_sol_cases_sum,max_time,timesteps,seas_case_threshold=5e2)
# subset for plot
plot_season_agedist=subset(season_peaks_AUC, agegroup<=7 & season>=min(3,npi_year-1))
selvar="auc_case_share_season"; ylimvals=c(min(plot_season_agedist[,selvar]),max(plot_season_agedist[,selvar]))
# if (length(unique(plot_season_agedist$pre_post_shtdn))==3) {linetype_vals=c("dashed","dotted","solid")} else {linetype_vals=c("dotted","solid")}
# PLOT
ggplot(plot_season_agedist,aes(x=agegroup_name,y=get(selvar),group=season)) +
  geom_line(aes(linetype=pre_post_shtdn,color=factor(season),size=pre_post_shtdn)) + 
  geom_point(aes(shape=pre_post_shtdn,color=factor(season)),size=3) + theme_bw() + standard_theme + 
  theme(axis.text=element_text(size=12),legend.text=element_text(size=12),legend.title=element_text(size=14)) +
  xlab("") + ylab("fraction of all symptomatic cases in season") + # facet_wrap(~season,ncol=3) + 
  scale_linetype_manual(values=c("solid","dotdash","dashed"))+scale_shape_manual(values=c(0,2,1)) + 
  scale_size_manual(values=c(1.75,1.25,1.25)) +
  scale_y_continuous(breaks=seq(round(ylimvals[1],2),round(ylimvals[2],2),by=0.02),limits=ylimvals) +
  labs(size="pre/post-NPI",shape="pre/post-NPI",color="season",linetype="pre/post-NPI",title="Age distribution before and after NPI in year 5",
  caption=paste0('NPI (-',preseas_npi_on,',',postseas_npi_off,") weeks from CI95 season, ",basal_rate*1e2,"% off-season activity"))
# SAVE
agedistr_filename=fcn_agedistrib_plot_tags(delta_susc_prop,delta_primary,plot_season_agedist,df_symptom_prop,
                                           preseas_npi_on,postseas_npi_off,seas_case_threshold=5e3)
ggsave(agedistr_filename,width=32,height=22,units="cm")

# mean age of infections
mean_age_perseason=season_peaks_AUC %>% group_by(season) %>% 
  summarise(season_mean_age=sum(auc_case_share_season*mean_age*agegr_size/sum(agegr_size)),pre_post_shtdn=unique(pre_post_shtdn))
meanage_max=ceiling(max(mean_age_perseason$season_mean_age))
# PLOT
ggplot(subset(mean_age_perseason,season>2),aes(x=season,y=season_mean_age,fill=pre_post_shtdn)) + geom_bar(stat="identity",color="black") + 
  scale_y_continuous(breaks=seq(0,meanage_max,0.5)) + scale_x_continuous(breaks=0:max(mean_age_perseason$season)) + 
  geom_rect(aes(xmin=ceiling(shutdwn_lims[1]/365)+0.6,xmax=ceiling(shutdwn_lims[1]/365)+1.4,ymin=0,ymax=meanage_max),
            fill="pink",color=NA,alpha=0.1) + theme_bw() + standard_theme + theme(legend.title=element_text("pre/post NPI"),
  axis.text.x=element_text(size=9,vjust=0.5)) + coord_cartesian(ylim=c(0.2,meanage_max)) + coord_flip() +
  geom_text(aes(label=round(season_mean_age,1)),hjust=-0.3,color="black",size=4) + ylab("mean age of cases") # xlab("season") + 
# SAVE
ggsave(gsub("age_distrib","mean_age",agedistr_filename),width=16,height=11,units="cm")

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