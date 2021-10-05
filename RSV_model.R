# install.packages("devtools"); library("devtools")
# contact data from https://bisaloo.github.io/contactdata/index.html (Prem 2017)
# functions
rm(list=ls()); currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
# library(contactdata); library(fitdistrplus);  library(bbmle); library(Rcpp); library(GillespieSSA)
lapply(c("tidyverse","deSolve","gtools","rstudioapi","wpp2019","plotly","Rcpp","zoo","lubridate","tsibble","pracma","qs","ungeviz"),
         library,character.only=TRUE)
source('fcns/RSV_model_functions.R')
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# SET PARAMETERS --------------------------------------------------------
### ### ### ### ### ### ### ###
# constant parameters
# selected country
country_sel="United Kingdom"
# time resolution (in days)
elem_time_step=1
# population data
standard_age_groups <- fun_cntr_agestr(country_sel,i_year="2020",seq(0,75,5),c(seq(4,74,5),99))
popul_struct=fcn_cntr_fullpop(n_year="2020",country_sel)
# RSV age groups (population data from wpp2019)
# rsv_age_groups<-fun_rsv_agegroups(standard_age_groups,popul_struct,rsv_age_groups_low=c(0,0.5,1,1.5, 2,3,4, 5,10,15, 20),
#                                  rsv_age_group_sizes=c(rep(0.4,4),rep(0.9,3),rep(4,3),79))
rsv_age_groups<-fun_rsv_agegroups(standard_age_groups,popul_struct,rsv_age_groups_low=c(0,0.5,1,1.5, 2,3,4, 5,15, 45, 65),
                                 rsv_age_group_sizes=c(rep(0.4,4),rep(0.9,3), 9, 29, 19, 34))
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
# if UK -> England's contact matrix # check: cm_parameters_SEI3R(cm_uk_locations("UK", 1))$pop[[1]]$matrices 
list_contmatrs=fun_covidm_contactmatrix(country_sel,currentdir_path,cm_path=cm_path) 
# make matrix reciprocal
C_m_polymod=Reduce('+',list_contmatrs) # fun_recipr_contmatr(Reduce('+',list_contmatrs),age_group_sizes=standard_age_groups$values)
# create for our age groups
C_m_merged_nonrecipr=fun_create_red_C_m(C_m_polymod,rsv_age_groups,
      orig_age_groups_duration=standard_age_groups$duration,orig_age_groups_sizes=standard_age_groups$values)
# make it reciprocal for the larger group
C_m=fun_recipr_contmatr(C_m_merged_nonrecipr,age_group_sizes=rsv_age_groups$value)
# bc of reinfections we need to input contact matrix repeatedly
contmatr_rowvector=t(do.call(cbind, lapply(1:nrow(C_m), function(x){diag(C_m[x,]) %*% matrix(1,n_age,n_inf)})))
# build kinetic matrix
# WANING (immunity) terms: R_i_j -> S_min(i+1,n_inf)_j
omega=1/350 # 1/runif(1,60,200)
# RECOVERY
rho=1/7 # 1/rho=rweibull(1, shape=4.1,scale=8.3)
# KINETIC MATRIX (aging terms need to be scaled by duration of age groups!)
K_m=fun_K_m_sirs_multiage(dim_sys,n_age,n_inf,n_compartment,rho,omega,varname_list,rsv_age_groups)
# BIRTH RATE into S_1_1 (Germany 2019: 778e3 births)
birth_rates=matrix(c(713e3/365,rep(0,dim_sys-1)),dim_sys,1)
# DEATHS (2019: 530841 deaths [England and Wales!]) # "uk_death_rate_byage_rsv_agegroups.csv" is for 1000 population!
uk_death_rate=read_csv("data/uk_death_rate_byage_rsv_agegroups.csv")
# we want population to be stationary (at 2019 or 2020 value), so deaths = births
rsv_age_groups <- rsv_age_groups %>% mutate(deaths_per_personyear=uk_death_rate$deaths_per_person/1e3,
                             rect_popul=sum(value)*read_csv("data/final_pop.csv")$final/sum(read_csv("data/final_pop.csv")$final))
death_corr_ratios=(rsv_age_groups$duration/sum(rsv_age_groups$duration))/rsv_age_groups$fraction
# 
init_ratio_birth_death <- sum(death_corr_ratios*uk_death_rate$deaths_per_1000person_peryear*rsv_age_groups$value/1000)/(
  sum(birth_rates)*365)
# sum(birth_rates)*365/sum(uk_death_rate$deaths_per_1000person_peryear*rsv_age_groups$value/1e3)
  # read_csv("data/final_pop.csv")$final/rsv_age_groups$value
death_rates=matrix(unlist(lapply(( (death_corr_ratios/(init_ratio_birth_death*1.241))*
  uk_death_rate$deaths_per_1000person_peryear/(365*1000)), function(x) rep(x,n_compartment*n_inf))),ncol=1)
# g(rsv_age_groups,death_rates) %=% fun_death_rates(rsv_age_groups,uk_death_rate,nage=n_age,ninf=n_inf,dimsys=dim_sys)
# estimated attack rates 
estim_attack_rates <- data.frame(agegroup_name=paste0("age=",rsv_age_groups$agegroup_name,"yr"),
                               median_est=c(rep(65,4),rep(40,4),10,8,5)) %>% mutate(min_est=median_est*0.5,max_est=median_est*1.5)
### ### ### ### ### ### ### ###
# variable parameters  --------------------------------------------------------
# SUSCEPTIBILITY (normalised by age group sizes, for infection terms ~ delta*(I1+I2+...+In)*S_i/N_i)
# agedep_fact determines strength of age dependence, if agedep_fact>1, decreasing susceptibility with age
agedep_fact=1; delta_primary=c(0.08,0.07,0.06)/2.8 # c(0.09,0.07,0.05)/2.5
# c(0.27,0.03,0.01)/5 # c(0.21,0.11,0.01)/5 # c(0.09,0.07,0.05)/1.6 # rep(0.15,3) # rep(0.09,3)
delta_susc <- sapply(1:n_age, function(x) {delta_primary/((agedep_fact^(x-1))*rsv_age_groups$value[x])})
delta_susc_prop <- delta_susc*matrix(rep(rsv_age_groups$value,3),nrow=3,byrow=T)
dep_subfolder_name<-fun_subfld(delta_primary,delta_susc_prop)
### PLOT susceptibility: 
# fcn_plot_suscept_table(fcn_suscept_agedeptable(rsv_age_groups,delta_susc,n_inf)) 
# ggsave(paste0("simul_output/suscept_age_dep/",agedep_fact,"delta_prim",unique(delta_primary),".png"),width=32,height=22,units="cm")
# calculate R0 (at max seasonal forcing=1)
R0_calc_SIRS(C_m,delta_susc_prop,rho,n_inf)
####
# DURATION of SIMULATION
# seasonal forcing (baseline level=1, forcing_strength=2 means 200% above baseline) | npi_reduc_strength: reduction from baseline 
# set seas lims from UK data: peak is weeks 49/50, on/off is 41,11
npi_dates=as.Date(c("2020-03-26","2021-04-01")); seaspeakval=1; seasforc_width_wks=4
g(n_years,timesteps,simul_start_end,forcing_vector_npi) %=% fun_shutdown_seasforc(npi_dates,years_pre_post_npi=c(10,50),
        season_width_wks=seasforc_width_wks,init_mt_day="06-01",peak_week=44,forcing_above_baseline=seaspeakval,npireduc_strength=0.5)
# R0, tags
R0val=round(R0_calc_SIRS(C_m,delta_susc_prop,rho,n_inf),2); filetag <- ifelse(grepl("exp_dep",dep_subfolder_name),
  paste0("susc_",paste0(round(delta_primary,2),collapse="_"),"_R0_",R0val,"_seasforc_",seaspeakval,"_seaswidth_wks",seasforc_width_wks),"")
# plot seasonal forcing
fcn_plot_seas_forc(simul_start_end,forcing_vector_npi,seas_lims_wks=c(7,42),npi_dates,date_resol="3 month")
# SAVE: ggsave(paste0("simul_output/NPI_y",npi_year,"_on",preseas_npi_on,"w_off",postseas_npi_off,"w.png"),units="cm",height=10,width=20)

# INITIAL CONDITIONS. # introduce stationary state as init state? | init_set: "previous" or anything else
if (!exists("ode_solution")) {ode_solution<-NA}
initvals_sirs_model <- fcn_set_initconds(rsv_age_groups,init_set=c("previous","fromscratch")[1],init_cond_src=c("output","file")[1],
    ode_solution,init_seed=10,seed_vars="all",filename="simul_output/df_ode_solution_UK_multiple_adult_grps.RDS")
# OR manually choose a timepoint from prev simul
# t_sel=(df_cases_infs %>% filter(date==ymd("2019-07-01")) %>% group_by(date) %>% summarise(t=unique(t)))$t
# initvals_sirs_model=matrix(as.numeric(round(ode_solution[t_sel,2:ncol(ode_solution)])))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### integrate ODE --------------------------------------------------------
params<-list(list(birth_rates,death_rates),K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds,delta_susc)
# interpolation fcns for seas forcing & extern introds
approx_seas_forc <- approxfun(data.frame(t=timesteps,seas_force=forcing_vector_npi))
# how many introductions every 30 days?
approx_introd <- approxfun(data.frame(t=timesteps,as.numeric(timesteps %% 30==0)*5))
tm<-proc.time(); ode_sol<-lsoda(initvals_sirs_model,timesteps,func=sirs_seasonal_forc,parms=params); round(proc.time()-tm,2)
# reshape data | # check size of objs: fcn_objs_mem_use(1)
g(final_pop,ode_solution,df_cases_infs) %=% fun_process_simul_output(ode_sol,varname_list,incidvar="newinf",incid_only=F,
                                        init_date=simul_start_end[1],n_age,n_inf,rsv_age_groups,neg_thresh=-0.01)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# PLOT how age group totals change (for this "incid_only" should be FALSE)
fcn_plotagegroup_totals(df_cases_infs,incidvar="newinf",scale_val=c('free_y','fixed')[1])
# final popul:
# sum((df_cases_infs %>% filter(!grepl("newinf",name) & t==max(t)))$value)
# ggsave(paste("simul_output/agegroup_totals_",scale_val,"yscale.png",sep=""),width=28,height=16,units="cm")
# check if % change between initial and final age group sizes smaller than 0.1 | fun_agegroup_init_final_pop(df_cases_infs)
# all(abs(1-(final_pop$final/final_pop$init))<0.001)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# simple lineplot
sel_weeks <- df_cases_infs %>% mutate(week=week(date),year=year(date)) %>% filter(week %in% c(9,41,49)) %>% group_by(year,agegroup,week) %>%
  filter(date==min(date) & infection=="infection #1"); varname=c("value_fract","value")[2]; sel_age_lim=9
# plot
df_cases_infs %>% filter(t %% 7==0 & agegroup<=sel_age_lim & date>as.Date("2015-07-01") & date<as.Date("2022-04-01")) %>% # 
  ggplot() + geom_area(aes(x=date,y=get(varname)*ifelse(grepl("fract",varname),1e2,1),fill=infection),position=position_stack(reverse=T),
            alpha=0.3,color="black",size=0.25) + facet_wrap(~agegroup_name,scales="free_y") + xlab("") + 
 geom_vline(data=sel_weeks%>%filter(agegroup<=sel_age_lim),aes(xintercept=date,linetype=ifelse(week==49,"solid","dashed")),
    size=1/4,show.legend=F) + scale_x_date(date_breaks="3 months",expand=expansion(0.01,0)) +
  geom_rect(xmin=npi_dates[1],xmax=npi_dates[2],ymin=-Inf,ymax=Inf,fill="pink",alpha=0.01) + scale_y_continuous(expand=expansion(0.01,0)) +
  labs(caption=paste0("seas peak=",max(forcing_vector_npi),"x baseline, NPI contact red.=-",round((1-min(forcing_vector_npi))*1e2),
     "%, R0 (baseline)=",R0val)) + ggtitle("",subtitle=paste0("suscept=[",paste0(round(delta_primary,3),collapse=","),"]") ) +
  theme_bw() + standard_theme + ylab(ifelse(varname=="value","# cases","% age group"))
# save
ggsave(paste0("simul_output/",dep_subfolder_name,"timecourse_",filetag,".png"),width=30,height=15,units="cm")

# phase plot
# df_cases_infs %>% filter(agegroup<=2 & grepl("1|2",infection) & date>as.Date("2007-07-01") & date<as.Date("2019-07-02")) %>% #
#   ungroup() %>% mutate(infection=paste0("inf",as.numeric(factor(infection))),year=year(date)) %>% select(t,value,infection,agegroup,year) %>% 
#   pivot_wider(names_from=c(infection,agegroup)) %>% ggplot() + geom_path(aes(x=inf1_1,y=inf1_2,color=factor(year)),size=1.2) + 
#   labs(caption=paste0("seas peak=",max(forcing_vector_npi),"x baseline, NPI contact red.=-",round((1-min(forcing_vector_npi))*1e2),
#       "%, R0 (baseline)=",R0val)) + ggtitle("",subtitle=paste0("suscept=[",paste0(round(delta_primary,3),collapse=","),"]") ) +
#   theme_bw() + standard_theme
# ggsave(paste0("simul_output/",dep_subfolder_name,"phaseplot_",filetag,".png"),width=30,height=15,units="cm")

# calculate attack rates by epi-year and age group # 81-97% of cases should be within weeks 41 and 9
sum_inf_epiyear_age=left_join(df_cases_infs %>% mutate(year=year(date),epi_year=ifelse(date>ymd(paste(year(date),"-07-01")),
      year(date),year(date)-1),in_out_season=ifelse(week(date)<=9 | week(date)>=41,"in","out")) %>% 
      group_by(epi_year,agegroup_name,agegroup) %>% summarise(inf_tot=sum(value),inf_in_seas=sum(value[in_out_season=="in"])) %>% 
  group_by(agegroup) %>% filter(epi_year>min(epi_year)),final_pop,by="agegroup") %>% 
  mutate(attack_rate_perc=100*inf_tot/final,seas_share=inf_in_seas/inf_tot)
# plot attack rates
ggplot(sum_inf_epiyear_age %>% mutate(attack_rate_perc=ifelse(epi_year==2020,NA,attack_rate_perc))) + 
         geom_hpline(aes(x=factor(epi_year),y=attack_rate_perc),width=0.9,size=1) + 
 facet_wrap(~agegroup_name,scales="free_y") + theme_bw() + standard_theme + xlab("epi_year") + ylab("attack rate % age group") +
  geom_rect(xmin=which(unique(sum_inf_epiyear_age$epi_year)==2020)-0.5,xmax=which(unique(sum_inf_epiyear_age$epi_year)==2020)+1/2,
      ymin=-Inf,ymax=Inf,fill="pink",alpha=0.1) + 
  geom_hline(data=estim_attack_rates %>% pivot_longer(!agegroup_name) %>% filter(name!="median_est"),
    aes(yintercept=value),size=1/3,color="red") + labs(caption=paste0("seas peak=",max(forcing_vector_npi),"x baseline, NPI contact red.=-",
    round((1-min(forcing_vector_npi))*1e2),"%, R0 (baseline)=",R0val)) + ggtitle("attack rates by age group",
    subtitle=paste0("suscept=[",paste0(round(delta_primary,3),collapse=","),"]"))
# save
ggsave(paste0("simul_output/",dep_subfolder_name,"attack_rates_",filetag,".png"),width=20,height=15,units="cm")

# plot concentration of infections within seasons
ggplot(sum_inf_epiyear_age %>% mutate(seas_share=ifelse(epi_year==2020,NA,seas_share))) + 
  geom_hpline(aes(x=factor(epi_year),y=seas_share*100),width=0.9,size=1) + facet_wrap(~agegroup_name,scales="free") + 
  geom_rect(xmin=which(unique(sum_inf_epiyear_age$epi_year)==2020)-0.5,xmax=which(unique(sum_inf_epiyear_age$epi_year)==2020)+1/2,
            ymin=-Inf,ymax=Inf,fill="pink",alpha=0.1) + theme_bw() + standard_theme + xlab("epiyear") + ylab("infections within season") +
  labs(caption=paste0("seas peak=",max(forcing_vector_npi),"x baseline, NPI contact red.=-",
    round((1-min(forcing_vector_npi))*1e2),"%, R0 (baseline)=",round(R0_calc_SIRS(C_m,delta_susc_prop,rho,n_inf),2) )) + 
  ggtitle("% infections within season",subtitle=paste0("suscept=[",paste0(round(delta_primary,3),collapse=","),"]")) + 
  geom_hline(yintercept=c(80,100),size=1/3,color="red")
# save
ggsave(paste0("simul_output/",dep_subfolder_name,"seas_share_infs.png"),width=20,height=15,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Plot time course --------------------------------------------------------
xval_lims=c(npi_year-2.23,npi_year+4); xval_breaks=seq(0,max_time/365,by=1/4)
# PLOT
# tags: y-axis fixed/free | facet by AGE/(age&INFECTION) | abs values/fraction
# k=7 --> (free scale, infects separate, absval). k=8 --> (free, infects separate, fraction)
g(scale_val,facet2tag,value_type,y_axis_tag,ncol_val,facet_formula,foldername,caption_txt,subtitle_str,timecourse_filename) %=% 
    fun_tcourse_plottags(k=8,nval=2,rval=3,n_inf,n_age,colvar="age",agegr_lim=7,delta_susc_prop,delta_primary,agedep_fact,
          preseas_npi_on,postseas_npi_off,npi_reduc_strength,forcing_strength)
# PLOT time course absolute/fract values by age group  --------------------------------------------------------
fcn_plot_allcases_absval_stackedinfs(df_cases_infs,by_date=T,sel_var="newinf",value_type,x_lims=xval_lims,t_subset=7,agegrlim=9,ncolval=3,
    y_axis_tag,scale_val,vertline_x=subset(seas_lims,on>xval_lims[1]&off<xval_lims[2]),shutdwn_lims,xval_breaks,subtitle_str,caption_txt)
## SAVE
ggsave(timecourse_filename, width=32,height=20,units="cm")

### % 1st infections (of all infections) by season -----------------
# for the whole 'epi-year': 0.78-1.78 = epiyear1 # for seasons only: 0.78-1.21 = season 1
inf_fracts_season <- df_cases_infs %>% filter(compartment=="newinf") %>% group_by(t,agegroup,agegroup_name) %>% 
  mutate(value_fract_infs=value/sum(value)) %>% mutate(epiyear=ceiling(t/365-unique(round(seas_lims$on-floor(seas_lims$on),3))),#+1
    season=ifelse(t/365-floor(t/365)>unique(round(seas_lims$on-floor(seas_lims$on),3)) | 
                    t/365-floor(t/365)<unique(round(seas_lims$off-floor(seas_lims$off),3)),"in","off"),
        infection_bin=ifelse(grepl("#1",infection),"first","reinfection") ) %>% group_by(epiyear,infection_bin,agegroup_name) %>% 
  summarise(max_fract_epiyear=max(value_fract_infs),mean_fract_epiyear=mean(value_fract_infs),
            max_fract_season=max(value_fract_infs[season=="in"]),mean_fract_season=mean(value_fract_infs[season=="in"]))
# plot
# mean_fract_epiyear max_fract_epiyear # max_fract_season mean_fract_season
plot_xlims=c(npi_year-2,npi_year+4); selval="max_fract_season" 
ggplot(subset(inf_fracts_season,epiyear>plot_xlims[1]&epiyear<=plot_xlims[2]&as.numeric(agegroup_name)<=6 &infection_bin=="first"))+
  geom_hpline(aes(x=epiyear,y=get(selval)*1e2,color=ifelse(epiyear==npi_year+2,"red",NA)),width=0.9,show.legend=F) +
  facet_wrap(~agegroup_name,scales="free") + theme_bw() + standard_theme + # labs(color="") +  # 
  geom_rect(aes(xmin=npi_year+0.5,xmax=npi_year+1.5,ymin=-Inf,ymax=Inf),fill="pink",alpha=0.2) + 
  geom_vline(xintercept=(plot_xlims[1]:plot_xlims[2])+0.5,size=0.2,linetype="dashed",color="black") + 
  theme(panel.grid.major.x=element_blank(),axis.text.x=element_text(vjust=0.5)) + 
  xlab("RSV season (epi-year)") + ylab("% first infections") + labs(subtitle=subtitle_str,caption=caption_txt)
# SAVE
ggsave(paste0(paste0(strsplit(timecourse_filename,"/")[[1]][1:2],collapse = "/"), "/",paste0("proportion_firstinf_",selval),
              ifelse(agedep_fact>1,paste0("_agedep_",agedep_fact),paste0("_expdep_",paste0(delta_primary,collapse="_")) ),
              "_peakforc",(1+forcing_strength)*1e2,"_NPIred",npi_reduc_strength*1e2,".png"), width=32,height=18,units="cm")

### mean age of 1st/2nd/3rd infection --------------------------------------------------------
mean_age_byinf_at_t=fcn_calc_mean_age_byinfs(rsv_age_groups,df_cases_infs,sel_var="newinf",seas_lims,low_thresh=1e3,n_aver=30)
# mean age per season
season_means_byinf=subset(mean_age_byinf_at_t, grepl("in",on_off)) %>% mutate(ageweight_aver=mean_age_at_t*suminfs) %>% 
  group_by(epi_year,infection) %>% summarise(weighted_sum=ifelse(sum(suminfs)>5,sum(ageweight_aver)/sum(suminfs),NA))
# check timecourse
# ggplot(subset(x,t>10.7*365 & t<13.3*365 & agegroup<=6),aes(x=t/365,y=value_fract,color=infection)) +
# geom_line() + facet_wrap(~agegroup_name) + theme_bw()
# plot means per season by infection type
plot_xlims=c(npi_year-2,npi_year+3)
ggplot(subset(season_means_byinf %>% mutate(weighted_sum=ifelse(epi_year==npi_year+1,NA,weighted_sum)),epi_year>plot_xlims[1]&epi_year<=plot_xlims[2])) + 
  geom_hpline(aes(x=epi_year,y=weighted_sum,color=ifelse(epi_year==npi_year+2,"red",NA)),width=0.9,show.legend=F) +
  facet_wrap(~infection,scales="free") + theme_bw() + standard_theme + theme(panel.grid.major.x=element_blank()) +
  geom_rect(aes(xmin=npi_year+0.5,xmax=npi_year+1.5,ymin=-Inf,ymax=Inf),fill="pink",alpha=0.2) + xlab("RSV season")+ylab("mean age") +
  geom_vline(xintercept=(plot_xlims[1]:plot_xlims[2])+0.5,size=0.2,linetype="dashed",color="black") + 
  theme(axis.text.x=element_text(vjust=0.5)) + labs(subtitle=subtitle_str,caption=caption_txt)
# save
ggsave(paste0(paste0(strsplit(timecourse_filename,"/")[[1]][1:2],collapse = "/"),
        "/meanage_perseas_byinf",ifelse(agedep_fact>1,paste0("_agedep_",agedep_fact),paste0("_expdep_",paste0(delta_primary,collapse="_")) ),
        "_peakforc",(1+forcing_strength)*1e2,"_NPIred",npi_reduc_strength*1e2,".png"), width=32,height=18,units="cm")

# PLOT time course normalised by max per age group  --------------------------------------------------------
xval_lims=c(npi_year-1.27,npi_year+3.25); seas_lims_plot=subset(seas_lims,on>xval_lims[1] & off<xval_lims[2])
fcn_plot_norm_max_byage(df_cases_infs,xval_lims,seas_lims_plot,shutdwn_lims,xval_breaks,subtitle_str,caption_txt)
# SAVE
ggsave(gsub(".png","_byage_maxnorm.png",timecourse_filename), width=30,height=20,units="cm")

### PLOT time course faceted by infection --------------------------------------------------------
xval_lims=c(npi_year-1.27,npi_year+3.25); seas_lims_plot=subset(seas_lims,on>xval_lims[1] & off<xval_lims[2])
fcn_plot_norm_max_byage(df_cases_infs,xval_lims,seas_lims_plot,shutdwn_lims,xval_breaks,subtitle_str,caption_txt)
# SAVE
ggsave(gsub(".png","_byinf_maxnorm.png",timecourse_filename), width=30,height=20,units="cm")

### distribution of 1st/2nd/3rd infections faceted by age --------------------------------------------------------
xval_lims=c(npi_year-1.27,npi_year+3.25); seas_lims_plot=subset(seas_lims,on>xval_lims[1] & off<xval_lims[2])
fcn_plot_share_infs_agefaceted(df_cases_infs,xval_lims,seas_lims_plot,n_aver=40,agegrlim=8,ncolval=4,
                               shutdwn_lims,xval_breaks,subtitle_str,caption_txt)
# SAVE
ggsave(gsub(".png","_byinf_share.png",timecourse_filename), width=32,height=18,units="cm")

### age distribution of infections faceted by 1st/2nd/3rd infections --------------------------------------------------------
xval_lims=c(npi_year-2.27,npi_year+3.25); seas_lims_plot=subset(seas_lims,on>xval_lims[1] & off<xval_lims[2])
# (df,x_lims,vertline_x,highl_lims,xvalbreaks,n_aver,agegrlim,ncolval,subtitle_str,caption_txt)
fcn_plot_share_infs_inf_faceted(df_cases_infs,xval_lims,seas_lims_plot,shutdwn_lims,xval_breaks,n_aver=35,agegrlim=n_age,ncolval=1,
                                subtitle_str,caption_txt)
# SAVE
# ggsave(gsub(".png","_byinf_share.png",timecourse_filename), width=32,height=18,units="cm")

### plot timecourse of mean age ----------------------------
xval_lims=c(npi_year-2.27,npi_year+2.5); seas_lims_plot=subset(seas_lims,on>xval_lims[1] & off<xval_lims[2])
# mean_age_at_t_smooth | season_mean. if plot_flag="seas" -> per season average
fcn_plot_mean_age_byinfs(subset(mean_age_byinf_at_t,grepl("#1",infection)),
                         plot_flag="",plotvar="mean_age_at_t_smooth",low_thresh=1e3,xval_lims,seas_lims_plot,
                         shutdwn_lims,xval_breaks,ncolval=1,shutdwn_lims,subtitle_str,caption_txt)
# SAVE
# ggsave(gsub(".png","_meanagebyinf_t.png",timecourse_filename), width=32,height=18,units="cm")
# ggsave(gsub(".png","_meanageperseasbyinf.png",timecourse_filename), width=32,height=18,units="cm")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### CALCULATE SYMPTOMATIC CASES & sum of 1,2,3rd infections -------------------------------
# dependence of SEVERITY on age (>1) | dependence on prev exposure
sever_agedep_fact=1; clinfract_expos=rep(1,n_inf) # c(0.8,0.4,0.2) # rep(0.75,3)
clin_fract_age_exp=fcn_clin_fract_age_exp(sever_agedep_fact,clinfract_expos,rsv_age_groups,n_age,n_inf)
# plot: fcn_plot_clinfract_table(clin_fract_age_exp)
df_sympt_cases=fun_symptomcases_table(df_cases_infs,clin_fract_age_exp,c("infection","agegroup"))

# PLOT time course SYMPTOMATIC CASES -----------------------------------
# k_plot: free/fixed | fraction/cases
g(scale_val,y_axis_tag,plotvar,foldername,symptcases_filename,caption_txt,subtitle_str) %=% fun_sumcase_plot_tags(n_val=2,r_val=2,
      k_plot=1,clin_fract_age_exp,delta_susc_prop,delta_primary,preseas_npi_on,postseas_npi_off,npi_reduc_strength,forcing_strength)
xval_lims=c(npi_year-2.5,npi_year+3.5); seas_lims_plot=subset(seas_lims,on>xval_lims[1]&off<xval_lims[2]) %>% pivot_longer(col=!season)
# plot
fcn_plot_symptomcases_agefacet(df_sympt_cases %>% filter(agegroup<=6),xval_lims,y_plotvar=plotvar,y_axis_tag,scale_val,seas_lims_plot,
                               shutdwn_lims,xval_breaks,n_col=3,n_agelim=11,subtitle_str,caption_txt)
# SAVE
# full_filename=gsub('.png','_reducedresurg.png',full_filename)
ggsave(symptcases_filename,width=31,height=22,units="cm")

### mean age per season ---
low_thr=1e3
mean_age_symptcases_at_t <- left_join(rsv_age_groups %>% mutate(agegroup=as.numeric(factor(agegroup_name,levels=agegroup_name))) %>% select(agegroup,mean_age_weighted),
            subset(df_sympt_cases,compartment %in% "I"),by="agegroup") %>% group_by(t) %>% # t,
  summarise(suminfs=sum(symptom_cases),mean_age_at_t=sum(symptom_cases*mean_age_weighted/sum(symptom_cases))) %>% 
  mutate(mean_age_at_t_smooth=ifelse(suminfs>low_thr,rollmean(mean_age_at_t,k=30,align="center",fill=NA),NA)) %>%
  mutate(on_off=ifelse(mod(findInterval(t/365,array(matrix(t(seas_lims[,c("on","off")])))),2)==1 & t/365>=seas_lims$on[1],
                       "in-season","off-season"),epi_year=findInterval(t/365,seas_lims$on)) %>% group_by(epi_year,on_off) %>% # ,infection
  mutate(season_mean=ifelse(mean(suminfs)>low_thr,mean(mean_age_at_t),NA))
# mean age per season
season_means_symptcases=subset(mean_age_symptcases_at_t, grepl("in",on_off)) %>% mutate(ageweight_aver=mean_age_at_t*suminfs) %>% 
  group_by(epi_year) %>% summarise(weighted_sum=ifelse(sum(suminfs)>5e4,sum(ageweight_aver)/sum(suminfs),NA))
# PLOT
plot_xlims=c(npi_year-2,npi_year+4)
ggplot(subset(season_means_symptcases %>% mutate(weighted_sum=ifelse(epi_year==npi_year+1,NA,weighted_sum)),epi_year>plot_xlims[1]&epi_year<=plot_xlims[2])) + 
  geom_hpline(aes(x=epi_year,y=weighted_sum,color=ifelse(epi_year==npi_year+2,"red",NA)),width=0.9,show.legend=F) +
  theme_bw() + standard_theme + theme(panel.grid.major.x=element_blank()) + scale_x_continuous(breaks=1:round(n_years)) +
  geom_rect(aes(xmin=npi_year+0.5,xmax=npi_year+1.5,ymin=-Inf,ymax=Inf),fill="pink",alpha=0.2) + xlab("RSV season")+ylab("mean age") +
  geom_vline(xintercept=(plot_xlims[1]:plot_xlims[2])+0.5,size=0.2,linetype="dashed",color="black") + 
  theme(axis.text.x=element_text(vjust=0.5,size=13),axis.text.y=element_text(size=13),axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15)) + labs(subtitle=subtitle_str,caption=caption_txt)
# SAVE
ggsave(paste0(paste0(unlist(strsplit(symptcases_filename,"/"))[1:3],collapse="/"),"/mean_age_symptcases_peak_forc",
        (1+forcing_strength)*1e2,"_NPIred",npi_reduc_strength*1e2,"pct",
        ifelse(agedep_fact>1,paste0("_agedep",agedep_fact),""),"_delta",paste0(unique(delta_primary),collapse="_"),".png"),
       width=22,height=14,units="cm")

### plot symptom cases stacked as area plot -----------
xval_lims=c(npi_year-0.5,npi_year+4.5); seas_lims_plot=subset(seas_lims,on>xval_lims[1]&off<xval_lims[2]) %>% pivot_longer(col=!season)
fcn_plot_symptomcases_agestacked(df_sympt_cases,xval_lims,plotvar,y_axis_tag,seas_lims_plot,shutdwn_lims,
                                 xval_breaks,subtitle_str,caption_txt)
# SAVE
ggsave(gsub("symptomcases","symptomcases_areaplot",symptcases_filename),width=31,height=22,units="cm")

### plot fraction symptom cases BY AGE as stacked area plot -----------
fcn_plot_share_symptcases_byage(df_sympt_cases,n_aver=28,xval_lims,seas_lims_plot,shutdwn_lims,xval_breaks,subtitle_str,caption_txt)
# fcn_plot_share_symptcases_byage(df_plot,n_aver,x_lims,vertline_x,highl_lims,xvalbreaks,subtitle_str,caption_txt)
  
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Age distrib before and after shutdown ----------------------------
season_peaks_AUC=fcn_seas_agedistrib(df_sympt_cases,max_time,timesteps,seas_case_threshold=1e3) %>%
  pivot_longer(cols=!c(agegroup,agegroup_name,season,agegr_size,mean_age,pre_post_shtdn))
# PLOT
plot_season_agedist=fcn_agedistrib_calc_season(season_peaks_AUC,seaslims=c(npi_year-1,npi_year+4),selvar="auc_case",
                                                 agegr_merge_min=which(rsv_age_groups$age_low==5)-1,rsv_age_groups)
if (length(unique(plot_season_agedist$pre_post_shtdn))==3) {colorvals=c(2,1,3)} else {colorvals=c(2,3)}
# bar plot
fcn_plot_agedistrib_perseas(plot_season_agedist,nrowval=2,yexpval=c(0,0.02),xval_lims,subtitle_str,caption_txt,textsize=3.5)
# SAVE
agedistr_filename=fcn_agedistrib_plot_tags(delta_susc_prop,delta_primary,plot_season_agedist,clin_fract_age_exp,
                                           preseas_npi_on,postseas_npi_off,npi_reduc_strength,seas_case_threshold=1e2)
ggsave(agedistr_filename,width=32,height=16,units="cm")

### ### ### ### ### 
### mean age of infections ----------------------------
mean_age_perseason=fcn_calc_mean_age(season_peaks_AUC,season_min=4,dep_name = "age")
# segment plot
var_sel="value"; if(grepl("norm",var_sel)) {ylab_tag=" (normalised)"} else {ylab_tag=" (year)"}
fcn_plot_mean_age(mean_age_perseason,seas_lims,npi_year,shutdwn_lims,yexpand=c(0.2,0),ylab_tag)
# SAVE
filenm=gsub("pct","pct_geomsegment",gsub("age_distrib","mean_age",agedistr_filename))