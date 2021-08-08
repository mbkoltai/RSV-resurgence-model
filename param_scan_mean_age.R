# load parameters
source("load_params.R")
# parallelised
library(doParallel)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# compare shifts in aver age of infection: age-dependent vs exposure-dependent
list_subtitle=list(); r0_vals=array(); suscept_vals=array()
for (k_dep in 1:2){
  for (k in 1:10) {
    if (k_dep==1){
      agedep_fact=1+k/10; delta_primary=rep(0.45*k/10,3) # 1,c(0.09,0.07,0.05)/3.8 | 1.2, rep(0.09,3)
      r0_try=R0_calc_SIRS(C_m,sapply(1:n_age, function(x) 
        delta_primary/((agedep_fact^(x-1))*rsv_age_groups$value[x]))*matrix(rep(rsv_age_groups$value,3),nrow=3,byrow=T),rho,n_inf)
      while (r0_try>1.5 | r0_try<1.4) {
        delta_primary=delta_primary*ifelse(r0_try>1.5,0.99,1.01)
        r0_try=R0_calc_SIRS(C_m,sapply(1:n_age, function(x) 
          delta_primary/((agedep_fact^(x-1))*rsv_age_groups$value[x]))*matrix(rep(rsv_age_groups$value,3),nrow=3,byrow=T),rho,n_inf) }
    } else {
      agedep_fact=1; delta_primary=c(0.11+k/100,0.11,0.11-k/100)
      r0_try=R0_calc_SIRS(C_m,sapply(1:n_age, function(x) 
        delta_primary/((agedep_fact^(x-1))*rsv_age_groups$value[x]))*matrix(rep(rsv_age_groups$value,3),nrow=3,byrow=T),rho,n_inf)
      while (r0_try>1.5 | r0_try<1.4) {
        delta_primary=delta_primary*ifelse(r0_try>1.5,0.99,1.01); r0_try=R0_calc_SIRS(C_m,sapply(1:n_age, function(x) 
          delta_primary/((agedep_fact^(x-1))*rsv_age_groups$value[x]))*matrix(rep(rsv_age_groups$value,3),nrow=3,byrow=T),rho,n_inf)
        r0_try=R0_calc_SIRS(C_m,sapply(1:n_age, function(x) 
          delta_primary/((agedep_fact^(x-1))*rsv_age_groups$value[x]))*matrix(rep(rsv_age_groups$value,3),nrow=3,byrow=T),rho,n_inf)
      }
    }
    # assign
    delta_susc=sapply(1:n_age, function(x) {delta_primary/((agedep_fact^(x-1))*rsv_age_groups$value[x])})
    delta_susc_prop=delta_susc*matrix(rep(rsv_age_groups$value,3),nrow=3,byrow=T)
    r0_vals[k]=R0_calc_SIRS(C_m,delta_susc_prop,rho,n_inf); suscept_vals[k]=paste0(unique(delta_primary),collapse=",")
    print(paste0(k_dep,",",k,", R0=",round(r0_vals[k],2),", suscept=",paste0(round(delta_primary,3),collapse = ",")) )

    # set params
    params=list(list(birth_rates,death_rates),K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds,elem_time_step,delta_susc)
    # interpolation fcns for seas forcing & extern introds
    # approx_seas_forc <- approxfun(data.frame(t=timesteps,seas_force=forcing_vector_npi))
    # approx_introd <- approxfun(data.frame(t=timesteps,as.numeric(timesteps %% 30==0)*10))
    tm<-proc.time(); ode_solution <- lsoda(initvals_sirs_model,timesteps,func=sirs_seasonal_forc,parms=params); round(proc.time()-tm,2)
    # reshape data | # check size of objs: fcn_objs_mem_use(1)
    g(ode_solution,df_cases_infs) %=% fun_process_simul_output(ode_solution,varname_list,n_age,n_inf,rsv_age_groups,neg_thresh=-0.01)
    
    # labels for plots
    g(scale_val,facet2tag,value_type,y_axis_tag,ncol_val,facet_formula,foldername,caption_txt,subtitle_str,timecourse_filename) %=% 
      fun_tcourse_plottags(k=8,nval=2,rval=3,n_inf,n_age,colvar="age",agegr_lim=7,delta_susc_prop,delta_primary,agedep_fact,
                           preseas_npi_on,postseas_npi_off,npi_reduc_strength,forcing_strength)
    # fun_tcourse_plottags(k=8,nval=2,rval=3,n_inf,n_age,colvar="age",agegr_lim=7,delta_susc_prop,agedep_fact,delta_primary,
    # preseas_npi_on,postseas_npi_off,npi_reduc_strength,forcing_strength)
    if (!grepl("expdep|agedep",subtitle_str)){
      list_subtitle[[k]]=c(subtitle=paste0(subtitle_str, ifelse(agedep_fact==1,"",paste0(", agedep=",agedep_fact)), 
                    ifelse(length(unique(delta_primary))==1,"",paste0(", expdep=[",paste0(round(delta_primary,3),collapse = ","),"]"))),
                           caption=gsub("from CI95 ","wrt ",caption_txt) ) }
    # PLOT time course absolute/fract values by age group  --------------------------------------------------------
    # fcn_plot_allcases_absval_stackedinfs(df_cases_infs,value_type,x_lims=xval_lims,t_subset=7,agegrlim=6,ncolval=3,y_axis_tag,
    #         scale_val,vertline_x=subset(seas_lims,on>xval_lims[1]&off<xval_lims[2]),shutdwn_lims,xval_breaks,subtitle_str,caption_txt)
    dyn=df_cases_infs %>% filter(compartment=="I") %>% mutate(type=ifelse(agedep_fact>1,"age","immunity"),dep_strength=k )
    
    # get peaks (or mean) of % of 1st infections by season
    seas_lims_uniq=c(on=unique(round(seas_lims$on-floor(seas_lims$on),3)),off=unique(round(seas_lims$off-floor(seas_lims$off),3)))
    x = df_cases_infs %>% filter(compartment=="I") %>% group_by(t,agegroup,agegroup_name) %>% 
      mutate(value_fract_infs=value/sum(value)) %>% mutate(epiyear=ceiling(t/365-as.numeric(seas_lims_uniq["on"])),
            season=ifelse(t/365-floor(t/365)> seas_lims_uniq["on"] | t/365-floor(t/365)<seas_lims_uniq["off"],"in","off"),
            infection_bin=ifelse(grepl("#1",infection),"first","reinfection") ) %>% group_by(epiyear,infection_bin,agegroup_name) %>% 
      summarise(max_fract_epiyear=max(value_fract_infs),mean_fract_epiyear=mean(value_fract_infs),
                max_fract_season=max(value_fract_infs[season=="in"]),mean_fract_season=mean(value_fract_infs[season=="in"])) %>% 
      mutate(type=ifelse(agedep_fact>1,"age","immunity"),agedepfact=agedep_fact,
             suscept=paste0(round(unique(delta_primary),3),collapse=","),dep_strength=k )
    if (k==1 & k_dep==1) { inf_fracts_season=x; dyn_parscan=dyn } else {
      inf_fracts_season=bind_rows(inf_fracts_season,x); dyn_parscan=bind_rows(dyn_parscan,dyn)}
    
    ### mean age of 1st/2nd/3rd infection --------------------------------------------------------
    mean_age_byinf_at_t=fcn_calc_mean_age_byinfs(rsv_age_groups,df_cases_infs,seas_lims,low_thresh=8e2,n_aver=30)
    # mean age per season
    y = subset(mean_age_byinf_at_t, grepl("in",on_off)) %>% mutate(ageweight_aver=mean_age_at_t*suminfs) %>% 
      group_by(epi_year,infection) %>% summarise(weighted_sum=ifelse(sum(suminfs)>1e4,sum(ageweight_aver)/sum(suminfs),NA)) %>%
      mutate(type=ifelse(agedep_fact>1,"age","immunity"),agedepfact=agedep_fact,
             suscept=paste0(round(unique(delta_primary),3),collapse=","),dep_strength=k)
    if (k==1 & k_dep==1){season_means_byinf=y} else {season_means_byinf=bind_rows(season_means_byinf,y)}
    
    # average age of all infections
    sever_agedep_fact=1; clinfract_expos=rep(1,n_inf) # c(0.8,0.4,0.2) # rep(0.75,3)
    clin_fract_age_exp=fcn_clin_fract_age_exp(sever_agedep_fact,clinfract_expos,rsv_age_groups,n_age,n_inf)
    df_sympt_cases=fun_symptomcases_table(df_cases_infs,clin_fract_age_exp,c("infection","agegroup"))
    low_thr=1e3
    mean_age_symptcases_at_t <- left_join(rsv_age_groups %>% mutate(agegroup=as.numeric(factor(agegroup_name,levels=agegroup_name))) %>% 
              select(agegroup,mean_age_weighted),subset(df_sympt_cases,compartment %in% "I"),by="agegroup") %>% group_by(t) %>% # t,
      summarise(suminfs=sum(symptom_cases),mean_age_at_t=sum(symptom_cases*mean_age_weighted/sum(symptom_cases))) %>% 
      mutate(mean_age_at_t_smooth=ifelse(suminfs>low_thr,rollmean(mean_age_at_t,k=30,align="center",fill=NA),NA)) %>%
      mutate(on_off=ifelse(mod(findInterval(t/365,array(matrix(t(seas_lims[,c("on","off")])))),2)==1 & t/365>=seas_lims$on[1],
                  "in-season","off-season"),epi_year=findInterval(t/365,seas_lims$on)) %>% group_by(epi_year,on_off) %>% # ,infection
      mutate(season_mean=ifelse(mean(suminfs)>low_thr,mean(mean_age_at_t),NA))
    # mean age per season
    zz=subset(mean_age_symptcases_at_t, grepl("in",on_off)) %>% mutate(ageweight_aver=mean_age_at_t*suminfs) %>% 
      group_by(epi_year) %>% summarise(weighted_sum=ifelse(sum(suminfs)>5e4,sum(ageweight_aver)/sum(suminfs),NA)) %>%
      mutate(type=ifelse(agedep_fact>1,"age","immunity"),agedepfact=agedep_fact,
             suscept=paste0(round(unique(delta_primary),3),collapse=","),dep_strength=k)
    if (k==1 & k_dep==1){season_means_symptcases=zz; df_sympt_cases_all=df_sympt_cases} else {
      season_means_symptcases=bind_rows(season_means_symptcases,zz); df_sympt_cases_all=bind_rows(df_sympt_cases_all,df_sympt_cases)}
  } # for (k_dep in 1:2)
} # for (k in 1:10)

# save outputs
saveRDS(list("inf_fracts_season"=inf_fracts_season,"dyn_parscan"=dyn_parscan,"season_means_byinf"=season_means_byinf,
             "season_means_symptcases"=season_means_symptcases,"df_sympt_cases_all"=df_sympt_cases_all),
        file="simul_output/parscan/parscan_outputs.rds")

### plots --------------------------------------------------------
# dynamics
xval_lims=c(npi_year-0.27,npi_year+3.35); n_type="immunity"
ggplot(subset(dyn_parscan,grepl('I',name) & type %in% n_type & agegroup %in% c(1,3,5) & t/365>xval_lims[1] & t/365<xval_lims[2]& t %% 7==0 & 
                dep_strength %in% c(1,4,7,10)), aes(x=t/365,y=value_fract,group=name)) + 
  geom_area(aes(fill=infection),position=position_stack(reverse=T),color="black",size=0.25) +
  facet_grid(dep_strength~agegroup_name,labeller=labeller(dep_strength=label_both)) + # ,scales=scaleval xvalbreaks vertline_x highl_lims
  theme_bw() + standard_theme + theme(axis.text.x=element_text(size=11,vjust=0.5),
        axis.text.y=element_text(size=12),legend.position="bottom",legend.title=element_blank(),strip.text=element_text(size=12)) +
  scale_x_continuous(breaks=(0:30)*0.5,expand=expansion(0.01,0)) + scale_y_continuous(expand=expansion(0.01,0)) +
  geom_rect(aes(xmin=shutdwn_lims[1]/365,xmax=shutdwn_lims[2]/365,ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
  geom_vline(data=subset(seas_lims,on>xval_lims[1]&off<xval_lims[2]) %>% pivot_longer(cols=!season),aes(xintercept=value),
    color="blue",linetype="dashed",size=0.3) + xlab('years') + ylab("proportion of age group") + 
  labs(subtitle=paste0("susceptibility~f(",n_type,")")) # ,caption=caption_txt
# library(cowplot)
ggsave(paste0("simul_output/parscan/dynamics_by_inf_",n_type,"dep.png"),width=32,height=22,units="cm")

# % first infections (out of all infs)
# ggplot(subset(inf_fracts_season %>% mutate(max_fract=ifelse(epiyear!=npi_year+1,max_fract,NA)),
#               epiyear>5 & epiyear<=11 & infection_bin=="first" & as.numeric(agegroup_name)<=6 & dep_strength<=3)) + 
#   geom_hpline(aes(x=factor(epiyear),y=max_fract*1e2,color=factor(dep_strength),linetype=type),width=0.9) +
#   facet_wrap(~agegroup_name,scales="free") +theme_bw() + standard_theme + labs(color="",caption=caption_txt) + # ,subtitle=subtitle_str
#   geom_rect(aes(xmin=npi_year-5+0.5,xmax=npi_year-5+1.5,ymin=-Inf,ymax=Inf),fill="grey",alpha=0.2) + 
#   geom_vline(xintercept=(1:12)+0.5,size=0.2,linetype="dashed",color="black") +
#   theme(panel.grid.major.x=element_blank()) + xlab("RSV season (epi-year)") + ylab("% first infections")
# # SAVE
# ggsave("simul_output/parscan/firstinf_share.png", width=32,height=18,units="cm")

# plot mean age per season by infection type
# pal <- wes_palette("Zissou1", 10, type = "continuous")
df_plot_byinf <- subset(season_means_byinf,epi_year>=npi_year & epi_year<=npi_year+3 & dep_strength %in% c(1,3,4,7)) %>% 
  mutate(weighted_sum=ifelse(epi_year==npi_year+1,NA,weighted_sum)) %>% group_by(type,dep_strength,infection) %>%
  mutate(pre_npi_mean_age=mean(weighted_sum[epi_year>8 & epi_year<=11]),weighted_sum_div=weighted_sum/pre_npi_mean_age,
         weighted_sum_div=weighted_sum-pre_npi_mean_age)
ggplot(df_plot_byinf) + geom_hpline(aes(x=epi_year,y=weighted_sum_div,color=dep_strength),width=0.9) + 
  facet_grid(infection~type,scales="free") + theme_bw() + standard_theme + theme(panel.grid.major.x=element_blank()) +
   geom_rect(aes(xmin=npi_year+0.5,xmax=npi_year+1.5,ymin=-Inf,ymax=Inf),fill="grey",alpha=0.2) + scale_color_gradientn(colours=pal) + 
  labs(color="",caption=caption_txt,subtitle=paste0("susceptibility~",unique(season_means_byinf$type))) + 
  geom_vline(xintercept=(10:14)+0.5,size=0.2,linetype="dashed",color="black") + xlab("RSV season (epi year)") + ylab("mean age")
# save
ggsave("simul_output/parscan/weighted_meanage_perseas_byinf.png",width=32,height=18,units="cm")

# plot means per season for ALL infections
# season_means_symptcases=subset(mean_age_symptcases_at_t, grepl("in",on_off)) %>% mutate(ageweight_aver=mean_age_at_t*suminfs) %>% 
#   group_by(epi_year) %>% summarise(weighted_sum=ifelse(sum(suminfs)>5e4,sum(ageweight_aver)/sum(suminfs),NA))
# PLOT
plot_xlims=c(npi_year-2,npi_year+4)
df_plot = subset(season_means_symptcases %>% mutate(weighted_sum=ifelse(epi_year==npi_year+1,NA,weighted_sum)),
        epi_year>plot_xlims[1] & epi_year<=plot_xlims[2] & dep_strength %in% c(1,4,7,10)) %>% group_by(type,dep_strength) %>% 
  mutate(pre_npi_mean_age=mean(weighted_sum[epi_year>8 & epi_year<=11]),weighted_sum_div=weighted_sum/pre_npi_mean_age,
         weighted_sum_subtr=weighted_sum-pre_npi_mean_age)
# plot
ggplot(df_plot) +  # geom_hpline(aes(x=epi_year,y=weighted_sum,color=ifelse(epi_year==npi_year+2,"red",NA)),width=0.9,show.legend=F) +
  geom_hpline(aes(x=epi_year,y=weighted_sum_div,color=factor(type)),width=0.9) + # scale_color_gradientn(colours=pal) + 
  facet_wrap(~dep_strength,scales="free",labeller=labeller(dep_strength=label_both)) + 
  theme_bw() + standard_theme + theme(panel.grid.major.x=element_blank()) + scale_x_continuous(breaks=1:round(n_years)) +
  geom_rect(aes(xmin=npi_year+0.5,xmax=npi_year+1.5,ymin=-Inf,ymax=Inf),fill="pink",alpha=0.15) + 
  geom_vline(xintercept=(plot_xlims[1]:plot_xlims[2])+0.5,size=0.2,linetype="dashed",color="black") + 
  geom_hline(yintercept=1,size=0.2,linetype="dashed",color="black") +
  theme(axis.text.x=element_text(vjust=0.5,size=13),axis.text.y=element_text(size=13),axis.title.x=element_text(size=15),
     axis.title.y=element_text(size=15),legend.text=element_text(size=16)) + xlab("RSV season (epiyear)") +
  ylab("average age of infection (compared to pre-pandemic baseline)") + labs(color="",caption=caption_txt)
# ggsave(paste0("simul_output/parscan/all_cases_abs_mean_age_deptype.png"),width=32,height=22,units="cm")
# ggsave(paste0("simul_output/parscan/all_cases_DIFF_mean_age_deptype.png"),width=32,height=22,units="cm")
ggsave(paste0("simul_output/parscan/all_cases_RATIO_mean_age_deptype.png"),width=32,height=22,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# initvals_sirs_model <- fcn_set_initconds(init_set="previous",init_cond_src="file",ode_solution,init_seed=10,
#                                          seed_vars="all","simul_output/df_ode_solution_UK_long.RDS")
# foldername <- paste0("simul_output/parscan/",foldertag,paste0(names(scan_params),sapply(scan_params, length),collapse = "_"))
# if (dir.exists(foldername)) {foldername=paste0(foldername,gsub("2021-| |-","_",format(Sys.time(), "%F %H-%M")))}
# dir.create(foldername); sink(paste0(foldername,"/timestamp.txt")); cat(format(Sys.time(), "%F %H-%M")); sink()
# ### # param table
# scan_params=list(susc=seq(0.06,0.11,length.out=4),expdep=1,agedep=c(1.1,1.3,1.6,2.1),t_npi=c(0,2,4),seasf=c(0.5,1,2))
# param_table=round(expand.grid(scan_params),3); colnames(param_table)=names(scan_params)
# if (length(unique(param_table$agedep))==1) {sel_col=which(grepl("age",colnames(param_table))); foldertag="expdep_"} else {
#   sel_col=which(grepl("exp",colnames(param_table))); foldertag="agedep_"}
# write_csv(param_table,paste0(foldername,"/param_table.csv"))
# ### ### ### season limits
# g(forcing_vector_npi,shutdwn_lims,seas_force,seas_lims) %=% fun_shutdown_seasforc(timesteps,elem_time_step=0.5,
#  forcing_above_baseline=2,npi_strength=0.5,npi_year=round(n_years-4),peak_week=46,season_width=4,npi_on=2,npi_off=2,n_prec=0.01,n_sd=2)
# seas_lims$on=floor(seas_lims$on)+41*7/365; seas_lims$off=round(seas_lims$off)+11*7/365
# parnames=c("suscept_abs_level","suscept_expdep","suscept_agedep","npi_lims","forcing_strength"); 
# y=data.frame(matrix(rep(NA,length(parnames)),ncol=length(parnames))); colnames(y)=parnames
# ##### start loop
# cores=detectCores(); cl<-makeCluster(cores[1]-1); registerDoParallel(cl)
# df_ode_sol_sympt_cases <- foreach(k_paramscan=1:nrow(param_table),.combine=rbind,
#                                   .packages=c("tidyr","deSolve","dplyr","RcppRoll")) %dopar% {
#                                     # assign params
#                                     y[,parnames]=array(param_table[k_paramscan,])
# # NPI onset/end
# npi_reduc_strength=0.5; npi_year=round(n_years-4); preseas_npi_on=y$npi_lims; postseas_npi_off=y$npi_lims-3
# list_seas_force_npi=fun_shutdown_seasforc(timesteps,elem_time_step=0.5,
#     y$forcing_strength,npi_reduc_strength,npi_year,peak_week=46,season_width=4,preseas_npi_on,postseas_npi_off,n_prec=0.01,n_sd=2)
# # g(forcing_vector_npi,shutdwn_lims,seas_force_no_npi,seas_lims)
# # susceptibility
# delta_susc=sapply(1:n_age, function(x) (y$suscept_abs_level/y$suscept_expdep^(0:2))/((y$suscept_agedep^(x-1))*rsv_age_groups$value[x]))
# # SIMULATE
# approx_seas_forc <- approxfun(data.frame(t=timesteps,seas_force=list_seas_force_npi[[1]]))
# approx_introd <- approxfun(data.frame(t=timesteps,as.numeric(timesteps %% 30==0)*100))
# params=list(birth_term,K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds,elem_time_step,delta_susc)
# ode_solution <- lsoda(initvals_sirs_model,timesteps,func=sirs_seasonal_forc,parms=params)
# df_ode_solution_tidy=fun_process_simul_output(ode_solution,varname_list,n_age,n_inf,rsv_age_groups,neg_thresh=-1e-3)[[2]]
# # calc symptomatic cases
# clin_fract_age_exp=fcn_clin_fract_age_exp(agedep_fact=1,clinfract_expos=rep(1,n_inf),rsv_age_groups,n_age,n_inf)
# ### CALCULATE SYMPTOMATIC CASES & sum of 1,2,3rd infections -------------------------------
# temp=left_join(
#   data.frame(subset(fun_symptomcases_table(df_ode_solution_tidy,clin_fract_age_exp,c("infection","agegroup")),t/365>npi_year-2),
#   parset_id=k_paramscan),
#   data.frame(round(param_table[k_paramscan,],2),parset_id=k_paramscan))
#                                   }
# stopCluster(cl)
# qsave(df_ode_sol_sympt_cases,paste0(foldername,"/df_ode_sol_sympt_cases_EXP_DEP",nrow(param_table),"_parset.qs"))