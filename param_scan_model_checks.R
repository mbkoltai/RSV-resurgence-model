rm(list=ls()); currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
# load constant parameters and functions
source("load_params.R"); library(wesanderson)
# estimated attack rates
estim_attack_rates <- data.frame(agegroup_name=rsv_age_groups$agegroup_name, # paste0("age=",,"yr")
      median_est=c(rep(65,4),rep(40,4),10,8,5)) %>% mutate(min_est=median_est*0.25,max_est=median_est*2.5,
      median_all_inf=c(rep(70,4),rep(60,4),50,30,20),min_est_all_inf=median_all_inf/2,max_est_all_inf=median_all_inf*2)
# write estim_attack_rates # write_csv(estim_attack_rates,file="data/estim_attack_rates.csv")
# % cases within season (filtering parameter sets)
# seas_conc_lim=0.8
npi_dates<-as.Date(c("2020-03-26","2021-05-17"))
partable <- bind_rows(expand.grid(list(exp_dep=seq(1/4,2,1/8),age_dep=seq(1/8,1,1/16),seasforc_width_wks=c(3,5,7),
      R0=1+(0:5)/10,peak_week=c(48),seasforce_peak=c(3/4,1,5/4),omega=c(1/250,1/350,1/450)) ) )
l_delta_susc <- lapply(1:nrow(partable), function(n_p) {sapply(1:n_age,
                function(x) {(1*exp(-partable$exp_dep[n_p]*(1:3)))/(exp(partable$age_dep[n_p]*x))})} ) 
partable <- partable %>% mutate(par_id=row_number(),
      const_delta=R0/unlist(lapply(l_delta_susc, function(x) R0_calc_SIRS(C_m,x,rho,n_inf))),seas_conc_lim=0.85,
      npi_start=npi_dates[1],npi_stop=npi_dates[2],seas_start_wk=42,seas_stop_wk=8)
# filtering param sets, from 1st fit: 
# selected parsets are along the line `age=-exp/3+5/6` (and the point (age,exp)=(1/8,1.75))
partable <- partable %>% mutate(age_dep_fit=5/6-exp_dep/3) %>% filter(abs(age_dep-age_dep_fit)/age_dep<1/3) %>%
  select(!age_dep_fit)
write_csv(partable,"partable.csv"); write_csv(partable,"partable.csv")
# agegroup indices for maternal immunity
mat_imm_flag <- TRUE; mat_imm_inds<-list(fun_sub2ind(i_inf=1,j_age=1,"R",c("S","I","R"),n_age,3),
                                         fun_sub2ind(i_inf=c(1,2,3),j_age=9,"R",c("S","I","R"),n_age,3),
                                         fun_sub2ind(i_inf=c(1,2,3),j_age=9,"S",c("S","I","R"),n_age,3))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# RUN SIMULATIONS 
# write file that'll run scripts
simul_length_yr<-15; n_post_npi_yr<-4; n_core<-16; memory_max <- 8
partable_filename <- "simul_output/parscan/parallel/partable.csv"; write_csv(partable,file=partable_filename); 
system(paste0(c("Rscript fcns/write_run_file.R",n_core,nrow(partable),simul_length_yr,n_post_npi_yr,
                partable_filename,"data/estim_attack_rates.csv nosave sep",memory_max),collapse=" "))
# run calculation
system("sh run_all_parallel_scan.sh")
# download results from cluster
system("scp lshmk17@hpclogin:RSV-model/simul_output/parscan/parallel/results_summ_all.csv simul_output/parscan/parallel/")
system("scp lshmk17@hpclogin:RSV-model/simul_output/parscan/parallel/results_dyn_all.csv simul_output/parscan/parallel/")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# READ IN RESULTS
peak_week_lims <- c(48,2)
foldername<-"simul_output/parscan/parallel/parsets_2809/" # partable <- read_csv(paste0(foldername,"partable.csv"))
# results_dyn_all <- read_csv("simul_output/parscan/parallel/results_dyn_all.csv")
results_summ_all <- read_csv(paste0(foldername,"results_summ_all.csv")) %>%
  mutate(max_incid_week_check=ifelse(max_incid_week>=peak_week_lims[1]|max_incid_week<=peak_week_lims[2],TRUE,FALSE))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###  
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# plot attack rates by age group and years
check_crit=11/11; sel_yrs<-2019; n_sel_yr=length(sel_yrs)
# unlist(lapply(1:length(l_delta_susc), function(x) R0_calc_SIRS(C_m,partable$const_delta[x]*l_delta_susc[[x]],rho,n_inf)))
all_sum_inf_epiyear_age_filtered <- left_join(results_summ_all %>% filter(epi_year %in% sel_yrs),
              partable %>% rename(forcing_peak_week=peak_week),by=c("par_id","seasforce_peak","R0")) %>% 
  group_by(seasforce_peak,exp_dep,age_dep,seasforc_width_wks,par_id) %>% 
 filter(sum(attack_rate_check)>=round(n_age*n_sel_yr*check_crit) & 
        sum(seas_share_check)>=round(n_age*n_sel_yr*check_crit) )

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# check DYNAMICS of SELECTED SIMUL
start_date=as.Date("2018-09-01")
sel_parsets=unique(all_sum_inf_epiyear_age_filtered$par_id)
ggplot(results_dyn_all %>% filter((par_id %in% sel_parsets[1:11]) & agegroup<=5) %>% mutate(date=start_date+t-min(t)) %>%
         filter(date<as.Date("2023-04-15") & date>as.Date("2019-08-01"))) + 
  geom_line(aes(x=date,y=value,color=factor(par_id))) + 
  facet_grid(infection~agegroup,scales="free_y",labeller=labeller(infection=label_both,agegroup=label_both)) +
  # scale_color_brewer(palette = "YlOrRd") + # scale_color_discrete() + 
  geom_rect(xmin=npi_dates[1],xmax=npi_dates[2],ymin=-Inf,ymax=Inf,fill="grey",alpha=0.01) +
  geom_vline(xintercept=as.Date(paste0(2018:2022,"-12-13"))-56,linetype="dashed",size=1/4) +
  geom_vline(xintercept=as.Date(paste0(2018:2022,"-12-13"))+56,linetype="dashed",size=1/4) +theme_bw() + standard_theme + 
  theme(legend.position="none") + scale_x_date(date_breaks="4 month") + xlab("") + ylab("") + labs(color="# par ID")

### ### ### ### ### ### ### ### ### ### ### ### ### ###
# PLOT dynamics by age groups (one simulation)
n_sel=sel_parsets[11]
sel_weeks <- results_dyn_all %>% filter(par_id==n_sel) %>% mutate(date=start_date+t-min(t),week=week(date),year=year(date)) %>% 
  filter(week %in% c(9,41,49)) %>% group_by(year,agegroup,week) %>% filter(date==min(date) & infection==1)
fcn_plot_timecourse_by_agegr(results_dyn_all %>% filter(par_id==n_sel) %>% mutate(date=start_date+t-min(t)) %>%
    filter(t %% 7==0 & agegroup<=9 & date>as.Date("2017-07-01") & date<as.Date("2022-04-01")),
    agegroup_name=rsv_age_groups$agegroup_name,sel_agelim=9,varname="value",npidates=npi_dates,date_break_val="2 month",
              selweeks=sel_weeks,alphaval=0.01,vline_w=c(1/4,1/8))
# # sum of all cases
p<-fcn_plot_timecourse_sum(results_dyn_all %>% filter(par_id==n_sel) %>% mutate(date=start_date+t-min(t)) %>%
  filter(agegroup<=9 & date>as.Date("2017-07-01") & date<as.Date("2023-04-01")) %>% 
  group_by(date,infection) %>% summarise(value=sum(value)) %>% mutate(infection=factor(infection)),npi_dates,n_peak_week=50)
p + scale_x_date(date_breaks = "month")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# PLOT attack rates, seasonal share, peak week
sel_var<-c("attack_rate_perc","seas_share","max_incid_week")
# plot_labels=c("attack rate % age group","seasonal share of infections","peak week")
estim_rates <- estim_attack_rates %>% select(agegroup_name,median_est,min_est,max_est) %>% 
  pivot_longer(!c(agegroup_name)) %>% rename(type=name) %>% mutate(name="attack_rate_perc")
estim_rates <- bind_rows(estim_rates, 
        estim_rates %>% filter(type!="median_est") %>% mutate(name="max_incid_week",value=ifelse(grepl("min",type),3,48)))
color_var<-"exp_dep" # R0 exp_dep age_dep
ggplot(all_sum_inf_epiyear_age_filtered %>% mutate(attack_rate_perc=ifelse(epi_year==2020,NA,attack_rate_perc),
  agegroup_name=factor(agegroup_name,levels=unique(agegroup_name))) %>% ungroup() %>% select(c(par_id,epi_year,
  agegroup_name,attack_rate_perc,seas_share,max_incid_week,exp_dep,age_dep,seasforc_width_wks,R0)) %>% 
    pivot_longer(!c(epi_year,agegroup_name,par_id,exp_dep,age_dep,seasforc_width_wks,R0)) ) +
  geom_hpline(aes(x=age_dep,y=value,color=get(color_var),group=par_id),width=0.1,size=1/2)+#position=position_dodge(width=1)
  facet_grid(name~agegroup_name,scales="free_y") + scale_y_continuous(expand=expansion(0.02,0))+
  scale_color_gradient2(midpoint=median(c(t(unique(all_sum_inf_epiyear_age_filtered[,color_var])))),low="blue",mid="white",high="red") +
  geom_hline(data=estim_rates,aes(yintercept=value),linetype="dashed",size=1/4)+ 
  xlab("age-dependence")+ylab("")+theme(legend.position="top")+theme_bw()+standard_theme # +labs(color="parameter ID")
# save
ggsave(paste0(foldername,"parscan_attack_rates_filtered_",color_var,".png"),width=32,height=20,units="cm")

######
# which parsets selected? # library("ggrepel")
partable_filtered <- partable %>% filter(par_id %in% unique(all_sum_inf_epiyear_age_filtered$par_id))
write_csv(partable_filtered,paste0(foldername,"partable_filtered.csv"))
# plot
ggplot(partable_filtered %>% mutate(sel_par=TRUE),aes(x=exp_dep,y=age_dep)) + 
   geom_point(aes(fill=sel_par,color=sel_par),size=1.5) + facet_grid(R0~seasforce_peak+seasforc_width_wks,
        labeller=labeller(seasforc_width_wks=label_both,seasforce_peak=label_both,R0=label_both)) +
  # geom_smooth(method="lm",color="black",size=1/2,se=F) + # geom_smooth(method="loess",se=F,size=1/2) + 
  geom_point(data=partable,aes(x=exp_dep,y=age_dep),color="grey",size=1/2) + scale_x_continuous(breaks=2*(1:8)/8) + 
  # geom_text(data=partable_filtered,aes(x=exp_dep,y=age_dep,label=par_id),position=position_jitter(width=1/9,height=1/9),
  #           size=2) 
   theme(legend.position="none",axis.text.x=element_text(size=9),axis.text.y=element_text(size=9),strip.text=
         element_text(size=7))+ labs(color="accepted") + theme_bw() + xlab("exposure") + ylab("age") + standard_theme
# selected parameter sets
ggsave(paste0(foldername,"sel_parsets_scatterplot.png"),width=40,height=20,units="cm")

# PCA on parameter sets
par_pca <- prcomp(partable_filtered %>% select(exp_dep,age_dep,seasforc_width_wks,R0,peak_week,seasforce_peak),
              center=TRUE,scale.=TRUE)
# library(devtools);install_github("vqv/ggbiplot"); library(ggbiplot)
ggbiplot(par_pca,groups=factor(partable_filtered$seasforce_peak),ellipse=TRUE) # ,labels=partable_filtered$par_id

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# % age groups satisfy criteria each year
# model checks: 1) attack rate 2) seasonal concentration 3) peak week 4) variation in attack rates across epi_years below a limit
ar_diff_crit=c(10,2)
parsets_filtered_allyears <- all_sum_inf_epiyear_age %>% filter(epi_year>=2016 & epi_year<=2019) %>% 
  group_by(par_id,epi_year,dep_val,seasforce_peak,R0,dep_type) %>%
  summarise(seas_share_check=sum(seas_share_check)/n_age, attack_rate_check=sum(attack_rate_check)/n_age,
            peak_week_check=sum(max_incid_week>=peak_week_lims[1]|max_incid_week<=peak_week_lims[2])/n_age,
            attack_rate_youngest=attack_rate_perc[agegroup==1],attack_rate_oldest=attack_rate_perc[agegroup==max(agegroup)]) %>% 
  group_by(par_id) %>% 
  mutate(diff_AR_young=max(attack_rate_youngest)-min(attack_rate_youngest),diff_AR_old=(max(attack_rate_oldest)-min(attack_rate_oldest))) %>%
  filter(all(attack_rate_check>=check_crit) & all(seas_share_check>=check_crit) & all(peak_week_check>=check_crit) &
           diff_AR_young<attack_rate_youngest*0.1 & diff_AR_old<attack_rate_oldest*0.1 ) 
## SAVE
write_csv(parsets_filtered_allyears,paste0(foldername,"parsets_filtered_allyears.csv"))
## COLORS
uniq_dep_vals=(parsets_filtered_allyears %>% group_by(dep_type) %>% summarise(n=length(unique(R0))))$n
colorpal=c(colorRampPalette(colors=c("orange","red"))(uniq_dep_vals[1]),colorRampPalette(colors=c("grey","black"))(uniq_dep_vals[2]))
# plot
ggplot(parsets_filtered_allyears %>% pivot_longer(c(seas_share_check,attack_rate_check,peak_week_check)) %>% group_by(par_id) %>% 
         ungroup() %>% mutate(value=mean(value),dep_type_R0=factor(paste0(dep_type,", R=",R0)),name=gsub("_"," ",gsub("_check","",name))) %>%
         rename(`dependence strength`=dep_val)) + 
  geom_point(aes(x=name,y=round(value*1e2,1),color=dep_type_R0,shape=factor(name),group=R0),position=position_dodge(width=1),size=3) +
  facet_grid(seasforce_peak~`dependence strength`,labeller=labeller(`dependence strength`=label_both)) +
  geom_rect(data=parsets_filtered_allyears %>% rename(`dependence strength`=dep_val) %>% group_by(seasforce_peak,`dependence strength`) %>% 
              summarise(`dependence strength`=unique(`dependence strength`),seasforce_peak=unique(seasforce_peak)),xmin=-Inf,xmax=Inf,
              ymin=-Inf,ymax=Inf,alpha=0.1,fill=NA,color="blue",size=1.5) + 
  scale_color_manual(values=colorpal) + scale_y_continuous(limits=c(95,100),breaks=c(96,98,100)) +
  geom_vline(xintercept=(0:4)+1/2,size=1/3,linetype="dashed") + xlab("") + ylab("% age groups satisfy criteria") + 
  scale_x_discrete(expand=expansion(0.1,0)) + theme_bw() + standard_theme + theme(legend.position="top",legend.text=element_text(size=15),
    axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),strip.text=element_text(size=12),axis.title.y=element_text(size=15)) + 
  guides(color=guide_legend(byrow=TRUE,ncol=5))+ labs(color="",linetype="",shape="") 
#caption=paste0("model check minimum=",round(check_crit,2))
# save
ggsave(paste0(foldername,"parscan_modelcheck_highlight.png"),width=35,height=30,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Generate dynamics of selected parsets for the 2020- period
# calc or read in
parsets_filtered <- parsets_filtered_allyears %>% filter(epi_year==2019) %>% group_by(dep_type,dep_val) %>% mutate(n_dep_val=n()) %>%
  filter(n_dep_val>1) %>% group_by(dep_type) %>% mutate(dep_val_sel=as.numeric(factor(dep_val)))

calc_flag=T; plot_flag=F
for (k_npi_str in 0:4){
npi_red_str=0+k_npi_str/n_step
if (calc_flag){
for (k_par_filt in 1:nrow(parsets_filtered)) {
  if (k_par_filt==1) {list_df_cases_infs=list(); list_df_cases_infs_broad_age_grp<-list()}
  delta_primary=as.numeric(partable[parsets_filtered$par_id[k_par_filt],] %>% select(contains("delta")))
  # assign age dependence factor
    agedep_fact=partable$agedep_val[parsets_filtered$par_id[k_par_filt]]
  g(delta_susc,delta_susc_prop) %=% fcn_delta_susc(delta_primary,n_age,n_inf,agedep_fact,rsv_age_groups$stationary_popul)
# seasonal forcing
  g(n_years,timesteps,simul_start_end,forcing_vector_npi) %=% fun_shutdown_seasforc(npi_dates,years_pre_post_npi=c(3,3),
    season_width_wks=seasforc_width_wks,init_mt_day="06-01",ifelse(grepl("exp",parsets_filtered$dep_type[k_par_filt]),45,49),
    forcing_above_baseline=parsets_filtered$seasforce_peak[k_par_filt],npireduc_strength=0)
  npi_inds=as.numeric(npi_dates[1]-simul_start_end[1]):as.numeric(npi_dates[2]-simul_start_end[1])
  forcing_vector_npi[npi_inds]=1 + (forcing_vector_npi[npi_inds]-1)*npi_red_str
# set params
  params<-list(list(birth_rates,matrix(unlist(lapply(uk_death_rate,function(x) rep(x,n_inf*n_compartment)))) ),
             K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds,delta_susc)
# interpolation fcns for seas forcing & extern introds
  approx_seas_forc <- approxfun(data.frame(t=timesteps,seas_force=forcing_vector_npi))
  approx_introd <- approxfun(data.frame(t=timesteps,as.numeric(timesteps %% 30==0)*10))
  if (!mat_imm_flag){
  ode_solution<-lsoda(stat_sol_allparsets[,parsets_filtered$par_id[k_par_filt]],timesteps,func=sirs_seasonal_forc,parms=params)} else {
    params[[7]] <- mat_imm_inds; ode_solution <- lsoda(initvals_sirs_model,timesteps,func=sirs_seasonal_forc_mat_immun,parms=params)     }
# process output
  broad_age_grp<-c("0-2y","2-5y","5-15y",">15y")
  list_df_cases_infs[[k_par_filt]] <- fcn_process_odesol_incid(ode_solution,n_age,n_inf,n_compartment,simul_start_end) %>%
    filter(date>as.Date("2018-10-01") & date<as.Date("2023-02-01") ) %>% 
    mutate(R0=parsets_filtered$R0[k_par_filt],par_id=parsets_filtered$par_id[k_par_filt],dep_type=parsets_filtered$dep_type[k_par_filt],
         dep_val=parsets_filtered$dep_val[k_par_filt], seasforce_peak=parsets_filtered$seasforce_peak[k_par_filt],npi_str=npi_red_str) 
  list_df_cases_infs_broad_age_grp[[k_par_filt]] <- list_df_cases_infs[[k_par_filt]] %>% 
    mutate(agegroup_large=factor(broad_age_grp[ifelse(agegroup<=4,1,ifelse(agegroup<=7,2,ifelse(agegroup==8,3,4)))],
         levels=broad_age_grp)) %>% group_by(date,dep_val,seasforce_peak,dep_type,npi_str,R0,agegroup_large) %>% 
    summarise(value=ifelse(sum(value)>0,sum(value),NA))
# print progress
print(paste0(k_par_filt,"/",nrow(parsets_filtered),", npi=",npi_red_str))
if (k_par_filt==nrow(parsets_filtered)) {list_df_cases_infs <- bind_rows(list_df_cases_infs)
                                        list_df_cases_infs_broad_age_grp <- bind_rows(list_df_cases_infs_broad_age_grp)}
} # end of loop
} else { list_df_cases_infs<-read_csv(paste0(foldername,"df_incid_seaswidth",seasforc_width_wks,"_npi_contacts",npi_red_str,".csv"))}
# save when done with parsets
write_csv(list_df_cases_infs,paste0(foldername,"df_cases_infs_npi_contactlevel_",npi_red_str*100,"pct.csv"))
write_csv(list_df_cases_infs_broad_age_grp,paste0(foldername,"df_cases_infs_broad_agegr_npi_contactlevel_",npi_red_str*100,"pct.csv"))
if (k_npi_str==0) {npi_str_scan_df_cases_infs<-list_df_cases_infs_broad_age_grp} else {
  npi_str_scan_df_cases_infs<-rbind(npi_str_scan_df_cases_infs,list_df_cases_infs_broad_age_grp)}

} # end of loop on NPI strength

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# summarise results
out_season_date<-as.Date("2021-08-01")
# plot
df_plot <- left_join(npi_str_scan_df_cases_infs %>% filter(!agegroup_large %in% ">15y") %>% filter(date>as.Date("2019-09-15")) %>% 
  filter(date<as.Date("2022-03-15")) %>% mutate(type_depval=paste0(dep_type,dep_val)) %>% 
    filter(type_depval %in% unique(paste0(parsets_filtered$dep_type,parsets_filtered$dep_val)) ) %>% 
    group_by(dep_type) %>% mutate(dep_val_sel=as.numeric(factor(dep_val))), 
  partable %>% select(par_id,dep_type,R0,dep_val,seasforce_peak),by=c("dep_type","R0","dep_val","seasforce_peak")) %>%
  mutate(epi_year=ifelse(date>ymd(paste(year(date),"-07-01")),year(date),year(date)-1)) %>%
  group_by(par_id,epi_year,agegroup_large,npi_str) %>% mutate(peak_week=week(date[value==max(value,na.rm=T)])) %>% ungroup()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# susceptibility parameters for selected param sets
for (k_par in 1:nrow(parsets_filtered)){
  if (k_par==1) {df_suscept=data.frame()}
  delta_primary <- as.numeric(left_join(parsets_filtered,partable,by=c("par_id","dep_val","dep_type","seasforce_peak","R0"))[k_par,] %>% 
                                ungroup() %>% select(contains("delta")) )
  agedep_fact<-left_join(parsets_filtered,partable,by=c("par_id","dep_val","dep_type","seasforce_peak","R0"))[k_par,]$agedep_val
  g(delta_susc,delta_susc_prop) %=% fcn_delta_susc(delta_primary,n_age,n_inf,agedep_fact,rsv_age_groups$stationary_popul)
  df <- data.frame(t(round(delta_susc_prop,5))); colnames(df) <-c("inf #1","inf #2","inf #3")
  df <- bind_cols(df,parsets_filtered[k_par,]) %>% rowid_to_column("agegr_no") %>% 
    mutate(agegroup=factor(rsv_age_groups$agegroup_name,levels=unique(rsv_age_groups$agegroup_name)))
  if (k_par==1) {df_suscept=df} else {df_suscept=rbind(df_suscept,df)
  } 
}
# dep_val values present in both dep types
joint_dep_vals=(df_suscept %>% group_by(dep_type) %>% summarise(dep_val=unique(dep_val)) %>% group_by(dep_val) %>% 
                  summarise(n=n()) %>% filter(n==2))$dep_val
# plot susceptibility age dependencies
ggplot(df_suscept %>% group_by(dep_type,dep_val_sel,agegroup) %>% filter(R0==min(R0)) %>% pivot_longer(c(`inf #1`,`inf #2`,`inf #3`)) %>%
         mutate(name=gsub("inf","infection",name)) %>% rename(dependence=dep_val_sel)) + 
  geom_hpline(aes(x=agegroup,y=value,group=interaction(dependence,name),color=factor(name)),
              width=0.95/3,size=1,position=position_dodge(width=1)) + 
  facet_grid(dependence~dep_type,labeller=labeller(dependence=label_both)) + 
  geom_vline(xintercept=0.5+(0:11),size=1/3,linetype="dashed") + 
  scale_y_log10(expand=expansion(0.04,0)) + scale_x_discrete(expand=expansion(0.02,0)) + 
  theme_bw() + standard_theme + ylab("susceptibility to infection") + labs(color="",linetype="") + 
  guides(color=guide_legend(byrow=TRUE))  + theme(legend.text=element_text(size=14),legend.title=element_text(size=14),
      axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),strip.text=element_text(size=15))
# save
ggsave(paste0(foldername,"susceptibility_params.png"),width=28,height=22,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# plot dynamics separately by NPI strength
library(rlang)
sel_var<-c("dep_val_sel","npi_str")[1]; color_var=ifelse(sel_var=="dep_val_sel","npi_str","dep_val_sel") 
# linetp_var="dep_val_sel" # seasforce_peak
color_lab <- ifelse(grepl("npi",color_var),"contact levels under NPI",ifelse(grepl("dep_val_sel",color_var),"suscept. ~ age/exposure","R0"))
# LOOP
for (sel_var_val in unlist(array(unique(df_plot[,sel_var]))) ){
  caption_str<-paste0(ifelse(grepl("npi",sel_var),paste0("contacts during NPI=",sel_var_val*1e2,"% (above baseline)"),
     ifelse(sel_var=="dep_val_sel",paste0("strengh of dep.=",sel_var_val),paste0("R0=",sel_var_val))),
     ", width seas forc=",seasforc_width_wks,"wks")
  # colors
  n_col<-nrow(unique(df_plot[,color_var]))
  if (n_col>=5) {pal<-wes_palette("Zissou1",n_col,type="discrete")} else {pal<-wes_palette("Zissou1",n_col+1,type="discrete")[2:n_col+1]}
    # average across seas forcing terms
  df_plot_dyn<-df_plot %>% filter(!!sym(sel_var)==sel_var_val) %>%  # dep_val_sel %in% joint_dep_val_sels & 
    group_by(date,agegroup_large,dep_type,!!sym(sel_var),!!sym(color_var))
  if (nrow(df_plot_dyn)>0){
    ggplot(df_plot_dyn %>% summarise(mean_val=mean(value),min_val=min(value),max_val=max(value)),
           aes(x=date,y=mean_val/1e3,group=interaction(dep_type,get(sel_var),get(color_var)))) + # R0,
      geom_line(aes(color=get(color_var))) + # ,linetype=factor(get(linetp_var)) par_id,get(color_var)
      geom_ribbon(aes(ymin=min_val/1e3,ymax=max_val/1e3,fill=get(color_var)),alpha=0.1) + # 
      scale_color_gradientn(colours=pal) + scale_fill_gradientn(colours=pal) +
      facet_grid(agegroup_large~dep_type,scales="free_y") + # ,nrow=length(unique(npi_str_scan_df_cases_infs$agegroup_large))
      scale_x_date(date_breaks="2 month",expand=expansion(0.001,0)) + scale_y_continuous(expand=expansion(0,0)) + 
      geom_rect(xmin=npi_dates[1],xmax=npi_dates[2],ymin=-Inf,ymax=Inf,fill="grey",alpha=0.01) +
      geom_vline(data=list_df_cases_infs %>% mutate(year=year(date),week=week(date)) %>% 
                   filter((date<npi_dates[1] | date>npi_dates[2]) & week %in% c(41,7)) %>% group_by(year,week) %>% summarise(date=min(date)),
                 aes(xintercept=date),size=1/4,linetype="dashed",show.legend=F) + geom_vline(xintercept=out_season_date,color="black",size=1/4)+
      theme_bw() + standard_theme + xlab("") + ylab("new cases (1000s)") + theme(legend.position="top") + 
      labs(color=color_lab,fill=color_lab,caption=caption_str)
# save
if (!dir.exists(paste0(foldername,"dynamics"))) {dir.create(paste0(foldername,"dynamics"))}
# save
ggsave(paste0(foldername,"dynamics/sel_parsets_agefacet_under15_",sel_var,"_",ifelse(grepl("npi",sel_var),
        paste0(sel_var_val*1e2,"pct"),sel_var_val),"_R0_seasfor_averaged",".png"),width=30,height=20,units="cm") }
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# calculate total infections per epi year, 2019 vs 2021, GROUPED into 0-2y, 2-5y, 5-15y
cumul_inf_pre_post_npi <- df_plot %>% filter(epi_year!=2020) %>% group_by(par_id,npi_str,epi_year,agegroup_large) %>% 
  summarise(sum_inf=sum(value),max_inf=max(value),peak_week=unique(peak_week),resurg_week_20pct=week(min(date[value>=0.2*max(value)])),
            dep_val_sel=unique(dep_val_sel),dep_val=unique(dep_val),R0=unique(R0),
    seas_forc_peak=unique(seasforce_peak),dep_type=unique(dep_type)) %>% group_by(par_id,npi_str,agegroup_large) %>% 
  mutate(sum_inf_norm=sum_inf/sum_inf[epi_year==2019],max_inf_norm=max_inf/max_inf[epi_year==2019],
      peak_week_shift=peak_week[epi_year==2021]-ifelse(peak_week[epi_year==2019]>40,peak_week[epi_year==2019],peak_week[epi_year==2019]+52),
      shift_seasonstart_week=ifelse(epi_year>=2021,ifelse(resurg_week_20pct[epi_year==2021]<10,
      resurg_week_20pct[epi_year==2021]+52,resurg_week_20pct[epi_year==2021])-
      mean(ifelse(resurg_week_20pct[epi_year<2020]<10,resurg_week_20pct[epi_year<2020]+52,resurg_week_20pct[epi_year<2020])),NA))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
## plot sum of infections / shift in peak week / shift in start week by large agegroups
# geom_point, NPI color coded, averaged over R0 and seas_forc_peak values (OR NOT)
uniq_R0=(cumul_inf_pre_post_npi %>% group_by(dep_type) %>% summarise(n=length(unique(R0))))$n; dodge_val=1
colorpal=c(colorRampPalette(colors=c("orange","red"))(uniq_R0[1]),colorRampPalette(colors=c("grey","black"))(uniq_R0[2]))
sel_var="sum_inf_norm" # sum_inf_norm | peak_week_shift | shift_seasonstart_week
ggplot(cumul_inf_pre_post_npi %>% mutate(dep_val=dep_val_sel) %>%
         filter(epi_year==2021) %>% pivot_longer(c(sum_inf_norm,max_inf_norm,peak_week_shift,shift_seasonstart_week)) %>%
     filter(grepl(sel_var,name)) %>% mutate(value=abs(value)) %>% group_by(agegroup_large,dep_type,dep_val,R0,name) %>% # ,seas_forc_peak
         summarise(mean_val=mean(value),min_val=min(value),max_val=max(value)) %>%
     mutate(name=ifelse(grepl("max",name),"PEAK","SUM"),dep_type=ifelse(grepl("age",dep_type),"~AGE","~IMMUNITY"),
    type_R0=paste0(dep_type,", R0=",R0)) %>% ungroup(), #  %>% mutate(dep_val=as.numeric(factor(dep_val)))
       aes(x=factor(dep_val),y=mean_val,color=type_R0,group=interaction(type_R0))) + # ,size=seas_forc_peak seas_forc_peak,
  geom_point(position=position_dodge(width=dodge_val)) + scale_size(range=c(0.7,2)/3) +
  geom_pointrange(aes(ymin=min_val,ymax=max_val),position=position_dodge(width=dodge_val)) +
  facet_wrap(~agegroup_large,scales="free",nrow=3) + scale_color_manual(values=colorpal) + # 
  geom_vline(xintercept=0.5+(0:8),linetype="dashed",size=1/3) + guides(color=guide_legend(ncol=5,byrow=TRUE)) + xlab("dependence on age/exp") + 
  ylab(ifelse(grepl("sum",sel_var),"ratio of infections 2021 to 2019",paste0("forward shift in ",
        ifelse(grepl("peak",sel_var),"peak","start")," week"))) + labs(size="seasonal forcing (>baseline)",color="") + theme_bw() + 
  standard_theme + theme(legend.position="top") + scale_x_discrete(expand=expansion(0.1,0))
# save
# ggsave(paste0(foldername,"resurgence_",ifelse(grepl("sum",sel_var),"cumul_inf",sel_var),"_rel_to_2019_large_agegroups.png"),
#        width=30,height=20,units="cm")
# ggsave(paste0(foldername,"resurgence_",ifelse(grepl("sum",sel_var),"cumul_inf",sel_var),"_rel_to_2019_large_agegroups_seasfor_sep.png"),
#         width=30,height=20,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# shift in mean age by large age groups
output_files <- list.files(foldername,pattern="df_cases_infs_npi_contactlevel_")
for (k_npi in 0:4){
full_output <- read_csv(paste0(foldername,output_files[grepl(paste0("_",k_npi*10,"pct"),output_files)])) %>% 
  mutate(npi_str=k_npi/10,epi_year=ifelse(date>ymd(paste(year(date),"-07-01")),year(date),year(date)-1),
  in_out_season=ifelse(week(date)<=9 | week(date)>=41,"in","out")) %>% filter(par_id %in% parsets_filtered$par_id) %>% 
  group_by(epi_year,par_id,dep_type,dep_val,R0,seasforce_peak,npi_str,agegroup) %>% 
  summarise(sum_inf=sum(value),max_inf=max(value),peak_week=week(date[value==max(value)]),
            resurg_week_20pct=week(min(date[value>=0.2*max(value)])) ) # %>% ungroup()
full_output <- left_join(full_output,rsv_age_groups %>% select(agegroup_name,mean_age_weighted) %>% mutate(agegroup=row_number()),by="agegroup")
if (k_npi==0) {npi_scan_sum_infs_by_age=full_output} else {npi_scan_sum_infs_by_age=rbind(npi_scan_sum_infs_by_age,full_output)} 
}
### end of read-in loop

# mean age before/after NPIs
full_output_mean_age = npi_scan_sum_infs_by_age %>% group_by(dep_type) %>% mutate(dep_val=as.numeric(factor(dep_val))) %>%
  group_by(epi_year,par_id,dep_type,dep_val,R0,seasforce_peak,npi_str) %>% 
  summarise(mean_age=sum((sum_inf/sum(sum_inf))*mean_age_weighted),
            mean_age_under_5=sum((sum_inf[agegroup<=7]/sum(sum_inf[agegroup<=7]))*mean_age_weighted[agegroup<=7]),
            mean_age_under_15=sum((sum_inf[agegroup<=8]/sum(sum_inf[agegroup<=8]))*mean_age_weighted[agegroup<=8])) %>%
  group_by(par_id,npi_str) %>% mutate(
   `shift in average age of infection <5y`=ifelse(epi_year>=2021,mean_age_under_5[epi_year==2021]-mean(mean_age_under_5[epi_year<2020]),NA),
   `shift in average age of infection <15y`=ifelse(epi_year>=2021,mean_age_under_15[epi_year==2021]-mean(mean_age_under_15[epi_year<2020]),NA))
   
# subset data
df_plot_mean_age_shift <- full_output_mean_age %>% pivot_longer(!c(epi_year,par_id,dep_type,dep_val,R0,seasforce_peak,npi_str)) %>%
  group_by(name,epi_year,par_id,dep_val,dep_type,R0) %>% # ,seasforce_peak
  summarise(mean_val=mean(value),min_npi_val=min(value),max_npi_val=max(value)) %>% filter(epi_year==2021 & grepl("shift",name)) %>% 
  mutate(dep_type=ifelse(grepl("age",dep_type),"~age","~immunity")) %>% ungroup() 
# plot params
ratio_diff_flag <- "shift"; dodge_val=1
title_str=ifelse(ratio_diff_flag=="shift","shift (in months) from 2019 to 2021","(mean age in 2021)/(mean age in 2019)")
# uniq_R0=(df_plot_mean_age_shift %>% group_by(dep_type) %>% summarise(n=length(unique(R0))))$n; dodge_val=1
# colorpal=c(colorRampPalette(colors=c("orange","red"))(uniq_R0[1]),colorRampPalette(colors=c("darkgrey","black"))(uniq_R0[2]))
# PLOT mean age shift by large age groups
ggplot(df_plot_mean_age_shift %>% filter(grepl("<",name)) %>% 
    mutate(name=factor(name,levels=c("shift in average age of infection <5y","shift in average age of infection <15y"))) %>% 
    mutate(type_R0=paste0(dep_type,", R0=",R0)), 
    aes(x=factor(dep_val),y=mean_val*12,group=interaction(type_R0))) + # ,size=seasforce_peak seasforce_peak,
  geom_point(aes(color=type_R0),position=position_dodge(width=dodge_val)) + scale_size(range=c(0.7,2)/2.5) +
  geom_pointrange(aes(color=type_R0,ymin=min_npi_val*12,ymax=max_npi_val*12),position=position_dodge(width=dodge_val)) +
  facet_wrap(~name,scales="free") + scale_color_manual(values=colorpal) + geom_vline(xintercept=0.5+0:8,linetype="dashed",size=1/3) +
  xlab("strength of dependence on age/exposure") + ylab(title_str) + labs(color=paste0(""),size="seasonal forcing (>baseline)") + 
  theme_bw() + standard_theme + scale_x_discrete(expand=expansion(0.12,0)) + guides(color=guide_legend(ncol=3)) + # ,byrow=TRUE
  theme(legend.position="top",axis.title.x=element_text(size=15),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        strip.text=element_text(size=15),legend.text=element_text(size=13),legend.title=element_text(size=14))
# save
# ggsave(paste0(foldername,"resurgence_mean_age_",ratio_diff_flag,"_2021_2019.png"),width=43,height=25,units="cm") 
# ggsave(paste0(foldername,"resurgence_mean_age_",ratio_diff_flag,"_2021_2019_seasfor_aver.png"),width=43,height=25,units="cm") 

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# resurgence features by YEARLY age groups
agegrp_yr<-paste0(c("0-1","1-2","2-3","3-4","4-5","5-15",">15"),"y")
sum_max_ratio_by_yearly_agegr <- npi_scan_sum_infs_by_age %>% group_by(dep_type) %>% mutate(dep_val=as.numeric(factor(dep_val))) %>%
  mutate(age_yr=agegrp_yr[findInterval(agegroup,c(0,3,5,6,7,8,9))]) %>%
  group_by(epi_year,par_id,dep_type,dep_val,R0,seasforce_peak,npi_str,age_yr) %>% 
  summarise(sum_inf=sum(sum_inf),max_inf=sum(max_inf),peak_week=mean(peak_week),resurg_week_20pct=mean(resurg_week_20pct)) %>% 
  group_by(par_id,npi_str,age_yr) %>% mutate(sum_ratio=ifelse(epi_year>=2021,sum_inf[epi_year==2021]/mean(sum_inf[epi_year<2020]),NA),
         max_ratio=ifelse(epi_year>=2021,max_inf[epi_year==2021]/mean(max_inf[epi_year<2020]),NA),
         shift_peak_week=ifelse(epi_year>=2021,
          ifelse(peak_week[epi_year==2021]<10,peak_week[epi_year==2021]+52,peak_week[epi_year==2021])-
            mean(ifelse(peak_week[epi_year<2020]<10,peak_week[epi_year<2020]+52,peak_week[epi_year<2020])),NA),
         shift_seasonstart_week=ifelse(epi_year>=2021,
                    ifelse(resurg_week_20pct[epi_year==2021]<10,resurg_week_20pct[epi_year==2021]+52,resurg_week_20pct[epi_year==2021])-
                    mean(ifelse(resurg_week_20pct[epi_year<2020]<10,resurg_week_20pct[epi_year<2020]+52,resurg_week_20pct[epi_year<2020])),NA))
#
# plot ratio of CUMUL infections 2021/2019
# subset data
df_plot_yr_cum_inf <- sum_max_ratio_by_yearly_agegr %>% filter(epi_year==2021 & !grepl(">15y",age_yr)) %>% 
  mutate(dep_type=ifelse(grepl("age",dep_type),"~age","~immunity")) %>% group_by(age_yr,dep_type,dep_val,R0) %>% # ,par_id 
  summarise(mean_val=mean(sum_ratio),min_npi_val=min(sum_ratio),max_npi_val=max(sum_ratio) ) %>% ungroup() %>%
  mutate(type_R0=factor(paste0(dep_type,", R0=",R0)) )
# filter(dep_val %in% joint_dep_vals) %>% mutate(dep_val=as.numeric(factor(dep_val))) %>% 
# plot params
# uni_dep_vals=(df_plot_yr_cum_inf %>% group_by(dep_type) %>% summarise(n=length(unique(R0))))$n
# colorpal=c(colorRampPalette(colors=c("orange","red"))(uni_dep_vals[1]),colorRampPalette(colors=c("grey","black"))(uni_dep_vals[2]))
varname <- "sum_ratio"; dodge_val=1
title_str <- ifelse(!grepl("ratio",varname),"shift in peak incidence week","ratio of cumulative infections 2021 to 2019")
# plot
ggplot(df_plot_yr_cum_inf,aes(x=factor(dep_val),y=mean_val,color=type_R0,group=interaction(seasforce_peak,type_R0),size=seasforce_peak)) + # 
  geom_point(position=position_dodge(width=dodge_val)) + scale_size(range=c(0.7,2)/3) +
  geom_pointrange(aes(ymin=min_npi_val,ymax=max_npi_val),position=position_dodge(width=dodge_val)) +
  facet_wrap(~age_yr,scales="free",nrow=3) + scale_color_manual(values=colorpal) + 
  geom_vline(xintercept=0.5+(0:8),linetype="dashed",size=1/3) + # guides(color=guide_legend(ncol=5,byrow=TRUE)) + 
  xlab("dependence on age/immunity") + ylab(title_str) + labs(size="seasonal forcing (+ baseline)",color="") + 
  theme_bw() + standard_theme  + scale_x_discrete(expand=expansion(0.1,0)) + theme(legend.position="top",
      axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x=element_text(size=16),axis.text.y=element_text(size=16),
      strip.text=element_text(size=15),legend.text=element_text(size=13),legend.title=element_text(size=14))
# save
# ggsave(paste0(foldername,"resurgence_",varname,"_2021_2019_agegr_yrly_seasfor_averaged.png"),width=35,height=25,units="cm")
# ggsave(paste0(foldername,"resurgence_",varname,"_2021_2019_agegr_yrly_seasfor_sep.png"),width=45,height=30,units="cm") 

### ### ### ### ### 
# season timing
# peak
sel_var<-"shift_seasonstart_week" # !!sym(sel_var)
df_plot_yr_seas_timing <- sum_max_ratio_by_yearly_agegr %>% filter(epi_year==2021 & !grepl(">15y",age_yr)) %>% 
  mutate(dep_type=ifelse(grepl("age",dep_type),"~age","~immunity")) %>% ungroup() %>%
  group_by(age_yr,dep_type,dep_val,R0,seasforce_peak) %>% # 
  summarise(mean_val=abs(mean(!!sym(sel_var))),min_npi_val=abs(min(!!sym(sel_var))),max_npi_val=abs(max(!!sym(sel_var))) ) %>% ungroup() %>%
  mutate(type_R0=factor(paste0(dep_type,", R0=",R0)) )
# seas timing
ggplot(df_plot_yr_seas_timing,aes(x=factor(dep_val),y=mean_val,color=type_R0,group=interaction(seasforce_peak,type_R0),size=seasforce_peak)) + #
  geom_point(position=position_dodge(width=dodge_val)) + scale_size(range=c(0.7,2)/3.5) +
  geom_pointrange(aes(ymin=min_npi_val,ymax=max_npi_val),position=position_dodge(width=dodge_val)) +
  facet_wrap(~age_yr,scales="free",nrow=3) + scale_color_manual(values=colorpal) + 
  geom_vline(xintercept=0.5+(0:8),linetype="dashed",size=1/3) + # guides(color=guide_legend(ncol=5,byrow=TRUE)) + 
  xlab("dependence on age/exp") + ylab(paste0("forward shift in season ",ifelse(grepl("start",sel_var),"start","peak")," (weeks)")) + 
  labs(size="seasonal forcing (+baseline)",color="") + theme_bw() + standard_theme  + scale_x_discrete(expand=expansion(0.1,0)) +
  theme(legend.position="top",axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),legend.title=element_text(size=14),
        axis.text.x=element_text(size=16),axis.text.y=element_text(size=16),strip.text=element_text(size=15),legend.text=element_text(size=13))
# SAVE
# ggsave(paste0(foldername,"resurgence_shift_",ifelse(grepl("start",sel_var),"start","peak"),"week_2021_2019_agegr_yrly_seasfor_sep.png"),
#       width=45,height=30,units="cm")
#
# ggsave(paste0(foldername,"resurgence_shift_",ifelse(grepl("start",sel_var),"start","peak"),"week_2021_2019_agegr_yrly_seasfor_aver.png"),
#      width=40,height=30,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# plot seas forcing
list_forcing_vector_npi=list()
for (npi_red_str in c(0:4)/10){
g(n_years,timesteps,simul_start_end,forcing_vector_npi) %=% fun_shutdown_seasforc(npi_dates,years_pre_post_npi=c(3,3),
              season_width_wks=seasforc_width_wks,init_mt_day="06-01",ifelse(grepl("exp",parsets_filtered$dep_type[k_par_filt]),45,49),
                                    forcing_above_baseline=1,npireduc_strength=0)
npi_inds=as.numeric(npi_dates[1]-simul_start_end[1]):as.numeric(npi_dates[2]-simul_start_end[1])
forcing_vector_npi[npi_inds]=1 + (forcing_vector_npi[npi_inds]-1)*npi_red_str
list_forcing_vector_npi[[1+npi_red_str*10]]=forcing_vector_npi }
# fcn_plot_seas_forc(simul_start_end,forcing_vector_npi,seas_lims_wks=c(7,42),npi_dates,date_resol="3 month")
xx=as.data.frame(list_forcing_vector_npi); colnames(xx)=paste0("NPI strength=",100*(0:4)/10)
df_seas_forc_npi=data.frame(date=seq(simul_start_end[1],simul_start_end[2],by=1),xx) %>% pivot_longer(!date) %>% 
  mutate(year=year(date),week=week(date),name=paste0(gsub("NPI.strength.","",name),"%"))
# plot
ggplot(df_seas_forc_npi %>% filter(date<=as.Date("2022-07-01")),aes(x=date,y=value,color=name)) + geom_line(size=1.05) + 
  geom_vline(data=df_seas_forc_npi %>% filter(week %in% c(7,42) & name=="1") %>% group_by(week,year) %>% filter(date==min(date)),
             aes(xintercept=date),linetype="dashed",color="black",size=1/4,show.legend=F) + labs(color="contacts during NPIs (% normal)") + 
  scale_x_date(date_breaks="2 month",expand=expansion(0.01,0)) + theme_bw() + standard_theme + xlab("") + ylab("strength of forcing") +
  geom_rect(xmin=npi_dates[1],xmax=npi_dates[2],ymin=-Inf,ymax=Inf,fill="grey",alpha=0.01,show.legend=F,color=NA) + 
  theme(axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),legend.text=element_text(size=14),legend.position="top",
        legend.title=element_text(size=16),axis.title.y=element_text(size=15))
# save
ggsave(paste0(foldername,"seasonal_forcing_NPI_contactlevel.png"),width=32,height=22,units="cm")

