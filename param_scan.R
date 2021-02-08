# this script is a 2D parscan to explore the spectrum of dynamical behaviors ~ f(basal_rate,onset/offset_NPI)
### define constant parameters ----
basal_rate_vals=c(0.08,0.1,0.15,0.2,0.21,0.225,0.25); npi_delta_vals=c(-2,0,1,2,3,6)
parscan_df_ode_sol_cases_sum=data.frame(); parscan_season_peaks_AUC=data.frame(); parscan_mean_age_perseason=data.frame()
shutdwn_lims_plot=data.frame()
list_symptom_agegroups=list(1:2,3:7,8:9,10:11)
prop_symptom=1-c(mean(rbeta(1e3,shape1=3,shape2=30)),mean(rbeta(1e3,shape1=9,shape2=43)),
                 mean(rbeta(1e3,shape1=38,shape2=35)),mean(rbeta(1e3,shape1=36,shape2=11)))
# EXPOSURE DEPENDENT or not? expos_dep_val=0 -> no dependence on exposure. if expos_dep_val>0 -> dependence on expos)
df_symptom_prop=fun_propsymptom_table(list_symptom_agegroups,expos_dep_val=1,rsv_age_groups$agegroup_name,n_inf=3)

### PARSCAN ----
for (k_basal in 1:length(basal_rate_vals)) {
  for (k_npi in 1:length(npi_delta_vals)) {
# basal rate of seas forcing
basal_rate=basal_rate_vals[k_basal]
# seasonal parameters
peak_week=49; season_width=3.8 # (in weeks)
# NPI onset/end
npi_year=4; preseas_npi_on=npi_delta_vals[k_npi]; postseas_npi_off=preseas_npi_on # shutdown_scale=0.2
shutdown_list=fun_shutdown_seasforc(timesteps,elem_time_step,basal_rate,npi_year,peak_week,season_width,
                                    preseas_npi_on,postseas_npi_off,n_prec=0.01,st_devs=2:3,n_sd=2)
forcing_vector_npi=shutdown_list[[1]]; shutdwn_lims=shutdown_list[[2]];seas_force=shutdown_list[[3]]; seas_lims=shutdown_list[[4]]; 
shutdown_scale=forcing_vector_npi[round(mean(shutdwn_lims))]
# season limits
shutdwn_lims=data.frame(t(shutdwn_lims)); shutdwn_lims[,"k_basal"]=basal_rate_vals[k_basal]; 
shutdwn_lims[,"k_npi"]=npi_delta_vals[k_npi]; colnames(shutdwn_lims)[1:2]=c("on","off")

# simulation
params=list(birth_term,K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds,forcing_vector_npi,elem_time_step)
# SIMULATE
ptm<-proc.time(); ode_solution<-lsoda(initvals_sirs_model,timesteps,func=sirs_seasonal_forc,parms=params); proc.time()-ptm
print(c(k_basal,k_npi))
# reshape data
list_simul_output=fun_process_simul_output(ode_solution,varname_list,n_age,n_inf,rsv_age_groups)
# df_ode_solution=list_simul_output[[1]]
df_ode_solution_tidy=list_simul_output[[2]]; rm(list_simul_output) 
# View(round(df_ode_solution[,grepl("I|t",colnames(df_ode_solution))],2))
# View(subset(df_ode_solution_tidy,compartment %in% "I"))

# symptomatic cases (currently fixed)
df_ode_sol_cases_sum=fun_symptomcases_table(df_ode_solution_tidy,df_symptom_prop, bindcolnames=c("infection","agegroup"),n_sd=2)
df_ode_sol_cases_sum[,"k_basal"]=basal_rate_vals[k_basal]; df_ode_sol_cases_sum[,"k_npi"]=npi_delta_vals[k_npi]
# age distribution
season_peaks_AUC=fcn_seas_agedistrib(df_ode_sol_cases_sum,max_time,timesteps,seas_case_threshold=5e2)
season_peaks_AUC[,"k_basal"]=basal_rate_vals[k_basal]; season_peaks_AUC[,"k_npi"]=npi_delta_vals[k_npi]

# mean age per season
mean_age_perseason=season_peaks_AUC %>% group_by(season) %>% 
  summarise(season_mean_age=sum(auc_case_share_season*mean_age*agegr_size/sum(agegr_size)),pre_post_shtdn=unique(pre_post_shtdn))
mean_age_perseason[,"k_basal"]=basal_rate_vals[k_basal]; mean_age_perseason[,"k_npi"]=npi_delta_vals[k_npi]

# concatenate
shutdwn_lims_plot=bind_rows(shutdwn_lims_plot,shutdwn_lims)
parscan_df_ode_sol_cases_sum=bind_rows(parscan_df_ode_sol_cases_sum,df_ode_sol_cases_sum)
parscan_season_peaks_AUC=bind_rows(parscan_season_peaks_AUC,season_peaks_AUC)
parscan_mean_age_perseason=bind_rows(parscan_mean_age_perseason,mean_age_perseason)
}
}
### end of parscan

# SAVE/LOAD
write_csv(parscan_df_ode_sol_cases_sum,"simul_output/parscan/parscan_df_ode_sol_cases_sum.csv")
write_csv(parscan_season_peaks_AUC,"simul_output/parscan/parscan_season_peaks_AUC.csv")
# parscan_season_peaks_AUC=read_csv("simul_output/parscan/parscan_season_peaks_AUC.csv")
write_csv(parscan_mean_age_perseason,"simul_output/parscan/parscan_mean_age_perseason.csv")
# parscan_mean_age_perseason=read_csv("simul_output/parscan/parscan_mean_age_perseason.csv")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### plot symptom cases for all param values -----
k_plot=2
g(scale_val,y_axis_tag,plotvar,foldername,full_filename,caption_txt) %=% fun_sumcase_plot_tags(n_val=2,r_val=2,k_plot=k_plot,
                                  df_symptom_prop,preseas_npi_on,postseas_npi_off,shutdown_scale,basal_rate)
# plot
ggplot(subset(parscan_df_ode_sol_cases_sum,t_years>xval_lims[1] & t_years<(max_time/365-0.05) & agegroup<=6 & agegroup>1),
       aes(x=t_years,y=get(plotvar),group=agegroup_name,color=agegroup_name)) + geom_line() + theme_bw() + standard_theme + 
  theme(axis.text.x=element_text(size=7,vjust=0.5),axis.text.y=element_text(size=6)) + 
  facet_grid(k_basal~k_npi,scales=scale_val,labeller=labeller(k_basal=label_both,k_npi=label_both)) + 
  scale_x_continuous(breaks=xval_breaks,minor_breaks=seq(0,max_time/365,by=1/12)) + # scale_y_continuous(breaks=(0:10)/10) + # 
  scale_color_manual(values = wes_palette(n=5, "Zissou1")) + # scale_color_gradientn(colours=wes_palette("Zissou1")) +
  geom_rect(data=shutdwn_lims_plot,aes(xmin=on/365,xmax=off/365,ymin=0,ymax=Inf),fill="grey",color=NA,inherit.aes=FALSE,alpha=0.5) +
  geom_vline(data=subset(seas_lims,sd==n_sd&on>xval_lims[1]&off<xval_lims[2]),aes(xintercept=on),color="blue",linetype="dashed",size=0.3) +
  geom_vline(data=subset(seas_lims,sd==n_sd&on>xval_lims[1]&off<xval_lims[2]),aes(xintercept=off),color="turquoise4",linetype="dashed",size=0.3) +
  xlab('years') + ylab(y_axis_tag) + ggtitle('Symptomatic cases by age group')
# save
ggsave(paste0("simul_output/parscan/parscan_cases_sum",c("_fraction","_absval")[k_plot],".png"),
       width=31,height=22,units="cm")

### plot symptom cases for selected param values -----
sel_scen_inds=t(data.frame(c(1,2),c(2,6),c(2,3),c(3,6),c(5,2),c(6,4))); rownames(sel_scen_inds)=c()
sel_scen_inds=data.frame(k_basal=basal_rate_vals[sel_scen_inds[,1]],k_npi=npi_delta_vals[sel_scen_inds[,2]])
sel_scenarios=right_join(subset(parscan_df_ode_sol_cases_sum,t_years>xval_lims[1] & t_years<(max_time/365-0.05) & agegroup<=6 & agegroup>1),
                         sel_scen_inds,by=c("k_basal","k_npi"))
sel_scenarios[,"scen_name"]=paste0("basal activity=",sel_scenarios$k_basal,", timing(NPI-season)=",sel_scenarios$k_npi,"w")
goodorder=match(array(t(sel_scen_inds %>% unite(x,everything()))), 
      array(t(unique(sel_scenarios[,c("k_basal","k_npi")]) %>% unite(x,everything()))) )
sel_scenarios$scen_name=factor(sel_scenarios$scen_name,levels=unique(sel_scenarios$scen_name)[goodorder])
# sel fraction/absval
k_plot=2
g(scale_val,y_axis_tag,plotvar,foldername,full_filename,caption_txt) %=% fun_sumcase_plot_tags(n_val=2,r_val=2,k_plot=k_plot,
                                                        df_symptom_prop,preseas_npi_on,postseas_npi_off,shutdown_scale,basal_rate)
ggplot(sel_scenarios,aes(x=t_years,y=get(plotvar),group=agegroup_name,color=agegroup_name)) + geom_line() + 
  theme_bw() + standard_theme + theme(axis.text.x=element_text(size=7,vjust=0.5),axis.text.y=element_text(size=6)) + 
  # facet_grid(k_basal~k_npi,scales=scale_val,labeller=labeller(k_basal=label_both,k_npi=label_both)) + 
  facet_wrap(~scen_name,scales=scale_val,ncol=2) + #,labeller=labeller(k_basal=label_both,k_npi=label_both)
  scale_x_continuous(breaks=xval_breaks,minor_breaks=seq(0,max_time/365,by=1/12)) + # scale_y_continuous(breaks=(0:10)/10) + # 
  scale_color_manual(values = wes_palette(n=5, "Zissou1")) + # scale_color_gradientn(colours=wes_palette("Zissou1")) +
  geom_rect(data=right_join(shutdwn_lims_plot, unique(sel_scenarios[,c("k_basal","k_npi")]),by=c("k_basal","k_npi")),
            aes(xmin=on/365,xmax=off/365,ymin=0,ymax=Inf),fill="grey",color=NA,inherit.aes=FALSE,alpha=0.05) +
  geom_vline(data=subset(seas_lims,sd==n_sd&on>xval_lims[1]&off<xval_lims[2]),aes(xintercept=on),color="blue",linetype="dashed",size=0.3) +
  geom_vline(data=subset(seas_lims,sd==n_sd&on>xval_lims[1]&off<xval_lims[2]),aes(xintercept=off),color="turquoise4",linetype="dashed",size=0.3) +
  xlab('years') + ylab(y_axis_tag) + ggtitle('Symptomatic cases by age group')
# SAVE
ggsave(paste0("simul_output/parscan/parscan_cases_sum",c("_fraction","_absval")[k_plot],"_SEL_SCEN.png"),
       width=31,height=22,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### plot age distributions ALL params ---
plot_season_agedist=subset(parscan_season_peaks_AUC, agegroup<=7 & season>=min(3,npi_year-1))
selvar="auc_case_share_season"; ylimvals=c(min(plot_season_agedist[,selvar]),max(plot_season_agedist[,selvar]))
# if (length(unique(plot_season_agedist$pre_post_shtdn))==3) {linetype_vals=c("dashed","dotted","solid")} else {linetype_vals=c("dotted","solid")}
# PLOT
ggplot(plot_season_agedist,aes(x=agegroup_name,y=get(selvar),group=season)) +
  geom_line(aes(linetype=pre_post_shtdn,color=factor(season)),size=1.05) + # ,size=pre_post_shtdn
  geom_point(aes(shape=pre_post_shtdn,color=factor(season)),size=1.5) + theme_bw() + standard_theme + 
  theme(axis.text=element_text(size=7),legend.text=element_text(size=12),legend.title=element_text(size=14)) +
  xlab("") + ylab("fraction of all symptomatic cases in season") + # facet_wrap(~season,ncol=3) + 
  scale_linetype_manual(values=c("solid","dotdash","dashed"))+scale_shape_manual(values=c(0,2,1)) + 
  facet_grid(k_basal~k_npi,labeller=labeller(k_basal=label_both,k_npi=label_both)) + # scale_size_manual(values=c(1.75,1.25,1.25)) +
  # scale_y_continuous(breaks=seq(round(ylimvals[1],2),round(ylimvals[2],2),by=0.02),limits=ylimvals) +
  labs(size="pre/post-NPI",shape="pre/post-NPI",color="season",linetype="pre/post-NPI",
       title="Age distribution before and after NPI in year 5")
# SAVE
ggsave("simul_output/parscan/parscan_peaks_AUC.png",width=31,height=22,units="cm")

### plot age distributions SELECTED params ---
plot_season_agedist=right_join(plot_season_agedist,sel_scen_inds,by=c("k_basal","k_npi"))
plot_season_agedist[,"scen_name"]=paste0("basal activity=",plot_season_agedist$k_basal,", timing(NPI-season)=",plot_season_agedist$k_npi,"w")
# goodorder=match(array(t(sel_scen_inds %>% unite(x,everything()))),
#                 array(t(unique(plot_season_agedist[,c("k_basal","k_npi")]) %>% unite(x,everything()))) )
plot_season_agedist$scen_name=factor(plot_season_agedist$scen_name,levels=unique(plot_season_agedist$scen_name)[goodorder])

ggplot(subset(plot_season_agedist,season>3),aes(x=agegroup_name,y=get(selvar),group=season)) +
  geom_line(aes(linetype=pre_post_shtdn,color=factor(season)),size=1.05) + # ,size=pre_post_shtdn
  geom_point(aes(shape=pre_post_shtdn,color=factor(season)),size=2.5) + theme_bw() + standard_theme + 
  # scale_color_manual(values=wes_palette(n=5, "Zissou1")) + # scale_color_gradientn(colours=wes_palette("Zissou1")) +
  theme(axis.text=element_text(size=11),legend.text=element_text(size=12),legend.title=element_text(size=14)) +
  xlab("") + ylab("fraction of all symptomatic cases in season") + # facet_wrap(~season,ncol=3) + 
  scale_linetype_manual(values=c("solid","dotdash","dashed"))+scale_shape_manual(values=c(0,2,1)) + 
  facet_wrap(~scen_name,ncol=2) + # scales=scale_val,scale_size_manual(values=c(1.75,1.25,1.25)) +
  # scale_y_continuous(breaks=seq(round(ylimvals[1],2),round(ylimvals[2],2),by=0.02),limits=ylimvals) +
  labs(size="pre/post-NPI",shape="pre/post-NPI",color="season",linetype="pre/post-NPI",
       title="Age distribution before and after NPI in year 5")
# SAVE
ggsave("simul_output/parscan/parscan_season_AUC_SELSCEN.png",width=31,height=22,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### plot mean age ALL ----------------
agemax=3; meanage_max=ceiling(max(subset(parscan_mean_age_perseason,season>agemax)[,"season_mean_age"]))+0.5
ggplot(subset(parscan_mean_age_perseason,season>agemax),aes(x=season,y=season_mean_age,fill=pre_post_shtdn)) + 
  geom_bar(stat="identity",color="black",width =0.7) + facet_grid(k_basal~k_npi,labeller=labeller(k_basal=label_both,k_npi=label_both)) +
  # scale_y_continuous(breaks=seq(0,meanage_max,0.5)) + # scale_x_continuous(breaks=0:max(mean_age_perseason$season)) + 
  geom_rect(aes(xmin=ceiling(shutdwn_lims[1]/365)+0.6,xmax=ceiling(shutdwn_lims[1]/365)+1.4,ymin=0,ymax=meanage_max),
            fill="pink",color=NA,alpha=0.1) + theme_bw() + standard_theme + theme(legend.title=element_text("pre/post NPI"),
            axis.text.x=element_text(size=9,vjust=0.5)) + ylab("mean age of cases") + # xlab("season") + 
  geom_text(aes(label=round(season_mean_age,1)),hjust=-0.3,color="black",size=3.5) + coord_cartesian(ylim=c(0.2,meanage_max-0.2)) +
  coord_flip() 
# SAVE
ggsave("simul_output/parscan/parscan_mean_age_perseason.png",width=31,height=22,units="cm")

### plot mean age SELECTED ----------------
parscan_mean_age_perseason[,"scen_name"]=paste0("basal activity=",parscan_mean_age_perseason$k_basal,
                                                ", timing(NPI-season)=",parscan_mean_age_perseason$k_npi,"w")
# goodorder=match(array(t(sel_scen_inds %>% unite(x,everything()))),
#                 array(t(unique(plot_season_agedist[,c("k_basal","k_npi")]) %>% unite(x,everything()))) )
parscan_mean_age_perseason=right_join(parscan_mean_age_perseason,sel_scen_inds,by=c("k_basal","k_npi"))
parscan_mean_age_perseason$scen_name=factor(parscan_mean_age_perseason$scen_name,
                                            levels=unique(parscan_mean_age_perseason$scen_name)[goodorder])

ggplot(subset(parscan_mean_age_perseason,season>agemax),aes(x=season,y=season_mean_age,fill=pre_post_shtdn)) + 
  geom_bar(stat="identity",color="black",width =0.7) + facet_wrap(~scen_name,ncol=2) +
  geom_rect(aes(xmin=ceiling(shutdwn_lims[1]/365)+0.6,xmax=ceiling(shutdwn_lims[1]/365)+1.4,ymin=0,ymax=meanage_max),
            fill="pink",color=NA,alpha=0.1) + theme_bw() + standard_theme + scale_y_continuous(breaks=0:7) +
  theme(legend.title=element_text("pre/post NPI"),axis.text.x=element_text(size=11,vjust=0.5)) + ylab("mean age of cases") + 
  geom_text(aes(label=round(season_mean_age,1)),hjust=-0.3,color="black",size=4) + coord_cartesian(ylim=c(0.2,meanage_max-0.2)) +
  coord_flip() 
# SAVE
ggsave("simul_output/parscan/parscan_mean_age_perseason_SEL_SCEN.png",width=31,height=22,units="cm")