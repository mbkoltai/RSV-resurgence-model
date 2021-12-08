### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Below are further plots/checks not included in the paper
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

## CHECK CORRELATIONS btwn PARAMS
ggplot(partable_filtered_AR_seasconc) +
  geom_jitter(aes(x=exp_dep,y=age_dep,color=factor(age_dep)),position=position_jitter(height=0.02,width=0.02)) +
  standard_theme + theme_bw()
ggsave(paste0(foldername,"age_exp_corr_partable_filtered_AR_seasconc.png"),width=25,height=20,units="cm")

ggplot(partable_filtered_AR_seasconc) + geom_jitter(aes(x=exp_dep,y=omega,color=factor(round(omega,4))),
                                                    position=position_jitter(height=0.0002,width=0.02)) + standard_theme + theme_bw()
ggsave(paste0(foldername,"waning_exp_corr_partable_filtered_AR_seasconc.png"),width=25,height=20,units="cm")


# PLOT attack rates, seasonal share, peak week (not a figure in manuscript)
sel_var<-c("attack_rate_perc","seas_share","max_incid_week")
estim_rates <- estim_attack_rates %>% select(agegroup_name,median_est,min_est,max_est) %>% 
  pivot_longer(!c(agegroup_name)) %>% rename(type=name) %>% mutate(name="attack_rate_perc")
estim_rates <- bind_rows(estim_rates, estim_rates %>% filter(type!="median_est") %>% 
                           mutate(name="max_incid_week",value=ifelse(grepl("min",type),3,48)))
color_var<-"exp_dep" # R0 exp_dep age_dep
ggplot(all_sum_inf_epiyear_age_filtered %>% mutate(attack_rate_perc=ifelse(epi_year==2020,NA,attack_rate_perc),
                                                   agegroup_name=factor(agegroup_name,levels=unique(agegroup_name))) %>% ungroup() %>% select(c(par_id,epi_year,
                                                                                                                                                agegroup_name,attack_rate_perc,seas_share,max_incid_week,exp_dep,age_dep,seasforc_width_wks,R0)) %>% 
         pivot_longer(!c(epi_year,agegroup_name,par_id,exp_dep,age_dep,seasforc_width_wks,R0)) ) +
  geom_hpline(aes(x=age_dep,y=value,color=get(color_var),group=par_id),width=0.1,size=1/2) +
  facet_grid(name~agegroup_name,scales="free_y") + scale_y_continuous(expand=expansion(0.02,0))+
  scale_color_gradient2(midpoint=median(c(t(unique(all_sum_inf_epiyear_age_filtered[,color_var])))),
                        low="blue",mid="white",high="red") +
  geom_hline(data=estim_rates %>% filter(!type %in% "median_est"),aes(yintercept=value),size=3/4)+
  geom_hline(data=estim_rates %>% filter(type %in% "median_est"),aes(yintercept=value),linetype="dashed",size=1/2)+ 
  xlab("age-dependence")+ylab("")+theme(legend.position="top")+theme_bw()+standard_theme+labs(color=color_var)
# save
# ggsave(paste0(foldername,"parscan_attack_rates_filtered_",color_var,".png"),width=32,height=20,units="cm")

# since the rate of waning and exposure-dependence (and inversely, age-dependence) are correlated, 
# we need to analyse the effect of waning at a given level of exp-dep
# need to find a value of exp-dep where there is a balanced number of param sets wrt waning rate
partable_regular_dyn %>% mutate(omega=1/omega) %>% group_by(exp_dep,omega) %>% summarise(n_parset=n()) %>% 
  group_by(exp_dep) %>% mutate(freq=n_parset/sum(n_parset)) %>% pivot_wider(names_from=omega,values_from=c(n_parset,freq))
# waning and R0 are also correlated!!
partable_regular_dyn %>% filter(seasforc_width_wks==3) %>% group_by(R0,omega,exp_dep) %>% 
  summarise(exp_dep=mean(exp_dep),waning=1/unique(omega),n_par=n())%>% arrange(R0,exp_dep,waning) %>% ungroup(omega) %>% 
  select(!omega)

# normalize time, calculate averages
summ_dyn_all_parsets_fixed_exp_dep <- right_join(dyn_all_parsets_broad_age, 
                                                 left_join(partable_regular_dyn %>% filter(R0==1 & exp_dep==1.625) %>% 
                                                             select(par_id,seasforc_width_wks,R0,seasforce_peak,omega),pred_pca %>% select(par_id,PC1),by="par_id"), by="par_id") %>%
  mutate(age_exp_par_bins=findInterval(PC1,seq(-1,1,by=1/5))) %>% group_by(age_exp_par_bins) %>% 
  mutate(age_exp_par_bins=round(mean(PC1),1)) %>% rename(incid_case=value) %>% pivot_longer(c(incid_case,incid_hosp)) %>% 
  rename(varname=name,varvalue=value) %>% select(!c(seasforc_width_wks,R0,seasforce_peak,PC1)) %>% 
  pivot_longer(c(omega,age_exp_par_bins)) %>% rename(parname=name,parvalue=value,value=varvalue) %>% 
  relocate(c(varname,value),.after=parvalue) %>%
  group_by(agegroup_broad,date,parname,parvalue,varname) %>% summarise(mean=mean(value),median=median(value),
                                                                       ci50_low=quantile(value,c(0.25,0.75))[1],ci50_up=quantile(value,c(0.25,0.75))[2],
                                                                       ci95_low=quantile(value,c(0.025,0.975))[1],ci95_up=quantile(value,c(0.025,0.975))[2]) %>% 
  pivot_longer(c(mean,median,ci50_low,ci50_up,ci95_low,ci95_up)) %>% 
  mutate(epi_year=ifelse(date>=paste0(year(date),"-07-01"),year(date),year(date)-1)) %>% 
  rename(metric=name) %>% relocate(epi_year,.after=date) %>% group_by(agegroup_broad,parname,parvalue,varname,metric) %>% 
  mutate(peak_2019=max(value[epi_year %in% 2019])) %>% ungroup() %>% mutate(value_norm=value/peak_2019) %>% 
  select(!peak_2019) %>%
  mutate(parname=ifelse(parname=="omega","waning","exposure (-1) <-> age (1)"),
         parvalue=ifelse(parname=="waning",1/parvalue,parvalue),
         value=round(value,1),value_norm=round(value_norm,3)) %>% filter(varname %in% "incid_hosp") %>% 
  group_by(agegroup_broad,parname,parvalue) %>% 
  mutate(max_val=max(value[(epi_year %in% 2019) & (metric %in% "median")])) %>% 
  group_by(agegroup_broad,parname,parvalue,epi_year) %>% 
  mutate(max_day_2019=as.Date(ifelse(epi_year %in% 2019,unique(date[value==max_val]),NA))) %>% 
  group_by(agegroup_broad,parname,parvalue) %>%
  mutate(max_day_2019_num=yday(as.Date(ifelse(epi_year==2019,max_day_2019,min(max_day_2019,na.rm=T)))),
         max_day_2019=as.Date(ifelse(epi_year==2019,max_day_2019,min(max_day_2019,na.rm=T)))) %>% ungroup() %>%
  mutate(max_day_year=as.Date(paste(max_day_2019_num,ifelse(year(max_day_2019)==2019,epi_year,epi_year+1)),format="%j %Y"),
         dist_peak_week=date-max_day_year,dist_peak_week_num=as.numeric(dist_peak_week)/7) %>% 
  select(!c(max_day_2019_num,max_day_2019,max_day_year))

# plot
df_plot <- summ_dyn_all_parsets_fixed_exp_dep %>% filter(agegroup_broad %in% c("<1y","1-2y","2-5y") & 
                                                           epi_year %in% c("2019","2021","2022","2023") & (metric %in% c("median","ci50_low","ci50_up")) & 
                                                           (parname %in% "waning")) %>% select(!c(varname,!!val_type_excl)) %>% 
  pivot_wider(names_from=metric,values_from=!!value_type) %>% rename(`epi-year`=epi_year)
n_par<-length(unique(df_plot$parvalue)); colorpal <- colorRampPalette(colors=c("blue","grey","red"))(n_par)
p <- ggplot(df_plot,aes(x=dist_peak_week_num,group=parvalue,color=factor(parvalue))) + geom_line(aes(y=median),size=1.1) +
  geom_ribbon(aes(ymin=ci50_low,ymax=ci50_up,fill=factor(parvalue)),color=NA,alpha=0.2) + 
  facet_grid(agegroup_broad~`epi-year`,labeller=labeller(`epi-year`=label_both),scales = "free_y") + 
  scale_x_continuous(limits=c(-23,15)) + theme_bw() + standard_theme + theme(legend.position="top",
                                                                             axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),
                                                                             axis.title.y=element_text(size=16),legend.text=element_text(size=15),legend.title=element_text(size=16),
                                                                             strip.text=element_text(size=18))+
  xlab("distance in weeks from (pre-pandemic) peak week") + labs(color=k_par,fill=k_par) + 
  ylab(paste0("weekly hospitalisations (% pre-pandemic)")) + # ",sel_agegr,ifelse(grepl("norm",value_type),", 
  geom_hline(yintercept=1,linetype="dashed",size=1/3) + geom_vline(xintercept=c(-9,9),linetype="dashed",size=1/3)
if (grepl("norm",value_type)) {p<-p+geom_hline(yintercept=1,linetype="dashed",size=1/3)}
if (n_par>3) {p <- p + scale_color_manual(values=colorpal) + scale_fill_manual(values=colorpal)}; p 
# save

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# correlations btwn params: partable_regular_dyn_corrs # partable_filtered
part_x <- partable_regular_dyn %>% select(c(exp_dep,age_dep,seasforc_width_wks,R0,seasforce_peak,omega)) %>% 
  rename(waning=omega) %>% mutate(waning=round(1/waning))
ggplot(part_x) + geom_jitter(aes(x=exp_dep,y=factor(waning),color=factor(waning)),
                             position=position_jitter(height=0.4,width=0.02)) + 
  geom_hline(yintercept = (0:3)+1/2) + scale_y_discrete(expand = expansion(0,0)) + labs(color="waning") +
  theme_bw() + standard_theme +theme(axis.title.y=element_text(size=16)) + ylab("waning period")
subfldr_name<-paste0(foldername,"median_interquant_by_param_value/param_corrs/")
if (!dir.exists(subfldr_name)) {dir.create(subfldr_name)}
ggsave(paste0(subfldr_name,"waning_exp_dep_correlation.png"),width=30,height=18,units="cm")
# all pairs
for (k in 1:sum(1:(ncol(part_x)-1))) {   col_ind<-combinations(n=ncol(part_x),r=2,v=1:ncol(part_x),repeats.allowed=F)[k,]
xx <- part_x[,col_ind] %>% mutate(pair_name=paste0(colnames(part_x[,col_ind]),collapse=" - "),
                                  x_var=colnames(part_x[,col_ind])[1],
                                  y_var=colnames(part_x[,col_ind])[2]); colnames(xx)[1:2] <- c("var1","var2")
if (k==1){param_pairs<-xx } else { param_pairs <- bind_rows(param_pairs,xx) } }
param_pairs <- param_pairs %>% group_by(pair_name) %>% mutate(var1_norm=var1/median(unique(var1)),
                                                              var2_norm=var2/median(unique(var2))) %>% group_by(pair_name) %>% 
  mutate(pair_name=paste0(pair_name," (corr=",round(cor(var1,var2),2),")"),corr=round(cor(var1,var2),2))
# plot all pairs
ggplot(param_pairs) +
  geom_jitter(aes(x=var1_norm,y=var2_norm,color=factor(var2)),position=position_jitter(height=0.05,width=0.025),
              alpha=1/2,size=3/4) +
  geom_jitter(data=param_pairs %>% group_by(pair_name,var1) %>% 
                summarise(var1_uni=unique(var1_norm),mean_var2=mean(var2_norm),corr=unique(corr)),aes(x=var1_uni,y=mean_var2)) + 
  facet_wrap(~pair_name,scales = "free") + 
  geom_rect(data=param_pairs %>% group_by(pair_name) %>% summarise(corr=unique(corr)) %>% filter(abs(corr)>0.2),
            xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,color="black",fill=NA,size=2) + theme_bw() + standard_theme + 
  theme(axis.title.y=element_text(size=16),legend.position="null") + 
  ggtitle("parameter correlations (normalised by median values)")
# plot
ggsave(paste0(subfldr_name,"param_correlations.png"),width=38,height=24,units="cm")

# distribution of param values
l_freq_table<-bind_rows(lapply(1:6, function(x) data.frame(varname=colnames(part_x)[x],t(t(table(part_x[,x])))) %>% 
                                 select(!Var2) )) %>% rename(n_par=Freq,varvalue=Var1) %>% 
  mutate(freq=round(n_par/nrow(part_x),3),varvalue=as.numeric(as.character(varvalue))) %>% 
  group_by(varname) %>% mutate(n_val=length(unique(varvalue)))
# plot frequencies
ggplot(l_freq_table) + geom_bar(aes(x=factor(varvalue),y=freq),stat="identity") + facet_wrap(~varname,scales="free") + 
  geom_hline(aes(yintercept=1/n_val),color="red") + theme_bw() + standard_theme + xlab("value of parameter") + 
  ylab("frequency")
# save
ggsave(paste0(subfldr_name,"param_freqs.png"),width=30,height=25,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Plotting individual trajectories from the parameter sampling

# selecting indiv parameter sets
ggplot(dyn_all_parsets_broad_age %>% filter(par_id %in% unique(par_id)[1:11]),
       aes(x=date,y=incid_hosp,group=par_id,color=factor(par_id))) + geom_line() + 
  facet_wrap(~agegroup_broad,scales="free_y") + scale_x_date(date_breaks="3 month",expand=expansion(0.01,0)) + 
  theme_bw() + standard_theme

# selected by parameter values
sel_inds<-c(1,7,12:15)
colorpal <- colorRampPalette(colors=c("blue","grey","red"))(length(unique(partable_regular_dyn$exp_dep)))[c(1,8,15)]
right_join(dyn_all_parsets_broad_age,
           left_join(partable_regular_dyn %>% select(par_id,seasforc_width_wks,R0,seasforce_peak,omega),pred_pca) %>% 
             mutate(omega=1/omega) %>% filter(seasforc_width_wks==5 & R0==1) ) %>% 
  filter(date>as.Date("2019-10-15") & date<as.Date("2022-04-15") & !(agegroup_broad %in% c("<1y")) & # ,"5+y"
           as.numeric(as.factor(K_exp)) %in% sel_inds) %>%
  ggplot(aes(x=date,y=incid_hosp,group=par_id,color=factor(omega))) + # ,linetype=factor(seasforce_peak)
  geom_line() + facet_grid(agegroup_broad~K_exp,scales="free_y",labeller = labeller(K_exp=label_both)) + 
  geom_vline(xintercept = as.Date("2021-09-01"),size=1/2,linetype="dashed") + labs(color="exp_dep") + xlab("") +
  scale_x_date(date_breaks="3 month",expand=expansion(0.01,0)) + theme_bw() + standard_theme + 
  theme(legend.position="top",legend.text=element_text(size=14),strip.text=element_text(size=14))
# + scale_color_manual(values=colorpal)
ggsave("repo_data/dynamics/waning_effect_grid_expdep_agegroup.png",width=30,height=16,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# BIN 
#
# norm_seas_length_wk <- round(length(as.Date("2020-10-07"):as.Date("2021-03-03"))/7,1)
# hosp_sum_prob_broad_agegr <- hosp_probabilities %>% mutate(agegroup_broad=c("<1y","1-2y","2-5y","5+y")[
#     findInterval(factor(agegroup_name,levels=unique(agegroup_name)),c(2,4,7)+1)+1]) %>%
#   group_by(agegroup_broad) %>% summarise(hosp_sum=sum(hosp_num_from_per_inf_prob)) %>% 
#   mutate(hosp_per_week_season=hosp_sum/norm_seas_length_wk)



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# cumul_peak_meanage_hosp_byvalue <- bind_rows(
#   # cumulative hospitalisations
#   parsets_broad_age_groups %>% filter(name %in% "hosp_tot" & epi_year>2020) %>% ungroup() %>% rename(varname=name) %>% 
#     select(!c(value,PC1,par_id)) %>% pivot_longer(c(age_exp_par_bins,seasforce_peak,waning,seasforc_width_wks,R0)) %>% 
#     rename(parname=name,parvalue=value) %>% relocate(value_norm,.after=parvalue) %>% 
#     group_by(agegroup_broad,epi_year,parname,varname) %>%
#     summarise(mean_val=mean(value_norm),med_val=median(value_norm),
#               min_val=mean(value_norm[parvalue==ifelse(parname %in% c("waning","age_exp_par_bins"),
#                                                        max(parvalue),min(parvalue))]),
#               max_val=mean(value_norm[parvalue==ifelse(parname %in% c("waning","age_exp_par_bins"),
#                                                        min(parvalue),max(parvalue))])),
#   # peak hospitalisations
#   parsets_max_incid_seas_length %>% filter(varname %in% "incid_hosp" & vartype %in% "max_value" & epi_year>2020) %>% 
#     ungroup() %>% 
#     select(!c(value,PC1,par_id)) %>% pivot_longer(c(age_exp_par_bins,seasforce_peak,waning,seasforc_width_wks,R0)) %>% 
#     rename(parname=name,parvalue=value) %>% relocate(value_norm,.after=parvalue) %>% 
#     group_by(agegroup_broad,epi_year,parname,varname) %>%
#     summarise(mean_val=mean(value_norm),med_val=median(value_norm),min_val=mean(value_norm[parvalue==min(parvalue)]),
#               max_val=mean(value_norm[parvalue==max(parvalue)])),
#   # mean age
#   parsets_mean_age_inf %>% filter(name %in% "mean_age_hosp_tot_under_5" & epi_year>2020) %>% ungroup() %>%
#     select(!c(value,PC1,par_id)) %>% rename(varname=name) %>% 
#     pivot_longer(c(age_exp_par_bins,seasforce_peak,waning,seasforc_width_wks,R0)) %>%
#     rename(parname=name,parvalue=value) %>% relocate(value_norm,.after=parvalue) %>% 
#     group_by(epi_year,parname,varname) %>%
#     summarise(mean_val=mean(value_norm),med_val=median(value_norm),min_val=mean(value_norm[parvalue==min(parvalue)]),
#               max_val=mean(value_norm[parvalue==max(parvalue)])) ) %>% filter(!agegroup_broad %in% "5+y")


####
# create merged dynamics from the component files
# left_join(bind_rows(
# lapply(list.files(dyn_folder,pattern="dyn_parsets.*csv"), 
#      function(x) read_csv(file=paste0(dyn_folder,x)) %>% # filter(par_id %in% parsets_regular_dyn$par_id) %>% 
#     mutate(date=t-min(t)+as.Date(start_date_dyn_save)) %>% select(!c(t,name)) %>% 
#     group_by(agegroup,date,par_id) %>% # filter(date>=as.Date("2018-10-01") & date<=as.Date("2024-07-01")) %>% 
#     summarise(value=sum(value)) %>% group_by(par_id,agegroup) %>% 
#     mutate(value=round(roll_sum(value,n=7,fill=NA,align="right",by=7))) %>% filter(!is.na(value))  )  ), 
#     hosp_probabilities %>% mutate(agegroup=as.numeric(factor(agegroup_name,levels=unique(agegroup_name)))) %>%
# select(agegroup,prob_hosp_per_infection_adj),by="agegroup") %>% mutate(incid_hosp=value*prob_hosp_per_infection_adj,
# agegroup_broad=c("<1y","1-2y","2-5y","5+y")[findInterval(agegroup,c(2,4,7)+1)+1]) %>% 
# group_by(agegroup_broad,date,par_id) %>% summarise(value=sum(value),incid_hosp=sum(incid_hosp))

# normalise time wrt peak week pre-NPI 
summ_dyn_all_parsets_norm_time <- summ_dyn_all_parsets_broad_age %>% filter(varname %in% "incid_hosp") %>% 
  group_by(agegroup_broad,parname,parvalue) %>% 
  mutate(max_val=max(value[(epi_year %in% 2019) & (metric %in% "median")])) %>% 
  group_by(agegroup_broad,parname,parvalue,epi_year) %>% 
  mutate(max_day_2019=as.Date(ifelse(epi_year %in% 2019,unique(date[value==max_val]),NA))) %>% 
  group_by(agegroup_broad,parname,parvalue) %>% 
  mutate(max_day_2019_num=yday(as.Date(ifelse(epi_year==2019,max_day_2019,min(max_day_2019,na.rm=T)))),
         max_day_2019=as.Date(ifelse(epi_year==2019,max_day_2019,min(max_day_2019,na.rm=T)))) %>% ungroup() %>%
  mutate(max_day_year=as.Date(paste(max_day_2019_num,
                                    ifelse(year(max_day_2019)==2019,epi_year,epi_year+1)),format="%j %Y"),
         dist_peak_week=date-max_day_year,dist_peak_week_num=as.numeric(dist_peak_week)/7) %>% filter(!is.na(parvalue))

# calculate statistics at abs time points
# summ_dyn_all_parsets_broad_age_abs_time <- dyn_all_parsets_broad_age_params %>% 
#   group_by(agegroup_broad,date,parname,parvalue,varname) %>%
#       summarise(mean=mean(value_norm),median=median(value_norm),
#             ci50_low=quantile(value_norm,c(0.25,0.75))[1],ci50_up=quantile(value_norm,c(0.25,0.75))[2],
#             ci95_low=quantile(value_norm,c(0.025,0.975))[1],ci95_up=quantile(value_norm,c(0.025,0.975))[2]) %>% 
#   mutate(epi_year=ifelse(week(date)>=epi_year_week_start,year(date),year(date)-1)) %>% 
#   pivot_longer(c(mean,median,ci50_low,ci50_up,ci95_low,ci95_up)) %>% rename(metric=name) %>% 
#   relocate(epi_year,.after=date) %>% group_by(agegroup_broad,parname,parvalue,varname,metric) %>% 
#   mutate(parname=ifelse(parname=="omega","waning","exposure (-1) <-> age (1)"),
#          parvalue=ifelse(parname=="waning",1/parvalue,parvalue),value=round(value,3))