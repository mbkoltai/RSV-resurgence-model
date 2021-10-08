rm(list=ls()); currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
# load constant parameters and functions
source("load_params.R")
# estimated attack rates
estim_attack_rates <- data.frame(agegroup_name=rsv_age_groups$agegroup_name, # paste0("age=",,"yr")
  median_est=c(rep(65,4),rep(40,4),10,8,5)) %>% mutate(min_est=median_est*0.25,max_est=median_est*2.5,
  median_all_inf=c(rep(70,4),rep(60,4),50,30,20),min_est_all_inf=median_all_inf*0.5,max_est_all_inf=median_all_inf*1.5)
# % cases within season (filtering parameter sets)
seas_conc_lim=0.8
# parameter sets to search through
selsets<-c(2,4:8) 
p_table <- bind_rows(expand.grid(list(dep_type="age",dep_val=seq(0.5,4,by=0.5)[selsets],R0=seq(12,14,0.5)/10,
                                      seasforce_peak=seq(1.1,1.5,by=0.1))),
              bind_rows(lapply(selsets, function(x) expand.grid(list(dep_type="exp",dep_val=x,
              seasforce_peak=list(c(0.8,1,1.2),c(1.125,1.25,1.375),c(1,1.25,1.5),c(1,1.25,1.5),
              c(1,1.125,1.25),c(1.25,1.375,1.5),c(1.25,1.375,1.5),c(1.25,1.375,1.5))[[x]],
              R0=list( (12:14)/10,(12:14)/10,seq(13,16,1.5)/10,seq(15,18,1.5)/10, seq(14,20,2)/10, 
              seq(15,18,1)/10,(18:20)/10,(24:26)/10)[[x]]))))) %>% arrange(dep_type,dep_val,R0,seasforce_peak) %>% 
  rowid_to_column("par_id")
# p_table <- read_csv("simul_output/parscan/sel_parsets/sel_parsets.csv") %>% rowid_to_column("par_id")
partable <- fcn_create_partable(p_table,nstep=10, scale_age_exp=c(0.35,0.29),pop_struct=rsv_age_groups$stationary_popul,
                                susc_denomin=100,susc_min=0.11,nage=11,ninf=3,rhoval=rho) %>% 
  mutate(dep_val=ifelse(dep_type=="age",dep_val*2,dep_val)) %>%
  group_by(dep_type,R0) %>% mutate(R0_no=row_number())
# save the stat sol of all param sets
stat_sol_allparsets=matrix(0,nrow=(n_compartment+1)*n_age*n_inf,ncol=nrow(partable))
# NPI dates
npi_dates=as.Date(c("2020-03-26","2021-05-17"))
# width of season (from peak)
seasforc_width_wks<-8
# length of simulations
simul_length_yr=7
# agegroup indices for maternal immunity
mat_imm_flag <- TRUE
mat_imm_inds<-list(fun_sub2ind(i_inf=1,j_age=1,"R",c("S","I","R"),11,3),
                   fun_sub2ind(i_inf=c(1,2,3),j_age=9,"R",c("S","I","R"),11,3),
                   fun_sub2ind(i_inf=c(1,2,3),j_age=9,"S",c("S","I","R"),11,3))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# parallelisation (write file that'll run scripts)
n_core=7; initcond_file <- "simul_output/parscan/parallel/initconds_all.csv"
system(paste0(c("Rscript write_run_file.R",n_core,nrow(partable),simul_length_yr,3,initcond_file),collapse=" "))
# run calculation
system("sh run_all_parallel_scan.sh")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
results_folder<-"simul_output/parscan/parallel/"
# collect all results
dyn_parsets_all <- bind_rows(lapply(list.files(path=results_folder,pattern="dyn_parsets*"),
                 function(x) read_csv(paste0(parall_foldername,x))))
summ_parsets_all <- bind_rows(lapply(list.files(path=results_folder,pattern="summ_parsets*"),
                                    function(x) read_csv(paste0(parall_foldername,x))))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# sneak peak of dynamics
ggplot(dyn_parsets_all %>% filter(par_id>=101 & par_id<=155 & date>as.Date("2019-10-01") & date<as.Date("2022-10-01") & agegroup<=3)) + 
  geom_line(aes(x=date,y=value,color=par_id)) + 
  facet_grid(infection~agegroup,scales="free_y",labeller=labeller(infection=label_both,agegroup=label_both)) +
  scale_color_gradientn(colours=wes_palette("Zissou1",55, type = "continuous")) +
  theme_bw() + standard_theme + scale_x_date(date_breaks="year") + xlab("") + ylab("") + labs(color="# infection")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# plot attack rates by age group and years
check_crit=10/11; n_sel_yr=4; peak_week_lims <- c(47,2)
all_sum_inf_epiyear_age_filtered <- all_sum_inf_epiyear_age %>% filter(epi_year>=2016 & epi_year<=2019) %>% 
  group_by(seasforce_peak,dep_val,R0,dep_type) %>% filter(sum(attack_rate_check)>=round(n_age*n_sel_yr*check_crit) & 
    sum(seas_share_check)>=round(n_age*n_sel_yr*check_crit) & (max_incid_week>=peak_week_lims[1]|max_incid_week<=peak_week_lims[2]) )
# plot
uni_dep_vals=(all_sum_inf_epiyear_age_filtered %>% group_by(dep_type) %>% summarise(n=length(unique(R0))))$n
colorpal=c(colorRampPalette(colors=c("orange","red"))(uni_dep_vals[1]),colorRampPalette(colors=c("grey","black"))(uni_dep_vals[2]))
# PLOT
ggplot(all_sum_inf_epiyear_age_filtered %>% mutate(attack_rate_perc=ifelse(epi_year==2020,NA,attack_rate_perc),
  agegroup_name=factor(agegroup_name,levels=unique(agegroup_name)),dep_type_R0=factor(paste0(dep_type,", R0=",R0),
  levels=unique(paste0(all_sum_inf_epiyear_age_filtered$dep_type,", R0=",all_sum_inf_epiyear_age_filtered$R0)))) ) +
  geom_hpline(aes(x=factor(epi_year),y=attack_rate_perc,color=dep_type_R0),width=0.9,size=1/2) + scale_color_manual(values=colorpal) + 
  facet_grid(dep_val~agegroup_name,scales="free_y",labeller=labeller(seasforce_peak=label_both)) + # ,linetype=dep_type
  geom_hline(data=estim_attack_rates %>% pivot_longer(!agegroup_name) %>% filter(name!="median_est"&!grepl("all_inf",name)),
             aes(yintercept=value),linetype="dashed",size=1/3) +
  geom_hline(data=estim_attack_rates %>% pivot_longer(!agegroup_name) %>% filter(name=="median_est"&!grepl("all_inf",name)),
     aes(yintercept=value),size=1/3) + ylab("attack rate % age group") + scale_y_log10() + theme(legend.position="top") +
  theme_bw() + standard_theme + labs(color="") + xlab("") 
# ,caption=paste0("model check minimum=",round(check_crit,2)*1e2,"%")

