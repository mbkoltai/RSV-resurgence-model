# plots of cumul/peak hosp normalised to pre-NPI by values of parameters, 
# as calculated from dynamics

sel_pars <- unique(summ_dyn_peak_cumul_meanage_byparvalue$parname)
sel_varnames <- unique(summ_dyn_peak_cumul_meanage_byparvalue$varname)

start_year=2021
all_plots=F
# save only peak hospitalisations plot

n_cnt=0 
for (k_epi_year_start in 1){ # 1:2
  for (k_plot_par in 2) { # 1:length(sel_pars)
    for (k_plot_var in 2) { # 1:length(sel_varnames)
      if(all(c(k_epi_year_start,k_plot_par,k_plot_var)==1)) {n_cnt=0} else {n_cnt=n_cnt+1}
      sel_var <- sel_varnames[k_plot_var]; sel_par <- sel_pars[k_plot_par]; dodge_val=1
      epi_year_week_start=c(23,40)[k_epi_year_start]
      # sel data to plot
      for (k_start_yr in 2021) { # start_year=2021
      df_plot <- summ_dyn_peak_cumul_meanage_byparvalue %>% 
        filter(epi_year>=start_year & epi_year<=2023 & (varname %in% sel_var) & (parname %in% sel_par) &
                 !(agegroup %in% 5) & epi_year_wk_start %in% epi_year_week_start) %>%
        mutate(varname=case_when(grepl("peak_hosp",varname) ~ "peak hospitalisations",
                                 grepl("peak_yday",varname) ~ "peak hospitalisation (day of year)", 
                                 grepl("cumul_hosp",varname) ~ "cumulative hospitalisation",
                                 grepl("mean_age_under5y",varname) ~ "shift in mean age")) %>% 
        mutate(parname=case_when(grepl("age_exp_PC1",parname) ~ "exposure (-1) <-> age (1)", 
                                 grepl("age_exp_ratio",parname) ~ "ratio of age- \nto exposure-dependence", 
                                 grepl("age_dep",parname) ~ "age-dependence (1/15 to 1/3)",
                                 grepl("exp_dep",parname) ~ "exposure-dependence (1/3 to 5/4)", 
                                 grepl("seasforc_width_wks",parname) ~ "season width (weeks)",
                                 grepl("seasforce_peak",parname) ~ "seasonal forcing (% above baseline)",
                                 grepl("peak_week",parname) ~ "peak week",
                                 grepl("R0",parname) ~ "R0 (baseline)",
                                 grepl("waning",parname) ~ "waning (1/day)"),
               agegroup=case_when(agegroup==1 ~ "0-1y", agegroup==2 ~ "1-2y",agegroup==3 ~ "2-5y"))
      # plot params
      ylab_tag <- paste0(unique(df_plot$varname), ifelse(grepl("yday|age",unique(df_plot$varname)),
                                                         " (change from pre-NPI)"," (relative to pre-NPI)"))
      n_par_value <- length(unique(df_plot$parvalue)); 
      median_width=1/100
      # colour palette
      if (!grepl("age|exp",sel_par)){
        colorpal=colorRampPalette(colors=c("orange","red"))(n_par_value) } else  {
          colorpal=colorRampPalette(colors=c("blue","grey","red"))(n_par_value); print("blue/red") }
      # create plot
      width_scale=1.32
      p <- ggplot(df_plot,aes(x=factor(epi_year),color=factor(signif(parvalue,2)),group=parvalue)) + 
        geom_linerange(aes(ymin=ci95_l,ymax=ci95_u),position=position_dodge(width=dodge_val),
                       alpha=0.3,size=width_scale*ifelse(any(grepl("mean",df_plot$varname)),72,24)/n_par_value) +
        geom_linerange(aes(ymin=ci50_l,ymax=ci50_u),position=position_dodge(width=dodge_val),
                       alpha=0.6,size=width_scale*ifelse(any(grepl("mean",df_plot$varname)),72,24)/n_par_value) +
        geom_linerange(aes(x=factor(epi_year),ymin=median-median_width,ymax=median+median_width),
                       position=position_dodge(width=dodge_val),color="black",
                       size=width_scale*ifelse(any(grepl("mean",df_plot$varname)),72,24)/n_par_value) + 
        geom_vline(xintercept=(0:4)+1/2,size=1/5) + # geom_hline(yintercept=1,linetype="dashed",size=1/2) +
        xlab("") + ylab(ylab_tag) + labs(color=unique(df_plot$parname)) +
        scale_x_discrete(expand=expansion(0.02,0)) + scale_color_manual(values=colorpal) + 
        theme_bw() + standard_theme + theme(strip.text=element_text(size=15),axis.text.x=element_text(size=13),
                          axis.text.y=element_text(size=12),legend.text=element_text(size=11),
                          legend.title=element_text(size=12),
                          legend.position="bottom") + guides(fill=guide_legend(ncol=3)) +
        manuscript_large_font_theme
      # ifelse(!grepl("exp|age|kappa",
      # unique(df_plot$parname)),"bottom","right")
      
      if (!grepl("mean_age",sel_var)) { p <- p  + facet_wrap(~agegroup,scales="free_y")}
      if (grepl("yday|age",unique(df_plot$varname))) {
        p <- p + geom_hline(yintercept=0,linetype="dashed",size=1/2) # + scale_y_continuous(breaks=break_vals) 
      } else { 
        p <- p + geom_hline(yintercept=1,linetype="dashed",size=1/2) }; p
      # save
      # folders
      if (all_plots) {
      subfldr_name <- paste0("median_interquant_by_param_value/from_dynamics/",sel_var,"/",start_year,"/")
      if ( !dir.exists(here(foldername,subfldr_name)) ) { dir.create(here(foldername,subfldr_name)) }
      # filename
      plot_filename <- gsub("//","/",here(foldername,paste0(subfldr_name,sel_par,
                                      ifelse(start_year==2020,"_incl2020",""),
                                      "_epiyear_startwk",epi_year_week_start,".png")))
      # save
      ggsave(plot_filename,width=28,height=16,units="cm") } else {
        ggsave(here(figs_folder,"SI_FIG8.png"),width=28,height=16,units="cm")  
      }
        
      # print progress of loop
      print(paste0(c(sel_varnames[k_plot_var],sel_pars[k_plot_par],
                     sel_var,epi_year_week_start,start_year),collapse=", "))
        }
      }
  }
} # epi_year_start week
