
# peak values are not generated from this dataframe because of the agegroup resolution being different,
# instead they are calculated from the dynamics below

sel_vars <- unique(summ_broad_age_groups_byvalue$varname)
sel_pars <- unique(summ_broad_age_groups_byvalue$parname) 

for (k_plot_var in 1:2) { # 1:length(sel_vars)
  for (k_plot_par in 2) { # 1:length(sel_pars)
    sel_par <- sel_pars[k_plot_par]; dodge_val=1
    ylab_tag <- ifelse(grepl("mean",sel_vars[k_plot_var])," (from pre-NPI)"," (relative to pre-NPI)" )
    # subset for plotting
    df_plot <- summ_broad_age_groups_byvalue %>% 
      filter(epi_year %in% c("2020-21","2022","2023") & 
               (varname %in% sel_vars[k_plot_var]) & 
               (parname %in% sel_par) & (!agegroup_broad %in% "5+y") ) %>%
      mutate(varname=case_when(
        grepl("hosp_tot_norm",varname) ~ "cumulative hospitalisations",
        grepl("peak_hosp_norm",varname) ~ "peak hospitalisations",
        grepl("mean",varname) ~ "change (months) in mean age <5y hospitalisations"),
        parname=case_when(grepl("age_dep",parname) ~ "age-dependence", 
                          grepl("exp_dep",parname) ~ "exposure-dependence",
                          grepl("kappa",parname) ~ "exposure (-1) <-> age (1)", 
                          grepl("age_exp_ratio",parname) ~ "ratio of age- to \nexposure-dependence", 
                          grepl("seasforc_width_wks",parname) ~ "season width (weeks)",
                          grepl("seasforce_peak",parname) ~ "seasonal forcing (above baseline)",
                          grepl("R0",parname) ~ "R0 (baseline)",
                          grepl("omega",parname) ~ "waning (days)")) %>%
      mutate(parvalue=ifelse(grepl("waning",parname),1/parvalue,parvalue),
             epi_year=gsub("2020-21","2021",epi_year))
      
    
    n_par_value <- length(unique(df_plot$par_bin)); median_width=0.5/100
    if (grepl("mean",sel_vars[k_plot_var])) {median_width=2/100}
    # colour palette
    if (!any(c("age_dep","exp_dep","kappa","age_exp_ratio") %in% sel_par)){ 
      colorpal=colorRampPalette(colors=c("orange","red"))(n_par_value)} else  {
        colorpal=colorRampPalette(colors=c("blue","grey","red"))(n_par_value)  }
    scale_size_val=1.32
    p <- ggplot(df_plot,aes(x=factor(epi_year),color=factor(signif(parvalue,2)),group=parvalue)) + 
      facet_wrap(~agegroup_broad,scales="free_y") + 
      geom_linerange(aes(ymin=ci95_l,ymax=ci95_u),position=position_dodge(width=dodge_val),
                     alpha=0.3,size=scale_size_val*ifelse(any(grepl("mean",df_plot$varname)),72,24)/n_par_value) +
      geom_linerange(aes(ymin=ci50_l,ymax=ci50_u),position=position_dodge(width=dodge_val),
                     alpha=0.6,size=scale_size_val*ifelse(any(grepl("mean",df_plot$varname)),72,24)/n_par_value) +
      geom_linerange(aes(x=factor(epi_year),ymin=median-median_width,ymax=median+median_width),
                     position=position_dodge(width=dodge_val),color="black",
                     size=scale_size_val*ifelse(any(grepl("mean",df_plot$varname)),72,24)/n_par_value) + 
      geom_vline(xintercept=(0:4)+1/2,size=1/5) + labs(color=unique(df_plot$parname)) +
      scale_x_discrete(expand=expansion(0.02,0)) +
      xlab("") + ylab(paste0(gsub("mean age","mean age\n",unique(df_plot$varname)),ylab_tag)) + 
      scale_color_manual(values=colorpal) + geom_hline(yintercept=1,linetype="dashed",size=1/2) +
      theme_bw() + standard_theme + 
      theme(strip.text=element_text(size=15),axis.text.x=element_text(size=13),
            axis.title.y=element_text(size=11),axis.text.y=element_text(size=13),
            legend.text=element_text(size=11),legend.title=element_text(size=12),
            legend.position="bottom") + guides(fill=guide_legend(ncol=3)) + 
      manuscript_large_font_theme
    # ifelse(!grepl("exp|age|kappa",unique(df_plot$parname)),"bottom","right")
    # 
    if (grepl("mean_age",sel_vars[k_plot_var])) {p <- p + scale_y_continuous(breaks=-5:5); p}
    
    # create filename and folders
    all_plots=F
    if (all_plots){
    sel_var_filename <- gsub("-|\\s","_",sel_vars[k_plot_var])
    if (!dir.exists(here(foldername,"median_interquant_by_param_value"))) {
      dir.create(here(foldername,"median_interquant_by_param_value"))}
    subfldr_name <- ifelse(grepl("peak",sel_vars[k_plot_var]),
                           "median_interquant_by_param_value/from_summ/peak/",
                           "median_interquant_by_param_value/from_summ/cumul/")
    if (grepl("mean",sel_vars[k_plot_var])) { 
      subfldr_name <- "median_interquant_by_param_value/from_summ/mean_age/" }
    if (!dir.exists(here(foldername,subfldr_name)) ) { dir.create(here(foldername,subfldr_name)) }
    # save
    ggsave(here(foldername,paste0(subfldr_name,"summ_stats_relative_2019_",
                                  paste0(sel_var_filename,collapse="_"),"_",sel_par,".png")),
           width=28,height=16,units="cm")
    } else {
      if (k_plot_var==1){
      ggsave(here(figs_folder,"FIG2B_summ.png"),width=28,height=16,units="cm") 
        } else {
      ggsave(here(figs_folder,"SI_FIG9.png"),width=28,height=16,units="cm")
      }
           }
    
    
    print(paste0(c(sel_vars[k_plot_var],sel_pars[k_plot_par]),collapse=", "))
  }
}