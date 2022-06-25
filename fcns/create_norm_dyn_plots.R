all_plots=F

sel_years <- c("2019","2021","2022","2023")
for (k_par in unique(summ_dyn_all_parsets_broad_age_relat_time$parname)[6] ){
  for (k_age in 2) {
    # up to which age group?
    sel_agegr<-c("<1y","1-2y","2-5y")[1:ifelse(k_age==1,2,3)]
    # subset data to plot
    df_plot <- summ_dyn_all_parsets_broad_age_relat_time %>% 
      filter(agegroup %in% sel_agegr & epi_year %in% sel_years & parname %in% k_par) %>% 
      rename(`epi-year`=epi_year)
    # color palette
    n_par<-length(unique(df_plot$par_bin))
    colorpal <- colorRampPalette(colors=c("blue","grey","red"))(n_par)
    ggplot(df_plot,
           aes(x=peak_week_distance,group=par_median,color=factor(signif(par_median,3)))) + 
      geom_line(aes(y=median)) + # geom_line(aes(y=ci50_up),linetype="dashed") +# ,size=1.1
      # geom_ribbon(aes(ymin=ci50_low,ymax=ci50_up,
      #   fill=factor(signif(par_median,3))),color=NA,alpha=0.05) + 
      facet_grid(agegroup~`epi-year`,
                 labeller=labeller(`epi-year`=label_both),scales="free_y") +
      scale_x_continuous(limits=c(-20,15),expand=expansion(0.02,0)) + theme_bw() + standard_theme + 
      scale_color_manual(values=colorpal) + scale_fill_manual(values=colorpal) +
      manuscript_large_font_theme +
      xlab("distance in weeks from (pre-pandemic) peak week") +  
      ylab(paste0("weekly hospitalisations (1=pre-pandemic peak)")) + labs(color=k_par,fill=k_par) +
      geom_vline(xintercept=c(-10,10),linetype="dashed",size=1/5) + geom_vline(xintercept=0,size=1/3) +
      geom_hline(yintercept=1,linetype="dashed",size=1/3)
    
    # p 
    # save
    # folder
    if (all_plots){
    if (!dir.exists(here(foldername,"dynamics"))) { dir.create(here(foldername,"dynamics"))}
    # filename
    plot_fn <- gsub("//","/",here(foldername, paste0("dynamics/weekly_norm_hosp_under5y_by_",k_par,".png") ) )
    # plot_fn <- gsub("//","/",here(foldername, 
    #   paste0("dynamics/early_offseason/weekly_norm_hosp_under5y_by_",k_par,".png") ) )
    # comp_year,"_peak_until",max(sel_years), ifelse(k_age>1,"_1_5y",""),
    # save
    message(gsub(here(foldername),"",plot_fn))
    ggsave(plot_fn,width=30,height=18,units="cm"); print(gsub(foldername,"",plot_fn)) 
    } else {
      
      ggsave(here(figs_folder,"FIG3.png"),width=30,height=18,units="cm")
    }
    
  }
}
