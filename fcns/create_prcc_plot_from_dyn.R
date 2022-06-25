# plot PRCC plots of cumulative/peak hospit, timing of peak relative to pre-NPI and mean age

for (k_epiyear_start_wk in c(23,40)[1]) {
  for (k_epiyear_plot in 2){ # 1:3
    for (k_output in 1){ # 1
    
    sel_years=list(2020,2021,2020:2021)[[k_epiyear_plot]]
    k_output_name = list("peak_hosp|cumul_hosp","peak_yday","mean_age_under5y")[[k_output]]
      
    df_plot <- bind_rows(list_prcc_dyn) %>% 
      filter(grepl(k_output_name,output) & 
               epi_year_wk_start==k_epiyear_start_wk & 
               epi_year %in% sel_years & 
               !(grepl("ratio",parname)) ) %>% 
      mutate(parname=factor(parname,levels=unique(parname)),
             agegroup=unique(results_fullscan_hosp$agegroup_broad)[agegroup],
             output=case_when(output %in% "peak_hosp" ~ "peak hospitalisation",
                              output %in% "cumul_hosp" ~ "cumulative hospitalisation",
                              output %in% "peak_yday" ~ "peak timing (day of year)",
                              output %in% "mean_age_under5y" ~ "mean age of hospitalisation (<5y)"),
             est=ifelse(p_below_0.05,est,NA))
      
    # create plot
    p <- ggplot(df_plot) + facet_wrap(~output,scales="free_x") +
      labs(fill="") + coord_flip() + xlab("") + ylab("PRCC") + 
      theme_bw() + standard_theme + manuscript_large_font_theme + 
      theme(panel.grid.major.y=element_blank())
    
    if (length(sel_years)>1) {
      p <- p + geom_col(aes(y=est,x=parname,
                            group=interaction(epi_year,agegroup),color=factor(epi_year),
                            fill=factor(agegroup)),size=2/3,
                            position=position_dodge(width=0.85),width=4/5) +
        scale_color_manual(values=c("black","red")) + labs(color="epi-year") + 
        geom_vline(xintercept=1/2+(1:4),size=1/3) + geom_hline(yintercept=0) + 
        guides(color=guide_legend(override.aes=list(fill=NA)))
    } else {
      p <- p + geom_col(aes(y=est,x=parname,group=agegroup,fill=factor(agegroup)),
                        size=1/3,color="black",position=position_dodge(width=0.85),width=4/5) +
        geom_vline(xintercept=1/2+(1:4),size=1/3) + geom_hline(yintercept=0)
    }
    if (length(unique(df_plot$agegroup))==1) { p <- p + scale_fill_manual(values="lightgrey")
    if (length(unique(df_plot$epi_year))==1) { p <- p + theme(legend.position=NULL)    }     }
    p
    
    # create folder/filenames
    subfldr_name="PRCC_from_dynamics/"; 
    if (!dir.exists(here(figs_folder,subfldr_name))) {
      dir.create(here(figs_folder,subfldr_name))
    }
    filename=paste0("PRCC_",paste0(gsub("peak_hosp","peak",gsub("\\|","",k_output_name)),collapse="_"),
                    paste0(sel_years,collapse="_"),"_startwk",k_epiyear_start_wk,".png")
    # SAVE
    ggsave(here(figs_folder,"FIG2A.png"),width=28,height=20,units="cm")
    message(filename)
    }
  }
}
