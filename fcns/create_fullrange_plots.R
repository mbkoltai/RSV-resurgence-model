# PLOT
# load tidybayes library
# if (!any(grepl("tidybayes",row.names(installed.packages())))) {
#   install.packages("tidybayes")}; library(tidybayes)

df_plot_fullrange <- output_ranges_full_scan %>% filter(!agegroup_broad %in% "5+y") %>%
  mutate(
    scan_param=case_when(
      grepl("age_dep",scan_param) ~ "age-dependence",
      grepl("exp_dep",scan_param) ~ "exposure-dependence",
      grepl("seasforce_peak",scan_param) ~ "seasonal forcing (strength)",
      grepl("seasforc_width_wks",scan_param) ~ "seasonal forcing (width)",
      grepl("R0",scan_param) ~ "baseline R0",
      grepl("omega",scan_param) ~ "waning rate"),
    epi_year=case_when(epi_year %in% "2021" ~ "2021-2022")) %>% 
  filter(epi_year %in% "2021-2022")

# plot change in CUMUL and PEAK HOSP
dodge_val=1
ggplot() +
  # geom_boxplot(data=df_plot_fullrange %>% filter(range %in% "full"),
  #              aes(y=agegroup_broad,xmiddle=norm_median,
  #                      xlower=norm_ci95_low,xupper=norm_ci95_up,xmin=NA,xmax=NA,
  #                  group=interaction(scan_param,agegroup_broad)),
  #                  color="grey",# color=factor(scan_param)),
  #                 position=position_dodge(width=dodge_val),stat="identity",width=1/4) +
  geom_boxplot(data=df_plot_fullrange %>% filter(range %in% "sel"),
               aes(y=agegroup_broad,xmiddle=norm_median,
                   xlower=norm_ci95_low,xupper=norm_ci95_up,xmin=norm_min,xmax=norm_max,
                   group=interaction(scan_param,agegroup_broad),
                   color=factor(scan_param)),fill=NA,
               position=position_dodge(width=dodge_val),stat="identity",width=3/4) +
  facet_grid(~vartype,scales="free_x") + 
  geom_hline(yintercept=(1:2)+1/2,size=1/2) + 
  geom_vline(xintercept=1,size=1/2,linetype="dashed") + 
  guides(color=guide_legend(nrow=2,byrow=TRUE)) +
  scale_x_log10(breaks=c(0.3,0.5,0.75,1,1.5,2,3,5,10)) + 
  scale_y_discrete(expand=expansion(0.1,1/3)) + 
  xlab("relative hospitalisation risk compared to pre-pandemic years") + 
  ylab("") + labs(color="") +
  theme_bw() + standard_theme + manuscript_large_font_theme

# save
subfldr_name <- here(figs_folder,"full_range_plots/")
# median_interquant_by_param_value/from_summ/summary_range/
if (!dir.exists(subfldr_name)) {dir.create(subfldr_name)}
# ggsave(paste0(subfldr_name,"cumul_peak_hosp_summary_plot_with_reject_parsets.png"),
#   width=28,height=24,units="cm")
ggsave(paste0(subfldr_name,"cumul_peak_hosp_summary_plot.png"),width=28,height=24,units="cm")

# check entire range across all params (percentage change from pre-pandemic)
df_plot_fullrange %>% 
  group_by(agegroup_broad,epi_year,range,vartype) %>% 
  summarise(norm_median=round(100*median(norm_median)-100),
            norm_min=round(100*min(norm_min)-100),
            norm_max=round(100*max(norm_max)-100)) %>% 
  filter(range %in% "sel" & vartype %in% "cumulative")

# shift in mean age
df_plot_mean_age <- mean_age_shift_ranges %>%
  mutate(scan_param=case_when(grepl("age_dep",scan_param) ~ "age-dependence",
                              grepl("exp_dep",scan_param) ~ "exposure-dependence",
                              grepl("seasforce_peak",scan_param) ~ "seasonal forcing (strength)",
                              grepl("seasforc_width_wks",scan_param) ~ "seasonal forcing (width)",
                              grepl("R0",scan_param) ~ "baseline R0",
                              grepl("omega",scan_param) ~ "waning rate"),
         epi_year=case_when(epi_year %in% "2021" ~ "2021-2022")) %>% 
  filter(epi_year %in% "2021-2022")

# plot shift in mean age
p_mean_age_shift <- ggplot() +
  geom_boxplot(data=df_plot_mean_age %>% filter(range %in% "sel"),
               aes(y=1,xmiddle=norm_median,xlower=norm_ci95_low,
                   xupper=norm_ci95_up,xmin=norm_min,xmax=norm_max,
                   group=scan_param,color=factor(scan_param)),
               position=position_dodge(width=dodge_val),stat="identity",width=3/4,size=1) +
  guides(color=guide_legend(nrow=2,byrow=TRUE)) + 
  scale_x_continuous(breaks=(-1:7)) + scale_y_discrete(expand=expansion(0,0.05)) + 
  xlab("shift in average age (months)") + ylab("") + labs(color="") +
  theme_bw() + standard_theme + manuscript_large_font_theme
# also show rejected params?
all_pars_show=F
if (all_pars_show) {
  p_mean_age_shift + geom_boxplot(data=df_plot_mean_age %>% filter(range %in% "full") %>% 
                                    mutate(norm_min=NaN,norm_max=NaN),
                                  aes(y=1,xmiddle=norm_median,xlower=norm_ci95_low,xupper=norm_ci95_up,
                                      xmin=norm_min,xmax=norm_max,group=scan_param),
                                  color="grey",position=position_dodge(width=dodge_val),
                                  stat="identity",width=3/4,fill=NA)
  # save
  ggsave(here(subfldr_name,"aver_age_hosp_summary_plot_all_par.png"),width=28,height=6,units="cm")
} else { 
  p_mean_age_shift
  ggsave(here(subfldr_name,"aver_age_hosp_summary_plot_accepted_par.png"),width=28,height=12,units="cm") 
}
