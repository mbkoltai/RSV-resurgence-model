# create the inset of Fig4
# comparing the param distribs of simuls that can replicate
# the early resurgence to others that were fitted only
# to pre-COVID-19 data

plot_partable_histogram_offseas <- partable_regular_dyn %>% filter(par_id %in% subsample_par) %>%
  mutate(`early off season`=par_id %in% early_off_season) %>% select(!const_delta) %>%
  mutate(`waning (days)`=1/omega,`maximal forcing (% above baseline)`=1e2*seasforce_peak) %>% 
  rename(`R0 (baseline)`=R0) %>%
  select(!c(omega,seasforce_peak)) %>% # `R0 peak`=R0*(1+seasforce_peak)
  rename(`age-dependence`=age_dep,`exposure-dependence`=exp_dep,
         `peak forcing (week)`=peak_week,`season width (weeks)`=seasforc_width_wks) %>%
  group_by(`early off season`) %>% pivot_longer(!c(par_id,`early off season`))

median_parvals = plot_partable_histogram_offseas %>% group_by(`early off season`,name) %>% 
  filter(!name %in% "peak forcing (week)") %>% summarise(median_parval=median(value))
KS_test = plot_partable_histogram_offseas %>% filter(!name %in% "peak forcing (week)") %>% 
  # select(!const_delta) %>% pivot_longer(!c(par_id,`early off season`)) %>%
  group_by(name) %>% summarise(min_val=min(value),max_val=max(value),
                               p_val=ks.test(x=value[`early off season`],y=value[!`early off season`])$p.value,
                               signif=p_val<0.01)
# plot
p <- ggplot(plot_partable_histogram_offseas %>% 
              filter(!name %in% "peak forcing (week)") %>%
              filter(!grepl("expos|waning",name)) %>%
              mutate(name=ifelse(grepl("maximal",name),"max. forcing (%)",name)),
            aes(x=value,color=`early off season`)) + stat_ecdf(geom="step") +
  facet_wrap(~name,scales="free_x",nrow=2) +
  scale_color_manual(values=c("grey","blue"),guide=guide_legend(override.aes=list(size=3))) + 
  geom_vline(data=median_parvals %>% filter(!grepl("expos|waning",name)) %>%
               mutate(name=ifelse(grepl("maximal",name),"max. forcing (%)",name)),
             aes(xintercept=median_parval,color=`early off season`),
             linetype="dashed",size=1/2,show.legend=F) +
  geom_text(data=KS_test %>% filter(!grepl("expos|waning",name)) %>%
              mutate(name=ifelse(grepl("maximal",name),"max. forcing (%)",name)),
            aes(x=max_val*0.95,y=0.32,label=paste0("p=",signif(p_val,3),ifelse(signif,"**","")) ),
            color="black",size=2.5,angle=90) +
  xlab("parameter values") + ylab("CDF") + theme_bw() + standard_theme + 
  theme(strip.text=element_text(size=15),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),legend.position="top",
        legend.text=element_text(size=14),legend.title=element_text(size=14))
# SAVE
# ggsave("simul_output/2e3_accepted_linear_relaxing/param_distrib_offseas_timing_ECDF.png",
#        width=28,height=18,units="cm")

p <- p + theme(legend.position="none") + 
  theme(strip.text=element_text(size=9),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=8),axis.text.y=element_text(size=8)); p
# ggsave("simul_output/2e3_accepted_linear_relaxing/param_distrib_offseas_timing_ECDF_largefont.png",
#        width=10,height=8,units="cm")
ggsave(here(figs_folder,"FIG4_inset.png"),width=10,height=8,units="cm")
