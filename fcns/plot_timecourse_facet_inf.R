ggplot(subset(df_ode_solution_tidy,grepl('I',name) & (t %% 7 ==0) & t/365>xval_lims[1] & t/365<xval_lims[2]),
       aes(x=t/365,y=get(value_type),group=name)) + geom_area(aes(fill=agegroup_name),color="black",position=position_stack(reverse=T)) + 
  facet_wrap(~infection,ncol=1,scales=scale_val) + theme_bw() + standard_theme + theme(axis.text.x=element_text(size=11,vjust=0.5),
                                                                                       axis.text.y=element_text(size=12),legend.position="top",legend.title=element_blank(),strip.text=element_text(size=12)) +
  scale_x_continuous(breaks=xval_breaks,expand=expansion(0,0)) + scale_y_continuous(expand=expansion(0.01,0)) +
  geom_rect(aes(xmin=shutdwn_lims[1]/365,xmax=shutdwn_lims[2]/365,ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
  geom_vline(data=seas_lims_plot %>% pivot_longer(cols=!season),aes(xintercept=value),color="blue",linetype="dashed",size=0.3) + 
  xlab('years') + ylab(y_axis_tag) + labs(subtitle=subtitle_str,caption=caption_txt)