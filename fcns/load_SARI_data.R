# load RSV hospitalisation data, estimate underreporting

SARIwatch_RSVhosp_under5_2018_2020_weekly_counts <- 
  read_csv("data/SARIwatch_RSVhosp_under5_2018_2020_weekly_counts.csv",col_types="ffddd") %>% 
  mutate(wk_n=gsub("\\.","-W",wk_n),wk_n=factor(wk_n,levels=unique(wk_n)),
         date=ISOweek2date(paste0(gsub("\\.","-W",wk_n),"-1")))
# eventually save this to `repo_data`
# compare SARI-Watch with literature estimates 
# annual RATE of hospitalisations per 100K population: 500/1e5 (18-19), 494/1e5
SARIwatch_RSVhosp_under5_2018_2020_weekly_counts %>% group_by(year) %>% 
  summarise(annual_cumul_rate=sum(rate_under5yrs))

# under-reporting for under-5 hospitalisations
# annual RATE of hospitalisations from Reeves 2017: <1y: 35.1/1000, 1-4y: 5.31/1000 = 1083.8/1e5
# annual RATE from Taylor 2016: <0.5y: 4184/1e5, 6-23mts: 1272/1e5, 2-4y: 114/1e5 = 822.7/1e5
reported_hosp_rate_per_100k <- c(
              reeves_2017=(35.1*sum(rsv_age_groups$value[1:2]) + 
                                  5.31*sum(rsv_age_groups$value[3:7]))*100/sum(rsv_age_groups$value[1:7]),
              taylor_2016=(4184*sum(rsv_age_groups$value[1]) + 1272*sum(rsv_age_groups$value[2:4]) + 
                                  114*sum(rsv_age_groups$value[5:7]))/sum(rsv_age_groups$value[1:7]) )
under_report_factor_under5 = mean( (SARIwatch_RSVhosp_under5_2018_2020_weekly_counts %>% group_by(year) %>% 
              summarise(annual_cumul_rate=
                          sum(rate_under5yrs)))$annual_cumul_rate)/mean(reported_hosp_rate_per_100k)

# hospitalisation counts for 65+
# for 65+y
SARIwatch_RSVhosp_over65_2018_2020_weekly_counts <- 
  read_csv("data/SARIwatch_RSVhosp_over65y_2018_2020_weekly_counts.csv") %>%
  mutate(wk_n=gsub("-","-W",wk_n),wk_n=factor(wk_n,levels=unique(wk_n)),
         date=ISOweek2date(paste0(wk_n,"-1")))

# rate per 100K from SARI-Watch
SARIwatch_RSVhosp_over65_2018_2020_weekly_counts %>% group_by(year) %>% summarise(annual_rate=sum(rate_65yplus))
# Estimates from literature 
# (Fleming 2015, https://bmcinfectdis.biomedcentral.com/articles/10.1186/s12879-015-1218-z/tables/3)
# 65-74: 86/100e3 (62-101); 75+: 234/100e3 (180-291)
over65_hosp_rate_100k_lit_estim = (
  sum(ons_2020_midyear_estimates_uk$value[ons_2020_midyear_estimates_uk$age %in% 65:74])*86 + 
    sum(ons_2020_midyear_estimates_uk$value[76:nrow(ons_2020_midyear_estimates_uk)])*234)/
  sum(ons_2020_midyear_estimates_uk$value[66:nrow(ons_2020_midyear_estimates_uk)])
# under-reporting rate
under_report_factor_over65y <- mean((SARIwatch_RSVhosp_over65_2018_2020_weekly_counts %>% 
                           group_by(year) %>% summarise(annual_rate=sum(rate_65yplus)))$annual_rate) / 
  over65_hosp_rate_100k_lit_estim