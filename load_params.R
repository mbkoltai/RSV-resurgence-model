# functions
rm(list=ls()); # currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
# library(contactdata); library(fitdistrplus);  library(bbmle); library(Rcpp); library(GillespieSSA)
x1<-c("tidyverse","deSolve","gtools","rstudioapi","devtools","wpp2019","Rcpp","lubridate",
      "tsibble","pracma","qs","zoo","RcppRoll") # ,"ungeviz"
x2 <- x1 %in% row.names(installed.packages()); if (any(x2 == FALSE)) { install.packages(x1[!x2]) }
# if (!any(grepl("ungeviz",row.names(installed.packages())))) {devtools::install_github("wilkelab/ungeviz")}
# Load all packages (unless already loaded) # as.Date <- zoo::as.Date
lapply(x1,library,character.only=TRUE)
select <- dplyr::select; # row_number <- dplyr::row_number; summarise <- dplyr::summarise
# source('fcns/RSV_model_functions.R')
source("fcns/essential_fcns.R")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# set plotting theme
standard_theme=theme(plot.title=element_text(hjust=0.5,size=16),axis.text.x=element_text(size=9,angle=90,vjust=1/2),
  axis.text.y=element_text(size=9),axis.title=element_text(size=14), panel.grid.minor=element_blank())
# text=element_text(family="Calibri") # some versions will throw errors if font set to Calibri
# panel.grid=element_line(linetype="solid",colour="black",size=0.1),
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# SET PARAMETERS --------------------------------------------------------
### ### ### ### ### ### ### ###
# constant parameters
# selected country
country_sel="United Kingdom"
# time resolution (in days)
elem_time_step=1
# population data
standard_age_groups <- fun_cntr_agestr(country_sel,i_year="2020",seq(0,75,5),c(seq(4,74,5),99))
popul_struct=fcn_cntr_fullpop(n_year="2020",country_sel)
# RSV age groups (population data from wpp2019)
rsv_age_groups<-fun_rsv_agegroups(standard_age_groups,popul_struct,rsv_age_groups_low=c(0,0.5,1,1.5, 2,3,4, 5,15, 45, 65),
                                  rsv_age_group_sizes=c(rep(0.4,4),rep(0.9,3), 9, 29, 19, 34))
# rsv_age_groups$value=rsv_age_groups$value*67e6/sum(rsv_age_groups$value)
ons_2020_midyear_estimates_uk <- read_csv(here::here("repo_data/ons_2020_midyear_estimates_uk.csv")) %>% 
  mutate(age_num=as.numeric(gsub("\\+","",age)))
low_inds<-findInterval(rsv_age_groups$age_low,ons_2020_midyear_estimates_uk$age_num)
high_inds <- findInterval(rsv_age_groups$age_low+rsv_age_groups$duration-0.1,ons_2020_midyear_estimates_uk$age_num)
rsv_age_groups$value <- unlist(lapply(1:length(low_inds), 
              function(x) sum(ons_2020_midyear_estimates_uk$value[low_inds[x]:high_inds[x]])*ifelse(rsv_age_groups$duration[x]<1,
                                            rsv_age_groups$duration[x],1) ))
# DEATHS (2019: 530841 deaths [England and Wales!]) # "uk_death_rate_byage_rsv_agegroups.csv" is for 1000 population!
# read_csv("data/uk_death_rate_byage_rsv_agegroups.csv")
# i slightly adjusted the age-specific deaths rates to get a stationary population at the 2019 total and age struct
uk_death_rate=c(rep(1e-5,2)*3,rep(1e-6,5),rep(0,2),1e-6,1.79e-4) 
# number of age groups, reinfections and variables (S,I,R)
n_age=nrow(rsv_age_groups); varname_list=c('S','I','R'); n_compartment=length(varname_list); n_inf=3
dim_sys=n_age*n_compartment*n_inf; n_days_year=365
# BIRTH RATE into S_1_1 (Germany 2019: 778e3 births)
daily_births=2314; birth_rates=matrix(c(daily_births,rep(0,dim_sys-1)),dim_sys,1)
# we want population to be stationary (at 2019 or 2020 value), so deaths = births
if (!any(grepl("death",colnames(rsv_age_groups)))){
  rsv_age_groups <- rsv_age_groups %>% mutate(deaths_per_person_per_day=uk_death_rate,
          stationary_popul=fcn_calc_stat_popul(rsv_age_groups,rsv_age_groups$duration,daily_births,uk_death_rate,output_type="")) 
}
# query variables: fun_sub2ind(1:3,11,"R",varname_list,n_age,n_inf)
# force of infection terms
# linear indices of the I & S variables
l_inf_susc=fun_inf_susc_index_lists(n_age,n_inf,varname_list); inf_vars_inds=l_inf_susc[[1]]; susc_vars_inds=l_inf_susc[[2]]
# CONTACT MATRIX
# contact matrix from covidm ("home","work","school","other")
# cm_path="~/Desktop/research/models/epid_models/covid_model/lmic_model/covidm/"
# if UK -> England's contact matrix # check: cm_parameters_SEI3R(cm_uk_locations("UK", 1))$pop[[1]]$matrices 
# list_contmatrs=fun_covidm_contactmatrix(country_sel,currentdir_path,cm_path=cm_path) 
# make matrix reciprocal
# C_m_polymod=Reduce('+',list_contmatrs) # fun_recipr_contmatr(Reduce('+',list_contmatrs),
#   age_group_sizes=standard_age_groups$values)
C_m_polymod<-readRDS(here::here("repo_data/UK_contact_matrix_sum.RDS"))

# create for our age groups
C_m_merged_nonrecipr=fun_create_red_C_m(C_m_polymod,rsv_age_groups,
        orig_age_groups_duration=standard_age_groups$duration,orig_age_groups_sizes=standard_age_groups$values)
# make it reciprocal for the larger group
C_m=fun_recipr_contmatr(C_m_merged_nonrecipr,age_group_sizes=rsv_age_groups$stationary_popul)
# bc of reinfections we need to input contact matrix repeatedly: 
# the force of infection is the same for the 1st, 2nd, 3rd infections, but in the full equation they are multiplied by a different
# susceptibility parameter (delta1>delta2>delta3). 
# Therefore we repeat the contact matrix's rows three times, corresponding to 1st, 2nd, 3rd infxs of each age group
# Also the columns of the matrix are normalised by the age group sizes, as they will multiply the infection terms by age groups
contmatr_rowvector=t(do.call(cbind, lapply(1:nrow(C_m), function(x){diag(C_m[x,]) %*% matrix(1,n_age,n_inf)}))) / 
                    matrix(rep(rsv_age_groups$stationary_popul,n_age*n_inf),nrow=n_age*n_inf,byrow=T)
# rsv_age_groups$stationary_popul[col(C_m)]
# matrix(1,nrow=33,ncol = 11)/matrix(rep(1:11,33),nrow=33,byrow=T)
every_nth = function(n) { return(function(x) {x[c(TRUE, rep(FALSE, n-1))]}) }
# build kinetic matrix
# WANING (immunity) terms: R_i_j -> S_min(i+1,n_inf)_j
omega=1/350 # 1/runif(1,60,200)
# RECOVERY
rho=1/7 # 1/rho=rweibull(1, shape=4.1,scale=8.3)
# KINETIC MATRIX (aging terms need to be scaled by duration of age groups!)
K_m=fun_K_m_sirs_multiage(dim_sys,n_age,n_inf,n_compartment,rho,omega,varname_list,
                          agegroup_durations=rsv_age_groups$duration)
# agegroup indices for maternal immunity
mat_imm_flag <- TRUE; mat_imm_inds<-list(fun_sub2ind(i_inf=1,j_age=1,"R",c("S","I","R"),n_age,3),
                                         fun_sub2ind(i_inf=c(1,2,3),j_age=9,"R",c("S","I","R"),n_age,3),
                                         fun_sub2ind(i_inf=c(1,2,3),j_age=9,"S",c("S","I","R"),n_age,3))
print("end of load_params.R")