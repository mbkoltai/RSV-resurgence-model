library(profvis)

g(n_years,timesteps,simul_start_end,forcing_vector_npi) %=% fun_shutdown_seasforc(npi_dates,years_pre_post_npi=c(0,1),
  season_width_wks=seasforc_width_wks,init_mt_day="06-01",peak_week=44,forcing_above_baseline=seaspeakval,npireduc_strength=0.5)
params<-list(list(birth_rates,death_rates),K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds,delta_susc)
# interpolation fcns for seas forcing & extern introds
approx_seas_forc <- approxfun(data.frame(t=timesteps,seas_force=forcing_vector_npi))
# how many introductions every 30 days?
approx_introd <- approxfun(data.frame(t=timesteps,as.numeric(timesteps %% 30==0)*5))

profvis({ ode_sol<-lsoda(initvals_sirs_model,timesteps,func=sirs_seasonal_forc,parms=params) })

###
# stationary demographics
uk_death_rate
list(sum(birth_rates),uk_death_rate$deaths_per_1000person_peryear/1e3)

ode_demogr <- function(t,X,parms){
  birthrates=parms[[1]]; deathrates=parms[[2]]; aging_terms=parms[[3]]; dimsys=nrow(birthrates)
  # print(dim(X))
  dXdt=birthrates- diag(deathrates) %*% t(t(X)) - diag(c(aging_terms[1:dimsys-1],0)) %*% t(t(X)) + 
    c(0,diag(aging_terms[1:dimsys-1]) %*% t(t(X[1:dimsys-1]))); list(dXdt) }

daily_births=round(sum(birth_rates))
deaths_balance_births=daily_births/(sum(birth_rates)*365*rsv_age_groups$duration[nrow(rsv_age_groups)])
# demogr params
demogr_par=list(birthrates=matrix(c(daily_births,rep(0,nrow(rsv_age_groups)-1))),
                deathrates=c(rep(0,nrow(rsv_age_groups)-1),deaths_balance_births), # 
                aging_terms=1/(rsv_age_groups$duration*365))
init_cond=(ode_demogr_sol %>% filter(t==max(t)))$value # sum(birth_rates)*365*rsv_age_groups$duration 
ode_demogr_sol <- lsoda(y=init_cond,times=1:1e4,func=ode_demogr,parms=demogr_par) %>% as.data.frame() %>%
   setNames(c("t",paste("X",1:11,sep=""))) %>% pivot_longer(!t) %>% mutate(name=factor(name,levels=unique(name)))
ggplot(ode_demogr_sol,aes(x=t,y=value,color=name)) + geom_line() + facet_wrap(~name,scales="free_y")+ theme_bw() + standard_theme

# steady state sol compared to rectangular popul struct
round((ode_demogr_sol %>% filter(t==max(t)))$value/(sum(birth_rates)*365*rsv_age_groups$duration),3)
# kinet_matrix = birthrates + 

# solve steady state problem in matrix form
daily_births=2314
# all rates PER DAY
D_m = -diag(c(rep(1e-5,2)*3,rep(1e-6,5),rep(0,2),1e-6,1.79e-4)) # # uk_death_rate$deaths_per_1000person_peryear/(365*1e3)
B_m = matrix(0,nrow = nrow(rsv_age_groups),ncol=ncol(rsv_age_groups))
diag(B_m) = -c(1/(rsv_age_groups$duration[1:(nrow(rsv_age_groups)-1)]*365),0)
# subdiag terms
subdiag_inds=(nrow(rsv_age_groups)+1)*(0:(nrow(rsv_age_groups)-2))+2
B_m[subdiag_inds]=1/(rsv_age_groups$duration[1:(nrow(rsv_age_groups)-1)]*365)
K_matr=D_m+B_m
# has almost zero determinant
x_ss=-round(inv(K_matr) %*% matrix(c(daily_births,rep(0,nrow(rsv_age_groups)-1))) )
# sum popul
sum(x_ss)
# compare to current UK population
round(c(x_ss/rsv_age_groups$value),2)
# popul struct
rbind(1e2*round(c(rsv_age_groups$value/sum(rsv_age_groups$value)),3),
  rbind(1e2*round(c(x_ss/sum(x_ss)),3),  1e2*(round(c(rsv_age_groups$value/sum(rsv_age_groups$value)),3) - round(c(x_ss/sum(x_ss)),3))))

# check if correct (this should be vectors of 0)
matrix(c(daily_births,rep(0,nrow(rsv_age_groups)-1))) + round(K_matr %*% x_ss)
# deaths = births
round(t(matrix(uk_death_rate$deaths_per_1000person_peryear/(365*1e3))) %*% t(t(x_ss)) - daily_births)

# rectangular
rect_popul=daily_births*365*rsv_age_groups$duration/(sum(daily_births*365*rsv_age_groups$duration)/sum(rsv_age_groups$value))