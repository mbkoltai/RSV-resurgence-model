# install.packages("devtools"); library("devtools")
# install_github("SineadMorris/shinySIR")

# install.packages("deSolve"); install.packages('odin')
library(dplyr); library(readr); library(reshape2); library(stringr); library(ggplot2); library(tidyr); 
library(deSolve)

######################################################################
# Function to compute the derivative of the ODE system
#
#  t - time
#  y - current state vector of the ODE at time t
#  parms - Parameter vector used by the ODE system
#
# Returns:
#  list with one component being a vector of length two containing
#  dS(t)/dt and dI(t)/dt
######################################################################

sir <- function(t, y, parms) {
  beta <- parms[1]; gamma <- parms[2]; S <- y[1]; I <- y[2]
  return(list(c(S = -beta*S*I, I = beta*S*I - gamma*I)))
}

######################################################################

# Population size 
N<-1e6
# Rate at which person stays in the infectious compartment (disease specific and tracing specific)
gamma<-0.5*1/5
# Infectious contact rate - beta = R0/N*gamma and when R0 \approx 2.25 then  2.25/N*gamma
beta<-0.5*4.5e-07
# R0 for the beta and gamma values
R0 <- beta*N/gamma

######################################################################
# Load package to numerically solve ODEs
# suppressPackageStartupMessages(library(deSolve))

# Grid where to evaluate
max_time<-300; times <- seq(0, max_time, by=0.01)

# Solve ODE system using Runge-Kutta numerical method.
initvals_S_I=c(N-10,10)
for (t_scale in seq(0.5,2.5,0.25)) {
  ode_solution <- rk4(y=initvals_S_I, times = times, func = sir, parms = t_scale*c(beta, gamma))
  df_ode_solution=ode_solution %>% as.data.frame() %>% setNames(c("t", "S", "I")) %>% 
    mutate(beta=beta,gamma=gamma,R0=N*beta/gamma,s=S/N,i=I/N,t_scale=t_scale, type="without_intervention")
  if (t_scale==0.5) {df_ode_solution_all=df_ode_solution} else {
    df_ode_solution_all=rbind(df_ode_solution_all,df_ode_solution)
    }
}

ggplot(df_ode_solution_all, aes(x=t)) + 
  geom_line(aes(y=s),color="darkred") + geom_line(aes(y=i), color="steelblue", linetype="twodash") +
  ggtitle('SIR model') + theme_bw() + 
  # facet_grid(rows=vars(period),cols=vars(region),switch="y") +
  facet_wrap(vars(t_scale),ncol=3) + 
  theme(axis.title=element_text(size=9),
        # axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=5),
        # axis.text.y=element_text(size=8), plot.title=element_text(size=11), 
        panel.grid=element_line(linetype="dashed",colour="black",size=0.1)) + 
  xlab('days') + ylab('S,I') + xlim(0,max(times)) + ylim(0,1) # + scale_y_log10()
