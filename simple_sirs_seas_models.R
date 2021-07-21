# simple SIRS model that shows oscillations
require(deSolve)
seirmod = function(t, y, parms) {
  S = y[1];   E = y[2];   I = y[3]; R = y[4]
  mu=parms["mu"]; N=parms["N"]; beta = parms["beta"]; sigma = parms["sigma"]; gamma = parms["gamma"]
  inf_term=beta*S*I/N
  dS = mu * (N - S) - inf_term; dE = beta * S * I/N - (mu + sigma) * E
  dI = sigma * E - (mu + gamma) * I; dR = gamma * I - mu * R # res = c(dS, dE, dI, dR)
  list(c(dS, dE, dI, dR)) }
# params
times = seq(0, 100, by=1/120)
paras = c(mu = 1/50, N = 1, beta = 1000,sigma = 365/8, gamma = 365/5)
start = c(S=0.06, E=0, I=0.001, R = 0.939)
# calc R0
R0 = expression(sigma/(sigma + mu) * beta/(gamma + mu)); with(as.list(paras), eval(R0))
# solve ODE
out = as.data.frame(ode(start, times, seirmod, paras))
par(mfrow = c(1,2))
ggplot(out,aes(x=S,y=I,color=time)) + geom_path() + theme_bw() + standard_theme
ggplot(out,aes(x=time,y=I)) + geom_line() + theme_bw() + standard_theme
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### chain SIR
chainSIR=function(t, logx, params){
  x=exp(logx)
  u=params["u"]
  S=x[1]; I=x[2:(u+1)]; R=x[u+2]
  with(as.list(params),{
    dS = mu * (N - S) - sum(beta * S * I) / N; dI = rep(0, u) 
    dI[1] = sum(beta * S * I) / N - (mu + u*gamma) * I[1]
    if(u>1){ for(i in 2:u){ dI[i]= u*gamma * I[i-1] - (mu+u*gamma)* I[i] } }
    dR = u*gamma * I[u] - mu * R #    res=c(dS/S, dI/I, dR/R)
    list(c(dS/S, dI/I, dR/R)) }) }
# chain SIR solve
times = seq(0, 20, by=1/52)
paras2 = c(mu = 1/75, N = 1, beta = 625,gamma = 365/14, u=1)
xstart2 = log(c(S=.06, I=c(0.001, rep(0.0001,paras2["u"]-1)), R = 0.0001))
out = as.data.frame(ode(xstart2, times, chainSIR,paras2))
plot(times, exp(out[,3]), ylab="Infected", xlab="Time", ylim=c(0, 0.01), type="l")
# u=2
paras2["u"]=2
xstart2 = log(c(S=.06, I=c(0.001, rep(0.0001/paras2["u"], paras2["u"]-1)), R = 0.0001))
out2 = as.data.frame(ode(xstart2, times, chainSIR,paras2))
lines(times, apply(exp(out2[,-c(1:2,length(out2))]),1 ,sum), col="blue")
# u=73
paras2["u"] =73
xstart2 = log(c(S=.06, I=c(0.001, rep(0.0001/paras2["u"], paras2["u"]-1)), R = 0.0001))
out3 = as.data.frame(ode(xstart2, times, chainSIR,paras2))
lines(times, apply(exp(out3[,-c(1:2,length(out3))]),1, sum), col="red", lwd=2, lty=2)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# algebr rootfinding
require(nleqslv)
rootfn=function(x, params){
  r=with(as.list(params), c(mu * (N - x[1]) - beta * x[1]* x[2] / N,
           beta * x[1] * x[2] / N - (mu + gamma) * x[2],
           gamma *x[2] - mu*x[3])); r }
parms = c(mu = 1/(50*52), N = 1, beta = 2.5,gamma = 1/2)
ans = nleqslv(x=c(0.1,0.5, 0.4), fn = rootfn,params=parms)
ans$x
# dis free equil?
ans=grid=expand.grid(seq(0,1, by=.25), seq(0,1, by=.25),seq(0,1, by=.25)); ans[,]=NA
for(i in 1:nrow(ans)){ ans[i,]=nleqslv(as.numeric(grid[i,]), fn=rootfn,params=parms)$x }
ans2=round(ans, 4)[!duplicated(ans2),]; ans2
# rootsolve package
sirmod=function(t, y, parameters){
  S=y[1];   I=y[2];   R=y[3]
  beta=parameters["beta"];   mu=parameters["mu"];   gamma=parameters["gamma"]; N=parameters["N"]
  dS = mu * (N - S) - beta * S * I / N
  dI = beta * S * I / N - (mu + gamma) * I
  dR = gamma*I - mu*R #   res=c(dS, dI, dR)
  list(c(dS, dI, dR)) }
# solve
require(rootSolve)
parms = c(mu=1/(50*52), N = 1,beta = 2.5, gamma = 1/2)
equil = runsteady(y = c(S=1-1E-4, I = 1E-4, R = 0),times = c(0, 1E05), func = sirmod, parms = parms)
round(equil$y, 4)
###
parms = c(mu = 1/(50*52), N = 1, beta = 2.5,gamma = 1/2)
N=parms["N"]; gamma=parms["gamma"]; beta=parms["beta"]; mu=parms["mu"]
Istar=as.numeric(mu*(beta/(gamma+mu)-1)/beta)
Sstar=as.numeric((gamma+mu)/beta)
# jacobian
dS = expression(mu * (1 - S) - beta * S * I/1)
dI = expression(beta * S * I/1 - (mu + gamma) * I)
vals = list(mu = 1/(50*52), N = 1, beta = 2.5,gamma = 1/2, S=Sstar, I=Istar)
J=with(vals, matrix(c(eval(D(dS, "S")), eval(D(dS, "I")), eval(D(dI,"S")),eval(D(dI,"I"))),ncol=2,byrow=T))
eigen(J, only.values=TRUE)$values
# resonant fewquency
2*pi/Im(eigen(J)$values[1])
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# SIRS model
N=1; gamma=7/3.8; omega = 1/(52 * 4); mu = 1/(52 * 70); R0 = 2.9
# R0=beta/(gamma+mu)
beta=R0*(gamma+mu); paras=c(beta=beta, gamma=gamma, mu=mu, omega=omega)
A=(omega+mu+gamma)/((omega+mu)*(beta-gamma-mu))
GI=1/(gamma+mu); GR = 1/(omega+mu)
Tper=4*pi/sqrt(4*(R0-1)/(GI*GR)-((1/GR)-(1/A))^2); Tper
# steady state
Sstar=1/R0; Istar=mu*(1-1/R0)/(gamma+mu-(omega*gamma)/(omega+mu)); Rstar=gamma*Istar/(omega+mu)
eq=list(S=Sstar, I=Istar, R=Rstar)
F_m = expression(mu * (1-S) - beta * S * I / N + omega * R)
G = expression(beta * S * I / N - (mu + gamma) * I)
H = expression(gamma*I -(mu +omega)*R)
j11 = D(F_m, "S");j12 = D(F_m, "I");j13 = D(F_m, "R")
j21 = D(G,"S"); j22=D(G, "I"); j23 = D(G, "R")
j31 = D(H,"S"); j32=D(H, "I"); j33 = D(H, "R")
J=with(eq, matrix(c(eval(j11),eval(j12),eval(j13),
                    eval(j21), eval(j22), eval(j23), eval(j31),
                    eval(j32), eval(j33)), nrow=3, byrow=T))
round(eigen(J)$values, 4)
sirs_flu <- function(t,y,parms){ S=y[1]; I=y[2]; R=y[3]
  with(as.list(parms), { inf_term=beta*I*S/N
    dS = mu*(N-S)-inf_term + omega*R; dI=inf_term - (mu+gamma)*I; dR=gamma*I - (mu+omega)*R   # res=c(dS, dE, dI, dR)
    list(c(dS,dI, dR))}) }
# solve
flu_sol <- as.data.frame(ode(c(S=0.299,I=0.001,R=0.7),seq(0,500,by=1/120),sirs_flu,paras))
# plot
ggplot(flu_sol,aes(x=time,y=I)) + geom_line() + scale_x_continuous(breaks=(0:20)*30,expand=expansion(0.01,0)) +
  theme_bw() + standard_theme + geom_vline(data=flu_sol %>% filter(time %in% c(44+(0:10)*49+(0:-10))), # cumsum(50:41)
             aes(xintercept=time),size=1/4,color="red",linetype="dashed")
### ### ### ### ### ### ### ### ###
# sirs seas forcing
sirs_flu_seasforc <- function(t,y,parms){ S=y[1]; I=y[2]; R=y[3]
with(as.list(parms), { inf_term=beta0*(1+beta1*cos((t-phi)/((t_period/2)/pi)))*I*S/N
dS = mu*(N-S)-inf_term + omega*R; dI=inf_term - (mu+gamma)*I; dR=gamma*I - (mu+omega)*R   # res=c(dS, dE, dI, dR)
list(c(dS,dI, dR))}) }
# no forcing: beta        gamma           mu        omega 
#             5.342       1.84          0.000275    0.00481 
paras=c(mu=0.000275,N=1,beta0=5,beta1=0.01,gamma=1.84,omega=1/(52*4),phi=55,t_period=55)
for (k in 1:9) { if (k==1) {list_flu_seasforc_sol=list()}
  paras["beta1"]=((k-1)^2)/250
  flu_seasforc_sol<-as.data.frame(ode(c(S=0.299,I=0.001,R=0.7),times=seq(0,1e4,by=1/5),func=sirs_flu_seasforc,parms=paras))
  list_flu_seasforc_sol[[k]] <- flu_seasforc_sol %>% mutate(beta1=round(paras["beta1"],3))}
# time course plot
ggplot(bind_rows(list_flu_seasforc_sol) %>% filter(time<=2e3),aes(x=time,y=I*100)) + #  %>% filter(time<80) 
  geom_line() + facet_wrap(~beta1,scales="free_y") + theme_bw() + standard_theme + ylab("new infs (% population)")
  # scale_x_continuous(breaks=(0:(max(flu_seasforc_sol$time)/20))*20,expand=expansion(0.01,0))
# phase plot
ggplot(bind_rows(list_flu_seasforc_sol),aes(x=S,y=I,color=time)) + #  %>% filter(time<80) 
  geom_path() + facet_wrap(~beta1,scales="free_y",labeller=label_both) + theme_bw() + standard_theme 
# + ylab("new infs (% population)")

findPeaks <- function (x, thresh = 0) {
  pks <- which(diff(sign(diff(x, na.pad = FALSE)), na.pad = FALSE) < 0) + 2
  if (!missing(thresh)) { pks[x[pks-1] - x[pks] > thresh] }
  else pks}

# bifurc plot
bind_rows(list_flu_seasforc_sol) %>% filter(time>1e3 & time %% 50==0) %>%
  ggplot(aes(x=beta1,y=I,color=time)) + geom_point() + theme_bw() + standard_theme

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### seasonal forcing
seirmod2=function(t, y, parms){
  S=y[1];   E=y[2];   I=y[3];   R=y[4]
  with(as.list(parms), { inf_term=beta0*(1+beta1*cos(2*pi*t))*S*I/N
  dS = mu * (N - S) - inf_term;   dE = inf_term - (mu + sigma) * E
  dI = sigma * E - (mu + gamma) * I;   dR = gamma * I - mu * R   # res=c(dS, dE, dI, dR)
  list(c(dS, dE, dI, dR)) }) }

times=seq(0,150,by=1/60)
paras = c(mu=1/50, N=1, beta0=1000, beta1 = 0.2, sigma = 365/8, gamma = 365/5)
start = c(S=0.16, E=0, I=0.001, R =0.839)
# out=as.data.frame(ode(start, times, seirmod2, paras))
for (k_beta in 0:5) { if (k_beta==0) {list_out=list()}; paras["beta1"]=k_beta/20
  list_out[[k_beta+1]]=as.data.frame(ode(start, times, seirmod2, paras)) %>% mutate(beta1=paras["beta1"]) }
# par(mfrow=c(1,2)) #Side-by-side plot
# plot(times, out$I, ylab="Infected", xlab="Time", xlim=c(90, 100), ylim=c(0,max(out$I[11001:12000])), type="l")
# plot(out$S[11001:12000], out$I[11001:12000], ylab="Infected", xlab="Susceptible", type="l")
ggplot(bind_rows(list_out) %>% filter(time>100),aes(x=time,y=I*100)) + geom_path() + 
  facet_wrap(~beta1,scales="free") + theme_bw() + standard_theme + ylab("new infs as % population")
# at beta=0.2 we get bi-annual oscillations

# steady state
mu = paras["mu"]; N = paras["N"]
beta0 = paras["beta0"]; beta1 = paras["beta1"]; sigma = paras["sigma"]; gamma = paras["gamma"]
R0 = sigma/(sigma + mu) * beta0/(gamma + mu)
Sstar = 1/R0; Istar = mu * (1 - 1/R0) * R0/beta0; Estar = (mu + gamma) * Istar/sigma
eq = list(S = Sstar, E = Estar, I = Istar)
# 
dS = expression(mu * (N - S) - beta0 * S * I / N)
dE= expression(beta0 * S * I / N - (mu + sigma) * E)
dI = expression(sigma*E - (mu + gamma) * I)
j11 = D(dS, "S"); j12 = D(dS, "E"); j13 = D(dS, "I")
j21 = D(dE, "S"); j22 = D(dE, "E"); j23 = D(dE, "I")
j31 = D(dI, "S"); j32 = D(dI, "E"); j33 = D(dI, "I")
J=with(eq, matrix(c(eval(j11),eval(j12),eval(j13),  eval(j21),eval(j22), eval(j23),  
                    eval(j31),eval(j32), eval(j33)), nrow=3, byrow=TRUE))
# eigenvalues
round(eigen(J)$values, 3)
# resonant period
2*pi/(Im(eigen(J)$values[2]))

# bifurc analysis
times = seq(0, 100, by = 1/120)
start = c(S=0.06, E=0, I=0.001, R=0.939)
beta1 = seq(0,0.25, length=101); sel = seq(7001, 12000, by = 120)
#Matrix to store infecteds
#Loop over beta1’s
for(i in 1:101){ if (i==1) {Imat = matrix(NA, ncol=12001, nrow=101); list_out=list()}
  paras = c(mu = 1/50, N = 1, beta0 = 1000, beta1=beta1[i], sigma = 365/8, gamma = 365/5)
  out = as.data.frame(ode(start, times, seirmod2, paras)); Imat[i,] = out$I 
  list_out[[i]]=as.data.frame(ode(start, times, seirmod2, paras)) %>% mutate(beta=beta1[i]) %>% select(c(time,I,beta)) }
df_out = bind_rows(list_out) %>% filter(time %in% times[sel]) %>% select(c(time,I,beta))
# plot
# plot(NA, xlim = range(beta1), ylim = c(1E-7,max(Imat[,sel])), log="y", xlab="beta1",ylab="prevalence")
# for(i in 1:101){   points(rep(beta1[i], length(sel)),Imat[i, sel], pch=20,shape=21) }
ggplot(df_out %>% filter(beta>0.11)) + geom_point(aes(x=beta,y=I,color=factor(round(time,1))),shape=21) + 
  theme_bw() + standard_theme

###
# long-term chaotic dynamics
times = seq(0, 1e3, by = 1/120); start = c(S = 0.06, E = 0, I = 0.001, R = 0.939)
paras = c(mu = 1/50, N = 1, beta0 = 1800, beta1=0.28,sigma = 35.84, gamma = 100)
out=as.data.frame(ode(start, times,seirmod2, paras)); sel=seq(7001/100, 12000000/100, by=120)
par(mfrow=c(1,2))
plot(out$time[7001:13001], out$I[7001:13001],type="l", xlab="Year", ylab="Prevalence")
plot(out$S[sel], out$I[sel], type="p", xlab="S",ylab="I", log="y", pch=20, cex=0.25)

### susceptible recruitment
times = seq(0, 100, by = 1/120)
start = c(S = 0.06, E = 0, I = 0.001, R = 0.939)
mu=seq(from = 0.005, to = 0.02, length = 101)
ImatF=ImatB=matrix(NA, ncol = 12001, nrow = 101)
# loop forward
for(i in 1:101){
  paras=c(mu=mu[i],N=1, beta0=2500,beta1=0.12, sigma = 365/8, gamma=365/5)
  out = as.data.frame(ode(start, times, seirmod2,paras)); if (i %% 10==0) {print(i)}
  ImatF[i,]=out$I; start = c(S=out$S[nrow(out)],E=out$E[nrow(out)],I=out$I[nrow(out)], R=out$R[12001]) }
# loop backward
start = c(S = 0.06, E = 0, I = 0.001, R = 0.939)
for(i in 101:1){
  paras = c(mu = mu[i], N = 1, beta0 = 2500,beta1=0.12, sigma = 365/8, gamma = 365/5)
  out = as.data.frame(ode(start, times, seirmod2,paras)); if (i %% 10==0) {print(i)}
  ImatB[i,]=out$I; start=c(S=out$S[nrow(out)], E = out$E[nrow(out)],I = out $I[nrow(out)], R = out$R[nrow(out)]) }
# selected tpoints
sel=seq(7001, 12000, by=120)
# plot
par(mfrow=c(1,1))
plot(NA, xlim=range(mu), ylim=range(ImatF[,sel]),log="y", xlab="mu", ylab="prevalence")
for(i in 1:101){
  points(rep(mu[i], dim(ImatF)[2]), ImatF[i, ],pch=20, cex=0.25)
  points(rep(mu[i], dim(ImatB)[2]), ImatB[i, ],pch=20, cex=0.25,col=2) }

df_Imat <- bind_rows(data.frame(t(ImatF)) %>% mutate(time=1:ncol(ImatF)) %>% pivot_longer(!time) %>%
  mutate(mu_val=mu[as.numeric(gsub("X","",name))],looptype="fwd"),
  data.frame(t(ImatB)) %>% mutate(time=1:ncol(ImatB)) %>% pivot_longer(!time) %>%
              mutate(mu_val=mu[as.numeric(gsub("X","",name))],looptype="bckwd"))
# timecourse
ggplot(df_Imat %>% filter(time>11e3 & mu_val %in% unique(mu)[(1:10)*10])) + geom_line(aes(x=time,y=value,color=looptype)) +
  facet_wrap(~mu_val) + theme_bw() + standard_theme
# slice
ggplot(df_Imat %>% filter(time %in% sel)) + geom_point(aes(x=mu_val,y=value,color=time),shape=21) + facet_wrap(~looptype) +
  scale_y_log10() + theme_bw()+ standard_theme # facet_wrap(~mu_val)

