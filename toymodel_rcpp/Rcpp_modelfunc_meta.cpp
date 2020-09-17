#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

// Sinusoidal function to describe temperature and biting rate

// [[Rcpp::export]]
double sinefunc(double t, double c, double d, double phi) {
  double pi = M_PI;
  return c + d*sin((2*pi/365.0)*(t-phi));
}


// Sigmoid and parabolic function to describe flea life cycle parameters

// [[Rcpp::export]]
double sigmoidfunc(double temp, double asym, double xmid, double scale) {
  return asym/((1+exp((xmid-temp)/scale)));
}

// [[Rcpp::export]]
double polyfunc(double temp, double a2, double a1, double a0) {
  return a2*temp*temp + a1*temp + a0;
  
}



// [[Rcpp::export]]
List diseasemodel_meta(double time, NumericVector logstate, NumericVector par, NumericMatrix M) {
  
  NumericVector state = exp(logstate); // exponentiate the log states
  List rparam(par);
  NumericVector dxdt(state.length());
  
  // temperature params
  int dintro = as<int>(rparam["dintro"]);
  double cT =  as<double>(rparam["cT"]);
  double dT =  as<double>(rparam["dT"]);
  double phiT =  as<double>(rparam["phiT"]);
  double u =  as<double>(rparam["u"]); // scaling factor for temperature
  double temp = sinefunc(dintro + time-1.0, cT, dT, phiT) + u;
  
  // params for flea life cycle functions
  double asymr = as<double>(rparam["asymr"]);
  double asymeps = as<double>(rparam["asymeps"]);
  double xmidr = as<double>(rparam["xmidr"]);
  double xmideps = as<double>(rparam["xmideps"]);
  double scalr = as<double>(rparam["scalr"]);
  double scaleps = as<double>(rparam["scaleps"]);
  double a2y = as<double>(rparam["a2y"]);
  double a1y = as<double>(rparam["a1y"]);
  double a0y = as<double>(rparam["a0y"]);
  double a2s = as<double>(rparam["a2s"]);
  double a1s = as<double>(rparam["a1s"]);
  double a0s = as<double>(rparam["a0s"]);
  
  // calculate reproduction, mortality and transition rates at current temperature T
  double r = sigmoidfunc(temp, asymr, xmidr, scalr);
  double epsilon = sigmoidfunc(temp, asymeps, xmideps, scaleps);
  double muvy = polyfunc(temp, a2y, a1y, a0y);
  double muvs = polyfunc(temp, a2s, a1s, a0s);
  
  // other params
  double muh = as<double>(rparam["muh"]); // host birth and mortality rate
  double f = as<double>(rparam["f"]); // proportion female vectors
  double K = as<double>(rparam["K"]); // pre-adult vector carrying capacity
  double bhv = as<double>(rparam["bhv"]); // transmission rate from vector to human
  double bvh = as<double>(rparam["bvh"]); // transmission rate from high-infectious human to vector
  double bhp = as<double>(rparam["bhp"]); // transmission rate from pneumonic cases
  double omega = as<double>(rparam["omega"]); // transition from exposed compartment to Recovery, Septicaemia or secondary pneumonic plague
  double gammas = as<double>(rparam["gammas"]); // transition from septicaemic compartment to death
  double gammap = as<double>(rparam["gammap"]); // transition from pneumonic compartment to death
  double delta = as<double>(rparam["delta"]); // recovery rate of vector
  double p = as<double>(rparam["p"]); // probabilty of host dying from bubonic plague
  double q = as<double>(rparam["q"]); // probability of host developing secondary pneumonic plague
  int npatch = as<int>(rparam["npatch"]); // Number of patches
  double tau = as<double>(rparam["tau"]); // probability of host leaving patch
  
  // sine-cosine biting rate
  double ca =  as<double>(rparam["ca"]);
  double da =  as<double>(rparam["da"]);
  double phia =  as<double>(rparam["phia"]);
  double a = sinefunc(dintro + time-1.0, ca+da, da, phia);
  
  // States in all patches
  enum  state_variables {Sh, Eh, Ihs, Rh, Ihp, Ihpsec, Dh, Ihcumb, Ihcump, Yv, Sv, Iv};
  NumericVector lambdav(npatch);
  NumericVector lambdahb(npatch);
  NumericVector lambdahp(npatch);
  
  for(int i=0 ; i<npatch ; i++) {
    
    // calculate number of living host population in patch
    int Nh = state[12*i+Sh]+ state[12*i+Eh] + state[12*i+Ihs] + state[12*i+Rh] + state[12*i+Ihp] + state[12*i+Ihpsec];
    
    // FOI to vectors
    lambdav[i] = a*bvh*state[12*i+Ihs]/Nh;
    
    // FOI to humans within same patch
    double lambdahbint = (a*bhv*state[12*i+Iv])/Nh;
    double lambdahpint = (bhp*(state[12*i+Ihp]+state[12*i+Ihpsec]))/Nh;
    
    // FOI to humans from other patches
    double lambdahbext = 0.0;
    double lambdahpext = 0.0;
    for(int ii=0 ; ii<npatch ; ii++) {
      
      lambdahbext += a*bhv*state[12*ii+Iv]*M(ii,i)/Nh;
      lambdahpext += bhp*(state[12*ii+Ihp] + state[12*ii+Ihpsec])*M(ii,i)/Nh;
  
    }
    
    // Total FOI to humans for patch
    lambdahb[i] = lambdahbint + tau*lambdahbext;
    lambdahp[i] = lambdahpint + tau*lambdahpext;
    

    // transitions
    dxdt[12*i+Sh] = -(lambdahb[i]+lambdahp[i])*state[12*i+Sh] + muh*Nh - muh*state[12*i+Sh] + gammas*state[12*i+Ihs] + gammap*(state[12*i+Ihp]+state[12*i+Ihpsec]); // Sh
    dxdt[12*i+Eh] = lambdahb[i]*state[12*i+Sh] - omega*state[12*i+Eh] - muh*state[12*i+Eh]; // Eh
    dxdt[12*i+Ihs] = (p-q)*omega*state[12*i+Eh] - gammas*state[12*i+Ihs] - muh*state[12*i+Ihs]; // Ihs
    dxdt[12*i+Rh] = (1-p)*omega*state[12*i+Eh]- muh*state[12*i+Rh]; // Rh
    
    dxdt[12*i+Ihp] = lambdahp[i]*state[12*i+Sh] - muh*state[12*i+Ihp] - gammap*state[12*i+Ihp]; // Ihp
    dxdt[12*i+Ihpsec] = q*omega*state[12*i+Eh] - muh*state[12*i+Ihpsec] - gammap*state[12*i+Ihpsec]; // Ihpsec
    
    dxdt[12*i+Dh] = gammas*state[12*i+Ihs] + gammap*(state[12*i+Ihp]+state[12*i+Ihpsec]); // Dh
    
    dxdt[12*i+Ihcumb] = lambdahb[i]*state[12*i+Sh]; // bubonic cumulative incidence in hosts
    dxdt[12*i+Ihcump] = lambdahp[i]*state[12*i+Sh]; // pneumonic cumulative incidence in hosts
    
    dxdt[12*i+Yv] = r*f*state[12*i+Sv]*(1-(state[12*i+Yv]/(K*Nh))) - epsilon*state[12*i+Yv] - muvy*state[12*i+Yv]; // Yv
    dxdt[12*i+Sv] = epsilon*state[12*i+Yv] + delta*state[12*i+Iv] - lambdav[i]*state[12*i+Sv] - muvs*state[12*i+Sv]; // Sv
    dxdt[12*i+Iv] = lambdav[i]*state[12*i+Sv] - muvs*state[12*i+Iv] - delta*state[12*i+Iv]; // Iv
    

  }
  
  return List::create(dxdt/state,
                      Named("lambdahb")=lambdahb, 
                      Named("lambdahp")=lambdahp, 
                      Named("lambdav")=lambdav,
                      Named("temp")=temp, 
                      Named("r")=r, 
                      Named("epsilon")=epsilon, 
                      Named("muvy")=muvy, 
                      Named("muvs")=muvs,
                      Named("a")=a);
  
}

