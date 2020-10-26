#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List toymodel(double time, NumericVector logstate, NumericVector par) {

  NumericVector state = exp(logstate); // exponentiate the log states
  List rparam(par);
  NumericVector dxdt(state.length());

  // params 
  double beta = as<double>(rparam["beta"]); 
  double gamma = as<double>(rparam["gamma"]);
  double mu = as<double>(rparam["mu"]);
  double sigma = as<double>(rparam["sigma"]);
  double eta = as<double>(rparam["eta"]);
  double phi = as<double>(rparam["phi"]);
  
  // states
  double S = state[0];   double I = state[1];  double R = state[2];  
  double N = S + I + R;
  
  // Seasonally forced FOI
  double beta_eff = beta * (1+eta*pow(sin(M_PI*(time-phi)/365.0),2)); // can also be another sinusoidal function
  double FOI = beta_eff*I/N;
 
  // change in states
  dxdt[0] = -FOI*S + mu*N - mu*S + sigma*R; 
  dxdt[1] = FOI*S - gamma*I - mu*I; 
  dxdt[2] = gamma*I - mu*R - sigma*R; 

  return List::create(dxdt/state); 
}