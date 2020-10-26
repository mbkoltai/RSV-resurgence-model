#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List SIS_ode_cpp(NumericVector times, vec state, List parms) {
  
  mat betaVT = as<mat>(parms["betaVT"]);
  mat betaNVT = as<mat>(parms["betaNVT"]);
  vec clearVT = as<vec>(parms["clearVT"]);
  vec clearNVT = as<vec>(parms["clearNVT"]);
  double comp = as<double>(parms["comp"]);
  vec pop = as<vec>(parms["Population"]);
  double noagegps = as<double>(parms["no.agegps"]);
  vec ageout = as<vec>(parms["ageout"]); 
  vec agein = as<vec>(parms["agein"]);
  
  vec S = state.subvec(0,noagegps-1);
  vec VT = state.subvec(noagegps,2*noagegps-1);
  vec NVT = state.subvec(2*noagegps,3*noagegps-1);
  vec B = state.subvec(3*noagegps,4*noagegps-1);
  vec N = S + VT + NVT + B;
  
  vec FOIVT = betaVT * ((VT+B)/N % pop); 
  vec FOINVT = betaNVT * ((NVT+B)/N % pop);
  
  vec SAge(noagegps), VTAge(noagegps), NVTAge(noagegps), BAge(noagegps);
  SAge(0)=sum(N.subvec(noagegps-1,noagegps-1)); VTAge(0)=0.0; NVTAge(0)=0.0; BAge(0)=0.0; /* .subvec(noagegps-1,noagegps-1) */
    SAge.subvec(1,noagegps-1) = S.subvec(0,noagegps-2);
    VTAge.subvec(1,noagegps-1) = VT.subvec(0,noagegps-2);
    NVTAge.subvec(1,noagegps-1) = NVT.subvec(0,noagegps-2);
    BAge.subvec(1,noagegps-1) = B.subvec(0,noagegps-2);
    
    mat res(noagegps,4);
    res.col(0) = -FOIVT%S - FOINVT%S +clearVT%VT +clearNVT%NVT -ageout%S +agein%SAge;
    res.col(1) = FOIVT%S -comp*FOINVT%VT  -clearVT%VT +clearNVT%B -ageout%VT +agein%VTAge;
    res.col(2) = FOINVT%S -comp*FOIVT%NVT -clearNVT%NVT +clearVT%B - ageout%NVT +agein%NVTAge;
    res.col(3) = comp*FOIVT%NVT +comp*FOINVT%VT -clearVT%B -clearNVT%B - ageout%B + agein%BAge;
    
    return List::create(res);
}