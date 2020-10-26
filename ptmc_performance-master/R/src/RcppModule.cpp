#include <Rcpp.h>

using namespace Rcpp;
#include "EvaluateLogLikelihood.h"

RCPP_MODULE(EvaluateLogLikelihoodModule) {
    class_<EvaluateLogLikelihood>( "EvaluateLogLikelihood" )
    .constructor<double, double, NumericVector>()
    .field( "contactMatrix", &EvaluateLogLikelihood::contactMatrix )
    .field( "observedData", &EvaluateLogLikelihood::observedData )
    .method( "evaluateLogLikelihoodCpp", &EvaluateLogLikelihood::evaluateLogLikelihoodCpp )
    ;
}
