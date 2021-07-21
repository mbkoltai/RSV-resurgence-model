#include <Rcpp.h>

using namespace Rcpp;
#include "EvaluateLogLikelihood.h"

RCPP_MODULE(EvaluateLogLikelihoodModule) {
    class_<EvaluateLogLikelihood>( "EvaluateLogLikelihood" )
    .constructor<double, double, NumericVector>()
    .field( "contactMatrixPhy", &EvaluateLogLikelihood::contactMatrixPhy )
    .field( "contactMatrixCon", &EvaluateLogLikelihood::contactMatrixCon )
    .field( "observedData", &EvaluateLogLikelihood::observedData )
    .field( "lowerParamSupport", &EvaluateLogLikelihood::lowerParamSupport )
    .field( "upperParamSupport", &EvaluateLogLikelihood::upperParamSupport )
    .field( "run_full", &EvaluateLogLikelihood::run_full )

    .method( "evaluateLogLikelihoodCpp", &EvaluateLogLikelihood::evaluateLogLikelihoodCpp )
    .method( "getWeeklySampleCpp", &EvaluateLogLikelihood::getWeeklySampleCpp )
    .method( "getAnnualIncidenceCpp", &EvaluateLogLikelihood::getAnnualIncidenceCpp )
    .method( "getProportionBornProtectedCpp", &EvaluateLogLikelihood::getProportionBornProtectedCpp )

    ;
}
