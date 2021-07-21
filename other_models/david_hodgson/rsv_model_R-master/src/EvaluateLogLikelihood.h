//
//  Created by David Hodgson on 03/02/2020.
//  Copyright © 2020 David Hodgson. All rights reserved.
//

#include <Rcpp.h>
#include <RcppEigen.h>
#include <random>
#include <boost/random.hpp>
#include <boost/math/distributions.hpp>
#include "ascent/Ascent.h"
//#include "epmgp.h"

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins("cpp14")]]
using namespace Rcpp;
using namespace std;
using namespace Eigen;
using namespace asc;

// Stuff for random number generation
std::random_device dev;
std::mt19937 engine(dev());
typedef boost::mt19937 PRNG_s;
PRNG_s rng(engine());


class EvaluateLogLikelihood
{
public:
    // Need to be defined in constructor
    int A;
    vector< double >  populationPerAgeGroup, eta, modelIncidencePerTime, pA, ep_t;
    double dailyBirthRate, totPopulation;
    double run_start, run_burn, run_oneyr, run_full, dt;
    NumericVector ageStratification;
    
    int dayNoAfterBurn,  weekNo;
    double valueLogLikelihood;
    
    // Defined later but of size A
    double currentODETime;
    bool loglikelihoodError;
    
    EvaluateLogLikelihood(double dailyBirthRate_t, double totPopulation_t, NumericVector ageStratification_t): dailyBirthRate(dailyBirthRate_t), totPopulation(totPopulation_t), ageStratification(ageStratification_t){
       
        A = ageStratification.size();
        eta.push_back(0);
        for (int i = 0; i < A-1; i++){
            populationPerAgeGroup.push_back(dailyBirthRate*365*(ageStratification[i+1]-ageStratification[i]));
            eta.push_back( 1.0/(365.0*(ageStratification[i+1] - ageStratification[i])) );
            modelIncidencePerTime.push_back(0);
        }            
        modelIncidencePerTime.push_back(0);
        populationPerAgeGroup.push_back(totPopulation - (dailyBirthRate*365)*ageStratification[A-1]);
        eta.push_back(dailyBirthRate/(totPopulation - (dailyBirthRate*365)*ageStratification[A-1]) );

        run_start = 0;
        run_burn = 52*7 + 1;
        run_oneyr = 52*7 + run_burn;
        dt = 1;
        currentODETime = 0;
        dayNoAfterBurn = 0;
        weekNo = 0;
        valueLogLikelihood = 0;
        
        loglikelihoodError = false;
    }
    
    double evaluateLogLikelihoodCpp(VectorXd currentParamValues);
    NumericMatrix getWeeklySampleCpp(VectorXd currentParamValues, bool epFlag);
    NumericVector getAnnualIncidenceCpp(VectorXd currentParamValues);
    NumericVector getProportionBornProtectedCpp(VectorXd currentParamValues);
    // define in R after construction
    NumericMatrix contactMatrixPhy;
    NumericMatrix contactMatrixCon;
    NumericMatrix observedData;
    vector< double > lowerParamSupport, upperParamSupport;

    // Related to parameter values
    NumericVector parameterValues;
    
    void ParameterValuesforODE(VectorXd currentParamValues){

        this->ep_t.clear();
        this->pA.clear();
        
        NumericVector parameterValuesTemp(25);
        
         parameterValuesTemp["xi"] = currentParamValues(0);
         parameterValuesTemp["si"] = currentParamValues(1);
         parameterValuesTemp["ga0"] = currentParamValues(2);
         parameterValuesTemp["g1"] = currentParamValues(3);
         parameterValuesTemp["g2"] = currentParamValues(4);
         parameterValuesTemp["om"] = currentParamValues(5);
         parameterValuesTemp["pA1"] = currentParamValues(6);
         parameterValuesTemp["pA2"] = currentParamValues(7);
         parameterValuesTemp["pA3"] = currentParamValues(8);
         parameterValuesTemp["pA4"] = currentParamValues(9);
         parameterValuesTemp["alpha_i"] =currentParamValues(10);
         parameterValuesTemp["d1"] = currentParamValues(11);
         parameterValuesTemp["d2"] = currentParamValues(12);
         parameterValuesTemp["d3"] = currentParamValues(13);
         parameterValuesTemp["phi"] = currentParamValues(14);
         parameterValuesTemp["qp"] = currentParamValues(15);
         parameterValuesTemp["qc"] = currentParamValues(16);
         parameterValuesTemp["b1"] = currentParamValues(17);
         parameterValuesTemp["psi"] = currentParamValues(18);
         parameterValuesTemp["c5ep1"] = currentParamValues(19);
         parameterValuesTemp["c5ep2"] = currentParamValues(20);
         parameterValuesTemp["ep5"] = currentParamValues(21);
         parameterValuesTemp["ep6"] = currentParamValues(22);
         parameterValuesTemp["I1"] = currentParamValues(23);
         parameterValuesTemp["I2"] = currentParamValues(24);
        this->parameterValues = parameterValuesTemp;
        
        // Define ep_t;
        for (int a = 0; a < 16; a++)
            this->ep_t.push_back(exp(parameterValues["c5ep1"] + a*parameterValues["c5ep2"]));
        for (int a = 16; a < 23; a++)
            this->ep_t.push_back(parameterValues["ep5"]);
        for (int a = 23; a < 25; a++)
            this->ep_t.push_back(parameterValues["ep6"]);
        // define pA;
        for (int a = 0; a < 12; a++)
            this->pA.push_back(parameterValues["pA1"]);
        for (int a = 12; a < 16; a++)
            this->pA.push_back(parameterValues["pA2"]);
        for (int a = 16; a < 18; a++)
            this->pA.push_back(parameterValues["pA3"]);
        for (int a = 18; a < 25; a++)
            this->pA.push_back(parameterValues["pA4"]);
    }
    
    vector< double >  generateInitialStates(){
        vector< double >  initialStates;
        
        vector<double > populationMatPro = initial_M();
        
        double a1, a2;
        double I1 = this->parameterValues["I1"];
        double I2 = this->parameterValues["I2"];
        double I3 = 0.5;
        
        double si = 1.0/parameterValues["si"];
        double g0 = 1.0/parameterValues["ga0"];
        double g1 = 1.0/((parameterValues["ga0"])*(parameterValues["g1"]));
        double g2 = 1.0/((parameterValues["ga0"])*(parameterValues["g1"])*(parameterValues["g2"]));
        double d1 = parameterValues["d1"];
        double d2 = parameterValues["d1"]*parameterValues["d2"];
        double d3 = parameterValues["d1"]*parameterValues["d2"]*parameterValues["d3"];
        
        for (int a = 0; a < this->A ; a++){
            if (a < A-1){
                a1 = this->ageStratification(a); a2 =  this->ageStratification(a+1);
            }
            else{
                a1 =  this->ageStratification(a); a2 = 90;
            }
            
            vector< double >  propEachExposureGroup = initialProportionExposure(I3, a1, a2);
            double populationNoMatPro = populationPerAgeGroup[a] - populationMatPro[a];

            initialStates.push_back(populationMatPro[a]);

            initialStates.push_back(populationNoMatPro*propEachExposureGroup[0]*(1-I1)*(1-I2));
            initialStates.push_back(populationNoMatPro*propEachExposureGroup[0]*I1*si/(si+g0));
            initialStates.push_back(populationNoMatPro*propEachExposureGroup[0]*I1*g0/(si+g0)*pA[a]);
            initialStates.push_back(populationNoMatPro*propEachExposureGroup[0]*I1*g0/(si+g0)*(1-pA[a]));
            initialStates.push_back(populationNoMatPro*propEachExposureGroup[0]*(1-I1)*I2);
            
            initialStates.push_back(populationNoMatPro*propEachExposureGroup[1]*(1-d1*I1)*(1-I2));
            initialStates.push_back(populationNoMatPro*propEachExposureGroup[1]*d1*I1*si/(si+g1));
            initialStates.push_back(populationNoMatPro*propEachExposureGroup[1]*d1*I1*g1/(si+g1)*pA[a]);
            initialStates.push_back(populationNoMatPro*propEachExposureGroup[1]*d1*I1*g1/(si+g1)*(1-pA[a]));
            initialStates.push_back(populationNoMatPro*propEachExposureGroup[1]*(1-d1*I1)*I2);
            
            initialStates.push_back(populationNoMatPro*propEachExposureGroup[2]*(1-d2*I1)*(1-I2));
            initialStates.push_back(populationNoMatPro*propEachExposureGroup[2]*d2*I1*si/(si+g2));
            initialStates.push_back(populationNoMatPro*propEachExposureGroup[2]*d2*I1*g2/(si+g2)*pA[a]);
            initialStates.push_back(populationNoMatPro*propEachExposureGroup[2]*d2*I1*g2/(si+g2)*(1-pA[a]));
            initialStates.push_back(populationNoMatPro*propEachExposureGroup[2]*(1-d2*I1)*I2);
            
            initialStates.push_back(populationNoMatPro*propEachExposureGroup[3]*(1-d3*I1)*(1-I2));
            initialStates.push_back(populationNoMatPro*propEachExposureGroup[3]*d3*I1*si/(si+g2));
            initialStates.push_back(populationNoMatPro*propEachExposureGroup[3]*d3*I1*g2/(si+g2)*pA[a]);
            initialStates.push_back(populationNoMatPro*propEachExposureGroup[3]*d3*I1*g2/(si+g2)*(1-pA[a]));
            initialStates.push_back(populationNoMatPro*propEachExposureGroup[3]*(1-d3*I1)*I2);
            initialStates.push_back(0);
            initialStates.push_back(0);

        }
        return initialStates;
    }
    
    vector<double > initial_M()
    {
        double xi = 1.0/this->parameterValues["xi"];
        boost::math::exponential_distribution <> exp ( xi );
        
        vector<double > init_con;
        for (int i = 0; i< this->A - 1; i++)
        {
            double init_con_temp = (cdf(exp, 365*ageStratification[i+1])-cdf(exp, 365*ageStratification[i]))/((365*ageStratification[i+1]-365*ageStratification[i])*xi);
            init_con.push_back(init_con_temp*populationPerAgeGroup[i]);
        }
        init_con.push_back( (cdf(exp, 365*90)-cdf(exp, 365*ageStratification[A-1]))/((365*90-365*ageStratification[A-1])*xi)*populationPerAgeGroup[A-1]);
        return init_con;
    }
    
    vector< double >  initialProportionExposure(double l, double a1, double a2){
        vector< double >  prop(this->A);
        prop[0] = abs(poisson_cdf(l,a2,0)-poisson_cdf(l,a1,0))/((a2-a1)*l);
        prop[1] = abs(poisson_cdf(l,a2,1)-poisson_cdf(l,a1,1))/((a2-a1)*l);
        prop[2] = abs(poisson_cdf(l,a2,2)-poisson_cdf(l,a1,2))/((a2-a1)*l);
        prop[3] = 1 - (prop[2]+prop[1]+prop[0]);
        return prop;
    }
    
    double poisson_cdf(double l, double a, double x){
      if( l == 0.0 || a == 0.0){
        boost::math::poisson_distribution<> p(0.000001); return cdf(p,x);
      }
      else{
        boost::math::poisson_distribution<> p(l*a); return cdf(p,x);
      }
    }
    
    double evaluateTimeStepLogLikelihood(int weekNo){
        double ll = 0;
        double estimatedLogBinomialCoeff;
        for (int a = 0; a < this->A; a++){
            double dataNewInfections = this->observedData(weekNo, a);
            if (dataNewInfections > this->modelIncidencePerTime[a]){
                this->loglikelihoodError=true;
                return log(0);
            }
            
            estimatedLogBinomialCoeff = stirlingApproximation(this->modelIncidencePerTime[a]) - stirlingApproximation(dataNewInfections) - stirlingApproximation(this->modelIncidencePerTime[a]-dataNewInfections);
            ll += estimatedLogBinomialCoeff + dataNewInfections*log(this->ep_t[a]) + (this->modelIncidencePerTime[a]-dataNewInfections)*log(1-this->ep_t[a]);
        }
        if (std::isinf(ll) || std::isnan(ll)){
          this->loglikelihoodError=true;
        }
        return ll;
    }
    
    long double stirlingApproximation(double n){
        if (n == 0)
            return 0;
        else {
            double x = n + 1;
            return (x - 0.5)*log(x) - x + 0.5*log(2*PI) + 1.0/(12*x) - 1.0/(360.0*pow(x,3)); // https://en.wikipedia.org/wiki/Stirling%27s_approximation#Speed_of_convergence_and_error_estimates
        }
    }
    
    inline void checkStability(vector< double >  x)
    {
        long double X_w = 0;
        // check state values are never negative
        for (int j = 0; j < this->A*23; j++){
            if (x[j] < 0) {
              this->loglikelihoodError=true;
              return; }
            else {X_w += x[j]; }
        }
        // Check state values are non-infinte numeric values
        if (std::isinf(X_w) || std::isnan(X_w)) {
            this->loglikelihoodError=true;
            return;
        }
    }
    
    void getWeeklyLikelihood(vector<double> &x0, VectorXd& inc_tot )
    {
        double ll = 0;
        if (this->dayNoAfterBurn == 0){
            for (int a = 0; a < this->A; a++)
                x0[22 + 23*a] = 0.0; //Incidence at t_d = 0;
        }
        if (this->dayNoAfterBurn%7 == 0 && this->dayNoAfterBurn > 0)
        {
            for (int a = 0; a < this->A; a++){
                this->modelIncidencePerTime[a] = x0[22 + 23*a]; //Incidence at t_d = 7;
                inc_tot(a) += this->modelIncidencePerTime[a];
                x0[22 + 23*a] = 0.0;
            }
            this->valueLogLikelihood += this->evaluateTimeStepLogLikelihood(this->weekNo);
            if (this->loglikelihoodError){return;}
            this->weekNo ++;
        }
        if (this->dayNoAfterBurn%365 == 0 && this->currentODETime > 0) {
          //  ll += check_incidence(inc_tot);
            if (std::isinf(ll)){
                if (this->loglikelihoodError){return;}
            }
            inc_tot = VectorXd::Zero(this->A);
        }
        this->dayNoAfterBurn++;
    }
    
    void getWeeklyIncidence(vector<double> &x0, NumericMatrix &sampleWeeklyIncidence, bool epFlag)
    {
      if (this->dayNoAfterBurn == 0){
        for (int a = 0; a < this->A; a++)
          x0[22 + 23*a] = 0.0; //Incidence at t_d = 0;
      }
      if (this->dayNoAfterBurn%7 == 0 && this->dayNoAfterBurn > 0)
      {
        for (int a = 0; a < this->A; a++){
          if (epFlag)
            sampleWeeklyIncidence(this->weekNo, a) = x0[22 + 23*a]*this->ep_t[a]; //Incidence at t_d = 7;
          else
            sampleWeeklyIncidence(this->weekNo, a) = x0[22 + 23*a]; //Incidence at t_d = 7;
          
          x0[22 + 23*a] = 0.0;
        }
        this->weekNo ++;
      }
      this->dayNoAfterBurn++;
    }
    
    void getAnnualIncidence(vector<double> &x0, NumericVector &sampleAnnualIncidence)
    {
        if (this->dayNoAfterBurn == 0){
            for (int a = 0; a < this->A; a++)
            x0[22 + 23*a] = 0.0; //Incidence at t_d = 0;
        }
        if (this->dayNoAfterBurn%7 == 0 && this->dayNoAfterBurn > 0)
        {
            for (int a = 0; a < this->A; a++){
                sampleAnnualIncidence[a] += x0[22 + 23*a]; //Incidence at t_d = 7;
                x0[22 + 23*a] = 0.0;
            }
            this->weekNo ++;
        }
        this->dayNoAfterBurn++;
    }
    
    double getProportionBornProtected(vector<double> &x)
    {
        vector<double > num_vec_wcba;
        double sum_wcb = 0.0;
        double CB2_temp = 0.0;
        for(int j = 18; j < 21 ; j++)
        {
            CB2_temp = (x[1 + j*23] + x[6 + j*23] + x[11 + j*23] + x[16 + j*23] + x[2 + j*23] + x[7 + j*23] + x[12 + j*23] + x[17 + j*23])/(double)this->populationPerAgeGroup[j];
            num_vec_wcba.push_back(CB2_temp/3.0);
        }
        for(int i=0; i < num_vec_wcba.size() ;i++){sum_wcb = sum_wcb + num_vec_wcba[i];}
        
        return sum_wcb;
    }
    
    inline double check_incidence(VectorXd inc_tot)
    {
        double pl = 0;
        for (int a = 0; a < 25; a++)
        {
            if (inc_tot(a) > populationPerAgeGroup[a]*0.8){
                return log(0);
            }
        }
        return pl;
    }
};

class ODE_desc
{
  EvaluateLogLikelihood* finELL;
  
public:
  double M, S0, S1, S2, S3, E0, E1, E2, E3, A0, A1, A2, A3, I0, I1, I2, I3, R1, R2, R3, N;
  double xi, si, ga0, ga1, ga2, ga3, d1, d2, d3, a1, a2, a3, alpha_i, rho, om, b1, qp, qc, psi, phi, beta;
  double dailyBirthRate;
  NumericMatrix contactMatrixPhy, contactMatrixCon;
  int A;
  vector< double >  populationPerAgeGroup, eta, pA;
  
  ODE_desc(EvaluateLogLikelihood* finELL_t): finELL(finELL_t){
    
    xi = 1.0/finELL->parameterValues["xi"];
    si = 1.0/finELL->parameterValues["si"];
    ga0 = 1.0/(finELL->parameterValues["ga0"]);
    ga1 = 1.0/(finELL->parameterValues["ga0"]*finELL->parameterValues["g1"]);
    ga2 = 1.0/(finELL->parameterValues["ga0"]*finELL->parameterValues["g1"]*finELL->parameterValues["g2"]);
    ga3 = ga2;
    om = 1.0/finELL->parameterValues["om"];
    
    rho = 1.0;
    alpha_i = finELL->parameterValues["alpha_i"];
    d1 = finELL->parameterValues["d1"];
    d2 = finELL->parameterValues["d1"]*finELL->parameterValues["d2"];
    d3 = finELL->parameterValues["d1"]*finELL->parameterValues["d2"]*finELL->parameterValues["d3"];
    a1 = 1.0, a2 = 1.0, a3 = 1.0;

    phi = finELL->parameterValues["phi"];
    qp = finELL->parameterValues["qp"];
    qc = finELL->parameterValues["qc"];
    b1 = finELL->parameterValues["b1"];
    psi = finELL->parameterValues["psi"];
      
    A = finELL->A;
    eta = finELL->eta;
    populationPerAgeGroup = finELL->populationPerAgeGroup;
    contactMatrixPhy = finELL->contactMatrixPhy;
    contactMatrixCon = finELL->contactMatrixCon;
    dailyBirthRate = finELL->dailyBirthRate;
    pA = finELL->pA;
    
   };
  
  
  void operator() (  vector< double >  &x , vector< double >  &dxdt , const double  t )
  {
      int t_d = (int)t%365;
      beta = 1 + b1*(1 + exp(-((t_d/365.0 - phi))*((t_d/365.0 - phi))/(2*psi*psi)));
      vector<double > num_vec_wcba;
      double sum_wcb = 0.0;
      for(int j = 18; j < 21 ; j++)
      {
          double CB2_temp = (x[1 + j*23] + x[6 + j*23] + x[11 + j*23] + x[16 + j*23] + x[2 + j*23] + x[7 + j*23] + x[12 + j*23] + x[17 + j*23])/(double)populationPerAgeGroup[j];
          num_vec_wcba.push_back(CB2_temp/3.0);
      }
      for(int i=0; i < num_vec_wcba.size() ;i++){sum_wcb = sum_wcb + num_vec_wcba[i];}
      double p_vul = sum_wcb;

      
    for (int a = 0; a < A ; a++)
    {
        double mu = a==0 ? dailyBirthRate : 0;
      
        double I_temp = 0;
        for (int b = 0; b < A ; b++)
        {
            I_temp += (x[3 + 23*b]*alpha_i+x[4+23*b]+a1*(x[8+23*b]*alpha_i+x[9+23*b])+a2*(x[13+ 23*b]*alpha_i+x[14+23*b])+a3*(x[18+23*b]*alpha_i+x[19+23*b]))*(qp*(contactMatrixPhy(b,a)+qc*contactMatrixCon(b,a)))/((double)populationPerAgeGroup[b]);
        }
        
        int pa = max(23*(a-1), 0);
        
        dxdt[0 + 23*a] = (1.0-p_vul)*mu - x[0 + 23*a]*xi                                   - x[0+23*a]*eta[a+1] + x[0+pa]*eta[a];
        
        dxdt[1 + 23*a] =  p_vul*mu + x[0 + 23*a]*xi       - x[1 + 23*a]*I_temp*beta        - x[1+23*a]*eta[a+1] + x[1+pa]*eta[a];
        
        dxdt[2 + 23*a] = x[1 + 23*a]*I_temp*beta               - x[2+23*a]*si                - x[2+23*a]*eta[a+1] + x[2+pa]*eta[a];
        dxdt[3 + 23*a] = x[2 + 23*a]*si*pA[a]                - x[3+23*a]*ga0*rho             - x[3+23*a]*eta[a+1] + x[3+pa]*eta[a];
        dxdt[4 + 23*a] = x[2 + 23*a]*si*(1.0-pA[a])            - x[4+23*a]*ga0               - x[4+23*a]*eta[a+1] + x[4+pa]*eta[a];
        dxdt[5 + 23*a] = x[4 + 23*a]*ga0 + x[3+23*a]*ga0*rho        - x[5+23*a]*om           - x[5+23*a]*eta[a+1] + x[5+pa]*eta[a];
        
        dxdt[6 + 23*a] = x[5+23*a]*om                       - d1*x[6+23*a]*I_temp*beta       - x[6+23*a]*eta[a+1] + x[6+pa]*eta[a];
        dxdt[7 + 23*a] = d1*x[6+23*a]*I_temp*beta            - x[7+23*a]*si                  - x[7+23*a]*eta[a+1] + x[7+pa]*eta[a];
        dxdt[8 + 23*a] = x[7+23*a]*si*pA[a]                - x[8+23*a]*ga1*rho               - x[8+23*a]*eta[a+1] + x[8+pa]*eta[a];
        dxdt[9 + 23*a] = x[7+23*a]*si*(1.0-pA[a])            - x[9+23*a]*ga1                 - x[9+23*a]*eta[a+1] + x[9+pa]*eta[a];
        dxdt[10 + 23*a] = x[9+23*a]*ga1 + x[8+23*a]*ga1*rho       - x[10+23*a]*om            - x[10+23*a]*eta[a+1] + x[10+pa]*eta[a];
        
        dxdt[11 + 23*a] = x[10+23*a]*om                     - d2*x[11+23*a]*I_temp*beta      - x[11+23*a]*eta[a+1] + x[11+pa]*eta[a];
        dxdt[12 + 23*a] = d2*x[11+23*a]*I_temp*beta           - x[12+23*a]*si                - x[12+23*a]*eta[a+1] + x[12+pa]*eta[a];
        dxdt[13 + 23*a] = x[12+23*a]*si*pA[a]              - x[13+23*a]*ga2*rho              - x[13+23*a]*eta[a+1] + x[13+pa]*eta[a];
        dxdt[14 + 23*a] = x[12+23*a]*si*(1.0-pA[a])           - x[14+23*a]*ga2               - x[14+23*a]*eta[a+1] + x[14+pa]*eta[a];
        dxdt[15 + 23*a] = x[14+23*a]*ga2 + x[13+23*a]*ga2*rho       - x[15+23*a]*om          - x[15+23*a]*eta[a+1] + x[15+pa]*eta[a];
        
        dxdt[16 + 23*a] = x[15+23*a]*om + x[20+23*a]*om     - d3*x[16+23*a]*I_temp*beta      - x[16+23*a]*eta[a+1] + x[16+pa]*eta[a];
        dxdt[17 + 23*a] = d3*x[16+23*a]*I_temp*beta           - x[17+23*a]*si                - x[17+23*a]*eta[a+1] + x[17+pa]*eta[a];
        dxdt[18 + 23*a] = x[17+23*a]*si*pA[a]               - x[18+23*a]*ga3*rho             - x[18+23*a]*eta[a+1] + x[18+pa]*eta[a];
        dxdt[19 + 23*a] = x[17+23*a]*si*(1.0-pA[a])           - x[19+23*a]*ga3               - x[19+23*a]*eta[a+1] + x[19+pa]*eta[a];
        dxdt[20 + 23*a] = x[19+23*a]*ga3 + x[18+23*a]*ga3*rho       - x[20+23*a]*om          - x[20+23*a]*eta[a+1] + x[20+pa]*eta[a];
        
        dxdt[21 + 23*a] = 0;
        
        dxdt[22 + 23*a] = si*(x[2 + 23*a] + x[7 + 23*a] + x[12 + 23*a] + x[17 + 23*a]);
    }
  }
};


double EvaluateLogLikelihood::evaluateLogLikelihoodCpp(VectorXd currentParamValues)
{
  // Restart values
  this->currentODETime = this->run_start;
  this->loglikelihoodError = false;
  this->dayNoAfterBurn = 0;
  this->weekNo = 0;
  this->valueLogLikelihood = 0;
  VectorXd inc_tot(this->A);
  inc_tot = VectorXd::Zero(this->A);
    
  ParameterValuesforODE(currentParamValues);
  
  // Set up and Run ODE solver
  EulerT<state_t> integrator;
  SystemT<state_t, system_t> System;
  asc::Recorder recorder;
  vector< double >  x0 = generateInitialStates();
  ODE_desc ODE_desc_inst(this);
  
  while (this->currentODETime < (this->run_full + this->run_burn)){
    integrator(ODE_desc_inst, x0, this->currentODETime, this->dt);
      
    checkStability(x0);
    if (this->loglikelihoodError){return log(0);}
    
    if (this->currentODETime > this->run_burn){
      getWeeklyLikelihood(x0, inc_tot);
      if (this->loglikelihoodError){return log(0);}
    }
  }
  // correction here
  return this->valueLogLikelihood; // RETURN THE LOG LIKELIHOOD
}

NumericMatrix EvaluateLogLikelihood::getWeeklySampleCpp(VectorXd currentParamValues, bool epFlag)
{
  // Restart values
  this->currentODETime = this->run_start;
  this->loglikelihoodError = false;
  this->dayNoAfterBurn = 0;
  this->weekNo = 0;
  this->valueLogLikelihood = 0;
  ParameterValuesforODE(currentParamValues);

  // Set up and Run ODE solver
  EulerT<state_t> integrator;
  SystemT<state_t, system_t> System;
  asc::Recorder recorder;
  vector< double >  x0 = generateInitialStates();
  ODE_desc ODE_desc_inst(this);
  
  NumericMatrix sampleWeeklyIncidence(52*7, this->A );
  
  while (this->currentODETime < (this->run_full + this->run_burn)){
    integrator(ODE_desc_inst, x0, this->currentODETime, this->dt);
    
    checkStability(x0);
    if (this->loglikelihoodError){return log(0);}
    
    if (this->currentODETime > this->run_burn){
      getWeeklyIncidence(x0, sampleWeeklyIncidence, epFlag);
      if (this->loglikelihoodError){return log(0);}
        
    }
  }
  return sampleWeeklyIncidence; // RETURN THE LOG LIKELIHOOD
}

NumericVector EvaluateLogLikelihood::getAnnualIncidenceCpp(VectorXd currentParamValues)
{
    // Restart values
    this->currentODETime = this->run_start;
    this->loglikelihoodError = false;
    this->dayNoAfterBurn = 0;
    this->weekNo = 0;
    this->valueLogLikelihood = 0;
    ParameterValuesforODE(currentParamValues);
    
    // Set up and Run ODE solver
    EulerT<state_t> integrator;
    SystemT<state_t, system_t> System;
    asc::Recorder recorder;
    vector< double >  x0 = generateInitialStates();
    ODE_desc ODE_desc_inst(this);
    
    NumericVector sampleAnnualIncidence(this->A);
    
    while (this->currentODETime < (this->run_full + this->run_burn)){
        integrator(ODE_desc_inst, x0, this->currentODETime, this->dt);
        
        checkStability(x0);
        if (this->loglikelihoodError){return log(0);}
        
        if ((this->currentODETime > this->run_burn) && (this->currentODETime < this->run_oneyr)){
            getAnnualIncidence(x0, sampleAnnualIncidence);
        }
    }
    for (int a = 0; a < 25; a++){
        sampleAnnualIncidence[a] = sampleAnnualIncidence[a]/populationPerAgeGroup[a];
    }
    return sampleAnnualIncidence; // RETURN THE LOG LIKELIHOOD
}
    
NumericVector EvaluateLogLikelihood::getProportionBornProtectedCpp(VectorXd currentParamValues)
{
    this->currentODETime = this->run_start;
    this->loglikelihoodError = false;
    this->dayNoAfterBurn = 0;
    this->weekNo = 0;
    this->valueLogLikelihood = 0;
    ParameterValuesforODE(currentParamValues);
    
    // Set up and Run ODE solver
    EulerT<state_t> integrator;
    SystemT<state_t, system_t> System;
    asc::Recorder recorder;
    vector< double >  x0 = generateInitialStates();
    ODE_desc ODE_desc_inst(this);
    
    NumericVector vectorProportionBornProtected(364);
    int i = 0;
    while (this->currentODETime < (this->run_full + this->run_burn))
    {
        integrator(ODE_desc_inst, x0, this->currentODETime, this->dt);
      
        if ((this->currentODETime > this->run_burn) && (this->currentODETime < (this->run_oneyr + 1))){
            vectorProportionBornProtected[i] = getProportionBornProtected(x0);
            i++;
        }
    }
    return vectorProportionBornProtected;
}



