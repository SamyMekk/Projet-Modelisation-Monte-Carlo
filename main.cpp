//
//  main.cpp
//  Pricing Black-Scholes Monte Carlo DIFFICILE ACTUEL C++
//
//  Created by GIOVANETTO on 19/01/2022.
//

#include <iostream>
#include <iterator>
#include <cmath>
#include <numeric>
#include <algorithm>    // Needed for the "max" function
#include <limits>
#include <string>
#include <vector>
#include <iomanip> // For rounding the numbers

double normalcdf(double x) // Phi(-∞, x) aka N(x)    // Funtion that computes the Cumulative Distribution function of a N(0,1)
{
    return std::erfc(-x/std::sqrt(2))/2;
}

double norm_pdf(const double x) {    // Funtion that computes the Probability Distribution function of a N(0,1)
  return (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x); // M_PI = PI
}

double gaussian_box_muller() {         // Algorithm generating random numbers from a N(0,1): Box-Muller Method
  double x = 0.0;
  double y = 0.0;
  double euclid_sq = 0.0;

  // Continue generating two uniform random variables
  // until the square of their "euclidean distance"
  // is less than unity
  do {
    x = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
    y = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
    euclid_sq = x*x + y*y;
  } while (euclid_sq >= 1.0);

  return x*sqrt(-2*log(euclid_sq)/euclid_sq);
}

std::vector<double> randn(int n) // Creating a vecotor of n numbers generated randomly from the previous algorithm
{
    std::vector<double> vector2(n, 0.0);
    double x = 0;
    for (int i=0;i<n;i++)
    {x = gaussian_box_muller();
        vector2[i]=x;}
return vector2; }

class BlackSholesModel{ // Class that stocks all the values common for all the options : it sets the environment in which we'll work : the Black-Scholes Model
    public:
std::vector<double> generateRiskNeutralPricePath(double toDate,int nSteps,double drift) const;
double stockPrice; // Prix initial
double volatility;
double riskFreeRate;
double drift;
double date;
void print()
    {
        std:: cout << stockPrice << ","<< volatility << ","<<riskFreeRate<<","<<date<< std:: endl;}
BlackSholesModel(double stockPrice, double volatility, double riskFreeRate, double date);// Constructor
BlackSholesModel(double stockPrice, double volatility, double riskFreeRate, double date,double drift); // Constructor
~BlackSholesModel();// Destructor
};

std::vector<double> BlackSholesModel::generateRiskNeutralPricePath(double toDate,int nSteps,double drift) const   // Generates a Risk Neutral Price Path for the Stock according to the Black-Scholes Model
    {
        std::vector<double> path(nSteps+1,0.0);
        path[0]=stockPrice;
        std::vector<double> epsilon = randn( nSteps );
        double dt = (toDate-date)/nSteps;
        double a = (drift - volatility*volatility*0.5)*dt;
        double b = volatility*sqrt(dt);
        double currentLogS = log( stockPrice );
        for (int i=0; i<nSteps; i++) {
            double dLogS = a + b*epsilon[i];
            double logS = currentLogS + dLogS;
            path[i+1] = exp( logS );
            currentLogS = logS;
        }
        return path;
    }
  

BlackSholesModel::BlackSholesModel(double stockPriceuser, double volatilityuser, double riskFreeRateuser, double dateuser) // Constructor
{
    stockPrice = stockPriceuser;
    volatility = volatilityuser;
    riskFreeRate = riskFreeRateuser;
    date = dateuser;
}

BlackSholesModel::BlackSholesModel(double stockPriceuser, double volatilityuser, double riskFreeRateuser, double dateuser, double driftuser)  // Constructor
{
    stockPrice = stockPriceuser;
    volatility = volatilityuser;
    riskFreeRate = riskFreeRateuser;
    date = dateuser;
    drift = driftuser;
}

BlackSholesModel::~BlackSholesModel()  // Destructor
{
}

class ContinuousTimeOption { // Class common for all the option types whether their pricing depends on the 'price path' of the stock (underlying) or not
    protected:
    double maturity;
    double strike;
    double barrier;
    public:
    virtual double priceclassic (const BlackSholesModel& bsm) const=0; // Pricing with explicit method if Possible (European Vanilla)
    // Const= 0 because this function will be used in Classes Inherited from this one, we have some polymorphism here --> this function is called a "Pure Virtual Function".
    virtual double priceMonteCarlo (const BlackSholesModel& bsm) const=0; // Pricing thanks to Monte-Carlo Method if it is not possible to use the Explicit Pricing
    virtual double payoff (std::vector<double>& stock) const=0; // Computation of the Payoff of the option
    virtual double computedeltaanalytic (const BlackSholesModel&bsm) const=0; // Computation of the Greeks with the Explicit Formulas
    virtual double computegammaanalytic (const BlackSholesModel&bsm) const=0; // Computation of the Greeks with the Explicit Formulas
    virtual double computethetaanalytic (const BlackSholesModel&bsm) const=0; // Computation of the Greeks with the Explicit Formulas
    virtual double computevegaanalytic (const BlackSholesModel&bsm) const=0; // Computation of the Greeks with the Explicit Formulas
    virtual double computerhoanalytic (const BlackSholesModel&bsm) const=0; // Computation of the Greeks with the Explicit Formulas
    virtual double computedeltaMonteCarlo (const BlackSholesModel&bsm) const=0; // Computation of the Greeks with the Finite Difference and Monte-Carlo Method
    virtual double computegammaMonteCarlo (const BlackSholesModel&bsm) const=0; // Computation of the Greeks with the Finite Difference and Monte-Carlo Method
    virtual double computethetaMonteCarlo (const BlackSholesModel&bsm) const=0; // Computation of the Greeks with the Finite Difference and Monte-Carlo Method
    virtual double computevegaMonteCarlo (const BlackSholesModel&bsm) const=0; // Computation of the Greeks with the Finite Difference and Monte-Carlo Method
    virtual double computerhoMonteCarlo (const BlackSholesModel&bsm) const=0; // Computation of the Greeks with the Finite Difference and Monte-Carlo Method
    virtual std::vector<double>  computeConfidenceInterval95 (double price,const BlackSholesModel&bsm) const=0; // Computation of the 95% Confidence Interval for the price of the Options
    virtual ~ContinuousTimeOption() {};
    virtual bool isPathDependent() const=0;
    void setMaturity(double maturityuser) { // Because this Variable is protected
        maturity = maturityuser;
        }
    int getMaturity() const // Because this Variable is protected
    {
            return maturity;
        }
    void setStrike(double strikeuser) { // Because this Variable is protected
        strike = strikeuser;
        }
    int getStrike() const // Because this Variable is protected
    {
            return strike;
        }
    void setBarrier(double barrieruser) { // Because this Variable is protected
        barrier = barrieruser;
        }
    int getBarrier() const // Because this Variable is protected
    {
            return barrier;
        }
};

class PathdependentOption : public ContinuousTimeOption { // Class common for all the option types that have their pricing that depends on the 'price path' of the stock (underlying)
    public:
    virtual ~PathdependentOption() {};
    bool isPathDependent() const {return true;};
    void setMaturity(double maturityuser) {
        maturity = maturityuser;
        }
    int getMaturity() const
    {
            return maturity;
        }
    void setStrike(double strikeuser) {
        strike = strikeuser;
        }
    int getStrike() const
    {
            return strike;
        }
    void setBarrier(double barrieruser) {
        barrier = barrieruser;
        }
    int getBarrier() const
    {
            return barrier;
        }
    virtual double payoff (std::vector<double>& stock) const=0;
    virtual double priceMonteCarlo(const BlackSholesModel&bsm)const=0;
};

class PathIndependentOption: public ContinuousTimeOption{ // Class common for all the option types that have their pricing that doesn't depend on the 'price path' of the stock (underlying)
    public:
    void setMaturity(double maturityuser) {
        maturity = maturityuser;
        }
    int getMaturity() const
    {
            return maturity;
        }
    void setStrike(double strikeuser) {
        strike = strikeuser;
        }
    int getStrike() const
    {
            return strike;
        }
    bool isPathDependent() const {return false;};
    virtual double payoff (std::vector<double>& stock) const=0;
    virtual double priceclassic(const BlackSholesModel&bsm)const=0; // Const pour indiquer que la fonction va être utilisée dans une classe dérivée, c'est une "Pure Virtual Function"
    virtual double priceMonteCarlo(const BlackSholesModel&bsm)const=0;
    virtual ~PathIndependentOption(){};
    
};

class MonteCarloPricer {  // Class related to the Monte-Carlo method
    public:
    /*  Constructor */
    MonteCarloPricer();
    /*  Number of scenarios */
    int nScenarios;
/* number of steps */
    int nSteps;
    double price(const ContinuousTimeOption & option,const BlackSholesModel& model) const; // Price of an option whether it is Path Dependent or Independent //
    double delta(const ContinuousTimeOption & option,const BlackSholesModel& model) const;
    double vega(const ContinuousTimeOption & option,const BlackSholesModel& model) const;
    double theta(const ContinuousTimeOption & option,const BlackSholesModel& model) const;
    double rho(const ContinuousTimeOption & option,const BlackSholesModel& model) const;
    double gamma(const ContinuousTimeOption & option,const BlackSholesModel& model) const;
    double standarderror (const ContinuousTimeOption& option,const BlackSholesModel& model )const;
    ~MonteCarloPricer();
};

MonteCarloPricer::MonteCarloPricer()
{
    nScenarios = 100000 ; // We chose 100000 simulations
    nSteps = 250; // We chose 250 Steps
}

double MonteCarloPricer::delta(const ContinuousTimeOption& option,const BlackSholesModel& model)const
    {
    double total = 0.0;
    double maturity = option.getMaturity();
    double r = model.riskFreeRate;
    double T = maturity - model.date;
    for (int i=0; i<nScenarios; i++) // We generate all our scenarios
    {
        std::vector<double> path1(nSteps+1,0.0);
        std::vector<double> path2(nSteps+1,0.0);
        path1[0]=model.stockPrice;
        path2[0]=model.stockPrice+0.001;
        std::vector<double> epsilon = randn( nSteps );
        double dt = (maturity-model.date)/nSteps;
        double a = (r - model.volatility*model.volatility*0.5)*dt;
        double b = model.volatility*sqrt(dt);
        double currentLogS1 = log(model.stockPrice);
        double currentLogS2 =log(model.stockPrice+0.001);
        for (int i=0; i<nSteps; i++) {
            double dLogS1 = a + b*epsilon[i];
            double logS1 = currentLogS1 + dLogS1;
            path1[i+1] = exp( logS1 );
            currentLogS1 = logS1;
            double dLogS2 = a + b*epsilon[i];
            double logS2 = currentLogS2 + dLogS2;
            path2[i+1] = exp( logS2 );
            currentLogS2 = logS2;
        }
        double payoff1 = option.payoff(path1); // We compute the payoff of the option according to this precise price path
        double payoff2 = option.payoff(path2); // We compute the payoff of the option according to this precise price path
        total+= (exp(-r*T)*(payoff2-payoff1))/(0.001); } // We sum all the payoffs are required in Monte Carlo
        double mean = total/nScenarios;
        return mean;
}

double MonteCarloPricer::gamma(const ContinuousTimeOption& option,const BlackSholesModel& model)const
    {
    double total = 0.0;
    double maturity = option.getMaturity();
    double r = model.riskFreeRate;
    double T = maturity - model.date;
    for (int i=0; i<nScenarios; i++) // We generate all our scenarios
    {
        std::vector<double> path1(nSteps+1,0.0);
        std::vector<double> path2(nSteps+1,0.0);
        std::vector<double> path3(nSteps+1,0.0);
        path1[0]=model.stockPrice;
        path2[0]=model.stockPrice+0.001;
        path3[0]=model.stockPrice-0.001;
        std::vector<double> epsilon = randn( nSteps );
        double dt = (maturity-model.date)/nSteps;
        double a = (r - model.volatility*model.volatility*0.5)*dt;
        double b = model.volatility*sqrt(dt);
        double currentLogS1 = log(model.stockPrice);
        double currentLogS2 =log(model.stockPrice+0.001);
        double currentLogS3 =log(model.stockPrice-0.001);
        for (int i=0; i<nSteps; i++) {
            double dLogS1 = a + b*epsilon[i];
            double logS1 = currentLogS1 + dLogS1;
            path1[i+1] = exp( logS1 );
            currentLogS1 = logS1;
            double dLogS2 = a + b*epsilon[i];
            double logS2 = currentLogS2 + dLogS2;
            path2[i+1] = exp( logS2 );
            currentLogS2 = logS2;
            double dLogS3 = a + b*epsilon[i];
            double logS3 = currentLogS3 + dLogS3;
            path3[i+1] = exp( logS3 );
            currentLogS3 = logS3;
        }
        double payoff1 = option.payoff(path1); // We compute the payoff of the option according to this precise price path
        double payoff2 = option.payoff(path2); // We compute the payoff of the option according to this precise price path
        double payoff3 = option.payoff(path3); // We compute the payoff of the option according to this precise price path
        total+= (exp(-r*T)*(payoff2 - 2*payoff1 + payoff3))/(0.001*0.001); } // We sum all the payoffs are required in Monte Carlo
        double mean = total/nScenarios;
        return mean;
}

double MonteCarloPricer::vega(const ContinuousTimeOption& option,const BlackSholesModel& model)const
    {
    double total = 0.0;
    double maturity = option.getMaturity();
    double r = model.riskFreeRate;
    double T = maturity - model.date;
    for (int i=0; i<nScenarios; i++) // We generate all our scenarios
    {  std::vector<double> path1(nSteps+1,0.0);
        std::vector<double> path2(nSteps+1,0.0);
        path1[0]=model.stockPrice;
        path2[0]=model.stockPrice;
        std::vector<double> epsilon = randn( nSteps );
        double dt = (maturity-model.date)/nSteps;
        double a1 = (r - model.volatility*model.volatility*0.5)*dt;
        double b1 = model.volatility*sqrt(dt);
        double a2 = (r -((model.volatility)+0.001)*((model.volatility)+0.001)*0.5)*dt;
        double b2 = ((model.volatility)+0.001)*sqrt(dt);
        double currentLogS1 = log(model.stockPrice);
        double currentLogS2 =log(model.stockPrice);
        for (int i=0; i<nSteps; i++) {
            double dLogS1 = a1 + b1*epsilon[i];
            double logS1 = currentLogS1 + dLogS1;
            path1[i+1] = exp( logS1 );
            currentLogS1 = logS1;
            double dLogS2 = a2 + b2*epsilon[i];
            double logS2 = currentLogS2 + dLogS2;
            path2[i+1] = exp( logS2 );
            currentLogS2 = logS2;
        }
    double payoff1 = option.payoff(path1); // We compute the payoff of the option according to this precise price path
    double payoff2 = option.payoff(path2); // We compute the payoff of the option according to this precise price path
    total+= (exp(-r*T)*(payoff2 - payoff1))/0.001; } // We sum all the payoffs are required in Monte Carlo
        double mean = total/nScenarios;
        return mean;
}

double MonteCarloPricer::theta(const ContinuousTimeOption& option,const BlackSholesModel& model)const
    {
    double total = 0.0;
    double maturity = option.getMaturity();
    double r = model.riskFreeRate;
    double T = maturity - model.date;
    for (int i=0; i<nScenarios; i++) // We generate all our scenarios
    {   std::vector<double> path1(nSteps+1,0.0);
        std::vector<double> path2(nSteps+1,0.0);
        path1[0]=model.stockPrice;
        path2[0]=model.stockPrice;
        std::vector<double> epsilon = randn( nSteps );
        double dt1 = (maturity-model.date)/nSteps;
        double dt2 = (maturity-model.date+0.001)/nSteps;
        double a1 = (r - (model.volatility)*(model.volatility)*0.5)*dt1;
        double b1 = model.volatility*sqrt(dt1);
        double a2 = (r - (model.volatility)*(model.volatility)*0.5)*dt2;
        double b2 = (model.volatility)*sqrt(dt2);
        double currentLogS1 = log(model.stockPrice);
        double currentLogS2 =log(model.stockPrice);
        for (int i=0; i<nSteps; i++) {
            double dLogS1 = a1 + b1*epsilon[i];
            double logS1 = currentLogS1 + dLogS1;
            path1[i+1] = exp( logS1 );
            currentLogS1 = logS1;
            double dLogS2 = a2 + b2*epsilon[i];
            double logS2 = currentLogS2 + dLogS2;
            path2[i+1] = exp( logS2 );
            currentLogS2 = logS2;
        }
            double payoff1 = option.payoff(path1); // We compute the payoff of the option according to this precise price path
            double payoff2 = option.payoff(path2); // We compute the payoff of the option according to this precise price path
            total+= -((exp(-r*(T+0.001))*payoff2) - (exp(-r*T)*payoff1))/0.001;
        
    } // We sum all the payoffs are required in Monte Carlo
        double mean = total/nScenarios;
        return mean;
    }

double MonteCarloPricer::rho(const ContinuousTimeOption& option,const BlackSholesModel& model)const
    {
    double total = 0.0;
    double maturity = option.getMaturity();
    double r = model.riskFreeRate;
    double T = maturity - model.date;
    for (int i=0; i<nScenarios; i++) // We generate all our scenarios
    {   std::vector<double> path1(nSteps+1,0.0);
        std::vector<double> path2(nSteps+1,0.0);
        path1[0]=model.stockPrice;
        path2[0]=model.stockPrice;
        std::vector<double> epsilon = randn( nSteps );
        double dt = (maturity-model.date)/nSteps;
        double a1 = (r - model.volatility*model.volatility*0.5)*dt;
        double b1 = model.volatility*sqrt(dt);
        double a2 = ((r+0.001) - model.volatility*model.volatility*0.5)*dt;
        double b2 = model.volatility*sqrt(dt);
        double currentLogS1 = log(model.stockPrice);
        double currentLogS2 =log(model.stockPrice);
        for (int i=0; i<nSteps; i++) {
            double dLogS1 = a1 + b1*epsilon[i];
            double logS1 = currentLogS1 + dLogS1;
            path1[i+1] = exp( logS1 );
            currentLogS1 = logS1;
            double dLogS2 = a2 + b2*epsilon[i];
            double logS2 = currentLogS2 + dLogS2;
            path2[i+1] = exp( logS2 );
            currentLogS2 = logS2;
        }
            double payoff1 = option.payoff(path1); // We compute the payoff of the option according to this precise price path
            double payoff2 = option.payoff(path2); // We compute the payoff of the option according to this precise price path
        total+= ((exp(-(r+0.001)*T)*payoff2 - exp(-r*T)*payoff1))/0.001;} // We sum all the payoffs are required in Monte Carlo
        double mean = total/nScenarios;
        return mean;
    }

double MonteCarloPricer::price(const ContinuousTimeOption& option,const BlackSholesModel& model )const
    {
    double total = 0.0;
    double maturity = option.getMaturity();
    for (int i=0; i<nScenarios; i++) // We generate all our scenarios
    {
std::vector<double> path=model.generateRiskNeutralPricePath(maturity,nSteps,model.drift); // We have the price Path here
            double payoff = option.payoff(path); // We compute the payoff of the option according to this precise price path
            total+= payoff;} // We sum all the payoffs are required in Monte Carlo
        double mean = total/nScenarios;
        double r = model.riskFreeRate;
        double T = maturity - model.date;
        return exp(-r*T)*mean;
}

double MonteCarloPricer::standarderror(const ContinuousTimeOption& option,const BlackSholesModel& model )const
    {
    double total = 0.0;
        double totalsquare = 0.0;
    double maturity = option.getMaturity();
    for (int i=0; i<nScenarios; i++) // We generate all our scenarios
    {
std::vector<double> path=model.generateRiskNeutralPricePath(maturity,nSteps,model.drift); // We have avec price Path here
            double payoff = option.payoff(path); // We compute the payoff of the option according to this precise price path
            total+= payoff;
        totalsquare += payoff*payoff;
    }
        double mean = total/nScenarios;
        double r = model.riskFreeRate;
        double T = maturity - model.date;
        double part1 = exp(-r*T)*mean;
        double meansquared = totalsquare/nScenarios;
        double part2 = exp(-2*r*T)*meansquared;
        return (nScenarios/(nScenarios-1))*(part2 - (part1*part1));
}

MonteCarloPricer::~MonteCarloPricer()
{
}

class EuropeanVanillaCallOption: public PathIndependentOption { // Type : European Vanilla Call --> Path Independent
    public:
    void setMaturity(double maturityuser) {
        maturity = maturityuser;
        }
    int getMaturity() const
    {
            return maturity;
        }
    void setStrike(double strikeuser) {
        strike = strikeuser;
        }
    int getStrike() const
    {
            return strike;
        }
    double payoff (std::vector<double>& stock) const;
    bool isPathDependent() const {return false;};
    double priceclassic (const BlackSholesModel& bsm) const;
    double  priceMonteCarlo (const BlackSholesModel& bsm) const;
    double computedeltaanalytic (const BlackSholesModel&bsm) const;
    double computegammaanalytic (const BlackSholesModel&bsm) const;
    double computethetaanalytic (const BlackSholesModel&bsm) const;
    double computevegaanalytic (const BlackSholesModel&bsm) const;
    double computerhoanalytic (const BlackSholesModel&bsm) const;
    double computedeltaMonteCarlo (const BlackSholesModel&bsm) const;
    double computegammaMonteCarlo (const BlackSholesModel&bsm) const;
    double computethetaMonteCarlo (const BlackSholesModel&bsm) const;
    double computevegaMonteCarlo (const BlackSholesModel&bsm) const;
    double computerhoMonteCarlo (const BlackSholesModel&bsm) const;
    std::vector<double> computeConfidenceInterval95 (double price,const BlackSholesModel&bsm) const;
    ~EuropeanVanillaCallOption(); };

double EuropeanVanillaCallOption :: payoff (std::vector<double>& stock) const {
    double SMaturity=stock.back();
    double K=strike;
    return fmax(SMaturity-K,0); } // Payoff of an European Vanilla Call

double EuropeanVanillaCallOption :: computedeltaanalytic (const BlackSholesModel& bsm) const
 {
double S=bsm.stockPrice;
double K=strike;
double sigma=bsm.volatility;
double r=bsm.riskFreeRate;
double T=maturity-bsm.date;
double numerator=log(S/K)+(r+sigma*sigma*0.5)*T;
double denominator=sigma*sqrt(T);
double d1=numerator/denominator;
return normalcdf(d1); }

double EuropeanVanillaCallOption :: computegammaanalytic (const BlackSholesModel& bsm) const
 {
double S=bsm.stockPrice;
double K=strike;
double sigma=bsm.volatility;
double r=bsm.riskFreeRate;
double T=maturity-bsm.date;
double numerator=log(S/K)+(r+sigma*sigma*0.5)*T;
double denominator=sigma*sqrt(T);
double d1=numerator/denominator;
double numeratortotal = norm_pdf(d1);
double denominatortotal = S*sigma*sqrt(T);
return numeratortotal/denominatortotal;
 }

double EuropeanVanillaCallOption :: computethetaanalytic (const BlackSholesModel& bsm) const
 {
double S=bsm.stockPrice;
double K=strike;
double sigma=bsm.volatility;
double r=bsm.riskFreeRate;
double T=maturity-bsm.date;
double numerator=log(S/K)+(r+sigma*sigma*0.5)*T;
double denominator=sigma*sqrt(T);
double d1=numerator/denominator;
double d2=d1-denominator;
double numeratortotal = S*norm_pdf(d1)*sigma;
double denominatortotal = 2*sqrt(T);
double part1 = numeratortotal/denominatortotal;
double part2 = r*K*exp(-(r)*T)*normalcdf(d2);
return - part1 - part2;
 }

double EuropeanVanillaCallOption :: computevegaanalytic (const BlackSholesModel& bsm) const
 {
double S=bsm.stockPrice;
double K=strike;
double sigma=bsm.volatility;
double r=bsm.riskFreeRate;
double T=maturity-bsm.date;
double numerator=log(S/K)+(r+sigma*sigma*0.5)*T;
double denominator=sigma*sqrt(T);
double d1=numerator/denominator;
return S*norm_pdf(d1)*sqrt(T);
 }

double EuropeanVanillaCallOption :: computerhoanalytic (const BlackSholesModel& bsm) const
 {
double S=bsm.stockPrice;
double K=strike;
double sigma=bsm.volatility;
double r=bsm.riskFreeRate;
double T=maturity-bsm.date;
double numerator=log(S/K)+(r+sigma*sigma*0.5)*T;
double denominator=sigma*sqrt(T);
double d1=numerator/denominator;
double d2=d1-denominator;
return K*T*exp(-(r)*T)*normalcdf(d2);
 }

double EuropeanVanillaCallOption :: priceclassic (const BlackSholesModel& bsm) const
 {
double S=bsm.stockPrice;
double K=strike;
double sigma=bsm.volatility;
double r=bsm.riskFreeRate;
double T=maturity-bsm.date;
double numerator=log(S/K)+(r+sigma*sigma*0.5)*T;
double denominator=sigma*sqrt(T);
double d1=numerator/denominator;
double d2=d1-denominator;
return S*normalcdf(d1)-exp(-r*T)*K*normalcdf(d2); }

double EuropeanVanillaCallOption :: computedeltaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.delta(*this, bsm); // On a ici l'utilisation de *this --> c'est un pointeur vers l'objet EuropeanVanillaCallOption sur lequel nous travaillons.//
 }

double EuropeanVanillaCallOption :: computegammaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.gamma(*this, bsm);
 }


double EuropeanVanillaCallOption :: computevegaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.vega(*this, bsm);
 }

double EuropeanVanillaCallOption :: computethetaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.theta(*this, bsm);
 }

double EuropeanVanillaCallOption :: computerhoMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.rho(*this, bsm);
 }

double EuropeanVanillaCallOption::priceMonteCarlo(const BlackSholesModel& bsm) const
{
    MonteCarloPricer pricer;
    return pricer.price( *this, bsm );
}

std::vector<double>  EuropeanVanillaCallOption :: computeConfidenceInterval95(double price, const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     double variance = pricer.standarderror(*this, bsm);
     std::vector<double> confidenceinderval(2);
     double lb = price - (1.960*variance/sqrt(100000));
     double ub = price + (1.960*variance/sqrt(100000));
     confidenceinderval[0] = lb;
     confidenceinderval[1] = ub;
     return confidenceinderval;
 }

EuropeanVanillaCallOption::~EuropeanVanillaCallOption()
{
}

class EuropeanVanillaPutOption: public PathIndependentOption {  // Type : European Vanilla Put --> Path Independent
    public:
    void setMaturity(double maturityuser) {
        maturity = maturityuser;
        }
    int getMaturity() const
    {
            return maturity;
        }
    void setStrike(double strikeuser) {
        strike = strikeuser;
        }
    int getStrike() const
    {
            return strike;
        }
    double payoff (std::vector<double>& stock) const;
    bool isPathDependent() const {return false;};
    double computedeltaanalytic (const BlackSholesModel&bsm) const;
    double computegammaanalytic (const BlackSholesModel&bsm) const;
    double computevegaanalytic (const BlackSholesModel&bsm) const;
    double computethetaanalytic (const BlackSholesModel&bsm) const;
    double computerhoanalytic (const BlackSholesModel&bsm) const;
    double computedeltaMonteCarlo (const BlackSholesModel&bsm) const;
    double computegammaMonteCarlo (const BlackSholesModel&bsm) const;
    double computethetaMonteCarlo (const BlackSholesModel&bsm) const;
    double computevegaMonteCarlo (const BlackSholesModel&bsm) const;
    double computerhoMonteCarlo (const BlackSholesModel&bsm) const;
    std::vector<double> computeConfidenceInterval95 (double price,const BlackSholesModel&bsm) const;
    double priceclassic (const BlackSholesModel& bsm) const;
    double  priceMonteCarlo (const BlackSholesModel& bsm) const;
    // We can use Both Monte Carlo and the Explicit Method here
    ~EuropeanVanillaPutOption(); };

double EuropeanVanillaPutOption:: payoff (std::vector<double>& stock) const
 {
double SMaturity=stock.back();
double K=strike;
return fmax(K-SMaturity,0); }

double EuropeanVanillaPutOption :: computedeltaanalytic (const BlackSholesModel& bsm) const
 {
double S=bsm.stockPrice;
double K=strike;
double sigma=bsm.volatility;
double r=bsm.riskFreeRate;
double T=maturity-bsm.date;
double numerator=log(S/K)+(r+sigma*sigma*0.5)*T;
double denominator=sigma*sqrt(T);
double d1=numerator/denominator;
return normalcdf(d1)-1; }

double EuropeanVanillaPutOption :: computegammaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double EuropeanVanillaPutOption :: computevegaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double EuropeanVanillaPutOption :: computethetaanalytic (const BlackSholesModel& bsm) const
 {
double S=bsm.stockPrice;
double K=strike;
double sigma=bsm.volatility;
double r=bsm.riskFreeRate;
double T=maturity-bsm.date;
double numerator=log(S/K)+(r+sigma*sigma*0.5)*T;
double denominator=sigma*sqrt(T);
double d1=numerator/denominator;
double d2=d1-denominator;
double numeratortotal = S*norm_pdf(d1)*sigma;
double denominatortotal = 2*sqrt(T);
double part1 = numeratortotal/denominatortotal;
double part2 = r*K*exp(-(r)*T)*normalcdf(-d2);
return - part1 + part2;
 }

double EuropeanVanillaPutOption :: computedeltaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.delta(*this, bsm);
 }

double EuropeanVanillaPutOption :: computegammaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.gamma(*this, bsm);
 }

double EuropeanVanillaPutOption :: computevegaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.vega(*this, bsm);
 }

double EuropeanVanillaPutOption :: computethetaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.theta(*this, bsm);
 }

double EuropeanVanillaPutOption :: computerhoMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.rho(*this, bsm);
 }

double EuropeanVanillaPutOption :: computerhoanalytic (const BlackSholesModel& bsm) const
 {
double S=bsm.stockPrice;
double K=strike;
double sigma=bsm.volatility;
double r=bsm.riskFreeRate;
double T=maturity-bsm.date;
double numerator=log(S/K)+(r+sigma*sigma*0.5)*T;
double denominator=sigma*sqrt(T);
double d1=numerator/denominator;
double d2=d1-denominator;
return -K*T*exp(-(r)*T)*normalcdf(-d2);
 }

double EuropeanVanillaPutOption:: priceclassic (const BlackSholesModel& bsm) const //We use the Call-Put Parity here, as seen in the Financial Instruments Course
 {
double S=bsm.stockPrice;
double K=strike;
double sigma=bsm.volatility;
double r=bsm.riskFreeRate;
double T=maturity-bsm.date;
double numerator=log(S/K)+(r+sigma*sigma*0.5)*T;
double denominator=sigma*sqrt(T);
double d1=numerator/denominator;
double d2=d1-denominator;
return K*exp(-r*T) * normalcdf(-d2)-S*normalcdf(-d1); }

double EuropeanVanillaPutOption::priceMonteCarlo(const BlackSholesModel& bsm) const
{
    MonteCarloPricer pricer;
    return pricer.price(*this,bsm );
}

std::vector<double>  EuropeanVanillaPutOption :: computeConfidenceInterval95(double price, const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     double variance = pricer.standarderror(*this, bsm);
     std::vector<double> confidenceinderval(2);
     double lb = price - (1.960*variance/sqrt(100000));
     double ub = price + (1.960*variance/sqrt(100000));
     confidenceinderval[0] = lb;
     confidenceinderval[1] = ub;
     return confidenceinderval;
 }

EuropeanVanillaPutOption::~EuropeanVanillaPutOption()
{

}

class AsianArithmeticCallOption: public PathdependentOption { // Type : Asian Arithmetic Call --> Path Dependent
    public:
    void setMaturity(double maturityuser) {
        maturity = maturityuser;
        }
    int getMaturity() const
    {
            return maturity;
        }
    void setStrike(double strikeuser) {
        strike = strikeuser;
        }
    int getStrike() const
    {
            return strike;
        }
    double payoff (std::vector<double>& stock) const;
    double computedeltaanalytic (const BlackSholesModel&bsm) const;
    double computegammaanalytic (const BlackSholesModel&bsm) const;
    double computethetaanalytic (const BlackSholesModel&bsm) const;
    double computevegaanalytic (const BlackSholesModel&bsm) const;
    double computerhoanalytic (const BlackSholesModel&bsm) const;
    double computedeltaMonteCarlo (const BlackSholesModel&bsm) const;
    double computegammaMonteCarlo (const BlackSholesModel&bsm) const;
    double computethetaMonteCarlo (const BlackSholesModel&bsm) const;
    double computevegaMonteCarlo (const BlackSholesModel&bsm) const;
    double computerhoMonteCarlo (const BlackSholesModel&bsm) const;
    std::vector<double> computeConfidenceInterval95 (double price,const BlackSholesModel&bsm) const;
    bool isPathDependent() const {return true;};
    double priceclassic (const BlackSholesModel& bsm) const;
    double  priceMonteCarlo (const BlackSholesModel& bsm) const;
    // We have to mention both of the methods here in order for this class not to be an abstract class but we can only use the Monte Carlo Method here
    ~AsianArithmeticCallOption(); };

double AsianArithmeticCallOption :: payoff (std::vector<double>& stock) const { // Computation of the Payoff according to this precise type of option
    double SMaturity= std::accumulate(stock.begin(), stock.end(), 0.0);
    double SMaturityMean = SMaturity/stock.size();
    double K=strike;
    return fmax(SMaturityMean-K,0); }

double AsianArithmeticCallOption :: computedeltaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double AsianArithmeticCallOption :: computegammaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double AsianArithmeticCallOption :: computethetaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double AsianArithmeticCallOption :: computevegaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double AsianArithmeticCallOption :: computerhoanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double AsianArithmeticCallOption :: computedeltaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.delta(*this, bsm);
 }

double AsianArithmeticCallOption :: computegammaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.gamma(*this, bsm);
 }

double AsianArithmeticCallOption :: computevegaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.vega(*this, bsm);
 }

double AsianArithmeticCallOption :: computethetaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.theta(*this, bsm);
 }

double AsianArithmeticCallOption :: computerhoMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.rho(*this, bsm);
 }

double AsianArithmeticCallOption::priceMonteCarlo(const BlackSholesModel& bsm) const
{
    MonteCarloPricer pricer;
    return pricer.price( *this, bsm );
}

double AsianArithmeticCallOption::priceclassic(const BlackSholesModel& bsm) const
{
    return 0; // We can't use this method here, therefore we chose to return always 0 --> we'll do our computations only by using Monte Carlo
}

std::vector<double>  AsianArithmeticCallOption :: computeConfidenceInterval95(double price, const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     double variance = pricer.standarderror(*this, bsm);
     std::vector<double> confidenceinderval(2);
     double lb = price - (1.960*variance/sqrt(100000));
     double ub = price + (1.960*variance/sqrt(100000));
     confidenceinderval[0] = lb;
     confidenceinderval[1] = ub;
     return confidenceinderval;
 }

AsianArithmeticCallOption::~AsianArithmeticCallOption()
{
}

class AsianArithmeticPutOption: public PathdependentOption { // Type : Asian Arithmetic Put --> Path Dependent
    public:
    void setMaturity(double maturityuser) {
        maturity = maturityuser;
        }
    int getMaturity() const
    {
            return maturity;
        }
    void setStrike(double strikeuser) {
        strike = strikeuser;
        }
    int getStrike() const
    {
            return strike;
        }
    double payoff (std::vector<double>& stock) const;
    double computedeltaanalytic (const BlackSholesModel&bsm) const;
    double computegammaanalytic (const BlackSholesModel&bsm) const;
    double computethetaanalytic (const BlackSholesModel&bsm) const;
    double computevegaanalytic (const BlackSholesModel&bsm) const;
    double computerhoanalytic (const BlackSholesModel&bsm) const;
    double computedeltaMonteCarlo (const BlackSholesModel&bsm) const;
    double computegammaMonteCarlo (const BlackSholesModel&bsm) const;
    double computethetaMonteCarlo (const BlackSholesModel&bsm) const;
    double computevegaMonteCarlo (const BlackSholesModel&bsm) const;
    double computerhoMonteCarlo (const BlackSholesModel&bsm) const;
    std::vector<double> computeConfidenceInterval95 (double price,const BlackSholesModel&bsm) const;
    bool isPathDependent() const {return true;};
    double priceclassic (const BlackSholesModel& bsm) const;
    double  priceMonteCarlo (const BlackSholesModel& bsm) const;
    ~AsianArithmeticPutOption(); };

double AsianArithmeticPutOption :: payoff (std::vector<double>& stock) const {
    // Computation of the Payoff according to this precise type of option
    double SMaturity= std::accumulate(stock.begin(), stock.end(), 0.0);
    double SMaturityMean = SMaturity/stock.size();
    double K=strike;
    return fmax(K-SMaturityMean,0); }

double AsianArithmeticPutOption :: computedeltaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double AsianArithmeticPutOption :: computegammaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double AsianArithmeticPutOption :: computethetaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double AsianArithmeticPutOption :: computevegaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double AsianArithmeticPutOption :: computerhoanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double AsianArithmeticPutOption :: computedeltaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.delta(*this, bsm);
 }

double AsianArithmeticPutOption :: computegammaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.gamma(*this, bsm);
 }

double AsianArithmeticPutOption :: computevegaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.vega(*this, bsm);
 }

double AsianArithmeticPutOption :: computethetaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.theta(*this, bsm);
 }

double AsianArithmeticPutOption :: computerhoMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.rho(*this, bsm);
 }

double AsianArithmeticPutOption::priceMonteCarlo(const BlackSholesModel& bsm) const
{
    MonteCarloPricer pricer;
    return pricer.price( *this, bsm );
}

double AsianArithmeticPutOption::priceclassic(const BlackSholesModel& bsm) const
{
    return 0;
}

std::vector<double>  AsianArithmeticPutOption :: computeConfidenceInterval95(double price, const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     double variance = pricer.standarderror(*this, bsm);
     std::vector<double> confidenceinderval(2);
     double lb = price - (1.960*variance/sqrt(100000));
     double ub = price + (1.960*variance/sqrt(100000));
     confidenceinderval[0] = lb;
     confidenceinderval[1] = ub;
     return confidenceinderval;
 }

AsianArithmeticPutOption::~AsianArithmeticPutOption()
{
}

class AsianGeometricCallOption: public PathdependentOption { // Type : Asian Geometric Call --> Path Dependent
    public:
    void setMaturity(double maturityuser) {
        maturity = maturityuser;
        }
    int getMaturity() const
    {
            return maturity;
        }
    void setStrike(double strikeuser) {
        strike = strikeuser;
        }
    int getStrike() const
    {
            return strike;
        }
    double payoff (std::vector<double>& stock) const;
    double computedeltaanalytic (const BlackSholesModel&bsm) const;
    double computegammaanalytic (const BlackSholesModel&bsm) const;
    double computethetaanalytic (const BlackSholesModel&bsm) const;
    double computevegaanalytic (const BlackSholesModel&bsm) const;
    double computerhoanalytic (const BlackSholesModel&bsm) const;
    double computedeltaMonteCarlo (const BlackSholesModel&bsm) const;
    double computegammaMonteCarlo (const BlackSholesModel&bsm) const;
    double computethetaMonteCarlo (const BlackSholesModel&bsm) const;
    double computevegaMonteCarlo (const BlackSholesModel&bsm) const;
    double computerhoMonteCarlo (const BlackSholesModel&bsm) const;
    std::vector<double> computeConfidenceInterval95 (double price,const BlackSholesModel&bsm) const;
    bool isPathDependent() const {return true;};
    double priceclassic (const BlackSholesModel& bsm) const;
    double  priceMonteCarlo (const BlackSholesModel& bsm) const;
    ~AsianGeometricCallOption(); };

double AsianGeometricCallOption :: payoff (std::vector<double>& stock) const {
    // Computation of the Payoff according to this precise type of option
    double log_sum= 0;
    for (int i=0; i<stock.size(); i++) {
      log_sum += log(stock[i]);
    }
    double geom_mean = exp(log_sum /stock.size());
    double K=strike;
    return fmax(geom_mean-K,0); }

double AsianGeometricCallOption :: computedeltaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double AsianGeometricCallOption :: computegammaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double AsianGeometricCallOption :: computethetaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double AsianGeometricCallOption :: computevegaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double AsianGeometricCallOption :: computerhoanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double AsianGeometricCallOption :: computedeltaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.delta(*this, bsm);
 }

double AsianGeometricCallOption :: computegammaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.gamma(*this, bsm);
 }

double AsianGeometricCallOption :: computevegaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.vega(*this, bsm);
 }

double AsianGeometricCallOption :: computethetaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.theta(*this, bsm);
 }

double AsianGeometricCallOption :: computerhoMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.rho(*this, bsm);
 }

double AsianGeometricCallOption::priceMonteCarlo(const BlackSholesModel& bsm) const
{
    MonteCarloPricer pricer;
    return pricer.price( *this, bsm );
}

double AsianGeometricCallOption::priceclassic(const BlackSholesModel& bsm) const
{
    return 0;
}

std::vector<double>  AsianGeometricCallOption :: computeConfidenceInterval95(double price, const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     double variance = pricer.standarderror(*this, bsm);
     std::vector<double> confidenceinderval(2);
     double lb = price - (1.960*variance/sqrt(100000));
     double ub = price + (1.960*variance/sqrt(100000));
     confidenceinderval[0] = lb;
     confidenceinderval[1] = ub;
     return confidenceinderval;
 }

AsianGeometricCallOption::~AsianGeometricCallOption()
{
}

class AsianGeometricPutOption: public PathdependentOption { // Type : Asian Geometric Put --> Path Dependent
    public:
    void setMaturity(double maturityuser) {
        maturity = maturityuser;
        }
    int getMaturity() const
    {
            return maturity;
        }
    void setStrike(double strikeuser) {
        strike = strikeuser;
        }
    int getStrike() const
    {
            return strike;
        }
    double payoff (std::vector<double>& stock) const;
    double computedeltaanalytic (const BlackSholesModel&bsm) const;
    double computegammaanalytic (const BlackSholesModel&bsm) const;
    double computethetaanalytic (const BlackSholesModel&bsm) const;
    double computevegaanalytic (const BlackSholesModel&bsm) const;
    double computerhoanalytic (const BlackSholesModel&bsm) const;
    double computedeltaMonteCarlo (const BlackSholesModel&bsm) const;
    double computegammaMonteCarlo (const BlackSholesModel&bsm) const;
    double computethetaMonteCarlo (const BlackSholesModel&bsm) const;
    double computevegaMonteCarlo (const BlackSholesModel&bsm) const;
    double computerhoMonteCarlo (const BlackSholesModel&bsm) const;
    std::vector<double> computeConfidenceInterval95 (double price,const BlackSholesModel&bsm) const;
    bool isPathDependent() const {return true;};
    double priceclassic (const BlackSholesModel& bsm) const;
    double  priceMonteCarlo (const BlackSholesModel& bsm) const;
    ~AsianGeometricPutOption(); };

double AsianGeometricPutOption :: payoff (std::vector<double>& stock) const {
    double log_sum= 0;
    for (int i=0; i<stock.size(); i++) {
      log_sum += log(stock[i]);
    }
    double geom_mean = exp(log_sum /stock.size());
    double K=strike;
    return fmax(K-geom_mean,0); }

double AsianGeometricPutOption :: computedeltaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double AsianGeometricPutOption :: computegammaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double AsianGeometricPutOption :: computethetaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double AsianGeometricPutOption :: computevegaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double AsianGeometricPutOption :: computerhoanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double AsianGeometricPutOption :: computedeltaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.delta(*this, bsm);
 }

double AsianGeometricPutOption :: computegammaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.gamma(*this, bsm);
 }

double AsianGeometricPutOption :: computevegaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.vega(*this, bsm);
 }

double AsianGeometricPutOption :: computethetaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.theta(*this, bsm);
 }

double AsianGeometricPutOption :: computerhoMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.rho(*this, bsm);
 }

double AsianGeometricPutOption::priceMonteCarlo(const BlackSholesModel& bsm) const
{
    MonteCarloPricer pricer;
    return pricer.price( *this, bsm );
}

double AsianGeometricPutOption::priceclassic(const BlackSholesModel& bsm) const
{
    return 0;
}

std::vector<double>  AsianGeometricPutOption :: computeConfidenceInterval95(double price, const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     double variance = pricer.standarderror(*this, bsm);
     std::vector<double> confidenceinderval(2);
     double lb = price - (1.960*variance/sqrt(100000));
     double ub = price + (1.960*variance/sqrt(100000));
     confidenceinderval[0] = lb;
     confidenceinderval[1] = ub;
     return confidenceinderval;
 }

AsianGeometricPutOption::~AsianGeometricPutOption()
{
}


class upandoutknockoutcalloption: public PathdependentOption { // Type : Up-and-Out Knock-out Call --> Path Dependent
    public:
    void setMaturity(double maturityuser) {
        maturity = maturityuser;
        }
    int getMaturity() const
    {
            return maturity;
        }
    void setStrike(double strikeuser) {
        strike = strikeuser;
        }
    int getStrike() const
    {
            return strike;
        }
    void setBarrier(double barrieruser) {
        barrier = barrieruser;
        }
    int getBarrier() const
    {
            return barrier;
        }
    double payoff (std::vector<double>& stock) const;
    double computedeltaanalytic (const BlackSholesModel&bsm) const;
    double computegammaanalytic (const BlackSholesModel&bsm) const;
    double computethetaanalytic (const BlackSholesModel&bsm) const;
    double computevegaanalytic (const BlackSholesModel&bsm) const;
    double computerhoanalytic (const BlackSholesModel&bsm) const;
    double computedeltaMonteCarlo (const BlackSholesModel&bsm) const;
    double computegammaMonteCarlo (const BlackSholesModel&bsm) const;
    double computethetaMonteCarlo (const BlackSholesModel&bsm) const;
    double computevegaMonteCarlo (const BlackSholesModel&bsm) const;
    double computerhoMonteCarlo (const BlackSholesModel&bsm) const;
    std::vector<double> computeConfidenceInterval95 (double price,const BlackSholesModel&bsm) const;
    bool isPathDependent() const {return true;};
    double priceclassic (const BlackSholesModel& bsm) const;
    double  priceMonteCarlo (const BlackSholesModel& bsm) const;
    ~upandoutknockoutcalloption(); };

double upandoutknockoutcalloption :: payoff (std::vector<double>& stock) const {
    double barrierofcall = barrier;
    for (int i=0; i<stock.size(); i++) {
        if (stock[i]>barrierofcall){
            return 0;}
        }
    double SMaturity= stock.back();
    double K=strike;
    return fmax(SMaturity-K,0); }

double upandoutknockoutcalloption :: computedeltaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double upandoutknockoutcalloption :: computegammaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double upandoutknockoutcalloption :: computethetaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double upandoutknockoutcalloption :: computevegaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double upandoutknockoutcalloption :: computerhoanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double upandoutknockoutcalloption :: computedeltaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.delta(*this, bsm);
 }

double upandoutknockoutcalloption :: computegammaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.gamma(*this, bsm);
 }

double upandoutknockoutcalloption :: computevegaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.vega(*this, bsm);
 }

double upandoutknockoutcalloption :: computethetaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.theta(*this, bsm);
 }

double upandoutknockoutcalloption :: computerhoMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.rho(*this, bsm);
 }


double upandoutknockoutcalloption::priceMonteCarlo(const BlackSholesModel& bsm) const
{
    MonteCarloPricer pricer;
    return pricer.price( *this, bsm );
}

double upandoutknockoutcalloption::priceclassic(const BlackSholesModel& bsm) const
{
    return 0;
}

std::vector<double>  upandoutknockoutcalloption :: computeConfidenceInterval95(double price, const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     double variance = pricer.standarderror(*this, bsm);
     std::vector<double> confidenceinderval(2);
     double lb = price - (1.960*variance/sqrt(100000));
     double ub = price + (1.960*variance/sqrt(100000));
     confidenceinderval[0] = lb;
     confidenceinderval[1] = ub;
     return confidenceinderval;
 }

upandoutknockoutcalloption::~upandoutknockoutcalloption()
{
}

class upandoutknockoutputoption: public PathdependentOption {
    public:
    void setMaturity(double maturityuser) {
        maturity = maturityuser;
        }
    int getMaturity() const
    {
            return maturity;
        }
    void setStrike(double strikeuser) {
        strike = strikeuser;
        }
    int getStrike() const
    {
            return strike;
        }
    void setBarrier(double barrieruser) {
        barrier = barrieruser;
        }
    int getBarrier() const
    {
            return barrier;
        }
    double payoff (std::vector<double>& stock) const;
    double computedeltaanalytic (const BlackSholesModel&bsm) const;
    double computegammaanalytic (const BlackSholesModel&bsm) const;
    double computethetaanalytic (const BlackSholesModel&bsm) const;
    double computevegaanalytic (const BlackSholesModel&bsm) const;
    double computerhoanalytic (const BlackSholesModel&bsm) const;
    double computedeltaMonteCarlo (const BlackSholesModel&bsm) const;
    double computegammaMonteCarlo (const BlackSholesModel&bsm) const;
    double computethetaMonteCarlo (const BlackSholesModel&bsm) const;
    double computevegaMonteCarlo (const BlackSholesModel&bsm) const;
    double computerhoMonteCarlo (const BlackSholesModel&bsm) const;
    std::vector<double> computeConfidenceInterval95 (double price,const BlackSholesModel&bsm) const;
    bool isPathDependent() const {return true;};
    double priceclassic (const BlackSholesModel& bsm) const;
    double  priceMonteCarlo (const BlackSholesModel& bsm) const;
    ~upandoutknockoutputoption(); };

double upandoutknockoutputoption :: payoff (std::vector<double>& stock) const {
    double barrierofcall = barrier;
    for (int i=0; i<stock.size(); i++) {
        if (stock[i]>barrierofcall){
            return 0;}
        }
    double SMaturity= stock.back();
    double K=strike;
    return fmax(K-SMaturity,0); }

double upandoutknockoutputoption :: computedeltaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double upandoutknockoutputoption :: computegammaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double upandoutknockoutputoption :: computethetaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double upandoutknockoutputoption :: computevegaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double upandoutknockoutputoption :: computerhoanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double upandoutknockoutputoption :: computedeltaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.delta(*this, bsm);
 }

double upandoutknockoutputoption :: computegammaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.gamma(*this, bsm);
 }

double upandoutknockoutputoption :: computevegaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.vega(*this, bsm);
 }

double upandoutknockoutputoption :: computethetaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.theta(*this, bsm);
 }

double upandoutknockoutputoption :: computerhoMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.rho(*this, bsm);
 }


double upandoutknockoutputoption::priceMonteCarlo(const BlackSholesModel& bsm) const
{
    MonteCarloPricer pricer;
    return pricer.price( *this, bsm );
}

double upandoutknockoutputoption::priceclassic(const BlackSholesModel& bsm) const
{
    return 0;
}

std::vector<double>  upandoutknockoutputoption :: computeConfidenceInterval95(double price, const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     double variance = pricer.standarderror(*this, bsm);
     std::vector<double> confidenceinderval(2);
     double lb = price - (1.960*variance/sqrt(100000));
     double ub = price + (1.960*variance/sqrt(100000));
     confidenceinderval[0] = lb;
     confidenceinderval[1] = ub;
     return confidenceinderval;
 }


upandoutknockoutputoption::~upandoutknockoutputoption()
{
}

class LookbackCallOption: public PathdependentOption { // Type : Asian Geometric Call --> Path Dependent
    public:
    void setMaturity(double maturityuser) {
        maturity = maturityuser;
        }
    int getMaturity() const
    {
            return maturity;
        }
    void setStrike(double strikeuser) {
        strike = strikeuser;
        }
    int getStrike() const
    {
            return strike;
        }
    double payoff (std::vector<double>& stock) const;
    double computedeltaanalytic (const BlackSholesModel&bsm) const;
    double computegammaanalytic (const BlackSholesModel&bsm) const;
    double computethetaanalytic (const BlackSholesModel&bsm) const;
    double computevegaanalytic (const BlackSholesModel&bsm) const;
    double computerhoanalytic (const BlackSholesModel&bsm) const;
    double computedeltaMonteCarlo (const BlackSholesModel&bsm) const;
    double computegammaMonteCarlo (const BlackSholesModel&bsm) const;
    double computethetaMonteCarlo (const BlackSholesModel&bsm) const;
    double computevegaMonteCarlo (const BlackSholesModel&bsm) const;
    double computerhoMonteCarlo (const BlackSholesModel&bsm) const;
    std::vector<double> computeConfidenceInterval95 (double price,const BlackSholesModel&bsm) const;
    bool isPathDependent() const {return true;};
    double priceclassic (const BlackSholesModel& bsm) const;
    double  priceMonteCarlo (const BlackSholesModel& bsm) const;
    ~LookbackCallOption(); };

double LookbackCallOption :: payoff (std::vector<double>& stock) const {
    // Computation of the Payoff according to this precise type of option
    double maximum_element = *max_element(stock.begin(), stock.end());
    double K=strike;
    return fmax(maximum_element-K,0); }

double LookbackCallOption :: computedeltaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double LookbackCallOption :: computegammaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double LookbackCallOption :: computethetaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double LookbackCallOption :: computevegaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double LookbackCallOption :: computerhoanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double LookbackCallOption :: computedeltaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.delta(*this, bsm);
 }

double LookbackCallOption :: computegammaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.gamma(*this, bsm);
 }

double LookbackCallOption :: computevegaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.vega(*this, bsm);
 }

double LookbackCallOption :: computethetaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.theta(*this, bsm);
 }

double LookbackCallOption :: computerhoMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.rho(*this, bsm);
 }

double LookbackCallOption::priceMonteCarlo(const BlackSholesModel& bsm) const
{
    MonteCarloPricer pricer;
    return pricer.price( *this, bsm );
}

double LookbackCallOption::priceclassic(const BlackSholesModel& bsm) const
{
    return 0;
}

std::vector<double>  LookbackCallOption :: computeConfidenceInterval95(double price, const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     double variance = pricer.standarderror(*this, bsm);
     std::vector<double> confidenceinderval(2);
     double lb = price - (1.960*variance/sqrt(100000));
     double ub = price + (1.960*variance/sqrt(100000));
     confidenceinderval[0] = lb;
     confidenceinderval[1] = ub;
     return confidenceinderval;
 }

LookbackCallOption::~LookbackCallOption()
{
}

class LookbackPutOption: public PathdependentOption { // Type : Asian Geometric Call --> Path Dependent
    public:
    void setMaturity(double maturityuser) {
        maturity = maturityuser;
        }
    int getMaturity() const
    {
            return maturity;
        }
    void setStrike(double strikeuser) {
        strike = strikeuser;
        }
    int getStrike() const
    {
            return strike;
        }
    double payoff (std::vector<double>& stock) const;
    double computedeltaanalytic (const BlackSholesModel&bsm) const;
    double computegammaanalytic (const BlackSholesModel&bsm) const;
    double computethetaanalytic (const BlackSholesModel&bsm) const;
    double computevegaanalytic (const BlackSholesModel&bsm) const;
    double computerhoanalytic (const BlackSholesModel&bsm) const;
    double computedeltaMonteCarlo (const BlackSholesModel&bsm) const;
    double computegammaMonteCarlo (const BlackSholesModel&bsm) const;
    double computethetaMonteCarlo (const BlackSholesModel&bsm) const;
    double computevegaMonteCarlo (const BlackSholesModel&bsm) const;
    double computerhoMonteCarlo (const BlackSholesModel&bsm) const;
    std::vector<double> computeConfidenceInterval95 (double price,const BlackSholesModel&bsm) const;
    bool isPathDependent() const {return true;};
    double priceclassic (const BlackSholesModel& bsm) const;
    double  priceMonteCarlo (const BlackSholesModel& bsm) const;
    ~LookbackPutOption(); };

double LookbackPutOption :: payoff (std::vector<double>& stock) const {
    // Computation of the Payoff according to this precise type of option
    double minimum_element = *min_element(stock.begin(), stock.end());
    double K=strike;
    return fmax(K-minimum_element,0); }

double LookbackPutOption :: computedeltaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double LookbackPutOption :: computegammaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double LookbackPutOption :: computethetaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double LookbackPutOption :: computevegaanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double LookbackPutOption :: computerhoanalytic (const BlackSholesModel& bsm) const
 {
return 0; } // On ne peut pas calculer analytiquement donc on met à 0

double LookbackPutOption :: computedeltaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.delta(*this, bsm);
 }

double LookbackPutOption :: computegammaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.gamma(*this, bsm);
 }

double LookbackPutOption :: computevegaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.vega(*this, bsm);
 }

double LookbackPutOption :: computethetaMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.theta(*this, bsm);
 }

double LookbackPutOption :: computerhoMonteCarlo(const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     return pricer.rho(*this, bsm);
 }

double LookbackPutOption::priceMonteCarlo(const BlackSholesModel& bsm) const
{
    MonteCarloPricer pricer;
    return pricer.price( *this, bsm );
}

double LookbackPutOption::priceclassic(const BlackSholesModel& bsm) const
{
    return 0;
}

std::vector<double>  LookbackPutOption :: computeConfidenceInterval95(double price, const BlackSholesModel& bsm) const
 {
     MonteCarloPricer pricer;
     double variance = pricer.standarderror(*this, bsm);
     std::vector<double> confidenceinderval(2);
     double lb = price - (1.960*variance/sqrt(100000));
     double ub = price + (1.960*variance/sqrt(100000));
     confidenceinderval[0] = lb;
     confidenceinderval[1] = ub;
     return confidenceinderval;
 }

LookbackPutOption::~LookbackPutOption()
{
}


int main()
{
    // We have decided to code an iterative program therefore you just have to lauch it and follow the instructions
    int numbermethod;
    int number;
    int numbergreeks;
    double underlying;
    double volatility;
    double riskfreerate;
    double date;
    double maturity;
    double strike;
    double barrier;
    double strikeput;
    double strikecall;
    double strikecall2;
    double strikecall3;
    double strikecall4;
    int numbermethod12;
    std :: cout << "1: No Strategy just 1 Computation"<< std::endl;
    std :: cout << "2: Strategy"<< std::endl;
    std :: cout << "Choose your number"<< std::endl;
    std:: cin >> numbermethod;
    if (numbermethod==1){
    std :: cout << "1: Call European"<< std::endl;
    std :: cout << "2: Put European"<< std::endl;
    std :: cout << "3: Call Asian Arithmetic"<< std::endl;
    std :: cout << "4: Put Asian Arithmetic"<< std::endl;
    std :: cout << "5: Call Asian Geometric"<< std::endl;
    std :: cout << "6: Put Asian Arithmetic"<< std::endl;
    std :: cout << "7: Up-and-Out knock-out Call"<< std::endl;
    std :: cout << "8: Up-and-Out knock-out Put"<< std::endl;
    std :: cout << "9: Lookback Call"<< std::endl;
    std :: cout << "10:Lookback Put"<< std::endl;
    std :: cout << "Choose your number"<< std::endl;
    std:: cin >> number;
    std :: cout << "Choose your underlying price"<< std::endl;
    std:: cin >> underlying;
    std :: cout << "Choose your volatility : Do not write in %, write in decimals"<< std::endl;
    std :: cout << "For example : do not write 30% but 0.30 (follow the same rule for the other variables), Thanks you!"<< std::endl;
    std:: cin >> volatility;
    std :: cout << "Choose your risk free rate"<< std::endl;
    std:: cin >> riskfreerate;
    std :: cout << "Choose your date (only in integer numbers) in years"<< std::endl;
    std:: cin >> date;
    std :: cout << "Choose your Maturity (only in interger numbers) in years"<< std::endl;
    std:: cin >> maturity;
    std :: cout << "Choose your Strike"<< std::endl;
    std:: cin >> strike;
    std :: cout << "Do you want the greeks ? (1=Yes/ 0=No) / Advice : The computation of the Greeks with Monte Carlo will be more precise as the strike is closer to the price of the underlying"<< std::endl;
    std:: cin >> numbergreeks;
    BlackSholesModel test(underlying,volatility,riskfreerate,date,riskfreerate); // Entrer les grandeurs de Black-Scholes --> The drift has to be equal to the risk-free rate to generate the Risk-Free Price Path of the Stock
    if (number==1){
        EuropeanVanillaCallOption test1;
        test1.setStrike(strike);
        test1.setMaturity(maturity);
        double price =0;
        // We enable the user to chose which method he wants to use
        std :: cout << "Choose your Method: 1. Classic/Explicit"<< std::endl;
        std :: cout << "2. Monte Carlo"<< std::endl;
        std:: cin >> numbermethod12;
        if (numbermethod12==1){
            price= test1.priceclassic(test);}
        if (numbermethod12==2){
            price= test1.priceMonteCarlo(test);}
        std :: cout <<" The price of the Stock is "<<test.stockPrice<< std::endl;
        std :: cout <<" The Volatility of the Stock is "<<test.volatility * 100<<" %"<<std::endl;
        std :: cout <<" The Risk-free rate is "<<test.riskFreeRate * 100<<" %"<< std::endl;
        std :: cout <<" The Date is "<<test.date<< std::endl;
        std :: cout <<" The Drift is "<<test.drift*100<<" %"<< std::endl;
        std :: cout <<" The Maturity is of the European Vanilla Call is "<<maturity<< std::endl;
        std :: cout <<" The Strike is of the European Vanilla Call is "<<strike<< std::endl;
        if (numbergreeks==1){
            if (numbermethod12==1){
                double delta = test1.computedeltaanalytic(test);
                double gamma = test1.computegammaanalytic(test);
                double theta = test1.computethetaanalytic(test);
                double vega = test1.computevegaanalytic(test);
                double rho = test1.computerhoanalytic(test);
                std :: cout <<"The Delta of the Call:"<<delta<< std::endl;
                std :: cout <<" The Gamma of the Call is "<<gamma<< std::endl;
                std :: cout <<" The Theta of the Call is "<<theta<< std::endl;
                std :: cout <<" The Vega of the Call is "<<vega<< std::endl;
                std :: cout <<" The Rho of the Call is "<<rho<< std::endl;
            }
            if (numbermethod12==2){
                double delta = test1.computedeltaMonteCarlo(test);
                double gamma = test1.computegammaMonteCarlo(test);
                double theta = test1.computethetaMonteCarlo(test);
                double vega = test1.computevegaMonteCarlo(test);
                double rho = test1.computerhoanalytic(test);
                std :: cout <<"The Delta of the Call:"<<delta<< std::endl;
                std :: cout <<"The Gamma of the Call:"<<gamma<< std::endl;
                std :: cout <<" The Theta of the Call is "<<theta<< std::endl;
                std :: cout <<" The Vega of the Call is "<<vega<< std::endl;
                std :: cout <<" The Rho of the Call is "<<rho<< std::endl;
            }
        }
        std :: cout <<" The price of the European Vanilla Call is "<<price<< std::endl;
        if (numbermethod12==2){
            std::vector<double> confidenceinterval = test1.computeConfidenceInterval95(price, test);
            double lb = confidenceinterval[0];
            double ub = confidenceinterval[1];
            std::cout << " Confidence Interval at level 95% : ";
            std::cout << "[" << std::setprecision(2) << std::fixed <<lb<< " , "<< std::setprecision(2) << std::fixed<< ub<< "]" << std::endl << std::endl;
        }
        std :: cout <<" Thank you for using this program "<< std::endl;
    }
    else if(number==2){
        EuropeanVanillaPutOption test1;
        test1.setStrike(strike);
        test1.setMaturity(maturity);
        double price =0;
        std :: cout << "Choose your Method: 1. Classic/Explicit"<< std::endl;
        std :: cout << "2. Monte Carlo"<< std::endl;
        std:: cin >> numbermethod12;
        if (numbermethod12==1){
            price= test1.priceclassic(test);}
        if (numbermethod12==2){
            price= test1.priceMonteCarlo(test);}
        std :: cout <<" The price of the Stock is "<<test.stockPrice<< std::endl;
        std :: cout <<" The Volatility of the Stock is "<<test.volatility * 100<<" %"<<std::endl;
        std :: cout <<" The Risk-free rate is "<<test.riskFreeRate * 100<<" %"<< std::endl;
        std :: cout <<" The Date is "<<test.date<< std::endl;
        std :: cout <<" The Drift is "<<test.drift*100<<" %"<< std::endl;
        std :: cout <<" The Maturity is of the European Vanilla Put is "<<maturity<< std::endl;
        std :: cout <<" The Strike is of the European Vanilla Put is "<<strike<< std::endl;
        if (numbergreeks==1){
            if (numbermethod12==1){
                double delta = test1.computedeltaanalytic(test);
                double gamma = test1.computegammaanalytic(test);
                double theta = test1.computethetaanalytic(test);
                double vega = test1.computevegaanalytic(test);
                double rho = test1.computerhoanalytic(test);
                std :: cout <<"The Delta of the Put:"<<delta<< std::endl;
                std :: cout <<" The Gamma of the Put is "<<gamma<< std::endl;
                std :: cout <<" The Theta of the Put is "<<theta<< std::endl;
                std :: cout <<" The Vega of the Put is "<<vega<< std::endl;
                std :: cout <<" The Rho of the Put is "<<rho<< std::endl;
            }
            if (numbermethod12==2){
                double delta = test1.computedeltaMonteCarlo(test);
                double gamma = test1.computegammaMonteCarlo(test);
                double theta = test1.computethetaMonteCarlo(test);
                double vega = test1.computevegaMonteCarlo(test);
                double rho = test1.computerhoanalytic(test);
                std :: cout <<"The Delta of the Put:"<<delta<< std::endl;
                std :: cout <<"The Gamma of the Put:"<<gamma<< std::endl;
                std :: cout <<" The Theta of the Put is "<<theta<< std::endl;
                std :: cout <<" The Vega of the Put is "<<vega<< std::endl;
                std :: cout <<" The Rho of the Put is "<<rho<< std::endl;
            }
        }
        std :: cout <<" The price of the European Vanilla Put is "<<price<< std::endl;
        if (numbermethod12==2){
            std::vector<double> confidenceinterval = test1.computeConfidenceInterval95(price, test);
            double lb = confidenceinterval[0];
            double ub = confidenceinterval[1];
            std::cout << " Confidence Interval at level 95% : ";
            std::cout << "[" << std::setprecision(2) << std::fixed <<lb<< " , "<< std::setprecision(2) << std::fixed<< ub<< "]" << std::endl << std::endl;
        }
        std :: cout <<" Thank you for using this program "<< std::endl;
    }
    else if (number==3){
        AsianArithmeticCallOption test1;
        test1.setStrike(strike);
        test1.setMaturity(maturity);
        double price =0;
        price= test1.priceMonteCarlo(test);
        std :: cout <<" The price of the Stock is "<<test.stockPrice<< std::endl;
        std :: cout <<" The Volatility of the Stock is "<<test.volatility * 100<<" %"<<std::endl;
        std :: cout <<" The Risk-free rate is "<<test.riskFreeRate * 100<<" %"<< std::endl;
        std :: cout <<" The Date is "<<test.date<< std::endl;
        std :: cout <<" The Drift is "<<test.drift*100<<" %"<< std::endl;
        std :: cout <<" The Maturity is of the Asian Arithmetic Call is "<<maturity<< std::endl;
        std :: cout <<" The Strike is of the Asian Arithmetic Call is "<<strike<< std::endl;
        if (numbergreeks==1){
                double delta = test1.computedeltaMonteCarlo(test);
                double gamma = test1.computegammaMonteCarlo(test);
                double theta = test1.computethetaMonteCarlo(test);
                double vega = test1.computevegaMonteCarlo(test);
                double rho = test1.computerhoMonteCarlo(test);
                std :: cout <<"The Delta of the Call:"<<delta<< std::endl;
                std :: cout <<"The Gamma of the Call:"<<gamma<< std::endl;
                std :: cout <<" The Theta of the Call is "<<theta<< std::endl;
                std :: cout <<" The Vega of the Call is "<<vega<< std::endl;
                std :: cout <<" The Rho of the Call is "<<rho<< std::endl;
        }
        std :: cout <<" The price of the Asian Arithmetic Call is "<<price<< std::endl;
        std::vector<double> confidenceinterval = test1.computeConfidenceInterval95(price, test);
        double lb = confidenceinterval[0];
        double ub = confidenceinterval[1];
        std::cout << " Confidence Interval at level 95% : ";
        std::cout << "[" << std::setprecision(2) << std::fixed <<lb<< " , "<< std::setprecision(2) << std::fixed<< ub<< "]" << std::endl << std::endl;
        std :: cout <<" Thank you for using this program "<< std::endl;
    }
    else if(number==4){
        AsianArithmeticPutOption test1;
        test1.setStrike(strike);
        test1.setMaturity(maturity);
        double price =0;
        price= test1.priceMonteCarlo(test);
        std :: cout <<" The price of the Stock is "<<test.stockPrice<< std::endl;
        std :: cout <<" The Volatility of the Stock is "<<test.volatility * 100<<" %"<<std::endl;
        std :: cout <<" The Risk-free rate is "<<test.riskFreeRate * 100<<" %"<< std::endl;
        std :: cout <<" The Date is "<<test.date<< std::endl;
        std :: cout <<" The Drift is "<<test.drift*100<<" %"<< std::endl;
        std :: cout <<" The Maturity is of the Asian Arithmetic Put is "<<maturity<< std::endl;
        std :: cout <<" The Strike is of the Asian Arithmetic Put is "<<strike<< std::endl;
        if (numbergreeks==1){
                double delta = test1.computedeltaMonteCarlo(test);
                double gamma = test1.computegammaMonteCarlo(test);
                double theta = test1.computethetaMonteCarlo(test);
                double vega = test1.computevegaMonteCarlo(test);
                double rho = test1.computerhoMonteCarlo(test);
                std :: cout <<"The Delta of the Put:"<<delta<< std::endl;
                std :: cout <<"The Gamma of the Put:"<<gamma<< std::endl;
                std :: cout <<" The Theta of the Put is "<<theta<< std::endl;
                std :: cout <<" The Vega of the Put is "<<vega<< std::endl;
                std :: cout <<" The Rho of the Put is "<<rho<< std::endl;
        }
        std :: cout <<" The price of the Asian Arithmetic Put is "<<price<< std::endl;
        std::vector<double> confidenceinterval = test1.computeConfidenceInterval95(price, test);
        double lb = confidenceinterval[0];
        double ub = confidenceinterval[1];
        std::cout << " Confidence Interval at level 95% : ";
        std::cout << "[" << std::setprecision(2) << std::fixed <<lb<< " , "<< std::setprecision(2) << std::fixed<< ub<< "]" << std::endl << std::endl;
        std :: cout <<" Thank you for using this program "<< std::endl;
    }
    else if(number==5){
        AsianGeometricCallOption test1;
        test1.setStrike(strike);
        test1.setMaturity(maturity);
        double price =0;
        price= test1.priceMonteCarlo(test);
        std :: cout <<" The price of the Stock is "<<test.stockPrice<< std::endl;
        std :: cout <<" The Volatility of the Stock is "<<test.volatility * 100<<" %"<<std::endl;
        std :: cout <<" The Risk-free rate is "<<test.riskFreeRate * 100<<" %"<< std::endl;
        std :: cout <<" The Date is "<<test.date<< std::endl;
        std :: cout <<" The Drift is "<<test.drift*100<<" %"<< std::endl;
        std :: cout <<" The Maturity is of the Asian Geometric Call is "<<maturity<< std::endl;
        std :: cout <<" The Strike is of the Asian Geometric Call is "<<strike<< std::endl;
        if (numbergreeks==1){
                double delta = test1.computedeltaMonteCarlo(test);
                double gamma = test1.computegammaMonteCarlo(test);
                double theta = test1.computethetaMonteCarlo(test);
                double vega = test1.computevegaMonteCarlo(test);
                double rho = test1.computerhoMonteCarlo(test);
                std :: cout <<"The Delta of the Call:"<<delta<< std::endl;
                std :: cout <<"The Gamma of the Call:"<<gamma<< std::endl;
                std :: cout <<" The Theta of the Call is "<<theta<< std::endl;
                std :: cout <<" The Vega of the Call is "<<vega<< std::endl;
                std :: cout <<" The Rho of the Call is "<<rho<< std::endl;
        }
        std :: cout <<" The price of the Asian Geometric Call is "<<price<< std::endl;
        std::vector<double> confidenceinterval = test1.computeConfidenceInterval95(price, test);
        double lb = confidenceinterval[0];
        double ub = confidenceinterval[1];
        std::cout << " Confidence Interval at level 95% : ";
        std::cout << "[" << std::setprecision(2) << std::fixed <<lb<< " , "<< std::setprecision(2) << std::fixed<< ub<< "]" << std::endl << std::endl;
        std :: cout <<" Thank you for using this program "<< std::endl;
    }
    else if(number==6){
        AsianGeometricPutOption test1;
        test1.setStrike(strike);
        test1.setMaturity(maturity);
        double price =0;
        price= test1.priceMonteCarlo(test);
        std :: cout <<" The price of the Stock is "<<test.stockPrice<< std::endl;
        std :: cout <<" The Volatility of the Stock is "<<test.volatility * 100<<" %"<<std::endl;
        std :: cout <<" The Risk-free rate is "<<test.riskFreeRate * 100<<" %"<< std::endl;
        std :: cout <<" The Date is "<<test.date<< std::endl;
        std :: cout <<" The Drift is "<<test.drift*100<<" %"<< std::endl;
        std :: cout <<" The Maturity is of the Asian Geometric Put is "<<maturity<< std::endl;
        std :: cout <<" The Strike is of the Asian Geometric Put is "<<strike<< std::endl;
        if (numbergreeks==1){
                double delta = test1.computedeltaMonteCarlo(test);
                double gamma = test1.computegammaMonteCarlo(test);
                double theta = test1.computethetaMonteCarlo(test);
                double vega = test1.computevegaMonteCarlo(test);
                double rho = test1.computerhoMonteCarlo(test);
                std :: cout <<"The Delta of the Put:"<<delta<< std::endl;
                std :: cout <<"The Gamma of the Put:"<<gamma<< std::endl;
                std :: cout <<" The Theta of the Put is "<<theta<< std::endl;
                std :: cout <<" The Vega of the Put is "<<vega<< std::endl;
                std :: cout <<" The Rho of the Put is "<<rho<< std::endl;
        }
        std :: cout <<" The price of the Asian Geometric Put is "<<price<< std::endl;
        std::vector<double> confidenceinterval = test1.computeConfidenceInterval95(price, test);
        double lb = confidenceinterval[0];
        double ub = confidenceinterval[1];
        std::cout << " Confidence Interval at level 95% : ";
        std::cout << "["<< std::setprecision(2) << std::fixed <<lb<< " , "<< std::setprecision(2) << std::fixed<< ub<<"]" << std::endl << std::endl;
        std :: cout <<" Thank you for using this program "<< std::endl;
    }
    else if(number==7){
        std :: cout << "Choose your Barrier"<< std::endl;
        std:: cin >> barrier;
        upandoutknockoutcalloption test1;
        test1.setStrike(strike);
        test1.setMaturity(maturity);
        test1.setBarrier(barrier);
        double price =0;
        price= test1.priceMonteCarlo(test);
        std :: cout <<" The price of the Stock is "<<test.stockPrice<< std::endl;
        std :: cout <<" The Volatility of the Stock is "<<test.volatility * 100<<" %"<<std::endl;
        std :: cout <<" The Risk-free rate is "<<test.riskFreeRate * 100<<" %"<< std::endl;
        std :: cout <<" The Date is "<<test.date<< std::endl;
        std :: cout <<" The Drift is "<<test.drift*100<<" %"<< std::endl;
        std :: cout <<" The Maturity of this Up-and-Out Knock-out Call is "<<maturity<< std::endl;
        std :: cout <<" The Strike of this Up-and-Out Knock-out Call is "<<strike<< std::endl;
        std :: cout <<" The Barrier of this Up-and-Out Knock-out Call is "<<barrier<< std::endl;
        if (numbergreeks==1){
                double delta = test1.computedeltaMonteCarlo(test);
                double gamma = test1.computegammaMonteCarlo(test);
                double theta = test1.computethetaMonteCarlo(test);
                double vega = test1.computevegaMonteCarlo(test);
                double rho = test1.computerhoMonteCarlo(test);
                std :: cout <<"The Delta of the Call:"<<delta<< std::endl;
                std :: cout <<"The Gamma of the Call:"<<gamma<< std::endl;
                std :: cout <<" The Theta of the Call is "<<theta<< std::endl;
                std :: cout <<" The Vega of the Call is "<<vega<< std::endl;
                std :: cout <<" The Rho of the Call is "<<rho<< std::endl;
        }
        std :: cout <<" The price of this Up-and-Out Knock-out Call is "<<price<< std::endl;
        std::vector<double> confidenceinterval = test1.computeConfidenceInterval95(price, test);
        double lb = confidenceinterval[0];
        double ub = confidenceinterval[1];
        std::cout << " Confidence Interval at level 95% : ";
        std::cout << "["<< std::setprecision(2) << std::fixed <<lb<< " , "<< std::setprecision(2) << std::fixed<< ub<<"]" << std::endl << std::endl;
        std :: cout <<" Thank you for using this program "<< std::endl;
    }
    else if (number==8){
            std :: cout << "Choose your Barrier"<< std::endl;
            std:: cin >> barrier;
            upandoutknockoutputoption test1;
            test1.setStrike(strike);
            test1.setMaturity(maturity);
            test1.setBarrier(barrier);
            double price =0;
            price= test1.priceMonteCarlo(test);
            std :: cout <<" The price of the Stock is "<<test.stockPrice<< std::endl;
            std :: cout <<" The Volatility of the Stock is "<<test.volatility * 100<<" %"<<std::endl;
            std :: cout <<" The Risk-free rate is "<<test.riskFreeRate * 100<<" %"<< std::endl;
            std :: cout <<" The Date is "<<test.date<< std::endl;
            std :: cout <<" The Drift is "<<test.drift*100<<" %"<< std::endl;
            std :: cout <<" The Maturity of this Up-and-Out Knock-out Put is "<<maturity<< std::endl;
            std :: cout <<" The Strike of this Up-and-Out Knock-out Put is "<<strike<< std::endl;
            std :: cout <<" The Barrier of this Up-and-Out Knock-out Put is "<<barrier<< std::endl;
        if (numbergreeks==1){
                double delta = test1.computedeltaMonteCarlo(test);
                double gamma = test1.computegammaMonteCarlo(test);
                double theta = test1.computethetaMonteCarlo(test);
                double vega = test1.computevegaMonteCarlo(test);
                double rho = test1.computerhoMonteCarlo(test);
                std :: cout <<"The Delta of the Put:"<<delta<< std::endl;
                std :: cout <<"The Gamma of the Put:"<<gamma<< std::endl;
                std :: cout <<" The Theta of the Put is "<<theta<< std::endl;
                std :: cout <<" The Vega of the Put is "<<vega<< std::endl;
                std :: cout <<" The Rho of the Put is "<<rho<< std::endl;
        }
        std :: cout <<" The price of this Up-and-Out Knock-out Put is "<<price<< std::endl;
        std::vector<double> confidenceinterval = test1.computeConfidenceInterval95(price, test);
        double lb = confidenceinterval[0];
        double ub = confidenceinterval[1];
        std::cout << " Confidence Interval at level 95% : ";
        std::cout << "["<< std::setprecision(2) << std::fixed <<lb<< " , "<< std::setprecision(2) << std::fixed<< ub<<"]" << std::endl << std::endl;
            std :: cout <<" Thank you for using this program "<< std::endl;
    }
    else if (number==9){
        LookbackCallOption test1;
        test1.setStrike(strike);
        test1.setMaturity(maturity);
        double price =0;
        price= test1.priceMonteCarlo(test);
        std :: cout <<" The price of the Stock is "<<test.stockPrice<< std::endl;
        std :: cout <<" The Volatility of the Stock is "<<test.volatility * 100<<" %"<<std::endl;
        std :: cout <<" The Risk-free rate is "<<test.riskFreeRate * 100<<" %"<< std::endl;
        std :: cout <<" The Date is "<<test.date<< std::endl;
        std :: cout <<" The Drift is "<<test.drift*100<<" %"<< std::endl;
        std :: cout <<" The Maturity is of the Lookback Call is "<<maturity<< std::endl;
        std :: cout <<" The Strike is of the Lookback Call is "<<strike<< std::endl;
        if (numbergreeks==1){
                double delta = test1.computedeltaMonteCarlo(test);
                double gamma = test1.computegammaMonteCarlo(test);
                double theta = test1.computethetaMonteCarlo(test);
                double vega = test1.computevegaMonteCarlo(test);
                double rho = test1.computerhoMonteCarlo(test);
                std :: cout <<"The Delta of the Call:"<<delta<< std::endl;
                std :: cout <<"The Gamma of the Call:"<<gamma<< std::endl;
                std :: cout <<" The Theta of the Call is "<<theta<< std::endl;
                std :: cout <<" The Vega of the Call is "<<vega<< std::endl;
                std :: cout <<" The Rho of the Call is "<<rho<< std::endl;
        }
        std :: cout <<" The price of the Lookback Call is "<<price<< std::endl;
        std::vector<double> confidenceinterval = test1.computeConfidenceInterval95(price, test);
        double lb = confidenceinterval[0];
        double ub = confidenceinterval[1];
        std::cout << " Confidence Interval at level 95% : ";
        std::cout << "[" << std::setprecision(2) << std::fixed <<lb<< " , "<< std::setprecision(2) << std::fixed<< ub<< "]" << std::endl << std::endl;
        std :: cout <<" Thank you for using this program "<< std::endl;
    }
    else if (number==10){
        LookbackPutOption test1;
        test1.setStrike(strike);
        test1.setMaturity(maturity);
        double price =0;
        price= test1.priceMonteCarlo(test);
        std :: cout <<" The price of the Stock is "<<test.stockPrice<< std::endl;
        std :: cout <<" The Volatility of the Stock is "<<test.volatility * 100<<" %"<<std::endl;
        std :: cout <<" The Risk-free rate is "<<test.riskFreeRate * 100<<" %"<< std::endl;
        std :: cout <<" The Date is "<<test.date<< std::endl;
        std :: cout <<" The Drift is "<<test.drift*100<<" %"<< std::endl;
        std :: cout <<" The Maturity is of the Lookback Put is "<<maturity<< std::endl;
        std :: cout <<" The Strike is of the Lookback Put is "<<strike<< std::endl;
        if (numbergreeks==1){
                double delta = test1.computedeltaMonteCarlo(test);
                double gamma = test1.computegammaMonteCarlo(test);
                double theta = test1.computethetaMonteCarlo(test);
                double vega = test1.computevegaMonteCarlo(test);
                double rho = test1.computerhoMonteCarlo(test);
                std :: cout <<"The Delta of the Put:"<<delta<< std::endl;
                std :: cout <<"The Gamma of the Put:"<<gamma<< std::endl;
                std :: cout <<" The Theta of the Put is "<<theta<< std::endl;
                std :: cout <<" The Vega of the Put is "<<vega<< std::endl;
                std :: cout <<" The Rho of the Put is "<<rho<< std::endl;
        }
        std :: cout <<" The price of the Lookback Put is "<<price<< std::endl;
        std::vector<double> confidenceinterval = test1.computeConfidenceInterval95(price, test);
        double lb = confidenceinterval[0];
        double ub = confidenceinterval[1];
        std::cout << " Confidence Interval at level 95% : ";
        std::cout << "[" << std::setprecision(2) << std::fixed <<lb<< " , "<< std::setprecision(2) << std::fixed<< ub<< "]" << std::endl << std::endl;
        std :: cout <<" Thank you for using this program "<< std::endl;
    }
    }
    if (numbermethod==2){
        std :: cout << "1: Stellage Bought : Buy 1 European Vanilla Call and 1 European Vanilla Put both with the same Strike"<< std::endl;
        std :: cout << "2: Strangle Bought: Buy 1 European Vanilla Call and 1 European Vanilla Put with two different Strikes (The Strike of the Put is the smaller one) "<< std::endl;
        std :: cout << "3: Bull Spread : Buy 1 European Vanilla Call and sell 1 European Vanilla Call with two different Strikes (The Strike of the Call you buy is the smaller one) "<< std::endl;
        std :: cout << "4: Bear Spread: Buy 1 European Vanilla Call and sell 1 European Vanilla Call with two different Strikes (The Strike of the Call you sell is the smaller one) "<< std::endl;
        std :: cout << "5: Butterfly Bought: Buy 1 European Vanilla Call, sell 2 European Vanilla Calls and finally Buy 1 European Vanilla Call, with three different Strikes for each 'move' (The Strike of the first Call you buy is the smaller one, the Strike of the 2nd Call you buy is the larger one and finally the Strike of the 2 Calls you'll sell is between the 2 previous ones ) "<< std::endl;
        std :: cout << "6: Condor Sold: Sell 1 European Vanilla Call, Buy 1 European Vanilla Call, Buy 1 European Vanilla Call and finally Sell 1 European Vanilla Call, with four different Strikes for each 'move' (The Strike of the first Call you sell is the smaller one, then you have the Strike of the 1st Call you buy, then the Strike of the 2nd Call you buy and finally the Strike of the 2nd Call you sell) "<< std::endl;
        std :: cout << "Choose your number"<< std::endl;
        std:: cin >> number;
        std :: cout << "Choose your underlying price"<< std::endl;
        std:: cin >> underlying;
        std :: cout << "Choose your volatility : Do not write in %, write in decimals"<< std::endl;
        std :: cout << "For example : do not write 30% but 0.30 (follow the same rule for the other variables), Thanks you!"<< std::endl;
        std:: cin >> volatility;
        std :: cout << "Choose your risk free rate"<< std::endl;
        std:: cin >> riskfreerate;
        std :: cout << "Choose your date (only in integer numbers) in years"<< std::endl;
        std:: cin >> date;
        BlackSholesModel test(underlying,volatility,riskfreerate,date,riskfreerate); /// Entrer les grandeurs de Black-Scholes
        if (number==1){
            std :: cout << "Choose your Maturity (only in interger numbers) in years"<< std::endl;
            std:: cin >> maturity;
            std :: cout << "Choose your Strike"<< std::endl;
            std:: cin >> strike;
            EuropeanVanillaCallOption test1;
            test1.setStrike(strike);
            test1.setMaturity(maturity);
            double price =0;
            std :: cout << "Choose your Method: 1. Classic/Explicit"<< std::endl;
            std :: cout << "2. Monte Carlo"<< std::endl;
            std:: cin >> numbermethod12;
            if (numbermethod12==1){
                price= test1.priceclassic(test);}
            if (numbermethod12==2){
                price= test1.priceMonteCarlo(test);}
            EuropeanVanillaPutOption test2;
            test2.setStrike(strike);
            test2.setMaturity(maturity);
            double price2 =0;
            if (numbermethod12==1){
                price2= test2.priceclassic(test);}
            if (numbermethod12==2){
                price2= test2.priceMonteCarlo(test);}
            std :: cout <<" The price of the Stock is "<<test.stockPrice<< std::endl;
            std :: cout <<" The Volatility of the Stock is "<<test.volatility * 100<<" %"<<std::endl;
            std :: cout <<" The Risk-free rate is "<<test.riskFreeRate * 100<<" %"<< std::endl;
            std :: cout <<" The Date is "<<test.date<< std::endl;
            std :: cout <<" The Drift is "<<test.drift*100<<" %"<< std::endl;
            std :: cout <<" The Maturity of this Strategy is "<<maturity<< std::endl;
            std :: cout <<" The shared Strike of this Strategy is "<<strike<< std::endl;
            std :: cout <<" The price of the Call you'll buy is "<<price<< std::endl;
            std :: cout <<" The price of the Put you'll buy is "<<price2<< std::endl;
            std :: cout <<" The price of this Strategy is "<<price+price2<< std::endl;
            if (price+price2>=0)
            {
            std :: cout <<" You have to spend this money to use this strategy "<< std::endl;
            }
            if (price+price2<0)
            {
            std :: cout <<" You earn this money (the absolute value of the result just above) to use this strategy "<<std::endl;
            }
            std :: cout <<" Thank you for using this program "<< std::endl;
        }
        else if(number==2){
            std :: cout << "Choose your Maturity (only in interger numbers) in years"<< std::endl;
            std:: cin >> maturity;
            std :: cout << "Choose your Strike for the Put, it has to be smaller than the Call's one"<< std::endl;
            std:: cin >> strikeput;
            std :: cout << "Choose your Strike for the Call"<< std::endl;
            std:: cin >> strikecall;
            EuropeanVanillaPutOption test1;
            test1.setStrike(strikeput);
            test1.setMaturity(maturity);
            double price =0;
            std :: cout << "Choose your Method: 1. Classic/Explicit"<< std::endl;
            std :: cout << "2. Monte Carlo"<< std::endl;
            std:: cin >> numbermethod12;
            if (numbermethod12==1){
                price= test1.priceclassic(test);}
            if (numbermethod12==2){
                price= test1.priceMonteCarlo(test);}
            EuropeanVanillaCallOption test2;
            test2.setStrike(strikecall);
            test2.setMaturity(maturity);
            double price2 =0;
            if (numbermethod12==1){
                price2= test2.priceclassic(test);}
            if (numbermethod12==2){
                price2= test2.priceMonteCarlo(test);}
            std :: cout <<" The price of the Stock is "<<test.stockPrice<< std::endl;
            std :: cout <<" The Volatility of the Stock is "<<test.volatility * 100<<" %"<<std::endl;
            std :: cout <<" The Risk-free rate is "<<test.riskFreeRate * 100<<" %"<< std::endl;
            std :: cout <<" The Date is "<<test.date<< std::endl;
            std :: cout <<" The Drift is "<<test.drift*100<<" %"<< std::endl;
            std :: cout <<" The Maturity of this Strategy is "<<maturity<< std::endl;
            std :: cout <<" The Strike of the Put in this Strategy is "<<strikeput<<" and the Strike of the Call in this Strategy is"<<strikecall<< std::endl;
            std :: cout <<" The price of the Call you'll buy is "<<price2<< std::endl;
            std :: cout <<" The price of the Put you'll buy is "<<price<< std::endl;
            std :: cout <<" The price of this Strategy is "<<price+price2<< std::endl;
            if (price+price2>=0)
            {
            std :: cout <<" You have to spend this money to use this strategy "<< std::endl;
            }
            if (price+price2<0)
            {
            std :: cout <<" You earn this money (the absolute value of the result just above) to use this strategy "<<std::endl;
            }
            std :: cout <<" Thank you for using this program "<< std::endl;
        }
        else if (number==3){
            std :: cout << "Choose your Maturity (only in interger numbers) in years"<< std::endl;
            std:: cin >> maturity;
            std :: cout << "Choose your Strike for the Call you will buy, it has to be smaller than the one of the Call you'll sell"<< std::endl;
            std:: cin >> strikecall;
            std :: cout << "Choose your Strike for the Call you'll sell"<< std::endl;
            std:: cin >> strikecall2;
            EuropeanVanillaCallOption test1;
            test1.setStrike(strikecall);
            test1.setMaturity(maturity);
            double price =0;
            std :: cout << "Choose your Method: 1. Classic/Explicit"<< std::endl;
            std :: cout << "2. Monte Carlo"<< std::endl;
            std:: cin >> numbermethod12;
            if (numbermethod12==1){
                price= test1.priceclassic(test);}
            if (numbermethod12==2){
                price= test1.priceMonteCarlo(test);}
            EuropeanVanillaCallOption test2;
            test2.setStrike(strikecall2);
            test2.setMaturity(maturity);
            double price2 =0;
            if (numbermethod12==1){
                price2= test2.priceclassic(test);}
            if (numbermethod12==2){
                price2= test2.priceMonteCarlo(test);}
            std :: cout <<" The price of the Stock is "<<test.stockPrice<< std::endl;
            std :: cout <<" The Volatility of the Stock is "<<test.volatility * 100<<" %"<<std::endl;
            std :: cout <<" The Risk-free rate is "<<test.riskFreeRate * 100<<" %"<< std::endl;
            std :: cout <<" The Date is "<<test.date<< std::endl;
            std :: cout <<" The Drift is "<<test.drift*100<<" %"<< std::endl;
            std :: cout <<" The Maturity of this Strategy is "<<maturity<< std::endl;
            std :: cout <<" The Strike of the Call you'll buy in this Strategy is "<<strikecall<< std::endl;
            std :: cout <<" The Strike of the Call you'll sell in this Strategy is "<<strikecall2<< std::endl;
            std :: cout <<" The price of the Call you'll buy is "<<price<< std::endl;
            std :: cout <<" The price of the Call you'll sell is "<<price2<< std::endl;
            std :: cout <<" The price of this Strategy is (you earn the money of the call you'll sell and you spend the money of the call you'll buy) : "<<price-price2<< std::endl;
            if (price-price2>=0)
            {
            std :: cout <<" You have to spend this money to use this strategy "<< std::endl;
            }
            if (price-price2<0)
            {
            std :: cout <<" You earn this money (the absolute value of the result just above) to use this strategy "<<std::endl;
            }
            std :: cout <<" Thank you for using this program "<<std::endl;
        }
        else if(number==4){
            std :: cout << "Choose your Maturity (only in interger numbers) in years"<< std::endl;
            std:: cin >> maturity;
            std :: cout << "Choose your Strike for the Call you will sell, it has to be smaller than the one of the Call you'll sell"<< std::endl;
            std:: cin >> strikecall;
            std :: cout << "Choose your Strike for the Call you'll buy"<< std::endl;
            std:: cin >> strikecall2;
            EuropeanVanillaCallOption test1;
            test1.setStrike(strikecall);
            test1.setMaturity(maturity);
            double price =0;
            std :: cout << "Choose your Method: 1. Classic/Explicit"<< std::endl;
            std :: cout << "2. Monte Carlo"<< std::endl;
            std:: cin >> numbermethod12;
            if (numbermethod12==1){
                price= test1.priceclassic(test);}
            if (numbermethod12==2){
                price= test1.priceMonteCarlo(test);}
            EuropeanVanillaCallOption test2;
            test2.setStrike(strikecall2);
            test2.setMaturity(maturity);
            double price2 =0;
            if (numbermethod12==1){
                price2= test2.priceclassic(test);}
            if (numbermethod12==2){
                price2= test2.priceMonteCarlo(test);}
            std :: cout <<" The price of the Stock is "<<test.stockPrice<< std::endl;
            std :: cout <<" The Volatility of the Stock is "<<test.volatility * 100<<" %"<<std::endl;
            std :: cout <<" The Risk-free rate is "<<test.riskFreeRate * 100<<" %"<< std::endl;
            std :: cout <<" The Date is "<<test.date<< std::endl;
            std :: cout <<" The Drift is "<<test.drift*100<<" %"<< std::endl;
            std :: cout <<" The Maturity of this Strategy is "<<maturity<< std::endl;
            std :: cout <<" The Strike of the Call you'll buy in this Strategy is "<<strikecall<< std::endl;
            std :: cout <<" The Strike of the Call you'll sell in this Strategy is "<<strikecall2<< std::endl;
            std :: cout <<" The price of the Call you'll buy is "<<price2<< std::endl;
            std :: cout <<" The price of the Call you'll sell is "<<price<< std::endl;
            std :: cout <<" The price of this Strategy is (you earn the money of the call you'll sell and you spend the money of the call you'll buy) : "<<price2-price<< std::endl;
            if (price2-price>=0)
            {
            std :: cout <<" You have to spend this money to use this strategy "<< std::endl;
            }
            if (price2-price<0)
            {
            std :: cout <<" You earn this money (the absolute value of the result just above) to use this strategy "<<std::endl;
            }
            std :: cout <<" Thank you for using this program "<<std::endl;
            }
        else if(number==5){
            std :: cout << "Choose your Maturity (only in interger numbers) in years"<< std::endl;
            std:: cin >> maturity;
            std :: cout << "Choose your Strike for the first Call you will buy, it has to be the smaller"<< std::endl;
            std:: cin >> strikecall;
            std :: cout << "Choose your Strike for the 2 Calls you'll sell"<< std::endl;
            std:: cin >> strikecall2;
            std :: cout << "Choose your Strike for the 2nd Call you'll buy, it has to be the bigger one"<< std::endl;
            std:: cin >> strikecall3;
            EuropeanVanillaCallOption test1;
            test1.setStrike(strikecall);
            test1.setMaturity(maturity);
            double price =0;
            std :: cout << "Choose your Method: 1. Classic/Explicit"<< std::endl;
            std :: cout << "2. Monte Carlo"<< std::endl;
            std:: cin >> numbermethod12;
            if (numbermethod12==1){
                price= test1.priceclassic(test);}
            if (numbermethod12==2){
                price= test1.priceMonteCarlo(test);}
            EuropeanVanillaCallOption test2;
            test2.setStrike(strikecall2);
            test2.setMaturity(maturity);
            double price2 =0;
            if (numbermethod12==1){
                price2= test2.priceclassic(test);}
            if (numbermethod12==2){
                price2= test2.priceMonteCarlo(test);}
            EuropeanVanillaCallOption test3;
            test3.setStrike(strikecall3);
            test3.setMaturity(maturity);
            double price3 =0;
            if (numbermethod12==1){
                price3= test3.priceclassic(test);}
            if (numbermethod12==2){
                price3= test3.priceMonteCarlo(test);}
            std :: cout <<" The price of the Stock is "<<test.stockPrice<< std::endl;
            std :: cout <<" The Volatility of the Stock is "<<test.volatility * 100<<" %"<<std::endl;
            std :: cout <<" The Risk-free rate is "<<test.riskFreeRate * 100<<" %"<< std::endl;
            std :: cout <<" The Date is "<<test.date<< std::endl;
            std :: cout <<" The Drift is "<<test.drift*100<<" %"<< std::endl;
            std :: cout <<" The Maturity of this Strategy is "<<maturity<< std::endl;
            std :: cout <<" The Strike of the 1st Call you'll buy in this Strategy is "<<strikecall<< std::endl;
            std :: cout <<" The Strike of the 2 Calls you'll sell in this Strategy is "<<strikecall2<< std::endl;
            std :: cout <<" The Strike of the 2nd Call you'll buy in this Strategy is "<<strikecall3<< std::endl;
            std :: cout <<" The price of the 1st Call you'll buy is "<<price<< std::endl;
            std :: cout <<" The price of the 2 Calls you'll sell is (Price of 1 Call) "<<price2<< std::endl;
            std :: cout <<" The price of the 2nd Call you'll buy is "<<price3<< std::endl;
            std :: cout <<" The price of this Strategy is (you earn the money of the calls you'll sell and you spend/'lose' the money of the call you'll buy) : "<<price+price3-2*price2<< std::endl;
            if (price+price3-2*price2>=0)
            {
            std :: cout <<" You have to spend this money to use this strategy "<< std::endl;
            }
            if (price+price3-2*price2<0)
            {
            std :: cout <<" You earn this money (the absolute value of the result just above) to use this strategy "<<std::endl;
            }
            std :: cout <<" Thank you for using this program "<<std::endl;
            }
        else if(number==6){
            std :: cout << "Choose your Maturity (only in interger numbers) in years"<< std::endl;
            std:: cin >> maturity;
            std :: cout << "Choose your Strike for the first Call you will sell, it has to be the smaller"<< std::endl;
            std:: cin >> strikecall;
            std :: cout << "Choose your Strike for the 1st Call you'll buy"<< std::endl;
            std:: cin >> strikecall2;
            std :: cout << "Choose your Strike for the 2nd Call you'll buy"<< std::endl;
            std:: cin >> strikecall3;
            std :: cout << "Choose your Strike for the 2nd Call you'll sell"<< std::endl;
            std:: cin >> strikecall4;
            EuropeanVanillaCallOption test1;
            test1.setStrike(strikecall);
            test1.setMaturity(maturity);
            double price =0;
            std :: cout << "Choose your Method: 1. Classic/Explicit"<< std::endl;
            std :: cout << "2. Monte Carlo"<< std::endl;
            std:: cin >> numbermethod12;
            if (numbermethod12==1){
                price= test1.priceclassic(test);}
            if (numbermethod12==2){
                price= test1.priceMonteCarlo(test);}
            EuropeanVanillaCallOption test2;
            test2.setStrike(strikecall2);
            test2.setMaturity(maturity);
            double price2 =0;
            if (numbermethod12==1){
                price2= test2.priceclassic(test);}
            if (numbermethod12==2){
                price2= test2.priceMonteCarlo(test);}
            EuropeanVanillaCallOption test3;
            test3.setStrike(strikecall3);
            test3.setMaturity(maturity);
            double price3 =0;
            if (numbermethod12==1){
                price3= test3.priceclassic(test);}
            if (numbermethod12==2){
                price3= test3.priceMonteCarlo(test);}
            EuropeanVanillaCallOption test4;
            test4.setStrike(strikecall4);
            test4.setMaturity(maturity);
            double price4 =0;
            if (numbermethod12==1){
                price4= test4.priceclassic(test);}
            if (numbermethod12==2){
                price4= test4.priceMonteCarlo(test);}
            std :: cout <<" The price of the Stock is "<<test.stockPrice<< std::endl;
            std :: cout <<" The Volatility of the Stock is "<<test.volatility * 100<<" %"<<std::endl;
            std :: cout <<" The Risk-free rate is "<<test.riskFreeRate * 100<<" %"<< std::endl;
            std :: cout <<" The Date is "<<test.date<< std::endl;
            std :: cout <<" The Drift is "<<test.drift*100<<" %"<< std::endl;
            std :: cout <<" The Maturity of this Strategy is "<<maturity<< std::endl;
            std :: cout <<" The Strike of the 1st Call you'll sell in this Strategy is "<<strikecall<< std::endl;
            std :: cout <<" The Strike of the 1st Call you'll buy in this Strategy is "<<strikecall2<< std::endl;
            std :: cout <<" The Strike of the 2nd Call you'll buy in this Strategy is "<<strikecall3<< std::endl;
            std :: cout <<" The Strike of the 2nd Call you'll sell in this Strategy is "<<strikecall4<< std::endl;
            std :: cout <<" The price of the 1st Call you'll sell is "<<price<< std::endl;
            std :: cout <<" The price of the 1st Call you'll buy is "<<price2<< std::endl;
            std :: cout <<" The price of the 2nd Call you'll buy is "<<price3<< std::endl;
            std :: cout <<" The price of the 2nd Call you'll sell is "<<price4<< std::endl;
            std :: cout <<" The price of this Strategy is (you earn the money of the calls you'll sell and you spend/'lose' the money of the call you'll buy) : "<<price2+price3-price-price4<< std::endl;
            if (price2+price3-price-price4>=0)
            {
            std :: cout <<" You have to spend this money to use this strategy "<< std::endl;
            }
            if (price2+price3-price-price4<0)
            {
            std :: cout <<" You earn this money (the absolute value of the result just above) to use this strategy "<<std::endl;
            }
            std :: cout <<" Thank you for using this program "<<std::endl;
            }
    }
    return 0;
    
}


