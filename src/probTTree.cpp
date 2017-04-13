/* Rewrite probTTree with allowTransPostSamp = TRUE
   Calculates the log-probability of a transmission tree.
   @param ttree Transmission tree
   @param off.r First parameter of the negative binomial distribution for offspring number
   @param off.p Second parameter of the negative binomial distribution for offspring number
   @param pi probability of sampling an infected individual
   @param w.shape Shape parameter of the Gamma probability density function representing the generation time
   @param w.scale Scale parameter of the Gamma probability density function representing the generation time 
   @param ws.shape Shape parameter of the Gamma probability density function representing the sampling time
   @param ws.scale Scale parameter of the Gamma probability density function representing the sampling time 
   @param dateT Date when process stops (this can be Inf for fully simulated outbreaks)
   @param allowTransPostSamp Whether or not to allow transmission after sampling of a host
   @return Probability of the transmission tree */


// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <boost/math/tools/roots.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <limits>
using namespace Rcpp;

struct wstar_functor
{ // Functor also returning 1st derivative.
  wstar_functor(double const& pi, double const& p, double const& r) : pi(pi), p(p), r(r){};
  
  std::pair<double, double> operator()(double const& x)
  {
    // Return both f(x) and f'(x).
    double temp = pow((1-p)/(1-p*x),r);
    double fx = x - (1-pi)*temp; 
    double dx = 1 - (1-pi)*p*r/(1-p*x)*temp;
    return std::make_pair(fx, dx);   
  }
private:
  double pi;
  double p;
  double r;
};


double wstar_rootFinder(double pi, double p, double r)
{
  using namespace boost::math::tools;
  
  const int digits = std::numeric_limits<double>::digits;  // Maximum possible binary digits accuracy for type T.
  int get_digits = static_cast<int>(digits * 0.6);    // Accuracy doubles with each step, so stop when we have
                                                      // just over half the digits correct.
  const boost::uintmax_t maxit = 20;
  boost::uintmax_t it = maxit;
  double result = newton_raphson_iterate(wstar_functor(pi, p, r), 0.5, 0.0, 1.0, get_digits, it);
  return result;
}


double alphastar(int d, double p, double r, double wstar)
{
  if(std::abs(r-1.0)<1e-6) // Exact solution available
    return (1-p)/(1-p*wstar)*pow(p/(1-p*wstar), d);

  int k = d;
  std::vector<double> toSum;
  
  boost::math::negative_binomial_distribution<double> nbinom(r,p);
  while(true){
    
    double dnb = pdf(nbinom,k);
    double term = dnb * pow(wstar,k);

    toSum.push_back(term);
    if(term < 1e-8) break; // Achieved desired accuracy
    
    k++;
  }

  NumericVector toSumR = wrap(toSum); // Convert toSum to Rcpp NumericVector
  NumericVector v(k-d+1);
  for(int i=0; i<v.size(); i++) v[i] = i+d;

  return sum(choose(v, d)*toSumR)/pow(wstar,d);
}
    

/* alpha() computes equation (10) in TransPhylo paper
   wbar0 --- wbar computed from wbar() using the oldest infection time. */
double alpha(double tinf, int d, double p, double r, NumericVector wbar0, double gridStart, double delta_t)
{
  double wbar_tinf = wbar0[std::round((tinf - gridStart)/delta_t)];
  if(std::abs(r-1.0)<1e-6) // Exact solution available
    return (1-p)/(1-p*wbar_tinf)*pow(p/(1-p*wbar_tinf), d);


  int k = d;
  std::vector<double> toSum;
  
  while(true){
    
    double dnb = R::dnbinom(k,r,p,0);
    double term = dnb * pow(wbar_tinf,k);

    toSum.push_back(term);
    if(Rf_choose(k,d)*term < 1e-8) break; // Achieved desired accuracy
    
    k++;
  }

  NumericVector toSumR = wrap(toSum); // Convert toSum to Rcpp NumericVector
  NumericVector v(k-d+1);
  for(int i=0; i<v.size(); i++) v[i] = i+d;

  return sum(choose(v, d)*toSumR)/pow(wbar_tinf,d);
}

  

// [[Rcpp::export]]
NumericVector wbar(double tinf, double dateT, double rOff, double pOff, double pi, double shGen, double scGen, double shSam, double scSam, double delta_t=0.05)
{
  int n = std::round((dateT-tinf)/delta_t); 
  NumericVector grid(n);
  for(int i=0; i<n; ++i) // use the left point of each subinterval
    grid[i] = dateT-n*delta_t+i*delta_t;

  NumericVector pi2 = pi*pgamma(dateT-grid, shSam, scSam);
  NumericVector F = 1-pgamma(dateT-grid, shGen, scGen);

  NumericVector w(n+1), out(n+1);
  out[n] = w[n] = 1.0;

  IntegerVector seq = seq_len(n);
  NumericVector gam = dgamma(as<NumericVector>(seq)*delta_t,shGen,scGen);
  double sumPrev = 0.5 * gam[0];
  for(int i=n-1; i>=0; --i){

    w[i] = (1-pi2[i]) * pow((1-pOff)/(1-pOff*F[i]-pOff*delta_t*sumPrev), rOff);
    out[i] = F[i] + sumPrev*delta_t;
    
    sumPrev = 0.0;
    for(int j=0; j<n-i; ++j)
      sumPrev += gam[j]*w[i+j];
    sumPrev += 0.5 * gam[n-i];
  }
  return out;
}
      


// [[Rcpp::export]]
double probTTree(NumericMatrix ttree, double rOff, double pOff, double pi, double shGen, double scGen, double shSam, double scSam, double dateT, double delta_t=0.05){

  int numCases = ttree.nrow();
  boost::math::gamma_distribution<double> genGamma(shGen, scGen);
  
  if(dateT == INFINITY){ // finished outbreak
    double wstar = wstar_rootFinder(pi, pOff, rOff);

    NumericVector sstatus = ifelse(is_na(ttree(_,1)), 1-pi, pi*dgamma(ttree(_,1)-ttree(_,0),shSam,scSam));
    NumericVector lsstatus = log(sstatus);

    std::map<int, std::vector<int> > infMap; // Map from infector to infected
    std::vector<std::vector<int> > progeny(numCases);
    for(int i=0; i<numCases; ++i){
      if(ttree(i,2) == 0) continue; // Found root node i 

      progeny[ttree(i,2)-1].push_back(i); // C++ index starts from 0
      infMap[ttree(i,2)-1] = progeny[ttree(i,2)-1]; 
    }

    double accum = 0.0;
    for(int i=0; i<numCases; ++i){
      accum += log(alphastar(progeny[i].size(), pOff, rOff, wstar));
     
      for(int j=0; j<progeny[i].size(); ++j)
	accum += log(pdf(genGamma, ttree(progeny[i][j],0) - ttree(i,0)));
    }
      
    return sum(lsstatus) + accum;
  }
  else{
    // Ongoing outbreak -- observation ends at finite dateT
    NumericVector probSam = pi*pgamma(dateT-ttree(_,0),shSam,scSam);
    NumericVector sstatus = ifelse(is_na(ttree(_,1)), 1-probSam, pi*dgamma(ttree(_,1)-ttree(_,0),shSam,scSam));
    NumericVector lsstatus = log(sstatus);

    std::map<int, std::vector<int> > infMap; // Map from infector to infected
    std::vector<std::vector<int> > progeny(numCases);
    for(int i=0; i<numCases; ++i){
      if(ttree(i,2) == 0) continue; // Found root node i 

      progeny[ttree(i,2)-1].push_back(i); // C++ index starts from 0
      infMap[ttree(i,2)-1] = progeny[ttree(i,2)-1]; 
    }

    double accum = 0.0;
    double tinfmin = min(ttree(_,0));
    NumericVector wbar0 = wbar(tinfmin, dateT, rOff, pOff, pi, shGen, scGen, shSam, scSam, delta_t);
    // wbar0.size = grid size +1 
    double gridStart = dateT-(wbar0.size()-1)*delta_t;

    for(int i=0; i<numCases; ++i){

      accum += log(alpha(ttree(i,0), progeny[i].size(), pOff, rOff, wbar0, gridStart, delta_t));

      for(int j=0; j<progeny[i].size(); ++j)
	accum += (R::dgamma(ttree(progeny[i][j],0)-ttree(i,0), shGen, scGen, 1) - R::pgamma(dateT-ttree(i,0), shGen, scGen, 1, 1));
    }

    return sum(lsstatus) + accum;
  }
}

