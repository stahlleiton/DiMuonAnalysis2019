/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

 /** \class RooModExtCBShape
     \ingroup Roofit
    
     P.d.f implementing the Modified Extended Crystal Ball line shape
 **/

#include "Riostream.h"

#include "RooModExtCBShape.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooMath.h"

ClassImp(RooModExtCBShape);

////////////////////////////////////////////////////////////////////////////////
    
Double_t RooModExtCBShape::ApproxErf(Double_t arg) const
{
  static const double erflim = 5.0;
  if( arg > erflim )
    return 1.0;
  if( arg < -erflim )
    return -1.0;
    
  return RooMath::erf(arg);
}
    
////////////////////////////////////////////////////////////////////////////////
    
RooModExtCBShape::RooModExtCBShape(const char *name, const char *title,
                                   RooAbsReal& _m, RooAbsReal& _m0, RooAbsReal& _sigma,
                                   RooAbsReal& _alpha, RooAbsReal& _n, RooAbsReal& _alpha2) :
RooAbsPdf(name, title),
  m("m", "Dependent", this, _m),
  m0("m0", "M0", this, _m0),
  sigma("sigma", "Sigma", this, _sigma),
  alpha("alpha", "Alpha", this, _alpha),
  n("n", "Order", this, _n),
  alpha2("alpha2", "Alpha2", this, _alpha2)
{
}

RooModExtCBShape::RooModExtCBShape()
{
   RooRealVar mRV("mRV","x", 0.0, 1.0);
   RooRealVar m0RV("m0RV", "m0", 0.0, 1.0);
   RooRealVar sigmaRV("sigmaRV", "sigma", 0.0, 1.0);
   RooRealVar alphaRV("alphaRV", "alpha", 0.0, 1.0);
   RooRealVar nRV("nRV","n", 0.0, 1.0);
   RooRealVar alpha2RV("alpha2RV", "alpha2", 0.0, 1.0);
   RooModExtCBShape("RooModExtCBShape", "RooModExtCBShape", mRV, m0RV, sigmaRV, alphaRV, nRV, alpha2RV);
}
    
////////////////////////////////////////////////////////////////////////////////
    
RooModExtCBShape::RooModExtCBShape(const RooModExtCBShape& other, const char* name) :
  RooAbsPdf(other, name), m("m", this, other.m), m0("m0", this, other.m0),
  sigma("sigma", this, other.sigma), alpha("alpha", this, other.alpha),
  n("n", this, other.n), alpha2("alpha2", this, other.alpha2)
{
}
    
////////////////////////////////////////////////////////////////////////////////
    
Double_t RooModExtCBShape::evaluate() const {
  Double_t t = (m-m0)/sigma;
  if (alpha < 0) t = -t;
  
  Double_t absAlpha = fabs((Double_t)alpha);
  Double_t absAlpha2 = fabs((Double_t)alpha2);
  
  if (t >= -absAlpha && t < absAlpha2) {
    return exp(-0.5*t*t);
  }
  if (t < -absAlpha) {
    Double_t a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
    Double_t b = n/absAlpha - absAlpha;
    return a/TMath::Power(b - t, n);
  }
  if (t >= absAlpha2) {
    return exp(0.5*absAlpha2*absAlpha2 - absAlpha2*t);
  }
  return 0.0;
}
    
////////////////////////////////////////////////////////////////////////////////
    
Int_t RooModExtCBShape::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  if( matchArgs(allVars,analVars,m) )
    return 1 ;
    
  return 0;
}
    
////////////////////////////////////////////////////////////////////////////////
    
Double_t RooModExtCBShape::analyticalIntegral(Int_t code, const char* rangeName) const
{
  static const double sqrtPiOver2 = 1.2533141373;
  static const double sqrt2 = 1.4142135624;
    
  R__ASSERT(code==1);
  double result = 0.0;
  bool useLog = false;
    
  if( fabs(n-1.0) < 1.0e-05 ) useLog = true;
    
  double sig = fabs((Double_t)sigma);
    
  double tmin = (m.min(rangeName)-m0)/sig;
  double tmax = (m.max(rangeName)-m0)/sig;
    
  if(alpha < 0) {
    double tmp = tmin;
    tmin = -tmax;
    tmax = -tmp;
  }
   
  double absAlpha = fabs((Double_t)alpha);
  double absAlpha2 = fabs((Double_t)alpha2);
  
  if( tmin >= -absAlpha && tmax < absAlpha2) {
    result += sig*sqrtPiOver2*(   ApproxErf(tmax/sqrt2)
                                  - ApproxErf(tmin/sqrt2) );
  }
  else if( tmax < -absAlpha ) {
    double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
    double b = n/absAlpha - absAlpha;
    
    if(useLog) {
      result += a*sig*( log(b-tmin) - log(b-tmax) );
    }
    else {
      result += a*sig/(1.0-n)*(   1.0/(TMath::Power(b-tmin,n-1.0))
                                  - 1.0/(TMath::Power(b-tmax,n-1.0)) );
    }
  }
  else if( tmin >= absAlpha2 ) {
    result += (sig/absAlpha2)*exp(0.5*absAlpha2*absAlpha2)*( exp(-absAlpha2*tmin) - exp(-absAlpha2*tmax) );
  }
  else if( tmax < absAlpha2 ) {
    double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
    double b = n/absAlpha - absAlpha;
    
    double term1 = 0.0;
    if(useLog) {
      term1 = a*sig*(  log(b-tmin) - log(n/absAlpha));
    }
    else {
      term1 = a*sig/(1.0-n)*(   1.0/(TMath::Power(b-tmin,n-1.0))
                                - 1.0/(TMath::Power(n/absAlpha,n-1.0)) );
    }
    
    double term2 = sig*sqrtPiOver2*(   ApproxErf(tmax/sqrt2)
                                       - ApproxErf(-absAlpha/sqrt2) );
    
    result += term1 + term2;
  }
  else if( tmin >= -absAlpha ) {
    double term1 = sig*sqrtPiOver2*( ApproxErf(absAlpha2/sqrt2) - ApproxErf(tmin/sqrt2) );
    double term2 = (sig/absAlpha2)*exp(0.5*absAlpha2*absAlpha2)*( exp(-absAlpha2*absAlpha2) - exp(-absAlpha2*tmax) );
    result += term1 + term2;
  }
  else {
    double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
    double b = n/absAlpha - absAlpha;
    
    double term1 = 0.0;
    if(useLog) {
      term1 = a*sig*(  log(b-tmin) - log(n/absAlpha));
    }
    else {
      term1 = a*sig/(1.0-n)*(   1.0/(TMath::Power(b-tmin,n-1.0))
                                - 1.0/(TMath::Power(n/absAlpha,n-1.0)) );
    }
    
    double term2 = sig*sqrtPiOver2*(   ApproxErf(absAlpha2/sqrt2)
                                       - ApproxErf(-absAlpha/sqrt2) );

    double term3 = (sig/absAlpha2)*exp(0.5*absAlpha2*absAlpha2)*( exp(-absAlpha2*absAlpha2) - exp(-absAlpha2*tmax) );
    
    result += term1 + term2 + term3;
  }
    
  return result;
}
    
////////////////////////////////////////////////////////////////////////////////
/// Advertise that we know the maximum of self for given (m0,alpha,n,sigma)
    
Int_t RooModExtCBShape::getMaxVal(const RooArgSet& vars) const
{
  RooArgSet dummy ;
    
  if (matchArgs(vars,dummy,m)) {
    return 1 ;
  }
  return 0 ;
}
    
////////////////////////////////////////////////////////////////////////////////
    
Double_t RooModExtCBShape::maxVal(Int_t code) const
{
  R__ASSERT(code==1) ;
    
  // The maximum value for given (m0,alpha,n,alpha2,sigma)
  return 1.0 ;
}
