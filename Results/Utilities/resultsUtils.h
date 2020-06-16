#ifndef Utilities_resultUtils_h
#define Utilities_resultUtils_h
// Auxiliary Headers
#include "bins.h"
#include "../../Utilities/dataUtils.h"
#include "../../Utilities/CMS/tdrstyle.C"
#include "../../Utilities/CMS/CMS_lumi.C"
// ROOT headers
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TPaletteAxis.h"
#include "TF1.h"
#include "TMatrixD.h"
// RooFit headers
// c++ headers
#include <iostream>
#include <string>
#include <map>
// CMS headers


// ------------------ TYPE -------------------------------
typedef std::vector< double                          > DoubleVec_t;
typedef std::pair< double        , double            > DoublePair_t;
typedef std::pair< AnaBin_t        , int               > BinPair_t;
typedef std::map< AnaBin_t         , DoubleDiMap_t     > DoubleBinDiMap_t;
typedef std::map< std::string    , DoubleBinDiMap_t  > VarBinMap_t;
typedef std::map< std::string    , VarBinMap_t       > VarBinDiMap_t;
typedef std::map< std::string    , VarBinDiMap_t     > VarBinTriMap_t;
typedef std::map< std::string    , VarBinTriMap_t    > VarBinQuadMap_t;
typedef std::map< AnaBin_t         , double            > BinMap_t;
typedef std::map< std::string    , BinMap_t          > BinDiMap_t;
typedef std::map< std::string    , BinDiMap_t        > BinTriMap_t;
typedef std::map< std::string    , BinTriMap_t       > BinQuadMap_t;
typedef std::map< std::string    , BinQuadMap_t      > BinPentaMap_t;
typedef std::map< std::string    , BinPentaMap_t     > BinSextaMap_t;
typedef std::vector< BinSextaMap_t                   > BinSextaMapVec_t;
typedef std::map< std::string    , BinSextaMapVec_t  > BinSeptaMapVec_t;
typedef std::map< std::string    , BinSeptaMapVec_t  > BinOctaMapVec_t;
typedef std::map< std::string    , BinSextaMap_t     > BinSeptaMap_t;
typedef std::map< std::string    , BinSeptaMap_t     > BinOctaMap_t;
typedef std::map< BinPair_t      , std::set<AnaBin_t>  > BinSetMap_t;
typedef std::map< AnaBin_t         , BinSetMap_t       > BinSetDiMap_t;
typedef std::map< std::string    , BinSetDiMap_t     > BinSetTriMap_t;
typedef std::map< std::string    , BinSetTriMap_t    > BinSetQuadMap_t;
typedef std::map< std::string    , BinSetQuadMap_t   > BinSetPentaMap_t;
typedef std::map< std::string    , BinSetPentaMap_t  > BinCont_t;
typedef std::map< std::string    , TGraphAsymmErrors > GraphMap_t;
typedef std::map< BinPair_t      , GraphMap_t        > GraphDiMap_t;
typedef std::map< std::string    , GraphDiMap_t      > GraphTriMap_t;
typedef std::map< AnaBin_t         , GraphTriMap_t     > GraphQuadMap_t;
typedef std::map< std::string    , GraphQuadMap_t    > GraphPentaMap_t;
typedef std::map< std::string    , GraphPentaMap_t   > GraphSextaMap_t;
typedef std::map< std::string    , GraphSextaMap_t   > GraphSeptaMap_t;
typedef std::pair< StringVector_t, IntMap_t          > StringVecIntPair_t;
typedef std::map< std::string    , StringVecIntPair_t > StringVecIntPairMap_t;
typedef std::map< std::string , StringVecIntPairMap_t > WSDirMap_t;

// TO REMOVE
typedef std::pair< AnaBin_t        , AnaBin_t            > BinPair2_t;
typedef std::map< BinPair2_t      , std::set<AnaBin_t>  > BinSetMap2_t;
typedef std::map< std::string    , BinSetMap2_t       > BinSetDiMap2_t;
typedef std::map< std::string    , BinSetDiMap2_t     > BinSetTriMap2_t;
typedef std::map< std::string    , BinSetTriMap2_t    > BinSetQuadMap2_t;
typedef std::map< std::string    , BinSetQuadMap2_t   > BinCont2_t;
typedef std::map< std::string    , TGraphAsymmErrors > GraphMap2_t;
typedef std::map< std::string      , GraphMap2_t        > GraphDiMap2_t;
typedef std::map< BinPair2_t      , GraphDiMap2_t       > GraphTriMap2_t;
typedef std::map< std::string    , GraphTriMap2_t     > GraphQuadMap2_t;
typedef std::map< std::string    , GraphQuadMap2_t    > GraphPentaMap2_t;
typedef std::map< std::string    , GraphPentaMap2_t   > GraphSextaMap2_t;


// ------------------ GLOBAL -------------------------------
const std::vector<int> COLOR({kBlack, kRed, kGreen+2, kBlue+2, kAzure-7, kYellow-3, kMagenta+3, kCyan+3, kOrange+9, kSpring+4, kTeal+4, kViolet+9, kPink+8});
const StringVector_t   STAT({ "Err_Stat_High" , "Err_Stat_Low" });
const StringVector_t   SYST({ "Err_Syst_High" , "Err_Syst_Low" });


// ------------------ FUNCTION -------------------------------
double symError(const double& errLow, const double& errHigh)
{
  // check inputs
  if (errLow<0. || isnan(errLow)) { throw std::logic_error(Form("[ERROR] symErrors: Invalid lower uncertainty ( %.6f )", errLow)); }
  if (errHigh<0. || isnan(errHigh)) { throw std::logic_error(Form("[ERROR] symErrors: Invalid higher uncertainty ( %.6f )", errHigh)); }
  // symmetrize uncertainties
  const auto res = std::max(errLow, errHigh);
  // check result
  if (res<0. || isnan(res)) { throw std::logic_error(Form("[ERROR] Invalid symmetric error %.6f ( %.6f , %.6f )", res, errLow, errHigh)); }
  // return result
  return res;
};


void checkVec(const DoubleVec_t& vec, const bool& testNeg=false)
{
  for (const auto& v : vec) {
    if ((testNeg && v<0.) || isnan(v)) { throw std::logic_error(Form("[ERROR] Invalid value %.2f", v)); }
  }
}

  
double sumErrors(const DoubleVec_t& errV, const TMatrixD& corrM=TMatrixD())
{
  //if (errV==std::vector<double>(errV.size())) { return 0.0; }
  // check inputs
  checkVec(errV);
  // compute quadratic sum of uncertainties
  double sumErrSqr = 0.0;
  if (corrM.GetNoElements()>0) {
    if (corrM.GetNcols()!=int(errV.size()) || !corrM.IsSymmetric()) { corrM.Print(); throw std::logic_error(Form("[ERROR] corrM(%d, %d) for errV(%lu) is not valid", corrM.GetNrows(), corrM.GetNcols(), errV.size())); }
    for (size_t i1=0; i1<errV.size(); i1++) {
      for (size_t i2=0; i2<errV.size(); i2++) {
	sumErrSqr += errV[i1]*errV[i2]*corrM(i1,i2);
      }
    }
  }
  else {
    for (const auto& e : errV) {
      sumErrSqr += e*e;
    }
  }
  const auto res = std::sqrt(sumErrSqr);
  // check result
  const bool isUnCorr = (corrM.GetNoElements()==0 || corrM.Determinant()==1.);
  if ((isUnCorr && res==0.) || res<0. || isnan(res)) { throw std::logic_error(Form("[ERROR] Invalid sum of errors %.2f ( sqrt(%.2f) )", res, sumErrSqr)); }
  // return result
  return res;
};


double divide(const DoubleVec_t& numV, const DoubleVec_t& denV)
{
  // check inputs
  checkVec(numV);
  checkVec(denV);
  // compute division
  double num=1.0, den=1.0;
  for (const auto& n : numV) { num *= n; }
  for (const auto& n : denV) { den *= n; }
  //if (den==0.0) { den = 1.0E-10; }
  const auto res = (num / den);
  // check division
  if (isnan(res)) { throw std::logic_error(Form("[ERROR] Invalid divide arguments %.2f ( %.2f / %.2f )", res, num, den)); }
  // return division
  return res;
};


double divide(const double& num, const double& den)
{
  return divide(DoubleVec_t({num}), DoubleVec_t({den}));
};

  
double divideError(const DoubleVec_t& numV, const DoubleVec_t& denV,
		   const DoubleVec_t& numErrV, const DoubleVec_t& denErrV,
		   const TMatrixD& corrIM=TMatrixD())
{
  // check inputs
  checkVec(numV);
  checkVec(denV, true);
  checkVec(numErrV, true);
  checkVec(denErrV, true);
  // check input vectors
  if (numErrV.size()!=numV.size()) { throw std::logic_error(Form("[ERROR] numErrV.size(%lu) != numV.size(%lu)", numErrV.size(), numV.size())); }
  if (denErrV.size()!=denV.size()) { throw std::logic_error(Form("[ERROR] denErrV.size(%lu) != denV.size(%lu)", denErrV.size(), denV.size())); }
  // compute relative errors
  DoubleVec_t relErrV;
  for (size_t i=0; i<numV.size(); i++) { relErrV.push_back(divide(numErrV[i], numV[i])); }
  for (size_t i=0; i<denV.size(); i++) { relErrV.push_back(divide(denErrV[i], denV[i])); }
  // compute correlation matrix for division
  auto corrM = corrIM;
  if (corrM.GetNoElements()>0) {
    if (corrM.GetNcols()!=int(relErrV.size()) || !corrM.IsSymmetric()) { throw std::logic_error(Form("[ERROR] corrM(%d, %d) for relErrV(%lu) is not valid", corrM.GetNrows(), corrM.GetNcols(), relErrV.size())); }
    for (int i1=0; i1<corrM.GetNrows(); i1++) {
      for (int i2=0; i2<corrM.GetNcols(); i2++) {
	corrM(i1, i2) = (i1<int(numV.size()) ? 1 : -1)*(i2<int(numV.size()) ? 1 : -1)*corrM(i1, i2);
      }
    }
  }
  // compute total uncertainty
  const auto res = (std::abs(divide(numV, denV)) * sumErrors(relErrV, corrM));
  // check total uncertainty
  const bool isUnCorr = (corrM.GetNoElements()==0 || corrM.Determinant()==1.);
  if ((isUnCorr && res==0.) || res<0. || isnan(res)) { throw std::logic_error(Form("[ERROR] Invalid divide uncertainty %.2f", res)); }
  // return total uncertainty
  return res;
};

  
double divideError(const DoubleVec_t& numV, const DoubleVec_t& denV,
		   const DoubleVec_t& numErrV, const DoubleVec_t& denErrV,
		   const bool& isCorrelated)
{
  const auto n = numV.size() + denV.size();
  const auto corrM = (isCorrelated ? TMatrixD(n,n,DoubleVec_t(n*n,1).data()) : TMatrixD());
  return divideError(numV, denV, numErrV, denErrV, corrM);
};

  
double divideError(const double& num, const double& den, const double& numErr, const double& denErr, const TMatrixD& corrM=TMatrixD())
{
  return divideError(DoubleVec_t({num}), DoubleVec_t({den}), DoubleVec_t({numErr}), DoubleVec_t({denErr}), corrM);
};

  
double divideError(const double& num, const double& den, const double& numErr, const double& denErr, const bool& isCorrelated)
{
  return divideError(DoubleVec_t({num}), DoubleVec_t({den}), DoubleVec_t({numErr}), DoubleVec_t({denErr}), isCorrelated);
};



void getBins(BinSetMap2_t& binS, const std::set<AnaBin_t>& allBin, const std::string& obs, const StringVector_t& vars)
{
  if (contain(vars, obs)) return;
  std::set<BinPair2_t> bSet;
  for (const auto& b : allBin) {
    if (b.getbin(obs).name()=="") continue;
    AnaBin_t varBin, incBin;
    bool keepBin = true;
    for (const auto& p : b) {
      if (p.name()==obs) continue;
      if (contain(vars, p.name())) { varBin.setbin(p); } else { incBin.setbin(p); }
      bool isInc = false;
      if      (p.name()=="Cand_Rap"   ) { isInc = (isEqual(p.high(),   2.4, 1) && isEqual(p.low(),  -2.4, 1)); }
      else if (p.name()=="Cand_AbsRap") { isInc = (isEqual(p.high(),   2.4, 1) && isEqual(p.low(),   0.0, 1)); }
      else if (p.name()=="Cand_RapCM" ) { isInc = (isEqual(p.high(),  1.93, 2) && isEqual(p.low(), -2.86, 2)); }
      else if (p.name()=="Centrality" ) { isInc = (isEqual(p.high(), 100.0, 1) && isEqual(p.low(),   0.0, 1)); }
      keepBin = (keepBin && (contain(vars, p.name()) ? !isInc : isInc));
    }
    const auto& binP = std::make_pair(varBin, incBin);
    if (keepBin) { binS[binP].insert(b); bSet.insert(binP); }
  }
  // Check if we only have one main bin
  for (const auto& binP : bSet) {
    std::set<BinF_t> mainBin; for (const auto& b : binS.at(binP)) { mainBin.insert(b.getbin(obs)); }
    if (mainBin.size()<=2) { binS.erase(binP); }
  }
};


void combineBins(std::vector<StringVector_t>& cmbObs, const StringVector_t& allObs, const uint& K=0)
{
  const auto& N = allObs.size(); if (K>=N) return;
  std::string bitmask(K, 1); bitmask.resize(N, 0);
  do {
    StringVector_t obsV;
    for (uint i=0; i<N; ++i) { if (bitmask[i]) obsV.push_back(allObs[i]); }
    cmbObs.push_back(obsV);
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
  combineBins(cmbObs, allObs, K+1);
};


void defineGeneralBins(BinCont2_t& binMap, const BinSextaMap_t& var)
{
  for (const auto& o : var) {
    for (const auto& c : o.second) {
      for (const auto& vv : c.second.begin()->second) {
	// Get all bins
	std::set<AnaBin_t> allBin;
	for (const auto& pd : c.second) { for (const auto& b : pd.second.at(vv.first).at("Val")) {  allBin.insert(b.first); } }
	// Get all observables
	StringVector_t allObs;
	for (const auto& b : allBin) { for (const auto& v : b) { if (!contain(allObs, v.name())) { allObs.push_back(v.name()); } } }
	// Combine all observables
	std::vector<StringVector_t> cmbObs;
	combineBins(cmbObs, allObs);
	// Loop over each observable
	for (const auto& p : allObs) {
	  auto& binS = binMap[o.first][c.first][vv.first][p];
	  // Loop over each combination of observables
	  for (const auto& v : cmbObs) { getBins(binS, allBin, p, v); }
	}
      }
    }
  }
};


std::string obsUnit(const std::string& obs)
{
  const auto& obsL = StringMap_t({{"Cand_Mass", "GeV/c^{2}"}, {"Cand_Pt", "GeV/c"}, {"Centrality", "%"}});
  return (contain(obsL, obs) ? obsL.at(obs) : "");
};


std::string formatObsTag(const std::string& obs)
{
  const auto& obsL = StringMap_t({{"Cand_Mass", "M"}, {"Cand_Pt", "p_{T}"}, {"Cand_Rap", "y"}, {"Cand_RapCM", "y_{CM}"},
				  {"Cand_AbsRap", "|y|"}, {"Centrality", "Cent"}, {"NTrack", "Ntrk"}});
  return (contain(obsL, obs) ? obsL.at(obs) : obs);
};


std::string formatObsName(const std::string& obs)
{
  const auto& obsL = StringMap_t({{"Cand_Mass", "Mass [GeV/c^{2}]"}, {"Cand_Pt", "p_{T} [GeV/c]"}, {"Cand_Rap", "y"}, {"Cand_RapCM", "y_{CM}"},
				  {"Cand_AbsRap", "|y|"}, {"Centrality", "Centrality [%]"}, {"NTrack", "Track multiplicity"}});
  std::string obsF = (obs.rfind("Cand_",0)==0 ? "#mu^{+}#mu^{#font[122]{\55}} " : "");
  obsF += (contain(obsL, obs) ? obsL.at(obs) : "");
  return obsF;
};


std::string formatObsRange(const BinF_t& bin)
{
  const auto& bName = bin.name();
  const auto& bN = formatObsTag(bName);
  const auto& bU = obsUnit(bName);
  const auto& obsMin = DoubleMap_t({{"Cand_Mass",   0.0}, {"Cand_Pt",   0.0}, {"Cand_Rap", -2.5}, {"Cand_RapCM", -2.86}, {"Cand_AbsRap", 0.0}, {"Centrality",   0.0}, {"NTrack",   0.0}});
  const auto& obsMax = DoubleMap_t({{"Cand_Mass", 500.0}, {"Cand_Pt", 500.0}, {"Cand_Rap",  2.5}, {"Cand_RapCM",  1.93}, {"Cand_AbsRap", 2.5}, {"Centrality", 100.0}, {"NTrack", 500.0}});
  const auto& min = (contain(obsMin, bName) ? obsMin.at(bName) : -999.0);
  const auto& max = (contain(obsMax, bName) ? obsMax.at(bName) :  999.0);
  std::string label = "";
  if (bin.low()>min && bin.high()>=max) { label = Form("%g %s #leq %s", bin.low(), bU.c_str(), bN.c_str()); }
  else if (bin.low()<=min && bin.high()<max) { label = Form("%s < %g %s", bN.c_str(), bin.high(), bU.c_str()); }
  else if (bin.low()>min && bin.high()<max) { label = Form("%g #leq %s < %g %s", bin.low(), bN.c_str(), bin.high(), bU.c_str()); }
  return label;
};


std::string formatObjName(const std::string& obj)
{
  const auto& objL = StringMap_t({{"DY", "Z/#gamma*"}, {"TTbar", "t#bar{t}"}, {"JPsi", "J/#psi"}, {"Psi2S", "#psi(2S)"},
				  {"Ups1S", "#Upsilon(1S)"}, {"Ups2S", "#Upsilon(2S)"}, {"Ups3S", "#Upsilon(3S)"}});
  auto ob = obj;
  if (ob.find("NoPR")!=std::string::npos) { ob.erase(ob.find("NoPR"), 4); }
  else if (ob.find("PR")!=std::string::npos) { ob.erase(ob.find("PR"), 2); }
  for (const auto& o : objL) { stringReplace(ob, o.first, o.second); }
  return ob;
};


std::string formatParName(const std::string& par)
{
  const auto& parL = StringMap_t({{"Alpha", "#alpha_{"}, {"m", "Mean"}, {"N", "Yield"}, {"Lambda", "#lambda_{"}, {"Sigma", "#sigma_{"}, {"rSigma21","#sigma2/#sigma1"}});
  auto parF = par; for (const auto& pL : parL) { if (parF.rfind(pL.first, 0)==0) { stringReplace(parF, pL.first, pL.second); break; } }
  if (parF.rfind("{")!=std::string::npos) { parF += "}"; }
  return parF;
};


std::string formatResultVarName(const std::string& var, const std::string& par, const std::string& obj, const std::string& col, const bool& isRatio=false)
{
  std::string label = "";
  // Format object name
  const auto& objF = formatObjName(obj);
  const std::string& ref = (obj.rfind("Ups",0)==0 ? "#Upsilon(1S)" : (obj.rfind("Psi",0)==0 ? "J/#psi" : ""));
  const auto& ratioLbl = objF+"/"+ref;
  // Add result variable
  if (var=="Cross_Section") {
    const auto& parF = formatObsTag(par);
    const std::string& unit = (col.rfind("PbPb",0)==0 ? "#mub" : (col.rfind("PP",0)==0 ? "pb" : "nb"));
    label = "B #times d#sigma/d"+parF+" ["+unit+"]";
  }
  else if (var=="RatioTo1S" || var=="R") {
    label = "N("+objF+") / N("+ref+")";
  }
  else if (var=="ForwardBackward_Ratio") {
    label = "N(+y_{CM}) / N(#font[122]{\55}y_{CM})";
  }
  // Add fit parameter
  else {
    const auto& varF = formatParName(var);
    label = (isRatio? ratioLbl : objF)+" "+varF;
  }
  return label;
};


std::string formatDecayLabel(const std::string& var, const std::string& obj)
{
  std::string label = formatObjName(obj);
  if (var=="R" || var=="RatioTo1S") { label = (obj.rfind("Ups",0)==0 ? "#Upsilon(nS)" : (obj.rfind("Psi",0)==0 ? "#psi(nS)" : label)); }
  if (obj.find("NoPR")!=std::string::npos) {
    label = "b#rightarrow "+ label +"#rightarrow #mu^{+}#mu^{#font[122]{\55}}";
  }
  else {
    label += "#rightarrow #mu^{+}#mu^{#font[122]{\55}}";
  }
  return label;
};


void setStyle()
{
  // Set the CMS style
  setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TGaxis::SetMaxDigits(3); // to display powers of 10
  //
  // Set Palette
  gStyle->SetPalette(55);
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  //
};


void formatLegendEntry(TLegendEntry& e, const double& size=0.040)
{
  e.SetTextSize(size);
};


void setRangeYAxisGraph(TGraphAsymmErrors& graph, const double& fracDo = 0.1, const double& fracUp = 0.5, const double& yMin=-99., const double& yMax=-99.)
{
  // Find maximum and minimum points of Plot to rescale Y axis
  double YMax = -1e99;
  if (yMax!=-99.) { YMax = yMax; }
  else { for (int i=0; i<=graph.GetN(); i++) { double x, y; graph.GetPoint(i, x, y); y += graph.GetErrorY(i); YMax = std::max(YMax, y); } }
  double YMin = 1e99;
  if (yMin!=-99.) { YMin = yMin; }
  else { for (int i=0; i<=graph.GetN(); i++) { double x, y; graph.GetPoint(i, x, y); y -= graph.GetErrorY(i); YMin = std::min(YMin, y); } }
  //
  const double YTot = (fracDo>=0. ? (YMax - YMin)/(1.0 - fracUp - fracDo) : YMax/(1.0 - fracUp));
  const double Yup   = YMax + fracUp*YTot;
  const double Ydown = (fracDo>=0. ? YMin - fracDo*YTot : 0.0);
  //
  graph.GetYaxis()->SetRangeUser(Ydown, Yup);
};


void formatResultsGraph(TGraphAsymmErrors& graph, const std::string& var, const std::string& par, const std::string& obj, const std::string& col, const int& color=-1, const bool& isRatio=false, const std::vector<double>& pRange={})
{
  //
  // Set the Axis Titles
  std::string xLabel = formatObsName(par);
  std::string yLabel = formatResultVarName(var, par, obj, col, isRatio);
  graph.SetTitle(Form(";%s;%s", xLabel.c_str(), yLabel.c_str()));
  //
  // General
  graph.SetMarkerColor((color>0 ? color : kBlack));
  graph.SetLineColor((color>0 ? color : kBlack)); 
  graph.SetMarkerStyle(20);
  graph.SetMarkerSize(1.5);
  graph.SetLineWidth(3);
  graph.SetLineStyle(1);
  graph.SetFillStyle(0);
  // X-axis
  graph.GetXaxis()->CenterTitle(kTRUE);
  graph.GetXaxis()->SetTitleOffset(0.75);
  graph.GetXaxis()->SetTitleSize(0.065);
  graph.GetXaxis()->SetLabelSize(0.035);
  double xMin = graph.GetXaxis()->GetXmin();
  double xMax = graph.GetXaxis()->GetXmax();
  if (!pRange.empty()) { xMin = pRange[0]; xMax = pRange[2]; }
  xMin -= 0.1; xMax += 0.1;
  if      (par=="Cand_Pt"    ) { xMin = -0.1; xMax = 66.0; }
  else if (par=="Cand_Rap"   ) { xMin = -2.6; xMax = 2.6; }
  else if (par=="Cand_RapCM" ) { xMin = -2.9; xMax = 2.0; }
  else if (par=="Cand_AbsRap") { xMin = -0.1; xMax = 2.6; }
  else if (par=="NTrack") { xMin = 0.0; xMax = 275.0; }
  graph.GetXaxis()->SetRangeUser(xMin, xMax);
  graph.GetXaxis()->SetLimits(xMin, xMax);
  if (par=="NTrack") graph.GetXaxis()->SetNdivisions(510);
  else graph.GetXaxis()->SetNdivisions(505);
  // Y-axis
  graph.GetYaxis()->CenterTitle(kTRUE);
  graph.GetYaxis()->SetTitleOffset(1.05);
  graph.GetYaxis()->SetTitleSize(0.065);
  graph.GetYaxis()->SetLabelSize(0.035);
  double fracD = ((var=="ForwardBackward_Ratio" || var=="RatioTo1S") ? -0.1 : 0.1); 
  if (pRange.empty()) { setRangeYAxisGraph(graph, fracD, 0.4); }
  else { setRangeYAxisGraph(graph, fracD, 0.4, pRange[1], pRange[3]); }
};


std::tuple<double, double, double> computeStats(TGraphAsymmErrors& gr)
{
  double mean = 0;
  double N = 0;
  for (int i=0; i<gr.GetN(); i++) {
    const auto& y = gr.GetY()[i];
    const auto& yErrLo = gr.GetEYlow()[i];
    const auto& yErrHi = gr.GetEYhigh()[i];
    auto yErr = std::max(yErrLo, yErrHi);
    if (yErr==0.) { yErr = 1E-9; }
    const auto yVar = yErr*yErr;
    const auto w = 1./yVar;
    mean += w*y;
    N += w;
  }
  mean /= N;
  const auto& meanErr = 1./std::sqrt(N);
  double RMS = 0;
  for (int i=0; i<gr.GetN(); i++) {
    const auto& y = gr.GetY()[i];
    const auto& yErrLo = gr.GetEYlow()[i];
    const auto& yErrHi = gr.GetEYhigh()[i];
    auto yErr = std::max(yErrLo, yErrHi);
    if (yErr==0.) { yErr = 1E-9; }
    const auto yVar = yErr*yErr;
    const auto w = 1./yVar;
    RMS += w*(y-mean)*(y-mean);
  }
  RMS = std::sqrt(RMS/N);
  return std::make_tuple(mean, meanErr, RMS);
};


void drawParametersGraph(GraphSextaMap2_t& graphMap, const std::string& outDir, const bool& isMC)
{
  //
  // Set Style
  setStyle();
  //
  std::cout << "[INFO] Drawing the output graphs" << std::endl;
  //
  // Draw all graphs
  for (auto& o : graphMap) {
    for (auto& c : o.second) {
      for (auto& p : c.second) {
	for (auto& pp : p.second) {
	  for (auto& v : pp.second) {
	    //
	    const std::string& obj = o.first;
	    const std::string& col = c.first;
	    const std::string& var = v.first;
	    const std::string& par = p.first;
	    //
	    auto& ggr = graphMap.at(obj).at(col).at(par).at(pp.first).at(var);
	    //
	    // Sort points
	    for (auto& gr : ggr) { gr.second.Sort(); }
	    //
	    // Extract inclusive bin (largest width bin)
	    int incBin = -1;
	    double xMaxErr = -9999999999.;
	    TGraphAsymmErrors incGraph;
	    for (int i=0; i<ggr.at("Err_Stat").GetN(); i++) {
	      const auto& xErr = ggr.at("Err_Stat").GetEXlow()[i];
	      if (xErr > xMaxErr) {
		incBin = i;
		xMaxErr = xErr;
	      }
	    }
	    if (incBin>=0) {
	      incGraph.Set(1);
	      incGraph.SetPoint(0, 0.0, ggr.at("Err_Stat").GetY()[incBin]);
	      incGraph.SetPointEYlow(0, ggr.at("Err_Stat").GetEYlow()[incBin]);
	      incGraph.SetPointEYhigh(0, ggr.at("Err_Stat").GetEYhigh()[incBin]);
	      for (auto& gr : ggr) { gr.second.RemovePoint(incBin); }
	    }
	    //
	    // Exclude intermediate bins
	    std::vector<int> excBin;
	    for (int i=0; i<ggr.at("Err_Stat").GetN(); i++) {
	      const auto& xVal1 = ggr.at("Err_Stat").GetX()[i];
	      const auto& xErr1 = ggr.at("Err_Stat").GetEXlow()[i];
	      for (int j=1; j<ggr.at("Err_Stat").GetN(); j++) {
		const auto& xVal2 = ggr.at("Err_Stat").GetX()[j];
		const auto& xErr2 = ggr.at("Err_Stat").GetEXlow()[j];
		if (xVal1>xVal2 && xVal1-xErr1+0.001 < xVal2+xErr2) {
		  const auto& excB = (xErr1>xErr2 ? i : j);
		  if (xErr1!=xErr2 && !contain(excBin, excB)) { excBin.push_back(excB); }
		}
	      }
	    }
	    TGraphAsymmErrors midGraph; midGraph.Set(excBin.size());
	    for (uint i=0; i<excBin.size(); i++) {
	      const auto& excB = excBin[i];
	      midGraph.SetPoint(i, ggr.at("Err_Stat").GetX()[excB], ggr.at("Err_Stat").GetY()[excB]);
	      midGraph.SetPointEYlow(i, ggr.at("Err_Stat").GetEYlow()[excB]);
	      midGraph.SetPointEYhigh(i, ggr.at("Err_Stat").GetEYhigh()[excB]);
	      midGraph.SetPointEXlow(i, ggr.at("Err_Stat").GetEXlow()[excB]);
	      midGraph.SetPointEXhigh(i, ggr.at("Err_Stat").GetEXhigh()[excB]);
	    }
	    for (int i=0; i<midGraph.GetN(); i++) {
	      const auto& xVal1 = midGraph.GetX()[i];
	      const auto& xErr1 = midGraph.GetEXlow()[i];
	      for (int j=1; j<ggr.at("Err_Stat").GetN(); j++) {
		const auto& xVal2 = ggr.at("Err_Stat").GetX()[j];
		const auto& xErr2 = ggr.at("Err_Stat").GetEXlow()[j];
		if (xVal1==xVal2 && xErr1==xErr2) {
		  for (auto& gr : ggr) { gr.second.RemovePoint(j); }
		  break;
		}
	      }
	    }
	    //
	    // Create Canvas
	    TCanvas c("c", "c", 1000, 1000); c.cd();
	    //
	    // Create the Text Info
	    TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035); float dy = 0;
	    std::vector< std::string > textToPrint;
	    textToPrint.push_back(formatDecayLabel(var, obj));
	    for (const auto& ppN : pp.first.first) { textToPrint.push_back(formatObsRange(ppN)); }
	    for (const auto& ppN : pp.first.second) { textToPrint.push_back(formatObsRange(ppN)); }
	    //
	    // Declare the graph vector (for drawing with markers)
	    uint iGr=0;
	    std::vector<std::string> legLblV;
	    std::vector< std::vector<TGraphAsymmErrors> > grDiVec;
	    //
	    auto& graph = pp.second.at(var);
	    // Draw graph
	    std::vector<TGraphAsymmErrors> grV;
	    grV.push_back(graph.at("Err_Stat"));
	    if (contain(graph, "Err_Syst")) {
	      grV.push_back(graph.at("Err_Tot"));
	      grV.push_back(graph.at("Err_Tot"));
	    }
	    grDiVec.push_back(grV);
	    auto& grVec = grDiVec.back();
	    if (contain(graph, "Err_Syst")) {
	      for (int i=0; i<grVec[0].GetN(); i++) { double x, y; grVec[0].GetPoint(i, x, y); grVec[2].SetPoint(i, x, y+grVec[2].GetErrorYhigh(i)); grVec[3].SetPoint(i, x, y-grVec[3].GetErrorYlow(i)); }
	    }
	    for (auto& gr : grVec) { formatResultsGraph(gr, var, par, obj, col, COLOR[iGr]); }; iGr++;
	    if (contain(graph, "Err_Syst")) {
	      for (uint j=1; j<grVec.size(); j++) {
		grVec[j].SetMarkerSize(0);
		for (int i=0; i<grVec[j].GetN(); i++) { grVec[j].SetPointEYhigh(i, 0.0); grVec[j].SetPointEYlow(i, 0.0); }
		for (int i=0; i<grVec[j].GetN(); i++) { grVec[j].SetPointEXhigh(i, 0.5*grVec[j].GetErrorXhigh(i)); grVec[j].SetPointEXlow(i, 0.5*grVec[j].GetErrorXlow(i)); }
	      }
	    }
	    if (incBin>=0) {
	      incGraph.SetPoint(0, grVec[0].GetXaxis()->GetXmax()-0.1, incGraph.GetY()[0]);
	      incGraph.SetMarkerColor(kRed);
	      incGraph.SetMarkerSize(2.0);
	      incGraph.SetMarkerStyle(21);
	    }
	    if (excBin.size()>0) {
	      midGraph.SetMarkerColor(kBlue);
	      midGraph.SetMarkerSize(2.0);
	      midGraph.SetMarkerStyle(21);
	    }
	    //for (int i=0; i<grVec[0].GetN(); i++) { grVec[0].SetPointEXhigh(i, 0.0); grVec[0].SetPointEXlow(i, 0.0); }
	    // Compute the weighted mean and RMS
	    const auto& grStat = computeStats(grVec[0]);
	    const auto meanVal = std::get<0>(grStat);
	    const auto meanErr = std::get<1>(grStat);
	    const auto RMS = std::get<2>(grStat);
	    TLine line(grVec[0].GetXaxis()->GetXmin(), meanVal, grVec[0].GetXaxis()->GetXmax(), meanVal);
	    line.SetLineStyle(7);
	    line.SetLineColor(kBlack);
	    line.SetLineWidth(4);
	    // Draw the graphs
	    grVec[0].Draw("ap");
	    line.Draw("same");
	    incGraph.Draw("samep");
	    midGraph.Draw("samep");
	    if (contain(graph, "Err_Syst")) {grVec[1].Draw("samep"); grVec[2].Draw("samep"); }
	    grVec[0].Draw("samep");
	    // Add legend text
	    std::string legLbl = ""; for (const auto& ppN : pp.first.first) { legLbl += formatObsRange(ppN)+" , "; };
	    if (legLbl.find(" , ")!=std::string::npos) { legLbl.erase(legLbl.find(" , "), 3); }
	    legLblV.push_back(legLbl);
	    //
	    // Initialize the Legend
	    std::unique_ptr<TLegend> leg;
	    if (legLblV.size()>1) {
	      // Define the position and size of the legend
	      double xmin = 0.49 , xmax = 0.66 , ymin = 0.58 , ymax = 0.69 , legSize = 0.047;
	      ymin = std::max(ymax-legLblV.size()*legSize*1.2, ymin);
	      leg.reset(new TLegend(xmin, ymin, xmax, ymax));
	      for (uint i=0; i<legLblV.size(); i++) { formatLegendEntry(*leg->AddEntry(&grDiVec[i][0], legLblV[i].c_str(), "p"), legSize); }
	      // Draw the Legend
	      leg->Draw("same");
	    }
	    // Update
	    c.Modified(); c.Update();
	    // Draw the text
	    tex.SetTextSize(0.055); tex.DrawLatex(0.22, 0.84, textToPrint[0].c_str());
	    tex.SetTextSize(0.060); tex.SetTextFont(61); tex.DrawLatex(0.78, 0.84, "CMS"); tex.SetTextFont(62);
	    tex.SetTextSize(0.046); tex.SetTextFont(52); tex.DrawLatex(0.69, 0.79, "Preliminary"); tex.SetTextFont(62);
	    for (uint i=1; i<textToPrint.size(); i++) { tex.SetTextSize(0.045); tex.DrawLatex(0.22, 0.76-dy, textToPrint[i].c_str()); dy+=0.060; }
	    if (var=="Cross_Section") { tex.SetTextSize(0.030); tex.DrawLatex(0.25, 0.17, "Lumi. uncertainty not shown"); }
	    tex.SetTextSize(0.030); tex.DrawLatex(0.66, 0.74, Form("Mean: %.3f #pm %.3f", meanVal, meanErr));
	    tex.SetTextSize(0.030); tex.DrawLatex(0.66, 0.69, Form("RMS: %.3f", RMS));
	    if (incBin>=0) { tex.SetTextSize(0.030); tex.DrawLatex(0.66, 0.64, Form("Incl.: %.3f #pm %.3f", incGraph.GetY()[0], std::max(incGraph.GetEYlow()[0], incGraph.GetEYhigh()[0]))); }
	    // Update
	    c.Modified(); c.Update(); // Pure paranoia
	    //
	    // set the CMS style
	    StringVector_t lumiLabels; getLumiLabels(lumiLabels, "", col, isMC);
	    CMS_lumi(&c, 33, (" "+lumiLabels[0]), lumiLabels[1], false, 0.65, false);
	    // Update
	    c.Modified(); c.Update(); // Pure paranoia
	    //
	    // Create Output Directory
	    const std::string& plotDir = outDir+"/Plots/"+par;
	    makeDir(plotDir + "/png/");
	    makeDir(plotDir + "/pdf/");
	    makeDir(plotDir + "/root/");
	    makeDir(plotDir + "/C/");
	    //
	    // Save Canvas
	    const std::string& grName = v.second.at("Err_Stat").GetName();
	    const auto& name = grName.substr(0, grName.rfind("_Err_Stat"));
	    c.SaveAs((plotDir + "/png/"  + name + ".png" ).c_str());
	    c.SaveAs((plotDir + "/pdf/"  + name + ".pdf" ).c_str());
	    c.SaveAs((plotDir + "/root/" + name + ".root").c_str());
	    //c.SaveAs((plotDir + "/C/"    + name + ".C"   ).c_str());
	    //
	    // Clean up memory
	    c.Clear(); c.Close();
	  }
	}
      }
    }
  }
};


void compareObjectGraph(GraphSextaMap2_t& graphMap, const std::string& outDir, const bool& isMC, const StringVector_t& objTags)
{
  //
  // Set Style
  setStyle();
  //
  std::cout << "[INFO] Drawing the output graphs" << std::endl;
  //
  // Main object
  const auto& mainObj = objTags[0];
  //
  // Draw all graphs
  for (auto& c : graphMap.at(mainObj)) {
    for (auto& p : c.second) {
      for (auto& pp : p.second) {
	for (auto& v : pp.second) {
	  //
	  const std::string& col = c.first;
	  const std::string& var = v.first;
	  const std::string& par = p.first;
	  //
	  const auto mainObjGraph = graphMap.at(mainObj).at(col).at(par).at(pp.first).at(var);
	  //
	  for (const auto& obj : objTags) {
	    if (obj==mainObj) continue;
	    //
	    // Create Canvas
	    TCanvas c("c", "c", 1000, 1000); c.cd();
	    //
	    // Create the Text Info
	    TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035); float dy = 0;
	    std::vector< std::string > textToPrint;
	    std::string sampleLabel = formatObjName(obj);
	    sampleLabel = (obj.rfind("Ups",0)==0 ? "#Upsilon(nS)" : (obj.rfind("Psi",0)==0 ? "#psi(nS)" : sampleLabel));
	    sampleLabel += "#rightarrow #mu^{+}#mu^{#font[122]{\55}}";
	    textToPrint.push_back(sampleLabel);
	    for (const auto& ppN : pp.first.first) { textToPrint.push_back(formatObsRange(ppN)); }
	    for (const auto& ppN : pp.first.second) { textToPrint.push_back(formatObsRange(ppN)); }
	    //
	    // Declare the graph vector (for drawing with markers)
	    uint iGr=0;
	    std::vector<std::string> legLblV;
	    std::vector< std::vector<TGraphAsymmErrors> > grDiVec;
	    //
	    auto& objGraph = graphMap.at(obj).at(col).at(par).at(pp.first).at(var);
	    //
	    auto graph = objGraph;
	    //
	    for (const auto& gr : mainObjGraph) {
	      const auto& mObjGr = gr.second;
	      auto& objGr = graph.at(gr.first);
	      for (int i=0; i<mObjGr.GetN(); i++) {
		const auto& mVal = (mObjGr.GetY()[i]!=0.0 ? mObjGr.GetY()[i] : 1E-9);
		const auto& oVal = objGr.GetY()[i]/mVal;
		objGr.SetPoint(i, objGr.GetX()[i], oVal);
		const auto& mErrUp = mObjGr.GetEYhigh()[i];
		const auto& mErrDw = mObjGr.GetEYlow()[i];
		const auto& ooVal = (objGr.GetY()[i]!=0.0 ? objGr.GetY()[i] : 1E-9);
		const auto& oErrUp = oVal*std::sqrt( std::pow(mErrUp/mVal, 2.0) + std::pow(objGr.GetEYhigh()[i]/ooVal, 2.0) );
		const auto& oErrDw = oVal*std::sqrt( std::pow(mErrDw/mVal, 2.0) + std::pow(objGr.GetEYlow()[i]/ooVal, 2.0) );
		objGr.SetPointEYhigh(i, oErrUp);
		objGr.SetPointEYlow (i, oErrDw);
	      }
	    }
	    //
	    // Sort points
	    for (auto& gr : graph) { gr.second.Sort(); }
	    //
	    // Extract inclusive bin (largest width bin)
	    int incBin = -1;
	    double xMaxErr = -9999999999.;
	    TGraphAsymmErrors incGraph;
	    for (int i=0; i<graph.at("Err_Stat").GetN(); i++) {
	      const auto& xErr = graph.at("Err_Stat").GetEXlow()[i];
	      if (xErr > xMaxErr) {
		incBin = i;
		xMaxErr = xErr;
	      }
	    }
	    if (incBin>=0) {
	      incGraph.Set(1);
	      incGraph.SetPoint(0, 0.0, graph.at("Err_Stat").GetY()[incBin]);
	      incGraph.SetPointEYlow(0, graph.at("Err_Stat").GetEYlow()[incBin]);
	      incGraph.SetPointEYhigh(0, graph.at("Err_Stat").GetEYhigh()[incBin]);
	      for (auto& gr : graph) { gr.second.RemovePoint(incBin); }
	    }
	    //
	    // Exclude intermediate bins
	    std::vector<int> excBin;
	    for (int i=0; i<graph.at("Err_Stat").GetN(); i++) {
	      const auto& xVal1 = graph.at("Err_Stat").GetX()[i];
	      const auto& xErr1 = graph.at("Err_Stat").GetEXlow()[i];
	      for (int j=1; j<graph.at("Err_Stat").GetN(); j++) {
		const auto& xVal2 = graph.at("Err_Stat").GetX()[j];
		const auto& xErr2 = graph.at("Err_Stat").GetEXlow()[j];
		if (xVal1>xVal2 && xVal1-xErr1+0.001 < xVal2+xErr2) {
		  const auto& excB = (xErr1>xErr2 ? i : j);
		  if (xErr1!=xErr2 && !contain(excBin, excB)) { excBin.push_back(excB); }
		}
	      }
	    }
	    TGraphAsymmErrors midGraph; midGraph.Set(excBin.size());
	    for (uint i=0; i<excBin.size(); i++) {
	      const auto& excB = excBin[i];
	      midGraph.SetPoint(i, graph.at("Err_Stat").GetX()[excB], graph.at("Err_Stat").GetY()[excB]);
	      midGraph.SetPointEYlow(i, graph.at("Err_Stat").GetEYlow()[excB]);
	      midGraph.SetPointEYhigh(i, graph.at("Err_Stat").GetEYhigh()[excB]);
	      midGraph.SetPointEXlow(i, graph.at("Err_Stat").GetEXlow()[excB]);
	      midGraph.SetPointEXhigh(i, graph.at("Err_Stat").GetEXhigh()[excB]);
	    }
	    for (int i=0; i<midGraph.GetN(); i++) {
	      const auto& xVal1 = midGraph.GetX()[i];
	      const auto& xErr1 = midGraph.GetEXlow()[i];
	      for (int j=1; j<graph.at("Err_Stat").GetN(); j++) {
		const auto& xVal2 = graph.at("Err_Stat").GetX()[j];
		const auto& xErr2 = graph.at("Err_Stat").GetEXlow()[j];
		if (xVal1==xVal2 && xErr1==xErr2) {
		  for (auto& gr : graph) { gr.second.RemovePoint(j); }
		  break;
		}
	      }
	    }
	    // Draw graph
	    std::vector<TGraphAsymmErrors> grV;
	    grV.push_back(graph.at("Err_Stat"));
	    if (contain(graph, "Err_Syst")) {
	      grV.push_back(graph.at("Err_Tot"));
	      grV.push_back(graph.at("Err_Tot"));
	    }
	    grDiVec.push_back(grV);
	    auto& grVec = grDiVec.back();
	    if (contain(graph, "Err_Syst")) {
	      for (int i=0; i<grVec[0].GetN(); i++) { double x, y; grVec[0].GetPoint(i, x, y); grVec[2].SetPoint(i, x, y+grVec[2].GetErrorYhigh(i)); grVec[3].SetPoint(i, x, y-grVec[3].GetErrorYlow(i)); }
	    }
	    for (auto& gr : grVec) { formatResultsGraph(gr, var, par, obj, col, COLOR[iGr], true); }; iGr++;
	    if (contain(graph, "Err_Syst")) {
	      for (uint j=1; j<grVec.size(); j++) {
		grVec[j].SetMarkerSize(0);
		for (int i=0; i<grVec[j].GetN(); i++) { grVec[j].SetPointEYhigh(i, 0.0); grVec[j].SetPointEYlow(i, 0.0); }
		for (int i=0; i<grVec[j].GetN(); i++) { grVec[j].SetPointEXhigh(i, 0.5*grVec[j].GetErrorXhigh(i)); grVec[j].SetPointEXlow(i, 0.5*grVec[j].GetErrorXlow(i)); }
	      }
	    }
	    if (incBin>=0) {
	      incGraph.SetPoint(0, grVec[0].GetXaxis()->GetXmax()-0.1, incGraph.GetY()[0]);
	      incGraph.SetMarkerColor(kRed);
	      incGraph.SetMarkerSize(2.0);
	      incGraph.SetMarkerStyle(21);
	    }
	    if (excBin.size()>0) {
	      midGraph.SetMarkerColor(kBlue);
	      midGraph.SetMarkerSize(2.0);
	      midGraph.SetMarkerStyle(21);
	    }
	    //for (int i=0; i<grVec[0].GetN(); i++) { grVec[0].SetPointEXhigh(i, 0.0); grVec[0].SetPointEXlow(i, 0.0); }
	    // Compute the weighted mean and RMS
	    const auto& grStat = computeStats(grVec[0]);
	    const auto meanVal = std::get<0>(grStat);
	    const auto meanErr = std::get<1>(grStat);
	    const auto RMS = std::get<2>(grStat);
	    TLine line(grVec[0].GetXaxis()->GetXmin(), meanVal, grVec[0].GetXaxis()->GetXmax(), meanVal);
	    line.SetLineStyle(7);
	    line.SetLineColor(kBlack);
	    line.SetLineWidth(4);
	    // Draw the graphs
	    grVec[0].Draw("ap");
	    incGraph.Draw("samep");
	    midGraph.Draw("samep");
	    line.Draw("same");
	    if (contain(graph, "Err_Syst")) {grVec[1].Draw("samep"); grVec[2].Draw("samep"); }
	    grVec[0].Draw("samep");
	    // Add legend text
	    std::string legLbl = ""; for (const auto& ppN : pp.first.first) { legLbl += formatObsRange(ppN)+" , "; };
	    if (legLbl.find(" , ")!=std::string::npos) { legLbl.erase(legLbl.find(" , "), 3); }
	    legLblV.push_back(legLbl);
	    //
	    // Initialize the Legend
	    std::unique_ptr<TLegend> leg;
	    if (legLblV.size()>1) {
	      // Define the position and size of the legend
	      double xmin = 0.49 , xmax = 0.66 , ymin = 0.58 , ymax = 0.69 , legSize = 0.047;
	      ymin = std::max(ymax-legLblV.size()*legSize*1.2, ymin);
	      leg.reset(new TLegend(xmin, ymin, xmax, ymax));
	      for (uint i=0; i<legLblV.size(); i++) { formatLegendEntry(*leg->AddEntry(&grDiVec[i][0], legLblV[i].c_str(), "p"), legSize); }
	      // Draw the Legend
	      leg->Draw("same");
	    }
	    // Update
	    c.Modified(); c.Update();
	    // Draw the text
	    tex.SetTextSize(0.055); tex.DrawLatex(0.22, 0.84, textToPrint[0].c_str());
	    tex.SetTextSize(0.060); tex.SetTextFont(61); tex.DrawLatex(0.78, 0.84, "CMS"); tex.SetTextFont(62);
	    tex.SetTextSize(0.046); tex.SetTextFont(52); tex.DrawLatex(0.69, 0.79, "Preliminary"); tex.SetTextFont(62);
	    for (uint i=1; i<textToPrint.size(); i++) { tex.SetTextSize(0.045); tex.DrawLatex(0.22, 0.76-dy, textToPrint[i].c_str()); dy+=0.060; }
	    if (var=="Cross_Section") { tex.SetTextSize(0.030); tex.DrawLatex(0.25, 0.17, "Lumi. uncertainty not shown"); }
	    tex.SetTextSize(0.030); tex.DrawLatex(0.66, 0.74, Form("Mean: %.3f #pm %.3f", meanVal, meanErr));
	    tex.SetTextSize(0.030); tex.DrawLatex(0.66, 0.69, Form("RMS: %.3f", RMS));
	    if (incBin>=0) { tex.SetTextSize(0.030); tex.DrawLatex(0.66, 0.64, Form("Incl.: %.3f #pm %.3f", incGraph.GetY()[0], std::max(incGraph.GetEYlow()[0], incGraph.GetEYhigh()[0]))); }
	    if (var.rfind("Sigma",0)==0 || var=="m") { tex.SetTextSize(0.030); tex.DrawLatex(0.62, 0.58, Form("PDG mass ratio: %.2f", ANA::MASS.at(obj).at("Val")/ANA::MASS.at(mainObj).at("Val"))); }
	    // Update
	    c.Modified(); c.Update(); // Pure paranoia
	    //
	    // set the CMS style
	    StringVector_t lumiLabels; getLumiLabels(lumiLabels, "", col, isMC);
	    CMS_lumi(&c, 33, (" "+lumiLabels[0]), lumiLabels[1], false, 0.65, false);
	    // Update
	    c.Modified(); c.Update(); // Pure paranoia
	    //
	    // Create Output Directory
	    const std::string& plotDir = outDir+"/Plots/"+par;
	    makeDir(plotDir + "/png/");
	    makeDir(plotDir + "/pdf/");
	    makeDir(plotDir + "/root/");
	    makeDir(plotDir + "/C/");
	    //
	    // Save Canvas
	    const std::string& grName = v.second.at("Err_Stat").GetName();
	    const auto& name = grName.substr(0, grName.rfind("_Err_Stat"));
	    c.SaveAs((plotDir + "/png/"  + name + ".png" ).c_str());
	    c.SaveAs((plotDir + "/pdf/"  + name + ".pdf" ).c_str());
	    c.SaveAs((plotDir + "/root/" + name + ".root").c_str());
	    //c.SaveAs((plotDir + "/C/"    + name + ".C"   ).c_str());
	    //
	    // Clean up memory
	    c.Clear(); c.Close();
	  }
	}
      }
    }
  }
};


#endif // ifndef resultUtils_h
