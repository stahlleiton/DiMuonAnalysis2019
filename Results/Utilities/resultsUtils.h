#ifndef Utilities_resultUtils_h
#define Utilities_resultUtils_h
// Auxiliary Headers
#include "bin.h"
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
// RooFit headers
// c++ headers
#include <iostream>
#include <string>
#include <map>
// CMS headers


// ------------------ TYPE -------------------------------
typedef std::pair< anabin        , anabin            > BinPair_t;
typedef std::map< anabin         , DoubleDiMap_t     > DoubleBinDiMap_t;
typedef std::map< std::string    , DoubleBinDiMap_t  > VarBinMap_t;
typedef std::map< std::string    , VarBinMap_t       > VarBinDiMap_t;
typedef std::map< std::string    , VarBinDiMap_t     > VarBinTriMap_t;
typedef std::map< anabin         , double            > BinMap_t;
typedef std::map< std::string    , BinMap_t          > BinDiMap_t;
typedef std::map< std::string    , BinDiMap_t        > BinTriMap_t;
typedef std::map< std::string    , BinTriMap_t       > BinQuadMap_t;
typedef std::map< std::string    , BinQuadMap_t      > BinPentaMap_t;
typedef std::map< std::string    , BinPentaMap_t     > BinSextaMap_t;
typedef std::map< std::string    , BinSextaMap_t     > BinSeptaMap_t;
typedef std::map< std::string    , BinSeptaMap_t     > BinOctaMap_t;
typedef std::map< BinPair_t      , std::set<anabin>  > BinSetMap_t;
typedef std::map< std::string    , BinSetMap_t       > BinSetDiMap_t;
typedef std::map< std::string    , BinSetDiMap_t     > BinSetTriMap_t;
typedef std::map< std::string    , BinSetTriMap_t    > BinCont_t;
typedef std::map< std::string    , TGraphAsymmErrors > GraphMap_t;
typedef std::map< std::string    , GraphMap_t        > GraphDiMap_t;
typedef std::map< BinPair_t      , GraphDiMap_t      > GraphTriMap_t;
typedef std::map< std::string    , GraphTriMap_t     > GraphQuadMap_t;
typedef std::map< std::string    , GraphQuadMap_t    > GraphPentaMap_t;
typedef std::map< std::string    , GraphPentaMap_t   > GraphSextaMap_t;


// ------------------ GLOBAL -------------------------------
const auto& COLOR = std::vector<int>({kBlack, kRed, kGreen+2, kBlue+2, kAzure-7, kYellow-3, kMagenta+3, kCyan+3, kOrange+9, kSpring+4, kTeal+4, kViolet+9, kPink+8});


// ------------------ FUNCTION -------------------------------


double sumErrors(const std::vector<double>& errV)
{
  double sumErrSqr = 0.;
  for (const auto& e : errV) { sumErrSqr += std::pow(e, 2.0); }
  return std::sqrt(sumErrSqr);
};


void getBins(BinSetMap_t& binS, const std::set<anabin>& allBin, const std::string& obs, const StringVector_t& vars)
{
  if (contain(vars, obs)) return;
  std::set<BinPair_t> bSet;
  for (const auto& b : allBin) {
    if (b.getbin(obs).name()=="") continue;
    anabin varBin, incBin;
    bool keepBin = true;
    for (const auto& p : b) {
      if (p.name()==obs) continue;
      if (contain(vars, p.name())) { varBin.setbin(p); } else { incBin.setbin(p); }
      bool isInc = false;
      if      (p.name()=="Cand_Rap"   ) { isInc = (isEqual(p.high(),   2.4, 1) && isEqual(p.low(), -2.4, 1)); }
      else if (p.name()=="Cand_AbsRap") { isInc = (isEqual(p.high(),   2.4, 1) && isEqual(p.low(),  0.0, 1)); }
      else if (p.name()=="Centrality" ) { isInc = (isEqual(p.high(), 100.0, 1) && isEqual(p.low(),  0.0, 1)); }
      keepBin = (keepBin && (contain(vars, p.name()) ? !isInc : isInc));
    }
    const auto& binP = std::make_pair(varBin, incBin);
    if (keepBin) { binS[binP].insert(b); bSet.insert(binP); }
  }
  // Check if we only have one main bin
  for (const auto& binP : bSet) {
    std::set<binF> mainBin; for (const auto& b : binS.at(binP)) { mainBin.insert(b.getbin(obs)); }
    if (mainBin.size()<=1) { binS.erase(binP); }
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


void defineBins(BinCont_t& binMap, const BinSextaMap_t& var)
{
  for (const auto& o : var) {
    for (const auto& c : o.second) {
      // Get all bins
      std::set<anabin> allBin;
      for (const auto& pd : c.second) { for (const auto& b : pd.second.begin()->second.at("Val")) {  allBin.insert(b.first); } }
      // Get all observables
      StringVector_t allObs;
      for (const auto& b : allBin) { for (const auto& v : b) { if (!contain(allObs, v.name())) { allObs.push_back(v.name()); } } }
      // Combine all observables
      std::vector<StringVector_t> cmbObs;
      combineBins(cmbObs, allObs);
      // Loop over each observable
      for (const auto& p : allObs) {
	auto& binS = binMap[o.first][c.first][p];
	// Loop over each combination of observables
	for (const auto& v : cmbObs) { getBins(binS, allBin, p, v); }
      }
    }
  }
};


void iniResultsGraph(GraphSextaMap_t& graphMap, const BinCont_t& binMap, const BinSextaMap_t& var)
{
  //
  const StringVector_t gType = {"Err_Tot", "Err_Stat"};
  //
  for (const auto& o : binMap) {
    for (const auto& c : o.second) {
      for (const auto& p : c.second) {
	for (const auto& pp : p.second) {
	  for (const auto& v : var.at(o.first).at(c.first).begin()->second) {
	    for (const auto& t : gType) {
	      auto& graph = graphMap[o.first][c.first][p.first][pp.first][v.first][t];
	      graph.Set(pp.second.size());
	      // Set Graph Name
	      std::string ppLbl = "_"; if (pp.first.first.size()>0) { ppLbl += "_"; for (const auto& ppN : pp.first.first) { ppLbl += Form("%s_%.0f_%.0f_", ppN.name().c_str(), ppN.low()*100., ppN.high()*100.); } }
	      const auto& name = ("gr_"+o.first+"_"+c.first+"_"+v.first+"_"+p.first+ppLbl+t);
	      graph.SetName(name.c_str());
	    }
	  }
	}
      }
    }
  }
};


void fillResultsGraph(GraphSextaMap_t& graphMap, const BinCont_t& binMap, const BinSextaMap_t& iVar)
{
  //
  std::cout << "[INFO] Filling the output graphs" << std::endl;
  //
  for (auto& o : graphMap) {
    for (auto& c : o.second) {
      for (auto& p : c.second) {
	for (auto& pp : p.second) {
	  //
	  // Determine the bin index
	  const auto& bins = binMap.at(o.first).at(c.first).at(p.first).at(pp.first);
	  std::map<anabin , uint> binIdx; uint iBin = 0; for (const auto& b : bins) { binIdx[b] = iBin; iBin++; }
	  //
	  for (auto& v : pp.second) {
	    for (auto& gr : v.second) {
	      auto& graph = gr.second;
	      for (const auto& b : bins) {
		//
		// Extract the parameters needed for each axis
		//
		// X Value
		const auto& binX = b.getbin(p.first);
		const auto& X = binX.mean(); // Mean value of bin
		// X Error
		const auto& Err_X = binX.width(); ; // Width of bin
		const auto& Err_X_High = Err_X;
		const auto& Err_X_Low  = Err_X;
		//
		double norm = 1.0;
		if (v.first=="Cross_Section") {
		  norm *= (binX.high()-binX.low());
		  if (binX.name()=="Centrality") { norm *= 0.01; }
		}
		//
		for (const auto& pd : iVar.at(o.first).at(c.first)) {
		  const auto& var = pd.second.at(v.first);
		  if (!contain(var.at("Val"), b)) continue;
		  //
		  // Y Value
		  const auto& Y = var.at("Val").at(b)/norm;
		  //
		  // Compute total systematic error
		  const auto& Err_Y_Syst_High = var.at("Err_Syst_High").at(b)/norm;
		  const auto& Err_Y_Syst_Low  = var.at("Err_Syst_Low" ).at(b)/norm;
		  //
		  // Compute total statistic error
		  const auto& Err_Y_Stat_High = var.at("Err_Stat_High").at(b)/norm;
		  const auto& Err_Y_Stat_Low  = var.at("Err_Stat_Low" ).at(b)/norm;
		  //
		  // Y Error
		  double Err_Y_High = 0.0 , Err_Y_Low  = 0.0;
		  if (gr.first=="Err_Tot") {
		    Err_Y_High = sumErrors({Err_Y_Stat_High , Err_Y_Syst_High});
		    Err_Y_Low  = sumErrors({Err_Y_Stat_Low  , Err_Y_Syst_Low });
		  }
		  else if (gr.first=="Err_Stat") {
		    Err_Y_High = Err_Y_Stat_High;
		    Err_Y_Low  = Err_Y_Stat_Low;
		  }
		  else if (gr.first=="Err_Syst") {
		    Err_Y_High = Err_Y_Syst_High;
		    Err_Y_Low  = Err_Y_Syst_Low;
		  }
		  //
		  // Fill the nominal graph
		  //
		  const auto& iBin = binIdx.at(b);
		  graph.SetPoint(iBin, X, Y);
		  graph.SetPointError(iBin, Err_X_Low, Err_X_High, Err_Y_Low, Err_Y_High);
		}
	      }
	    }
	  }
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
  const auto& obsL = StringMap_t({{"Cand_Mass", "M"}, {"Cand_Pt", "p_{T}"}, {"Cand_Rap", "y"},
				  {"Cand_AbsRap", "|y|"}, {"Centrality", "Cent"}, {"NTrack", "Ntrk"}});
  return (contain(obsL, obs) ? obsL.at(obs) : obs);
};


std::string formatObsName(const std::string& obs)
{
  const auto& obsL = StringMap_t({{"Cand_Mass", "Mass [GeV/c^{2}]"}, {"Cand_Pt", "p_{T} [GeV/c]"}, {"Cand_Rap", "y"},
				  {"Cand_AbsRap", "|y|"}, {"Centrality", "Centrality [%]"}, {"NTrack", "Track Multiplicity"}});
  std::string obsF = (obs.rfind("Cand_",0)==0 ? "#mu^{+}#mu^{#font[122]{\55}} " : "");
  obsF += (contain(obsL, obs) ? obsL.at(obs) : "");
  return obsF;
};


std::string formatObsRange(const binF& bin)
{
  const auto& bName = bin.name();
  const auto& bN = formatObsTag(bName);
  const auto& bU = obsUnit(bName);
  const auto& obsMin = DoubleMap_t({{"Cand_Mass",   0.0}, {"Cand_Pt",   0.0}, {"Cand_Rap", -2.5}, {"Cand_AbsRap", 0.0}, {"Centrality",   0.0}, {"NTrack",   0.0}});
  const auto& obsMax = DoubleMap_t({{"Cand_Mass", 500.0}, {"Cand_Pt", 500.0}, {"Cand_Rap",  2.5}, {"Cand_AbsRap", 2.5}, {"Centrality", 100.0}, {"NTrack", 500.0}});
  const auto& min = (contain(obsMin, bName) ? obsMin.at(bName) : -999.0);
  const auto& max = (contain(obsMax, bName) ? obsMax.at(bName) :  999.0);
  std::string label = "";
  if (bin.low()>min && bin.high()>=max) { label = Form("%g%s #leq %s", bin.low(), bU.c_str(), bN.c_str()); }
  else if (bin.low()<=min && bin.high()<max) { label = Form("%s < %g%s", bN.c_str(), bin.high(), bU.c_str()); }
  else if (bin.low()>min && bin.high()<max) { label = Form("%g #leq %s < %g%s", bin.low(), bN.c_str(), bin.high(), bU.c_str()); }
  return label;
};


std::string formatObjName(const std::string& obj)
{
  const auto& objL = StringMap_t({{"DY", "Z/#gamma*"}, {"TTbar", "t#bar{t}"}, {"JPsi", "J/#psi"}, {"Psi2S", "#psi(2S)"},
				  {"Ups1S", "#Upsilon(1S)"}, {"Ups2S", "#Upsilon(2S)"}, {"Ups3S", "#Upsilon(3S)"}});
  auto ob = obj;
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


std::string formatResultVarName(const std::string& var, const std::string& par, const std::string& obj, const std::string& col)
{
  std::string label = "";
  // Format object name
  // Add result variable
  if (var=="Cross_Section") {
    const auto& parF = formatObsTag(par);
    const std::string& unit = (col.rfind("PbPb",0)==0 ? "#mub" : (col.rfind("PP",0)==0 ? "pb" : "nb"));
    label = "B #times d#sigma/d"+parF+" ["+unit+"]";
  }
  else if (var=="RatioTo1S" || var=="R") {
    const auto& objF = formatObjName(obj);
    const std::string& ref = (obj.rfind("Ups",0)==0 ? "#Upsilon(1S)" : (obj.rfind("Psi",0)==0 ? "J/#psi" : ""));
    label = objF+"/"+ref+" Ratio";
  }
  // Add fit parameter
  else {
    const auto& objF = formatObjName(obj);
    const auto& varF = formatParName(var);
    label = objF+" "+varF;
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


void setRangeYAxisGraph(TGraphAsymmErrors& graph, const double& fracDo = 0.1, const double& fracUp = 0.5)
{
  // Find maximum and minimum points of Plot to rescale Y axis
  double YMax = -1e99;
  for (int i=0; i<=graph.GetN(); i++) { double x, y; graph.GetPoint(i, x, y); y += graph.GetErrorY(i); YMax = std::max(YMax, y); }
  double YMin = 1e99;
  for (int i=0; i<=graph.GetN(); i++) { double x, y; graph.GetPoint(i, x, y); y -= graph.GetErrorY(i); YMin = std::min(YMin, y); }
  //
  const double YTot = (YMax - YMin)/(1.0 - fracUp - fracDo);
  const double Yup   = YMax + fracUp*YTot;
  const double Ydown = YMin - fracDo*YTot;
  //
  graph.GetYaxis()->SetRangeUser(Ydown, Yup);
};


void formatResultsGraph(TGraphAsymmErrors& graph, const std::string& var, const std::string& par, const std::string& obj, const std::string& col, const int& color=-1)
{
  //
  // Set the Axis Titles
  std::string xLabel = formatObsName(par);
  std::string yLabel = formatResultVarName(var, par, obj, col);
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
  double xMin = graph.GetXaxis()->GetXmin()-0.1;
  double xMax = graph.GetXaxis()->GetXmax()+0.1;
  if      (par=="Cand_Pt"    ) { xMin = -0.1; xMax = 20.0; }
  else if (par=="Cand_Rap"   ) { xMin = -2.5; xMax = 2.5; }
  else if (par=="Cand_AbsRap") { xMin = -0.1; xMax = 2.5; }
  else if (par=="NTrack") { xMin = 0.0; xMax = 260.0; }
  graph.GetXaxis()->SetRangeUser(xMin, xMax);
  graph.GetXaxis()->SetLimits(xMin, xMax);
  if (par=="NTrack") graph.GetXaxis()->SetNdivisions(510);
  else graph.GetXaxis()->SetNdivisions(505);
  // Y-axis
  graph.GetYaxis()->CenterTitle(kTRUE);
  graph.GetYaxis()->SetTitleOffset(1.05);
  graph.GetYaxis()->SetTitleSize(0.065);
  graph.GetYaxis()->SetLabelSize(0.035);
  setRangeYAxisGraph(graph, 0.1, 0.4);
};


void drawResultsGraph(GraphSextaMap_t& graphMap, const std::string& outDir, const bool& isMC)
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
	    // Create Canvas
	    TCanvas c("c", "c", 1000, 1000); c.cd();
	    //
	    // Create the Text Info
	    TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035); float dy = 0;
	    std::vector< std::string > textToPrint;
	    std::string sampleLabel = formatObjName(obj);
	    if (var=="R" || var=="RatioTo1S") { sampleLabel = (obj.rfind("Ups",0)==0 ? "#Upsilon(nS)" : (obj.rfind("Psi",0)==0 ? "#psi(nS)" : sampleLabel)); }
	    sampleLabel += " #rightarrow #mu^{+} + #mu^{#font[122]{\55}}";
	    textToPrint.push_back(sampleLabel);
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
	    //for (int i=0; i<grVec[0].GetN(); i++) { grVec[0].SetPointEXhigh(i, 0.0); grVec[0].SetPointEXlow(i, 0.0); }
	    // Draw the graphs
	    grVec[0].Draw("ap");
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
	    c.SaveAs((plotDir + "/C/"    + name + ".C"   ).c_str());
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
