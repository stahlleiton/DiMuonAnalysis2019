#ifndef drawUtils_h
#define drawUtils_h

#include "TROOT.h"
#include "TSystem.h"
#include "TIterator.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGaxis.h"

#include "RooWorkspace.h"
#include "RooClassFactory.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooPlot.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "RooStringVar.h"
#include "RooCurve.h"
#include "RooHist.h"
#include "RooHistPdf.h"
#include "RooFitResult.h"

#include "../../../Utilities/CMS/tdrstyle.C"
#include "../../../Utilities/CMS/CMS_lumi.C"
#include "../../../Utilities/RooGoF.C"

#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>


typedef std::unordered_map< std::string , std::string > StringMap_d;
typedef std::unordered_map< std::string , StringMap_d > StringDiMap_d;
typedef std::unordered_map< std::string , std::unique_ptr<RooPlot> > RooPlotPtrMap_d;
typedef std::unordered_map< std::string , TPad*       > PadPtrMap_d; // Unique Pointer does produce Segmentation Fault, so don't use it


void formatLegendEntry(TLegendEntry& e, const double& size=0.060)
{
  e.SetTextSize(size);
};


void getPoint(const RooHist& rH, const uint& i, double& x, double& y, double& exl, double& exh, double& eyl, double& eyh)
{
  //
  rH.GetPoint(i, x, y);
  //
  eyl = rH.GetErrorYlow(i);
  eyh = rH.GetErrorYhigh(i);
  //
  exl = rH.GetErrorXlow(i);
  exh = rH.GetErrorXhigh(i);
  if (exl<=0.0 ) { exl = rH.GetErrorX(i); }
  if (exh<=0.0 ) { exh = rH.GetErrorX(i); }
  if (exl<=0.0 ) { exl = 0.5*rH.getNominalBinWidth(); }
  if (exh<=0.0 ) { exh = 0.5*rH.getNominalBinWidth(); }
};


bool rooPlotToTH1(TH1D& hData, TH1D& hFit, const RooPlot& frame, const bool& useAverage = true)
{
  // Find curve object
  const auto& rFit = dynamic_cast<RooCurve*>(frame.findObject(0, RooCurve::Class()));
  if (!rFit) { std::cout << "[ERROR] The latest RooCurve was not found" << std::endl; return false; }
  // Find histogram object
  const auto& rData = dynamic_cast<RooHist*>(frame.findObject(0, RooHist::Class()));
  if (!rData) { std::cout << "[ERROR] The latest RooHist was not found" << std::endl; return false; }
  // Determine range of curve
  double xstart, xstop, yDummy;
  rFit->GetPoint(0, xstart, yDummy);
  rFit->GetPoint((rFit->GetN()-1), xstop, yDummy);
  // Get Binning
  std::vector<double> binV;
  for (int i = 0; i < rData->GetN(); i++) {
    double x, y; rData->GetPoint(i, x, y);
    // Only consider bins inside curve range
    if (x<xstart || x>xstop) continue;
    const double& binW = rData->getNominalBinWidth();
    binV.push_back(x - (binW*0.5));
    if (i==(rData->GetN()-1)) { binV.push_back(x + (binW*0.5)); }
  }
  const uint& nBin = (binV.size()-1);
  double bin[nBin+1];
  for (uint i = 0; i < binV.size(); i++) { bin[i] = binV[i]; }
  //
  hData.Reset(); hData = TH1D(Form("hData_%s", rData->GetName()), rData->GetTitle(), nBin, bin); hData.Sumw2();
  hFit.Reset();  hFit  = TH1D(Form("hFit_%s" , rFit->GetName()) , rFit->GetTitle() , nBin, bin); hFit.Sumw2();
  // Set Histogram entries
  for (uint i = 0; i < nBin; i++) {
    double x, dataVal, exl, exh, eyl, eyh;
    getPoint(*rData, i, x, dataVal, exl, exh, eyl, eyh);
    double fitVal = 0.0;
    if (useAverage) { fitVal = rFit->average(x-exl, x+exh); }
    else            { fitVal = rFit->interpolate(x);        }
    hData.SetBinContent((i+1), dataVal);
    hData.SetBinError((i+1), std::sqrt((eyl*eyl + eyh*eyh)/2.0));
    hFit.SetBinContent((i+1), fitVal);
    hFit.SetBinError((i+1), 0.0);
  }
  return true;
};


void setPlotRange(RooPlot& frame, const RooWorkspace& ws, const std::string& varName, const std::string& dsName, const bool& setLogScale, const int& nBins)
{
  // Find maximum and minimum points of Plot to rescale Y axis
  TH1D hData, hFit;
  if (ws.data(dsName.c_str())!=NULL) {
    auto h = std::unique_ptr<TH1>(ws.data(dsName.c_str())->createHistogram("hist", *ws.var(varName.c_str()), RooFit::Binning(nBins, frame.GetXaxis()->GetXmin(), frame.GetXaxis()->GetXmax())));
    const auto& h1D = dynamic_cast<TH1D*>(h.get());
    if (h1D) { hData = *h1D; }
  }
  else {
    if (!rooPlotToTH1(hData, hFit, frame)) { std::cout << "[ERROR] Could not find the RooHist from the frame!" << std::endl; return; }
  }
  double YMax = hData.GetBinContent(hData.GetMaximumBin());
  double YMin = 1e99;
  for (int i=1; i<=hData.GetNbinsX(); i++) if (hData.GetBinContent(i)>0) YMin = min(YMin, hData.GetBinContent(i));
  double Yup(0.),Ydown(0.);
  if(setLogScale)
  {
    Yup = YMax*pow((YMax), 0.55);
    Ydown = 1.0;
  }
  else
  {
    Yup = YMax+(YMax-0.0)*0.55;
    Ydown = 0.0;
  }
  frame.GetYaxis()->SetRangeUser(Ydown,Yup);
  //
  // Draw Lines for the var range if cut
  if (varName=="MET" && ws.data((dsName+"_FIT").c_str())!=NULL) {
    const double& varMin = ws.var(varName.c_str())->getMin(Form("%sWindow", varName.c_str()));
    if (varMin > 0.0) {
      auto minline = new TLine(varMin, 0.0, varMin, (setLogScale?(Ydown*TMath::Power((Yup/Ydown),0.5)):(Ydown + (Yup-Ydown)*0.5)));
      minline->SetLineStyle(2); minline->SetLineColor(1); minline->SetLineWidth(3);
      frame.addObject(minline);
    }
    double varMax = ws.var(varName.c_str())->getMax(Form("%sWindow", varName.c_str()));
    auto maxline = new TLine(varMax, 0.0, varMax, (setLogScale?(Ydown*TMath::Power((Yup/Ydown),0.5)):(Ydown + (Yup-Ydown)*0.5)));
    maxline->SetLineStyle(2); maxline->SetLineColor(1); maxline->SetLineWidth(3);
    frame.addObject(maxline);
  }
  //
  return;
};


void addGoFToWS(RooWorkspace& ws, const std::string& parName, const double& parValue)
{
  if (ws.var(parName.c_str())) { ws.var(parName.c_str())->setVal(parValue);               }
  else                         { ws.factory(Form("%s[%.6f]", parName.c_str(), parValue)); }
};


bool printGoF(TPad& pad, RooWorkspace& ws, const RooPlot& frame, const string& varLabel, const string& dataLabel, const string& pdfLabel)
{
  //
  if (ws.data(dataLabel.c_str())==NULL && ws.var(Form("testStat_BCChi2_%s", varLabel.c_str()))!=NULL) {
    const int& ndof_BCChi2 = int(ws.var(Form("ndofc_BCChi2_%s", varLabel.c_str()))->getVal());
    const double& testStat_BCChi2 = ws.var(Form("testStat_BCChi2_%s", varLabel.c_str()))->getVal();
    TLatex t = TLatex(); t.SetNDC(); t.SetTextSize(0.15);
    t.DrawLatex(0.72, 0.83, Form("#chi^{2}/ndof = %.0f / %d ", testStat_BCChi2, ndof_BCChi2));
    return true;
  }
  //
  if (ws.data(dataLabel.c_str())==NULL) { std::cout << "[ERROR] Dataset " << dataLabel << " was not found!" << std::endl; return false; }
  if (ws.pdf (pdfLabel.c_str() )==NULL) { std::cout << "[ERROR] PDF "     << pdfLabel  << " was not found!" << std::endl; return false; }
  const auto& dataP = dynamic_cast<RooDataSet*>(ws.data(dataLabel.c_str()));
  const auto& pdfP  = dynamic_cast<RooAbsPdf* >(ws.pdf (pdfLabel.c_str() ));
  const auto& varP  = dynamic_cast<RooRealVar*>(ws.var (varLabel.c_str() ));
  //
  // Find curve object
  RooCurve* rFit = frame.getCurve(Form("plot_%s", pdfLabel.c_str()));
  if (!rFit) { std::cout << "[ERROR] The latest RooCurve was not found" << std::endl; return false; }
  // Find histogram object
  RooHist* rData = frame.getHist(Form("plot_Tot%s", dataLabel.c_str()));
  if (!rData) { std::cout << "[ERROR] The latest RooHist was not found" << std::endl; return false; }
  //
  pad.cd();
  //
  // Unbinned Goodness-of-Fit tests
  //
  RooFit::RooGoF GoF_Unbinned(dataP, pdfP, varP);
  GoF_Unbinned.setRange(varP->getMin(), varP->getMax());
  //GoF_Unbinned.setNtoys(100, true, RooFit::Extended(kTRUE), RooFit::Strategy(2), RooFit::NumCPU(32));
  //
  // Kolmogorov-Smirnov test
  double pvalue_KS = -1., testStat_KS = -1.;
  GoF_Unbinned.KSTest(pvalue_KS, testStat_KS);
  if (testStat_KS>=0.0) {
    std::cout << "[INFO] Using Kolmogorov-Smirnov test gives result " << testStat_KS << " ( " << pvalue_KS << " ) " << std::endl;
    addGoFToWS(ws, Form("testStat_KS_%s", varLabel.c_str()), testStat_KS);
    addGoFToWS(ws, Form("pvalue_KS_%s"  , varLabel.c_str()), pvalue_KS  );
  }
  //
  // Anderson-Darling test
  double pvalue_AD = -1., testStat_AD = -1.;
  GoF_Unbinned.ADTest(pvalue_AD, testStat_AD);
  if (testStat_AD>=0.0) {
    std::cout << "[INFO] Using Anderson-Darling test gives result " << testStat_AD << " ( " << pvalue_AD << " ) " << std::endl;
    addGoFToWS(ws, Form("testStat_AD_%s", varLabel.c_str()), testStat_AD);
    addGoFToWS(ws, Form("pvalue_AD_%s"  , varLabel.c_str()), pvalue_AD  );
  }
  //
  // Binned Goodness-of-Fit tests
  //
  RooFit::RooGoF GoF_Binned(rData, rFit);
  GoF_Binned.setRange(varP->getMin(), varP->getMax());
  GoF_Binned.setRebin(5, false); // We use 5 in approval plots
  //
  // Determine the number of free parameters
  auto parList = std::unique_ptr<RooArgSet>(ws.pdf(pdfLabel.c_str())->getParameters(*ws.data(dataLabel.c_str())));
  const int& nFitPar = parList->selectByAttrib("Constant", kFALSE)->getSize();
  //
  // Baker-Cousins chi2 test
  int ndof_BCChi2 = -1;
  double pvalue_BCChi2 = -1., testStat_BCChi2 = -1.;
  GoF_Binned.BCChi2Test(pvalue_BCChi2, testStat_BCChi2, ndof_BCChi2, nFitPar);
  if (ndof_BCChi2>=0.0) {
    std::cout << "[INFO] Using Baker-Cousins chi2 test gives result " << testStat_BCChi2 << "/" << ndof_BCChi2 << " ( " << pvalue_BCChi2 << " ) " << std::endl;
    addGoFToWS(ws, Form("testStat_BCChi2_%s", varLabel.c_str()), testStat_BCChi2);
    addGoFToWS(ws, Form("pvalue_BCChi2_%s"  , varLabel.c_str()), pvalue_BCChi2  );
    addGoFToWS(ws, Form("ndofc_BCChi2_%s"   , varLabel.c_str()), ndof_BCChi2    );
  }
  //
  // Pearson chi2 test
  int ndof_PChi2 = -1;
  double pvalue_PChi2 = -1., testStat_PChi2 = -1.;
  GoF_Binned.PearsonChi2Test(pvalue_PChi2, testStat_PChi2, ndof_PChi2, nFitPar);
  if (ndof_PChi2>=0.0) {
    std::cout << "[INFO] Using Pearson chi2 test gives result " << testStat_PChi2 << "/" << ndof_PChi2 << " ( " << pvalue_PChi2 << " ) " << std::endl;
    addGoFToWS(ws, Form("testStat_PChi2_%s", varLabel.c_str()), testStat_PChi2);
    addGoFToWS(ws, Form("pvalue_PChi2_%s"  , varLabel.c_str()), pvalue_PChi2  );
    addGoFToWS(ws, Form("ndofc_PChi2_%s"   , varLabel.c_str()), ndof_PChi2    );
  }
  //
  // Default RooFit chi2 test (NOT RECOMMENDED)
  int ndof_RooFitChi2 = -1;
  double pvalue_RooFitChi2 = -1., testStat_RooFitChi2 = -1.;
  GoF_Binned.RooFitChi2Test(pvalue_RooFitChi2, testStat_RooFitChi2, ndof_RooFitChi2, nFitPar);
  if (ndof_RooFitChi2>=0.0) {
    std::cout << "[INFO] Using RooFit chi2 test gives result " << testStat_RooFitChi2 << "/" << ndof_RooFitChi2 << " ( " << pvalue_RooFitChi2 << " ) " << std::endl;
    addGoFToWS(ws, Form("testStat_RooFitChi2_%s", varLabel.c_str()), testStat_RooFitChi2);
    addGoFToWS(ws, Form("pvalue_RooFitChi2_%s"  , varLabel.c_str()), pvalue_RooFitChi2  );
    addGoFToWS(ws, Form("ndofc_RoofitChi2_%s"   , varLabel.c_str()), ndof_RooFitChi2    );
  }
  //
  TLatex t = TLatex(); t.SetNDC(); t.SetTextSize(0.15);
  if (ndof_BCChi2>0.) { t.DrawLatex(0.72, 0.83, Form("#chi^{2}/ndof = %.0f / %d ", testStat_BCChi2, ndof_BCChi2)); }
  //
  return true;
};


void divideBand(RooCurve& band, const RooCurve& cDen)
{
  double x, y;
  const auto& cNum = band;
  for (int i = 0; i < band.GetN(); i++) {
    band.GetPoint(i, x, y);
    if (cDen.interpolate(x) > 0.) {
      band.SetPoint(i, x, (y/cDen.interpolate(x)) );
      std::cout << i << "  " << x << "  " << y << "  " << cDen.interpolate(x) << std::endl;
    }
  }
};


bool addRatioBand(RooPlot& outFrame, RooPlot& frame, const RooWorkspace& ws, const std::string& pdfName, const std::string& dsName,
		  const std::string& varRange, const std::string& cName = "cenCurve", const double& Z = 1.0)
{
  // Get normalization
  const double& norm = ws.data(dsName.c_str())->sumEntries();
  // Find central curve
  ws.pdf(pdfName.c_str())->plotOn(&frame, RooFit::Name(cName.c_str()), RooFit::Range(varRange.c_str()), RooFit::NormRange(varRange.c_str()),
                                  RooFit::Normalization(norm, RooAbsReal::NumEvent)
                                  );
  const auto& rFit = dynamic_cast<RooCurve*>(frame.findObject(cName.c_str()));
  if (!rFit) { std::cout << "[ERROR] addRatioBand("<<cName<<") cannot find curve" << std::endl; return false; }
  RooCurve cenCurve = *rFit;
  frame.remove(0, kFALSE);
  // Get the fit results
  const auto& fitResult = dynamic_cast<RooFitResult*>(ws.obj(Form("fitResult_%s", pdfName.c_str())));
  // Create the Error Band
  ws.pdf(pdfName.c_str())->plotOn(&frame, RooFit::Name("band"), RooFit::Range(varRange.c_str()), RooFit::NormRange(varRange.c_str()),
                                  RooFit::Normalization(norm, RooAbsReal::NumEvent), RooFit::VisualizeError(*fitResult, Z, true)
                                  );
  // Find band
  const auto& band = dynamic_cast<RooCurve*>(frame.getCurve("band"));
  if (!band) { std::cout << "[ERROR] addRatioBand(band) cannot find curve" << std::endl; return false; }
  frame.remove(0, kFALSE);
  // Normalize the cenCurve
  divideBand(*band, cenCurve);
  // Add to frame
  outFrame.addPlotable(band, "CF2");
  outFrame.getAttFill()->SetFillColor(kOrange);
  // Return
  return true;
};


bool makePullHist(RooHist& pHist, const RooPlot& frame, const std::string& histname, const std::string& curvename, const bool& useAverage)
{
  // Find curve object
  const auto& rFit = dynamic_cast<RooCurve*>(frame.findObject(((curvename=="") ? 0 : curvename.c_str()), RooCurve::Class()));
  if (!rFit) { std::cout << "[ERROR] makePullHist(" << curvename << ") cannot find curve" << std::endl; return false; }
  // Find histogram object
  const auto& rData = dynamic_cast<RooHist*>(frame.findObject(((histname=="") ? 0 : histname.c_str()), RooHist::Class()));
  if (!rData) { std::cout << "[ERROR] makePullHist(" << histname  << ") cannot find histogram" << std::endl; return false; }
  // Determine range of curve
  double xstart, xstop, yDummy;
  rFit->GetPoint(0, xstart, yDummy);
  rFit->GetPoint((rFit->GetN()-1), xstop, yDummy);
  // Add histograms, calculate Poisson confidence interval on sum value
  for (int i = 0; i < rData->GetN(); i++) {
    // Get Data Value and Error
    double x, dataVal, exl, exh, eyl, eyh;
    getPoint(*rData, i, x, dataVal, exl, exh, eyl, eyh);
    // Only calculate pull for bins inside curve range
    if (x<xstart || x>xstop) continue;
    // Get Fit Value
    double fitVal = 0.0;
    if (useAverage) { fitVal = rFit->average(x-exl, x+exh); }
    else            { fitVal = rFit->interpolate(x);        }
    // Compute the Residual
    double y = (dataVal - fitVal);
    // Get the Norm factor
    const double& norm = ( (y>0) ? eyl : eyh );
    // Set Values
    if (norm>0.0) {
      y   /= norm;
      eyl /= norm;
      eyh /= norm;
      pHist.addBinWithXYError(x, y, exl, exh, eyl, eyh);
    }
    else {
      pHist.addBinWithXYError(x, 0, exl, exh, 0, 0);
    }
  }
  return true;
};


bool makeRatioHist(RooHist& rHist, const RooPlot& frame, const std::string& histname, const std::string& curvename, const bool& useAverage)
{
  // Find curve object
  const auto& rFit = dynamic_cast<RooCurve*>(frame.findObject(((curvename=="") ? 0 : curvename.c_str()), RooCurve::Class()));
  if (!rFit) { std::cout << "[ERROR] makeRatioHist(" << curvename << ") cannot find curve" << std::endl; return false; }
  // Find histogram object
  const auto& rData = dynamic_cast<RooHist*>(frame.findObject(((histname=="") ? 0 : histname.c_str()), RooHist::Class()));
  if (!rData) { std::cout << "[ERROR] makeRatioHist(" << histname  << ") cannot find histogram" << std::endl; return false; }
  // Determine range of curve
  double xstart, xstop, yDummy;
  rFit->GetPoint(0, xstart, yDummy);
  rFit->GetPoint((rFit->GetN()-1), xstop, yDummy);
  // Add histograms, calculate Poisson confidence interval on sum value
  for (int i = 0; i < rData->GetN(); i++) {
    // Get Data Value and Error
    double x, dataVal, exl, exh, eyl, eyh;
    getPoint(*rData, i, x, dataVal, exl, exh, eyl, eyh);
    // Only calculate pull for bins inside curve range
    if (x<xstart || x>xstop) continue;
    // Get Fit Value
    double fitVal = 0.0;
    if (useAverage) { fitVal = rFit->average(x-exl, x+exh); }
    else            { fitVal = rFit->interpolate(x);        }
    // Set Values
    if (fitVal>0.0 && dataVal>0.0) {
      const double& y = (dataVal / fitVal);
      eyl /= fitVal;
      eyh /= fitVal;
      rHist.addBinWithXYError(x, y, exl, exh, eyl, eyh);
    }
    else {
      rHist.addBinWithXYError(x, 0, exl, exh, 0, 0);
    }
  }
  return true;
};


double getErrorHi(const RooRealVar& var)
{
  if (var.getErrorLo()==0.0 && var.getErrorHi()==0.0) { return var.getError(); }
  return std::abs(var.getErrorHi());
};


double getErrorLo(const RooRealVar& var)
{
  if (var.getErrorLo()==0.0 && var.getErrorHi()==0.0) { return var.getError(); }
  return std::abs(var.getErrorLo());
};


bool isParAtLimit(const RooRealVar& var)
{
  if (
      ( (std::abs(var.getValV() - var.getMin())/getErrorLo(var)) <= 3.0 ) ||
      ( (std::abs(var.getValV() - var.getMax())/getErrorHi(var)) <= 3.0 )
      )
    { return true; }
  return false;
};


bool getVar(std::vector<RooRealVar>& varVec, const RooWorkspace& ws, const std::string& name, const std::string& pdfName)
{
  varVec.clear();
  auto parIt = std::unique_ptr<TIterator>(ws.allVars().selectByAttrib("Constant", kFALSE)->createIterator());
  for (auto itp = parIt->Next(); itp!=NULL; itp = parIt->Next() ) {
    const auto& it = dynamic_cast<RooRealVar*>(itp); if (!it) continue;
    std::string s(it->GetName());
    if((s.find("Pl")!=std::string::npos)!=(pdfName.find("Pl")!=std::string::npos)){ continue; }
    if((s.find("OS")!=std::string::npos)!=(pdfName.find("OS")!=std::string::npos)){ continue; }
    if (s.find(name)!=std::string::npos) { varVec.push_back(*it);}
  }
  if (varVec.size()>0) return true;
  auto fncIt = std::unique_ptr<TIterator>(ws.allFunctions().selectByAttrib("Constant", kFALSE)->createIterator());
  for (auto itp = fncIt->Next(); itp!=NULL; itp = fncIt->Next() ) {
    const auto& it = dynamic_cast<RooRealVar*>(itp); if (!it) continue;
    std::string s(it->GetName());
    if((s.find("Pl")!=std::string::npos)!=(pdfName.find("Pl")!=std::string::npos)){ continue; }
    if((s.find("OS")!=std::string::npos)!=(pdfName.find("OS")!=std::string::npos)){ continue; }
    if (s.find(name)!=std::string::npos) { varVec.push_back(*it); }
  }
  return (varVec.size()>0);
};


void parseVarName(const std::string& name, std::string& label)
{
  label = "";
  // Parse the parameter's labels
  stringstream ss(name); std::string s1, s2, s3;
  getline(ss, s1, '_'); getline(ss, s2, '_'); getline(ss, s3, '_');
  // Format model parameters
  if (s1=="Alpha"){ s1="#alpha"; } else if (s1=="Beta"){ s1="#beta"; }
  else if (s1=="Sigma0"){ s1="#sigma0"; } else if (s1=="Sigma1"){ s1="#sigma1"; } else if (s1=="Sigma2"){ s1="#sigma2"; } else if (s1=="rSigma21"){ s1="#frac{#sigma2}{#sigma1}"; }
  else if (s1=="Lambda1"){ s1="#lambda1"; } else if (s1=="Lambda2"){ s1="#lambda2"; } else if (s1=="Lambda3"){ s1="#lambda3"; }
  else if (s1=="Lambda4"){ s1="#lambda4"; } else if (s1=="Lambda5"){ s1="#lambda5"; } else if (s1=="Lambda6"){ s1="#lambda6"; }
  else if (s1=="XSection"){ s1="#sigma"; } else if (s1=="AccXEff"){ s1="#alphax#epsilon"; }
  // Format Object name
  std::string chg = ""; if (s2.find("Pl")!=std::string::npos) { chg = "+"; } else if (s2.find("Mi")!=std::string::npos) { chg = "-"; }
  if (s2.find("WToTau")!=std::string::npos) { s2 = "W#rightarrow#tau"; } else if (s2.find("WW")!=std::string::npos) { s2 = "WW"; }
  else if (s2.find("WZ")!=std::string::npos) { s2 = "WZ"; } else if (s2.find("W")!=std::string::npos) { s2 = "W#rightarrow#mu"; }
  else if (s2.find("DYToTau")!=std::string::npos) { s2 = "Z/#gamma*#rightarrow#tau"; }
  else if (s2.find("DY")!=std::string::npos) { s2 = "Z/#gamma*"; } else if (s2.find("Z")!=std::string::npos) { s2 = "Z"; }
  else if (s2.find("QCD")!=std::string::npos) { s2 = "QCD"; } else if (s2.find("TTbar")!=std::string::npos) { s2 = "t#bar{t}"; }
  else if (s2.find("JPsi")!=std::string::npos) { s2 = "J/#psi"; } else if (s2.find("Psi2S")!=std::string::npos) { s2 = "#psi(2S)"; }
  else if (s2.find("Ups1S")!=std::string::npos) { s2 = "#Upsilon(1S)"; } else if (s2.find("Ups2S")!=std::string::npos) { s2 = "#Upsilon(2S)"; }
  else if (s2.find("Ups3S")!=std::string::npos) { s2 = "#Upsilon(3S)"; } else if (s2.find("Bkg")!=std::string::npos) { s2 = "Bkg"; }
  s2 = ( s2 + chg );
  if(s3!=""){ label = Form("%s_{%s}^{%s}", s1.c_str(), s2.c_str(), s3.c_str()); } else { label = Form("%s^{%s}", s1.c_str(), s2.c_str()); }
  return;
};


#endif // #ifndef drawUtils_h
