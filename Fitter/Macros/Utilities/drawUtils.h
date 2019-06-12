#ifndef drawUtils_h
#define drawUtils_h


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
#include "RooExtendPdf.h"
#include "RooRealVar.h"
#include "RooStringVar.h"
#include "RooCurve.h"
#include "RooHist.h"
#include "RooHistPdf.h"
#include "RooFitResult.h"

#include "../../../Utilities/CMS/tdrstyle.C"
#include "../../../Utilities/CMS/CMS_lumi.C"
#include "../../../Utilities/RooGoF.C"
#include "../../../Utilities/dataUtils.h"

#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <vector>


typedef std::vector<std::string > StringVector_d;
typedef std::map< std::string , std::string > StringMap_d;
typedef std::map< std::string , StringMap_d > StringDiMap_d;
typedef std::map< std::string , TPad*       > PadPtrMap_d; // Unique Pointer does produce Segmentation Fault, so don't use it
typedef std::map< std::string , std::unique_ptr<RooPlot> > RooPlotPtrMap_d;


const std::map< std::string , std::vector<int> > PDFMAP_ =
  {
   {"JPsi", {11, kRed+1}}, {"Ups1S", {10, kGreen+3}}, {"Ups2S", {9, kGreen-5}}, {"Psi2S", {8, kRed-7}},
   {"Ups3S", {7, kGreen-6}}, {"DY", {6, kBlue+2}}, {"Z", {5, kBlue-4}}, {"TTbar", {4, kOrange+1}},
   {"D0", {3, kViolet+9}}, {"D0Swap", {2, kViolet-6}}, {"Bkg", {1, kAzure-9}}
  };


bool changeToHistPdf(RooWorkspace& ws, std::string& pdfName, const std::string& varName)
{
  if (pdfName.find("pdfHIST")!=std::string::npos) { return true; }
  const auto& oPDF = ws.pdf(pdfName.c_str());
  if (!oPDF) { std::cout << "[ERROR] changeToHistPdf: PDF " << pdfName << " was not found!" << std::endl; return false; }
  const auto& aPDF = dynamic_cast<RooAddPdf*>(oPDF);
  if (!aPDF) { std::cout << "[ERROR] changeToHistPdf: PDF " << pdfName << " does not derive from RooAddPdf!" << std::endl; return false; }
  const auto& var = ws.var(varName.c_str());
  if (!var) { std::cout << "[ERROR] changeToHistPdf: Variable " << varName << " was not found!" << std::endl; return false; }
  // Loop over the sub PDFs
  RooArgList pdfs;
  auto pdfIt = std::unique_ptr<TIterator>(aPDF->pdfList().createIterator());
  for (auto itp = pdfIt->Next(); itp!=NULL; itp = pdfIt->Next()) {
    const auto& it = dynamic_cast<RooAbsPdf*>(itp); if (!it) continue;
    // Check if PDF is already a RooHistPdf
    if (dynamic_cast<RooHistPdf*>(itp)){ pdfs.add(*it); continue; }
    // Define the Hist PDF name
    std::string hPDFName = it->GetName(); stringReplace(hPDFName, "pdf", "pdfHIST");
    std::string hDSName = hPDFName; stringReplace(hDSName, "pdf", "ds");
    // Create the Hist PDF
    if (!ws.data(hDSName.c_str())) {
      auto data = std::unique_ptr<RooDataHist>(it->generateBinned(RooArgSet(*var), 1000000, RooFit::Name(hDSName.c_str()), RooFit::ExpectedData()));
      if (!data) { std::cout << "[ERROR] changeToHistPdf: RooDataHist " << hDSName << " was not created!" << std::endl; return false; }
      if (ws.import(*data)) { return false; }
    }
    if (!ws.pdf(hPDFName.c_str())) {
      auto pdf = std::unique_ptr<RooHistPdf>(new RooHistPdf(hPDFName.c_str(), hPDFName.c_str(), *var, *dynamic_cast<RooDataHist*>(ws.data(hDSName.c_str()))));
      if (!pdf) { std::cout << "[ERROR] changeToHistPdf: RooHistPdf " << hPDFName << " was not created!" << std::endl; return false; }
      if (ws.import(*pdf)) { return false; }
      pdfs.add(*pdf);
    }
  }
  stringReplace(pdfName, "pdf", "pdfHIST");
  const auto& addPDF = RooAddPdf(pdfName.c_str(), pdfName.c_str(), pdfs);
  if (!ws.pdf(pdfName.c_str()) && ws.import(addPDF)) { return false; }
  return true;
};


bool addPdfToList(RooArgList& pdfs, RooArgList& yields, const RooWorkspace& ws, const RooArgList& pdfList)
{
  auto pdfIt = std::unique_ptr<TIterator>(pdfList.createIterator());
  for (auto it = pdfIt->Next(); it!=NULL; it = pdfIt->Next()) {
    std::string label = it->GetName(); label = label.substr(label.find("_")+1);
    std::string pTag = it->GetName(); pTag = pTag.substr(0, pTag.find("Tot_"));
    // Add the PDF
    const auto& pName = (ws.pdf((pTag+"BIN_"+label).c_str()) ? (pTag+"BIN_"+label) : (pTag+"_"+label));
    const auto& pdf = ws.pdf(pName.c_str());
    if (!pdf) { std::cout << "[ERROR] addPdfToList: PDF " << pName << " was not found!" << std::endl; return false; }
    pdfs.add(*pdf);
    // Add the yield
    const auto& yield = (ws.var(("N_"+label).c_str()) ? ws.var(("N_"+label).c_str()) : ws.function(("N_"+label).c_str()));
    if (!yield) { std::cout << "[ERROR] addPdfToList: Parameter N_" << label << " was not found!" << std::endl; return false; }
    yields.add(*yield);
  }
  return true;
};


void getLumiLabels(StringVector_d& labels, const std::string& PD, const std::string& col, const bool& isMC)
{
  std::string lumiLabel="";
  const auto& colV = parseColStr(col);
  if (!colV.empty()) {
    auto colN = colV[0];
    if (colN=="PP") { colN = "pp"; } else if (colN=="PA") { colN = "pPb"; }
    lumiLabel += colN;
  }
  if (isMC) { lumiLabel += " Simulation"; }
  else if (col=="PbPb5Y18") { lumiLabel += Form(" %.0f #mub^{-1}", PbPb::R5TeV::Y2018::LumiFromPD(PD)); }
  else if (col=="PP13Y18" ) { lumiLabel += Form(" %.1f pb^{-1}", pp::R13TeV::Y2018::LumiFromPD(PD)); }
  else if (col=="PP5Y17"  ) { lumiLabel += Form(" %.1f pb^{-1}", pp::R5TeV::Y2017::LumiFromPD(PD)); }
  else if (col.rfind("8Y16")!=std::string::npos) { lumiLabel += Form(" %.1f nb^{-1}", pPb::R8TeV::Y2016::LumiFromPD(PD, col)); }
  else if (col=="PbPb5Y15") { lumiLabel += Form(" %.0f #mub^{-1}", PbPb::R5TeV::Y2015::LumiFromPD(PD)); }
  labels.push_back(lumiLabel);
  //
  std::string lumiLabel2="";
  if (colV.size()>1) {
    auto colE = colV[1];
    if (colE=="5") { colE = "5.02 TeV"; } else if (colE=="8") { colE = "8.16 TeV"; } else if (colE=="13") { colE = "13 TeV"; } 
    colE = ((colV[0]=="PP") ? "#sqrt{s} = " : "#sqrt{s_{NN}} = ") +colE;
    lumiLabel2 += colE;
  }
  labels.push_back(lumiLabel2);
};


void formatLegendEntry(TLegendEntry& e, const double& size=0.060)
{
  e.SetTextSize(size);
};


void setVarToTag(RooWorkspace& ws, const std::string& varName, const std::string& tag="")
{
  const auto& var = ws.var(varName.c_str());
  if (!var || tag=="") return;
  // Extract the info
  const auto& varMin = var->getMin(tag.c_str());
  const auto& varMax = var->getMax(tag.c_str());
  const auto& varBin = var->getBins(tag.c_str());
  // Set the info
  var->setBins(varBin);
  var->setRange(varMin, varMax);
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


bool rooPlotToTH1(TH1D& hData, TH1D& hFit, const RooPlot& frame, const bool& useAverage = true, const int& sBinFit=1)
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
  std::vector<double> binDataV, binFitV;
  for (int i = 0; i < rData->GetN(); i++) {
    double x, y; rData->GetPoint(i, x, y);
    // Only consider bins inside curve range
    if (x<xstart || x>xstop) continue;
    const double& binW = rData->getNominalBinWidth();
    for (int j = 0; j < sBinFit; j++) { binFitV.push_back(x - (binW*0.5) + j*(binW/sBinFit)); }
    binDataV.push_back(x - (binW*0.5));
    if (i==(rData->GetN()-1)) { binFitV.push_back(x + (binW*0.5)); }
    if (i==(rData->GetN()-1)) { binDataV.push_back(x + (binW*0.5)); }
  }
  const uint& nBinFit = (binFitV.size()-1), nBinData = (binDataV.size()-1);
  double binFit[nBinFit+1], binData[nBinData+1];
  for (uint i = 0; i < binFitV.size(); i++) { binFit[i] = binFitV[i]; }
  for (uint i = 0; i < binDataV.size(); i++) { binData[i] = binDataV[i]; }
  //
  hData.Reset(); hData = TH1D(Form("hData_%s", rData->GetName()), rData->GetTitle(), nBinData, binData); hData.Sumw2();
  hFit.Reset();  hFit  = TH1D(Form("hFit_%s" , rFit->GetName()) , rFit->GetTitle() , nBinFit,  binFit);  hFit.Sumw2();
  // Set Histogram entries
  for (uint i = 0; i < nBinData; i++) {
    double x, dataVal, exl, exh, eyl, eyh;
    getPoint(*rData, i, x, dataVal, exl, exh, eyl, eyh);
    hData.SetBinContent((i+1), dataVal);
    hData.SetBinError((i+1), std::sqrt((eyl*eyl + eyh*eyh)/2.0));
    for (int j = 0; j < sBinFit; j++) {
      double fitVal = 0.0, binW = (exh+exl), x1 = x-exl+j*(binW/sBinFit), x2 = x-exl+(j+1)*(binW/sBinFit);
      if (useAverage) { fitVal = rFit->average(x1, x2);         }
      else            { fitVal = rFit->interpolate((x1+x2)/2.); }
      hFit.SetBinContent((sBinFit*i+j+1), fitVal);
      hFit.SetBinError((sBinFit*i+j+1), 0.0);
    }
  }
  return true;
};


void setPlotRange(RooPlot& frame, const RooWorkspace& ws, const std::string& varName, const std::string& dsName, const bool& setLogScale)
{
  // Find maximum and minimum points of Plot to rescale Y axis
  TH1D hData, hFit;
  if (!rooPlotToTH1(hData, hFit, frame, true, 4)) { std::cout << "[ERROR] Could not find the RooHist from the frame!" << std::endl; return; }
  double YMax = std::max(hData.GetBinContent(hData.GetMaximumBin()), hFit.GetBinContent(hFit.GetMaximumBin()));
  double YMin = 1e99;
  for (int i=1; i<=hData.GetNbinsX(); i++) if (hData.GetBinContent(i)>0) YMin = min(YMin, hData.GetBinContent(i));
  double Yup(0.), Ydown(0.), rDown=0.05, rUp=0.4;
  if(setLogScale)
  {
    YMin = std::max(YMin, 0.1); YMax = std::max(YMax, 0.1);
    Ydown = std::max(YMin/(std::pow(YMax/YMin, (rDown/(1.0-rUp-rDown)))), 0.1);
    Yup = Ydown*std::pow(YMax/Ydown, (1.0/(1.0 - rUp)));
  }
  else
  {
    Ydown = std::max(YMin - (rDown/(1.0-rUp-rDown))*(YMax-YMin), 0.0);
    Yup = Ydown + (YMax - Ydown)/(1.0 - rUp);
  }
  frame.GetYaxis()->SetRangeUser(Ydown, Yup);
  std::cout << "[INFO] Setting plot y-axis range to: " << Ydown << " - " << Yup << std::endl;
  //
  // Draw Lines for the var range if cut
  if (dsName.find("_FIT")!=std::string::npos) {
    const auto& varMin = ws.var(varName.c_str())->getMin("FitWindow");
    if (varMin > 0.0) {
      auto minline = new TLine(varMin, 0.0, varMin, (setLogScale?(Ydown*TMath::Power((Yup/Ydown),0.5)):(Ydown + (Yup-Ydown)*0.5)));
      minline->SetLineStyle(2); minline->SetLineColor(1); minline->SetLineWidth(3);
      frame.addObject(minline);
    }
    const auto& varMax = ws.var(varName.c_str())->getMax("FitWindow");
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
  RooHist* rData = frame.getHist(Form("plot_%s", dataLabel.c_str()));
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
  auto varList = std::unique_ptr<RooAbsCollection>(parList->selectByAttrib("Constant", false));
  const int& nFitPar = varList->getSize();
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
    const double norm = ( (y>0) ? eyl : eyh );
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


RooAbsCollection* skimSet(const RooArgSet& set, const std::string& varStr="*", const StringVector_d& excV={}, const std::string& state="Constant")
{
  // Select variables passing varStr
  auto vSet = std::unique_ptr<RooAbsCollection>(set.selectByName(varStr.c_str(), true));
  // Remove variables passing state
  const auto& sSet = vSet->selectByAttrib(state.c_str(), false);
  // Remove variables in excV
  for (const auto& o : excV) { const auto& p = sSet->find(o.c_str()); if (p) { sSet->remove(*p); } }
  return sSet;
};
  

std::vector<RooRealVar> getModelVar(const RooWorkspace& ws, const std::string& name="*")
{
  // Get information from workspace
  const std::string& cha     = (ws.obj("channel") ? dynamic_cast<RooStringVar*>(ws.obj("channel"))->getVal() : "");
  const std::string& pdfName = (ws.obj("pdfName") ? dynamic_cast<RooStringVar*>(ws.obj("pdfName"))->getVal() : "");
  const std::string& dsName  = (ws.obj("dsName")  ? dynamic_cast<RooStringVar*>(ws.obj("dsName") )->getVal() : "");
  // Define the var name
  const auto& label = (pdfName.find("_")!=std::string::npos ? pdfName.substr(pdfName.find(cha)) : "");
  const auto& varStr = name + "_*" + (label=="" ? "" : label);
  // Define the list of observables
  StringVector_d obsV;
  if (ws.data(dsName.c_str())) {
    const  auto& listObs = const_cast<RooWorkspace*>(&ws)->set(("SET_"+dsName).c_str());
    auto parIt = std::unique_ptr<TIterator>(listObs->createIterator());
    for (auto it = parIt->Next(); it!=NULL; it = parIt->Next() ) { obsV.push_back(it->GetName()); }
  }
  // Get the var set
  auto parSet = std::unique_ptr<RooArgSet>(ws.pdf(pdfName.c_str()) ? ws.pdf(pdfName.c_str())->getParameters(RooArgSet()) : NULL);
  auto varSet = std::unique_ptr<RooAbsCollection>(skimSet((parSet ? *parSet : ws.allVars()), varStr, obsV));
  if (varSet->getSize()==0) { varSet.reset(skimSet(ws.allFunctions(), varStr, obsV)); }
  // Fill the vector of variables
  std::vector<RooRealVar> varV;
  auto parIt = std::unique_ptr<TIterator>(varSet->createIterator());
  for (auto itp = parIt->Next(); itp!=NULL; itp = parIt->Next()) {
    const auto& it = dynamic_cast<RooRealVar*>(itp); if (!it) continue;
    varV.push_back(*it);
  }
  // Set yields at the beginning and order based on object
  std::vector<RooRealVar> varVec;
  for (const auto& v : varV) { if (std::string(v.GetName()).rfind("N_",0)==0) { varVec.push_back(v); } }
  for (const auto& v : varV) { if (std::string(v.GetName()).rfind("R_",0)==0) { varVec.push_back(v); } }
  for (const auto& p : PDFMAP_) {
    if (p.first=="Bkg") continue;
    for (const auto& v : varV) { const std::string& vv = v.GetName(); if (vv.rfind("N_",0)!=0 && vv.rfind("R_",0)!=0 && vv.find(p.first)!=std::string::npos) { varVec.push_back(v); } }
  }
  for (const auto& v : varV) { const std::string& vv = v.GetName(); if (vv.rfind("N_",0)!=0 && vv.find("Bkg")!=std::string::npos) { varVec.push_back(v); } }
  return varVec;
};


std::string formatCut(const std::string& cut)
{
  if (cut=="") return cut;
  std::string str = cut;
  str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
  // Format the variable string
  const auto& varL = StringMap_d({{"Cand_Pt", "p_{T}"}, {"Cand_Mass", "M"}});
  for (const auto& vL : varL) { stringReplace(str, vL.first, vL.second); }
  // Format the observables
  for (const auto& elem : VARLABEL_) { stringReplace(str, elem.first, elem.second); }
  // Format the object string
  const auto& objL = StringMap_d({{"Pl", "^{+}+x"}, {"Mi", "^{#font[122]{\55}}+x"}, {"To", "#rightarrow"},
				  {"El", "e"}, {"Mu", "#mu"}, {"Tau","#tau"},
				  {"DY", "Z/#gamma*"}, {"TTbar", "t#bar{t}"},
				  {"JPsi", "J/#psi"}, {"Psi(2S)", "#psi(2S)"},
				  {"Ups(1S)", "#Upsilon(1S)"}, {"Ups(2S)", "#Upsilon(2S)"}, {"Ups(3S)", "#Upsilon(3S)"},
				  {"&&", " & "}, {"(", ""}, {")", ""}, {"<=", "#leq"}, {"<", " < "}, {">=", " #geq "}, {">", " > "}});
  for (const auto& oL : objL) { stringReplace(str, oL.first, oL.second); }
  return str;
};


std::string formatPar(const std::string& name, const std::string& cha)
{
  std::string label = cha;
  if (label=="ToMuMu") { label = "#mu^{+}#mu^{#font[122]{\55}} "; }
  const auto& parL = StringMap_d({{"Cand_Mass", "Mass GeV/c^{2}"}, {"Cand_Pt", "p_{T} GeV/c"}});
  for (const auto& pL : parL) { if (name==pL.first) { label += pL.second; break; } }
  return label;
};


std::string parseVarName(const std::string& name)
{
  std::string label = "";
  // Parse the parameter's labels
  stringstream ss(name); std::string s1, s2, s3;
  getline(ss, s1, '_'); getline(ss, s2, '_'); getline(ss, s3, '_');
  // Format model parameters
  const auto& varL = StringMap_d({{"Alpha", "#alpha"}, {"Beta", "#beta"}, {"Lambda", "#lambda"}, {"Sigma", "#sigma"}, {"rSigma21","#sigma2/#sigma1"},
				  {"XSection","#sigma"}, {"AccXEff","#alphax#epsilon"}});
  for (const auto& vL : varL) { if (s1.rfind(vL.first, 0)==0) { stringReplace(s1, vL.first, vL.second); break; } }
  // Format Object name
  const auto& objL = StringMap_d({{"Bkg", "Bkg"}, {"Z", "Z"}, {"DY", "Z/#gamma*"}, {"TTbar", "t#bar{t}"},
				  {"JPsi", "J/#psi"}, {"Psi2S", "#psi(2S)"},
				  {"Ups1S", "#Upsilon(1S)"}, {"Ups2S", "#Upsilon(2S)"}, {"Ups3S", "#Upsilon(3S)"}});
  for (const auto& oL : objL) { if (s2.find(oL.first)!=std::string::npos) { s2 = oL.second; break; } }
  // Format System name
  const auto& s3V = parseColStr(s3);
  auto col = (!s3V.empty() ? s3V[0] : ""); if (col=="PA") { col = "pA"; } else if (col=="PP") { col = "pp"; }
  if (col!="") { label = Form("%s_{%s}^{%s}", s1.c_str(), s2.c_str(), col.c_str()); } else { label = Form("%s^{%s}", s1.c_str(), s2.c_str()); }
  return label;
};


std::string parseObject(const std::string& inObj)
{
  // Format Object name
  const auto& objL = StringMap_d({{"DY", "Z/#gamma*"}, {"TTbar", "t#bar{t}"},
				  {"PsinS", "#psi(nS)"}, {"JPsi", "J/#psi"}, {"Psi2S", "#psi(2S)"},
				  {"UpsnS", "#Upsilon(nS)"}, {"Ups1S", "#Upsilon(1S)"}, {"Ups2S", "#Upsilon(2S)"}, {"Ups3S", "#Upsilon(3S)"}});
  return (objL.find(inObj)!=objL.end() ? objL.at(inObj) : inObj);
};


std::string parseProcess(const std::string& obj, const std::string& cha)
{
  // Format Object name
  std::string proc = obj;
  if (proc=="") return proc;
  const auto& objL = StringMap_d({{"To", " #rightarrow "}, {"El", "e"}, {"Mu", "#mu"}, {"Tau","#tau"},
				  {"DY", "Z/#gamma*"}, {"TTbar", "t#bar{t}"},
				  {"PsinS", "#psi(nS)"}, {"JPsi", "J/#psi"}, {"Psi2S", "#psi(2S)"},
				  {"UpsnS", "#Upsilon(nS)"}, {"Ups1S", "#Upsilon(1S)"}, {"Ups2S", "#Upsilon(2S)"}, {"Ups3S", "#Upsilon(3S)"}});
  for (const auto& oL : objL) { stringReplace(proc, oL.first, oL.second); }
  // Format the channel
  proc += cha;
  const auto& chaL = StringMap_d({{"ToMuMu", " #rightarrow #mu^{+} + #mu^{#font[122]{\55}}"}});
  for (const auto& cL : chaL) { stringReplace(proc, cL.first, cL.second); }
  proc = ("#font[62]{#scale[1.1]{"+proc+"}}");
  return proc;
};


std::string parseProcess(const std::set<std::string>& objV, const std::string& cha)
{
  std::string obj="";
  if (objV.size()==1) { obj = *objV.begin(); }
  else if (objV.find("Ups1S")!=objV.end()) { obj = "UpsnS"; }
  else if (objV.find("JPsi")!=objV.end()) { obj = "PsinS"; }
  return parseProcess(obj, cha);
};


#endif // #ifndef drawUtils_h
