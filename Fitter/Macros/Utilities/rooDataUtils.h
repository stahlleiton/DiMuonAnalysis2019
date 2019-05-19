#ifndef rooDataUtils_h
#define rooDataUtils_h

#include "TSystem.h"
#include "TFile.h"
#include "TObject.h"
#include "TObjString.h"
#include "TIterator.h"
#include "TH1.h"

#include "RooFit.h"
#include "RooMsgService.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooAbsPdf.h"
#include "RooStringVar.h"

#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <map>

#include "initClasses.h"
#include "../../../Utilities/dataUtils.h"


std::string formatCut(const std::string& cut, const StringMap_t& map = StringMap_t())
{
  std::string str = cut;
  str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
  stringReplace( str, "<=Cand_Pt", " GeV/c<=DiMuon_Pt" ); stringReplace( str, "<Cand_Pt", " GeV/c<Cand_Pt" );
  stringReplace( str, "<=Cand_Mass", " GeV/c^{2}<=DiMuon_Mass" ); stringReplace( str, "<Cand_Mass", " GeV/c^{2}<Cand_Mass" );
  for (auto& elem : map) { stringReplace( str, elem.first, elem.second ); }
  stringReplace( str, "Pl", "^{+}+x" ); stringReplace( str, "Mi", "^{-}+x" );
  stringReplace( str, "Mu", "#mu" ); stringReplace( str, "Tau", "#tau" ); stringReplace( str, "DY", "Z/#gamma*" ); stringReplace( str, "TTbar", "t#bar{t}" );
  stringReplace( str, "Ups(1S)", "#Upsilon(1S)" ); stringReplace( str, "Ups(2S)", "#Upsilon(2S)" ); stringReplace( str, "Ups(3S)", "#Upsilon(3S)" );
  stringReplace( str, "JPsi", "J/#psi" ); stringReplace( str, "Psi(2S)", "#psi(2S)" );
  stringReplace( str, "To", "#rightarrow" ); stringReplace( str, "&&", " & " ); stringReplace( str, "(", "" ); stringReplace( str, ")", "" );
  stringReplace( str, "<=", " #leq " ); stringReplace( str, "<", " < " ); stringReplace( str, ">=", " #geq " ); stringReplace( str, ">", " > " );
  return str;
};

void getRange(std::vector<double>& range, const TH1D& hist, const int& nMaxBins)
{
  // 1) Find the bin with the maximum Y value
  const auto& binMaximum = hist.GetMaximumBin();
  // 2) Loop backward and find the first bin
  int binWithContent = -1;
  int firstBin = 1;
  for (int i = binMaximum; i > 0; i--) {
    if (hist.GetBinContent(i) > 0.0) {
      if ( (binWithContent > 0) && ((binWithContent-i) > nMaxBins) && (hist.GetBinContent(i) < 5.0) ) { firstBin = binWithContent; break; }
      else { binWithContent = i; }
    }
  }
  // 3) Loop forward and find the last bin
  binWithContent = -1;
  int lastBin = hist.GetNbinsX();
  for (int i = binMaximum; i < hist.GetNbinsX(); i++) {
    if (hist.GetBinContent(i) > 0.0) {
      if ( ( binWithContent > 0) && ((i - binWithContent) > nMaxBins) && (hist.GetBinContent(i) < 5.0) ) { lastBin = binWithContent+1; break; }
      else { binWithContent = i; }
    }
  }
  // 4) Build the set of bins
  const auto& startBin = ( (firstBin > 1) ? (firstBin - 1) : firstBin );
  const auto& nNewBins = lastBin - startBin + 1;
  double binning[nNewBins+2];
  binning[0] = hist.GetXaxis()->GetXmin();
  binning[nNewBins+1] = hist.GetXaxis()->GetXmax();
  for (int i = 1; i <= nNewBins; i++) {
    int iBin = startBin + i;
    binning[i] = hist.GetBinLowEdge(iBin);
  }
  // 5) Save the bin range
  range.push_back(binning[(firstBin>1)?1:0]);
  range.push_back(binning[nNewBins]);
  //
  return;
};


void setFixedVarsToContantVars(RooWorkspace& ws)
{
  auto listVar = ws.allVars();
  std::unique_ptr<TIterator> parIt = std::unique_ptr<TIterator>(listVar.createIterator());
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    if ( it->getMin()==it->getMax() && !it->isConstant() ) {
      std::cout << "[INFO] Setting " << it->GetName() << " constant!" << std::endl;
      it->setConstant(kTRUE);
    }
  }
};


bool setConstant(RooWorkspace& myws, const string& parName, const bool& CONST)
{
  if (myws.var(parName.c_str())) { 
    myws.var(parName.c_str())->setConstant(CONST);
    if (CONST) { std::cout << "[INFO] Setting parameter " << parName << " : " << myws.var(parName.c_str())->getVal() << " to constant value!" << std::endl; }
  }
  else if (!myws.function(parName.c_str())) { 
    std::cout << "[ERROR] Parameter " << parName << " was not found!" << std::endl;
    return false;
  }
  return true;
};


void copyWorkspace(RooWorkspace& outWS, const RooWorkspace& inWS, const std::string& smp="All", const bool& addVar=true, const bool& addSnap=false)
{
  //
  const auto& level = RooMsgService::instance().globalKillBelow();
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  // Copy all category variables
  if (inWS.allCats().getSize()>0 && addVar) {
    auto listCat = inWS.allCats();
    std::unique_ptr<TIterator> catIt = std::unique_ptr<TIterator>(listCat.createIterator());
    for (RooCategory* it = (RooCategory*)catIt->Next(); it!=NULL; it = (RooCategory*)catIt->Next() ) { if (!outWS.cat(it->GetName())) { outWS.import(*it); } }
  }
  // Copy all category functions
  if (inWS.allCatFunctions().getSize()>0 && addVar) {
    auto listFnc = inWS.allCatFunctions();
    std::unique_ptr<TIterator> fncIt = std::unique_ptr<TIterator>(listFnc.createIterator());
    for (RooCategory* it = (RooCategory*)fncIt->Next(); it!=NULL; it = (RooCategory*)fncIt->Next() ) { if (!outWS.catfunc(it->GetName())) { outWS.import(*it); } }
  }
  // Copy all variables
  if (inWS.allVars().getSize()>0 && addVar) {
    auto listVar = inWS.allVars();
    std::unique_ptr<TIterator> parIt = std::unique_ptr<TIterator>(listVar.createIterator());
    for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) { if (!outWS.var(it->GetName())) { outWS.import(*it); } }
  }
  // Copy all functions
  if (inWS.allFunctions().getSize()>0 && addVar) {
    auto listFnc = inWS.allFunctions();
    std::unique_ptr<TIterator> fncIt = std::unique_ptr<TIterator>(listFnc.createIterator());
    for (RooRealVar* it = (RooRealVar*)fncIt->Next(); it!=NULL; it = (RooRealVar*)fncIt->Next() ) { if (!outWS.function(it->GetName())) { outWS.import(*it); } }
  }
  // Copy all PDFs
  if (inWS.allPdfs().getSize()>0) {
    auto listPdf = inWS.allPdfs();
    std::unique_ptr<TIterator> pdfIt = std::unique_ptr<TIterator>(listPdf.createIterator());
    for (RooAbsPdf* it = (RooAbsPdf*)pdfIt->Next(); it!=NULL; it = (RooAbsPdf*)pdfIt->Next() ) { if (!outWS.pdf(it->GetName())) { outWS.import(*it, RooFit::RecycleConflictNodes()); } }
  }
  // Copy all Datasets
  if (inWS.allData().size()>0 && smp!="") {
    auto listData = inWS.allData();
    for (const auto& it : listData) { if (!outWS.data(it->GetName()) && (smp=="All" || std::string(it->GetName()).find(smp)!=std::string::npos)) { outWS.import(*it); } }
  }
  // Copy all Embedded Datasets
  if (inWS.allEmbeddedData().size()>0 && smp!="") {
    auto listData = inWS.allEmbeddedData();
    for (const auto& it : listData) { if (!outWS.embeddedData(it->GetName()) && (smp=="All" || std::string(it->GetName()).find(smp)!=std::string::npos)) { outWS.import(*it, RooFit::Embedded(1)); } }
  }
  // Copy all Generic Objects
  if (inWS.allGenericObjects().size()>0) {
    auto listObj = inWS.allGenericObjects();
    for (const auto& it : listObj) { if (!outWS.genobj(it->GetTitle())) { outWS.import(*it, it->GetTitle()); } }
  }
  // Copy all snapshots variables
  if (addSnap) {
    std::vector<std::string> snapNames = { "initialParameters", "fittedParameters" };
    for (const auto& snapName : snapNames) {
      const auto& snap = inWS.getSnapshot(snapName.c_str());
      if (snap) { outWS.saveSnapshot(snapName.c_str(), *snap, kTRUE); }
    }
  }
  // Return the RooMessenger Level
  RooMsgService::instance().setGlobalKillBelow(level);
};


bool saveWorkSpace(const RooWorkspace& ws, const string& outputDir, const string& fileName, const bool& saveAll=true)
{
  // Save the workspace
  gSystem->mkdir(outputDir.c_str(), kTRUE);
  auto file = std::unique_ptr<TFile>(new TFile((outputDir+fileName).c_str(), "RECREATE"));
  if (!file || !file->IsOpen() || file->IsZombie()) {
    file->Close();
    std::cout << "[ERROR] Output root file with fit results could not be created!" << std::endl; return false;
  }
  else {
    file->cd();
    if (saveAll) { ws.Write("workspace"); }
    else {
      RooWorkspace tmpWS;
      copyWorkspace(tmpWS, ws, "", true, true);
      tmpWS.Write("workspace");
    }
    file->Write(); file->Close();
  }
  return true;
};


bool isCompatibleDataset(const RooDataSet& ds, const RooDataSet& ref, const bool& checkRange=true)
{
  // Check that the DataSets have the same number of events
  if (ds.numEntries()!=ref.numEntries()){ std::cout << "[ERROR] DataSet disagreement in number of events : "
                                                    << ds.GetName() << " (" << ds.numEntries() << ") and "
                                                    << ref.GetName() << " (" << ref.numEntries() << ") !" << std::endl; return false; }
  if (ds.sumEntries()!=ref.sumEntries()){ std::cout << "[ERROR] DataSet disagreement in sum of weights : "
                                                    << ds.GetName() << " (" << ds.sumEntries() << ") and "
                                                    << ref.GetName() << " (" << ref.sumEntries() << ") !" << std::endl; return false; }
  // Check that the input DataSet have the same variables and distributions as the reference
  const auto& listVar = ref.get();
  std::unique_ptr<TIterator> parIt = std::unique_ptr<TIterator>(listVar->createIterator());
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    if ( ds.get()->find(it->GetName()) == NULL ) { std::cout << "[ERROR] DataSet " << ds.GetName() << " does not contain the variable " << it->GetName() << " !" << std::endl; return false; }
    if (checkRange) {
      if ( it->getMin() != ((RooRealVar*)ds.get()->find(it->GetName()))->getMin() ) { std::cout << "[ERROR] " << it->GetName() << " Min Range disagreement : "
                                                                                                << ds.GetName() << " ( " << ((RooRealVar*)ds.get()->find(it->GetName()))->getMin() << " ) " << " and "
                                                                                                << ref.GetName() << " ( " << it->getMin() << " ) ! " << std::endl; return false; }
      if ( it->getMax() != ((RooRealVar*)ds.get()->find(it->GetName()))->getMax() ) { std::cout << "[ERROR] " << it->GetName() << " Max Range disagreement : "
                                                                                                << ds.GetName() << " ( " << ((RooRealVar*)ds.get()->find(it->GetName()))->getMax() << " ) " << " and "
                                                                                                << ref.GetName() << " ( " << it->getMax() << " ) ! " << std::endl; return false; }
    }
    for (uint i = 1; i < 5; i++) {
      if ( ref.moment(*it, i) != ds.moment(*it, i) ) { std::cout << "[ERROR] " << it->GetName() << " " << i << " Moment Value disagreement : "
                                                                 << ds.GetName() << " ( " << ds.moment(*it, i) << " ) " << " and "
                                                                 << ref.GetName() << " ( " << ref.moment(*it, i) << " ) " << " ! " << std::endl; return false; }
    }
  }
  // DataSets are compatible if they passed all tests
  std::cout << "[INFO] DataSets " << ds.GetName() << " and " << ref.GetName() << " are compatible!" << std::endl; return true;
};


bool compareSnapshots(const RooArgSet& pars1, const RooArgSet& pars2)
{
  std::unique_ptr<TIterator> parIt = std::unique_ptr<TIterator>(pars1.createIterator());
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    double val = pars2.getRealValue(it->GetName(),-1e99);
    const std::string& name = it->GetName();
    if ( (name=="Cand_Pt") || (name=="Cand_AbsRap") || (name=="Cand_Rap") || (name=="Cand_Mass") || (name=="Centrality") ) continue;
    if ( (name=="Dau1_InvBeta") || (name=="Dau1_P") || (name=="Dau2_InvBeta") || (name=="Dau1_P") ) continue;
    if (val==-1e99) return false;          // the parameter was not found!
    if (val != it->getVal()) return false;  // the parameter was found, but with a different value!
    if ( ((RooRealVar&)pars2[it->GetName()]).getMin() != it->getMin() ) return false;  // the parameter has different lower limit
    if ( ((RooRealVar&)pars2[it->GetName()]).getMax() != it->getMax() ) return false;  // the parameter has different upper limit
  }
  return true;
};


bool isFitAlreadyFound(const RooArgSet& newpars, const string& fileName, const string& pdfName)
{
  std::cout << "[INFO] Checking if fit was already done!" << std::endl;
  if (gSystem->AccessPathName(fileName.c_str())) {
    std::cout << "[INFO] FileName: " << fileName << " was not found" << std::endl;
    return false; // File was not found
  }
  auto file = std::unique_ptr<TFile>(new TFile(fileName.c_str()));
  if (!file || !file->IsOpen() || file->IsZombie()) return false;
  RooWorkspace* ws = (RooWorkspace*) file->Get("workspace");
  if (!ws) {
    std::cout << "[INFO] Workspace was not found" << std::endl;
    file->Close();
    return false;
  }
  const auto& params = ( (ws->getSnapshot("initialParameters")) ? ws->getSnapshot("initialParameters") :  ws->getSnapshot(Form("%s_parIni", pdfName.c_str())) );
  if (!params) {
    std::cout << "[INFO] Snapshot of initial parameters was not found!" << std::endl;
    file->Close();
    return false;
  }
  bool result = compareSnapshots(newpars, *params);
  file->Close();
  return result;
};


RooRealVar getVar(const RooArgSet& set, const std::string& varName)
{
  if (set.find(varName.c_str()) != NULL) {
    return *(RooRealVar*)set.find(varName.c_str());
  }
  std::cout << "[WARNING] Variable " << varName << " was not found in the set, return empty RooRealVar" << std::endl;
  return RooRealVar();
};


void updateParameterRange(RooWorkspace& myws, GlobalInfo&  info, const std::string& chg, const std::string& DSTAG, const std::string& var="Cand_Mass", const double& maxRange=-1.0)
{
  // Check if maxRange is the same as current range, otherwise return
  if (maxRange==myws.var(var.c_str())->getMax()) { return; }
  //
  const std::string& dsName = ( "d" + chg + "_" + DSTAG );
  //
  double varMin , varMax;
  varMin = myws.var(var.c_str())->getMin(); varMax = myws.var(var.c_str())->getMax();
  if (maxRange > 0.0) { varMax = maxRange; }
  else {
    myws.data(dsName.c_str())->getRange(*myws.var(var.c_str()), varMin, varMax);
  }
  const std::string& varFitRange = Form("(%g <= %s && %s < %g)", varMin, var.c_str(), var.c_str(), varMax);
  //
  if (myws.data(dsName.c_str())->reduce(varFitRange.c_str())->numEntries() <= myws.data(dsName.c_str())->numEntries()) {
    const auto& nBins = int( (varMax - varMin)/myws.var(var.c_str())->getBinWidth(0) );
    myws.var(var.c_str())->setRange(Form("%sWindow", var.c_str()), varMin, varMax);
    myws.var(var.c_str())->setBins(nBins, Form("%sWindow", var.c_str()));
    //
    auto dataToFit = std::unique_ptr<RooDataSet>((RooDataSet*)(myws.data(dsName.c_str())->reduce(varFitRange.c_str())->Clone((dsName+"_FIT").c_str())));
    myws.import(*dataToFit);
    info.Var.at("numEntries").at(chg) = myws.data((dsName+"_FIT").c_str())->sumEntries();
    //
    std::cout << Form("[INFO] %s range was updated to : %g <= %s < %g", var.c_str(), varMin, var.c_str(), varMax) << std::endl;
  }
  //
  return;
};


bool loadParameterRange(GlobalInfo&  info, const std::string& var, const string& fileName)
{
  if (gSystem->AccessPathName(fileName.c_str())) {
    std::cout << "[ERROR] File " << fileName << " was not found!" << std::endl;
    return false; // File was not found
  }
  auto file = std::unique_ptr<TFile>(new TFile(fileName.c_str()));
  if (!file || !file->IsOpen() || file->IsZombie()) return false;
  auto ws = (RooWorkspace*)file->Get("workspace");
  if (!ws) {
    std::cout << "[ERROR] Workspace was not found in: " << fileName << std::endl;
    file->Close();
    return false;
  }
  if (ws->var(var.c_str())) {
    info.Var.at(var).at("Min") = ws->var(var.c_str())->getMin();
    info.Var.at(var).at("Max") = ws->var(var.c_str())->getMax();
  }
  else {
    std::cout << "[ERROR] " << var << " ctauErr was not found!" << std::endl;
    file->Close();
    return false;
  }
  file->Close();
  return true;
};


bool loadYields(RooWorkspace& myws, const std::string& fileName, const std::string& dsName, const std::string& pdfName, const std::string& col="PbPb")
{
  if (gSystem->AccessPathName(fileName.c_str())) {
    std::cout << "[INFO] File " << fileName << " was not found!" << std::endl;
    return false; // File was not found
  }
  auto file = std::unique_ptr<TFile>(new TFile(fileName.c_str()));
  if (!file || !file->IsOpen() || file->IsZombie()) return false;
  auto ws = (RooWorkspace*)file->Get("workspace");
  if (!ws) {
    std::cout << "[INFO] Workspace was not found in: " << fileName << std::endl;
    file->Close();
    return false;
  }
  bool compDS = true;
  if (ws->data(dsName.c_str()) && myws.data(dsName.c_str())) {
    if (isCompatibleDataset(*(RooDataSet*)myws.data(dsName.c_str()), *(RooDataSet*)ws->data(dsName.c_str()))) {
      const auto& params = ws->getSnapshot(Form("%s_parIni", pdfName.c_str()));
      if (!params) {
        std::cout << "[INFO] Snapshot " << pdfName << "_parIni was not found!" << std::endl;
        file->Close();
        return false;
      }
      std::unique_ptr<TIterator> parIt = std::unique_ptr<TIterator>(params->createIterator());
      for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
        std::string name = it->GetName();
        if (name.rfind("N_",0)==0 && name.rfind(col)!=std::string::npos) {
          const std::string& par = Form("%s[ %.2f, %.2f, %.2f ]", name.c_str(), it->getVal(), it->getMin(), it->getMax());
          if (myws.var(name.c_str())) {
            myws.var(name.c_str())->setVal(it->getVal());
            myws.var(name.c_str())->setMin(it->getMin());
            myws.var(name.c_str())->setMax(it->getMax());
          }
          else { myws.factory(par.c_str()); }
          std::cout << "[INFO] Yield loaded : " << name << std::endl;
        }
      }
    }
    else { std::cout << "[INFO] RooDatasets used to extract the Yields are not compatible!" << std::endl; compDS = false; }
  }
  else { std::cout << "[INFO] RooDatasets used to extract the Yields were not found!" << std::endl; compDS = false; }
  file->Close();
  return compDS;
};


bool loadSPlotDS(RooWorkspace& myws, const string& fileName, const string& dsName)
{
  if (gSystem->AccessPathName(fileName.c_str())) {
    std::cout << "[ERROR] File " << fileName << " was not found!" << std::endl;
    return false; // File was not found
  }
  auto file = std::unique_ptr<TFile>(new TFile(fileName.c_str()));
  if (!file || !file->IsOpen() || file->IsZombie()) return false;
  auto ws = (RooWorkspace*)file->Get("workspace");
  if (!ws) {
    std::cout << "[ERROR] Workspace was not found in: " << fileName << std::endl;
    file->Close();
    return false;
  }
  if (ws->data(dsName.c_str())) {
    myws.import(*ws->data(dsName.c_str()), RooFit::Rename((dsName+"_INPUT").c_str()));
    if (myws.data((dsName+"_INPUT").c_str())) { std::cout << "[INFO] RooDataset " << (dsName+"_INPUT") << " was imported!" << std::endl; }
    else { std::cout << "[ERROR] Importing RooDataset " << (dsName+"_INPUT") << " failed!" << std::endl; }
  }
  else {
    std::cout << "[ERROR] RooDataset " << dsName << " was not found!" << std::endl;
    file->Close();
    return false;
  }
  file->Close();
  return true;
};


bool createBinnedDataset(RooWorkspace& ws, const std::string& var="Cand_Mass")
{
  //
  const std::string& DSTAG = (ws.obj("DSTAG"))     ? ((RooStringVar*)ws.obj("DSTAG")    )->getVal() : "";
  const std::string& chg   = (ws.obj("fitCharge")) ? ((RooStringVar*)ws.obj("fitCharge"))->getVal() : "";
  const std::string& dsName = ( "d" + chg + "_" + DSTAG );
  //
  if (ws.data(dsName.c_str())==NULL) { std::cout << "[WARNING] DataSet " << dsName << " was not found!" << std::endl; return false; }
  if (ws.data(dsName.c_str())->numEntries()<=10.0) { std::cout << "[WARNING] DataSet " << dsName << " has too few events!" << std::endl; return false; }
  if (ws.var(var.c_str())==NULL) { std::cout << "[WARNING] Variable " << var << " was not found!" << std::endl; return false; }
  //
  const auto& min  = ws.var(var.c_str())->getMin();
  const auto& max  = ws.var(var.c_str())->getMax();
  const uint& nBin = ws.var(var.c_str())->getBins();
  //
  // Create the histogram
  string histName = dsName + "_" + var;
  histName.replace(histName.find("d"), std::string("d").length(), "h");
  std::unique_ptr<TH1D> hist = std::unique_ptr<TH1D>((TH1D*)ws.data(dsName.c_str())->createHistogram(histName.c_str(), *ws.var(var.c_str()), RooFit::Binning(nBin, min, max)));
  if (hist==NULL) { std::cout << "[WARNING] Histogram " << histName << " is NULL!" << std::endl; return false; }
  // Cleaning the input histogram
  // 1) Remove the Under and Overflow bins
  hist->ClearUnderflowAndOverflow();
  // 2) Set negative bin content to zero
  for (int i=0; i<=hist->GetNbinsX(); i++) { if (hist->GetBinContent(i)<0.0) { hist->SetBinContent(i, 0.0); } }
  // 2) Reduce the range of histogram and rebin it
  //hist.reset((TH1D*)rebinhist(*hist, range[1], range[2]));
  if (hist==NULL) { std::cout << "[WARNING] Cleaned Histogram of " << histName << " is NULL!" << std::endl; return false; }
  string dataName = (dsName +"_"+var+"_FIT");
  std::unique_ptr<RooDataHist> dataHist = std::unique_ptr<RooDataHist>(new RooDataHist(dataName.c_str(), "", *ws.var(var.c_str()), hist.get()));
  if (dataHist==NULL) { std::cout << "[WARNING] DataHist used to create " << dsName << " failed!" << std::endl; return false; }
  if (dataHist->sumEntries()==0) { std::cout << "[WARNING] DataHist used to create " << dsName << " is empty!" << std::endl; return false; }
  if (std::abs(dataHist->sumEntries() - hist->GetSumOfWeights())>0.001) { std::cout << "[ERROR] DataHist used to create " << dsName << "  " << " is invalid!  " << std::endl; return false; }
  ws.import(*dataHist);
  ws.var(var.c_str())->setBins(nBin); // Bug Fix
  return true;
};

//
//---------------------------------------------------------------------------------------------
//


bool setModel(StringDiMap_t& model, GlobalInfo&  info, const std::string& type="Cand_Mass")
{
  const auto& cha = info.Par.at("channel");
  for (const auto& col : info.StrV.at("fitSystem")) {
    for (const auto& obj : info.StrV.at("fitObject")) {
      for (const auto& chg : info.StrV.at("fitCharge")) {
        const std::string& label = Form("Model_%s_%s", (obj+cha+chg).c_str(), col.c_str());
        std::string inputLabel = label;
        info.StrV["tags"].push_back(obj+cha+chg+"_"+col);
        const std::vector<std::string> tryChannel = { cha , "" };
        std::vector<std::string> trySystem = { col };
        const std::vector<std::string> tryCharge = { chg , "" };
        for (const auto& tryCha : tryChannel) {
          bool trySuccess = false;
          for (const auto& tryCol : trySystem) {
            for (const auto& tryChg : tryCharge) {
              if (info.Par.count(inputLabel)==0) { inputLabel = ("Model_" + obj + tryCha + tryChg + "_" + tryCol); } else { trySuccess = true; break; }
            }
            if (trySuccess) break;
          }
          if (trySuccess) break;
        }
        if (info.Par.count(inputLabel)>0) {
          std::string value = info.Par.at(inputLabel);
          info.Par[label] = value;
          std::vector<std::string> k;
          if (value.find("+")!=std::string::npos) { if (!splitString(k, value, "+")) { return false; } }
          else { k.push_back(value); }
          for (auto& kk : k) {
            std::string modelName = kk;
            std::vector<std::string> p;
            if (kk.find("[")!=std::string::npos) { 
              modelName = kk.substr(0, kk.find("["));
              kk.erase(0, kk.find("[")+std::string("[").length());
              if (kk.rfind("]")==std::string::npos) { std::cout << "[ERROR] Missing ']' in model: " << value << std::endl; return false; }
              kk.erase(kk.rfind("]"),kk.length());
              if (kk.find(";")!=std::string::npos) { if (!splitString(p, kk, ";")) { return false; } }
              else if (kk.find(",")!=std::string::npos) { if (!splitString(p, kk, ",")) { return false; } }
              else { p.push_back(kk); }
            }
            else { p.push_back(obj); }
            for (const auto& ll : p) {
              std::string objectName  = Form("%s", (ll+cha+chg).c_str());;
              if (info.Flag.at("fitMC") && ll!=obj && modelName=="TEMP") continue;
              if (modelName=="TEMP") { modelName = "Template";    }
              if (modelName=="MJET") { modelName = "MultiJetBkg"; }
              if (ModelDictionary.at(modelName)==0) {
                std::cout << "[ERROR] The " << (ll+cha+chg) << " " << type << " model: " << modelName << " is invalid" << std::endl; return false;
              }
              model[label][ll] = modelName;
              if (modelName == "Template") {
                std::string dsTag = ( "MC_" + ll + "_" + info.Par.at("channelDS") + "_" + col );
                info.StrV[Form("TEMPDS_%s_%s", (obj+cha+chg).c_str(), col.c_str())].push_back(dsTag);
              }
              info.StrV["addObjectModel_"+(obj+cha+chg)+"_"+col].push_back(ll);
            }
          }
        } else {
          std::cout << "[ERROR] " << (obj+cha+chg) << " " << type << " model for " << col << " was not found in the initial parameters!" << std::endl; return false;
        }
      }
    }
  }
  return true;
};


void setGlobalParameterRange(RooWorkspace& myws, const GlobalInfo& info, const std::string& var="Cand_Mass")
{
  if (!myws.var(var.c_str())) { std::cout << "[ERROR] Parameter " << var << " does not exist, failed to set global parameter range!" << std::endl; return; }
  myws.var(var.c_str())->setRange((var+"Window").c_str(), info.Var.at(var).at("Min"), info.Var.at(var).at("Max"));
  const auto& nBins = std::min(int( std::round((info.Var.at(var).at("Max") - info.Var.at(var).at("Min"))/info.Var.at(var).at("binWidth")) ), 2000);
  myws.var(var.c_str())->setBins(nBins, (var+"Window").c_str());
  myws.var(var.c_str())->setBins(nBins);
  return;
};


int importDataset(RooWorkspace& myws, const std::map<string, RooWorkspace>& inputWS, const GlobalInfo& info, const std::string& chg)
{
  // Check info container
  if (info.StrV.count("dsList")==0) { std::cout << "[ERROR] DSList was not found while importing dataset!" << std::endl; return -1; }
  //
  // Define the selection string
  std::string cutDS = "";
  if (info.Par.count("Event_Type")>0 && info.Par.at("Event_Type")!="") { cutDS += Form("(Event_Type==Event_Type::%s)&&", info.Par.at("Event_Type").c_str()); }
  for (const auto& var : info.Var) {
    bool addCut = false;
    // For Mass Resonance analysis
    if (var.first=="Cand_Mass" || var.first=="Cand_Pt" || var.first=="Cand_Rap" || var.first=="Cand_AbsRap" || var.first=="Dau1_InvBeta" || var.first=="Dau2_InvBeta") { addCut = true; }
    // For PbPb analysis
    if (info.Flag.at("fitPbPb") && var.first=="Centrality") { addCut = true; }
    if (addCut) {
      if (var.second.at("Min")==var.second.at("Max")) { cutDS += Form("(%s == %g)", var.first.c_str(), var.second.at("Max")); }
      else { cutDS += Form("(%g <= %s && %s < %g)", var.second.at("Min"), var.first.c_str(), var.first.c_str(), var.second.at("Max")); }
      cutDS += "&&";
    }
  }
  cutDS.erase(cutDS.size()-string("&&").length(), cutDS.size());
  TObjString tmp; tmp.SetString(cutDS.c_str()); myws.import(*((TObject*)&tmp), "Cut_DataSet"); // Save the cut expression for bookkeeping
  std::cout << "[INFO] Importing local RooDataSets with cuts: " << cutDS << std::endl;
  // Reduce and import the datasets
  for (const auto& labelT : info.StrV.at("dsList")) {
    const auto& label = (labelT!="MC_D0Swap_PIONKAON_PbPb" ? labelT : "MC_D0_PIONKAON_PbPb");
    const bool& isSwapDS = (labelT.find("Swap_")!=std::string::npos);
    // Extract the RooDatasets
    std::string dsType = "";
    for (const auto& type : std::vector<std::string>({ "COR", "LUM", "SET", "RAW" })) {
      if (inputWS.at(label).data(Form("d%s_%s_%s", chg.c_str(), type.c_str(), label.c_str()))) { dsType = type; break; }
    }
    if (dsType=="") { std::cout << "[ERROR] Sample RAW_" << label << " was not found!" << std::endl; return -1; }
    const std::string& extLabel = dsType + "_" + label;
    std::cout << "[INFO] Importing local RooDataSet " << extLabel << std::endl;
    //
    if (chg!="") {
      const std::string& dsName = Form("d%s_%s", chg.c_str(), labelT.c_str());
      if (!myws.data(dsName.c_str())) {
        if ( inputWS.count(label)==0 || !(inputWS.at(label).data(Form("d%s_%s", chg.c_str(), extLabel.c_str())))){ 
          std::cout << "[ERROR] The dataset " <<  Form("d%s_%s", chg.c_str(), extLabel.c_str()) << " was not found!" << std::endl;
          return -1;
        }
        const auto& cutDST = cutDS + (isSwapDS ? "&&(Cand_IsSwap==1)" : "");
        auto data = std::unique_ptr<RooDataSet>((RooDataSet*)inputWS.at(label).data(Form("d%s_%s", chg.c_str(), extLabel.c_str()))->reduce(cutDST.c_str()));
        if (data==NULL || data->sumEntries()==0){
          if (extLabel.rfind("MC_",0)==0 || chg=="SS") {
            std::cout << "[WARNING] No events from dataset " <<  Form("d%s_%s", chg.c_str(), extLabel.c_str()) << " passed the kinematic cuts!" << std::endl;
          }
          else { std::cout << "[ERROR] No events from dataset " <<  Form("d%s_%s", chg.c_str(), extLabel.c_str()) << " passed the kinematic cuts!" << std::endl; return -1; }
        }
        else {
          data->SetName(dsName.c_str());
          myws.import(*data);
        }
        std::cout << "[INFO] " << Form("%.0f", data->sumEntries()) << " weighted entries imported from local RooDataSet " << dsName << std::endl;
        //
        // Set the range of each global parameter in the local roodataset
        if (myws.data(dsName.c_str())!=NULL) {
          const auto& row = myws.data(dsName.c_str())->get();
          for (const auto& var : info.Var) {
            if (row->find(var.first.c_str()) && (var.second.count("Min")>0)) {
              ((RooRealVar*)row->find(var.first.c_str()))->setMin(var.second.at("Min"));
              ((RooRealVar*)row->find(var.first.c_str()))->setMax(var.second.at("Max"));
            }
          }
        }
      }
      //
      const auto& dsSPLOTInputName = (dsName+"_SPLOT_INPUT");
      if (myws.data(dsSPLOTInputName.c_str())) {
        auto data = std::unique_ptr<RooDataSet>((RooDataSet*)myws.data(dsSPLOTInputName.c_str())->reduce(cutDS.c_str()));
        // Set the range of each global parameter in the local roodataset
        if (data!=NULL) {
          const auto& row = data->get();
          for (const auto& var : info.Var) {
            if (row->find(var.first.c_str()) && (var.second.count("Min")>0)) {
              ((RooRealVar*)row->find(var.first.c_str()))->setMin(var.second.at("Min"));
              ((RooRealVar*)row->find(var.first.c_str()))->setMax(var.second.at("Max"));
            }
          }
        }
        if (data->sumEntries()==0){ std::cout << "[ERROR] No events from dataset " <<  dsSPLOTInputName << " passed the kinematic cuts!" << std::endl; }
        else if (!isCompatibleDataset(*data, *(RooDataSet*)myws.data(dsName.c_str()))) { cout << "[ERROR] sPlot and Original Datasets are inconsistent!" << std::endl; return -1; }
        else {
          const auto& dsSPLOTName = (dsName+"_SPLOT");
          data->SetName(dsSPLOTName.c_str());
          myws.import(*data, RooFit::Rename(dsSPLOTName.c_str()));
          if (myws.data(dsSPLOTName.c_str())) { std::cout << "[INFO] RooDataset " << dsSPLOTName << " was imported!" << std::endl; }
          else { std::cout << "[ERROR] Importing RooDataset " << dsSPLOTName << " failed!" << std::endl; return -1; }
          std::cout << "[INFO] SPlotDS Events: " << data->sumEntries() << " , origDS Events: " << myws.data(Form("dOS_%s", label.c_str()))->sumEntries() << std::endl;
        }
      }
    }
  }
  // Set the range of each global parameter in the local workspace
  for (const auto& var : info.Var) {
    if ( myws.var(var.first.c_str()) && (var.second.count("Min")>0) ) {
      myws.var(var.first.c_str())->setMin(var.second.at("Min"));
      myws.var(var.first.c_str())->setMax(var.second.at("Max"));
    }
    else if ( (myws.var(var.first.c_str())==NULL) && (var.second.count("Val")>0) ) {
      myws.factory(Form("%s[%.10f]", var.first.c_str(), var.second.at("Val")));
    }
  }
  // Print bin information
  std::string binInfo = "[INFO] Analyzing bin:";
  for (const auto& var : std::vector<std::string>({ "Cand_Pt", "Cand_AbsRap", "Centrality" })) {
    if (myws.var(var.c_str())) {
      std::string varName = var;
      double varMin = myws.var(var.c_str())->getMin();
      double varMax = myws.var(var.c_str())->getMax();
      binInfo += Form(" %g < %s < %g ,", varMin, varName.c_str(), varMax);
    }
  }
  binInfo.erase(binInfo.size()-string(" ,").length(), binInfo.size());
  std::cout << binInfo << std::endl;
  return 1;
};
 

#endif // #ifndef rooDataUtils_h
