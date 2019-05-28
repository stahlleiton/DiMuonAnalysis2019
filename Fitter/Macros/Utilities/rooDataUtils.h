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
#include "RooFitResult.h"
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
  for (const auto& elem : map) { stringReplace( str, elem.first, elem.second ); }
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
  const auto& listVar = ws.allVars();
  auto parIt = std::unique_ptr<TIterator>(listVar.createIterator());
  for (auto itp = parIt->Next(); itp!=NULL; itp = parIt->Next() ) {
    const auto& it = dynamic_cast<RooRealVar*>(itp); if (!it) continue;
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


void clearWorkspace(RooWorkspace& inWS, const std::string& smp="All", const bool& delVar=true)
{
  //
  const auto& level = RooMsgService::instance().globalKillBelow();
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  // Copy all category variables
  if (inWS.allCats().getSize()>0 && delVar) {
    const auto& listCat = inWS.allCats();
    auto catIt = std::unique_ptr<TIterator>(listCat.createIterator());
    for (auto it = catIt->Next(); it!=NULL; it = catIt->Next()) { if (inWS.cat(it->GetName())) { inWS.RecursiveRemove(it); if(it) delete it; } }
  }
  // Copy all category functions
  if (inWS.allCatFunctions().getSize()>0 && delVar) {
    const auto& listFnc = inWS.allCatFunctions();
    auto fncIt = std::unique_ptr<TIterator>(listFnc.createIterator());
    for (auto it = fncIt->Next(); it!=NULL; it = fncIt->Next()) { if (inWS.catfunc(it->GetName())) { inWS.RecursiveRemove(it); if(it) delete it; } }
  }
  // Copy all variables
  if (inWS.allVars().getSize()>0 && delVar) {
    const auto& listVar = inWS.allVars();
    auto parIt = std::unique_ptr<TIterator>(listVar.createIterator());
    for (auto it = parIt->Next(); it!=NULL; it = parIt->Next()) { if (inWS.var(it->GetName())) { inWS.RecursiveRemove(it); if(it) delete it; } }
  }
  // Copy all functions
  if (inWS.allFunctions().getSize()>0 && delVar) {
    const auto& listFnc = inWS.allFunctions();
   auto fncIt = std::unique_ptr<TIterator>(listFnc.createIterator());
    for (auto it = fncIt->Next(); it!=NULL; it = fncIt->Next()) { if (inWS.function(it->GetName())) { inWS.RecursiveRemove(it); if(it) delete it; } }
  }
  // Copy all PDFs
  if (inWS.allPdfs().getSize()>0) {
    const auto& listPdf = inWS.allPdfs();
    auto pdfIt = std::unique_ptr<TIterator>(listPdf.createIterator());
    for (auto it = pdfIt->Next(); it!=NULL; it = pdfIt->Next()) { if (inWS.pdf(it->GetName())) { inWS.RecursiveRemove(it); if(it) delete it; } }
  }
  // Copy all Datasets
  if (inWS.allData().size()>0 && smp!="") {
    const auto& listData = inWS.allData();
    for (const auto& it : listData) { if (inWS.data(it->GetName()) && (smp=="All" || std::string(it->GetName()).find(smp)!=std::string::npos)) { inWS.RecursiveRemove(it); if(it) delete it; } }
  }
  // Copy all Embedded Datasets
  if (inWS.allEmbeddedData().size()>0 && smp!="") {
    const auto& listData = inWS.allEmbeddedData();
    for (const auto& it : listData) { if (inWS.embeddedData(it->GetName()) && (smp=="All" || std::string(it->GetName()).find(smp)!=std::string::npos)) { inWS.RecursiveRemove(it); if(it) delete it; } }
  }
  // Copy all Generic Objects
  if (inWS.allGenericObjects().size()>0) {
    const auto& listObj = inWS.allGenericObjects();
    for (const auto& it : listObj) { if (inWS.genobj(it->GetTitle())) { inWS.RecursiveRemove(it); if(it) delete it; } }
  }
  // Return the RooMessenger Level
  RooMsgService::instance().setGlobalKillBelow(level);
};


void copyWorkspace(RooWorkspace& outWS, const RooWorkspace& inWS, const std::string& smp="All", const bool& addVar=true, const bool& addSnap=false)
{
  //
  const auto& level = RooMsgService::instance().globalKillBelow();
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  // Copy all category variables
  if (inWS.allCats().getSize()>0 && addVar) {
    const auto& listCat = inWS.allCats();
    auto catIt = std::unique_ptr<TIterator>(listCat.createIterator());
    for (auto itp = catIt->Next(); itp!=NULL; itp = catIt->Next() ) { const auto& it = dynamic_cast<RooCategory*>(itp); if (it && !outWS.cat(it->GetName())) { outWS.import(*it); } }
  }
  // Copy all category functions
  if (inWS.allCatFunctions().getSize()>0 && addVar) {
    const auto& listFnc = inWS.allCatFunctions();
    auto fncIt = std::unique_ptr<TIterator>(listFnc.createIterator());
    for (auto itp = fncIt->Next(); itp!=NULL; itp = fncIt->Next() ) { const auto& it = dynamic_cast<RooCategory*>(itp); if (it && !outWS.catfunc(it->GetName())) { outWS.import(*it); } }
  }
  // Copy all variables
  if (inWS.allVars().getSize()>0 && addVar) {
    const auto& listVar = inWS.allVars();
    auto parIt = std::unique_ptr<TIterator>(listVar.createIterator());
    for (auto itp = parIt->Next(); itp!=NULL; itp = parIt->Next() ) { const auto& it = dynamic_cast<RooRealVar*>(itp); if (it && !outWS.var(it->GetName())) { outWS.import(*it); } }
  }
  // Copy all functions
  if (inWS.allFunctions().getSize()>0 && addVar) {
    const auto& listFnc = inWS.allFunctions();
    auto fncIt = std::unique_ptr<TIterator>(listFnc.createIterator());
    for (auto itp = fncIt->Next(); itp!=NULL; itp = fncIt->Next() ) { const auto& it = dynamic_cast<RooRealVar*>(itp); if (it && !outWS.function(it->GetName())) { outWS.import(*it); } }
  }
  // Copy all PDFs
  if (inWS.allPdfs().getSize()>0) {
    const auto& listPdf = inWS.allPdfs();
    auto pdfIt = std::unique_ptr<TIterator>(listPdf.createIterator());
    for (auto itp = pdfIt->Next(); itp!=NULL; itp = pdfIt->Next() ) { const auto& it = dynamic_cast<RooAbsPdf*>(itp); if (it && !outWS.pdf(it->GetName())) { outWS.import(*it, RooFit::RecycleConflictNodes()); } }
  }
  // Copy all Datasets
  if (inWS.allData().size()>0 && smp!="") {
    const auto& listData = inWS.allData();
    for (const auto& it : listData) { if (!outWS.data(it->GetName()) && (smp=="All" || std::string(it->GetName()).find(smp)!=std::string::npos)) { outWS.import(*it); } }
  }
  // Copy all Embedded Datasets
  if (inWS.allEmbeddedData().size()>0 && smp!="") {
    const auto& listData = inWS.allEmbeddedData();
    for (const auto& it : listData) { if (!outWS.embeddedData(it->GetName()) && (smp=="All" || std::string(it->GetName()).find(smp)!=std::string::npos)) { outWS.import(*it, RooFit::Embedded(1)); } }
  }
  // Copy all Generic Objects
  if (inWS.allGenericObjects().size()>0) {
    const auto& listObj = inWS.allGenericObjects();
    for (const auto& it : listObj) { if (!outWS.genobj(it->GetTitle())) { outWS.import(*it, it->GetTitle()); } }
  }
  // Copy all snapshots variables
  if (addSnap) {
    const auto& snapNames = StringVector_t({ "initialParameters", "fittedParameters" });
    for (const auto& snapName : snapNames) {
      const auto& snap = inWS.getSnapshot(snapName.c_str());
      if (snap) { outWS.saveSnapshot(snapName.c_str(), *snap, kTRUE); }
    }
  }
  // Return the RooMessenger Level
  RooMsgService::instance().setGlobalKillBelow(level);
};


bool saveWorkSpace(const RooWorkspace& ws, RooFitResult* fitResult, const string& outputDir, const string& fileName, const bool& saveAll=true)
{
  // Save the workspace
  gSystem->mkdir(outputDir.c_str(), kTRUE);
  auto file = std::unique_ptr<TFile>(new TFile((outputDir+fileName).c_str(), "RECREATE"));
  if (!file || !file->IsOpen() || file->IsZombie()) {
    std::cout << "[ERROR] Output root file with fit results could not be created!" << std::endl; if (file) { file->Close(); }; return false;
  }
  else {
    file->cd();
    if (saveAll) { ws.Write("workspace"); }
    else {
      RooWorkspace tmpWS;
      copyWorkspace(tmpWS, ws, "", true, true);
      tmpWS.Write("workspace");
    }
    if (fitResult) { fitResult->Write("fitResult"); }
    file->Write(); file->Close();
  }
  std::cout << "[INFO] RooWorkspace saved in: " << (outputDir+fileName) << std::endl;
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
  auto parIt = std::unique_ptr<TIterator>(listVar->createIterator());
  for (auto itp = parIt->Next(); itp!=NULL; itp = parIt->Next()) {
    const auto& it = dynamic_cast<RooRealVar*>(itp); if (!it) continue;
    if ( ds.get()->find(it->GetName()) == NULL ) { std::cout << "[ERROR] DataSet " << ds.GetName() << " does not contain the variable " << it->GetName() << " !" << std::endl; return false; }
    if (checkRange) {
      const auto& var = static_cast<RooRealVar*>(ds.get()->find(it->GetName()));
      if ( it->getMin() != var->getMin() ) { std::cout << "[ERROR] " << it->GetName() << " Min Range disagreement : "
						       << ds.GetName() << " ( " << var->getMin() << " ) " << " and "
						       << ref.GetName() << " ( " << it->getMin() << " ) ! " << std::endl; return false; }
      if ( it->getMax() != var->getMax() ) { std::cout << "[ERROR] " << it->GetName() << " Max Range disagreement : "
						       << ds.GetName() << " ( " << var->getMax() << " ) " << " and "
						       << ref.GetName() << " ( " << it->getMax() << " ) ! " << std::endl; return false; }
    }
    for (uint i = 1; i < 5; i++) {
      const auto& rMom = ref.moment(*it, i);
      const auto& dMom = ds.moment(*it, i);
      if (rMom!=dMom) { std::cout << "[ERROR] " << it->GetName() << " " << i << " Moment Value disagreement : "
				  << ds.GetName() << " ( " << dMom << " ) " << " and "
				  << ref.GetName() << " ( " << rMom << " ) " << " ! " << std::endl; return false; }
    }
  }
  // DataSets are compatible if they passed all tests
  std::cout << "[INFO] DataSets " << ds.GetName() << " and " << ref.GetName() << " are compatible!" << std::endl; return true;
};


bool compareSnapshots(const RooArgSet& pars1, const RooArgSet& pars2)
{
  auto parIt = std::unique_ptr<TIterator>(pars1.createIterator());
  for (auto itp = parIt->Next(); itp!=NULL; itp = parIt->Next() ) {
    const auto& it = dynamic_cast<RooRealVar*>(itp); if (!it) continue;
    const auto& val = pars2.getRealValue(it->GetName(),-1e99);
    if (val==-1e99) return false;           // the parameter was not found!
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
  if (!file || !file->IsOpen() || file->IsZombie()) { if(file) { file->Close(); }; return false; }
  const auto& ws = dynamic_cast<RooWorkspace*>(file->Get("workspace"));
  if (!ws) { std::cout << "[INFO] Workspace was not found" << std::endl; file->Close(); return false; }
  const auto& params = ( ws->getSnapshot("initialParameters") ? ws->getSnapshot("initialParameters") :  ws->getSnapshot(Form("%s_parIni", pdfName.c_str())) );
  if (!params) { std::cout << "[INFO] Snapshot of initial parameters was not found!" << std::endl; file->Close(); return false; }
  bool result = compareSnapshots(newpars, *params);
  file->Close();
  return result;
};


RooRealVar getVar(const RooArgSet& set, const std::string& varName)
{
  if (set.find(varName.c_str()) != NULL) {
    const auto& it = dynamic_cast<RooRealVar*>(set.find(varName.c_str()));
    if (it) { return *it; }
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
    myws.var(var.c_str())->setRange("FitWindow", varMin, varMax);
    myws.var(var.c_str())->setBins(nBins, "FitWindow");
    //
    auto dataToFit = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(myws.data(dsName.c_str())->reduce(varFitRange.c_str())->Clone((dsName+"_FIT").c_str())));
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
  const auto& ws = dynamic_cast<RooWorkspace*>(file->Get("workspace"));
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
  const auto& ws = dynamic_cast<RooWorkspace*>(file->Get("workspace"));
  if (!ws) {
    std::cout << "[INFO] Workspace was not found in: " << fileName << std::endl;
    file->Close();
    return false;
  }
  bool compDS = true;
  if (ws->data(dsName.c_str()) && myws.data(dsName.c_str())) {
    if (isCompatibleDataset(*dynamic_cast<RooDataSet*>(myws.data(dsName.c_str())), *dynamic_cast<RooDataSet*>(ws->data(dsName.c_str())))) {
      const auto& params = ws->getSnapshot(Form("%s_parIni", pdfName.c_str()));
      if (!params) {
        std::cout << "[INFO] Snapshot " << pdfName << "_parIni was not found!" << std::endl;
        file->Close();
        return false;
      }
      auto parIt = std::unique_ptr<TIterator>(params->createIterator());
      for (auto itp = parIt->Next(); itp!=NULL; itp = parIt->Next() ) {
	const auto& it = dynamic_cast<RooRealVar*>(itp); if (!it) continue;
        const std::string& name = it->GetName();
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
  const auto& ws = dynamic_cast<RooWorkspace*>(file->Get("workspace"));
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
  const std::string& DSTAG = (ws.obj("DSTAG"))     ? dynamic_cast<RooStringVar*>(ws.obj("DSTAG")    )->getVal() : "";
  const std::string& chg   = (ws.obj("fitCharge")) ? dynamic_cast<RooStringVar*>(ws.obj("fitCharge"))->getVal() : "";
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
  auto histName = dsName + "_" + var;
  histName.replace(histName.find("d"), std::string("d").length(), "h");
  std::unique_ptr<TH1D> hist = std::unique_ptr<TH1D>(dynamic_cast<TH1D*>(ws.data(dsName.c_str())->createHistogram(histName.c_str(), *ws.var(var.c_str()), RooFit::Binning(nBin, min, max))));
  if (hist==NULL) { std::cout << "[WARNING] Histogram " << histName << " is NULL!" << std::endl; return false; }
  // Cleaning the input histogram
  // 1) Remove the Under and Overflow bins
  hist->ClearUnderflowAndOverflow();
  // 2) Set negative bin content to zero
  for (int i=0; i<=hist->GetNbinsX(); i++) { if (hist->GetBinContent(i)<0.0) { hist->SetBinContent(i, 0.0); } }
  // 2) Reduce the range of histogram and rebin it
  //hist.reset(dynamic_cast<TH1D*>(rebinhist(*hist, range[1], range[2])));
  if (hist==NULL) { std::cout << "[WARNING] Cleaned Histogram of " << histName << " is NULL!" << std::endl; return false; }
  const auto& dataName = (dsName +"_"+var+"_FIT");
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


std::string findLabel(const std::string& par, const std::string& obj, const std::string& chg,
		      const std::string& col, const std::string& cha, const GlobalInfo& info)
{
  std::string tryLabel="";
  const StringVector_t tryChannel = { cha , "" };
  StringVector_t trySystem  = { col , "" };
  if (info.Flag.at("doPA8Y16")) { trySystem.push_back("PA8Y16"); }
  const StringVector_t tryCharge  = { chg , "" };
  for (const auto& tryCha : tryChannel) {
    bool trySuccess = false;
    for (const auto& tryCol : trySystem) {
      for (const auto& tryChg : tryCharge) {
	tryLabel = obj + tryCha + tryChg + (tryCol!="" ? "_"+tryCol : "");
	if (info.Par.count(par+"_"+tryLabel)>0) { trySuccess = true; break; }
      }
      if (trySuccess) break;
    }
    if (trySuccess) break;
  }
  return tryLabel;
};


bool setModel(StringDiMap_t& model, GlobalInfo&  info, const std::string& type="Cand_Mass")
{
  const auto& cha = info.Par.at("channel");
  for (const auto& col : info.StrS.at("fitSystem")) {
    for (const auto& obj : info.StrS.at("fitObject")) {
      for (const auto& chg : info.StrS.at("fitCharge")) {
        const std::string& label = Form("Model_%s_%s", (obj+cha+chg).c_str(), col.c_str());
        info.StrS["tags"].insert(obj+cha+chg+"_"+col);
	const std::string inputLabel = "Model_"+findLabel("Model", obj, chg, col, cha, info);
        if (info.Par.count(inputLabel)>0) {
          const auto& value = info.Par.at(inputLabel);
          info.Par[label] = value;
          StringVector_t k;
          if (value.find("+")!=std::string::npos) { splitString(k, value, "+"); }
          else { k.push_back(value); }
          for (auto& kk : k) {
            std::string modelName = kk;
            StringVector_t p;
            if (kk.find("[")!=std::string::npos) { 
              modelName = kk.substr(0, kk.find("["));
              kk.erase(0, kk.find("[")+std::string("[").length());
              if (kk.rfind("]")==std::string::npos) { std::cout << "[ERROR] Missing ']' in model: " << value << std::endl; return false; }
              kk = kk.substr(0, kk.rfind("]"));
              if (kk.find(";")!=std::string::npos) { splitString(p, kk, ";"); }
              else if (kk.find(",")!=std::string::npos) { splitString(p, kk, ","); }
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
                info.StrS[Form("TEMPDS_%s_%s", (obj+cha+chg).c_str(), col.c_str())].insert(dsTag);
              }
              info.StrS["addObjectModel_"+(obj+cha+chg)+"_"+col].insert(ll);
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


void setDSParamaterRange(const RooDataSet& ds, const GlobalInfo& info)
{
  const auto& row = ds.get();
  for (const auto& var : info.Var) {
    const bool& isAbs = (var.first.find("Abs")!=std::string::npos);
    auto varN = var.first; if (isAbs) { varN.erase(varN.find("Abs"), 3); }
    if (row->find(varN.c_str()) && (var.second.count("Min")>0)) {
      const auto& v = dynamic_cast<RooRealVar*>(row->find(varN.c_str()));
      v->setMin(isAbs ? -var.second.at("Max") : var.second.at("Min"));
      v->setMax(var.second.at("Max"));
    }
  }
};


bool setFitParameterRange(RooWorkspace& myws, const GlobalInfo& info)
{
  for (const auto& var : info.StrS.at("fitVariable")) {
    if (!myws.var(var.c_str())) { std::cout << "[ERROR] Parameter " << var << " does not exist, failed to set fit parameter range!" << std::endl; return false; }
    myws.var(var.c_str())->setRange("FitWindow", info.Var.at(var).at("Min"), info.Var.at(var).at("Max"));
    const auto& nBins = std::min(int( std::round((info.Var.at(var).at("Max") - info.Var.at(var).at("Min"))/info.Var.at(var).at("binWidth")) ), 2000);
    myws.var(var.c_str())->setBins(nBins, "FitWindow");
    myws.var(var.c_str())->setBins(nBins);
  }
  return true;
};


int importDataset(RooWorkspace& myws, GlobalInfo& info, const RooWorkspaceMap_t& inputWS, const std::string& chg)
{
  // Check info container
  if (info.StrS.count("dsList")==0) { std::cout << "[ERROR] DSList was not found while importing dataset!" << std::endl; return -1; }
  for (const auto& inWS : inputWS) { if (inWS.second.allData().size()==0) { std::cout << "[ERROR] Input workspace " << inWS.first << " is empty!" << std::endl; return -1; } }
  //
  // Define the selection string
  std::string cutDS = "";
  info.StrS["cutPars"].clear();
  for (const auto& it : inputWS.at(*info.StrS.at("dsList").begin()).allData()) {
    const std::string& inDS = it->GetName();
    if (inDS.rfind("d"+chg+"_",0)!=0) continue;
    const auto& dsVars = it->get();
    for (const auto& p : info.Par) {
      if (dsVars->find(p.first.c_str()) && p.second!="") { cutDS += Form("(%s == %s::%s) && ", p.first.c_str(), p.first.c_str(), p.second.c_str()); info.StrS.at("cutPars").insert(p.first); }
    }
    for (const auto& v : info.Var) {
      const bool& isAbs = (v.first.find("Abs")!=std::string::npos);
      auto varN = v.first; if (isAbs) { varN.erase(varN.find("Abs"), 3); }
      if (!dsVars->find(varN.c_str())) continue;
      if (v.second.at("Min")==v.second.at("Default_Min") && v.second.at("Max")==v.second.at("Default_Max")) continue;
      if (v.first=="Centrality" && (inDS.rfind("PbPb")==std::string::npos || info.Par.at("PD")=="UPC")) continue;
      if (isAbs) { varN = "abs("+varN+")"; }
      if (v.second.at("Min")==v.second.at("Max")) { cutDS += Form("(%s == %g) && ", varN.c_str(), v.second.at("Max")); }
      else if (v.second.at("Min")==v.second.at("Default_Min")) { cutDS += Form("(%s < %g) && ", varN.c_str(), v.second.at("Max")); }
      else if (v.second.at("Max")==v.second.at("Default_Max")) { cutDS += Form("(%g <= %s) && ", v.second.at("Min"), varN.c_str()); }
      else { cutDS += Form("(%g <= %s && %s < %g) && ", v.second.at("Min"), varN.c_str(), varN.c_str(), v.second.at("Max")); }
      info.StrS.at("cutPars").insert(v.first);
    }
    break;
  }
  cutDS = cutDS.substr(0, cutDS.rfind(" && "));
  TObjString tmp; tmp.SetString(cutDS.c_str()); myws.import(*dynamic_cast<TObject*>(&tmp), "Cut_DataSet"); // Save the cut expression for bookkeeping
  std::cout << "[INFO] Importing local RooDataSets with cuts: " << cutDS << std::endl;
  //
  // Reduce and import the datasets
  for (const auto& labelT : info.StrS.at("dsList")) {
    auto label = labelT;
    const bool& isSwapDS = (labelT.find("Swap_")!=std::string::npos);
    if (isSwapDS) { stringReplace(label, "Swap", ""); }
    // Extract the RooDatasets
    std::string dsType = "";
    for (auto type = info.StrV.at("dsType").rbegin(); type != info.StrV.at("dsType").rend(); ++type) {
      if (inputWS.at(label).data(Form("d%s_%s_%s", chg.c_str(), type->c_str(), label.c_str()))) { dsType = *type; break; }
    }
    if (dsType=="") { std::cout << "[ERROR] Sample " << label << " was not found!" << std::endl; return -1; }
    const auto& extLabel = dsType + "_" + label;
    std::cout << "[INFO] Importing local RooDataSet " << extLabel << std::endl;
    //
    if (chg!="") {
      const auto& dsName = "d"+chg+"_"+labelT;
      const auto& dsExtName = "d"+chg+"_"+extLabel;
      if (!myws.data(dsName.c_str())) {
        if (inputWS.count(label)==0 || !inputWS.at(label).data(dsExtName.c_str())){ 
          std::cout << "[ERROR] The dataset " <<  dsExtName << " was not found!" << std::endl; return -1;
        }
        const auto& cutDST = cutDS + (isSwapDS ? "&&(Cand_IsSwap==Cand_IsSwap::Yes)" : "");
        auto data = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(inputWS.at(label).data(dsExtName.c_str())->reduce(RooFit::Cut(cutDST.c_str()), RooFit::Name(dsName.c_str()), RooFit::Title(dsName.c_str()))));
        if (!data || data->sumEntries()==0){
          if (extLabel.rfind("MC_",0)==0 || chg=="SS") {
            std::cout << "[WARNING] No events from dataset " <<  dsExtName << " passed the kinematic cuts!" << std::endl;
          }
          else { std::cout << "[ERROR] No events from dataset " <<  dsExtName << " passed the kinematic cuts!" << std::endl; return -1; }
        }
        else {
          myws.import(*data);
	  if (!myws.data(dsName.c_str())) { std::cout << "[ERROR] Importing RooDataSet " <<  dsName << " failed!" << std::endl; return -1; }
        }
        std::cout << "[INFO] " << data->numEntries() << " entries imported from local RooDataSet " << dsName << std::endl;
        // Set the range of each global parameter in the local roodataset
	if (myws.data(dsName.c_str())) { setDSParamaterRange(*dynamic_cast<RooDataSet*>(myws.data(dsName.c_str())), info); }
      }
      //
      const auto& dsSPLOTInputName = (dsName+"_SPLOT_INPUT");
      if (myws.data(dsSPLOTInputName.c_str())) {
        auto data = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(myws.data(dsSPLOTInputName.c_str())->reduce(cutDS.c_str())));
        // Set the range of each global parameter in the local roodataset
	setDSParamaterRange(*data, info);
        if (data->sumEntries()==0){ std::cout << "[ERROR] No events from dataset " <<  dsSPLOTInputName << " passed the kinematic cuts!" << std::endl; }
        else if (!isCompatibleDataset(*data, *dynamic_cast<RooDataSet*>(myws.data(dsName.c_str())))) { cout << "[ERROR] sPlot and Original Datasets are inconsistent!" << std::endl; return -1; }
        else {
          const auto& dsSPLOTName = (dsName+"_SPLOT");
          data->SetName(dsSPLOTName.c_str());
          myws.import(*data, RooFit::Rename(dsSPLOTName.c_str()));
          if (myws.data(dsSPLOTName.c_str())) { std::cout << "[INFO] RooDataSet " << dsSPLOTName << " was imported!" << std::endl; }
          else { std::cout << "[ERROR] Importing RooDataSet " << dsSPLOTName << " failed!" << std::endl; return -1; }
          std::cout << "[INFO] SPlotDS Events: " << data->sumEntries() << " , origDS Events: " << myws.data(Form("dOS_%s", label.c_str()))->sumEntries() << std::endl;
        }
      }
    }
  }
  // Check if the user wants to use the Center-of-Mass frame
  for (auto& v : info.Var) {
    if (info.Flag.count("use"+v.first+"CM")>0 && info.Flag.at("use"+v.first+"CM")) {
      myws.factory(Form("use%sCM[1.0]", v.first.c_str()));
    }
  }
  // Set the range of each global parameter in the local workspace
  for (const auto& var : info.Var) {
    const bool& isAbs = (var.first.find("Abs")!=std::string::npos);
    auto varN = var.first; if (isAbs) { varN.erase(varN.find("Abs"), 3); }
    if (myws.var(varN.c_str()) && var.second.count("Min")>0) {
      myws.var(varN.c_str())->setMin(isAbs ? -var.second.at("Max") : var.second.at("Min"));
      myws.var(varN.c_str())->setMax(var.second.at("Max"));
    }
    else if (!myws.var(var.first.c_str()) && var.second.count("Val")>0) {
      myws.factory(Form("%s[%.10f]", var.first.c_str(), var.second.at("Val")));
    }
  }
  // Print bin information
  std::string binInfo = "[INFO] Analyzing bin:";
  for (const auto& var : info.StrS.at("cutPars")) {
    const bool& isAbs = (var.find("Abs")!=std::string::npos);
    auto varN = var; if (isAbs) { varN.erase(varN.find("Abs"), 3); }
    if (myws.var(varN.c_str())) {
      const auto& varMin = (isAbs ? info.Var.at(var).at("Min") : myws.var(varN.c_str())->getMin());
      const auto& varMax = myws.var(varN.c_str())->getMax();
      if (info.Flag.count("use"+var+"CM")>0  && info.Flag.at("use"+var+"CM")) {
	const bool& ispPb = (info.Flag.at("fitpPb8Y16") || info.Flag.at("fitPA8Y16"));
	binInfo += Form(" %g <= %s < %g ,", pPb::EtaLABtoCM(varMin, ispPb), (var+"CM").c_str(), pPb::EtaLABtoCM(varMax, ispPb));
      }
      else {
	binInfo += Form(" %g <= %s < %g ,", varMin, var.c_str(), varMax);
      }
    }
    else if (myws.cat(var.c_str())) {
      const auto& catLbl = myws.cat(var.c_str())->getLabel();
      binInfo += Form(" %s == %s ,", var.c_str(), catLbl);
    }
  }
  binInfo = binInfo.substr(0, binInfo.rfind(" ,"));
  std::cout << binInfo << std::endl;
  return 1;
};


void setFileName(std::string& fileName, std::string& outputDir, const StringSet_t& fitV, const std::string& DSTAG, const std::string& plotLabel, const GlobalInfo& info)
{
  const auto& dsTag  = DSTAG.substr(0, DSTAG.rfind("_DIMUON"));
  const auto& colTag = DSTAG.substr(DSTAG.find_last_of("_")+1);
  const auto& objTag = *info.StrS.at("fitObject").begin();
  std::string fitVar = "";
  for (const auto& v : fitV) { auto s = v; stringReplace(s, "_", ""); fitVar += s+"_"; }
  fitVar = fitVar.substr(0, fitVar.rfind("_"));
  outputDir = Form("%s%s/%s/%s/%s/", outputDir.c_str(), fitVar.c_str(), dsTag.c_str(), objTag.c_str(), colTag.c_str());
  std::string varLbl = "";
  for (const auto& var : info.StrS.at("cutPars")) {
    bool incVar = true;
    for (const auto& v : fitV) { if (var==v) { incVar = false; break; } }
    if (!incVar || info.Var.count(var)==0 || info.Var.at(var).count("Min")==0) continue;
    if (info.Var.at(var).at("Min")==info.Var.at(var).at("Default_Min") && info.Var.at(var).at("Max")==info.Var.at(var).at("Default_Max")) continue;
    if (var=="Centrality" && (DSTAG.rfind("PbPb")==std::string::npos || info.Par.at("PD")=="UPC")) continue;
    auto varN = var; stringReplace(varN, "_", "");
    auto varMin = info.Var.at(var).at("Min");
    auto varMax = info.Var.at(var).at("Max");
    if (info.Flag.count("use"+var+"CM")>0 && info.Flag.at("use"+var+"CM")) {
      varN += "CM";
      const bool& ispPb = (info.Flag.at("fitpPb8Y16") || info.Flag.at("fitPA8Y16"));
      varMin = pPb::EtaLABtoCM(varMin, ispPb);
      varMax = pPb::EtaLABtoCM(varMax, ispPb);
    }
    if (varMin==varMax) { varLbl += Form("%s_%.0f_", varN.c_str(), varMax*100.); }
    else { varLbl += Form("%s_%.0f_%.0f_", varN.c_str(), varMin*100., varMax*100.); }
  }
  varLbl = varLbl.substr(0, varLbl.rfind("_"));
  fileName = Form("%s_%s_%s%s", fitVar.c_str(), dsTag.c_str(), plotLabel.c_str(), varLbl.c_str());
};
 

#endif // #ifndef rooDataUtils_h
