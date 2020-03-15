#ifndef rooDataUtils_h
#define rooDataUtils_h

#include "TSystem.h"
#include "TFile.h"
#include "TObject.h"
#include "TIterator.h"
#include "TH1.h"
#include "TKey.h"

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
#include "RooList.h"

#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <map>

#include "initClasses.h"
#include "../../../Utilities/dataUtils.h"


void setFileName(std::string& fileName, std::string& outputDir, const std::string& lbl, const GlobalInfo& info, const StringSet_t& varV={})
{
  const auto& cha = info.Par.at("channel");
  const auto& label = lbl.substr(lbl.rfind(cha));
  // Define plot label
  auto varS = info.StrS.at("fitVariable"); if (!varV.empty()) { varS = varV; }
  auto varT = info.StrS.at("fitVarName"); if (!varV.empty()) { varT.clear(); for (const auto& v : varV) { auto t = v; stringReplace(t, "_", ""); varT.insert(t); } }
  std::string plotLabel = "";
  for (const auto& var : varT) {
    for (const auto& o : info.StrS.at("fitObject")) {
      const auto& obj = info.Par.at("objTag_"+var+"_"+o);
      if (!contain(info.Par, "Model"+var+"_"+obj+label)) { std::cout << "[ERROR] Model tag: " << var+"_"+obj+label << " does not exist!" << std::endl; }
      std::string modelN = info.Par.at("Model"+var+"_"+obj+label);
      modelN.erase(std::remove(modelN.begin(), modelN.end(), ' '), modelN.end());
      stringReplace( modelN, "[", "" ); stringReplace( modelN, "]", "" ); stringReplace( modelN, "+", "_" ); stringReplace( modelN, ",", "_" ); stringReplace( modelN, ";", "_" );
      plotLabel += "Mod"+var+"_" + modelN+"_";
    }
  }
  for (const auto& m : ModelDictionary) { stringReplace(plotLabel, m.first, m.second.second); }; stringReplace(plotLabel, "Cat", ""); stringReplace(plotLabel, "Resolution", "Res");
  // Set file name
  const auto& DSTAG = info.Par.at("DSTAG");
  auto dsTag = DSTAG.substr(0, DSTAG.rfind("_"+info.Par.at("channelDS")));
  if (info.Flag.at("fitMC")) { dsTag += "_"+info.Par.at("PD"); }
  const auto& colTag = DSTAG.substr(DSTAG.find_last_of("_")+1);
  std::string objTag = "";
  for (const auto& o : info.StrS.at("fitObject")) { objTag += o+"_"; }; objTag = objTag.substr(0, objTag.rfind("_"));
  if (contain(info.StrS.at("fitObject"), "Bkg")) { objTag = info.Par.at("objTag_"+*varT.begin()+"_Bkg"); }
  std::string fitVar = "";
  for (const auto& var : varT) { fitVar += var+"_"; }; fitVar = fitVar.substr(0, fitVar.rfind("_"));
  outputDir = Form("%s%s/%s/%s/%s/", outputDir.c_str(), fitVar.c_str(), dsTag.c_str(), objTag.c_str(), colTag.c_str());
  std::string varLbl = "";
  for (const auto& var : info.StrS.at("cutPars")) {
    if (var=="Cand_Mass" || var=="Cand_DLenErr" || var=="Cand_DLenRes" || var=="Cand_DLenGen" || var=="Cand_DLen") continue;
    bool incVar = true; for (const auto& v : varS) { if (var==v) { incVar = false; break; } }
    if (!incVar || !contain(info.Var, var) || !contain(info.Var.at(var), "Min")) continue;
    if (info.Var.at(var).at("Min")==info.Var.at(var).at("Default_Min") && info.Var.at(var).at("Max")==info.Var.at(var).at("Default_Max")) continue;
    if (var=="Centrality" && (DSTAG.rfind("PbPb")==std::string::npos || info.Par.at("PD")=="UPC")) continue;
    auto varN = var; stringReplace(varN, "_", "");
    auto varMin = info.Var.at(var).at("Min");
    auto varMax = info.Var.at(var).at("Max");
    if (contain(info.Flag, "use"+var+"CM") && info.Flag.at("use"+var+"CM")) {
      varN += "CM";
      const bool& ispPb = (info.Flag.at("fitpPb8Y16") || info.Flag.at("fitPA8Y16"));
      varMin = pPb::EtaLABtoCM(varMin, ispPb);
      varMax = pPb::EtaLABtoCM(varMax, ispPb);
    }
    if (varMin==varMax) { varLbl += Form("%s_%.0f_", varN.c_str(), varMax*100.); }
    else { varLbl += Form("%s_%.0f_%.0f_", varN.c_str(), varMin*100., varMax*100.); }
  }
  varLbl = varLbl.substr(0, varLbl.rfind("_"));
  fileName = Form("%s_%s%s", dsTag.c_str(), plotLabel.c_str(), varLbl.c_str());
};


void getRange(double& varMin, double& varMax, const TH1D& hist, const bool& ignoreNegW, const int& nMaxBins, const double& minEntries)
{
  // 1) Find the bin with the maximum Y value
  const auto& binMaximum = hist.GetMaximumBin();
  // 2) Find the start bin
  int startBin = 1;
  if (ignoreNegW) { for (int i = 1; i < hist.GetNbinsX(); i++) { if (hist.GetBinContent(i) > 0.0) { startBin = i; break; } } }
  // 3) Loop backward from maximum Y bin and find the first bin
  int firstBin = startBin;
  for (int i = binMaximum; i >= startBin; i--) {
    if (hist.GetBinContent(i) > 0.1) { firstBin = i; }
    else if (((firstBin-i) > nMaxBins && (hist.GetBinContent(firstBin) < minEntries || hist.GetBinContent(i) < 0)) || (firstBin-i) > 10) break;
  }
  firstBin -= (ignoreNegW ? 0 : 1);
  // 4) Find the end bin
  int endBin = hist.GetNbinsX();
  if (ignoreNegW) { for (int i = hist.GetNbinsX(); i > 0; i--) { if (hist.GetBinContent(i) > 0.0) { endBin = i; break; } } }
  // 5) Loop forward from maximum Y bin and find the last bin
  int lastBin = endBin;
  for (int i = binMaximum; i <= endBin; i++) {
    if (hist.GetBinContent(i) > 0.1) { lastBin = i; }
    else if (((i-lastBin) > nMaxBins && (hist.GetBinContent(lastBin) < minEntries || hist.GetBinContent(i) < 0)) || (i-lastBin) > 10) break;
  }
  // 6) Set the range
  varMin = hist.GetBinLowEdge(firstBin);
  varMax = hist.GetBinLowEdge(lastBin+1)+hist.GetBinWidth(lastBin+1);
};


bool getRange(RooAbsData* ds, RooRealVar& var, double& varMin, double& varMax, const int& nMaxBins, const double& minEntries=5.)
{
  if (!ds) { std::cout << "[ERROR] getRange: DataSet is NULL!" << std::endl; return false; }
  if (ds->getRange(var, varMin, varMax)) { std::cout << "[ERROR] getRange: Error in dataset range!" << std::endl; return false; }
  if (nMaxBins>=0) {
    const auto& varMaxPlot = var.getMax("PlotWindow");
    varMin -= std::abs(varMin)*0.01; varMax += std::abs(varMax)*0.01;
    varMin = varMaxPlot - std::ceil((varMaxPlot-varMin)/var.getBinWidth(0))*var.getBinWidth(0);
    if (std::string(var.GetName())=="Cand_DLenErr" && varMin<=0.0) { varMin = var.getBinWidth(0); }
    varMax = varMin + std::ceil((varMax-varMin)/var.getBinWidth(0))*var.getBinWidth(0);
  }
  if (nMaxBins>=1) {
    const auto& nBins = getNBins(varMin, varMax, var.getBinWidth(0));
    auto hTot = std::unique_ptr<TH1D>(static_cast<TH1D*>(ds->createHistogram("TMP", var, RooFit::Binning(nBins, varMin, varMax))));
    if (!hTot) { std::cout << "[ERROR] getRange: Histogram " << ds->GetName() << " was not created!" << std::endl; return false; }
    getRange(varMin, varMax, *hTot, ds->isWeighted(), nMaxBins, minEntries);
  }
  if (nMaxBins>=0) { std::cout << "[INFO] " << var.GetName() << " range set from data: [ " << varMin << " , " << varMax << " ]" << std::endl; }
  return true; 
};


void addString(RooWorkspace& ws, const std::string& strName, const std::string& strVal, const bool& asObj=true)
{
  if (ws.obj(strName.c_str())) {
    const auto& strVar = dynamic_cast<RooStringVar*>(ws.obj(strName.c_str()));
    if (strVar) { strVar->setVal(strVal.c_str()); }
  }
  else if (ws.arg(strName.c_str())) {
    const auto& strVar = dynamic_cast<RooStringVar*>(ws.arg(strName.c_str()));
    if (strVar) { strVar->setVal(strVal.c_str()); }
  }
  else {
    RooStringVar strVar(strName.c_str(), strName.c_str(), strVal.c_str(), strVal.length()+100);
    if (asObj) { ws.import(*dynamic_cast<TObject*>(&strVar), strVar.GetTitle()); }
    else { ws.import(strVar); }
  }
};


std::string getString(const RooWorkspace& ws, const std::string& strName)
{
  std::string strVal = "";
  if (ws.obj(strName.c_str())) {
    const auto& strVar = dynamic_cast<RooStringVar*>(ws.obj(strName.c_str()));
    if (strVar) { strVal = strVar->getVal(); }
  }
  else if (ws.arg(strName.c_str())) {
    const auto& strVar = dynamic_cast<RooStringVar*>(ws.arg(strName.c_str()));
    if (strVar) { strVal = strVar->getVal(); }
  }
  return strVal;
};


bool setConstant(RooWorkspace& ws, const std::string& parName, const bool& CONST)
{
  if (ws.var(parName.c_str())) { 
    ws.var(parName.c_str())->setConstant(CONST);
    if (CONST) { std::cout << "[INFO] Setting parameter " << parName << " : " << ws.var(parName.c_str())->getVal() << " to constant!" << std::endl; }
  }
  else if (!ws.function(parName.c_str())) { 
    std::cout << "[ERROR] Parameter " << parName << " was not found!" << std::endl;
    return false;
  }
  return true;
};


void setFixedVarsToContantVars(RooWorkspace& ws)
{
  const auto& listVar = ws.allVars();
  auto parIt = std::unique_ptr<TIterator>(listVar.createIterator());
  for (auto itp = parIt->Next(); itp!=NULL; itp = parIt->Next() ) {
    const auto& it = dynamic_cast<RooRealVar*>(itp); if (!it) continue;
    if (it->getMin()==it->getMax() && !it->isConstant()) {
      setConstant(ws, it->GetName(), kTRUE);
    }
  }
};


RooArgSet defineSet(RooWorkspace& ws, const std::string& setName, const RooArgSet& inSet, const std::string& setGroup="SET")
{
  auto set = inSet;
  if (ws.set(setName.c_str())) { set.add(*ws.set(setName.c_str())); }
  if (ws.defineSet(setName.c_str(), set, kTRUE)) { std::cout << "[ERROR] Failed to define set " << setName << std::endl; std::exit(0); }
  if (!ws.arg(setName.c_str())) {
    addString(ws, setName, setName, false);
    if (ws.extendSet(setGroup.c_str(), setName.c_str())) { std::cout << "[ERROR] Failed to extend set " << setGroup << std::endl; std::exit(0); }
  }
  if (setGroup!="SET" && !ws.arg(setGroup.c_str())) {
    addString(ws, setGroup.c_str(), setGroup.c_str(), false);
    if (ws.extendSet("SET", setGroup.c_str())) { std::cout << "[ERROR] Failed to extend set SET" << std::endl; std::exit(0); }
  }
  return set;
};


void defineSet(RooWorkspace& ws, const std::string& setName, const StringSet_t& setStr, const std::string& setGroup="SET")
{
  RooArgSet set;
  for (const auto& s : setStr) { if (ws.arg(s.c_str())) { set.add(*ws.arg(s.c_str())); } }
  defineSet(ws, setName, set, setGroup);
};


void saveSnapshot(RooWorkspace& ws, const RooArgSet& inSet, const std::string& snapName)
{
  if (snapName=="") { std::cout << "[ERROR] Snapshot name is empty!" << std::endl; std::exit(0); }
  // Register the snapshot in WS
  const auto& set = defineSet(ws, snapName, inSet, "snapshots");
  // Save snapshot
  ws.saveSnapshot(snapName.c_str(), set, kTRUE);
};


void saveSnapshot(RooWorkspace& ws, const std::string& snapName, const std::string& dsName="")
{
  if (snapName=="") { std::cout << "[ERROR] Snapshot name is empty!" << std::endl; std::exit(0); }
  // Add the variables and categories
  RooArgSet set(ws.allVars(), snapName.c_str()); set.add(ws.allCats());
  // Remove those from the dataset
  if (dsName!="") { set.remove(*ws.data(dsName.c_str())->get(), false, true); }
  // Save the snapshot
  saveSnapshot(ws, set, snapName);
};


std::set<std::string> MODELCLASS_;
bool importModelClass(RooWorkspace& ws, const std::string& className)
{
  const std::string& CWD = getcwd(NULL, 0);
  const auto& classDir = CWD+"/Macros/Utilities/Models/";
  TInterpreter::EErrorCode ecode;
  if (!contain(MODELCLASS_, className)) {
    gInterpreter->ProcessLineSynch(Form(".L %s%s.cxx+",classDir.c_str(), className.c_str()), &ecode);
    if (ecode!=TInterpreter::kNoError) { std::cout << "[ERROR] Class " << className << " did not compile!" << std::endl; return false; }
    MODELCLASS_.insert(className);
  }
  auto classPdf = std::unique_ptr<RooAbsPdf>((RooAbsPdf*)gInterpreter->ProcessLineSynch(Form("new %s()", className.c_str()), &ecode));
  if (ecode!=TInterpreter::kNoError) { std::cout << "[ERROR] Class " << className << " was not created!" << std::endl; return false; }
  return ws.importClassCode(classPdf->IsA());
};


bool importModelClass(RooWorkspace& ws, RooAbsPdf* pdf)
{
  if (!pdf) return false;
  const std::string& className = pdf->ClassName();
  if (!contain(MODELCLASS_, className)) return true;
  return importModelClass(ws, pdf->ClassName());
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
    for (auto it = catIt->Next(); it!=NULL; it = catIt->Next()) { inWS.RecursiveRemove(it); if(it) delete it; }
  }
  // Clear all category functions
  if (inWS.allCatFunctions().getSize()>0 && delVar) {
    const auto& listFnc = inWS.allCatFunctions();
    auto fncIt = std::unique_ptr<TIterator>(listFnc.createIterator());
    for (auto it = fncIt->Next(); it!=NULL; it = fncIt->Next()) { inWS.RecursiveRemove(it); if(it) delete it; }
  }
  // Clear all variables
  if (inWS.allVars().getSize()>0 && delVar) {
    const auto& listVar = inWS.allVars();
    auto parIt = std::unique_ptr<TIterator>(listVar.createIterator());
    for (auto it = parIt->Next(); it!=NULL; it = parIt->Next()) { inWS.RecursiveRemove(it); if(it) delete it; }
  }
  // Clear all functions
  if (inWS.allFunctions().getSize()>0 && delVar) {
    const auto& listFnc = inWS.allFunctions();
    auto fncIt = std::unique_ptr<TIterator>(listFnc.createIterator());
    for (auto it = fncIt->Next(); it!=NULL; it = fncIt->Next()) { inWS.RecursiveRemove(it); if(it) delete it; }
  }
  // Clear all PDFs
  if (inWS.allPdfs().getSize()>0) {
    const auto& listPdf = inWS.allPdfs();
    auto pdfIt = std::unique_ptr<TIterator>(listPdf.createIterator());
    for (auto it = pdfIt->Next(); it!=NULL; it = pdfIt->Next()) { inWS.RecursiveRemove(it); if(it) delete it; }
  }
  // Clear all Datasets
  if (inWS.allData().size()>0 && smp!="") {
    const auto& listData = inWS.allData();
    for (const auto& it : listData) { if (it && (smp=="All" || std::string(it->GetName()).find(smp)!=std::string::npos)) { inWS.RecursiveRemove(it); if(it) delete it; } }
  }
  // Clear all Embedded Datasets
  if (inWS.allEmbeddedData().size()>0 && smp!="") {
    const auto& listData = inWS.allEmbeddedData();
    for (const auto& it : listData) { if (it && (smp=="All" || std::string(it->GetName()).find(smp)!=std::string::npos)) { inWS.RecursiveRemove(it); if(it) delete it; } }
  }
  // Clear all Generic Objects
  if (inWS.allGenericObjects().size()>0) {
    const auto& listObj = inWS.allGenericObjects();
    for (const auto& it : listObj) { if (it) { inWS.RecursiveRemove(it); } }
  }
  // Clear all sets
  if (const_cast<RooWorkspace*>(&inWS)->set("SET")) {
    const auto& listSet = const_cast<RooWorkspace*>(&inWS)->set("SET");
    auto setIt = std::unique_ptr<TIterator>(listSet->createIterator());
    for (auto it = setIt->Next(); it!=NULL; it = setIt->Next()) { inWS.removeSet(it->GetName()); inWS.RecursiveRemove(it); }
    inWS.removeSet("SET");
  }
  // Clear RooStudyManager
  inWS.clearStudies();
  // Clear any remaining component
  if (inWS.components().getSize()>0 && delVar) {
    auto cmpIt = std::unique_ptr<TIterator>(inWS.componentIterator());
    for (auto it = cmpIt->Next(); it!=NULL; it = cmpIt->Next()) { inWS.RecursiveRemove(it); if(it) delete it; }
  }
  // Return the RooMessenger Level
  RooMsgService::instance().setGlobalKillBelow(level);
};


void copyWorkspace(RooWorkspace& outWS, const RooWorkspace& inWS, const std::string& smp="All", const bool& addVar=true, const bool& addSnap=true, const bool& addPDF=true, const bool& addObj=true)
{
  //
  const auto& level = RooMsgService::instance().globalKillBelow();
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  // Copy all category variables
  if (inWS.allCats().getSize()>0 && addVar) {
    const auto& listCat = inWS.allCats();
    auto catIt = std::unique_ptr<TIterator>(listCat.createIterator());
    for (auto itp = catIt->Next(); itp!=NULL; itp = catIt->Next()) { const auto& it = dynamic_cast<RooCategory*>(itp); if (it && !outWS.cat(it->GetName())) { outWS.import(*it); } }
  }
  // Copy all category functions
  if (inWS.allCatFunctions().getSize()>0 && addVar) {
    const auto& listFnc = inWS.allCatFunctions();
    auto fncIt = std::unique_ptr<TIterator>(listFnc.createIterator());
    for (auto itp = fncIt->Next(); itp!=NULL; itp = fncIt->Next()) { const auto& it = dynamic_cast<RooCategory*>(itp); if (it && !outWS.catfunc(it->GetName())) { outWS.import(*it); } }
  }
  // Copy all variables
  if (inWS.allVars().getSize()>0 && addVar) {
    const auto& listVar = inWS.allVars();
    auto parIt = std::unique_ptr<TIterator>(listVar.createIterator());
    for (auto itp = parIt->Next(); itp!=NULL; itp = parIt->Next()) { const auto& it = dynamic_cast<RooRealVar*>(itp); if (it && !outWS.var(it->GetName())) { outWS.import(*it); } }
  }
  // Copy all functions
  if (inWS.allFunctions().getSize()>0 && addVar) {
    const auto& listFnc = inWS.allFunctions();
    auto fncIt = std::unique_ptr<TIterator>(listFnc.createIterator());
    for (auto itp = fncIt->Next(); itp!=NULL; itp = fncIt->Next()) { const auto& it = dynamic_cast<RooRealVar*>(itp); if (it && !outWS.function(it->GetName())) { outWS.import(*it); } }
  }
  // Copy all PDFs
  if (inWS.allPdfs().getSize()>0 && addPDF) {
    const auto& listPdf = inWS.allPdfs();
    auto pdfIt = std::unique_ptr<TIterator>(listPdf.createIterator());
    for (auto itp = pdfIt->Next(); itp!=NULL; itp = pdfIt->Next() ) { const auto& it = dynamic_cast<RooAbsPdf*>(itp); if (it && !outWS.pdf(it->GetName()) && importModelClass(outWS, it)) {
	outWS.import(*it, RooFit::RecycleConflictNodes()); }
    }
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
  if (inWS.allGenericObjects().size()>0 && addObj) {
    const auto& listObj = inWS.allGenericObjects();
    for (const auto& it : listObj) { if (!outWS.genobj(it->GetTitle())) { outWS.import(*it, it->GetTitle()); } }
  }
  // Copy all snapshots variables
  if (const_cast<RooWorkspace*>(&inWS)->set("snapshots") && addSnap) {
    const auto& listSnap = const_cast<RooWorkspace*>(&inWS)->set("snapshots");
    auto snapIt = std::unique_ptr<TIterator>(listSnap->createIterator());
    for (auto it = snapIt->Next(); it!=NULL; it = snapIt->Next()) { if (!outWS.set(it->GetName())) { saveSnapshot(outWS, *inWS.getSnapshot(it->GetName()), it->GetName()); } }
  }
  // Copy all sets
  if (const_cast<RooWorkspace*>(&inWS)->set("SET")) {
    const auto& listSet = const_cast<RooWorkspace*>(&inWS)->set("SET");
    auto setIt = std::unique_ptr<TIterator>(listSet->createIterator());
    for (auto it = setIt->Next(); it!=NULL; it = setIt->Next()) { if (!outWS.arg(it->GetName())) { addString(outWS, it->GetName(), it->GetName(), false); outWS.extendSet("SET", it->GetName()); } }
  }
  // Return the RooMessenger Level
  RooMsgService::instance().setGlobalKillBelow(level);
};


bool saveDataSet(const RooDataSet& ds, const std::string& fileName)
{
  // Create the output file
  const auto& outputDir = fileName.substr(0, fileName.rfind("/")+1);
  gSystem->mkdir(outputDir.c_str(), kTRUE);
  auto file = std::unique_ptr<TFile>(new TFile(fileName.c_str(), "RECREATE"));
  if (!file || !file->IsOpen() || file->IsZombie()) {
    std::cout << "[ERROR] Output root file with fit results could not be created!" << std::endl; if (file) { file->Close(); }; return false;
  }
  file->cd();
  // Save the dataset
  ds.Write(ds.GetName());
  // Write and close output file
  file->Write(); file->Close();
  std::cout << "[INFO] RooDataSet " << ds.GetName() << " stored in: " << fileName << std::endl;
  return true;
};


bool getDataSet(RooWorkspace& ws, const std::string& dsName, const std::string& fileName)
{
  std::cout << "[INFO] Loading RooDataSet " << dsName << " from: " << fileName << std::endl;
  // Check if dataset already exist in workspace
  if (ws.data(dsName.c_str())) { std::cout << "[INFO] DataSet " << dsName << " already exist, will not load it!" << std::endl; return true; }
  // Check if input file exists
  if (gSystem->AccessPathName(fileName.c_str())) {
    std::cout << "[INFO] File: " << fileName << " not found!" << std::endl; return false;
  }
  // Open the input file
  auto file = std::unique_ptr<TFile>(TFile::Open(fileName.c_str()));
  if (!file || !file->IsOpen() || file->IsZombie()) {
    std::cout << "[ERROR] File: " << fileName << " is corrupted!" << std::endl; if(file) { file->Close(); }; return false;
  }
  file->cd();
  TH1::AddDirectory(kFALSE); // Avoid adding objects to memory
  const auto& inDS = dynamic_cast<RooDataSet*>(file->Get(dsName.c_str()));
  if (!inDS) { std::cout << "[ERROR] DataSet " << dsName << " in: " << fileName << " not found!" << std::endl; file->Close(); return false; }
  if (ws.import(*inDS)) { std::cout << "[ERROR] DataSet " << dsName << " in: " << fileName << " was not imported!" << std::endl; file->Close(); return false; }
  // Close input file
  file->Close();
  return true;
};


bool getPDFData(RooWorkspace& ws, const std::string& pdfName, const std::string& fileName)
{
  auto dataName = pdfName;
  dataName.replace(dataName.find("pdf"), std::string("pdf").length(), "dh");
  std::cout << "[INFO] Extracting RooDataHist " << dataName << " from: " << fileName << std::endl;
  // Check if dataset already exist in workspace
  if (ws.data(dataName.c_str())) { std::cout << "[INFO] RooDataHist " << dataName << " already exist, will not load it!" << std::endl; return true; }
  // Check if input file exists
  if (gSystem->AccessPathName(fileName.c_str())) {
    std::cout << "[INFO] File: " << fileName << " not found!" << std::endl; return false;
  }
  // Open the input file
  auto file = std::unique_ptr<TFile>(TFile::Open(fileName.c_str()));
  if (!file || !file->IsOpen() || file->IsZombie()) {
    std::cout << "[ERROR] File: " << fileName << " is corrupted!" << std::endl; if(file) { file->Close(); }; return false;
  }
  file->cd();
  TH1::AddDirectory(kFALSE); // Avoid adding objects to memory
  const auto& inData = dynamic_cast<RooDataHist*>(file->Get(dataName.c_str()));
  if (!inData) { std::cout << "[ERROR] RooDataHist " << dataName << " in: " << fileName << " not found!" << std::endl; file->Close(); return false; }
  if (ws.import(*inData)) { std::cout << "[ERROR] RooDataHist " << dataName << " in: " << fileName << " was not imported!" << std::endl; file->Close(); return false; }
  // Close input file
  file->Close();
  return true;
};


bool getPDFData(RooWorkspace& ws, const std::string& pdfName, const GlobalInfo& info, const std::string& label, const std::string& var)
{
  std::string fileName = "";
  auto dir = info.Par.at("outputDir");
  setFileName(fileName, dir, label, info, {var});
  auto varT = var; stringReplace(varT, "_", "");
  const auto& filePath = (dir+"result/FIT_"+varT+"_"+fileName+".root");
  return getPDFData(ws, pdfName, filePath);
};


bool saveWorkSpace(const RooWorkspace& ws, const std::string& outputDir, const std::string& fileName, const bool& saveDS=true)
{
  // Create the output file
  gSystem->mkdir(outputDir.c_str(), kTRUE);
  auto file = std::unique_ptr<TFile>(new TFile((outputDir+fileName).c_str(), "RECREATE"));
  if (!file || !file->IsOpen() || file->IsZombie()) {
    std::cout << "[ERROR] Output root file with fit results could not be created!" << std::endl; if (file) { file->Close(); }; return false;
  }
  file->cd();
  // Save the workspace
  if (saveDS) { ws.Write("workspace"); }
  else {
    RooWorkspace tmpWS;
    copyWorkspace(tmpWS, ws, (saveDS ? "All" : ""));
    tmpWS.Write("workspace");
  }
  // Save the snapshots for later checks (BUG FIX)
  if (const_cast<RooWorkspace*>(&ws)->set("snapshots")) {
    const auto& listSnap = const_cast<RooWorkspace*>(&ws)->set("snapshots");
    auto snapIt = std::unique_ptr<TIterator>(listSnap->createIterator());
    for (auto its = snapIt->Next(); its!=NULL; its = snapIt->Next()) {
      const auto& it = dynamic_cast<RooStringVar*>(its); if (!it) continue;
      const auto& snap = ws.getSnapshot(it->GetName());
      if (snap) { snap->Write(it->GetName()); }
    }
  }
  // Save the RooDataHists
  const auto& listData = ws.allData();
  for (const auto& itp : listData) {
    const auto& it = dynamic_cast<RooDataHist*>(itp);
    if (it) { it->Write(it->GetName()); }
  }
  // Write and close output file
  file->Write(); file->Close();
  std::cout << "[INFO] RooWorkspace saved in: " << (outputDir+fileName) << std::endl;
  return true;
};


bool getWorkSpace(RooWorkspace& ws, const std::string& fileName, const std::string& getDS="All", const bool& onlySnap=false)
{
  // Check if input file exists
  if (gSystem->AccessPathName(fileName.c_str())) {
    std::cout << "[INFO] File: " << fileName << " not found!" << std::endl; return false;
  }
  // Open the input file
  auto file = std::unique_ptr<TFile>(TFile::Open(fileName.c_str()));
  if (!file || !file->IsOpen() || file->IsZombie()) {
    std::cout << "[ERROR] File: " << fileName << " is corrupted!" << std::endl; if(file) { file->Close(); }; return false;
  }
  file->cd();
  TH1::AddDirectory(kFALSE); // Avoid adding objects to memory
  if (onlySnap) {
    TIter keyIt(file->GetListOfKeys());
    for (auto itk = keyIt(); itk!=NULL; itk = keyIt() ) {
      const auto& key = dynamic_cast<TKey*>(itk); if (!key) continue;
      if (std::string(key->GetClassName())!="RooArgSet") continue;
      const auto& snap = dynamic_cast<RooArgSet*>(key->ReadObj());
      if (snap) { saveSnapshot(ws, *snap, key->GetName()); }
    }
  }
  else {
    const auto& inWS = dynamic_cast<RooWorkspace*>(file->Get("workspace"));
    if (!inWS) { std::cout << "[ERROR] Workspace in: " << fileName << " not found!" << std::endl; file->Close(); return false; }
    copyWorkspace(ws, *inWS, getDS);
  }
  // Close input file
  file->Close();
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
  auto parIt2 = std::unique_ptr<TIterator>(pars1.createIterator());
  for (auto itp = parIt2->Next(); itp!=NULL; itp = parIt2->Next() ) {
    const auto& it = dynamic_cast<RooRealVar*>(itp); if (!it) continue;
    if (!pars2.find(it->GetName())) continue;
    const auto& val = pars2.getRealValue(it->GetName(),-1e99);
    if (val != it->getVal() || ((RooRealVar&)pars2[it->GetName()]).getMin() != it->getMin() || ((RooRealVar&)pars2[it->GetName()]).getMax() != it->getMax()){
      std::cout << "OLD: " << it->GetName() << "  " << ((RooRealVar&)pars2[it->GetName()]).getVal() << "  " << ((RooRealVar&)pars2[it->GetName()]).getMin() << "  " << ((RooRealVar&)pars2[it->GetName()]).getMax() << std::endl;
      std::cout << "NEW: " << it->GetName() << "  " << it->getVal() << "  " << it->getMin() << "  " << it->getMax() << std::endl;
    }
  }
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


bool isFitAlreadyFound(const RooArgSet& newpars, const std::string& fileName, const std::string& snapName="initialParameters")
{
  RooWorkspace ws; if (!getWorkSpace(ws, fileName, "", true)) { return false; }
  const RooArgSet* params = (ws.set(snapName.c_str()) ? ws.getSnapshot(snapName.c_str()) : NULL);
  if (!params) { std::cout << "[INFO] Snapshot " << snapName << " was not found!" << std::endl; return false; }
  return compareSnapshots(newpars, *params);
};


RooRealVar getVar(const RooArgSet& set, const std::string& varName)
{
  if (set.find(varName.c_str())) {
    const auto& it = dynamic_cast<RooRealVar*>(set.find(varName.c_str()));
    if (it) { return *it; }
  }
  std::cout << "[WARNING] Variable " << varName << " was not found in the set, return empty RooRealVar" << std::endl;
  return RooRealVar();
};


bool loadParameterRange(GlobalInfo& info, const std::string& var, const std::string& fileName, const std::string& snap="fittedParameters")
{
  RooWorkspace ws; if (!getWorkSpace(ws, fileName, "", true)) { return false; }
  ws.loadSnapshot(snap.c_str());
  if (ws.var(var.c_str())) {
    info.Var.at(var).at("Min") = ws.var(var.c_str())->getMin();
    info.Var.at(var).at("Max") = ws.var(var.c_str())->getMax();
  }
  else { std::cout << "[ERROR] Variable " << var << " was not found!" << std::endl; return false; }
  return true;
};


bool loadVarRange(GlobalInfo& info, const std::string& var, const std::string& label, const std::string& snap="fittedParameters")
{
  std::string fileName = "";
  auto dir = info.Par.at("outputDir");
  setFileName(fileName, dir, label, info, {var});
  auto varT = var; stringReplace(varT, "_", "");
  const auto& filePath = (dir+"result/FIT_"+varT+"_"+fileName+".root");
  std::cout << "[INFO] Loading " << var << " parameter range from: " << fileName << std::endl;
  const auto& status = loadParameterRange(info, var, filePath, snap);
  if (status) {
    std::cout << "[INFO] Parameter range loaded: " << var << " > [ " << info.Var.at(var).at("Min") << " , " << info.Var.at(var).at("Max") << " ]" << std::endl;
  }
  return status;
};


bool loadParameters(RooWorkspace& ws, const std::string& fileName, const std::string& varName="*", const std::string& col="", 
		    const std::string& dsName="", const std::string& snap="fittedParameters")
{
  std::cout << "[INFO] Loading model parameters from snapshot " << snap << " in " << fileName << std::endl;
  RooWorkspace inWS; if (!getWorkSpace(inWS, fileName, dsName, (dsName==""))) { return false; }
  const auto& inDS = dynamic_cast<RooDataSet*>(inWS.data(dsName.c_str()));
  const auto& ds = dynamic_cast<RooDataSet*>(ws.data(dsName.c_str()));
  if (dsName!="" && (!inDS || !ds)) {
    std::cout << "[INFO] loadParameters: RooDataset " << dsName << " was not found!" << std::endl; return false;
  }
  else if (dsName!="" && !isCompatibleDataset(*ds, *inDS)) {
    std::cout << "[INFO] loadParameters: RooDataset " << dsName << " is different!" << std::endl; return false;
  }
  const auto& params = inWS.getSnapshot(snap.c_str());
  if (!params) { std::cout << "[INFO] Snapshot " << snap << " was not found!" << std::endl; return false; }
  auto parIt = std::unique_ptr<TIterator>(params->createIterator());
  for (auto itp = parIt->Next(); itp!=NULL; itp = parIt->Next() ) {
    const auto& it = dynamic_cast<RooRealVar*>(itp); if (!it) continue;
    const std::string& name = it->GetName();
    if ((varName!="*" && name.rfind(varName,0)!=0) || (col!="" && name.rfind(col)==std::string::npos)) continue;
    if (ws.var(name.c_str())) {
      ws.var(name.c_str())->setVal(it->getVal());
      ws.var(name.c_str())->setError(it->getError());
      ws.var(name.c_str())->setMin(it->getMin());
      ws.var(name.c_str())->setMax(it->getMax());
    }
    else {
      ws.factory((name+Form("[ %.2f, %.2f, %.2f ]", it->getVal(), it->getMin(), it->getMax())).c_str());
    }
    std::cout << "[INFO] Parameter " << name << " set to " << it->getVal() << "+-" << it->getError() << std::endl;
  }
  return true;
};


bool loadFitResult(RooWorkspace& ws, const std::string& fileName, const std::string& pdfName, const std::string& var,
		   const bool& loadYield=true, const bool& fixPar=true, const std::string& snapName="fittedParameters")
{
  std::cout << "[INFO] Loading parameters for " << pdfName << " from " << fileName << std::endl;
  const auto& loadedFileName = getString(ws, "loadedFile_"+var);
  if (loadedFileName==fileName) {
    std::cout << "[INFO] Parameters from " << fileName << " has already been loaded!" << std::endl;
    return true;
  }
  // Extract the workspace
  RooWorkspace inWS;
  if (!ws.pdf(pdfName.c_str())) {
    if (!getWorkSpace(inWS, fileName, "")) { return false; }
    const auto& pdf = inWS.pdf(pdfName.c_str());
    if (!pdf) { std::cout << "[ERROR] PDF " << pdfName << " was not found in " << fileName << std::endl; return false; }
    if (ws.import(*pdf, RooFit::RecycleConflictNodes())) { std::cout << "[ERROR] PDF " << pdfName << " was not imported!" << std::endl; return false; }
  }
  else if (!getWorkSpace(inWS, fileName, "", true)) { return false; }
  // Extract PDF
  const auto& pdf = ws.pdf(pdfName.c_str());
  if (!pdf) { std::cout << "[INFO] PDF " << pdfName << " was not found!" << std::endl; return false; }
  // Extract the snapshot
  const RooArgSet* snapPar = (inWS.set(snapName.c_str()) ? inWS.getSnapshot(snapName.c_str()) : NULL);
  if (!snapPar) { std::cout << "[INFO] Snapshot " << snapName << " was not found!" << std::endl; return false; }
  // Loop over PDF variables and load from snapshot
  std::string res = ""; RooArgSet loadSet;
  auto params = std::unique_ptr<RooArgSet>(pdf->getVariables());
  auto parIt = std::unique_ptr<TIterator>(params->createIterator());
  for (auto it = parIt->Next(); it!=NULL; it = parIt->Next()) {
    const std::string& name = it->GetName();
    if (!ws.var(name.c_str())) continue;
    if (!loadYield && (name.rfind("N_",0)==0 || name.rfind("R_",0)==0)) continue;
    const auto& sPar = dynamic_cast<RooRealVar*>(snapPar->find(name.c_str()));
    double val = -99999.;
    double unc = -99999.;
    if (name.rfind("N_",0)==0 && !sPar && snapPar->find(("R_"+name.substr(2)).c_str())) {
      const auto& vN = name.substr(2);
      const auto& obj = vN.substr(0, vN.find("To"));
      auto rN = vN; stringReplace(rN, obj, (obj=="Psi2S"?"JPsi":(obj=="Ups2S"?"Ups1S":(obj=="Ups3S"?"Ups1S":obj))));
      const auto& rPar = dynamic_cast<RooRealVar*>(snapPar->find(("R_"+vN).c_str()));
      const auto& nPar = dynamic_cast<RooRealVar*>(snapPar->find(("N_"+rN).c_str()));
      if (!nPar) { std::cout << "[ERROR] Parameter " << ("N_"+rN) << " was not found in snapshot " << snapName << std::endl; return false; }
      val = rPar->getVal()*nPar->getVal();
      unc = std::sqrt(std::pow(rPar->getVal()*nPar->getError(), 2.0) + std::pow(rPar->getError()*nPar->getVal(), 2.0));
    }
    else if (sPar) {
      val = sPar->getVal(); if (name.rfind("N_",0)==0) { val = std::max(val, 0.0); }
      unc = sPar->getError();
    }
    if (val!=-99999.) {
      ws.var(name.c_str())->setVal(val);
      if (unc!=-99999.) { ws.var(name.c_str())->setError(unc); }
      if (name.rfind("N_",0)!=0 && name.rfind("R_",0)!=0 && fixPar) ws.var(name.c_str())->setConstant(true);
      std::string cnst = ""; if (ws.var(name.c_str())->isConstant()) { cnst = "[C]"; }
      res += Form("%s = %g %s, ", name.c_str(), val, cnst.c_str());
      loadSet.add(*ws.var(name.c_str()));
    }
  }
  saveSnapshot(ws, loadSet, "loadedParameters");
  addString(ws, "loadedFile_"+var, fileName);
  std::cout << "[INFO] Loaded " << pdfName << " parameters: " << res << std::endl;
  return true;
};


bool createBinnedDataset(RooWorkspace& ws, const std::string& dsName, const std::string& var)
{
  if (ws.data(dsName.c_str())==NULL) { std::cout << "[ERROR] createBinnedDataset: DataSet " << dsName << " was not found!" << std::endl; return false; }
  if (ws.data(dsName.c_str())->numEntries()<=10.0) { std::cout << "[ERROR] createBinnedDataset: DataSet " << dsName << " has too few events!" << std::endl; return false; }
  if (ws.var(var.c_str())==NULL) { std::cout << "[ERROR] createBinnedDataset: Variable " << var << " was not found!" << std::endl; return false; }
  //
  const auto& min  = ws.var(var.c_str())->getMin();
  const auto& max  = ws.var(var.c_str())->getMax();
  const uint& nBin = ws.var(var.c_str())->getBins();
  //
  // Create the histogram
  auto histName = "h" + dsName + "_" + var;
  std::unique_ptr<TH1D> hist = std::unique_ptr<TH1D>(static_cast<TH1D*>(ws.data(dsName.c_str())->createHistogram(histName.c_str(), *ws.var(var.c_str()), RooFit::Binning(nBin, min, max))));
  if (hist==NULL) { std::cout << "[ERROR] createBinnedDataset: Histogram " << histName << " is NULL!" << std::endl; return false; }
  // Cleaning the input histogram
  // 1) Remove the Under and Overflow bins
  hist->ClearUnderflowAndOverflow();
  // 2) Set negative bin content to zero
  for (int i=0; i<=hist->GetNbinsX(); i++) { if (hist->GetBinContent(i)<0.0) { hist->SetBinContent(i, 0.0); } }
  // 2) Reduce the range of histogram and rebin it
  //hist.reset(static_cast<TH1D*>(rebinhist(*hist, range[1], range[2])));
  if (hist==NULL) { std::cout << "[ERROR] createBinnedDataset: Cleaned Histogram of " << histName << " is NULL!" << std::endl; return false; }
  const auto& dataName = (dsName +"_"+var+"_FIT");
  std::unique_ptr<RooDataHist> dataHist = std::unique_ptr<RooDataHist>(new RooDataHist(dataName.c_str(), "", *ws.var(var.c_str()), hist.get()));
  if (dataHist==NULL) { std::cout << "[ERROR] createBinnedDataset: DataHist used to create " << dsName << " failed!" << std::endl; return false; }
  if (dataHist->sumEntries()==0) { std::cout << "[ERROR] createBinnedDataset: DataHist used to create " << dsName << " is empty!" << std::endl; return false; }
  if (std::abs(dataHist->sumEntries() - hist->GetSumOfWeights())>0.001) { std::cout << "[ERROR] createBinnedDataset: DataHist used to create " << dsName << "  " << " is invalid!  " << std::endl; return false; }
  if (ws.import(*dataHist)) { std::cout << "[ERROR] createBinnedDataset: DataHist " << dataName << " was not imported!" << std::endl; return false; }
  ws.var(var.c_str())->setBins(nBin); // Bug Fix
  return true;
};

//
//---------------------------------------------------------------------------------------------
//


bool setModel(GlobalInfo&  info)
{
  const auto& cha = info.Par.at("channel");
  for (const auto& col : info.StrS.at("fitSystem")) {
    for (const auto& oo : info.StrS.at("fitObject")) {
      for (const auto& chg : info.StrS.at("fitCharge")) {
	for (const auto& var : info.StrS.at("incVarName")) {
	  const auto& obj = info.Par.at("objTag_"+var+"_"+oo);
	  const std::string& label = Form("Model%s_%s_%s", var.c_str(), (obj+cha+chg).c_str(), col.c_str());
	  info.StrS["tags"].insert(obj+cha+chg+"_"+col);
	  const std::string inputLabel = "Model"+findLabel("Model", var, obj, chg, col, cha, info);
	  if (contain(info.Par, inputLabel)) {
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
		if (info.Flag.at("fitMC") && ll!=obj && modelName=="TEMP") continue;
		if (modelName=="TEMP") { modelName = "Template";    }
		if (modelName=="MJET") { modelName = "MultiJetBkg"; }
		if (!contain(ModelDictionary, modelName) || ModelDictionary.at(modelName).first==0) {
		  std::cout << "[ERROR] The " << (ll+cha+chg) << " " << var << " model: " << modelName << " is invalid" << std::endl; return false;
		}
		info.Par[label+"_"+ll] = modelName;
		if (modelName == "Template") {
		  const auto& dsTag = ( "MC_" + ll + "_" + info.Par.at("channelDS") + "_" + col );
		  info.StrS["TEMPDS_"+var+"_"+(obj+cha+chg)+"_"+col].insert(dsTag);
		}
		info.StrS["addObjectModel_"+var+"_"+(obj+cha+chg)+"_"+col].insert(ll);
	      }
	    }
	  } else {
	    std::cout << "[ERROR] " << (obj+cha+chg) << " " << var << " model for " << col << " was not found in the initial parameters!" << std::endl; return false;
	  }
	}
      }
    }
  }
  return true;
};


bool storeDSRange(RooWorkspace& ws, RooAbsData* ds, const std::string& label="DSWindow")
{
  auto varIt = std::unique_ptr<TIterator>(ds->get()->createIterator());
  for (auto itp = varIt->Next(); itp!=NULL; itp = varIt->Next()) {
    const auto& it = dynamic_cast<RooRealVar*>(itp); if (!it) continue;
    double varMin, varMax; if (!getRange(ds, *it, varMin, varMax, -1)) { return false; }
    const auto& var = ws.var(it->GetName());
    if (!var) { std::cout << "[ERROR] Variable " << it->GetName() << " was not found!" << std::endl; return false; }
    var->setRange(label.c_str(), varMin, varMax);
  }
  return true;
};


bool storeDSStat(RooWorkspace& ws, const std::string& dsName)
{
  const auto& ds = dynamic_cast<RooDataSet*>(ws.data(dsName.c_str()));
  if (!ds) { std::cout << "[ERROR] Dataset " << dsName << " was not found!" << std::endl; return false; }
  auto varIt = std::unique_ptr<TIterator>(ds->get()->createIterator());
  for (auto itp = varIt->Next(); itp!=NULL; itp = varIt->Next()) {
    const auto& it = dynamic_cast<RooRealVar*>(itp); if (!it) continue;
    auto mean = std::unique_ptr<RooRealVar>(ds->meanVar(*it));
    if (!ws.var(mean->GetName()) && ws.import(*mean)) { std::cout << "[ERROR] Failed to import " << mean->GetName() << std::endl; return false; }
    auto rms = std::unique_ptr<RooRealVar>(ds->rmsVar(*it));
    if (!ws.var(rms->GetName()) && ws.import(*rms)) { std::cout << "[ERROR] Failed to import " << rms->GetName() << std::endl; return false; }
  }
  if (ws.var("Cand_Rap")) {
    auto rapDS = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(ds->reduce(RooArgSet(*ws.var("Cand_Rap")))));
    for (const auto& v : StringMap_t({{"Cand_AbsRap", "abs(Cand_Rap)"}, {"Cand_RapCM", "Cand_Rap-0.465"}})) {
      if (ws.var(v.first.c_str())) {
	RooFormulaVar rapVar(v.first.c_str(), v.second.c_str(), RooArgList(*ws.var("Cand_Rap")));
	rapDS->addColumn(rapVar);
	if (rapDS->get()->find(v.first.c_str())) {
	  const auto& var = dynamic_cast<RooRealVar*>(ws.var(v.first.c_str()));
	  ws.factory(Form("%sMean[%g]", v.first.c_str(), rapDS->meanVar(*var)->getVal()));
	  ws.factory(Form("%sRMS[%g]", v.first.c_str(), rapDS->rmsVar(*var)->getVal()));
	}
      }
    }
  }
  return true;
};


void setDSParamaterRange(RooDataSet& ds, const GlobalInfo& info, const std::string& tag="")
{
  const std::string& dsName = ds.GetName();
  const auto& row = ds.get();
  for (const auto& var : info.Var) {
    const bool& isAbs = (var.first.find("Abs")!=std::string::npos);
    auto varN = var.first; if (isAbs) { varN.erase(varN.find("Abs"), 3); }
    if (row->find(varN.c_str()) && (contain(var.second, tag+"Min"))) {
      const auto& v = dynamic_cast<RooRealVar*>(row->find(varN.c_str()));
      v->setMin(isAbs ? -var.second.at(tag+"Max") : var.second.at(tag+"Min"));
      v->setMax(var.second.at(tag+"Max"));
    }
  }
};


void setDSParamaterRange(RooWorkspace& ws, const std::string& dsName, const GlobalInfo& info, const std::string& tag="")
{
  const auto& ds = dynamic_cast<RooDataSet*>(ws.data(dsName.c_str()));
  if (!ds) { throw std::runtime_error("[ERROR] Dataset "+dsName+" was not found!"); }
  defineSet(ws, "SET_"+dsName, *ds->get());
  setDSParamaterRange(*ds, info, tag);
};


void setDefaultRange(RooWorkspace& ws, const GlobalInfo& info)
{
  for (const auto& v : info.Var) {
    if (!ws.var(v.first.c_str()) || !contain(v.second, "Default_Min")) continue;
    ws.var(v.first.c_str())->setRange("DEFAULT", v.second.at("Default_Min"), v.second.at("Default_Max"));
    if (contain(v.second, "binWidth")) {
      const auto& nBins = getNBins(v.second.at("Default_Min"), v.second.at("Default_Max"), v.second.at("binWidth"));
      ws.var(v.first.c_str())->setBins(nBins, "DEFAULT");
    }
  }
};


bool setFitParameterRange(RooWorkspace& ws, GlobalInfo& info)
{
  // Store the default range
  setDefaultRange(ws, info);
  // Set the fit parameter range
  for (const auto& var : info.StrS.at("fitCondVariable")) {
    if (!ws.var(var.c_str())) { std::cout << "[ERROR] Parameter " << var << " does not exist, failed to set fit parameter range!" << std::endl; return false; }
    const auto& varMin = info.Var.at(var).at("Min");
    const auto& varMax = info.Var.at(var).at("Max");
    const auto& varNBins = getNBins(varMin, varMax, info.Var.at(var).at("binWidth"));
    ws.var(var.c_str())->setRange("FitWindow", varMin, varMax);
    ws.var(var.c_str())->setBins(varNBins, "FitWindow");
    const auto& varFullMin = info.Var.at(var).at("Full_Min");
    const auto& varFullMax = info.Var.at(var).at("Full_Max");
    const auto& varFullNBins = getNBins(varFullMin, varFullMax, info.Var.at(var).at("binWidth"));
    ws.var(var.c_str())->setRange("FullWindow", varFullMin, varFullMax);
    ws.var(var.c_str())->setBins(varFullNBins, "FullWindow");
    const auto& varPlotMin = info.Var.at(var).at("Plot_Min");
    const auto& varPlotMax = info.Var.at(var).at("Plot_Max");
    const auto& varPlotNBins = getNBins(varPlotMin, varPlotMax, info.Var.at(var).at("binWidth"));
    ws.var(var.c_str())->setRange("PlotWindow", varPlotMin, varPlotMax);
    ws.var(var.c_str())->setBins(varPlotNBins, "PlotWindow");
    const auto& nBins = getNBins(ws.var(var.c_str())->getMin(), ws.var(var.c_str())->getMax(), info.Var.at(var).at("binWidth"));
    ws.var(var.c_str())->setBins(nBins);
  }
  return true;
};


bool updateFitRange(RooWorkspace& ws, GlobalInfo&  info, const std::string& chg)
{
  DoubleDiMap_t inVarList;
  for (const auto& varN : info.StrS.at("fitVariable")) {
    int emptyBins = -1;
    double varMin = -99999.9, varMax = -99999.9;
    if (info.Flag.at("fitData")) {
      emptyBins = 5;
      if (varN.rfind("Cand_DLen", 0)==0) {
	varMin = -999.9;
	varMax = -999.9;
      }
      if (varN=="Cand_DLenErr") {
	emptyBins = 0;
      }
      if (varN=="Cand_DLenRes") {
	varMin = (info.Flag.at("fitPsi2S") ? -999.9 : info.Var.at(varN).at("Plot_Min"));
	varMax = 0.0;
      }
      if (varN=="Cand_DLen" && info.Flag.at("fitBkg")) {
	varMin = info.Var.at(varN).at("Plot_Min");
	varMax = info.Var.at(varN).at("Plot_Max");
      }
    }
    else if (info.Flag.at("fitMC")) {
      emptyBins = 5;
      if (varN=="Cand_DLenErr" || varN=="Cand_DLen") {
	varMin = -999.9;
	varMax = -999.9;
      }
    }
    inVarList[varN] = {{"Min", varMin}, {"Max", varMax}, {"EmptyBins", emptyBins}};
  }
  //
  std::string cutDS = "";
  DoubleDiMap_t varList;
  const auto& dsName = info.Par.at("dsName"+chg);
  const auto& ds = ws.data(dsName.c_str());
  const auto& dsVars = ds->get();
  for (const auto& v : inVarList) {
    const auto& varN = v.first;
    if (!dsVars->find(varN.c_str())) continue;
    if (v.second.at("Min")==-99999.9 && v.second.at("Max")==-99999.9) continue;
    double vMin=0., vMax=0.;
    if ((v.second.at("Min")==-999.9 || v.second.at("Max")==-999.9) && !getRange(ds, *ws.var(varN.c_str()), vMin, vMax, v.second.at("EmptyBins"))) { return false; }
    const double& varMin = ((v.second.at("Min") != -999.9) ? ((v.second.at("Min") != -99999.9) ? v.second.at("Min") : ws.var(varN.c_str())->getMin()) : vMin);
    const double& varMax = ((v.second.at("Max") != -999.9) ? ((v.second.at("Max") != -99999.9) ? v.second.at("Max") : ws.var(varN.c_str())->getMax()) : vMax);
    if (varMin==varMax) { cutDS += Form("(%s == %g) && ", varN.c_str(), varMax); }
    else { cutDS += Form("(%g <= %s && %s < %g) && ", varMin, varN.c_str(), varN.c_str(), varMax); }
    info.StrS.at("cutPars").insert(varN);
    varList[varN] = {{"Min", varMin}, {"Max", varMax}};
  }
  if (varList.empty()) return true;
  cutDS = cutDS.substr(0, cutDS.rfind(" && "));
  addString(ws, "cutDSFit", cutDS); // Save the cut expression for bookkeeping
  //
  for (const auto& v : varList) {
    const auto& nBins = getNBins(v.second.at("Min"), v.second.at("Max"), ws.var(v.first.c_str())->getBinWidth(0));
    ws.var(v.first.c_str())->setRange("FitWindow", v.second.at("Min"), v.second.at("Max"));
    ws.var(v.first.c_str())->setBins(nBins, "FitWindow");
    ws.var(v.first.c_str())->setRange(v.second.at("Min"), v.second.at("Max"));
    ws.var(v.first.c_str())->setBins(nBins);
    info.Var.at(v.first.c_str()).at("Min") = v.second.at("Min");
    info.Var.at(v.first.c_str()).at("Max") = v.second.at("Max");
    std::cout << "[INFO] " << v.first << " range was updated to : " << ws.var(v.first.c_str())->getMin() << " <= " << v.first << " < " << ws.var(v.first.c_str())->getMax() << std::endl;
  }
  //
  const auto& listData = ws.allData();
  for (const auto& ds : listData) {
    std::string dsN = ds->GetName();
    const auto& tmpDS = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(ds->reduce(RooFit::Cut(cutDS.c_str()))));
    if (!tmpDS) { std::cout << "[ERROR] Reduced RooDataSet " << dsN << " was not created!" << std::endl; return false; }
    if (tmpDS->numEntries() > 0 && tmpDS->numEntries() < ds->numEntries()) {
      std::cout << "[INFO] Applying cut: " << cutDS << " on dataset: " << dsN << std::endl;
      if (dsN==dsName) {
	dsN += "_FIT";
	info.Par.at("dsNameFit"+chg) = dsN;
	if (contain(info.Par, "dsSPlotName"+chg) && info.Par.at("dsSPlotName"+chg)==dsName) { info.Par.at("dsSPlotNameFit"+chg) = dsN; }
	if (ws.import(*tmpDS, RooFit::Rename(dsN.c_str()))) { std::cout << "[ERROR] RooDataSet " << dsN << " was not created!" << std::endl; return false; }
	info.Var.at("numEntries")[dsN] = tmpDS->sumEntries();
	ws.factory(Form("numEntries_%s[%.0f]", dsN.c_str(), tmpDS->sumEntries()));
      }
      else {
	ws.RecursiveRemove(ds); if(ds) delete ds;
	if (ws.import(*tmpDS, RooFit::Rename(dsN.c_str()))) { std::cout << "[ERROR] RooDataSet " << dsN << " was not created!" << std::endl; return false; }
      }
      std::cout << "[INFO] RooDataSet " << dsN << " reduced from " << ds->numEntries() << " to " << tmpDS->numEntries() << " events" << std::endl;
      setDSParamaterRange(ws, dsN, info);
    }
    else if (tmpDS->numEntries() <= 0) {  std::cout << "[ERROR] Reduced RooDataSet " << dsN << " has 0 entries!" << std::endl; return false; }
    else if (tmpDS->numEntries() > ds->numEntries()) { std::cout << "[ERROR] Reduced RooDataSet " << dsN << " has more entries than original!" << std::endl; return false; }
  }
  //
  return true;
};


bool fitPDF(std::unique_ptr<RooFitResult>& fitResult, RooWorkspace& ws, const std::vector<RooCmdArg>& cmdList, const std::string& pdfName, const std::string& dsName, const std::string& snapshot="")
{
  if (snapshot!="") { ws.loadSnapshot(snapshot.c_str()); }
  RooLinkedList fitConf; auto cmdL = cmdList; for (auto& cmd : cmdL) { fitConf.Add(dynamic_cast<TObject*>(&cmd)); };
  const auto& tmp = ws.pdf(pdfName.c_str())->fitTo(*ws.data(dsName.c_str()), fitConf);
  fitResult.reset(tmp);
  bool fitPass = true; for (uint iSt = 0; iSt < fitResult->numStatusHistory(); iSt++) { if (fitResult->statusCodeHistory(iSt)!=0) { fitPass = false; break; } }
  fitResult->Print("v");
  return fitPass;
};


int importDataset(RooWorkspace& myws, GlobalInfo& info, const RooWorkspaceMap_t& inputWS, const std::string& chg)
{
  // Check info container
  if (!contain(info.StrS, "dsList")) { std::cout << "[ERROR] DSList was not found while importing dataset!" << std::endl; return -1; }
  for (const auto& inWS : inputWS) { if (inWS.second.allData().empty()) { std::cout << "[ERROR] Input workspace " << inWS.first << " is empty!" << std::endl; return -1; } }
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
  //
  // Add the decay length input selection
  if (contain(info.Par, "Cut") && info.Par.at("Cut")!="") {
    std::string cutSel = "";
    const auto& cutLbl = info.Par.at("Cut");
    const auto& PD = info.Par.at("PD");
    const bool isMuonTrigger = (PD.rfind("MUON")!=std::string::npos);
    std::string threshold = "0.90";
    if (cutLbl.rfind("[")!=std::string::npos) {
      const auto& tmp = cutLbl.substr(cutLbl.rfind("[")+1);
      threshold = tmp.substr(0, tmp.rfind("]"));
    }
    if (cutLbl.rfind("PromptDecay", 0)==0) {
      cutSel = ANA::CHARMONIA::decayLenCut("Cand_DLen", "Cand_Pt", "Cand_Rap", threshold, true);
      std::cout << "[INFO] Cutting on candidate decay length to select PROMPT decays" << std::endl;
    }
    else if (cutLbl.rfind("NonPromptDecay", 0)==0) {
      cutSel = ANA::CHARMONIA::decayLenCut("Cand_DLen", "Cand_Pt", "Cand_Rap", threshold, false);
      std::cout << "[INFO] Cutting on candidate decay length to select NON-PROMPT decays" << std::endl;
    }
    if (cutSel!="") {
      cutDS += " && "+cutSel;
      addString(myws, "cutSelExp", cutSel); // Save the cut expression for bookkeeping
      addString(myws, "cutSelStr", cutLbl); // Save the cut label for bookkeeping
    }
  }
  addString(myws, "cutDS", cutDS); // Save the cut expression for bookkeeping
  //
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
        if (!contain(inputWS, label) || !inputWS.at(label).data(dsExtName.c_str())){ 
          std::cout << "[ERROR] The dataset " <<  dsExtName << " was not found!" << std::endl; return -1;
        }
        const auto& cutDST = cutDS + (isSwapDS ? "&&(Cand_IsSwap==Cand_IsSwap::Yes)" : "");
        auto data = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(inputWS.at(label).data(dsExtName.c_str())->reduce(RooFit::Cut(cutDST.c_str()), RooFit::Name(dsName.c_str()), RooFit::Title(dsName.c_str()))));
        if (!data) { std::cout << "[ERROR] Dataset " <<  dsExtName << " failed to reduce!" << std::endl; return -1; }
        else if (data->sumEntries()==0){
          if (extLabel.rfind("MC_",0)==0 || chg=="SS") {
            std::cout << "[WARNING] No events from dataset " <<  dsExtName << " passed the kinematic cuts!" << std::endl;
          }
          else { std::cout << "[ERROR] No events from dataset " <<  dsExtName << " passed the kinematic cuts!" << std::endl; return -1; }
        }
        else {
          if (myws.import(*data)) { std::cout << "[ERROR] DataSet " << dsName << " was not imported!" << std::endl; return -1; }
	  if (!myws.data(dsName.c_str())) { std::cout << "[ERROR] Importing RooDataSet " <<  dsName << " failed!" << std::endl; return -1; }
	  copyWorkspace(myws, inputWS.at(label), "", false, false, true); // Copy the generic objects
	  storeDSRange(myws, data.get());
        }
        std::cout << "[INFO] " << data->numEntries() << " entries imported from local RooDataSet " << dsName << std::endl;
        // Set the range of each global parameter in the local roodataset
	if (myws.data(dsName.c_str())) { setDSParamaterRange(myws, dsName, info); }
      }
    }
  }
  // Check if user wants to use the Center-of-Mass frame
  for (auto& v : info.Var) {
    if (contain(info.Flag, "use"+v.first+"CM") && info.Flag.at("use"+v.first+"CM")) {
      myws.factory(Form("use%sCM[1.0]", v.first.c_str()));
      const auto& var = info.Var.at(v.first+"CM");
      myws.factory(Form("%sCM[%.10f,%.10f,%.10f]", v.first.c_str(), 0.0, var.at("Min"), var.at("Max")));
      myws.var((v.first+"CM").c_str())->setRange("DEFAULT", var.at("Default_Min"), var.at("Default_Max"));
    }
  }
  // Set the range of each global parameter in the local workspace
  StringSet_t varSet;
  for (const auto& var : info.Var) {
    const bool& isAbs = (var.first.find("Abs")!=std::string::npos && contain(info.StrS.at("cutPars"), var.first));
    auto varN = var.first; if (isAbs) { varN.erase(varN.find("Abs"), 3); }
    if (contain(varSet, varN)) continue;
    if (myws.var(varN.c_str()) && contain(var.second, "Min")) {
      myws.var(varN.c_str())->setMin(isAbs ? -var.second.at("Max") : var.second.at("Min"));
      myws.var(varN.c_str())->setMax(var.second.at("Max"));
      varSet.insert(varN);
    }
    else if (!myws.var(var.first.c_str()) && contain(var.second, "Val")) {
      myws.factory(Form("%s[%.10f]", var.first.c_str(), var.second.at("Val")));
    }
    if (contain(info.StrS.at("cutPars"), var.first) && !myws.var(var.first.c_str()) && contain(var.second, "Min")) {
      myws.factory(Form("%s[%.10f,%.10f,%.10f]", var.first.c_str(), 0.0, var.second.at("Min"), var.second.at("Max")));
      myws.var(var.first.c_str())->setRange("DEFAULT", var.second.at("Default_Min"), var.second.at("Default_Max"));
    }
  }
  // Print bin information
  std::string binInfo = "[INFO] Analyzing bin:";
  for (const auto& var : info.StrS.at("cutPars")) {
    // Ignore if it is the fitting variable
    if (contain(info.StrS.at("fitVariable"), var)) continue;
    // Check if is an absolute value variable
    // Print the binning
    if (myws.var(var.c_str())) {
      const auto& varMin = myws.var(var.c_str())->getMin();
      const auto& varMax = myws.var(var.c_str())->getMax();
      if (contain(info.Flag, "use"+var+"CM")  && info.Flag.at("use"+var+"CM")) {
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


#endif // #ifndef rooDataUtils_h
