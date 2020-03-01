#ifndef storeWS2ResultsTree_C
#define storeWS2ResultsTree_C

// Auxiliary Headers
#include "../Utilities/dataUtils.h"
#include "../Utilities/RunInfo/eventUtils.h"
// ROOT headers
#include "TFile.h"
#include "TTree.h"
// RooFit headers
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooRealVar.h"
#include "RooAddPdf.h"
#include "RooFormulaVar.h"
#include "RooFitResult.h"
#include "RooStringVar.h"
#include "RooCategory.h"
// c++ headers
#include <iostream>
#include <string>


double getErrorFromWS     ( const RooRealVar& var , const std::string& type="" );
void   iniResultsTreeInfo ( GlobalInfo& info      , const RooWorkspace& ws     );
void   setBranches        ( TTree& tree           , GlobalInfo& info           );
double getLumiFromPD      ( const std::string& PD , const std::string& col     );


bool storeWS2ResultsTree(
			 const std::string& workDirName = "Test",
			 const std::string& trgTag      = "DIMUON",
			 const std::string& colTag      = "PA8Y16",
			 const std::string& objTag      = "JPsi",
			 const std::string& dataTag     = "DATA",
			 const std::string& varTag      = "Cand_Mass"
			 )
{
  //
  // --------------------------------------------------------------------------------- //
  //
  // Define the output file info
  const std::string& CWD = getcwd(NULL, 0);
  const std::string& outputFileName = "tree_allvars.root";
  const auto& dsTag = (dataTag+"_"+(dataTag=="DATA" ? "" : (objTag+(trgTag.rfind("Cat",0)==0?"":"_")))+trgTag);
  auto vTag = varTag; if (vTag.find("_")!=std::string::npos) { vTag.erase(vTag.find("_"), 1); }
  const auto& outputDirPath  = (CWD+"/Tree/"+workDirName+"/"+vTag+"/"+dsTag+"/"+objTag+"/"+colTag);
  const auto& outputFilePath = (outputDirPath+"/"+outputFileName);
  //
  // --------------------------------------------------------------------------------- //
  //
  // Define the tree info container
  GlobalInfo info;
  //
  // Initialize the tree
  TTree tree("fitResults", "Fit Results");
  //
  // --------------------------------------------------------------------------------- //
  //
  // Get the list of input files
  //
  StringVector_t inputFileNames;
  std::string preCWD = CWD; preCWD.erase(preCWD.find_last_of("/"), 10);
  const auto& inputDirPath = (preCWD+"/Fitter/Output/"+workDirName+"/"+vTag+"/"+dsTag+"/"+objTag+"/"+colTag+"/result");
  if (!existDir(inputDirPath)) { std::cout << "[ERROR] Workspace directory " << inputDirPath << " was not found!" << std::endl; return false; }
  if (!fileList(inputFileNames, inputDirPath)) { return false; };
  //
  // --------------------------------------------------------------------------------- //
  //
  // Loop over the input files
  //
  bool isFirstFile = true;
  for (const auto& inputFileName : inputFileNames) {
    //
    //std::cout << "Processing file: " << inputFileName << std::endl;
    //
    // Open input file
    const auto& inputFilePath = (inputDirPath+"/"+inputFileName);
    TFile inputFile(inputFilePath.c_str(), "READ");
    //
    if (inputFile.IsOpen()==false || inputFile.IsZombie()==true) {
      std::cout << "[ERROR] The input file " << inputFilePath << " could not be created!" << std::endl; return false;
    }
    inputFile.cd();
    //
    // Extract the Workspace
    const auto& ws = dynamic_cast<RooWorkspace*>(inputFile.Get("workspace"));
    if (!ws) { std::cout << "[ERROR] Workspace not found in " << inputFilePath << std::endl; inputFile.Close(); return false; }
    //
    // Initialize the tree
    if (isFirstFile) {
      iniResultsTreeInfo(info, *ws);
      setBranches(tree, info);
      isFirstFile = false;
    }
    //
    // Fill the string information
    for (auto& o : info.Par) {
      const auto& obj = dynamic_cast<RooStringVar*>(ws->obj(o.first.c_str()));
      if (obj) { o.second = obj->getVal(); }
    }
    // Check the information
    const StringSet_t checkStr = { "DSTAG" , "channel", "fitSystem", "PD", "dsName", "pdfName" };
    for (const auto& s : checkStr) { if (!contain(info.Par, s)) { std::cout << "[ERROR] " << s << " was not found in workspace!" << std::endl; inputFile.Close(); return false; } }
    if (dsTag.rfind("DATA_",0)==0 && info.Par.at("DSTAG").find(dsTag)==std::string::npos) {
      std::cout << "[ERROR] Workspace DSTAG " << info.Par.at("DSTAG") << " is not consistent with input dsTag " << dsTag << std::endl; inputFile.Close(); return false;
    }
    if (info.Par.at("fitSystem")!=colTag) { std::cout << "[ERROR] Workspace COL " << info.Par.at("fitSystem") << " is not consistent with input colTag " << colTag << std::endl; inputFile.Close(); return false; }
    if (info.Par.at("fitSystem")!=colTag) { std::cout << "[ERROR] Workspace COL " << info.Par.at("fitSystem") << " is not consistent with input colTag " << colTag << std::endl; inputFile.Close(); return false; }
    if (info.Par.at("channel")!="ToMuMu") { std::cout << "[ERROR] Only dimuon channel is currently supported in result macros!" << std::endl; return false; }
    //
    // Fill the flag information
    for (auto& f : info.Flag) {
      const auto& var = dynamic_cast<RooRealVar*>(ws->var(f.first.c_str()));
      if (var) { f.second = (var->getVal()==1.0); }
    }
    //
    // Get the dataset
    auto ds = dynamic_cast<RooAbsData*>(ws->data(("CutAndCount_"+info.Par.at("dsName")).c_str()));
    if (!ds) { ds = dynamic_cast<RooAbsData*>(ws->data(info.Par.at("dsName").c_str())); }
    //
    // Fill the observables
    for (auto& v : info.Var) {
      if (v.first.find("OBS_")==std::string::npos) continue;
      auto obsName = v.first; obsName.erase(obsName.find("OBS_"), 4);
      //
      const auto& var  = dynamic_cast<RooRealVar*>(ws->var(obsName.c_str()));
      const auto& inDS  = (ds ? (ds->get()->find(obsName.c_str())!=NULL) : false);
      const auto& mean = ((var && inDS) ? dynamic_cast<RooRealVar*>(ds->meanVar(*var)) : NULL);
      const auto& rms  = ((var && inDS) ? dynamic_cast<RooRealVar*>(ds->rmsVar(*var))  : NULL);
      //
      if (contain(v.second, "Min")) { v.second.at("Min") = var  ? var->getMin()   : -99.0; }
      if (contain(v.second, "Max")) { v.second.at("Max") = var  ? var->getMax()   : -99.0; }
      if (contain(v.second, "Val")) { v.second.at("Val") = mean ? mean->getVal()  : -99.0; }
      if (contain(v.second, "Err")) { v.second.at("Err") = rms  ? rms->getVal() : -99.0; }
      if (contain(v.second, "DefaultMin")) { v.second.at("DefaultMin") = var ? var->getMin("DEFAULT") : -99.0; }
      if (contain(v.second, "DefaultMax")) { v.second.at("DefaultMax") = var ? var->getMax("DEFAULT") : -99.0; }
    }
    //
    // Get the snapshots
    const auto& parIni = ws->getSnapshot("initialParameters");
    //
    // Get the fit result
    RooFitResult* fitResult=0;
    const auto& listObj = ws->allGenericObjects();
    for (const auto& ito : listObj) { const auto& it = dynamic_cast<RooFitResult*>(ito); if (it) { fitResult = it; break; } }
    //
    // Fill the parameters
    for (auto& p : info.Var) {
      if (p.first.find("PAR_")==std::string::npos) continue;
      auto parName = p.first; parName.erase(parName.find("PAR_"), 4);
      //
      if (ws->var(parName.c_str())) {
        const auto& par = dynamic_cast<RooRealVar*>(ws->var(parName.c_str()));
        const auto& par_Ini = (parIni ? dynamic_cast<RooRealVar*>(parIni->find(parName.c_str())) : NULL);
        if (contain(p.second, "Min")  ) { p.second.at("Min")   = par ? par->getMin()    : -99.0; }
        if (contain(p.second, "Max")  ) { p.second.at("Max")   = par ? par->getMax()    : -99.0; }
        if (contain(p.second, "Val")  ) { p.second.at("Val")   = par ? par->getVal()    : -99.0; }
        if (contain(p.second, "ErrLo")) { p.second.at("Err")   = par ? getErrorFromWS(*par) : -99.0; }
        if (contain(p.second, "ErrLo")) { p.second.at("ErrLo") = par ? getErrorFromWS(*par, "Lo") : -99.0; }
        if (contain(p.second, "ErrHi")) { p.second.at("ErrHi") = par ? getErrorFromWS(*par, "Hi") : -99.0; }
        if (contain(p.second, "iniVal")) { p.second.at("iniVal") = par_Ini ? par_Ini->getVal()   : -99.0; }
        if (contain(p.second, "iniErr")) { p.second.at("iniErr") = par_Ini ? par_Ini->getError() : -99.0; }
      }
      else if (ws->function(parName.c_str())) {
        const auto& par = dynamic_cast<RooFormulaVar*>(ws->function(parName.c_str()));
        const double error = ((par && fitResult) ? par->getPropagatedError(*fitResult) : -99.0);
        if (contain(p.second, "Min")  ) { p.second.at("Min")   = -99.0; }
        if (contain(p.second, "Max")  ) { p.second.at("Max")   = -99.0; }
        if (contain(p.second, "Val")  ) { p.second.at("Val")   = par ? par->getVal() : -99.0; }
        if (contain(p.second, "ErrLo")) { p.second.at("Err")   = error; }
        if (contain(p.second, "ErrLo")) { p.second.at("ErrLo") = error; }
        if (contain(p.second, "ErrHi")) { p.second.at("ErrHi") = error; }
        if (contain(p.second, "iniVal")) { p.second.at("iniVal") = (par && parIni) ? par->getValV(parIni) : -99.0; }
        if (contain(p.second, "iniErr")) { p.second.at("iniErr") = -99.0; }
      }
    }
    //
    // Fill the remaining Variable Information
    for (auto& v : info.Var) {
      if (v.first=="Luminosity") {
	v.second.at("Val") = getLumiFromPD(info.Par.at("PD"), info.Par.at("fitSystem"));
      }
      else if (v.first=="N_DS_Entries") {
	const auto& varName = "numEntries_"+info.Par.at("dsName");
        v.second.at("Val") = (ds ? ds->sumEntries() : (contain(info.Var, varName) ? info.Var.at(varName).at("Val") : -99.0));
      }
      else if (v.first=="N_FIT_Entries") {
        const auto& pdf = dynamic_cast<RooAddPdf*>(ws->pdf(info.Par.at("pdfName").c_str()));
        v.second.at("Val") = (pdf ? pdf->expectedEvents(pdf->coefList()) : -99.0);
      }
      else if (v.first=="TEST_FIT") {
        v.second.at("Val")  = (contain(info.Var, "PAR_pvalue_BCChi2_"+varTag  ) ? info.Var.at("PAR_pvalue_BCChi2_"+varTag).at("Val")  : -99.0);
        v.second.at("Chi2") = (contain(info.Var, "PAR_testStat_BCChi2_"+varTag) ? info.Var.at("PAR_testStat_BCChi2_"+varTag).at("Val") : -99.0);
        v.second.at("NDoF") = (contain(info.Var, "PAR_ndofc_BCChi2_"+varTag   ) ? info.Var.at("PAR_ndofc_BCChi2_"+varTag).at("Val")   : -99.0);
      }
    }
    //
    // Close the input file
    inputFile.Close();
    //
    // Fill the tree
    tree.Fill();
  } // loop on the files
  //
  // Create the output file
  gSystem->mkdir(outputDirPath.c_str(), kTRUE);
  TFile outputFile(outputFilePath.c_str(), "RECREATE");
  if (outputFile.IsOpen()==false || outputFile.IsZombie()==true) {
    std::cout << "[ERROR] The output file " << outputFilePath << " could not be created!" << std::endl; return false;
  }
  outputFile.cd();
  //
  // Write the output tree
  tree.Write();
  //
  // Write the output file
  outputFile.Write();
  // Close the output file
  outputFile.Close();
  //
  // return
  return true;
};


double getErrorFromWS(const RooRealVar& var, const std::string& type)
{
  const bool hasAsymErrors = (var.getErrorLo()!=0.0 || var.getErrorHi()!=0.0);
  if (type=="Hi" && hasAsymErrors) { return std::abs(var.getErrorHi()); }
  if (type=="Lo" && hasAsymErrors) { return std::abs(var.getErrorLo()); }
  if (hasAsymErrors) { return 0.5*(std::abs(var.getErrorLo()) + std::abs(var.getErrorHi())); }
  return var.getError();
};


void iniResultsTreeInfo(GlobalInfo& info , const RooWorkspace& ws)
{
  //
  // Find the dataset name
  const auto& dsStr = dynamic_cast<RooStringVar*>(ws.obj("dsName"));
  const std::string& dsName = (dsStr ? dsStr->getVal() : "");
  // Find the observable names
  StringSet_t obsNames, catNames;
  const auto& listObs = const_cast<RooWorkspace*>(&ws)->set(("SET_"+dsName).c_str());
  if (listObs) {
    auto obsIt = std::unique_ptr<TIterator>(listObs->createIterator());
    for (auto it = obsIt->Next(); it!=NULL; it = obsIt->Next()) {
      if (dynamic_cast<RooRealVar*>(it)) { obsNames.insert(it->GetName()); }
      else if (dynamic_cast<RooCategory*>(it)) { catNames.insert(it->GetName()); }
    }
  }
  else { obsNames = StringSet_t({"Cand_Mass", "Cand_APhi", "Cand_Pt", "Cand_Rap", "Cand_Len", "Centrality", "NTrack"}); } // BUG FIX
  // Add CM and abs variables
  for (const auto& obs : obsNames) {
    bool found = false;
    const auto varAbs = ("Cand_Abs")+obs.substr(obs.rfind("_")+1);
    if (ws.var((obs+"CM").c_str())) {  obsNames.insert(obs+"CM"); found = true; }
    else if (ws.var(varAbs.c_str())) { obsNames.insert(varAbs); obsNames.insert(obs+"CM"); found = true; }
    if (found) { obsNames.erase(obsNames.find(obs)); }
  }
  // Initialize the observables
  const StringSet_t obsType = { "Min" , "Max" , "DefaultMin" , "DefaultMax" , "Val" , "Err" };
  for (const auto& o : obsNames) { for (const auto& t : obsType) { info.Var["OBS_"+o][t] = -99.0; } }
  //
  // Initialize the categories
  const StringSet_t catType = { "Val" };
  for (const auto& o : catNames) { for (const auto& t : catType) { info.Var["CAT_"+o][t] = -99.0; } }
  //
  // Define the parameter names
  StringSet_t parNames;
  auto parIt = std::unique_ptr<TIterator>(ws.componentIterator());
  if (parIt) {
    for (auto it = parIt->Next(); it!=NULL; it = parIt->Next()) {
      if (contain(obsNames, it->GetName())) continue;
      if (dynamic_cast<RooRealVar*>(it) || dynamic_cast<RooFormulaVar*>(it)) { parNames.insert(it->GetName()); }
    }
  }
  const StringSet_t parType = {"Min" , "Max" , "Val" , "Err" , "ErrLo" , "ErrHi" , "iniVal" , "iniErr"};
  // Initialize the parameters
  for (const auto& p : parNames) {
    const auto& par = dynamic_cast<RooRealVar*>(ws.arg(p.c_str()));
    const bool isConst = (par && ((par->getMin()==par->getMax()) || par->isConstant()));
    for (const auto& t : parType) { if (!isConst || t=="Val") { info.Var["PAR_"+p][t] = -99.0; } }
  }
  //
  // Define the names of the remaining variables
  info.Var["Luminosity"]["Val"] = -99.0;
  info.Var["N_DS_Entries"]["Val"] = -99.0;
  info.Var["N_FIT_Entries"]["Val"] = -99.0;
  info.Var["TEST_FIT"]["Val"] = -99.0;
  info.Var["TEST_FIT"]["Chi2"] = -99.0;
  info.Var["TEST_FIT"]["NDoF"] = -99.0;
  //
  // Initialize the information strings
  const auto& listObj = ws.allGenericObjects();
  for (const auto& its : listObj) { const auto& it = dynamic_cast<RooStringVar*>(its); if (it) { info.Par[it->GetTitle()] = ""; } }
  //
  // Boolean flags
  info.Flag["useCand_RapCM"] = false;
  info.Flag["useCand_AbsRap"] = false;
  //

};


void setBranches(TTree& tree , GlobalInfo& info)
{
  //
  // Create Tree branches for variables
  for (auto& v : info.Var) {
    for (auto& p : v.second) {
      tree.Branch(Form("%s_%s", v.first.c_str(), p.first.c_str()) , &(p.second) , Form("%s_%s/D", v.first.c_str(), p.first.c_str()));
    }
  }
  //
  // Create Tree branches for flags
  for (auto& f : info.Flag) {
    tree.Branch(Form("%s", f.first.c_str()) , &(f.second) , Form("%s/O", f.first.c_str()));
  }
  //
  // Create Tree branches for char arrays
  for (auto& c : info.Par) {
    tree.Branch(Form("%s", c.first.c_str()) , &(c.second));
  }
  //
};


double getLumiFromPD(const std::string& PD, const std::string& col)
{
  double lumi = -99.0;
  if      (col=="PbPb5Y18") { lumi = PbPb::R5TeV::Y2018::LumiFromPD(PD); }
  else if (col=="PP13Y18" ) { lumi = pp::R13TeV::Y2018::LumiFromPD(PD);  }
  else if (col=="PP5Y17"  ) { lumi = pp::R5TeV::Y2017::LumiFromPD(PD);   }
  else if (col.rfind("8Y16")!=std::string::npos) { lumi = pPb::R8TeV::Y2016::LumiFromPD(PD, col); }
  else if (col=="PbPb5Y15") { lumi =  PbPb::R5TeV::Y2015::LumiFromPD(PD); }
  return lumi;
};


#endif // #ifndef storeWS2ResultsTree_C
