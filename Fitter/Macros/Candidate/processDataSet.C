// -*- C++ -*-
//
// Package:    Fitter
//
/*
 Description: Process RooDataSet.
 Implementation:
 This program process RooDataSets.
 */
// Original Author:  Andre Stahl,
//         Created:  May 23 19:08 CET 2019
//
//
#ifndef Candidate_processDataSet_C
#define Candidate_processDataSet_C


#include "RooWorkspace.h"
#include "RooArgSet.h"

#include <iostream>
#include <string>

#include "../Utilities/rooDataUtils.h"
#include "../Utilities/initClasses.h"


bool skimDataSet        ( RooWorkspaceMap_t& workspaces , const GlobalInfo& info );
bool combineDataSet     ( RooWorkspaceMap_t& workspaces , GlobalInfo& info , const std::string& run="8Y16" );
bool processDecayLength ( RooWorkspaceMap_t& workspaces , const GlobalInfo& info );


bool processDataSet(RooWorkspaceMap_t& workspaces, GlobalInfo& info)
{
  info.StrV["dsType"].push_back("RAW");
  //
  // -------------------------------------------------------------------------------
  // STEP 1: SKIM THE DATASETS
  //
  if (!skimDataSet(workspaces, info)) { return false; }
  //
  // -------------------------------------------------------------------------------
  // STEP 2: COMBINE THE DATASETS FOR pPb
  //
  if (!combineDataSet(workspaces, info)) { return false; }
  //
  // -------------------------------------------------------------------------------
  // STEP 3: PROCESS THE CANDIDATE DECAY LENGTH INFORMATION
  //
  if (!processDecayLength(workspaces, info)) { return false; }
  //
  return true;
};


bool skimDataSet(RooWorkspaceMap_t& workspaces, const GlobalInfo& info)
{
  //
  std::cout << "[INFO] Proceed to skim the input datasets" << std::endl;
  //
  // Define the objecst needed for the mass range
  StringSet_t objS;
  for (const auto& v : info.StrS.at("incVarName")) {
    for (const auto& o : info.StrS.at("allObject_"+v)) {
      if (o.rfind("Bkg",0)==0) continue;
      auto obj = o; if (!contain(ANA::MASS, obj)) { for (const auto& m : ANA::MASS) { if (obj.rfind(m.first,0)==0) { obj = m.first; break; } } }
      objS.insert(obj);
    }
  }
  //
  // Variables to delete
  const std::vector<std::string> delVarNames = {"Cand_Qual", "Cand_Trig", "Cand_VtxP", "Dau1_Eta", "Dau1_Pt", "Dau2_Eta", "Dau2_Pt", "Cand_DLen2D", "Cand_DLenErr2D", "Cand_DLenGen2D", "Event_Sel"};
  //
  // Loop over the RooDataSets
  for (auto& ws : workspaces) {
    //
    // Determine the collision system of the sample
    std::string evtCol = "";
    if (ws.first.rfind("_")!=std::string::npos) { evtCol = ws.first.substr(ws.first.rfind("_")+1); }
    if (evtCol=="") { std::cout << "[ERROR] Could not determine the collision system in the sample" << std::endl; return false; }
    //
    // Determine the cut string
    const std::string& PD = (ws.second.obj("PD") ? dynamic_cast<RooStringVar*>(ws.second.obj("PD"))->getVal() : "");
    auto cutStr = ANA::analysisSelection("Dau1_Pt", "Dau1_Eta", "Dau2_Pt", "Dau2_Eta", "Cand_Mass", "Cand_VtxP", "Cand_Trig", "Cand_Qual", PD, evtCol, objS, true);
    if (ws.second.var("Event_Sel")) { cutStr += " && (int(Event_Sel) & 1)"; }
    //
    // Skim the datasets
    RooWorkspace tmpWS; copyWorkspace(tmpWS, ws.second, "", false);
    const auto& listData = ws.second.allData();
    for (const auto& ds : listData) {
      // Define the skimmed dataset variables
      auto vars = *ds->get();
      for (const auto& n : delVarNames) { const auto& var = vars.find(n.c_str()); if (var) { vars.remove(*var); } }
      // Reduce the dataset
      std::cout << "[INFO] Applying skim selection: " << cutStr << std::endl;
      auto data = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(ds->reduce(RooFit::Cut(cutStr.c_str()), RooFit::SelectVars(vars), RooFit::Name(ds->GetName()), RooFit::Title(ds->GetTitle()))));
      if (!data) { std::cout << "[ERROR] Skimmed dataset " << ds->GetName() << " is NULL!" << std::endl; return false; }
      else if (data->sumEntries()==0){ std::cout << "[ERROR] Skimmed dataset " <<  ds->GetName() << " is empty!" << std::endl; return false; }
      std::cout << "[INFO] Reduced dataset " << ds->GetName() << " from " << ds->sumEntries() << " to " << data->sumEntries() << " entries" << std::endl;
      addString(tmpWS, "skimDS", cutStr); // Save the cut expression for bookkeeping
      if (tmpWS.import(*data)) { std::cout << "[ERROR] skimDataSet: Failed to import " << data->GetName() << std::endl; return false; }
    }
    clearWorkspace(ws.second);
    copyWorkspace(ws.second, tmpWS);
  }
  //
  return true;
};


bool invertEtaAndFill(RooDataSet& dsPA, const RooDataSet& dsPbp, const double& lumiW=1.0)
{
  // Find the rapidity related variables
  StringSet_t varNames;
  const auto& iniVars = dsPbp.get();
  auto parIt = std::unique_ptr<TIterator>(iniVars->createIterator());
  for (auto itp = parIt->Next(); itp!=NULL; itp = parIt->Next()) {
    const std::string& vName = itp->GetName();
    if (vName.find("_Rap")!=std::string::npos || vName.find("_Eta")!=std::string::npos) {
      varNames.insert(vName);
    }
  }
  // Fill the PA dataset
  for(int i = 0; i < dsPbp.numEntries(); i++){
    const auto& dsSet = *dsPbp.get(i);
    // Invert the rapidities
    for (const auto& vName : varNames) {
      const auto& var = dynamic_cast<RooRealVar*>(dsSet.find(vName.c_str()));
      var->setVal(-1.0*var->getVal());
    }
    // Add the new event
    const auto& weight = lumiW*dsPbp.weight();
    if (weight <= 0.0) { std::cout << "[ERROR] invertEtaAndFill: Pbp weight is negative ( " << weight << " )" << std::endl; return false; }
    dsPA.add(dsSet, weight);
  }
  //
  return true;
};


bool combineDataSet(RooWorkspaceMap_t& workspaces, GlobalInfo& info, const std::string& run)
{
  if (!contain(info.Flag, "fitPA"+run) || info.Flag.at("fitPA"+run)==false) { return true; }
  //
  std::cout << "[INFO] Proceed to combine the pPb" << run << " RooDatasets" << std::endl;
  //
  // Determine the sample tags (sample name without the collision tag)
  StringSet_t sampleTags;
  for (const auto& ws : workspaces) { sampleTags.insert(ws.first.substr(0, ws.first.rfind("_"))); }
  //
  // Loop over the sample tags
  for (const auto& smTag : sampleTags) {
    const auto& sample_pPb = smTag+"_pPb"+run;
    const auto& sample_Pbp = smTag+"_Pbp"+run;
    const auto& sample_PA  = smTag+"_PA"+run;
    //
    // Check input datasets
    if (!contain(workspaces, sample_pPb)) { std::cout << "[ERROR] RooWorkspace for sample " << sample_pPb << " does not exist!" << std::endl; return false; }
    if (!contain(workspaces, sample_Pbp)) { std::cout << "[ERROR] RooWorkspace for sample " << sample_Pbp << " does not exist!" << std::endl; return false; }
    const auto& ws_pPb = workspaces.at(sample_pPb);
    const auto& ws_Pbp = workspaces.at(sample_Pbp);
    auto& ws_PA = workspaces[sample_PA];
    copyWorkspace(ws_PA, ws_pPb, "");
    const std::string& PD = (ws_PA.obj("PD") ? dynamic_cast<RooStringVar*>(ws_PA.obj("PD"))->getVal() : "");
    if (PD=="") { std::cout << "[ERROR] PD was not defined!" << std::endl; return false; }
    //
    // Determine the dataset tags
    StringSet_t dsTags;
    const auto& listData = ws_pPb.allData();
    for (const auto& ds : listData) {
      const std::string dsName = ds->GetName();
      dsTags.insert(dsName.substr(0, dsName.rfind(sample_pPb)));
    }
    //
    // Loop over the dataset tags
    for (const auto& dsTag : dsTags) {
      const auto& dsTag_pPb = dsTag+sample_pPb;
      const auto& dsTag_Pbp = dsTag+sample_Pbp;
      const auto& dsTag_PA  = dsTag+sample_PA;
      //
      // Extract the datasets
      const auto& ds_pPb = dynamic_cast<RooDataSet*>(ws_pPb.data(dsTag_pPb.c_str()));
      const auto& ds_Pbp = dynamic_cast<RooDataSet*>(ws_Pbp.data(dsTag_Pbp.c_str()));
      if (!ds_pPb) { std::cout << "[ERROR] RooDataSet " << dsTag_pPb << " does not exist!" << std::endl; return false; }
      if (!ds_Pbp) { std::cout << "[ERROR] RooDataSet " << dsTag_Pbp << " does not exist!" << std::endl; return false; }
      //
      std::cout << "[INFO] Creating RooDataSet " << dsTag_PA << std::endl;
      // Copy the pPb datasets
      std::unique_ptr<RooDataSet> ds_PA;
      if (sample_pPb.find("DATA",0)==0) { ds_PA.reset(dynamic_cast<RooDataSet*>(ds_pPb->Clone(dsTag_PA.c_str()))); }
      else if (sample_pPb.find("MC",0)==0) {
	ds_PA.reset(dynamic_cast<RooDataSet*>(ds_pPb->emptyClone(dsTag_PA.c_str())));
	const auto& lumiW = pPb::R8TeV::Y2016::LumiWeightFromPD(PD, "pPb8Y16", sample_pPb);
	for(int i = 0; i < ds_pPb->numEntries(); i++){
	  const auto& dsSet = *ds_pPb->get(i);
	  const auto& weight = lumiW*ds_pPb->weight();
	  if (weight <= 0.0) { std::cout << "[ERROR] combineDataSet: pPb weight is negative ( " << weight << " )" << std::endl; return false; }
	  ds_PA->add(dsSet, weight);
	}
      }
      if (!ds_PA || ds_PA->sumEntries()==0) { std::cout << "[ERROR] RooDataSet " << dsTag_PA << " was not created!" << std::endl; return false; }
      // Invert the rapidity of Pbp dataset and fill the PA dataset
      if (sample_Pbp.find("DATA",0)==0) { if (!invertEtaAndFill(*ds_PA, *ds_Pbp)) { return false; } }
      else if (sample_Pbp.find("MC",0)==0) {
	const auto& lumiW = pPb::R8TeV::Y2016::LumiWeightFromPD(PD, "Pbp8Y16", sample_Pbp);
	if (!invertEtaAndFill(*ds_PA, *ds_Pbp, lumiW)) { return false; }
      }
      // Check the consistency of the combined dataset
      if (ds_PA->numEntries()!=(ds_pPb->numEntries()+ds_Pbp->numEntries())) { std::cout << "[ERROR] Number of entries for the combined " << dsTag_PA << " dataset is inconsistent!" << std::endl; return false; }
      if (std::abs(ds_PA->sumEntries()-(ds_pPb->sumEntries()+ds_Pbp->sumEntries()))>0.05*(ds_pPb->sumEntries()+ds_Pbp->sumEntries())) {
        std::cout << "[ERROR] Number of weighted entries for the combined " << dsTag_PA << " dataset ( " << ds_PA->sumEntries() << " , " << ds_pPb->sumEntries() << " , " << ds_Pbp->sumEntries() << " ) is inconsistent!" << std::endl; return false;
      }
      // Import the new PA dataset
      if (ws_PA.import(*ds_PA)) { std::cout << "[ERROR] combineDataSet: Failed to import " << ds_PA->GetName() << std::endl; return false; }
      else { std::cout << "[INFO] RooDataSet " << dsTag_PA << " was created with " << ds_PA->numEntries() << " entries!" << std::endl; }
      if (contain(info.StrS.at("DSTAG"), sample_pPb)) { info.StrS.at("DSTAG").insert(sample_PA); }
    }
    // Delete old datasets if not used
    if (!info.Flag.at("fitpPb"+run)) { clearWorkspace(workspaces.at(sample_pPb)); workspaces.erase(sample_pPb); info.StrS.at("DSTAG").erase(sample_pPb); }
    if (!info.Flag.at("fitPbp"+run)) { clearWorkspace(workspaces.at(sample_Pbp)); workspaces.erase(sample_Pbp); info.StrS.at("DSTAG").erase(sample_Pbp); }
  }
  //
  return true;
};


bool processDecayLength(RooWorkspaceMap_t& workspaces, const GlobalInfo& info)
{
  //
  // Check if user wants to fit the decay length information
  if (!info.Flag.at("fitCand_DLen") && !info.Flag.at("fitCand_DLenRes") && !info.Flag.at("fitCand_DLenErr") && !info.Flag.at("fitCand_DLenGen")) return true;
  //
  // Find the mass reference
  std::string ref = "";
  double minM = 9999999.;
  for (const auto& v : info.StrS.at("incVarName")) {
    for (const auto& o : info.StrS.at("allObject_"+v)) {
      if (o.rfind("Bkg",0)==0 || o=="DLenRes") continue;
      auto obj = o; if (!contain(ANA::MASS, obj)) { for (const auto& m : ANA::MASS) { if (obj.rfind(m.first,0)==0) { obj = m.first; break; } } }
      if (contain(ANA::MASS, obj)) {
	if (minM > ANA::MASS.at(obj).at("Loose_Min")) { ref = obj; minM = ANA::MASS.at(obj).at("Loose_Min"); }
      }
    }
  }
  // Keep mass reference of excited states to mass of 1S state in decay length
  if (ref=="Psi2S") { ref = "JPsi"; } else if (ref=="Ups2S") { ref = "Ups1S"; } else if (ref=="Ups3S") { ref = "Ups1S"; }
  //
  // Return if reference is JPsi and we dont want to produce the decay length resolution
  if (ref=="JPsi" && !info.Flag.at("fitCand_DLen") && !info.Flag.at("fitCand_DLenRes") && !info.Flag.at("fitCand_DLenErr")) return true;
  //
  std::cout << "[INFO] Proceed to process the decay length information using " << ref << " mass" << std::endl;
  //
  const auto& massRef   = ANA::MASS.at(ref).at("Val");
  const auto& massJPsi  = ANA::MASS.at("JPsi").at("Val");
  const auto& massRatio = massRef/massJPsi;
  //
  // Loop over the RooWorkspaces
  for (auto& ws : workspaces) {
    //
    // Loop over the datasets
    const auto& listData = ws.second.allData();
    for (const auto& ds : listData) {
      //
      // Define dataset variables
      auto vars = *ds->get();
      const auto& candDLen = dynamic_cast<RooRealVar*>(vars.find("Cand_DLen"));
      const auto& candDLenErr = dynamic_cast<RooRealVar*>(vars.find("Cand_DLenErr"));
      const auto& candDLenGen = dynamic_cast<RooRealVar*>(vars.find("Cand_DLenGen"));
      if (!candDLen || !candDLenErr) continue;
      auto candDLenRes = RooRealVar("Cand_DLenRes", "Candidate c#tau resolution", -100000.0, 100000.0, "");
      vars.add(candDLenRes);
      //
      // Create temporary dataset
      auto tmpDS = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(ds->emptyClone(0, 0, &vars, 0)));
      //
      // Loop over the entries//
      for (int i = 0; i < ds->numEntries(); i++) {
	ds->get(i);
	//
	// Process the decay length information
	candDLen->setVal(candDLen->getVal()*massRatio);
	candDLenErr->setVal(candDLenErr->getVal()*massRatio);
	if (candDLenGen) { candDLenGen->setVal(candDLenGen->getVal()*massRatio); }
	const auto& dLenDiff = candDLen->getVal() - (candDLenGen ? candDLenGen->getVal() : 0.0);
	candDLenRes.setVal(dLenDiff/candDLenErr->getVal());
	//
	// Fill temporary dataset
	tmpDS->addFast(vars, ds->weight());
      }
      // Import to RooWorkspace
      ws.second.RecursiveRemove(ds); if(ds) delete ds;
      if (ws.second.import(*tmpDS)) { std::cout << "[ERROR] processDecayLength: Failed to import " << tmpDS->GetName() << std::endl; return false; }
      else { std::cout << "[INFO] Processed " << tmpDS->numEntries() << " entries in RooDataSet " << tmpDS->GetName() << std::endl; }
      // Store dataset range
      storeDSRange(ws.second, tmpDS.get(), "DSFullWindow");
    }
    addString(ws.second, "dLenMassRef", ref); // Save the reference for bookkeeping
  }
  return true;
};


#endif
