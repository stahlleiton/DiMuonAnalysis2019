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


bool skimDataSet    ( RooWorkspaceMap_t& workspaces , const GlobalInfo& info );
bool combineDataSet ( RooWorkspaceMap_t& workspaces , GlobalInfo& info , const std::string& run="8Y16" );


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
  return true;
};


bool skimDataSet(RooWorkspaceMap_t& workspaces, const GlobalInfo& info)
{
  //
  std::cout << "[INFO] Proceed to skim the input datasets" << std::endl;
  //
  // Define the invariant mass range
  double minM = 999999., maxM = -999999.;
  for (const auto& obj : info.StrS.at("incObject")) {
    if (contain(MASS, obj)) {
      if (minM > MASS.at(obj).at("Min")) { minM = MASS.at(obj).at("Min"); }
      if (maxM < MASS.at(obj).at("Max")) { maxM = MASS.at(obj).at("Max"); }
    }
  }
  const std::string massCut = Form("(Cand_Mass > %g && Cand_Mass < %g)", minM, maxM);
  //
  // Define the muon kinematic cuts
  std::string PsiAcceptance_2015_Hybrid = "((abs(Dau1_Eta)<1.2 && Dau1_Pt>=3.5) || (1.2<=abs(Dau1_Eta) && abs(Dau1_Eta)<2.1 && Dau1_Pt>=5.77-1.89*abs(Dau1_Eta)) || (2.1<=abs(Dau1_Eta) && abs(Dau1_Eta)<2.4 && Dau1_Pt>=1.8))";
  PsiAcceptance_2015_Hybrid += " && ((abs(Dau2_Eta)<1.2 && Dau2_Pt>=3.5) || (1.2<=abs(Dau2_Eta) && abs(Dau2_Eta)<2.1 && Dau2_Pt>=5.77-1.89*abs(Dau2_Eta)) || (2.1<=abs(Dau2_Eta) && abs(Dau2_Eta)<2.4 && Dau2_Pt>=1.8))";
  std::string PsiAcceptance_2015_Soft = "((abs(Dau1_Eta)<1.0 && Dau1_Pt>=3.3) || (1.0<=abs(Dau1_Eta) && abs(Dau1_Eta)<2.2 && Dau1_Pt>=2.9) || (2.2<=abs(Dau1_Eta) && abs(Dau1_Eta)<2.4 && Dau1_Pt>=0.8))";
  PsiAcceptance_2015_Soft += " && ((abs(Dau2_Eta)<1.0 && Dau2_Pt>=3.3) || (1.0<=abs(Dau2_Eta) && abs(Dau2_Eta)<2.2 && Dau2_Pt>=2.9) || (2.2<=abs(Dau2_Eta) && abs(Dau2_Eta)<2.4 && Dau2_Pt>=0.8))";
  //
  std::string PsiAcceptance_2018_Hybrid = "((abs(Dau1_Eta)<1.2 && Dau1_Pt>=3.5) || (1.2<=abs(Dau1_Eta) && abs(Dau1_Eta)<2.1 && Dau1_Pt>=5.47-1.89*abs(Dau1_Eta)) || (2.1<=abs(Dau1_Eta) && abs(Dau1_Eta)<2.4 && Dau1_Pt>=1.5))";
  PsiAcceptance_2018_Hybrid += " && ((abs(Dau2_Eta)<1.2 && Dau2_Pt>=3.5) || (1.2<=abs(Dau2_Eta) && abs(Dau2_Eta)<2.1 && Dau2_Pt>=5.47-1.89*abs(Dau2_Eta)) || (2.1<=abs(Dau2_Eta) && abs(Dau2_Eta)<2.4 && Dau2_Pt>=1.5))";
  std::string PsiAcceptance_2018_Loose = "((abs(Dau1_Eta)<0.3 && Dau1_Pt>=3.4) || (0.3<=abs(Dau1_Eta) && abs(Dau1_Eta)<1.1 && Dau1_Pt>=3.3) || (1.1<=abs(Dau1_Eta) && abs(Dau1_Eta)<1.4 && Dau1_Pt>=7.7-4.0*abs(Dau1_Eta)) || (1.4<=abs(Dau1_Eta) && abs(Dau1_Eta)<1.55 && Dau1_Pt>=2.1) || (1.55<=abs(Dau1_Eta) && abs(Dau1_Eta)<2.2 && Dau1_Pt>=4.25-1.39*abs(Dau1_Eta)) || (2.2<=abs(Dau1_Eta) && abs(Dau1_Eta)<2.4 && Dau1_Pt>=1.2))";
  PsiAcceptance_2018_Loose += " && ((abs(Dau2_Eta)<0.3 && Dau2_Pt>=3.4) || (0.3<=abs(Dau2_Eta) && abs(Dau2_Eta)<1.1 && Dau2_Pt>=3.3) || (1.1<=abs(Dau2_Eta) && abs(Dau2_Eta)<1.4 && Dau2_Pt>=7.7-4.0*abs(Dau2_Eta)) || (1.4<=abs(Dau2_Eta) && abs(Dau2_Eta)<1.55 && Dau2_Pt>=2.1) || (1.55<=abs(Dau2_Eta) && abs(Dau2_Eta)<2.2 && Dau2_Pt>=4.25-1.39*abs(Dau2_Eta)) || (2.2<=abs(Dau2_Eta) && abs(Dau2_Eta)<2.4 && Dau2_Pt>=1.2))";
  //
  const std::string PsiAcceptance = "((abs(Dau1_Eta)<2.4 && Dau1_Pt*cosh(Dau1_Eta)>=3.5) && (abs(Dau2_Eta)<2.4 && Dau2_Pt*cosh(Dau2_Eta)>=3.5))";
  const std::string UpsAcceptance = "((abs(Dau1_Eta)<2.4 && Dau1_Pt>=3.4) && (abs(Dau2_Eta)<2.4 && Dau2_Pt>=3.4))";
  const std::string ZAcceptance = "((abs(Dau1_Eta)<2.4 && Dau1_Pt>=15.0) && (abs(Dau2_Eta)<2.4 && Dau2_Pt>=15.0))";
  //
  // Define the muon quality cuts
  const std::string useSoftMuons = "(Cand_Qual & 1)";
  const std::string useHybridMuons = "(Cand_Qual & 2)";
  const std::string useTightMuons = "(Cand_Qual & 4)";
  //
  // Variables to delete
  const std::vector<std::string> delVarNames = {"Cand_Qual", "Cand_Trig", "Cand_VtxP", "Dau1_Eta", "Dau1_Pt", "Dau2_Eta", "Dau2_Pt"}; 
  //
  // Loop over the RooDataSets
  for (auto& ws : workspaces) {
    //
    // Determine the collision system of the sample
    std::string evtCol = "";
    if (ws.first.rfind("_")!=std::string::npos) { evtCol = ws.first.substr(ws.first.rfind("_")+1); }
    if (evtCol=="") { std::cout << "[ERROR] Could not determine the collision system in the sample" << std::endl; return false; }
    //
    // Define the trigger selection
    std::vector<uint> trigIdx;
    const auto& PD = info.Par.at("PD");
    if      (evtCol=="PP13Y18" ) { trigIdx = pp::R13TeV::Y2018::HLTBitsFromPD(PD); }
    else if (evtCol=="PP5Y17"  ) { trigIdx = pp::R5TeV::Y2017::HLTBitsFromPD(PD); }
    else if (evtCol=="PbPb5Y18") { trigIdx = PbPb::R5TeV::Y2018::HLTBitsFromPD(PD); }
    else if (evtCol=="PbPb5Y15") { trigIdx = PbPb::R5TeV::Y2015::HLTBitsFromPD(PD); }
    else if (evtCol.rfind("8Y16")!=std::string::npos) { trigIdx = pPb::R8TeV::Y2016::HLTBitsFromPD(PD); }
    if (trigIdx.empty()) { std::cout << "[ERROR] Could not determine the trigger index for the sample" << std::endl; return false; }
    std::string trigCut = "";
    for (const auto& idx : trigIdx) { trigCut += Form("(Cand_Trig & %.0f) ||", std::pow(2.0, idx)); }
    trigCut = trigCut.substr(0, trigCut.rfind(" ||")); if (trigIdx.size()>1) { trigCut = "("+trigCut+")"; }
    //
    // Determine the cut string
    std::string cutStr = massCut+" && "+trigCut;
    if (minM > 30.) { cutStr += " && "+ZAcceptance+" && "+useTightMuons; }
    else if (minM > 5.0) {
      const bool isSoft = (info.Par.at("PD")=="UPC" || evtCol.find("5Y1")==std::string::npos);
      cutStr += " && "+UpsAcceptance+" && "+(isSoft ? useSoftMuons : useHybridMuons);
    }
    else {
      const bool isSoft = (info.Par.at("PD")=="UPC" || evtCol.find("5Y1")==std::string::npos);
      const bool isAcc2018 = (evtCol.rfind("Y18")!=std::string::npos || evtCol.rfind("Y17")!=std::string::npos);
      cutStr += " && "+(isAcc2018 ? PsiAcceptance_2018_Hybrid : (isSoft ? PsiAcceptance : PsiAcceptance_2015_Hybrid));
      cutStr += " && "+(isSoft ? useSoftMuons : useHybridMuons);
    }
    //
    // Skim the datasets
    RooWorkspace tmpWS; copyWorkspace(tmpWS, ws.second, "", false);
    const auto& listData = ws.second.allData();
    for (const auto& ds : listData) {
      // Define the skimmed dataset variables
      auto vars = *ds->get();
      for (const auto& n : delVarNames) { const auto& var = vars.find(n.c_str()); vars.remove(*var); }
      // Reduce the dataset
      std::cout << "[INFO] Applying skim selection: " << cutStr << std::endl;
      auto data = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(ds->reduce(RooFit::Cut(cutStr.c_str()), RooFit::SelectVars(vars), RooFit::Name(ds->GetName()), RooFit::Title(ds->GetTitle()))));
      if (!data) { std::cout << "[ERROR] Skimmed dataset " << ds->GetName() << " is NULL!" << std::endl; return false; }
      if (data->sumEntries()==0){ std::cout << "[WARNING] Skimmed dataset " <<  ds->GetName() << " is empty!" << std::endl; }
      std::cout << "[INFO] Reduced dataset " << ds->GetName() << " from " << ds->sumEntries() << " to " << data->sumEntries() << " entries" << std::endl;
      addString(tmpWS, "skimDS", cutStr); // Save the cut expression for bookkeeping
      tmpWS.import(*data);
    }
    clearWorkspace(ws.second);
    copyWorkspace(ws.second, tmpWS);
  }
  //
  return true;
};


bool invertEtaAndFill(RooDataSet& dsPA, RooDataSet& dsPbp)
{
  const auto& ds = dsPbp;
  // Find the rapidity related variables
  StringSet_t varNames;
  const auto& iniVars = ds.get();
  auto parIt = std::unique_ptr<TIterator>(iniVars->createIterator());
  for (auto itp = parIt->Next(); itp!=NULL; itp = parIt->Next()) {
    const std::string& vName = itp->GetName();
    if (vName.find("_Rap")!=std::string::npos || vName.find("_Eta")!=std::string::npos) {
      varNames.insert(vName);
    }
  }
  // Fill the PA dataset
  if (!varNames.empty()) {
    for(int i = 0; i < ds.numEntries(); i++){
      auto set = *ds.get(i);
      // Invert the rapidities
      for (const auto& vName : varNames) {
        const auto& var = dynamic_cast<RooRealVar*>(set.find(vName.c_str()));
        var->setVal(-1.0*var->getVal());
      }
      // Add the new event
      if (ds.weight() <= 0.0) { std::cout << "[ERROR] invertEtaAndFill: Weight is negative ( " << ds.weight() << " )" << std::endl; return false; }
      dsPA.add(set, ds.weight());
    }
  }
  else { dsPA.append(dsPbp); }
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
      auto ds_PA = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(ds_pPb->Clone(dsTag_PA.c_str())));
      if (!ds_PA || ds_PA->sumEntries()==0) { std::cout << "[ERROR] RooDataSet " << dsTag_PA << " was not created!" << std::endl; return false; }
      // Invert the rapidity of Pbp dataset and fill the PA dataset
      if (!invertEtaAndFill(*ds_PA, *ds_Pbp)) { return false; }
      // Check the consistency of the combined dataset
      if (ds_PA->numEntries()!=(ds_pPb->numEntries()+ds_Pbp->numEntries())) { std::cout << "[ERROR] Number of entries for the combined " << dsTag_PA << " dataset is inconsistent!" << std::endl; return false; }
      if (std::abs(ds_PA->sumEntries()-(ds_pPb->sumEntries()+ds_Pbp->sumEntries()))>0.05*(ds_pPb->sumEntries()+ds_Pbp->sumEntries())) {
        std::cout << "[ERROR] Number of weighted entries for the combined " << dsTag_PA << " dataset ( " << ds_PA->sumEntries() << " , " << ds_pPb->sumEntries() << " , " << ds_Pbp->sumEntries() << " ) is inconsistent!" << std::endl; return false;
      }
      // Import the new PA dataset
      ws_PA.import(*ds_PA);
      if (ws_PA.data(dsTag_PA.c_str())==NULL) { std::cout << "[ERROR] RooDataSet " << dsTag_PA << " was not imported!" << std::endl; return false; }
      else { std::cout << "[INFO] RooDataSet " << dsTag_PA << " was created with " << ds_PA->numEntries() << " entries!" << std::endl; }
      if (contain(info.StrS.at("DSTAG"), sample_pPb)) { info.StrS.at("DSTAG").insert(sample_PA); }
    }
  }
  //
  return true;
};


#endif
