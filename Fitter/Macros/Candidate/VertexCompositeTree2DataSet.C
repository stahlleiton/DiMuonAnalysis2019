// -*- C++ -*-
//
// Package:    Fitter
//
/*
 Description: VertexCompositeTree to RooDataSet converter.
 Implementation:
 This program create RooDataSets from VertexCompositeTrees.
 */
// Original Author:  Andre Stahl,
//         Created:  Feb 17 19:08 CET 2019
//
//
#ifndef Candidate_VertexCompositeTree2DataSet_C
#define Candidate_VertexCompositeTree2DataSet_C


#include "TROOT.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TMessageHandler.h"

#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooArgSet.h"

#include <string>
#include <memory>
#include <vector>

#include "../../../Utilities/Ntuple/VertexCompositeTree.h"
#include "../../../Utilities/RunInfo/eventUtils.h"
#include "../../../Utilities/dataUtils.h"
#include "../Utilities/initClasses.h"


bool checkVertexCompositeDS ( const RooDataSet& ds , const std::string& analysis );


bool VertexCompositeTree2DataSet(RooWorkspaceMap_t& workspaces, const StringVectorMap_t& fileInfo, const GlobalInfo& info, const bool& updateDS)
{
  StringVector_t outputFileNames;
  const auto& chaDir = info.Par.at("channelDir");
  const auto& outputFileDir  = fileInfo.at("outputFileDir");
  const auto& inputFileNames = fileInfo.at("inputFileNames");
  const auto& dsNames        = fileInfo.at("dsNames");
  const bool& isData = (dsNames[0].rfind("DATA", 0)==0);
  const bool& isMC   = (dsNames[0].rfind("MC", 0)==0);
  for (const auto& tag : dsNames) {
    std::string o = (outputFileDir[0] + chaDir + "/") + "DATASET_" + tag + ".root"; 
    if (gSystem->AccessPathName(o.c_str())) { makeDir(outputFileDir[1] + chaDir + "/"); o = (outputFileDir[1] + chaDir + "/") + "DATASET_" + tag + ".root"; }
    outputFileNames.push_back(o);
  }
  // Extract Input Information
  const auto& PD = info.Par.at("PD");
  const auto& type = info.Par.at("analysis");
  if (type.rfind("CandTo", 0)!=0) { std::cout << "[ERROR] The analysis: " << type << " is not supported!" << std::endl; return false; }
  // Create RooDataSets
  std::vector< std::unique_ptr<RooDataSet> > dataOS, dataSS;
  bool createDS = updateDS;
  bool doSS = true;
  // Check if RooDataSets exist and are not corrupt
  for (uint i=0; i<outputFileNames.size(); i++) {
    if ( !gSystem->AccessPathName(outputFileNames[i].c_str()) ) {
      std::cout << "[INFO] Loading RooDataSets from " << outputFileNames[i] << std::endl;
      auto dbFile = std::unique_ptr<TFile>(TFile::Open(outputFileNames[i].c_str(),"READ"));
      if (!dbFile || !dbFile->IsOpen() || dbFile->IsZombie()) { std::cout << "[ERROR] File: " << outputFileNames[i] << " is corrupted!" << std::endl; return false; }
      dataOS.push_back( std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dbFile->Get(Form("dOS_RAW_%s", dsNames[i].c_str())))) );
      dataSS.push_back( std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dbFile->Get(Form("dSS_RAW_%s", dsNames[i].c_str())))) );
      if (dataOS[i]==NULL || checkVertexCompositeDS(*dataOS[i], type)==false) { createDS = true; }
      if (dataSS[i]==NULL || checkVertexCompositeDS(*dataSS[i], type)==false) { doSS = false; }
      dbFile->Close();
    }
    else { createDS = true; break; }
  }
  if (createDS) {
    ///// Input Forest
    //
    // Find directory in ROOT file
    std::string dirName = "";
    findDirInFile(dirName, inputFileNames[0]);
    if (dirName=="") { std::cout << "[ERROR] Failed to find the directory in: " << inputFileNames[0] << std::endl; return false; }
    const auto& dirNameSS = (dirName=="dimucontana_mc" ? "dimucontana_wrongsign_mc" : (dirName+"_wrongsign"));
    //
    auto candOSTree = std::unique_ptr<VertexCompositeTree>(new VertexCompositeTree());
    if (!candOSTree->GetTree(inputFileNames, dirName)) return false;
    const auto& nentries = candOSTree->GetEntries();
    auto candSSTree = std::unique_ptr<VertexCompositeTree>(new VertexCompositeTree());
    doSS = candSSTree->GetTree(inputFileNames, dirNameSS);
    if (doSS==false) { std::cout << "[INFO] Tree: " << dirName+"_wrongsign not found, will be ignored!" << std::endl; }
    if (doSS && candSSTree->GetEntries() != nentries) { std::cout << "[ERROR] Inconsistent number of entries in candTreeSS!" << std::endl; return false; }
    //
    ///// RooDataSet Variables
    auto candMass = RooRealVar ( "Cand_Mass"  , "Candidate Mass"  ,    -1.0 ,     100.0 ,  "GeV/c^{2}" );
    auto candPt   = RooRealVar ( "Cand_Pt"    , "Candidate p_{T}" ,    -1.0 ,    1000.0 ,  "GeV/c"     );
    auto candRap  = RooRealVar ( "Cand_Rap"   , "Candidate y"     ,   -10.0 ,      10.0 ,  ""          );
    auto candAPhi = RooRealVar ( "Cand_APhi"  , "Candidate #phi_{#mu}" , -10.0 ,   10.0 ,  ""          );
    auto candDLen = RooRealVar ( "Cand_DLen"  , "Candidate c#tau" , -1000.0 ,    1000.0 ,  "mm"        );
    auto candDLenErr = RooRealVar ( "Cand_DLenErr" , "Candidate c#tau error" , -1.0 , 1000.0 , "mm"    );
    auto candDLenGen = RooRealVar ( "Cand_DLenGen" , "GEN candidate c#tau" , -1000.0 , 1000.0 , "mm"   );
    auto dau1Pt   = RooRealVar ( "Dau1_Pt"    , "Daughter1 p_{T}" ,  -1.0 ,       100.0 ,  "GeV/c"     );
    auto dau2Pt   = RooRealVar ( "Dau2_Pt"    , "Daughter2 p_{T}" ,  -1.0 ,       100.0 ,  "GeV/c"     );
    auto dau1Eta  = RooRealVar ( "Dau1_Eta"   , "Daughter1 #eta"  , -10.0 ,        10.0 ,  ""          );
    auto dau2Eta  = RooRealVar ( "Dau2_Eta"   , "Daughter2 #eta"  , -10.0 ,        10.0 ,  ""          );
    auto cent     = RooRealVar ( "Centrality" , "Centrality"      ,  -1.0 ,      1000.0 ,  ""          );
    auto nTrk     = RooRealVar ( "NTrack"     , "Number of Tracks",  -1.0 ,   1000000.0 ,  ""          );
    auto candQual = RooRealVar ( "Cand_Qual"  , "Candidate Quality", -1.0 ,        10.0 ,  ""          );
    auto candTrig = RooRealVar ( "Cand_Trig"  , "Candidate Trigger", -1.0 ,      1000.0 ,  ""          );
    auto candVtxP = RooRealVar ( "Cand_VtxP"  , "Cand. Vertex Prob." , -1.0 ,      10.0 ,  ""          );
    auto weight   = RooRealVar ( "Weight"     , "Weight"          ,  -1.0 , 100000000.0 ,  ""          );
    auto isSwap   = RooCategory( "Cand_IsSwap" , "Candidate IsSwap");
    isSwap.defineType("None", -1); isSwap.defineType("No", 0); isSwap.defineType("Yes", 1);
    //
    auto cols = RooArgSet(candMass, candPt, candRap, candDLen, candDLenErr, cent, nTrk);
    cols.add(dau1Pt); cols.add(dau1Eta); cols.add(dau2Pt); cols.add(dau2Eta);
    cols.add(candAPhi); cols.add(candQual); cols.add(candTrig); cols.add(candVtxP);
    if (isMC) { cols.add(isSwap); cols.add(candDLenGen); }
    //
    ///// Initiliaze RooDataSets
    dataOS.clear(); dataSS.clear();
    for (uint i=0; i<dsNames.size(); i++) {
      if (isMC) {
        cols.add(weight);
        dataOS.push_back( std::unique_ptr<RooDataSet>(new RooDataSet(Form("dOS_RAW_%s", dsNames[i].c_str()), "dOS", cols, RooFit::WeightVar(weight))) );
        if (doSS) { dataSS.push_back( std::unique_ptr<RooDataSet>(new RooDataSet(Form("dSS_RAW_%s", dsNames[i].c_str()), "dSS", cols, RooFit::WeightVar(weight))) ); }
      }
      else {
        dataOS.push_back( std::unique_ptr<RooDataSet>(new RooDataSet(Form("dOS_RAW_%s", dsNames[i].c_str()), "dOS", cols)) );
        if (doSS) { dataSS.push_back( std::unique_ptr<RooDataSet>(new RooDataSet(Form("dSS_RAW_%s", dsNames[i].c_str()), "dSS", cols)) ); }
      }
    }
    //
    // Determine the collision system of the sample
    std::string evtCol = "";
    if (dsNames[0].rfind("_")!=std::string::npos) { evtCol = dsNames[0].substr(dsNames[0].rfind("_")+1); }
    if (evtCol=="") { std::cout << "[ERROR] Could not determine the collision system in the sample" << std::endl; return false; }
    //
    // Determine the event selections
    std::vector<uint> evtSelIdx;
    if      (PD=="UPC" && evtCol=="PbPb5Y18") { evtSelIdx.push_back(PbPb::R5TeV::Y2018::primaryVertexFilter); evtSelIdx.push_back(PbPb::R5TeV::Y2018::clusterCompatibilityFilter); }
    else if (PD=="UPC" && evtCol=="PbPb5Y15") { evtSelIdx.push_back(PbPb::R5TeV::Y2015::primaryVertexFilter); evtSelIdx.push_back(PbPb::R5TeV::Y2018::clusterCompatibilityFilter); }
    else if (evtCol=="PP13Y18" ) { evtSelIdx.push_back(pp::R13TeV::Y2018::colEvtSel); }
    else if (evtCol=="PP5Y17"  ) { evtSelIdx.push_back(pp::R5TeV::Y2017::colEvtSel); }
    else if (evtCol=="PbPb5Y18") { evtSelIdx.push_back(PbPb::R5TeV::Y2018::colEvtSel); }
    else if (evtCol=="PbPb5Y15") { evtSelIdx.push_back(PbPb::R5TeV::Y2015::colEvtSel); }
    else if (evtCol.rfind("8Y16")!=std::string::npos) { evtSelIdx.push_back(pPb::R8TeV::Y2016::colEvtSel); }
    if (evtSelIdx.empty()) { std::cout << "[ERROR] Could not determine the event selection index for the sample" << std::endl; return false; }
    //
    // Determine the trigger paths
    std::vector<uint> trigIdx;
    std::map<uint, std::string> allTrig;
    if      (evtCol=="PP13Y18" ) { trigIdx = pp::R13TeV::Y2018::HLTBitsFromPD(PD);  allTrig = pp::R13TeV::Y2018::HLTBits();  }
    else if (evtCol=="PP5Y17"  ) { trigIdx = pp::R5TeV::Y2017::HLTBitsFromPD(PD);   allTrig = pp::R5TeV::Y2017::HLTBits();   }
    else if (evtCol=="PbPb5Y18") { trigIdx = PbPb::R5TeV::Y2018::HLTBitsFromPD(PD); allTrig = PbPb::R5TeV::Y2018::HLTBits(); }
    else if (evtCol=="PbPb5Y15") { trigIdx = PbPb::R5TeV::Y2015::HLTBitsFromPD(PD); allTrig = PbPb::R5TeV::Y2015::HLTBits(); }
    else if (evtCol.rfind("8Y16")!=std::string::npos) { trigIdx = pPb::R8TeV::Y2016::HLTBitsFromPD(PD); allTrig = pPb::R8TeV::Y2016::HLTBits(); }
    if (trigIdx.empty()) { std::cout << "[ERROR] Could not determine the trigger index for the sample" << std::endl; return false; }
    //
    ///// Iterate over the Input Ntuple
    int treeIdx = -1;
    std::cout << "[INFO] Starting to process " << nentries << " nentries" << std::endl;
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
      //
      // Get the entry in the trees
      if (candOSTree->GetEntry(jentry)<0) { std::cout << "[ERROR] CandOS Tree invalid entry!"  << std::endl; return false; }
      if (doSS && candSSTree->GetEntry(jentry)<0) { std::cout << "[ERROR] CandSS Tree invalid entry!"  << std::endl; return false; }
      // 
      // Check if different trees are synchronized
      if (doSS && candSSTree->RunNb()!=candOSTree->RunNb()    ) { std::cout << "[ERROR] Run number ("<<candSSTree->RunNb()<<") in wrong-sign tree is invalid!" << std::endl; return false; }
      if (doSS && candSSTree->EventNb()!=candOSTree->EventNb()) { std::cout << "[ERROR] Event number ("<<candSSTree->EventNb()<<") in wrong-sign tree is invvalid!" << std::endl; return false; }
      // 
      if (candOSTree->GetTreeNumber()!=treeIdx) {
        treeIdx = candOSTree->GetTreeNumber();
        std::cout << "[INFO] Processing Root File: " << inputFileNames[treeIdx] << std::endl;
      }
      //
      loadBar(jentry, nentries);
      //
      // For pPb, find out the run on data
      if (isData && evtCol.rfind("8Y16")!=std::string::npos) {
        if      (candOSTree->RunNb() >= 285410 && candOSTree->RunNb() <= 285951) evtCol = "Pbp8Y16"; // for Pbp8Y16
        else if (candOSTree->RunNb() >= 285952 && candOSTree->RunNb() <= 286504) evtCol = "pPb8Y16"; // for pPb8Y16
      }
      //
      // Event Based Information
      //
      // Apply Event Filters
      bool evtSel = true;
      for (const auto& idx : evtSelIdx) { evtSel = (evtSel && candOSTree->evtSel()[idx]); }
      if (evtSel==false) continue;
      //
      // Check Trigger Decisions
      std::map<uint, bool> trigMap;
      for (const auto& idx : allTrig) { trigMap[idx.first] = candOSTree->trigHLT()[idx.first]; }
      uint evtTrig = 0;
      for (const auto& idx : trigMap) { if (idx.second) { evtTrig += std::pow(2.0, idx.first); } }
      if (evtTrig==0) continue;
      //
      // Apply Trigger Selection
      bool trigSel = false;
      for (const auto& idx : trigIdx) { trigSel = (trigSel || trigMap.at(idx)); }
      if (isData && trigSel==false) continue;
      //
      // Candidate Based Information
      //
      for (uint iC = 0; iC < candOSTree->candSize(); iC++) {
        //
        // Check if we are doing dimuon analysis
        if (type=="CandToMuMu") {
          //
          // Apply loose muon acceptance
          const auto& d1Pt  = candOSTree->pTD1()[iC];
          const auto& d2Pt  = candOSTree->pTD2()[iC];
          const auto& d1Eta = candOSTree->EtaD1()[iC];
          const auto& d2Eta = candOSTree->EtaD2()[iC];
          const auto& d1P   = d1Pt*std::cosh(d1Eta);
          const auto& d2P   = d2Pt*std::cosh(d2Eta);
          if ( (std::abs(d1Eta) > 2.4 || d1P < 2.5) || (std::abs(d2Eta) > 2.4 || d2P < 2.5) ) continue;
          //
          // Apply loose muon quality cuts
          const auto& centV = candOSTree->centrality();
          const auto& isTightCand  = candOSTree->tightCand(iC);
          const auto& isHybridCand = candOSTree->hybridCand(iC);
          const auto& isSoftCand   = candOSTree->softCand(iC);
          int candQ = 0; if (isSoftCand) { candQ += 1; }; if (isHybridCand) { candQ += 2; }; if (isTightCand) { candQ += 4; }
          if (candQ==0) continue;
          if (isData  && evtCol.rfind("PbPb", 0)==0 && centV<50 && candQ==1) continue;
          //
          // Apply muon trigger matching
          bool matchTrig = false;
          for (const auto& idx : trigIdx) {
	    bool isTrigMatch = false;
	    if      (allTrig.at(idx)=="Muon"  ) { isTrigMatch = (trigMap.at(idx) && candOSTree->trigCand(idx, iC, true )); }
            else if (allTrig.at(idx)=="DiMuon") { isTrigMatch = (trigMap.at(idx) && candOSTree->trigCand(idx, iC, false)); }
	    else { matchTrig = true; break; }
	    matchTrig = (matchTrig || isTrigMatch);
          }
          if (isData && matchTrig==false) continue;
          //
          // Apply double muon acceptance
          const auto& pT   = candOSTree->pT()[iC];
          const auto& mass = candOSTree->mass()[iC];
          const bool& massCut = ( (mass > 0.0 && mass < 1.6) || (mass > 4.5 && mass < 6.0) || (mass > 15.0 && mass < 55.0) );
          if (isData && PD!="UPC" && pT>1.0 && massCut) continue;
          //
          // Apply MC cuts
          if (isMC) {
            // Check that candidate is matched to gen
            if (candOSTree->matchGEN()[iC]==false) continue;
            // Check the PID of the matched gen particle
            const auto& tmp = dsNames[0].substr(dsNames[0].find("_")+1);
            auto par = tmp.substr(0, tmp.find("_"));
            for (const auto& p : MASS) { if (par.find(p.first)!=std::string::npos) { par = p.first; break; } }
            if (!contain(MASS, par)) { std::cout << "[ERROR] MC particle "<<par<<" is not valid!" << std::endl; return false; }
            if (fabs(candOSTree->idmom_reco()[iC])!=int(MASS.at(par).at("PID"))) continue;
          }
          //
          // Store muon trigger matching info
          uint trigM = 0;
	  for (const auto& idx : allTrig) {
            bool isTrigMatch = true;
            if      (idx.second=="Muon"  ) { isTrigMatch = (trigMap.at(idx.first) && candOSTree->trigCand(idx.first, iC, true )); }
            else if (idx.second=="DiMuon") { isTrigMatch = (trigMap.at(idx.first) && candOSTree->trigCand(idx.first, iC, false)); }
            if (isTrigMatch) { trigM += std::pow(2.0, idx.first); }
          }
	  if (trigM==0) continue;
          //
          // Compute the pseudo-proper-decay length
          const auto& p    = pT*std::cosh(candOSTree->eta()[iC]);
          const auto& rap  = candOSTree->y()[iC];
          const auto& massJPsi = MASS.at("JPsi").at("Val");
          const auto& dLen = (candOSTree->V3DDecayLength()[iC] * candOSTree->V3DCosPointingAngle()[iC])*(massJPsi/p)*10.0;
          const auto& dLenErr = (candOSTree->V3DDecayLengthError()[iC])*(massJPsi/p)*10.0;
          //
          // Compute the azimuthal asymmetry
          const auto& aPhi = candOSTree->phiAsym(iC);
          //
          // Set the variables
          candMass.setVal( mass  );
          candPt.setVal  ( pT    );
          candRap.setVal ( rap   );
          candDLen.setVal( dLen  );
          candDLenErr.setVal( dLenErr );
          candAPhi.setVal( aPhi  );
          candQual.setVal( candQ );
	  candTrig.setVal( trigM );
          candVtxP.setVal( candOSTree->VtxProb()[iC] );
          dau1Pt.setVal  ( d1Pt  );
          dau2Pt.setVal  ( d2Pt  );
          dau1Eta.setVal ( d1Eta );
          dau2Eta.setVal ( d2Eta );
          cent.setVal    ( centV );
          nTrk.setVal    ( candOSTree->Ntrkoffline() );
          weight.setVal  ( 1.0   );
          isSwap.setLabel( "None");
	  candDLenGen.setVal( -1.0 );
          if (isMC) {
            weight.setVal( candOSTree->weight_gen() );
            isSwap.setLabel( candOSTree->isSwap() ? "Yes" : "No" );
	    const auto& iGen = candOSTree->GenIdx(iC);
	    if (iGen>=0) {
	      const auto& pGen = candOSTree->pT_gen()[iGen]*std::cosh(candOSTree->eta_gen()[iGen]);
	      const auto& dLenGen = (candOSTree->V3DDecayLength_gen()[iGen] * std::cos(candOSTree->V3DPointingAngle_gen()[iGen]))*(massJPsi/pGen)*10.0;
	      candDLenGen.setVal( dLenGen );
	    }
          }
          //
          // Fill the RooDataSets
          for (uint i=0; i<dsNames.size(); i++) {
            if (dsNames[i].rfind(evtCol)!=std::string::npos) { dataOS[i]->addFast(cols); }
          }
        }
      }
      //
      for (uint iC = 0; iC < (doSS ? candSSTree->candSize() : 0); iC++) {
        //
        // Check if we are doing dimuon analysis
        if (type=="CandToMuMu") {
          //
          // Apply loose muon acceptance
          const auto& d1Pt  = candSSTree->pTD1()[iC];
          const auto& d2Pt  = candSSTree->pTD2()[iC];
          const auto& d1Eta = candSSTree->EtaD1()[iC];
          const auto& d2Eta = candSSTree->EtaD2()[iC];
          const auto& d1P   = d1Pt*std::cosh(d1Eta);
          const auto& d2P   = d2Pt*std::cosh(d2Eta);
          if ( (std::abs(d1Eta) > 2.4 || d1P < 2.5) || (std::abs(d2Eta) > 2.4 || d2P < 2.5) ) continue;
          //
          // Apply loose muon quality cuts
          const auto& centV = candSSTree->centrality();
          const auto& isTightCand  = candSSTree->tightCand(iC);
          const auto& isHybridCand = candSSTree->hybridCand(iC);
          const auto& isSoftCand   = candSSTree->softCand(iC);
          int candQ = 0; if (isSoftCand) { candQ += 1; }; if (isHybridCand) { candQ += 2; }; if (isTightCand) { candQ += 4; }
          if (candQ==0) continue;
          if (isData && evtCol.rfind("PbPb", 0)==0 && centV<50 && candQ==1) continue;
          //
          // Apply muon trigger matching
          bool matchTrig = false;
          for (const auto& idx : trigIdx) {
            bool isTrigMatch = false;
            if      (allTrig.at(idx)=="Muon"  ) { isTrigMatch = (trigMap.at(idx) && candSSTree->trigCand(idx, iC, true )); }
            else if (allTrig.at(idx)=="DiMuon") { isTrigMatch = (trigMap.at(idx) && candSSTree->trigCand(idx, iC, false)); }
            else { matchTrig = true; break; }
            matchTrig = (matchTrig || isTrigMatch);
          }
          if (isData && matchTrig==false) continue;
          //
          // Apply double muon acceptance
          const auto& pT   = candSSTree->pT()[iC];
          const auto& mass = candSSTree->mass()[iC];
          const bool& massCut = ( (mass > 0.0 && mass < 1.6) || (mass > 4.5 && mass < 6.0) || (mass > 15.0 && mass < 55.0) );
          if (isData && PD!="UPC" && pT>1.0 && massCut) continue;
          //
          // Store muon trigger matching info
          uint trigM = 0;
          for (const auto& idx : allTrig) {
            bool isTrigMatch = true;
            if      (idx.second=="Muon"  ) { isTrigMatch = (trigMap.at(idx.first) && candSSTree->trigCand(idx.first, iC, true )); }
            else if (idx.second=="DiMuon") { isTrigMatch = (trigMap.at(idx.first) && candSSTree->trigCand(idx.first, iC, false)); }
            if (isTrigMatch) { trigM += std::pow(2.0, idx.first); }
          }
          if (trigM==0) continue;
          //
          // Compute the pseudo-proper-decay length
          const auto& p    = pT*std::cosh(candSSTree->eta()[iC]);
          const auto& rap  = candSSTree->y()[iC];
          const auto& massJPsi = MASS.at("JPsi").at("Val");
          const auto& dLen = (candSSTree->V3DDecayLength()[iC] * candSSTree->V3DCosPointingAngle()[iC])*(massJPsi/p)*10.0;
          const auto& dLenErr = (candSSTree->V3DDecayLengthError()[iC])*(massJPsi/p)*10.0;
          //
          // Compute the azimuthal asymmetry  
          const auto& aPhi = candSSTree->phiAsym(iC);
          //
          // Set the variables
          candMass.setVal( mass  );
          candPt.setVal  ( pT    );
          candRap.setVal ( rap   );
          candDLen.setVal ( dLen  );
          candDLenErr.setVal( dLenErr );
          candAPhi.setVal( aPhi  );
          candQual.setVal( candQ );
	  candTrig.setVal( trigM );
      	  candVtxP.setVal( candSSTree->VtxProb()[iC] );
          dau1Pt.setVal  ( d1Pt  );
          dau2Pt.setVal  ( d2Pt  );
          dau1Eta.setVal ( d1Eta );
          dau2Eta.setVal ( d2Eta );
          cent.setVal    ( centV );
          nTrk.setVal    ( candSSTree->Ntrkoffline() );
          weight.setVal  ( 1.0 );
          isSwap.setLabel("None");
	  candDLenGen.setVal( -1.0 );
          //
          // Fill the RooDataSets
          for (uint i=0; i<dsNames.size(); i++) {
            if (dsNames[i].rfind(evtCol)!=std::string::npos) { dataSS[i]->addFast(cols); }
          }
        }
      }
    }
    //// Save the RooDataSets
    for (uint i=0; i<dsNames.size(); i++) {
      std::cout << "[INFO] Creating output file: " << outputFileNames[i] << std::endl;
      // Write the datasets
      auto dbFile = std::unique_ptr<TFile>(TFile::Open(outputFileNames[i].c_str(),"RECREATE"));
      dbFile->cd();
      std::cout << "[INFO] Converting datasets " << dsNames[i] << " to tree store" << std::endl;
      dataOS[i]->convertToTreeStore();
      if (doSS) { dataSS[i]->convertToTreeStore(); }
      std::cout << "[INFO] Saving datasets " << dsNames[i] << " in " << outputFileNames[i] << std::endl;
      dataOS[i]->Write(Form("dOS_RAW_%s", dsNames[i].c_str()));
      if (doSS) { dataSS[i]->Write(Form("dSS_RAW_%s", dsNames[i].c_str())); }
      std::cout << "[INFO] Closing output file " << outputFileNames[i] << std::endl;
      dbFile->Write(); dbFile->Close();
      std::cout << "[INFO] Converting datasets " << dsNames[i] << " back to vector store" << std::endl;
      dataOS[i]->convertToVectorStore();
      if (doSS) { dataSS[i]->convertToVectorStore(); }
    }
  }
  // Import datasets to the workspaces
  for (uint i=0; i<dsNames.size(); i++) {
    if (!dataOS[i]) { std::cout << "[ERROR] " << dsNames[i] << " OS dataset was not found" << std::endl; return false; }
    if (dataOS[i]->numEntries()==0) { std::cout << "[ERROR] " << dsNames[i] << " OS dataset is empty!" << std::endl; return false; }
    workspaces[dsNames[i]].import(*dataOS[i]);
    if (doSS) {
      if (!dataSS[i]) { std::cout << "[ERROR] " << dsNames[i] << " SS dataset was not found" << std::endl; return false; }
      if (dataSS[i]->numEntries()==0) { std::cout << "[WARNING] " << dsNames[i] << " SS dataset is empty!" << std::endl; }
      workspaces[dsNames[i]].import(*dataSS[i]);
    }
  }
  dataOS.clear(); dataSS.clear();
  return true;
};


bool checkVertexCompositeDS(const RooDataSet& ds, const std::string& analysis)
{
  if (ds.numEntries()==0 || ds.sumEntries()==0) { std::cout << "[WARNING] Original dataset: " << ds.GetName() << " is empty, will remake it!" << std::endl; return false; }
  const bool& isMC = (std::string(ds.GetName()).find("_MC")!=std::string::npos);
  const auto& row = ds.get();
  if (analysis.rfind("CandTo", 0)==0) {
    if (
        ( row->find("Cand_Mass") !=0 ) &&
        ( row->find("Cand_Pt") !=0   ) &&
        ( row->find("Cand_Rap") !=0  ) &&
        ( row->find("Cand_DLen") !=0 ) &&
        ( row->find("Cand_DLenErr") !=0 ) &&
        ( row->find("Cand_APhi") !=0 ) &&
        ( row->find("Cand_Qual") !=0 ) &&
	( row->find("Cand_Trig") !=0 ) &&
	( row->find("Cand_VtxP") !=0 ) &&
        ( !isMC || row->find("Cand_IsSwap") !=0  ) &&
        ( !isMC || row->find("Cand_DLenGen") !=0 ) &&
        ( row->find("Dau1_Pt") !=0  ) &&
        ( row->find("Dau1_Eta") !=0 ) &&
        ( row->find("Dau2_Pt") !=0  ) &&
        ( row->find("Dau2_Eta") !=0 ) &&
        ( row->find("Centrality") !=0 ) &&
        ( row->find("NTrack") !=0 )
        )
      { return true; }
    else { std::cout << "[WARNING] Original dataset: " << ds.GetName() << " is corrupted, will remake it!" << std::endl; }
  }
  return false;
};


#endif
