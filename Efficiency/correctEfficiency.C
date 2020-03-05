#if !defined(__CINT__) || defined(__MAKECINT__)
// Auxiliary Headers
#include "../Utilities/Ntuple/VertexCompositeTree.h"
#include "../Utilities/RunInfo/eventUtils.h"
#include "../Utilities/dataUtils.h"
#include "util.h"
#include "tnp_weight_lowPt.h"
// ROOT headers
#include "TSystem.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TVectorD.h"
#include "TMath.h"
#include "TH1D.h"
// c++ headers
#include <dirent.h>
#include <iostream>
#include <map>
#include <vector>
#include <tuple>
#include <string>

#endif


// ------------------ TYPE -------------------------------
using TnPVec_t     =  std::map< std::string , std::vector< double > >;
using AnaBinPair_t =  std::pair< AnaBin_t , std::string >;
using AnaVarMap_t  =  std::map< AnaBinPair_t, Var_t >;
using Unc1DVec_t   =  std::map< std::string , TVectorD >;
using Unc1DMap_t   =  std::map< std::string , std::map< std::string , std::map< std::string , std::map< AnaBinPair_t , Unc1DVec_t > > > >;
using TH1DVec_t    =  std::map< std::string , std::vector< std::tuple< TH1D , TH1D , double > > >;
using TH1DMap_t    =  std::map< std::string , std::map< std::string , std::map< std::string , std::map< AnaBinPair_t , TH1DVec_t > > > >;
using EffVec_t     =  std::map< std::string , std::vector< TEfficiency > >;
using EffMap_t     =  std::map< std::string , std::map< std::string , std::map< std::string , std::map< AnaBinPair_t , EffVec_t > > > >;
using VarMap_t     =  std::map< std::string , double >;
using CorrMap_t    =  std::map< std::string , uint >;


// ------------------ FUNCTION -------------------------------
TnPVec_t getTnPScaleFactors  ( const double& ptD1 , const double& etaD1 , const double& ptD2 , const double& etaD2 , const CorrMap_t& corrType , const bool& incMuTrig );
bool     getTnPUncertainties ( Unc1DVec_t& unc , const EffVec_t& eff );
bool     getTnPUncertainties ( Unc1DMap_t& unc , const EffMap_t& eff );
void     initEff1D           ( TH1DMap_t& h , const AnaVarMap_t& binMap , const CorrMap_t& corrType );
bool     fillEff1D           ( TH1DVec_t& h , const bool& den_pass , const bool& num_pass , const double& xVar , const TnPVec_t& sfTnP , const double& evtWeight );
bool     fillEff1D           ( TH1DMap_t& h , const bool& den_pass , const bool& num_pass, const std::string& sample, const std::string& col, const std::string& type,
			       const VarMap_t& var, const TnPVec_t& sfTnP, const double& evtWeight );
bool     loadEff1D           ( EffMap_t& eff, const TH1DMap_t& h );
void     mergeEff            ( EffMap_t& eff );
void     writeEff            ( TFile& file , const EffMap_t& eff , const Unc1DMap_t& unc , const std::string& mainDirName );
void     saveEff             ( const std::string& outDir , const EffMap_t& eff1D , const Unc1DMap_t& unc );
void     setGlobalWeight     ( TH1DMap_t& h , const double& weight , const std::string& sample , const std::string& col );
bool     checkWeights        ( const TH1& pass , const TH1& total );
const char* clStr            ( const std::string& in );


// ------------------ GLOBAL ------------------------------- 
//
// Collision System
const std::vector<std::string> COLL_ = { "pPb8Y16", "Pbp8Y16", "PA8Y16" };
//
// Efficiency Categories
const std::vector< std::string > EFFTYPE_ = {"Acceptance", "Efficiency_Total", "Efficiency_DecayCut_90"}; //"Efficiency_DecayCut_85", "Efficiency_DecayCut_95", "Efficiency_DecayCut_99"};
//
// Correction Categories
const CorrMap_t corrType_ = {
  { "NoCorr"           , 1   },
  { "TnP_Nominal"      , 1   },
  /*
  { "TnP_Stat_TrkM"    , 2   },
  { "TnP_Stat_MuID"    , 2   },
  { "TnP_Stat_Trig"    , 2   },
  { "TnP_Syst_TrkM"    , 2   },
  { "TnP_Syst_MuID"    , 2   },
  { "TnP_Syst_Trig"    , 2   },
  { "TnP_Syst_BinTrkM" , 1   },
  */
};
//
// Input Files for analysis
const std::string path_MC = "/Users/andre/DiMuonAnalysis2019/Tree";
const std::map< std::string , std::string > inputFileMap_ =
  {
   {"MC_JPsiPR_pPb"       , Form("%s/%s", path_MC.c_str(), "VertexCompositeTree_JPsiToMuMu_pPb-Bst_pPb816Summer16_DiMuMC.root")   },
   {"MC_JPsiPR_Pbp"       , Form("%s/%s", path_MC.c_str(), "VertexCompositeTree_JPsiToMuMu_Pbp-Bst_pPb816Summer16_DiMuMC.root")   },
   {"MC_Psi2SPR_pPb"      , Form("%s/%s", path_MC.c_str(), "VertexCompositeTree_Psi2SToMuMu_pPb-Bst_pPb816Summer16_DiMuMC.root")  },
   {"MC_Psi2SPR_Pbp"      , Form("%s/%s", path_MC.c_str(), "VertexCompositeTree_Psi2SToMuMu_Pbp-Bst_pPb816Summer16_DiMuMC.root")  },
   {"MC_JPsiNoPR_pPb"     , Form("%s/%s", path_MC.c_str(), "VertexCompositeTree_BToPsiToMuMu_pPb-Bst_pPb816Summer16_DiMuMC.root") },
   {"MC_JPsiNoPR_Pbp"     , Form("%s/%s", path_MC.c_str(), "VertexCompositeTree_BToPsiToMuMu_Pbp-Bst_pPb816Summer16_DiMuMC.root") },
   {"MC_Psi2SNoPR_pPb"    , Form("%s/%s", path_MC.c_str(), "VertexCompositeTree_BToPsiToMuMu_pPb-Bst_pPb816Summer16_DiMuMC.root") },
   {"MC_Psi2SNoPR_Pbp"    , Form("%s/%s", path_MC.c_str(), "VertexCompositeTree_BToPsiToMuMu_Pbp-Bst_pPb816Summer16_DiMuMC.root") },
   {"MC_JPsiPR_pPbGEN"    , Form("%s/%s", path_MC.c_str(), "VertexCompositeTree_JPsiToMuMu_pPb-Bst_GENonly_pPb816Summer16_DiMuGENONLY.root")   },
   {"MC_Psi2SPR_pPbGEN"   , Form("%s/%s", path_MC.c_str(), "VertexCompositeTree_Psi2SToMuMu_pPb-Bst_GENonly_pPb816Summer16_DiMuGENONLY.root")  },
   {"MC_JPsiNoPR_pPbGEN"  , Form("%s/%s", path_MC.c_str(), "VertexCompositeTree_BToPsiToMuMu_pPb-Bst_GENonly_pPb816Summer16_DiMuGENONLY.root") },
   {"MC_Psi2SNoPR_pPbGEN" , Form("%s/%s", path_MC.c_str(), "VertexCompositeTree_BToPsiToMuMu_pPb-Bst_GENonly_pPb816Summer16_DiMuGENONLY.root") }
  };
std::map< std::string , std::vector< std::string > > sampleType_;


void correctEfficiency(const std::string workDirName = "Nominal", const std::string PD = "DIMUON")
{
  //
  std::cout << "[INFO] Starting to compute efficiencies" << std::endl;
  //
  // Initialize the Kinematic Bin info
  BinMapMap_t  ANA_BIN;
  if (workDirName == "Nominal") {
    ANA_BIN = BINMAP_Psi2S;
  }
  else if (workDirName == "General") {
    ANA_BIN = BINMAP_General;
  }
  else {
    std::cout << "[ERROR] WorkDirName " << workDirName << " has not been defined" << std::endl; return;
  }
  AnaVarMap_t ANA_BIN_MAP;
  for (const auto& v1 : ANA_BIN) {
    const auto& var1 = v1.first;
    for (const auto& v2 : v1.second) {
      const auto& var2 = v2.first;
      const auto& anaBin = AnaBin_t({var1, var2});
      const auto& v3 = v2.second;
      for (uint i=0; i<v3.size(); i++) {
	const auto& anaBinPair = std::make_pair(anaBin, Form("%s_%d", std::get<0>(v3[i]).c_str(), i));
	ANA_BIN_MAP[anaBinPair] = v3[i];
      }
    }
  }
  //
  // Change the working directory
  const std::string CWD = getcwd(NULL, 0);
  const std::string mainDir = Form("%s/Output/", CWD.c_str());
  gSystem->mkdir(mainDir.c_str(), kTRUE);
  gSystem->ChangeDirectory(mainDir.c_str());
  //
  // Create list of samples
  sampleType_["sample"] = {};
  for (const auto& f : inputFileMap_) { std::string tmp = f.first; tmp = tmp.substr(0, tmp.find_last_of("_"));
    if (std::find(sampleType_.at("sample").begin(), sampleType_.at("sample").end(), tmp)==sampleType_.at("sample").end()) { sampleType_.at("sample").push_back(tmp); }
  }
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Initialize the Correction Type Map
  CorrMap_t corrType;
  for (const auto& cor : corrType_) {
    corrType[cor.first] = cor.second;
  }
  corrType["NoCorr"] = 1;
  corrType["TnP_Nominal"] = 1;
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Declare the histograms for the efficiency
  TH1DMap_t h1D;   // Stores the total and passing histograms separately
  //
  // Initialize the efficiencies
  initEff1D(h1D , ANA_BIN_MAP, corrType);
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Extract all the samples
  std::map< std::string , Long64_t > nentries;
  std::map< std::string , std::unique_ptr< VertexCompositeTree > > trees;
  for (const auto & inputFile : inputFileMap_) {
    const auto& sample = inputFile.first;
    const auto& fileInfo = inputFile.second;
    //
    trees[sample] = std::unique_ptr<VertexCompositeTree>(new VertexCompositeTree());
    const std::string dir = (sample.rfind("GEN")!=std::string::npos ? "dimuana_mc" : "dimucontana_mc");
    if (!trees.at(sample)->GetTree(fileInfo, dir)) return;
    nentries[sample] = trees.at(sample)->GetEntries();
  }
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Loop over the samples
  //
  for (auto & inputFile : inputFileMap_) {
    const auto& sample = inputFile.first;
    auto& tree = trees.at(sample);
    // Loop over the events
    int treeIdx = -1;
    std::cout << "[INFO] Starting to process " << nentries.at(sample) << " nentries" << std::endl;
    for (Long64_t jentry = 0; jentry < nentries.at(sample); jentry++) {
      //
      // Get the entry in the trees
      if (tree->GetEntry(jentry)<0) { std::cout << "[ERROR] Muon Tree invalid entry!"  << std::endl; return; }
      //
      if (tree->GetTreeNumber()!=treeIdx) {
        treeIdx = tree->GetTreeNumber();
        std::cout << "[INFO] Processing " << sample << " using file: " << inputFile.second << std::endl;
      }
      //
      loadBar(jentry, nentries.at(sample));
      //
      // Determine the collision system of the sample
      std::string col = "";
      if (sample.find("Pbp")!=std::string::npos) col = "Pbp8Y16"; // for Pbp
      if (sample.find("pPb")!=std::string::npos) col = "pPb8Y16"; // for pPb
      if (col=="") { std::cout << "[ERROR] Could not determine the collision system in the sample" << std::endl; return; }
      //
      // Determine the type of sample : i.e. MC_DYToMuMu
      auto sampleType = sample;
      sampleType = sampleType.substr(0, (sampleType.find(col)-1));
      //
      // Get the Lumi re-weight for MC (global weight)
      const auto lumi = pPb::R8TeV::Y2016::LumiFromPD(PD, col);
      const double mcWeight = (lumi / tree->GetTreeEntries());
      // Set the global weight only in the first event (i.e. once per sample)
      if (jentry==0) { setGlobalWeight(h1D, mcWeight, sampleType, col); }
      //
      // Define the event weight (set to 1.0 by default)
      double evtWeight = 1.0;
      //
      // Check Event Conditions
      //
      // Determine if the event pass the event filters
      const bool passEventSelection = (sample.find("pPbGEN")==std::string::npos ? tree->evtSel()[0] : true);
      //
      // Check Muon Conditions
      //
      // Loop over the generated muons
      for (ushort iGen = 0; iGen < tree->candSize_gen(); iGen++) {
	// Check the type of particle
        const auto pid = (sample.rfind("MC_JPsi",0)==0 ? 443 : (sample.rfind("MC_Psi2S",0)==0 ? 100443 : -1));
        if (fabs(tree->PID_gen()[iGen])!=pid) continue;

	// Check that both generated muons are inside the single muon acceptance kinematic region
        const auto& pTD1 = tree->pTD1_gen()[iGen];
        const auto& etaD1 = tree->EtaD1_gen()[iGen];
        const bool mu1InAccep = (PD=="DIMUON" ? pPb::R8TeV::Y2016::triggerMuonAcceptance(pTD1, etaD1) : pPb::R8TeV::Y2016::muonAcceptance(pTD1, etaD1));
        const auto& pTD2 = tree->pTD2_gen()[iGen];
        const auto& etaD2 = tree->EtaD2_gen()[iGen];
        const bool mu2InAccep = (PD=="DIMUON" ? pPb::R8TeV::Y2016::triggerMuonAcceptance(pTD2, etaD2) : pPb::R8TeV::Y2016::muonAcceptance(pTD2, etaD2));
	const bool passMuAccep = mu1InAccep && mu2InAccep;

	// Check that candidate is within analysis kinematic range
	const auto& pT = tree->pT_gen()[iGen];
	const auto& rap = (col=="Pbp8Y16" ? -1.0 : 1.0) * tree->y_gen()[iGen];
	const bool passCandAccep = (fabs(rap)<2.4 && (pT>6.5 || fabs(rap)>1.4)) && passMuAccep;
	
        // Fill the VarMap with the kinematic information
	const auto nTrack = (sample.find("pPbGEN")==std::string::npos ? tree->Ntrkoffline() : 0);
	const auto rapCM = pPb::EtaLABtoCM(rap, (col=="pPb8Y16"));
	VarMap_t  varInfo =
	  {
	   { "Cand_Pt",     pT        },
	   { "Cand_Rap",    rap       },
	   { "Cand_AbsRap", fabs(rap) },
	   { "Cand_RapCM" , rapCM     },
	   { "NTrack",      nTrack    }
	  };
	
	if (sample.find("pPbGEN")!=std::string::npos) {
	  //
	  // Total Acceptance (Based on Generated muons)
	  //
	  TnPVec_t sfMC = {{"NoCorr", {1.0}}};
	  if (!fillEff1D(h1D, true, passCandAccep, sample, "pPb8Y16", "Acceptance", varInfo, sfMC, evtWeight)) { return; }
	  if (!fillEff1D(h1D, true, passCandAccep, sample, "Pbp8Y16", "Acceptance", varInfo, sfMC, evtWeight)) { return; }
	  //
	  continue;
	}
	
        // Initialize the boolean flags
	bool passRecoCandAccep  = false;
	bool passVertexProbCut  = false;
        bool passIdentification = false;
        bool passTrigger        = false;
	bool passDecayCut_85    = false;
	bool passDecayCut_90    = false;
	bool passDecayCut_95    = false;
	bool passDecayCut_99    = false;
	
        // Initialize the Tag-And-Probe scale factos
        TnPVec_t sfTnP = {};

        // Check that the generated candidate is within the analysis kinematic range
        if (passCandAccep) {
          // Find the reconstructed candidate matched to gen
          const short iReco = tree->RecIdx_gen()[iGen];
          if (iReco >= 0) {
            //
            // Candidate was matched to generated candidate
	    
            // Extract the kinematic information of reconstructed candidate
            const auto& cand_Pt    = tree->pT()[iReco];
            const auto& cand_Rap   = tree->y()[iReco];
	    const auto& cand_Eta   = tree->eta()[iReco];
	    const auto  cand_P     = cand_Pt*std::cosh(cand_Eta);
	    const auto  decayLen   = (tree->V3DDecayLength()[iReco]*tree->V3DCosPointingAngle()[iReco])*(3.0969/cand_P)*10.0;
            const auto& cand_PtD1  = tree->pTD1()[iReco];
            const auto& cand_EtaD1 = tree->EtaD1()[iReco];
            const auto& cand_PtD2  = tree->pTD2()[iReco];
            const auto& cand_EtaD2 = tree->EtaD2()[iReco];
	    
            // Determine the Tag-And-Probe scale factors
            sfTnP = getTnPScaleFactors(cand_PtD1, cand_EtaD1, cand_PtD2, cand_EtaD2, corrType, PD=="DIMUON");

	    // Check if the reconstructed candidate is within analysis kinematic range
	    const bool candMu1InAccep = (PD=="DIMUON" ? pPb::R8TeV::Y2016::triggerMuonAcceptance(cand_PtD1, cand_EtaD1) : pPb::R8TeV::Y2016::muonAcceptance(cand_PtD1, cand_EtaD1));
	    const bool candMu2InAccep = (PD=="DIMUON" ? pPb::R8TeV::Y2016::triggerMuonAcceptance(cand_PtD2, cand_EtaD2) : pPb::R8TeV::Y2016::muonAcceptance(cand_PtD2, cand_EtaD2));
	    passRecoCandAccep = (fabs(cand_Rap)<2.4 && (cand_Pt>6.5 || fabs(cand_Rap)>1.4)) && candMu1InAccep && candMu2InAccep;
	    
            // Check if the reconstructed muons pass muon ID
            passIdentification = tree->softCand(iReco);
	    
            // Check if the reconstructed candidate is matched to the trigger
	    const auto trigIdx = pPb::R8TeV::Y2016::HLTBitsFromPD(PD)[0];
	    passTrigger = tree->trigHLT()[trigIdx];
	    if (PD=="DIMUON") { passTrigger = passTrigger && tree->trigCand(trigIdx, iReco); }

	    // Check if the reconstructed candidate pass vertex probability cut
	    passVertexProbCut = tree->VtxProb()[iReco] > 0.0001;

	    // Check if the reconstructed candidate pass decay lenght cut
	    passDecayCut_85 = passRecoCandAccep && (decayLen < pPb::R8TeV::Y2016::decayLenCut(cand_Pt, cand_Rap, (PD=="DIMUON"), 0.85));
	    passDecayCut_90 = passRecoCandAccep && (decayLen < pPb::R8TeV::Y2016::decayLenCut(cand_Pt, cand_Rap, (PD=="DIMUON"), 0.90));
	    passDecayCut_95 = passRecoCandAccep && (decayLen < pPb::R8TeV::Y2016::decayLenCut(cand_Pt, cand_Rap, (PD=="DIMUON"), 0.95));
	    passDecayCut_99 = passRecoCandAccep && (decayLen < pPb::R8TeV::Y2016::decayLenCut(cand_Pt, cand_Rap, (PD=="DIMUON"), 0.99));
          }
        }
	//
	// Total Efficiency (Based on Generated muons)
	//
	const bool& isAnaCand = (passCandAccep && passRecoCandAccep && passIdentification && passTrigger && passVertexProbCut && passEventSelection);
	if (!fillEff1D(h1D, passCandAccep, isAnaCand, sample, col, "Efficiency_Total", varInfo, sfTnP, evtWeight)) { return; }
	//
	// Decay cut Efficiency (Based on Generated muons)
	//
	if (!fillEff1D(h1D, isAnaCand, isAnaCand && passDecayCut_85, sample, col, "Efficiency_DecayCut_85", varInfo, sfTnP, evtWeight)) { return; }
	if (!fillEff1D(h1D, isAnaCand, isAnaCand && passDecayCut_90, sample, col, "Efficiency_DecayCut_90", varInfo, sfTnP, evtWeight)) { return; }
	if (!fillEff1D(h1D, isAnaCand, isAnaCand && passDecayCut_95, sample, col, "Efficiency_DecayCut_95", varInfo, sfTnP, evtWeight)) { return; }
	if (!fillEff1D(h1D, isAnaCand, isAnaCand && passDecayCut_99, sample, col, "Efficiency_DecayCut_99", varInfo, sfTnP, evtWeight)) { return; }
      }
    }
  }
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Declare the efficiencies
  EffMap_t eff1D; // Stores the efficiency
  //
  // Load the efficiencies with the histograms
  if (!loadEff1D(eff1D, h1D)) { return; };
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Merge Efficiencies
  mergeEff(eff1D);
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Declare the uncertainties container
  Unc1DMap_t unc1D; // Stores the efficiency
  //
  // Calculate Uncertainties
  if (!getTnPUncertainties(unc1D, eff1D)) { return; };
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Define the output directory
  //
  // Store the Efficiencies
  //
  std::string outDir = mainDir + workDirName + "/" + PD + "/";
  gSystem->mkdir(outDir.c_str(), kTRUE);
  saveEff(outDir, eff1D, unc1D);
};


double getTnPScaleFactor(const double& pt, const double& eta, const std::string& cor, const int& i, const bool& incMuTrig)
{
  //
  // - TrkM: (tnp_weight_trkM_ppb)
  //   * idx = 0:  nominal
  //   * idx = 1..100: toy variations, stat. only
  //   * idx = -1: syst variation, "new_MAX", +1 sigma NOT YET IMPLEMENTED
  //   * idx = -2: syst variation, "new_MAX", -1 sigma NOT YET IMPLEMENTED
  //   * idx = -10: binned

  // - MuID: (tnp_weight_muid_ppb)
  //   * idx = -10: binned
  //   * idx = -11: binned, stat up
  //   * idx = -12: binned, stat down
  //   * idx = -13: binned, syst up
  //   * idx = -14: binned, syst down

  // - Trigger: (tnp_weight_trg_ppb)
  //   * idx = 0: nominal
  //   * idx = -1: syst variation,  +1 sigma
  //   * idx = -2: syst variation,  -1 sigma
  //   * idx = -11: stat variation,  +1 sigma
  //   * idx = -12: stat variation,  -1 sigma
  //
  double sf_TrkM = tnp_weight_trk_ppb ( pt , eta ,   0 );
  double sf_MuID = tnp_weight_muid_ppb( pt , eta , -10 );
  double sf_Trig = tnp_weight_trg_ppb ( pt , eta ,   0 );
  //
  if (cor=="NoCorr") { sf_TrkM = 1.0; sf_MuID = 1.0; sf_Trig = 1.0; }
  else if (cor=="TnP_Stat_TrkM"   ) { sf_TrkM = tnp_weight_trk_ppb ( pt , eta ,  i    ); }
  else if (cor=="TnP_Stat_MuID"   ) { sf_MuID = tnp_weight_muid_ppb( pt , eta , -10-i ); }
  else if (cor=="TnP_Stat_Trig"   ) { sf_Trig = tnp_weight_trg_ppb ( pt , eta , -10-i ); }
  else if (cor=="TnP_Syst_TrkM"   ) { sf_TrkM = tnp_weight_trk_ppb ( pt , eta ,  -i   ); }
  else if (cor=="TnP_Syst_MuID"   ) { sf_MuID = tnp_weight_muid_ppb( pt , eta ,   0   ); }
  else if (cor=="TnP_Syst_Trig"   ) { sf_Trig = tnp_weight_trg_ppb ( pt , eta ,   0   ); }
  else if (cor=="TnP_Syst_BinTrkM") { sf_TrkM = tnp_weight_trk_ppb ( pt , eta , -10   ); }
  if (!incMuTrig) { sf_Trig = 1.0; }
  //
  return ( sf_TrkM * sf_MuID * sf_Trig );
};


TnPVec_t getTnPScaleFactors(const double& ptD1, const double& etaD1, const double& ptD2, const double& etaD2, const CorrMap_t& corrType, const bool& incMuTrig)
{
  TnPVec_t sfTnP;
  for (const auto& cor : corrType) {
    sfTnP[cor.first].clear();
    for (uint i = 1; i <= cor.second; i++) {
      const auto sf_TnP_D1 = getTnPScaleFactor(ptD1, etaD1, cor.first, i, incMuTrig);
      const auto sf_TnP_D2 = getTnPScaleFactor(ptD2, etaD2, cor.first, i, incMuTrig);
      sfTnP.at(cor.first).push_back( sf_TnP_D1 * sf_TnP_D2 );
    }
  }
  return sfTnP;
};


bool getTnPUncertainties(Unc1DVec_t& unc, const EffVec_t& eff)
{
  if (eff.count("TnP_Nominal") == 0) { return true; }
  // Compute individual uncertainties
  const TEfficiency& nom = eff.at("TnP_Nominal")[0];
  const uint& nBin = nom.GetCopyTotalHisto()->GetNbinsX();
  for (const auto& co : eff) {
    if (co.first=="TnP_Nominal" || co.first.find("TnP_")==std::string::npos) continue;
    if (eff.count(co.first+"_")>0) continue;
    unc[co.first].ResizeTo(nBin);
    for (uint iBin = 1; iBin <= nBin; iBin++) {
      double uncVal = 0.0;
      if (co.second.size() == 100) {
        double sum = 0.0;
        for(uint i = 0; i < co.second.size(); i++) {
          const double diff = ( co.second[i].GetEfficiency(iBin) - nom.GetEfficiency(iBin) );
          sum += ( diff * diff );
        }
        uncVal = std::sqrt( sum / 100.0 );
      }
      else if (co.second.size() == 2) {
        uncVal = std::max(std::abs(co.second[0].GetEfficiency(iBin) - nom.GetEfficiency(iBin)) , std::abs(co.second[1].GetEfficiency(iBin) - nom.GetEfficiency(iBin)));
      }
      else if (co.second.size() == 1) {
        uncVal = std::abs(co.second[0].GetEfficiency(iBin) - nom.GetEfficiency(iBin));
      }
      else {
        std::cout << "[ERROR] Number of variations " << co.second.size() << " for " << co.first << " is wrong!" << std::endl; return false;
      }
      unc.at(co.first)[iBin-1] = uncVal;
    }
  }
  // Compute total statistical uncertainty
  unc["TnP_Stat"].ResizeTo(nBin);
  for (uint iBin = 0; iBin < nBin; iBin++) {
    double uncVal = 0.0;
    for (const auto& co : eff) {
      if (co.first.find("TnP_Stat")!=std::string::npos) {
        uncVal += std::pow( unc.at(co.first)[iBin] , 2.0 );
      }
    }
    unc.at("TnP_Stat")[iBin] = std::sqrt( uncVal );
  }
  // Compute total systematic uncertainty
  unc["TnP_Syst"].ResizeTo(nBin);
  for (uint iBin = 0; iBin < nBin; iBin++) {
    double uncVal = 0.0;
    for (const auto& co : eff) {
      if (co.first.find("TnP_Syst")!=std::string::npos) {
        uncVal += std::pow( unc.at(co.first)[iBin] , 2.0 );
      }
    }
    unc.at("TnP_Syst")[iBin] = std::sqrt( uncVal );
  }
  // Compute total TnP uncertainty
  unc["TnP_Tot"].ResizeTo(nBin);
  for (uint iBin = 0; iBin < nBin; iBin++) {
    double uncVal = 0.0;
    uncVal += std::pow( unc.at("TnP_Syst")[iBin] , 2.0 );
    uncVal += std::pow( unc.at("TnP_Stat")[iBin] , 2.0 );
    unc.at("TnP_Tot")[iBin] = std::sqrt( uncVal );
  }
  //
  return true;
};


bool getTnPUncertainties(Unc1DMap_t& unc, const EffMap_t& eff)
{
  for (const auto& s : eff) {
    for (const auto& c : s.second) {
      for (const auto& t : c.second) {
        for (const auto& b : t.second) {
	  if (!getTnPUncertainties(unc[s.first][c.first][t.first][b.first], b.second)) { return false; };
	}
      }
    }
  }
  return true;
};


void initEff1D(TH1DMap_t& h, const AnaVarMap_t& binMap, const CorrMap_t& corrType)
{
  for (const auto& sample : sampleType_.at("sample")) {
    for (const auto& col : COLL_) {
      for (const auto& effType : EFFTYPE_) {
	for (const auto& b : binMap) {
	  const auto& bin = b.first;
	  const auto& var = b.second;
	  const auto& hVar = std::get<0>(var);
	  const auto& varType = std::get<1>(var);
	  const auto& hBin = std::get<2>(var);
	  std::string hVarN;
	  for (const auto& v : bin.first) { hVarN += Form("%s_%.0f_%.0f_", v.name().c_str(), v.low()*10., v.high()*10.); }
	  hVarN += bin.second;
	  for (const auto& cor : corrType) {
	    if (effType=="Acceptance" && (cor.first.find("TnP_")!=std::string::npos)) continue;
	    for (uint i = 0; i < cor.second; i++) {
	      h[sample][col][effType][bin][cor.first].push_back( std::make_tuple(TH1D() , TH1D() , 1.0) );
	      auto& hist = h.at(sample).at(col).at(effType).at(bin).at(cor.first)[i];
	      if (varType == "FIX") {
		std::get<0>(hist) = TH1D("Passed", "Passed", hBin[0], hBin[1], hBin[2]);
		std::get<1>(hist) = TH1D("Total" , "Total" , hBin[0], hBin[1], hBin[2]);
	      }
	      else if (varType == "VAR") {
		std::get<0>(hist) = TH1D("Passed", "Passed", hBin.size()-1, hBin.data());
		std::get<1>(hist) = TH1D("Total" , "Total" , hBin.size()-1, hBin.data());
	      }
	      std::get<2>(hist) = 1.0;
	      //
	      const std::string& name = "h1D_" + sample +"_"+ col +"_"+ effType +"_"+ hVarN +"_"+ cor.first + ((cor.second>1) ? Form("_%d",(i+1)) : "");
	      std::get<0>(hist).SetName((name+"_Passed").c_str());
	      std::get<1>(hist).SetName((name+"_Total").c_str());
	      std::get<0>(hist).SetTitle(hVar.c_str());
	      std::get<1>(hist).SetTitle(hVar.c_str());
	      std::get<0>(hist).Sumw2();
	      std::get<1>(hist).Sumw2();
	    }
	  }
	}
      }
    }
  }
};


bool fillEff1D(TH1DVec_t& h, const bool& den_pass, const bool& num_pass, const double& xVar, const TnPVec_t& sfTnP, const double& evtWeight, const std::string& type)
{
  for (auto& cor : h) {
    for (uint i = 0; i < cor.second.size(); i++) {
      //
      bool found = false;
      //
      double sf = 1.0;
      if (sfTnP.count(cor.first)>0 && sfTnP.at(cor.first).size()>i) { sf = sfTnP.at(cor.first)[i]; found = true; }
      else if (sfTnP.size()==0 || sfTnP.count(cor.first)==0) { sf = 1.0; }
      else { std::cout << "[ERROR] Correction " << cor.first << " has invalid number of entries: " << cor.second.size() << "  " << sfTnP.at(cor.first).size() << " !" << std::endl; return false; }
      if ((cor.first.find("TnP_")!=std::string::npos) && num_pass && sfTnP.size()==0) { std::cout << "[ERROR] TnP scale factor vector is empty!" << std::endl; return false; }
      //
      if (found==false && sfTnP.size()>0) {
        std::cout << "[ERROR] Correction " << cor.first << " was not found!" << std::endl; return false;
      }
      //
      // Fill the passing histogram (numerator)
      if (num_pass) {
	if      (cor.first=="NoCorr") { std::get<0>(cor.second[i]).Fill(xVar , 1.0          ); }
	else if (type=="Acceptance" ) { std::get<0>(cor.second[i]).Fill(xVar , evtWeight    ); }
	else                          { std::get<0>(cor.second[i]).Fill(xVar , evtWeight*sf ); }
      }
      // Fill the total histogram (denominator)
      if (den_pass) {
	const bool applySF = (type.rfind("Efficiency_DecayCut",0)==0);
	if      (cor.first=="NoCorr") { std::get<1>(cor.second[i]).Fill(xVar , 1.0);           }
	else if (applySF)             { std::get<1>(cor.second[i]).Fill(xVar , evtWeight*sf ); }
	else                          { std::get<1>(cor.second[i]).Fill(xVar , evtWeight );    }
      }
    }
  }
  return true;
};


bool fillEff1D(TH1DMap_t& h, const bool& den_pass, const bool& num_pass, const std::string& sample, const std::string& col, const std::string& type,
	       const VarMap_t& var, const TnPVec_t& sfTnP, const double& evtWeight)
{
  if (sample.find(col)!=std::string::npos) { std::cout << "[ERROR] Sample name " << sample << " has wrong format!" << std::endl; return false; }
  for (auto& s : h) {
    //
    if (sample.rfind(s.first,0)==0) continue;
    //
    for (auto& c : s.second) {
      //
      if ((c.first=="PA8Y16" && col!="Pbp8Y16") || (c.first!="PA8Y16" && c.first!=col)) continue;
      //
      if (c.second.count(type)==0) { return true; }
      //
      for (auto& b : c.second.at(type)) {
	//
	const auto& bin1 = b.first.first.getbin(0);
	const auto& bin2 = b.first.first.getbin(1);
	const auto& bin3 = b.first.second;
	if (var.count(bin1.name())==0) { std::cout << "[ERROR] Variable " << bin1.name() << " was not loaded properly" << std::endl; return false; }
	if (var.count(bin2.name())==0) { std::cout << "[ERROR] Variable " << bin2.name() << " was not loaded properly" << std::endl; return false; }
	const auto& val1 = var.at(bin1.name());
	const auto& val2 = var.at(bin2.name());
	const auto& xVarName = bin3.substr(0, bin3.rfind("_"));
	if (var.count(xVarName)==0) { std::cout << "[ERROR] Variable " << xVarName << " was not loaded properly" << std::endl; return false; }
	const auto& xVar = var.at(xVarName);
	//
	bool inBin = true;
	if (type!="Acceptance" || bin1.name()!="NTrack") { inBin = (val1>=bin1.low() && val1<bin1.high()); }
	if (type!="Acceptance" || bin2.name()!="NTrack") { inBin = inBin && (val2>=bin2.low() && val2<bin2.high()); }
	if (inBin) { // Don't include values outside of range
	  // Fill histograms
	  if (!fillEff1D(b.second, den_pass, num_pass, xVar, sfTnP, evtWeight, type)) { return false; }
	}
      }
    }
  }
  return true;
};


bool loadEff1D(EffMap_t& ef, const TH1DMap_t& h)
{
  for (const auto& s : h) {
    for (const auto& c : s.second) {
      for (const auto& t : c.second) {
	for (const auto& b : t.second) {
	  for (const auto& co : b.second) {
	    for (uint i = 0; i < co.second.size(); i++) {
	      //
	      const auto& hPassed = std::get<0>(co.second[i]);
	      const auto& hTotal  = std::get<1>(co.second[i]);
	      const auto& weight  = std::get<2>(co.second[i]);
	      std::string hVarN;
	      for (const auto& v : b.first.first) { hVarN += Form("%s_%.0f_%.0f_", v.name().c_str(), v.low()*10., v.high()*10.); }
	      hVarN += b.first.second;
	      //
	      const std::string& name = "eff1D_" + s.first +"_"+ c.first +"_"+ t.first +"_"+ hVarN +"_"+ co.first + ((co.second.size()>1) ? Form("_%d",(i+1)) : "");
	      //
	      if ( TEfficiency::CheckConsistency(hPassed, hTotal,"w") ) {
		ef[s.first][c.first][t.first][b.first][co.first].push_back( TEfficiency(hPassed , hTotal) );
		auto& eff = ef.at(s.first).at(c.first).at(t.first).at(b.first).at(co.first)[i];
		eff.SetName(name.c_str());
		// Set Global Weight
                eff.SetWeight(weight);
		if ( checkWeights(hPassed, hTotal) ) {
		  eff.SetStatisticOption(TEfficiency::kBJeffrey);
		}
	      }
	      else {
		std::cout << "[ERROR] 1D Histograms for " << name << " are not consistent!" << std::endl;
		for(Int_t i = 0; i < (hPassed.GetNbinsX()+2); ++i) {
		  if(hPassed.GetBinContent(i) > hTotal.GetBinContent(i)) {
		    std::cout << "[ERROR] Bin " << i << " has Pass: " << hPassed.GetBinContent(i) << " Total: " << hTotal.GetBinContent(i) << std::endl;
		  }
		}
		ef[s.first][c.first][t.first][b.first][co.first].push_back( TEfficiency() );
		auto& eff = ef.at(s.first).at(c.first).at(t.first).at(b.first).at(co.first)[i];
                eff.SetName(name.c_str());
                eff.SetTotalHistogram(hTotal, "f");
                eff.SetPassedHistogram(hPassed, "f");
		// Set Global Weight
                eff.SetWeight(weight);
		if ( checkWeights(hPassed, hTotal) ) {
		  eff.SetStatisticOption(TEfficiency::kBJeffrey);
		}
	      }
	    }
	  }
	}
      }
    }
  }
  return true;
};
  

void mergeEff(EffMap_t& ef)
{
  // Merge pPb and Pbp (inverted) -> PA
  if(std::find(COLL_.begin(), COLL_.end(), "PA") != COLL_.end()) {
    for (auto& s : ef) {
      if (s.second.count("PA8Y16")>0 && s.second.count("pPb8Y16")>0) {
	for (auto& t : s.second.at("PA8Y16")) {
	  for (auto& b : t.second) {
	    for (auto& co : b.second) {
	      for (uint i = 0; i < co.second.size(); i++) {
		// Just add the pPb, no need to combine
		const auto& eff_pPb = s.second.at("pPb8Y16").at(t.first).at(b.first).at(co.first)[i];
		const auto& eff_Pbp = s.second.at("Pbp8Y16").at(t.first).at(b.first).at(co.first)[i]; // Only used for the weight
		// Passed Histogram
		auto hPassed = *dynamic_cast<const TH1D*>(eff_pPb.GetPassedHistogram());
		hPassed.Add(eff_pPb.GetPassedHistogram(), co.second[i].GetPassedHistogram(), eff_pPb.GetWeight(), eff_Pbp.GetWeight());
		co.second[i].SetPassedHistogram(hPassed, "f");
		// Total Histogram
		auto hTotal = *dynamic_cast<const TH1D*>(eff_pPb.GetTotalHistogram());
		hTotal.Add(eff_pPb.GetTotalHistogram(), co.second[i].GetTotalHistogram(), eff_pPb.GetWeight(), eff_Pbp.GetWeight());
		co.second[i].SetTotalHistogram(hTotal, "f");
		// Set the statistics in case of weighted histograms
		if ( checkWeights(hPassed, hTotal) ) { co.second[i].SetStatisticOption(TEfficiency::kBJeffrey); }
	      }
	    }
	  }
	}
      }
      else {
	std::cout << "[WARNING] Can't merge pPb and Pbp in " << s.first << " ! "<< std::endl;
	if (s.second.count("PA8Y16")>0) { s.second.erase("PA8Y16"); }
      }
    }
  }
};


void writeEff(TFile& file, const EffMap_t& eff, const Unc1DMap_t& unc, const std::string& sample, const std::string& mainDirName)
{
  if (eff.size()>0 && eff.count(sample)>0) {
    auto mainDir = file.mkdir(mainDirName.c_str());
    mainDir->cd();
    for (auto& c : eff.at(sample)) {
      auto colDir = mainDir->GetDirectory(c.first.c_str());
      if (colDir==NULL) { colDir = mainDir->mkdir(c.first.c_str()); }
      colDir->cd();
      for (auto& t : c.second) {
	auto typeDir = colDir->GetDirectory(t.first.c_str());
	if (typeDir==NULL) { typeDir = colDir->mkdir(t.first.c_str()); }
	typeDir->cd();
	for (auto& b : t.second) {
	  const std::map<std::string , TVectorD>& u = unc.at(sample).at(c.first).at(t.first).at(b.first);
	  std::string binN;
	  for (const auto& v : b.first.first) { binN += Form("%s_%.0f_%.0f_", v.name().c_str(), v.low()*10., v.high()*10.); }
	  binN += b.first.second;
	  const auto& name = "unc1D_" + sample +"_"+ c.first +"_"+ t.first +"_"+ binN +"_";
	  auto binDir = typeDir->GetDirectory(binN.c_str());
	  if (binDir==NULL) { binDir = typeDir->mkdir(binN.c_str()); }
	  binDir->cd();
	  for (auto& co : b.second) {
	    auto corDir = binDir->GetDirectory(co.first.c_str());
	    if (corDir==NULL) { corDir = binDir->mkdir(co.first.c_str()); }
	    corDir->cd();
	    for (uint i = 0; i < co.second.size(); i++) {
	      std::cout << co.second[i].GetName() << "  " << co.second[i].GetTotalHistogram()->GetEntries() << "  " << co.second[i].GetPassedHistogram()->GetEntries() << std::endl;
	      co.second[i].Write(clStr(co.second[i].GetName()));
	    }
	    if (u.count(co.first)>0) {
	      if (co.first.find("TnP_S")!=std::string::npos) { u.at(co.first).Write((name + co.first).c_str()); }
	    }
	    binDir->cd();
	  }
	  if (u.count("TnP_Tot") > 0) { 
	    u.at("TnP_Stat").Write((name + "TnP_Stat").c_str());
	    u.at("TnP_Syst").Write((name + "TnP_Syst").c_str());
	    u.at("TnP_Tot").Write((name + "TnP_Tot").c_str());
	  }
	  typeDir->cd();
	}
	colDir->cd();
      }
      mainDir->cd();
    }
  }
};


void saveEff(const std::string& outDir, const EffMap_t& eff1D, const Unc1DMap_t& unc1D)
{
  for (const auto& s : eff1D) {
    // Step 0: Create the output file
    const std::string fileName = "effContainer_"+s.first+".root";
    TFile file((outDir + fileName).c_str(), "RECREATE");
    file.cd();
    // Step 1: Store all the 1D Efficiency objects
    writeEff(file, eff1D, unc1D, s.first, "Eff1D");
    // Step 3: Write and Close the file
    file.Write();
    file.Close("R");
    std::cout << "[INFO] Information store in file " << outDir << fileName << std::endl;
  }
};


void setGlobalWeight(TH1DMap_t& h, const double& weight, const std::string& sample, const std::string& col)
{
  if (h.count(sample) && h.at(sample).count(col)>0) {
    for (auto& t : h.at(sample).at(col)) {
      for (auto& b : t.second) {
	for (auto& co : b.second) {
	  for (auto& p : co.second) {
	    std::get<2>(p) = weight;
	  }
	}
      }
    }
  }
};


bool checkWeights(const TH1& pass , const TH1& total)
{
  if (pass.GetSumw2N() == 0 && total.GetSumw2N() == 0) return false;
 
  // check also that the total sum of weight and weight squares are consistent
  Double_t statpass[TH1::kNstat];
  Double_t stattotal[TH1::kNstat];
 
  pass.GetStats(statpass);
  total.GetStats(stattotal);
 
  double tolerance = (total.IsA() == TH1F::Class() ) ? 1.E-5 : 1.E-12; 
    
  //require: sum of weights == sum of weights^2
  if(!TMath::AreEqualRel(statpass[0],statpass[1],tolerance) ||
     !TMath::AreEqualRel(stattotal[0],stattotal[1],tolerance) ) {
    return true;
  }

  // histograms are not weighted 
  return false;  
};


const char* clStr(const std::string& in)
{
  std::string out = in;
  while (out.find("_copy")!=std::string::npos) { out.erase(out.find("_copy"), 5); }
  return out.c_str();
};


using KeyPVec_t  =  std::vector< TKey* >;


KeyPVec_t getK(TList* list)
{
  TIter iter(list);
  KeyPVec_t out;
  while (TKey* key = (TKey*)iter()) { out.push_back(key); }
  return out;
};


bool getEffObjectsFromFile(EffMap_t& ef, Unc1DMap_t& un, const std::string& filePath)
{
  // Open the input file
  TFile file(filePath.c_str(), "READ");
  if (file.IsZombie()) { std::cout << "[ERROR] The input file " << filePath << " was not found!" << std::endl; return false; }
  file.cd();
  //
  // Get sample
  std::string sample;
  auto tmp = filePath.substr(filePath.rfind("_MC")+1);
  sample = tmp.substr(0, tmp.rfind(".root"));
  //
  //gain time, do not add the objects in the list in memory
  Bool_t status = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  //
  const std::string mainDirName = gDirectory->GetListOfKeys()->First()->GetName();
  TDirectory* mainDir = file.GetDirectory(mainDirName.c_str());
  if (mainDir==NULL) { std::cout << "[ERROR] Dir: " << mainDirName << " was not found!" << std::endl; return false; }
  mainDir->cd();
  // loop over all keys in this directory
  for (const auto& c : getK(mainDir->GetListOfKeys())) {
    if (c->ReadObj()->IsA()->InheritsFrom(TDirectory::Class()) == false) continue;
    auto colDir = mainDir->GetDirectory(c->GetName());
    colDir->cd();
    for (const auto& t : getK(colDir->GetListOfKeys())) {
      if (t->ReadObj()->IsA()->InheritsFrom(TDirectoryFile::Class()) == false) continue;
      auto typeDir = colDir->GetDirectory(t->GetName());
      typeDir->cd();
      for (const auto& b : getK(typeDir->GetListOfKeys())) {
        if (b->ReadObj()->IsA()->InheritsFrom(TDirectoryFile::Class()) == false) continue;
        auto binDir = typeDir->GetDirectory(b->GetName());
        binDir->cd();
	const auto& binN = std::string(b->GetName());
	//
	const auto varA1 = binN.substr(0, binN.rfind("NTrack")-1);
	const auto varB1 = varA1.substr(varA1.rfind("_")+1);
	const auto varC1 = varA1.substr(0, varA1.rfind("_"));
	const auto varD1 = varC1.substr(varC1.rfind("_")+1);
	const auto varE1 = varC1.substr(0, varC1.rfind("_"));
	const auto& var1 = BinF_t({varE1, float(std::stof(varD1)/10.), float(std::stof(varB1)/10.)});
	//
	const auto tmp1 = binN.substr(binN.rfind("NTrack"));
	//
	const auto varA2 = tmp1.substr(0, tmp1.rfind("Cand")-1);
	const auto varB2 = varA2.substr(varA2.rfind("_")+1);
	const auto varC2 = varA2.substr(0, varA2.rfind("_"));
	const auto varD2 = varC2.substr(varC2.rfind("_")+1);
	const auto varE2 = varC2.substr(0, varC2.rfind("_"));
	const auto& var2 = BinF_t({varE2, float(std::stof(varD2)/10.), float(std::stof(varB2)/10.)});
	//
	const auto var3 = tmp1.substr(tmp1.rfind("Cand"));
	//
	const auto anaBin = AnaBin_t({var1, var2});
	const auto bin = AnaBinPair_t({anaBin, var3});
	//
        for (const auto& co : getK(binDir->GetListOfKeys())) {
          if (co->ReadObj()->IsA()->InheritsFrom(TDirectoryFile::Class())) {
	    auto corDir = binDir->GetDirectory(co->GetName());
	    corDir->cd();
	    for (const auto& key : getK(corDir->GetListOfKeys())) {
	      if (key->ReadObj()->IsA()->InheritsFrom(TEfficiency::Class())) {
		const auto& effP = static_cast<TEfficiency*>(key->ReadObj());
		if (effP) {
		  ef[sample][c->GetName()][t->GetName()][bin][co->GetName()].push_back(*effP);
		  auto& eff = ef.at(sample).at(c->GetName()).at(t->GetName()).at(bin).at(co->GetName());
		  eff.back().SetName(key->GetName());
		}
	      }
	      if (key->ReadObj()->IsA()->InheritsFrom(TVectorD::Class())) {
		const auto& uncP = static_cast<TVectorD*>(key->ReadObj());
		if (uncP) {
		  un[sample][c->GetName()][t->GetName()][bin][co->GetName()].ResizeTo(uncP->GetNrows());
		  auto& unc = un.at(sample).at(c->GetName()).at(t->GetName()).at(bin).at(co->GetName());
		  unc = *uncP;
		}
	      }
	    }
	    binDir->cd();
	  }
	  if (co->ReadObj()->IsA()->InheritsFrom(TVectorD::Class())) {
	    std::string corrName = "";
	    if (std::string(co->GetName()).find("TnP_Stat")!=std::string::npos) { corrName = "TnP_Stat"; }
	    if (std::string(co->GetName()).find("TnP_Syst")!=std::string::npos) { corrName = "TnP_Syst"; }
	    if (std::string(co->GetName()).find("TnP_Tot" )!=std::string::npos) { corrName = "TnP_Tot";  }
	    if (std::string(co->GetName()).find("MC_Syst" )!=std::string::npos) { corrName = "MC_Syst";  }
	    if (corrName!="") {
	      const auto& uncP = static_cast<TVectorD*>(co->ReadObj());
	      if (uncP) {
		un[sample][c->GetName()][t->GetName()][bin][co->GetName()].ResizeTo(uncP->GetNrows());
		auto& unc = un.at(sample).at(c->GetName()).at(t->GetName()).at(bin).at(co->GetName());

	      }
	    }
	  }
	}
	typeDir->cd();
      }
      colDir->cd();
    }
    mainDir->cd();
  }
  // Close File
  file.Close("R");
  // return
  return true;
};
