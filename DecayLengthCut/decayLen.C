#if !defined(__CINT__) || defined(__MAKECINT__)
// Auxiliary Headers
#include "../Utilities/Ntuple/VertexCompositeTree.h"
#include "../Utilities/dataUtils.h"
#include "../Efficiency/util.h"
#include "../Utilities/CMS/tdrstyle.C"
#include "../Utilities/CMS/CMS_lumi.C"
// ROOT headers
#include "TSystem.h"
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TEfficiency.h"
#include "TMatrixD.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TGaxis.h"
// c++ headers
#include <dirent.h>
#include <iostream>
#include <map>
#include <vector>
#include <tuple>
#include <string>

#endif


// ------------------ TYPE -------------------------------
using BinHistMap_t    = std::map< BinF_t, std::map< BinF_t, std::map< std::string , std::map< BinF_t , std::tuple< TH1D , int , double > > > > >;
using BinHistDiMap_t  = std::map< std::string, BinHistMap_t >;
using BinHistTriMap_t = std::map< std::string, BinHistDiMap_t >;
using BinHistQuadMap_t = std::map< std::string, BinHistTriMap_t >;
using BinValMap_t     = std::map< BinF_t, std::map< BinF_t, std::map< std::string , std::map< double , TMatrixD > > > >;
using BinValDiMap_t   = std::map< std::string, BinValMap_t >;
using BinValTriMap_t  = std::map< std::string, BinValDiMap_t >;
using BinValQuadMap_t  = std::map< std::string, BinValTriMap_t >;


// ------------------ FUNCTION -------------------------------
bool fillHist      ( BinHistTriMap_t& histMap   , const StringVector_t& inputFiles , const std::string& sample );
bool storeHist     ( const BinHistQuadMap_t& histMap , const std::string& outDir );
bool extractHist   ( BinHistQuadMap_t& histMap  , const std::string& outDir      );
void getThreshold  ( BinValQuadMap_t& thrMap    , const BinHistQuadMap_t& histMap );
void findThreshold ( double& val , double& unc_low , double& unc_high , const TH1D& hist , const double& eff_thr);
bool storeThreshold   ( const BinValQuadMap_t& thrMap , const std::string& outDir );
bool extractThreshold ( BinValQuadMap_t& thrMap       , const std::string& outDir );
void plotThreshold    ( const BinValQuadMap_t& thrMap , const std::string& outDir , const bool& altFunc );
  
  
// ------------------ GLOBAL -------------------------------
//
// Efficiency thresholds
const std::vector<double> EFF_THR_({ 0.85, 0.90, 0.95, 0.99 });
//
// Collision System
const std::vector<std::string> COLL_ = { "PA8Y16" };
//
// Datasets
const std::vector<std::string> PD_ = { "DIMUON", "MINBIAS" };
//
// Binning
BinMapMap_t BINMAP =
  {
   {
    {"Cand_Pt", 6.5, 50.0},
    {
     {
      {"Cand_AbsRap", 1.4, 2.4},
      {
       {"NTrack", "VAR",   {15.,30.,40.,50.,60.,70.,75.,80.,85.,90.,95.,100.,105.,110.,115.,120.,130.,140.,150.,170.,200.,250.}}
      }
     },
     {
      {"NTrack", 15., 250.},
      {
       {"Cand_Rap", "VAR", {-2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4}},
      }
     }
    }
   },
   {
    {"Cand_Pt", 3.0, 50.0},
    {
     {
      {"Cand_AbsRap", 0.0, 1.4},
      {
       {"NTrack", "VAR",   {15.,30.,40.,50.,60.,70.,75.,80.,85.,90.,95.,100.,105.,110.,115.,120.,130.,140.,150.,170.,200.,250.}}
      }
     }
    }
   },
   {
    {"NTrack", 15., 250.},
    {
     {
      {"Cand_AbsRap", 0.0, 0.4},
      {
       {"Cand_Pt",  "VAR", {6.5,8,9,10,11,12,13,15,17,19,21,25,30,35,40,50}},
      }
     },
     {
      {"Cand_AbsRap", 0.4, 1.2},
      {
       {"Cand_Pt",  "VAR", {6.5,8,9,10,11,12,13,15,17,19,21,25,30,35,40,50}},
      }
     },
     {
      {"Cand_AbsRap", 1.2, 1.6},
      {
       {"Cand_Pt",  "VAR", {6.5,8,9,10,11,12,13,15,17,19,21,25,30,35,40,50}},
      }
     },
     {
      {"Cand_AbsRap", 1.6, 2.4},
      {
       {"Cand_Pt",  "VAR", {0,1,2,3,4,5,6.5,8,9,10,11,12,13,15,17,19,21,25,30,35,40,50}},
      }
     }
    }
   }
  };
//
// Input Files for analysis
const std::string path_MC = "/Users/andre/DiMuonAnalysis2019/Tree";
const std::map< std::string , std::string > inputFileMap_ =
  {
   {"MC_JPsiPR_pPb"    , Form("%s/%s", path_MC.c_str(), "VertexCompositeTree_JPsiToMuMu_pPb-Bst_pPb816Summer16_DiMuMC.root")   },
   {"MC_JPsiPR_Pbp"    , Form("%s/%s", path_MC.c_str(), "VertexCompositeTree_JPsiToMuMu_Pbp-Bst_pPb816Summer16_DiMuMC.root")   },
   {"MC_Psi2SPR_pPb"   , Form("%s/%s", path_MC.c_str(), "VertexCompositeTree_Psi2SToMuMu_pPb-Bst_pPb816Summer16_DiMuMC.root")  },
   {"MC_Psi2SPR_Pbp"   , Form("%s/%s", path_MC.c_str(), "VertexCompositeTree_Psi2SToMuMu_Pbp-Bst_pPb816Summer16_DiMuMC.root")  }
  };
std::map< std::string , std::vector< std::string > > inputSamples_;


void decayLen(const bool& altFunc=false)
{
  //
  std::cout << "[INFO] Starting to derive the decay length thresholds" << std::endl;
  //
  // Change the working directory
  const std::string CWD = getcwd(NULL, 0);
  const std::string mainDir = Form("%s/Output/", CWD.c_str());
  gSystem->mkdir(mainDir.c_str(), kTRUE);
  gSystem->ChangeDirectory(mainDir.c_str());
  //
  // Create list of samples
  inputSamples_ = {};
  for (const auto& f : inputFileMap_) {
    std::string tmp = f.first; tmp = tmp.substr(0, tmp.find_last_of("_"));
    inputSamples_[tmp].push_back(f.second);
  }
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Initialize the thresholds
  BinValQuadMap_t thrMap;
  for (const auto& s : inputSamples_) {
    for (const auto& col : COLL_) {
      for (const auto& pd : PD_) {
	const auto& smp = s.first;
	for (const auto& v1: BINMAP) {
	  for (const auto& v2: v1.second) {
	    for (const auto& v3: v2.second) {
	      const auto& var1 = BinF_t(v1.first);
	      const auto& var2 = BinF_t(v2.first);
	      const auto& var3N = std::get<0>(v3);
	      for (const auto& effThr : EFF_THR_) {
		thrMap[smp][col][pd][var1][var2][var3N][effThr] = TMatrixD();
		auto& thr = thrMap.at(smp).at(col).at(pd).at(var1).at(var2).at(var3N).at(effThr);
		thr.ResizeTo(std::get<2>(v3).size()-1, 6);
	      }
	    }
	  }
	}
      }
    }
  }
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Extract the thresholds
  //
  //
  const std::string& outTHRDir = mainDir + "Thresholds";
  const auto redoTHR = !extractThreshold(thrMap, outTHRDir);
  //
  if (redoTHR) {
    //
    // ------------------------------------------------------------------------------------------------------------------------
    //
    // Initialize the histograms
    BinHistQuadMap_t histMap;
    for (const auto& s : inputSamples_) {
      for (const auto& col : COLL_) {
	for (const auto& pd : PD_) {
	  const auto& smp = s.first;
	  for (const auto& v1: BINMAP) {
	    for (const auto& v2: v1.second) {
	      for (const auto& v3: v2.second) {
		const auto& var1 = BinF_t(v1.first);
		const auto& var2 = BinF_t(v2.first);
		const auto& var3N = std::get<0>(v3);
		const auto& var3V = std::get<2>(v3);
		for (uint i=1; i<var3V.size(); i++) {
		  const auto& var3 = BinF_t({var3N, float(var3V[i-1]), float(var3V[i])});
		  const std::string name = Form("hDLC_%s_%s_%s_%s_%.0f_%.0f_%s_%.0f_%.0f_%s_%.0f_%.0f",
						smp.c_str(), col.c_str(), pd.c_str(),
						var1.name().c_str(), var1.low()*10., var1.high()*10.,
						var2.name().c_str(), var2.low()*10., var2.high()*10.,
						var3.name().c_str(), var3.low()*10., var3.high()*10.);
		  histMap[smp][col][pd][var1][var2][var3N][var3] = std::make_tuple(TH1D(name.c_str(), "", 60000, -4.0, 8.0), 0, 0.0);
		  auto& hist = std::get<0>(histMap.at(smp).at(col).at(pd).at(var1).at(var2).at(var3N).at(var3));
		  hist.Sumw2();
		}
	      }
	    }
	  }
	}
      }
    }
    //
    // ------------------------------------------------------------------------------------------------------------------------
    //
    // Extract the histograms
    //
    //
    const std::string& outHISTDir = mainDir + "Histograms";
    const auto redoHist = !extractHist(histMap, outHISTDir);
    //
    if (redoHist) {
      //
      // ------------------------------------------------------------------------------------------------------------------------
      //
      // Fill the histograms
      //
      //
      for (auto& s : histMap) {
	const auto& smp = s.first;
	const auto& inputFiles = inputSamples_.at(smp);
	auto& histM = s.second;
	if (!fillHist(histM, inputFiles, smp)) { return; }
      }
      std::cout << "[INFO] Completed to fill the histograms" << std::endl;
      //
      // ------------------------------------------------------------------------------------------------------------------------
      //
      // Store the histograms
      //
      //
      if (!storeHist(histMap , outHISTDir)) { return; }
    }
    //
    // ------------------------------------------------------------------------------------------------------------------------
    //
    // Derive the decay length thresholds
    //
    //
    getThreshold(thrMap, histMap);
    //
    // ------------------------------------------------------------------------------------------------------------------------
    //
    // Store the decay length thresholds
    //
    //
    if (!storeThreshold(thrMap , outTHRDir)) { return; }
  }
  //
  plotThreshold(thrMap, outTHRDir, altFunc);
  //
};


bool fillHist(BinHistTriMap_t& histMap, const StringVector_t& inputFiles, const std::string& sample)
{
  // Loop over input files
  for(const auto& file : inputFiles) {
    
    VertexCompositeTree tree;
    if (!tree.GetTree(file, "dimucontana_mc")) { std::cout << "Invalid tree for: " << file << "!" << std::endl; return false; }
    
    StringSet_t objS;
    if (sample.rfind("MC_JPsi",0)==0) { objS.insert("JPsi"); }
    else if (sample.rfind("MC_Psi2S",0)==0) { objS.insert("Psi2S"); }
      
    // Determine the collision system of the sample
    std::string col = "";
    if (file.find("_Pbp-")!=std::string::npos) col = "Pbp8Y16"; // for Pbp
    if (file.find("_pPb-")!=std::string::npos) col = "pPb8Y16"; // for pPb
    if (col=="") { std::cout << "[ERROR] Could not determine the collision system in the sample" << std::endl; return false; }
    
    // Loop over the events
    int treeIdx = -1;
    const auto nentries = tree.GetEntries();
    std::cout << "[INFO] Starting to process " << nentries << " nentries" << std::endl;
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
      
      // Get the entry in the trees
      if (tree.GetEntry(jentry)<0) { std::cout << "Invalid entry for: " << file << "!" << std::endl; return 0.0; }
      
      if (tree.GetTreeNumber()!=treeIdx) {
        treeIdx = tree.GetTreeNumber();
        std::cout << "[INFO] Processing " << sample << " using file: " << file << std::endl;
      }
      loadBar(jentry, nentries);

      // Loop over candidates
      for(uint iReco=0; iReco<tree.candSize(); iReco++) {
	
	// Check that candidate is matched to gen
	if (tree.matchGEN()[iReco]==false) continue;
	// Check the PID of the matched gen particle
	int pid = 443;
	if (sample.rfind("MC_Psi2S",0)==0) { pid = 100443; }
	if (fabs(tree.idmom_reco()[iReco])!=pid) continue;

	const auto pT = tree.pT()[iReco];
	const auto eta = tree.eta()[iReco];
        const auto p = pT*std::cosh(eta);
        const auto decayLen = (tree.V3DDecayLength()[iReco]*tree.V3DCosPointingAngle()[iReco])*(3.0969/p)*10.0;

	// Extract variables
	std::map<std::string , double> varInfo =
	  {
	   { "Cand_Rap"  , tree.y()[iReco] },
	   { "Cand_AbsRap"  , std::abs(tree.y()[iReco]) },
	   { "Cand_Pt"   , tree.pT()[iReco] },
	   { "NTrack"    , tree.NTracks()[iReco] },
	   { "Cand_DLen" , decayLen }
	  };

	// Fill the histogram
	for (auto& c : histMap) {
	  if (c.first!="PA8Y16" && c.first!=col) continue;
	  for (auto& pd : c.second) {
	    // Apply analysis cuts as in data
	    const auto& COL = pd.first;
	    const auto& PD = pd.first;
	    // Apply kinematic range
	    if (PD=="DIMUON" && tree.pT()[iReco]<3.0) continue;
	    if (std::abs(tree.y()[iReco])<1.4 && tree.pT()[iReco]<6.5) continue;
	    // Derive the luminosity weight
	    const auto& weight = pPb::R8TeV::Y2016::LumiWeightFromPD(PD, col, sample);
	    if (ANA::analysisSelection(tree, iReco, PD, col, objS, false)) {
	      for (auto& v1 : pd.second) {
		const auto& bin1 = v1.first;
		const auto& val1 = varInfo.at(bin1.name());
		if (val1<bin1.low() || val1>=bin1.high()) continue;
		for (auto& v2 : v1.second) {
		  const auto& bin2 = v2.first;
		  const auto& val2 = varInfo.at(bin2.name());
		  if (val2<bin2.low() || val2>=bin2.high()) continue;
		  for (auto& v3 : v2.second) {
		    for (auto& b3 : v3.second) {
		      const auto& bin3 = b3.first;
		      const auto& val3 = varInfo.at(v3.first);
		      if (val3<bin3.low() || val3>=bin3.high()) continue;
		      auto& hist = std::get<0>(b3.second);
		      auto& n = std::get<1>(b3.second);
		      auto& m = std::get<2>(b3.second);
		      hist.Fill(varInfo.at("Cand_DLen"), weight);
		      n += 1;
		      m += val3;
		    }
		  }
		}
	      }
	    }
	  }
	}
	//
      }
    }
  }
  std::cout << "[INFO] Filling histograms done!" << std::endl;
  return true;
};


bool storeHist(const BinHistQuadMap_t& histMap, const std::string& outDir)
{
  std::cout << "[INFO] Storing the histograms" << std::endl;
  bool notStored = false;
  for (const auto& s : histMap) {
    for (const auto& c : s.second) {
      // Step 0: Create the output file
      const std::string fileDir = outDir + "/" + s.first + "/" + c.first + "/";
      makeDir(fileDir);
      for (const auto& pd : c.second) {
	const std::string fileName = fileDir + "effContainer_"+s.first+"_"+c.first+"_"+pd.first+".root";
	TFile file(fileName.c_str(), "RECREATE");
	file.cd();
	if (file.IsOpen() && !file.IsZombie()) {
	  // Step 1: Store all the 1D Efficiency objects
	  for (const auto& v1 : pd.second) {
	    for (const auto& v2 : v1.second) {
	      for (const auto& v3 : v2.second) {
		for (const auto& b3 : v3.second) {
		  auto& hist = std::get<0>(b3.second);
		  auto& n = std::get<1>(b3.second);
		  auto& m = std::get<2>(b3.second);
		  const_cast<TH1D*>(&hist)->SetTitle(Form("%.3f_%d", m, n));
		  hist.Write(hist.GetName());
		}
	      }
	    }
	  }
	  // Step 3: Write and Close the file
	  file.Write();
	  file.Close("R");
	}
	else { std::cout << "[INFO] File " << fileName << " could not be created" << std::endl; notStored = true; }
	std::cout << "[INFO] Information store in file " << outDir << fileName << std::endl;
      }
    }
  }
  std::cout << "[INFO] Storing histograms done!" << std::endl;
  return !notStored;
};


bool extractHist(BinHistQuadMap_t& histMap, const std::string& outDir)
{
  std::cout << "[INFO] Extracting the histograms" << std::endl;
  bool redoHist = false;
  for (auto& s : histMap) {
    for (auto& c : s.second) {
      const std::string fileDir = outDir + "/" + s.first + "/" + c.first + "/";
      if (!existDir(fileDir)) { std::cout << "[INFO] Directory " << fileDir << "does not exist!" << std::endl; redoHist = true; return false; }
      for (auto& pd : c.second) {
	// Step 0: Create the output file
	const std::string fileName = fileDir + "effContainer_"+s.first+"_"+c.first+"_"+pd.first+".root";
	if (!existFile(fileName)) { std::cout << "[INFO] File " << fileName << "does not exist!" << std::endl; redoHist = true; return false; }
	TFile file(fileName.c_str(), "READ");
	file.cd();
	if (file.IsOpen() && !file.IsZombie()) {
	  // Step 1: Store all the 1D Efficiency objects
	  for (auto& v1 : pd.second) {
	    for (auto& v2 : v1.second) {
	      for (auto& v3 : v2.second) {
		for (auto& b3 : v3.second) {
		  auto& hist = std::get<0>(b3.second);
		  auto& n = std::get<1>(b3.second);
		  auto& m = std::get<2>(b3.second);
		  auto objP = dynamic_cast<TH1D*>(file.Get(hist.GetName()));
		  if (objP) { hist = *objP; }
		  else { std::cout << "[ERROR] Object stored in " << fileName << " from " << v3.first << " does not exist" << std::endl; redoHist = true; return false; }
		  const std::string lbl = hist.GetTitle();
		  m = stod(lbl.substr(0,lbl.find("_")));
		  n = stod(lbl.substr(lbl.find("_")+1));
		}
	      }
	    }
	  }
	  // Step 3: Close the file
	  file.Close("R");
	}
	else { std::cout << "[INFO] File " << fileName << " could not be opened" << std::endl; redoHist = true; return false; }
	std::cout << "[INFO] Extracted information from file " << fileName << std::endl;
      }
    }
  }
  std::cout << "[INFO] Extracting histograms done!" << std::endl;
  return !redoHist;
};
 

void findThreshold(double& val, double& unc_low, double& unc_high, const TH1D& hist, const double& eff_thr)
{
  //
  val = -99.; unc_low = -99.; unc_high = -99.;
  const double total = hist.Integral();
  if (total==0.0) { return; }
  //
  const int pass = int(total*eff_thr);
  const int tot = int(total);
  const auto eff_high = TEfficiency::ClopperPearson(tot, pass, 0.682689492137, true);
  const auto eff_low  = TEfficiency::ClopperPearson(tot, pass, 0.682689492137, false);
  //
  double thr = -999.0, thr_low = -999.0, thr_high = -999.0, binError=-999.0;
  const auto& nBins = hist.GetNbinsX();
  const auto nStart = hist.GetXaxis()->FindBin(0.0);
  for(int i=nStart; i<nBins; i++) {
    const auto x = hist.GetBinCenter(i);
    const auto ratio = hist.Integral(0, i)/total;
    if(ratio>=eff_low && thr_low==-999.0) { thr_low = x; }
    if(ratio>=eff_thr && thr==-999.0) { thr = x; binError = hist.GetBinWidth(i)/2.0; }
    if(ratio>=eff_high && thr_high==-999.0) { thr_high = x; break; }
  }
  if (thr==-999.0) { throw std::runtime_error("[ERROR] Threshold is -999"); }
  if (thr_low==-999.0) { throw std::runtime_error("[ERROR] Threshold low error is -999"); }
  if (thr_high==-999.0) { throw std::runtime_error("[ERROR] Threshold high error is -999"); }
  //
  val = thr;
  unc_low = std::sqrt(std::pow(std::abs(thr - thr_low), 2.0) + std::pow(binError, 2.0));
  unc_high = std::sqrt(std::pow(std::abs(thr - thr_high), 2.0) + std::pow(binError, 2.0));
};


void getThreshold(BinValQuadMap_t& thrMap, const BinHistQuadMap_t& histMap)
{
  //
  // Loop over histograms
  std::cout << "[INFO] Computing the thresholds" << std::endl;
  for (const auto& s : histMap) {
    std::cout << "[INFO] Processing sample: " << s.first << std::endl;
    for (const auto& c : s.second) {
      for (const auto& pd : c.second) {
	for (const auto& v1 : pd.second) {
	  for (const auto& v2 : v1.second) {
	    for (const auto& v3 : v2.second) {
	      for (const auto& effThr : EFF_THR_) {
		auto& thrM = thrMap[s.first][c.first][pd.first][v1.first][v2.first][v3.first][effThr];
		uint iRow = 0;
		for (const auto& b3 : v3.second) {
		  auto& hist = std::get<0>(b3.second);
		  auto& n = std::get<1>(b3.second);
		  auto& m = std::get<2>(b3.second);
		  double val=0., unc_low=0., unc_high=0.;
		  findThreshold(val, unc_low, unc_high, hist, effThr);
		  thrM(iRow, 0) = m/n;
		  thrM(iRow, 1) = b3.first.low();
		  thrM(iRow, 2) = b3.first.high();
		  thrM(iRow, 3) = val;
		  thrM(iRow, 4) = unc_low;
		  thrM(iRow, 5) = unc_high;
		  iRow += 1;
		}
	      }
	    }
	  }
	}
      }
    }
  }
};


bool extractThreshold(BinValQuadMap_t& thrMap, const std::string& outDir)
{
  std::cout << "[INFO] Extracting the thresholds" << std::endl;
  bool redoTHR = false;
  for (auto& s : thrMap) {
    for (auto& c : s.second) {
      const std::string fileDir = outDir + "/" + s.first + "/" + c.first + "/";
      if (!existDir(fileDir)) { std::cout << "[INFO] Directory " << fileDir << "does not exist!" << std::endl; redoTHR = true; return false; }
      for (auto& pd : c.second) {
	// Step 0: Define the input file
	const std::string fileName = fileDir + "thrContainer_"+s.first+"_"+c.first+"_"+pd.first+".root";
	if (!existFile(fileName)) { std::cout << "[INFO] File " << fileName << "does not exist!" << std::endl; redoTHR = true; return false; }
	TFile file(fileName.c_str(), "READ");
	file.cd();
	if (file.IsOpen() && !file.IsZombie()) {
	  // Step 1: Store all the treshold objects
	  for (auto& v1 : pd.second) {
	    for (auto& v2 : v1.second) {
	      for (auto& v3 : v2.second) {
		for (auto& th : v3.second) {
		  const std::string name = Form("matrixTHR_%s_%s_%s_%s_%.0f_%.0f_%s_%.0f_%.0f_%s_Eff_%.0f",
						s.first.c_str(), c.first.c_str(), pd.first.c_str(),
						v1.first.name().c_str(), v1.first.low()*10., v1.first.high()*10.,
						v2.first.name().c_str(), v2.first.low()*10., v2.first.high()*10.,
						v3.first.c_str(), th.first*100.);
		  auto& thr = th.second;
		  auto objP = dynamic_cast<TMatrixD*>(file.Get(name.c_str()));
		  if (objP) { thr = *objP; }
		  else { std::cout << "[ERROR] Object " <<  name << " searched in " << fileName << " does not exist" << std::endl; redoTHR = true; return false; }
		  std::cout << th.first << std::endl; thr.Print();
		}
	      }
	    }
	  }
	  // Step 3: Close the file
	  file.Close("R");
	}
	else { std::cout << "[INFO] File " << fileName << " could not be opened" << std::endl; redoTHR = true; return false; }
	std::cout << "[INFO] Extracted information from file " << fileName << std::endl;
      }
    }
  }
  std::cout << "[INFO] Extracting thresholds done!" << std::endl;
  return !redoTHR;
};


bool storeThreshold(const BinValQuadMap_t& thrMap, const std::string& outDir)
{
  std::cout << "[INFO] Storing the thresholds" << std::endl;
  bool notStored = false;
  for (const auto& s : thrMap) {
    for (const auto& c : s.second) {
      // Step 0: Create the output file
      const std::string fileDir = outDir + "/" + s.first + "/" + c.first + "/";
      makeDir(fileDir);
      for (const auto& pd : c.second) {
	const std::string fileName = fileDir + "thrContainer_"+s.first+"_"+c.first+"_"+pd.first+".root";
	TFile file(fileName.c_str(), "RECREATE");
	file.cd();
	if (file.IsOpen() && !file.IsZombie()) {
	  // Step 1: Store all the threshold objects
	  for (const auto& v1 : pd.second) {
	    for (const auto& v2 : v1.second) {
	      for (const auto& v3 : v2.second) {
		for (const auto& th : v3.second) {
		  const std::string name = Form("matrixTHR_%s_%s_%s_%s_%.0f_%.0f_%s_%.0f_%.0f_%s_Eff_%.0f",
						s.first.c_str(), c.first.c_str(), pd.first.c_str(),
						v1.first.name().c_str(), v1.first.low()*10., v1.first.high()*10.,
						v2.first.name().c_str(), v2.first.low()*10., v2.first.high()*10.,
						v3.first.c_str(), th.first*100.);
		  th.second.Write(name.c_str());
		}
	      }
	    }
	  }
	  // Step 3: Write and Close the file
	  file.Write();
	  file.Close("R");
	}
	else { std::cout << "[INFO] File " << fileName << " could not be created" << std::endl; notStored = true; }
	std::cout << "[INFO] Information store in file " << outDir << fileName << std::endl;
      }
    }
  }
  std::cout << "[INFO] Storing thresholds done!" << std::endl;
  return !notStored;
};


void setRangeYAxisGraph(TGraphAsymmErrors& graph, const double& fracDo = 0.1, const double& fracUp = 0.5)
{
  // Find maximum and minimum points of Plot to rescale Y axis
  double YMax = -1e99;
  for (int i=0; i<=graph.GetN(); i++) { double x, y; graph.GetPoint(i, x, y); y += graph.GetErrorY(i); YMax = std::max(YMax, y); }
  double YMin = 1e99;
  for (int i=0; i<=graph.GetN(); i++) { double x, y; graph.GetPoint(i, x, y); y -= graph.GetErrorY(i); YMin = std::min(YMin, y); }
  //
  const double YTot = (YMax - YMin)/(1.0 - fracUp - fracDo);
  const double Yup   = YMax + fracUp*YTot;
  const double Ydown = YMin - fracDo*YTot;
  //
  graph.GetYaxis()->SetRangeUser(Ydown, Yup);
};


std::map<std::string, std::string> VARMAP1 = {
					      {"Cand_Pt", "p_{T}"},
					      {"Cand_Rap", "y"},
					      {"Cand_AbsRap", "|y|"},
					      {"NTrack", "N^{Trk}"}
};

std::map<std::string, std::string> VARMAP2 = {
					      {"Cand_Pt", "#mu^{+}#mu^{-} p_{T} (GeV/c)"},
					      {"Cand_Rap", "#mu^{+}#mu^{-} y"},
					      {"Cand_AbsRap", "#mu^{+}#mu^{-} |y|"},
					      {"Cand_RapCM", "#mu^{+}#mu^{-} y_{CM}"},
					      {"NTrack", "Track multiplicity"}
};


void formatGraph(TGraphAsymmErrors& graph, const std::string& var, const double& effThr, const std::string yVar="")
{
  //
  // Set the Axis Titles
  std::string xLabel = VARMAP2.at(var);
  std::string yLabel = Form("#font[12]{l}^{J/#psi}_{xyz} threshold (%.0f%%)", effThr*100.);
  if (yVar!="") { yLabel = yVar; }
  graph.SetTitle(Form(";%s;%s", xLabel.c_str(), yLabel.c_str()));
  //
  // General
  graph.SetMarkerColor(kBlack);
  graph.SetLineColor(kBlack); 
  graph.SetMarkerStyle(20);
  graph.SetMarkerSize(1.5);
  graph.SetLineWidth(3);
  graph.SetLineStyle(1);
  graph.SetFillStyle(0);
  // X-axis
  graph.GetXaxis()->CenterTitle(kTRUE);
  graph.GetXaxis()->SetTitleOffset(0.75);
  graph.GetXaxis()->SetTitleSize(0.065);
  graph.GetXaxis()->SetLabelSize(0.035);
  double xMin = graph.GetXaxis()->GetXmin()-0.1;
  double xMax = graph.GetXaxis()->GetXmax()+0.1;
  if      (var=="Cand_Pt"    ) { xMin = -0.1; xMax = 66.0; }
  else if (var=="Cand_Rap"   ) { xMin = -2.6; xMax = 2.6; }
  else if (var=="Cand_RapCM" ) { xMin = -2.9; xMax = 2.0; }
  else if (var=="Cand_AbsRap") { xMin = -0.1; xMax = 2.6; }
  else if (var=="NTrack") { xMin = 0.0; xMax = 275.0; }
  graph.GetXaxis()->SetRangeUser(xMin, xMax);
  if (var=="NTrack") graph.GetXaxis()->SetNdivisions(510);
  else graph.GetXaxis()->SetNdivisions(505);
  // Y-axis
  graph.GetYaxis()->CenterTitle(kTRUE);
  graph.GetYaxis()->SetTitleOffset(1.05);
  graph.GetYaxis()->SetTitleSize(0.065);
  graph.GetYaxis()->SetLabelSize(0.035);
  setRangeYAxisGraph(graph, 0.1, 0.4);
};


void plotThreshold(const BinValQuadMap_t& thrMap, const std::string& outDir, const bool& altFunc)
{
  // Set the CMS style
  setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TGaxis::SetMaxDigits(3); // to display powers of 10
  //
  std::cout << "[INFO] Plotting the thresholds" << std::endl;
  for (const auto& s : thrMap) {
    for (const auto& cc : s.second) {
      for (const auto& pd : cc.second) {
	for (const auto& v1 : pd.second) {
	  //
	  std::map< double , std::map< std::string , std::map< std::string , std::map< BinF_t , std::pair<double, double> > > > > fitResV;
	  //
	  for (const auto& v2 : v1.second) {
	    for (const auto& v3 : v2.second) {
	      for (const auto& th : v3.second) {
		//
		// Create Canvas
		TCanvas c("c", "c", 1000, 1000); c.cd();
		//
		// Create the Text Info
		TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035); float dy = 0;
		std::vector< std::string > textToPrint;
		std::string sampleLabel = (s.first.rfind("JPsi")!=std::string::npos ? "J/#psi" : "#psi(2S)");
		sampleLabel += "#rightarrow #mu^{+}#mu^{#font[122]{\55}}";
		textToPrint.push_back(sampleLabel);
		if (v1.first.name()!="NONE") { textToPrint.push_back(Form("%g #leq %s < %g", v1.first.low(), VARMAP1.at(v1.first.name()).c_str(), v1.first.high())); }
		if (v2.first.name()!="NONE") { textToPrint.push_back(Form("%g #leq %s < %g", v2.first.low(), VARMAP1.at(v2.first.name()).c_str(), v2.first.high())); }
		//
		const auto& thr = th.second;
		const auto& nRows = thr.GetNrows();
		// 
		// Create graph
		TGraphAsymmErrors gr(nRows);
		double xMin = 99999., xMax = -99999.;
		for (int i=0; i<nRows; i++) {
		  const auto& xRngLo = thr(i, 1);
		  const auto& xRngHi = thr(i, 2);
		  const auto xVal = thr(i, 0);
		  const auto xErrLo = (xVal - xRngLo);
		  const auto xErrHi = (xRngHi - xVal);
		  const auto yVal = thr(i, 3);
		  const auto yErrLo = thr(i, 4);
		  const auto yErrHi = thr(i, 5);
		  //
		  if (yVal!=-99.) {
		    xMin = std::min(xMin, xRngLo);
		    xMax = std::max(xMax, xRngHi);
		    gr.SetPoint(i, xVal, yVal);
		    gr.SetPointError(i, xErrLo, xErrHi, yErrLo, yErrHi);
		  }
		}
		formatGraph(gr, v3.first, th.first);
		//
		// Fit graph
		std::string fitFunc = "[0] + [1]/pow(x,[2])";
		if (altFunc) { fitFunc = "[0]*(1 - [1]*exp(-1/x))"; }
		if (v3.first=="NTrack") { fitFunc = "[0]"; }
		double xFMin = xMin, xFMax = xMax;
		if (v3.first=="Cand_Pt") { xFMin = std::max(3.0, xFMin); }
		TF1 func("function", fitFunc.c_str(), xFMin, xFMax);
		std::string fitResult = "";
		if (v3.first=="Cand_Pt" && altFunc) {
		  func.SetParameter(0, 0.4);
		  func.SetParameter(1, 0.4);
		  gr.Fit(&func, "QMRN");
		  const auto p1 = func.GetParameter(0);
		  const auto p2 = func.GetParameter(1);
		  const auto pE1 = func.GetParError(0);
		  const auto pE2 = func.GetParError(1);
		  fitResult = Form("%.2f*(1 - %.2f*e^{-1/p_{T}})", p1, p2);
		  fitResV[th.first][v3.first]["p1"][v2.first] = std::make_pair(p1, pE1);
		  fitResV[th.first][v3.first]["p2"][v2.first] = std::make_pair(p2, pE2);
		}
		else if (v3.first=="Cand_Pt") {
		  func.SetParameter(0, 0.02);
		  func.SetParameter(1, 0.2);
		  func.SetParameter(2, 1.0);
		  gr.Fit(&func, "QMRN");
		  const auto p1 = func.GetParameter(0);
		  const auto p2 = func.GetParameter(1);
		  const auto p3 = func.GetParameter(2);
		  const auto pE1 = func.GetParError(0);
		  const auto pE2 = func.GetParError(1);
		  const auto pE3 = func.GetParError(2);
		  std::cout << p1 << "  " << p2 << std::endl;
		  fitResult = Form("%.3f + #frac{%.2f}{p_{T}^{%.2f}}", p1, p2, p3);
		  fitResV[th.first][v3.first]["p1"][v2.first] = std::make_pair(p1, pE1);
		  fitResV[th.first][v3.first]["p2"][v2.first] = std::make_pair(p2, pE2);
		  fitResV[th.first][v3.first]["p3"][v2.first] = std::make_pair(p3, pE3);
		}
		else if (v3.first=="NTrack") {
		  gr.Fit(&func, "QMRNF");
		  const auto& p1 = func.GetParameter(0);
		  fitResult = Form("Mean: %.3f", p1);
		}
		func.SetLineColor(kRed);
		func.SetLineWidth(2.0);
		//
		// Draw the graph
		gr.Draw("ap");
		if (v3.first=="Cand_Pt" || v3.first=="NTrack") {
		  func.Draw("same");
		}
		c.Modified(); c.Update();
		//
		// Draw the text
		tex.SetTextSize(0.055); tex.DrawLatex(0.22, 0.84, textToPrint[0].c_str());
		tex.SetTextSize(0.060); tex.SetTextFont(61); tex.DrawLatex(0.78, 0.84, "CMS"); tex.SetTextFont(62);
		tex.SetTextSize(0.046); tex.SetTextFont(52); tex.DrawLatex(0.69, 0.79, "Preliminary"); tex.SetTextFont(62);
		for (uint i=1; i<textToPrint.size(); i++) { tex.SetTextSize(0.045); tex.DrawLatex(0.22, 0.76-dy, textToPrint[i].c_str()); dy+=0.060; }
		if (fitResult!="") { tex.SetTextSize(0.040); tex.DrawLatex(0.62, 0.70, fitResult.c_str()); }
		c.Modified(); c.Update(); // Pure paranoia
		//
		// set the CMS style
		StringVector_t lumiLabels; getLumiLabels(lumiLabels, "", cc.first, true);
		CMS_lumi(&c, 33, (" "+lumiLabels[0]), lumiLabels[1], false, 0.65, false);
		c.Modified(); c.Update(); // Pure paranoia
		//
		// Create Output Directory
		const std::string& plotDir = outDir+"/Plots"+(altFunc?"_AltFunc/":"/")+s.first+"/"+cc.first+"/"+Form("THR_%.0f",th.first*100.);
		makeDir(plotDir + "/png/");
		makeDir(plotDir + "/pdf/");
		//
		// Save Canvas
		const std::string name = Form("plotTHR_%s_%s_%s_%s_%.0f_%.0f_%s_%.0f_%.0f_%s_Eff_%.0f",
					      s.first.c_str(), cc.first.c_str(), pd.first.c_str(),
					      v1.first.name().c_str(), v1.first.low()*10., v1.first.high()*10.,
					      v2.first.name().c_str(), v2.first.low()*10., v2.first.high()*10.,
					      v3.first.c_str(), th.first*100.);
		c.SaveAs((plotDir + "/png/"  + name + ".png" ).c_str());
		c.SaveAs((plotDir + "/pdf/"  + name + ".pdf" ).c_str());
		//
		// Clean up memory
		c.Clear(); c.Close();
	      }
	    }
	  }
	  //
	  for (const auto& th : fitResV) {
	    for (const auto& v3 : th.second) {
	      for (const auto& p : v3.second) {
		//
		// Create Canvas
		TCanvas c("c", "c", 1000, 1000); c.cd();
		//
		// Create the Text Info
		TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035); float dy = 0;
		std::vector< std::string > textToPrint;
		std::string sampleLabel = (s.first.rfind("JPsi")!=std::string::npos ? "J/#psi" : "#psi(2S)");
		sampleLabel += "#rightarrow #mu^{+}#mu^{#font[122]{\55}}";
		textToPrint.push_back(sampleLabel);
		if (v1.first.name()!="NONE") { textToPrint.push_back(Form("%g #leq %s < %g", v1.first.low(), VARMAP1.at(v1.first.name()).c_str(), v1.first.high())); }
		//
		const auto& nBins = p.second.size();
		// 
		// Create graph
		TGraphAsymmErrors gr(nBins);
		double xMin = 99999., xMax = -99999.;
		uint iBin = 0;
		for (const auto& v2 : p.second) {
		  const auto& xRngLo = v2.first.low();
		  const auto& xRngHi = v2.first.high();
		  const auto xVal = (xRngHi + xRngLo)/2.0;
		  const auto xErr = (xRngHi - xRngLo)/2.0;
		  const auto yVal = v2.second.first;
		  const auto yErrLo = v2.second.second;
		  const auto yErrHi = yErrLo;
		  //
		  if (yVal!=-99.) {
		    xMin = std::min(xMin, xVal-xErr);
		    xMax = std::max(xMax, xVal+xErr);
		    gr.SetPoint(iBin, xVal, yVal);
		    gr.SetPointError(iBin, xErr, xErr, yErrLo, yErrHi);
		  }
		  iBin += 1;
		}
		formatGraph(gr, p.second.begin()->first.name(), th.first, p.first);
		//
		// Draw the graph
		gr.Draw("ap");
		c.Modified(); c.Update();
		//
		// Draw the text
		tex.SetTextSize(0.055); tex.DrawLatex(0.22, 0.84, textToPrint[0].c_str());
		tex.SetTextSize(0.060); tex.SetTextFont(61); tex.DrawLatex(0.78, 0.84, "CMS"); tex.SetTextFont(62);
		tex.SetTextSize(0.046); tex.SetTextFont(52); tex.DrawLatex(0.69, 0.79, "Preliminary"); tex.SetTextFont(62);
		for (uint i=1; i<textToPrint.size(); i++) { tex.SetTextSize(0.045); tex.DrawLatex(0.22, 0.76-dy, textToPrint[i].c_str()); dy+=0.060; }
		c.Modified(); c.Update(); // Pure paranoia
		//
		// set the CMS style
		StringVector_t lumiLabels; getLumiLabels(lumiLabels, "", cc.first, true);
		CMS_lumi(&c, 33, (" "+lumiLabels[0]), lumiLabels[1], false, 0.65, false);
		c.Modified(); c.Update(); // Pure paranoia
		//
		// Create Output Directory
		const std::string& plotDir = outDir+"/Plots/"+s.first+"/"+cc.first+"/"+Form("THR_%.0f",th.first*100.);
		makeDir(plotDir + "/png/");
		makeDir(plotDir + "/pdf/");
		//
		// Save Canvas
		const std::string name = Form("parTHR_%s_%s_%s_%s_%s_%.0f_%.0f_Eff_%.0f",
					      p.first.c_str(), s.first.c_str(), cc.first.c_str(), pd.first.c_str(),
					      v1.first.name().c_str(), v1.first.low()*10., v1.first.high()*10., th.first*100.);
		c.SaveAs((plotDir + "/png/"  + name + ".png" ).c_str());
		c.SaveAs((plotDir + "/pdf/"  + name + ".pdf" ).c_str());
		//
		// Clean up memory
		c.Clear(); c.Close();
	      }
	    }
	  }
	}
      }
    }
  }
};
	     
	     
