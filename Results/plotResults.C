#ifndef plotResults_C
#define plotResults_C

// Auxiliary Headers
#include "Utilities2/ResultManager.h"


void plotResults(
		 const std::string& anaDirName = "Psi2S/Charmonia_Fit",
		 const std::string& nominalWorkDirName = "Nominal",
		 const bool& doSyst = true,
		 const StringVector_t& colTags  = {"PA8Y16"} //"PP5Y17", "PP13Y18", "pPb8Y16", "Pbp8Y16", "PA8Y16", "PbPb5Y15"
               )
{
  //
  // Define the Systematic varations
  // Variation   Method: +1: Fully Uncorrelated , +2: Fully Correlated
  // Propagation Method: -1: Fully Uncorrelated , -2: Fully Correlated
  WSDirMap_t workDirInfo = {
			    { "Nominal" ,
			      {
			       { "Nominal"  , { { nominalWorkDirName } , {{"Rap", 1}, {"Obj", 1}} } }
			      }
			    },
			    // Background
			    { "Background" ,
			      {
			       { "Shape" , { { "Systematic_Background_ExpCheb" } , {{"Rap", 1}, {"Obj", 2}} } },
			       { "LLR"   , { { "Systematic_Background_LLR25",
					       "Systematic_Background_LLR100"  } , {{"Rap", 1}, {"Obj", 2}} } }
			      }
			    },
			    // Decay Length cut
			    { "DecayCut" ,
			      {
			       { "Threshold" , { { "Systematic_DecayCut_85",
						   "Systematic_DecayCut_95"    } , {{"Rap", 2}, {"Obj", 2}} } },
			       { "Function"  , { { "Systematic_DecayCut_Alt"   } , {{"Rap", 2}, {"Obj", 2}} } },
			       { "Object"    , { { "Systematic_DecayCut_Psi2S" } , {{"Rap", 2}, {"Obj", 2}} } }
			      }
			    }
  };
  //
  const std::string    dataTag  = "DATA";
  const StringVector_t trgTags  = {"DIMUON", "MINBIAS"};
  const StringVector_t objTags  = {"JPsi"};
  const std::string&   varTag   = "Cand_Mass";
  const StringVector_t dLenDirs = {"Prompt_Cut", "NonPrompt_Cut"};
  const StringVector_t subDirs = {"Inclusive", "NTrack_15_250", "NTrack_15_50", "NTrack_120_250", "NTrack_50_80", "NTrack_80_120"};
  //
  // initialize the result manager
  ResultManager result;
  result.setFitDir(anaDirName, subDirs, dLenDirs);
  result.setFitInfo(dataTag, varTag, trgTags, colTags, objTags);
  result.setWorkDirInfo(workDirInfo);
  result.doSyst(doSyst);
  //
  // extract and process the information
  if (!result.extractAndProcess()) { return; }
  //
  // define binning
  result.defineBins();
  //
  // plot results
  result.plot();
  /*
  //
  // Create the main plots
  //
  // Define the output directory
  const std::string& CWD = getcwd(NULL, 0);
  auto vTag = varTag; if (vTag.find("_")!=std::string::npos) { vTag.erase(vTag.find("_"), 1); }
  auto tTag = trgTags[0]; for (uint i=1; i<trgTags.size(); i++) { tTag += "_"+trgTags[i]; }
  auto cTag = colTags[0]; for (uint i=1; i<colTags.size(); i++) { cTag += "_"+colTags[i]; }
  auto oTag = objTags[0]; for (uint i=1; i<objTags.size(); i++) { oTag += "_"+objTags[i]; }
  const std::string& outDir = CWD + "/Output/Results/" + nominalWorkDirName+"/" + vTag+"/" + dataTag+"_"+tTag+"/" + oTag+"/" + cTag;
  //
  // Draw the plots
  const bool& isMC = (dataTag!="DATA");
  drawResultsGraph(graphMap, outDir, isMC);
  */
  //
};


#endif // ifndef plotResults_C
