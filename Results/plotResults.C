#ifndef plotResults_C
#define plotResults_C

// Auxiliary Headers
#include "Utilities/resultsUtils.h"
#include "extractResultsTree.C"
#include "processResultsTree.C"


void plotResults(
		 const std::string& workDirName = "Test",
		 const StringVector_t& trgTags  = {"DIMUON", "MINBIAS", "HIGHMULT"}, //"MUON", "DIMUON", "DIMUONPERI", "HIGHMULT", "HIGHMULT2", "MINBIAS", "UPC"
		 const StringVector_t& colTags  = {"PA8Y16"}, //"PP5Y17", "PP13Y18", "pPb8Y16", "Pbp8Y16", "PA8Y16", "PbPb5Y15"
		 const StringVector_t& objTags  = {"JPsi"}, //"Bkg", "JPsi", "Psi2S", "Ups1S", "Ups2S", "Ups3S", "Z", "D0"
		 const std::string&    dataTag  = "DATA",
		 const std::string&    varTag   = "Cand_Mass",
		 const std::string&    effType  = "",
		 const std::string&    accType  = ""
               )
{
  //
  // Initialize the Result
  BinSextaMap_t var;
  //
  // Get the Result
  VarBinTriMap_t inputVar;
  if (workDirName!="") {
    std::cout << "[INFO] Adding results for: " << workDirName << std::endl;
    for (const auto& trgTag : trgTags) {
      for (const auto& colTag : colTags) {
	for (const auto& objTag : objTags) {
	  if (!extractResultsTree(inputVar, workDirName, trgTag, colTag, objTag, dataTag, varTag)) { return; }
	}
      }
    }
    if (!processResults(var, inputVar, effType, accType)) { return; }
  }
  //
  // Define the bins
  BinCont_t binMap;
  defineBins(binMap, var);
  //
  // Create the main plots
  GraphSextaMap_t graphMap;
  iniResultsGraph(graphMap, binMap, var);
  fillResultsGraph(graphMap, binMap, var);
  //
  // Define the output directory
  const std::string& CWD = getcwd(NULL, 0);
  auto vTag = varTag; if (vTag.find("_")!=std::string::npos) { vTag.erase(vTag.find("_"), 1); }
  auto tTag = trgTags[0]; for (uint i=1; i<trgTags.size(); i++) { tTag += "_"+trgTags[i]; }
  auto cTag = colTags[0]; for (uint i=1; i<colTags.size(); i++) { cTag += "_"+colTags[i]; }
  auto oTag = objTags[0]; for (uint i=1; i<objTags.size(); i++) { oTag += "_"+objTags[i]; }
  const std::string& outDir = CWD + "/Output/Results/" + workDirName+"/" + vTag+"/" + dataTag+"_"+tTag+"/" + oTag+"/" + cTag;
  //
  // Draw the plots
  const bool& isMC = (dataTag!="DATA");
  drawResultsGraph(graphMap, outDir, isMC);
  //
};


#endif // ifndef plotResults_C
