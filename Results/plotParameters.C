#ifndef plotParameters_C
#define plotParameters_C

// Auxiliary Headers
#include "Utilities/resultsUtils.h"
#include "extractResultsTree.C"
#include "processResultsTree.C"


void plotParameters(
		    const std::string& workDirName = "Test",
		    const StringVector_t& objTags  = { "JPsi" },    //"Bkg", "JPsi", "Psi2S", "Ups1S", "Ups2S", "Ups3S", "Z", "D0"
		    const StringVector_t& trgTags  = { "CatPR_DIMUON" },   //"MUON", "DIMUON", "DIMUONPERI", "HIGHMULT", "HIGHMULT2", "MINBIAS", "UPC"
		    const StringVector_t& colTags  = { "PA8Y16" }, //"PP5Y17", "PP13Y18", "pPb8Y16", "Pbp8Y16", "PA8Y16", "PbPb5Y15"
		    const std::string&    dataTag  = "MC",         //"DATA", "MC" 
		    const std::string&    varTag   = "Cand_Mass"     //"Cand_Mass"
		    )
{
  //
  // Check input
  if (workDirName=="") { std::cout << "[ERROR] Name of working directory was not set!" << std::endl; return; }
  if (trgTags.empty()) { std::cout << "[ERROR] PD tags were not set!" << std::endl; return; }
  if (colTags.empty()) { std::cout << "[ERROR] Collision tags were not set!" << std::endl; return; }
  if (objTags.empty()) { std::cout << "[ERROR] Object tags were not set!" << std::endl; return; }
  if (dataTag==""    ) { std::cout << "[ERROR] Data tag was not set!" << std::endl; return; }
  if (varTag==""     ) { std::cout << "[ERROR] Fit variable tag was not set!" << std::endl; return; }
  //
  // Get the result
  VarBinTriMap_t inputVar;
  std::cout << "[INFO] Adding results for: " << workDirName << std::endl;
  for (const auto& trgTag : trgTags) {
    for (const auto& colTag : colTags) {
      for (const auto& objTag : objTags) {
	if (!extractResultsTree(inputVar, workDirName, trgTag, colTag, objTag, dataTag, varTag)) { return; }
      }
    }
  }
  // Extract the parameters
  BinSextaMap_t par;
  addParameters(par, inputVar);
  //
  // Define the bins
  BinCont_t binMap;
  defineBins(binMap, par);
  //
  // Create the plots
  GraphSextaMap_t graphMap;
  iniResultsGraph(graphMap, binMap, par);
  fillResultsGraph(graphMap, binMap, par);
  //
  // Define the output directory
  const std::string& CWD = getcwd(NULL, 0);
  auto vTag = varTag; if (vTag.find("_")!=std::string::npos) { vTag.erase(vTag.find("_"), 1); }
  auto tTag = trgTags[0]; for (uint i=1; i<trgTags.size(); i++) { tTag += "_"+trgTags[i]; }
  auto cTag = colTags[0]; for (uint i=1; i<colTags.size(); i++) { cTag += "_"+colTags[i]; }
  auto oTag = objTags[0]; for (uint i=1; i<objTags.size(); i++) { oTag += "_"+objTags[i]; }
  const std::string& outDir = CWD + "/Output/Parameters/" + workDirName+"/" + vTag+"/" + dataTag+"_"+tTag+"/" + oTag+"/" + cTag;
  //
  // Draw the plots
  const bool& isMC = (dataTag!="DATA");
  drawResultsGraph(graphMap, outDir, isMC);
  //
};


#endif // ifndef plotParameters_C
