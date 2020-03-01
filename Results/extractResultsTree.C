#ifndef extractResultsTree_C
#define extractResultsTree_C

// Auxiliary Files
#include "../Utilities/dataUtils.h"
#include "Utilities/resultsUtils.h"
#include "storeWS2ResultsTree.C"
// ROOT headers
#include "TFile.h"
#include "TTree.h"
// c++ headers
#include <iostream>
#include <string>


void getInfoFromTree ( GlobalInfo& info , TTree& tree );


bool extractResultsTree(
			VarBinTriMap_t& inputVar,
			const std::string& workDirName,
			const std::string& trgTag  = "DIMUON",
			const std::string& colTag  = "PA8Y16",
			const std::string& objTag  = "JPsi",
			const std::string& dataTag = "DATA",
			const std::string& varTag  = "Cand_Mass"
			)
{
  //
  std::cout << "[INFO] Extracting results from " << workDirName << " " << dataTag << " " << trgTag << " " << colTag << " " << objTag << " fits" << std::endl;
  //
  // --------------------------------------------------------------------------------- //
  //
  // Define the input file info
  const std::string& CWD = getcwd(NULL, 0);
  const std::string& inputFileName = "tree_allvars.root";
  const auto& dsTag = (dataTag+"_"+(dataTag=="DATA" ? "" : (objTag+(trgTag.rfind("Cat",0)==0?"":"_")))+trgTag);
  auto vTag = varTag; if (vTag.find("_")!=std::string::npos) { vTag.erase(vTag.find("_"), 1); }
  const auto& inputDirPath  = (CWD+"/Tree/"+workDirName+"/"+vTag+"/"+dsTag+"/"+objTag+"/"+colTag);
  const auto& inputFilePath = (inputDirPath+"/"+inputFileName);
  //
  // Define the tree info container
  GlobalInfo info;
  //
  // --------------------------------------------------------------------------------- //
  //
  // Open the input file
  std::unique_ptr<TFile> inputFile;
  if (existFile(inputFilePath)) { inputFile.reset(TFile::Open(inputFilePath.c_str(), "READ")); }
  if (!inputFile || !inputFile->IsOpen() || inputFile->IsZombie()) {
    std::cout << "[WARNING] The input file " << inputFilePath << " was not found, will create it!" << std::endl; if (inputFile) { inputFile->Close(); }
    if (!storeWS2ResultsTree(workDirName, trgTag, colTag, objTag, dataTag, varTag)) { return false; };
    inputFile.reset(TFile::Open(inputFilePath.c_str(), "READ"));
    if (!inputFile || !inputFile->IsOpen() || inputFile->IsZombie()) {
      std::cout << "[ERROR] The input file " << inputFilePath << " could not be re-created!" << std::endl; return false;
    }
  }
  //
  // Extract the input tree
  const auto& tree = dynamic_cast<TTree*>(inputFile->Get("fitResults"));
  if (tree==NULL) { std::cout << "[ERROR] The input tree fitResults was not found in " << inputFilePath << "" << std::endl; inputFile->Close(); return false; }
  // Set the address of the tree branches
  getInfoFromTree(info, *tree);
  //
  // Extract the input variables
  for (uint i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    //
    // Extract the string information
    const auto& col = *info.StrP.at("fitSystem");
    const auto& PD  = *info.StrP.at("PD");
    const auto& obj = objTag;
    //
    // Determine the analysis bin
    anabin bin;
    for (const auto& v : info.Var) {
      if (v.first.rfind("OBS_",0)!=0) continue;
      auto name = v.first; name.erase(name.find("OBS_"), 4);
      auto var = v.second;
      auto min = var.at("Min");
      auto max = var.at("Max");
      // Check if bin was properly set
      if (min==-99. || max==-99.) { std::cout << "[WARNING] The bin of " << v.first << " was not set properly!" << std::endl; continue; }
      // Ignore fit variable
      if (name=="Cand_Mass") continue;
      // Check if bin is set to default values
      const auto& defMin = v.second.at("DefaultMin");
      const auto& defMax = v.second.at("DefaultMax");
      if (min==defMin && max==defMax) continue;
      // Divide by two if centrality
      if (name=="Centrality") { min /= 2.0; max /= 2.0; }
      // Round up the bin boundaries to 3rd decimal
      roundValue(min, 3); roundValue(max, 3);
      // Get mean and width
      auto mean = var.at("Val");
      auto width = var.at("Err");
      roundValue(mean, 3); roundValue(width, 3);
      // Store the bin
      bin.setbin(name, min, max, mean, width);
    }
    //
    // Store the parameter information
    for (const auto& v : info.Var) {
      if (v.first.rfind("PAR_",0)==0) {
	auto name = v.first; name.erase(0, 4);
        if (contain(v.second, "Val")) { inputVar[obj][col][PD][bin][name]["Val"] = v.second.at("Val"); }
        inputVar[obj][col][PD][bin][name]["Err_Stat_High"] = contain(v.second, "ErrHi") ? v.second.at("ErrHi") : 0.0;
        inputVar[obj][col][PD][bin][name]["Err_Stat_Low" ] = contain(v.second, "ErrLo") ? v.second.at("ErrLo") : 0.0;
        inputVar[obj][col][PD][bin][name]["Err_Syst_High"] = 0.0;
	inputVar[obj][col][PD][bin][name]["Err_Syst_Low" ] = 0.0;
      }
      else {
	auto name = v.first;
	if (v.first.rfind("OBS_",0)==0) { name.erase(0, 4); }
	for (const auto& t : v.second) {
	  inputVar[obj][col][PD][bin][name][t.first] = t.second;
	}
      }
    }
  }
  //
  // Close the input file
  inputFile->Close();
  //
  // return
  return true;
};


void getInfoFromTree(GlobalInfo& info , TTree& tree)
{
  //
  // Loop over the tree branches
  const auto& branchList = tree.GetListOfBranches();
  for (int i=0; i<branchList->GetEntries(); i++) {
    const std::string& brName = branchList->At(i)->GetName();
    std::string brType = branchList->At(i)->GetTitle();
    if (brType.rfind("/")!=std::string::npos) { brType = brType.substr(brType.rfind("/")+1); }
    else { brType = tree.GetBranch(brName.c_str())->GetClassName(); }
    const auto& vName = brName.substr(0, brName.rfind("_"));
    const auto& vType = ((brName.rfind("_")!=std::string::npos) ? brName.substr(brName.rfind("_")+1) : "");
    //
    if      (brType=="D"     ) { tree.SetBranchAddress(brName.c_str(), &(info.Var[vName][vType])); }
    else if (brType=="O"     ) { tree.SetBranchAddress(brName.c_str(), &(info.Flag[brName]));      }
    else if (brType=="string") { tree.SetBranchAddress(brName.c_str(), &(info.StrP[brName]));      }
  }
};


#endif // #ifndef extractResultsTree_C
