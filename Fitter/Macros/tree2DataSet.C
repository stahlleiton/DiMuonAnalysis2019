// -*- C++ -*-
//
// Package:    Fitter
//
/*
 Description: TTree to RooDataSet converter.
 Implementation:
 This program create RooDataSets from TTrees.
 */
// Original Author:  Andre Stahl,
//         Created:  Feb 28 17:35 CET 2019
//
//
#ifndef tree2DataSet_C
#define tree2DataSet_C

#include "Utilities/initClasses.h"
#include "Candidate/VertexCompositeTree2DataSet.C"


bool checkFileInfo    ( const StringVectorMap_t& fileInfo );
bool checkAnalysis    ( const GlobalInfo& info );


bool tree2DataSet(RooWorkspaceMap_t& workspaces, const StringVectorMap_t& fileInfo, const GlobalInfo& userInfo, const bool& updateDS=false)
{
  if (!checkFileInfo(fileInfo)) return false;
  // Check if Analysis type is supported
  if (!checkAnalysis(userInfo)) return false;
  // Make RooDatasets
  if (userInfo.Par.at("treeType") == "VertexCompositeTree") {
    if (!VertexCompositeTree2DataSet(workspaces, fileInfo, userInfo, updateDS)) return false;
  }
  return true;
};


bool checkAnalysis(const GlobalInfo& info)
{
  const auto& analysis = info.Par.at("analysis");
  if (analysis.rfind("CandTo", 0)==0) {
    if (info.Par.at("treeType") == "VertexCompositeTree") {
      std::cout << "[INFO] Proceed to make "<<info.Par.at("PD")<<" datasets for Mass Resonance analysis using VertexComposite Ntuple" << std::endl; return true;
    }
    else if (info.Par.at("treeType") == "OniaTree") {
      std::cout << "[INFO] Proceed to make "<<info.Par.at("PD")<<" datasets for Mass Resonance analysis using Onia Tree" << std::endl; return true;
    }
    else {
      std::cout << "[INFO] Tree type " << info.Par.at("treeType") << " is invalid for Mass Resonance analysis!" << std::endl; return false;
    }
  }
  std::cout << "[ERROR] The input analysis: " << analysis << " is not supported by this fitter!" << std::endl;
  return false;
};


bool checkFileInfo(const StringVectorMap_t& fileInfo)
{
  if ( fileInfo.at("outputFileDir").size()==0  ) { std::cout << "[ERROR] OutputFileDir is empty!" << std::endl; return false;  }
  if ( fileInfo.at("inputFileNames").size()==0 ) { std::cout << "[ERROR] InputFileNames is empty!" << std::endl; return false; }
  if ( fileInfo.at("dsNames").size()==0        ) { std::cout << "[ERROR] DSNames is empty!" << std::endl; return false;        }
  return true;
};


#endif // #ifndef tree2DataSet_C
