#ifndef makeInputForLLR_C
#define makeInputForLLR_C

#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <string>
#include <vector>
#include <algorithm>
#include <future>

#include "../Utilities/dataUtils.h"


bool readFile(std::vector<std::string>& rows, const std::string& fileName);


void makeInputForLLR(const std::string& workDirName,
		     const bool& isNominal = true,
		     const std::string& par="JPsi",
		     const std::string& col="PA8Y16")
{
  //
  std::vector<std::string> BKGMODELS;
  if (isNominal) { BKGMODELS = {"Uniform", "Chebychev1", "Chebychev2", "Chebychev3", "Chebychev4", "Chebychev5", "Chebychev6"}; }
  else { BKGMODELS = {"Uniform", "ExpChebychev1", "ExpChebychev2", "ExpChebychev3", "ExpChebychev4", "ExpChebychev5", "ExpChebychev6"}; }
  //
  const std::string& CWD = getcwd(NULL, 0);
  const auto& inputDir = CWD + "/Input/" + workDirName + "/";
  if (existDir(inputDir)==false){ 
    std::cout << "[ERROR] Input directory: " << inputDir << " doesn't exist!" << std::endl;
    return; 
  }
  std::vector<std::string> inputDirList;
  findSubDir(inputDirList, inputDir);
  //
  std::vector<std::string> inputFileList;
  for (const auto& dir : inputDirList) {
    inputFileList.push_back(dir+"InitialParam_Cand_Mass_"+par+"_"+col+".csv");
  }
  //
  std::map<std::string, std::vector<std::string> > inputFileRows;
  for (const auto& file : inputFileList) {
    readFile(inputFileRows[file], file);
  }
  //
  auto outputDir = inputDir;
  stringReplace(outputDir, workDirName, workDirName+"_LLR");
  std::vector<std::string> outputDirList;
  for (const auto& dir : inputDirList) {
    outputDirList.push_back(dir);
    stringReplace(outputDirList.back(), workDirName, workDirName+"_LLR");
  }
  //
  std::vector<std::string> outputFileList;
  std::map<std::string, std::vector<std::string> > outputFileRows;
  for (const auto& file : inputFileList) {
    outputFileList.push_back(file);
    stringReplace(outputFileList.back(), workDirName, workDirName+"_LLR");
    for (const auto& row : inputFileRows.at(file)) {
      for (const auto& bkgModel : BKGMODELS) {
	auto outputRow = row;
	std::string model = "";
	if (row.rfind("[Bkg]")!=std::string::npos) { model = row.substr(0, row.rfind("[Bkg]")); }
	if (model.rfind("+")!=std::string::npos) { model = model.substr(model.rfind("+")+1); }
	if (model!="") {
	  stringReplace(outputRow, model, bkgModel);
	  outputFileRows[outputFileList.back()].push_back(outputRow);
	}
	else {
	  outputFileRows[outputFileList.back()].push_back(outputRow);
	  break;
	}
      }
    }
  }
  //
  makeDir(outputDir);
  for (const auto& dir : outputDirList) {
    makeDir(dir);
  }
  //
  for (const auto& file : outputFileList) {
    ofstream outputFile;
    outputFile.open(file.c_str());
    for (const auto& row : outputFileRows.at(file)) {
      outputFile << row << std::endl;
    }
    outputFile.close();
  } 
};

bool readFile(std::vector<std::string>& rows, const std::string& fileName)
{
  std::cout << "[INFO] Reading file: " << fileName << std::endl;
  ifstream myfile(fileName.c_str());
  if (myfile.is_open()){ 
    std::string line, CHAR;
    while (getlineSafe(myfile, line)) {
      std::stringstream row(line), tmp(line); tmp >> CHAR;
      if ((!tmp) || (CHAR.find('#')!=std::string::npos) || (CHAR.find("//")!=std::string::npos)) continue;
      rows.push_back(line);
    }
  }
  else {
    std::cout << "[ERROR] File: " << fileName << " was not found!" << std::endl; return false;
  }
  return true;
};

#endif // #ifndef makeInputForLLR_C
