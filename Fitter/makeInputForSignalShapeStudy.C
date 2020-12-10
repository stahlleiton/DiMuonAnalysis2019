#ifndef makeInputForSignalShapeStudy_C
#define makeInputForSignalShapeStudy_C

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


void makeInputForSignalShapeStudy(const std::string& workDirName,
				  const std::string& par="JPsi",
				  const std::string& col="PA8Y16")
{
  //
  const std::vector<std::string> SIGMODELS = {"SingleCrystalBall", "SingleExtCrystalBall", "DoubleCrystalBall", "DoubleExtCrystalBall", "GaussianAndCrystalBall", "GaussianAndExtCrystalBall"};
  //
  const std::string& CWD = getcwd(NULL, 0);
  const auto& inputDir = CWD + "/Input/" + workDirName + "/";
  if (existDir(inputDir)==false){ 
    std::cout << "[ERROR] Input directory: " << inputDir << " doesn't exist!" << std::endl;
    return; 
  }
  std::vector<std::string> inputDirList({inputDir});
  findSubDir(inputDirList, inputDir);
  //
  std::vector<std::string> inputFileList;
  for (const auto& dir : inputDirList) {
    inputFileList.push_back(dir+"InitialParam_Cand_Mass_"+par+"_"+col+".csv");
  }
  //
  std::map<std::string, std::vector<std::string> > inputFileRows;
  for (const auto& file : inputFileList) {
    std::cout << file << std::endl;
    readFile(inputFileRows[file], file);
  }
  //
  auto outputDir = inputDir;
  stringReplace(outputDir, workDirName, workDirName+"_SignalShapeStudy");
  std::vector<std::string> outputDirList;
  for (const auto& dir : inputDirList) {
    outputDirList.push_back(dir);
    stringReplace(outputDirList.back(), workDirName, workDirName+"_SignalShapeStudy");
  }
  //
  std::vector<std::string> outputFileList;
  std::map<std::string, std::vector<std::string> > outputFileRows;
  for (const auto& file : inputFileList) {
    outputFileList.push_back(file);
    stringReplace(outputFileList.back(), workDirName, workDirName+"_SignalShapeStudy");
    std::cout << "A " << file << std::endl;
    for (const auto& row : inputFileRows.at(file)) {
      for (const auto& sigModel : SIGMODELS) {
	auto outputRow = row;
	std::string model = "";
	if (row.rfind("SingleExtCrystalBall")!=std::string::npos) { model = "SingleExtCrystalBall"; }
	if (model!="") {
	  stringReplace(outputRow, model, sigModel);
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
    std::cout << "B " << file << "  " << contain(outputFileRows, file) <<  std::endl;
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

#endif // #ifndef makeInputForSignalShapeStudy_C
