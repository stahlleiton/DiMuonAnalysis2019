#ifndef makeInputForBestModel_C
#define makeInputForBestModel_C

#include "fitter.C"
#include "../Results/Utilities/bin.h"


typedef std::map< anabin , GlobalInfo > GlobalInfoMap_t;
typedef std::vector<std::tuple<std::string, double, double>> VarMap_t;

const StringVector_t OBSLIST_({"Cand_Rap", "Cand_RapCM", "Cand_AbsRap", "Cand_Pt"});
const VarMap_t FIXVAR_({{"nR", 1.0, 49.0}, {"Alpha", 1.0, 2.9}, {"n", 1.0, 3.4}});
const VarMap_t FREEVAR_;

bool findDir         ( StringVectorMap_t& DIR , const std::string& workDirName , const std::string& workDirName_Fit , const std::string& dirLbl );
bool readWSFiles     ( GlobalInfoMap_t& content , const StringVector_t& fileNames , const std::string& dirPath );
bool inputFileForBestModel ( const std::string& outputFile , const std::string& inputFile , const GlobalInfoMap_t& infoMap , const std::string& objLbl , const bool& fixPar );


void makeInputForBestModel(const std::string& workDirName,
			   const std::string& workDirName_Fit,
			   const bool& fixPar = true,
			   const std::string& objTag="JPsi",
			   const StringVector_t& PDList = {"MINBIAS", "DIMUON"},
			   const std::string& colTag="PA8Y16")
{
  //
  const std::string dirLbl = "_BestModel";
  const std::string objLbl = objTag+"ToMuMuOS_"+colTag;
  //
  //
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Minimization);
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  //
  // Step 0: Set all input and output directories
  StringVectorMap_t DIR;
  if(!findDir(DIR, workDirName, workDirName_Fit, dirLbl)){ return; }
  //
  // Step 1: Read all fit output files
  for(uint j = 0; j < DIR.at("output").size(); j++) {
    if (DIR.at("output").size()>1 && j==0) continue; // First entry is always the main input directory
    if (DIR.at("output")[j]=="") continue;
    std::cout << DIR.at("output")[j] << std::endl;
    //
    GlobalInfoMap_t infoMap;
    for(const auto& PD : PDList) {
      // Set the fit output directory
      const std::string& dirPath = Form("%sCandMass/DATA_%s/%s/%s/", DIR.at("output")[j].c_str(), PD.c_str(), objTag.c_str(), colTag.c_str());
      if (!existDir(dirPath)) continue;
      // Find the fit output files
      StringVector_t fileNames;
      if (!fileList(fileNames, dirPath+"result/", false)) { return; }
      // Read the fit output files
      if (!readWSFiles(infoMap, fileNames, dirPath+"result/")) { return; }
      //
    }
    //
    // Create the input files for LLR
    const auto& dir = DIR.at("input")[j];
    std::string inputFile = dir+"InitialParam_Cand_Mass_"+objTag+"_"+colTag+".csv";
    auto outputFile = inputFile; stringReplace(outputFile, DIR["input"][j], DIR["inputLLR"][j]);
    inputFileForBestModel(outputFile, inputFile, infoMap, objLbl, fixPar);
  }
};


bool findDir(StringVectorMap_t& DIR, const std::string& workDirName, const std::string& workDirName_Fit, const std::string& dirLbl)
{
  std::cout << "[INFO] Searching for input fit directories in: " << workDirName << std::endl;
  const std::string& CWD = getcwd(NULL, 0);
  DIR["main"].push_back(CWD);
  DIR["input"].push_back(DIR.at("main")[0] + "/Input/" + workDirName + "/");
  if (!existDir(DIR.at("input")[0])){
    std::cout << "[ERROR] Input directory: " << DIR.at("input")[0] << " doesn't exist!" << std::endl; return false;
  }
  else {
    findSubDir(DIR.at("input"), DIR.at("input")[0]);
  }
  for (const auto& inDir : DIR.at("input")) {
    auto dir = inDir;
    stringReplace(dir, workDirName, workDirName+dirLbl);
    DIR["inputLLR"].push_back(dir);
    makeDir(DIR.at("inputLLR").back());
  }
  for (const auto& inDir : DIR.at("input")) {
    auto dir = inDir;
    if (dir.find("/NTrack_")!=std::string::npos) {  dir = dir.substr(0, dir.rfind("/NTrack_"))+"/NTrack_15_250/"; }
    stringReplace(dir, workDirName, workDirName_Fit);
    stringReplace(dir, "/Input/", "/Output/");
    DIR["output"].push_back(existDir(dir) ? dir : "");
  }
  return true;
};


bool extractInfo(GlobalInfo& info, anabin& bin, const std::string& fileName)
{
  // Open the input file
  TFile inputFile(fileName.c_str(), "READ");
  if (!inputFile.IsOpen() || inputFile.IsZombie()) { std::cout << "[ERROR] The input file " << fileName << " was not opened!" << std::endl; return false; }
  inputFile.cd();
  // Extract the workspace
  const auto& ws = dynamic_cast<RooWorkspace*>(inputFile.Get("workspace"));
  if (!ws) { std::cout << "[ERROR] Workspace not found in " << fileName << std::endl; inputFile.Close(); return false; }
  // Get variables
  auto vars = ws->allVars();
  auto varIt = std::unique_ptr<TIterator>(vars.createIterator());
  for (auto itp = varIt->Next(); itp!=NULL; itp = varIt->Next()) {
    const auto& it = dynamic_cast<RooRealVar*>(itp); if (!it) continue;
    const std::string& name = it->GetName();
    info.Var[name]["Min"] = it->getMin();
    info.Var[name]["Max"] = it->getMax();
    info.Var[name]["Val"] = it->getVal();
  }
  // Define bin
  for (const auto& obs : OBSLIST_) {
    if (obs=="Cand_Rap" && (contain(info.Var, "Cand_RapCM") || contain(info.Var, "Cand_AbsRap"))) continue;
    if (contain(info.Var, obs)) { bin.setbin(obs, info.Var.at(obs).at("Min"), info.Var.at(obs).at("Max")); }
  }
  // Return
  inputFile.Close();
  return true;
};


bool readWSFiles(GlobalInfoMap_t& content, const StringVector_t& fileNames, const std::string& dirPath)
{
  for (const auto& fileName : fileNames) {
    std::cout << "[INFO] Processing file: " << fileName << std::endl;
    //
    anabin bin;
    GlobalInfo info;
    if (extractInfo(info, bin, dirPath+fileName)) { content[bin] = info; }
    else { return false; }
  }
  if (content.size()==0) {
    std::cout << "[ERROR] readFiles: Workspace was not found in the input files" << std::endl; return false;
  }
  return true;
};


bool inputFileForBestModel(const std::string& outputFileName, const std::string& inputFile, const GlobalInfoMap_t& infoMap, const std::string& objLbl, const bool& fixPar)
{
  std::cout << "[INFO] Processing input file: " << inputFile << std::endl;
  //
  // Extract information from input files
  StringMapVector_t inputData;
  if(!parseFile(inputData, inputFile, true)) { return false; }
  //
  // Add extra parameters
  StringVector_t tmpV;
  for (const auto& v : FREEVAR_) { tmpV.push_back(std::get<0>(v)); }
  for (const auto& v : FIXVAR_) { tmpV.push_back(std::get<0>(v)); }
  for (const auto& v : tmpV) {
    const bool isAlt = v.rfind("_")!=std::string::npos;
    const auto var = v + (isAlt ? objLbl.substr(objLbl.rfind("To")) : "_"+objLbl);
    auto& header = inputData.at(0).at("HEADER");
    bool notFound = false;
    if (header.rfind(var)==std::string::npos) {
      header += ","+var;
      notFound = true;
    }
    if (notFound) {
      for (auto& row : inputData) {
	if (!contain(row, "HEADER")) {
	  row[var] = "DUMMY";
	  row.at("ROW") += ","+row.at(var);
	}
      }
    }
  }
  //
  // Define output config file
  StringVector_t outputFileRows;
  outputFileRows.push_back(inputData.at(0).at("HEADER"));
  //
  for (const auto& row : inputData) {
    if (contain(row, "HEADER")) continue;
    anabin bin;
    for (const auto& obs : OBSLIST_) {
      if (contain(row, obs) && row.at(obs)!="NONE") {
        std::vector<double> v;
        if (!parseString(v, row.at(obs))) { return false; }
	bin.setbin(obs, v.at(0), v.at(1));
      }
    }
    //
    // Get the output file information
    if (!contain(infoMap, bin)) { std::cout << "[ERROR] Bin not found! "; bin.print(); return false; }
    const auto& info = infoMap.at(bin);
    //
    // Define the output row
    auto oRow = row.at("ROW");
    for (const auto& vv : FREEVAR_) {
      const auto& v = std::get<0>(vv);
      const bool isAlt = v.rfind("_")!=std::string::npos;
      const auto var = v + (isAlt ? objLbl.substr(objLbl.rfind("To")) : "_"+objLbl);
      if (contain(info.Var, var) && contain(row, var)) { 
	const auto min = std::max(std::get<1>(vv), info.Var.at(var).at("Min"));
	const auto max = std::min(std::get<2>(vv), info.Var.at(var).at("Max"));
	stringReplace(oRow, row.at(var), Form("[%g;%g;%g]", info.Var.at(var).at("Val"), min, max), 1);
      }
    }
    if (fixPar) {
      for (const auto& vv : FIXVAR_) {
	const auto& v = std::get<0>(vv);
	const auto min = std::get<1>(vv);
	const auto max = std::get<2>(vv);
	const bool isAlt = v.rfind("_")!=std::string::npos;
	const auto var = v + (isAlt ? objLbl.substr(objLbl.rfind("To")) : "_"+objLbl);
	const bool fitFail = contain(info.Var, "FIT_FAILED");
	if (contain(info.Var, var) && contain(row, var) && !fitFail) {
	  const bool isGood = (info.Var.at(var).at("Val")>min && info.Var.at(var).at("Val")<max);
	  if (isGood) { stringReplace(oRow, row.at(var), Form("[%.2f]", info.Var.at(var).at("Val")), 1); }
	}
      }
      const auto fVar = "f_"+objLbl;
      const auto fVal = (contain(row,fVar) ? row.at(fVar).substr(0, row.at(fVar).find(";"))+";0.0;1.0]" : "");
      if (contain(info.Var, fVar) && contain(row, fVar)) { stringReplace(oRow, row.at(fVar), fVal, 1); }
    }
    stringReplace(oRow, "DUMMY", "");
    outputFileRows.push_back(oRow);
  }
  //
  // Save output config file
  ofstream outputFile;
  std::cout << "[INFO] Creating input file: " << outputFileName << std::endl;
  outputFile.open(outputFileName.c_str());
  for (const auto& row : outputFileRows) {
    outputFile << row << std::endl;
  }
  outputFile.close();
  //
  return true;
};


#endif // #ifndef makeInputForBestModel_C
