#ifndef makeInputForLLR_C
#define makeInputForLLR_C

#include "fitter.C"
#include "../Results/Utilities/bin.h"


typedef std::map< anabin , GlobalInfo > GlobalInfoMap_t;

const StringVector_t OBSLIST_({"Cand_Rap", "Cand_RapCM", "Cand_AbsRap", "Cand_Pt", "NTrack"});
const StringVector_t SETVAR_({"AlphaR", "rSigma21"});


bool findDir         ( StringVectorMap_t& DIR , const std::string& workDirName , const std::string& dirLbl );
bool readWSFiles     ( GlobalInfoMap_t& content , const StringVector_t& fileNames , const std::string& dirPath );
bool inputFileForLLR ( const std::string& outputFile , const std::string& inputFile , const GlobalInfoMap_t& infoMap , const std::string& objLbl , const StringVector_t& bkgModels , const bool& fixPar );


void makeInputForLLR(const std::string& workDirName,
		      const bool& isNominal = true,
		      const bool& fixPar = true,
		      const std::string& objTag="JPsi",
		      const StringVector_t& PDList = {"MINBIAS", "DIMUON"},
		      const std::string& colTag="PA8Y16")
{
  //
  StringVector_t BKGMODELS;
  if (isNominal) { BKGMODELS = {"Uniform", "Chebychev1", "Chebychev2", "Chebychev3", "Chebychev4", "Chebychev5", "Chebychev6"}; }
  else { BKGMODELS = {"Uniform", "ExpChebychev1", "ExpChebychev2", "ExpChebychev3", "ExpChebychev4", "ExpChebychev5", "ExpChebychev6"}; }
  const std::string dirLbl = (isNominal?"_LLR":"_ExpLLR");
  const std::string objLbl = objTag+"ToMuMuOS_"+colTag;
  //
  //
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Minimization);
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  //
  // Step 0: Set all input and output directories
  StringVectorMap_t DIR;
  if(!findDir(DIR, workDirName, dirLbl)){ return; }
  //
  // Step 1: Read all fit output files
  for(uint j = 0; j < DIR.at("output").size(); j++) {
    if (DIR.at("output").size()>1 && j==0) continue; // First entry is always the main input directory
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
    inputFileForLLR(outputFile, inputFile, infoMap, objLbl, BKGMODELS, fixPar);
  }
};


bool findDir(StringVectorMap_t& DIR, const std::string& workDirName, const std::string& dirLbl)
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
    stringReplace(dir, "/Input/", "/Output/");
    DIR["output"].push_back(dir);
    if (!existDir(dir)) { std::cout << "[ERROR] Output fit directory: " << dir << " doesn't exist!" << std::endl; return false; }
    dir = inDir;
    stringReplace(dir, workDirName, workDirName+dirLbl);
    DIR["inputLLR"].push_back(dir);
    makeDir(DIR.at("inputLLR").back());
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


bool inputFileForLLR(const std::string& outputFileName, const std::string& inputFile, const GlobalInfoMap_t& infoMap, const std::string& objLbl, const StringVector_t& bkgModels, const bool& fixPar)
{
  std::cout << "[INFO] Processing input file: " << inputFile << std::endl;
  //
  // Extract information from input files
  StringMapVector_t inputData;
  if(!parseFile(inputData, inputFile, true)) { return false; }
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
    if (fixPar) {
      for (const auto& v : SETVAR_) {
	const auto var = v + "_" + objLbl;
	if (contain(info.Var, var) && contain(row, var)) { stringReplace(oRow, row.at(var), Form("[%.4f]", info.Var.at(var).at("Val"))); }
      }
      const auto fVar = "f_"+objLbl;
      const auto fVal = (contain(row,fVar) ? row.at(fVar).substr(0, row.at(fVar).find(";"))+";0.0;1.0]" : "");
      if (contain(info.Var, fVar) && contain(row, fVar)) { stringReplace(oRow, row.at(fVar), fVal); }
    }
    std::string model = "";
    if (oRow.rfind("[Bkg]")!=std::string::npos) { model = oRow.substr(0, oRow.rfind("[Bkg]")); }
    if (model.rfind("+")!=std::string::npos) { model = model.substr(model.rfind("+")+1); }
    if (model!="") {
      for (const auto& bkgModel : bkgModels) {
	auto tmpRow = oRow;
	stringReplace(tmpRow, model, bkgModel);
	outputFileRows.push_back(tmpRow);
      }
    }
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


#endif // #ifndef makeInputForLLR_C
