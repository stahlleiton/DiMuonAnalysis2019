#ifndef printLLRStudy_C
#define printLLRStudy_C

#include "fitter.C"
#include "../Results/Utilities/bin.h"


typedef std::map< std::string , GlobalInfo > GlobalInfoMap_t;
typedef std::map< anabin , GlobalInfoMap_t > GlobalInfoDiMap_t;
typedef std::map< anabin , std::string > BinModelMap_t;


const StringVector_t OBSLIST_({"Cand_Rap", "Cand_RapCM", "Cand_AbsRap", "Cand_Pt", "NTrack"});


bool findDir            ( StringVectorMap_t& DIR , const std::string& workDirName , const double& pvalcut);
bool readWSFiles        ( GlobalInfoDiMap_t& content , const StringVector_t& fileNames , const std::string& dirPath );
bool inputFileFromLLR   ( const std::string& outputFile , const std::string& inputFile , const BinModelMap_t& winnerLabels );
StringVector_t printNLL ( BinModelMap_t& winnerLabels , const GlobalInfoDiMap_t& content , const std::string& outputDir , const double& pvalcut );


void printLLRStudy(
		   const std::string workDirName = "Nominal_LLR", // Working directory
                   const double pvalcut = 5. // cut pvalue, in %  
                   )
{
  //
  const std::string& objTag = "JPsi";
  const std::string& colTag = "PA8Y16";
  const StringVector_t& PDList = {"MINBIAS", "DIMUON"};
  const std::string& channel = "ToMuMu";
  //
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Minimization);
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  //
  // Step 0: Set all input and output directories
  StringVectorMap_t DIR;
  if(!findDir(DIR, workDirName, pvalcut)){ return; }
  //
  // Step 1: Read all fit output files
  for(uint j = 0; j < DIR.at("output").size(); j++) {
    if (DIR.at("output").size()>1 && j==0) continue; // First entry is always the main input directory
    //
    BinModelMap_t winnerLabels;
    for(const auto& PD : PDList) {
      // Set the fit output directory
      const std::string& dirPath = Form("%sCandMass/DATA_%s/%s/%s/", DIR.at("output")[j].c_str(), PD.c_str(), objTag.c_str(), colTag.c_str());
      if (!existDir(dirPath)) continue;
      // Find the fit output files
      StringVector_t fileNames;
      if (!fileList(fileNames, dirPath+"result/", false)) { return; }
      // Read the fit output files
      GlobalInfoDiMap_t content;
      if (!readWSFiles(content, fileNames, dirPath+"result/")) { return; }
      //
      // Perform the LLR test
      const std::string& outputDir = Form("%sCandMass/DATA_%s/%s/%s/", DIR.at("outputLLR")[j].c_str(), PD.c_str(), objTag.c_str(), colTag.c_str());
      for (const auto& s : StringVector_t({"result", "plot/CandMass/png", "plot/CandMass/pdf", "plot/CandMass/root", "plot/CandMass/C", "tex"})) {
	makeDir(outputDir+s);
      }
      const auto& bestModelFiles = printNLL(winnerLabels, content, outputDir, pvalcut); 
      std::cout << "[INFO] " << "Background study summary file done for PD: " << PD << "!" << std::endl;
      //
      // Copy all the fit output files for the best models
      std::cout << "[INFO] The files for the best models are: " << std::endl;
      for (auto bestModelFile : bestModelFiles) {
	std::cout << bestModelFile << std::endl;
	stringReplace(bestModelFile, ".root", "");
	// Copy the fit output files
	gSystem->CopyFile((dirPath+"result/"+bestModelFile+".root").c_str(), (outputDir+"result/"+bestModelFile+".root").c_str());
	// Copy the fit output plots
	stringReplace(bestModelFile, "FIT_", "PLOT_");
	for (const auto& s : StringVector_t({"png", "pdf", "root", "C"})) {
	  std::cout << "[INFO] Copying file: " << (dirPath+"plot/CandMass/"+s+"/"+bestModelFile+"."+s) << " to " << (outputDir+"plot/CandMass/"+s+"/"+bestModelFile+"."+s) << std::endl;
	  gSystem->CopyFile((dirPath+"plot/CandMass/"+s+"/"+bestModelFile+"."+s).c_str(), (outputDir+"plot/CandMass/"+s+"/"+bestModelFile+"."+s).c_str());
	}
      }
    }
    //
    // Create the input files for LLR
    const auto& dir = DIR.at("input")[j];
    std::string inputFile = dir+"InitialParam_Cand_Mass_"+objTag+"_"+colTag+".csv";
    auto outputFile = inputFile; stringReplace(outputFile, DIR["input"][j], DIR["inputLLR"][j]);
    inputFileFromLLR(outputFile, inputFile, winnerLabels);
  }
  return;
};


bool findDir(StringVectorMap_t& DIR, const std::string& workDirName, const double& pvalcut)
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
    stringReplace(dir, "/Input/", Form("/Input/LLR_%.0f/", pvalcut*10.));
    DIR["inputLLR"].push_back(dir);
    makeDir(DIR.at("inputLLR").back());
    dir = inDir;
    stringReplace(dir, "/Input/", Form("/Output/LLR_%.0f/", pvalcut*10.));
    DIR["outputLLR"].push_back(dir);
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
  // Get PDF info
  const std::string& pdfName = (ws->obj("pdfName") ? dynamic_cast<RooStringVar*>(ws->obj("pdfName"))->getVal() : "");
  if (!ws->pdf(pdfName.c_str())) { std::cout << "[ERROR] extractInfo: PDF " << pdfName << " was not found!" << std::endl; inputFile.Close(); return false; }
  auto pars = std::unique_ptr<RooArgSet>(ws->pdf(pdfName.c_str())->getParameters(*ws->var("Cand_Mass")));
  if (!pars) { std::cout << "[ERROR] extractNLL: Parameters for " << pdfName << " was not found in " << fileName << std::endl; inputFile.Close(); return false; }
  // Extract the NLL info
  if (!contain(info.Var, "NLL_"+pdfName)) { std::cout << "[ERROR] extractNLL: NLL was not found in " << fileName << std::endl; inputFile.Close(); return false; }
  info.Var["npar"]["Val"] = pars->getSize();
  info.Var["nll"]["Val"] = info.Var.at("NLL_"+pdfName).at("Val");
  // Define bin
  for (const auto& obs : OBSLIST_) {
    if (obs=="Cand_Rap" && (contain(info.Var, "Cand_RapCM") || contain(info.Var, "Cand_AbsRap"))) continue;
    if (contain(info.Var, obs)) { bin.setbin(obs, info.Var.at(obs).at("Min"), info.Var.at(obs).at("Max")); }
  }
  // Return
  inputFile.Close();
  return true;
};


bool readWSFiles(GlobalInfoDiMap_t& content, const StringVector_t& fileNames, const std::string& dirPath)
{
  for (const auto& fileName : fileNames) {
    std::cout << "[INFO] Processing file: " << fileName << std::endl;
    //
    anabin bin;
    GlobalInfo info;
    const auto& str0 = fileName.substr(0, fileName.rfind("Bkg_"));
    const auto& str1 = fileName.substr(fileName.rfind("Bkg_")+4);
    auto modelName = str0.substr(str0.rfind("_")+1);
    if (modelName=="Uniform") { modelName = "Cheb0"; } // Temporary lazy solution
    // Extract the info
    info.Par["fileName"] = fileName;
    info.Par["modelName"] = modelName;
    info.Var["cnt"]["Val"] = 0;
    if (info.Par.at("modelName")=="") continue;
    if (extractInfo(info, bin, dirPath+fileName)) { content[bin][modelName] = info; }
    else { return false; }
  }
  if (content.size()==0) {
    std::cout << "[ERROR] readFiles: Workspace was not found in the input files" << std::endl; return false;
  }
  return true;
};


void setLines(StringVector_t& strLin, const StringVector_t& lin) 
{
  const std::string& empty  = "                          ";
  if (strLin.size() < lin.size()) { strLin = lin; }
  else {
    for (uint i = 0; i < lin.size(); i++) { strLin[i] += lin[i]; }
  }
  for (uint i = 0; i < lin.size(); i++) { strLin[i].append(empty, 0, (25-lin[i].length())); }
};


StringVector_t printNLL(BinModelMap_t& winnerLabels, const GlobalInfoDiMap_t& content, const std::string& outputDir, const double& pvalcut)
{
  StringVector_t ans;
  //
  const std::string& outputFile = Form("%s/LLRTest_Bkg.txt", outputDir.c_str());
  ofstream fout(outputFile);
  const std::string& outputFileTexTable = Form("%s/tex/LLRTest_Bkg_TexTables.txt", outputDir.c_str());
  ofstream foutTexTable(outputFileTexTable);
  // Loop over the bins
  for (const auto& b : content) {
    const auto& binName = b.first.binInfo();
    const auto& binCont = b.second;
    std::cout << "Analzing Kinematic Bin: " << binName << std::endl; std::cout << " " << std::endl;
    fout << "Analzing Kinematic Bin: " << binName << endl; fout << " " << endl;
    // Loop over the models
    auto tmpCont = binCont;
    for (auto& modRow : tmpCont) {
      auto& modelRow = modRow.second;
      const auto& modelNameB = modelRow.Par.at("modelName");
      const auto& nParB      = modelRow.Var.at("npar").at("Val");
      const auto& NLLB       = modelRow.Var.at("nll").at("Val");
      const auto& AICB       = 2.0*(nParB + NLLB);
      //
      StringVector_t strLin;
      for (const auto& modCol : binCont) {
	const auto& modelCol   = modCol.second;
        if (modelCol.Var.at("npar").at("Val")>=nParB) {
	  const auto& modelNameA = modelCol.Par.at("modelName");
	  const auto& nParA      = modelCol.Var.at("npar").at("Val");
	  const auto& NLLA       = modelCol.Var.at("nll").at("Val");
	  const auto& AICA       = 2.0*(nParA + NLLA);
	  //
	  StringVector_t lin;
          if (modelNameA==modelNameB) {
            lin.push_back("|  TEST: "+modelNameA);
            lin.push_back(Form("|    NPar: %.0f  ", nParA));
            lin.push_back(Form("|    NLL: %.2f  ", NLLA));
            lin.push_back(Form("|    AIC: %.2f  ", AICA));
            lin.push_back("|     ");
            lin.push_back("|     ");
            lin.push_back("|------------------------");
          }
	  else if (nParA >= nParB) {
            const double& diffNLL  = -2.0*(NLLA - NLLB);
            const double& diffNPar =  2.0*(nParA-nParB);
            double probChi2 = 100.0*TMath::Prob(diffNLL, diffNPar);
            if (diffNLL<0) probChi2 = 100.0;
            if (probChi2>pvalcut && (nParA-nParB)<=2.0) modelRow.Var.at("cnt").at("Val")++;
            lin.push_back("|   "+modelNameA);
            lin.push_back(Form("|    NPar: %.0f  ", nParA));
            lin.push_back(Form("|    NLL: %.2f  ", NLLA));
            lin.push_back(Form("|    Diff: %.2f  ", diffNLL));
            lin.push_back(Form("|    Prob: %.1f%s   ", probChi2, "%"));
            lin.push_back(Form("|    AIC: %.2f  ", -(AICA-AICB)));
            lin.push_back("|------------------------");
          }
          setLines(strLin, lin);
        }
      }
      for (auto& l : strLin) { l += "|"; }
      for (const auto& l : strLin) { std::cout << l << std::endl; fout << l << endl; }
    }
    // which is the best model for this bin?
    BinModelMap_t bestModelLabels;
    std::string bestModelFile="NOTFOUND", bestModelName=""; double minok=999., maxpar=0.;
    for (const auto& modRow : tmpCont) {
      const auto& modelRow = modRow.second;
      maxpar = std::max(maxpar, modelRow.Var.at("npar").at("Val"));
      if (modelRow.Var.at("cnt").at("Val")>=2 && modelRow.Var.at("npar").at("Val")<minok) {
	bestModelFile = modelRow.Par.at("fileName");
	bestModelName = modelRow.Par.at("modelName");
	winnerLabels[b.first] = bestModelName;
	minok = modelRow.Var.at("npar").at("Val");
      }
    }
    if (minok==999.) { // sometimes the best model is one of the two highest orders...
      for (const auto& modRow : tmpCont) {
	const auto& modelRow = modRow.second;
	const auto& npar = modelRow.Var.at("npar").at("Val");
	if (modelRow.Var.at("cnt").at("Val")>=maxpar-npar && npar<minok) {
	  bestModelFile = modelRow.Par.at("fileName")+" WARNING, HIGH ORDER";
	  bestModelName = modelRow.Par.at("modelName");
	  winnerLabels[b.first] = bestModelName;
	  minok = npar;
	}
      }
    }
    //
    std::cout << std::endl << " And the winner is... " << bestModelFile << std::endl << std::endl << std::endl;
    fout << std::endl << " And the winner is... " << bestModelFile << std::endl << std::endl << std::endl;
    ans.push_back(bestModelFile);
  }
  // Return
  return ans;
};


bool inputFileFromLLR(const std::string& outputFileName, const std::string& inputFile, const BinModelMap_t& winnerLabels)
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
    // Define the output row
    auto oRow = row.at("ROW");
    if (!contain(winnerLabels, bin)) { std::cout << "[ERROR] Bin not found in winnerLabels! "; bin.print(); return false; }
    auto model = winnerLabels.at(bin);
    if (model=="Cheb0") { model = "Uniform"; }
    else { stringReplace(model, "Cheb", "Chebychev"); }
    if (oRow.find(model)!=std::string::npos) { outputFileRows.push_back(oRow); }
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


#endif // #ifndef printLLRStudy_C
