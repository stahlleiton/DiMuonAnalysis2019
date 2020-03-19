#include "fitter.C"


typedef std::map< std::string , GlobalInfo      > GlobalInfoMap_t;
typedef std::map< std::string , GlobalInfoMap_t > GlobalInfoDiMap_t;


bool           findDir   ( StringVectorMap_t& DIR , const std::string& workDirName , const double& pvalcut );
bool           readFiles ( GlobalInfoDiMap_t& content , const StringVector_t& fileNames , const std::string& dirPath );
StringVector_t printNLL  ( StringDiVector_t& winnerLabels , const GlobalInfoDiMap_t& content , const std::string& outputDir , const double& pvalcut );
bool           reduceInputFile ( const std::string& outputFile , const std::string& inputFile , const StringDiVector_t& winnerLabels );


void printLLRStudy(
		   const std::string workDirName = "Nominal_Inclusive_AllBkg", // Working directory
                   const double pvalcut = 5. // cut pvalue, in %  
                   )
{
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
    StringDiVector_t winnerLabels;
    for(const auto& PD : PDList) {
      // Set the fit output directory
      const std::string& dirPath = Form("%sCandMass/DATA_%s/%s/%s/", DIR.at("output")[j].c_str(), PD.c_str(), objTag.c_str(), colTag.c_str());
      if (!existDir(dirPath)) continue;
      // Find the fit output files
      StringVector_t fileNames;
      if (!fileList(fileNames, dirPath+"result/", false)) { return; }
      // Read the fit output files
      GlobalInfoDiMap_t content;
      if (!readFiles(content, fileNames, dirPath+"result/")) { return; }
      // Perform the LLR test
      const std::string& outputDir = Form("%sCandMass/DATA_%s/%s/%s/", DIR.at("outputLLR")[j].c_str(), PD.c_str(), objTag.c_str(), colTag.c_str());
      for (const auto& s : StringVector_t({"result", "plot/CandMass/png", "plot/CandMass/pdf", "plot/CandMass/root", "plot/CandMass/C", "tex"})) {
	makeDir(outputDir+s);
      }
      makeDir(DIR.at("inputLLR")[j]);
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
    // Reduce the input files
    BoolDiMap_t VARMAP; VARMAP["Cand_Mass"][objTag] = true;
    BoolMap_t COLMAP; COLMAP[colTag] = true;
    for (const auto& VAR : VARMAP) {
      auto PARMAP = VAR.second;
      for (const auto& PAR : PARMAP) {
	if (PAR.second) {
	  for (const auto& COL : COLMAP) {
	    if(COL.second) {
	      const auto& dir = DIR.at("input")[j];
	      std::string inputFile = "", name = (dir + "InitialParam_" + VAR.first + "_" + PAR.first);
	      StringVector_t tryChannel = { channel , "" };
	      StringVector_t trySystem  = { COL.first };
	      if (colTag.rfind("8Y16")!=std::string::npos) { trySystem.push_back("PA8Y16"); }; trySystem.push_back("");
	      for (const auto& tryCha : tryChannel) {
		bool trySuccess = false;
		for (const auto& tryCol : trySystem) {
		  if (ifstream(inputFile).good()==false) { inputFile = (name + tryCha + (tryCol!="" ? "_"+tryCol : "") + ".csv"); } else { trySuccess = true; break; }
		}
		if (trySuccess) break;
	      }
	      auto outputFile = inputFile; stringReplace(outputFile, DIR["input"][j], DIR["inputLLR"][j]);
	      reduceInputFile(outputFile, inputFile, winnerLabels);
	    }
	  }
	}
      }
    }
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
    dir = inDir;
    stringReplace(dir, "/Input/", Form("/Output/LLR_%.0f/", pvalcut*10.));
    DIR["outputLLR"].push_back(dir);
  }
  return true;
};


bool extractNLL(GlobalInfo& info, const std::string& fileName)
{
  // Open the input file
  TFile inputFile(fileName.c_str(), "READ");
  if (!inputFile.IsOpen() || inputFile.IsZombie()) { std::cout << "[ERROR] The input file " << fileName << " was not opened!" << std::endl; return false; }
  inputFile.cd();
  // Extract the workspace
  const auto& ws = dynamic_cast<RooWorkspace*>(inputFile.Get("workspace"));
  if (!ws) { std::cout << "[ERROR] Workspace not found in " << fileName << std::endl; inputFile.Close(); return false; }
  // Get the dataset variables and PDF names
  const std::string& dsName  = (ws->obj("dsName")  ? dynamic_cast<RooStringVar*>(ws->obj("dsName"))->getVal() : "");
  const std::string& pdfName = (ws->obj("pdfName") ? dynamic_cast<RooStringVar*>(ws->obj("pdfName"))->getVal() : "");
  if (!ws->pdf(pdfName.c_str())) { std::cout << "[ERROR] extractNLL: DataSet " << pdfName << " was not found!" << std::endl; inputFile.Close(); return false; }
  // Extract the NLL info
  auto nll = ws->var(("NLL_"+pdfName).c_str());
  if (!nll) { std::cout << "[ERROR] extractNLL: NLL was not found in " << fileName << std::endl; inputFile.Close(); return false; }
  const auto& valNLL = nll->getVal();
  auto pars = std::unique_ptr<RooArgSet>(ws->pdf(pdfName.c_str())->getParameters(*ws->var("Cand_Mass")));
  if (!pars) { std::cout << "[ERROR] extractNLL: Parameters for " << pdfName << " was not found in " << fileName << std::endl; inputFile.Close(); return false; }
  const auto& npar = pars->getSize();
  // Store the NLL info
  info.Var["npar"]["Val"] = npar;
  info.Var["nll"]["Val"] = valNLL;
  // Extract the observable's info
  StringSet_t obsNames;
  const auto& listObs = ws->set(("SET_"+dsName).c_str());
  if (listObs) {
    auto obsIt = std::unique_ptr<TIterator>(listObs->createIterator());
    for (auto it = obsIt->Next(); it!=NULL; it = obsIt->Next()) { obsNames.insert(it->GetName()); }
  }
  else { obsNames = StringSet_t({"Cand_Mass", "Cand_Pt", "Cand_Rap", "Cand_DLen", "Centrality", "NTrack"}); } // BUG FIX
  //
  for (const auto& obs : obsNames) {
    auto name = obs;
    if (name=="Cand_Rap" && ws->var("Cand_RapCM")) { name = "Cand_RapCM"; }
    else if (name=="Cand_Rap" && ws->var("Cand_AbsRap")) { name = "Cand_AbsRap"; }
    else if (name=="Cand_Mass") continue;
    auto var = ws->var(name.c_str());
    if (var->getMin()!=var->getMin("DEFAULT") || var->getMax()!=var->getMax("DEFAULT")) {
      const auto vMin = var->getMin();
      const auto vMax = var->getMax();
      std::string lbl = ""; 
      lbl += fmod(vMin*10., 1.)==0 ? Form("%.1f", vMin) : Form("%g", vMin);
      lbl += (name.rfind("Cand_Rap",0)==0 ? "_" : "-");
      lbl += fmod(vMax*10., 1.)==0 ? Form("%.1f", vMax) : Form("%g", vMax);
      info.StrV["OBS"].push_back(lbl);
    }
  }
  // Return
  inputFile.Close();
  return true;
};


bool readFiles(GlobalInfoDiMap_t& content, const StringVector_t& fileNames, const std::string& dirPath)
{
  for (const auto& fileName : fileNames) {
    std::cout << "[INFO] Processing file: " << fileName << std::endl;
    //
    GlobalInfo info;
    // Get info from the fileName
    const auto& str0 = fileName.substr(0, fileName.rfind("Bkg_"));
    const auto& str1 = fileName.substr(fileName.rfind("Bkg_")+4);
    // Store the info
    info.Par["fileName"] = fileName;
    info.Par["modelName"] = str0.substr(str0.rfind("_")+1);
    if (info.Par.at("modelName")=="Uniform") { info.Par.at("modelName") = "Cheb0"; } // Temporary lazy solution
    info.Par["binName"] = str1.substr(0, str1.rfind(".root"));
    info.Var["cnt"]["Val"] = 0;
    //
    if (info.Par.at("modelName")=="") continue;
    if (extractNLL(info, dirPath+fileName)) { content[info.Par.at("binName")][info.Par.at("modelName")] = info; }
    else { return false; }
  }
  if (content.size()==0) {
    std::cout << "[ERROR] readFiles: No NLL values were found in the input files" << std::endl; return false;
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


StringVector_t printNLL(StringDiVector_t& winnerLabels, const GlobalInfoDiMap_t& content, const std::string& outputDir, const double& pvalcut)
{
  StringVector_t ans;
  //
  const std::string& outputFile = Form("%s/LLRTest_Bkg.txt", outputDir.c_str());
  ofstream fout(outputFile);
  const std::string& outputFileTexTable = Form("%s/tex/LLRTest_Bkg_TexTables.txt", outputDir.c_str());
  ofstream foutTexTable(outputFileTexTable);
  // Loop over the bins
  for (const auto& b : content) {
    const auto& binName = b.first;
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
    StringVector_t bestModelLabels;
    std::string bestModelFile="NOTFOUND", bestModelName=""; double minok=999., maxpar=0.;
    for (const auto& modRow : tmpCont) {
      const auto& modelRow = modRow.second;
      maxpar = std::max(maxpar, modelRow.Var.at("npar").at("Val"));
      if (modelRow.Var.at("cnt").at("Val")>=2 && modelRow.Var.at("npar").at("Val")<minok) {
	bestModelFile = modelRow.Par.at("fileName");
	bestModelName = modelRow.Par.at("modelName");
	bestModelLabels = modelRow.StrV.at("OBS");
	bestModelLabels.push_back(bestModelName);
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
	  bestModelLabels = modelRow.StrV.at("OBS");
	  bestModelLabels.push_back(bestModelName);
	  minok = npar;
	}
      }
    }
    //
    std::cout << std::endl << " And the winner is... " << bestModelFile << std::endl << std::endl << std::endl;
    fout << std::endl << " And the winner is... " << bestModelFile << std::endl << std::endl << std::endl;
    ans.push_back(bestModelFile);
    winnerLabels.push_back(bestModelLabels);
  }
  // Return
  return ans;
};


bool readInputFile(StringVector_t& content, const std::string& fileName)
{
  ifstream myfile(fileName.c_str());
  if (myfile.is_open()){ 
    std::string line;
    while (getline(myfile, line)){ content.push_back(line); }
  } else {
    std::cout << "[ERROR] File: " << fileName << " was not found!" << std::endl; return false;
  }
  return true;
};


bool reduceInputFile(const std::string& outputFile, const std::string& inputFile, const StringDiVector_t& winnerLabels)
{
  std::cout << "[INFO] Processing input file: " << inputFile << std::endl;
  // Open the input file and extract its content
  StringVector_t inputContent;
  if(!readInputFile(inputContent, inputFile)){ return false; }
  // Loop over the input file and determine the lines to keep
  StringVector_t outputContent;
  for(const auto& row : inputContent) {
    bool found = false;
    if (&row == &inputContent.front()) { found = true; }
    else {
      for(const auto& winLbl : winnerLabels) {
	bool foundLbl = true;
	for (const auto& lbl : winLbl) {
	  auto label = lbl;
	  stringReplace(label, "Cheb", "Chebychev");
	  if (row.find(label)==std::string::npos) { foundLbl = false; break; }
	}
	if (foundLbl) { found = true; break; }
      }
    }
    if (found) { outputContent.push_back(row); }
  }
  // Save the reduced input file
  if (outputContent.size()>0) {
    ofstream myfile(outputFile.c_str());
    if (myfile.is_open()){ for (const auto& line : outputContent) { myfile << line << std::endl; } }
    std::cout << "[INFO] Reduced input file stored in: " << outputFile << std::endl;
  }
  return true;
};
