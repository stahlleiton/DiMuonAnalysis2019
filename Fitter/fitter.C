#include "TROOT.h"
#include "TSystem.h"
#include "TSystemDirectory.h"

#include "RooFit.h"
#include "RooMsgService.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <string>
#include <bitset>
#include <algorithm>
#include <future>

#include "../Utilities/dataUtils.h"
#include "Macros/Utilities/initClasses.h"
#include "Macros/tree2DataSet.C"
#include "Macros/Candidate/processDataSet.C"
#include "Macros/Candidate/fitCandidateModel.C"


bool checkSettings     ( const GlobalInfo& userInput );
bool parseFile         ( StringMapVector_t& data , const std::string& fileName );
bool parseString       ( std::vector<double>& output , std::string input );
bool iniWorkEnv        ( StringVectorMap_t& DIR , const std::string& workDirName );
void iniFileDir        ( StringMapVector_t& inputFitDirs , StringDiMapVector_t& inputInitialFilesDirs ,
                         const StringMap_t& inputFitDir  , const StringDiMap_t& inputInitialFilesDir , const StringVectorMap_t& DIR );
bool readFile          ( StringDiVector_t& content, const std::string& fileName , const int& nCol=-1, int nRow=-1 );
bool getInputFileNames ( StringDiVectorMap_t& inputFileCollection , const std::string& inputTrees );
bool setParameters     ( GlobalInfo& info , GlobalInfo& userInfo , const StringMap_t& row );
bool addParameters     ( GlobalInfoVector_t& infoVector , GlobalInfo& userInfo , const std::string& inputFile );
bool createDataSets    ( RooWorkspaceMap_t& Workspace , GlobalInfo& userInput , const StringVectorMap_t& DIR );


void fitter(
            const std::string workDirName = "Test",       // Working directory
            const std::bitset<1> useExt   = 1,            // Use external: (bit 0 (1)) Input DataSets
            // Select the type of datasets to fit
            const std::bitset<2> fitData  = 1,            // Fit Sample: (bit 0 (1)) Data , (bit 1 (2)) MC
            const std::bitset<7> fitColl  = 32,           // Fit System: (bit 0  (1)) PbPb5Y18 ,
                                                          //             (bit 1  (2)) PP5Y17   , (bit 2  (4)) PP13Y18 ,  
                                                          //             (bit 3  (8)) pPb8Y16  , (bit 4 (16)) Pbp8Y16 , (bit 5 (32)) PA8Y16  ,
                                                          //             (bit 6 (64)) PbPb5Y15
            const std::bitset<3> fitChg   = 1,            // Fit Charge: (bit 0 (1)) OS , (bit 1 (2)) SS , (bit 2 (4)) Inclusive
            const unsigned int   usePD    = 0,            // Use PD: (0) MUON , (1) DIMUON , (2) DIMUONPERI ,
                                                          //         (3) HIGHMULT , (4) HIGHMULT2 (100/255) ,
                                                          //         (5) MINBIAS , (6) UPC
            // Select the type of objects to fit
            const std::bitset<8> fitObj   = 8,            // Fit Objects: (bit 0  (1)) Bkg   ,
                                                          //              (bit 1  (2)) JPsi  , (bit 2   (4)) Psi2S ,
                                                          //              (bit 3  (8)) Ups1S , (bit 4  (16)) Ups2S , (bit 5 (32)) Ups3S ,
                                                          //              (bit 6 (64)) Z     , (bit 7 (128)) D0
            // Select the fitting options
            const std::bitset<5> fitVar   = 1,            // Fit Variable: (bit 0 (1)) Cand_Mass    , (bit 1 (2)) Cand_DLen ,
	                                                  //               (bit 2 (4)) Cand_DLenErr , (bit 3 (8)) Cand_DLenRes , (bit 4 (16)) Cand_DLenGen
            const std::bitset<2> fitCat   = 3,            // Fit Objects: (bit 0  (1)) PR , (bit 1  (2)) NoPR
            const unsigned int   numCores = 24,           // Number of cores used for fitting
            const std::string    analysis = "CandToMuMu", // Type of analysis: CandToXX (Mass Resonance)
            // Select the drawing options
            const bool setLogScale  = true                // Draw plot with log scale
            )
{
  //
  // -------------------------------------------------------------------------------
  // STEP 0: INITIALIZE THE FITTER WORK ENVIROMENT
  // The work enviroment is divided as follows:
  /*
    main |-> Macros: Contain all the macros
    |-> Input   |-> <WorkDir> : Contain Input File, Bin and Parameter List for a given work directory (e.g. 20160201)
    |-> Output  |-> <WorkDir> : Contain Output Plots and Results for a given work directory (e.g. 20160201)
    |-> DataSet : Contain all the datasets (MC and Data)
  */
  //
  // Suppress Messages for RooFit
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Caching);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Integration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Minimization);
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  //
  const std::string& CWD = getcwd(NULL, 0);
  GlobalInfo userInput;
  bool saveAll = true;
  //
  // Set all the Parameters from the input settings
  //
  userInput.Par["analysis"] = analysis;
  userInput.Int["numCores"] = numCores;
  userInput.Flag["setLogScale"] = setLogScale;
  //
  // Store more information for fitting
  userInput.Par["extTreesFileDir"] = Form("%s/Input/", CWD.c_str());
  userInput.Par["extDSDir_DATA"]   = "/storage1/users/ags9/DataSet/";
  userInput.Par["extDSDir_MC"]     = "/storage1/users/ags9/DataSet/";
  //
  // Settings for Mass Resonance Analysis
  if (userInput.Par.at("analysis").rfind("CandTo",0)==0) {
    userInput.Par["treeType"] = "VertexCompositeTree";
  }
  //
  // Set all the Boolean Flags from the input settings
  //
  // Check if use external information
  userInput.StrV["setExt"] = {"ExtDS"};
  for (uint i=0; i<userInput.StrV.at("setExt").size(); i++) { userInput.Flag["use"+userInput.StrV.at("setExt")[i]] = useExt[i]; }
  //
  // Type of data to fit
  userInput.StrV["sample"] = {"Data", "MC"};
  for (uint i=0; i<userInput.StrV.at("sample").size(); i++) { userInput.Flag["fit"+userInput.StrV.at("sample")[i]] = fitData[i]; }
  //
  // Dataset to use (ignored if fitting MC)
  userInput.StrV["PD"] = {"MUON", "DIMUON", "DIMUONPERI", "HIGHMULT", "HIGHMULT2", "MINBIAS", "UPC"};
  userInput.Par["PD"] = userInput.StrV.at("PD")[usePD];
  //
  // Collision system
  userInput.StrV["system"] = {"PbPb5Y18", "PP5Y17", "PP13Y18", "pPb8Y16", "Pbp8Y16", "PA8Y16", "PbPb5Y15"};
  for (uint i=0; i<userInput.StrV.at("system").size(); i++) { userInput.Flag["fit"+userInput.StrV.at("system")[i]] = fitColl[i]; }
  for (uint i=0; i<userInput.StrV.at("system").size(); i++) { userInput.Flag[ "do"+userInput.StrV.at("system")[i]] = fitColl[i]; }
  if (userInput.Flag.at("fitpPb8Y16") || userInput.Flag.at("fitPbp8Y16")) { userInput.Flag.at("doPA8Y16") = true; }
  if (userInput.Flag.at("fitPA8Y16")) { userInput.Flag.at("dopPb8Y16") = true; userInput.Flag.at("doPbp8Y16") = true; }
  //
  // Decay channel
  userInput.StrV["channel"] = { "DiMuon", "PionKaon" };
  const StringVector_t chaLabel = { "MuMu", "PiKa" };
  for (uint iCha = 0; iCha < userInput.StrV.at("channel").size(); iCha++) {
    const auto& chaV = userInput.StrV.at("channel")[iCha];
    const auto& chaL = chaLabel[iCha];
    if (userInput.Par.at("analysis").rfind(chaL)!=std::string::npos) {
      auto CHA = chaV; std::transform(CHA.begin(), CHA.end(), CHA.begin(), ::toupper);
      userInput.Flag["do"+chaV]   = true;
      userInput.Par["channelDS"]  = CHA;
      userInput.Par["channelDir"] = chaV;
      userInput.Par["channel"]    = ("To" + chaL);
    }
    else { userInput.Flag["do"+chaV] = false; }
  }
  //
  // Decay charge
  if (userInput.Par.at("analysis").rfind("CandTo",0)==0) {
    userInput.StrV["charge"] = {"OS", "SS", "ChgInc"};
  }
  for (uint i=0; i<userInput.StrV.at("charge").size(); i++) { userInput.Flag["fit"+userInput.StrV.at("charge")[i]] = fitChg[i]; }
  for (const auto& chg : userInput.StrV.at("charge")) {
    std::string label = ""; if (chg=="ChgInc") { label = ""; } else { label = chg.substr(0,2); }
    if (userInput.Flag.at("fit"+chg)) { userInput.StrS["fitCharge"].insert(label); }
  }
  //
  // Fit variable
  userInput.StrV["variable"] = {"Cand_Mass", "Cand_DLen", "Cand_DLenErr", "Cand_DLenRes", "Cand_DLenGen"};
  for (uint i=0; i<userInput.StrV.at("variable").size(); i++) { userInput.Flag["fit"+userInput.StrV.at("variable")[i]] = fitVar[i]; }
  for (const auto& v : userInput.StrV.at("variable")) {
    if (userInput.Flag.at("fit"+v)) {
      userInput.StrS["fitVariable"].insert(v);
      userInput.StrS["incVariable"].insert(v);
      auto var = v; stringReplace(var, "_", "");
      userInput.StrS["fitVarName"].insert(var);
      userInput.StrS["incVarName"].insert(var);
      if (v=="Cand_DLenErr" || v=="Cand_DLenRes") {
	userInput.StrS["incVariable"].insert("Cand_Mass");
	userInput.StrS["incVarName"].insert("CandMass");
      }
    }
  }
  //
  // Type of candidates
  userInput.StrV["object"] = {"Bkg", "JPsi", "Psi2S", "Ups1S", "Ups2S", "Ups3S", "Z", "D0"};
  for (uint i=0; i<userInput.StrV.at("object").size(); i++) {
    userInput.Flag["fit"+userInput.StrV.at("object")[i]] = fitObj[i];
    userInput.Flag["inc"+userInput.StrV.at("object")[i]] = fitObj[i];
  }
  for (const auto& obj : userInput.StrV.at("object")) { if (userInput.Flag.at("fit"+obj)) { userInput.StrS["fitObject"].insert(obj); } }
  for (const auto& var : userInput.StrS.at("fitVarName")) {
    userInput.StrS["incObject_"+var].clear();
    userInput.StrS["template_"+var].clear();
  }
  userInput.StrV["category"] = {"PR", "NoPR"};
  for (uint i=0; i<userInput.StrV.at("category").size(); i++) {
    userInput.Flag["fit"+userInput.StrV.at("category")[i]] = fitCat[i];
    userInput.Flag["inc"+userInput.StrV.at("category")[i]] = fitCat[i];
  }
  //  
  if (userInput.Par.at("analysis").rfind("CandTo",0)==0) {
    double massWidth = 0.1;
    if (userInput.Flag.at("fitJPsi") || userInput.Flag.at("fitPsi2S")) { massWidth = 0.01; }
    else if (userInput.Flag.at("fitUps1S") || userInput.Flag.at("fitUps2S") || userInput.Flag.at("fitUps3S")) { massWidth = 0.05; }
    else if (userInput.Flag.at("fitZ")) { massWidth = 1.0; }
    userInput.Var["Cand_Mass"]["binWidth"] = massWidth;
    userInput.Var["Cand_DLen"]["binWidth"] = 0.10;
    userInput.Var["Cand_DLenErr"]["binWidth"] = 0.005;
    userInput.Var["Cand_DLenRes"]["binWidth"] = 0.25;
    userInput.Var["Cand_DLenGen"]["binWidth"] = 0.025;
  }
  //
  // Clear extra variables if missing
  if (!contain(userInput.Par, "anaType")) { userInput.Par["anaType"]  = ""; }
  if (!contain(userInput.Int, "triggerIndex")) { userInput.Int["triggerIndex"] = -1; }
  if (!contain(userInput.Par, "treeType")) { userInput.Par["treeType"] = ""; }
  const bool& useExtDS = (userInput.Flag.at("useExtDS") && (workDirName.find("Test")==std::string::npos));
  for (const auto& var : userInput.StrV.at("variable")) { if (userInput.Flag.at("fit"+var) || !useExtDS) { userInput.Par["extFitDir_"+var] = ""; userInput.Par["extInitFileDir_"+var] = "";} }
  for (const auto& var : userInput.StrV.at("variable")) { for (const auto& obj : userInput.StrV.at("object")) { if (userInput.Flag.at("fit"+var) || !useExtDS) { userInput.Par["extInitFileDir_"+var+"_"+obj] = ""; } } }
  //
  // Check the User Input Settings
  if (!checkSettings(userInput)){ return; }
  
  // Set the Local Work Enviroment
  StringVectorMap_t DIR;
  if(!iniWorkEnv(DIR, workDirName)){ return; }
  /////////////////////
  StringMap_t inputFitDir;
  for (const auto& var : userInput.StrV.at("variable")) { inputFitDir[var] = userInput.Par["extFitDir_"+var]; }
  StringDiMap_t inputInitialFilesDir;
  for (const auto& var : userInput.StrV.at("variable")) {
    for (const auto& obj : userInput.StrV.at("object")) { inputInitialFilesDir[var][obj] = userInput.Par["extInitFileDir_"+var+"_"+obj]; }
  }

  // Initiliaze all the input Fit and Initial File Directories
  StringMapVector_t inputFitDirs;
  StringDiMapVector_t inputInitialFilesDirs;
  iniFileDir(inputFitDirs, inputInitialFilesDirs, inputFitDir, inputInitialFilesDir, DIR);


  // -------------------------------------------------------------------------------
  // STEP 1: LOAD THE INITIAL PARAMETERS
  /*
    Input : List of initial parameters with format PT <tab> RAP <tab> CEN <tab> iniPar ... 
    Output: two vectors with one entry per kinematic bin filled with the cuts and initial parameters
  */

  std::vector< GlobalInfoVectorMap_t > infoMapVectors;
  BoolDiMap_t VARMAP;
  for (const auto& var : userInput.StrS.at("incVariable")) {
    for (const auto& obj : userInput.StrV.at("object")) {
      VARMAP[var][obj] = userInput.Flag.at("fit"+obj);
    }
  }
  BoolMap_t COLMAP;
  for (const auto& col : userInput.StrV.at("system")) {
    COLMAP[col] = userInput.Flag.at("fit"+col);
  }
  for(uint j = 0; j < DIR.at("input").size(); j++) {
    if (DIR.at("input").size()>1 && j==0) continue; // First entry is always the main input directory
    GlobalInfoVectorMap_t infoMapVector;
    for (const auto& VAR : VARMAP) {
      auto PARMAP = VAR.second;
      for (const auto& PAR : PARMAP) {
        if (PAR.second) {
          for (const auto& COL : COLMAP) {
            if(COL.second) {
              std::string dir = DIR.at("input")[j];
              if (inputInitialFilesDirs[j].at(VAR.first).at(PAR.first)!="") { dir = inputInitialFilesDirs[j].at(VAR.first).at(PAR.first); }
              std::string inputFile = "", name = (dir + "InitialParam_" + VAR.first + "_" + PAR.first);
              StringVector_t tryChannel = { userInput.Par.at("channel") , "" };
              StringVector_t trySystem  = { COL.first };
              if (userInput.Flag.at("doPA8Y16")) { trySystem.push_back("PA8Y16"); }; trySystem.push_back("");
              for (const auto& tryCha : tryChannel) {
                bool trySuccess = false;
                for (const auto& tryCol : trySystem) {
                  if (ifstream(inputFile).good()==false) { inputFile = (name + tryCha + (tryCol!="" ? "_"+tryCol : "") + ".csv"); } else { trySuccess = true; break; }
                }
                if (trySuccess) break;
              }
              if (!addParameters(infoMapVector[COL.first], userInput, inputFile)) { return; }
            }
          }
        }
      }
    }
    infoMapVectors.push_back(infoMapVector);
  }

  // -------------------------------------------------------------------------------
  // STEP 2: CREATE/LOAD THE ROODATASETS
  /*
    Input : List of TTrees with format:  TAG <tab> FILE_NAME
    Output: Collection of RooDataSets splitted by tag name
  */
  RooWorkspaceMap_t iniWorkspaces;
  if (!createDataSets(iniWorkspaces, userInput, DIR)) { return; }

  // -------------------------------------------------------------------------------
  // STEP 3: PROCESS THE ROODATASETS
  /*
    Input : Collection of RooWorkspaces containing the original RooDataSets
    Output: Collection of RooWorkspaces containing the processed RooDatasets
  */

  // STEP 3.1: Skim the input RooDataSets
  if (!processDataSet(iniWorkspaces, userInput)) { return; }

  // -------------------------------------------------------------------------------  
  // STEP 4: FIT THE DATASETS
  /*
    Input : 
    -> The cuts and initial parameters per kinematic bin
    -> The workspace with the full datasets included.
    Output: 
    -> Plots (png, pdf and C format) of each fit.
    -> The local workspace used for each fit.
  */

  for(uint j = 0; j < infoMapVectors.size(); j++) {
    const auto& index = (DIR.at("output").size()>1 ? j+1 : j); // First entry is always the main output directory
    const auto& outputDir = DIR.at("output")[index];
    //
    for (const auto& DSTAG : userInput.StrS.at("DSTAG")) {
      const auto& dsCol = DSTAG.substr(DSTAG.rfind("_")+1);
      //
      if (contain(iniWorkspaces, DSTAG)) {
        //
        for (const auto& infoMapVector : infoMapVectors[j]) {
          const auto& col = infoMapVector.first;
          if (userInput.Flag.at("fit"+col) && (col==dsCol)) {
            //
            for (const auto& infoVector : infoMapVector.second) {
              //
	      std::cout << "[INFO] Proceed to fit the dataset " << DSTAG << std::endl;
              if (userInput.Par.at("analysis").rfind("CandTo", 0)==0) {
                if (!fitCandidateModel( iniWorkspaces, infoVector,
					userInput,
					// Select the type of datasets to fit
					outputDir,
					DSTAG,
					saveAll
					)
                    ) { return; }
              }
            }
          }
        }
      }
      else {
        std::cout << "[ERROR] The workspace for " << DSTAG << " was not found!" << std::endl; return;
      }
    }
  }
  std::cout << "[INFO] All fits done!" << std::endl;
};


bool createDataSets(RooWorkspaceMap_t& workspace, GlobalInfo& userInput, const StringVectorMap_t& DIR)
{
  //
  std::string inputTrees = DIR.at("input")[0] + "InputTrees.txt";
  if (existFile(inputTrees)==false && userInput.Par.at("extTreesFileDir")!="") { inputTrees = userInput.Par.at("extTreesFileDir") + "InputTrees.txt"; }
  StringDiVectorMap_t inputFileCollection;
  if(!getInputFileNames(inputFileCollection, inputTrees)){ return false; }
  //
  gSystem->Exec(Form("rm -f %s/cpp/AutoDict*", getcwd(NULL,0)));
  userInput.StrS["DSTAG"].clear(); // Array to store the different tags in the list of trees
  for (const auto& fileCollection : inputFileCollection) {
    // Get the file tag which has the following format: DSTAG_PD_CHAN_COLL , i.e. DATA_MUON_DIMUON_pPb5Y16
    const auto& FILETAG = fileCollection.first;
    if (!FILETAG.size()) { std::cout << "[ERROR] FILETAG is empty!" << std::endl; return false; }
    userInput.Par["localDSDir"] = DIR.at("dataset")[0];
    // Check file tag
    bool ignore = false;
    for (const auto& cha : userInput.StrV.at("channel")) {
      auto CHA = cha; std::transform(CHA.begin(), CHA.end(), CHA.begin(), ::toupper);
      if ((FILETAG.rfind(CHA)!=std::string::npos) && !userInput.Flag.at("do"+cha)) { ignore = true; break; }
    }
    if (ignore) continue;
    const auto& col = FILETAG.substr(FILETAG.rfind("_")+1);
    if (!contain(userInput.Flag, "do"+col) || !userInput.Flag.at("do"+col)) continue;
    // Extract the filenames
    std::string dir = "";
    bool fitDS = false;
    // If we have data, check if the user wants to fit data
    if (FILETAG.rfind("DATA_"+userInput.Par.at("PD")+"_", 0)==0) {
      if (userInput.Flag.at("fitData")) {
        dir = userInput.Par.at("localDSDir");
        if (userInput.Flag.at("useExtDS") && userInput.Par.at("extDSDir_DATA")!="" && existDir(userInput.Par.at("extDSDir_DATA"))) { dir = userInput.Par.at("extDSDir_DATA"); }
        fitDS = true;
      }
    }
    // If we find MC, check if the user wants to fit MC
    if (FILETAG.rfind("MC", 0)==0) {
      bool keep = false;
      if (userInput.Flag.at("fitMC")) {
        for (const auto& obj : userInput.StrS.at("fitObject")) { if (userInput.Flag.at("fit"+obj) && FILETAG.find(obj)!=std::string::npos) { keep = true; fitDS = true; break; } }
      }
      for(const auto& v : userInput.StrS.at("fitVarName")) {
	for(const auto& s : userInput.StrS.at("fitObject")) {
	  if (userInput.Flag.at("incMCTemp_"+v+"_"+s)) {
	    for (const auto& tmp : userInput.StrS.at("template_"+v)) {
	      auto n = tmp; if (n.rfind("Swap")!=std::string::npos) { n.erase(n.rfind("Swap"), 4); }
	      if (userInput.Flag.at("incMCTemp_"+v+"_"+s+"_"+tmp) && FILETAG.find(n)!=std::string::npos) { keep = true; fitDS = false; break; }
	    }
	  }
	}
      }
      if (keep) {
        dir = userInput.Par.at("localDSDir");
        if (userInput.Flag.at("useExtDS") && userInput.Par.at("extDSDir_MC")!="" && (existDir(userInput.Par.at("extDSDir_MC"))==true)) { dir = userInput.Par.at("extDSDir_MC"); }
      }
    }
    if (dir!="") {
      StringVectorMap_t fileInfo;
      fileInfo["inputFileNames"].clear(); fileInfo["treeTags"].clear();
      for (const auto& row : fileCollection.second) {
        fileInfo.at("inputFileNames").push_back(row[0]);
        if (row.size()>1) { fileInfo.at("treeTags").push_back(row[1]); }
      }
      fileInfo["outputFileDir"].push_back(dir);
      fileInfo["outputFileDir"].push_back(DIR.at("dataset")[0]);
      if (FILETAG.rfind("_PA8Y16")==std::string::npos) { fileInfo["dsNames"].push_back(FILETAG); }
      else {
        std::string NAMETAG = FILETAG; NAMETAG.erase(NAMETAG.rfind("_"), 10);
        if (userInput.Flag.at("dopPb8Y16")) fileInfo["dsNames"].push_back(NAMETAG+"_pPb8Y16");
        if (userInput.Flag.at("doPbp8Y16")) fileInfo["dsNames"].push_back(NAMETAG+"_Pbp8Y16");
      }
      // Produce the output datasets
      if(!tree2DataSet(workspace, fileInfo, userInput)){ return false; }
      if (fitDS) { for (const auto& DSTAG : fileInfo.at("dsNames")) { userInput.StrS.at("DSTAG").insert(DSTAG); } }
      else { userInput.Flag.at("fit"+col) = false; }
    }
  }
  if (workspace.empty()) {
    std::cout << "[ERROR] No tree files were found matching the user's input settings!" << std::endl; return false;
  }
  //
  return true;
};


bool addParameters(GlobalInfoVector_t& infoVector , GlobalInfo& userInfo , const std::string& inputFile)
{
  StringMapVector_t  data;
  if(!parseFile(data, inputFile)) { return false; }
  if (infoVector.empty()) {
    for (const auto& row : data) {
      GlobalInfo info = GlobalInfo();
      if(!setParameters(info, userInfo, row)) { return false; }
      infoVector.push_back(info);
    }
  }
  else {
    if (data.size()!=infoVector.size()) {
      std::cout << "[ERROR] The initial parameters in file " << inputFile << " ( " << data.size() << " ) are not consistent with previous files ( " << infoVector.size() << " ) !" << std::endl; return false;
    }
    for (unsigned int i=0; i<data.size(); i++) {
      GlobalInfo info = GlobalInfo();
      if (!setParameters(info, userInfo, data.at(i))) { return false; };
      if (info.Var != infoVector.at(i).Var) { std::cout << "[ERROR] The bins in file " << inputFile << " are not consistent with previous files!" << std::endl; return false; }
      infoVector.at(i).Copy(info, true);
    }
  }
  return true;
};


bool setParameters(GlobalInfo& info, GlobalInfo& userInfo, const StringMap_t& row)
{
  //
  const auto& analysis = userInfo.Par.at("analysis");
  // set initial values of variables
  if (analysis.rfind("CandTo", 0)==0) {
    info.Var["Cand_Mass"]["Min"]    = 0.0;
    info.Var["Cand_Mass"]["Max"]    = 100000.0;
    info.Var["Cand_Pt"]["Min"]      = 0.0;
    info.Var["Cand_Pt"]["Max"]      = 100000.0;
    info.Var["Cand_DLen"]["Min"]    = -100000.0;
    info.Var["Cand_DLen"]["Max"]    = 100000.0;
    info.Var["Cand_DLenErr"]["Min"] = -1.0;
    info.Var["Cand_DLenErr"]["Max"] = 100000.0;
    info.Var["Cand_DLenRes"]["Min"] = -100000.0;
    info.Var["Cand_DLenRes"]["Max"] = 100000.0;
    info.Var["Cand_DLenGen"]["Min"] = -100000.0;
    info.Var["Cand_DLenGen"]["Max"] = 100000.0;
    info.Var["Cand_Rap"]["Min"]     = -2.5;
    info.Var["Cand_Rap"]["Max"]     = 2.5;
    info.Var["Cand_AbsRap"]["Min"]  = 0.0;
    info.Var["Cand_AbsRap"]["Max"]  = 2.5;
    info.Var["Cand_APhi"]["Min"]    = -3.5;
    info.Var["Cand_APhi"]["Max"]    = 3.5;
    info.Var["Dau1_Pt"]["Min"]      = -1.;
    info.Var["Dau1_Pt"]["Max"]      = 100000.0;
    info.Var["Dau2_Pt"]["Min"]      = -1.;
    info.Var["Dau2_Pt"]["Max"]      = 100000.0;
    info.Var["Dau1_Eta"]["Min"]     = -2.5;
    info.Var["Dau1_Eta"]["Max"]     = 2.5;
    info.Var["Dau2_Eta"]["Min"]     = -2.5;
    info.Var["Dau2_Eta"]["Max"]     = 2.5;
    info.Var["Centrality"]["Min"]   = -1.0;
    info.Var["Centrality"]["Max"]   = 100000.0;
    info.Var["NTrack"]["Min"]       = -1.0;
    info.Var["NTrack"]["Max"]       = 1000000.0;
    info.Par["Cand_isSwap"] = "";
    info.Par["Cand_Qual"] = "";
  }
  for (auto& v : info.Var) {
    v.second["Default_Min"] = v.second.at("Min");
    v.second["Default_Max"] = v.second.at("Max");
  }
  info.Par["Cut"]   = "";
  for (const auto& var : userInfo.StrS.at("incVarName")) { info.Par["Model"+var] = ""; }
  for (const auto& v : info.Var) { if (contain(row, v.first+"CM")) { info.Flag["use"+v.first+"CM"] = true; } }
  // set parameters from file
  for (const auto& col : row) {
    std::string colName = col.first; stringReplace(colName, "CM", "");
    if (contain(info.Flag, "use"+col.first+"CM") && info.Flag.at("use"+col.first+"CM")) continue;
    bool found = false;
    for (const auto& var : info.Var) {
      const auto& varName = var.first;
      if (colName==varName) {
        if (col.second=="") {
          std::cout << "[ERROR] Input column " << varName << " has invalid value: " << col.second << std::endl; return false;
        }
	else if (col.second=="NONE") { found = true; break; }
        std::vector<double> v;
        if (!parseString(v, col.second)) { return false; }
        if (v.size()!=2 && v.size()!=1) {
          std::cout << "[ERROR] Input column " << varName << " has incorrect number of values, it should have 1 or 2 values but has: " << v.size() << std::endl; return false;
        }
        info.Var.at(varName).at("Min") = v.at(0);
        info.Var.at(varName).at("Max") = v.at(v.size()-1);
        found = true; break;
      }
    }
    if (found==false) {
      if (colName.rfind("Model",0)==0) {
	if (col.second=="") {
	  std::cout << "[ERROR] Input column " << colName << " has empty value" << std::endl; return false;
	}
	auto var = colName.substr(0, colName.find("_")).substr(5);
	if (var=="" && userInfo.StrS.at("incVarName").size()>1) {
	  std::cout << "[ERROR] Model name " << colName << " is not valid for multi-variable fits!" << std::endl; return false;
	}
	else if (var!="" && !contain(userInfo.StrS.at("incVarName"), var)) continue;
	var = (var=="" ? *userInfo.StrS.at("incVarName").begin() : var);
	const auto& modelName = "Model"+var+colName.substr(colName.find("_"));
	info.Par[modelName] = col.second;
	for(const auto& s : userInfo.StrS.at("fitObject")) {
	  if (colName.find(s)!=std::string::npos) {
	    const auto& mcTemp = "incMCTemp_"+var+"_"+s;
	    const bool& incTemp = (userInfo.Flag.at("fitMC")==false && col.second.find("TEMP")!=std::string::npos);
	    if (!contain(userInfo.Flag, mcTemp)) { userInfo.Flag[mcTemp] = incTemp; }
	    StringVector_t modV;
	    splitString(modV, col.second, "+");
	    for (const auto& mod : modV) {
	      if (mod.find("[")==std::string::npos) continue;
	      if (mod.rfind("]")==std::string::npos) { std::cout << "[ERROR] Missing ']' in " << colName << std::endl; return false; }
	      const auto& modT = mod.substr(0, mod.find("["));
	      auto modR = mod; modR.erase(0, modR.find("[")+1); modR.erase(modR.rfind("]"), modR.length());
	      StringVector_t modP;
	      if (modR.find(";")!=std::string::npos) { splitString(modP, modR, ";"); }
	      else { splitString(modP, modR, ","); }
	      for (const auto& o : modP) {
		bool incObj = (modT!="TEMP");
		if (modT=="TEMP" && userInfo.Flag.at(mcTemp)) {
		  userInfo.StrS["template_"+var].insert(o);
		  userInfo.Flag["incMCTemp_"+var+"_"+s+"_"+o] = true;
		  incObj = true;
		}
		if (incObj) {
		  userInfo.StrS["incObject_"+var].insert(o);
		  userInfo.Flag["inc"+o] = true;
		}
	      }
	    }
	  }
	}
	found = true;
      }
    }
    if (found==false) {
      for (const auto& par : info.Par) {
        std::string parName = par.first;
        if (colName.find(parName)!=std::string::npos) {
          if (col.second=="") {
            std::cout << "[ERROR] Input column " << colName << " has empty value" << std::endl; return false;
          }
	  info.Par[colName] = col.second;
          found = true;
        }
      }
    }
    if (found==false) {
      if (col.second != "") {
        std::string value = col.second;
        // check that initial parameters format is correct: [ num, num, num ]
        if ((value.find("[")==std::string::npos)||(value.rfind("]")==std::string::npos)) {
          // Special cases like parameter constrains could be set here but for now, let's keep it simple
          std::cout << "[ERROR] Either ']' or '[' are missing in the initial parameter values!" << std::endl; return false;
        }
        else {
          value.erase(0, value.find("[")+std::string("[").length());
          value.erase(value.rfind("]"), value.length());
        }
        std::vector<double> v;
        if (!parseString(v, value)) { return false; }
        if (v.size()>3 || v.size()<1) {
          std::cout << "[ERROR] Initial parameter " << col.first << " has incorrect number of values, it has: " << v.size() << std::endl; return false;
        }
        // everything seems alright, then proceed to save the values
        if (v.size()==1) {
          // if only one value is given i.e. [ num ], consider it a constant value
          info.Par[col.first] = Form("%s[%.10f]", col.first.c_str(), v.at(0));
        }
        else if (v.size()==2) {
          if (col.second.find("RNG")!=std::string::npos) {
            info.Par[col.first] = Form("%s[%.10f, %.10f]", col.first.c_str(), v.at(0), v.at(1));
          } 
          else { // For Constrained Fits
            info.Par[col.first] = Form("%s[%.10f, %.10f, %.10f]", col.first.c_str(), v.at(0), (v.at(0)-10.*v.at(1)), (v.at(0)+10.*v.at(1)));
            info.Par["val"+col.first] = Form("%s[%.10f]", ("val"+col.first).c_str(), v.at(0));
            info.Par["sig"+col.first] = Form("%s[%.10f]", ("sig"+col.first).c_str(), v.at(1));
          }
        }
        else if (v.size()==3) {
          info.Par[col.first] = Form("%s[%.10f, %.10f, %.10f]", col.first.c_str(), v.at(0), v.at(1), v.at(2));
        }
      }
      else {
        info.Par[col.first] = "";
      }
    }
  }
  return true;
};


bool parseString(std::vector<double>& output, std::string input)
{
  // remove spaces from input string 
  input.erase(std::remove(input.begin(), input.end(), ' '), input.end());
  // proceed to parse input string
  char *end;
  while(input!="") {
    double d = strtod(input.c_str(), &end);
    if (end != input) {
      output.push_back(d);
    }
    else {
      std::cout << "[ERROR] The conversion from string " << input << " to double failed!" << std::endl; return false;
    }
    input = end;
    if (input!="") { input.erase(input.begin()); } // Delete the delimiter, should have length = 1
  }
  return true;
};


bool parseFile(StringMapVector_t& data, const std::string& fileName)
{
  StringDiVector_t content, tmp;
  if(!readFile(tmp, fileName, -1, 1)){ return false; }
  StringVector_t header = tmp.at(0);
  if (header.empty()) { std::cout << "[ERROR] The header is null!" << std::endl; return false; }
  if(!readFile(content, fileName, header.size())){ return false; }
  for (const auto& rHeader : header) {
    if (rHeader=="") { std::cout << "[ERROR] A column has no label!" << std::endl; return false; }
  }
  content.erase(content.begin()); // remove header
  for (const auto& row : content) {
    StringMap_t col;
    for (unsigned int i=0; i<header.size(); i++) {
      if (i<row.size()) { col[header.at(i)] = row.at(i); }
      else { col[header.at(i)] = ""; }
    }
    data.push_back(col);
  }
  return true;
};


bool getInputFileNames(StringDiVectorMap_t& inputFileCollection, const std::string& inputTrees)
{
  StringDiVector_t content;
  if(!readFile(content, inputTrees, 3)){ return false; }
  for(auto& row : content) {
    for (auto& col : row) {
      col.erase(std::remove(col.begin(), col.end(), ' '), col.end());  // remove spaces
      col.erase(std::remove(col.begin(), col.end(), '\t'), col.end()); // remove tabs
    }
    if (row.at(0)!="" && row.at(1)=="") { std::cout << "[ERROR] There is an empty file name in your InputTrees.txt, please fix it" << std::endl; return false; }
    if (row.at(0)=="" && row.at(1)!="") { std::cout << "[ERROR] There is an empty file tag in your InputTrees.txt, please fix it" << std::endl; return false; }
    if (row.at(0)!="" && row.at(1)!="") {
      // store the filenames mapped by the tag
      StringVector_t tmp;
      for (uint i = 1; i < row.size(); i++) {
        tmp.push_back(row.at(i));
      }
      inputFileCollection[row.at(0)].push_back(tmp);
    }
  }
  return true;
};


bool readFile(StringDiVector_t& content, const std::string& fileName, const int& nCol, int nRow)
{
  if (nCol==0 || nRow==0) { 
    std::cout << "[WARNING] Ignoring content of File: " << fileName << std::endl; return true; 
  }
  if (nRow!=1) { std::cout << "[INFO] Reading file: " << fileName << std::endl; }
  ifstream myfile(fileName.c_str());
  char delimiter = ' ';
  if (myfile.is_open()){ 
    std::string line, CHAR;
    while (getline(myfile, line)) {
      std::stringstream row(line), tmp(line); tmp >> CHAR;
      if ((!tmp) || (CHAR.find('#')!=std::string::npos) || (CHAR.find("//")!=std::string::npos)) continue;
      if (delimiter == ' ' && line.find(',')!=std::string::npos) { delimiter = ','; }
      if (delimiter == ' ' && line.find(';')!=std::string::npos) { delimiter = ';'; }
      if (delimiter == ' ') { std::cout << "[ERROR] File: " << fileName << " has unknown delimiter!" << std::endl; return false; }
      if (nRow==0) break; else { nRow = nRow-1; }
      StringVector_t cols; int i=0;
      while (true){
        std::string col; getline(row, col, delimiter);
	if ((nCol>=0) ? (i>=nCol) : (col=="")) { break; }
	cols.push_back(col);
	i++;
      }
      content.push_back(cols);
    }
  }
  else {
    std::cout << "[ERROR] File: " << fileName << " was not found!" << std::endl; return false;
  }
  return true;
};


bool iniWorkEnv(StringVectorMap_t& DIR, const std::string& workDirName)
{
  std::cout << "[INFO] Initializing the work enviroment" << std::endl;
  DIR["main"].push_back(gSystem->ExpandPathName(gSystem->pwd()));
  DIR["macros"].push_back(DIR.at("main")[0] + "/Macros/");
  if (existDir(DIR.at("macros")[0].c_str())==false){ 
    std::cout << "[ERROR] Input directory: " << DIR.at("macros")[0] << " doesn't exist!" << std::endl;
    return false; 
  }
  DIR["input"].push_back(DIR.at("main")[0] + "/Input/" + workDirName + "/");
  if (existDir(DIR.at("input")[0])==false){ 
    std::cout << "[ERROR] Input directory: " << DIR.at("input")[0] << " doesn't exist!" << std::endl;
    return false; 
  }
  else {
    findSubDir(DIR.at("input"), DIR.at("input")[0]);
  }
  DIR["output"].push_back(DIR.at("main")[0] + "/Output/" + workDirName + "/");
  makeDir(DIR.at("output")[0]);
  for(uint j = 1; j < DIR.at("input").size(); j++) {
    auto subdir = DIR.at("input")[j];
    stringReplace(subdir, "/Input/", "/Output/");
    makeDir(subdir);
    DIR.at("output").push_back(subdir);
  }
  DIR["dataset"].push_back(DIR.at("main")[0] + "/DataSet/");
  makeDir(DIR.at("dataset")[0]);
  return true;
};


void iniFileDir(StringMapVector_t& inputFitDirs, StringDiMapVector_t& inputInitialFilesDirs, const StringMap_t& inputFitDir, const StringDiMap_t& inputInitialFilesDir, const StringVectorMap_t& DIR)
{
  inputFitDirs.push_back(inputFitDir);
  for(uint i=1; i<DIR.at("input").size(); i++) {
    inputFitDirs.push_back(inputFitDir);
    for (const auto& iter : inputFitDirs[i]) {
      const auto& key = iter.first;
      if (inputFitDirs[i][key]!="") {
        inputFitDirs[i][key] = DIR.at("input")[i];
        inputFitDirs[i][key].replace(inputFitDirs[i][key].find(DIR.at("input")[0]), DIR.at("input")[0].length(), inputFitDirs[0][key]);
      }
    }
  }
  inputInitialFilesDirs.push_back(inputInitialFilesDir);
  for(uint i=1; i<DIR.at("input").size(); i++) {
    inputInitialFilesDirs.push_back(inputInitialFilesDir);
    for (const auto& iter : inputInitialFilesDirs[i]) {
      for (const auto& iter2 : iter.second) {
        const auto& var  = iter.first;
        const auto& type = iter2.first;
        if (inputInitialFilesDirs[i].at(var).at(type)!="") {
          inputInitialFilesDirs[i].at(var).at(type) = DIR.at("input")[i];
          inputInitialFilesDirs[i].at(var).at(type).replace(inputInitialFilesDirs[i].at(var).at(type).find(DIR.at("input")[0]), DIR.at("input")[0].length(), inputInitialFilesDirs[0].at(var).at(type));
        }
      }
    }
  }
};


bool checkSettings(const GlobalInfo& userInput)
{ 
  std::cout << "[INFO] Checking user settings " << std::endl;
  const auto& analysis = userInput.Par.at("analysis");
  //
  if (analysis.rfind("CandTo", 0)==0) { /*NEED TO IMPLEMENT CHECK*/ }
  //
  // Check data to fit
  bool foundSample = false;
  for (const auto& s : userInput.StrV.at("sample")) { if (userInput.Flag.at("fit"+s)) { foundSample = true; break; } }
  if (!foundSample) { std::cout << "[ERROR] Neither data or MC were selected!" << std::endl; return false; }
  // Check dataset
  if (userInput.Par.at("PD")=="") { std::cout << "[ERROR] No trigger dataset was selected!" << std::endl; return false; }
  // Check collision system
  bool foundSystem = false;
  for (const auto& s : userInput.StrV.at("system")) { if (userInput.Flag.at("fit"+s)) { foundSystem = true; break; } }
  if (!foundSystem) { std::cout << "[ERROR] No collision system was selected!" << std::endl; return false; }
  //
  std::cout << "[INFO] All user setting are correct " << std::endl;
  return true;
};
