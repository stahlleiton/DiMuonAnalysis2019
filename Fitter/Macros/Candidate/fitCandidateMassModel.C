#ifndef fitCandidateMassModel_C
#define fitCandidateMassModel_C


#include "RooFitResult.h"

#include "../../../Utilities/RunInfo/eventUtils.h"
#include "../Utilities/initClasses.h"
#include "../Utilities/rooDataUtils.h"
#include "buildCandidateMassModel.C"
#include "drawCandidateMassPlot.C"


void setCandidateMassRange ( GlobalInfo& info );


bool fitCandidateMassModel( const RooWorkspaceMap_t& inputWorkspaces, // Workspace with all the input RooDatasets
                            const GlobalInfo&  inputInfo,             // Contains information on initial Parameters, cut values, flags, ...
                            const GlobalInfo&  userInput,             // Contains information on initial Parameters, cut values, flags, ...
                            const std::string& outputDir,             // Path to output directory
                            // Select the type of datasets to fit
                            const std::string& DSTAG,                 // Specifies the name of the dataset to fit
                            const bool& saveAll=true
                            )
{
  //
  // Set up the local workspace and the input information
  RooWorkspaceMap_t myws;
  GlobalInfo info(userInput);
  info.Copy(inputInfo, true); // Copy the user input information (avoid duplicating information in fitter)

  // Check the input settings
  // Check if is MC
  info.Flag["fitMC"] = (DSTAG.rfind("MC", 0)==0);
  // Figure out the collision system to fit
  info.StrV["fitSystem"].clear();
  for (const auto& col : info.StrV.at("system")) { 
    info.Flag["fit"+col] = (DSTAG.rfind(col)!=std::string::npos);
    if (info.Flag.at("fit"+col)) { info.StrV.at("fitSystem").push_back(col); } 
  }
  //

  // Set the range of all the parameters
  setCandidateMassRange(info);

  // Set models based on input files
  StringDiMap_t model;
  if (!setModel(model, info)) { return false; }
  
  // Import all the datasets needed for the fit
  bool doFit = false;
  for (const auto& chg : info.StrV.at("fitCharge")) { info.Par["dsName"+chg] = ("d"+chg+"_"+DSTAG); }
  // Add the main dataset to the list
  info.StrV["dsList"].push_back(DSTAG);

  // Check the datasets needed for templates
  bool addTemp = false;
  for (const auto& tag : info.StrV.at("tags")) {
    for (const auto& obj : info.StrV.at("fitObject")) {
      if (tag.find(obj)!=std::string::npos && (info.StrV.count("TEMPDS_"+tag)>0)) {
        if (!info.Flag.at("incMCTemp_"+obj)) { std::cout << "[ERROR] The input file for " << tag << " include templates but the input flag " << ("incMCTemp_"+obj) << " is false" << std::endl; return false; }
        addTemp = true;
        break;
      }
    }
  }
  // Add the datasets needed for the template fits
  if (addTemp) {
    for (const auto& tag : info.StrV.at("tags")) {
      if (tag.find("Pl_")==std::string::npos && tag.find("OS_")==std::string::npos) continue; // Only look at one charge since they are symmetric
      if (info.StrV.count("TEMPDS_"+tag)>0) {
        for (const auto& tempDS : info.StrV.at("TEMPDS_"+tag)) {
          if (std::find(info.StrV.at("dsList").begin(), info.StrV.at("dsList").end(), tempDS) == info.StrV.at("dsList").end()) { info.StrV.at("dsList").push_back(tempDS); }
        }
      }
    }
  }

  // Proceed to import the list of datasets
  for (const auto& chg : info.StrV.at("fitCharge")) {
    const std::string& dsName = ("d"+chg+"_"+DSTAG);
    info.Par[Form("dsName_%s", chg.c_str())] = dsName;
    if ( !(myws[chg].data(dsName.c_str())) ) {
      const auto& importID = importDataset(myws.at(chg), inputWorkspaces, info, chg);
      if (importID<0) { return false; }
      else if (importID==0) { doFit = false; }
    }
    info.Var["numEntries"][chg] = myws.at(chg).data(dsName.c_str())->sumEntries(); if (info.Var.at("numEntries").at(chg)<=0) { doFit = false; }
  }

  // Store the number of MC ontries passing analysis cuts
  for (const auto& mainObj : info.StrV.at("fitObject")) {
    if (info.Flag.at("incMCTemp_"+mainObj)) {
      for (const auto& col : info.StrV.at("fitSystem")) {
        for (const auto& obj : info.StrV.at("template")) {
          std::string dsLabel = "MC_" + obj + "_" + info.Par.at("channelDS") + "_" + col;
          for (const auto& chg : info.StrV.at("fitCharge")) {
            if (myws.at(chg).data(Form("d%s_%s", chg.c_str(), dsLabel.c_str()))) {
              std::string label = obj + info.Par.at("channel") + chg + "_" + col;
              info.Var["recoMCEntries"][label] = myws.at(chg).data(Form("d%s_%s", chg.c_str(), dsLabel.c_str()))->sumEntries();
            }
          }
        }
      }
    }
  }

  // Set global parameters
  for (const auto& chg : info.StrV.at("fitCharge")) { setGlobalParameterRange(myws.at(chg), info, "Cand_Mass"); }

  // Build the fit model
  for (const auto& chg : info.StrV.at("fitCharge")) { if (!buildCandidateMassModel(myws.at(chg), model, info, chg))  { return false; } }

  // Proceed to Fit and Save the results
  std::string cha = info.Par.at("channel");
  for (const auto& col : info.StrV.at("fitSystem")) {
    for (const auto& chg : info.StrV.at("fitCharge")) {
      // Save the info in the workspace
      RooStringVar tmp = RooStringVar();
      if (myws.at(chg).obj("DSTAG")) { ((RooStringVar*)myws.at(chg).obj("DSTAG"))->setVal(DSTAG.c_str()); }
      else { tmp.setVal(DSTAG.c_str()); tmp.SetTitle("DSTAG"); myws.at(chg).import(*((TObject*)&tmp), tmp.GetTitle()); }
      if (myws.at(chg).obj("channel")) { ((RooStringVar*)myws.at(chg).obj("channel"))->setVal(cha.c_str()); }
      else { tmp.setVal(cha.c_str()); tmp.SetTitle("channel"); myws.at(chg).import(*((TObject*)&tmp), tmp.GetTitle()); }
      if (myws.at(chg).obj("fitSystem")) { ((RooStringVar*)myws.at(chg).obj("fitSystem"))->setVal(col.c_str()); }
      else { tmp.setVal(col.c_str()); tmp.SetTitle("fitSystem"); myws.at(chg).import(*((TObject*)&tmp), tmp.GetTitle()); }
      if (myws.at(chg).obj("fitCharge")) { ((RooStringVar*)myws.at(chg).obj("fitCharge"))->setVal(chg.c_str()); }
      else { tmp.setVal(chg.c_str()); tmp.SetTitle("fitCharge"); myws.at(chg).import(*((TObject*)&tmp), tmp.GetTitle()); }
      // Total PDF Name
      const auto& label = cha + chg + "_" + col;
      const auto& pdfName = ( "pdfMass_Tot" + label );
      // Plot Name
      std::string plotLabel = "";
      for (const auto& obj : info.StrV.at("fitObject")) {
        const auto& objLabel = obj + cha + chg + "_" + col;
        std::string modelN = info.Par.at("Model_"+objLabel);
        modelN.erase(std::remove(modelN.begin(), modelN.end(), ' '), modelN.end());
        stringReplace( modelN, "[", "_" ); stringReplace( modelN, "]", "" ); stringReplace( modelN, "+", "_" ); stringReplace( modelN, ",", "" ); stringReplace( modelN, ";", "" );
        if (myws.at(chg).obj("modelName")) { ((RooStringVar*)myws.at(chg).obj("modelName"))->setVal(modelN.c_str()); }
        else { tmp.setVal(modelN.c_str()); tmp.SetTitle("modelName"); myws.at(chg).import(*((TObject*)&tmp), tmp.GetTitle()); }
        plotLabel += obj + "Model_" + modelN+"_";
      }
      // Output File Name
      std::string fileName = "";
      std::string outDir = outputDir;
      //setFileName(fileName, outDir, DSTAG, plotLabel, info);
      // Dataset Name
      const auto& dsName = ( "d" + chg + "_" + DSTAG );
      // Check if the user wants to do binned fits and proceed to bin the data
      if (info.Flag.count("doBinnedFit")>0 && info.Flag.at("doBinnedFit")) { if (!createBinnedDataset(myws.at(chg))) { return false; } }
      // Set the name of the dataset to fit
      const auto& dsNameFit = ( (myws.at(chg).data((dsName+"_FIT").c_str())!=NULL) ? (dsName+"_FIT") : dsName );
      // check if we have already done this fit. If yes, do nothing and return true.
      bool found =  true; bool skipFit = false;
      std::unique_ptr<RooArgSet> newpars;
      if (myws.at(chg).pdf(pdfName.c_str())) { newpars = std::unique_ptr<RooArgSet>(myws.at(chg).pdf(pdfName.c_str())->getParameters(*myws.at(chg).data(dsName.c_str()))); }
      else { newpars = std::unique_ptr<RooArgSet>(new RooArgSet(myws.at(chg).allVars())); }
      found = found && isFitAlreadyFound(*newpars, Form("%sresult/%s.root", outDir.c_str(), ("FIT_"+fileName).c_str()), pdfName.c_str());
      if (found) {
        std::cout << "[INFO] This fit for " << label << " was already done, so I'll just go to the next one." << std::endl;
        continue;
      }
      // Fit the Datasets
      if (skipFit==false) {
        if (myws.at(chg).pdf(pdfName.c_str())) {
          bool isWeighted = myws.at(chg).data(dsNameFit.c_str())->isWeighted();
          if (dsName.rfind("DATA", 0)==0) { isWeighted = false; } // BUG FIX
          const int& numCores = info.Int.at("numCores");
          RooArgList* pdfConstrains = (RooArgList*)myws.at(chg).genobj(Form("pdfConstr%s", label.c_str()));
          std::unique_ptr<RooFitResult> fitResult;
          bool fitFailed = false;
          if (pdfConstrains!=NULL && pdfConstrains->getSize()>0) {
            std::cout << "[INFO] Fitting with constrain PDFs" << std::endl;
            auto tmp = myws.at(chg).pdf(pdfName.c_str())->fitTo(*myws.at(chg).data(dsNameFit.c_str()), RooFit::Extended(kTRUE), RooFit::SumW2Error(isWeighted), RooFit::Strategy(2),
                                                                RooFit::Range("Cand_MassWindow"), RooFit::ExternalConstraints(*pdfConstrains), RooFit::NumCPU(numCores), RooFit::Save());
            fitResult.reset(tmp);
            bool fitFailed = false; for (uint iSt = 0; iSt < fitResult->numStatusHistory(); iSt++) { if (fitResult->statusCodeHistory(iSt)!=0) { fitFailed = true; break; } }
          }
          else {
            auto tmp = myws.at(chg).pdf(pdfName.c_str())->fitTo(*myws.at(chg).data(dsNameFit.c_str()), RooFit::Extended(kTRUE), RooFit::SumW2Error(isWeighted), RooFit::Strategy(2),
                                                                RooFit::Range("Cand_MassWindow"), RooFit::NumCPU(numCores), RooFit::Save());
            fitResult.reset(tmp);
            bool fitFailed = false; for (uint iSt = 0; iSt < fitResult->numStatusHistory(); iSt++) { if (fitResult->statusCodeHistory(iSt)!=0) { fitFailed = true; break; } }
          }
          if (fitResult!=NULL) {
            fitResult->Print("v");
            fitResult->SetTitle(Form("fitResult_%s", pdfName.c_str()));
            myws.at(chg).import(*fitResult, fitResult->GetTitle());
          }
          else { std::cout << "[ERROR] Fit Result returned by the PDF is NULL!" << std::endl; return false; }
          for (uint iSt = 0; iSt < fitResult->numStatusHistory(); iSt++) {
            if (fitResult->statusCodeHistory(iSt)!=0) {
              std::cout << "[ERROR] Fit failed in " << fitResult->statusLabelHistory(iSt) << " with status " << fitResult->statusCodeHistory(iSt) << " !" << std::endl; //return false;
            }
          }
        }
        else if ( myws.at(chg).obj(("CutAndCount_"+label).c_str()) ) {
          // cut and count
          std::cout << Form("[INFO] Using the CutAndCount Method with the following cut: %s", info.Par.at("Cut_"+label).c_str()) << std::endl;
          auto ds = std::unique_ptr<RooDataSet>((RooDataSet*)myws.at(chg).data(dsName.c_str())->reduce(RooFit::Cut(info.Par.at("Cut_"+label).c_str()), RooFit::Name(("CutAndCount_"+dsName).c_str())));
          myws.at(chg).var( ("N_"+label).c_str() )->setVal( ds->sumEntries() );
          myws.at(chg).var( ("N_"+label).c_str() )->setError( sqrt(ds->sumEntries()) );
          myws.at(chg).import(*ds);
          cout << Form("[INFO] Number of events that passed the cut: %.1f", myws.at(chg).var( ("N_"+label).c_str() )->getValV() ) << std::endl;
        }
        else {
          std::cout << "[ERROR] The PDF " << pdfName << " was not found!" << std::endl; return false;
        }
        if (!drawCandidateMassPlot(myws.at(chg), ("PLOT_"+fileName), outDir, info.Flag.at("setLogScale"), -1., true, false)) { return false; }
        myws.at(chg).saveSnapshot("fittedParameters", myws.at(chg).allVars(), kTRUE);
  myws.at("OS").Print("v");
  info.Print();
  return false;
        // Save the results
        //if (!saveWorkSpace(myws.at(chg), Form("%sresult/", outDir.c_str()), Form("%s.root", ("FIT_"+fileName).c_str()), saveAll)) { return false; }
      }
    }
  }


  myws.at("OS").Print("v"); return false;
  //
  return true;
  //
};


void setCandidateMassRange(GlobalInfo& info)
{
  // Define the Candidate Mass range
  double massMin = 100000000., massMax = -1.;
  if (info.Flag.count("fitD0") && info.Flag.at("fitD0")) {
    if (massMin > 1.75) { massMin = 1.75; }
    if (massMax < 2.00) { massMax = 2.00; }
  }
  if (info.Flag.count("fitJPsi") && info.Flag.at("fitJPsi")) {
    if (massMin > 2.6) { massMin = 2.6; }
    if (massMax < 3.5) { massMax = 3.5; }
  }
  if (info.Flag.count("fitPsi2S") && info.Flag.at("fitPsi2S")) {
    if (massMin > 3.4) { massMin = 3.4; }
    if (massMax < 4.0) { massMax = 4.0; }
  }
  if (info.Flag.count("fitUps1") && info.Flag.at("fitUps1S")) {
    if (massMin >  8.5) { massMin =  8.5; }
    if (massMax < 10.5) { massMax = 10.5; }
  }
  if (info.Flag.count("fitUps2S") && info.Flag.at("fitUps2S")) {
    if (massMin >  9.0) { massMin =  9.0; }
    if (massMax < 11.0) { massMax = 11.0; }
  }
  if (info.Flag.count("fitUps3S") && info.Flag.at("fitUps3S")) {
    if (massMin >  9.5) { massMin =  9.5; }
    if (massMax < 12.5) { massMax = 12.5; }
  }
  if (info.Flag.count("fitZ") && info.Flag.at("fitZ")) {
    if (massMin >  70.0) { massMin =  70.0; }
    if (massMax < 110.0) { massMax = 110.0; }
  }
  if (massMax > massMin) {
    if (info.Var.at("Cand_Mass").at("Min") <=     0.0) { info.Var.at("Cand_Mass").at("Min") = massMin; }
    if (info.Var.at("Cand_Mass").at("Max") >= 10000.0) { info.Var.at("Cand_Mass").at("Max") = massMax; }
  }
  // Print Information
  std::cout << "[INFO] Setting candidate mass range to min: " << info.Var.at("Cand_Mass").at("Min") << " and max: " << info.Var.at("Cand_Mass").at("Max") << endl;
  return;
};


#endif // #ifndef fitCandidateMassModel_C
