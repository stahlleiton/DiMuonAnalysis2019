#ifndef Candidate_fitCandidateModel_C
#define Candidate_fitCandidateModel_C


#include "RooFitResult.h"

#include "../../../Utilities/RunInfo/eventUtils.h"
#include "../Utilities/initClasses.h"
#include "../Utilities/rooDataUtils.h"
#include "buildCandidateModel.C"
#include "drawCandidatePlot.C"


void defineFitParameterRange ( GlobalInfo& info );


bool fitCandidateModel( const RooWorkspaceMap_t& inputWorkspaces, // Workspace with all the input RooDatasets
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
  info.Par["DSTAG"] = DSTAG;
  info.Par["outputDir"] = outputDir;
  info.Flag["saveAll"] = saveAll;

  // Check the input settings
  // Check if is MC
  info.Flag["fitMC"] = (DSTAG.rfind("MC", 0)==0);
  // Figure out the collision system to fit
  info.StrS["fitSystem"].clear();
  for (const auto& col : info.StrV.at("system")) { 
    info.Flag["fit"+col] = (DSTAG.rfind(col)!=std::string::npos);
    if (info.Flag.at("fit"+col)) { info.StrS.at("fitSystem").insert(col); } 
  }

  // If we use Center-of-Mass frame, then go back to LAB frame because dataset variables are in LAB frame
  for (auto& v : info.Var) {
    if (contain(info.Flag, "use"+v.first+"CM") && info.Flag.at("use"+v.first+"CM")) {
      const bool& ispPb = (info.Flag.at("fitpPb8Y16") || info.Flag.at("fitPA8Y16"));
      std::cout << "[INFO] Using " << v.first << " at Centre of Mass from " << (ispPb ? "p-Pb" : "Pb-p") << " LAB system"  << std::endl;
      const auto& varCM = info.Var.at(v.first+"CM");
      v.second.at("Min") = pPb::EtaCMtoLAB(varCM.at("Min"), ispPb);
      v.second.at("Max") = pPb::EtaCMtoLAB(varCM.at("Max"), ispPb);
      v.second["Full_Min"] = v.second.at("Min");
      v.second["Full_Max"] = v.second.at("Max");
      v.second["Plot_Min"] = v.second.at("Min");
      v.second["Plot_Max"] = v.second.at("Max");
      std::cout << "LAB: " << v.second.at("Min") << "  " << v.second.at("Max") << std::endl;
    }
  }

  // Define the range of all the fit parameters
  defineFitParameterRange(info);
  
  // Set models based on input files
  if (!setModel(info)) { return false; }

  // Define all the datasets needed for the fit
  for (const auto& chg : info.StrS.at("fitCharge")) {
    info.Par["dsName"+chg] = ("d"+chg+"_"+DSTAG);
    info.Par["dsNameFit"+chg] = ("d"+chg+"_"+DSTAG);
  }
  // Add the main dataset to the list
  info.StrS["dsList"].insert(DSTAG);

  // Add the datasets needed for the template fits
  for (const auto& tag : info.StrS.at("tags")) {
    if (!contain(info.StrS, "TEMPDS_"+tag)) continue;
    const auto& obj = tag.substr(0, tag.find(info.Par.at("channel")));
    if (!info.Flag.at("incMCTemp_"+obj)) { std::cout << "[ERROR] The input file for " << tag << " include templates but the input flag " << ("incMCTemp_"+obj) << " is false" << std::endl; return false; }
    for (const auto& tempDS : info.StrS.at("TEMPDS_"+tag)) { info.StrS.at("dsList").insert(tempDS); }
  }

  // Proceed to import the list of datasets
  bool doFit = true;
  for (const auto& chg : info.StrS.at("fitCharge")) {
    const auto& dsName = info.Par.at("dsName"+chg);
    if ( !(myws[chg].data(dsName.c_str())) ) {
      const auto& importID = importDataset(myws.at(chg), info, inputWorkspaces, chg);
      if (importID<0) { return false; }
      else if (importID==0) { doFit = false; }
    }
    info.Var["numEntries"][dsName] = myws.at(chg).data(dsName.c_str())->sumEntries();
    if (info.Var.at("numEntries").at(dsName)<=0) { doFit = false; }
  }
  if (!doFit) { std::cout << "[ERROR] No entries to fit!" << std::endl; return false; }

  // Proceed to import the sPlot datasets
  if (userInput.Par.at("fitSampleType")=="SPLOT") {
    for (const auto& chg : info.StrS.at("fitCharge")) { if (!importSPlotDataset(myws.at(chg), info, chg, true)) { return false; } }
  }

  // Set fit parameter range in workspace
  for (const auto& chg : info.StrS.at("fitCharge")) { if (!setFitParameterRange(myws.at(chg), info)) { return false; } }

  // Update fit parameter range
  for (const auto& chg : info.StrS.at("fitCharge")) { if (!updateFitRange(myws.at(chg), info, chg)) { return false; } }

  // Store the number of MC entries passing analysis cuts
  for (const auto& col : info.StrS.at("fitSystem")) {
    for (const auto& var : info.StrS.at("fitCondVarName")) {
      for (const auto& obj : info.StrS.at("template_"+var)) {
	const auto& dsLabel = "MC_" + obj + "_" + info.Par.at("channelDS") + "_" + col;
	for (const auto& chg : info.StrS.at("fitCharge")) {
	  const auto& dsName = "d"+chg+"_"+dsLabel;
	  if (myws.at(chg).data(dsName.c_str())) {
	    const auto& label = obj + info.Par.at("channel") + chg + "_" + col;
	    info.Var["recoMCEntries"][label] = myws.at(chg).data(dsName.c_str())->sumEntries();
	  }
	}
      }
    }
  }

  // Build the fit model
  for (const auto& chg : info.StrS.at("fitCharge")) { if (!buildCandidateModel(myws.at(chg), info, chg))  { return false; } }

  // Proceed to Fit and Save the results
  const auto& cha = info.Par.at("channel");
  for (const auto& col : info.StrS.at("fitSystem")) {
    for (const auto& chg : info.StrS.at("fitCharge")) {
      // Save the info in the workspace
      addString(myws.at(chg), "DSTAG", DSTAG);
      addString(myws.at(chg), "channel", cha);
      addString(myws.at(chg), "fitSystem", col);
      addString(myws.at(chg), "fitCharge", chg);
      addString(myws.at(chg), "PD", info.Par.at("PD"));
      defineSet(myws.at(chg), "fitVariable", info.StrS.at("fitVariable"));
      defineSet(myws.at(chg), "condVariable", info.StrS.at("condVariable"));
      //
      // Define output file name
      std::string fileName = "";
      std::string outDir = info.Par.at("outputDir");
      const auto& label = cha + chg + "_" + col;
      setFileName(fileName, outDir, label, info);
      std::string fitVar = ""; for (const auto& var : info.StrS.at("fitVarName")) { fitVar += var+"_"; };
      //
      // Get dataset and PDF names
      auto dsName = info.Par.at("dsName"+chg);
      auto dsNameFit = info.Par.at("dsNameFit"+chg);
      const auto& pdfName = info.Par.at("pdfName"+chg);
      addString(myws.at(chg), "dsName", dsName);
      addString(myws.at(chg), "dsNameFit", dsNameFit);
      addString(myws.at(chg), "pdfName", pdfName);
      const bool& isData = (dsName.find("_DATA_")!=std::string::npos);
      //
      // Skip fit for Cand_DLenErr
      const bool& skipFit = contain(info.StrS.at("fitVariable"), "Cand_DLenErr");
      //
      if (!skipFit) {
	// Check if we have already done this fit. If yes, continue.
	bool found =  true;
	std::unique_ptr<RooArgSet> newpars;
	if (myws.at(chg).pdf(pdfName.c_str())) { newpars = std::unique_ptr<RooArgSet>(myws.at(chg).pdf(pdfName.c_str())->getParameters(*myws.at(chg).data(dsName.c_str()))); }
	else { newpars = std::unique_ptr<RooArgSet>(new RooArgSet(myws.at(chg).allVars())); }
	const auto& outFileName = (outDir+"result/FIT_"+fitVar+fileName+".root");
	found = found && isFitAlreadyFound(*newpars, outFileName);
	if (found) {
	  std::cout << "[INFO] This fit for " << pdfName << " was already done, so I'll just go to the next one." << std::endl;
	  continue;
	}
      }
      //
      // Check if the user wants to do binned fits and proceed to bin the data
      const auto& fitVars = info.StrS.at("fitVariable");
      if (info.Flag["doBinnedFit"] && !createBinnedDataset(myws.at(chg), dsNameFit, *fitVars.begin())) { return false; }
      //
      // Store the mean and RMS of each dataset variable
      if (!storeDSStat(myws.at(chg), dsNameFit)) { return false; }
      // Store number of entries
      for (const auto& v : info.Var.at("numEntries")) { myws.at(chg).factory(Form("numEntries_%s[%.0f]", v.first.c_str(), v.second)); }
      //
      if (!skipFit) {
	// Fit the datasets
	bool fitFailed = false;
	std::unique_ptr<RooFitResult> fitResult;
        if (myws.at(chg).pdf(pdfName.c_str())) {
          const bool& isWeighted = myws.at(chg).data(dsNameFit.c_str())->isWeighted();
          const auto& numCores = info.Int.at("numCores");
          const auto& pdfConstrains = dynamic_cast<RooArgList*>(myws.at(chg).genobj(("pdfConstr"+label).c_str()));
	  uint opt = 2; if (info.Flag.at("fitBkg")) { opt = 0; }
          std::vector<RooCmdArg> cmdList = { RooFit::Extended(true), RooFit::AsymptoticError(false), RooFit::InitialHesse(true), RooFit::Minos(false), RooFit::Strategy(opt), RooFit::Minimizer("Minuit2"),
					     RooFit::Optimize(false), RooFit::NumCPU(numCores, 1), RooFit::Save(true), RooFit::Timer(true), RooFit::PrintLevel(1), RooFit::BatchMode(true)};
          if (pdfConstrains!=NULL && pdfConstrains->getSize()>0) {
            std::cout << "[INFO] Using constrain PDFs to fit " << pdfName << " on " << dsNameFit << std::endl;
            cmdList.push_back(RooFit::ExternalConstraints(*pdfConstrains));
            fitFailed = !fitPDF(fitResult, myws.at(chg), cmdList, pdfName, dsNameFit);
          }
          else {
            std::cout << "[INFO] Fitting " << pdfName << " on " << dsNameFit << std::endl;
            fitFailed = !fitPDF(fitResult, myws.at(chg), cmdList, pdfName, dsNameFit);
	    if (fitFailed) {
	      for (uint iTry=1; iTry<=opt; iTry++) {
		cmdList[4] = RooFit::Strategy(opt-iTry);
		fitFailed = !fitPDF(fitResult, myws.at(chg), cmdList, pdfName, dsNameFit, "initialParameters");
		if (!fitFailed) break;
	      }
	    }
	    if (false && isData && !fitFailed && contain(fitVars, "Cand_Mass")) {
              cmdList[3] = RooFit::Minos(true);
              fitFailed = !fitPDF(fitResult, myws.at(chg), cmdList, pdfName, dsNameFit);
	    }
          }
          if (fitResult) {
	    fitResult->SetTitle(fitResult->GetName());
	    myws.at(chg).import(*fitResult, fitResult->GetName());
	  }
          else { std::cout << "[ERROR] Fit Result returned by the PDF is NULL!" << std::endl; return false; }
	  if (fitFailed) {
	    for (uint iSt = 0; iSt < fitResult->numStatusHistory(); iSt++) {
	      if (fitResult->statusCodeHistory(iSt)!=0) {
		std::cout << "[ERROR] Fit failed in " << fitResult->statusLabelHistory(iSt) << " with status " << fitResult->statusCodeHistory(iSt) << " !" << std::endl; break;
	      }
	    }
	    myws.at(chg).factory("FAILED[1.0]");
	  }
        }
        else if ( myws.at(chg).obj(pdfName.c_str()) ) {
	  const auto& lbl = *info.StrS.at("fitObject").begin()+label;
          std::cout << "[INFO] Using the CutAndCount method with the following cut: " << info.Par.at("Cut_"+lbl) << std::endl;
          auto ds = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(myws.at(chg).data(dsName.c_str())->reduce(RooFit::Cut(info.Par.at("Cut_"+lbl).c_str()), RooFit::Name(("CutAndCount_"+dsName).c_str()))));
          myws.at(chg).var(("N_"+lbl).c_str())->setVal( ds->sumEntries() );
          myws.at(chg).var(("N_"+lbl).c_str())->setError( std::sqrt(ds->sumEntries()) );
          myws.at(chg).import(*ds);
          cout << Form("[INFO] Number of events that passed the cut: %.1f", myws.at(chg).var( ("N_"+lbl).c_str() )->getValV() ) << std::endl;
        }
        else {
          std::cout << "[ERROR] The PDF " << pdfName << " was not found!" << std::endl; return false;
        }
      }
      //
      // Draw the plot
      if (!drawCandidatePlot(myws.at(chg), fileName, outDir, info.Flag.at("setLogScale"), -1., true, false)) { return false; }
      //
      // Save the fit results
      saveSnapshot(myws.at(chg), "fittedParameters", info.Par.at("dsName"+chg));
      if (!saveWorkSpace(myws.at(chg), Form("%sresult/", outDir.c_str()), Form("%s.root", ("FIT_"+fitVar+fileName).c_str()), saveAll)) { return false; }
    }
  }
  //
  std::cout << "[INFO] Fit done, go to next bin!" << std::endl;  
  //
  return true;
  //
};


void defineFitParameterRange(GlobalInfo& info)
{
  for (const auto& var : info.StrV.at("variable")) {
    if (var!="Cand_Mass" && !info.Flag.at("fit"+var) && !info.Flag.at("cond"+var)) { continue; }
    double varMin = 100000000., varMax = -100000000.;
    double plotVarMin = 100000000., plotVarMax = -100000000.;
    if (var=="Cand_Mass") {
      const auto& vars = (contain(info.StrS, "incObject_CandMass") ? StringSet_t({"CandMass"}) : info.StrS.at("fitVarName"));
      StringSet_t objS;
      for (const auto& var : vars) {
	for (const auto& obj : info.StrS.at("incObject_"+var)) { objS.insert(obj); }
      }
      const auto& massRange = ANA::getMassRange(objS, info.Flag.at("fitMC"));
      varMin = massRange.first; varMax = massRange.second;
    }
    else if (var=="Cand_DLenErr") {
      if (info.Flag.at("fit"+var)==false) {
	const auto& cha = info.Par.at("channel");
	const auto& col = *info.StrS.at("fitSystem").begin();
	const auto& label = "JPsi" + cha + "OS" + "_" + col;
	if (!loadVarRange(info, var, label)) { assert("[ERROR] Could not load Cand_DLenErr range!"); }
      }
      else {
	varMin = 0.0;
	varMax = 10.0;
      }
      plotVarMin = 0.0;
      plotVarMax = 0.6;
    }
    else if (var=="Cand_DLenRes") {
      plotVarMin = -12.0;
      plotVarMax = 12.0;
    }
    else if (var=="Cand_DLenGen") {
      plotVarMin = -1.0;
      plotVarMax =  8.0;
    }
    else if (var=="Cand_DLen") {
      varMin = -30.0;
      varMax = 100.0;
      plotVarMin = -6.0;
      plotVarMax =  8.0;
    }
    if (varMax > varMin) {
      if (info.Var.at(var).at("Min") <= info.Var.at(var).at("Default_Min")) { info.Var.at(var).at("Min") = varMin; }
      if (info.Var.at(var).at("Max") >= info.Var.at(var).at("Default_Max")) { info.Var.at(var).at("Max") = varMax; }
    }
    info.Var.at(var).at("Full_Min") = info.Var.at(var).at("Min");
    info.Var.at(var).at("Full_Max") = info.Var.at(var).at("Max");
    info.Var.at(var).at("Plot_Min") = ((plotVarMax > plotVarMin) ? plotVarMin : info.Var.at(var).at("Min"));
    info.Var.at(var).at("Plot_Max") = ((plotVarMax > plotVarMin) ? plotVarMax : info.Var.at(var).at("Max"));
    // Print Information
    std::cout << "[INFO] Setting " << var << " fit range to min: " << info.Var.at(var).at("Min") << " and max: " << info.Var.at(var).at("Max");
    std::cout << " and plot range to min: " << info.Var.at(var).at("Plot_Min") << " and max: " << info.Var.at(var).at("Plot_Max") << endl;
  }
};


#endif // #ifndef Candidate_fitCandidateModel_C
