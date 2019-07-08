#ifndef Candidate_fitCandidateModel_C
#define Candidate_fitCandidateModel_C


#include "RooFitResult.h"

#include "../../../Utilities/RunInfo/eventUtils.h"
#include "../Utilities/initClasses.h"
#include "../Utilities/rooDataUtils.h"
#include "buildCandidateModel.C"
#include "drawCandidateMassPlot.C"


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
      std::cout << "CM: " << v.second.at("Min") << "  " << v.second.at("Max") << std::endl;
      v.second.at("Min") = pPb::EtaCMtoLAB(v.second.at("Min"), ispPb);
      v.second.at("Max") = pPb::EtaCMtoLAB(v.second.at("Max"), ispPb);
      std::cout << "LAB: " << v.second.at("Min") << "  " << v.second.at("Max") << std::endl;
    }
  }

  // Define the range of all the fit parameters
  defineFitParameterRange(info);

  // Set models based on input files
  if (!setModel(info)) { return false; }

  // Define all the datasets needed for the fit
  bool doFit = false;
  for (const auto& chg : info.StrS.at("fitCharge")) { info.Par["dsName"+chg] = ("d"+chg+"_"+DSTAG); }
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
  for (const auto& chg : info.StrS.at("fitCharge")) {
    const auto& dsName = info.Par.at("dsName"+chg);
    if ( !(myws[chg].data(dsName.c_str())) ) {
      const auto& importID = importDataset(myws.at(chg), info, inputWorkspaces, chg);
      if (importID<0) { return false; }
      else if (importID==0) { doFit = false; }
    }
    info.Var["numEntries"][chg] = myws.at(chg).data(dsName.c_str())->sumEntries();
    myws.at(chg).factory(Form("numEntries_%s[%.0f]", dsName.c_str(), info.Var.at("numEntries").at(chg)));
    if (info.Var.at("numEntries").at(chg)<=0) { doFit = false; }
  }

  // Store the number of MC ontries passing analysis cuts
  for (const auto& col : info.StrS.at("fitSystem")) {
    for (const auto& var : info.StrS.at("fitVarName")) {
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

  // Set fit parameter range in workspace
  for (const auto& chg : info.StrS.at("fitCharge")) { if (!setFitParameterRange(myws.at(chg), info)) { return false; } }

  // Build the fit model
  for (const auto& chg : info.StrS.at("fitCharge")) { if (!buildCandidateModel(myws.at(chg), info, chg))  { return false; } }

  // Proceed to Fit and Save the results
  std::string cha = info.Par.at("channel");
  for (const auto& col : info.StrS.at("fitSystem")) {
    for (const auto& chg : info.StrS.at("fitCharge")) {
      // Save the info in the workspace
      addString(myws.at(chg), "DSTAG", DSTAG);
      addString(myws.at(chg), "channel", cha);
      addString(myws.at(chg), "fitSystem", col);
      addString(myws.at(chg), "fitCharge", chg);
      addString(myws.at(chg), "PD", info.Par.at("PD"));
      defineSet(myws.at(chg), "fitVariable", info.StrS.at("fitVariable"));
      //
      // Define output file name
      std::string fileName = "";
      std::string outDir = outputDir;
      const auto& label = cha + chg + "_" + col;
      setFileName(fileName, outDir, label, info);
      //
      // Get dataset and PDF names
      const auto& dsName  = info.Par.at("dsName"+chg);
      const auto& pdfName = info.Par.at("pdfName"+chg);
      addString(myws.at(chg), "dsName", dsName);
      addString(myws.at(chg), "pdfName", pdfName);
      //
      // Check if the user wants to do binned fits and proceed to bin the data
      if (contain(info.Flag, "doBinnedFit") && info.Flag.at("doBinnedFit")) { if (!createBinnedDataset(myws.at(chg))) { return false; } }
      // Set the name of the dataset to fit
      const auto& dsNameFit = ((myws.at(chg).data((dsName+"_FIT").c_str())!=NULL) ? (dsName+"_FIT") : dsName);
      // Store the mean and RMS of each dataset variable
      if (!storeDSStat(myws.at(chg), dsNameFit)) { return false; }
      //
      // Check if we have already done this fit. If yes, do nothing and return true.
      bool found =  true; bool skipFit = false;
      std::unique_ptr<RooArgSet> newpars;
      if (myws.at(chg).pdf(pdfName.c_str())) { newpars = std::unique_ptr<RooArgSet>(myws.at(chg).pdf(pdfName.c_str())->getParameters(*myws.at(chg).data(dsName.c_str()))); }
      else { newpars = std::unique_ptr<RooArgSet>(new RooArgSet(myws.at(chg).allVars())); }
      const auto& outFileName = (outDir+"result/FIT_"+fileName+".root");
      found = found && isFitAlreadyFound(*newpars, outFileName);
      if (found) {
        std::cout << "[INFO] This fit for " << pdfName << " was already done, so I'll just go to the next one." << std::endl;
	continue;
      }
      //
      // Fit the datasets
      if (skipFit==false) {
	bool fitFailed = false;
	std::unique_ptr<RooFitResult> fitResult;
        if (myws.at(chg).pdf(pdfName.c_str())) {
          bool isWeighted = myws.at(chg).data(dsNameFit.c_str())->isWeighted();
          if (dsName.find("_DATA_")!=std::string::npos) { isWeighted = false; } // BUG FIX
          const auto& numCores = info.Int.at("numCores");
          const auto& pdfConstrains = dynamic_cast<RooArgList*>(myws.at(chg).genobj(("pdfConstr"+label).c_str()));
          if (pdfConstrains!=NULL && pdfConstrains->getSize()>0) {
            std::cout << "[INFO] Using constrain PDFs to fit " << pdfName << " on " << dsNameFit << std::endl;
            const auto& tmp = myws.at(chg).pdf(pdfName.c_str())->fitTo(*myws.at(chg).data(dsNameFit.c_str()), RooFit::Extended(kTRUE), RooFit::SumW2Error(isWeighted), RooFit::Strategy(2),
								       RooFit::Range("FitWindow"), RooFit::ExternalConstraints(*pdfConstrains), RooFit::NumCPU(numCores), RooFit::Save());
            fitResult.reset(tmp);
            fitFailed = false; for (uint iSt = 0; iSt < fitResult->numStatusHistory(); iSt++) { if (fitResult->statusCodeHistory(iSt)!=0) { fitFailed = true; break; } }
          }
          else {
            std::cout << "[INFO] Fitting " << pdfName << " on " << dsNameFit << std::endl;
            const auto& tmp = myws.at(chg).pdf(pdfName.c_str())->fitTo(*myws.at(chg).data(dsNameFit.c_str()), RooFit::Extended(kTRUE), RooFit::SumW2Error(isWeighted), RooFit::Strategy(2),
								       RooFit::Range("FitWindow"), RooFit::NumCPU(numCores), RooFit::Save());
            fitResult.reset(tmp);
            fitFailed = false; for (uint iSt = 0; iSt < fitResult->numStatusHistory(); iSt++) { if (fitResult->statusCodeHistory(iSt)!=0) { fitFailed = true; break; } }
	    if (fitFailed) {
	      std::cout << std::endl; std::cout << "[WARNING] Fit failed, trying again with Minuit2 and Minimizer" << std::endl; std::cout << std::endl;
	      const auto& tmp = myws.at(chg).pdf(pdfName.c_str())->fitTo(*myws.at(chg).data(dsNameFit.c_str()), RooFit::Extended(kTRUE), RooFit::SumW2Error(isWeighted),
									 RooFit::Minimizer("Minuit2","minimize"), RooFit::Strategy(2),
									 RooFit::Range("FitWindow"), RooFit::NumCPU(numCores), RooFit::Save());
	      fitResult.reset(tmp);
	      fitFailed = false; for (uint iSt = 0; iSt < fitResult->numStatusHistory(); iSt++) { if (fitResult->statusCodeHistory(iSt)!=0) { fitFailed = true; break; } }
	    }
	    if (fitFailed) {
	      std::cout << std::endl; std::cout << "[WARNING] Fit failed, trying again with Minuit2 and Scan" << std::endl; std::cout << std::endl;
	      const auto& tmp = myws.at(chg).pdf(pdfName.c_str())->fitTo(*myws.at(chg).data(dsNameFit.c_str()), RooFit::Extended(kTRUE), RooFit::SumW2Error(isWeighted),
									 RooFit::Minimizer("Minuit2","scan"), RooFit::Strategy(2),
									 RooFit::Range("FitWindow"), RooFit::NumCPU(numCores), RooFit::Save());
	      fitResult.reset(tmp);
	      fitFailed = false; for (uint iSt = 0; iSt < fitResult->numStatusHistory(); iSt++) { if (fitResult->statusCodeHistory(iSt)!=0) { fitFailed = true; break; } }
	    }
          }
          if (fitResult) {
	    fitResult->Print("v");
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
	//
	// Draw the plot
        if (!drawCandidateMassPlot(myws.at(chg), ("PLOT_"+fileName), outDir, info.Flag.at("setLogScale"), -1., true, false)) { return false; }
	//
	// Save the fit results
	if (!fitFailed) {
	  saveSnapshot(myws.at(chg), "fittedParameters", info.Par.at("dsName"+chg));
	  if (!saveWorkSpace(myws.at(chg), Form("%sresult/", outDir.c_str()), Form("%s.root", ("FIT_"+fileName).c_str()), saveAll)) { return false; }
	}
      }
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
  for (const auto& var : info.StrS.at("fitVariable")) {
    double varMin = 100000000., varMax = -100000000.;
    if (var=="Cand_Mass") {
      if (info.Flag.at("fitMC")) {
  	for (const auto& obj : info.StrS.at("incObject")) {
    	  if (contain(MASS, obj)) {
      	    if (varMin > MASS.at(obj).at("Min")) { varMin = MASS.at(obj).at("Min"); }
      	    if (varMax < MASS.at(obj).at("Max")) { varMax = MASS.at(obj).at("Max"); }
    	  }
  	}
      }
      else {
        // Define the Candidate Mass range
        if (contain(info.Flag, "incD0") && info.Flag.at("incD0")) {
	  if (varMin > 1.75) { varMin = 1.75; }
	  if (varMax < 2.00) { varMax = 2.00; }
        }
        if (contain(info.Flag, "incJPsi") && info.Flag.at("incJPsi")) {
	  if (varMin > 2.7) { varMin = 2.7; }
	  if (varMax < 3.5) { varMax = 3.5; }
        }
        if (contain(info.Flag, "incPsi2S") && info.Flag.at("incPsi2S")) {
  	  if (varMin > 3.4) { varMin = 3.4; }
  	  if (varMax < 4.0) { varMax = 4.0; }
        }
        if (contain(info.Flag, "incUps1S") && info.Flag.at("incUps1S")) {
	  if (varMin >  8.0) { varMin =  8.0; }
	  if (varMax < 10.5) { varMax = 10.5; }
        }
        if (contain(info.Flag, "incUps2S") && info.Flag.at("incUps2S")) {
	  if (varMin >  9.0) { varMin =  9.0; }
	  if (varMax < 11.0) { varMax = 11.0; }
        }
        if (contain(info.Flag, "incUps3S") && info.Flag.at("incUps3S")) {
	  if (varMin >  9.5) { varMin =  9.5; }
	  if (varMax < 13.5) { varMax = 14.0; }
        }
        if (contain(info.Flag, "incZ") && info.Flag.at("incZ")) {
	  if (varMin >  70.0) { varMin =  70.0; }
	  if (varMax < 110.0) { varMax = 110.0; }
        }
      }
    }
    if (varMax > varMin) {
      if (info.Var.at(var).at("Min") <= info.Var.at(var).at("Default_Min")) { info.Var.at(var).at("Min") = varMin; }
      if (info.Var.at(var).at("Max") >= info.Var.at(var).at("Default_Max")) { info.Var.at(var).at("Max") = varMax; }
    }
    // Print Information
    std::cout << "[INFO] Setting " << var << " range to min: " << info.Var.at(var).at("Min") << " and max: " << info.Var.at(var).at("Max") << endl;
  }
};


#endif // #ifndef Candidate_fitCandidateModel_C
