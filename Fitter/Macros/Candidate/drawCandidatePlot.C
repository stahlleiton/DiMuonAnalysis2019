#ifndef Candidate_drawCandidatePlot_C
#define Candidate_drawCandidatePlot_C


#include "../Utilities/drawUtils.h"


void     printCandidateParameters  ( TPad& pad , const RooWorkspace& ws , const std::string& pdfName , const uint& drawMode );
void     printCandidateTextInfo    ( TPad& pad , const RooWorkspace& ws , const std::string& varN , const uint& drawMode , const int& plotStyle );
TLegend* printCandidateLegend      ( TPad& pad , const RooPlot& frame   , const StringMap_t& legInfo , const int& plotStyle );


bool drawCandidatePlot( RooWorkspace& ws,  // Local Workspace
			const std::string& varName,
			const std::string& fileName,
			const std::string& outputDir,
			const bool& yLogScale,
			const double& maxRng = -1.0,
			const bool& doGoF = true,
			const int&  plotStyle = 0, // 2: Paper , 1: PAS , 0: AN
			const bool& redoFrame = false,
			const bool& doHistPDF = false
			)
{
  //
  const auto& setLogScale = yLogScale;
  auto varTag = varName; stringReplace(varTag, "_", "");
  //
  // set the CMS style
  setTDRStyle();
  //
  const std::string& DSTAG     = (ws.obj("DSTAG")     ? dynamic_cast<RooStringVar*>(ws.obj("DSTAG")    )->getVal() : "");
  const std::string& cha       = (ws.obj("channel")   ? dynamic_cast<RooStringVar*>(ws.obj("channel")  )->getVal() : "");
  const std::string& col       = (ws.obj("fitSystem") ? dynamic_cast<RooStringVar*>(ws.obj("fitSystem"))->getVal() : "");
  const std::string& chg       = (ws.obj("fitCharge") ? dynamic_cast<RooStringVar*>(ws.obj("fitCharge"))->getVal() : "");
  const std::string& PDFName   = (ws.obj("pdfName")   ? dynamic_cast<RooStringVar*>(ws.obj("pdfName")  )->getVal() : "");
  const std::string& dsName    = (ws.obj("dsName")    ? dynamic_cast<RooStringVar*>(ws.obj("dsName")   )->getVal() : "");
  const std::string& DSNameFit = (ws.obj("dsNameFit") ? dynamic_cast<RooStringVar*>(ws.obj("dsNameFit"))->getVal() : "");
  const std::string& PD        = (ws.obj("PD")        ? dynamic_cast<RooStringVar*>(ws.obj("PD")       )->getVal() : "");
  //
  const auto& dsNameFit = (ws.data(("CutAndCount_"+dsName).c_str()) ? ("CutAndCount_"+dsName) : DSNameFit);
  const auto& ds = ws.data(dsName.c_str());
  if (!ds) { std::cout << "[ERROR] Fit dataset " << dsName << " was not found!" << std::endl;  return false; }
  const auto& fitDS = ws.data(dsNameFit.c_str());
  if (!fitDS) { std::cout << "[ERROR] Fit dataset " << dsNameFit << " was not found!" << std::endl;  return false; }
  auto pdfName = "pdf"+varTag+PDFName.substr(PDFName.find("_")); if (!ws.pdf(pdfName.c_str())) { pdfName = PDFName; }
  if (doHistPDF) { stringReplace(pdfName, "pdf", "pdfHIST"); if (!ws.pdf(pdfName.c_str())) { pdfName = PDFName; } }
  const auto& fitPDF = dynamic_cast<RooAddPdf*>(ws.pdf(pdfName.c_str()));
  const auto& label = (cha + chg + "_" + col);
  //
  // Create the Range for Plotting
  const auto& fitVar = ws.var(varName.c_str());
  if (!fitVar) { std::cout << "[ERROR] Fit variable " << varName << " was not found!" << std::endl;  return false; }
  const auto& fitSet = RooArgSet(*fitVar);
  const auto& binWidth = fitVar->getBinWidth(0, "FitWindow");
  const auto& minDSRange = fitVar->getMin("PlotWindow");
  const auto& maxDSRange = ( (maxRng>0.0) ? maxRng : fitVar->getMax("PlotWindow") );
  const auto& nDSBins    = getNBins(minDSRange, maxDSRange, binWidth);
  fitVar->setRange("DSPlotWindow", minDSRange, maxDSRange);
  fitVar->setBins(nDSBins, "DSPlotWindow");
  setVarToTag(ws, varName, "DSPlotWindow");
  const auto& minPDFRange = std::max(fitVar->getMin("FitWindow"), minDSRange);
  const auto& maxPDFRange = std::min(fitVar->getMax("FitWindow"), maxDSRange);
  const auto& nPDFBins    = getNBins(minPDFRange, maxPDFRange, binWidth);
  fitVar->setRange("PDFPlotWindow", minPDFRange, maxPDFRange);
  fitVar->setBins(nPDFBins, "PDFPlotWindow");
  setVarToTag(ws, varName, "PDFPlotWindow");
  //
  const bool& useDS = (fitDS!=NULL);
  const bool& isMC = (DSTAG.rfind("MC",0)==0);
  const bool& isWeighted = (useDS ? ws.data(dsNameFit.c_str())->isWeighted() : false);
  //
  auto projSet = *fitPDF->getObservables(*fitDS);
  projSet.remove(*fitVar, true, true);
  if (projSet.find("Cand_Mass")!=NULL) { projSet.remove(*ws.var("Cand_Mass"), true, true); } // Remove mass variable
  const bool& useProjWData = (projSet.getSize() > 0);
  auto projDS = std::unique_ptr<RooAbsData>(useProjWData ? ws.data(dsNameFit.c_str())->reduce(RooFit::Name("projDS"), RooFit::SelectVars(projSet)) : new RooDataSet("projDS", "", RooArgSet()));
  const auto& NORM = (useProjWData ? projDS->sumEntries() : 1.0);
  if (useProjWData) { std::cout << "[INFO] Using projected dataset for " << fitPDF->GetName() << std::endl; }
  //
  int drawMode = 0;
  bool drawPull = true;  // false : Draw DATA/FIT , true : Draw the Pull
  //
  StringMap_t legInfo;
  RooPlotPtrMap_d frame;
  PadPtrMap_d pad; // Unique Pointer does produce Segmentation Fault, so don't use it
  //
  // Create the main plot of the fit
  //
  std::cout << "[INFO] Drawing " << varName << " plot for " << pdfName << " fitted on " << dsNameFit << std::endl;
  //
  const auto& frameName = "frame"+varTag+"_Tot"+label;
  const bool& reFrame = (ws.obj(frameName.c_str()) ? redoFrame : false);
  //
  if (!reFrame && ws.obj(frameName.c_str())) {
    frame["MAIN"] = std::unique_ptr<RooPlot>(dynamic_cast<RooPlot*>(ws.obj(frameName.c_str())));
  }
  else {
    // Create the frame
    if (useDS==false) { std::cout << "[ERROR] Dataset " << dsNameFit << " was not found!" << std::endl; return false; }
    frame["MAIN"] = std::unique_ptr<RooPlot>(fitVar->frame(RooFit::Range("DSPlotWindow")));
    // Check if we want to change PDFs to Histograms
    if (ws.pdf(pdfName.c_str()) && doHistPDF && !changeToHistPdf(ws, pdfName, varName)) { return false; }
    // Plot the data
    if (fitDS) {
      fitDS->plotOn(frame.at("MAIN").get(), RooFit::Name(("plot_"+dsNameFit).c_str()), RooFit::Binning("PDFPlotWindow"),
		    RooFit::MarkerColor(kBlack), RooFit::LineColor(kBlack), RooFit::MarkerSize(1.2));
    }
    // Plot the signal and background PDFs
    if (fitPDF) {
      // Check number of PDFs used
      const auto& fitPDFList = fitPDF->pdfList();
      // One PDF Case: Draw PDF subcomponents
      if (fitPDFList.getSize()==1) {
	std::string label = fitPDFList.at(0)->GetName(); label = label.substr(label.find("_")+1);
	const auto& pdf = dynamic_cast<RooAddPdf*>(ws.pdf(("pdf"+varTag+"_"+label).c_str()));
	if (pdf && pdf->pdfList().getSize()>1) {
	  const auto& norm = fitPDF->expectedEvents(fitSet);
	  const auto& pdfColor = std::vector<int>({kRed+1, kBlue+2, kGreen+3, kViolet+2, kAzure-7});
	  auto pdfIt = std::unique_ptr<TIterator>(pdf->pdfList().createIterator()); int iPdf = 0;
	  for (auto it = pdfIt->Next(); it!=NULL; it = pdfIt->Next(), iPdf++) {
	    RooArgList pdfL; pdfL.add(*dynamic_cast<RooAbsPdf*>(it));
	    ws.pdf(pdfName.c_str())->plotOn(frame.at("MAIN").get(), RooFit::Name(Form("plot_%s", it->GetName())), RooFit::Components(pdfL),
					    RooFit::Range("PDFPlotWindow"), RooFit::NormRange("PDFPlotWindow"), RooFit::Normalization(norm/NORM, RooAbsReal::NumEvent),
					    RooFit::LineColor(pdfColor[iPdf]), RooFit::LineStyle(2),
					    RooFit::ProjWData(projSet, *projDS, true), RooFit::NumCPU(32, 1));
	  }
	}
      }
      // Case: Multiple PDFs
      else {	
	// Extract the PDFs
	std::map< int , RooAbsPdf* > pdfStackMap, pdfDrawMap;
	auto pdfIt = std::unique_ptr<TIterator>(fitPDFList.createIterator());
	for (auto itp = pdfIt->Next(); itp!=NULL; itp = pdfIt->Next()) {
 	  const auto& it = dynamic_cast<RooAbsPdf*>(itp); if (!it) continue;
	  std::string obj = it->GetName(); obj = obj.substr(obj.find("_")+1); obj = obj.substr(0, obj.find(cha));
	  if (!contain(PDFMAP_, obj)) { std::cout << "[ERROR] Object " << obj << " is not defined in PDFMAP!" << std::endl; return false; }
	  // Store the PDFs
	  const bool& doStack = ((varName!="Cand_Mass" && varName!="Cand_DLen" && varName!="Cand_DLenErr") || obj.find("Swap")!=std::string::npos || obj=="Bkg");
	  if  (doStack) { pdfStackMap[PDFMAP_.at(obj)[0]] = it; }
	  else { pdfDrawMap[PDFMAP_.at(obj)[0]] = it; }
	}
	// Loop over the stacked PDFs
	if (pdfStackMap.size()>1) {
	  double norm = 99999999.999;
	  RooArgList pdfList; for (const auto& p : pdfStackMap) { pdfList.add(*p.second); }
	  for (const auto& elem : pdfStackMap) {
	    const auto& pdf = elem.second;
	    const std::string& pName = pdf->GetName();
	    std::string label = pName; label = label.substr(label.find("_")+1);
	    const auto& obj = label.substr(0, label.find(cha));
	    const auto& yield = (ws.var(("N_"+label).c_str()) ? ws.var(("N_"+label).c_str()) : ws.function(("N_"+label).c_str()));
	    if (!yield) { std::cout << "[ERROR] Yield N_" << label << " was not found!" << std::endl; return false; }
	    if (norm > 0.0) {
	      // Add the PDFs from the list
	      RooArgList coef, pdfs;
	      if(!addPdfToList(pdfs, coef, ws, pdfList)) { return false; }
	      const auto& pdfPlotName = std::string("pdfPlot")+(reFrame?"RE_":"_")+pName;
	      const auto& addPDF = RooAddPdf(pdfPlotName.c_str(), pdfPlotName.c_str(), pdfs, coef);
	      if (norm==99999999.999) { norm = addPDF.expectedEvents(fitSet); }
	      if (!ws.pdf(pdfPlotName.c_str()) && ws.import(addPDF)) { return false; }
	      // Plot the sum of PDFs
	      ws.pdf(pdfPlotName.c_str())->plotOn(frame.at("MAIN").get(), RooFit::Name(("plot_"+pName).c_str()), RooFit::Range("PDFPlotWindow"),
						  RooFit::Normalization(norm/NORM, RooAbsReal::NumEvent), RooFit::Precision(1e-7),
						  RooFit::FillStyle(1001), RooFit::FillColor(PDFMAP_.at(obj)[1]), RooFit::VLines(), RooFit::DrawOption("B"),
						  RooFit::ProjWData(projSet, *projDS, true), RooFit::NumCPU(32, 1));
	    }
	    norm -= yield->getVal();
	    pdfList.remove(*pdf);
	  }
	}
	else if (pdfStackMap.size()==1) {
	  const auto& pdf = pdfStackMap.begin()->second;
	  const auto& norm = dynamic_cast<RooExtendPdf*>(pdf)->expectedEvents(fitSet);
	  const std::string& pName = pdf->GetName();
	  std::string obj = pName; obj = obj.substr(obj.find("_")+1); obj = obj.substr(0, obj.find("To"));
	  // Plot the PDF
	  ws.pdf(pName.c_str())->plotOn(frame.at("MAIN").get(), RooFit::Name(("plot_"+pName).c_str()), RooFit::Range("PDFPlotWindow"),
					RooFit::Normalization(norm/NORM, RooAbsReal::NumEvent), RooFit::Precision(1e-7),
					RooFit::FillStyle(1001), RooFit::FillColor(PDFMAP_.at(obj)[1]), RooFit::VLines(), RooFit::DrawOption("B"),
					RooFit::ProjWData(projSet, *projDS, true), RooFit::NumCPU(32, 1));
	}
	// Loop over the non-stacked PDFs
	if (!pdfDrawMap.empty()) {
	  const auto& norm = fitPDF->expectedEvents(fitSet);
	  RooArgList stackPdfs; for (const auto& p : pdfStackMap) { stackPdfs.add(*p.second); }
	  for (const auto& elem : pdfDrawMap) {
	    const auto& pdf = elem.second;
	    const std::string& pName = pdf->GetName();
	    std::string obj = pName; obj = obj.substr(obj.find("_")+1); obj = obj.substr(0, obj.find(cha));
	    // Add the components
	    RooArgList pdfs; pdfs.add(*pdf);
	    // Plot the PDF
	    ws.pdf(pdfName.c_str())->plotOn(frame.at("MAIN").get(), RooFit::Name(("plot_"+pName).c_str()), RooFit::Components(pdfs),
					    RooFit::Range("PDFPlotWindow"), RooFit::NormRange("PDFPlotWindow"), RooFit::Normalization(norm/NORM, RooAbsReal::NumEvent),
					    RooFit::LineColor(PDFMAP_.at(obj)[1]), RooFit::LineStyle(2),
					    RooFit::ProjWData(projSet, *projDS, true), RooFit::NumCPU(32, 1));
	  }
	}
      }
    }
    // Plot the data
    if (ds) {
      ds->plotOn(frame.at("MAIN").get(), RooFit::Name(("plot_"+dsName).c_str()), RooFit::Binning("DSPlotWindow"),
		 RooFit::MarkerColor(kBlack), RooFit::LineColor(kBlack), RooFit::MarkerSize(1.2));
    }
    // Plot total PDF
    if (fitPDF) {
      const auto& norm = fitDS->sumEntries();
      fitPDF->plotOn(frame.at("MAIN").get(), RooFit::Name(("plot_"+pdfName).c_str()), RooFit::Range("PDFPlotWindow"), RooFit::NormRange("PDFPlotWindow"),
		     RooFit::Normalization(norm/NORM, RooAbsReal::NumEvent), RooFit::Precision(1e-7), RooFit::LineColor(kBlack), RooFit::LineStyle(1),
		     RooFit::ProjWData(projSet, *projDS, true), RooFit::NumCPU(32, 1));
    }
    // Store the frame
    frame.at("MAIN")->SetTitle(frameName.c_str());
    if (!ws.obj(frameName.c_str()) && ws.import(*frame.at("MAIN"), frame.at("MAIN")->GetTitle())) { return false; }
  }
  //
  // Create the extra frames and legend information
  if (fitDS) { legInfo["plot_"+dsNameFit] = (isMC ? "MC" : "Data"); }
  if (fitPDF) {
    legInfo["plot_"+pdfName] = "Fit";
    const auto& fitPDFList = fitPDF->pdfList();
    if (fitPDFList.getSize()==1) { drawPull = true; }
    else {
      auto pdfIt = std::unique_ptr<TIterator>(fitPDFList.createIterator());
      for (auto it = pdfIt->Next(); it!=NULL; it = pdfIt->Next()) {
	const std::string& pName = it->GetName();
	auto obj = pName; obj = obj.substr(obj.find("_")+1); obj = obj.substr(0, obj.find(cha));
        legInfo["plot_"+pName] = formatCut(obj);
      }
    }
    frame["EXTRA"] = std::unique_ptr<RooPlot>(dynamic_cast<RooPlot*>(frame.at("MAIN")->emptyClone("EXTRA")));
    RooHist* hExtra = new RooHist(2.0); // !!!DONT USE UNIQUE POINTER!!!!, 2 represents the binWidth of var
    if (drawPull) { if (!makePullHist (*hExtra, *frame.at("MAIN"), "", "", true)) { return false; } }
    else          { if (!makeRatioHist(*hExtra, *frame.at("MAIN"), "", "", true)) { return false; } }
    hExtra->SetName("hExtra");
    frame.at("EXTRA")->addPlotable(hExtra, "EP");
    drawMode = 1;
  }
  //
  // Create the canvas
  std::unique_ptr<TCanvas> cFig  = std::unique_ptr<TCanvas>(new TCanvas(("c"+varTag+"Fig_Tot"+label).c_str(), "cFig", 800, 800));
  cFig->cd();
  //
  // Edit the main frame
  if (contain(frame, "MAIN")) {
    //TGaxis::SetMaxDigits(4);
    frame.at("MAIN")->SetTitle("");
    frame.at("MAIN")->SetMarkerSize(1.5);
    frame.at("MAIN")->GetYaxis()->CenterTitle(kTRUE);
    frame.at("MAIN")->GetXaxis()->CenterTitle(kTRUE);
    frame.at("MAIN")->GetXaxis()->SetTitleOffset(1.1);
    frame.at("MAIN")->GetXaxis()->SetTitleSize(0.036);
    frame.at("MAIN")->GetXaxis()->SetLabelSize(0.033);
    if (drawMode==0) {
      frame.at("MAIN")->GetYaxis()->SetTitleOffset(1.5);
      frame.at("MAIN")->GetYaxis()->SetTitleSize(0.036);
      frame.at("MAIN")->GetYaxis()->SetLabelSize(0.033);
      frame.at("MAIN")->GetXaxis()->SetTitle(formatPar(varName, cha).c_str());
      pad["MAIN"] = new TPad(("padMAIN_Tot"+label).c_str(), "", 0, 0, 1, 1 );
    }
    else if (drawMode>0) {
      frame.at("MAIN")->GetYaxis()->SetTitleOffset(0.9);
      frame.at("MAIN")->GetYaxis()->SetTitleSize(0.060*(1./0.8));
      frame.at("MAIN")->GetYaxis()->SetLabelSize(0.033*(1./0.8));
      frame.at("MAIN")->GetXaxis()->SetTitleOffset(3);
      frame.at("MAIN")->GetXaxis()->SetLabelOffset(3);
      frame.at("MAIN")->GetXaxis()->SetTitle("");
      pad["MAIN"] = new TPad(("padMAIN_Tot"+label).c_str(), "", 0, 0.2, 1, 1 );
      pad.at("MAIN")->SetFixedAspectRatio(kTRUE);
      pad.at("MAIN")->SetBottomMargin(0.015);
      if (drawPull==false && (plotStyle==1 || plotStyle==2 || plotStyle==3)) { pad.at("MAIN")->SetBottomMargin(0.0); }
    }
  }
  // Edit the extra frame
  std::unique_ptr<TLine> pLine;
  if (contain(frame, "EXTRA")) {
    frame.at("EXTRA")->SetTitle("");
    frame.at("EXTRA")->SetMarkerSize(1.5);
    frame.at("EXTRA")->GetYaxis()->CenterTitle(kTRUE);
    frame.at("EXTRA")->GetXaxis()->CenterTitle(kTRUE);
    frame.at("EXTRA")->GetYaxis()->SetTitleOffset(0.25);
    frame.at("EXTRA")->GetYaxis()->SetTitleSize(0.25);
    frame.at("EXTRA")->GetYaxis()->SetLabelSize(0.15);
    frame.at("EXTRA")->GetYaxis()->SetNdivisions(204);
    if (drawPull) { frame.at("EXTRA")->GetYaxis()->SetTitle("Pull"); }
    else          { frame.at("EXTRA")->GetYaxis()->SetTitle("#frac{Data}{Fit}"); }
    frame.at("EXTRA")->GetXaxis()->SetTitleOffset(0.7);
    frame.at("EXTRA")->GetXaxis()->SetTitleSize(0.25);
    frame.at("EXTRA")->GetXaxis()->SetLabelSize(0.15);
    frame.at("EXTRA")->GetXaxis()->SetTitle(formatPar(varName, cha).c_str());
    if (drawPull) { frame.at("EXTRA")->GetYaxis()->SetRangeUser(-6.0, 6.0); }
    else if (plotStyle==1 || plotStyle==2 || plotStyle==3) { frame.at("EXTRA")->GetYaxis()->SetRangeUser(0., 2.5); }
    else { frame.at("EXTRA")->GetYaxis()->SetRangeUser(-0.01, 2.1); }
    pad["EXTRA"] = new TPad(("padEXTRA_Tot"+label).c_str(), "", 0, 0, 1, 0.2 );
    pad.at("EXTRA")->SetFixedAspectRatio(kTRUE);
    pad.at("EXTRA")->SetTopMargin(0.02);
    if (drawPull==false && (plotStyle==1 || plotStyle==2 || plotStyle==3)) { pad.at("EXTRA")->SetTopMargin(0.0); }
    pad.at("EXTRA")->SetBottomMargin(0.5);
    pad.at("EXTRA")->SetFillStyle(4000);
    pad.at("EXTRA")->SetFrameFillStyle(4000);
    // Draw the Pull
    pad.at("EXTRA")->Draw();
    pad.at("EXTRA")->cd();
    frame.at("EXTRA")->Draw();
    if (doGoF) { printGoF(*pad.at("EXTRA"), ws, *frame.at("MAIN"), varName.c_str(), dsNameFit, pdfName); }
    if (drawPull) { pLine = std::unique_ptr<TLine>(new TLine(frame.at("EXTRA")->GetXaxis()->GetXmin(), 0.0, frame.at("EXTRA")->GetXaxis()->GetXmax(), 0.0)); }
    else          { pLine = std::unique_ptr<TLine>(new TLine(frame.at("EXTRA")->GetXaxis()->GetXmin(), 1.0, frame.at("EXTRA")->GetXaxis()->GetXmax(), 1.0)); }
    pLine->Draw("same");
    pad.at("EXTRA")->Update();
  }
  //
  setPlotRange(*frame.at("MAIN"), ws, varName, dsName, setLogScale);
  //
  cFig->cd();
  pad.at("MAIN")->Draw();
  pad.at("MAIN")->cd();
  frame.at("MAIN")->Draw();
  //
  // Set the CMS style and add text info
  StringVector_t lumiLabels; getLumiLabels(lumiLabels, PD, col, isMC);
  CMS_lumi(pad.at("MAIN"), 33, lumiLabels[0], lumiLabels[1], false, 0.8, false);
  printCandidateTextInfo(*pad.at("MAIN"), ws, varName, drawMode, plotStyle);
  //
  // Draw the parameter results of the fit
  if (plotStyle==0) { printCandidateParameters(*pad.at("MAIN"), ws, pdfName, drawMode); }
  //
  // Draw the legend
  auto leg = std::unique_ptr<TLegend>(printCandidateLegend(*pad.at("MAIN"), *frame.at("MAIN"), legInfo, plotStyle));
  //
  // Set log scale if requested
  pad.at("MAIN")->SetLogy(setLogScale);
  pad.at("MAIN")->Update();
   //
  // Save the plot in different formats
  StringVector_t formats = {"png", "pdf", "root", "C"};
  for (const auto& f : formats) {
    makeDir(outputDir+"plot/"+varTag+"/"+f+"/");
    cFig->SaveAs((outputDir+"plot/"+varTag+"/"+f+"/"+fileName+"."+f).c_str());
  }
  // Close canvas
  cFig->Clear();
  cFig->Close();
  //
  // Return to fit values as default
  setVarToTag(ws, varName, "FitWindow");
  //
  return true;
};


bool drawCandidatePlot( RooWorkspace& ws,
			const std::string& fileName,
			const std::string& outputDir,
			const bool& yLogScale,
			const double& maxRng = -1.0,
			const bool& doGoF = true,
			const int&  plotStyle = 0,
			const bool& redoFrame = false,
			const bool& doHistPDF = false
			)
{
  const auto& fitVars = ws.set("fitVariable");
  if (!fitVars) { std::cout << "[ERROR] Set fitVariable was not found!" << std::endl; return false; }
  const auto& condVars = ws.set("condVariable");
  if (!condVars) { std::cout << "[ERROR] Set condVariable was not found!" << std::endl; return false; }
  auto set = *fitVars; set.add(*condVars);
  auto varIt = std::unique_ptr<TIterator>(set.createIterator());
  for (auto var = varIt->Next(); var!=NULL; var = varIt->Next() ) {
    std::string varTag = var->GetName(); stringReplace(varTag, "_", "");
    const auto& fileN = ("PLOT_"+varTag+"_"+fileName);
    if (!drawCandidatePlot(ws, var->GetName(), fileN, outputDir, yLogScale, maxRng, doGoF, plotStyle, redoFrame, doHistPDF)) { return false; }
  }
  return true;
};


void printCandidateParameters(TPad& pad, const RooWorkspace& ws, const std::string& pdfName, const uint& drawMode)
{
  pad.cd();
  float xPos = 0.69, yPos = 0.76, dYPos = 0.042, dy = 0.025;
  TLatex t = TLatex(); t.SetNDC(); t.SetTextSize(0.013);
  if (drawMode>0) { dy = 0.065; dYPos *= (1./0.8); t.SetTextSize(0.027*(1./0.8)); }
  // Get the parameters
  const auto& vars = getModelVar(ws, pdfName, "*");
  // Print the parameters
  for (const auto& v : vars) {
    const std::string& s = v.GetName();
    // Parse the parameter's labels
    const auto& parLbl = parseVarName(s); if (parLbl=="") continue;
    // Get number of decimals
    const int& n = std::max(-std::floor(std::log10(v.getError()>0. ? v.getError() : 1.)), 0.);
    // Print the parameter's results
    std::string txtLbl;
    if (s.rfind("N_",0)==0) { txtLbl = Form("%s = %.0f#pm%.0f", parLbl.c_str(), v.getValV(), v.getError()); }
    else if (n==0) { txtLbl = Form("%s = %.0f#pm%.0f", parLbl.c_str(), v.getValV(), v.getError()); }
    else if (n==1) { txtLbl = Form("%s = %.1f#pm%.1f", parLbl.c_str(), v.getValV(), v.getError()); }
    else if (n==2) { txtLbl = Form("%s = %.2f#pm%.2f", parLbl.c_str(), v.getValV(), v.getError()); }
    else           { txtLbl = Form("%s = %.3f#pm%.3f", parLbl.c_str(), v.getValV(), v.getError()); }
    const bool& isAtLimit = isParAtLimit(v);
    if (isAtLimit) { txtLbl += " (!)"; }
    const bool printTxt = (isAtLimit || (s.rfind("b_",0)==0) ||  (s.rfind("N_",0)==0) || (s.find("BkgTo")==std::string::npos));
    if (printTxt) { t.DrawLatex(xPos, yPos-dy, txtLbl.c_str()); dy+=dYPos; }
  }
  pad.Update();
  return;
};


void printCandidateTextInfo(TPad& pad, const RooWorkspace& ws, const std::string& fitVar, const uint& drawMode, const int& plotStyle)
{
  pad.cd();
  // Get information from workspace
  const std::string& cha       = (ws.obj("channel")   ? dynamic_cast<RooStringVar*>(ws.obj("channel")  )->getVal() : "");
  const std::string& col       = (ws.obj("fitSystem") ? dynamic_cast<RooStringVar*>(ws.obj("fitSystem"))->getVal() : "");
  const std::string& chg       = (ws.obj("fitCharge") ? dynamic_cast<RooStringVar*>(ws.obj("fitCharge"))->getVal() : "");
  const std::string& dsName    = (ws.obj("dsName")    ? dynamic_cast<RooStringVar*>(ws.obj("dsName")   )->getVal() : "");
  const std::string& dsNameFit = (ws.obj("dsNameFit") ? dynamic_cast<RooStringVar*>(ws.obj("dsNameFit"))->getVal() : "");
  //
  // Include the CMS labels
  TLatex t = TLatex(); t.SetNDC(); t.SetTextSize(0.030);
  double xPos = 0.20, yPos = 0.89, dYPos = 0.050, dy = 0.035;
  if (plotStyle==1 || plotStyle==2 || plotStyle==3) { xPos = 0.22; yPos = 0.87; dYPos = 0.055; dy = 0.035; }
  t.SetTextSize(0.058*1.25); t.SetTextFont(61); t.DrawLatex(0.78, 0.82, "CMS"); t.SetTextFont(62);
  if (plotStyle!=2) { t.SetTextSize(0.044*1.25); t.SetTextFont(52); t.DrawLatex(0.69, 0.75, "Preliminary"); t.SetTextFont(62); }
  //
  // Include the process
  std::set<std::string> objS;
  for (const auto& p : PDFMAP_) {
    if (p.first.rfind("Swap")!=std::string::npos) continue;
    if (ws.arg(("N_"+p.first+cha+chg+"_"+col).c_str())) { objS.insert(p.first); }
  }
  const auto& process = parseProcess(objS, cha, dsNameFit);
  if (drawMode>0) { dy *= (1./0.8); dYPos *= (1./0.8); t.SetTextSize(0.040*(1./0.8)); }
  if (plotStyle==1 || plotStyle==2 || plotStyle==3) { t.SetTextSize(0.055*1.25); t.DrawLatex(xPos, 0.82, process.c_str()); dy+=dYPos; }
  else { t.SetTextSize(0.040); t.DrawLatex(xPos, yPos-dy, process.c_str()); dy+=dYPos; }
  //
  // Extract the dataset variables
  std::vector<RooRealVar> dsVar;
  const auto& dsSet = const_cast<RooWorkspace*>(&ws)->set(("SET_"+dsNameFit).c_str());
  auto parIt = std::unique_ptr<TIterator>(dsSet->createIterator());
  for (auto itp = parIt->Next(); itp!=NULL; itp = parIt->Next()) {
    const auto& it = dynamic_cast<RooRealVar*>(itp); if (!it) continue;
    const std::string& varN = it->GetName();
    if (varN==fitVar || varN=="Cand_Mass" || varN.rfind("Cand_DLen",0)==0) continue;
    if (varN=="Centrality" && col.rfind("PbPb",0)!=0) continue;
    auto absVarN = varN; if (absVarN.find("_")!=std::string::npos) { absVarN.insert(absVarN.find("_")+1, "Abs"); }
    auto varCM = varN+"CM";
    if (ws.var(varCM.c_str())) { dsVar.push_back(*ws.var(varCM.c_str())); }
    if (ws.var(absVarN.c_str())) { dsVar.push_back(*ws.var(absVarN.c_str())); }
    if (!ws.var(varCM.c_str()) && !ws.var(absVarN.c_str())) { dsVar.push_back(*it); }
  }
  // Draw the dataset binning
  for (const auto& var : dsVar) {
    std::string varName = var.GetName();
    if (varName.rfind("N_", 0)==0 || varName.rfind("L_N_", 0)==0) continue;
    // Get the default range
    const auto& defaultMin = var.getMin("DEFAULT");
    const auto& defaultMax = var.getMax("DEFAULT");
    // Get the dataset range
    const auto& dsMin = var.getMin("DSWindow");
    const auto& dsMax = var.getMax("DSWindow");
    // Draw the binning
    if (!contain(VARLABEL_, varName)) { std::cout << "[WARNING] Parameter " << varName << " is not in VARLABEL!" << std::endl; continue; }
    if (ws.var(varName.c_str())) {
      const auto& varMin = ws.var(varName.c_str())->getMin();
      const auto& varMax = ws.var(varName.c_str())->getMax();
      const auto& varLbl = VARLABEL_.at(varName);
      std::string vUnit  = ws.var(varName.c_str())->getUnit(); if (vUnit!="") { vUnit = " "+vUnit; }
      if (varName.rfind("Cand_DLen",0)==0) {
	if (varMin>dsMin && varMax>=dsMax) { t.DrawLatex(xPos, yPos-dy, Form("%g%s #leq %s", varMin, vUnit.c_str(), varLbl.c_str())); dy+=dYPos; }
	else if (varMin<=dsMin && varMax<dsMax) { t.DrawLatex(xPos, yPos-dy, Form("%s < %g%s", varLbl.c_str(), varMax, vUnit.c_str())); dy+=dYPos; }
	else if (varMin>dsMin && varMax<dsMax) { t.DrawLatex(xPos, yPos-dy, Form("%g #leq %s < %g%s", varMin, varLbl.c_str(), varMax, vUnit.c_str())); dy+=dYPos; }
      }
      else if (varMin!=defaultMin && varMax==defaultMax) { t.DrawLatex(xPos, yPos-dy, Form("%g%s #leq %s", varMin, vUnit.c_str(), varLbl.c_str())); dy+=dYPos; }
      else if (varMin==defaultMin && varMax!=defaultMax) { t.DrawLatex(xPos, yPos-dy, Form("%s < %g%s", varLbl.c_str(), varMax, vUnit.c_str())); dy+=dYPos; }
      else if (varMin!=defaultMin && varMax!=defaultMax) { t.DrawLatex(xPos, yPos-dy, Form("%g #leq %s < %g%s", varMin, varLbl.c_str(), varMax, vUnit.c_str())); dy+=dYPos; }
    }
  }
  // Draw the extra information
  const auto& cutDS = dynamic_cast<RooStringVar*>(ws.obj(("CutAndCount_Tot"+cha+chg+"_"+col).c_str()));
  if (cutDS) { t.DrawLatex(xPos, yPos-dy, formatCut(cutDS->getVal()).c_str()); dy+=dYPos; }
  const std::string& cutSel = (ws.obj("cutSelStr") ? dynamic_cast<RooStringVar*>(ws.obj("cutSelStr"))->getVal() : "");
  if (cutSel!="") { t.DrawLatex(xPos, yPos-dy, (formatCutLbl(cutSel)+" cut").c_str()); dy+=dYPos; }
  if (ws.var("FAILED")) { t.DrawLatex(xPos, yPos-dy, "FIT FAILED"); dy+=dYPos; }
  //
  // Display the number of events lost if the dataset was reduced before fitting
  const auto& dsEntries = ws.var(("numEntries_"+dsName).c_str());
  const auto& fitDSEntries = ws.var(("numEntries_"+dsNameFit).c_str());
  if (fitVar!="Cand_DLenRes" && dsEntries && fitDSEntries) {
    const double& outTot = dsEntries->getVal();
    const double& outCut = fitDSEntries->getVal();
    if ((outTot-outCut) >= 0.5) { t.DrawLatex(xPos, yPos-dy, Form("Cut: (%.2f%%) %.0f evts", ((outTot-outCut)*100.0/outTot), (outTot-outCut))); }
  }
  // Update pad and return
  pad.Update();
  return;
};


TLegend* printCandidateLegend(TPad& pad, const RooPlot& frame, const StringMap_t& legInfo, const int& plotStyle)
{
  pad.cd();
  // Define the position and size of the legend
  double xmin = 0.49 , xmax = 0.66 , ymin = 0.58 , ymax = 0.89 , legSize = 0.047;
  if (plotStyle==1 || plotStyle==2 || plotStyle==3) {
    xmin = 0.74; xmax = 0.90; ymin = 0.25; ymax = 0.71; legSize = 0.06;
  }
  ymin = std::max(ymax-legInfo.size()*legSize*1.2, ymin);
  auto leg = new TLegend(xmin, ymin, xmax, ymax);
  // Define the legend order
  auto objMap = PDFMAP_; objMap["DATA"].push_back(-2); objMap["_TotTo"].push_back(-1);
  StringVector_t legOrder(objMap.size()+2, "");
  for (const auto& p : objMap) {
    for (const auto& l : legInfo) { if (l.first.find(p.first)!=std::string::npos) { legOrder[p.second[0]+2] = l.first; break; } }
  }
  // Fill the legend
  for (const auto& plotN : legOrder) {
    if (plotN=="") continue;
    const auto& objFrame = frame.findObject(plotN.c_str()); if (!objFrame) continue;
    // Find the draw option
    const std::string& drawOpt = frame.getDrawOptions(plotN.c_str()).Data();
    std::string legOpt = "pe"; if (drawOpt=="L") { legOpt = "l"; } else if (drawOpt=="F" || drawOpt=="B") { legOpt = "f"; }
    // Add the legend
    formatLegendEntry(*leg->AddEntry(objFrame, parseObject(legInfo.at(plotN)).c_str(), legOpt.c_str()), legSize);
  }
  // Draw the legend
  leg->Draw("same");
  // Update pad and return
  pad.Update();
  return leg;
};


#endif // #ifndef Candidate_drawCandidatePlot_C
