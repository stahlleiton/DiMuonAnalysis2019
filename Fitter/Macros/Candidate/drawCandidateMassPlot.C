#ifndef Candidate_drawCandidateMassPlot_C
#define Candidate_drawCandidateMassPlot_C

#include "RooWorkspace.h"

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "../Utilities/drawUtils.h"


void printCandidateMassParameters  ( TPad& pad , const RooWorkspace& ws , const std::string& pdfName, const std::string& varName , const uint& drawMode );
void printCandidateMassBinning     ( TPad& pad , const RooWorkspace& ws , const std::string& dsName , const std::vector< std::string >& text , const uint& drawMode , const int& plotStyle );
void printCandidateMassLegend      ( TPad& pad , TLegend& leg , const RooPlot& frame , const StringDiMap_d& legInfo , const double& size=0.05 );


bool drawCandidateMassPlot( RooWorkspace& ws,  // Local Workspace
                            // Select the type of datasets to fit
                            const std::string& fileName,
                            const std::string& outputDir,
                            const bool& yLogScale,
                            const double& maxRng = -1.0,
                            const bool& doGoF = true,
                            const int&  plotStyle = 0, // 3: Thesis , 2: Paper , 1: PAS , 0: AN
                                  bool redoFrame = false
                            )
{
  //
  // set the CMS style
  setTDRStyle();
  //
  const std::string& varName = "Cand_Mass";
  auto varType = varName; stringReplace(varType, "_", "");
  //
  const std::string& DSTAG = (ws.obj("DSTAG")    ) ? dynamic_cast<RooStringVar*>(ws.obj("DSTAG")   )->getVal() : "";
  const std::string& cha   = (ws.obj("channel")  ) ? dynamic_cast<RooStringVar*>(ws.obj("channel")  )->getVal() : "";
  const std::string& col   = (ws.obj("fitSystem")) ? dynamic_cast<RooStringVar*>(ws.obj("fitSystem"))->getVal() : "";
  const std::string& chg   = (ws.obj("fitCharge")) ? dynamic_cast<RooStringVar*>(ws.obj("fitCharge"))->getVal() : "";
  const std::string& obj   = (ws.obj("fitObject")) ? dynamic_cast<RooStringVar*>(ws.obj("fitObject"))->getVal() : "";
  //
  const std::string& tag = ( obj + cha + chg + "_" + col );
  const std::string& dsName = ( "d" + chg + "_" + DSTAG );
  const std::string& dsNameFit = ( (ws.data((dsName+"_FIT").c_str())!=NULL) ? (dsName+"_FIT") : dsName );
  const std::string& pdfName = Form("pdf%s_Tot%s", varType.c_str(), tag.c_str());
  const auto& setLogScale = yLogScale;
  //
  // Create the Range for Plotting
  const auto& binWidth = ws.var(varName.c_str())->getBinWidth(0);
  const auto& minRange = ws.var(varName.c_str())->getMin();
  const auto& maxRange = ( (maxRng>0.0) ? maxRng : ws.var(varName.c_str())->getMax() );
  const auto& nBins    = int(std::round((maxRange - minRange)/binWidth));
  ws.var(varName.c_str())->setRange("PlotWindow", minRange, maxRange);
  ws.var(varName.c_str())->setBins(nBins, "PlotWindow");
  // BUG FIX
  const int& oNBins = ws.var(varName.c_str())->getBins(); const double& oMinRange = ws.var(varName.c_str())->getMin(); const double& oMaxRange = ws.var(varName.c_str())->getMax();
  ws.var(varName.c_str())->setBins(nBins); ws.var(varName.c_str())->setRange(minRange, maxRange);
  //
  const bool& useDS = (ws.data(dsName.c_str())!=NULL);
  const bool& isMC = (DSTAG.rfind("MC",0)==0);
  const bool& isWeighted = (useDS ? ws.data(dsName.c_str())->isWeighted() : false);
  int drawMode = 0;
  bool drawPull = false;  // false : Draw DATA/FIT , true : Draw the Pull
  //
  // Format Object name
  std::string process = "";
  std::string chgL = " "; if (chg=="Pl") { chgL = "+"; } else if (chg=="Mi") { chgL = "#font[122]{\55}"; }
  if      (obj=="WToTau" ) { process = Form("W^{%s}#rightarrow#tau^{%s}", chgL.c_str(), chgL.c_str()); }
  else if (obj=="DYToTau") { process = Form("Z/#gamma*#rightarrow#tau^{%s}", chgL.c_str()); }
  else if (obj=="W"      ) { process = Form("W^{%s}" , chgL.c_str()); }
  else if (obj=="WZ"     ) { process = Form("W^{%s}Z", chgL.c_str()); }
  else if (obj=="WW"     ) { process = "W^{+}W^{-}"; }
  else if (obj=="DY"     ) { process = "Z/#gamma*";  }
  else if (obj=="Z"      ) { process = "Z";          }
  else if (obj=="ZZ"     ) { process = "ZZ";         }
  else if (obj=="QCD"    ) { process = "QCD";        }
  else if (obj=="TTbar"  ) { process = "t#bar{t}";   }
  else if (obj=="JPsi"   ) { process = "J/#psi";     }
  else if (obj=="Psi2S"  ) { process = "#psi(2S)";   }
  else if (obj=="Ups1S"  ) { process = "#Upsilon(1S)"; }
  else if (obj=="Ups2S"  ) { process = "#Upsilon(2S)"; }
  else if (obj=="Ups3S"  ) { process = "#Upsilon(3S)"; }
  else if (obj=="Bkg"    ) { process = "Bkg";        }
  else { process = "#Upsilon"; }
  if (cha=="ToMuMu") { process += " #rightarrow #mu^{+} + #mu^{-}"; }
  else if (cha=="ToMu") {
    if (obj=="W") {
      if (plotStyle < 3) { process += Form("#kern[0.2]{#rightarrow}#kern[0.2]{#mu^{%s}}#kern[0.2]{%s}", chgL.c_str(), (chg=="Pl"?"#nu_{#mu}":"#bar{#nu}_{#mu}")); }
      else { process += Form(" #rightarrow #mu^{%s} + %s", chgL.c_str(), (chg=="Pl"?"#nu_{#mu}":"#bar{#nu}_{#mu}")); }
    }
    else { process += Form(" #rightarrow #mu^{%s} + x", chgL.c_str()); }
  }
  if (obj=="DY" ) { process = "Z/#gamma* #rightarrow #mu^{+} + #mu^{-}"; }
  if (obj=="Z"  ) { process = "Z #rightarrow #mu^{+} + #mu^{-}"; }
  process = Form("#font[62]{#scale[1.1]{%s}}", process.c_str());
  //
  StringDiMap_d legInfo;
  //
  RooPlotPtrMap_d frame;
  PadPtrMap_d pad; // Unique Pointer does produce Segmentation Fault, so don't use it
  //
  // Create the main plot of the fit
  //
  const std::string& frameName = Form("frame_Tot%s", tag.c_str());
  if (ws.obj(frameName.c_str())==NULL) { redoFrame = false; }
  //
  if (!redoFrame && ws.obj(frameName.c_str())!=NULL) { frame["MAIN"] = std::unique_ptr<RooPlot>(dynamic_cast<RooPlot*>(ws.obj(frameName.c_str()))); }
  else {
    if (useDS==false) { std::cout << "[ERROR] Dataset " << dsName << " was not found!" << std::endl; return false; }
    frame["MAIN"] = std::unique_ptr<RooPlot>(ws.var(varName.c_str())->frame( RooFit::Range("PlotWindow") ));
    if (ws.data(("CutAndCount_"+dsName).c_str())) {
      ws.data(("CutAndCount_"+dsName).c_str())->plotOn(frame.at("MAIN").get(), RooFit::Name(Form("plot_Tot%s", dsName.c_str())), RooFit::Binning("PlotWindow"),
                                                       RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0),
                                                       RooFit::MarkerColor(kBlack), RooFit::LineColor(kBlack), RooFit::MarkerSize(1.2));
    }
    else {
      ws.data(dsName.c_str())->plotOn(frame.at("MAIN").get(), RooFit::Name(Form("plot_Tot%s", dsName.c_str())), RooFit::Binning("PlotWindow"),
                                      RooFit::MarkerColor(kBlack), RooFit::LineColor(kBlack), RooFit::MarkerSize(1.2));
    }
    //
    if (ws.pdf(pdfName.c_str())) {
      RooArgList pdfList = dynamic_cast<RooAddPdf*>(ws.pdf(pdfName.c_str()))->pdfList();
      if (pdfList.getSize()==1) {
        const double& norm = dynamic_cast<RooAddPdf*>(ws.pdf(pdfName.c_str()))->expectedEvents(RooArgSet(*ws.var(varName.c_str())));
        ws.pdf(pdfName.c_str())->plotOn(frame.at("MAIN").get(), RooFit::Name(Form("plot_%s", pdfName.c_str())), RooFit::Range("PlotWindow"), RooFit::NormRange("PlotWindow"),
                                        RooFit::Normalization(norm, RooAbsReal::NumEvent), RooFit::Precision(1e-7),
                                        RooFit::LineColor(kBlack), RooFit::LineStyle(1)
                                        );
      }
      else {
        //
        const double& norm = ws.data(dsNameFit.c_str())->sumEntries();
        //
        std::string pdfNameTot = Form("pdfPlot%s_pdf%sTot_%s%s%s_%s", (redoFrame?"RE":""), varName.c_str(), obj.c_str(), cha.c_str(), chg.c_str(), col.c_str());
        if (ws.pdf(pdfNameTot.c_str())==NULL) { pdfNameTot = pdfName; }
        ws.pdf(pdfNameTot.c_str())->plotOn(frame.at("MAIN").get(), RooFit::Name(Form("plot_%s", pdfName.c_str())), RooFit::Range("PlotWindow"), RooFit::NormRange("PlotWindow"),
                                           RooFit::Normalization(norm, RooAbsReal::NumEvent), RooFit::Precision(1e-7),
                                           RooFit::LineColor(kBlack), RooFit::LineStyle(1)
                                           );
      }
    }
    // Store the frame
    frame.at("MAIN")->SetTitle(frameName.c_str());
    if (ws.obj(frameName.c_str())==NULL) { ws.import(*frame.at("MAIN"), frame.at("MAIN")->GetTitle()); }
  }
  return false;
  //
  legInfo["DATA"][Form("plot_Tot%s", dsName.c_str())] = ( isMC ? "Simulation" : "Data" );
  //
  if (ws.pdf(pdfName.c_str())) {
    RooArgList pdfList = dynamic_cast<RooAddPdf*>(ws.pdf(pdfName.c_str()))->pdfList();
    if (pdfList.getSize()==1) {
      legInfo["PDF"][Form("plot_%s", pdfName.c_str())] = "Fit";
      frame["EXTRA"] = std::unique_ptr<RooPlot>(dynamic_cast<RooPlot*>(frame.at("MAIN")->emptyClone("EXTRA")));
      const auto& hPull = new RooHist(2.0); // !!!DONT USE UNIQUE POINTER!!!!, 2 represents the binWidth of var
      if (!makePullHist (*hPull, *frame.at("MAIN").get(), "", "", true)) { return false; }; drawPull = true;
      hPull->SetName("hPull");
      frame.at("EXTRA")->addPlotable(hPull, "EP");
      drawMode = 1;
    }
    else {
      const std::vector< std::string > pdfOrder = { "W" , "QCD" , "DY" , "WToTau" , "DYToTau" , "TTbar" };
      for (const auto& pdfT : pdfOrder) {
        const auto& obj = pdfT;
        const std::string& name = Form("pdf%sTot_%s%s%s_%s", varName.c_str(), obj.c_str(), cha.c_str(), chg.c_str(), col.c_str());
        //legInfo["TEMP"][Form("plot_%s", name.c_str())] = formatCut(obj);
      }
      frame["EXTRA"] = std::unique_ptr<RooPlot>(dynamic_cast<RooPlot*>(frame.at("MAIN")->emptyClone("EXTRA")));
      RooHist* hExtra = new RooHist(2.0); // !!!DONT USE UNIQUE POINTER!!!!, 2 represents the binWidth of var
      if (drawPull) { if (!makePullHist (*hExtra, *frame.at("MAIN"), "", "", true)) { return false; } }
      else          { if (!makeRatioHist(*hExtra, *frame.at("MAIN"), "", "", true)) { return false; } }
      hExtra->SetName("hExtra");
      frame.at("EXTRA")->addPlotable(hExtra, "EP");
      drawMode = 1;
    }
  }
  //
  std::unique_ptr<TCanvas> cFig  = std::unique_ptr<TCanvas>(new TCanvas( Form("c%sFig_Tot%s", varType.c_str(), tag.c_str()), "cFig", 800, 800 ));
  cFig->cd();
  //
  std::unique_ptr<TLine> pLine;
  if (drawMode==0) {
    //TGaxis::SetMaxDigits(4);
    // Main Frame
    frame.at("MAIN")->SetTitle("");
    frame.at("MAIN")->SetMarkerSize(1.5);
    frame.at("MAIN")->GetYaxis()->CenterTitle(kTRUE);
    frame.at("MAIN")->GetXaxis()->CenterTitle(kTRUE);
    frame.at("MAIN")->GetYaxis()->SetTitleOffset(1.5);
    frame.at("MAIN")->GetYaxis()->SetTitleSize(0.036);
    frame.at("MAIN")->GetYaxis()->SetLabelSize(0.033);
    frame.at("MAIN")->GetXaxis()->SetTitleOffset(1.1);
    frame.at("MAIN")->GetXaxis()->SetTitleSize(0.036);
    frame.at("MAIN")->GetXaxis()->SetLabelSize(0.033);
    frame.at("MAIN")->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
    pad["MAIN"] = new TPad( Form("padMAIN_Tot%s", tag.c_str()), "", 0, 0, 1, 1 );
  }
  else if (drawMode>0) {
    //TGaxis::SetMaxDigits(4);
    // Main Frame
    frame.at("MAIN")->SetTitle("");
    frame.at("MAIN")->SetMarkerSize(1.5);
    frame.at("MAIN")->GetYaxis()->CenterTitle(kTRUE);
    frame.at("MAIN")->GetXaxis()->CenterTitle(kTRUE);
    frame.at("MAIN")->GetYaxis()->SetTitleOffset(0.9);
    frame.at("MAIN")->GetYaxis()->SetTitleSize(0.060*(1./0.8));
    frame.at("MAIN")->GetYaxis()->SetLabelSize(0.033*(1./0.8));
    frame.at("MAIN")->GetXaxis()->SetTitleOffset(1.1);
    frame.at("MAIN")->GetXaxis()->SetTitleSize(0.036);
    frame.at("MAIN")->GetXaxis()->SetLabelSize(0.033);
    frame.at("MAIN")->GetXaxis()->SetTitleOffset(3);
    frame.at("MAIN")->GetXaxis()->SetLabelOffset(3);
    frame.at("MAIN")->GetXaxis()->SetTitle("");
    pad["MAIN"] = new TPad( Form("padMAIN_Tot%s", tag.c_str()), "", 0, 0.2, 1, 1 );
    pad.at("MAIN")->SetFixedAspectRatio(kTRUE);
    pad.at("MAIN")->SetBottomMargin(0.015);
    if (drawPull==false && (plotStyle==1 || plotStyle==2 || plotStyle==3)) { pad.at("MAIN")->SetBottomMargin(0.0); }
  }
  if (drawMode==1) {
    // Pull Frame
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
    frame.at("EXTRA")->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
    if (drawPull) { frame.at("EXTRA")->GetYaxis()->SetRangeUser(-6.0, 6.0); }
    else if (plotStyle==1 || plotStyle==2 || plotStyle==3) { frame.at("EXTRA")->GetYaxis()->SetRangeUser(0., 2.5); }
    else { frame.at("EXTRA")->GetYaxis()->SetRangeUser(-0.01, 2.1); }
    pad["EXTRA"] = new TPad( Form("padEXTRA_Tot%s", tag.c_str()), "", 0, 0, 1, 0.2 );
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
    if (doGoF) { printGoF(*pad.at("EXTRA"), ws, *frame.at("MAIN"), varName.c_str(), dsName, pdfName); }
    if (drawPull) { pLine = std::unique_ptr<TLine>(new TLine(frame.at("EXTRA")->GetXaxis()->GetXmin(), 0.0, frame.at("EXTRA")->GetXaxis()->GetXmax(), 0.0)); }
    else          { pLine = std::unique_ptr<TLine>(new TLine(frame.at("EXTRA")->GetXaxis()->GetXmin(), 1.0, frame.at("EXTRA")->GetXaxis()->GetXmax(), 1.0)); }
    pLine->Draw("same");
    pad.at("EXTRA")->Update();
  }
  //
  setPlotRange(*frame.at("MAIN"), ws, varName, dsName, setLogScale, ws.var(varName.c_str())->getBins("PlotWindow"));
  //
  cFig->cd();
  pad.at("MAIN")->Draw();
  pad.at("MAIN")->cd();
  frame.at("MAIN")->Draw();
  //
  int lumiId = 213;
  CMS_lumi(pad.at("MAIN"), lumiId, 33, "", false, 0.8, false);
  //
  if (plotStyle==0) { printCandidateMassParameters(*pad.at("MAIN"), ws, pdfName, varName.c_str(), drawMode); }
  std::vector< std::string > text = { process };
  std::cout << "A 4" << std::endl;
  //
  printCandidateMassBinning(*pad.at("MAIN"), ws, dsName, text, drawMode, plotStyle);
  //
  double xmin = 0.49 , xmax = 0.66 , ymin = 0.58 , ymax = 0.89;
  if (maxRng>0. && maxRng<120.) { ymax = 0.72; xmax = 0.40; ymin = 0.49; }
  double legSize = 0.047 , dy = (ymax-ymin);
  if (plotStyle==1 || plotStyle==2 || plotStyle==3) {
    if (legInfo.size()>2) { xmin = 0.74; ymin = 0.25; xmax = 0.90; ymax = 0.71; legSize = 0.06; }
    if (legInfo.size()<3) { xmin = 0.74; ymin = 0.35; xmax = 0.90; ymax = 0.71; legSize = 0.06; }
  }
  std::cout << "A 5" << std::endl;
  auto leg = std::unique_ptr<TLegend>(new TLegend(xmin, ymin, xmax, ymax));
  if (drawMode>0) { dy *= (1./0.8); }
  printCandidateMassLegend(*pad.at("MAIN"), *leg, *frame.at("MAIN"), legInfo, legSize);
  std::cout << "A 6" << std::endl;
  //
  pad.at("MAIN")->SetLogy(setLogScale);
  pad.at("MAIN")->Update();
  std::cout << "A 7" << std::endl;
  //
  // Save the plot in different formats
  gSystem->mkdir(Form("%splot/C/", outputDir.c_str()), kTRUE);
  cFig->SaveAs(Form("%splot/C/%s.C", outputDir.c_str(), fileName.c_str()));
  gSystem->mkdir(Form("%splot/pdf/", outputDir.c_str()), kTRUE);
  cFig->SaveAs(Form("%splot/pdf/%s.pdf", outputDir.c_str(), fileName.c_str()));
  gSystem->mkdir(Form("%splot/png/", outputDir.c_str()), kTRUE);
  cFig->SaveAs(Form("%splot/png/%s.png", outputDir.c_str(), fileName.c_str()));
  std::cout << "A 8" << std::endl;
  //
  cFig->Clear();
  cFig->Close();
  //
  // Undo the changes
  ws.var(varName.c_str())->setBins(oNBins); ws.var(varName.c_str())->setRange(oMinRange, oMaxRange);
  //
  return true;
};


void printCandidateMassParameters(TPad& pad, const RooWorkspace& ws, const std::string& pdfName, const std::string& varName, const uint& drawMode)
{
  pad.cd();
  float xPos = 0.69, yPos = 0.74, dYPos = 0.050, dy = 0.025;
  TLatex t = TLatex(); t.SetNDC(); t.SetTextSize(0.015);
  if (drawMode>0) { dy = 0.065; dYPos *= (1./0.8); t.SetTextSize(0.030*(1./0.8)); }
  std::vector<RooRealVar> vars; std::string label;
  if (vars.size()==0) {
    if (getVar(vars, ws, "N_", pdfName)) {
      for (const auto& v : vars) { parseVarName(v.GetName(), label); if(label!="") { t.DrawLatex(xPos, yPos-dy, Form("%s = %.0f#pm%.0f", label.c_str(), v.getValV(), v.getError())); dy+=dYPos; } }
    }
  }
  std::unique_ptr<TIterator> parIt;
  if (ws.pdf(pdfName.c_str())) {
    auto parList = std::unique_ptr<RooArgSet>(ws.pdf(pdfName.c_str())->getParameters(RooArgSet(*ws.var(varName.c_str()))));
    parIt = std::unique_ptr<TIterator>(parList->selectByAttrib("Constant", kFALSE)->createIterator());
  }
  else { parIt = std::unique_ptr<TIterator>(ws.allVars().selectByAttrib("Constant", kFALSE)->createIterator()); }
  for (auto itp = parIt->Next(); itp!=NULL; itp = parIt->Next()) {
    const auto& it = dynamic_cast<RooRealVar*>(itp); if (!it) continue;
    // Parse the parameter's labels
    std::string label="", s(it->GetName());
    // Ignore dataset variables
    if(s=="MET" || s=="Muon_Pt" || s=="Muon_Eta" || s=="Muon_Iso" || s=="Muon_MT" || s=="Event_Type" || s=="Centrality") continue;
    if(s=="Cand_Mass" || s=="Cand_Pt" || s=="Cand_Rap" || s=="Cand_AbsRap" || s=="Cand_Ctau") continue;
    if((s.find("Pl")!=std::string::npos)!=(pdfName.find("Pl")!=std::string::npos)) continue;
    if(s.rfind("N_",0)==0) continue;
    parseVarName(it->GetName(), label); if (label=="") continue;
    // Print the parameter's results
    std::string txtLbl;
    txtLbl = Form("%s = %.3f#pm%.3f", label.c_str(), it->getValV(), it->getError());
    if (isParAtLimit(*it)) { txtLbl += " (!)"; }
    t.DrawLatex(xPos, yPos-dy, txtLbl.c_str()); dy+=dYPos;
  }
  pad.Update();
  return;
};


void printCandidateMassBinning(TPad& pad, const RooWorkspace& ws, const std::string& dsName, const std::vector< std::string >& text, const uint& drawMode, const int& plotStyle)
{
  pad.cd();
  TLatex t = TLatex(); t.SetNDC(); t.SetTextSize(0.030);
  double xPos = 0.20, yPos = 0.89, dYPos = 0.050, dy = 0.035;
  if (plotStyle==1 || plotStyle==2 || plotStyle==3) { xPos = 0.22; yPos = 0.87; dYPos = 0.055; dy = 0.035; }
  t.SetTextSize(0.058*1.25); t.SetTextFont(61); t.DrawLatex(0.78, 0.82, "CMS"); t.SetTextFont(62);
  if (plotStyle!=2) { t.SetTextSize(0.044*1.25); t.SetTextFont(52); t.DrawLatex(0.69, 0.75, "Preliminary"); t.SetTextFont(62); }
  if (drawMode>0) { dy *= (1./0.8); dYPos *= (1./0.8); t.SetTextSize(0.040*(1./0.8)); }
  if (plotStyle==1 || plotStyle==2 || plotStyle==3) {
    t.SetTextSize(0.055*1.25); t.DrawLatex(xPos, 0.82, Form("%s", text[0].c_str())); dy+=dYPos;
  }
  else {
    t.SetTextSize(0.040); t.DrawLatex(xPos, yPos-dy, Form("%s", text[0].c_str())); dy+=dYPos;
  }
  std::vector<std::string> varNameList;
  if (ws.data(dsName.c_str())!=NULL) {
    auto parIt = std::unique_ptr<TIterator>(dynamic_cast<RooDataSet*>(ws.data(dsName.c_str()))->get()->createIterator());
    for (auto it = parIt->Next(); it!=NULL; it = parIt->Next()) {
      if (std::string(it->GetName())=="Cand_Mass") continue;
      varNameList.push_back(it->GetName());
    }
  }
  else { varNameList = std::vector<std::string>({ "Cand_Rap" , "Cand_AbsRap" , "Cand_Pt" , "Centrality"  }); }
  for (const auto& varName : varNameList) {
    double defaultMin = 0.0 , defaultMax = 100000.0;
    if (varName=="Muon_Eta") { defaultMin = -2.5; defaultMax = 2.5; }
    if (varName=="Muon_Iso") { defaultMin = 0.0; }
    if (varName=="Cand_Rap") { defaultMin = -2.5; defaultMax = 2.5; }
    if (varName=="Cand_AbsRap") { defaultMin = 0.0; defaultMax = 2.5; }
    if (varName=="Cand_Ctau") { defaultMin = -100000.0; }
    if (varName=="Centrality") { defaultMax = 100.0; }
    if (ws.var(varName.c_str())) {
      double minVal = ws.var(varName.c_str())->getMin();
      double maxVal = ws.var(varName.c_str())->getMax();
      string fVarName = "";//varLabel.at(varName);
      //
      if (minVal!=defaultMin && maxVal==defaultMax) {
        t.DrawLatex(xPos, yPos-dy, Form("%g #leq %s", minVal, fVarName.c_str())); dy+=dYPos;
      }
      if (minVal==defaultMin && maxVal!=defaultMax) {
        t.DrawLatex(xPos, yPos-dy, Form("%s < %g", fVarName.c_str(), maxVal)); dy+=dYPos;
      }
      if (minVal!=defaultMin && maxVal!=defaultMax) {
        t.DrawLatex(xPos, yPos-dy, Form("%g #leq %s < %g", minVal, fVarName.c_str(), maxVal)); dy+=dYPos;
      }
    }
  }
  for (const auto& txt : text) { if (text[0]!=txt) { t.DrawLatex(xPos, yPos-dy, Form("%s", txt.c_str())); dy+=dYPos; } }
  //
  if (ws.data(dsName.c_str())!=NULL) {
    const double& outTot = ws.data(dsName.c_str())->sumEntries();
    const double& outCut = ( (ws.data((dsName+"_FIT").c_str())!=NULL) ? ws.data((dsName+"_FIT").c_str())->sumEntries() : outTot );
    if (outCut != outTot) { t.DrawLatex(xPos, yPos-dy, Form("Loss: (%.4f%%) %.0f evts", ((outTot-outCut)*100.0/outTot), (outTot-outCut))); }
  }
  //
  pad.Update();
  return;
};


void printCandidateMassLegend(TPad& pad, TLegend& leg, const RooPlot& frame, const StringDiMap_d& legInfo, const double& size)
{
  pad.cd();
  StringMap_d drawOption = { { "DATA" , "pe" } , { "PDF" , "l" } , { "TEMP" , "f" } };
  const std::vector< std::string > pdfMapOrder = { "WToMu" , "QCD" , "DY" , "WToTau" , "DYToTau" , "TTbar" };
  for (const auto& map : legInfo) {
    if (map.first=="TEMP") {
      for (const auto& pdfM : pdfMapOrder) {
        std::pair< std::string , std::string > elem;
        for (const auto& el : map.second) { if (el.first.find(pdfM)!=std::string::npos) { elem = std::make_pair( el.first , el.second ); break; } }
        if (frame.findObject(elem.first.c_str())) { formatLegendEntry(*leg.AddEntry(frame.findObject(elem.first.c_str()), elem.second.c_str(), drawOption[map.first].c_str()), size); }
      }
    }
    else {
      for (const auto& elem : map.second) {
        if (frame.findObject(elem.first.c_str())) { formatLegendEntry(*leg.AddEntry(frame.findObject(elem.first.c_str()), elem.second.c_str(), drawOption[map.first].c_str()), size); }
      }
    }
  }
  leg.Draw("same");
  pad.Update();
  return;
};


#endif // #ifndef Candidate_drawCandidateMassPlot_C
