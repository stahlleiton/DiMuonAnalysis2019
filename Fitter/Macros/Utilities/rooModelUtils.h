#ifndef rooModelUtils_h
#define rooModelUtils_h

#include "TH1.h"
#include "TInterpreter.h"

#include "RooFit.h"
#include "RooWorkspace.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooHistPdf.h"

#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <set>

#include "rooDataUtils.h"
#include "RooStats/SPlot.h"


std::set<std::string> MODELCLASS_;


bool importModelClass(RooWorkspace& ws, const std::string& className)
{
  const std::string& CWD = getcwd(NULL, 0);
  const auto& classDir = CWD+"/Macros/Utilities/Models/";
  TInterpreter::EErrorCode ecode;
  if (MODELCLASS_.find(className)==MODELCLASS_.end()) {
    gInterpreter->ProcessLineSynch(Form(".L %s%s.cxx+",classDir.c_str(), className.c_str()), &ecode);
    if (ecode!=TInterpreter::kNoError) return false;
    MODELCLASS_.insert(className);
  }
  auto classPdf = std::unique_ptr<RooAbsPdf>((RooAbsPdf*)gInterpreter->ProcessLineSynch(Form("new %s()", className.c_str()), &ecode));
  if (ecode!=TInterpreter::kNoError) return false;
  return ws.importClassCode(classPdf->IsA());
};


void constrainQuarkoniumMassParameters(GlobalInfo& info, const StringVector_t& varList, const std::string& excLabel)
{
  const auto& cha = info.Par.at("channel");
  if (excLabel.find(cha)==std::string::npos) { std::cout << "[ERROR] Invalid channel " << cha << " in label " << excLabel << std::endl; return; }
  const auto& obj = excLabel.substr(0, excLabel.find(cha));
  const auto& col = excLabel.substr(excLabel.find("_")+1);
  const auto& chg = excLabel.substr(excLabel.find(cha)+cha.size(), 2);
  if (obj!="Psi2S" && obj!="Ups2S" && obj!="Ups3S") return;
  std::cout << "[INFO] Constraining " << obj << " mass model parameters to 1S state using PDG mass ratio!" << std::endl;
  const auto& refState = StringMap_t({{"Psi2S" , "JPsi"}, {"Ups2S" , "Ups1S"}, {"Ups3S" , "Ups1S"}});
  const auto& refLabel = (refState.at(obj)+cha+chg+"_"+col);
  for (const auto& v : varList) {
    if (contain(info.Par, v+"_"+excLabel)) {
      if (v=="m" || v=="Sigma1" || v=="Sigma2") {
	const double& massRatioValue = (MASS.at(obj).at("Val")/MASS.at(refState.at(obj)).at("Val"));
	const auto& massRatioLabel = ("MassRatio_" + obj + "Over" + refState.at(obj));
	info.Par[massRatioLabel] = Form("%s[%.10f]", massRatioLabel.c_str(), massRatioValue);
	info.Par.at(v+"_"+excLabel) = Form("RooFormulaVar::%s('@0*@1',{%s,%s})", (v+"_"+excLabel).c_str(), info.Par.at(massRatioLabel).c_str(), (v+"_"+refLabel).c_str());
      }
      else if (v=="N") {
	const auto& rN = "R_"+excLabel;
	info.Par[rN] = rN+"[0.4,-0.1,1.0]";
	info.Par.at(v+"_"+excLabel) = Form("RooFormulaVar::%s('@0*@1',{%s,%s})", (v+"_"+excLabel).c_str(), info.Par.at(rN).c_str(), (v+"_"+refLabel).c_str());
      }
      else { info.Par.at(v+"_"+excLabel) = Form("RooFormulaVar::%s('@0',{%s})", (v+"_"+excLabel).c_str(), (v+"_"+refLabel).c_str()); }
    }
  }
};


bool setModelPar(GlobalInfo& info, const StringVector_t& parNames, const std::string& label, const std::string& fitVar, const std::string& modelN)
{
  if (modelN!="") std::cout << "[INFO] Setting " << fitVar << " model parameters for " << modelN << std::endl;
  //
  const auto& cha = info.Par.at("channel");
  if (label.find(cha)==std::string::npos) { std::cout << "[ERROR] Invalid channel " << cha << " in label " << label << std::endl; return false; }
  const auto& obj = label.substr(0, label.find(cha));
  const auto& col = label.substr(label.find("_")+1);
  const auto& chg = label.substr(label.find(cha)+cha.size(), 2);
  //
  const auto& objLabel = obj + cha  + chg  + (col!="" ? "_"+col : "");
  auto vTmp = fitVar; stringReplace(vTmp, "_", "");
  auto objFoundLabel = findLabel("Model", vTmp, obj, chg, col, cha, info);
  if (vTmp!="") { objFoundLabel.erase(0, vTmp.size()+1); } else { objFoundLabel.erase(0, 1); }
  //
  for (const auto& v : parNames) {
    const bool& found = (contain(info.Par, v+"_"+objLabel) || contain(info.Par, v+"_"+objFoundLabel));
    if (!contain(info.Par, v+"_"+objLabel) || info.Par.at(v+"_"+objLabel)=="") {
      if (!contain(info.Par, v+"_"+objFoundLabel) || info.Par.at(v+"_"+objFoundLabel)=="") {
	if (v=="Cut") { info.Par[v+"_"+objLabel] = ""; }
	else if (v=="N") {
	  const auto& numEntries = info.Var.at("numEntries").at(chg);
	  info.Par[v+"_"+objLabel] = Form("%s[%.10f,%.10f,%.10f]", ("N_"+objLabel).c_str(), numEntries, -200.0, 2.0*numEntries);
	}
	else if (v.rfind("rSigma", 0)==0) { info.Par[v+"_"+objLabel] = Form("%s[%.6f,%.6f,%.6f]", (v+"_"+objLabel).c_str(),  2.000,    1.000,   8.000); }
	else if (v=="Alpha2"   && !found) { info.Par[v+"_"+objLabel] = Form("RooFormulaVar::%s('@0',{%s})", (v+"_"+objLabel).c_str(), ("Alpha_"+objLabel).c_str());  }
	else if (v=="AlphaR2"  && !found) { info.Par[v+"_"+objLabel] = Form("RooFormulaVar::%s('@0',{%s})", (v+"_"+objLabel).c_str(), ("AlphaR_"+objLabel).c_str()); }
	else if (v.rfind("Alpha" , 0)==0) { info.Par[v+"_"+objLabel] = Form("%s[%.6f,%.6f,%.6f]", (v+"_"+objLabel).c_str(),  2.000,    0.500,   8.000); }
	else if (v=="nR2"      && !found) { info.Par[v+"_"+objLabel] = Form("RooFormulaVar::%s('@0',{%s})", (v+"_"+objLabel).c_str(), ("nR_"+objLabel).c_str()); }
	else if (v.rfind("nR"    , 0)==0) { info.Par[v+"_"+objLabel] = Form("%s[%.6f,%.6f,%.6f]", (v+"_"+objLabel).c_str(),  6.000,    0.500,  25.000); }
	else if (v=="n2"       && !found) { info.Par[v+"_"+objLabel] = Form("RooFormulaVar::%s('@0',{%s})", (v+"_"+objLabel).c_str(), ("n_"+objLabel).c_str()); }
	else if (v.rfind("n"     , 0)==0) { info.Par[v+"_"+objLabel] = Form("%s[%.6f,%.6f,%.6f]", (v+"_"+objLabel).c_str(),  2.000,    0.500,  10.000); }
	else if (v.rfind("f"     , 0)==0) { info.Par[v+"_"+objLabel] = Form("%s[%.6f,%.6f,%.6f]", (v+"_"+objLabel).c_str(),  0.500,    0.000,   1.000); }
	else if (v.rfind("b"     , 0)==0) { info.Par[v+"_"+objLabel] = Form("%s[%.6f,%.6f,%.6f]", (v+"_"+objLabel).c_str(),  0.500,    0.000,   1.000); }
	else if (v=="Lambda"            ) { info.Par[v+"_"+objLabel] = Form("%s[%.6f,%.6f,%.6f]", (v+"_"+objLabel).c_str(),  0.200,   -2.000,   2.000); }
	else if (v=="LambdaF"           ) { info.Par[v+"_"+objLabel] = Form("%s[%.6f,%.6f,%.6f]", (v+"_"+objLabel).c_str(),  0.300,  0.00001,   5.000); }
	else if (v=="LambdaDS"          ) { info.Par[v+"_"+objLabel] = Form("%s[%.6f,%.6f,%.6f]", (v+"_"+objLabel).c_str(),  0.060,   0.0001,   5.000); }
	else if (v.rfind("LambdaS",0)==0) { info.Par[v+"_"+objLabel] = Form("%s[%.6f,%.6f,%.6f]", (v+"_"+objLabel).c_str(),  0.450,   0.0001,   5.000); }
	else if (v.rfind("Lambda", 0)==0) { info.Par[v+"_"+objLabel] = Form("%s[%.6f,%.6f,%.6f]", (v+"_"+objLabel).c_str(),  0.000, -100.000, 100.000); }
	else if (v=="Sigma2"   && !found) { info.Par[v+"_"+objLabel] = Form("RooFormulaVar::%s('@0*@1',{%s,%s})", ("Sigma2_"+objLabel).c_str(), ("rSigma21_"+objLabel).c_str(), ("Sigma1_"+objLabel).c_str()); }
	else if (v=="Sigma3"   && !found) { info.Par[v+"_"+objLabel] = Form("RooFormulaVar::%s('@0*@1',{%s,%s})", ("Sigma3_"+objLabel).c_str(), ("rSigma32_"+objLabel).c_str(), ("Sigma2_"+objLabel).c_str()); }
	else if (v=="Sigma4"   && !found) { info.Par[v+"_"+objLabel] = Form("RooFormulaVar::%s('@0*@1',{%s,%s})", ("Sigma4_"+objLabel).c_str(), ("rSigma43_"+objLabel).c_str(), ("Sigma3_"+objLabel).c_str()); }
	else if (v=="Sigma"             ) { info.Par[v+"_"+objLabel] = Form("%s[%.6f,%.6f,%.6f]", (v+"_"+objLabel).c_str(),  1.000,    0.010,  40.000); }
	else if (v=="xb"                ) { info.Par[v+"_"+objLabel] = Form("%s[%.6f,%.6f,%.6f]", (v+"_"+objLabel).c_str(), 10.000,    0.000,  40.000); }
	else if (v.rfind("Sigma",  0)==0) {
	  const auto& varR = ((fitVar=="Cand_Mass" && contain(MASS, obj)) ? MASS.at(obj).at("Width") : ((info.Var.at(fitVar).at("Max")-info.Var.at(fitVar).at("Min"))/6.0));
	  info.Par[v+"_"+objLabel] = Form("%s[%.10f,%.10f,%.10f]", (v+"_"+objLabel).c_str(), varR, 0.2*varR, 2.0*varR);
	}
	else if (v=="m") {
	  const auto& varV = ((fitVar=="Cand_Mass" && contain(MASS, obj)) ? MASS.at(obj).at("Val")   : ((info.Var.at(fitVar).at("Max")+info.Var.at(fitVar).at("Min"))/2.0));
	  const auto& varR = ((fitVar=="Cand_Mass" && contain(MASS, obj)) ? MASS.at(obj).at("Width") : ((info.Var.at(fitVar).at("Max")-info.Var.at(fitVar).at("Min"))/6.0));
	  if (fitVar=="Cand_Mass" && !contain(MASS, obj) && (obj!="Bkg" || obj.find("Swap")!=std::string::npos)) {
	    std::cout << "[WARNING] Initial value for " << (v+"_"+objLabel) << " was not found!" << std::endl;
	  }
	  info.Par[v+"_"+objLabel] = Form("%s[%.10f,%.10f,%.10f]", (v+"_"+objLabel).c_str(), varV, (varV - 2.0*varR), (varV + 2.0*varR));
	}
      }
      else {
	std::string content = info.Par.at(v+"_"+objFoundLabel);
	if (content.find("[")!=std::string::npos) {
	  content = content.substr(content.find("["));
	  info.Par[v+"_"+objLabel] = (v+"_"+objLabel+content);
	}
	else { info.Par[v+"_"+objLabel] = content; }
      }
    }
    // Check parameters for constrained fits
    if (v!="N" && v!="Cut") {
      for (const auto& con : StringVector_t({"val", "sig"})) {
	const auto& name = (con+v+"_"+objFoundLabel);
	if (contain(info.Par, name) && info.Par.at(name)!="") {
	  auto content = info.Par.at(name); content = content.substr( content.find("[") );
	  info.Par[con+v+"_"+objLabel] = (con+v+"_"+objLabel+content);
	}
	else break;
      }
    }
  }
  //
  return true;
};


bool addModelPar(RooWorkspace& ws, GlobalInfo& info, const StringVector_t& parNames, const std::string& fitVar,
		 const std::string& label, const std::string& modelN)
{
  // Check if parameters are already added
  bool addPar = false; for (const auto& v : parNames) { addPar = addPar || !ws.arg((v+"_"+label).c_str()); };
  if (!addPar) { return true; }
  // initialize all input parameters
  if (!setModelPar(info, parNames, label, fitVar, modelN)) { return false; }
  // Constrain Quarkonium excited states to 1S state
  if (fitVar=="Cand_Mass") constrainQuarkoniumMassParameters(info, parNames, label);
  // check that all input parameters are defined
  for (const auto& v : parNames) {
    if (!contain(info.Par, v+"_"+label)) {
      std::cout << "[ERROR] Initial parameter " << v << " was not found for " << modelN << " of " << label << std::endl; return false;
    }
  }
  // create the variables for this model
  RooArgList pdfConstrains;
  for (const auto& v : parNames) {
    if (ws.arg((v+"_"+label).c_str()) || v=="Cut") continue;
    if (contain(info.Par, v+"_"+label)) {
      if (!ws.factory(info.Par.at(v+"_"+label).c_str())) { std::cout << "[ERROR] Failed to create variable " << v+"_"+label << " defined as " << info.Par.at(v+"_"+label) << std::endl; return false; }
    }
    // create the Gaussian PDFs for Constrain fits
    if (contain(info.Par, "val"+v+"_"+label) && contain(info.Par, "sig"+v+"_"+label)) {
      if (!ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", (v+"_"+label).c_str(), (v+"_"+label).c_str(),
			   info.Par.at("val"+v+"_"+label).c_str(),
			   info.Par.at("sig"+v+"_"+label).c_str()
			   ))) { std::cout << "[ERROR] Failed to create PDF Constr" << v+"_"+label << std::endl; return false; }
      pdfConstrains.add(*ws.pdf(("Constr"+v+"_"+label).c_str()));
    }
  }
  if (pdfConstrains.getSize()>0) {
    const auto& pdfConsName = ("pdfConstr"+label.substr(label.find("To")));
    if (ws.genobj(pdfConsName.c_str())) { dynamic_cast<RooArgList*>(ws.genobj(pdfConsName.c_str()))->add(pdfConstrains); }
    else if (ws.import(pdfConstrains, pdfConsName.c_str())) { std::cout << "[ERROR] Failed to create constrain PDF list" << pdfConsName << std::endl; return false; }
  }
  return true;
};


TH1* rebinhist(const TH1& hist, const double& xmin, const double& xmax, const std::string& type="ReBin")
{
  auto hcopy = std::unique_ptr<TH1>(dynamic_cast<TH1*>(hist.Clone("hcopy")));
  // range of the new hist
  int imin = hcopy->FindBin(xmin);
  if (imin>=hcopy->GetNbinsX()) imin = 1;
  int imax = hcopy->FindBin(0.999999*xmax);
  if (imax<=1) imax = hcopy->GetNbinsX();
  std::vector<double> newbins;
  newbins.push_back(hcopy->GetBinLowEdge(imin));
  for (int i=imin; i<=imax; i++) {
    if (hcopy->GetBinContent(i)>0.0) {
      newbins.push_back(hcopy->GetBinLowEdge(i)+hcopy->GetBinWidth(i));
    } else {
      int nrebin=2;
      for (i++; i<=imax; i++) {
        if (hcopy->GetBinContent(i)>0.0) {
          newbins.push_back(hcopy->GetBinLowEdge(i)+hcopy->GetBinWidth(i));
          const double& newval = (hcopy->GetBinContent(i)/nrebin);
          hcopy->SetBinContent(i, newval);
          if (type=="FixBin") { for (int j=1; j<nrebin; j++) { hcopy->SetBinContent(i-j, newval); } }
          break;
        }
        nrebin++;
      }
    }
  }
  if (type=="ReBin") {
    if (xmin < newbins[1]) newbins[0] = xmin;
    if (xmax > newbins[newbins.size()-1]) { newbins.push_back(xmax); }
    return hcopy->Rebin(newbins.size()-1, "hnew", newbins.data());
  }
  else if (type=="FixBin") {
    return dynamic_cast<TH1*>(hcopy->Clone("hnew"));
  }
  return NULL;
};

  
bool histToPdf(RooWorkspace& ws, const std::string& pdfName, const RooDataSet& ds, const std::string& var, const std::vector< double >& range)
{
  //
  if (ws.pdf(pdfName.c_str())) { std::cout << "[INFO] The " << pdfName << " Template has already been created!" << std::endl; return true; }
  if (!ws.var(var.c_str())) { std::cout << "[ERROR] Variable " << var << " was not found!" << std::endl; return false; }
  std::cout << "[INFO] Implementing " << pdfName << " Template" << std::endl;
  //  
  // Create the histogram
  auto histName = pdfName;
  histName.replace(histName.find("pdf"), std::string("pdf").length(), "h");
  auto hist = std::unique_ptr<TH1D>(static_cast<TH1D*>(ds.createHistogram(histName.c_str(), *ws.var(var.c_str()), RooFit::Binning(int(range[0]), range[1], range[2]))));
  if (!hist) { std::cout << "[ERROR] Histogram " << histName << " is NULL!" << std::endl; return false; }
  // Cleaning the input histogram
  // 1) Remove the Under and Overflow bins
  hist->ClearUnderflowAndOverflow();
  // 2) Set negative bin content to zero
  for (int i=0; i<=hist->GetNbinsX(); i++) { if (hist->GetBinContent(i)<0.0) { hist->SetBinContent(i, 0.0); } }
  // 2) Reduce the range of histogram and rebin it
  hist.reset(static_cast<TH1D*>(rebinhist(*hist, range[1], range[2])));
  if (!hist) { std::cout << "[ERROR] Cleaned Histogram of " << histName << " is NULL!" << std::endl; return false; }
  auto dataName = pdfName;
  dataName.replace(dataName.find("pdf"), std::string("pdf").length(), "dh");
  std::unique_ptr<RooDataHist> dataHist = std::unique_ptr<RooDataHist>(new RooDataHist(dataName.c_str(), "", *ws.var(var.c_str()), hist.get()));
  if (!dataHist) { std::cout << "[ERROR] DataHist used to create " << pdfName << " failed!" << std::endl; return false; }
  if (dataHist->sumEntries()==0) { std::cout << "[WARNING] DataHist used to create " << pdfName << " is empty!" << std::endl; return false; } 
  if (std::abs(dataHist->sumEntries() - hist->GetSumOfWeights())>0.001) {
    std::cout << "[ERROR] DataHist (" << dataHist->sumEntries() << ")  used to create histogram (" << hist->GetSumOfWeights() << ") for PDF " << pdfName << "  " << " is invalid!  " << std::endl; return false;
  }
  const auto nBinTmp = ws.var(var.c_str())->getBins(); // Bug Fix
  if (ws.import(*dataHist)) { std::cout << "[ERROR] RooDataHist " << dataName << " was not imported!" << std::endl; return false; }
  ws.var(var.c_str())->setBins(nBinTmp); // Bug Fix
  auto pdf = std::unique_ptr<RooHistPdf>(new RooHistPdf(pdfName.c_str(), pdfName.c_str(), *ws.var(var.c_str()), *dynamic_cast<RooDataHist*>(ws.data(dataName.c_str()))));
  //std::unique_ptr<RooKeysPdf> pdf = std::unique_ptr<RooKeysPdf>(new RooKeysPdf(pdfName.c_str(), pdfName.c_str(), *ws.var(var.c_str()), *dynamic_cast<RooDataSet*>(ws.data(dsName.c_str())), RooKeysPdf::NoMirror, 0.4));
  if (!pdf) { std::cout << "[WARNING] RooHistPDF " << pdfName << " is NULL!" << std::endl; return false; }
  if (ws.import(*pdf)) { std::cout << "[ERROR] PDF " << pdfName << " was not imported!" << std::endl; return false; }
  // Return
  return true;
};


bool histToPdf(RooWorkspace& ws, const std::string& pdfName, const std::string& dsName, const std::string& var, const std::vector< double >& range)
{
  const auto& pd = dynamic_cast<RooDataSet*>(ws.data(dsName.c_str()));
  if (!pd) { std::cout << "[ERROR] DataSet " << dsName << " was not found!" << std::endl; return false; }
  if (pd->numEntries()<=2.0) { std::cout << "[WARNING] DataSet " << dsName << " has too few events!" << std::endl; return false; }
  return histToPdf(ws, pdfName, *pd, var, range);
};


#endif // #ifndef rooModelUtils_h
