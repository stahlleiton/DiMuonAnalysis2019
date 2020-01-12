#ifndef Candidate_addModel_C
#define Candidate_addModel_C


#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooRealVar.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooStringVar.h"

#include <iostream>
#include <string>
#include <memory>

#include "../Utilities/rooModelUtils.h"
#include "../Utilities/rooDataUtils.h"


RooAbsPdf*  getTotalPDF ( RooWorkspace& ws , const std::string& var , const std::string& label , const GlobalInfo& info   );
bool        makeSPlotDS ( RooWorkspace& ws , GlobalInfo& info       , const std::string& var   , const std::string& label );


bool addModel(RooWorkspace& ws, GlobalInfo& info, const std::string& chg, const StringSet_t& varV={})
{
  //
  const auto& cha = info.Par.at("channel");
  auto varS = info.StrS.at("fitVariable"); if (!varV.empty()) { varS = varV; }
  auto varT = info.StrS.at("fitVarName"); if (!varV.empty()) { varT.clear(); for (const auto& v : varV) { auto t = v; stringReplace(t, "_", ""); varT.insert(t); } }
  std::string varTot = ""; for (const auto& v : varT) { varTot += v; }
  for (const auto& col : info.StrS.at("fitSystem")) {
    const auto& lbl = cha + chg + "_" + col;
    //
    RooArgList pdfListTot;
    for (const auto& mainObj : info.StrS.at("fitObject")) {
      const auto& mainTag = mainObj + cha + chg;
      const auto& mainLabel = mainTag + "_" + col;
      //
      std::map<std::string, RooArgList> pdfListVar;
      for (const auto& varName : varS) {
	const std::string& varWindow   = "FitWindow";
	const auto& varNormName = (varName+"Norm");
	auto varType = varName; stringReplace(varType, "_", "");
	std::cout << "[INFO] Implementing " << mainTag << " " << varName << " Model for " << col << std::endl;
	// Make sure that DLenRes and main object are done first
	StringVector_t objV;
	if (contain(info.StrS.at("addObjectModel_"+varType+"_"+mainLabel), "DLenRes")) { objV.push_back("DLenRes"); }
	if (contain(info.StrS.at("addObjectModel_"+varType+"_"+mainLabel), mainObj)) { objV.push_back(mainObj); }
	for (const auto& m : info.StrS.at("addObjectModel_"+varType+"_"+mainLabel)) { if (m!="DLenRes" && m!=mainObj) { objV.push_back(m); } }
	std::map<std::string, RooArgList> pdfList;
	for (const auto& obj : objV) {
	  const auto& modelN = info.Par.at("Model"+varType+"_"+mainLabel+"_"+obj);
	  const auto& tag = obj + cha + chg;
	  const auto& label = tag + "_" + col;
	  auto objI = obj;
	  if (objI!="DLenRes" && !contain(info.StrV.at("variable"), objI)) {
	    for (const auto& c : info.StrV.at("category")) { if (objI.rfind(c)!=std::string::npos) { stringReplace(objI, c, ""); break; } }
	  }
	  addString(ws, ("Model"+varType+"_"+label), modelN); // Save the model name for bookkeeping
	  //
	  const std::string& pdfName    = Form("pdf%s_%s",    varType.c_str(), label.c_str());
	  const std::string& pdf1Name   = Form("pdf%s1_%s",   varType.c_str(), label.c_str());
	  const std::string& pdf2Name   = Form("pdf%s2_%s",   varType.c_str(), label.c_str());
	  const std::string& pdf3Name   = Form("pdf%s3_%s",   varType.c_str(), label.c_str());
	  const std::string& pdf4Name   = Form("pdf%s4_%s",   varType.c_str(), label.c_str());
	  const std::string& pdfPolName = Form("pdf%sPol_%s", varType.c_str(), label.c_str());
	  const std::string& pdfResName = Form("pdf%s_DLenRes%s", varType.c_str(), lbl.c_str());
	  //
	  // Create Models
	  switch(ModelDictionary.at(modelN))
	    {
	      //-------------------------------------------
	      //
	      // General Models
	      //
	      //-------------------------------------------
	    case (int(Model::CutAndCount)):
	      {
		// input variables
		const StringVector_t parNames = {"Cut"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		addString(ws, "CutAndCount_"+label, info.Par.at("Cut_"+label));
		info.Par["pdfName"+chg] = "CutAndCount_"+label;
		std::cout << "[INFO] " << tag << " in " << col << " added for CutAndCount!" << std::endl; break;
	      }
	    case (int(Model::Template)):
	      {
		// create the Template
		std::string dsName = ( "d" + chg + "_MC_" + obj + "_" + info.Par.at("channelDS") + "_" + col );
		const std::vector< double > range = { double(ws.var(varName.c_str())->getBins(varWindow.c_str())) , ws.var(varName.c_str())->getMin(varWindow.c_str()) , ws.var(varName.c_str())->getMax(varWindow.c_str()) };
		const auto& proceed = histToPdf(ws, pdfName, dsName, varName, range);
		if (!contain(info.Var, "recoMCEntries") || !contain(info.Var.at("recoMCEntries"), label)) {
		  if (proceed && (std::abs(info.Var.at("recoMCEntries").at(label)-ws.data(Form("dh%s_%s", varName.c_str(), label.c_str()))->sumEntries())>0.5)) {
		    std::cout << "[WARNING] The number of events in " << Form("dh%s_%s", varName.c_str(), label.c_str()) << " changed from (" << info.Var.at("recoMCEntries").at(label) << ") to (" <<
		      ws.data(Form("dh%s_%s", varName.c_str(), label.c_str()))->sumEntries() << ")" << std::endl;
		  }
		}
		//
		if (proceed) {
		  ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		  pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		  std::cout << "[INFO] " << tag << " Template " << varName << " PDF in " << col << " added!" << std::endl; break;
		}
		else {
		  std::cout << "[INFO] %s Template " << varName << " PDF in " << col << " was NOT created!" << std::endl; break;
		}
	      }
	    case (int(Model::SPLOT)):
	      {
		// create the sPlot DataSet
		if (!makeSPlotDS(ws, info, varName, label)) { return false; }
		// create the PDF
		const auto& sPlotDS = info.Par.at("dsName"+chg)+"_SPLOT";
		const auto& sPlotN = ("N_"+label+"_sw");
		const auto& range = std::vector<double>({((double)getNBins(varName, info)), info.Var.at(varName).at("Min"), info.Var.at(varName).at("Max")});
		auto dataw = std::unique_ptr<RooDataSet>(new RooDataSet("TMP","TMP", dynamic_cast<RooDataSet*>(ws.data(sPlotDS.c_str())), RooArgSet(*ws.var(varName.c_str()), *ws.var(sPlotN.c_str())), 0, sPlotN.c_str()));
		if (!histToPdf(ws, pdfName, *dataw, varName, range)) { return false; }
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " sPlot " << varName << " Template in " << col << " added!" << std::endl; break;
	      }
	      //-------------------------------------------
	      //
	      // Signal Candidate Mass Models
	      //
	      //-------------------------------------------
	    case (int(Model::SingleGaussian)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the PDF
		if (!ws.factory(Form("Gaussian::%s(%s, %s, %s)", pdfName.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma1_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Single Gaussian " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::DoubleGaussian)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "f"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the two PDFs
		if (!ws.factory(Form("Gaussian::%s(%s, %s, %s)", pdf1Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma1_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf1Name << std::endl; return false; }
		if (!ws.factory(Form("Gaussian::%s(%s, %s, %s)", pdf2Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma2_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf2Name << std::endl; return false; }
		// Sum the PDFs to get the signal PDF
		if (!ws.factory(Form("SUM::%s(%s*%s, %s)", pdfName.c_str(),
				     ("f_"+label).c_str(),
				     pdf1Name.c_str(),
				     pdf2Name.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Double Gaussian " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::SingleCrystalBall)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "Alpha", "n"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the PDF
		if (!ws.factory(Form("CBShape::%s(%s, %s, %s, %s, %s)", pdfName.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma1_"+label).c_str(),
				     ("Alpha_"+label).c_str(),
				     ("n_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Single Crystal Ball " << varName << " PDF in " << col << " included" << std::endl; break;
	      }
	    case (int(Model::DoubleCrystalBall)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "Alpha", "Alpha2", "n", "n2", "f"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the two PDFs
		if (!ws.factory(Form("CBShape::%s(%s, %s, %s, %s, %s)", pdf1Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma1_"+label).c_str(),
				     ("Alpha_"+label).c_str(),
				     ("n_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf1Name << std::endl; return false; }
		if (!ws.factory(Form("CBShape::%s(%s, %s, %s, %s, %s)", pdf2Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma2_"+label).c_str(),
				     ("Alpha2_"+label).c_str(),
				     ("n2_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf2Name << std::endl; return false; }
		// Sum the PDFs to get the signal PDF
		if (!ws.factory(Form("SUM::%s(%s*%s, %s)", pdfName.c_str(),
				     ("f_"+label).c_str(),
				     pdf1Name.c_str(),
				     pdf2Name.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Double Crystal Ball " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::GaussianAndCrystalBall)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "Alpha", "n", "f"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the two PDFs
		if (!ws.factory(Form("CBShape::%s(%s, %s, %s, %s, %s)", pdf1Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma1_"+label).c_str(),
				     ("Alpha_"+label).c_str(),
				     ("n_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf1Name << std::endl; return false; }
		if (!ws.factory(Form("Gaussian::%s(%s, %s, %s)", pdf2Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma2_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf2Name << std::endl; return false; }
		// Sum the PDFs to get the signal PDF
		if (!ws.factory(Form("SUM::%s(%s*%s, %s)", pdfName.c_str(),
				     ("f_"+label).c_str(),
				     pdf1Name.c_str(),
				     pdf2Name.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Gaussian and Crystal Ball " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::SingleExtCrystalBall)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "Alpha", "n", "AlphaR", "nR"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// import the model class
		if (!importModelClass(ws, "RooExtCBShape")) { std::cout << "[ERROR] Could not import the class RooExtCBShape!" << std::endl; return false; }
		// create the PDF
		if (!ws.factory(Form("RooExtCBShape::%s(%s, %s, %s, %s, %s, %s, %s)", pdfName.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma1_"+label).c_str(),
				     ("Alpha_"+label).c_str(),
				     ("n_"+label).c_str(),
				     ("AlphaR_"+label).c_str(),
				     ("nR_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Single Extended Crystal Ball " << varName << " PDF in " << col << " included" << std::endl; break;
	      }
	    case (int(Model::DoubleExtCrystalBall)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "Alpha", "n", "AlphaR", "nR", "Alpha2", "n2", "AlphaR2", "nR2", "f"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// import the model class
		if (!importModelClass(ws, "RooExtCBShape")) { std::cout << "[ERROR] Could not import the class RooExtCBShape!" << std::endl; return false; }
		// create the two PDFs
		if (!ws.factory(Form("RooExtCBShape::%s(%s, %s, %s, %s, %s, %s, %s)", pdf1Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma1_"+label).c_str(),
				     ("Alpha_"+label).c_str(),
				     ("n_"+label).c_str(),
				     ("AlphaR_"+label).c_str(),
				     ("nR_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf1Name << std::endl; return false; }
		if (!ws.factory(Form("RooExtCBShape::%s(%s, %s, %s, %s, %s, %s, %s)", pdf2Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma2_"+label).c_str(),
				     ("Alpha2_"+label).c_str(),
				     ("n2_"+label).c_str(),
				     ("AlphaR2_"+label).c_str(),
				     ("nR2_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf2Name << std::endl; return false; }
		// sum the PDFs to get the signal PDF
		if (!ws.factory(Form("SUM::%s(%s*%s, %s)", pdfName.c_str(),
				     ("f_"+label).c_str(),
				     pdf1Name.c_str(),
				     pdf2Name.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Double Extended Crystal Ball " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::GaussianAndExtCrystalBall)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "Alpha", "n", "AlphaR", "nR", "f"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// import the model class
		if (!importModelClass(ws, "RooExtCBShape")) { std::cout << "[ERROR] Could not import the class RooExtCBShape!" << std::endl; return false; }
		// create the two PDFs
		if (!ws.factory(Form("RooExtCBShape::%s(%s, %s, %s, %s, %s, %s, %s)", pdf1Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma1_"+label).c_str(),
				     ("Alpha_"+label).c_str(),
				     ("n_"+label).c_str(),
				     ("AlphaR_"+label).c_str(),
				     ("nR_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf1Name << std::endl; return false; }
		if (!ws.factory(Form("Gaussian::%s(%s, %s, %s)", pdf2Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma2_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf2Name << std::endl; return false; }
		// sum the PDFs to get the signal PDF
		if (!ws.factory(Form("SUM::%s(%s*%s, %s)", pdfName.c_str(),
				     ("f_"+label).c_str(),
				     pdf1Name.c_str(),
				     pdf2Name.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Gaussian and Extended Crystal Ball " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::SingleModCrystalBall)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "Alpha"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// import the model class
		if (!importModelClass(ws, "RooModCBShape")) { std::cout << "[ERROR] Could not import the class RooModCBShape!" << std::endl; return false; }
		// create the PDF
		if (!ws.factory(Form("RooModCBShape::%s(%s,%s,%s,%s)", pdfName.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma1_"+label).c_str(),
				     ("Alpha_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Single Modified Crystal Ball " << varName << " PDF in " << col << " included" << std::endl; break;
	      }
	    case (int(Model::DoubleModCrystalBall)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "Alpha", "Alpha2", "f"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// import the model class
		if (!importModelClass(ws, "RooModCBShape")) { std::cout << "[ERROR] Could not import the class RooModCBShape!" << std::endl; return false; }
		// create the two PDFs
		if (!ws.factory(Form("RooModCBShape::%s(%s,%s,%s,%s)", pdf1Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma1_"+label).c_str(),
				     ("Alpha_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf1Name << std::endl; return false; }
		if (!ws.factory(Form("RooModCBShape::%s(%s,%s,%s,%s)", pdf2Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma2_"+label).c_str(),
				     ("Alpha2_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf2Name << std::endl; return false; }
		// sum the PDFs to get the signal PDF
		if (!ws.factory(Form("SUM::%s(%s*%s, %s)", pdfName.c_str(),
				     ("f_"+label).c_str(),
				     pdf1Name.c_str(),
				     pdf2Name.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Double Modified Crystal Ball " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::GaussianAndModCrystalBall)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "Alpha", "f"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// import the model class
		if (!importModelClass(ws, "RooModCBShape")) { std::cout << "[ERROR] Could not import the class RooModCBShape!" << std::endl; return false; }
		// create the two PDFs
		if (!ws.factory(Form("RooModCBShape::%s(%s,%s,%s,%s)", pdf1Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma1_"+label).c_str(),
				     ("Alpha_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf1Name << std::endl; return false; }
		if (!ws.factory(Form("Gaussian::%s(%s, %s, %s)", pdf2Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma2_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf2Name << std::endl; return false; }
		// sum the PDFs to get the signal PDF
		if (!ws.factory(Form("SUM::%s(%s*%s, %s)", pdfName.c_str(),
				     ("f_"+label).c_str(),
				     pdf1Name.c_str(),
				     pdf2Name.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Gaussian and Modified Crystal Ball " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::SingleModExtCrystalBall)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "Alpha", "n", "AlphaR"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// import the model class
		if (!importModelClass(ws, "RooModExtCBShape")) { std::cout << "[ERROR] Could not import the class RooModExtCBShape!" << std::endl; return false; }
		// create the PDF
		if (!ws.factory(Form("RooModExtCBShape::%s(%s, %s, %s, %s, %s, %s)", pdfName.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma1_"+label).c_str(),
				     ("Alpha_"+label).c_str(),
				     ("n_"+label).c_str(),
				     ("AlphaR_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Single Modified Extended Crystal Ball " << varName << " PDF in " << col << " included" << std::endl; break;
	      }
	    case (int(Model::DoubleModExtCrystalBall)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "Alpha", "n", "AlphaR", "Alpha2", "n2", "AlphaR2", "f"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// import the model class
		if (!importModelClass(ws, "RooModExtCBShape")) { std::cout << "[ERROR] Could not import the class RooModExtCBShape!" << std::endl; return false; }
		// create the two PDFs
		if (!ws.factory(Form("RooModExtCBShape::%s(%s, %s, %s, %s, %s, %s)", pdf1Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma1_"+label).c_str(),
				     ("Alpha_"+label).c_str(),
				     ("n_"+label).c_str(),
				     ("AlphaR_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf1Name << std::endl; return false; }
		if (!ws.factory(Form("RooModExtCBShape::%s(%s, %s, %s, %s, %s, %s)", pdf2Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma2_"+label).c_str(),
				     ("Alpha2_"+label).c_str(),
				     ("n2_"+label).c_str(),
				     ("AlphaR2_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf2Name << std::endl; return false; }
		// sum the PDFs to get the signal PDF
		if (!ws.factory(Form("SUM::%s(%s*%s, %s)", pdfName.c_str(),
				     ("f_"+label).c_str(),
				     pdf1Name.c_str(),
				     pdf2Name.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Double Modified Extended Crystal Ball " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::GaussianAndModExtCrystalBall)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "Alpha", "n", "AlphaR", "f"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// import the model class
		if (!importModelClass(ws, "RooModExtCBShape")) { std::cout << "[ERROR] Could not import the class RooModExtCBShape!" << std::endl; return false; }
		// create the two PDFs
		if (!ws.factory(Form("RooModExtCBShape::%s(%s, %s, %s, %s, %s, %s)", pdf1Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma1_"+label).c_str(),
				     ("Alpha_"+label).c_str(),
				     ("n_"+label).c_str(),
				     ("AlphaR_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf1Name << std::endl; return false; }
		if (!ws.factory(Form("Gaussian::%s(%s, %s, %s)", pdf2Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma2_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf2Name << std::endl; return false; }
		// sum the PDFs to get the signal PDF
		if (!ws.factory(Form("SUM::%s(%s*%s, %s)", pdfName.c_str(),
				     ("f_"+label).c_str(),
				     pdf1Name.c_str(),
				     pdf2Name.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Gaussian and Modified Extended Crystal Ball " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::Voigtian)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the PDF
		if (!ws.factory(Form("Voigtian::%s(%s, %s, %s, %s)", pdfName.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma1_"+label).c_str(),
				     ("Sigma2_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Voigtian " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::BWCrystalBall)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "Alpha", "n"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the two PDFs
		if (!ws.factory(Form("BreitWigner::%s(%s, %s, %s)", pdf1Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma1_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf1Name << std::endl; return false; }
		if (!ws.factory(Form("CBShape::%s(%s, %s, %s, %s, %s)", pdf2Name.c_str(), varName.c_str(),
				     ("Peak_"+label+"[0]").c_str(),
				     ("Sigma2_"+label).c_str(),
				     ("Alpha_"+label).c_str(),
				     ("n_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf2Name << std::endl; return false; }
		// covolve the PDFs to get the signal PDF
		if (!ws.factory(Form("FCONV::%s(%s, %s, %s)", pdfName.c_str(), varName.c_str(),
				     pdf1Name.c_str(),
				     pdf2Name.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Breit-Wigner-CrystallBall " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	      //-------------------------------------------
	      //
	      // Background Candidate Mass Models
	      //
	      //-------------------------------------------
	    case (int(Model::Uniform)):
	      {
		// create the PDF
		if (!ws.factory(Form("Uniform::%s(%s)", pdfName.c_str(), varName.c_str()))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Uniform " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::Chebychev1)):
	      {
		// input variables
		const StringVector_t parNames = {"Lambda1"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the PDF
		if (!ws.factory(Form("Chebychev::%s(%s, {%s})", pdfName.c_str(), varName.c_str(),
				     ("Lambda1_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " 1st Order Chebychev " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::Chebychev2)):
	      {
		// input variables
		const StringVector_t parNames = {"Lambda1", "Lambda2"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the PDF
		if (!ws.factory(Form("Chebychev::%s(%s, {%s, %s})", pdfName.c_str(), varName.c_str(),
				     ("Lambda1_"+label).c_str(),
				     ("Lambda2_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " 2nd Order Chebychev " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::Chebychev3)):
	      {
		// input variables
		const StringVector_t parNames = {"Lambda1", "Lambda2", "Lambda3"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the PDF
		if (!ws.factory(Form("Chebychev::%s(%s, {%s, %s, %s})", pdfName.c_str(), varName.c_str(),
				     ("Lambda1_"+label).c_str(),
				     ("Lambda2_"+label).c_str(),
				     ("Lambda3_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " 3rd Order Chebychev " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::Chebychev4)):
	      {
		// input variables
		const StringVector_t parNames = {"Lambda1", "Lambda2", "Lambda3", "Lambda4"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the PDF
		if (!ws.factory(Form("Chebychev::%s(%s, {%s, %s, %s, %s})", pdfName.c_str(), varName.c_str(),
				     ("Lambda1_"+label).c_str(),
				     ("Lambda2_"+label).c_str(),
				     ("Lambda3_"+label).c_str(),
				     ("Lambda4_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " 4th Order Chebychev " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::Chebychev5)):
	      {
		// input variables
		const StringVector_t parNames = {"Lambda1", "Lambda2", "Lambda3", "Lambda4", "Lambda5"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the PDF
		if (!ws.factory(Form("Chebychev::%s(%s, {%s, %s, %s, %s, %s})", pdfName.c_str(), varName.c_str(),
				     ("Lambda1_"+label).c_str(),
				     ("Lambda2_"+label).c_str(),
				     ("Lambda3_"+label).c_str(),
				     ("Lambda4_"+label).c_str(),
				     ("Lambda5_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " 5th Order Chebychev " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::Chebychev6)):
	      {
		// input variables
		const StringVector_t parNames = {"Lambda1", "Lambda2", "Lambda3", "Lambda4", "Lambda5", "Lambda6"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the PDF
		if (!ws.factory(Form("Chebychev::%s(%s, {%s, %s, %s, %s, %s, %s})", pdfName.c_str(), varName.c_str(),
				     ("Lambda1_"+label).c_str(),
				     ("Lambda2_"+label).c_str(),
				     ("Lambda3_"+label).c_str(),
				     ("Lambda4_"+label).c_str(),
				     ("Lambda5_"+label).c_str(),
				     ("Lambda6_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " 6th Order Chebychev " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::ExpChebychev1)):
	      {
		// input variables
		const StringVector_t parNames = {"Lambda1"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the PDF
		if (!ws.factory(Form("RooFormulaVar::%s('( 2.0*@0 - @2 - @1 )/( @2 - @1 )', {%s, vMin[%.6f], vMax[%.6f]})",
				     varNormName.c_str(), varName.c_str(),
				     info.Var.at(varName).at("Min"),
				     info.Var.at(varName).at("Max")
				     ))) { std::cout << "[ERROR] Failed to create variable " << varNormName << std::endl; return false; }
		if (!ws.factory(Form("RooFormulaVar::%s('@1*(@0) + 1.0', {%s, %s})",
				     pdfPolName.c_str(), varNormName.c_str(),
				     ("Lambda1_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfPolName << std::endl; return false; }
		if (!ws.factory(Form("Exponential::%s(%s, One[1.0])", pdfName.c_str(),
				     pdfPolName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " 1st Order Exponential Chebychev " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::ExpChebychev2)):
	      {
		// input variables
		const StringVector_t parNames = {"Lambda1", "Lambda2"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the PDF
		if (!ws.factory(Form("RooFormulaVar::%s('( 2.0*@0 - @2 - @1 )/( @2 - @1 )', {%s, vMin[%.6f], vMax[%.6f]})",
				     varNormName.c_str(), varName.c_str(),
				     info.Var.at(varName).at("Min"),
				     info.Var.at(varName).at("Max")
				     ))) { std::cout << "[ERROR] Failed to create variable " << varNormName << std::endl; return false; }
		if (!ws.factory(Form("RooFormulaVar::%s('@2*(2.0*@0*@0 - 1.0) + @1*(@0) + 1.0', {%s, %s, %s})",
				     pdfPolName.c_str(), varNormName.c_str(),
				     ("Lambda1_"+label).c_str(),
				     ("Lambda2_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfPolName << std::endl; return false; }
		if (!ws.factory(Form("Exponential::%s(%s, One[1.0])", pdfName.c_str(),
				     pdfPolName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " 2nd Order Exponential Chebychev " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::ExpChebychev3)):
	      {
		// input variables
		const StringVector_t parNames = {"Lambda1", "Lambda2", "Lambda3"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the PDF
		if (!ws.factory(Form("RooFormulaVar::%s('( 2.0*@0 - @2 - @1 )/( @2 - @1 )', {%s, vMin[%.6f], vMax[%.6f]})",
				     varNormName.c_str(), varName.c_str(),
				     info.Var.at(varName).at("Min"),
				     info.Var.at(varName).at("Max")
				     ))) { std::cout << "[ERROR] Failed to create variable " << varNormName << std::endl; return false; }
		if (!ws.factory(Form("RooFormulaVar::%s('@3*(4.0*@0*@0*@0 - 3.0*@0) + @2*(2.0*@0*@0 - 1.0) + @1*(@0) + 1.0', {%s, %s, %s, %s})",
				     pdfPolName.c_str(), varNormName.c_str(),
				     ("Lambda1_"+label).c_str(),
				     ("Lambda2_"+label).c_str(),
				     ("Lambda3_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfPolName << std::endl; return false; }
		if (!ws.factory(Form("Exponential::%s(%s, One[1.0])", pdfName.c_str(),
				     pdfPolName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " 3rd Order Exponential Chebychev " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::ExpChebychev4)):
	      {
		// input variables
		const StringVector_t parNames = {"Lambda1", "Lambda2", "Lambda3", "Lambda4"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the PDF
		if (!ws.factory(Form("RooFormulaVar::%s('( 2.0*@0 - @2 - @1 )/( @2 - @1 )', {%s, vMin[%.6f], vMax[%.6f]})",
				     varNormName.c_str(), varName.c_str(),
				     info.Var.at(varName).at("Min"),
				     info.Var.at(varName).at("Max")
				     ))) { std::cout << "[ERROR] Failed to create variable " << varNormName << std::endl; return false; }
		if (!ws.factory(Form("RooFormulaVar::%s('@4*(8.0*@0*@0*@0*@0 - 8.0*@0*@0 + 1.0) + @3*(4.0*@0*@0*@0 - 3.0*@0) + @2*(2.0*@0*@0 - 1.0) + @1*(@0) + 1.0', {%s, %s, %s, %s, %s})",
				     pdfPolName.c_str(), varNormName.c_str(),
				     ("Lambda1_"+label).c_str(),
				     ("Lambda2_"+label).c_str(),
				     ("Lambda3_"+label).c_str(),
				     ("Lambda4_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfPolName << std::endl; return false; }
		if (!ws.factory(Form("Exponential::%s(%s, One[1.0])", pdfName.c_str(),
				     pdfPolName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " 4th Order Exponential Chebychev " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::ExpChebychev5)):
	      {
		// input variables
		const StringVector_t parNames = {"Lambda1", "Lambda2", "Lambda3", "Lambda4", "Lambda5"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the PDF
		if (!ws.factory(Form("RooFormulaVar::%s('( 2.0*@0 - @2 - @1 )/( @2 - @1 )', {%s, vMin[%.6f], vMax[%.6f]})",
				     varNormName.c_str(), varName.c_str(),
				     info.Var.at(varName).at("Min"),
				     info.Var.at(varName).at("Max")
				     ))) { std::cout << "[ERROR] Failed to create variable " << varNormName << std::endl; return false; }
		if (!ws.factory(Form("RooFormulaVar::%s('@5*(16.0*@0*@0*@0*@0*@0 - 20.0*@0*@0*@0 + 5.0*@0) + @4*(8.0*@0*@0*@0*@0 - 8.0*@0*@0 + 1.0) + @3*(4.0*@0*@0*@0 - 3.0*@0) + @2*(2.0*@0*@0 - 1.0) + @1*(@0) + 1.0', {%s, %s, %s, %s, %s, %s})",
				     pdfPolName.c_str(), varNormName.c_str(),
				     ("Lambda1_"+label).c_str(),
				     ("Lambda2_"+label).c_str(),
				     ("Lambda3_"+label).c_str(),
				     ("Lambda4_"+label).c_str(),
				     ("Lambda5_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfPolName << std::endl; return false; }
		if (!ws.factory(Form("Exponential::%s(%s, One[1.0])", pdfName.c_str(),
				     pdfPolName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " 5th Order Exponential Chebychev " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::ExpChebychev6)):
	      {
		// input variables
		const StringVector_t parNames = {"Lambda1", "Lambda2", "Lambda3", "Lambda4", "Lambda5", "Lambda6"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the PDF
		if (!ws.factory(Form("RooFormulaVar::%s('( 2.0*@0 - @2 - @1 )/( @2 - @1 )', {%s, vMin[%.6f], vMax[%.6f]})",
				     varNormName.c_str(), varName.c_str(),
				     info.Var.at(varName).at("Min"),
				     info.Var.at(varName).at("Max")
				     ))) { std::cout << "[ERROR] Failed to create variable " << varNormName << std::endl; return false; }
		if (!ws.factory(Form("RooFormulaVar::%s('@6*(32.0*@0*@0*@0*@0*@0*@0 - 48.0*@0*@0*@0*@0 + 18.0*@0*@0 - 1.0) + @5*(16.0*@0*@0*@0*@0*@0 - 20.0*@0*@0*@0 + 5.0*@0) + @4*(8.0*@0*@0*@0*@0 - 8.0*@0*@0 + 1.0) + @3*(4.0*@0*@0*@0 - 3.0*@0) + @2*(2.0*@0*@0 - 1.0) + @1*(@0) + 1.0', {%s, %s, %s, %s, %s, %s, %s})",
				     pdfPolName.c_str(), varNormName.c_str(),
				     ("Lambda1_"+label).c_str(),
				     ("Lambda2_"+label).c_str(),
				     ("Lambda3_"+label).c_str(),
				     ("Lambda4_"+label).c_str(),
				     ("Lambda5_"+label).c_str(),
				     ("Lambda6_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfPolName << std::endl; return false; }
		if (!ws.factory(Form("Exponential::%s(%s, One[1.0])", pdfName.c_str(),
				     pdfPolName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " 6th Order Exponential Chebychev " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::Exponential)):
	      {
		// input variables
		const StringVector_t parNames = {"Lambda1"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the PDF
		if (!ws.factory(Form("Exponential::%s(%s, %s)", pdfName.c_str(), varName.c_str(),
				     ("Lambda1_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Exponential " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::ExpError)):
	      {
		// input variables
		const StringVector_t parNames = {"Sigma", "xb", "Lambda"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the PDF
		if (!ws.factory(Form("RooGenericPdf::%s('TMath::Exp(-@0*@1)*(1.0+TMath::Erf((@0-@2)/@3))', {%s, %s, %s, %s})", pdfName.c_str(), varName.c_str(),
				     ("Lambda_"+label).c_str(),
				     ("xb_"+label).c_str(),
				     ("Sigma_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " ExpError " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	      //-------------------------------------------
	      //
	      // Candidate Decay Length Resolution Models
	      //
	      //-------------------------------------------
	    case (int(Model::DeltaResolution)):
	      {
		// create the PDF
		if (!ws.factory(Form("TruthModel::%s(%s)", pdfName.c_str(), varName.c_str()))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		if (varName=="Cand_DLenRes") pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Truth " << varName << " Model in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::SingleGaussianResolution)):
	      {
		// input variables
		StringVector_t parNames = {"m", "Sigma1"};
		// create the variables for this model
		if (!ws.var("One")) { ws.factory("One[1.0]"); }
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the PDF
		if (!ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", pdfName.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma1_"+label).c_str(),
				     (varName=="Cand_DLenRes" ? "One" : "Cand_DLenErr")
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		if (varName=="Cand_DLenRes") pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Single Gaussian " << varName << " Model in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::DoubleGaussianResolution)):
	      {
		// input variables
		StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "f"};
		// create the variables for this model
		if (!ws.var("One")) { ws.factory("One[1.0]"); }
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the two PDFs
		if (!ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", pdf1Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma1_"+label).c_str(),
				     (varName=="Cand_DLenRes" ? "One" : "Cand_DLenErr")
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf1Name << std::endl; return false; }
		if (!ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", pdf2Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma2_"+label).c_str(),
				     (varName=="Cand_DLenRes" ? "One" : "Cand_DLenErr")
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf2Name << std::endl; return false; }
		// sum the PDFs to get the total PDF
		if (!ws.factory(Form("AddModel::%s({%s, %s}, %s)", pdfName.c_str(),
				     pdf1Name.c_str(),
				     pdf2Name.c_str(),
				     ("f_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		if (varName=="Cand_DLenRes") pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Double Gaussian " << varName << " Model in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::TripleGaussianResolution)):
	      {
		// input variables
		StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "f", "rSigma32", "Sigma3", "f2"};
		// create the variables for this model
		if (!ws.var("One")) { ws.factory("One[1.0]"); }
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the two PDFs
		if (!ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", pdf1Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma1_"+label).c_str(),
				     (varName=="Cand_DLenRes" ? "One" : "Cand_DLenErr")
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf1Name << std::endl; return false; }
		if (!ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", pdf2Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma2_"+label).c_str(),
				     (varName=="Cand_DLenRes" ? "One" : "Cand_DLenErr")
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf2Name << std::endl; return false; }
		if (!ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", pdf3Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma3_"+label).c_str(),
				     (varName=="Cand_DLenRes" ? "One" : "Cand_DLenErr")
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf3Name << std::endl; return false; }
		// sum the PDFs to get the total PDF
		if (!ws.factory(Form("AddModel::%s({%s, %s, %s}, {%s, %s})", pdfName.c_str(),
				     pdf1Name.c_str(),
				     pdf2Name.c_str(),
				     pdf3Name.c_str(),
				     ("f_"+label).c_str(),
				     ("f2_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		if (varName=="Cand_DLenRes") pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Triple Gaussian " << varName << " Model in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::QuadrupleGaussianResolution)):
	      {
		// input variables
		StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "f", "rSigma32", "Sigma3", "f2", "rSigma43", "Sigma4", "f3"};
		// create the variables for this model
		if (!ws.var("One")) { ws.factory("One[1.0]"); }
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the two PDFs
		if (!ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", pdf1Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma1_"+label).c_str(),
				     (varName=="Cand_DLenRes" ? "One" : "Cand_DLenErr")
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf1Name << std::endl; return false; }
		if (!ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", pdf2Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma2_"+label).c_str(),
				     (varName=="Cand_DLenRes" ? "One" : "Cand_DLenErr")
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf2Name << std::endl; return false; }
		if (!ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", pdf3Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma3_"+label).c_str(),
				     (varName=="Cand_DLenRes" ? "One" : "Cand_DLenErr")
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf3Name << std::endl; return false; }
		if (!ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", pdf4Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma4_"+label).c_str(),
				     (varName=="Cand_DLenRes" ? "One" : "Cand_DLenErr")
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf4Name << std::endl; return false; }
		// sum the PDFs to get the total PDF
		if (!ws.factory(Form("AddModel::%s({%s, %s, %s, %s}, {%s, %s, %s})", pdfName.c_str(),
				     pdf1Name.c_str(),
				     pdf2Name.c_str(),
				     pdf3Name.c_str(),
				     pdf4Name.c_str(),
				     ("f_"+label).c_str(),
				     ("f2_"+label).c_str(),
				     ("f3_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		if (varName=="Cand_DLenRes") pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Quadruple Gaussian " << varName << " Model in " << col << " added!" << std::endl; break;
	      }
	      //-------------------------------------------
	      //
	      // Candidate Decay Length Models
	      //
	      //-------------------------------------------
	    case (int(Model::Delta)):
	      {
		// create the PDF
		if (!ws.factory(Form("SUM::%s(%s)", pdfName.c_str(),
				     pdfResName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Delta " << varName << " Model in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::SingleSidedDecay)):
	      {
		// input variables
		StringVector_t parNames = {"LambdaSS"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the PDF
		if (!ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", pdfName.c_str(), varName.c_str(),
				     ("LambdaSS_"+label).c_str(),
				     pdfResName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Single Sided Decay " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::TripleDecay)):
	      {
		// input variables
		StringVector_t parNames = {"LambdaSS", "LambdaF", "LambdaDS", "f", "f2"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the PDF
		if (!ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", pdf1Name.c_str(), varName.c_str(),
				     ("LambdaSS_"+label).c_str(),
				     pdfResName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf1Name << std::endl; return false; }
		if (!ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", pdf2Name.c_str(), varName.c_str(),
				     ("LambdaF_"+label).c_str(),
				     pdfResName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf2Name << std::endl; return false; }
		if (!ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::DoubleSided)", pdfName.c_str(), varName.c_str(),
				     ("LambdaDS_"+label).c_str(),
				     pdfResName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf3Name << std::endl; return false; }
		// sum the PDFs to get the total PDF
		if (!ws.factory(Form("SUM::%s({%s, %s, %s}, {%s, %s})", pdfName.c_str(),
				     pdf1Name.c_str(),
				     pdf2Name.c_str(),
				     pdf3Name.c_str(),
				     ("f_"+label).c_str(),
				     ("f2_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Triple Decay " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::QuadrupleDecay)):
	      {
		// input variables
		StringVector_t parNames = {"LambdaSS", "LambdaSS2", "LambdaF", "LambdaDS", "f", "f2", "f3"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the PDF
		if (!ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", pdf1Name.c_str(), varName.c_str(),
				     ("LambdaSS_"+label).c_str(),
				     pdfResName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf1Name << std::endl; return false; }
		if (!ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", pdf2Name.c_str(), varName.c_str(),
				     ("LambdaSS2_"+label).c_str(),
				     pdfResName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf2Name << std::endl; return false; }
		if (!ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", pdf3Name.c_str(), varName.c_str(),
				     ("LambdaF_"+label).c_str(),
				     pdfResName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf3Name << std::endl; return false; }
		if (!ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::DoubleSided)", pdf4Name.c_str(), varName.c_str(),
				     ("LambdaDS_"+label).c_str(),
				     pdfResName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf4Name << std::endl; return false; }
		// sum the PDFs to get the total PDF
		if (!ws.factory(Form("SUM::%s({%s, %s, %s, %s}, {%s, %s, %s})", pdfName.c_str(),
				     pdf1Name.c_str(),
				     pdf2Name.c_str(),
				     pdf3Name.c_str(),
				     pdf4Name.c_str(),
				     ("f_"+label).c_str(),
				     ("f2_"+label).c_str(),
				     ("f3_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Quadruple Decay " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    default :
	      {
		if (modelN=="") { std::cout << "[ERROR] Candidate Mass Model for " << modelN << " was not defined (is empty)!" << std::endl; return false; }
		else { std::cout << "[ERROR] Selected Candidate Mass Model: " << modelN << " has not been implemented" << std::endl; return false; }
	      }
	    }
	}
	for (const auto& p : pdfList) {
	  if (p.second.getSize()>1) {
	    const auto& pdfName = ("pdf"+varType+"_"+p.first+lbl);
	    // create the parameters for the PDF sum
	    StringVector_t parNames = {"b"}; for (int i=3; i<=p.second.getSize(); i++) { parNames.push_back(Form("b%d",i-1)); }
	    if (!addModelPar(ws, info, parNames, varName, (p.first+lbl), ("SUM_"+p.first))) { return false; }
	    RooArgList coefList; for (const auto& par : parNames) { if (ws.var(par.c_str())) { coefList.add(*ws.var(par.c_str())); } }
	    // Sum the PDFs
	    auto themodel = std::unique_ptr<RooAddPdf>(new RooAddPdf(pdfName.c_str(), pdfName.c_str(), p.second, coefList, true));
	    if (ws.import(*themodel)) { std::cout << "[ERROR] Failed to import PDF " << pdfName << std::endl; return false; }
	    pdfListVar[p.first].add(*ws.pdf(pdfName.c_str()));
	  }
	  else { pdfListVar[p.first].add(p.second); }
	}
      }
      RooArgList pdfListVarTot;
      for (const auto& p : pdfListVar) {
	// Multiply the PDFs
	std::string pdfName = p.second.at(0)->GetName();
	if (p.second.getSize()>1) {
	  pdfName = ("pdf"+varTot+"_"+p.first+lbl);
	  auto themodel = std::unique_ptr<RooProdPdf>(new RooProdPdf(pdfName.c_str(), pdfName.c_str(), p.second));
	  if (ws.import(*themodel)) { std::cout << "[ERROR] Failed to import PDF " << pdfName << std::endl; return false; }
	}
	// create the yield
	if (!addModelPar(ws, info, {"N"}, "Cand_Mass", (p.first+lbl), "")) { return false; }
	// extend the PDF
	const auto& pdfTotName = ("pdf"+varTot+"Tot_"+p.first+lbl);
	if (!ws.factory(Form("RooExtendPdf::%s(%s,%s)", pdfTotName.c_str(),
			     pdfName.c_str(),
			     ("N_"+p.first+lbl).c_str()
			     ))) { std::cout << "[ERROR] Failed to create extended PDF " << pdfTotName << std::endl; return false; }
	pdfListVarTot.add(*ws.pdf(pdfTotName.c_str()));
      }
      if (pdfListVarTot.getSize()>0) {
	const auto& pdfName = ("pdf"+varTot+"_Tot"+mainLabel);
	auto themodel = std::unique_ptr<RooAddPdf>(new RooAddPdf(pdfName.c_str(), pdfName.c_str(), pdfListVarTot));
	ws.import(*themodel);
      }
      pdfListTot.add(pdfListVarTot);
    }
    if (pdfListTot.getSize()>0) {
      const auto& pdfName = ("pdf"+varTot+"_Tot"+lbl);
      info.Par["pdfName"+chg] = pdfName;
      auto themodel = std::unique_ptr<RooAddPdf>(new RooAddPdf(pdfName.c_str(), pdfName.c_str(), pdfListTot));
      ws.import(*themodel);
    }
  }
  //
  return true;
};


RooAbsPdf* getTotalPDF(RooWorkspace& ws, const std::string& var, const std::string& label, const GlobalInfo& info)
{
  // Define input file name
  std::string fileName = "";
  auto dir = info.Par.at("outputDir");
  setFileName(fileName, dir, label, info, {var});
  const auto& inFileName = (dir+"result/FIT_"+fileName+".root");
  // Initialize the previous PDF
  const auto& cha = info.Par.at("channel");
  const auto& chg = label.substr(label.find(cha)+cha.size(), 2);
  GlobalInfo infoTMP(info);
  if (!addModel(ws, infoTMP, chg, {var})) { return NULL; }
  // Extract the previous PDF
  auto varT = var; stringReplace(varT, "_", "");
  const auto& pdfName = ("pdf"+varT+"_Tot"+label);
  if (!loadFitResult(ws, inFileName, pdfName)) { return NULL; }
  return ws.pdf(pdfName.c_str());
};


bool makeSPlotDS(RooWorkspace& ws, GlobalInfo& info, const std::string& var, const std::string& label)
{
  //
  const auto& cha = info.Par.at("channel");
  if (label.find(cha)==std::string::npos) { std::cout << "[ERROR] makeSPlotDS: Invalid channel " << cha << " in label " << label << std::endl; return false; }
  const auto& chg = label.substr(label.find(cha)+cha.size(), 2);
  // Check if sPlot DataSet has already been created
  const auto& dsName = info.Par.at("dsName"+chg);
  if (ws.data((dsName+"_SPLOT").c_str())) { return true; }
  const auto& obj = label.substr(0, label.find(cha));
  const auto& col = label.substr(label.find("_")+1);
  auto varT = var; stringReplace(varT, "_", "");
  //
  // Load the sPlot DataSet if already done
  std::string fileName = "";
  auto dir = info.Par.at("outputDir");
  setFileName(fileName, dir, label, info, {"Cand_Mass"});
  const auto& inFileName = (dir+"dataset/SPLOT_"+fileName+".root");
  const auto& foundDS = getDataSet(ws, dsName+"_SPLOT", inFileName);
  //
  if (!foundDS) {
    std::cout << "[INFO] Building the sPlot DataSets using the fitted mass PDF results!" << std::endl;
    //
    // Initialize the yields
    std::cout << "[INFO] Initializing " << var << " model yields for SPLOT" << std::endl;
    for (const auto& o : info.StrS.at("addObjectModel_"+varT+"_"+label)) {
      const auto& lbl = (o + cha + chg + "_" + col);
      if (!addModelPar(ws, info, {"N"}, var, lbl, "")) { return false; }
    }
    // Extract the yields
    RooArgList yieldList;
    if (ws.components().getSize()>0) {
      auto cmpIt = std::unique_ptr<TIterator>(ws.componentIterator());
      for (auto itp = cmpIt->Next(); itp!=NULL; itp = cmpIt->Next()) {
	const auto& it = dynamic_cast<RooAbsArg*>(itp); if (!it) continue;
	if (std::string(it->GetName()).rfind("N_",0)==0) { yieldList.add(*it); }
      }
    }
    if (yieldList.getSize()==0) { std::cout << "[ERROR] makeSPlotDS: Workspace has no yields!" << endl; return false; }
    //
    // Extract the RooDataSet
    auto data = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(ws.data(dsName.c_str())->Clone("TMP_DATA")));
    if (!data) { std::cout << "[ERROR] RooDataSet " << dsName << " was not found!" << endl; return false; }
    //
    // Extract the mass PDF
    const auto& massPDF = getTotalPDF(ws, "Cand_Mass", label, info);
    if (!massPDF) { std::cout << "[ERROR] makeSPlotDS: Cand_Mass PDF was not loaded!" << endl; return false; }
    auto cloneSet = std::unique_ptr<RooArgSet>(dynamic_cast<RooArgSet*>(RooArgSet(*massPDF, massPDF->GetName()).snapshot(kTRUE)));
    if (!cloneSet) { std::cout << "[ERROR] Couldn't deep-clone " << massPDF->GetName() << std::endl; return false; }
    const auto& clonePDF = dynamic_cast<RooAbsPdf*>(cloneSet->find(massPDF->GetName()));
    if (!clonePDF) { cout << "[ERROR] Couldn't deep-clone " << massPDF->GetName() << endl; return false; }
    clonePDF->setOperMode(RooAbsArg::ADirty, kTRUE);
    //
    // Create the sPlot DataSet
    std::cout << "[INFO] Creating the sPlot DataSets!" << std::endl;
    const auto& sData = RooStats::SPlot("sData", "An SPlot", *data, clonePDF, yieldList);
    if (ws.import(*data, RooFit::Rename((dsName+"_SPLOT").c_str()))) { std::cout << "[ERROR] sPlot DataSets were not imported!" << std::endl; return false; }
    else { std::cout << "[INFO] sPlot DataSets created succesfully!" << std::endl; }
    ws.loadSnapshot("loadedParameters");
    auto yIt = std::unique_ptr<TIterator>(yieldList.createIterator());
    for (auto it = yIt->Next(); it!=NULL; it = yIt->Next()) {
      const std::string& name = it->GetName();
      const auto& fitVal = ws.var(name.c_str())->getVal();
      const auto& sVal = sData.GetYieldFromSWeight(name.c_str());
      if (std::abs(fitVal - sVal)>0.1) { std::cout << "[ERROR] Variable " << name << " has different fitted (" << fitVal << ") and sPlot (" << sVal << ") results!" << std::endl; return false; }
    }
    // Store the sPlot DataSet
    const auto& dsSPlot = dynamic_cast<RooDataSet*>(ws.data((dsName+"_SPLOT").c_str()));
    if (!saveDataSet(*dsSPlot, dir, fileName)) { return false; }
  }
  //
  return true;
};


#endif // #ifndef Candidate_addModel_C
