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


bool makeSPlotDS ( RooWorkspace& ws , GlobalInfo& info       , const std::string& label );


bool addModel(RooWorkspace& ws, GlobalInfo& info, const std::string& chg, const StringSet_t& varV={})
{
  //
  const auto& cha = info.Par.at("channel");
  auto varS = info.StrS.at("fitCondVariable"); if (!varV.empty()) { varS = varV; }
  auto varT = info.StrS.at("fitCondVarName"); if (!varV.empty()) { varT.clear(); for (const auto& v : varV) { auto t = v; stringReplace(t, "_", ""); varT.insert(t); } }
  std::string varTot = ""; for (const auto& v : varT) { varTot += v; }
  for (const auto& col : info.StrS.at("fitSystem")) {
    const auto& lbl = cha + chg + "_" + col;
    //
    std::map<std::string, std::pair<RooArgList, RooArgList>> pdfMapTot;
    for (const auto& mainObj : info.StrS.at("fitObject")) {
      const auto& mainTag = mainObj + cha + chg;
      const auto& mainLabel = mainTag + "_" + col;
      //
      std::map<std::string, RooArgList> pdfListVar, condPdfListVar;
      std::map<std::string,  std::map<std::string, RooArgList> > pdfMapVar;
      for (const auto& varName : varS) {
	const std::string& varWindow   = "FitWindow";
	const auto& varNormName = (varName+"Norm");
	auto varType = varName; stringReplace(varType, "_", "");
	std::cout << "[INFO] Implementing " << mainTag << " " << varName << " Model for " << col << std::endl;
	// Make sure that DLenRes and main object are done first
	StringVector_t objV;
	if (contain(info.StrS.at("addObjectModel_"+varType+"_"+mainLabel), "DLenRes")) { objV.push_back("DLenRes"); }
	if (contain(info.StrS.at("addObjectModel_"+varType+"_"+mainLabel), mainObj)) { objV.push_back(mainObj); }
	for (const auto& m : info.StrS.at("addObjectModel_"+varType+"_"+mainLabel)) { if (!contain(objV, m)) { objV.push_back(m); } }
	std::map<std::string, RooArgList> pdfList;
	for (const auto& obj : objV) {
	  const auto& modelN = info.Par.at("Model"+varType+"_"+mainLabel+"_"+obj);
	  const auto& tag = obj + cha + chg;
	  const auto& label = tag + "_" + col;
	  auto objI = obj.substr(0, obj.find("Cat"));
	  const auto& labelI = objI + cha + chg + "_" + col;
	  addString(ws, ("Model"+varType+"_"+label), modelN); // Save the model name for bookkeeping
	  //
	  const std::string& pdfName    = Form("pdf%s_%s",    varType.c_str(), label.c_str());
	  const std::string& pdf1Name   = Form("pdf%s1_%s",   varType.c_str(), label.c_str());
	  const std::string& pdf2Name   = Form("pdf%s2_%s",   varType.c_str(), label.c_str());
	  const std::string& pdf3Name   = Form("pdf%s3_%s",   varType.c_str(), label.c_str());
	  const std::string& pdf4Name   = Form("pdf%s4_%s",   varType.c_str(), label.c_str());
	  const std::string& pdf5Name   = Form("pdf%s5_%s",   varType.c_str(), label.c_str());
	  const std::string& pdfPolName = Form("pdf%sPol_%s", varType.c_str(), label.c_str());
	  const std::string& pdfResName = Form("pdf%s_DLenRes%s", varType.c_str(), lbl.c_str());
	  //
	  // Create Models
	  switch(ModelDictionary.at(modelN).first)
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
		const auto& proceed = histToPdf(ws, pdfName, dsName, varName, range, "HIST");
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
		// load model from previous fit, if done
		if (getPDFData(ws, pdfName, info, label, varName)) {
		  if (!histToPdf(ws, pdfName, varName, "HIST")) { return false; }
		}
		else {
		  // create the sPlot DataSet
		  if (!makeSPlotDS(ws, info, label)) { return false; }
		  // create the PDF
		  const auto& sPlotDS = info.Par.at("dsSPlotNameFit"+chg);
		  const auto& sPlotN = ("N_"+labelI+"_sw");
		  const auto& range = std::vector<double>({double(getNBins(varName, info)), info.Var.at(varName).at("Min"), info.Var.at(varName).at("Max")});
		  std::cout << "[INFO] Using " << sPlotN << " of " << sPlotDS << " to create " << pdfName << " SPlot template" << std::endl;
		  auto dataw = std::unique_ptr<RooDataSet>(new RooDataSet("TMP","TMP", dynamic_cast<RooDataSet*>(ws.data(sPlotDS.c_str())), RooArgSet(*ws.var(varName.c_str()), *ws.var(sPlotN.c_str())), 0, sPlotN.c_str()));
		  if (!histToPdf(ws, pdfName, *dataw, varName, range, "HIST")) { return false; }
		}
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " " << " binned SPLot " << " " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::KEYS)):
	      {
		// load model from previous fit, if done
		if (getPDFData(ws, pdfName, info, label, varName)) {
		  if (!histToPdf(ws, pdfName, varName, "KEYS")) { return false; }
		}
		else {
		  // create the PDF using RooKeys from binned dataset
		  const auto& dsName = info.Par.at("dsNameFit"+chg);
		  const auto& range = std::vector<double>({((double)getNBins(varName, info)), info.Var.at(varName).at("Min"), info.Var.at(varName).at("Max")});
		  std::cout << "[INFO] Using " << dsName << " to create " << pdfName << " binned RooKeysPdf template" << std::endl;
		  auto data = std::unique_ptr<RooDataSet>(new RooDataSet("TMP","TMP", dynamic_cast<RooDataSet*>(ws.data(dsName.c_str())), RooArgSet(*ws.var(varName.c_str()))));
		  if (!histToPdf(ws, pdfName, *data, varName, range, "KEYS")) { return false; }
		}
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " " << " binned RooKeysPdf " << " " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::DSKEYS)):
	      {
		// TODO: implement loading of previous RooKeysPdf
		if (true) {
		  // create the PDF using RooKeys from unbinned dataset
		  const auto& dsName = info.Par.at("dsNameFit"+chg);
		  std::cout << "[INFO] Using " << dsName << " to create " << pdfName << " unbinned RooKeysPdf template" << std::endl;
		  auto data = dynamic_cast<RooDataSet*>(ws.data(dsName.c_str()));
		  auto pdf = std::unique_ptr<RooKeysPdf>(new RooKeysPdf(pdfName.c_str(), pdfName.c_str(), *ws.var(varName.c_str()), *data, RooKeysPdf::MirrorAsymBoth));
		  if (ws.import(*pdf)) { std::cout << "[ERROR] RooKeysPdf " << pdfName << " failed to import!" << std::endl; }
		}
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " " << " unbinned RooKeysPdf " << " " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::HIST)):
	      {
		// load model from previous fit, if done
		if (getPDFData(ws, pdfName, info, label, varName)) {
		  if (!histToPdf(ws, pdfName, varName, "HIST")) { return false; }
		}
		else {
		  // create the PDF using RooHistPdf
		  const auto& dsName = info.Par.at("dsNameFit"+chg);
		  const auto& range = std::vector<double>({((double)getNBins(varName, info)), info.Var.at(varName).at("Min"), info.Var.at(varName).at("Max")});
		  std::cout << "[INFO] Using " << dsName << " to create " << pdfName << " RooHistPdf template" << std::endl;
		  auto data = std::unique_ptr<RooDataSet>(new RooDataSet("TMP","TMP", dynamic_cast<RooDataSet*>(ws.data(dsName.c_str())), RooArgSet(*ws.var(varName.c_str()))));
		  if (!histToPdf(ws, pdfName, *data, varName, range, "HIST")) { return false; }
		}
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " " << " RooHistPdf " << " " << varName << " PDF in " << col << " added!" << std::endl; break;
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
	    case (int(Model::TripleGaussian)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "f", "rSigma32", "Sigma3", "f2", "recf2"};
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
		if (!ws.factory(Form("Gaussian::%s(%s, %s, %s)", pdf3Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma3_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf3Name << std::endl; return false; }
		// Sum the PDFs to get the signal PDF
		if (!ws.factory(Form("SUM::%s(%s*%s, %s*%s, %s)", pdfName.c_str(),
				     ("f_"+label).c_str(), 
				     pdf1Name.c_str(),
				     ("recf2_"+label).c_str(),
				     pdf2Name.c_str(),
				     pdf3Name.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Triple Gaussian " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::QuadrupleGaussian)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "f", "rSigma32", "Sigma3", "f2", "recf2", "rSigma43", "Sigma4", "f3", "recf3"};
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
		if (!ws.factory(Form("Gaussian::%s(%s, %s, %s)", pdf3Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma3_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf3Name << std::endl; return false; }
		if (!ws.factory(Form("Gaussian::%s(%s, %s, %s)", pdf4Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma4_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf4Name << std::endl; return false; }
		// Sum the PDFs to get the signal PDF
		if (!ws.factory(Form("SUM::%s(%s*%s, %s*%s, %s*%s, %s)", pdfName.c_str(),
				     ("f_"+label).c_str(),
				     pdf1Name.c_str(),
				     ("recf2_"+label).c_str(),
				     pdf2Name.c_str(),
				     ("recf3_"+label).c_str(),
				     pdf3Name.c_str(),
				     pdf4Name.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Quadruple Gaussian " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::PentaGaussian)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "f", "rSigma32", "Sigma3", "f2", "recf2", "rSigma43", "Sigma4", "f3", "recf3", "rSigma54", "Sigma5", "f4", "recf4"};
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
		if (!ws.factory(Form("Gaussian::%s(%s, %s, %s)", pdf3Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma3_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf3Name << std::endl; return false; }
		if (!ws.factory(Form("Gaussian::%s(%s, %s, %s)", pdf4Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma4_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf4Name << std::endl; return false; }
		if (!ws.factory(Form("Gaussian::%s(%s, %s, %s)", pdf5Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma5_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf5Name << std::endl; return false; }
		// Sum the PDFs to get the signal PDF
		if (!ws.factory(Form("SUM::%s(%s*%s, %s*%s, %s*%s, %s*%s, %s)", pdfName.c_str(),
				     ("f_"+label).c_str(),
				     pdf1Name.c_str(),
				     ("recf2_"+label).c_str(),
				     pdf2Name.c_str(),
				     ("recf3_"+label).c_str(),
				     pdf3Name.c_str(),
				     ("recf4_"+label).c_str(),
				     pdf4Name.c_str(),
				     pdf5Name.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Penta Gaussian " << varName << " PDF in " << col << " added!" << std::endl; break;
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
	    case (int(Model::DoubleGaussianAndCrystalBall)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "rSigma32", "Sigma3", "Alpha", "n", "f", "f2", "recf2"};
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
		if (!ws.factory(Form("Gaussian::%s(%s, %s, %s)", pdf3Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma3_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf3Name << std::endl; return false; }
		// Sum the PDFs to get the signal PDF
		if (!ws.factory(Form("SUM::%s(%s*%s, %s*%s, %s)", pdfName.c_str(),
				     ("f_"+label).c_str(), 
				     pdf1Name.c_str(),
				     ("recf2_"+label).c_str(),
				     pdf2Name.c_str(),
				     pdf3Name.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Double Gaussian and Crystal Ball " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::SingleExtCrystalBall)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "Alpha", "n", "rAlphaR", "AlphaR", "nR"};
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
		const StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "Alpha", "n", "rAlphaR", "AlphaR", "nR", "Alpha2", "n2", "AlphaR2", "nR2", "f"};
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
		const StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "Alpha", "n", "rAlphaR", "AlphaR", "rnR", "nR", "f"};
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
	    case (int(Model::DoubleGaussianAndExtCrystalBall)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "rSigma32", "Sigma3", "Alpha", "n", "rAlphaR", "AlphaR", "nR", "f", "f2", "recf2"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
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
		if (!ws.factory(Form("Gaussian::%s(%s, %s, %s)", pdf3Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma3_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf3Name << std::endl; return false; }
		// Sum the PDFs to get the signal PDF
		if (!ws.factory(Form("SUM::%s(%s*%s, %s*%s, %s)", pdfName.c_str(),
				     ("f_"+label).c_str(), 
				     pdf1Name.c_str(),
				     ("recf2_"+label).c_str(),
				     pdf2Name.c_str(),
				     pdf3Name.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Double Gaussian and Extended Crystal Ball " << varName << " PDF in " << col << " added!" << std::endl; break;
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
	    case (int(Model::SingleModExtCrystalBall)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "Alpha", "n", "rAlphaR", "AlphaR"};
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
		const StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "Alpha", "n", "rAlphaR", "AlphaR", "Alpha2", "n2", "AlphaR2", "f"};
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
		const StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "Alpha", "n", "rAlphaR", "AlphaR", "f"};
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
		std::cout << "[INFO] " << tag << " Truth " << varName << " Resolution Model in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::SingleGaussianResolution)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1"};
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
		std::cout << "[INFO] " << tag << " Single Gaussian " << varName << " Resolution Model in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::DoubleGaussianResolution)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "f"};
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
		std::cout << "[INFO] " << tag << " Double Gaussian " << varName << " Resolution Model in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::TripleGaussianResolution)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "f", "rSigma32", "Sigma3", "f2", "recf2"};
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
				     ("recf2_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		std::cout << "[INFO] " << tag << " Triple Gaussian " << varName << " Resolution Model in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::QuadrupleGaussianResolution)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "f", "rSigma32", "Sigma3", "f2", "recf2", "rSigma43", "Sigma4", "f3", "recf3"};
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
				     ("recf2_"+label).c_str(),
				     ("recf3_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		std::cout << "[INFO] " << tag << " Quadruple Gaussian " << varName << " Resolution Model in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::PentaGaussianResolution)):
	      {
		// input variables
		const StringVector_t parNames = {"m", "Sigma1", "rSigma21", "Sigma2", "f", "rSigma32", "Sigma3", "f2", "recf2", "rSigma43", "Sigma4", "f3", "recf3", "rSigma54", "Sigma5", "f4", "recf4"};
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
		if (!ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", pdf5Name.c_str(), varName.c_str(),
				     ("m_"+label).c_str(),
				     ("Sigma5_"+label).c_str(),
				     (varName=="Cand_DLenRes" ? "One" : "Cand_DLenErr")
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf5Name << std::endl; return false; }
		// sum the PDFs to get the total PDF
		if (!ws.factory(Form("AddModel::%s({%s, %s, %s, %s, %s}, {%s, %s, %s, %s})", pdfName.c_str(),
				     pdf1Name.c_str(),
				     pdf2Name.c_str(),
				     pdf3Name.c_str(),
				     pdf4Name.c_str(),
				     pdf5Name.c_str(),
				     ("f_"+label).c_str(),
				     ("recf2_"+label).c_str(),
				     ("recf3_"+label).c_str(),
				     ("recf4_"+label).c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		std::cout << "[INFO] " << tag << " Penta Gaussian " << varName << " Resolution Model in " << col << " added!" << std::endl; break;
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
		const StringVector_t parNames = {"LambdaSS"};
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
            case (int(Model::DoubleSingleSidedDecay)): 
              { 
                // input variables
                const StringVector_t parNames = {"LambdaSS", "rLambdaSS21", "LambdaSS2", "f"};
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
                // sum the PDFs to get the total PDF
                if (!ws.factory(Form("SUM::%s(%s*%s, %s)", pdfName.c_str(),
                                     ("f_"+label).c_str(),
                                     pdf1Name.c_str(),
                                     pdf2Name.c_str()
                                     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
                ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
                // add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
                std::cout << "[INFO] " << tag << " Double Single Sided Decay " << varName << " PDF in " << col << " added!" << std::endl; break;
              }
	    case (int(Model::TripleDecay)):
	      {
		// input variables
		const StringVector_t parNames = {"LambdaSS", "LambdaF", "LambdaDS", "f", "f2", "recf2"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the PDF
		if (!ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", pdf1Name.c_str(), varName.c_str(),
				     ("LambdaSS_"+label).c_str(),
				     pdfResName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf1Name << std::endl; return false; }
		if (!ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::DoubleSided)", pdf2Name.c_str(), varName.c_str(),
				     ("LambdaDS_"+label).c_str(),
				     pdfResName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf2Name << std::endl; return false; }
		if (!ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", pdf3Name.c_str(), varName.c_str(),
				     ("LambdaF_"+label).c_str(),
				     pdfResName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf3Name << std::endl; return false; }
		// sum the PDFs to get the total PDF
		if (!ws.factory(Form("SUM::%s(%s*%s, %s*%s, %s)", pdfName.c_str(),
				     ("f_"+label).c_str(),
				     pdf1Name.c_str(),
				     ("recf2_"+label).c_str(),
				     pdf2Name.c_str(),
				     pdf3Name.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Triple Decay " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::QuadrupleDecay)):
	      {
		// input variables
		const StringVector_t parNames = {"LambdaSS", "rLambdaSS21", "LambdaSS2", "LambdaF", "LambdaDS", "f", "f2", "recf2", "f3", "recf3"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the PDF
		if (!ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::DoubleSided)", pdf1Name.c_str(), varName.c_str(),
				     ("LambdaDS_"+label).c_str(),
				     pdfResName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf1Name << std::endl; return false; }
		if (!ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", pdf2Name.c_str(), varName.c_str(),
				     ("LambdaSS2_"+label).c_str(),
				     pdfResName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf2Name << std::endl; return false; }
		if (!ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", pdf3Name.c_str(), varName.c_str(),
				     ("LambdaSS_"+label).c_str(),
				     pdfResName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf3Name << std::endl; return false; }
		if (!ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", pdf4Name.c_str(), varName.c_str(),
				     ("LambdaF_"+label).c_str(),
				     pdfResName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf4Name << std::endl; return false; }
		// sum the PDFs to get the total PDF
		if (!ws.factory(Form("SUM::%s(%s*%s, %s*%s, %s*%s, %s)", pdfName.c_str(),
				     ("f_"+label).c_str(),
				     pdf1Name.c_str(),
				     ("recf2_"+label).c_str(),
				     pdf2Name.c_str(),
				     ("recf3_"+label).c_str(),
				     pdf3Name.c_str(),
				     pdf4Name.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Quadruple Decay " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::QuadrupleDecay2)):
	      {
		// input variables
		const StringVector_t parNames = {"LambdaSS", "LambdaF", "rLambdaF21", "LambdaF2", "LambdaDS", "f", "f2", "recf2", "f3", "recf3"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the PDF
		if (!ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", pdf1Name.c_str(), varName.c_str(),
				     ("LambdaSS_"+label).c_str(),
				     pdfResName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf1Name << std::endl; return false; }
		if (!ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::DoubleSided)", pdf2Name.c_str(), varName.c_str(),
				     ("LambdaDS_"+label).c_str(),
				     pdfResName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf2Name << std::endl; return false; }
		if (!ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", pdf3Name.c_str(), varName.c_str(),
				     ("LambdaF2_"+label).c_str(),
				     pdfResName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf3Name << std::endl; return false; }
		if (!ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", pdf4Name.c_str(), varName.c_str(),
				     ("LambdaF_"+label).c_str(),
				     pdfResName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf4Name << std::endl; return false; }
		// sum the PDFs to get the total PDF
		if (!ws.factory(Form("SUM::%s(%s*%s, %s*%s, %s*%s, %s)", pdfName.c_str(),
				     ("f_"+label).c_str(),
				     pdf1Name.c_str(),
				     ("recf2_"+label).c_str(),
				     pdf2Name.c_str(),
				     ("recf3_"+label).c_str(),
				     pdf3Name.c_str(),
				     pdf4Name.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdfName << std::endl; return false; }
		ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
		// add PDF to list
		pdfList[objI].add(*ws.pdf(pdfName.c_str()));
		std::cout << "[INFO] " << tag << " Quadruple Decay " << varName << " PDF in " << col << " added!" << std::endl; break;
	      }
	    case (int(Model::QuadrupleDecay3)):
	      {
		// input variables
		const StringVector_t parNames = {"LambdaSS", "LambdaF", "LambdaDS", "rLambdaDS21", "LambdaDS2", "f", "f2", "recf2", "f3", "recf3"};
		// create the variables for this model
		if (!addModelPar(ws, info, parNames, varName, label, modelN)) { return false; }
		// create the PDF
		if (!ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", pdf1Name.c_str(), varName.c_str(),
				     ("LambdaSS_"+label).c_str(),
				     pdfResName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf1Name << std::endl; return false; }
		if (!ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::DoubleSided)", pdf2Name.c_str(), varName.c_str(),
				     ("LambdaDS_"+label).c_str(),
				     pdfResName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf2Name << std::endl; return false; }
		if (!ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::DoubleSided)", pdf3Name.c_str(), varName.c_str(),
				     ("LambdaDS2_"+label).c_str(),
				     pdfResName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf3Name << std::endl; return false; }
		if (!ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", pdf4Name.c_str(), varName.c_str(),
				     ("LambdaF_"+label).c_str(),
				     pdfResName.c_str()
				     ))) { std::cout << "[ERROR] Failed to create PDF " << pdf4Name << std::endl; return false; }
		// sum the PDFs to get the total PDF
		if (!ws.factory(Form("SUM::%s(%s*%s, %s*%s, %s*%s, %s)", pdfName.c_str(),
				     ("f_"+label).c_str(),
				     pdf1Name.c_str(),
				     ("recf2_"+label).c_str(),
				     pdf2Name.c_str(),
				     ("recf3_"+label).c_str(),
				     pdf3Name.c_str(),
				     pdf4Name.c_str()
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
	    StringVector_t parFullNames, parNames = {"b"};
	    for (int i=3; i<=p.second.getSize(); i++) { parNames.push_back(Form("b%d",i-1)); }
	    if (!addModelPar(ws, info, parNames, varName, (p.first+lbl), ("SUM_"+p.first), parFullNames)) { return false; }
	    RooArgList coefList; for (const auto& par : parFullNames) { if (ws.var(par.c_str())) { coefList.add(*ws.var(par.c_str())); } }
	    // Sum the PDFs
	    auto pdf = std::unique_ptr<RooAddPdf>(new RooAddPdf(pdfName.c_str(), pdfName.c_str(), p.second, coefList));
	    if (!pdf) { std::cout << "[ERROR] RooAddPdf " << pdfName << " is NULL!" << std::endl; return false; }
	    if (ws.import(*pdf)) { std::cout << "[ERROR] RooAddPdf " << pdfName << " was not imported!" << std::endl; return false; }
	    const bool& isCondPdf = (varType=="CandDLen" && info.Flag.at("condCand_DLenErr"));
	    if (isCondPdf) { condPdfListVar[p.first].add(*ws.pdf(pdfName.c_str())); } else { pdfListVar[p.first].add(*ws.pdf(pdfName.c_str())); }
	    if (!isCondPdf) { pdfMapVar[varType][p.first].add(*ws.pdf(pdfName.c_str())); }
	  }
	  else { pdfListVar[p.first].add(p.second); pdfMapVar[varType][p.first].add(p.second); }
	}
      }
      StringSet_t varList = {varTot}; for (const auto& v : pdfMapVar) { varList.insert(v.first); }
      for (const auto& var : varList) {
	auto pdfListVarTot = std::pair<RooArgList, RooArgList>(RooArgList(), RooArgList());
	for (const auto& p : (var==varTot ? pdfListVar : pdfMapVar[var])) {
	  // Multiply the PDFs
	  std::string pdfName = p.second.at(0)->GetName();
	  const auto& condPdfSize = (var==varTot ? condPdfListVar[p.first].getSize() : 0);
	  if ((p.second.getSize()+condPdfSize)>1) {
	    pdfName = ("pdf"+var+"_"+p.first+lbl);
	    std::unique_ptr<RooProdPdf> pdf;
	    if (condPdfSize>0) {
	      pdf.reset(new RooProdPdf(pdfName.c_str(), pdfName.c_str(), p.second, RooFit::Conditional(condPdfListVar[p.first], RooArgSet(*ws.var("Cand_DLen")))));
	    }
	    else { pdf.reset(new RooProdPdf(pdfName.c_str(), pdfName.c_str(), p.second)); }
	    if (!pdf) { std::cout << "[ERROR] RooProdPdf " << pdfName << " is NULL!" << std::endl; return false; }
	    if (ws.import(*pdf)) { std::cout << "[ERROR] RooProdPdf " << pdfName << " was not imported!" << std::endl; return false; }
	  }
	  // create the yield
	  if (!addModelPar(ws, info, {"N"}, "Cand_Mass", (p.first+lbl), "", true)) { return false; }
	  // extend the PDF
	  const auto& pdfTotName = ("pdf"+var+"Tot_"+p.first+lbl);
	  if (!ws.factory(Form("RooExtendPdf::%s(%s,%s)", pdfTotName.c_str(),
			       pdfName.c_str(),
			       ("N_"+p.first+lbl).c_str()
			       ))) { std::cout << "[ERROR] Failed to create extended PDF " << pdfTotName << std::endl; return false; }
	  pdfListVarTot.first.add(*ws.pdf(pdfName.c_str()));
	  pdfListVarTot.second.add(*ws.arg(("N_"+p.first+lbl).c_str()));
	}
	if (pdfListVarTot.first.getSize()>0) {
	  const auto& pdfName = ("pdf"+var+"_Tot"+mainLabel);
	  auto pdf = std::unique_ptr<RooAddPdf>(new RooAddPdf(pdfName.c_str(), pdfName.c_str(), pdfListVarTot.first, pdfListVarTot.second));
	  if (!pdf) { std::cout << "[ERROR] RooAddPdf " << pdfName << " is NULL!" << std::endl; return false; }
	  if (ws.import(*pdf)) { std::cout << "[ERROR] RooAddPdf " << pdfName << " was not imported!" << std::endl; return false; }
	}
	pdfMapTot[var].first.add(pdfListVarTot.first);
	pdfMapTot[var].second.add(pdfListVarTot.second);
      }
    }
    for (const auto& v : pdfMapTot) {
      if (v.second.first.getSize()>0) {
	const auto& pdfName = ("pdf"+v.first+"_Tot"+lbl);
	if (v.first==varTot) { info.Par["pdfName"+chg] = pdfName; }
	auto pdf = std::unique_ptr<RooAddPdf>(new RooAddPdf(pdfName.c_str(), pdfName.c_str(), v.second.first, v.second.second));
	if (!pdf) { std::cout << "[ERROR] RooAddPdf " << pdfName << " is NULL!" << std::endl; return false; }
	if (ws.import(*pdf)) { std::cout << "[ERROR] RooAddPdf " << pdfName << " was not imported!" << std::endl; return false; }
      }
    }
    // Add N for other objects not used in the PDF fit if mass PDF info has been loaded
    if (contain(info.StrS, "incObject_CandMass")) {
      for (const auto& o : info.StrS.at("incObject_CandMass")) {
	if (!addModelPar(ws, info, {"N"}, "", o+lbl, "")) { return false; }
      }
    }
  }
  //
  return true;
};


bool loadPDFParameters(RooWorkspace& ws, const std::string& var, const std::string& label, const GlobalInfo& info,
		       const bool& loadYield=false, const bool& fixPar=true)
{
  // Define input file name
  std::string fileName = "";
  auto dir = info.Par.at("outputDir");
  setFileName(fileName, dir, label, info, {var});
  auto varT = var; stringReplace(varT, "_", "");
  const auto& inFileName = (dir+"result/FIT_"+varT+"_"+fileName+".root");
  // Extract the previous PDF results
  const auto& cha = info.Par.at("channel");
  const auto& chg = label.substr(label.find(cha)+cha.size(), 2);
  const auto& pdfName = info.Par.at("pdfName"+chg);
  if (!loadFitResult(ws, inFileName, pdfName, var, loadYield, fixPar)) { return false; }
  return true;
};


RooAbsPdf* getTotalPDF(RooWorkspace& ws, const std::string& var, const std::string& label, const GlobalInfo& info)
{
  // Define input file name
  std::string fileName = "";
  auto dir = info.Par.at("outputDir");
  setFileName(fileName, dir, label, info, {var});
  auto varT = var; stringReplace(varT, "_", "");
  const auto& inFileName = (dir+"result/FIT_"+varT+"_"+fileName+".root");
  // Initialize the previous PDF
  const auto& cha = info.Par.at("channel");
  const auto& chg = label.substr(label.find(cha)+cha.size(), 2);
  GlobalInfo infoTMP(info);
  if (!addModel(ws, infoTMP, chg, {var})) { return NULL; }
  // Extract the previous PDF
  const auto& pdfName = ("pdf"+varT+"_Tot"+label);
  if (!loadFitResult(ws, inFileName, pdfName, var)) { return NULL; }
  return ws.pdf(pdfName.c_str());
};


bool makeSPlotDS(RooWorkspace& ws, GlobalInfo& info, const std::string& label)
{
  //
  const auto& cha = info.Par.at("channel");
  if (label.find(cha)==std::string::npos) { std::cout << "[ERROR] makeSPlotDS: Invalid channel " << cha << " in label " << label << std::endl; return false; }
  const auto& chg = label.substr(label.find(cha)+cha.size(), 2);
  // Check if sPlot DataSet has already been created
  const auto& dsName = info.Par.at("dsName"+chg);
  const auto& dsSPlotName = (dsName.rfind("_sPlot")==std::string::npos ? (dsName+"_sPlot") : dsName);
  info.Par["dsSPlotName"+chg] = dsSPlotName;
  if (ws.data(dsSPlotName.c_str())) { std::cout << "[INFO] sPlot dataset " << dsSPlotName << " already made!" << std::endl; return true; }
  const auto& dsNameFit = info.Par.at("dsNameFit"+chg);
  const auto& dsSPlotNameFit = (dsNameFit.rfind("_sPlot")==std::string::npos ? (dsNameFit!=dsName ? (dsName+"_sPlot_FIT") : (dsName+"_sPlot")) : dsNameFit);
  info.Par["dsSPlotNameFit"+chg] = dsSPlotNameFit;
  if (ws.data(dsSPlotNameFit.c_str())) { std::cout << "[INFO] sPlot dataset " << dsSPlotNameFit << " already made!" << std::endl; return true; }
  const auto& col = label.substr(label.rfind("_")+1);
  const auto& sumEntries = ws.data(dsName.c_str())->sumEntries();
  const auto& numEntries = ws.data(dsName.c_str())->numEntries();
  //
  // Load the sPlot DataSet if already done
  std::string fileName = "";
  auto dir = info.Par.at("outputDir");
  setFileName(fileName, dir, label, info);
  const auto& outDir = (dir+"dataset/");
  const auto& inFileName = (outDir+"SPlot_"+fileName+".root");
  const auto& foundDS = getDataSet(ws, dsSPlotName, inFileName);
  //
  if (!foundDS) {
    std::cout << "[INFO] Building the sPlot datasets using the fitted mass PDF results!" << std::endl;
    //
    // Extract the RooDataSet
    if (!ws.data(dsName.c_str())) { std::cout << "[ERROR] RooDataSet " << dsName << " was not found!" << endl; return false; }
    auto data = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(ws.data(dsName.c_str())->Clone("TMP_DATA")));
    if (!data) { std::cout << "[ERROR] RooDataSet " << dsName << " was not cloned!" << endl; return false; }
    //
    // Extract the mass PDF
    info.Flag["notConstrainYields"] = true;
    const auto& massPDF = getTotalPDF(ws, "Cand_Mass", label, info);
    if (!massPDF) { std::cout << "[ERROR] makeSPlotDS: Cand_Mass PDF was not loaded!" << endl; return false; }
    auto cloneSet = std::unique_ptr<RooArgSet>(dynamic_cast<RooArgSet*>(RooArgSet(*massPDF, massPDF->GetName()).snapshot(kTRUE)));
    if (!cloneSet) { std::cout << "[ERROR] Couldn't deep-clone " << massPDF->GetName() << std::endl; return false; }
    const auto& clonePDF = dynamic_cast<RooAbsPdf*>(cloneSet->find(massPDF->GetName()));
    if (!clonePDF) { cout << "[ERROR] Couldn't deep-clone " << massPDF->GetName() << endl; return false; }
    clonePDF->setOperMode(RooAbsArg::ADirty, kTRUE);
    //
    // Extract the yields
    RooArgList yieldList;
    auto parSet = std::unique_ptr<RooArgSet>(clonePDF->getParameters(RooArgSet()));
    if (parSet->getSize()>0) {
      auto parIt = std::unique_ptr<TIterator>(parSet->createIterator());
      for (auto itp = parIt->Next(); itp!=NULL; itp = parIt->Next()) {
	const auto& it = dynamic_cast<RooRealVar*>(itp); if (!it) continue;
	if (std::string(it->GetName()).rfind("N_",0)==0) {
	  it->setMin(0.0); // Set minimum range of yields to zero
	  yieldList.add(*it);
	}
      }
    }
    if (yieldList.getSize()==0) { std::cout << "[ERROR] makeSPlotDS: Workspace has no yields!" << endl; return false; }
    //
    // Create sPlot dataset
    std::cout << "[INFO] Creating the sPlot datasets!" << std::endl;
    const auto& sData = RooStats::SPlot("sData", "An SPlot", *data, clonePDF, yieldList);
    // Remove extra variables
    auto skimVars = *data->get();
    auto varIt = std::unique_ptr<TIterator>(data->get()->createIterator());
    for (auto itp = varIt->Next(); itp!=NULL; itp = varIt->Next()) {
      const std::string& name = itp->GetName();
      if (name.rfind("L_N_", 0)==0) { const auto& var = skimVars.find(name.c_str()); if (var) { skimVars.remove(*var); } }
    }
    auto skimData = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(data->reduce(RooFit::SelectVars(skimVars))));
    auto skimDataFit = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(skimData->reduce(RooFit::Cut(getString(ws, "cutDSFit").c_str()))));
    // Import sPlot dataset
    if (ws.import(*skimData, RooFit::Rename(dsSPlotName.c_str()))) { std::cout << "[ERROR] sPlot datasets were not imported!" << std::endl; return false; }
    else {
      if (ws.data(dsName.c_str())) { ws.RecursiveRemove(ws.data(dsName.c_str())); }
      info.Var.at("numEntries")[dsSPlotName] = ws.data(dsSPlotName.c_str())->sumEntries();
      info.Par.at("dsName"+chg) = dsSPlotName;
      setDSParamaterRange(ws, dsSPlotName, info, "Full_");
      std::cout << "[INFO] sPlot datasets created succesfully!" << std::endl;
    }
    if (skimDataFit->numEntries()>0 && skimDataFit->numEntries()<skimData->numEntries()) {
      if (ws.import(*skimDataFit, RooFit::Rename(dsSPlotNameFit.c_str()))) { std::cout << "[ERROR] sPlot fit datasets were not imported!" << std::endl; return false; }
      else {
	if (ws.data(dsNameFit.c_str())) { ws.RecursiveRemove(ws.data(dsNameFit.c_str())); }
	info.Var.at("numEntries")[dsSPlotNameFit] = ws.data(dsSPlotNameFit.c_str())->sumEntries();
	info.Par.at("dsNameFit"+chg) = dsSPlotNameFit;
	setDSParamaterRange(ws, dsSPlotNameFit, info);
	std::cout << "[INFO] sPlot fit datasets created succesfully!" << std::endl;
      }
    }
    else if (skimDataFit->numEntries()<=0) { std::cout << "[ERROR] No events from dataset " <<  dsSPlotNameFit << " passed kinematic cuts!" << std::endl; return false; }
    else { info.Par.at("dsNameFit"+chg) = dsSPlotName; }
    // Check sPlot results
    ws.loadSnapshot("loadedParameters");
    auto yIt = std::unique_ptr<TIterator>(yieldList.createIterator());
    for (auto it = yIt->Next(); it!=NULL; it = yIt->Next()) {
      const std::string& name = it->GetName();
      const auto& fitVal = ws.var(name.c_str())->getVal();
      const auto& fitUnc = ws.var(name.c_str())->getError();
      const auto& sVal = sData.GetYieldFromSWeight(name.c_str());
      if (std::abs(fitVal - sVal)>0.005*fitUnc) { std::cout << "[ERROR] Variable " << name << " has different fitted (" << fitVal << " +- " << fitUnc << ") and sPlot (" << sVal << ") results!" << std::endl; return false; }
    }
    // Store the sPlot dataset
    const auto& sPlotDS = dynamic_cast<RooDataSet*>(ws.data(dsSPlotName.c_str()));
    if (!saveDataSet(*sPlotDS, inFileName)) { return false; }
    else { std::cout << "[INFO] Created " << dsSPlotName << " with " << sPlotDS->numEntries() << " events (" << numEntries << " origDS events)" << " and " << sPlotDS->sumEntries() << " wevents (" << sumEntries << " origDS wevents)" << std::endl; }
  }
  else {
    std::cout << "[INFO] SPlot dataset " << dsSPlotName << " found!" << std::endl;
    auto data = dynamic_cast<RooDataSet*>(ws.data(dsSPlotName.c_str()));
    setDSParamaterRange(*data, info, "Full_");
    if (data->numEntries()<=0) { std::cout << "[ERROR] No events from dataset " <<  dsSPlotName << " passed kinematic cuts!" << std::endl; return false; }
    else if (!isCompatibleDataset(*data, *dynamic_cast<RooDataSet*>(ws.data(dsName.c_str())))) { std::cout << "[ERROR] sPlot and original datasets are inconsistent!" << std::endl; return false; }
    else {
      if (ws.data(dsName.c_str())) { ws.RecursiveRemove(ws.data(dsName.c_str())); }
      info.Var.at("numEntries")[dsSPlotName] = ws.data(dsSPlotName.c_str())->sumEntries();
      info.Par.at("dsName"+chg) = dsSPlotName;
      defineSet(ws, "SET_"+dsSPlotName, *data->get());
      std::cout << "[INFO] Imported " << dsSPlotName << " with " << data->numEntries() << " events (" << ws.data(dsName.c_str())->numEntries() << " origDS events)" << " and " << data->sumEntries() << " wevents (" << ws.data(dsName.c_str())->sumEntries() << " origDS wevents)" << std::endl;
    }
    auto dataFit = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(data->reduce(RooFit::Cut(getString(ws, "cutDSFit").c_str()), RooFit::Name(dsSPlotNameFit.c_str()))));
    setDSParamaterRange(*dataFit, info);
    if (dataFit->numEntries()<=0) { std::cout << "[ERROR] No events from dataset " <<  dsSPlotNameFit << " passed kinematic cuts!" << std::endl; return false; }
    else if (dataFit->numEntries()>data->numEntries()) { std::cout << "[ERROR] Dataset " <<  dsSPlotNameFit << " has more events than original!" << std::endl; return false; }
    else if (dataFit->numEntries()==data->numEntries()) { info.Par.at("dsNameFit"+chg) = dsSPlotName; }
    else if (!isCompatibleDataset(*dataFit, *dynamic_cast<RooDataSet*>(ws.data(dsNameFit.c_str())))) { std::cout << "[ERROR] sPlot and original fit datasets are inconsistent!" << std::endl; return false; }
    else {
      if (ws.import(*dataFit, RooFit::Rename(dsSPlotNameFit.c_str()))) { std::cout << "[ERROR] RooDataSet " << dsSPlotNameFit << " was not imported!" << std::endl; return false; }
      else {
	if (ws.data(dsNameFit.c_str())) { ws.RecursiveRemove(ws.data(dsNameFit.c_str())); }
	info.Var.at("numEntries")[dsSPlotNameFit] = ws.data(dsSPlotNameFit.c_str())->sumEntries();
	info.Par.at("dsNameFit"+chg) = dsSPlotNameFit;
	defineSet(ws, "SET_"+dsSPlotNameFit, *dataFit->get());
	std::cout << "[INFO] Imported " << dsSPlotNameFit << " with " << dataFit->numEntries() << " events (" << ws.data(dsNameFit.c_str())->numEntries() << " origDS events)" << " with " << dataFit->sumEntries() << " wevents (" << ws.data(dsNameFit.c_str())->sumEntries() << " origDS wevents)" << std::endl;
      }
    }
  }
  //
  return true;
};


bool addSPlotWeight(RooWorkspace& ws, GlobalInfo& info, const std::string& chg, const std::string& dsLbl)
{
  const auto& dsName = info.Par.at(dsLbl+chg);
  if (dsLbl=="dsSPlotName" || dsName!=info.Par.at("dsSPlotName"+chg)) {
    const auto& ds = ws.data(dsName.c_str());
    if (!ds) { std::cout << "[ERROR] RooDataSet " << dsName << " was not found!" << endl; return false; }
    //
    // Create temporary dataset
    auto vars = *ds->get();
    const auto& wVar = RooRealVar("Weight", "Weight", -100000000.0, 100000000.0, "");
    if (!vars.find("Weight")) { vars.add(wVar); }
    const auto& weight = dynamic_cast<RooRealVar*>(vars.find("Weight"));
    auto tmpDS = std::unique_ptr<RooDataSet>(new RooDataSet(ds->GetName(), ds->GetTitle(), vars, RooFit::WeightVar(*weight)));
    //
    // Extract the sPlot weight variables
    StringSet_t objList;
    for (const auto& varN : info.StrS.at("fitCondVarName")) {
      for (const auto& obj : info.StrS.at("incObject_"+varN)) { objList.insert(obj); }
    }
    if (contain(objList, "DLenRes") && objList.size()==1) { objList.insert("JPsi"); }
    StringSet_t wVars;
    auto varIt = std::unique_ptr<TIterator>(vars.createIterator());
    for (auto itp = varIt->Next(); itp!=NULL; itp = varIt->Next()) {
      const auto& it = dynamic_cast<RooRealVar*>(itp); if (!it) continue;
      const std::string& name = it->GetName();
      for (const auto& obj : objList) {
	auto objI = obj.substr(0, obj.find("Cat"));
	if (name.rfind("N_"+objI, 0)==0 && name.rfind("_sw")!=std::string::npos) { wVars.insert(name); }
      }
    }
    if (wVars.empty()) { std::cout << "[ERROR] SPlot weight variables were not found in " << dsName << " !" << endl; return false; }
    std::cout << "[INFO] Using sPlot weight variables on " << dsName << " : "; for (const auto& wVar : wVars) { std::cout << wVar << " , "; }; std::cout << std::endl;
    //
    // Loop over the entries
    for (int i = 0; i < ds->numEntries(); i++) {
      ds->get(i);
      double w = 0.0;
      for (const auto& wVarN : wVars) { w += vars.getRealValue(wVarN.c_str(), 0.0); }
      weight->setVal(w * ds->weight());
      const auto& N = dynamic_cast<RooRealVar*>(vars.find("N_JPsiToMuMuOS_PA8Y16_sw"));
      tmpDS->add(vars, weight->getVal());
    }
    //
    // Import to RooWorkspace
    ws.RecursiveRemove(ds); if(ds) delete ds;
    if (ws.import(*tmpDS)) { std::cout << "[ERROR] addSPlotWeight: Failed to import " << tmpDS->GetName() << std::endl; }
  }
  //
  info.Var.at("numEntries")[dsName] = ws.data(dsName.c_str())->sumEntries();
  //
  return true;
};


bool importSPlotDataset(RooWorkspace& ws, GlobalInfo& info, const std::string& chg, const bool& setWeight=false)
{
  //
  std::cout << "[INFO] Importing sPlot datasets" << std::endl;
  //
  const auto& cha = info.Par.at("channel");
  for (const auto& col : info.StrS.at("fitSystem")) {
    if (chg!="" && !makeSPlotDS(ws, info, (cha+chg+"_"+col))) { return false; }
    if (setWeight && !addSPlotWeight(ws, info, chg, "dsSPlotName")) { return false; }
    if (setWeight && !addSPlotWeight(ws, info, chg, "dsSPlotNameFit")) { return false; }
  }
  return true;
};


bool loadFitResults(RooWorkspace& ws, const GlobalInfo& info, const std::string& chg)
{
  //
  const auto& cha = info.Par.at("channel");
  for (const auto& col : info.StrS.at("fitSystem")) {
    const auto& label = (cha+chg+"_"+col);
    // Load mass fit
    if (info.Flag.at("fitCand_DLenErr") && info.StrS.at("incObject_CandDLenErr").size()>1) {
      if (!loadPDFParameters(ws, "Cand_Mass", label, info, true, true)) { std::cout << "[ERROR] loadFitResults: Cand_Mass PDF was not loaded!" << endl; return false; }
    }
    // Load mass fit
    if (info.Flag.at("fitCand_DLen") && info.Flag.at("fitCand_Mass")) {
      if (!loadPDFParameters(ws, "Cand_Mass", label, info, true, true)) { std::cout << "[ERROR] loadFitResults: Cand_Mass PDF was not loaded!" << endl; return false; }
      //
      const std::string& pdfName = "pdfCandMass_TotToMuMuOS_PA8Y16";
      const auto& dsNameFit = info.Par.at("dsNameFit"+chg);
      std::vector<RooCmdArg> cmdList = { RooFit::Extended(true), RooFit::AsymptoticError(false), RooFit::InitialHesse(true), RooFit::Minos(false), RooFit::Strategy(2), RooFit::Minimizer("Minuit2"),
					 RooFit::Optimize(false), RooFit::NumCPU(32, 1), RooFit::Save(true), RooFit::Timer(true), RooFit::PrintLevel(1), RooFit::BatchMode(true)};
      std::cout << "[INFO] Fitting " << pdfName << " on " << dsNameFit << std::endl;
      std::unique_ptr<RooFitResult> fitResult;
      if (fitPDF(fitResult, ws, cmdList, pdfName, dsNameFit)) {
	for (const auto& vv : StringVector_t({"N_JPsi"+label, "R_Psi2S"+label, "N_Bkg"+label})) {
	  if (ws.var(vv.c_str())) {
	    ws.var(vv.c_str())->setConstant(true);
	    //ws.var(vv.c_str())->setMin(ws.var(vv.c_str())->getVal() - ws.var(vv.c_str())->getError()*10.0);
	    //ws.var(vv.c_str())->setMax(ws.var(vv.c_str())->getVal() + ws.var(vv.c_str())->getError()*10.0);
	  }
	}
      }
    }
    // Load decay length resolution fit
    if (info.Flag.at("fitCand_DLen")) {
      GlobalInfo infoTMP(info);
      const auto& o = *info.StrS.at("fitObject").begin();
      const auto& obj = info.Par.at("objTag_CandDLenRes_"+o);
      infoTMP.Par["objTag_CandDLenRes_"+o] = obj;
      const bool& setConst = !(info.Flag.at("fitMC") && info.Par.at("MC_CAT")=="PR");
      infoTMP.Par["ModelCandDLenRes_"+obj+label] = info.Par.at("ModelCandDLen_"+o+label+"_DLenRes")+"[DLenRes] + Delta["+obj+"CatPR]";
      if (!loadPDFParameters(ws, "Cand_DLenRes", label, infoTMP, false, setConst)) { std::cout << "[ERROR] loadFitResults: Cand_DLenRes PDF was not loaded!" << endl; return false; }
    }
    // Load background decay length fit
    if (info.Flag.at("fitCand_DLen") && !info.Flag.at("fitBkg")) {
      for (const auto& objN : info.StrS.at("incObject_CandDLen")) {
	if (objN=="BkgCatNoPR") {
	  GlobalInfo infoTMP(info);
	  infoTMP.StrS.at("fitObject") = {"Bkg"};
	  infoTMP.Par["objTag_CandDLen_Bkg"] = "Bkg";
	  const auto& o = *info.StrS.at("fitObject").begin();
	  infoTMP.Par["ModelCandDLen_Bkg"+label] = info.Par.at("ModelCandDLen_"+o+label+"_DLenRes")+"[DLenRes] + Delta[BkgCatPR] + "+info.Par.at("ModelCandDLen_"+o+label+"_"+objN)+"["+objN+"]";
	  if (!loadPDFParameters(ws, "Cand_DLen", label, infoTMP, false, true)) { std::cout << "[ERROR] loadFitResults: " << objN << " Cand_DLen PDF was not loaded!" << endl; return false; }
	}
      }
    }
    // Load signal truth decay length fit
    if (info.Flag.at("fitCand_DLen")) {
      for (const auto& objN : info.StrS.at("incObject_CandDLen")) {
	if (objN.rfind("Bkg", 0)!=0 && objN.rfind("CatNoPR")!=std::string::npos) {
	  GlobalInfo infoTMP(info);
	  infoTMP.Flag.at("fitMC") = true;
	  infoTMP.Par.at("DSTAG") = "MC_"+objN+"_"+info.Par.at("channelDS")+"_"+col;
	  const auto& objI = objN.substr(0, objN.find("Cat"));
	  infoTMP.StrS.at("fitObject") = {objI};
	  infoTMP.Par["objTag_CandDLenGen_"+objI] = objI;
	  const auto& o = *info.StrS.at("fitObject").begin();
	  infoTMP.Par["ModelCandDLenGen_"+objI+label] = "DeltaResolution[DLenRes] + "+info.Par.at("ModelCandDLen_"+o+label+"_"+objN)+"["+objN+"]";
	  const bool& setConst = true;//!info.Flag.at("fit"+objI);
	  if (!loadPDFParameters(ws, "Cand_DLenGen", label, infoTMP, false, setConst)) { std::cout << "[ERROR] loadFitResults: " << objN << " Cand_DLenGen PDF was not loaded!" << endl; return false; }
	  if (ws.var(("f_"+objN+label).c_str())) { ws.var(("f_"+objN+label).c_str())->setConstant(true); }
	  if (ws.var(("rLambdaSS21_"+objN+label).c_str())) { ws.var(("rLambdaSS21_"+objN+label).c_str())->setConstant(true); }
	}
      }
    }
  }
  //
  return true;
};


#endif // #ifndef Candidate_addModel_C
