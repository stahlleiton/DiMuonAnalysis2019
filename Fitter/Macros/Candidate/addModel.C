#ifndef Candidate_addModel_C
#define Candidate_addModel_C


#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooRealVar.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooStringVar.h"

#include <iostream>
#include <string>
#include <memory>

#include "../Utilities/initClasses.h"
#include "../Utilities/rooModelUtils.h"


bool addModel(RooWorkspace& ws, GlobalInfo& info, const StringDiMap_t& models, const std::string& chg, const std::string& varName)
{
  //
  const auto& cha         = info.Par.at("channel");
  const auto& varWindow   = (varName + "Window");
  const auto& varNormName = (varName+"Norm");
  std::string varType = varName; stringReplace(varType, "_", "");
  //
  for (const auto& col : info.StrS.at("fitSystem")) {
    RooArgList pdfListTot;
    for (const auto& mainObj : info.StrS.at("fitObject")) {;
      const auto& mainTag = mainObj + cha + chg;
      const auto& mainLabel = mainTag + "_" + col;
      std::cout << "[INFO] Implementing " << mainTag << " " << varName << " Model for " << col << std::endl;
      RooArgList pdfList;
      // Make sure that the main object is done first
      StringVector_t objV = { mainObj }; for (const auto& m : models.at("Model_"+mainLabel)) { if (m.first!=mainObj) { objV.push_back(m.first); } }
      for (const auto& obj : objV) {
	const auto& modelN = models.at("Model_"+mainLabel).at(obj);  
	const auto& tag = obj + cha + chg;
	const auto& label = tag + "_" + col;
	RooStringVar tmp; tmp.setVal(modelN.c_str()); tmp.SetTitle(("Model_"+label).c_str());
        ws.import(*dynamic_cast<TObject*>(&tmp), tmp.GetTitle()); // Save the model name for bookkeeping
        //
        const std::string& pdfName    = Form("pdf%s_%s",    varType.c_str(), label.c_str());
        const std::string& pdf1Name   = Form("pdf%s1_%s",   varType.c_str(), label.c_str());
        const std::string& pdf2Name   = Form("pdf%s2_%s",   varType.c_str(), label.c_str());
        const std::string& pdfPolName = Form("pdf%sPol_%s", varType.c_str(), label.c_str());
        const std::string& pdfTotName = Form("pdf%sTot_%s", varType.c_str(), label.c_str());
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
	      // check that all input parameters are defined
	      if (!( 
		    info.Par.count("N_"+label) &&
		    info.Par.count("Cut_"+label)
		     )) {
		std::cout << "[ERROR] Initial parameters where not found for " << tag << " CutAndCount Model in " << col << std::endl; return false;
	      }
	      // create the variables for this model
	      ws.factory( info.Par.at("N_"+label).c_str() );
	      std::string cut = info.Par.at("Cut_"+label);
	      RooStringVar tmp; tmp.setVal(cut.c_str()); tmp.SetTitle(("CutAndCount_"+label).c_str()); ws.import(*dynamic_cast<TObject*>(&tmp), tmp.GetTitle());
	      std::cout << "[INFO] " << tag << " in " << col << " added for CutAndCount!" << std::endl; break;
	    }
          case (int(Model::Template)):
	    {
	      // check that all input parameters are defined
	      if (!( 
		    info.Par.count("N_"+label)
		     )) {
		std::cout << "[ERROR] Initial parameters where not found for " << tag << " Template Model in " << col << std::endl; return false;
	      }
	      // create the Template
	      std::string dsName = ( "d" + chg + "_MC_" + obj + "_" + info.Par.at("channelDS") + "_" + col );
	      const std::vector< double > range = { double(ws.var(varName.c_str())->getBins(varWindow.c_str())) , ws.var(varName.c_str())->getMin(varWindow.c_str()) , ws.var(varName.c_str())->getMax(varWindow.c_str()) };
	      const auto& proceed = histToPdf(ws, pdfName, dsName, varName, range);
              if (info.Var.count("recoMCEntries")==0 || info.Var.at("recoMCEntries").count(label)==0) {
                if (proceed && (std::abs(info.Var.at("recoMCEntries").at(label)-ws.data(Form("dh%s_%s", varName.c_str(), label.c_str()))->sumEntries())>0.5)) {
                  std::cout << "[WARNING] The number of events in " << Form("dh%s_%s", varName.c_str(), label.c_str()) << " changed from (" << info.Var.at("recoMCEntries").at(label) << ") to (" <<
                    ws.data(Form("dh%s_%s", varName.c_str(), label.c_str()))->sumEntries() << ")" << std::endl;
                }
              }
	      //
	      if (proceed) {
                // create the variables for this model
                ws.factory( info.Par.at("N_"+label).c_str() );
                ws.factory(Form("RooExtendPdf::%s(%s,%s)", pdfTotName.c_str(),
                                pdfName.c_str(),
                                Form("N_%s", label.c_str())
                                ));
                ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
                pdfList.add( *ws.pdf(pdfTotName.c_str()) );
                std::cout << "[INFO] " << tag << " Template " << varName << " PDF in " << col << " added!" << std::endl; break;
              }
	      else {
		std::cout << "[INFO] %s Template " << varName << " PDF in " << col << " was NOT created!" << std::endl; break;
	      }
            }
            //-------------------------------------------
            //
            // Signal Candidate Mass Models
            //
            //-------------------------------------------
          case (int(Model::SingleGaussian)):
            {
      
              // check that all input parameters are defined
              if (!( 
		    info.Par.count("N_"+label) &&
                    info.Par.count("m_"+label) &&
                    info.Par.count("Sigma1_"+label)
                     )) { 
                std::cout << "[ERROR] Initial parameters where not found for " << tag << " Single Gaussian Model in " << col << std::endl; return false;
              }
              // create the variables for this model
              ws.factory( info.Par.at("N_"+label).c_str() );
              RooArgList pdfConstrains;
              StringVector_t varNames = {"m", "Sigma1"};
              for (const auto& v : varNames) {
		if (info.Par.count(v+"_"+label)) { ws.factory( info.Par.at(v+"_"+label).c_str() ); }
                // create the Gaussian PDFs for Constrain fits
                if (info.Par.count("val"+v+"_"+label) && info.Par.count("sig"+v+"_"+label)) {
                  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", (v+"_"+label).c_str(), (v+"_"+label).c_str(), info.Par.at("val"+v+"_"+label).c_str(), info.Par.at("sig"+v+"_"+label).c_str()));
                  pdfConstrains.add( *ws.pdf(("Constr"+v+"_"+label).c_str()) );
                }
              }
              if (pdfConstrains.getSize()>0) { ws.import(pdfConstrains, Form("pdfConstr%s", mainLabel.c_str())); }
	      //
	      // create the PDF
              ws.factory(Form("Gaussian::%s(%s, %s, %s)", pdfName.c_str(), varName.c_str(),
                              Form("m_%s", label.c_str()),
                              Form("Sigma1_%s", label.c_str())
                              ));
              ws.factory(Form("RooExtendPdf::%s(%s,%s)", pdfTotName.c_str(),
                              pdfName.c_str(),
                              Form("N_%s", label.c_str())
                              ));
	      ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
	      pdfList.add( *ws.pdf(pdfTotName.c_str()) );
              std::cout << "[INFO] " << tag << " Single Gaussian " << varName << " PDF in " << col << " added!" << std::endl; break;
            }
          case (int(Model::DoubleGaussian)):
            {
              // check that all input parameters are defined
              if (!( 
		    info.Par.count("N_"+label) &&
                    info.Par.count("m_"+label) &&
                    info.Par.count("Sigma1_"+label) &&
                    info.Par.count("Sigma2_"+label) &&
                    info.Par.count("f_"+label)
                     )) { 
                std::cout << "[ERROR] Initial parameters where not found for " << tag << " Double Gaussian Model in " << col << std::endl; return false;
              }
              // create the variables for this model
              ws.factory( info.Par.at("N_"+label).c_str() );
              RooArgList pdfConstrains;
              StringVector_t varNames = {"m", "Sigma1", "rSigma21", "Sigma2", "f"};
              for (const auto& v : varNames) {
		if (info.Par.count(v+"_"+label)) { ws.factory( info.Par.at(v+"_"+label).c_str() ); }
                // create the Gaussian PDFs for Constrain fits
                if (info.Par.count("val"+v+"_"+label) && info.Par.count("sig"+v+"_"+label)) {
                  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", (v+"_"+label).c_str(), (v+"_"+label).c_str(), info.Par.at("val"+v+"_"+label).c_str(), info.Par.at("sig"+v+"_"+label).c_str()));
                  pdfConstrains.add( *ws.pdf(("Constr"+v+"_"+label).c_str()) );
                }
              }
              if (pdfConstrains.getSize()>0) { ws.import(pdfConstrains, Form("pdfConstr%s", mainLabel.c_str())); }
              //
              // create the two PDFs             
              ws.factory(Form("Gaussian::%s(%s, %s, %s)", pdf1Name.c_str(), varName.c_str(), 
                              Form("m_%s", label.c_str()),
                              Form("Sigma1_%s", label.c_str())
                              ));
              ws.factory(Form("Gaussian::%s(%s, %s, %s)", pdf2Name.c_str(), varName.c_str(), 
                              Form("m_%s", label.c_str()),
                              Form("Sigma2_%s", label.c_str())
                              ));
              // Sum the PDFs to get the signal PDF
              ws.factory(Form("SUM::%s(%s*%s, %s)", pdfName.c_str(),
                              Form("f_%s", label.c_str()),
                              pdf1Name.c_str(),
                              pdf2Name.c_str()
                              ));
              ws.factory(Form("RooExtendPdf::%s(%s,%s)", pdfTotName.c_str(),
                              pdfName.c_str(),
                              Form("N_%s", label.c_str())
                              ));
	      ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
	      pdfList.add( *ws.pdf(pdfTotName.c_str()) );
              std::cout << "[INFO] " << tag << " Double Gaussian " << varName << " PDF in " << col << " added!" << std::endl; break;
            }
          case (int(Model::SingleCrystalBall)):
            {
              // check that all input parameters are defined
              if (!( 
		    info.Par.count("N_"+label) &&
                    info.Par.count("m_"+label) &&
                    info.Par.count("Sigma1_"+label) &&
                    info.Par.count("Alpha_"+label) &&
                    info.Par.count("n_"+label)
                     )) {
                std::cout << "[ERROR] Initial parameters where not found for " << tag << " Single Crystal Ball Model in " << col << std::endl; return false;
              }
              // create the variables for this model
              ws.factory( info.Par.at("N_"+label).c_str() );
              RooArgList pdfConstrains;
              StringVector_t varNames = {"m", "Sigma1", "Alpha", "n"};
              for (const auto& v : varNames) {
		if (info.Par.count(v+"_"+label)) { ws.factory( info.Par.at(v+"_"+label).c_str() ); }
                // create the Gaussian PDFs for Constrain fits
                if (info.Par.count("val"+v+"_"+label) && info.Par.count("sig"+v+"_"+label)) {
                  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", (v+"_"+label).c_str(), (v+"_"+label).c_str(), info.Par.at("val"+v+"_"+label).c_str(), info.Par.at("sig"+v+"_"+label).c_str()));
                  pdfConstrains.add( *ws.pdf(("Constr"+v+"_"+label).c_str()) );
                }
              }
              if (pdfConstrains.getSize()>0) { ws.import(pdfConstrains, Form("pdfConstr%s", mainLabel.c_str())); }
              //
              // create the PDF
              ws.factory(Form("CBShape::%s(%s, %s, %s, %s, %s)", pdfName.c_str(), varName.c_str(), 
                              Form("m_%s", label.c_str()),
                              Form("Sigma1_%s", label.c_str()),
                              Form("Alpha_%s", label.c_str()),
                              Form("n_%s", label.c_str())
                              ));
              ws.factory(Form("RooExtendPdf::%s(%s,%s)", pdfTotName.c_str(),
                              pdfName.c_str(),
                              Form("N_%s", label.c_str())
                              ));
	      ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
	      pdfList.add( *ws.pdf(pdfTotName.c_str()) );
              std::cout << "[INFO] " << tag << " Single Crystal Ball " << varName << " PDF in " << col << " included" << std::endl; break;
            }
          case (int(Model::DoubleCrystalBall)):
            {
              // check that all input parameters are defined
              if (!(
		    info.Par.count("N_"+label) &&
                    info.Par.count("m_"+label) &&
                    info.Par.count("Sigma1_"+label) &&
                    info.Par.count("Sigma2_"+label) &&
                    info.Par.count("Alpha_"+label) &&
                    info.Par.count("Alpha2_"+label) &&
                    info.Par.count("n_"+label) &&
                    info.Par.count("n2_"+label) &&
                    info.Par.count("f_"+label)
                    )) {
                std::cout << "[ERROR] Initial parameters where not found for " << tag << " Double Crystal Ball Model in " << col << std::endl; return false;
              }
              // create the variables for this model
              ws.factory( info.Par.at("N_"+label).c_str() );
              RooArgList pdfConstrains;
              StringVector_t varNames = {"m", "Sigma1", "rSigma21", "Sigma2", "Alpha", "Alpha2", "n", "n2", "f"};
              for (const auto& v : varNames) {
		if (info.Par.count(v+"_"+label)) { ws.factory( info.Par.at(v+"_"+label).c_str() ); }
                // create the Gaussian PDFs for Constrain fits
                if (info.Par.count("val"+v+"_"+label) && info.Par.count("sig"+v+"_"+label)) {
                  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", (v+"_"+label).c_str(), (v+"_"+label).c_str(), info.Par.at("val"+v+"_"+label).c_str(), info.Par.at("sig"+v+"_"+label).c_str()));
                  pdfConstrains.add( *ws.pdf(("Constr"+v+"_"+label).c_str()) );
                }
              }
              if (pdfConstrains.getSize()>0) { ws.import(pdfConstrains, Form("pdfConstr%s", mainLabel.c_str())); }
              //
              // create the two PDFs
              ws.factory(Form("CBShape::%s(%s, %s, %s, %s, %s)", pdf1Name.c_str(), varName.c_str(), 
                              Form("m_%s", label.c_str()),
                              Form("Sigma1_%s", label.c_str()),
                              Form("Alpha_%s", label.c_str()),
                              Form("n_%s", label.c_str())
                              ));
              ws.factory(Form("CBShape::%s(%s, %s, %s, %s, %s)", pdf2Name.c_str(), varName.c_str(), 
                              Form("m_%s", label.c_str()), 
                              Form("Sigma2_%s", label.c_str()),
                              Form("Alpha2_%s", label.c_str()),
                              Form("n2_%s", label.c_str())
                              ));
              // Sum the PDFs to get the signal PDF
              ws.factory(Form("SUM::%s(%s*%s, %s)", pdfName.c_str(),
                              Form("f_%s", label.c_str()),
                              pdf1Name.c_str(),
                              pdf2Name.c_str()
                              ));
              ws.factory(Form("RooExtendPdf::%s(%s,%s)", pdfTotName.c_str(),
                              pdfName.c_str(),
                              Form("N_%s", label.c_str())
                              ));
	      ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
	      pdfList.add( *ws.pdf(pdfTotName.c_str()) );
              std::cout << "[INFO] " << tag << " Double Crystal Ball " << varName << " PDF in " << col << " added!" << std::endl; break;
            }
          case (int(Model::GaussianAndCrystalBall)):
            {
              // check that all input parameters are defined
              if (!(
		    info.Par.count("N_"+label) &&
                    info.Par.count("m_"+label) &&
                    info.Par.count("Sigma1_"+label) &&
                    info.Par.count("Sigma2_"+label) &&
                    info.Par.count("Alpha_"+label) &&
                    info.Par.count("n_"+label) &&
                    info.Par.count("f_"+label)
                    )) {
                std::cout << "[ERROR] Initial parameters where not found for " << tag << " Gaussian and Crystal Ball Model in " << col << std::endl; return false;
              }
              // create the variables for this model
              ws.factory( info.Par.at("N_"+label).c_str() );
              RooArgList pdfConstrains;
              StringVector_t varNames = {"m", "Sigma1", "rSigma21", "Sigma2", "Alpha", "n", "f"};
              for (const auto& v : varNames) {
		if (info.Par.count(v+"_"+label)) { ws.factory( info.Par.at(v+"_"+label).c_str() ); }
                // create the Gaussian PDFs for Constrain fits
                if (info.Par.count("val"+v+"_"+label) && info.Par.count("sig"+v+"_"+label)) {
                  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", (v+"_"+label).c_str(), (v+"_"+label).c_str(), info.Par.at("val"+v+"_"+label).c_str(), info.Par.at("sig"+v+"_"+label).c_str()));
                  pdfConstrains.add( *ws.pdf(("Constr"+v+"_"+label).c_str()) );
                }
              }
              if (pdfConstrains.getSize()>0) { ws.import(pdfConstrains, Form("pdfConstr%s", mainLabel.c_str())); }
              //
              // create the two PDFs
              ws.factory(Form("CBShape::%s(%s, %s, %s, %s, %s)", pdf1Name.c_str(), varName.c_str(),
                              Form("m_%s", label.c_str()),
                              Form("Sigma1_%s", label.c_str()),
                              Form("Alpha_%s", label.c_str()),
                              Form("n_%s", label.c_str())
                              ));
              ws.factory(Form("Gaussian::%s(%s, %s, %s)", pdf2Name.c_str(), varName.c_str(),
                              Form("m_%s", label.c_str()), 
                              Form("Sigma2_%s", label.c_str())
                              ));
              // Sum the PDFs to get the signal PDF 
              ws.factory(Form("SUM::%s(%s*%s, %s)", pdfName.c_str(),
                              Form("f_%s", label.c_str()),
                              pdf1Name.c_str(),
                              pdf2Name.c_str()
                              ));
              ws.factory(Form("RooExtendPdf::%s(%s,%s)", pdfTotName.c_str(),
                              pdfName.c_str(),
                              Form("N_%s", label.c_str())
                              ));
	      ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
	      pdfList.add( *ws.pdf(pdfTotName.c_str()) );
              std::cout << "[INFO] " << tag << " Gaussian and Crystal Ball " << varName << " PDF in " << col << " added!" << std::endl; break;
            }
	    case (int(Model::Voigtian)):
            {
              // check that all input parameters are defined
              if (!(
		    info.Par.count("N_"+label) &&
                    info.Par.count("m_"+label) &&
                    info.Par.count("Sigma1_"+label) &&
                    info.Par.count("Sigma2_"+label)
                    )) {
                std::cout << "[ERROR] Initial parameters where not found for " << tag << " Voigtian Model in " << col << std::endl; return false;
              }
              // create the variables for this model
              ws.factory( info.Par.at("N_"+label).c_str() );
              RooArgList pdfConstrains;
              StringVector_t varNames = {"m", "Sigma1", "rSigma21", "Sigma2"};
              for (const auto& v : varNames) {
		if (info.Par.count(v+"_"+label)) { ws.factory( info.Par.at(v+"_"+label).c_str() ); }
                // create the Gaussian PDFs for Constrain fits
                if (info.Par.count("val"+v+"_"+label) && info.Par.count("sig"+v+"_"+label)) {
                  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", (v+"_"+label).c_str(), (v+"_"+label).c_str(), info.Par.at("val"+v+"_"+label).c_str(), info.Par.at("sig"+v+"_"+label).c_str()));
                  pdfConstrains.add( *ws.pdf(("Constr"+v+"_"+label).c_str()) );
                }
              }
              if (pdfConstrains.getSize()>0) { ws.import(pdfConstrains, Form("pdfConstr%s", mainLabel.c_str())); }
              //
              // create the PDF
              ws.factory(Form("Voigtian::%s(%s, %s, %s, %s)", pdfName.c_str(), varName.c_str(),
                              Form("m_%s", label.c_str()),
                              Form("Sigma1_%s", label.c_str()),
                              Form("Sigma2_%s", label.c_str())
                              ));
              ws.factory(Form("RooExtendPdf::%s(%s,%s)", pdfTotName.c_str(),
                              pdfName.c_str(),
                              Form("N_%s", label.c_str())
                              ));
	      ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
	      pdfList.add( *ws.pdf(pdfTotName.c_str()) );
              std::cout << "[INFO] " << tag << " Voigtian " << varName << " PDF in " << col << " added!" << std::endl; break;
            }
          case (int(Model::BWCrystalBall)):
            {
              // check that all input parameters are defined
              if (!(
		    info.Par.count("N_"+label) &&
                    info.Par.count("m_"+label) &&
                    info.Par.count("Sigma1_"+label) &&
                    info.Par.count("Sigma2_"+label) &&
                    info.Par.count("Alpha_"+label) &&
                    info.Par.count("n_"+label)
                    )) {
                std::cout << "[ERROR] Initial parameters where not found for " << tag << " Breit-Wigner-CrystalBall Model in " << col << std::endl; return false;
              }
              // create the variables for this model
              ws.factory( info.Par.at("N_"+label).c_str() );
              RooArgList pdfConstrains;
              StringVector_t varNames = {"m", "Sigma1", "rSigma21", "Sigma2", "Alpha", "n"};
              for (const auto& v : varNames) {
		if (info.Par.count(v+"_"+label)) { ws.factory( info.Par.at(v+"_"+label).c_str() ); }
                // create the Gaussian PDFs for Constrain fits
                if (info.Par.count("val"+v+"_"+label) && info.Par.count("sig"+v+"_"+label)) {
                  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", (v+"_"+label).c_str(), (v+"_"+label).c_str(), info.Par.at("val"+v+"_"+label).c_str(), info.Par.at("sig"+v+"_"+label).c_str()));
                  pdfConstrains.add( *ws.pdf(("Constr"+v+"_"+label).c_str()) );
                }
              }
              if (pdfConstrains.getSize()>0) { ws.import(pdfConstrains, Form("pdfConstr%s", mainLabel.c_str())); }
              //
              // create the two PDFs
	       ws.factory(Form("BreitWigner::%s(%s, %s, %s)", pdf1Name.c_str(), varName.c_str(),
                              Form("m_%s", label.c_str()), 
                              Form("Sigma1_%s", label.c_str())
                              ));
              ws.factory(Form("CBShape::%s(%s, %s, %s, %s, %s)", pdf2Name.c_str(), varName.c_str(),
			      Form("Peak_%s[0]", label.c_str()),
                              Form("Sigma2_%s", label.c_str()),
                              Form("Alpha_%s", label.c_str()),
                              Form("n_%s", label.c_str())
                              ));
              // covolve the PDFs to get the signal PDF 
              ws.factory(Form("FCONV::%s(%s, %s, %s)", pdfName.c_str(), varName.c_str(),
                              pdf1Name.c_str(),
                              pdf2Name.c_str()
                              ));
              ws.factory(Form("RooExtendPdf::%s(%s,%s)", pdfTotName.c_str(),
                              pdfName.c_str(),
                              Form("N_%s", label.c_str())
                              ));
	      ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
	      pdfList.add( *ws.pdf(pdfTotName.c_str()) );
              std::cout << "[INFO] " << tag << " Breit-Wigner-CrystallBall " << varName << " PDF in " << col << " added!" << std::endl; break;
            }
            //-------------------------------------------
            //
            // Background Candidate Mass Models
            //
            //-------------------------------------------
          case (int(Model::Uniform)):
            {
              // check that all input parameters are defined
              if (!info.Par.count("N_"+label)) {
                std::cout << "[ERROR] Initial parameters where not found for " << tag << " Uniform Model in " << col << std::endl; return false;
              }
	      // create the variables for this model
	      ws.factory( info.Par.at("N_"+label).c_str() );
	      //
              // create the PDF
              ws.factory(Form("Uniform::%s(%s)", pdfName.c_str(), varName.c_str()));
              //
              ws.factory(Form("RooExtendPdf::%s(%s,%s)", pdfTotName.c_str(),
                              pdfName.c_str(),
                              Form("N_%s", label.c_str())
                              ));
	      ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
	      pdfList.add( *ws.pdf(pdfTotName.c_str()) );
              std::cout << "[INFO] %s Uniform " << varName << " PDF in " << col << " added!" << std::endl; break;
            }
          case (int(Model::Chebychev1)):
            {
              // check that all input parameters are defined 
              if (!(
		    info.Par.count("N_"+label) &&
                    info.Par.count("Lambda1_"+label)
                    )) {
                std::cout << "[ERROR] Initial parameters where not found for " << tag << " First Order Chebychev in " << col << std::endl; return false;
              }
              // create the variables for this model
              ws.factory( info.Par.at("N_"+label).c_str() );
              RooArgList pdfConstrains;
	      StringVector_t varNames = {"Lambda1"};
	      for (const auto& v : varNames) {
		if (info.Par.count(v+"_"+label)) { ws.factory( info.Par.at(v+"_"+label).c_str() ); }
		// create the Gaussian PDFs for Constrain fits
		if (info.Par.count("val"+v+"_"+label) && info.Par.count("sig"+v+"_"+label)) {
		  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", (v+"_"+label).c_str(), (v+"_"+label).c_str(), info.Par.at("val"+v+"_"+label).c_str(), info.Par.at("sig"+v+"_"+label).c_str()));
		  pdfConstrains.add( *ws.pdf(("Constr"+v+"_"+label).c_str()) );
		}
	      }
	      if (pdfConstrains.getSize()>0) { ws.import(pdfConstrains, Form("pdfConstr%s", mainLabel.c_str())); }
	      //
              // create the PDF
              ws.factory(Form("Chebychev::%s(%s, {%s})", pdfName.c_str(), varName.c_str(), 
                              Form("Lambda1_%s", label.c_str())
                              ));
              ws.factory(Form("RooExtendPdf::%s(%s,%s)", pdfTotName.c_str(),
                              pdfName.c_str(),
                              Form("N_%s", label.c_str())
                              ));
	      ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
	      pdfList.add( *ws.pdf(pdfTotName.c_str()) );
              std::cout << "[INFO] " << tag << " 1st Order Chebychev " << varName << " PDF in " << col << " added!" << std::endl; break;
            }
          case (int(Model::Chebychev2)):
            {
              // check that all input parameters are defined
              if (!(
		    info.Par.count("N_"+label) &&
                    info.Par.count("Lambda1_"+label) &&
                    info.Par.count("Lambda2_"+label)
                    )) {
                std::cout << "[ERROR] Initial parameters where not found for " << tag << " Second Order Chebychev in " << col << std::endl; return false;
              }
	      // create the variables for this model
	      ws.factory( info.Par.at("N_"+label).c_str() );
	      RooArgList pdfConstrains;
	      StringVector_t varNames = {"Lambda1", "Lambda2"};
	      for (const auto& v : varNames) {
		if (info.Par.count(v+"_"+label)) { ws.factory( info.Par.at(v+"_"+label).c_str() ); }
		// create the Gaussian PDFs for Constrain fits
		if (info.Par.count("val"+v+"_"+label) && info.Par.count("sig"+v+"_"+label)) {
		  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", (v+"_"+label).c_str(), (v+"_"+label).c_str(), info.Par.at("val"+v+"_"+label).c_str(), info.Par.at("sig"+v+"_"+label).c_str()));
		  pdfConstrains.add( *ws.pdf(("Constr"+v+"_"+label).c_str()) );
		}
	      }
	      if (pdfConstrains.getSize()>0) { ws.import(pdfConstrains, Form("pdfConstr%s", mainLabel.c_str())); }
	      //
              // create the PDF
              ws.factory(Form("Chebychev::%s(%s, {%s, %s})", pdfName.c_str(), varName.c_str(), 
                              Form("Lambda1_%s", label.c_str()),
                              Form("Lambda2_%s", label.c_str())
                              ));
              ws.factory(Form("RooExtendPdf::%s(%s,%s)", pdfTotName.c_str(),
                              pdfName.c_str(),
                              Form("N_%s", label.c_str())
                              ));
	      ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
	      pdfList.add( *ws.pdf(pdfTotName.c_str()) );
              std::cout << "[INFO] " << tag << " 2nd Order Chebychev " << varName << " PDF in " << col << " added!" << std::endl; break;
            }
          case (int(Model::Chebychev3)):
            {
              // check that all input parameters are defined 
              if (!(
		    info.Par.count("N_"+label) &&
                    info.Par.count("Lambda1_"+label) &&
                    info.Par.count("Lambda2_"+label) &&
                    info.Par.count("Lambda3_"+label)
                    )) {
                std::cout << "[ERROR] Initial parameters where not found for " << tag << " Third Order Chebychev in " << col << std::endl; return false;
              }
	      // create the variables for this model
	      ws.factory( info.Par.at("N_"+label).c_str() );
	      RooArgList pdfConstrains;
	      StringVector_t varNames = {"Lambda1", "Lambda2", "Lambda3"};
	      for (const auto& v : varNames) {
		if (info.Par.count(v+"_"+label)) { ws.factory( info.Par.at(v+"_"+label).c_str() ); }
		// create the Gaussian PDFs for Constrain fits
		if (info.Par.count("val"+v+"_"+label) && info.Par.count("sig"+v+"_"+label)) {
		  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", (v+"_"+label).c_str(), (v+"_"+label).c_str(), info.Par.at("val"+v+"_"+label).c_str(), info.Par.at("sig"+v+"_"+label).c_str()));
		  pdfConstrains.add( *ws.pdf(("Constr"+v+"_"+label).c_str()) );
		}
	      }
	      if (pdfConstrains.getSize()>0) { ws.import(pdfConstrains, Form("pdfConstr%s", mainLabel.c_str())); }
	      //
              // create the PDF
              ws.factory(Form("Chebychev::%s(%s, {%s, %s, %s})", pdfName.c_str(), varName.c_str(), 
                              Form("Lambda1_%s", label.c_str()),
                              Form("Lambda2_%s", label.c_str()),
                              Form("Lambda3_%s", label.c_str())
                              ));
              ws.factory(Form("RooExtendPdf::%s(%s,%s)", pdfTotName.c_str(),
                              pdfName.c_str(),
                              Form("N_%s", label.c_str())
                              ));
	      ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
	      pdfList.add( *ws.pdf(pdfTotName.c_str()) );
              std::cout << "[INFO] " << tag << " 3rd Order Chebychev " << varName << " PDF in " << col << " added!" << std::endl; break;
            }
          case (int(Model::Chebychev4)):
            {
              // check that all input parameters are defined 
              if (!(
		    info.Par.count("N_"+label) &&
                    info.Par.count("Lambda1_"+label) &&
                    info.Par.count("Lambda2_"+label) &&
                    info.Par.count("Lambda3_"+label) &&
                    info.Par.count("Lambda4_"+label)
                    )) {
                std::cout << "[ERROR] Initial parameters where not found for " << tag << " Fourth Order Chebychev in " << col << std::endl; return false;
              }
	      // create the variables for this model
	      ws.factory( info.Par.at("N_"+label).c_str() );
	      RooArgList pdfConstrains;
	      StringVector_t varNames = {"Lambda1", "Lambda2", "Lambda3", "Lambda4"};
	      for (const auto& v : varNames) {
		if (info.Par.count(v+"_"+label)) { ws.factory( info.Par.at(v+"_"+label).c_str() ); }
		// create the Gaussian PDFs for Constrain fits
		if (info.Par.count("val"+v+"_"+label) && info.Par.count("sig"+v+"_"+label)) {
		  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", (v+"_"+label).c_str(), (v+"_"+label).c_str(), info.Par.at("val"+v+"_"+label).c_str(), info.Par.at("sig"+v+"_"+label).c_str()));
		  pdfConstrains.add( *ws.pdf(("Constr"+v+"_"+label).c_str()) );
		}
	      }
	      if (pdfConstrains.getSize()>0) { ws.import(pdfConstrains, Form("pdfConstr%s", mainLabel.c_str())); }
	      //
              // create the PDF
              ws.factory(Form("Chebychev::%s(%s, {%s, %s, %s, %s})", pdfName.c_str(), varName.c_str(), 
                              Form("Lambda1_%s", label.c_str()),
                              Form("Lambda2_%s", label.c_str()),
                              Form("Lambda3_%s", label.c_str()),
                              Form("Lambda4_%s", label.c_str())
                              ));
              ws.factory(Form("RooExtendPdf::%s(%s,%s)", pdfTotName.c_str(),
                              pdfName.c_str(),
                              Form("N_%s", label.c_str())
                              ));
	      ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
	      pdfList.add( *ws.pdf(pdfTotName.c_str()) );
              std::cout << "[INFO] " << tag << " 4th Order Chebychev " << varName << " PDF in " << col << " added!" << std::endl; break;
            }
          case (int(Model::Chebychev5)):
            {
              // check that all input parameters are defined 
              if (!(
		    info.Par.count("N_"+label) &&
                    info.Par.count("Lambda1_"+label) &&
                    info.Par.count("Lambda2_"+label) &&
                    info.Par.count("Lambda3_"+label) &&
                    info.Par.count("Lambda4_"+label) &&
                    info.Par.count("Lambda5_"+label)
                    )) {
                std::cout << "[ERROR] Initial parameters where not found for " << tag << " Fifth Order Chebychev in " << col << std::endl; return false;
              }
	      // create the variables for this model
	      ws.factory( info.Par.at("N_"+label).c_str() );
	      RooArgList pdfConstrains;
	      StringVector_t varNames = {"Lambda1", "Lambda2", "Lambda3", "Lambda4", "Lambda5"};
	      for (const auto& v : varNames) {
		if (info.Par.count(v+"_"+label)) { ws.factory( info.Par.at(v+"_"+label).c_str() ); }
		// create the Gaussian PDFs for Constrain fits
		if (info.Par.count("val"+v+"_"+label) && info.Par.count("sig"+v+"_"+label)) {
		  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", (v+"_"+label).c_str(), (v+"_"+label).c_str(), info.Par.at("val"+v+"_"+label).c_str(), info.Par.at("sig"+v+"_"+label).c_str()));
		  pdfConstrains.add( *ws.pdf(("Constr"+v+"_"+label).c_str()) );
		}
	      }
	      if (pdfConstrains.getSize()>0) { ws.import(pdfConstrains, Form("pdfConstr%s", mainLabel.c_str())); }
	      //
              // create the PDF
              ws.factory(Form("Chebychev::%s(%s, {%s, %s, %s, %s, %s})", pdfName.c_str(), varName.c_str(), 
                              Form("Lambda1_%s", label.c_str()),
                              Form("Lambda2_%s", label.c_str()),
                              Form("Lambda3_%s", label.c_str()),
                              Form("Lambda4_%s", label.c_str()),
                              Form("Lambda5_%s", label.c_str())
                              ));
              ws.factory(Form("RooExtendPdf::%s(%s,%s)", pdfTotName.c_str(),
                              pdfName.c_str(),
                              Form("N_%s", label.c_str())
                              ));
	      ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
	      pdfList.add( *ws.pdf(pdfTotName.c_str()) );
              std::cout << "[INFO] " << tag << " 5th Order Chebychev " << varName << " PDF in " << col << " added!" << std::endl; break;
            }
          case (int(Model::Chebychev6)):
            {
              // check that all input parameters are defined 
              if (!(
		    info.Par.count("N_"+label) &&
                    info.Par.count("Lambda1_"+label) &&
                    info.Par.count("Lambda2_"+label) &&
                    info.Par.count("Lambda3_"+label) &&
                    info.Par.count("Lambda4_"+label) &&
                    info.Par.count("Lambda5_"+label) &&
                    info.Par.count("Lambda6_"+label)
                    )) {
                std::cout << "[ERROR] Initial parameters where not found for " << tag << " Sixth Order Chebychev in " << col << std::endl; return false;
              }
	      // create the variables for this model
	      ws.factory( info.Par.at("N_"+label).c_str() );
	      RooArgList pdfConstrains;
	      StringVector_t varNames = {"Lambda1", "Lambda2", "Lambda3", "Lambda4", "Lambda5", "Lambda6"};
	      for (const auto& v : varNames) {
		if (info.Par.count(v+"_"+label)) { ws.factory( info.Par.at(v+"_"+label).c_str() ); }
		// create the Gaussian PDFs for Constrain fits
		if (info.Par.count("val"+v+"_"+label) && info.Par.count("sig"+v+"_"+label)) {
		  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", (v+"_"+label).c_str(), (v+"_"+label).c_str(), info.Par.at("val"+v+"_"+label).c_str(), info.Par.at("sig"+v+"_"+label).c_str()));
		  pdfConstrains.add( *ws.pdf(("Constr"+v+"_"+label).c_str()) );
		}
	      }
	      if (pdfConstrains.getSize()>0) { ws.import(pdfConstrains, Form("pdfConstr%s", mainLabel.c_str())); }
	      //
              // create the PDF
              ws.factory(Form("Chebychev::%s(%s, {%s, %s, %s, %s, %s, %s})", pdfName.c_str(), varName.c_str(), 
                              Form("Lambda1_%s", label.c_str()),
                              Form("Lambda2_%s", label.c_str()),
                              Form("Lambda3_%s", label.c_str()),
                              Form("Lambda4_%s", label.c_str()),
                              Form("Lambda5_%s", label.c_str()),
                              Form("Lambda6_%s", label.c_str())
                              ));
              ws.factory(Form("RooExtendPdf::%s(%s,%s)", pdfTotName.c_str(),
                              pdfName.c_str(),
                              Form("N_%s", label.c_str())
                              ));
	      ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
	      pdfList.add( *ws.pdf(pdfTotName.c_str()) );
              std::cout << "[INFO] " << tag << " 6th Order Chebychev " << varName << " PDF in " << col << " added!" << std::endl; break;
            }
          case (int(Model::ExpChebychev1)):
            {
              // check that all input parameters are defined 
              if (!( 
		    info.Par.count("N_"+label) &&
                    info.Par.count("Lambda1_"+label)
                     )) { 
                std::cout << "[ERROR] Initial parameters where not found for " << tag << " First Order Exponential Chebychev " << col << std::endl; return false;
              }
	      // create the variables for this model
	      ws.factory( info.Par.at("N_"+label).c_str() );
	      RooArgList pdfConstrains;
	      StringVector_t varNames = {"Lambda1"};
	      for (const auto& v : varNames) {
		if (info.Par.count(v+"_"+label)) { ws.factory( info.Par.at(v+"_"+label).c_str() ); }
		// create the Gaussian PDFs for Constrain fits
		if (info.Par.count("val"+v+"_"+label) && info.Par.count("sig"+v+"_"+label)) {
		  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", (v+"_"+label).c_str(), (v+"_"+label).c_str(), info.Par.at("val"+v+"_"+label).c_str(), info.Par.at("sig"+v+"_"+label).c_str()));
		  pdfConstrains.add( *ws.pdf(("Constr"+v+"_"+label).c_str()) );
		}
	      }
	      if (pdfConstrains.getSize()>0) { ws.import(pdfConstrains, Form("pdfConstr%s", mainLabel.c_str())); }
	      //
              // create the PDF
              ws.factory(Form("RooFormulaVar::%s('( 2.0*@0 - @2 - @1 )/( @2 - @1 )', {%s, vMin[%.6f], vMax[%.6f]})",
                              varNormName.c_str(), varName.c_str(),
                              info.Var.at(varName).at("Min"),
                              info.Var.at(varName).at("Max")
                              ));
              ws.factory(Form("RooFormulaVar::%s('@1*(@0) + 1.0', {%s, %s})", 
                              pdfPolName.c_str(), varNormName.c_str(), 
                              Form("Lambda1_%s", label.c_str())
                              ));
              ws.factory(Form("Exponential::%s(%s, One[1.0])", pdfName.c_str(),
                              pdfPolName.c_str()
                              ));
              ws.factory(Form("RooExtendPdf::%s(%s,%s)", pdfTotName.c_str(),
                              pdfName.c_str(),
                              Form("N_%s", label.c_str())
                              ));
	      ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
	      pdfList.add( *ws.pdf(pdfTotName.c_str()) );
              std::cout << "[INFO] " << tag << " 1st Order Exponential Chebychev " << varName << " PDF in " << col << " added!" << std::endl; break;
            }
          case (int(Model::ExpChebychev2)):
            {
              // check that all input parameters are defined 
              if (!(
		    info.Par.count("N_"+label) &&
                    info.Par.count("Lambda1_"+label) &&
                    info.Par.count("Lambda2_"+label)
                    )) {
                std::cout << "[ERROR] Initial parameters where not found for " << tag << " Second Order Exponential Chebychev in " << col << std::endl; return false;
              }
	      // create the variables for this model
	      ws.factory( info.Par.at("N_"+label).c_str() );
	      RooArgList pdfConstrains;
	      StringVector_t varNames = {"Lambda1", "Lambda2"};
	      for (const auto& v : varNames) {
		if (info.Par.count(v+"_"+label)) { ws.factory( info.Par.at(v+"_"+label).c_str() ); }
		// create the Gaussian PDFs for Constrain fits
		if (info.Par.count("val"+v+"_"+label) && info.Par.count("sig"+v+"_"+label)) {
		  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", (v+"_"+label).c_str(), (v+"_"+label).c_str(), info.Par.at("val"+v+"_"+label).c_str(), info.Par.at("sig"+v+"_"+label).c_str()));
		  pdfConstrains.add( *ws.pdf(("Constr"+v+"_"+label).c_str()) );
		}
	      }
	      if (pdfConstrains.getSize()>0) { ws.import(pdfConstrains, Form("pdfConstr%s", mainLabel.c_str())); }
	      //
              // create the PDF
              ws.factory(Form("RooFormulaVar::%s('( 2.0*@0 - @2 - @1 )/( @2 - @1 )', {%s, vMin[%.6f], vMax[%.6f]})",
                              varNormName.c_str(), varName.c_str(),
                              info.Var.at(varName).at("Min"),
                              info.Var.at(varName).at("Max")
                              ));            
              ws.factory(Form("RooFormulaVar::%s('@2*(2.0*@0*@0 - 1.0) + @1*(@0) + 1.0', {%s, %s, %s})", 
                              pdfPolName.c_str(), varNormName.c_str(),
                              Form("Lambda1_%s", label.c_str()),
                              Form("Lambda2_%s", label.c_str())
                              ));
              ws.factory(Form("Exponential::%s(%s, One[1.0])", pdfName.c_str(),
                              pdfPolName.c_str()
                              ));
              ws.factory(Form("RooExtendPdf::%s(%s,%s)", pdfTotName.c_str(),
                              pdfName.c_str(),
                              Form("N_%s", label.c_str())
                              ));
	      ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
	      pdfList.add( *ws.pdf(pdfTotName.c_str()) );
              std::cout << "[INFO] " << tag << " 2nd Order Exponential Chebychev " << varName << " PDF in " << col << " added!" << std::endl; break;
            }
          case (int(Model::ExpChebychev3)):
            {
              // check that all input parameters are defined 
              if (!(
		    info.Par.count("N_"+label) &&
                    info.Par.count("Lambda1_"+label) &&
                    info.Par.count("Lambda2_"+label) &&
                    info.Par.count("Lambda3_"+label)
                    )) {
                std::cout << "[ERROR] Initial parameters where not found for " << tag << " Third Order Exponential Chebychev in " << col << std::endl; return false;
              }
	      // create the variables for this model
	      ws.factory( info.Par.at("N_"+label).c_str() );
	      RooArgList pdfConstrains;
	      StringVector_t varNames = {"Lambda1", "Lambda2", "Lambda3"};
	      for (const auto& v : varNames) {
		if (info.Par.count(v+"_"+label)) { ws.factory( info.Par.at(v+"_"+label).c_str() ); }
		// create the Gaussian PDFs for Constrain fits
		if (info.Par.count("val"+v+"_"+label) && info.Par.count("sig"+v+"_"+label)) {
		  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", (v+"_"+label).c_str(), (v+"_"+label).c_str(), info.Par.at("val"+v+"_"+label).c_str(), info.Par.at("sig"+v+"_"+label).c_str()));
		  pdfConstrains.add( *ws.pdf(("Constr"+v+"_"+label).c_str()) );
		}
	      }
	      if (pdfConstrains.getSize()>0) { ws.import(pdfConstrains, Form("pdfConstr%s", mainLabel.c_str())); }
	      //
              // create the PDF
              ws.factory(Form("RooFormulaVar::%s('( 2.0*@0 - @2 - @1 )/( @2 - @1 )', {%s, vMin[%.6f], vMax[%.6f]})",
                              varNormName.c_str(), varName.c_str(),
                              info.Var.at(varName).at("Min"),
                              info.Var.at(varName).at("Max")
                              ));
              ws.factory(Form("RooFormulaVar::%s('@3*(4.0*@0*@0*@0 - 3.0*@0) + @2*(2.0*@0*@0 - 1.0) + @1*(@0) + 1.0', {%s, %s, %s, %s})", 
                              pdfPolName.c_str(), varNormName.c_str(),
                              Form("Lambda1_%s", label.c_str()),
                              Form("Lambda2_%s", label.c_str()),
                              Form("Lambda3_%s", label.c_str())
                              ));
              ws.factory(Form("Exponential::%s(%s, One[1.0])", pdfName.c_str(),
                              pdfPolName.c_str()
                              ));
              ws.factory(Form("RooExtendPdf::%s(%s,%s)", pdfTotName.c_str(),
                              pdfName.c_str(),
                              Form("N_%s", label.c_str())
                              ));
	      ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
	      pdfList.add( *ws.pdf(pdfTotName.c_str()) );
              std::cout << "[INFO] " << tag << " 3rd Order Exponential Chebychev " << varName << " PDF in " << col << " added!" << std::endl; break;
            }
          case (int(Model::ExpChebychev4)):
            {
              // check that all input parameters are defined 
              if (!(
		    info.Par.count("N_"+label) &&
                    info.Par.count("Lambda1_"+label) &&
                    info.Par.count("Lambda2_"+label) &&
                    info.Par.count("Lambda3_"+label) &&
                    info.Par.count("Lambda4_"+label)
                    )) {
                std::cout << "[ERROR] Initial parameters where not found for " << tag << " Fourth Order Exponential Chebychev in " << col << std::endl; return false;
              }
	      // create the variables for this model
	      ws.factory( info.Par.at("N_"+label).c_str() );
	      RooArgList pdfConstrains;
	      StringVector_t varNames = {"Lambda1", "Lambda2", "Lambda3", "Lambda4"};
	      for (const auto& v : varNames) {
		if (info.Par.count(v+"_"+label)) { ws.factory( info.Par.at(v+"_"+label).c_str() ); }
		// create the Gaussian PDFs for Constrain fits
		if (info.Par.count("val"+v+"_"+label) && info.Par.count("sig"+v+"_"+label)) {
		  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", (v+"_"+label).c_str(), (v+"_"+label).c_str(), info.Par.at("val"+v+"_"+label).c_str(), info.Par.at("sig"+v+"_"+label).c_str()));
		  pdfConstrains.add( *ws.pdf(("Constr"+v+"_"+label).c_str()) );
		}
	      }
	      if (pdfConstrains.getSize()>0) { ws.import(pdfConstrains, Form("pdfConstr%s", mainLabel.c_str())); }
	      //
              // create the PDF
              ws.factory(Form("RooFormulaVar::%s('( 2.0*@0 - @2 - @1 )/( @2 - @1 )', {%s, vMin[%.6f], vMax[%.6f]})",
                              varNormName.c_str(), varName.c_str(),
                              info.Var.at(varName).at("Min"),
                              info.Var.at(varName).at("Max")
                              ));
              ws.factory(Form("RooFormulaVar::%s('@4*(8.0*@0*@0*@0*@0 - 8.0*@0*@0 + 1.0) + @3*(4.0*@0*@0*@0 - 3.0*@0) + @2*(2.0*@0*@0 - 1.0) + @1*(@0) + 1.0', {%s, %s, %s, %s, %s})", 
                              pdfPolName.c_str(), varNormName.c_str(),
                              Form("Lambda1_%s", label.c_str()),
                              Form("Lambda2_%s", label.c_str()),
                              Form("Lambda3_%s", label.c_str()),
                              Form("Lambda4_%s", label.c_str())
                              ));
              ws.factory(Form("Exponential::%s(%s, One[1.0])", pdfName.c_str(),
                              pdfPolName.c_str()
                              ));
              ws.factory(Form("RooExtendPdf::%s(%s,%s)", pdfTotName.c_str(),
                              pdfName.c_str(),
                              Form("N_%s", label.c_str())
                              ));
	      ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
	      pdfList.add( *ws.pdf(pdfTotName.c_str()) );
              std::cout << "[INFO] " << tag << " 4th Order Exponential Chebychev " << varName << " PDF in " << col << " added!" << std::endl; break;
            }
          case (int(Model::ExpChebychev5)):
            {
              // check that all input parameters are defined 
              if (!(
		    info.Par.count("N_"+label) &&
                    info.Par.count("Lambda1_"+label) &&
                    info.Par.count("Lambda2_"+label) &&
                    info.Par.count("Lambda3_"+label) &&
                    info.Par.count("Lambda4_"+label) &&
                    info.Par.count("Lambda5_"+label)
                    )) {
                std::cout << "[ERROR] Initial parameters where not found for " << tag << " Fifth Order Exponential Chebychev in " << col << std::endl; return false;
              }
	      // create the variables for this model
	      ws.factory( info.Par.at("N_"+label).c_str() );
	      RooArgList pdfConstrains;
	      StringVector_t varNames = {"Lambda1", "Lambda2", "Lambda3", "Lambda4", "Lambda5"};
	      for (const auto& v : varNames) {
		if (info.Par.count(v+"_"+label)) { ws.factory( info.Par.at(v+"_"+label).c_str() ); }
		// create the Gaussian PDFs for Constrain fits
		if (info.Par.count("val"+v+"_"+label) && info.Par.count("sig"+v+"_"+label)) {
		  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", (v+"_"+label).c_str(), (v+"_"+label).c_str(), info.Par.at("val"+v+"_"+label).c_str(), info.Par.at("sig"+v+"_"+label).c_str()));
		  pdfConstrains.add( *ws.pdf(("Constr"+v+"_"+label).c_str()) );
		}
	      }
	      if (pdfConstrains.getSize()>0) { ws.import(pdfConstrains, Form("pdfConstr%s", mainLabel.c_str())); }
	      //
              // create the PDF
              ws.factory(Form("RooFormulaVar::%s('( 2.0*@0 - @2 - @1 )/( @2 - @1 )', {%s, vMin[%.6f], vMax[%.6f]})",
                              varNormName.c_str(), varName.c_str(),
                              info.Var.at(varName).at("Min"),
                              info.Var.at(varName).at("Max")
                              ));
              ws.factory(Form("RooFormulaVar::%s('@5*(16.0*@0*@0*@0*@0*@0 - 20.0*@0*@0*@0 + 5.0*@0) + @4*(8.0*@0*@0*@0*@0 - 8.0*@0*@0 + 1.0) + @3*(4.0*@0*@0*@0 - 3.0*@0) + @2*(2.0*@0*@0 - 1.0) + @1*(@0) + 1.0', {%s, %s, %s, %s, %s, %s})", 
                              pdfPolName.c_str(), varNormName.c_str(), 
                              Form("Lambda1_%s", label.c_str()),
                              Form("Lambda2_%s", label.c_str()),
                              Form("Lambda3_%s", label.c_str()),
                              Form("Lambda4_%s", label.c_str()),
                              Form("Lambda5_%s", label.c_str())
                              ));
              ws.factory(Form("Exponential::%s(%s, One[1.0])", pdfName.c_str(),
                              pdfPolName.c_str()
                              ));
              ws.factory(Form("RooExtendPdf::%s(%s,%s)", pdfTotName.c_str(),
                              pdfName.c_str(),
                              Form("N_%s", label.c_str())
                              ));
	      ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
	      pdfList.add( *ws.pdf(pdfTotName.c_str()) );
              std::cout << "[INFO] " << tag << " 5th Order Exponential Chebychev " << varName << " PDF in " << col << " added!" << std::endl; break;
            }
          case (int(Model::ExpChebychev6)):
            {
              // check that all input parameters are defined 
              if (!(
		    info.Par.count("N_"+label) &&
                    info.Par.count("Lambda1_"+label) &&
                    info.Par.count("Lambda2_"+label) &&
                    info.Par.count("Lambda3_"+label) &&
                    info.Par.count("Lambda4_"+label) &&
                    info.Par.count("Lambda5_"+label) &&
                    info.Par.count("Lambda6_"+label)
                    )) {
                std::cout << "[ERROR] Initial parameters where not found for " << tag << " Sixth Order Exponential Chebychev in " << col << std::endl; return false;
              }
	      // create the variables for this model
	      ws.factory( info.Par.at("N_"+label).c_str() );
	      RooArgList pdfConstrains;
	      StringVector_t varNames = {"Lambda1", "Lambda2", "Lambda3", "Lambda4", "Lambda5", "Lambda6"};
	      for (const auto& v : varNames) {
		if (info.Par.count(v+"_"+label)) { ws.factory( info.Par.at(v+"_"+label).c_str() ); }
		// create the Gaussian PDFs for Constrain fits
		if (info.Par.count("val"+v+"_"+label) && info.Par.count("sig"+v+"_"+label)) {
		  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", (v+"_"+label).c_str(), (v+"_"+label).c_str(), info.Par.at("val"+v+"_"+label).c_str(), info.Par.at("sig"+v+"_"+label).c_str()));
		  pdfConstrains.add( *ws.pdf(("Constr"+v+"_"+label).c_str()) );
		}
	      }
	      if (pdfConstrains.getSize()>0) { ws.import(pdfConstrains, Form("pdfConstr%s", mainLabel.c_str())); }
	      //
              // create the PDF
              ws.factory(Form("RooFormulaVar::%s('( 2.0*@0 - @2 - @1 )/( @2 - @1 )', {%s, vMin[%.6f], vMax[%.6f]})",
                              varNormName.c_str(), varName.c_str(),
                              info.Var.at(varName).at("Min"),
                              info.Var.at(varName).at("Max")
                              ));
              ws.factory(Form("RooFormulaVar::%s('@6*(32.0*@0*@0*@0*@0*@0*@0 - 48.0*@0*@0*@0*@0 + 18.0*@0*@0 - 1.0) + @5*(16.0*@0*@0*@0*@0*@0 - 20.0*@0*@0*@0 + 5.0*@0) + @4*(8.0*@0*@0*@0*@0 - 8.0*@0*@0 + 1.0) + @3*(4.0*@0*@0*@0 - 3.0*@0) + @2*(2.0*@0*@0 - 1.0) + @1*(@0) + 1.0', {%s, %s, %s, %s, %s, %s, %s})", 
                              pdfPolName.c_str(), varNormName.c_str(),
                              Form("Lambda1_%s", label.c_str()),
                              Form("Lambda2_%s", label.c_str()),
                              Form("Lambda3_%s", label.c_str()),
                              Form("Lambda4_%s", label.c_str()),
                              Form("Lambda5_%s", label.c_str()),
                              Form("Lambda6_%s", label.c_str())
                              ));
              ws.factory(Form("Exponential::%s(%s, One[1.0])", pdfName.c_str(),
                              pdfPolName.c_str()
                              ));
              ws.factory(Form("RooExtendPdf::%s(%s,%s)", pdfTotName.c_str(),
                              pdfName.c_str(),
                              Form("N_%s", label.c_str())
                              ));
	      ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
	      pdfList.add( *ws.pdf(pdfTotName.c_str()) );
              std::cout << "[INFO] " << tag << " 6th Order Exponential Chebychev " << varName << " PDF in " << col << " added!" << std::endl; break;
            }
          case (int(Model::Exponential)):
            {
              // check that all input parameters are defined 
              if (!(
		    info.Par.count("N_"+label) &&
                    info.Par.count("Lambda1_"+label)
                    )) { 
                std::cout << "[ERROR] Initial parameters where not found for " << tag << " Exponential in " << col << std::endl; return false;
              }
	      // create the variables for this model
	      ws.factory( info.Par.at("N_"+label).c_str() );
	      RooArgList pdfConstrains;
	      StringVector_t varNames = {"Lambda1"};
	      for (const auto& v : varNames) {
		if (info.Par.count(v+"_"+label)) { ws.factory( info.Par.at(v+"_"+label).c_str() ); }
		// create the Gaussian PDFs for Constrain fits
		if (info.Par.count("val"+v+"_"+label) && info.Par.count("sig"+v+"_"+label)) {
		  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", (v+"_"+label).c_str(), (v+"_"+label).c_str(), info.Par.at("val"+v+"_"+label).c_str(), info.Par.at("sig"+v+"_"+label).c_str()));
		  pdfConstrains.add( *ws.pdf(("Constr"+v+"_"+label).c_str()) );
		}
	      }
	      if (pdfConstrains.getSize()>0) { ws.import(pdfConstrains, Form("pdfConstr%s", mainLabel.c_str())); }
              // create the PDF
              ws.factory(Form("Exponential::%s(%s, %s)", pdfName.c_str(), varName.c_str(), 
                              Form("Lambda1_%s", label.c_str())
                              ));
              ws.factory(Form("RooExtendPdf::%s(%s,%s)", pdfTotName.c_str(),
                              pdfName.c_str(),
                              Form("N_%s", label.c_str())
                              ));
	      ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
	      pdfList.add( *ws.pdf(pdfTotName.c_str()) );
              std::cout << "[INFO] " << tag << " Exponential " << varName << " PDF in " << col << " added!" << std::endl; break;
            }
          case (int(Model::ExpError)):
            {
              // check that all input parameters are defined 
              if (!(
		    info.Par.count("N_"+label) &&
                    info.Par.count("Sigma_"+label) &&
                    info.Par.count("xb_"+label) &&
                    info.Par.count("Lambda_"+label)
                    )) {
                std::cout << "[ERROR] Initial parameters where not found for " << tag << " ExpError in " << col << std::endl; return false;
              }
	      // create the variables for this model
	      ws.factory( info.Par.at("N_"+label).c_str() );
	      RooArgList pdfConstrains;
	      StringVector_t varNames = {"Sigma", "xb", "Lambda"};
	      for (const auto& v : varNames) {
		if (info.Par.count(v+"_"+label)) { ws.factory( info.Par.at(v+"_"+label).c_str() ); }
		// create the Gaussian PDFs for Constrain fits
		if (info.Par.count("val"+v+"_"+label) && info.Par.count("sig"+v+"_"+label)) {
		  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", (v+"_"+label).c_str(), (v+"_"+label).c_str(), info.Par.at("val"+v+"_"+label).c_str(), info.Par.at("sig"+v+"_"+label).c_str()));
		  pdfConstrains.add( *ws.pdf(("Constr"+v+"_"+label).c_str()) );
		}
	      }
	      if (pdfConstrains.getSize()>0) { ws.import(pdfConstrains, Form("pdfConstr%s", mainLabel.c_str())); }
              // create the PDF                 
              ws.factory(Form("RooGenericPdf::%s('TMath::Exp(-@0*@1)*(1.0+TMath::Erf((@0-@2)/@3))', {%s, %s, %s, %s})", pdfName.c_str(), varName.c_str(),
                              Form("Lambda_%s", label.c_str()),
                              Form("xb_%s", label.c_str()),
                              Form("Sigma_%s", label.c_str())
                              ));
              ws.factory(Form("RooExtendPdf::%s(%s,%s)", pdfTotName.c_str(),
                              pdfName.c_str(),
                              Form("N_%s", label.c_str())
                              ));
	      ws.pdf(pdfName.c_str())->setNormRange(varWindow.c_str());
	      pdfList.add( *ws.pdf(pdfTotName.c_str()) );
              std::cout << "[INFO] " << tag << " ExpError " << varName << " PDF in " << col << " added!" << std::endl; break;
            }
          default :
            {
	      if (modelN=="") { std::cout << "[ERROR] Candidate Mass Model for " << modelN << " was not defined (is empty)!" << std::endl; return false; }
	      else { std::cout << "[ERROR] Selected Candidate Mass Model: " << modelN << " has not been implemented" << std::endl; return false; }
	    }
          }
      }
      if (pdfList.getSize()>0) {
        const auto& pdfName = ( "pdf" + varType + "_Tot" + mainLabel );
        auto themodel = std::unique_ptr<RooAddPdf>(new RooAddPdf(pdfName.c_str(), pdfName.c_str(), pdfList));
        ws.import(*themodel);
      }
      pdfListTot.add(pdfList);
    }
    if (pdfListTot.getSize()>0) {
      const auto& pdfName = ( "pdf" + varType + "_Tot" + (cha+chg+"_"+col) );
      info.Par["pdfName"+chg] = pdfName;
      auto themodel = std::unique_ptr<RooAddPdf>(new RooAddPdf(pdfName.c_str(), pdfName.c_str(), pdfListTot));
      ws.import(*themodel);
    }
  }
  return true;
};


#endif // #ifndef Candidate_addModel_C
