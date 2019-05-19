#ifndef buildCandidateMassModel_C
#define buildCandidateMassModel_C


#include "TObject.h"

#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooStringVar.h"

#include <iostream>
#include <string>
#include <memory>

#include "../Utilities/initClasses.h"
#include "../Utilities/rooDataUtils.h"
#include "../Utilities/rooModelUtils.h"


void  constrainQuarkoniumMassParameters ( GlobalInfo& info , const std::string&   chg );
void  setCandidateMassModelParameters   ( GlobalInfo& info , const std::string&   chg );
bool  addCandidateModel                 ( RooWorkspace& ws , const StringDiMap_t& models , const GlobalInfo& info , const std::string& chg , const std::string& varName );


bool buildCandidateMassModel(RooWorkspace& ws, const StringDiMap_t& models, GlobalInfo&  info, const std::string& chg)
{
  //
  // Initialize all the Candidate Mass Model parameters needed for fitting
  setCandidateMassModelParameters(info, chg);
  //
  // Constrain Quarkonium excited states to 1S state
  constrainQuarkoniumMassParameters(info, chg);
  //
  // Import the Candidate Mass Models to the local workspace
  if (!addCandidateModel(ws, models, info, chg, "Cand_Mass")) { return false; }
  //
  // Set Fixed parameters to constant (clean up)
  setFixedVarsToContantVars(ws);
  //
  // save the initial values of the model we've just created
  ws.saveSnapshot("initialParameters", ws.allVars(), true);
  //
  return true;
};


bool addCandidateModel(RooWorkspace& ws, const StringDiMap_t& models, const GlobalInfo& info, const std::string& chg, const std::string& varName)
{
  //
  const auto& cha         = info.Par.at("channel");
  const auto& varWindow   = (varName + "Window");
  const auto& varNormName = (varName+"Norm");
  std::string varType = varName; if (varType.find("_")!=std::string::npos) { varType = varType.substr(varType.find("_")+1, varType.size()); }
  //
  for (const auto& col : info.StrV.at("fitSystem")) {
    RooArgList pdfListTot;
    for (const auto& obj : info.StrV.at("fitObject")) {;
      const auto& mainTag = obj + cha + chg;
      const auto& mainLabel = mainTag + "_" + col;
      std::cout << "[INFO] Implementing " << mainTag << " " << varName << " Model for " << col << std::endl;
      RooArgList pdfList;
      for (const auto& model : models.at("Model_"+mainLabel)) {
	const auto& tag = model.first + cha + chg;
	const auto& label = tag + "_" + col;
	RooStringVar tmp; tmp.setVal(model.second.c_str()); tmp.SetTitle(("Model_"+label).c_str());
        ws.import(*((TObject*)&tmp), tmp.GetTitle()); // Save the model name for bookkeeping
        //
        const std::string& pdfName    = Form("pdf%s_%s",    varType.c_str(), label.c_str());
        const std::string& pdf1Name   = Form("pdf%s1_%s",   varType.c_str(), label.c_str());
        const std::string& pdf2Name   = Form("pdf%s2_%s",   varType.c_str(), label.c_str());
        const std::string& pdfPolName = Form("pdf%sPol_%s", varType.c_str(), label.c_str());
        const std::string& pdfTotName = Form("pdf%sTot_%s", varType.c_str(), label.c_str());
        //
	// Create Models
	switch(ModelDictionary.at(model.second))
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
	      RooStringVar tmp; tmp.setVal(cut.c_str()); tmp.SetTitle(("CutAndCount_"+label).c_str()); ws.import(*((TObject*)&tmp), tmp.GetTitle());
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
	      std::string dsName = ( "d" + chg + "_MC_" + model.first + "_" + info.Par.at("channelDS") + "_" + col );
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
                ws.factory( info.Par.at(v+"_"+label).c_str() );
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
                ws.factory( info.Par.at(v+"_"+label).c_str() );
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
		ws.factory( info.Par.at(v+"_"+label).c_str() );
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
		ws.factory( info.Par.at(v+"_"+label).c_str() );
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
		ws.factory( info.Par.at(v+"_"+label).c_str() );
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
		ws.factory( info.Par.at(v+"_"+label).c_str() );
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
		ws.factory( info.Par.at(v+"_"+label).c_str() );
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
		ws.factory( info.Par.at(v+"_"+label).c_str() );
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
		ws.factory( info.Par.at(v+"_"+label).c_str() );
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
		ws.factory( info.Par.at(v+"_"+label).c_str() );
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
		ws.factory( info.Par.at(v+"_"+label).c_str() );
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
		ws.factory( info.Par.at(v+"_"+label).c_str() );
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
		ws.factory( info.Par.at(v+"_"+label).c_str() );
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
		ws.factory( info.Par.at(v+"_"+label).c_str() );
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
		ws.factory( info.Par.at(v+"_"+label).c_str() );
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
          default :
            {
	      if (model.second=="") { std::cout << "[ERROR] Candidate Mass Model for " << model.second << " was not defined (is empty)!" << std::endl; return false; }
	      else { std::cout << "[ERROR] Selected Candidate Mass Model: " << model.second << " has not been implemented" << std::endl; return false; }
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
      auto themodel = std::unique_ptr<RooAddPdf>(new RooAddPdf(pdfName.c_str(), pdfName.c_str(), pdfListTot));
      ws.import(*themodel);
    }
  }
  return true;
};


void setCandidateMassModelParameters(GlobalInfo& info, const std::string& chg)
{
  std::cout << "[INFO] Initializing Candidate Mass Model parameters based on user input!" << std::endl;
  std::string cha = info.Par.at("channel");
  for (const auto& col : info.StrV.at("fitSystem")) {
    for (const auto& mainObj : info.StrV.at("fitObject")) {
      for (const auto& obj : info.StrV.at("addObjectModel_"+(mainObj+cha+chg)+"_"+col)) {
	const auto& mainCha = cha , mainChg = chg , mainCol = col;
	std::string foundCha = cha , foundChg = chg , foundCol = col;
	const StringVector_t tryChannel = { cha , "" };
	StringVector_t trySystem  = { col , "" };
	const StringVector_t tryCharge  = { chg , "" };
	for (const auto& tryCha : tryChannel) {
	  bool trySuccess = false;
	  for (const auto& tryCol : trySystem) {
	    for (const auto& tryChg : tryCharge) {
	      foundCha = tryCha; foundChg = tryChg; foundCol = tryCol;
	      const auto& tryLabel = obj + foundCha + foundChg + ( (tryCol!="") ? ("_" + foundCol) : "" );
	      if (info.Par.count(("Model_" + tryLabel).c_str())>0) { trySuccess = true; break; }
	    }
	    if (trySuccess) break;
	  }
	  if (trySuccess) break;
	}
	//
	const auto& objLabel      = obj + mainCha  + mainChg  + ( (mainCol !="") ? ("_" + mainCol ) : "" );
	const auto& objFoundLabel = obj + foundCha + foundChg + ( (foundCol!="") ? ("_" + foundCol) : "" );
	//
	// NUMBER OF EVENTS
	if (info.Par.count("N_"+objLabel)==0 || info.Par.at("N_"+objLabel)=="") {
	  if (info.Par.count("N_"+objFoundLabel)==0 || info.Par.at("N_"+objFoundLabel)=="") {
            const auto& numEntries = info.Var.at("numEntries").at(chg);
            info.Par["N_"+objLabel] = Form("%s[%.10f,%.10f,%.10f]", ("N_"+objLabel).c_str(), numEntries, -2.0*numEntries, 2.0*numEntries);
          }
	  else {
	    std::string content = info.Par.at("N_"+objFoundLabel); content = content.substr( content.find("[") );
	    info.Par["N_"+objLabel] = Form("%s%s", ("N_"+objLabel).c_str(), content.c_str());
	  }
	}
        //
        // Check if it is for template, in which case continue
        if (info.Flag.count("incMCTemp_"+mainObj+"_"+obj)>0 && info.Flag.at("incMCTemp_"+mainObj+"_"+obj)) continue;
        //
        // CUTS FOR CUT AND COUNT ALGO
        if (info.Par.count("Cut_"+objLabel)==0 || info.Par.at("Cut_"+objLabel)=="") {
          if (info.Par.count("Cut_"+objFoundLabel)==0 || info.Par.at("Cut_"+objFoundLabel)=="") {
            info.Par["Cut_"+objLabel] = "";
          }
          else {
            info.Par["Cut_"+objLabel] = info.Par.at("Cut_"+objFoundLabel);
          }
        }
        // MASS MODEL PARAMETERS
        const StringVector_t varNames = {"m", "Sigma1", "rSigma21", "Sigma2", "Alpha", "Alpha2", "n", "n2", "f", "Lambda1", "Lambda2", "Lambda3", "Lambda4", "Lambda5", "Lambda6"};
        for (const auto v : varNames) {
          if (obj!="Bkg" && v.rfind("Lambda",0)==0) continue;
          if (info.Par.count(v+"_"+objLabel)==0 || info.Par.at(v+"_"+objLabel)=="") {
            if (info.Par.count(v+"_"+objFoundLabel)==0 || info.Par.at(v+"_"+objFoundLabel)=="") {
              if (v=="Sigma1"  ) { info.Par[v+"_"+objLabel] = Form("%s[%.4f,%.4f,%.4f]", (v+"_"+objLabel).c_str(), 0.120, 0.005,  0.240); }
              if (v=="rSigma21") { info.Par[v+"_"+objLabel] = Form("%s[%.4f,%.4f,%.4f]", (v+"_"+objLabel).c_str(), 2.000, 1.000,  4.000); }
              if (v=="Alpha"   ) { info.Par[v+"_"+objLabel] = Form("%s[%.4f,%.4f,%.4f]", (v+"_"+objLabel).c_str(), 2.000, 0.500, 30.000); }
              if (v=="n"       ) { info.Par[v+"_"+objLabel] = Form("%s[%.4f,%.4f,%.4f]", (v+"_"+objLabel).c_str(), 1.800, 0.500, 10.000); }
              if (v=="f"       ) { info.Par[v+"_"+objLabel] = Form("%s[%.4f,%.4f,%.4f]", (v+"_"+objLabel).c_str(), 0.500, 0.000,  1.000); }
              if (v=="m") {
                if (MASS.count(obj)>0) { info.Par[v+"_"+objLabel] = Form("%s[%.10f,%.10f,%.10f]", (v+"_"+objLabel).c_str(),
                                                                         MASS.at(obj).at("Val"),
                                                                         (MASS.at(obj).at("Val") - 2.0*MASS.at(obj).at("Width")),
                                                                         (MASS.at(obj).at("Val") + 2.0*MASS.at(obj).at("Width"))
                                                                         ); }
                else if (obj!="Bkg" || obj.find("Swap")!=std::string::npos) {
                  info.Par[v+"_"+objLabel] = (v+"_"+objLabel+"[10.0,0.0,1000.0]"); std::cout << "[WARNING] Initial value for " << (v+"_"+objLabel) << " was not found!" << std::endl;
                }
              }
              if (info.Par.count(v+"_"+objLabel)>0 || info.Par.count(v+"_"+objFoundLabel)>0) {
                if (v=="Sigma2") { info.Par[v+"_"+objLabel] = Form("%s[%.4f,%.4f,%.4f]", (v+"_"+objLabel).c_str(), 0.040, 0.010,  0.100); }
                if (v=="Alpha2") { info.Par[v+"_"+objLabel] = Form("%s[%.4f,%.4f,%.4f]", (v+"_"+objLabel).c_str(), 2.000, 0.500, 30.000); }
                if (v=="n2"    ) { info.Par[v+"_"+objLabel] = Form("%s[%.4f,%.4f,%.4f]", (v+"_"+objLabel).c_str(), 1.800, 0.500, 10.000); }
              }
              else {
                if (v=="Sigma2") { info.Par[v+"_"+objLabel] = Form("RooFormulaVar::%s('@0*@1',{%s,%s})", ("Sigma2_"+objLabel).c_str(), ("rSigma21_"+objLabel).c_str(), ("Sigma1_"+objLabel).c_str()); }
                if (v=="Alpha2") { info.Par[v+"_"+objLabel] = Form("RooFormulaVar::%s('@0',{%s})", (v+"_"+objLabel).c_str(), ("Alpha_"+objLabel).c_str()); }
                if (v=="n2"    ) { info.Par[v+"_"+objLabel] = Form("RooFormulaVar::%s('@0',{%s})", (v+"_"+objLabel).c_str(), ("n_"+objLabel).c_str());     }
              }
              if (v.rfind("Lambda",0)==0) { info.Par[v+"_"+objLabel] = Form("%s[%.4f,%.4f,%.4f]", (v+"_"+objLabel).c_str(), 0.0, -100.0, 100.0); }
            }
            else {
              std::string content = info.Par.at(v+"_"+objFoundLabel); content = content.substr( content.find("[") );
              info.Par[v+"_"+objLabel] = Form("%s%s", (v+"_"+objLabel).c_str(), content.c_str());
            }
          }
          // Check Parameters for Constrain Fits
          StringVector_t constrainLabel = { "val" , "sig" };
          for (const auto& con : constrainLabel) {
            const std::string& name = Form("%s%s_%s", con.c_str(), v.c_str(), objLabel.c_str());
            if (info.Par.count(name) && info.Par.at(name)!="") {
              std::string content = info.Par.at(name); content = content.substr( content.find("[") );
              info.Par[Form("%s%s_%s", con.c_str(), v.c_str(), objLabel.c_str())] = Form("%s%s_%s%s", con.c_str(), v.c_str(), objLabel.c_str(), content.c_str());
            }
            else break;
          }
        }
      }
    }
  }
};


void constrainQuarkoniumMassParameters(GlobalInfo& info, const std::string& chg)
{
  std::cout << "[INFO] Constraining excited state mass parameters to reference state using PDF Mass Ratio!" << std::endl;
  std::string cha = info.Par.at("channel");
  for (const auto& col : info.StrV.at("fitSystem")) {
    const StringVector_t refStates = { "JPsi" , "Ups1S" };
    for (const auto& refState : refStates) {
      const auto& refLabel = refState + cha + chg  + "_" + col;
      if (info.Par.count(("m_"+refLabel).c_str())>0) {
        StringVector_t excStates;
        if (refState=="JPsi" ) { excStates.push_back("Psi2S"); }
        if (refState=="Ups1S") { excStates.push_back("Ups2S"); excStates.push_back("Ups3S"); }
        for (const auto& excState : excStates) {
          const auto& excLabel = excState + cha + chg  + "_" + col;
          if (info.Par.count(("m_"+excLabel).c_str())>0) {
            const StringVector_t varList = { "m", "Sigma1", "Sigma2", "Alpha", "Alpha2", "n", "n2", "f"};
            for (const auto& v : varList) {
              if (v=="m" || v=="Sigma1" || v=="Sigma2") {
                const double& massRatioValue = (MASS.at(excState).at("Val")/MASS.at(refState).at("Val"));
                const auto& massRatioLabel = ("MassRatio_" + excState + "Over" + refState);
                info.Par[massRatioLabel] = Form("%s[%.6f]", massRatioLabel.c_str(), massRatioValue);
                info.Par.at(v+"_"+excLabel) = Form("RooFormulaVar::%s('@0*@1',{%s,%s})", (v+"_"+excLabel).c_str(), info.Par.at(massRatioLabel).c_str(), (v+"_"+refLabel).c_str());
              }
              else { info.Par.at(v+"_"+excLabel) = Form("RooFormulaVar::%s('@0',{%s})", (v+"_"+excLabel).c_str(), (v+"_"+refLabel).c_str()); }
            }
          }
        }
      }
    }
  }
};


#endif // #ifndef buildElectroWeakMETModel_C
