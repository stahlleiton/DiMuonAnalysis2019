#ifndef Candidate_buildCandidateMassModel_C
#define Candidate_buildCandidateMassModel_C


#include "RooWorkspace.h"
#include "RooArgSet.h"

#include <iostream>
#include <string>

#include "addModel.C"
#include "../Utilities/initClasses.h"
#include "../Utilities/rooDataUtils.h"


void  constrainQuarkoniumMassParameters ( GlobalInfo& info , const std::string&   chg );
void  setCandidateMassModelParameters   ( GlobalInfo& info , const std::string&   chg , const std::string& var="Cand_Mass" );


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
  if (!addModel(ws, info, models, chg, "Cand_Mass")) { return false; }
  //
  // Set Fixed parameters to constant (clean up)
  setFixedVarsToContantVars(ws);
  //
  // save the initial values of the model we've just created
  RooArgSet set = ws.allVars(); set.add(ws.allCats()); set.remove(*ws.data(info.Par.at("dsName"+chg).c_str())->get(), false, true);
  ws.saveSnapshot("initialParameters", set, true);
  //
  return true;
};


void setCandidateMassModelParameters(GlobalInfo& info, const std::string& chg, const std::string& var)
{
  std::cout << "[INFO] Initializing " << var << " Model parameters based on user input!" << std::endl;
  std::string cha = info.Par.at("channel");
  for (const auto& col : info.StrS.at("fitSystem")) {
    for (const auto& mainObj : info.StrS.at("fitObject")) {
      for (const auto& obj : info.StrS.at("addObjectModel_"+(mainObj+cha+chg)+"_"+col)) {
	const auto& objLabel      = obj + cha  + chg  + (col!="" ? "_"+col : "");
	const auto& objFoundLabel = findLabel("Model", obj, chg, col, cha, info);
	const bool& isSwap = (obj.find("Swap")!=std::string::npos);
	//
	// NUMBER OF EVENTS
	if (info.Par.count("N_"+objLabel)==0 || info.Par.at("N_"+objLabel)=="") {
	  if (info.Par.count("N_"+objFoundLabel)==0 || info.Par.at("N_"+objFoundLabel)=="") {
            const auto& numEntries = info.Var.at("numEntries").at(chg);
            info.Par["N_"+objLabel] = Form("%s[%.10f,%.10f,%.10f]", ("N_"+objLabel).c_str(), numEntries, -0.001*numEntries, 2.0*numEntries);
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
        const StringVector_t varNames = {"m", "Sigma1", "rSigma21", "Sigma2", "Alpha", "Alpha2", "n", "n2", "f", "Lambda1", "Lambda2", "Lambda3", "Lambda4", "Lambda5", "Lambda6", "Lambda", "Sigma", "xb"};
        for (const auto& v : varNames) {
          if (obj!="Bkg" && v.rfind("Lambda",0)==0) continue;
          if (info.Par.count(v+"_"+objLabel)==0 || info.Par.at(v+"_"+objLabel)=="") {
            if (info.Par.count(v+"_"+objFoundLabel)==0 || info.Par.at(v+"_"+objFoundLabel)=="") {
	      const auto& varV = ((var=="Cand_Mass" && MASS.count(obj)>0) ? MASS.at(obj).at("Val")   : ((info.Var.at(var).at("Max")+info.Var.at(var).at("Min"))/2.0));
	      const auto& varR = ((var=="Cand_Mass" && MASS.count(obj)>0) ? MASS.at(obj).at("Width") : ((info.Var.at(var).at("Max")-info.Var.at(var).at("Min"))/6.0));
	      if      (v=="Sigma1"  ) { info.Par[v+"_"+objLabel] = Form("%s[%.4f,%.4f,%.4f]", (v+"_"+objLabel).c_str(), varR, 0.05*varR, 3.0*varR); }
              else if (v=="rSigma21") { info.Par[v+"_"+objLabel] = Form("%s[%.4f,%.4f,%.4f]", (v+"_"+objLabel).c_str(), 2.000, 1.000,  4.000); }
              else if (v=="Alpha"   ) { info.Par[v+"_"+objLabel] = Form("%s[%.4f,%.4f,%.4f]", (v+"_"+objLabel).c_str(), 2.000, 0.500, 30.000); }
              else if (v=="n"       ) { info.Par[v+"_"+objLabel] = Form("%s[%.4f,%.4f,%.4f]", (v+"_"+objLabel).c_str(), 1.800, 0.500, 10.000); }
              else if (v=="f"       ) { info.Par[v+"_"+objLabel] = Form("%s[%.4f,%.4f,%.4f]", (v+"_"+objLabel).c_str(), 0.500, 0.000,  1.000); }
              else if (v=="m"       ) { info.Par[v+"_"+objLabel] = Form("%s[%.9f,%.9f,%.9f]", (v+"_"+objLabel).c_str(), varV, (varV - 3.0*varR), (varV + 3.0*varR)); }
	      if (v=="m" && var=="Cand_Mass" && MASS.count(obj)==0 && (obj!="Bkg" || isSwap)) {
		std::cout << "[WARNING] Initial value for " << (v+"_"+objLabel) << " was not found!" << std::endl;
	      }
              if (info.Par.count(v+"_"+objLabel)>0 || info.Par.count(v+"_"+objFoundLabel)>0) {
                if (v=="Sigma2") { info.Par[v+"_"+objLabel] = Form("%s[%.4f,%.4f,%.4f]", (v+"_"+objLabel).c_str(), varV, (varV - 3.0*varR), (varV + 3.0*varR)); }
                if (v=="Alpha2") { info.Par[v+"_"+objLabel] = Form("%s[%.4f,%.4f,%.4f]", (v+"_"+objLabel).c_str(), 2.000, 0.500, 30.000); }
                if (v=="n2"    ) { info.Par[v+"_"+objLabel] = Form("%s[%.4f,%.4f,%.4f]", (v+"_"+objLabel).c_str(), 1.800, 0.500, 10.000); }
              }
              else {
                if (v=="Sigma2") { info.Par[v+"_"+objLabel] = Form("RooFormulaVar::%s('@0*@1',{%s,%s})", ("Sigma2_"+objLabel).c_str(), ("rSigma21_"+objLabel).c_str(), ("Sigma1_"+objLabel).c_str()); }
                if (v=="Alpha2") { info.Par[v+"_"+objLabel] = Form("RooFormulaVar::%s('@0',{%s})", (v+"_"+objLabel).c_str(), ("Alpha_"+objLabel).c_str()); }
                if (v=="n2"    ) { info.Par[v+"_"+objLabel] = Form("RooFormulaVar::%s('@0',{%s})", (v+"_"+objLabel).c_str(), ("n_"+objLabel).c_str());     }
              }
              if (v=="Lambda") { info.Par[v+"_"+objLabel] = Form("%s[%.4f,%.4f,%.4f]", (v+"_"+objLabel).c_str(),  0.2, -2.00,  2.00); }
              else if (v.rfind("Lambda",0)==0) { info.Par[v+"_"+objLabel] = Form("%s[%.4f,%.4f,%.4f]", (v+"_"+objLabel).c_str(), 0.0, -100.0, 100.0); }
              else if (v=="Sigma") { info.Par[v+"_"+objLabel] = Form("%s[%.4f,%.4f,%.4f]", (v+"_"+objLabel).c_str(),  0.4,  0.01, 40.00); }
              else if (v=="xb"   ) { info.Par[v+"_"+objLabel] = Form("%s[%.4f,%.4f,%.4f]", (v+"_"+objLabel).c_str(), 10.0,  0.00, 40.00); }
            }
            else {
              std::string content = info.Par.at(v+"_"+objFoundLabel); content = content.substr( content.find("[") );
              info.Par[v+"_"+objLabel] = Form("%s%s", (v+"_"+objLabel).c_str(), content.c_str());
            }
          }
          // Check parameters for constrained fits
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
  const StringVector_t varList = { "m", "Sigma1", "Sigma2", "Alpha", "Alpha2", "n", "n2", "f" };
  const StringVectorMap_t excStates = { {"JPsi", {"Psi2S"}} , {"Ups1S", {"Ups2S", "Ups3S"}} };
  std::string cha = info.Par.at("channel");
  for (const auto& col : info.StrS.at("fitSystem")) {
    for (const auto& refSt : excStates) {
      const auto& refLabel = refSt.first + cha + chg  + "_" + col;
      if (info.Par.count(("m_"+refLabel).c_str())>0) {
        for (const auto& excState : refSt.second) {
          const auto& excLabel = excState + cha + chg  + "_" + col;
          if (info.Par.count(("m_"+excLabel).c_str())>0) {
            for (const auto& v : varList) {
              if (v=="m" || v=="Sigma1" || v=="Sigma2") {
                const double& massRatioValue = (MASS.at(excState).at("Val")/MASS.at(refSt.first).at("Val"));
                const auto& massRatioLabel = ("MassRatio_" + excState + "Over" + refSt.first);
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


#endif // #ifndef Candidate_buildCandidateMassModel_C
