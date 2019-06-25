#ifndef processResults_C
#define processResults_C

// Auxiliary Files
#include "extractResultsTree.C"
// ROOT headers
// c++ headers
#include <iostream>
#include <string>


void iniAcceptanceAndEfficiency ( BinSextaMap_t& eff , const VarBinTriMap_t& inputVar );
bool correctRawYields           ( VarBinTriMap_t& inputVar , const BinSextaMap_t& effMap , const std::string& accType = "" , const std::string& effType = "" );
bool computeRatioTo1S           ( BinSextaMap_t& var , const VarBinTriMap_t& inputVar );
bool computeCrossSection        ( BinSextaMap_t& var , const VarBinTriMap_t& inputVar );



bool processResults(
		    BinSextaMap_t& var,
		    VarBinTriMap_t& inputVar,
		    const std::string& effType = "",
		    const std::string& accType = ""
		    )
{
  //
  // --------------------------------------------------------------------------------- //
  //
  if (effType!="" || accType!="") {
    //
    // Initialize the Acceptance and Efficiency container
    BinSextaMap_t eff;
    iniAcceptanceAndEfficiency(eff, inputVar);
    //
    // Proceed to correct the Raw Yields
    //
    if (!correctRawYields(inputVar, eff, accType, effType)) { return false; }
    //
  }
  //
  // Compute the Ratio
  std::cout << "[INFO] Computing ratios to 1S" << std::endl;
  if (!computeRatioTo1S(var, inputVar)) { return false; }
  //
  // Compute the Cross-Section
  std::cout << "[INFO] Computing cross sections" << std::endl;
  if (!computeCrossSection(var, inputVar)) { return false; }
  //
  // return
  return true;
};


void addParameters(
		   BinSextaMap_t& var,
		   VarBinTriMap_t& inputVar
		   )
{
  //
  for (const auto& o : inputVar) {
    for (const auto& c : o.second) {
      for (const auto& pd : c.second) {
	for (const auto& b : pd.second) {
	  for (const auto& v : b.second) {
	    if (v.first.find("To")==std::string::npos || v.first.find("SS")!=std::string::npos) continue;
	    const auto& vN = v.first.substr(0, v.first.find("To"));
	    const auto& oN = (vN.find("_")!=std::string::npos ? vN.substr(vN.find("_")+1) : vN);
	    const auto& pN = vN.substr(0, vN.find("_"));
	    // Initialize the paramater
	    auto& oVar = var[oN][c.first][pd.first][pN];
	    // Store the paramater result
	    for (const auto& l : v.second) { oVar[l.first][b.first] = l.second; }
	  }
	}
      }
    }
  }
};


void iniAcceptanceAndEfficiency(BinSextaMap_t& eff, const VarBinTriMap_t& inputVar)
{
  //
  const StringVector_t lbl = {"Val", "Err_Stat_High", "Err_Stat_Low", "Err_Syst_High", "Err_Syst_Low", "Err_Tot_High", "Err_Tot_Low"};
  const StringVector_t obj = {"JPsi", "Psi2S", "Ups1S", "Ups2S", "Ups3S", "Z", "D0"};
  //
  for (const auto& o : obj) {
    for (const auto& c : (contain(inputVar, o) ? inputVar.at(o) : inputVar.begin()->second)) {
      for (const auto& pd : c.second) {
	for (const auto& b : pd.second) {
	  for (const auto& l : lbl) {
	    // For MC Acceptance
	    eff[o][c.first][pd.first]["Acceptance_MC"][l][b.first] = ((l=="Val") ? 1.0 : 0.0);
	    // For MC Efficiency
	    eff[o][c.first][pd.first]["Efficiency_MC"][l][b.first] = ((l=="Val") ? 1.0 : 0.0);
	  }
	}
      }
    }
  }
};


double getCorrectedYieldValue(const double& rawN, const double& acc, const double& eff)
{
  return (rawN / (acc * eff));
};

  
double getCorrectedYieldError(const double& rawN, const double& acc, const double& eff, const double& errRawN, const double& errAcc, const double& errEff)
{
  const auto& common = getCorrectedYieldValue(rawN, acc, eff);
  const auto& relErrRawN = (errRawN / rawN);
  const auto& relErrAcc  = (errAcc / acc);
  const auto& relErrEff  = (errEff / eff);
  return (std::abs(common) * sumErrors({relErrRawN, relErrAcc, relErrEff}));
};


bool correctRawYields(VarBinTriMap_t& inputVar, const BinSextaMap_t& effMap, const std::string& accType, const std::string& effType)
{
  //
  const std::string accName = ((accType!="") ? Form("Acceptance_%s", accType.c_str()) : "");
  const std::string effName = ((effType!="") ? Form("Efficiency_%s", effType.c_str()) : "");
  //
  for (const auto& o : inputVar) {
    for (const auto& c : o.second) {
      for (const auto& pd : c.second) {
	for (const auto& b : pd.second) {
	  auto& var = inputVar.at(o.first).at(c.first).at(pd.first).at(b.first);
	  // Loop over the variables
	  for (const auto& v : var) {
	    if (v.first.rfind("N_",0)!=0 || v.first.find("To")==std::string::npos) continue;
	    if (v.first.rfind("N_Bkg",0)==0 || v.first.find("SS")!=std::string::npos) continue;
	    const auto& vN = v.first;
	    const auto& obj = vN.substr(0, vN.find("To")).substr(2);
	    // Check the efficiency container
	    if (!contain(effMap, obj)) { std::cout << "[ERROR] Efficiency container does not have the object " << obj << std::endl; return false; }
	    if (!contain(effMap.at(obj), c.first)) { std::cout << "[ERROR] Efficiency container does not have the system " << c.first << std::endl; return false; }
	    if (!contain(effMap.at(obj).at(c.first), pd.first)) { std::cout << "[ERROR] Efficiency container does not have the PD " << pd.first << std::endl; return false; }
	    const auto& eff = effMap.at(obj).at(c.first).at(pd.first);
	    if (accName!="" && !contain(eff, accName)) { std::cout << "[ERROR] Efficiency container does not have the variable " << accName << std::endl; return false; }
	    if (effName!="" && !contain(eff, effName)) { std::cout << "[ERROR] Efficiency container does not have the variable " << effName << std::endl; return false; }
	    if ((accName!="" && !contain(eff.at(accName).at("Val"), b.first)) || (effName!="" && !contain(eff.at(effName).at("Val"), b.first))) {
	      std::cout << "[ERROR] Efficiency container does not have the bin: "; b.first.print(); return false;
	    }
	    // Extract the acceptance and efficiency
	    DoubleMap_t Acceptance , Efficiency;
	    for (const auto& t : eff.at(accName)) { Acceptance[t.first] = ((accName!="") ? t.second.at(b.first) : ((t.first=="Val") ? 1.0 : 0.0)); }
	    for (const auto& t : eff.at(effName)) { Efficiency[t.first] = ((effName!="") ? t.second.at(b.first) : ((t.first=="Val") ? 1.0 : 0.0)); }
	    //
	    // Create a backup
	    for (const auto& t : v.second) { var[vN+"_RAW"][t.first] = t.second; }
	    for (const auto& l : eff) {
	      if (l.first.rfind("_MC_")!=std::string::npos) continue;
	      for (const auto& t : l.second) { var[accName][t.first] = t.second.at(b.first); }
	    }
	    //
	    auto&       N_Corr = var.at(vN);
	    const auto& N_Raw  = var.at(vN+"_RAW");
	    //
	    // Fill with the corrected values;
	    N_Corr.at("Val") = getCorrectedYieldValue(N_Raw.at("Val"), Acceptance.at("Val"), Efficiency.at("Val"));
	    N_Corr.at("Err_Stat_Low" ) = getCorrectedYieldError(N_Raw.at("Val"), Acceptance.at("Val"), Efficiency.at("Val"), N_Raw.at("Err_Stat_Low" ), 0.0 , 0.0);
	    N_Corr.at("Err_Stat_High") = getCorrectedYieldError(N_Raw.at("Val"), Acceptance.at("Val"), Efficiency.at("Val"), N_Raw.at("Err_Stat_High"), 0.0 , 0.0);
	    N_Corr.at("Err_Syst_Low" ) = getCorrectedYieldError(N_Raw.at("Val"), Acceptance.at("Val"), Efficiency.at("Val"), N_Raw.at("Err_Syst_Low" ), Acceptance.at("Err_Tot_Low" ), Efficiency.at("Err_Tot_Low" ));
	    N_Corr.at("Err_Syst_High") = getCorrectedYieldError(N_Raw.at("Val"), Acceptance.at("Val"), Efficiency.at("Val"), N_Raw.at("Err_Syst_High"), Acceptance.at("Err_Tot_High"), Efficiency.at("Err_Tot_High"));
	  }
	}
      }
    }
  }
  return true;
};


bool computeRatioTo1S(BinSextaMap_t& var, const VarBinTriMap_t& inputVar)
{
  //
  for (const auto& o : inputVar) {
    for (const auto& c : o.second) {
      for (const auto& pd : c.second) {
	for (const auto& b : pd.second) {
	  for (const auto& v : b.second) {
	    if (v.first.rfind("R_",0)!=0 || v.first.find("To")==std::string::npos) continue;
	    if (v.first.rfind("N_Bkg",0)==0 || v.first.find("SS")!=std::string::npos) continue;
	    const auto& obj = v.first.substr(0, v.first.find("To")).substr(2);
	    // Initialize ratio to 1S
	    const auto& iVar = v.second;
	    auto& oVar = var[obj][c.first][pd.first]["RatioTo1S"];
	    // Copy value and errors of ratio to 1S
	    for (const auto& p : iVar) { oVar[p.first][b.first] = p.second; }
	  }
	}
      }
    }
  }
  //
  return true;
};


double getCrossSectionValue(const double& N, const double& lumi)
{
  return (N / lumi);
};

  
double getCrossSectionError(const double& N, const double& lumi, const double& errN, const double& errLumi)
{
  const auto& common = getCrossSectionValue(N, lumi);
  const auto& relErrN = (errN / N);
  const auto& relErrL = (errLumi / lumi); 
  return (std::abs(common) * sumErrors({relErrN, relErrL}));
};


bool computeCrossSection(BinSextaMap_t& var, const VarBinTriMap_t& inputVar)
{
  //
  for (const auto& o : inputVar) {
    for (const auto& c : o.second) {
      for (const auto& pd : c.second) {
	for (const auto& b : pd.second) {
	  for (const auto& v : b.second) {
	    if (v.first.rfind("N_",0)!=0 || v.first.find("To")==std::string::npos) continue;
	    if (v.first.rfind("N_Bkg",0)==0 || v.first.find("SS")!=std::string::npos) continue;
	    const auto& obj = v.first.substr(0, v.first.find("To")).substr(2);
	    // Initialize cross section
	    const auto& iVar = v.second;
	    auto& oVar = var[obj][c.first][pd.first]["Cross_Section"];
	    // Extract the luminosity
	    const auto& lumi = b.second.at("Luminosity").at("Val");
	    // Compute the cross-section value and errors
	    for (const auto& t : iVar) {
	      if (t.first=="Val") { oVar[t.first][b.first] = getCrossSectionValue(t.second, lumi); }
	      else if (t.first.rfind("Err_",0)==0) { oVar[t.first][b.first] = getCrossSectionError(iVar.at("Val"), lumi, t.second, 0.0); }
	    }
	  }
	}
      }
    }
  }
  //
  return true;
};


#endif // #ifndef extractResultsTree_C
