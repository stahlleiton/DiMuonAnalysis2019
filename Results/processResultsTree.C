#ifndef processResults_C
#define processResults_C

// Auxiliary Files
#include "extractResultsTree.C"
// ROOT headers
#include "TEfficiency.h"
// c++ headers
#include <iostream>
#include <string>


bool getAcceptanceAndEfficiency ( BinSextaMap_t& eff , const VarBinTriMap_t& inputVar );
bool correctRawYields           ( VarBinTriMap_t& inputVar , const BinSextaMap_t& effMap);
bool computeRatioTo1S           ( BinSextaMap_t& var , const VarBinTriMap_t& inputVar );
bool computeCrossSection        ( BinSextaMap_t& var , const VarBinTriMap_t& inputVar );



bool processResults(
		    BinSextaMap_t& var,
		    VarBinTriMap_t& inputVar
		    )
{
  //
  // Extract the Acceptance and Efficiency container
  BinSextaMap_t eff;
  std::cout << "Extracting the acceptance and efficiency" << std::endl;
  if (!getAcceptanceAndEfficiency(eff, inputVar)) { return false; }
  //
  // Proceed to correct the Raw Yields
  std::cout << "Correcting the raw yields" << std::endl;
  if (!correctRawYields(inputVar, eff)) { return false; }
  //
  // Compute the Ratio
  std::cout << "[INFO] Computing ratios to 1S" << std::endl;
  if (!computeRatioTo1S(var, inputVar)) { return false; }
  //
  // Compute the Cross-Section
  //std::cout << "[INFO] Computing cross sections" << std::endl;
  //if (!computeCrossSection(var, inputVar)) { return false; }
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


bool getAcceptanceAndEfficiency(BinSextaMap_t& eff, const VarBinTriMap_t& inputVar)
{
  //
  const StringVector_t lbl = {"Val", "Err_Stat_High", "Err_Stat_Low", "Err_Syst_High", "Err_Syst_Low", "Err_Tot_High", "Err_Tot_Low"};
  StringVector_t effObj = { "JPsiNoPR", "JPsiPR", "Psi2SNoPR", "Psi2SPR" };
  StringVector_t effTyp = { "Efficiency", "Acceptance" };
  //
  for (const auto& t : effTyp) {
    for (const auto& o : effObj) {
      for (const auto& c : inputVar.begin()->second) {
	for (const auto& pd : c.second) {
	  //
	  const std::string& pdT = (pd.first.rfind("MUON")!=std::string::npos ? "WT" : "NT");
	  std::map< binF , TEfficiency > effM;
	  for (const auto& r : std::vector< binF >({binF("Cand_AbsRap", 0.0, 1.4),binF("Cand_AbsRap", 1.4, 2.4)})) {
	    const std::string& rapS = (r==binF("Cand_AbsRap", 0.0, 1.4) ? "0_1.4" : (r==binF("Cand_AbsRap", 1.4, 2.4) ? "1.4_2.4" : ""));
	    const std::string& objT = std::string(o.rfind("NoPR")!=std::string::npos ? "NonPrompt" : "Prompt") + (pd.first.rfind("MUON")!=std::string::npos ? "WT" : "NT");
	    const std::string& objP = (o.rfind("Psi2S",0)==0 ? "Psi2S" : "JPsi");
	    const std::string& effL = (t=="Acceptance" ? "acc" : "eff_comb");
	    const std::string& dirName = ("Efficiency/MC/" + objT + "/" + objP + "/" + rapS + "/");
	    const auto& fileName = (objT+objP+rapS+"_"+effL+".root");
	    if (!existFile(dirName+fileName)) { std::cout << "[ERROR] File " << (dirName+fileName) << " does not exist!" << std::endl; return false; }
	    TFile file((dirName+fileName).c_str(), "READ");
	    if (!file.IsOpen() || file.IsZombie()) { std::cout << "[ERROR] File " << (dirName+fileName) << " was not found!" << std::endl; return false; }
	    const auto& effP = dynamic_cast<TEfficiency*>(file.Get(effL.c_str()));
	    if (!effP) { std::cout << "[ERROR] Efficiency " << effL << " in " << (dirName+fileName) << " was not found!" << std::endl; file.Close(); return false; }
            effM[r] = *effP;
          }
	  //
	  // Extract efficiency
          for (const auto& b : pd.second) {
	    const auto& rap = b.first.getbin("Cand_Rap"); // x-axis
	    const auto& pt  = b.first.getbin("Cand_Pt"); // y-axis
	    const binF absrap("Cand_AbsRap", std::min(fabs(rap.low()), fabs(rap.high())), std::max(fabs(rap.low()), fabs(rap.high())));
	    //
            const auto& effP = effM.at(absrap);
	    const auto& iBin = effP.FindFixBin(pt.mean());
	    if (effP.GetEfficiency(iBin)<=0.) {
	      b.first.print();
	      std::cout << "[ERROR] Negative efficiency in: bin: " << iBin << ", pT: " << pt.mean() << ", rap: " << rap.mean() << ", eff: " << effP.GetEfficiency(iBin) << ", obj: " << o << std::endl;
	      return false;
	    }
	    eff[o][c.first][pd.first][t+"_MC"]["Val"][b.first] = effP.GetEfficiency(iBin);
	    eff[o][c.first][pd.first][t+"_MC"]["Err_Stat_High"][b.first] = effP.GetEfficiencyErrorUp(iBin);
	    eff[o][c.first][pd.first][t+"_MC"]["Err_Stat_Low"][b.first] = effP.GetEfficiencyErrorLow(iBin);
	    eff[o][c.first][pd.first][t+"_MC"]["Err_Syst_High"][b.first] = 0.0;
	    eff[o][c.first][pd.first][t+"_MC"]["Err_Syst_Low"][b.first] = 0.0;
	    eff[o][c.first][pd.first][t+"_MC"]["Err_Tot_High"][b.first] = effP.GetEfficiencyErrorUp(iBin);
	    eff[o][c.first][pd.first][t+"_MC"]["Err_Tot_Low"][b.first] = effP.GetEfficiencyErrorLow(iBin);
	  }
	}
      }
    }
  }
  return true;
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


bool correctRawYields(VarBinTriMap_t& inputVar, const BinSextaMap_t& effMap)
{
  //
  const std::string& accName = "Acceptance_MC";
  const std::string& effName = "Efficiency_MC";
  //
  const auto tmpVar = inputVar;
  for (const auto& o : tmpVar) {
    for (const auto& c : o.second) {
      for (const auto& pd : c.second) {
	for (const auto& b : pd.second) {
	  auto& var = inputVar.at(o.first).at(c.first).at(pd.first).at(b.first);
	  // Loop over the variables
	  for (const auto& v : b.second) {
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
	  std::vector<std::string> rVec;
	  for (const auto& v : b.second) { if (v.first.rfind("R_",0)==0) { rVec.push_back(v.first); } }
	  if (rVec.empty()) { for (const auto& v : b.second) { if (v.first.rfind("N_",0)==0) { rVec.push_back(v.first); } } }
	  if (rVec.empty()) { std::cout << "[ERROR] Info for RatioTo1S not found!" << std::endl; return false; }
	  for (const auto& v : rVec) {
	    if (v.find("To")==std::string::npos || v.rfind("N_Bkg",0)==0 || v.find("SS")!=std::string::npos) continue;
	    const auto& obj = v.substr(0, v.find("To")).substr(2);
	    if (obj.rfind("JPsi",0)==0 || obj.rfind("Ups1S",0)==0) continue;
	    // Initialize ratio to 1S
	    const auto& iVar = b.second.at(v);
	    auto& oVar = var[obj][c.first][pd.first]["RatioTo1S"];
	    // Copy value and errors of ratio to 1S
	    if (v.rfind("R_",0)==0) { for (const auto& p : iVar) { oVar[p.first][b.first] = p.second; } }
	    else if (v.rfind("N_",0)==0) {
	      auto r = v;
	      if (obj.rfind("Psi2S",0)==0) { stringReplace(r, "Psi2S", "JPsi"); }
	      else if (obj.rfind("Ups2S",0)==0) { stringReplace(r, "Ups2S", "Ups1S"); }
	      else if (obj.rfind("Ups3S",0)==0) { stringReplace(r, "Ups3S", "Ups1S"); }
	      const auto& rVar = b.second.at(r);
	      oVar["Val"][b.first] = iVar.at("Val")/rVar.at("Val");
	      oVar["Err_Stat_High"][b.first] = fabs(oVar.at("Val").at(b.first))*sumErrors({(iVar.at("Err_Stat_High")/iVar.at("Val")), (rVar.at("Err_Stat_High")/rVar.at("Val"))});
	      oVar["Err_Stat_Low"][b.first] = fabs(oVar.at("Val").at(b.first))*sumErrors({(iVar.at("Err_Stat_Low")/iVar.at("Val")), (rVar.at("Err_Stat_Low")/rVar.at("Val"))});
	      oVar["Err_Syst_High"][b.first] = fabs(oVar.at("Val").at(b.first))*sumErrors({(iVar.at("Err_Syst_High")/iVar.at("Val")), (rVar.at("Err_Syst_High")/rVar.at("Val"))});
	      oVar["Err_Syst_Low"][b.first] = fabs(oVar.at("Val").at(b.first))*sumErrors({(iVar.at("Err_Syst_Low")/iVar.at("Val")), (rVar.at("Err_Syst_Low")/rVar.at("Val"))});
	    }
	  }
	}
      }
    }
  }
  //
  return true;
};


double getCrossSectionValue(const double& N, const double& lumi, const double& binWidth)
{
  return (N / (lumi*binWidth));
};

  
double getCrossSectionError(const double& N, const double& lumi, const double& binWidth, const double& errN, const double& errLumi)
{
  const auto& common = getCrossSectionValue(N, lumi, binWidth);
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
	      if (t.first=="Val") { oVar[t.first][b.first] = getCrossSectionValue(t.second, lumi, 1.0); }
	      else if (t.first.rfind("Err_",0)==0) { oVar[t.first][b.first] = getCrossSectionError(iVar.at("Val"), lumi, 1.0, t.second, 0.0); }
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
