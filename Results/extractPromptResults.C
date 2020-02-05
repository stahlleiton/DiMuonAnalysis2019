#ifndef extractPromptResults_C
#define extractPromptResults_C

// Auxiliary Files
// ROOT headers
#include "TEfficiency.h"
// c++ headers
#include <iostream>
#include <string>


typedef std::map< binF        , std::string > StrBinMap_t;
typedef std::map< std::string , StrBinMap_t > StrBinDiMap_t;
typedef std::map< std::string , StrBinDiMap_t > StrBinTriMap_t;
typedef std::map< binF        , TEfficiency > EffMap_t;
typedef std::map< std::string , EffMap_t    > EffDiMap_t;
typedef std::map< std::string , EffDiMap_t  > EffTriMap_t;


bool extractEfficiency   ( EffTriMap_t& effM , const StrBinTriMap_t& effFiles );
bool extractPromptYields ( VarBinTriMap_t& nom , const VarBinTriMap_t& prCut , const VarBinTriMap_t& nonPrCut , const EffTriMap_t& effM );


bool extractPromptResults(
			  VarBinTriMap_t& inputVar_Nominal,
			  const VarBinTriMap_t& inputVar_Prompt,
			  const VarBinTriMap_t& inputVar_NonPrompt
			  )
{
  //
  // Extract the prompt and non-prompt decayLen-cut efficiencies
  StrBinTriMap_t effFile =
    {
     {"NT" ,
      {
       {"JPsiNoPR" , {{binF("Cand_AbsRap", 0.0, 1.4), "0_1.4_0.9NT_NonpromptJPsi"},
		      {binF("Cand_AbsRap", 1.4, 2.4), "1.4_2.4_0.9NT_NonpromptJPsi"}}},
       {"JPsiPR"   , {{binF("Cand_AbsRap", 0.0, 1.4), "0_1.4_0.9NT_PromptJPsi"},
		      {binF("Cand_AbsRap", 1.4, 2.4), "1.4_2.4_0.9NT_PromptJPsi"}}},
       {"Psi2SNoPR", {{binF("Cand_AbsRap", 0.0, 1.4), "0_1.4_0.9NT_NonpromptPsi2S"},
		      {binF("Cand_AbsRap", 1.4, 2.4), "1.4_2.4_0.9NT_NonpromptPsi2S"}}},
       {"Psi2SPR"  , {{binF("Cand_AbsRap", 0.0, 1.4), "0_1.4_0.9NT_PromptPsi2S"},
		      {binF("Cand_AbsRap", 1.4, 2.4), "1.4_2.4_0.9NT_PromptPsi2S"}}}
      }
     },
     {"WT" ,
      {
       {"JPsiNoPR" , {{binF("Cand_AbsRap", 0.0, 1.4), "0_1.4_0.9WT_NonpromptJPsi"},
		      {binF("Cand_AbsRap", 1.4, 2.4), "1.4_2.4_0.9WT_NonpromptJPsi"}}},
       {"JPsiPR"   , {{binF("Cand_AbsRap", 0.0, 1.4), "0_1.4_0.9WT_PromptJPsi"},
		      {binF("Cand_AbsRap", 1.4, 2.4), "1.4_2.4_0.9WT_PromptJPsi"}}},
       {"Psi2SNoPR", {{binF("Cand_AbsRap", 0.0, 1.4), "0_1.4_0.9WT_NonpromptPsi2S"},
		      {binF("Cand_AbsRap", 1.4, 2.4), "1.4_2.4_0.9WT_NonpromptPsi2S"}}},
       {"Psi2SPR"  , {{binF("Cand_AbsRap", 0.0, 1.4), "0_1.4_0.9WT_PromptPsi2S"},
		      {binF("Cand_AbsRap", 1.4, 2.4), "1.4_2.4_0.9WT_PromptPsi2S"}}}
      }
     }
    };
  
  EffTriMap_t eff;
  std::cout << "Extracting the Decay Length cut efficiency" << std::endl;
  if (!extractEfficiency(eff, effFile)) { return false; }
  //
  // Compute the prompt and non-prompt yields
  std::cout << "Extracting the prompt and non-prompt yields" << std::endl;
  if (!extractPromptYields(inputVar_Nominal, inputVar_Prompt, inputVar_NonPrompt, eff)) { return false; }
  //
  // return
  return true;
};


double derivePromptYield(const bool& doPr, const double& yPrCut, const double& yNonPrCut, const double& effPr, const double& effNonPr)
{
  const auto& effD = (effPr - effNonPr);
  if (doPr) { return ((1.0 - effNonPr)/effD)*yPrCut - (effNonPr/effD)*yNonPrCut; }
  return (effPr/effD)*yNonPrCut - ((1.0 - effPr)/effD)*yPrCut;
};


double derivePromptUnc(const bool& doPr, const double& yPrCut, const double& yNonPrCut, const double& effPr, const double& effNonPr,
		       const double& uPrCut, const double& uNonPrCut, const double& uEffPr, const double& uEffNonPr)
{
  const auto& effI = 1.0/(effPr - effNonPr);
  const auto& uEffD = sumErrors({uEffPr, uEffNonPr});
  const auto& uEffI = uEffD*effI*effI;
  if (doPr) {
    const auto& u1 = fabs((1.0 - effNonPr)*effI*yPrCut)*sumErrors({(uEffNonPr/(1.0-effNonPr)), (uEffI/effI), (uPrCut/yPrCut)});
    const auto& u2 = fabs(effNonPr*effI*yNonPrCut)*sumErrors({(uEffNonPr/effNonPr), (uEffI/effI), (uNonPrCut/yNonPrCut)});
    return sumErrors({u1, u2});
  }
  const auto& u1 = fabs(effPr*effI*yNonPrCut)*sumErrors({(uEffPr/effPr), (uEffI/effI), (uNonPrCut/yNonPrCut)});
  const auto& u2 = fabs((1.0 - effPr)*effI*yPrCut)*sumErrors({(uEffPr/(1.0-effPr)), (uEffI/effI), (uPrCut/yPrCut)});
  return sumErrors({u1, u2});
};


bool extractPromptYields(VarBinTriMap_t& nom, const VarBinTriMap_t& prCut, const VarBinTriMap_t& nonPrCut, const EffTriMap_t& effM)
{
  for (const auto& o : prCut) {
    for (const auto& c : o.second) {
      for (const auto& pd : c.second) {
	for (const auto& b : pd.second) {
	  // Loop over the variables
	  for (const auto& v : b.second) {
	    const auto& vN = v.first;
	    const auto& obj = vN.substr(0, vN.find("To")).substr(2);
	    if (vN.rfind("N_",0)!=0 || vN.find("To")==std::string::npos) continue;
	    if (vN.rfind("N_Bkg",0)==0 || vN.find("SS")!=std::string::npos) continue;
	    // Extract the efficiencies
	    const auto& rap = b.first.getbin("Cand_Rap");
	    const auto& pT = b.first.getbin("Cand_Pt");
	    const binF absrap("Cand_AbsRap", std::min(fabs(rap.low()), fabs(rap.high())), std::max(fabs(rap.low()), fabs(rap.high())));
	    const std::string& pdT = (pd.first.rfind("MUON")!=std::string::npos ? "WT" : "NT");
	    if (!contain(effM, pdT)) { std::cout << "[ERROR] Efficiency for " << pdT << " was not found!" << std::endl; return false; }
	    if (!contain(effM.at(pdT), obj+"PR") && !contain(effM.at(pdT), obj+"NoPR")) { std::cout << "[ERROR] Efficiency for " << obj << " was not found!" << std::endl; return false; }
	    if (!contain(effM.at(pdT).at(obj+"PR"), absrap)) { std::cout << "[ERROR] Efficiency does not contain bin "; absrap.print(); return false; }
	    const auto& effPr = effM.at(pdT).at(obj+"PR").at(absrap);
	    const auto& effNonPr = effM.at(pdT).at(obj+"NoPR").at(absrap);
	    const auto& iBinPr = effPr.FindFixBin(pT.mean());
	    const auto& iBinNonPr = effNonPr.FindFixBin(pT.mean());
	    // Extract the yields with decayLen cuts
	    if (!contain(nonPrCut, o.first)) { std::cout << "[ERROR] NonPR does not contain: " << o.first << std::endl; return false; }
	    if (!contain(nonPrCut.at(o.first), c.first)) { std::cout << "[ERROR] NonPR does not contain: " << o.first << std::endl; return false; }
	    if (!contain(nonPrCut.at(o.first).at(c.first), pd.first)) { std::cout << "[ERROR] NonPR does not contain: " << pd.first << std::endl; return false; }
	    if (!contain(nonPrCut.at(o.first).at(c.first).at(pd.first), b.first)) { std::cout << "[ERROR] NonPR does not contain "; b.first.print(); return false; }
	    if (!contain(nonPrCut.at(o.first).at(c.first).at(pd.first).at(b.first), vN)) { std::cout << "[ERROR] NonPR does not contain: " << vN << std::endl; return false; }
	    const auto& yPrCut = prCut.at(o.first).at(c.first).at(pd.first).at(b.first).at(vN);
	    const auto& yNonPrCut = nonPrCut.at(o.first).at(c.first).at(pd.first).at(b.first).at(vN);
	    // Define the prompt and non-prompt yields
	    auto vPr = vN; stringReplace(vPr, obj, (obj+"PR"));
	    auto vNonPr = vN; stringReplace(vNonPr, obj, (obj+"NoPR"));
	    auto& yPr = nom[o.first][c.first][pd.first][b.first][vPr];
	    auto& yNonPr = nom[o.first][c.first][pd.first][b.first][vNonPr];
	    // Compute the prompt and non-prompt yield
	    yPr["Val"] = derivePromptYield(true, yPrCut.at("Val"), yNonPrCut.at("Val"), effPr.GetEfficiency(iBinPr), effNonPr.GetEfficiency(iBinNonPr));
	    yNonPr["Val"] = derivePromptYield(false, yPrCut.at("Val"), yNonPrCut.at("Val"), effPr.GetEfficiency(iBinPr), effNonPr.GetEfficiency(iBinNonPr));
	    // Compute the prompt and non-prompt uncertainties
	    yPr["Err_Stat_Low"] = derivePromptUnc(true, yPrCut.at("Val"), yNonPrCut.at("Val"), effPr.GetEfficiency(iBinPr), effNonPr.GetEfficiency(iBinNonPr), 
						  yPrCut.at("Err_Stat_Low"), yNonPrCut.at("Err_Stat_Low"), 0.0, 0.0);
	    yPr["Err_Syst_Low"] = derivePromptUnc(true, yPrCut.at("Val"), yNonPrCut.at("Val"), effPr.GetEfficiency(iBinPr), effNonPr.GetEfficiency(iBinNonPr), 
						  0.0, 0.0, effPr.GetEfficiencyErrorLow(iBinPr), effNonPr.GetEfficiencyErrorLow(iBinNonPr));
	    yPr["Err_Stat_High"] = derivePromptUnc(true, yPrCut.at("Val"), yNonPrCut.at("Val"), effPr.GetEfficiency(iBinPr), effNonPr.GetEfficiency(iBinNonPr), 
						  yPrCut.at("Err_Stat_High"), yNonPrCut.at("Err_Stat_High"), 0.0, 0.0);
	    yPr["Err_Syst_High"] = derivePromptUnc(true, yPrCut.at("Val"), yNonPrCut.at("Val"), effPr.GetEfficiency(iBinPr), effNonPr.GetEfficiency(iBinNonPr), 
	    					  0.0, 0.0, effPr.GetEfficiencyErrorLow(iBinPr), effNonPr.GetEfficiencyErrorLow(iBinNonPr));
	    yNonPr["Err_Stat_Low"] = derivePromptUnc(false, yPrCut.at("Val"), yNonPrCut.at("Val"), effPr.GetEfficiency(iBinPr), effNonPr.GetEfficiency(iBinNonPr), 
						  yPrCut.at("Err_Stat_Low"), yNonPrCut.at("Err_Stat_Low"), 0.0, 0.0);
	    yNonPr["Err_Syst_Low"] = derivePromptUnc(false, yPrCut.at("Val"), yNonPrCut.at("Val"), effPr.GetEfficiency(iBinPr), effNonPr.GetEfficiency(iBinNonPr), 
						  0.0, 0.0, effPr.GetEfficiencyErrorLow(iBinPr), effNonPr.GetEfficiencyErrorLow(iBinNonPr));
	    yNonPr["Err_Stat_High"] = derivePromptUnc(false, yPrCut.at("Val"), yNonPrCut.at("Val"), effPr.GetEfficiency(iBinPr), effNonPr.GetEfficiency(iBinNonPr), 
						  yPrCut.at("Err_Stat_High"), yNonPrCut.at("Err_Stat_High"), 0.0, 0.0);
	    yNonPr["Err_Syst_High"] = derivePromptUnc(false, yPrCut.at("Val"), yNonPrCut.at("Val"), effPr.GetEfficiency(iBinPr), effNonPr.GetEfficiency(iBinNonPr), 
						  0.0, 0.0, effPr.GetEfficiencyErrorLow(iBinPr), effNonPr.GetEfficiencyErrorLow(iBinNonPr));
	  }
	}
      }
    }
  }
  return true;
};


bool extractEfficiency(EffTriMap_t& effM, const StrBinTriMap_t& effFiles)
{
  for (const auto& t : effFiles) {
    for (const auto& o : t.second) {
      for (const auto& fN : o.second) {
	const auto& bin = fN.first;
	const auto& fileName = (fN.second+".root");
	const std::string& dirName = "Efficiency/DecayLenCut/";
	TFile file((dirName+fileName).c_str(), "READ");
	if (!file.IsOpen() || file.IsZombie()) { std::cout << "[ERROR] File " << (dirName+fileName) << " was not found!" << std::endl; return false; }
	const auto& effName = (fN.second.substr(fN.second.rfind("t")+1)+"_Efficiency");
	const auto& effP = dynamic_cast<TEfficiency*>(file.Get(effName.c_str()));
	if (!effP) { std::cout << "[ERROR] Efficiency " << effName << " in " << (dirName+fileName) << " was not found!" << std::endl; file.Close(); return false; }
	effM[t.first][o.first][bin] = *effP;
	file.Close();
      }
    }
  }
  return true;
};


#endif // #ifndef extractPromptResults_C
