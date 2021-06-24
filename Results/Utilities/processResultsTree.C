#ifndef processResultsTree_C
#define processResultsTree_C

// Auxiliary Files
#include "extractResultsTree.C"
#include "resultsUtils.h"
#include "../../Efficiency/correctEfficiency.C"
// ROOT headers
#include "TEfficiency.h"
// c++ headers
#include <iostream>
#include <string>


bool addDecayCutYield            ( VarBinTriMap_t& inputVar , const VarBinTriMap_t& inputVar_Cut   , const VarBinTriMap_t& inputVar_Inc );
bool extractPromptYields         ( VarBinTriMap_t& inputVar , const VarBinTriMap_t& inputVar_PrCut , const VarBinTriMap_t& inputVar_NonPrCut , const BinSextaMap_t& effM , const std::string& effType );
bool getEfficiency               ( BinSextaMap_t& eff , const VarBinTriMap_t& inputVar , const std::string& effType , const std::string& effDir , const bool& isNominal=true );
bool correctRawYields            ( VarBinTriMap_t& inputVar , const BinSextaMap_t& effMap , const bool& isNominal=true , const std::string& accType="MC" , const std::string& effType="TnP");
bool computeCrossSection         ( BinSextaMap_t& var , const VarBinTriMap_t& inputVar , const bool& doSyst , const VarBinTriMap_t& nomVar );
bool computeRatioTo1S            ( BinSextaMap_t& var , const VarBinTriMap_t& inputVar , const bool& doSyst , const VarBinTriMap_t& nomVar , const IntMap_t& systCorr );
bool computeForwardBackwardRatio ( BinSextaMap_t& var , const VarBinTriMap_t& inputVar , const bool& doSyst , const VarBinTriMap_t& nomVar , const IntMap_t& systCorr );


bool processResultsTree(
			BinSextaMap_t& var,
			BinSeptaMapVec_t& systVar,
			VarBinTriMap_t& inputVar,
			const bool& doSyst,
			const bool& isNominal,
			const VarBinTriMap_t& nomVar,
			const IntMap_t& systCorr,
			const std::string& effDir
			)
{
  //
  // Extract the Acceptance and Efficiency container
  BinSextaMap_t eff;
  std::cout << "[INFO] Extracting the acceptance" << std::endl;
  if (!getEfficiency(eff, inputVar, "Acceptance", effDir)) { return false; }
  std::cout << "[INFO] Extracting the total efficiency" << std::endl;
  if (!getEfficiency(eff, inputVar, "Efficiency_Total", effDir)) { return false; }
  //
  // Proceed to correct the Raw Yields
  std::cout << "[INFO] Correcting the raw yields" << std::endl;
  if (!correctRawYields(inputVar, eff, isNominal)) { return false; }
  //
  // Compute the cross section
  std::cout << "[INFO] Computing cross sections" << std::endl;
  if (!computeCrossSection(var, inputVar, doSyst, nomVar)) { return false; }
  //
  // Compute the ratio to 1S
  std::cout << "[INFO] Computing ratios to 1S" << std::endl;
  if (!computeRatioTo1S(var, inputVar, doSyst, nomVar, systCorr)) { return false; }
  //
  // Compute the Forward-Backward ratio
  std::cout << "[INFO] Computing forward-backward ratios" << std::endl;
  if (!computeForwardBackwardRatio(var, inputVar, doSyst, nomVar, systCorr)) { return false; }
  //
  return true;
};


bool processDecayCutResults(VarBinQuadMap_t& inputVar, const std::string& effLbl, const std::string& effDir)
{
  // Check decay cut results. If missing, add using inclusive results
  if (!contain(inputVar, "NonPrompt_Cut") || !contain(inputVar, "Prompt_Cut")) {
    if (!contain(inputVar, "Inclusive")) { std::logic_error("[ERROR] Missing inclusive result"); }
    if (contain(inputVar, "NonPrompt_Cut")) {
      if (!addDecayCutYield(inputVar["Prompt_Cut"], inputVar.at("NonPrompt_Cut"), inputVar.at("Inclusive"))) { return false; }
    }
    else if (contain(inputVar, "Prompt_Cut")) {
      if (!addDecayCutYield(inputVar["NonPrompt_Cut"], inputVar.at("Prompt_Cut"), inputVar.at("Inclusive"))) { return false; }
    }
    else { std::logic_error("[ERROR] Missing decay cut results"); }
  }
  // Process decay cut results
  auto& inputVar_Inc = inputVar["Inclusive"];
  const auto& inputVar_PrCut = inputVar.at("Prompt_Cut");
  const auto& inputVar_NonPrCut = inputVar.at("NonPrompt_Cut");
  //
  // Extract the Decay Cut Efficiency container
  BinSextaMap_t eff;
  std::cout << "[INFO] Extracting the decay cut efficiency for: " << effLbl << std::endl;
  if (!getEfficiency(eff, inputVar_PrCut, effLbl, effDir)) { return false; }
  //
  // Compute the prompt and non-prompt yields
  std::cout << "[INFO] Extracting the prompt and non-prompt yields" << std::endl;
  if (!extractPromptYields(inputVar_Inc, inputVar_PrCut, inputVar_NonPrCut, eff, effLbl+"_TnP")) { return false; }
  //
  return true; 
};


bool addDecayCutYield(VarBinTriMap_t& inputVar, const VarBinTriMap_t& inputVar_Cut, const VarBinTriMap_t& inputVar_Inc)
{
  for (const auto& o : inputVar_Cut) {
    for (const auto& c : o.second) {
      for (const auto& pd : c.second) {
	for (const auto& b : pd.second) {
	  // Loop over the variables
	  for (const auto& v : b.second) {
	    //
	    const auto& vN  = v.first;
	    if (vN=="Luminosity") { inputVar[o.first][c.first][pd.first][b.first][vN] = v.second; }
	    if (vN.rfind("N_",0)!=0 || vN.find("To")==std::string::npos) continue;
	    if (vN.rfind("N_Bkg",0)==0 || vN.find("SS")!=std::string::npos) continue;
	    //
	    // Extract the yields with decayLen cuts
	    if (!contain(inputVar_Inc, o.first)) { std::cout << "[ERROR] inputVar_Inc does not contain: " << o.first << std::endl; return false; }
	    if (!contain(inputVar_Inc.at(o.first), c.first)) { std::cout << "[ERROR] inputVar_Inc does not contain: " << c.first << std::endl; return false; }
	    if (!contain(inputVar_Inc.at(o.first).at(c.first), pd.first)) { std::cout << "[ERROR] inputVar_Inc does not contain: " << pd.first << std::endl; return false; }
	    if (!contain(inputVar_Inc.at(o.first).at(c.first).at(pd.first), b.first)) { std::cout << "[ERROR] inputVar_Inc does not contain "; b.first.print(); return false; }
	    if (!contain(inputVar_Inc.at(o.first).at(c.first).at(pd.first).at(b.first), vN)) { std::cout << "[ERROR] inputVar_Inc does not contain: " << vN << std::endl; b.first.print(); return false; }
	    //
	    auto& yield = inputVar[o.first][c.first][pd.first][b.first][vN];
	    const auto& yield_Inc = inputVar_Inc.at(o.first).at(c.first).at(pd.first).at(b.first).at(vN);
	    const auto& yield_Cut = inputVar_Cut.at(o.first).at(c.first).at(pd.first).at(b.first).at(vN);
	    //
	    // Add the yield
	    yield["Val"] = yield_Inc.at("Val") - yield_Cut.at("Val");
	    //
	    // Add the uncertainties
	    for (const auto& u : StringVector_t({"Err_Stat_Low", "Err_Stat_High", "Err_Syst_Low", "Err_Syst_High"})) {
	      yield[u] = std::sqrt(std::abs(yield_Inc.at(u)*yield_Inc.at(u) - yield_Cut.at(u)*yield_Cut.at(u)));
	      if (yield_Inc.at(u)<yield_Cut.at(u)) {
		b.first.print();
		//throw std::logic_error(Form("[ERROR] Invaid uncertainty for inclusive (%.1f +- %.1f) and cut (%.1f +- %.1f)", yield_Inc.at("Val"), yield_Inc.at(u), yield_Cut.at("Val"), yield_Cut.at(u)));
		//std::cout << Form("[WARNING] Invalid uncertainty for inclusive (%.1f +- %.1f) - cut (%.1f +- %.1f) = (%.1f +- %.1f)", yield_Inc.at("Val"), yield_Inc.at(u), yield_Cut.at("Val"), yield_Cut.at(u), yield.at("Val"), yield.at(u)) << std::endl;
	      }
	    }
	    // Check the yield information
	    for (const auto& v : yield) {
	      if (v.second<0. || isnan(v.second)) {
		std::cout << yield_Inc.at(v.first) << " " << yield_Cut.at(v.first) << std::endl;
		std::cout << "[ERROR] Invalid decayLen cut yield " << v.first << " ( " << v.second << " )  in pd: " << pd.first << " , obj: " << o.first << std::endl; b.first.print(); return false;
	      }
	    }
	  }
	}
      }
    }
  }
  return true;
};


DoublePair_t derivePromptYield(const double& yield_PrCut, const double& yield_NonPrCut,
			       const double& eff_PrCut_Prompt, const double& eff_PrCut_NonPrompt)
{
  const auto eff_PrCut_Diff  = (eff_PrCut_Prompt - eff_PrCut_NonPrompt);
  // Check inputs
  if (yield_PrCut<0.0 || isnan(yield_PrCut)) { throw std::logic_error(Form("[ERROR] Invalid PrCut yield value ( %.2f )", yield_PrCut)); }
  if (yield_NonPrCut<-300.0 || isnan(yield_NonPrCut)) { throw std::logic_error(Form("[ERROR] Invalid NonPrCut yield value ( %.2f )", yield_NonPrCut)); }
  if (eff_PrCut_Prompt<=0.0 || isnan(eff_PrCut_Prompt)) { throw std::logic_error(Form("[ERROR] Invalid PrCut_Prompt efficiency ( %.2f )", eff_PrCut_Prompt)); }
  if (eff_PrCut_NonPrompt<=0.0 || isnan(eff_PrCut_NonPrompt)) { throw std::logic_error(Form("[ERROR] Invalid PrCut_NonPrompt efficiency ( %.2f )", eff_PrCut_NonPrompt)); }
  if (eff_PrCut_Diff<=0.0 || isnan(eff_PrCut_Diff)) { throw std::logic_error(Form("[ERROR] Invalid PrCut_Diff efficiency ( %.2f )", eff_PrCut_Diff)); }
  // Compute prompt and non-prompt yields 
  const auto yield_Prompt    = divide(+yield_PrCut - eff_PrCut_NonPrompt*(yield_NonPrCut + yield_PrCut), eff_PrCut_Diff);
  const auto yield_NonPrompt = divide(-yield_PrCut + eff_PrCut_Prompt   *(yield_NonPrCut + yield_PrCut), eff_PrCut_Diff);
  // Check prompt and non-prompt yields
  if (yield_Prompt<0.0 || isnan(yield_Prompt)) { throw std::logic_error(Form("[ERROR] Invalid prompt yield value ( %g )", yield_Prompt)); }
  if (yield_NonPrompt<-300.0 || isnan(yield_NonPrompt)) { throw std::logic_error(Form("[ERROR] Invalid non-prompt yield value ( %g )", yield_NonPrompt)); }
  // Return results
  return std::make_pair(yield_Prompt, yield_NonPrompt);
};


DoublePair_t derivePromptUnc(const double& yield_PrCut, const double& yield_NonPrCut,
			     const double& eff_PrCut_Prompt, const double& eff_PrCut_NonPrompt,
			     const double& unc_yield_PrCut, const double& unc_yield_NonPrCut,
			     const double& unc_eff_PrCut_Prompt, const double& unc_eff_PrCut_NonPrompt,
			     const bool& isCorrelatedYield=false)
{
  // Check input uncertainties
  if (unc_yield_PrCut<0.0 || isnan(unc_yield_PrCut)) { throw std::logic_error(Form("[ERROR] Invalid PrCut yield uncertainty ( %.2f )", unc_yield_PrCut)); }
  if (unc_yield_NonPrCut<0.0 || isnan(unc_yield_NonPrCut)) { throw std::logic_error(Form("[ERROR] Invalid NonPrCut yield uncertainty ( %.2f )", unc_yield_NonPrCut)); }
  if (unc_yield_PrCut!=unc_yield_NonPrCut &&
      (unc_yield_PrCut==0. || unc_yield_NonPrCut==0.)) { throw std::logic_error(Form("[ERROR] Invalid PrCut yeild uncertainty ( %.2f , %.2f )", unc_yield_PrCut, unc_yield_NonPrCut)); }
  if (unc_eff_PrCut_Prompt<0.0 || isnan(unc_eff_PrCut_Prompt)) { throw std::logic_error(Form("[ERROR] Invalid PrCut_Prompt efficiency uncertainty ( %.2f )", unc_eff_PrCut_Prompt)); }
  if (unc_eff_PrCut_NonPrompt<0.0 || isnan(unc_eff_PrCut_NonPrompt)) { throw std::logic_error(Form("[ERROR] Invalid PrCut_NonPrompt efficiency uncertainty ( %.2f )", unc_eff_PrCut_NonPrompt)); }
  if (unc_eff_PrCut_NonPrompt!=unc_eff_PrCut_Prompt &&
      (unc_eff_PrCut_NonPrompt==0. || unc_eff_PrCut_Prompt==0.)) { throw std::logic_error(Form("[ERROR] Invalid PrCut efficiency uncertainty ( %.2f , %.2f )", unc_eff_PrCut_Prompt, unc_eff_PrCut_NonPrompt)); }
  // Compute prompt and non-prompt yield values
  const auto value = derivePromptYield(yield_PrCut, yield_NonPrCut, eff_PrCut_Prompt, eff_PrCut_NonPrompt);
  const auto eff_PrCut_Diff  = (eff_PrCut_Prompt - eff_PrCut_NonPrompt);
  // Compute prompt yield uncertainty parts
  const auto unc_part1_Prompt = unc_yield_PrCut         * divide(1.0 - eff_PrCut_NonPrompt, eff_PrCut_Diff);
  const auto unc_part2_Prompt = unc_yield_NonPrCut      * divide(-eff_PrCut_NonPrompt, eff_PrCut_Diff);
  const auto unc_part3_Prompt = unc_eff_PrCut_Prompt    * divide(-value.first, eff_PrCut_Diff);
  const auto unc_part4_Prompt = unc_eff_PrCut_NonPrompt * divide(-value.second, eff_PrCut_Diff);
  // Compute non-prompt yield uncertainty parts
  const auto unc_part1_NonPrompt = unc_yield_PrCut         * divide(-1.0 + eff_PrCut_Prompt, eff_PrCut_Diff);
  const auto unc_part2_NonPrompt = unc_yield_NonPrCut      * divide(eff_PrCut_Prompt, eff_PrCut_Diff);
  const auto unc_part3_NonPrompt = unc_eff_PrCut_Prompt    * divide(value.first, eff_PrCut_Diff);
  const auto unc_part4_NonPrompt = unc_eff_PrCut_NonPrompt * divide(value.second, eff_PrCut_Diff);
  // Compute correlation matrix
  const auto corrM = (isCorrelatedYield ? TMatrixD(4,4,DoubleVec_t({1,1,0,0,1,1,0,0,0,0,1,0,0,0,0,1}).data()) : TMatrixD());
  // Compute prompt and non-prompt uncertainties
  const auto unc_yield_Prompt    = sumErrors({unc_part1_Prompt   , unc_part2_Prompt   , unc_part3_Prompt   , unc_part4_Prompt   }, corrM);
  const auto unc_yield_NonPrompt = sumErrors({unc_part1_NonPrompt, unc_part2_NonPrompt, unc_part3_NonPrompt, unc_part4_NonPrompt}, corrM);
  // Check prompt and non-prompt uncertainties
  if (unc_yield_Prompt<0. || isnan(unc_yield_Prompt)) { throw std::logic_error(Form("[ERROR] Invalid prompt yield uncertainty ( %.2f )", unc_yield_Prompt)); }
  if (unc_yield_NonPrompt<0. || isnan(unc_yield_NonPrompt)) { throw std::logic_error(Form("[ERROR] Invalid non-prompt yield uncertainty ( %.2f )", unc_yield_NonPrompt)); }
  if (unc_yield_Prompt!=unc_yield_NonPrompt &&
      (unc_yield_Prompt==0. || unc_yield_NonPrompt==0.)) { throw std::logic_error(Form("[ERROR] Invalid prompt or non-prompt yield uncertainties ( %.2f , %.2f )", unc_yield_Prompt, unc_yield_NonPrompt)); }
  // Return results
  return std::make_pair(unc_yield_Prompt, unc_yield_NonPrompt);
};


bool extractPromptYields(VarBinTriMap_t& inputVar, const VarBinTriMap_t& inputVar_PrCut, const VarBinTriMap_t& inputVar_NonPrCut,
			 const BinSextaMap_t& effM, const std::string& effType)
{
  for (const auto& o : inputVar_PrCut) {
    for (const auto& c : o.second) {
      for (const auto& pd : c.second) {
	for (const auto& b : pd.second) {
	  // Loop over the variables
	  for (const auto& v : b.second) {
	    //
	    const auto& vN  = v.first;
	    if (vN=="Luminosity") { inputVar[o.first][c.first][pd.first][b.first][vN] = v.second; }
	    if (vN.rfind("N_",0)!=0 || vN.find("To")==std::string::npos) continue;
	    if (vN.rfind("N_Bkg",0)==0 || vN.find("SS")!=std::string::npos) continue;
	    //
	    const auto& obj = vN.substr(0, vN.find("To")).substr(2);
	    const auto& objPR = obj+"PR";
	    const auto& objNoPR = obj+"NoPR";
	    //
	    // Extract the efficiencies
	    const auto& col = c.first;
	    const auto& PD = pd.first;
	    const auto& anaBin = b.first;
	    for (const auto& oo : StringVector_t({objPR, objNoPR})) {
	      if (!contain(effM, oo)) { std::cout << "[ERROR] Efficiency_DecayCut for " << oo << " was not found!" << std::endl; return false; }
	      if (!contain(effM.at(oo), col)) { std::cout << "[ERROR] Efficiency_DecayCut of " << oo << " for " << col << " was not found!" << std::endl; return false; }
	      if (!contain(effM.at(oo).at(col), PD)) { std::cout << "[ERROR] Efficiency_DecayCut of " << oo << " , " << col << " for " << PD << " was not found!" << std::endl; return false; }
	      if (!contain(effM.at(oo).at(col).at(PD), effType)) { std::cout << "[ERROR] " << effType << " of " << oo << " , " << col << " , " << PD << " was not found!" << std::endl; return false; }
	      if (!contain(effM.at(oo).at(col).at(PD).at(effType), "Val")) { std::cout << "[ERROR] " << effType << " of " << oo << " , " << col << " , " << PD << " has no values!" << std::endl; return false; }
	      if (!contain(effM.at(oo).at(col).at(PD).at(effType).at("Val"), anaBin)) {
		std::cout << "[ERROR] " << effType << " of " << oo << " , " << col << " , " << PD  << " does not contain "; anaBin.print();
		for (const auto& bb : effM.at(oo).at(col).at(PD).at(effType).at("Val")) { bb.first.print(); }
		return false;
	      }
	    }
	    const auto& effM_PrCut_Prompt = effM.at(objPR).at(col).at(PD).at(effType);
	    const auto& effM_PrCut_NonPrompt = effM.at(objNoPR).at(col).at(PD).at(effType);
	    //
	    DoubleMap_t eff_PrCut_Prompt , eff_PrCut_NonPrompt;
	    for (const auto& t : effM_PrCut_Prompt  ) { eff_PrCut_Prompt[t.first] = t.second.at(b.first); }
	    for (const auto& t : effM_PrCut_NonPrompt) { eff_PrCut_NonPrompt[t.first] = t.second.at(b.first); }
	    if (eff_PrCut_Prompt.at("Val")==0.0) { std::cout << "[ERROR] Null efficiency for: " << objPR << " at "; anaBin.print(); return false; }
	    if (eff_PrCut_NonPrompt.at("Val")==0.0) { std::cout << "[ERROR] Null efficiency for: " << objNoPR << " at "; anaBin.print(); return false; }
	    //
	    // Extract the yields with decayLen cuts
	    if (!contain(inputVar_NonPrCut, o.first)) { std::cout << "[ERROR] inputVar_NonPrCut does not contain: " << o.first << std::endl; return false; }
	    if (!contain(inputVar_NonPrCut.at(o.first), col)) { std::cout << "[ERROR] inputVar_NonPrCut does not contain: " << col << std::endl; return false; }
	    if (!contain(inputVar_NonPrCut.at(o.first).at(col), PD)) { std::cout << "[ERROR] inputVar_NonPrCut does not contain: " << PD << std::endl; return false; }
	    if (!contain(inputVar_NonPrCut.at(o.first).at(col).at(PD), anaBin)) { std::cout << "[ERROR] inputVar_NonPrCut does not contain "; anaBin.print();
	      for (const auto& pp : inputVar_NonPrCut.at(o.first).at(col).at(PD)) { pp.first.print(); }
	      return false; }
	    if (!contain(inputVar_NonPrCut.at(o.first).at(col).at(PD).at(anaBin), vN)) {
	      std::cout << "[ERROR] inputVar_NonPrCut does not contain: " << vN << std::endl;
	      anaBin.print();
	      std::cout << o.first << "  " << col << "  " << PD << std::endl;
	      for (const auto& pp : inputVar_NonPrCut.at(o.first).at(col).at(PD).at(anaBin)) { std::cout << pp.first << std::endl; }
	      return false;
	    }
	    //
	    const auto& yield_PrCut = inputVar_PrCut.at(o.first).at(col).at(PD).at(anaBin).at(vN);
	    const auto& yield_NonPrCut = inputVar_NonPrCut.at(o.first).at(col).at(PD).at(anaBin).at(vN);
	    //
	    // Define the prompt and non-prompt yields
	    auto vPr = vN; stringReplace(vPr, obj, objPR);
	    auto vNonPr = vN; stringReplace(vNonPr, obj, objNoPR);
	    auto& yield_Prompt = inputVar[o.first][c.first][pd.first][b.first][vPr];
	    auto& yield_NonPrompt = inputVar[o.first][c.first][pd.first][b.first][vNonPr];
	    //
	    // Compute the prompt and non-prompt yields
	    const auto yields = derivePromptYield(yield_PrCut.at("Val"), yield_NonPrCut.at("Val"), eff_PrCut_Prompt.at("Val"), eff_PrCut_NonPrompt.at("Val"));
	    yield_Prompt["Val"] = yields.first;
	    yield_NonPrompt["Val"] = yields.second;
	    //
	    // Compute the uncertainties
	    for (const auto& u : StringVector_t({"Err_Stat_Low", "Err_Stat_High", "Err_Syst_Low", "Err_Syst_High"})) {
	      //
	      const auto unc_eff_Prompt = (u=="Err_Syst_Low" ? eff_PrCut_Prompt.at("Err_Tot_Low") : (u=="Err_Syst_High" ? eff_PrCut_Prompt.at("Err_Tot_High") : 0.0)); 
	      const auto unc_eff_NonPrompt = (u=="Err_Syst_Low" ? eff_PrCut_NonPrompt.at("Err_Tot_Low") : (u=="Err_Syst_High" ? eff_PrCut_NonPrompt.at("Err_Tot_High") : 0.0));
	      //
	      const auto unc = derivePromptUnc(yield_PrCut.at("Val"), yield_NonPrCut.at("Val"),
					       eff_PrCut_Prompt.at("Val"), eff_PrCut_NonPrompt.at("Val"),
					       yield_PrCut.at(u), yield_NonPrCut.at(u),
					       unc_eff_Prompt, unc_eff_NonPrompt);
	      //
	      yield_Prompt[u] = unc.first;
	      yield_NonPrompt[u] = unc.second;
	    }
	  }
	}
      }
    }
  }
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


std::map<std::string, BinSextaMap_t> EFF_; 
bool getEfficiency(BinSextaMap_t& eff, const VarBinTriMap_t& inputVar, const std::string& effType, const std::string& effDir, const bool& isNominal)
{
  if (contain(EFF_, effType)) { eff = EFF_.at(effType); return true; }
  //
  const StringVector_t lbl = {"Val", "Err_Stat_High", "Err_Stat_Low", "Err_Syst_High", "Err_Syst_Low", "Err_Tot_High", "Err_Tot_Low"};
  const StringVector_t effObj = { "JPsiNoPR", "JPsiPR", "Psi2SNoPR", "Psi2SPR" };
  //
  std::map<std::string, EffMap_t> eff1D;
  std::map<std::string, Unc1DMap_t> unc1D;
  for (const auto& o : effObj) {
    for (const auto& pd : inputVar.begin()->second.begin()->second) {
      auto PD = (pd.first.rfind("HIGHMULT",0)==0 ? "MINBIAS" : pd.first);
      const auto fileName = "../Efficiency/Output/"+effDir+"/" + PD + "/effContainer_MC_" + o + ".root";
      if (!getEffObjectsFromFile(eff1D[pd.first], unc1D[pd.first], fileName)) { return false; }
    }
  }
  //
  for (const auto& o : effObj) {
    for (const auto& c : inputVar.begin()->second) {
      for (const auto& pd : c.second) {
	//
	const auto& PD = pd.first;
	const auto& sample = "MC_" + o;
	const auto& col = c.first;
	//
	if (!contain(eff1D, PD)) { std::cout << "[ERROR] The extracted efficiency does not contain: " << PD << std::endl; return false; }
	if (!contain(eff1D.at(PD), sample)) { std::cout << "[ERROR] The extracted efficiency does not contain: " << sample << std::endl; return false; }
	if (!contain(eff1D.at(PD).at(sample), col)) { std::cout << "[ERROR] The extracted efficiency does not contain: " << col << std::endl; return false; }
	if (!contain(eff1D.at(PD).at(sample).at(col), effType)) { std::cout << "[ERROR] The extracted efficiency does not contain: " << effType << std::endl; return false; }
	const auto& effM = eff1D.at(PD).at(sample).at(col).at(effType);
	//
	std::map< AnaBinPair_t , Unc1DVec_t > uncM;
	if (effType!="Acceptance") {
	  if (!contain(unc1D, PD)) { std::cout << "[ERROR] The extracted efficiency uncertainty does not contain: " << PD << std::endl; return false; }
	  if (!contain(unc1D.at(PD), sample)) { std::cout << "[ERROR] The extracted efficiency uncertainty does not contain: " << sample << std::endl; return false; }
	  if (!contain(unc1D.at(PD).at(sample), col)) { std::cout << "[ERROR] The extracted efficiency uncertainty does not contain: " << col << std::endl; return false; }
	  if (!contain(unc1D.at(PD).at(sample).at(col), effType)) { std::cout << "[ERROR] The extracted efficiency uncertainty does not contain: " << effType << std::endl; return false; }
	  uncM = unc1D.at(PD).at(sample).at(col).at(effType);
	}
	//
	// Extract efficiency
	for (const auto& b : effM) {
	  const auto& bin = b.first;
	  auto eBin = bin;
	  if (eBin.first.getbin("NTrack").low()>150.) { eBin.first.setbin("NTrack", 150., 185.); } // HARDCODED WARNING
	  if (!contain(effM, eBin)) { eBin.first.print(); throw std::logic_error("[ERROR] Efficiency container does not cotain bin"); }
	  const auto& effV = effM.at(eBin);
	  const auto& effB = bin.first;
	  AnaBin_t anaB(effB);
	  const auto& varN = bin.second;
	  const auto& xVar = varN.substr(0, varN.rfind("_"));
	  if (!contain(effV, "NoCorr") || effV.at("NoCorr").empty()) { std::cout << "[ERROR] The extracted " << effType << " does not contain NoCorr" << std::endl; return false; }
	  const auto& hTotal = effV.at("NoCorr")[0].GetTotalHistogram();
	  for (int iBin=1; iBin<=hTotal->GetNbinsX(); iBin++) {
	    const double xMin = hTotal->GetXaxis()->GetBinLowEdge(iBin);
	    const double xMax = hTotal->GetXaxis()->GetBinUpEdge(iBin);
	    auto anaBin = anaB;
	    anaBin.setbin(xVar, xMin, xMax);
	    //
	    // Consider only bins that will be used
	    if (!contain(pd.second, anaBin)) continue;
	    //
	    for (const auto& co : StringMap_t({{"MC", "NoCorr"}, {"TnP", "TnP_Nominal"}})) {
	      if (contain(effV, co.second)) {
		if (effV.at(co.second).empty()) { std::cout << "[ERROR] The extracted " << effType << " does not contain " << co.second << std::endl; return false; }
		const auto& effP = effV.at(co.second)[0];
		const std::string effLbl = (effType=="Efficiency_Total" ? "Efficiency" : effType);
		const auto& type = effLbl + "_" + co.first;
		eff[o][col][PD][type]["Val"][anaBin] = effP.GetEfficiency(iBin);
		const bool dropUnc = ( (isNominal==false) || (type.find("TnP_")!=std::string::npos) || (type.find("MC_Syst")!=std::string::npos) );
		eff[o][col][PD][type]["Err_Stat_High"][anaBin] = (dropUnc ? 0.0 : effP.GetEfficiencyErrorUp(iBin));
		eff[o][col][PD][type]["Err_Stat_Low"][anaBin] = (dropUnc ? 0.0 : effP.GetEfficiencyErrorLow(iBin));
		eff[o][col][PD][type]["Err_Syst_High"][anaBin] = 0.0;
		eff[o][col][PD][type]["Err_Syst_Low"][anaBin] = 0.0;
		eff[o][col][PD][type]["Err_Tot_High"][anaBin] = 0.0;
		eff[o][col][PD][type]["Err_Tot_Low"][anaBin] = 0.0;
		for (const auto& e : eff.at(o).at(col).at(PD).at(type)) {
		  if (!contain(e.second, anaBin)) { anaBin.print(); std::cout << "EFF TYPE NOT FOUND FOR " << e.first << " " << PD << " " << col << " " << o << std::endl; }
		  if (contain(e.second, anaBin) && (e.second.at(anaBin)<0. || isnan(e.second.at(anaBin)))) {
		    std::cout << "[ERROR] Invalid " << type << " " << e.first << " ( " << e.second.at(anaBin) << " )  in pd: " << PD << " , obj: " << o << std::endl; anaBin.print(); return false;
		  }
		}
	      }
	    }
	    //
	    // Fill the total uncertainties
	    for (auto& e :eff[o][col][PD]) {
	      e.second["Err_Tot_Low" ][anaBin] = sumErrors({e.second["Err_Stat_Low" ][anaBin], e.second["Err_Syst_Low" ][anaBin]});
	      e.second["Err_Tot_High"][anaBin] = sumErrors({e.second["Err_Stat_High"][anaBin], e.second["Err_Syst_High"][anaBin]});
	    }
	  }
	}
      }
    }
  }
  EFF_[effType] = eff;
  return true;
};


double getCorrectedYieldValue(const double& rawN, const double& acc, const double& eff)
{
  // Check inputs
  if (rawN<-300.0 || isnan(rawN)) { throw std::logic_error(Form("[ERROR] Invalid raw yield value ( %.2f )", rawN)); }
  if (acc<=0.0 || isnan(acc)) { throw std::logic_error(Form("[ERROR] Invalid acceptance value ( %.2f )", acc)); }
  if (eff<=0.0 || isnan(eff)) { throw std::logic_error(Form("[ERROR] Invalid efficiency value ( %.2f )", eff)); }
  // Return result
  return divide({rawN}, {acc, eff});
};

  
double getCorrectedYieldError(const double& rawN, const double& acc, const double& eff, const double& errRawN, const double& errAcc, const double& errEff)
{
  // Check inputs
  if (errRawN<0.0 || isnan(errRawN)) { throw std::logic_error(Form("[ERROR] Invalid raw yield uncertainty ( %.2f )", errRawN)); }
  if (errAcc<0.0 || isnan(errAcc)) { throw std::logic_error(Form("[ERROR] Invalid acceptance uncertainty ( %.2f )", errAcc)); }
  if (errEff<0.0 || isnan(errEff)) { throw std::logic_error(Form("[ERROR] Invalid efficiency uncertainty ( %.2f )", errEff)); }
  // Return result
  return divideError({rawN}, {acc, eff}, {errRawN}, {errAcc, errEff});
};


bool correctRawYields(VarBinTriMap_t& inputVar, const BinSextaMap_t& effMap, const bool& isNominal, const std::string& accType, const std::string& effType)
{
  //
  const std::string accName = ( (accType!="") ? Form("Acceptance_%s", accType.c_str()) : "" );
  const std::string effName = ( (effType!="") ? Form("Efficiency_%s", effType.c_str()) : "" );
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
	    if (Acceptance.at("Val")==0.0) { std::cout << "[ERROR] Null acceptance for: " << obj << " at "; b.first.print(); return false; }
	    if (Efficiency.at("Val")==0.0) { std::cout << "[ERROR] Null efficiency for: " << obj << " at "; b.first.print(); return false; }
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
	    // Check raw yield, acceptance and efficiency
	    for (const auto& t : N_Raw) {
	      const auto& type = t.first;
	      const auto& n = t.second;
	      if (((type.rfind("Err_Tot",0)==0 || type.rfind("Err_Stat",0)==0) && n<=0.) || (type.rfind("Err_Syst",0)==0 && n<0.) || (type=="Val" && (pd.first!="DIMUON" ? n<-300. : n<=0.)) || isnan(n)) {
	        std::cout << "[ERROR] Invalid raw yield " << t.first << " ( " << t.second << " )  in pd: " << pd.first << " , obj: " << obj << std::endl; b.first.print(); return false;
	      }
	      const auto& acc = Acceptance.at(t.first);
	      if (((type.rfind("Err_Tot",0)==0 || type.rfind("Err_Stat",0)==0) && acc<=0.) || (type.rfind("Err_Syst",0)==0 && acc<0.) || (type=="Val" && acc<=0.) || isnan(acc)) {
	        std::cout << "[ERROR] Invalid acceptance " << t.first << " ( " << acc << " )  in pd: " << pd.first << " , obj: " << obj << std::endl; b.first.print(); return false;
	      }
	      const auto& eff = Efficiency.at(t.first);
	      if (((type.rfind("Err_Tot",0)==0 || type.rfind("Err_Stat",0)==0) && eff<=0.) || (type.rfind("Err_Syst",0)==0 && eff<0.) || (type=="Val" && eff<=0.) || isnan(eff)) {
	        std::cout << "[ERROR] Invalid efficiency " << t.first << " ( " << eff << " )  in pd: " << pd.first << " , obj: " << obj << std::endl; b.first.print(); return false;
	      }
	    }
	    // Fill with the corrected values
	    N_Corr.at("Val") = getCorrectedYieldValue(N_Raw.at("Val"), Acceptance.at("Val"), Efficiency.at("Val"));
	    N_Corr.at("Err_Stat_Low" ) = getCorrectedYieldError(N_Raw.at("Val"), Acceptance.at("Val"), Efficiency.at("Val"), N_Raw.at("Err_Stat_Low" ), 0.0 , 0.0);
	    N_Corr.at("Err_Stat_High") = getCorrectedYieldError(N_Raw.at("Val"), Acceptance.at("Val"), Efficiency.at("Val"), N_Raw.at("Err_Stat_High"), 0.0 , 0.0);
	    N_Corr.at("Err_Syst_Low" ) = getCorrectedYieldError(N_Raw.at("Val"), Acceptance.at("Val"), Efficiency.at("Val"), N_Raw.at("Err_Syst_Low" ), Acceptance.at("Err_Tot_Low" ), Efficiency.at("Err_Tot_Low" ));
	    N_Corr.at("Err_Syst_High") = getCorrectedYieldError(N_Raw.at("Val"), Acceptance.at("Val"), Efficiency.at("Val"), N_Raw.at("Err_Syst_High"), Acceptance.at("Err_Tot_High"), Efficiency.at("Err_Tot_High"));
	    //
	    if (isNominal) {
	      for (const auto& effT : eff) {
		if (
		    (accType == "MC"  && effT.first=="Acceptance_MC") ||
		    (effType == "MC"  && effT.first=="Efficiency_MC") ||
		    (effType == "TnP" && effT.first.find("Efficiency")!=std::string::npos)
		    ) {
		  auto& N_Corr = inputVar.at(o.first).at(c.first).at(pd.first).at(b.first)[vN+"_"+effT.first];
		  const auto& effVal = (effT.first=="Acceptance_MC" ? 1.0 : effT.second.at("Val").at(b.first));
		  N_Corr["Val"] = getCorrectedYieldValue(N_Raw.at("Val") , Acceptance.at("Val") , effVal);
		  for (const auto& errLbL : StringVector_t({"Err_Stat_Low", "Err_Stat_High", "Err_Syst_Low", "Err_Syst_High"})) {
		    const auto& nErr = (errLbL.find("_Low")!=std::string::npos ? N_Raw.at("Err_Syst_Low") : N_Raw.at("Err_Syst_High"));
		    const auto& effErr = (effT.first=="Acceptance_MC" ? 0.0 : effT.second.at(errLbL).at(b.first));
		    N_Corr[errLbL] = getCorrectedYieldError(N_Raw.at("Val"), Acceptance.at("Val"), effVal, nErr, Acceptance.at(errLbL), effErr);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  return true;
};


double getCrossSectionValue(const double& N, const double& lumi, const double& binWidth)
{
  // Check inputs
  if (N<-1700. || isnan(N)) { throw std::logic_error(Form("[ERROR] Invalid yield value ( %.2f )", N)); }
  if (lumi<=0. || isnan(lumi)) { throw std::logic_error(Form("[ERROR] Invalid luminosity value ( %.2f )", lumi)); }
  if (binWidth<=0.0 || isnan(binWidth)) { throw std::logic_error(Form("[ERROR] Invalid bin width ( %.2f )", binWidth)); }
  // Return result
  return divide({N}, {lumi, binWidth});
};

  
double getCrossSectionError(const double& N, const double& lumi, const double& binWidth, const double& errN, const double& errLumi)
{
  // Check inputs
  if (errN<=0.0 || isnan(errN)) { throw std::logic_error(Form("[ERROR] Invalid yield uncertainty ( %.2f )", errN)); }
  if (errLumi<0.0 || isnan(errLumi)) { throw std::logic_error(Form("[ERROR] Invalid luminosity uncertainty ( %.2f )", errLumi)); }
  // Return result
  return divideError({N}, {lumi, binWidth}, {errN}, {errLumi, 0.0});
};


bool computeCrossSection(BinSextaMap_t& var, const VarBinTriMap_t& inputVar, const bool& doSyst, const VarBinTriMap_t& nomVar)
{
  // Loop over data
  for (const auto& o : inputVar) {
    for (const auto& c : o.second) {
      for (const auto& pd : c.second) {
	for (const auto& b : pd.second) {
	  for (const auto& v : b.second) {
	    if (v.first.rfind("N_",0)!=0 || v.first.find("To")==std::string::npos) continue;
	    if (v.first.rfind("N_Bkg",0)==0 || v.first.find("SS")!=std::string::npos) continue;
	    const auto& obj = v.first.substr(0, v.first.find("To")).substr(2);
	    // Extract yield and luminosity
	    const auto& iVar = v.second;
	    const auto& lumi = b.second.at("Luminosity").at("Val");
	    // Check yield and luminosity
	    for (const auto& t : iVar) {
	      if ((t.first.rfind("Err_",0)==0 && t.second<=0.) || (pd.first!="DIMUON" ? t.second<-1700. : t.second<=0.) || isnan(t.second)) {
	        std::cout << "[ERROR] Invalid yield " << t.first << " ( " << t.second << " )  in pd: " << pd.first << " , obj: " << obj << std::endl; b.first.print(); return false;
	      }
	    }
	    if (lumi<=0. || isnan(lumi)) {
	      std::cout << "[ERROR] Invalid luminosity value ( " << lumi << " )  in pd: " << pd.first << " , obj: " << obj << std::endl; b.first.print(); return false;
	    }
	    // Compute the cross-section value and errors
	    auto& oVar = var[obj][c.first][pd.first]["Cross_Section"];
	    oVar["Val"][b.first] = getCrossSectionValue(iVar.at("Val"), lumi, 1.0);
	    // Compute the statistical error
	    for (const auto& t : STAT) { oVar[t][b.first] = getCrossSectionError(iVar.at("Val"), lumi, 1.0, iVar.at(t), 0.0); }
	    // Compute the systematic error
	    for (const auto& t : SYST) { oVar[t][b.first] = 0.0; }
	    if (!doSyst) continue;
	    else if (nomVar.size()==0) {
	      for (const auto& t : SYST) { oVar[t][b.first] = getCrossSectionError(iVar.at("Val"), lumi, 1.0, iVar.at(t), 0.0); }
	    }
	    else if (nomVar.size()>0) {
	      const auto& nVar = nomVar.at(o.first).at(c.first).at(pd.first).at(b.first).at(v.first);
	      const auto nomVal = getCrossSectionValue(nVar.at("Val"), lumi, 1.0);
	      const auto varVal = getCrossSectionValue(iVar.at("Val"), lumi, 1.0);
	      for (const auto& t : SYST) { oVar[t][b.first] = std::abs(varVal - nomVal); }
	      oVar.at("Val").at(b.first) = varVal;
	    }
	  }
	}
      }
    }
  }
  //
  return true;
};


bool computeRatioTo1S(BinSextaMap_t& var, const VarBinTriMap_t& inputVar, const bool& doSyst, const VarBinTriMap_t& nomVar, const IntMap_t& systCorr)
{
  // Check inputs
  if (doSyst) {
    if (!contain(systCorr, "Obj")) { throw std::runtime_error("computeRatioTo1S: systCorr has no Obj key"); }
    if (std::abs(systCorr.at("Obj"))!=1 && std::abs(systCorr.at("Obj"))!=2) { throw std::runtime_error(Form("computeRatioTo1S: systCorr has wrong value (%d)", systCorr.at("Obj"))); }
  }
  // Loop over data
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
	    if (v.rfind("R_",0)==0) {
	      for (const auto& p : iVar) { oVar[p.first][b.first] = p.second; }
	    }
	    else if (v.rfind("N_",0)==0) {
	      auto r = v;
	      if (obj.rfind("Psi2S",0)==0) { stringReplace(r, "Psi2S", "JPsi"); }
	      else if (obj.rfind("Ups2S",0)==0) { stringReplace(r, "Ups2S", "Ups1S"); }
	      else if (obj.rfind("Ups3S",0)==0) { stringReplace(r, "Ups3S", "Ups1S"); }
	      if (r==v) continue;
	      const auto& rVar = b.second.at(r);
	      // Check nominal yields
	      for (const auto& t : iVar) {
	        if ((t.first.rfind("Err_",0)==0 && t.second<=0.) || (pd.first!="DIMUON" ? t.second<-1700. : t.second<=0.) || isnan(t.second)) {
	          std::cout << "[ERROR] Invalid " << v << " " << t.first << " ( " << t.second << " )  in pd: " << pd.first << " , obj: " << obj << std::endl; b.first.print(); return false;
	        }
	      }
	      for (const auto& t : rVar) {
	        if (t.second<=0. || isnan(t.second)) {
	          std::cout << "[ERROR] Invalid " << r << " " << t.first << " ( " << t.second << " )  in pd: " << pd.first << " , obj: " << obj << std::endl; b.first.print(); return false;
	        }
	      }
	      // Compute the ratio to 1S
	      oVar["Val"][b.first] = divide(iVar.at("Val"), rVar.at("Val"));
	      // Compute the statistical error
	      for (const auto& t : STAT) { oVar[t][b.first] = divideError(iVar.at("Val"), rVar.at("Val"), iVar.at(t), rVar.at(t), r==v); }
	      // Compute the systematic error
	      for (const auto& t : SYST) { oVar[t][b.first] = 0.0; }
	      if (!doSyst) continue;
	      else if (nomVar.size()==0) {
		const bool isCorrelated = (r==v || std::abs(systCorr.at("Obj"))==2);
		for (const auto& t : SYST) { oVar[t][b.first] = divideError(iVar.at("Val"), rVar.at("Val"), iVar.at(t), rVar.at(t), isCorrelated); }
	      }
	      else if (nomVar.size()>0) {
		const auto& iNomVar = nomVar.at(o.first).at(c.first).at(pd.first).at(b.first).at(v);
		const auto& rNomVar = nomVar.at(o.first).at(c.first).at(pd.first).at(b.first).at(r);
		const auto nomVal = divide(iNomVar.at("Val"), rNomVar.at("Val"));
		oVar.at("Val").at(b.first) = nomVal;
		// Variation method: object correlated
		if (r==v || systCorr.at("Obj")==2) {
		  const auto varVal = divide(iVar.at("Val"), rVar.at("Val"));
		  for (const auto& t : SYST) { oVar[t][b.first] = std::abs(varVal - nomVal); }
		  oVar.at("Val").at(b.first) = varVal;
		}
		// Variation method: object uncorrelated
		else if (systCorr.at("Obj")==1) {
		  const auto vVar_Obj = divide(iVar.at("Val"), rNomVar.at("Val"));
		  const auto vVar_Ref = divide(iNomVar.at("Val"), rVar.at("Val"));
		  for (const auto& t : SYST) {
		    var.at(obj).at(c.first).at(pd.first)["RatioTo1S_Obj"][t][b.first] = std::abs(vVar_Obj - nomVal);
		    var.at(obj).at(c.first).at(pd.first)["RatioTo1S_Ref"][t][b.first] = std::abs(vVar_Ref - nomVal);
		  }
		}
		// Propagation method
		else if (systCorr.at("Obj")<0) {
		  const auto iVarErr = std::abs(iVar.at("Val") - iNomVar.at("Val"));
		  const auto rVarErr = std::abs(rVar.at("Val") - rNomVar.at("Val"));
		  const bool isCorrelated = (r==v || systCorr.at("Obj")==-2);
		  for (const auto& t : SYST) { oVar[t][b.first] = divideError(iNomVar.at("Val"), rNomVar.at("Val"), iVarErr, rVarErr, isCorrelated); }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  //
  return true;
};


bool computeForwardBackwardRatio(BinSextaMap_t& var, const VarBinTriMap_t& inputVar, const bool& doSyst, const VarBinTriMap_t& nomVar, const IntMap_t& systCorr)
{
  // check inputs
  if (doSyst) {
    if (!contain(systCorr, "Rap")) { throw std::runtime_error("computeForwardBackwardRatio: systCorr has no Obj key"); }
    if (std::abs(systCorr.at("Rap"))!=1 && std::abs(systCorr.at("Rap"))!=2) { throw std::runtime_error(Form("computeForwardBackwardRatio: systCorr has wrong value (%d)", systCorr.at("Rap"))); }
  }
  // loop over data
  for (const auto& o : inputVar) {
    for (const auto& c : o.second) {
      for (const auto& pd : c.second) {
	for (const auto& b : pd.second) {
	  const auto& anaBin = b.first;
	  const std::string rapLbl = (anaBin.hasbin("Cand_Rap") ? "Cand_Rap" : (anaBin.hasbin("Cand_RapCM") ? "Cand_RapCM" : ""));
	  if (rapLbl=="") continue; // ignore those that do not have rapidity bins
	  if (anaBin.getbin(rapLbl).low() < 0.0) continue; // Ignore backward bins when looping
	  // Build the backward and forward bin
	  auto binFw = anaBin; binFw.setbin(rapLbl,    anaBin.getbin(rapLbl).low() ,    anaBin.getbin(rapLbl).high());
	  auto binBw = anaBin; binBw.setbin(rapLbl, -1*anaBin.getbin(rapLbl).high(), -1*anaBin.getbin(rapLbl).low() );
	  // Check that everything is fine
	  if (pd.second.count(binBw)==0) continue; // Ignore if backward bin is not present
	  if (pd.second.count(binFw)==0) { std::cout << "[ERROR] Not found forward bin : ";  binFw.print(); return false; }
	  //
	  for (const auto& v : b.second) {
	    if (v.first.rfind("N_",0)!=0 || v.first.find("To")==std::string::npos) continue;
	    if (v.first.rfind("N_Bkg",0)==0 || v.first.find("SS")!=std::string::npos) continue;
	    const auto& obj = v.first.substr(0, v.first.find("To")).substr(2);
	    // Initialize the forward and backward yields
	    const auto& iVar_Fw = pd.second.at(binFw).at(v.first);
	    const auto& iVar_Bw = pd.second.at(binBw).at(v.first);
	    // Check forward and backward yields
	    for (const auto& t : iVar_Fw) {
	      if (t.second<=0. || isnan(t.second)) {
	        std::cout << "[ERROR] Invalid forward yield " << t.first << " ( " << t.second << " )  in pd: " << pd.first << " , obj: " << obj << std::endl; b.first.print(); return false;
	      }
	    }
	    for (const auto& t : iVar_Fw) {
	      if (t.second<=0. || isnan(t.second)) {
	        std::cout << "[ERROR] Invalid backward yield " << t.first << " ( " << t.second << " )  in pd: " << pd.first << " , obj: " << obj << std::endl; b.first.print(); return false;
	      }
	    }
	    // Compute the forward and backward ratio
	    auto& oVar = var[obj][c.first][pd.first]["ForwardBackward_Ratio"];
	    oVar["Val"][b.first] = divide(iVar_Fw.at("Val"), iVar_Bw.at("Val"));
	    // Compute the statistical error
	    for (const auto& t : STAT) { oVar[t][b.first] = divideError(iVar_Fw.at("Val"), iVar_Bw.at("Val"), iVar_Fw.at(t), iVar_Bw.at(t)); }
	    // Compute the systematic error
	    for (const auto& t : SYST) { oVar[t][b.first] = 0.0; }
	    if (!doSyst) continue;
	    else if (nomVar.size()==0) {
	      const bool isCorrelated = std::abs(systCorr.at("Rap"))==2;
	      for (const auto& t : SYST) { oVar[t][b.first] = divideError(iVar_Fw.at("Val"), iVar_Bw.at("Val"), iVar_Fw.at(t), iVar_Bw.at(t), isCorrelated); }
	    }
	    else if (nomVar.size()>0) {
	      //
	      const auto& nVar_Fw = nomVar.at(o.first).at(c.first).at(pd.first).at(binFw).at(v.first);
	      const auto& nVar_Bw = nomVar.at(o.first).at(c.first).at(pd.first).at(binBw).at(v.first);
	      const auto nomVal = divide(nVar_Fw.at("Val"), nVar_Bw.at("Val"));
	      oVar.at("Val").at(b.first) = nomVal;
	      // Variation method: object correlated
	      if (systCorr.at("Rap")==2) {
		const auto varVal = divide(iVar_Fw.at("Val"), iVar_Bw.at("Val"));
		for (const auto& t : SYST) { oVar[t][b.first] = std::abs(varVal - nomVal); }
		oVar.at("Val").at(b.first) = varVal;
	      }
	      // Variation method: object uncorrelated
	      else if (systCorr.at("Rap")==1) {
		const auto vVar_Fw = divide(iVar_Fw.at("Val"), nVar_Bw.at("Val"));
		const auto vVar_Bw = divide(nVar_Fw.at("Val"), iVar_Bw.at("Val"));
		for (const auto& t : SYST) {
		  var.at(obj).at(c.first).at(pd.first)["ForwardBackward_Ratio_Fw"][t][b.first] = std::abs(vVar_Fw - nomVal);
		  var.at(obj).at(c.first).at(pd.first)["ForwardBackward_Ratio_Bw"][t][b.first] = std::abs(vVar_Bw - nomVal);
		}
	      }
	      // Propagation method
	      else if (systCorr.at("Rap")<0) {
		const auto iVarErr_Fw = std::abs(iVar_Fw.at("Val") - nVar_Fw.at("Val"));
		const auto iVarErr_Bw = std::abs(iVar_Bw.at("Val") - nVar_Bw.at("Val"));
		const bool isCorrelated = (systCorr.at("Obj")==-2);
		for (const auto& t : SYST) { oVar[t][b.first] = divideError(iVar_Fw.at("Val"), iVar_Bw.at("Val"), iVarErr_Fw, iVarErr_Bw, isCorrelated); }
	      }
	    }
	  }
	}
      }
    }
  }
  //
  return true;
};


void computeSystematic(BinSeptaMap_t& varMap, BinOctaMapVec_t& systVarMap)
{
  // Extract nominal result
  auto& nomVar = varMap.at("Nominal");
  // Loop over systematic categories
  for (const auto& cat : systVarMap) {
    auto& var = varMap[cat.first];
    auto& systVarVec = systVarMap.at(cat.first);
    const auto origVarVec = systVarVec;
    var.clear(); systVarVec.clear();
    // Loop over systematic sub-categories
    for (const auto& lbl : origVarVec) {
      systVarVec[lbl.first].clear();
      const auto& nVariation = lbl.second.size();
      std::cout << "[INFO] Processing " << nVariation << " systematic variations for: " << cat.first << " " << lbl.first << std::endl;
      bool isEmpty = true;
      BinSextaMap_t systVarTot;
      // Loop over systematic variations
      for (size_t iVr=0; iVr<lbl.second.size(); iVr++) {
	BinSextaMap_t systVar;
	for (const auto& o : lbl.second[iVr]) {
	  for (const auto& c : o.second) {
	    for (const auto& pd : c.second) {
	      for (const auto& v : pd.second) {
		auto vLbl = v.first;
		if (std::count(v.first.begin(), v.first.end(), '_')>1) { vLbl = vLbl.substr(0, vLbl.rfind("_")); }
		auto& valMap = systVar[o.first][c.first][pd.first][vLbl];
		// Add systematic uncertainties
		for (const auto& t : StringVector_t({"Err_Syst_Low", "Err_Syst_High"})) {
		  for (const auto& b : v.second.at(t)) {
		    // Combine the varied uncertainties
		    valMap.at(t).at(b.first) = (b.second!=0. ? sumErrors({valMap[t][b.first], b.second}) : valMap[t][b.first]);
		    if (b.second!=0.) { isEmpty = false; }
		    // Compute total systematic uncertainty
		    if (iVr==0) {
		      double uncVal = 0.0;
		      if (nVariation == 1) {
			uncVal = std::abs(lbl.second[0].at(o.first).at(c.first).at(pd.first).at(v.first).at(t).at(b.first));
		      }
		      else if (nVariation == 2) {
			const auto diff_0 = std::abs(lbl.second[0].at(o.first).at(c.first).at(pd.first).at(v.first).at(t).at(b.first));
			const auto diff_1 = std::abs(lbl.second[1].at(o.first).at(c.first).at(pd.first).at(v.first).at(t).at(b.first));
			uncVal = symError(diff_0, diff_1);
		      }
		      else if (nVariation > 2) {
			double sum = 0.0;
			for (size_t i = 0; i < nVariation; i++) {
			  sum += std::pow(lbl.second[i].at(o.first).at(c.first).at(pd.first).at(v.first).at(t).at(b.first), 2.0);
			}
			uncVal = std::sqrt(sum / nVariation);
		      }
		      auto& val = systVarTot[o.first][c.first][pd.first][vLbl][t][b.first];
		      if (vLbl==v.first) { val = uncVal; }
		      else if (uncVal!=0.) { val = sumErrors({val, uncVal}); }
		    }
		  }
		}
	      }
	      for (const auto& v : pd.second) {
		if (std::count(v.first.begin(), v.first.end(), '_')>1) continue;
		auto& valMap = systVar[o.first][c.first][pd.first][v.first];
		// Add value and statistical uncertainty
		valMap["Nom"] = nomVar.at(o.first).at(c.first).at(pd.first).at(v.first).at("Val");
		for (const auto& t : StringVector_t({"Val", "Err_Stat_Low", "Err_Stat_High"})) {
		  valMap[t] = pd.second.at(v.first).at(t);
		}
		if (iVr==0) {
		  for (const auto& t : StringVector_t({"Nom", "Val", "Err_Stat_Low", "Err_Stat_High"})) {
		    systVarTot[o.first][c.first][pd.first][v.first][t] = valMap.at(t);
		  }
		}
		for (const auto& b : valMap.at("Val")) {
		  if (valMap.at("Val").at(b.first)==valMap.at("Nom").at(b.first)) {
		    valMap.at("Nom").at(b.first) += symError(valMap.at("Err_Syst_Low").at(b.first), valMap.at("Err_Syst_High").at(b.first));
		    if (iVr==0) {
		      auto& valMapTot = systVarTot.at(o.first).at(c.first).at(pd.first).at(v.first);
		      valMapTot.at("Val").at(b.first) += symError(valMapTot.at("Err_Syst_Low").at(b.first), valMapTot.at("Err_Syst_High").at(b.first));
		    }
		  }
		}
	      }
	    }
	  }
	}
        if (isEmpty) { std::cout << "[WARNING] Variation " << lbl.first << " index " << iVr << " is empty. Ignoring it!" << std::endl; break; }
	systVarVec.at(lbl.first).push_back(systVar);
      }
      if (isEmpty) { systVarVec.erase(lbl.first); continue; }
      // Add the total systematic uncertainty
      for (const auto& o : systVarTot) {
	for (const auto& c : o.second) {
	  for (const auto& pd : c.second) {
	    for (const auto& v : pd.second) {
	      for (const auto& t : v.second) {
		for (const auto& b : t.second) {
		  auto& val = var[o.first][c.first][pd.first][v.first][t.first][b.first];
		  auto& nomVal = nomVar.at(o.first).at(c.first).at(pd.first).at(v.first).at((t.first=="Nom")?"Val":t.first).at(b.first);
		  if (t.first.rfind("Err_Syst_",0)==0) {
		    val = (b.second!=0. ? sumErrors({val, b.second}) : val);
		    if (cat.first!="Efficiency") {// Avoid double counting of TnP eff unc
		      nomVal = (b.second!=0. ? sumErrors({nomVal, b.second}) : nomVal);
		    }
		  }
		  else if (t.first=="Nom") { val = nomVal; }
		  else if (t.first=="Val" && origVarVec.size()>1) {
		    val = nomVal + symError(var.at(o.first).at(c.first).at(pd.first).at(v.first).at("Err_Syst_Low").at(b.first),
					    var.at(o.first).at(c.first).at(pd.first).at(v.first).at("Err_Syst_High").at(b.first));
		  }
		  else if (t.first=="Val") { val = b.second; }
		}
	      }
	    }
	  }
	}
      }
      systVarVec.at(lbl.first).push_back(systVarTot);
    }
  }
};


#endif // #ifndef extractResultsTree_C
