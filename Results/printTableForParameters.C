#ifndef printTableForParameters_C
#define printTableForParameters_C

// Auxiliary Headers
#include "Utilities/resultsUtils.h"
#include "extractResultsTree.C"
#include "processResultsTree.C"
#include <fstream>


void printTableForParameters(
				   const std::string& workDirName = "Test",
				   const StringVector_t& objTags  = { "JPsi" },    //"Bkg", "JPsi", "Psi2S", "Ups1S", "Ups2S", "Ups3S", "Z", "D0"
				   const StringVector_t& trgTags  = { "CatPR_DIMUON",  "CatPR_MINBIAS"},
				   const std::string& dataTag = "MC"
				   )
{
  //
  // Check input
  if (workDirName=="") { std::cout << "[ERROR] Name of working directory was not set!" << std::endl; return; }
  if (trgTags.empty()) { std::cout << "[ERROR] PD tags were not set!" << std::endl; return; }
  if (objTags.empty()) { std::cout << "[ERROR] Object tags were not set!" << std::endl; return; }
  const std::string& varTag  = "Cand_Mass";
  const StringVector_t colTags = {"PA8Y16"};
  //
  const std::string& CWD = getcwd(NULL, 0);
  const std::string& outDir = CWD + "/Tables/Parameters/" + workDirName;
  makeDir(outDir.c_str());
  //
  // Get the result
  VarBinTriMap_t inputVar;
  std::cout << "[INFO] Adding results for: " << workDirName << std::endl;
  for (const auto& trgTag : trgTags) {
    for (const auto& colTag : colTags) {
      for (const auto& objTag : objTags) {
	if (!extractResultsTree(inputVar, workDirName, trgTag, colTag, objTag, dataTag, varTag, false)) { return; }
      }
    }
  }
  //
  std::vector<std::string> latexTable;
  //
  std::map<std::string, std::string> ParMap = {{"f", "f"}, {"Alpha", "\\alpha"}, {"n", "n"},
						 {"AlphaR", "\\AlphaR"}, {"nR", "\\nR"}, 
						 {"rSigma21", "\\rSigma"}, {"Sigma1", "\\sigma"}};
  //
  for (const auto& o : inputVar) {
    for (const auto& c : o.second) {
      //
      const auto& obj = o.first;
      const auto& col = c.first;
      //
      latexTable.push_back("\\begin{tabular}{|ccc|c");
      latexTable.push_back("  \\hline");
      //
      std::map<std::string, std::string> parMap;
      for (const auto& p : c.second.begin()->second.begin()->second) {
	auto par = p.first;
	par = par.substr(0,par.find("_"));
	if (contain(ParMap, par)) { parMap[p.first] = par; }
      }
      const uint& nPars = parMap.size();
      const auto& model = c.second.begin()->second.begin()->second.at("Model_Chi2").begin()->first;
      //
      for (uint i=1; i<nPars; i++) { latexTable.back() += "|c"; }; latexTable.back() += "}";
      //
      latexTable.push_back(Form("  \\multicolumn{3}{c}{Bin} & \\multicolumn{%d}{c}{Parameters}\\\\", nPars));
      latexTable.push_back("  Rapidity & \\pt & \\NTrk");
      for (const auto& p : parMap) { latexTable.back() += " & "+ParMap.at(p.second); }; latexTable.back() += "\\\\";
      latexTable.push_back("  \\hline\\hline");
      //
      for (const auto& pd : c.second) {
	for (const auto& b : pd.second) {
	  const auto& anabin = b.first;
	  const auto& bin_NTrack = anabin.getbin("NTrack");
	  const auto& bin_Rap = anabin.getbin("Cand_Rap");
	  const auto& bin_AbsRap = anabin.getbin("Cand_AbsRap");
	  const auto& bin_RapCM = anabin.getbin("Cand_RapCM");
	  const auto& bin_Pt = anabin.getbin("Cand_Pt");
	  //
	  std::string lbl_NTrack = Form("%g - %g", bin_NTrack.low(), bin_NTrack.high());
	  std::string lbl_Pt = Form("%g - %g", bin_Pt.low(), bin_Pt.high());
	  std::string lbl_Rap = "";
	  if (anabin.hasbin("Cand_Rap")) { lbl_Rap = Form("%g \\geq \\rap < %g", bin_Rap.low(), bin_Rap.high()); }
	  else if (anabin.hasbin("Cand_AbsRap")) { lbl_Rap = Form("%g \\geq \\absRap < %g", bin_AbsRap.low(), bin_AbsRap.high()); }
	  else if (anabin.hasbin("Cand_RapCM")) { lbl_Rap = Form("%g \\geq \\rapCM} < %g", bin_RapCM.low(), bin_RapCM.high()); }
	  latexTable.push_back(Form("    %s & %s & %s", lbl_Rap.c_str(), lbl_Pt.c_str(), lbl_NTrack.c_str()));
	  //
	  // Extract the parameter values
	  std::vector<double> parValV;
	  for (const auto& p : parMap) {
	    if (contain(b.second, p.first)) {
	      auto val = b.second.at(p.first).at("Val");
	      if (p.first.rfind("Sigma",0)==0) { val *= 1000.; }
	      parValV.push_back(val);
	    }
	  }
	  //
	  for (const auto& val : parValV) {
	    latexTable.back() += Form(" & %.2f", val);
	  }
	  latexTable.back() += "\\\\";
	}
      }
      //
      latexTable.push_back("  \\hline");
      latexTable.push_back("\\end{tabular}");
      //
      const std::string& outFileName = Form("%s/Table_SignalParameterStudy_%s_%s_%s_%s.tex", outDir.c_str(), dataTag.c_str(), model.c_str(), obj.c_str(), col.c_str());
      ofstream outFile(outFileName.c_str());
      if (outFile.is_open()) {
	for (const auto& line : latexTable) {
	  outFile << line << std::endl;
	}
	outFile.close();
      }
      else { std::cout << "[ERROR] File: " << outFileName << " failed to open!" << std::endl; }
    }
  }
  //
  for (const auto& line : latexTable) {
    std::cout << line << std::endl;
  }
};


#endif // ifndef printTableForParameters_C
