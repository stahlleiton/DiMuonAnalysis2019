#if !defined(__CINT__) || defined(__MAKECINT__)
// Auxiliary Headers
#include "../Utilities/Ntuple/VertexCompositeTree.h"
#include "../Utilities/RunInfo/eventUtils.h"
#include "util.h"
//
#include "TEfficiency.h"
#include "TH1D.h"
#include "tnp_weight_lowPt.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include <iostream>

#endif


void decayCutEfficiency(const std::string& type="PromptJPsi", const bool& muonTrig=true, const double& thr=0.95, const bool& corrEff=true)
{
  //
  std::cout << "[INFO] Computing efficiency for " << type << " including muon trigger " << muonTrig <<  " with decay cut threshold " << thr << " and TnP corrections " << corrEff << std::endl; 
  //
  const std::string& fType = (type=="PromptJPsi" ? "JPsiToMuMu" : (type=="PromptPsi2S" ? "Psi2SToMuMu" : ((type=="NonPromptJPsi" || type=="NonPromptPsi2S") ? "BToPsiToMuMu" : "")));
  const std::map<std::string, std::string> inputFiles = {
    { "pPb" , "/Users/andre/DiMuonAnalysis2019/Tree/VertexCompositeTree_"+fType+"_pPb-Bst_pPb816Summer16_DiMuMC.root" },
    { "Pbp" , "/Users/andre/DiMuonAnalysis2019/Tree/VertexCompositeTree_"+fType+"_PbP-Bst_pPb816Summer16_DiMuMC.root" }
  };
  const auto& treeDir = "dimucontana_mc"; // For MC use dimucontana_mc

  // Initialize the efficiencies
  std::map<std::string, long int> nevents;
  std::map<std::string, std::map< anabin_t, std::vector< std::unique_ptr<TEfficiency> > > > eff;

  for(const auto& t : inputFiles) {
    const auto& s = t.first;
    VertexCompositeTree tree;
    if (!tree.GetTree(t.second, treeDir)) { std::cout << "Invalid tree in: " << t.second << "!" << std::endl; return; }
    
    // Initialize the efficiency histograms
    std::map< anabin_t, std::vector< std::pair<std::unique_ptr<TH1D>, std::unique_ptr<TH1D>> > > hist;
    for(const auto& v1 : BINMAP) {
      const auto& var1 = v1.first;
      for(const auto& v2 : v1.second) {
	const auto& var2 = v2.first;
	const auto binPair = anabin_t({var1, var2});
	uint iHist = 0;
	for(const auto& v3 : v2.second) {
	  //
	  const auto& hBins = v3.second.data();
	  const auto& nBins = v3.second.size()-1;
	  const std::string& name = Form("Thr%.0f_%s_%.0f_%.0f_%s_%.0f_%.0f_%s_%d",
					 thr*100.,
					 var1.name().c_str(), var1.low()*10., var1.high()*10.,
					 var2.name().c_str(), var2.low()*10., var2.high()*10.,
					 v3.first.c_str(), iHist);
	  iHist++;
	  //
	  hist[binPair].push_back({std::unique_ptr<TH1D>(), std::unique_ptr<TH1D>()});
	  hist[binPair].back().first.reset(new TH1D(Form("hPassed_%s",name.c_str()), Form("hPassed_%s",name.c_str()), nBins, hBins));
	  hist[binPair].back().second.reset(new TH1D(Form("hTotal_%s",name.c_str()), Form("hTotal_%s",name.c_str()), nBins, hBins));
	  hist[binPair].back().first->Sumw2();
	  hist[binPair].back().second->Sumw2();
	}
      }
    }

    // Loop over the tree
    nevents[s] = tree.GetEntries();
    for(Long64_t jentry=0; jentry<nevents[s]; jentry++) { 
      if(!(jentry % 1000000)) std::cout<<"Processed "<<jentry<<" events out of "<<nevents[s]<<std::endl;
      if (tree.GetEntry(jentry)<0) { std::cout << "Invalid entry in: " << t.second << "!" << std::endl; return; }
      
      // Loop over the generated candidates
      for(uint iGen=0; iGen<tree.candSize(); iGen++){
	// Check the type of particle
        const auto& pid = (type.rfind("JPsi",0)==0 ? 443 : (type.rfind("Psi2S",0)==0 ? 100443 : -1));
        if (fabs(tree.PID_gen()[iGen])!=id) continue;

	// Check that both reconstructed muons are inside the single muon acceptance kinematic region
        const auto& pTD1 = tree.pTD1()[iReco];
        const auto& etaD1 = (s=="Pbp" ? -1.0 : 1.0) * tree.EtaD1()[iReco];
        const bool mu1InAccep = (muonTrig ? triggerMuonAcceptance(pTD1, etaD1) : muonAcceptance(pTD1, etaD1));
        const auto& pTD2 = tree.pTD2()[iReco];
        const auto& etaD2 = (s=="Pbp" ? -1.0 : 1.0) * tree.EtaD2()[iReco];
        const bool mu2InAccep = (muonTrig ? triggerMuonAcceptance(pTD2, etaD2) : muonAcceptance(pTD2, etaD2));
        if(!mu1InAccep || !mu2InAccep) continue;

	// Check that candidate is within acceptance
	const auto& pT = tree.pT()[iReco];
	

	// Check the passing condition
        const bool softCand = tree.softCand(iReco);
        const bool goodEvt = tree.evtSel()[0];
        const bool trigCand = (muonTrig ? (tree.trigHLT()[0] && tree.trigCand(0, iReco)) : true);
        const bool passed = (goodEvt && softCand && trigCand);
	if(!passed) continue;

	// Derive the efficiency correction
	double effCorr = 1.0;
	if(corrEff)
  	{
	  const auto& trgCorr = (muonTrig ? (tnp_weight_trg_ppb(pTD1, etaD1, 0) * tnp_weight_trg_ppb(pTD2, etaD2, 0)) : 1.0);
	  const auto& trkCorr = 1.0;//tnp_weight_trk_ppb(pTD1, etaD1, 0) * tnp_weight_trk_ppb(pTD2, etaD2, 0);
	  const auto& muidCorr = tnp_weight_muid_ppb(pTD1, etaD1, -10) * tnp_weight_muid_ppb(pTD2, etaD2, -10);
	  effCorr = trgCorr*trkCorr*muidCorr;
	}

	for(const auto& b : bins)
	{
	  const auto& yMin = b.first.first;
	  const auto& yMax = b.first.second;
	  const auto& bin = std::make_pair(yMin, yMax);
	  // Check if pass decay length cut
          const auto& y = (s=="Pbp" ? -1.0 : 1.0) * tree.y()[iReco];
	  const auto& pT = tree.pT[iReco];
	  const auto& eta = (s.first=="Pbp" ? -1.0 : 1.0)*tree.eta()[iReco];
          const auto& p = pT*std::cosh(eta);
          const auto& decayLen = (tree.V3DDecayLength()[iReco]*tree.V3DCosPointingAngle()[iReco])*(3.0969/p)*10.0;
	  const auto& passDL = (decayLen < decayLenCut(pT, y, muonTrig, thr));
	  // Fill histograms
	  hTotal.at(bin)->Fill(pT, effCorr);
	  if(passDL) { hPasseddat(bin)->Fill(pT, effCorr); }
 	}
      }
    }
    // Set the histograms to the efficiency
    for(const auto& b : bins)
    {
      const auto& yMin = b.first.first;
      const auto& yMax = b.first.second;
      const auto& bin = std::make_pair(yMin, yMax);
      const std::string& name = Form("effDL_Thr%.0f_Rap_%.0f_%.0f_%s", thr*100., yMin*10., yMax*10., s.c_str());
      eff[s][bin].reset(new TEfficiency(name.c_str(), name.c_str(), size, pTBin));
      eff.at(s).at(bin)->SetTotalHistogram(*hTotal.at(bin), "f");
      eff.at(s).at(bin)->SetPassedHistogram(*hPassed.at(bin), "f");
      if(TEfficiency::CheckWeights(*hPassed.at(bin), *hTotal.at(bin))) { eff.at(s).at(bin)->SetStatisticOption(TEfficiency::kBJeffrey); }
    }
  }
  
  // Combine the pPb and Pbp efficiencies
  for(const auto& b : bins)
  {
    const auto& yMin = b.first.first;
    const auto& yMax = b.first.second;
    const auto& bin = std::make_pair(yMin, yMax);
    const std::string& name = Form("effDL_Thr%.0f_Rap_%.0f_%.0f_pA", thr*100., yMin*10., yMax*10.);
    eff["pA"][bin].reset(new TEfficiency(name.c_str(), name.c_str(), size, pTBin));

    auto hPassed = *((TH1D*)eff.at("pPb").at(bin)->GetPassedHistogram());
    hPassed.Add(eff.at("Pbp").at(bin)->GetPassedHistogram(), eff.at("Pbp").at(bin)->GetPassedHistogram(), 110.8/nevents.at("pPb"), 62.6/nevents.at("Pbp"));
    auto hTotal = *((TH1D*)eff.at("pPb").at(bin)->GetTotalHistogram());
    hTotal.Add(eff.at("pPb").at(bin)->GetTotalHistogram(), eff.at("Pbp").at(bin)->GetTotalHistogram(), 110.8/nevents.at("pPb"), 62.6/nevents.at("Pbp"));
    eff.at("pA").at(bin)->SetTotalHistogram(hTotal, "f");
    eff.at("pA").at(bin)->SetPassedHistogram(hPassed, "f");
    if(TEfficiency::CheckWeights(hPassed, hTotal)) { eff.at("pA").at(bin)->SetStatisticOption(TEfficiency::kBJeffrey); }
  }
  
  // Draw and save efficiencies
  const std::string& fName = Form("rooFile/effDL_Thr*%.0f_%s%s.root", thr*100., type.c_str(), (muonTrig?"WT":"NT"));
  TFile outFile(fName.c_str(), "RECREATE");
  if(!outFile.isOpen() || outFile.isZombie()) { std::cout << "[ERROR] Output file: " << fName << " coud not be open!" << std::endl; return; }
  for(auto& c : eff)
  {
    for(auto& b : c.second)
    {
      const auto& col = c.first;
      const auto& bin = b.first;
      auto& effDL = b.second;
      const std::string& name = Form("Thr%.0f_%s%s_Rap_%.0f_%.0f_%s", thr*100., type.c_str(), (muonTrig?"WT":"NT"), yMin*10., yMax*10., s.c_str());
      TCanvas c; c.cd();
      effDL->Draw();
      c.Print(("fig_"+name+".png").c_str());
      outFile.cd();
      eff->Write(("effDL_"+name).c_str());
    }
  }
  outFile.Write();
  outFile.Close();
}
*/
