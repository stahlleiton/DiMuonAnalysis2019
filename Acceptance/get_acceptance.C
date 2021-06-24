#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TFormula.h"
#include "TStyle.h"
#include "../Utilities/Ntuple/VertexCompositeTree.h"


typedef std::map<std::string, std::string> StrToStr;

double min(double a, double b){
  if(a<b) {return a;}
  else{return b;}
}

double GetZVal(TLorentzVector* mu, TH2D *h){
  int b = h->FindBin(mu->Eta() , fabs(mu->Pt()) ) ;
  return h->GetBinContent(b);
}

bool MuonKinAcceptance(float pt, float eta){
  if (fabs(eta)>2.4) {return false;}
  else if (fabs(eta)>2.1) {return (pt>1.8);}
  else if (fabs(eta)<1.2) {return (pt>3.5);}
  else {return (pt>(5.77-1.89*fabs(eta)));}
}

bool MuonKinAcceptance2(float pt, float eta){
  if (fabs(eta)>2.4) {return false;}
  else if (fabs(eta)>2.2) {return (pt>0.8);}
  else if (fabs(eta)<1.0) {return (pt>3.3);}
  else {return (pt>(2.9/std::cosh(eta)));}
}
/*
bool MuonKinAcceptance2(float pt, float eta){
  if (fabs(eta)>2.4) {return false;}
  else if (fabs(eta)>2.1) {return (pt>1.5);}
  else if (fabs(eta)<1.2) {return (pt>3.5);}
  else {return (pt>(5.47-1.89*fabs(eta)));}
}
*/

bool IsAcceptedMuon(TLorentzVector* mu, int ST, int SelType){
  bool res = false;
  double feta = fabs(mu->Eta());
  double pt = mu->Pt();

  if(SelType==0){
    //0 is: harsh global muons
    if((ST&2)>0 && feta<2.4){
      if(feta>=2.1 && pt>1.8){res = true;}
      else if(feta<1.2 && pt>3.5 ){res = true;}
      else if(feta>=1.2 && feta<2.1 && (pt > (5.77-1.8*feta)) ){res = true;}
    }
  }
  
  else if(SelType==1){
    //1 is: Loose, global muons
    if((ST&2)>0 && feta<2.4){
      if(feta<=1.0 && pt>3.3){res = true;}
      if(feta<=1.35 && feta>1.0 && (pt>6.73-3.43*feta) ){res = true;}
      if(feta>1.35 && (pt > (3.52-1.05*feta)) ){res = true;}
    }
  }

  else if(SelType==2){
    //2 is: Loose, tracker muons
    if((ST&8)>0 && feta<2.4){
      if(feta<=0.8 && pt>3.3){res = true;}
      if(feta<=1.5 && feta>0.8 && (pt>5.81-3.14*feta) ){res = true;}
      if(feta>1.5 && (pt > (1.89-0.526*feta)) ){res = true;}
    }
  }
  /*
  else if(SelType==2){
    //2 is: Loose, tracker muons
    if((ST&8)>0 && feta<2.4){
      if(feta<=0.8 && pt>3.3){res = true;}
      if(feta<=1.5 && feta>0.8 && (pt>5.81-3.14*feta) ){res = true;}
      if(feta>1.5 && (pt > (1.89-0.526*feta)) ){res = true;}
    }
  }
  */

  else{
    std::cout<<"Wrong muon selection type: accepted int's are 0 (harsh global), 1 (loose global), and 2 (loose tracker)"<<std::endl;}

  return res;
}

void get_acceptance()
{

  bool ispp = true;
  bool withTrig = false;
  bool L1trig = false;
  bool doTnP = true;
  bool doJpsiEff = false;
  bool doJpsiAcc = false;

  auto h_test = new TH1D();
  h_test->SetDefaultSumw2(true);

  std::string obj = "JPsi";
  //****** Open the tree and make it scan branches one by one (SetMakeClass, to study one branch at a time)
  const auto& inputFile = "/eos/cms/store/group/phys_heavyions/anstahll/RiceHIN/pPb2016/Tree/VertexCompositeTree_"+obj+"ToMuMu_pPb-Bst_pPb816Summer16_DiMuMC.root"
  const auto& treeDir = "dimucontana_mc"; // For MC use dimucontana_mc

  VertexCompositeTree tree;
  if (!tree.GetTree(inputFile, treeDir)) { std::cout << "Invalid tree!" << std::endl; return; } 
  int nevents= (int)tree.GetEntries();;
  std::cout<<"nevents = "<<nevents<<"\n";
  
  //****** Some parameters
  // double ScaleToXS = 208*208 * 460 * 2.54e-3 * 0.668 / (double)nevents; // A^2 * Lumi_PbPb[mub-1] * (XS_Bc_pp * BF((J/psi -> mu mu) mu nu))[mub] * (XS(5.02 TeV) / XS(7 TeV))
  // std::cout<<"Scaling the number of generated events by : "<<ScaleToXS<<std::endl;
  std::string MuTypeForEfficiency = "tracker"; //"all", "global", "nonglobal", "tracker"
  std::string QQSelectionAcc = "all";
  int trigbit = 1;

  //***** Some histograms
  int AEnbins = 27;
  double acceffBins[AEnbins+1]; 
  for(int l=0;l<19;l++) acceffBins[l] = l*0.5;
  for(int l=19;l<25;l++) acceffBins[l] = 9. + (l-19);
  acceffBins[25] = 16.;  acceffBins[26] = 19;  acceffBins[27] = 25;

  TH1D *h_JpsiNb = new TH1D("JpsiNb","Jpsi number",2,0,2);
  TH1D *h_recJpsi_Pt = new TH1D("recJpsiPt","transverse momentum of reconstructed J/#psi;P_{t}(J/#psi) [GeV]",30,2,20);
  TH1D *h_rec2Jpsi_Pt = new TH1D("rec2JpsiPt","transverse momentum of reconstructed J/#psi + new cuts;P_{t}(J/#psi) [GeV]",30,2,20);
  TH2D *h_recJpsi_YPt = new TH2D("recJpsiYPt","Reconstructed J/#psi;|Rapidity|;P_{t} [GeV]",13,0,2.6,AEnbins,acceffBins);
  TH1D *h_genJpsi_Pt = new TH1D("genJpsiPt","transverse momentum of observable J/#psi;P_{t}(J/#psi) [GeV]",30,2,20);
  TH1D *h_genJpsi_Y = new TH1D("genJpsiY","rapidity of observable J/#psi;|Rapidity|",50,0,3);
  TH2D *h_genJpsi_YPt = new TH2D("genJpsiYPt","All generated J/#psi;|Rapidity|;P_{t} [GeV]",13,0,2.6,AEnbins,acceffBins);
  int muptbins = doTnP?45:60;
  int muetabins = doTnP?52:52;
  double ptlow = doTnP?0.03:0; double pthigh = doTnP?6.03:6;
  TH2D *h_genmu_EtaPt = new TH2D("genmuEtaPt","All generated #mu;|#eta|;P_{t} [GeV]",muetabins,0,2.6,muptbins,ptlow,pthigh);
  TH2D *h_recmu_EtaPt = new TH2D("recmuEtaPt","All reconstructed #mu;|#eta|;P_{t} [GeV]",muetabins,0,2.6,muptbins,ptlow,pthigh);

  int NQQ_recoAccepted =0;
  int NQQ_accepted =0;
  int NQQ_reco =0;
  int genmuNb_goodTag = 0;
  int recmuNb_goodTag = 0;

  //***** Loop on events
  for(Long64_t i=0; i<nevents; i++)
  {
    if(!(i % 100000)) std::cout<<"Processed "<<i<<" events out of "<<nevents<<std::endl;
    
    if (tree.GetEntry(i)<0) { std::cout << "Invalid entry!" << std::endl; return; }

    const auto& Gen_weight = 1.0;

    //Only muons from a gen Jpsi
    if(doTnP || doJpsiEff){
      for(uint genQQNb=0; genQQNb<tree.candSize_gen(); genQQNb++){
	const auto& genPt = tree.pT_gen()[genQQNb];
	const auto& genRap = tree.y_gen()[genQQNb];
	
	const auto& PID = tree.PID_gen()[genQQNb];
	if (fabs(PID)!=443) continue;
	
	uint muPlI_gen = (tree.chargeD1_gen()[genQQNb]>0 ? 1 : 2);
	uint muMiI_gen = (tree.chargeD1_gen()[genQQNb]<0 ? 1 : 2);
	
	const auto& genMuPlPt = (muPlI_gen==1 ? tree.pTD1_gen()[genQQNb] : tree.pTD2_gen()[genQQNb]);
	const auto& genMuPlEta = (muPlI_gen==1 ? tree.EtaD1_gen()[genQQNb] : tree.EtaD2_gen()[genQQNb]);
	const auto& genMuPlPhi = (muPlI_gen==1 ? tree.PhiD1_gen()[genQQNb] : tree.PhiD2_gen()[genQQNb]);
	TLorentzVector genMuPlVec; genMuPlVec.SetPtEtaPhiM(genMuPlPt, genMuPlEta, genMuPlPhi, 0.0);
	
	const auto& genMuMiPt = (muMiI_gen==1 ? tree.pTD1_gen()[genQQNb] : tree.pTD2_gen()[genQQNb]);
	const auto& genMuMiEta = (muMiI_gen==1 ? tree.EtaD1_gen()[genQQNb] : tree.EtaD2_gen()[genQQNb]);
	const auto& genMuMiPhi = (muMiI_gen==1 ? tree.PhiD1_gen()[genQQNb] : tree.PhiD2_gen()[genQQNb]);
	TLorentzVector genMuMiVec; genMuMiVec.SetPtEtaPhiM(genMuMiPt, genMuMiEta, genMuMiPhi, 0.0);

	if(doTnP){
	  int tagIdx, probeIdx, tagGenIdx, probeGenIdx;
	  for (int k=0;k<2;k++){ //try both muons for the tag
	    
	    if(k==0){
	      tagGenIdx = muMiI_gen;
	      probeGenIdx = muPlI_gen;
	    } else{
	      tagGenIdx = muPlI_gen;
	      probeGenIdx = muMiI_gen;
	    }

	    const auto& genMuTagPt = (tagGenIdx==1 ? tree.pTD1_gen()[genQQNb] : tree.pTD2_gen()[genQQNb]);
	    const auto& genMuTagEta = (tagGenIdx==1 ? tree.EtaD1_gen()[genQQNb] : tree.EtaD2_gen()[genQQNb]);
	    const auto& genMuTagPhi = (tagGenIdx==1 ? tree.PhiD1_gen()[genQQNb] : tree.PhiD2_gen()[genQQNb]);
	    TLorentzVector genMuTagVec; genMuTagVec.SetPtEtaPhiM(genMuTagPt, genMuTagEta, genMuTagPhi, 0.0);
	    if (!MuonKinAcceptance(genMuTagPt, genMuTagEta)) continue;

	    // Step 1: Match the tag and probe to one of the reco muons
	    int tagIdx = -1;
	    for(uint recMuNb=0; recMuNb<tree.candSize_mu(); recMuNb++){

	      const auto& recMuPt = tree.pT_mu()[recMuNb];
	      const auto& recMuEta = tree.eta_mu()[recMuNb];
	      const auto& recMuPhi = tree.phi_mu()[recMuNb];
	      TLorentzVector recMuVec; recMuVec.SetPtEtaPhiM(recMuPt, recMuEta, recMuPhi, 0.0);

	      const bool isTagMu = (recMuVec.DeltaR(genMuTagVec)<0.05 && ((fabs(genMuTagPt-recMuPt)/genMuTagPt)<0.5));
	      if (isTagMu) { tagIdx = recMuNb; break; }
	    }
	    // Step2: Check the results
	    if(1==1
	       && (tagIdx>=0 ? tree.softMuon_mu()[tagIdx] : false)
	       && (tree.trigHLT()[7] && tree.trigMuon_mu()[7][tagIdx])
	       ){
	      
		//Fill Pt,Eta for gen muons
		const auto& genMuProbePt = (probeGenIdx==1 ? tree.pTD1_gen()[genQQNb] : tree.pTD2_gen()[genQQNb]);
		const auto& genMuProbeEta = (probeGenIdx==1 ? tree.EtaD1_gen()[genQQNb] : tree.EtaD2_gen()[genQQNb]);
		const auto& genMuProbePhi = (probeGenIdx==1 ? tree.PhiD1_gen()[genQQNb] : tree.PhiD2_gen()[genQQNb]);
		TLorentzVector genMuProbeVec; genMuProbeVec.SetPtEtaPhiM(genMuProbePt, genMuProbeEta, genMuProbePhi, 0.0);
		h_genmu_EtaPt->Fill(fabs(genMuProbeEta),genMuProbePt, Gen_weight);
		genmuNb_goodTag+=1;
		
		int probeIdx = -1;
		for(uint recMuNb=0; recMuNb<tree.candSize_mu(); recMuNb++){

		  const auto& recMuPt = tree.pT_mu()[recMuNb];
		  const auto& recMuEta = tree.eta_mu()[recMuNb];
		  const auto& recMuPhi = tree.phi_mu()[recMuNb];
		  TLorentzVector recMuVec; recMuVec.SetPtEtaPhiM(recMuPt, recMuEta, recMuPhi, 0.0);

		  const bool isProbeMu = (recMuVec.DeltaR(genMuProbeVec)<0.05 && ((fabs(genMuProbePt-recMuPt)/genMuProbePt)<0.5));
		  if (isProbeMu) { probeIdx = recMuNb; break; }
		}
	      
		if(// ((recmu_SelType[probeIdx])&(int)pow(2,12))>0 &&
		   (probeIdx>=0 ? tree.softMuon_mu()[probeIdx] : false)
		   && (withTrig ? (tree.trigHLT()[0] && (probeIdx>=0 ? tree.trigMuon_mu()[0][probeIdx] : false)) : true)
		   ){
		  h_recmu_EtaPt->Fill(fabs(genMuProbeEta),genMuProbePt, Gen_weight);
		  recmuNb_goodTag+=1;
		}
	    }
	    
	  }
	}


	if(doJpsiEff){
	  const int recQQIdx = 0;
	  if(MuonKinAcceptance2(tree.pTD1_gen()[genQQNb],tree.EtaD1_gen()[genQQNb]) && MuonKinAcceptance2(tree.pTD2_gen()[genQQNb],tree.EtaD2_gen()[genQQNb])
	     ){
	    h_genJpsi_YPt->Fill(genRap, genPt, Gen_weight);

	    if(
	       tree.softMuon1()[recQQIdx]
	       && tree.trigMuon1()[0][recQQIdx]
	       && tree.softMuon2()[recQQIdx]
	       && tree.trigMuon2()[0][recQQIdx]
	       ){
	      h_recJpsi_YPt->Fill(genRap, genPt, Gen_weight);
	    }
	  }
	}

      }
      ///////////////end Jpsi loop
    }
  }

  //****************************************************************
  //Jpsi acceptance
  //****************************************************************

  //***** Some histograms
  TH2D *h_genJpsi_YPt_all = new TH2D("genJpsiYPtAccAll","Gen J/#psi;|Rapidity|;P_{t} [GeV]",13,0,2.6,AEnbins,acceffBins);
  TH2D *h_genJpsi_YPt_acc = new TH2D("genJpsiYPtAcc","Gen J/#psi with accepted muons;|Rapidity|;P_{t} [GeV]",13,0,2.6,AEnbins,acceffBins);

  
  //****************************************************************
  //Drawing all
  //****************************************************************
  
  // std::cout<<"Number of reconstructed QQ: "<<NQQ_reco<<std::endl;
  // std::cout<<"Number of reconstructed QQ with accepted muons: "<<NQQ_recoAccepted<<std::endl;
  // std::cout<<"Number of QQ with accepted muons: "<<NQQ_accepted<<std::endl;

  std::cout<<"Number of generated probe muons from Jpsi, with a passing tag = "<<genmuNb_goodTag<<std::endl;
  std::cout<<"Number of reconstructed(+cuts) probe muons from Jpsi, with a passing tag = "<<recmuNb_goodTag<<std::endl;
  
  //  gStyle->SetPalette(kRainBow);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  gStyle->SetStatW(0.21);
  gStyle->SetStatH(0.16);
  gStyle->SetOptStat("nime");

  // Set Palette
  const Int_t NCont = 999;
  // const Int_t NRGBs = 7;
  // Double_t stops[NRGBs] = { 0.00, 0.0999, 0.1, 0.34, 0.61, 0.84, 1.00 };
  // Double_t red[NRGBs]   = { 0.99, 0.99, 0.0, 0.00, 0.87, 1.00, 0.51 };
  // Double_t green[NRGBs] = { 0.86, 0.9, 0.0, 0.81, 1.00, 0.20, 0.00 };
  // Double_t blue[NRGBs]  = { 0.99, 0.99, 0.0, 1.00, 0.12, 0.00, 0.00 };
  const Int_t NRGBs = 7;
  Double_t stops[NRGBs] = { 0.00, 0.025, 0.1, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.96, 0.99, 0.0, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.85, 0.0, 0.0, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.96, 0.99, 0.0, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  // TCanvas * c2 = new TCanvas("c2","pt acceptance",800,700);
  // h_genQQ_Pt->Draw("B");
  // TCanvas * c3 = new TCanvas("c3","Y acceptance",800,700);
  // h_genQQ_Y->Draw("B");

  gStyle->SetOptStat(0);
  StrToStr leg;
  leg["2global1tracker"] = "2 global + 1 tracker";
  leg["2harsh1tracker"] = "2 harsh + 1 tracker";
  leg["3global"] = "3 global";
  leg["3tracker"] = "3 tracker";

  // TCanvas * c1 = new TCanvas("c1","acceptance",800,700);
  // c1->SetRightMargin(0.115);
  // TH2D *h_genQQ_YPt_accRatio = (TH2D*)h_genQQ_YPt_acc->Clone("genQQYPtAccRatio");
  // h_genQQ_YPt_accRatio->SetTitle(("Generated B_{c} acceptance x efficiency ("+leg[QQSelectionAcc]+" muons);Gen Y;Gen P_{t} [GeV]").c_str()); //x efficiency
  // h_genQQ_YPt_accRatio->Divide(h_genQQ_YPt);
  // h_genQQ_YPt_accRatio->GetZaxis()->SetRangeUser(0,0.25);
  // h_genQQ_YPt_accRatio->Draw("COLZ");
    
  //***** Lines to draw geometrical acceptance
  TLine *l1 = new TLine(0, 3.3, 1., 3.3);
  TLine *l2 = new TLine(1., 3.3, 2.1, 1.2);
  TLine *l3 = new TLine(2.1, 1.2, 2.4, 1.2);
  TLine *l4 = new TLine(2.4, 1.2, 2.4, 6);
  if(MuTypeForEfficiency=="nonglobal" || MuTypeForEfficiency=="tracker"){
    l1->SetX1(0); l1->SetY1(3.3); l1->SetX2(0.8); l1->SetY2(3.3);
    l2->SetX1(0.8); l2->SetY1(3.3); l2->SetX2(1.5); l2->SetY2(1.1);
    l3->SetX1(1.5); l3->SetY1(1.1); l3->SetX2(2.4); l3->SetY2(0.6);
    l4->SetX1(2.4); l4->SetY1(5); l4->SetX2(2.4); l4->SetY2(0.6);
  }
  l1->SetLineColor(kGreen+1);
  l1->SetLineWidth(6);
  l2->SetLineColor(kGreen+1);
  l2->SetLineWidth(6);
  l3->SetLineColor(kGreen+1);
  l3->SetLineWidth(6);
  l4->SetLineColor(kGreen+1);
  l4->SetLineWidth(6);

  TLine *l13 = new TLine(0, 3.4, 0.3, 3.4);
  TLine *l13b = new TLine(0.3, 3.4, 0.3, 3.3);
  TLine *l13c = new TLine(0.3, 3.3, 1.1, 3.3);
  TLine *l14 = new TLine(1.1, 3.3, 1.4, 2.1);
  TLine *l15 = new TLine(1.4, 2.1, 1.55, 2.1);
  TLine *l16 = new TLine(1.55, 2.1, 2.2, 1.2);
  TLine *l17 = new TLine(2.2, 1.2, 2.4, 1.2);
  TLine *l18 = new TLine(2.4, 1.2, 2.4, 6);
  l13->SetLineColor(kGreen+1);
  l13->SetLineWidth(6);
  l13b->SetLineColor(kGreen+1);
  l13b->SetLineWidth(6);
  l13c->SetLineColor(kGreen+1);
  l13c->SetLineWidth(6);
  l14->SetLineColor(kGreen+1);
  l14->SetLineWidth(6);
  l15->SetLineColor(kGreen+1);
  l15->SetLineWidth(6);
  l16->SetLineColor(kGreen+1);
  l16->SetLineWidth(6);
  l17->SetLineColor(kGreen+1);
  l17->SetLineWidth(6);
  l18->SetLineColor(kGreen+1);
  l18->SetLineWidth(6);

  TLine *l9 = new TLine(0, 3.5, 1.2, 3.5);
  TLine *l10 = new TLine(1.2, 3.2, 2.1, 1.5);
  TLine *l11 = new TLine(2.1, 1.5, 2.4, 1.5);
  TLine *l12 = new TLine(2.4, 1.5, 2.4, 6);
  l9->SetLineColor(kGreen+1);
  l9->SetLineWidth(6);
  l10->SetLineColor(kGreen+1);
  l10->SetLineWidth(6);
  l11->SetLineColor(kGreen+1);
  l11->SetLineWidth(6);
  l12->SetLineColor(kGreen+1);
  l12->SetLineWidth(6);

  bool drawHarsh=true;
  TLine *l5 = new TLine(0, 3.5, 1.2, 3.5);
  TLine *l6 = new TLine(1.2, 3.5, 2.1, 1.8);
  TLine *l7 = new TLine(2.1, 1.8, 2.4, 1.8);
  TLine *l8 = new TLine(2.4, 1.8, 2.4, 6);
  if(drawHarsh){
    l5->SetLineColor(kBlack);
    l5->SetLineWidth(6);
    l6->SetLineColor(kBlack);
    l6->SetLineWidth(6);
    l7->SetLineColor(kBlack);
    l7->SetLineWidth(6);
    l8->SetLineColor(kBlack);
    l8->SetLineWidth(6);
  }

  TCanvas * c5 = new TCanvas("c5","c5",800,700);
  //h_recmu_EtaPt->Scale(ScaleToXS);
  h_recmu_EtaPt->Draw("COLZ");
  // l1->DrawClone("same");
  // l2->DrawClone("same");
  // l3->DrawClone("same");
  // l4->DrawClone("same");
  // if(drawHarsh){
  //   l5->DrawClone("same");
  //   l6->DrawClone("same");
  //   l7->DrawClone("same");
  //   l8->DrawClone("same");
  //  }

  TCanvas * c6 = new TCanvas("c6","c6",800,700);
  //  h_genmu_EtaPt->Scale(ScaleToXS);
  h_genmu_EtaPt->Draw("COLZ");

  TCanvas * c4 = new TCanvas("c4","muon efficiency",1600,1400);
  c4->SetRightMargin(0.115);
  TH2D *h_muAccEff_EtaPt = (TH2D*)h_recmu_EtaPt->Clone("muAccEffEtaPt");
  h_muAccEff_EtaPt->SetTitle("(Reco+Soft ID"+(TString)(withTrig?"+"+(TString)(L1trig?"L1 ":"")+"Trigger":"")+")/Generated muons;Gen |#eta|;Gen P_{t} [GeV]");// ("Reconstructed/Generated for "+MuTypeForEfficiency+" muons;Gen |#eta|;Gen P_{t} [GeV]"+"").c_str() 
  h_muAccEff_EtaPt->Divide(h_genmu_EtaPt);
  h_muAccEff_EtaPt->GetZaxis()->SetRangeUser(0,1);

  // for (int i=1;i<h_muAccEff_EtaPt->GetNbinsX();i++){
  //   for (int j=1;j<h_muAccEff_EtaPt->GetNbinsY();j++){
  //     if(h_muAccEff_EtaPt->GetBinContent(i,j) < 0.1) h_muAccEff_EtaPt->SetBinContent(i,j,0);
  //   }
  // }

  h_muAccEff_EtaPt->Draw("COLZ");

    /*
  if(drawHarsh){
    l5->DrawClone("same");
    l6->DrawClone("same");
    l7->DrawClone("same");
    l8->DrawClone("same");
  }
    */
  // l1->DrawClone("same");
  // l2->DrawClone("same");
  // l3->DrawClone("same");
  // l4->DrawClone("same");
  /*
  TF1 *pcut = new TF1("pcut","(x<1.0 ? 3.3 : (x<2.2 ? 2.9/cosh(x) : 0.8))",0,3);
  pcut->SetLineColor(kRed);
  pcut->SetLineWidth(6);
  pcut->Draw("same");
  */
  TF1 *pcut3 = new TF1("pcut3","(x<1.2 ? 3.3 : (x<2.1 ? 3.93-1.11*abs(x) : (x<2.4 ? 1.3 : 99.)))",0,3);
  pcut3->SetLineColor(kBlue);
  pcut3->SetLineWidth(6);
  pcut3->Draw("same");
  TF1 *pcut2 = new TF1("pcut2","(x<1.0 ? 3.3 : (x<1.5 ? 7.5-4.2*abs(x) : (x<2.4 ? ((2.4-0.8*abs(x))>0.8 ? (2.4-0.8*abs(x)) : 0.8) : 99.)))",0,3);
  pcut2->SetLineColor(kBlack);
  pcut2->SetLineWidth(6);
  pcut2->Draw("same");
  /*
  TF1 *pcut2 = new TF1("pcut2","(x<0.8 ? 3.3 : (x<1.5 ? 5.81-3.14*abs(x) : (x<2.4 ? ((1.89-0.526*abs(x))>0.8 ? (1.89-0.526*abs(x)) : 0.8) : 99.)))",0,3);
  pcut2->SetLineColor(kBlack);
  pcut2->SetLineWidth(6);
  pcut2->Draw("same");
  */
  
  if(doTnP){
    c4->SaveAs("SingleMuAcceptance_"+(TString)(withTrig?(L1trig?"L1step":"L3step"):"HybridSoftMuId")+"_TagPassedL3Mu3"+(TString)((ispp)?"_pPb2016":"")+"_"+obj+".png");
    c4->SaveAs("SingleMuAcceptance_"+(TString)(withTrig?(L1trig?"L1step":"L3step"):"HybridSoftMuId")+"_TagPassedL3Mu3"+(TString)((ispp)?"_pPb2016":"")+"_"+obj+".pdf");
  } else{
    c4->SaveAs("SingleMuAcceptance_ALLMUONS_HybridSoftMuId"+(TString)((ispp)?"_pPb2016":"")+".pdf");
    c4->SaveAs("SingleMuAcceptance_ALLMUONS_HybridSoftMuId"+(TString)((ispp)?"_pPb2016":"")+".png");
  }


  // gStyle->SetOptStat(0);
  if(doJpsiEff){
    TCanvas * c8 = new TCanvas("c8","Jpsi Pt Y efficiency",1600,1400);
    c8->SetRightMargin(0.115);
    TH1D *h_JpsiEff_YPt = (TH1D*)h_recJpsi_YPt->Clone("JpsiEffYPt");
    h_JpsiEff_YPt->SetTitle("Efficiency of J/#psi's with accepted muons;Gen |Rapidity|; Gen J/#psi P_{t} [GeV]");// ("Reconstructed/Generated for "+MuTypeForEfficiency+" muons;Gen |#eta|;Gen P_{t} [GeV]"+"").c_str() 
    h_JpsiEff_YPt->Divide(h_genJpsi_YPt);
    h_JpsiEff_YPt->GetZaxis()->SetRangeUser(0,1);
    h_JpsiEff_YPt->Draw("COLZ");
    c8->SaveAs("JpsiEfficiency_NewMuonKinCuts.png");
    c8->SaveAs("JpsiEfficiency_NewMuonKinCuts.pdf");

    if(doJpsiAcc){
      TCanvas * c12 = new TCanvas("c12","Jpsi Pt Y acceptance",1600,1400);
      c12->SetRightMargin(0.115);
      TH1D *h_JpsiAcc_YPt = (TH1D*)h_genJpsi_YPt_acc->Clone("JpsiAccYPt");
      h_JpsiAcc_YPt->SetTitle("Acceptance of generated J/#psi's;Gen |Rapidity|; Gen J/#psi P_{t} [GeV]");
      h_JpsiAcc_YPt->Divide(h_genJpsi_YPt_all);
      h_JpsiAcc_YPt->GetZaxis()->SetRangeUser(0,1);
      h_JpsiAcc_YPt->Draw("COLZ");
      c12->SaveAs("JpsiAcceptance_NewMuonKinCuts.png");
      c12->SaveAs("JpsiAcceptance_NewMuonKinCuts.pdf");

      TH1D *h_JpsiAccEff_YPt = (TH1D*)h_JpsiAcc_YPt->Clone("JpsiAccEffYPt");
      TCanvas * c13 = new TCanvas("c13","Jpsi Pt Y acceptance",1600,1400);
      h_JpsiAccEff_YPt->SetTitle("Acceptance #times Efficiency for generated J/#psi's;Gen |Rapidity|; Gen J/#psi P_{t} [GeV]");
      h_JpsiAccEff_YPt->Multiply(h_JpsiEff_YPt);
      h_JpsiAccEff_YPt->GetZaxis()->SetRangeUser(0,1);
      h_JpsiAccEff_YPt->Draw("COLZ");
      c13->SaveAs("JpsiAccEff_NewMuonKinCuts.png");
      c13->SaveAs("JpsiAccEff_NewMuonKinCuts.pdf");

    }
  }

  // TCanvas * c8 = new TCanvas("c8","Jpsi Pt efficiency",1600,1400);
  // c8->SetRightMargin(0.115);
  // TH1D *h_JpsiAccEff1_Pt = (TH1D*)h_recJpsi_Pt->Clone("JpsiAccEffPt");
  // h_JpsiAccEff1_Pt->SetTitle("Efficiency of good J/#psi's to pass dR3p5 trigger;Gen J/#psi P_{t} [GeV]");// ("Reconstructed/Generated for "+MuTypeForEfficiency+" muons;Gen |#eta|;Gen P_{t} [GeV]"+"").c_str() 
  // h_JpsiAccEff1_Pt->Divide(h_genJpsi_Pt);
  // h_JpsiAccEff1_Pt->GetZaxis()->SetRangeUser(0,1);
  // h_JpsiAccEff1_Pt->SetLineWidth(3);
  // h_JpsiAccEff1_Pt->Draw("hist");

  // TH1D *h_JpsiAccEff2_Pt = (TH1D*)h_rec2Jpsi_Pt->Clone("JpsiAccEff2Pt");
  // h_JpsiAccEff2_Pt->SetTitle("Efficiency of J/#psi passing all cuts + old/new muon acceptance;Gen J/#psi P_{t} [GeV]");// ("Reconstructed/Generated for "+MuTypeForEfficiency+" muons;Gen |#eta|;Gen P_{t} [GeV]"+"").c_str() 
  // h_JpsiAccEff2_Pt->Divide(h_genJpsi_Pt);
  // h_JpsiAccEff2_Pt->GetZaxis()->SetRangeUser(0,1);
  // h_JpsiAccEff2_Pt->SetLineWidth(3);
  // h_JpsiAccEff2_Pt->SetLineColor(kRed);
  // h_JpsiAccEff2_Pt->Draw("histsame");

  // auto legend = new TLegend(0.1,0.7,0.48,0.9);
  // legend->AddEntry(h_JpsiAccEff1_Pt,"old muon acceptance","l");
  // legend->AddEntry(h_JpsiAccEff2_Pt,"new muon acceptance","l");
  // legend->Draw();
  //  c8->SaveAs("JpsiEfficiency_dR3p5.png");

  // for(int bini=1;bini<=muptbins;bini++){
  //   for(int binj=1;binj<=muetabins;binj++){
  // 	if(h_muAccEff_EtaPt->GetBinContent(bini,binj)>1){
  // 	  std::cout<<"bin "<<bini<<", "<<binj<<" contains ratio >1"<<std::endl;
  // 	  std::cout<<"Number of reco muons in this bin : "<<h_recmu_EtaPt->GetBinContent(bini,binj)<<std::endl;
  // 	  std::cout<<"Number of gen muons in this bin : "<<h_genmu_EtaPt->GetBinContent(bini,binj)<<std::endl;
  // 	}
  //   }
  // }

  // TFile out_file("SingleMuonEfficiency_new.root","NEW");
  // h_muAccEff_EtaPt->Write();

}
