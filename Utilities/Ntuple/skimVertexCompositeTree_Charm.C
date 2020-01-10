#include "TTree.h"
#include "TFile.h"
#include <iostream>
#include <memory>


TTree* skimTree(TTree* told, std::vector<int>& evtV);


void skimVertexCompositeTree_Charm(const std::string& filein, const std::string& fileout)
{
  auto fin = std::unique_ptr<TFile>(TFile::Open(filein.c_str(), "READ"));
  auto fout = std::unique_ptr<TFile>(TFile::Open(fileout.c_str(),"RECREATE"));
  if (!fin || !fin->IsOpen() || fin->IsZombie()) return;
  if (!fout || !fout->IsOpen() || fout->IsZombie()) return;

  fout->cd();
  const auto& tdir = fout->mkdir("dimucontana");
  tdir->cd();
  std::vector<int> evtVec;
  const auto& tree = skimTree(dynamic_cast<TTree*>(fin->Get("dimucontana/VertexCompositeNtuple")), evtVec);

  fout->Write();
  fout->Close();
  fin->Close();
};


TTree* skimTree(TTree* told, std::vector<int>& evtV)
{
  const auto& tnew = told->CloneTree(0);
  tnew->SetAutoSave(0);
  tnew->SetAutoFlush(0);
  uint candSize;
  float mass[500], pTD1[500], pTD2[500], EtaD1[500], EtaD2[500];
  bool evtSel[10], trigHLT[20]={0}, softMuon1[500]={0}, softMuon2[500]={0};
  std::vector< std::vector< UChar_t > > *trigMuon1=0, *trigMuon2=0;

  told->SetBranchAddress("candSize",&candSize);
  told->SetBranchAddress("trigHLT",trigHLT);
  told->SetBranchAddress("mass",mass);
  told->SetBranchAddress("pTD1",pTD1);
  told->SetBranchAddress("pTD2",pTD2);
  told->SetBranchAddress("EtaD1",EtaD1);
  told->SetBranchAddress("EtaD2",EtaD2);
  told->SetBranchAddress("softMuon1",softMuon1);
  told->SetBranchAddress("softMuon2",softMuon2);
  told->SetBranchAddress("trigMuon1",&trigMuon1);
  told->SetBranchAddress("trigMuon2",&trigMuon2);
  if (!trigMuon1 || !trigMuon2) { std::cout << "WTF" << std::endl; return NULL; }

  const auto& nentries = told->GetEntriesFast();
  const bool& fillV = (int(evtV.size())<nentries);
  if (fillV) {
    evtV.clear(); evtV.reserve(nentries);
    for (int i=0; i<nentries; i++) {
      told->GetEntry(i);
      bool hasGoodCandidate = false;
      for (uint j=0; j<std::min(candSize,uint(500)); j++) {
	const auto& inMassPeak = ((mass[j]>2.1 && mass[j]<4.8));
	const auto& isMuon1Accp = inMassPeak && (fabs(EtaD1[j])<1.2 ? pTD1[j]>3.3 : (fabs(EtaD1[j])<2.1 ? pTD1[j]>(3.93-1.11*fabs(EtaD1[j])) : (fabs(EtaD1[j])<2.4 ? pTD1[j]>1.3 : false)));
        const auto& isMuon2Accp = inMassPeak && (fabs(EtaD2[j])<1.2 ? pTD2[j]>3.3 : (fabs(EtaD2[j])<2.1 ? pTD2[j]>(3.93-1.11*fabs(EtaD2[j])) : (fabs(EtaD2[j])<2.4 ? pTD2[j]>1.3 : false)));
        if (inMassPeak && (softMuon1[j] && softMuon2[j])) { hasGoodCandidate = true; break; }
      }
      const bool& keepEvt = hasGoodCandidate;
      evtV.push_back(keepEvt);
      if (keepEvt) tnew->Fill();
    }
  }
  else {
    for (int i=0; i<nentries; i++) { if (evtV[i]) { told->GetEntry(i); tnew->Fill(); } }
  }
  return tnew;
};
