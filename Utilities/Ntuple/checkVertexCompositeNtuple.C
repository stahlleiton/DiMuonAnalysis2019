#include "TTree.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include <memory>


TTree* checkDuplicates(TTree* tree, std::map<int, std::set<int>>& lumiMap);
void saveJSON(const std::map<int, std::set<int>>& lumiMap, const std::string& fName);


void checkVertexCompositeNtuple(const std::string& filein)
{
  auto fin = std::unique_ptr<TFile>(TFile::Open(filein.c_str(), "READ"));
  if (!fin || !fin->IsOpen() || fin->IsZombie()) return;
  std::string oName = filein.substr(filein.rfind("/")+1); oName = oName.substr(0,oName.rfind(".root"))+"_TRIM.root";
  auto fout = std::unique_ptr<TFile>(TFile::Open(oName.c_str(),"RECREATE"));
  if (!fout || !fout->IsOpen() || fout->IsZombie()) return;

  std::map<int, std::set<int>> lumiMap_OS, lumiMap_SS;
  std::string name_OS = "dimucontana", name_SS = "dimucontana_wrongsign";
  if(!fin->Get((name_OS+"/VertexCompositeNtuple").c_str())) { name_OS = "dimucontana_mc"; }
  if(!fin->Get((name_SS+"/VertexCompositeNtuple").c_str())) { name_SS = "dimucontana_wrongsign_mc"; }

  fout->cd();
  const auto& tdir = fout->mkdir(name_OS.c_str());
  tdir->cd();
  const auto& tree = checkDuplicates(dynamic_cast<TTree*>(fin->Get((name_OS+"/VertexCompositeNtuple").c_str())), lumiMap_OS);
  
  fout->cd();
  const auto& tdir_ws = fout->mkdir(name_SS.c_str());
  tdir_ws->cd();
  const auto& tree_ws = checkDuplicates(dynamic_cast<TTree*>(fin->Get((name_OS+"/VertexCompositeNtuple").c_str())), lumiMap_SS);

  if (lumiMap_OS!=lumiMap_SS) { std::cout << "[ERROR] Luminosity maps are different between opposite and same sign vertex composite trees!" << std::endl; fout->Close(); fin->Close(); return; }

  fout->Write();
  fout->Close();
  fin->Close();

  std::string jName = filein.substr(filein.rfind("/")+1); jName = jName.substr(0,jName.rfind(".root"))+"_TREE.json";
  saveJSON(lumiMap_OS, jName);
};


void saveJSON(const std::map<int, std::set<int>>& lumiMap, const std::string& fName)
{
  ofstream jsonFile;
  jsonFile.open(fName);
  if (jsonFile.is_open()) {
    jsonFile << "{" << std::endl;
    bool last = false;
    for (const auto& ele : lumiMap) {
      jsonFile << "  \"" << ele.first << "\": [";
      bool first = true;
      for (auto i=ele.second.begin(); i!=ele.second.end(); i++) {
        if (i==ele.second.begin() || *i>*std::prev(i)+1) { jsonFile << (first ? "[" : ", [") << *i << ", "; if (first) { first=false; } }
        if (i==std::prev(ele.second.end()) || *i<*std::next(i)-1) { jsonFile << *i << "]"; }
      }
      if (ele.first==std::prev(lumiMap.end())->first) { last = true; }
      jsonFile << (last ? "]" : "],") << std::endl;
    }
    jsonFile << "}";
    jsonFile.close();
  }
  else { std::cout << "[ERROR] Unable to open JSON file!" << std::endl; }
};


TTree* checkDuplicates(TTree* tree, std::map<int, std::set<int>>& lumiMap)
{
  auto tnew = tree->CloneTree(0);
  tnew->SetAutoSave(0);
  tnew->SetAutoFlush(0);

  lumiMap.clear();
  uint RunNb, LSNb, EventNb;
  tree->SetBranchAddress("RunNb",&RunNb);
  tree->SetBranchAddress("LSNb",&LSNb);
  tree->SetBranchAddress("EventNb",&EventNb);
  typedef std::tuple<uint, uint, uint> EvtTuple;
  std::set<EvtTuple> evtInfo;

  const auto& nentries = tree->GetEntriesFast();
  for (int i=0; i<nentries; i++) {
    if (i%1000000==0) { std::cout << "Done: " << i << " / " << nentries << " (" << i*100./nentries << "%)" << std::endl; }
    tree->GetEntry(i);
    EvtTuple evtTuple(RunNb, LSNb, EventNb);
    if (evtInfo.find(evtTuple)!=evtInfo.end()) { std::cout << "Found duplicate in: " << RunNb << " , " << LSNb << " , " << EventNb << std::endl; }
    else { evtInfo.insert(evtTuple); tnew->Fill(); }
    lumiMap[RunNb].insert(LSNb);
  }

  return tnew;
};
