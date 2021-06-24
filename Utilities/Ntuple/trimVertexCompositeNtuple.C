#include "TTree.h"
#include "TFile.h"
#include <iostream>
#include <memory>


bool removeDuplicates(TTree& tree, std::map<int, std::set<int>>& lumiMap);

void trimVertexCompositeNtuple(const std::string& filein)
{
  auto fin = std::unique_ptr<TFile>(TFile::Open(filein.c_str(), "READ"));
  if (!fin || !fin->IsOpen() || fin->IsZombie()) return;
  auto fout = std::unique_ptr<TFile>(TFile::Open(("TRIM_"+filein).c_str(),"RECREATE"));
  if (!fout || !fout->IsOpen() || fout->IsZombie()) return;

  std::map<int, std::set<int>> lumiMap_OS, lumiMap_SS;
  std::string name_OS = "dimucontana", name_SS = "dimucontana_wrongsign";
  if(!fin->Get((name_OS+"/VertexCompositeNtuple").c_str())) { name_OS = "dimucontana_mc"; }
  if(!fin->Get((name_SS+"/VertexCompositeNtuple").c_str())) { name_SS = "dimucontana_wrongsign_mc"; }
  const auto& isDup_OS = checkDuplicates(*dynamic_cast<TTree*>(fin->Get((name_OS+"/VertexCompositeNtuple").c_str())), lumiMap_OS);
  if (isDup_OS) { std::cout << "[ERROR] Opposite sign vertex composite tree has duplicated events!" << std::endl; return; }
  const auto& isDup_SS = checkDuplicates(*dynamic_cast<TTree*>(fin->Get((name_OS+"/VertexCompositeNtuple").c_str())), lumiMap_SS);
  if (isDup_SS) { std::cout << "[ERROR] Same sign vertex composite tree has duplicated events!" << std::endl; return; }
  if (lumiMap_OS!=lumiMap_SS) { std::cout << "[ERROR] Luminosity maps are different between opposite and same sign vertex composite trees!" << std::endl; return; }
};


bool checkDuplicates(TTree& tree, std::map<int, std::set<int>>& lumiMap)
{
  lumiMap.clear();
  uint RunNb, LSNb, EventNb;
  tree.SetBranchStatus("*",0);
  tree.SetBranchStatus("RunNb",1);
  tree.SetBranchStatus("LSNb",1);
  tree.SetBranchStatus("EventNb",1);
  tree.SetBranchAddress("RunNb",&RunNb);
  tree.SetBranchAddress("LSNb",&LSNb);
  tree.SetBranchAddress("EventNb",&EventNb);
  typedef std::tuple<uint, uint, uint> EvtTuple;
  std::set<EvtTuple> evtInfo;

  const auto& nentries = tree.GetEntriesFast();
  for (int i=0; i<nentries; i++) {
    if (i%1000000==0) { std::cout << "Done: " << i << " / " << nentries << " (" << i*100./nentries << "%)" << std::endl; }
    tree.GetEntry(i);
    EvtTuple evtTuple(RunNb, LSNb, EventNb);
    if (evtInfo.find(evtTuple)!=evtInfo.end()) { std::cout << "Found duplicate in: " << RunNb << " , " << LSNb << " , " << EventNb << std::endl; return true; }
    else { evtInfo.insert(evtTuple); }
    lumiMap[RunNb].insert(LSNb);
  }
  return false;
};
