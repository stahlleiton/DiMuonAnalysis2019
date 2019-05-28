# DiMuonAnalysis2019


## How to read the VertexComposite trees:

 A wrapper has been developed to initialize and extract information from the VertexComposite trees. The code is located in this repository in: [Utilities/Ntuple/VertexCompositeTree.h](Utilities/Ntuple/VertexCompositeTree.h)
  
 An example on how to use it is shown below:
 
```
#include "Utilities/Ntuple/VertexCompositeTree.h"
#include <iostream>
  
void test()
{
  const auto& inputFile = "/storage1/users/wl33/DiMuTrees/pPb2016/Tree/VertexCompositeTree_PADoubleMuon_PARun2016C_DiMuMassMin2.root";
  const auto& treeDir = "dimucontana"; // For MC use dimucontana_mc
  // Extract the tree
  VertexCompositeTree tree;
  if (!tree.GetTree(inputFile, treeDir)) { std::cout << "Invalid tree!" << std::endl; return; }
  // Loop over the tree
  for(Long64_t jentry=0; jentry<tree.GetEntries(); jentry++)
  {
    // Get the entry
    if (tree.GetEntry(jentry)<0) { std::cout << "Invalid entry!" << std::endl; return; }
    // Loop over the candidates
    for(uint iCand=0; iCand<tree.candSize(); iCand++)
    {
      // Print the candidate mass
      std::cout << "Candidate " << iCand << " has mass: " << tree.mass()[iCand] << std::endl;
    }
    // Stop after 10 events
    if (jentry>10) break;
  }
}
```
  run this example in ROOT using the command: root -l -b -q test.C
