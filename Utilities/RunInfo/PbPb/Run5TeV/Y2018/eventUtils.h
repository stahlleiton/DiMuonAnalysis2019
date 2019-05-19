#ifndef PbPb_Run5TeV_Y2018_eventUtils_H_
#define PbPb_Run5TeV_Y2018_eventUtils_H_


namespace PbPb {
  namespace R5TeV {
    namespace Y2018 {
      // Trigger
      enum TRIGGERBIT {
        HLT_HIL1DoubleMuOpen_OS_Centrality_40_100 = 0,
        HLT_HIL1DoubleMuOpen_Centrality_50_100 = 1,
        HLT_HIL3Mu2p5NHitQ10_L2Mu2_M7toinf = 2,
        HLT_HIL1DoubleMu10 = 3,
        HLT_HIUPC_DoubleMu0_NotMBHF2AND = 4,
        HLT_HIL1MuOpen_Centrality_80_100 = 5,
        HLT_HIL3Mu12 = 6,
        HLT_HIUPC_SingleMuOpen_NotMBHF2AND = 7,
      };
      // Event Selection
      enum FILTERBIT {
        colEvtSel = 0,
        hfCoincFilter2Th4 = 1,
        primaryVertexFilter = 2,
        clusterCompatibilityFilter = 3,
      };
    };
  };  
};
  

#endif /* PbPb_Run5TeV_Y2018_eventUtils_H_ */
