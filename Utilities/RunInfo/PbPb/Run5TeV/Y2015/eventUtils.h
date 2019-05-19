#ifndef PbPb_Run5TeV_Y2015_eventUtils_H_
#define PbPb_Run5TeV_Y2015_eventUtils_H_


namespace PbPb {
  namespace R5TeV {
    namespace Y2015 {
      // Trigger
      enum TRIGGERBIT {
        HLT_HIL1DoubleMu0 = 0,
        HLT_HIL1DoubleMu0_part = 1,
        HLT_HIL1DoubleMu0_2HF = 2,
        HLT_HIL1DoubleMu0_2HF0 = 3,
        HLT_HIL1DoubleMu0_2HF_Cent30100 = 4,
        HLT_HIL1DoubleMu0_2HF0_Cent30100 = 5,
        HLT_HIL1DoubleMu10 = 6,
        HLT_HIUPCL1DoubleMuOpenNotZDCAND = 7,
        HLT_HIUPCL1DoubleMuOpenNotHF2 = 8,
        HLT_HIL3Mu15 = 9,
        HLT_HIUPCL1MuOpenNotZDCAND = 10,
        HLT_HIUPCSingleMuNotHF2Pixel_SingleTrack = 11,
      };
      // Event Selection
      enum FILTERBIT {
        colEvtSel = 0,
        hfCoincFilter3 = 1,
        primaryVertexFilter = 2,
        clusterCompatibilityFilter = 3,
      };
    };
  };  
};
  

#endif /* PbPb_Run5TeV_Y2015_eventUtils_H_ */
