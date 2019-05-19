#ifndef pPb_Run8TeV_Y2016_eventUtils_H_
#define pPb_Run8TeV_Y2016_eventUtils_H_



namespace pPb {
  namespace R8TeV {
    namespace Y2016 {
      // Trigger
      enum TRIGGERBIT {
        HLT_PAL1DoubleMuOpen = 0,
        HLT_PAL3Mu12 = 1,
        HLT_PAFullTracks_Multiplicity120 = 2,
        HLT_PAFullTracks_Multiplicity150  = 3,
        HLT_PAFullTracks_Multiplicity185  = 4,
        HLT_PAFullTracks_Multiplicity250 = 5,
        HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack = 6,
      };
      // Event Selection
      enum FILTERBIT {
        colEvtSel = 0,
        hfCoincFilter = 1,
        primaryVertexFilter = 2,
        NoScraping = 3,
        pileupVertexFilterCut = 4,
        pileupVertexFilterCutGplus =5,
      };
    };
  };
};


#endif /* pPb_Run8TeV_Y2016_eventUtils_H_ */
