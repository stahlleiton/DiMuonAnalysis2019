#ifndef pp_Run13TeV_Y2018_eventUtils_H_
#define pp_Run13TeV_Y2018_eventUtils_H_



namespace pp {
  namespace R13TeV {
    namespace Y2018 {
      // Trigger
      enum TRIGGERBIT {
        HLT_FullTrack_Multiplicity85 = 0,
        HLT_FullTrack_Multiplicity100 = 1,
        HLT_FullTrack_Multiplicity130 = 2,
        HLT_FullTrack_Multiplicity155  = 3,
        HLT_L1MinimumBiasHF_OR  = 4,
	HLT_L1DoubleMu0 = 5,
      };
      // Event Selection
      enum FILTERBIT {
        colEvtSel = 0,
        primaryVertexFilter = 1,
        NoScraping = 2,
        pileupVertexFilterCut = 3,
        pileupVertexFilterCutGplus = 4,
      };
    };
  };
};


#endif /* pp_Run13TeV_Y2017_eventUtils_H_ */
