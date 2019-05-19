#ifndef pp_Run5TeV_Y2017_eventUtils_H_
#define pp_Run5TeV_Y2017_eventUtils_H_



namespace pp {
  namespace R5TeV {
    namespace Y2017 {
      // Trigger
      enum TRIGGERBIT {
        HLT_HIL1DoubleMu0 = 0,
        HLT_HIL3Mu12 = 1,
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


#endif /* pp_Run5TeV_Y2017_eventUtils_H_ */
