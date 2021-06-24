#ifndef Utilities_bins_bin_Psi1S_h
#define Utilities_bins_bin_Psi1S_h

#include "../../../Utilities/bin.h"

// Initialize the binning
const std::vector<float> NTRACKBINS_JPsi_ = { 15.0, 40.0, 60.0, 80.0, 100.0, 120.0, 150.0, 185.0, 250.0, 350.0 };

BinInfoDiMap_t BINMAP_Pt_JPsi_ =
  {
   {
    {{"Cand_AbsRap", 0.0, 1.4}},
    {
     {
      {{"Cand_Pt", 6.5, 8.0}},
      {
       {"NTrack", NTRACKBINS_JPsi_}
      }
     },
     {
      {{"Cand_Pt", 8.0, 10.0}},
      {
       {"NTrack", NTRACKBINS_JPsi_}
      }
     },
     {
      {{"Cand_Pt", 10.0, 12.0}},
      {
       {"NTrack", NTRACKBINS_JPsi_}
      }
     },
     {
      {{"Cand_Pt", 12.0, 50.0}},
      {
       {"NTrack", NTRACKBINS_JPsi_}
      }
     }
    }
   },
   {
    {{"Cand_AbsRap", 1.4, 2.4}},
    {
     {
      {{"Cand_Pt", 0.0, 1.5}},
      {
       {"NTrack", NTRACKBINS_JPsi_}
      }
     },
     {
      {{"Cand_Pt", 1.5, 3.0}},
      {
       {"NTrack", NTRACKBINS_JPsi_}
      }
     },
     {
      {{"Cand_Pt", 3.0, 4.0}},
      {
       {"NTrack", NTRACKBINS_JPsi_}
      }
     },
     {
      {{"Cand_Pt", 4.0, 5.0}},
      {
       {"NTrack", NTRACKBINS_JPsi_}
      }
     },
     {
      {{"Cand_Pt", 5.0, 6.5}},
      {
       {"NTrack", NTRACKBINS_JPsi_}
      }
     },
     {
      {{"Cand_Pt", 6.5, 8.0}},
      {
       {"NTrack", NTRACKBINS_JPsi_}
      }
     },
     {
      {{"Cand_Pt", 8.0, 10.0}},
      {
       {"NTrack", NTRACKBINS_JPsi_}
      }
     },
     {
      {{"Cand_Pt", 10.0, 12.0}},
      {
       {"NTrack", NTRACKBINS_JPsi_}
      }
     },
     {
      {{"Cand_Pt", 12.0, 50.0}},
      {
       {"NTrack", NTRACKBINS_JPsi_}
      }
     }
    }
   }
  };


BinInfoDiMap_t BINMAP_Rap_JPsi_ =
  {
   {
    {{"Cand_Pt", 6.5, 50.0}},
    {
     {
      {{"Cand_Rap", -2.4, -1.4}},
      {
       {"NTrack", NTRACKBINS_JPsi_},
      }
     },
     {
      {{"Cand_Rap", -1.4, -0.8}},
      {
       {"NTrack", NTRACKBINS_JPsi_},
      }
     },
     {
      {{"Cand_Rap", -0.8, 0.0}},
      {
       {"NTrack", NTRACKBINS_JPsi_},
      }
     },
     {
      {{"Cand_Rap", 0.0, 0.8}},
      {
       {"NTrack", NTRACKBINS_JPsi_},
      }
     },
     {
      {{"Cand_Rap", 0.8, 1.4}},
      {
       {"NTrack", NTRACKBINS_JPsi_},
      }
     },
     {
      {{"Cand_Rap", 1.4, 2.4}},
      {
       {"NTrack", NTRACKBINS_JPsi_},
      }
     }
    }
   }
  };


BinInfoDiMap_t BINMAP_FBRatio_JPsi_ =
  {
   {
    {{"Cand_Pt", 6.5, 50.0}},
    {
     {
      {{"Cand_Rap", 0.0, 0.8}},
      {
       {"NTrack", NTRACKBINS_JPsi_},
      }
     },
     {
      {{"Cand_Rap", 0.8, 1.4}},
      {
       {"NTrack", NTRACKBINS_JPsi_},
      }
     },
     {
      {{"Cand_Rap", 1.4, 2.4}},
      {
       {"NTrack", NTRACKBINS_JPsi_},
      }
     }
    }
   },
   {
    {{"Cand_Rap", 1.4, 2.4}},
    {
     {
      {{"Cand_Pt", 0.0, 3.0}},
      {
       {"NTrack", NTRACKBINS_JPsi_},
      }
     },
     {
      {{"Cand_Pt", 3.0, 6.5}},
      {
       {"NTrack", NTRACKBINS_JPsi_},
      }
     },
     {
      {{"Cand_Rap", 6.5, 50.0}},
      {
       {"NTrack", NTRACKBINS_JPsi_},
      }
     }
    }
   }
  };


BinInfoTriMap_t BINDIMAP_JPsi_ =
  {
   {"ForwardBackward_Ratio", { BINMAP_FBRatio_JPsi_ } },
   {"Cross_Section",         { BINMAP_Pt_JPsi_ , BINMAP_Rap_JPsi_ } },
   {"RatioTo1S",             { BINMAP_Pt_JPsi_ , BINMAP_Rap_JPsi_ } }
  };


#endif // #ifndef bin_Psi1S_h
