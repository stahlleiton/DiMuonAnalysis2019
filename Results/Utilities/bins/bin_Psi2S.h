#ifndef Utilities_bins_bin_Psi2S_h
#define Utilities_bins_bin_Psi2S_h

#include "../../../Utilities/bin.h"

// Initialize the binning
BinInfoDiMap_t BINMAP_NTrack_ =
  {
   {
    {{"Cand_AbsRap", 0.0, 1.4}},
    {
     {
      {{"NTrack", 15.0, 50.0}},
      {
       {"Cand_Pt", {6.5, 9.0, 50.0}}
      }
     },
     {
      {{"NTrack", 50.0, 80.0}},
      {
       {"Cand_Pt", {6.5, 9.0, 50.0}}
      }
     },
     {
      {{"NTrack", 80.0, 120.0}},
      {
       {"Cand_Pt", {6.5, 9.0, 50.0}}
      }
     },
     {
      {{"NTrack", 120.0, 250.0}},
      {
       {"Cand_Pt", {6.5, 9.0, 50.0}}
      }
     }
    }
   },
   {
    {{"Cand_AbsRap", 1.4, 2.4}},
    {
     {
      {{"NTrack", 15.0, 50.0}},
      {
       {"Cand_Pt", {3.0, 5.0, 6.5, 9.0, 50.0}}
      }
     },
     {
      {{"NTrack", 50.0, 80.0}},
      {
       {"Cand_Pt", {3.0, 5.0, 6.5, 9.0, 50.0}}
      }
     },
     {
      {{"NTrack", 80.0, 120.0}},
      {
       {"Cand_Pt", {3.0, 5.0, 6.5, 9.0, 50.0}}
      }
     },
     {
      {{"NTrack", 120.0, 250.0}},
      {
       {"Cand_Pt", {3.0, 5.0, 6.5, 9.0, 50.0}}
      }
     }
    }
   }
  };

BinInfoDiMap_t BINMAP_Pt_ =
  {
   {
    {{"Cand_AbsRap", 0.0, 1.4}},
    {
     {
      {{"Cand_Pt", 6.5, 9.0}},
      {
       {"NTrack", {15.0, 50.0, 80.0, 120.0, 250.0}}
      }
     },
     {
      {{"Cand_Pt", 9.0, 50.0}},
      {
       {"NTrack", {15.0, 50.0, 80.0, 120.0, 250.0}}
      }
     }
    }
   },
   {
    {{"Cand_AbsRap", 1.4, 2.4}},
    {
     {
      {{"Cand_Pt", 3.0, 5.0}},
      {
       {"NTrack", {15.0, 50.0, 80.0, 120.0, 250.0}}
      }
     },
     {
      {{"Cand_Pt", 5.0, 6.5}},
      {
       {"NTrack", {15.0, 50.0, 80.0, 120.0, 250.0}}
      }
     },
     {
      {{"Cand_Pt", 6.5, 9.0}},
      {
       {"NTrack", {15.0, 50.0, 80.0, 120.0, 250.0}}
      }
     },
     {
      {{"Cand_Pt", 9.0, 50.0}},
      {
       {"NTrack", {15.0, 50.0, 80.0, 120.0, 250.0}}
      }
     }
    }
   }
  };


BinInfoDiMap_t BINMAP_RapCM_ =
  {
   {
    {{"Cand_Pt", 6.5, 50.0}},
    {
     {
      {{"NTrack", 15.0, 50.0}},
      {
       {"Cand_RapCM", {-2.86, -1.93, -1.2, -0.6, 0.0, 0.6, 1.2, 1.93}},
      }
     },
     {
      {{"NTrack", 50.0, 80.0}},
      {
       {"Cand_RapCM", {-2.86, -1.93, -1.2, -0.6, 0.0, 0.6, 1.2, 1.93}},
      }
     },
     {
      {{"NTrack", 80.0, 120.0}},
      {
       {"Cand_RapCM", {-2.86, -1.93, -1.2, -0.6, 0.0, 0.6, 1.2, 1.93}},
      }
     },
     {
      {{"NTrack", 120.0, 250.0}},
      {
       {"Cand_RapCM", {-2.86, -1.93, -1.2, -0.6, 0.0, 0.6, 1.2, 1.93}},
      }
     }
    }
   }
  };

BinInfoDiMap_t BINMAP_FBRatio_ =
  {
   {
    {{"Cand_Pt", 6.5, 50.0}},
    {
     {
      {{"Cand_RapCM", 0.0, 0.6}},
      {
       {"NTrack", {15.0, 50.0, 80.0, 120.0, 250.0}},
      }
     },
     {
      {{"Cand_RapCM", 0.6, 1.2}},
      {
       {"NTrack", {15.0, 50.0, 80.0, 120.0, 250.0}},
      }
     },
     {
      {{"Cand_RapCM", 1.2, 1.93}},
      {
       {"NTrack", {15.0, 50.0, 80.0, 120.0, 250.0}},
      }
     }
    }
   }
  };

BinInfoTriMap_t BINDIMAP_ =
  {
   {"ForwardBackward_Ratio", { BINMAP_FBRatio_ } },
   {"Cross_Section",         { BINMAP_Pt_ , BINMAP_RapCM_ } },
   {"RatioTo1S",             { BINMAP_Pt_ , BINMAP_RapCM_ } }
  };


#endif // #ifndef bin_h
