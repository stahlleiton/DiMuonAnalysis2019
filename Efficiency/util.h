#include <vector>
#include <tuple>
#include <string>
#include <set>
#include <map>
#include <iostream>
#include <cmath>

#include "../Utilities/bin.h"

typedef std::tuple<std::string, std::string, std::vector<double>> Var_t;
typedef std::vector<Var_t> VarVec_t;
typedef std::map< BinF_t , VarVec_t  > BinBMap_t;
typedef std::map< BinF_t , BinBMap_t > BinBDiMap_t;

// Initialize the binning

// FOR Psi(2S):
BinBMap_t BINMAP_Psi2S_NTrk =
  {
   {
    {"Cand_AbsRap", 0.0, 1.4},
    {
     {"Cand_Pt", "VAR", {6.5, 9.0, 50.0}}
    }
   },
   {
    {"Cand_AbsRap", 1.4, 2.4},
    {
     {"Cand_Pt", "VAR", {3.0, 5.0, 6.5, 9.0, 50.0}}
    }
   },
   {
    {"Cand_Pt", 3.0, 6.5},
    {
     {"Cand_Rap", "VAR", {-2.4, -1.4, 1.4, 2.4}}
    }
   },
   {
    {"Cand_Pt", 6.5, 50.0},
    {
     {"Cand_RapCM", "VAR", {-2.86, -1.93, -1.2, -0.6, 0.0, 0.6, 1.2, 1.93}},
     {"Cand_RapCM", "VAR", {-1.93, 0.0, 1.93}},
     {"Cand_Rap", "VAR", {-2.4, -1.4, 1.4, 2.4}},
    }
   }
  };

BinBDiMap_t BINMAP_Psi2S =
  {
   {
    {"NTrack", 15.0, 250.0},
    {
     {
      {"Cand_AbsRap", 0.0, 1.4},
      {
       {"Cand_Pt", "VAR", {6.5, 8.0, 10.0, 12.0, 50.0}},
       {"Cand_Pt", "VAR", {6.5, 9.0, 50.0}}
      }
     },
     {
      {"Cand_AbsRap", 1.4, 2.4},
      {
       {"Cand_Pt", "VAR", {0.0, 1.5, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0, 12.0, 50.0}},
       {"Cand_Pt", "VAR", {3.0, 50.0}},
       {"Cand_Pt", "VAR", {3.0, 5.0, 6.5, 9.0, 50.0}}
      }
     },
     {
      {"Cand_Rap", -2.4, -1.4},
      {
       {"Cand_Pt", "VAR", {3.0, 6.5, 50.0}}
      }
     },
     {
      {"Cand_Rap", 1.4, 2.4},
      {
       {"Cand_Pt", "VAR", {3.0, 6.5, 50.0}}
      }
     },
     {
      {"Cand_Pt", 6.5, 50.0},
      {
       {"Cand_RapCM", "VAR", {-2.86, -2.4, -1.93, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 1.93}},
       {"Cand_RapCM", "VAR", {-2.86, -1.93, -1.2, -0.6, 0.0, 0.6, 1.2, 1.93}},
       {"Cand_RapCM", "VAR", {-1.93, 0.0, 1.93}}
      }
     }
    }
   },
   {
    {"NTrack", 15.0, 50.0},
    BINMAP_Psi2S_NTrk
   },
   {
    {"NTrack", 50.0, 80.0},
    BINMAP_Psi2S_NTrk
   },
   {
    {"NTrack", 80.0, 120.0},
    BINMAP_Psi2S_NTrk
   },
   {
    {"NTrack", 120.0, 250.0},
    BINMAP_Psi2S_NTrk
   }
  };


// FOR Psi(1S):
BinBMap_t BINMAP_Psi1S_NTrk =
  {
   {
    {"Cand_AbsRap", 0.0, 1.4},
    {
     {"Cand_Pt", "VAR", {6.5, 8.0, 10.0, 12.0, 50.0}}
    }
   },
   {
    {"Cand_AbsRap", 1.4, 2.4},
    {
     {"Cand_Pt", "VAR", {0.0, 1.5, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0, 12.0, 50.0}},
     {"Cand_Pt", "VAR", {3.0, 50.0}}
    }
   },
   {
    {"Cand_Rap", -2.4, -1.4},
    {
     {"Cand_Pt", "VAR", {0.0, 3.0, 6.5}}
    }
   },
   {
    {"Cand_Rap", 1.4, 2.4},
    {
     {"Cand_Pt", "VAR", {0.0, 3.0, 6.5}}
    }
   },
   {
    {"Cand_Pt", 6.5, 50.0},
    {
     {"Cand_Rap", "VAR", {-2.4, -1.4, -0.8, 0.0, 0.8, 1.4, 2.4}},
    }
   }
  };

BinBDiMap_t BINMAP_Psi1S =
  {
   {
    {"NTrack", 15.0, 250.0},
    BINMAP_Psi1S_NTrk
   },
   {
    {"NTrack", 15.0, 40.0},
    BINMAP_Psi1S_NTrk
   },
   {
    {"NTrack", 40.0, 60.0},
    BINMAP_Psi1S_NTrk
   },
   {
    {"NTrack", 60.0, 80.0},
    BINMAP_Psi1S_NTrk
   },
   {
    {"NTrack", 80.0, 100.0},
    BINMAP_Psi1S_NTrk
   },
   {
    {"NTrack", 100.0, 120.0},
    BINMAP_Psi1S_NTrk
   },
   {
    {"NTrack", 120.0, 150.0},
    BINMAP_Psi1S_NTrk
   },
   {
    {"NTrack", 150.0, 350.0},
    BINMAP_Psi1S_NTrk
   }
  };

// FOR TEST:

BinBDiMap_t BINMAP_General =
  {
   {
    {"NONE", 0.0, 0.0},
    {
     {
      {"NONE", 0.0, 0.0},
      {
       {"Cand_Pt",  "VAR", {0.0, 1.5, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0, 12.0, 50.0}},
       {"Cand_Rap", "FIX", {24, -2.4, 2.4}},
       {"NTrack", "FIX", {20, 0.0, 20.0}},
       {"NTrack", "VAR", {200, 0.0, 400.0}}
      }
     }
    }
   }
  };

BinBDiMap_t BINMAP_TEST =
  {
   {
    {"Cand_AbsRap", 0.0, 1.4},
    {
     {
      {"NTrack", 0.0, 10000.0},
      {
       {"Cand_Pt",  "VAR", {6.5,7,8,9,10,11,12,13,15,17,19,21,25,30,35,40,50}},
      }
     }
    }
   },
   {
    {"Cand_AbsRap", 1.4, 2.4},
    {
     {
      {"NTrack", 0.0, 10000.0},
      {
       {"Cand_Pt",  "VAR", {3,4,5,6.5,7,8,9,10,11,12,13,15,17,19,21,25,30,35,40,50}},
      }
     }
    }
   }
  };
