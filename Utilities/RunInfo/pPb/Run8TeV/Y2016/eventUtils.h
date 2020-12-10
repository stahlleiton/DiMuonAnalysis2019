#ifndef pPb_Run8TeV_Y2016_eventUtils_H_
#define pPb_Run8TeV_Y2016_eventUtils_H_

#include<map>

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
        HLT_INVALID = 7,
      };
      static const std::map< TRIGGERBIT , std::string > TRIGNAME =
	{
	 { HLT_PAL1DoubleMuOpen                      , "HLT_PAL1DoubleMuOpen"                      },
	 { HLT_PAL3Mu12                              , "HLT_PAL3Mu12"                              },
	 { HLT_PAFullTracks_Multiplicity120          , "HLT_PAFullTracks_Multiplicity120"          },
	 { HLT_PAFullTracks_Multiplicity150          , "HLT_PAFullTracks_Multiplicity150"          },
	 { HLT_PAFullTracks_Multiplicity185          , "HLT_PAFullTracks_Multiplicity185"          },
	 { HLT_PAFullTracks_Multiplicity250          , "HLT_PAFullTracks_Multiplicity250"          },
	 { HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack , "HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack" },
	};
      std::string HLTPath(const TRIGGERBIT& bit)
      {
	if (TRIGNAME.count(bit)) { return TRIGNAME.at(bit); }
	throw std::runtime_error(Form("[ERROR] HLT bit %d is invalid", bit));
	return "HLT_INVALID";
      }; 
      TRIGGERBIT HLTBit(const std::string& path)
      {
	for (const auto& trg : TRIGNAME) { if (trg.second==path) { return trg.first; } }
	throw std::runtime_error(Form("[ERROR] HLT path %s is invalid", path.c_str()));
	return HLT_INVALID;
      };
      std::map<uint, std::string> HLTBits()
      {
        std::map<uint, std::string> trigBits;
        for (const auto& t : TRIGNAME) { trigBits[t.first] = ""; }
        trigBits.at(HLT_PAL1DoubleMuOpen) = "DiMuon";
	trigBits.at(HLT_PAL3Mu12) = "Muon";
        return trigBits;
      };
      static const std::map< std::string , std::vector<uint> > PDTRIG =
	{
	 { "MUON"      , { HLT_PAL3Mu12 } },
	 { "DIMUON"    , { HLT_PAL1DoubleMuOpen } },
	 { "HIGHMULT0" , { HLT_PAFullTracks_Multiplicity150 } },
	 { "HIGHMULT"  , { HLT_PAFullTracks_Multiplicity185 } },
	 { "HIGHMULT2" , { HLT_PAFullTracks_Multiplicity250 } },
	 { "MINBIAS"   , { HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack } },
	};
      std::vector<uint> HLTBitsFromPD(const std::string& PD)
      {
	if (PDTRIG.count(PD)>0) { return PDTRIG.at(PD); }
	throw std::runtime_error("[ERROR] PD name " + PD + " is invalid");
	return {};
      };
      // Event Selection
      enum FILTERBIT {
        colEvtSel = 0,
        hfCoincFilter = 1,
        primaryVertexFilter = 2,
        NoScraping = 3,
        pileupVertexFilterCut = 4,
        pileupVertexFilterCutGplus =5,
	allEvtSel =6,
      };
      // Luminosity
      static const std::map< TRIGGERBIT , std::vector<double> > TRIGLUMI =
	{
	 { HLT_PAL1DoubleMuOpen                      , { 110.78 , 62.64 , 173.42 } },
	 { HLT_PAL3Mu12                              , { 110.56 , 62.65 , 173.20 } },
	 { HLT_PAFullTracks_Multiplicity120          , {   1.67 ,  0.92 ,   2.59 } },
	 { HLT_PAFullTracks_Multiplicity150          , {   3.17 ,  2.47 ,   5.64 } },
	 { HLT_PAFullTracks_Multiplicity185          , {  65.80 , 28.04 ,  93.84 } },
	 { HLT_PAFullTracks_Multiplicity250          , { 110.77 , 61.36 , 172.13 } },
	 { HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack , {   2.95 ,  1.09 ,   4.04 } },
	};
      double LumiFromHLTBit(const TRIGGERBIT& bit, const std::string& col)
      {
	if (TRIGLUMI.count(bit)) {
	  if      (col=="pPb8Y16") { return TRIGLUMI.at(bit)[0]; }
	  else if (col=="Pbp8Y16") { return TRIGLUMI.at(bit)[1]; }
	  else if (col=="PA8Y16" ) { return TRIGLUMI.at(bit)[2]; }
	  throw std::runtime_error(Form("[ERROR] System %s is invalid", col.c_str()));
	}
	throw std::runtime_error(Form("[ERROR] HLT bit %d is invalid", bit));
	return 0.0;
      };
      static const std::map< std::string , std::vector<double> > PDLUMI =
	{
	 { "MUON"      , TRIGLUMI.at(HLT_PAL3Mu12) },
	 { "DIMUON"    , TRIGLUMI.at(HLT_PAL1DoubleMuOpen) },
	 { "HIGHMULT0" , TRIGLUMI.at(HLT_PAFullTracks_Multiplicity150) },
	 { "HIGHMULT"  , TRIGLUMI.at(HLT_PAFullTracks_Multiplicity185) },
	 { "HIGHMULT2" , TRIGLUMI.at(HLT_PAFullTracks_Multiplicity250) },
	 { "MINBIAS"   , TRIGLUMI.at(HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack) },
	};
      double LumiFromPD(const std::string& PD, const std::string& col)
      {
	if (PDLUMI.count(PD)) {
	  if      (col=="pPb8Y16") { return PDLUMI.at(PD)[0]; }
	  else if (col=="Pbp8Y16") { return PDLUMI.at(PD)[1]; }
	  else if (col=="PA8Y16" ) { return PDLUMI.at(PD)[2]; }
	  throw std::runtime_error(Form("[ERROR] System %s is invalid", col.c_str()));
	}
	throw std::runtime_error(Form("[ERROR] PD %s is invalid", PD.c_str()));
	return 0.0;
      };
      double LumiWeightFromPD(const std::string& PD, const std::string& col, const std::string& smp)
      {
	if (PDLUMI.count(PD)) {
	  if (col=="pPb8Y16" || col=="Pbp8Y16") {
	    const double lumiWeight = LumiFromPD(PD, col) / LumiFromPD(PD, "PA8Y16");
	    double entries_pPb=0.0, entries_Pbp=0.0;
	    if (smp.rfind("MC_JPsi",0)==0 || smp.rfind("MC_PRJPsi",0)==0) { entries_pPb = 10573573; entries_Pbp = 10561630; }
	    else if (smp.rfind("MC_Psi2S",0)==0 || smp.rfind("MC_PRPsi2S",0)==0) { entries_pPb = 9700119; entries_Pbp = 8547402; }
	    else if (smp.rfind("MC_NoPRJPsi",0)==0 || smp.rfind("MC_NoPRPsi2S",0)==0) { entries_pPb = 14867044; entries_Pbp = 14521819; }
	    else { throw std::runtime_error(Form("[ERROR] Sample %s is invalid", smp.c_str())); }
	    const auto entries_pA = entries_pPb + entries_Pbp;
	    const double statWeight = (col=="Pbp8Y16" ? entries_Pbp : entries_pPb) / entries_pA;
	    return lumiWeight / statWeight;
	  }
	  throw std::runtime_error(Form("[ERROR] System %s is invalid", col.c_str()));
	}
	throw std::runtime_error(Form("[ERROR] PD %s is invalid", PD.c_str()));
	return 0.0;
      };
      // Muon acceptance cuts
      bool triggerMuonAcceptance(const double& pt, const double& eta)
      {
	return ( fabs(eta) < 2.4 &&
		 (    ( fabs(eta) < 1.2 && pt >= 3.3 ) ||
		      ( 1.2 <= fabs(eta) && fabs(eta) < 2.1 && pt >= 3.93-1.11*fabs(eta) ) ||
		      ( 2.1 <= fabs(eta) && fabs(eta) < 2.4 && pt >= 1.3 )
		      )
		 );
      };
      bool muonAcceptance(const double& pt, const double& eta)
      {
	return ( fabs(eta) < 2.4 &&
		 (    ( fabs(eta) < 1.0 && pt >= 3.3 ) ||
		      ( 1.0 <= fabs(eta) && fabs(eta) < 1.5 && pt >= 7.5-4.2*fabs(eta) ) ||
		      ( 1.5 <= fabs(eta) && fabs(eta) < 2.4 && pt >= std::max(2.4-0.8*fabs(eta), 0.8) )
		      )
		 );
      };
      std::string triggerMuonAcceptance(const std::string& pT, const std::string& eta)
      {
	std::string cutStr;
	cutStr += "(";
	cutStr += "(abs("+eta+") < 1.2 && "+pT+" >= 3.3) || ";
	cutStr += "(1.2 <= abs("+eta+") && abs("+eta+") < 2.1 && "+pT+" >= 3.93-1.11*abs("+eta+")) || ";
	cutStr += "(2.1 <= abs("+eta+") && abs("+eta+") < 2.4 && "+pT+" >= 1.3)";
	cutStr += ")";
	return cutStr;
      };
      std::string muonAcceptance(const std::string& pT, const std::string& eta)
      {
	std::string cutStr;
	cutStr += "(";
	cutStr += "(abs("+eta+") < 1.0 && "+pT+" >= 3.3) || ";
	cutStr += "(1.0 <= abs("+eta+") && abs("+eta+") < 1.5 && "+pT+" >= 7.5-4.2*abs("+eta+")) || ";
	cutStr += "(1.5 <= abs("+eta+") && abs("+eta+") < 2.4 && "+pT+" >= max(2.4-0.8*abs("+eta+"), 0.8))";
	cutStr += ")";
	return cutStr;
      };
    };
  };
};


#endif /* pPb_Run8TeV_Y2016_eventUtils_H_ */
