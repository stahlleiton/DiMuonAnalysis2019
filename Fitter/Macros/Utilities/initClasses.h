#ifndef initClasses_h
#define initClasses_h

#include "RooWorkspace.h"

#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>


typedef std::unordered_map< std::string , RooWorkspace        > RooWorkspaceMap_t;
typedef std::unordered_set< std::string                       > StringSet_t;
typedef std::unordered_map< std::string , StringSet_t         > StringSetMap_t;
typedef std::vector<        std::string                       > StringVector_t;
typedef std::vector<        StringVector_t                    > StringDiVector_t;
typedef std::unordered_map< std::string , StringDiVector_t    > StringDiVectorMap_t;
typedef std::unordered_map< std::string , StringVector_t      > StringVectorMap_t;
typedef std::unordered_map< std::string , StringVectorMap_t   > StringVectorDiMap_t;
typedef std::unordered_map< std::string , StringVectorDiMap_t > StringVectorTriMap_t;
typedef std::unordered_map< std::string , double              > DoubleMap_t;
typedef std::unordered_map< std::string , DoubleMap_t         > DoubleDiMap_t;
typedef std::unordered_map< std::string , std::string         > StringMap_t;
typedef std::vector<        StringMap_t                       > StringMapVector_t;
typedef std::unordered_map< std::string , StringMap_t         > StringDiMap_t;
typedef std::vector<        StringDiMap_t                     > StringDiMapVector_t;
typedef std::unordered_map< std::string , int                 > IntMap_t;
typedef std::unordered_map< std::string , IntMap_t            > IntDiMap_t;
typedef std::unordered_map< std::string , bool                > BoolMap_t;
typedef std::unordered_map< std::string , BoolMap_t           > BoolDiMap_t;
typedef std::unordered_map< std::string , std::vector< int >  > IntVecMap_t;
typedef std::unordered_map< std::string , IntVecMap_t         > IntVecDiMap_t;
typedef std::unordered_map< std::string , IntDiMap_t          > ModelMap;


const DoubleDiMap_t MASS = {
  { "D0",    {{"Val",  1.865}, {"Width", 0.10}, {"PID",    421.0}, {"Min",  1.7}, {"Max",   2.1}} },
  { "LambC", {{"Val",  2.286}, {"Width", 0.10}, {"PID",   4122.0}, {"Min",  2.0}, {"Max",   2.5}} },
  { "JPsi",  {{"Val",  3.096}, {"Width", 0.05}, {"PID",    443.0}, {"Min",  2.5}, {"Max",   3.5}} },
  { "Psi2S", {{"Val",  3.686}, {"Width", 0.05}, {"PID", 100443.0}, {"Min",  3.3}, {"Max",   4.1}} },
  { "Ups1S", {{"Val",  9.460}, {"Width", 0.10}, {"PID",    553.0}, {"Min",  6.1}, {"Max",  10.6}} },
  { "Ups2S", {{"Val", 10.023}, {"Width", 0.10}, {"PID", 100553.0}, {"Min",  8.9}, {"Max",  11.1}} },
  { "Ups3S", {{"Val", 10.355}, {"Width", 0.10}, {"PID", 200553.0}, {"Min",  9.4}, {"Max",  14.9}} },
  { "Z",     {{"Val", 91.188}, {"Width", 3.00}, {"PID",     23.0}, {"Min", 60.0}, {"Max", 120.0}} },
};


// Fit Model
enum class Model 
{
    InvalidModel,
    // General Models
    CutAndCount,
    Template,
    // Candidate Mass Models
    SingleGaussian,
    DoubleGaussian,
    SingleCrystalBall,
    DoubleCrystalBall,
    GaussianAndCrystalBall,
    SingleModCrystalBall,
    DoubleModCrystalBall,
    GaussianAndModCrystalBall,
    Voigtian,
    BWCrystalBall,
    Uniform,
    Chebychev1,
    Chebychev2,
    Chebychev3,
    Chebychev4,
    Chebychev5,
    Chebychev6,
    Exponential,
    ExpChebychev1,
    ExpChebychev2,
    ExpChebychev3,
    ExpChebychev4,
    ExpChebychev5,
    ExpChebychev6,
    ExpError,
    // Candidate Decay Length Models
    SingleGaussianResolution,
    DoubleGaussianResolution,
    TripleGaussianResolution,
    QuadrupleGaussianResolution,
    Delta,
    SingleSidedDecay,
    DoubleSingleSidedDecay,
    TripleDecay,
    QuadrupleDecay,
    //
    Size
};
const IntMap_t ModelDictionary = {
  {"InvalidModel",                int(Model::InvalidModel)},
  // General Models
  {"CutAndCount",                 int(Model::CutAndCount)},
  {"Template",                    int(Model::Template)},
  // Candidate Mass Models
  {"SingleGaussian",              int(Model::SingleGaussian)},
  {"DoubleGaussian",              int(Model::DoubleGaussian)},
  {"SingleCrystalBall",           int(Model::SingleCrystalBall)},
  {"DoubleCrystalBall",           int(Model::DoubleCrystalBall)},
  {"GaussianAndCrystalBall",      int(Model::GaussianAndCrystalBall)},
  {"SingleModCrystalBall",        int(Model::SingleModCrystalBall)},
  {"DoubleModCrystalBall",        int(Model::DoubleModCrystalBall)},
  {"GaussianAndModCrystalBall",   int(Model::GaussianAndModCrystalBall)},
  {"Voigtian",                    int(Model::Voigtian)},
  {"BWCrystalBall",               int(Model::BWCrystalBall)},
  {"Uniform",                     int(Model::Uniform)},
  {"Chebychev1",                  int(Model::Chebychev1)},
  {"Chebychev2",                  int(Model::Chebychev2)},
  {"Chebychev3",                  int(Model::Chebychev3)},
  {"Chebychev4",                  int(Model::Chebychev4)},
  {"Chebychev5",                  int(Model::Chebychev5)},
  {"Chebychev6",                  int(Model::Chebychev6)},
  {"Exponential",                 int(Model::Exponential)},
  {"ExpChebychev1",               int(Model::ExpChebychev1)},
  {"ExpChebychev2",               int(Model::ExpChebychev2)},
  {"ExpChebychev3",               int(Model::ExpChebychev3)},
  {"ExpChebychev4",               int(Model::ExpChebychev4)},
  {"ExpChebychev5",               int(Model::ExpChebychev5)},
  {"ExpChebychev6",               int(Model::ExpChebychev6)},
  {"ExpError",                    int(Model::ExpError)},
  // Candidate Decay Length Models
  {"SingleGaussianResolution",    int(Model::SingleGaussianResolution)},
  {"DoubleGaussianResolution",    int(Model::DoubleGaussianResolution)},
  {"TripleGaussianResolution",    int(Model::TripleGaussianResolution)},
  {"QuadrupleGaussianResolution", int(Model::QuadrupleGaussianResolution)},
  {"Delta",                       int(Model::Delta)},
  {"SingleSidedDecay",            int(Model::SingleSidedDecay)},
  {"DoubleSingleSidedDecay",      int(Model::DoubleSingleSidedDecay)},
  {"TripleDecay",                 int(Model::TripleDecay)},
  {"QuadrupleDecay",              int(Model::QuadrupleDecay)}
};


const StringMap_t VARLABEL_ = {
  { "Cand_Mass"   , "M"         },
  { "Cand_Rap"    , "y"         },
  { "Cand_RapCM"  , "y_{CM}"    },
  { "Cand_AbsRap" , "|y|"       },
  { "Cand_Pt"     , "p_{T}"     },
  { "Cand_Len"    , "c#tau"     },
  { "Cand_APhi"   , "#phi_{#mu}"},
  { "Centrality"  , "Cent."     },
  { "NTrack"      , "Ntrk"      },
};


// Global Info Structure (wrapper to carry information around)
typedef struct GlobalInfo {
  DoubleDiMap_t     Var;
  StringMap_t       Par;
  IntMap_t          Int;
  StringVectorMap_t StrV;
  StringSetMap_t    StrS;
  BoolMap_t         Flag;
  void              Clear() { this->Var.clear(); this->Par.clear(); this->Int.clear(); this->StrV.clear(); this->StrS.clear(); this->Flag.clear(); }
  GlobalInfo() {}
  GlobalInfo(const GlobalInfo &ref, bool keep = false) {
    this->Copy(ref, keep);
  }
  ~GlobalInfo() {
    this->Clear();
  }
  void Copy(const DoubleDiMap_t &ref, bool keep = true) {
    if (!keep) this->Var.clear();
    for (const auto& var : ref) {
      for (const auto& ele : var.second) {
        this->Var[var.first][ele.first] = ele.second;
      }
    }
  }
  void Copy(const StringMap_t &ref, bool keep = true) {
    if (!keep) this->Par.clear();
    for (const auto& par : ref) {
      this->Par[par.first] = par.second;
    }
  }
  void Copy(const IntMap_t &ref, bool keep = true) {
    if (!keep) this->Int.clear();
    for (const auto& i : ref) {
      this->Int[i.first] = i.second;
    }
  }
  void Copy(const StringVectorMap_t &ref, bool keep = true) {
    if (!keep) this->StrV.clear();
    for (const auto& i : ref) {
      this->StrV[i.first] = i.second;
    }
  }
  void Copy(const StringSetMap_t &ref, bool keep = true) {
    if (!keep) this->StrS.clear();
    for (const auto& i : ref) {
      this->StrS[i.first] = i.second;
    }
  }
  void Copy(const BoolMap_t &ref, bool keep = true) {
    if (!keep) this->Flag.clear();
    for (const auto& flag : ref) {
      this->Flag[flag.first] = flag.second;
    }
  }
  void Copy(const GlobalInfo &ref, bool keep = true) {
    this->Copy(ref.Var, keep);
    this->Copy(ref.Par, keep);
    this->Copy(ref.Int, keep);
    this->Copy(ref.StrV, keep);
    this->Copy(ref.StrS, keep);
    this->Copy(ref.Flag, keep);
  }
  void Print(void) const 
  {
    for (const auto& var : this->Var) {
      for (const auto& ele : var.second) { std::cout << "VAR: " << var.first << " " << ele.first << " >> " << ele.second << std::endl; }
    }
    for (const auto& stV : this->StrV) {
      std::string n = "{ "; for (const auto& ele : stV.second) { n += ele+", "; } n += "}";
      std::cout << "STR: " << stV.first << " >> " << n << std::endl;
    }
    for (const auto& stS : this->StrS) {
      std::string n = "{ "; for (const auto& ele : stS.second) { n += ele+", "; } n += "}";
      std::cout << "STR: " << stS.first << " >> " << n << std::endl;
    }
    for (const auto& par : this->Par ) { std::cout << "PAR: "  << par.first << " >> " << par.second << std::endl; }
    for (const auto& inT : this->Int ) { std::cout << "INT: "  << inT.first << " >> " << inT.second << std::endl; }
    for (const auto& flg : this->Flag) { std::cout << "FLAG: " << flg.first << " >> " << flg.second << std::endl; }
  }
  bool operator == (const DoubleDiMap_t &ref) const
  {
    if (ref.size() != this->Var.size()) return false;
    for (const auto& var : this->Var) {
      if (ref.find(var.first)==ref.end() || ref.at(var.first).find("Min")==ref.at(var.first).end() || ref.at(var.first).find("Max")==ref.at(var.first).end()) return false;
      if (var.second.at("Min") != ref.at(var.first).at("Min")) return false;
      if (var.second.at("Max") != ref.at(var.first).at("Max")) return false;
    }
    return true;
  }
  bool operator == (const StringMap_t &ref) const
  {
    if (ref.size() != this->Par.size()) return false;
    for (const auto& par : this->Par) {
      if (ref.find(par.first)==ref.end()) return false;
      if (par.second != ref.at(par.first)) return false;
    }
    return true;
  }
  bool operator == (const IntMap_t &ref) const
  {
    if (ref.size() != this->Int.size()) return false;
    for (const auto& i : this->Int) {
      if (ref.find(i.first)==ref.end()) return false;
      if (i.second != ref.at(i.first)) return false;
    }
    return true;
  }
  bool operator == (const StringVectorMap_t &ref) const
  {
    if (ref.size() != this->StrV.size()) return false;
    for (const auto& i : this->StrV) {
      if (ref.find(i.first)==ref.end()) return false;
      if (i.second != ref.at(i.first)) return false;
    }
    return true;
  }
  bool operator == (const StringSetMap_t &ref) const
  {
    if (ref.size() != this->StrS.size()) return false;
    for (const auto& i : this->StrS) {
      if (ref.find(i.first)==ref.end()) return false;
      if (i.second != ref.at(i.first)) return false;
    }
    return true;
  }
  bool operator == (const BoolMap_t &ref) const
  {
    if (ref.size() != this->Flag.size()) return false;
    for (const auto& flag : this->Flag) {
      if (ref.find(flag.first)==ref.end()) return false;
      if (flag.second != ref.at(flag.first)) return false;
    }
    return true;
  }
  bool operator != ( const DoubleDiMap_t      &ref ) const { if (*this == ref) { return false; } return true; }
  bool operator != ( const StringMap_t        &ref ) const { if (*this == ref) { return false; } return true; }
  bool operator != ( const IntMap_t           &ref ) const { if (*this == ref) { return false; } return true; }
  bool operator != ( const StringVectorMap_t  &ref ) const { if (*this == ref) { return false; } return true; }
  bool operator != ( const StringSetMap_t     &ref ) const { if (*this == ref) { return false; } return true; }
  bool operator != ( const BoolMap_t          &ref ) const { if (*this == ref) { return false; } return true; }
  bool operator == ( const GlobalInfo         &ref) const
  {
    if ( *this != ref.Var  ) return false;
    if ( *this != ref.Par  ) return false;
    if ( *this != ref.Int  ) return false;
    if ( *this != ref.StrV ) return false;
    if ( *this != ref.StrS ) return false;
    if ( *this != ref.Flag ) return false;
    return true;
  }
} GlobalInfo;


typedef std::vector< GlobalInfo > GlobalInfoVector_t;
typedef std::unordered_map< std::string , GlobalInfoVector_t  > GlobalInfoVectorMap_t;


#endif // #ifndef initClasses_h
