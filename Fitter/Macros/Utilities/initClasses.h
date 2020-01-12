#ifndef initClasses_h
#define initClasses_h

#include "RooWorkspace.h"
#include "../../../Utilities/dataUtils.h"

#include <string>
#include <unordered_map>
#include <vector>


typedef std::vector< GlobalInfo > GlobalInfoVector_t;
typedef std::unordered_map< std::string , GlobalInfoVector_t  > GlobalInfoVectorMap_t;
typedef std::unordered_map< std::string , RooWorkspace        > RooWorkspaceMap_t;
typedef std::vector<        StringVector_t                    > StringDiVector_t;
typedef std::unordered_map< std::string , StringDiVector_t    > StringDiVectorMap_t;
typedef std::unordered_map< std::string , StringVectorMap_t   > StringVectorDiMap_t;
typedef std::unordered_map< std::string , StringVectorDiMap_t > StringVectorTriMap_t;
typedef std::vector<        StringMap_t                       > StringMapVector_t;
typedef std::unordered_map< std::string , StringMap_t         > StringDiMap_t;
typedef std::vector<        StringDiMap_t                     > StringDiMapVector_t;
typedef std::unordered_map< std::string , IntMap_t            > IntDiMap_t;
typedef std::unordered_map< std::string , BoolMap_t           > BoolDiMap_t;
typedef std::unordered_map< std::string , std::vector< int >  > IntVecMap_t;
typedef std::unordered_map< std::string , IntVecMap_t         > IntVecDiMap_t;
typedef std::unordered_map< std::string , IntDiMap_t          > ModelMap;


const DoubleDiMap_t MASS = {
  { "D0",    {{"Val",  1.865}, {"Width", 0.10}, {"PID",    421.0}, {"Min",  1.7}, {"Max",   2.1}} },
  { "LambC", {{"Val",  2.286}, {"Width", 0.10}, {"PID",   4122.0}, {"Min",  2.0}, {"Max",   2.5}} },
  { "JPsi",  {{"Val",  3.096}, {"Width", 0.03}, {"PID",    443.0}, {"Min",  2.2}, {"Max",   3.9}} },
  { "Psi2S", {{"Val",  3.686}, {"Width", 0.05}, {"PID", 100443.0}, {"Min",  2.6}, {"Max",   4.7}} },
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
    SPLOT,
    // Candidate Mass Models
    SingleGaussian,
    DoubleGaussian,
    SingleCrystalBall,
    DoubleCrystalBall,
    GaussianAndCrystalBall,
    SingleExtCrystalBall,
    DoubleExtCrystalBall,
    GaussianAndExtCrystalBall,
    SingleModCrystalBall,
    DoubleModCrystalBall,
    GaussianAndModCrystalBall,
    SingleModExtCrystalBall,
    DoubleModExtCrystalBall,
    GaussianAndModExtCrystalBall,
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
    DeltaResolution,
    SingleGaussianResolution,
    DoubleGaussianResolution,
    TripleGaussianResolution,
    QuadrupleGaussianResolution,
    Delta,
    SingleSidedDecay,
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
  {"SPLOT",                       int(Model::SPLOT)},
  // Candidate Mass Models
  {"SingleGaussian",              int(Model::SingleGaussian)},
  {"DoubleGaussian",              int(Model::DoubleGaussian)},
  {"SingleCrystalBall",           int(Model::SingleCrystalBall)},
  {"DoubleCrystalBall",           int(Model::DoubleCrystalBall)},
  {"GaussianAndCrystalBall",      int(Model::GaussianAndCrystalBall)},
  {"SingleExtCrystalBall",        int(Model::SingleExtCrystalBall)},
  {"DoubleExtCrystalBall",        int(Model::DoubleExtCrystalBall)},
  {"GaussianAndExtCrystalBall",   int(Model::GaussianAndExtCrystalBall)},
  {"SingleModCrystalBall",        int(Model::SingleModCrystalBall)},
  {"DoubleModCrystalBall",        int(Model::DoubleModCrystalBall)},
  {"GaussianAndModCrystalBall",   int(Model::GaussianAndModCrystalBall)},
  {"SingleModExtCrystalBall",     int(Model::SingleModExtCrystalBall)},
  {"DoubleModExtCrystalBall",     int(Model::DoubleModExtCrystalBall)},
  {"GaussianAndModExtCrystalBall",int(Model::GaussianAndModExtCrystalBall)},
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
  {"DeltaResolution",             int(Model::DeltaResolution)},
  {"SingleGaussianResolution",    int(Model::SingleGaussianResolution)},
  {"DoubleGaussianResolution",    int(Model::DoubleGaussianResolution)},
  {"TripleGaussianResolution",    int(Model::TripleGaussianResolution)},
  {"QuadrupleGaussianResolution", int(Model::QuadrupleGaussianResolution)},
  {"Delta",                       int(Model::Delta)},
  {"SingleSidedDecay",            int(Model::SingleSidedDecay)},
  {"TripleDecay",                 int(Model::TripleDecay)},
  {"QuadrupleDecay",              int(Model::QuadrupleDecay)}
};


const StringMap_t VARLABEL_ = {
  { "Cand_Mass"    , "M"          },
  { "Cand_Rap"     , "y"          },
  { "Cand_RapCM"   , "y_{CM}"     },
  { "Cand_AbsRap"  , "|y|"        },
  { "Cand_Pt"      , "p_{T}"      },
  { "Cand_DLen"    , "#font[12]{l}"},
  { "Cand_DLenErr" , "#sigma_{#font[12]{l}}"},
  { "Cand_DLenRes" , "#frac{#font[12]{l}}{#sigma_{#font[12]{l}}}"},
  { "Cand_DLenGen" , "#font[12]{l}^{GEN}" },
  { "Cand_APhi"    , "#phi_{#mu}" },
  { "Centrality"   , "Cent."      },
  { "NTrack"       , "Ntrk"       },
};


#endif // #ifndef initClasses_h
