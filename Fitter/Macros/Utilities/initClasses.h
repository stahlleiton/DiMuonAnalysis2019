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
typedef std::unordered_map< std::string , IntDiMap_t          > ModelMap_t;
typedef std::pair<                  int , std::string         > IntPair_t;
typedef std::unordered_map< std::string , IntPair_t           > IntPairMap_t;


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
    FULL,
    // Candidate Mass Models
    SingleGaussian,
    DoubleGaussian,
    TripleGaussian,
    QuadrupleGaussian,
    PentaGaussian,
    SingleCrystalBall,
    DoubleCrystalBall,
    GaussianAndCrystalBall,
    DoubleGaussianAndCrystalBall,
    SingleExtCrystalBall,
    DoubleExtCrystalBall,
    GaussianAndExtCrystalBall,
    DoubleGaussianAndExtCrystalBall,
    SingleModCrystalBall,
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
    PentaGaussianResolution,
    Delta,
    SingleSidedDecay,
    DoubleSingleSidedDecay,
    TripleDecay,
    QuadrupleDecay,
    QuadrupleDecay2,
    QuadrupleDecay3,
    //
    Size
};
const IntPairMap_t ModelDictionary =
{
  {"InvalidModel",                {int(Model::InvalidModel), "INV"}},
  // General Models
  {"CutAndCount",                 {int(Model::CutAndCount), "CutNCnt"}},
  {"Template",                    {int(Model::Template), "Temp"}},
  {"SPLOT",                       {int(Model::SPLOT), "SPLOT"}},
  {"FULL",                        {int(Model::FULL), "FULL"}},
  // Candidate Mass Models
  {"SingleGaussian",              {int(Model::SingleGaussian), "SngGaus"}},
  {"DoubleGaussian",              {int(Model::DoubleGaussian), "DobGaus"}},
  {"TripleGaussian",              {int(Model::TripleGaussian), "TriGaus"}},
  {"QuadrupleGaussian",           {int(Model::QuadrupleGaussian), "QuadGaus"}},
  {"PentaGaussian",               {int(Model::PentaGaussian), "PenGaus"}},
  {"SingleCrystalBall",           {int(Model::SingleCrystalBall), "SngCrysB"}},
  {"DoubleCrystalBall",           {int(Model::DoubleCrystalBall), "DobCrysB"}},
  {"GaussianAndCrystalBall",      {int(Model::GaussianAndCrystalBall), "GausNCrysB"}},
  {"DoubleGaussianAndCrystalBall",{int(Model::DoubleGaussianAndCrystalBall), "DobGausNCrysB"}},
  {"SingleExtCrystalBall",        {int(Model::SingleExtCrystalBall), "SngExtCrysB"}},
  {"DoubleExtCrystalBall",        {int(Model::DoubleExtCrystalBall), "DobExtCrysB"}},
  {"GaussianAndExtCrystalBall",   {int(Model::GaussianAndExtCrystalBall), "GausNExtCrysB"}},
  {"DoubleGaussianAndExtCrystalBall",{int(Model::DoubleGaussianAndExtCrystalBall), "DobGausNExtCrysB"}},
  {"SingleModCrystalBall",        {int(Model::SingleModCrystalBall), "SngModCrysB"}},
  {"SingleModExtCrystalBall",     {int(Model::SingleModExtCrystalBall), "SngModExtCrysB"}},
  {"DoubleModExtCrystalBall",     {int(Model::DoubleModExtCrystalBall), "DobModExtCrysB"}},
  {"GaussianAndModExtCrystalBall",{int(Model::GaussianAndModExtCrystalBall), "GausNModExtCrysB"}},
  {"Voigtian",                    {int(Model::Voigtian), "Voigt"}},
  {"BWCrystalBall",               {int(Model::BWCrystalBall), "BWCrysB"}},
  {"Uniform",                     {int(Model::Uniform), "Uniform"}},
  {"Chebychev1",                  {int(Model::Chebychev1), "Cheb1"}},
  {"Chebychev2",                  {int(Model::Chebychev2), "Cheb2"}},
  {"Chebychev3",                  {int(Model::Chebychev3), "Cheb3"}},
  {"Chebychev4",                  {int(Model::Chebychev4), "Cheb4"}},
  {"Chebychev5",                  {int(Model::Chebychev5), "Cheb5"}},
  {"Chebychev6",                  {int(Model::Chebychev6), "Cheb6"}},
  {"Exponential",                 {int(Model::Exponential), "Exp"}},
  {"ExpChebychev1",               {int(Model::ExpChebychev1), "ExpCheb1"}},
  {"ExpChebychev2",               {int(Model::ExpChebychev2), "ExpCheb2"}},
  {"ExpChebychev3",               {int(Model::ExpChebychev3), "ExpCheb3"}},
  {"ExpChebychev4",               {int(Model::ExpChebychev4), "ExpCheb4"}},
  {"ExpChebychev5",               {int(Model::ExpChebychev5), "ExpCheb5"}},
  {"ExpChebychev6",               {int(Model::ExpChebychev6), "ExpCheb6"}},
  {"ExpError",                    {int(Model::ExpError), "ExpError"}},
  // Candidate Decay Length Models
  {"DeltaResolution",             {int(Model::DeltaResolution), "DeltaRes"}},
  {"SingleGaussianResolution",    {int(Model::SingleGaussianResolution), "SngGaussRes"}},
  {"DoubleGaussianResolution",    {int(Model::DoubleGaussianResolution), "DobGaussRes"}},
  {"TripleGaussianResolution",    {int(Model::TripleGaussianResolution), "TriGaussRes"}},
  {"QuadrupleGaussianResolution", {int(Model::QuadrupleGaussianResolution), "QuadGaussRes"}},
  {"PentaGaussianResolution",     {int(Model::PentaGaussianResolution), "PenGaussRes"}},
  {"Delta",                       {int(Model::Delta), "Delta"}},
  {"SingleSidedDecay",            {int(Model::SingleSidedDecay), "SngSideDec"}},
  {"DoubleSingleSidedDecay",      {int(Model::DoubleSingleSidedDecay), "DobSngSideDec"}},
  {"TripleDecay",                 {int(Model::TripleDecay), "TriSideDec"}},
  {"QuadrupleDecay",              {int(Model::QuadrupleDecay), "QuadDec"}},
  {"QuadrupleDecay2",             {int(Model::QuadrupleDecay2), "QuadDec2"}},
  {"QuadrupleDecay3",             {int(Model::QuadrupleDecay3), "QuadDec3"}}
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
