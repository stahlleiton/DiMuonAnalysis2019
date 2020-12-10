#ifndef ResultManager_h
#define ResultManager_h

// Auxiliary Headers
#include "resultsUtils.h"
#include "extractResultsTree.C"
#include "processResultsTree.C"
#include "plotResultsTree.C"


class ResultManager
{
 public:
  ResultManager();
  ~ResultManager() {};

  // setters
  void doSyst(const bool& s) { doSyst_ = s; };
  void setWorkDirInfo(const WSDirMap_t& w) { workDirInfo_ = w; };
  void setWorkDirInfo(const std::string& w) { workDirInfo_.at("Nominal").at("Nominal").first[0] = w; }
  void setFitDir(const std::string& mD, const StringVector_t& sD={""}, const StringVector_t& dL={""});
  void setFitInfo(const std::string& s, const std::string& v, const StringVector_t& t, const StringVector_t& c, const StringVector_t& o);
  void setEffDir(const std::string& mD);

  // functions
  bool extractAndProcess();
  void defineBins(const BinInfoTriMap_t& mainBinMap=BINDIMAP_);
  void plot();

 private:
  // attributes
  bool doSyst_, isDLCut_;
  WSDirMap_t workDirInfo_;
  std::string sample_, fitVariable_, mainDirectory_, preCWD_, CWD_, effDirectory_;
  StringVector_t triggers_, objects_, collisionSystems_, decayLengthDirectories_, subDirectories_;

  // containers
  VarBinQuadMap_t inputVarNom_;
  BinSeptaMap_t var_;
  BinOctaMapVec_t systVar_;
  BinCont_t binMap_;
  GraphSeptaMap_t graphMap_;

  // functions
  bool extract(VarBinQuadMap_t& inputVar, const std::string& wkDir);
  bool process(BinSextaMap_t& var, BinSeptaMapVec_t& systVar, VarBinQuadMap_t& inputVar, const bool& isNominal, const std::string& fitLabel, const IntMap_t& systCorr);
  const std::string nDir(const std::string& s) const { return ((s!="" && strcmp(&s.back(),"/")) ? (s + "/") : s); }
};


ResultManager::ResultManager()
{
  doSyst_ = false;
  isDLCut_ = false;
  workDirInfo_ = WSDirMap_t({{"Nominal", {{"Nominal", {{""}, {{"Rap", 1}, {"Obj", 1}}}}}}});
  CWD_ = getcwd(NULL, 0);
  preCWD_ = CWD_;
  preCWD_.erase(preCWD_.find_last_of("/"), 10);
};


void ResultManager::setFitDir(const std::string& mD, const StringVector_t& sD, const StringVector_t& dL)
{
  mainDirectory_ = mD;
  subDirectories_ = sD;
  decayLengthDirectories_ = dL;
  if (dL.size()==2) { isDLCut_ = true; }
};


void ResultManager::setFitInfo(const std::string& s, const std::string& v, const StringVector_t& t, const StringVector_t& c, const StringVector_t& o)
{
  sample_ = s;
  fitVariable_ = v;
  if (fitVariable_.find("_")!=std::string::npos) { fitVariable_.erase(fitVariable_.find("_"), 1); }
  triggers_ = t;
  collisionSystems_ = c;
  objects_ = o;
};


void ResultManager::setEffDir(const std::string& mD)
{
  effDirectory_ = mD;
};


bool ResultManager::extract(VarBinQuadMap_t& inputVar, const std::string& wkDir)
{
  for (const auto& dLDir : decayLengthDirectories_) {
    for (const auto& sDir : subDirectories_) {
      const auto& wsN = nDir(mainDirectory_) + nDir(dLDir) + nDir(wkDir) + nDir(sDir);
      for (const auto& trgTag : triggers_) {
	const auto dir = (preCWD_+"/Fitter/Output/"+nDir(wsN)+nDir(fitVariable_)+sample_+"_"+nDir(trgTag));
	if (!existDir(dir)) { std::cout << "[WARNING] Directory: " << dir << " does not exist, skipping!" << std::endl; continue; }
	for (const auto& colTag : collisionSystems_) {
	  for (const auto& objTag : objects_) {
	    const auto& dL = (dLDir!="" ? dLDir : "Inclusive");
	    if (!extractResultsTree(inputVar[dL], wsN, trgTag, colTag, objTag, sample_, fitVariable_)) { return false; }
	  }
	}
      }
    }
  }
  return true;
};


bool ResultManager::process(BinSextaMap_t& var, BinSeptaMapVec_t& systVar, VarBinQuadMap_t& inputVar, const bool& isNominal, const std::string& fitLabel, const IntMap_t& systCorr)
{
  if (isDLCut_) {
    std::string decayCutEfficiency = "Efficiency_DecayCut_90";
    if      (fitLabel=="Systematic_DecayCut_85"   ) { decayCutEfficiency = "Efficiency_DecayCut_85";    }
    else if (fitLabel=="Systematic_DecayCut_95"   ) { decayCutEfficiency = "Efficiency_DecayCut_95";    }
    else if (fitLabel=="Systematic_DecayCut_Alt"  ) { decayCutEfficiency = "Efficiency_DecayCut_Alt";   }
    else if (fitLabel=="Systematic_DecayCut_Psi2S") { decayCutEfficiency = "Efficiency_DecayCut_Psi2S"; }
    if (!processDecayCutResults(inputVar, decayCutEfficiency, effDirectory_)) { return false; }
  }
  const auto& inputVarNom = (!inputVarNom_.empty() ? inputVarNom_.at("Inclusive") : VarBinTriMap_t());
  if (!contain(inputVar, "Inclusive")) { std::cout << "[ERROR] inputVar does not contain Inclusive" << std::endl; return false; }
  if (!processResultsTree(var, systVar, inputVar.at("Inclusive"), doSyst_, isNominal, inputVarNom, systCorr, effDirectory_)) { return false; }
  return true;
};


bool ResultManager::extractAndProcess()
{
  // set nominal directory first and then the others
  StringVector_t workDirs({"Nominal"});
  for (const auto& t1 : workDirInfo_) {
    if (t1.first!="Nominal") {
      workDirs.push_back(t1.first);
    }
  }
  // loop over fit directories
  for (const auto& t1 : workDirs) {
    const bool isNominal = (t1=="Nominal");
    if (!doSyst_ && !isNominal) continue;
    for (const auto& t2 : workDirInfo_.at(t1)) {
      for (const auto& wkDir : t2.second.first) {
	std::cout << "[INFO] Adding results for : " << nDir(mainDirectory_) + wkDir << std::endl;
	// extract the information
	VarBinQuadMap_t inputVar;
	if (!extract(inputVar, wkDir)) { return false; }
	// process the information
	BinSextaMap_t var;
	BinSeptaMapVec_t systVar;
	if (!process(var, systVar, inputVar, isNominal, wkDir, t2.second.second)) { return false; }
	if (isNominal) {
	  var_["Nominal"] = var;
	  systVar_["Efficiency"] = systVar;
	  inputVarNom_ = inputVar;
        }
	else {
	  systVar_[t1][t2.first].push_back(var);
	}
      }
    }
  }
  // compute systematics
  if (doSyst_) { computeSystematic(var_, systVar_); }
  // done
  return true;
};


void ResultManager::defineBins(const BinInfoTriMap_t& mainBinMap)
{
  for (const auto& o : var_.at("Nominal")) {
    for (const auto& c : o.second) {
      for (const auto& v : c.second.begin()->second) {
	if (!contain(mainBinMap, v.first)) continue;
	for (const auto& varBinMap : mainBinMap.at(v.first)) {
	  for (const auto& incB : varBinMap) {
	    for (const auto& subB : incB.second) {
	      AnaBin_t anaBin(std::vector<AnaBin_t>({incB.first, subB.first}));
	      for (size_t iB=0; iB<subB.second.size(); iB++) {
		const auto& obs = subB.second[iB].first;
		const auto& bVec = subB.second[iB].second;
		std::set<AnaBin_t> binS;
		for (size_t i=1; i<bVec.size(); i++) {
		  auto bin = anaBin; bin.setbin(BinF_t(obs, bVec[i-1], bVec[i]));
		  const auto& vBin = v.second.at("Val").find(bin);
		  if (vBin!=v.second.at("Val").end()) { binS.insert(vBin->first); }
		  else { bin.print(); throw std::logic_error("[ERROR] Bin not found for "+o.first+" "+c.first+" "+v.first+" !"); }
		}
		if (!binS.empty()) {
		  binMap_[o.first][c.first][v.first][obs][incB.first][{subB.first, iB}] = binS;
		}
	      }
	    }
	  }
	}
      }
    }
  }
};


void ResultManager::plot()
{
  // define the output directory
  auto tTag = triggers_[0]; for (uint i=1; i<triggers_.size(); i++) { tTag += "_"+triggers_[i]; }
  auto cTag = collisionSystems_[0]; for (uint i=1; i<collisionSystems_.size(); i++) { cTag += "_"+collisionSystems_[i]; }
  auto oTag = objects_[0]; for (uint i=1; i<objects_.size(); i++) { oTag += "_"+objects_[i]; }
  const std::string& outDir = CWD_ + "/Output/Results/" + mainDirectory_+"/" + fitVariable_+"/" + sample_+"_"+tTag+"/" + oTag+"/" + cTag;
  const bool isMC = (sample_!="DATA");
  // plot results
  plotResultsTree(graphMap_, binMap_, var_.at("Nominal"), outDir, isMC, doSyst_);
};


#endif // ifndef ResultManager_h
