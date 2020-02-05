//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

#ifndef DATAUTILS_H_
#define DATAUTILS_H_

// ROOT Headers
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"
#include "TFile.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"

// c++ headers
#include <dirent.h>
#include <string>
#include <sys/stat.h>
#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <chrono>
#include <memory>
#include <utility>
#include <algorithm>

// Auxiliary Headers
#include "RunInfo/eventUtils.h"


typedef std::set< std::string                  > StringSet_t;
typedef std::map< std::string , StringSet_t    > StringSetMap_t;
typedef std::vector< std::string               > StringVector_t;
typedef std::map< std::string , StringVector_t > StringVectorMap_t;
typedef std::map< std::string , double         > DoubleMap_t;
typedef std::map< std::string , DoubleMap_t    > DoubleDiMap_t;
typedef std::map< std::string , std::string    > StringMap_t;
typedef std::map< std::string , std::string*   > StringPMap_t;
typedef std::map< std::string , int            > IntMap_t;
typedef std::map< std::string , bool           > BoolMap_t;


// Global Info Structure (wrapper to carry information around)
typedef struct GlobalInfo {
  DoubleDiMap_t     Var;
  StringMap_t       Par;
  IntMap_t          Int;
  StringVectorMap_t StrV;
  StringSetMap_t    StrS;
  StringPMap_t      StrP;
  BoolMap_t         Flag;
  void              Clear() { this->Var.clear(); this->Par.clear(); this->Int.clear(); this->StrV.clear(); this->StrS.clear(); this->StrP.clear(); this->Flag.clear(); }
  GlobalInfo() {}
  GlobalInfo(const GlobalInfo &ref, bool keep = false) {
    this->Copy(ref, keep);
  }
  ~GlobalInfo() {
    for (auto& p : this->StrP) { if (p.second) { delete p.second; } }
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
    for (const auto& stP : this->StrP) {
      if (stP.second) { std::cout << "STR: " << stP.first << " >> "<< *stP.second << std::endl; }
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


// Utiliy Functions


template<class C, class T>
inline auto contain_impl(const C& c, const T& x, int)
-> decltype(c.find(x), true)
{ return std::end(c) != c.find(x); }

template<class C, class T>
inline bool contain_impl(const C& v, const T& x, long)
{ return std::end(v) != std::find(std::begin(v), std::end(v), x); }

template<class C, class T>
auto contain(const C& c, const T& x)
-> decltype(std::end(c), true)
{ return contain_impl(c, x, 0); }

template <typename  T>
std::vector<T> merge(const std::vector<T>& s, const std::vector<T>& e)
{
  auto v = s;
  v.insert(s.end(), e.begin(), e.end());
  return v;
}

template <typename  T>
std::set<T> merge(const std::set<T>& s, const std::set<T>& e)
{
  auto v = s;
  std::merge(s.begin(), s.end(), e.begin(), e.end(), std::inserter(v, v.begin()));
  return v;
}


bool fileList(std::vector< std::string >& fileNames, const std::string& dirPath, const bool& verb=false)
{
  // Open the directory
  DIR * dpdf = opendir(dirPath.c_str());
  // Search for all the files inside the directory
  if (dpdf != NULL){
    struct dirent *epdf;
    while ((epdf = readdir(dpdf))){
      if (strcmp(epdf->d_name,".")!=0 && strcmp(epdf->d_name,"..")!=0 ) {
        if (verb) { std::cout << "[INFO] Adding file: " << epdf->d_name << std::endl; }
        fileNames.push_back(epdf->d_name);
      }
    }
  } else {
    std::cout << "[ERROR] Working directory ( " << dirPath << " ) was not found!" << endl; return false;
  }
  return true;
};


bool existDir(const std::string& dir)
{
  bool exist = false;
  const auto& dirp = gSystem->OpenDirectory(dir.c_str());
  if (dirp) { gSystem->FreeDirectory(dirp); exist = true; }
  return exist;
};


bool existFile(const std::string& file)
{
  struct stat buffer;
  return (stat (file.c_str(), &buffer) == 0);
};


void makeDir(const std::string& dir)
{
  if (existDir(dir)==false){
    std::cout << "[INFO] DataSet directory: " << dir << " doesn't exist, will create it!" << std::endl;  
    gSystem->mkdir(dir.c_str(), kTRUE);
  }
};


void findSubDir(std::vector<std::string>& dirList, const std::string dirName)
{
  TSystemDirectory dir(dirName.c_str(), dirName.c_str());
  auto subdirs = std::unique_ptr<TList>(dir.GetListOfFiles());
  if (subdirs) {
    TSystemFile *subdir;
    TIter next(subdirs.get());
    while ( (subdir=dynamic_cast<TSystemFile*>(next())) ) {
      if (subdir->IsDirectory() && std::string(subdir->GetName())!="." && std::string(subdir->GetName())!="..") {
	const auto& subDirName = dirName + subdir->GetName() + "/";
        dirList.push_back(subDirName);
        std::cout << "[INFO] Input subdirectory: " << subDirName << " found!" << std::endl;
      }
    }
  }
};


void findDirInFile(std::string& dirName, const std::string& fileName)
{
  dirName = "";
  auto fName = fileName;
  if (fName.rfind("/store/", 0)==0) { fName = "root://cms-xrd-global.cern.ch/" + fName; }
  auto file = std::unique_ptr<TFile>(TFile::Open(fName.c_str()));
  if (file && file->IsOpen() && !file->IsZombie()) {
    const auto& dir = static_cast<TDirectory*>(file->GetListOfKeys()->First());
    if (dir) { dirName = dir->GetName(); }
  }
  if (file) file->Close();
};


void splitString(std::vector< std::string >& output, const std::string& instr, const std::string& delimiter)
{
  auto input = instr;
  // remove spaces from input string
  input.erase(std::remove(input.begin(), input.end(), ' '), input.end());
  // proceed to parse input string
  while(input!="") {
    std::string d = input.substr(0, input.find(delimiter));
    output.push_back(d);
    if(input.find(delimiter)!=std::string::npos){ input.erase(0, input.find(delimiter) + delimiter.length()); }
    else { input = ""; }
  }
};


void stringReplace(std::string& txt, const std::string& from, const std::string& to)
{
  std::string::size_type n = 0;
  while ( ( n = txt.find( from, n ) ) != std::string::npos ) {
    txt.replace( n, from.size(), to );
    n += to.size();
  }
};


std::vector<std::string> parseColStr(const std::string& colStr)
{
  std::vector<std::string> vecStr;
  if (colStr=="") return vecStr;
  const auto& cStr = (colStr.rfind("_")!=std::string::npos ? colStr.substr(colStr.rfind("_")+1) : colStr);
  const auto& posN = cStr.find(*std::find_if(cStr.begin(), cStr.end(), [](unsigned char c){ return std::isdigit(c); }));
  vecStr.push_back(cStr.substr(0, posN));
  const auto& posY = cStr.rfind("Y");
  if (posN!=std::string::npos) vecStr.push_back(cStr.substr(posN, posY-posN));
  if (posY!=std::string::npos) vecStr.push_back(cStr.substr(posY+1));
  return vecStr;
};


std::string findLabel(const std::string& par, const std::string& var, const std::string& obj,
		      const std::string& chg, const std::string& col, const std::string& cha, const GlobalInfo& info)
{
  std::string tryLabel="";
  const StringVector_t tryChannel = { cha , "" };
  StringVector_t trySystem  = { col , "" };
  if (col.rfind("pPb",0)==0 || col.rfind("Pbp",0)==0) { auto t = col; stringReplace(t, col, "PA"); trySystem.push_back(t); }
  const StringVector_t tryCharge  = { chg , "" };
  const StringVector_t tryVariable = { var , "" };
  for (const auto& tryCha : tryChannel) {
    for (const auto& tryChg : tryCharge) {
      for (const auto& tryCol : trySystem) {
	for (const auto& tryVar : tryVariable) {
	  tryLabel = tryVar +"_"+ obj + tryCha + tryChg + (tryCol!="" ? "_"+tryCol : "");
	  if (contain(info.Par, par+tryLabel)) { return tryLabel; }
	}
      }
    }
  }
  return "";
};


void getLumiLabels(StringVector_t& labels, const std::string& PD, const std::string& col, const bool& isMC)
{
  std::string lumiLabel="";
  const auto& colV = parseColStr(col);
  if (!colV.empty()) {
    auto colN = colV[0];
    if (colN=="PP") { colN = "pp"; } else if (colN=="PA") { colN = "pPb"; }
    lumiLabel += colN;
  }
  if (isMC) { lumiLabel += " Simulation"; }
  else if (PD=="") { lumiLabel += " Data"; }
  else if (col=="PbPb5Y18") { lumiLabel += Form(" %.0f #mub^{-1}", PbPb::R5TeV::Y2018::LumiFromPD(PD)); }
  else if (col=="PP13Y18" ) { lumiLabel += Form(" %.1f pb^{-1}", pp::R13TeV::Y2018::LumiFromPD(PD)); }
  else if (col=="PP5Y17"  ) { lumiLabel += Form(" %.1f pb^{-1}", pp::R5TeV::Y2017::LumiFromPD(PD)); }
  else if (col.rfind("8Y16")!=std::string::npos) { lumiLabel += Form(" %.1f nb^{-1}", pPb::R8TeV::Y2016::LumiFromPD(PD, col)); }
  else if (col=="PbPb5Y15") { lumiLabel += Form(" %.0f #mub^{-1}", PbPb::R5TeV::Y2015::LumiFromPD(PD)); }
  labels.push_back(lumiLabel);
  //
  std::string lumiLabel2="";
  if (colV.size()>1) {
    auto colE = colV[1];
    if (colE=="5") { colE = "5.02 TeV"; } else if (colE=="8") { colE = "8.16 TeV"; } else if (colE=="13") { colE = "13 TeV"; } 
    colE = ((colV[0]=="PP") ? "#sqrt{s} = " : "#sqrt{s_{NN}} = ") +colE;
    lumiLabel2 += colE;
  }
  labels.push_back(lumiLabel2);
};


void roundValue(double& value , const uint& nDecimals)
{
  double tmp = value;
  tmp *= std::pow(10.0, nDecimals);
  tmp = std::round(tmp);
  tmp /= std::pow(10.0, nDecimals);
  value = tmp;
};


int getNBins(const double& min, const double& max, const double& binW)
{
  const auto& nBins = std::round((max - min)/binW);
  return int(nBins);
};


int getNBins(const std::string& var, const GlobalInfo& info)
{
  if (!contain(info.Var, var)) { return -1; }
  return getNBins(info.Var.at(var).at("Min"), info.Var.at(var).at("Max"), info.Var.at(var).at("binWidth"));
};


bool isEqual(const double& inVal1 , const double& inVal2 , const uint& nDecimals)
{
  double val1 = inVal1; roundValue(val1, nDecimals);
  double val2 = inVal2; roundValue(val2, nDecimals);
  if (val1==val2) return true;
  return false;
};


TH1D graphToHist(const TGraphAsymmErrors& gr , const int& addOnly = 0)
{
  auto tmp = gr;
  // Remove all unwanted points
  if (addOnly!=0) {
    for (int i = 0; i < tmp.GetN(); ){
      double eta, xSec;
      tmp.GetPoint(i, eta, xSec);
      if      (addOnly<0 && eta>=0) { tmp.RemovePoint(i); }
      else if (addOnly>0 && eta<=0) { tmp.RemovePoint(i); }
      else { i++; }
    }
  }
  const uint& nBin = tmp.GetN();
  double bins[nBin+1];
  for (uint iBin = 0; iBin <= nBin; iBin++) {
    double edge;
    if (iBin < nBin) {
      double x, y; tmp.GetPoint(iBin, x, y);
      const auto& ex = tmp.GetErrorXlow(iBin);
      edge = x - ex;
    }
    else {
      double x, y; tmp.GetPoint((nBin-1), x, y);
      const auto& ex = tmp.GetErrorXhigh(nBin-1);
      edge = x + ex;
    }
    bins[iBin] = edge;
  }
  TH1D h(tmp.GetName(), tmp.GetTitle(), nBin, bins);
  for (uint iBin = 0; iBin < nBin; iBin++) {
    double x, y; tmp.GetPoint(iBin, x, y);
    const auto& ey = tmp.GetErrorY(iBin);
    h.SetBinContent(iBin+1, y);
    h.SetBinError(iBin+1, ey);
  }
  tmp.TAttLine::Copy(h);
  tmp.TAttFill::Copy(h);
  tmp.TAttMarker::Copy(h);
  return h;
};


auto TIME_START_ = std::chrono::high_resolution_clock::now();
static inline void loadBar(const int& iEvent, const int& nEvents, const int& r = 100, const int& w = 100)
{
  // Only update r times.
  if ( iEvent == (nEvents-1) ) { std::cout << std::endl; }
  if ( (iEvent % ((nEvents/r) + 1)) != 0 ) return;
  // Calculuate the ratio of complete-to-incomplete.
  const float& ratio = (iEvent / (float)nEvents);
  const int&   c     = (ratio * w);
  // Get Time Difference
  const auto& TIME_END = std::chrono::high_resolution_clock::now();
  const int&  TIME_DIF = std::chrono::duration_cast<std::chrono::seconds>(TIME_END - TIME_START_).count();
  TIME_START_ = std::chrono::high_resolution_clock::now();
  // Show the percentage complete.
  const int& sec  = int( double(TIME_DIF) * (1.0-ratio) *100. );
  const int& min  = (sec / 60);
  const int& hour = (min / 60);
  printf("[INFO] %3d%% (%02d:%02d:%02d) [", int(ratio*100), hour, int(min%60), int(sec%60));
  // Show the load bar.
  for (int i = 0; i < c; i++) { std::cout << "="; }
  for (int i = c; i < w; i++) { std::cout << " "; }
  // ANSI Control codes to go back to the
  // previous line and clear it.
  std::cout << "]\r" << std::flush;
};


void printCPUInfo()
{
  CpuInfo_t cpuInfo;
  gSystem->GetCpuInfo(&cpuInfo, 100);
  MemInfo_t memInfo;
  gSystem->GetMemInfo(&memInfo);
  ProcInfo_t sysInfo;
  gSystem->GetProcInfo(&sysInfo);
  std::cout << "[INFO] CPU memory usage: " << std::endl;
  std::cout << "Total RAM: " << memInfo.fMemTotal << std::endl;
  std::cout << "Free RAM: " << (memInfo.fMemFree*100./memInfo.fMemTotal) << "%" << std::endl;
  std::cout << "Used RAM: " << (memInfo.fMemUsed*100./memInfo.fMemTotal) << "%" << std::endl;
  std::cout << "Total SWAP: " << memInfo.fSwapTotal << std::endl;
  std::cout << "Free SWAP: " << (memInfo.fSwapFree*100./memInfo.fSwapTotal) << "%" << std::endl;
  std::cout << "Used SWAP: " << (memInfo.fSwapUsed*100./memInfo.fSwapTotal) << "%" << std::endl;
  std::cout << "Total CPU: " << cpuInfo.fTotal << std::endl;
  std::cout << "System CPU: " << (cpuInfo.fSys*100./cpuInfo.fTotal) << "%" << std::endl;
  std::cout << "User CPU: " << (cpuInfo.fUser*100./cpuInfo.fTotal) << "%" << std::endl;
  std::cout << "Process System CPU: " << sysInfo.fCpuSys << std::endl;
  std::cout << "Process User CPU: " << sysInfo.fCpuUser << std::endl;
  std::cout << "Process Resident Memory: " << sysInfo.fMemResident << std::endl;
  std::cout << "Process Resident Virtual: " << sysInfo.fMemVirtual << std::endl;
};


#endif /* DATAUTILS_H_ */
