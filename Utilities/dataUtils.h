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
#include <string>
#include <sys/stat.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <memory>
#include <utility>
#include <algorithm>


// Utiliy Functions


template<class C, class T>
inline auto contain_impl(const C& c, const T& x, int)
-> decltype(c.find(x), true)
{ return end(c) != c.find(x); }

template<class C, class T>
inline bool contain_impl(const C& v, const T& x, long)
{ return end(v) != std::find(begin(v), end(v), x); }

template<class C, class T>
auto contain(const C& c, const T& x)
-> decltype(end(c), true)
{ return contain_impl(c, x, 0); }


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


void findSubDir(std::vector<std::string>& dirList, const std::string& dirName)
{
  TSystemDirectory dir(dirName.c_str(), dirName.c_str());
  auto subdirs = std::unique_ptr<TList>(dir.GetListOfFiles());
  if (subdirs) {
    TSystemFile *subdir;
    TIter next(subdirs.get());
    while ( (subdir=dynamic_cast<TSystemFile*>(next())) ) {
      if (subdir->IsDirectory() && std::string(subdir->GetName())!="." && std::string(subdir->GetName())!="..") {
        dirList.push_back(dirName + subdir->GetName() + "/");
        std::cout << "[INFO] Input subdirectory: " << dirName + subdir->GetName() + "/" << " found!" << std::endl;
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


std::vector<std::string>
parseColStr(const std::string& colStr)
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
  return std::min(int(nBins), 2000);
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
