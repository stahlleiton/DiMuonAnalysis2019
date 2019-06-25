#ifndef rooModelUtils_h
#define rooModelUtils_h

#include "TH1.h"
#include "TInterpreter.h"

#include "RooFit.h"
#include "RooWorkspace.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooHistPdf.h"

#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <set>

#include "initClasses.h"


std::set<std::string> MODELCLASS_;


bool importModelClass(RooWorkspace& ws, const std::string& className)
{
  const std::string& CWD = getcwd(NULL, 0);
  const auto& classDir = CWD+"/Macros/Utilities/Models/";
  TInterpreter::EErrorCode ecode;
  if (MODELCLASS_.find(className)==MODELCLASS_.end()) {
    gInterpreter->ProcessLineSynch(Form(".L %s%s.cxx+",classDir.c_str(), className.c_str()), &ecode);
    if (ecode!=TInterpreter::kNoError) return false;
    MODELCLASS_.insert(className);
  }
  auto classPdf = std::unique_ptr<RooAbsPdf>((RooAbsPdf*)gInterpreter->ProcessLineSynch(Form("new %s()", className.c_str()), &ecode));
  if (ecode!=TInterpreter::kNoError) return false;
  return ws.importClassCode(classPdf->IsA());
}


TH1* rebinhist(const TH1& hist, const double& xmin, const double& xmax, const std::string& type="Old")
{
  auto hcopy = std::unique_ptr<TH1>(dynamic_cast<TH1*>(hist.Clone("hcopy")));
  // range of the new hist
  int imin = hcopy->FindBin(xmin);
  if (imin>=hcopy->GetNbinsX()) imin=1;
  int imax = hcopy->FindBin(0.999999*xmax);
  if (imax<=1) imax=hcopy->GetNbinsX();
  std::vector<double> newbins;
  newbins.push_back(hcopy->GetBinLowEdge(imin));
  for (int i=imin; i<=imax; i++) {
    if (hcopy->GetBinContent(i)>0.0) {
      newbins.push_back(hcopy->GetBinLowEdge(i)+hcopy->GetBinWidth(i));
    } else {
      int nrebin=2;
      for (i++; i<=imax; i++) {
        if (hcopy->GetBinContent(i)>0.0) {
          newbins.push_back(hcopy->GetBinLowEdge(i)+hcopy->GetBinWidth(i));
          const double& newval = (hcopy->GetBinContent(i)/nrebin);
          hcopy->SetBinContent(i, newval);
          if (type=="New") { for (int j=1; j<nrebin; j++) { hcopy->SetBinContent(i-j, newval); } }
          break;
        }
        nrebin++;
      }
    }
  }
  if (type=="Old") {
    if (xmin < newbins[1]) newbins[0] = xmin;
    if (xmax > newbins[newbins.size()-1]) { newbins.push_back(xmax); }
    return hcopy->Rebin(newbins.size()-1, "hnew", newbins.data());
  }
  if (type=="New") {
    return dynamic_cast<TH1*>(hcopy->Clone("hnew"));
  }
  return NULL;
};


bool histToPdf(RooWorkspace& ws, const std::string& pdfName, const std::string& dsName, const std::string& var, const std::vector< double >& range)
{
  //
  if (ws.pdf(pdfName.c_str())) { std::cout << Form("[INFO] The %s Template has already been created!", pdfName.c_str()) << std::endl; return true; }
  std::cout << Form("[INFO] Implementing %s Template", pdfName.c_str()) << std::endl;
  //
  if (ws.data(dsName.c_str())==NULL) { std::cout << "[WARNING] DataSet " << dsName << " was not found!" << std::endl; return false; }
  if (ws.data(dsName.c_str())->numEntries()<=2.0) { std::cout << "[WARNING] DataSet " << dsName << " has too few events!" << std::endl; return false; }
  if (ws.var(var.c_str())==NULL) { std::cout << "[WARNING] Variable " << var << " was not found!" << std::endl; return false; }
  // Create the histogram
  auto histName = pdfName;
  histName.replace(histName.find("pdf"), std::string("pdf").length(), "h");
  std::unique_ptr<TH1D> hist = std::unique_ptr<TH1D>(dynamic_cast<TH1D*>(ws.data(dsName.c_str())->createHistogram(histName.c_str(), *ws.var(var.c_str()), RooFit::Binning(int(range[0]), range[1], range[2]))));
  if (hist==NULL) { std::cout << "[WARNING] Histogram " << histName << " is NULL!" << std::endl; return false; }
  // Cleaning the input histogram
  // 1) Remove the Under and Overflow bins
  hist->ClearUnderflowAndOverflow();
  // 2) Set negative bin content to zero
  for (int i=0; i<=hist->GetNbinsX(); i++) { if (hist->GetBinContent(i)<0.0) { hist->SetBinContent(i, 0.0); } }
  // 2) Reduce the range of histogram and rebin it
  hist.reset(dynamic_cast<TH1D*>(rebinhist(*hist, range[1], range[2])));
  if (hist==NULL) { std::cout << "[WARNING] Cleaned Histogram of " << histName << " is NULL!" << std::endl; return false; }
  auto dataName = pdfName;
  dataName.replace(dataName.find("pdf"), std::string("pdf").length(), "dh");
  std::unique_ptr<RooDataHist> dataHist = std::unique_ptr<RooDataHist>(new RooDataHist(dataName.c_str(), "", *ws.var(var.c_str()), hist.get()));
  if (dataHist==NULL) { std::cout << "[WARNING] DataHist used to create " << pdfName << " failed!" << std::endl; return false; } 
  if (dataHist->sumEntries()==0) { std::cout << "[WARNING] DataHist used to create " << pdfName << " is empty!" << std::endl; return false; } 
  if (std::abs(dataHist->sumEntries() - hist->GetSumOfWeights())>0.001) {
    std::cout << "[ERROR] DataHist (" << dataHist->sumEntries() << ")  used to create histogram (" << hist->GetSumOfWeights() << ") for PDF " << pdfName << "  " << " is invalid!  " << std::endl; return false;
  }
  ws.import(*dataHist);
  ws.var(var.c_str())->setBins(int(range[0])); // Bug Fix
  std::unique_ptr<RooHistPdf> pdf = std::unique_ptr<RooHistPdf>(new RooHistPdf(pdfName.c_str(), pdfName.c_str(), *ws.var(var.c_str()), *dynamic_cast<RooDataHist*>(ws.data(dataName.c_str()))));
  //std::unique_ptr<RooKeysPdf> pdf = std::unique_ptr<RooKeysPdf>(new RooKeysPdf(pdfName.c_str(), pdfName.c_str(), *ws.var(var.c_str()), *dynamic_cast<RooDataSet*>(ws.data(dsName.c_str())), RooKeysPdf::NoMirror, 0.4));
  if (pdf==NULL) { std::cout << "[WARNING] RooHistPDF " << pdfName << " is NULL!" << std::endl; return false; }
  ws.import(*pdf);
  return true;
};


#endif // #ifndef rooModelUtils_h
