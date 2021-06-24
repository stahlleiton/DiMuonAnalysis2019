#ifndef plotResultsTree_C
#define plotResultsTree_C

// Auxiliary Files
#include "resultsUtils.h"
// c++ headers
#include <iostream>
#include <string>


void iniResultsGraph  ( GraphSeptaMap_t& graphMap , const BinCont_t& binMap   , const BinSextaMap_t& iVar , const bool& doSyst );
void fillResultsGraph ( GraphSeptaMap_t& graphMap , const BinCont_t& binMap   , const BinSextaMap_t& iVar );
void drawResultsGraph ( GraphSeptaMap_t& graphMap , const std::string& outDir , const bool& isMC );


void plotResultsTree(
		     GraphSeptaMap_t& graphMap,
		     const BinCont_t& binMap,
		     const BinSextaMap_t& var,
		     const std::string& outDir,
		     const bool& isMC,
		     const bool& doSyst
		     )
{
  //
  // Intialize the results graph
  iniResultsGraph(graphMap, binMap, var, doSyst);
  //
  // Fill the results graph
  fillResultsGraph(graphMap, binMap, var);
  //
  // Draw the results graph
  drawResultsGraph(graphMap, outDir, isMC);
  //
};


void iniResultsGraph(GraphSeptaMap_t& graphMap, const BinCont_t& binMap, const BinSextaMap_t& var, const bool& doSyst)
{
  //
  StringVector_t gType({"Err_Tot", "Err_Stat"});
  if (doSyst) { gType.push_back("Err_Syst"); }
  //
  for (const auto& o : binMap) {
    for (const auto& c : o.second) {
      for (const auto& v : c.second) {
	if (!contain(var.at(o.first).at(c.first).begin()->second, v.first)) continue;
	for (const auto& obs : v.second) {
	  for (const auto& incB : obs.second) {
	    for (const auto& subB : incB.second) {
	      for (const auto& t : gType) {
		auto& graph = graphMap[o.first][c.first][obs.first][incB.first][v.first][subB.first][t];
		graph.Set(subB.second.size());
		// Set Graph Name
		std::string label = "_";
		for (const auto& b : incB.first) { label += Form("%s_%.0f_%.0f_", b.name().c_str(), b.low()*100., b.high()*100.); }
		const auto& name = ("gr_"+o.first+"_"+c.first+"_"+v.first+"_"+obs.first+label+t);
		graph.SetName(name.c_str());
	      }
	    }
	  }
	}
      }
    }
  }
};


void fillResultsGraph(GraphSeptaMap_t& graphMap, const BinCont_t& binMap, const BinSextaMap_t& iVar)
{
  //
  std::cout << "[INFO] Filling the output graphs" << std::endl;
  //
  for (auto& o : graphMap) {
    for (auto& c : o.second) {
      for (auto& obs : c.second) {
	for (auto& incB : obs.second) {
	  for (auto& v : incB.second) {
	    for (auto& subB : v.second) {
	      const auto& bins = binMap.at(o.first).at(c.first).at(v.first).at(obs.first).at(incB.first).at(subB.first);
	      for (auto& gr : subB.second) {
		auto& graph = gr.second;
		for (size_t iBin=0; iBin<bins.size(); iBin++) {
		  const auto& b = *std::next(bins.begin(), iBin);
		  //
		  // Extract the parameters needed for each axis
		  //
		  // X Value
		  const bool useMean = (obs.first=="NTrack");
		  const auto& binX = b.getbin(obs.first);
		  const float& X = (useMean ? binX.mean() : (binX.high()+binX.low())/2.);
		  // X Error
		  const float& Err_X = (useMean ? binX.width() : (binX.high()-binX.low())/2.);
		  const auto& Err_X_High = std::min(X + Err_X, binX.high()) - X;
		  const auto& Err_X_Low  = X - std::max(X - Err_X, binX.low());
		  //
		  double norm = 1.0;
		  if (v.first=="Cross_Section") {
		    norm *= (binX.high()-binX.low());
		    if (binX.name()=="Centrality") { norm *= 0.01; }
		    if (norm==0.) { binX.print(); throw std::logic_error("[ERROR] fillResultsGraph: Cross section bin width is 0!"); }
		  }
		  //
		  for (const auto& pd : iVar.at(o.first).at(c.first)) {
		    if (!contain(pd.second, v.first)) continue;
		    const auto& var = pd.second.at(v.first);
		    if (!contain(var.at("Val"), b)) continue;
		    //
		    // Y Value
		    const auto& Y = var.at("Val").at(b)/norm;
		    //
		    // Compute total systematic error
		    const auto& Err_Y_Syst_High = var.at("Err_Syst_High").at(b)/norm;
		    const auto& Err_Y_Syst_Low  = var.at("Err_Syst_Low" ).at(b)/norm;
		    //
		    // Compute total statistic error
		    const auto& Err_Y_Stat_High = var.at("Err_Stat_High").at(b)/norm;
		    const auto& Err_Y_Stat_Low  = var.at("Err_Stat_Low" ).at(b)/norm;
		    //
		    // Y Error
		    double Err_Y_High = 0.0 , Err_Y_Low  = 0.0;
		    if (gr.first=="Err_Tot") {
		      Err_Y_High = sumErrors({Err_Y_Stat_High , Err_Y_Syst_High});
		      Err_Y_Low  = sumErrors({Err_Y_Stat_Low  , Err_Y_Syst_Low });
		    }
		    else if (gr.first=="Err_Stat") {
		      Err_Y_High = Err_Y_Stat_High;
		      Err_Y_Low  = Err_Y_Stat_Low;
		    }
		    else if (gr.first=="Err_Syst") {
		      Err_Y_High = Err_Y_Syst_High;
		      Err_Y_Low  = Err_Y_Syst_Low;
		    }
		    //
		    // Check results
		    //
		    if (Err_X_Low<=0. || Err_X_High<=0. || isnan(Err_X_Low) || isnan(Err_X_High) || isnan(X)) {
		      std::cout << o.first << " " << c.first << " " << obs.first << " " << v.first; subB.first.first.print();
		      throw std::logic_error(Form("[ERROR] Invalid x value: %.4f - %.4f + %.4f", X, Err_X_Low, Err_X_High));
		    }
		    if (Err_Y_Low<=0. || Err_Y_High<=0. || isnan(Err_Y_Low) || isnan(Err_Y_High) || isnan(Y)) {
		      std::cout << o.first << " " << c.first << " " << obs.first << " " << v.first; subB.first.first.print();
		      throw std::logic_error(Form("[ERROR] Invalid y value: %.4f - %.4f + %.4f", Y, Err_Y_Low, Err_Y_High));
		    }
		    //
		    // Fill the graph
		    //
		    graph.SetPoint(iBin, X, Y);
		    graph.SetPointError(iBin, Err_X_Low, Err_X_High, Err_Y_Low, Err_Y_High);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
};


void drawResultsGraph(GraphSeptaMap_t& graphMap, const std::string& outDir, const bool& isMC)
{
  // Set Style
  setStyle();
  std::cout << "[INFO] Drawing the output graphs" << std::endl;
  //
  // Draw all graphs
  for (const auto& o : graphMap) {
    for (const auto& c : o.second) {
      for (const auto& ob : c.second) {
	for (const auto& incB : ob.second) {
	  for (const auto& v : incB.second) {
	    const std::string& obj = o.first;
	    const std::string& col = c.first;
	    const std::string& var = v.first;
	    const std::string& obs = ob.first;
	    // Create Canvas
	    TCanvas c("c", "c", 1000, 1000); c.cd();
	    // Create the Text Info
	    TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035); float dy = 0;
	    std::vector< std::string > textToPrint;
	    textToPrint.push_back(formatDecayLabel(var, obj));
	    for (const auto& b : incB.first) { textToPrint.push_back(formatObsRange(b)); }
	    // Declare the graph vector (for drawing with markers)
	    std::map< BinPair_t , std::vector<TGraphAsymmErrors> > grDiMap;
	    // Format and add graphs
	    double xMin = 99999999., xMax = -99999999., yMin = 99999999., yMax = -99999999.;
	    for (const auto& subB : v.second) {
	      auto& graph = subB.second;
	      std::vector<TGraphAsymmErrors> grV;
	      grV.push_back(graph.at("Err_Stat"));
	      if (contain(graph, "Err_Syst")) {
		grV.push_back(graph.at("Err_Syst"));
		double xRng = 99999999.; for (int j=0; j<grV.back().GetN(); j++) { xRng = std::min(xRng, std::min(grV.back().GetEXlow()[j], grV.back().GetEXhigh()[j])); }
		for (int j=0; j<grV.back().GetN(); j++) { grV.back().GetEXlow()[j] = 0.4*xRng; grV.back().GetEXhigh()[j] = 0.4*xRng; }
		grV.push_back(graph.at("Err_Tot"));
		grV.push_back(graph.at("Err_Tot"));
		for (int i=0; i<grV[0].GetN(); i++) { grV[2].GetY()[i] += grV[2].GetErrorYhigh(i); grV[3].GetY()[i] -= grV[3].GetErrorYlow(i); }
		for (size_t j=2; j<=3; j++) { for (int i=0; i<grV[j].GetN(); i++) { grV[j].SetPointError(i, 0.5*xRng, 0.5*xRng, 0.0, 0.0); } }
	      }
	      grDiMap[subB.first] = grV;
	      //
	      for (int i=0; i<grV[0].GetN(); i++) {
		xMin = std::min(grV[0].GetX()[i]-grV[0].GetErrorX(i), xMin);  yMin = std::min(grV[0].GetY()[i]-grV[0].GetErrorY(i), yMin);
		xMax = std::max(grV[0].GetX()[i]+grV[0].GetErrorX(i), xMax);  yMax = std::max(grV[0].GetY()[i]+grV[0].GetErrorY(i), yMax);
	      }
	    }
	    int iGr = 0;
	    for (auto& grM : grDiMap) {
	      auto& grV = grM.second;
	      for (auto& gr : grV) {
		formatResultsGraph(gr, var, obs, obj, col, COLOR[iGr], false, {xMin, yMin, xMax, yMax});
		gr.SetFillStyle(1001);
	      }
	      iGr += 1;
	      //for (int i=0; i<grV[0].GetN(); i++) { grV[0].SetPointEXhigh(i, 0.0); grV[0].SetPointEXlow(i, 0.0); }
	      for (size_t j=2; j<=3; j++) { grV[j].SetMarkerSize(0); }
	      if (grDiMap.size()>1) { grV[1].SetFillColor(kGreen+2); }
	    }
	    // Draw the graphs
	    grDiMap.begin()->second[0].Draw("ap");
	    for (auto& grM : grDiMap) {
	      //if (grM.second.size()>1) { grM.second[1].Draw("same2"); }
	      if (grM.second.size()>3) { grM.second[2].Draw("samep"); grM.second[3].Draw("samep"); }
	      grM.second[0].Draw("samep");
	    }
	    const auto& graph = grDiMap.begin()->second[0];
	    // Draw the Line
	    TLine line(graph.GetXaxis()->GetXmin(), 1.0, graph.GetXaxis()->GetXmax(), 1.0); line.SetLineStyle(2);
	    if (var=="ForwardBackward_Ratio" || var=="RatioTo1S") { line.Draw("same"); }
	    // Initialize the legend
	    std::unique_ptr<TLegend> leg;
	    if (grDiMap.size()>1) {
	      // Define the position and size of the legend
	      double xmin = 0.61 , xmax = 0.78 , ymin = 0.52 , ymax = 0.76 , legSize = 0.034;
	      ymin = std::max(ymax-grDiMap.size()*legSize*1.25, ymin);
	      leg.reset(new TLegend(xmin, ymin, xmax, ymax));
	      for (const auto& grM : grDiMap) {
		// Add legend entry
		std::string legLbl = ""; for (const auto& sB : grM.first.first) { legLbl += formatObsRange(sB)+" , "; };
		if (legLbl.find(" , ")!=std::string::npos) { legLbl.erase(legLbl.find(" , "), 3); }
		formatLegendEntry(*leg->AddEntry(&grM.second[0], legLbl.c_str(), "pe"), legSize);
	      }
	      // Draw the Legend
	      leg->Draw("same");
	    }
	    // Update
	    c.Modified(); c.Update();
	    // Draw the text
	    tex.SetTextSize(0.055); tex.DrawLatex(0.22, 0.84, textToPrint[0].c_str());
	    tex.SetTextSize(0.060); tex.SetTextFont(61); tex.DrawLatex(0.78, 0.84, "CMS"); tex.SetTextFont(62);
	    tex.SetTextSize(0.046); tex.SetTextFont(52); tex.DrawLatex(0.69, 0.79, "Preliminary"); tex.SetTextFont(62);
	    for (size_t i=1; i<textToPrint.size(); i++) { tex.SetTextSize(0.045); tex.DrawLatex(0.22, 0.76-dy, textToPrint[i].c_str()); dy+=0.060; }
	    if (var=="Cross_Section") { tex.SetTextSize(0.030); tex.DrawLatex(0.25, 0.17, "Lumi. uncertainty not shown"); }
	    // Update
	    c.Modified(); c.Update(); // Pure paranoia
	    //
	    // set the CMS style
	    StringVector_t lumiLabels; getLumiLabels(lumiLabels, "DIMUON", col, isMC);
	    CMS_lumi(&c, 33, (" "+lumiLabels[0]), lumiLabels[1], false, 0.65, false);
	    // Update
	    c.Modified(); c.Update(); // Pure paranoia
	    //
	    // Create Output Directory
	    const std::string& plotDir = outDir+"/Plots/"+obs;
	    makeDir(plotDir + "/png/");
	    makeDir(plotDir + "/pdf/");
	    makeDir(plotDir + "/root/");
	    makeDir(plotDir + "/C/");
	    //
	    // Save Canvas
	    const std::string& grName = graph.GetName();
	    const auto& name = grName.substr(0, grName.rfind("_Err_Stat"));
	    c.SaveAs((plotDir + "/png/"  + name + ".png" ).c_str());
	    c.SaveAs((plotDir + "/pdf/"  + name + ".pdf" ).c_str());
	    c.SaveAs((plotDir + "/root/" + name + ".root").c_str());
	    //c.SaveAs((plotDir + "/C/"    + name + ".C"   ).c_str());
	    //
	    // Clean up memory
	    c.Clear(); c.Close();
	  }
	}
      }
    }
  }
};


#endif // #ifndef plotResultsTree_C
