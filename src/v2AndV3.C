#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

#include "Math/ProbFuncMathCore.h"

#include "include/mntToXRootdFileString.h"
#include "include/inToOutFileString.h"
#include "include/histDefUtility.h"
#include "include/kirchnerPalette.h"

int v2AndV3(const std::string inFileName, const Bool_t isMC = false)
{
  kirchnerPalette col;
  //points digitized from CMS-HIN-16-019, Fig.2 here
  //http://cms-results.web.cern.ch/cms-results/public-results/publications/HIN-16-019/index.html
  const Int_t nCent16019 = 11;
  /*  
  const Float_t centXVals16019[nCent16019] = {7.460595446584943, 
					    12.451838879159371,
					    17.443082311733797,
					    22.521891418563932,
					    27.425569176882668,
					    32.50437828371278, 
					    37.49562171628721, 
					    42.48686514886165, 
					    47.47810858143607, 
					    52.469352014010504,
					    57.46059544658491};
  */

					   
  const Float_t centYVals16019[nCent16019] = {0.04861047835990889,
					    0.06391799544419134,
					    0.07649202733485194,
					    0.08724373576309796,
					    0.09562642369020502,
					    0.10182232346241457,
					    0.10656036446469247,
					    0.10929384965831435,
					    0.11056947608200454,
					    0.11020501138952163,
					    0.10801822323462414};
  
  std::vector<std::string> fileList;
  if(inFileName.find(".root") != std::string::npos) fileList.push_back(inFileName);
  else if(inFileName.find(".txt") != std::string::npos){
    std::ifstream file(inFileName.c_str());
    std::string tempStr;
    while(std::getline(file, tempStr)){
      if(tempStr.size() == 0) continue;
      while(tempStr.substr(0,1).find(" ") != std::string::npos){tempStr.replace(0,1,"");}
      if(tempStr.find(".root") == std::string::npos) continue;
      fileList.push_back(tempStr);
    }
    file.close();
  }
  else{
    std::cout << "input \'" << inFileName << "\' produces no valid root files return 1" << std::endl;
    return 1;
  }

  if(fileList.size() == 0){
    std::cout << "input \'" << inFileName << "\' produces no valid root files return 1" << std::endl;
    return 1;
  }

  std::string mcStr = "_DATA";
  if(isMC) mcStr = "_MC";
  const std::string outFileName = "output/" + inToOutFileString(inFileName, "V2AndV3" + mcStr);
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");

  const Int_t nCentBins = 11;
  const Int_t centBinsLow[nCentBins] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55};
  const Int_t centBinsHi[nCentBins] = {10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60};
  const Float_t centBins[nCentBins+1] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60};

  TH1F* v2_h[nCentBins];
  TH1F* v3_h[nCentBins];

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
    std::string centStr2 = std::to_string(centBinsLow[cI]) + "-" + std::to_string(centBinsHi[cI]) + "%";

    std::string v2Name = "v2_" + centStr + "_h";
    std::string v3Name = "v3_" + centStr + "_h";

    v2_h[cI] = new TH1F(v2Name.c_str(), (";v_{2};Counts (" + centStr2 + ")").c_str(), 100, 0, 1.);
    v3_h[cI] = new TH1F(v3Name.c_str(), (";v_{3};Counts (" + centStr2 + ")").c_str(), 100, 0, 1.);

    centerTitles({v2_h[cI], v3_h[cI]});
  }

  TH1F* v2_Mean_h = new TH1F("v2_Mean_h", ";Centrality (%);#LTv_{2}#GT", nCentBins, centBins);
  TH1F* v3_Mean_h = new TH1F("v3_Mean_h", ";Centrality (%);#LTv_{3}#GT", nCentBins, centBins);

  TH1F* v2_HIN16019_h = new TH1F("v2_HIN16019_h", ";#LTv_{2}#GT;Centrality (%)", nCentBins, centBins);

  centerTitles({v2_Mean_h, v3_Mean_h, v2_HIN16019_h});

  Int_t hiBin_;

  std::vector<double>* rhoFlowFitParams_p=NULL;

  const Int_t nMaxEntries = 100000;
  Int_t totEntries = 0;

  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    TFile* inFile_p = TFile::Open(mntToXRootdFileString(fileList.at(fI)).c_str(), "READ");
    TTree* rhoTree_p = (TTree*)inFile_p->Get("hiPuRhoR3Analyzer/t");
    TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");

    rhoTree_p->SetBranchStatus("*", 0);
    rhoTree_p->SetBranchStatus("rhoFlowFitParams", 1);

    rhoTree_p->SetBranchAddress("rhoFlowFitParams", &rhoFlowFitParams_p);

    hiTree_p->SetBranchStatus("*", 0);
    hiTree_p->SetBranchStatus("hiBin", 1);

    hiTree_p->SetBranchAddress("hiBin", &hiBin_);

    const Int_t nEntries = rhoTree_p->GetEntries();

    for(Int_t entry = 0; entry < nEntries; ++entry){
      if(totEntries%10000 == 0) std::cout << " Entry: " << totEntries << "/" << nMaxEntries << std::endl;
      rhoTree_p->GetEntry(entry);
      hiTree_p->GetEntry(entry);

      Int_t centPos = -1;
      for(Int_t cI = 0; cI < nCentBins; ++cI){
	if(hiBin_/2 >= centBinsLow[cI] && hiBin_/2 < centBinsHi[cI]){
	  centPos = cI;
	  break;
	}
      }

      if(centPos < 0) continue;

      if(rhoFlowFitParams_p->size() != 0){
	double val = ROOT::Math::chisquared_cdf_c(rhoFlowFitParams_p->at(5), rhoFlowFitParams_p->at(6));
        bool minProb = val > .05;
        bool maxProb = val < .95;

	if(minProb && maxProb){
	  v2_h[centPos]->Fill(rhoFlowFitParams_p->at(1));
	  v3_h[centPos]->Fill(rhoFlowFitParams_p->at(3));
	}
      }

      ++totEntries;
      if(totEntries >= nMaxEntries) break;
    }

    inFile_p->Close();
    delete inFile_p;

    if(totEntries >= nMaxEntries) break;
  }

  outFile_p->cd();
  
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    v2_h[cI]->Write("", TObject::kOverwrite);
    v3_h[cI]->Write("", TObject::kOverwrite);

    v2_Mean_h->SetBinContent(cI+1, v2_h[cI]->GetMean());
    v2_Mean_h->SetBinError(cI+1, v2_h[cI]->GetMeanError());

    v3_Mean_h->SetBinContent(cI+1, v3_h[cI]->GetMean());
    v3_Mean_h->SetBinError(cI+1, v3_h[cI]->GetMeanError());
    
    v2_HIN16019_h->SetBinContent(cI+1, centYVals16019[cI]);
    v2_HIN16019_h->SetBinError(cI+1, 0.);

    delete v2_h[cI];
    delete v3_h[cI];
  }

  v2_Mean_h->Write("", TObject::kOverwrite);
  v3_Mean_h->Write("", TObject::kOverwrite);

  v2_HIN16019_h->Write("", TObject::kOverwrite);

  TCanvas* canvV2_p = new TCanvas("canvV2_p", "canvV2_p", 400, 400);
  canvV2_p->SetTopMargin(0.01);
  canvV2_p->SetRightMargin(0.01);
  canvV2_p->SetLeftMargin(canvV2_p->GetLeftMargin()*1.3);
  canvV2_p->SetBottomMargin(canvV2_p->GetLeftMargin());

  v2_Mean_h->SetMaximum(.2);
  v2_Mean_h->SetMinimum(.0);

  v2_Mean_h->SetMarkerStyle(20);
  v2_Mean_h->SetMarkerSize(0.8);
  v2_Mean_h->SetMarkerColor(col.getColor(0));
  v2_Mean_h->SetLineColor(col.getColor(0));

  gStyle->SetOptStat(0);
  v2_Mean_h->DrawCopy("HIST E1 P");

  v2_HIN16019_h->SetMarkerStyle(20);
  v2_HIN16019_h->SetMarkerSize(0.8);
  v2_HIN16019_h->SetMarkerColor(col.getColor(2));
  v2_HIN16019_h->SetLineColor(col.getColor(2));

  v2_HIN16019_h->DrawCopy("HIST E1 P SAME");

  TLegend* leg_p = new TLegend(0.2, .7, 0.4, .95);
  leg_p->SetFillStyle(0);
  leg_p->SetFillColor(0);
  leg_p->SetBorderSize(0);
  leg_p->SetTextFont(43);
  leg_p->SetTextSize(14);

  leg_p->AddEntry(v2_Mean_h, "Jet Extraction", "P L");
  leg_p->AddEntry(v2_HIN16019_h, "HIN-16-019 (UNFOLDED)", "P L");

  leg_p->Draw("SAME");

  TDatime* date = new TDatime();
  std::string saveStr = "pdfDir/v2Overlay" + mcStr + "_" + std::to_string(date->GetDate()) + ".pdf";
  delete date;
  canvV2_p->SaveAs(saveStr.c_str());

  delete canvV2_p;

  delete leg_p;

  delete v2_Mean_h;
  delete v3_Mean_h;

  delete v2_HIN16019_h;

  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2 && argc != 3){
    std::cout << "Usage: ./v2AndV3.exe <inFileName> <isMC>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 2) retVal += v2AndV3(argv[1]);
  else if(argc == 3) retVal += v2AndV3(argv[1], std::stoi(argv[2]));
  return retVal;
}
