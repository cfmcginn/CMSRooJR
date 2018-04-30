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
#include "TF1.h"
#include "TDatime.h"
#include "TObjArray.h"
#include "TSystem.h"

#include "Math/ProbFuncMathCore.h"

#include "include/mntToXRootdFileString.h"
#include "include/inToOutFileString.h"
#include "include/histDefUtility.h"
#include "include/kirchnerPalette.h"
#include "include/goodGlobalSelection.h"

#include "CustomCanvas.h"
#include "TexSlides.C"


int v2AndV3(const std::string inFileName, const Bool_t isMC = false)
{
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());

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

  const Int_t nFlowTypes = 6;
  const std::string flowTypes[nFlowTypes] = {"", "HFPlane", "FreePlane", "JettyExclude", "JettyExcludeHFPlane", "JettyExcludeFreePlane"};
  const std::string flowTypes2[nFlowTypes] = {"Default", "HFPlane", "FreePlane", "JettyExcludeDefault", "JettyExcludeHFPlane", "JettyExcludeFreePlane"};

  const Int_t nCentBins = 11;
  const Int_t centBinsLow[nCentBins] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55};
  const Int_t centBinsHi[nCentBins] = {10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60};
  const Float_t centBins[nCentBins+1] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60};

  TH1F* v2_h[nCentBins][nFlowTypes];
  TH1F* v3_h[nCentBins][nFlowTypes];

  for(Int_t fI = 0; fI < nFlowTypes; ++fI){
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      const std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
      const std::string centStr2 = std::to_string(centBinsLow[cI]) + "-" + std::to_string(centBinsHi[cI]) + "%";
      
      const std::string v2Name = "v2_" + flowTypes2[fI] + "_" + centStr + "_h";
      const std::string v3Name = "v3_" + flowTypes2[fI] + "_" + centStr + "_h";
      
      //      std::cout << "Cent: " << fI << "/" << nFlowTypes << ", " << cI << "/" << nCentBins << ", " << centStr << std::endl;

      v2_h[cI][fI] = new TH1F(v2Name.c_str(), (";v_{2};Counts (" + centStr2 + ")").c_str(), 150, 0, 0.6);
      v3_h[cI][fI] = new TH1F(v3Name.c_str(), (";v_{3};Counts (" + centStr2 + ")").c_str(), 150, 0, 0.6);
      
      centerTitles({v2_h[cI][fI], v3_h[cI][fI]});
    }
  }

  TH1F* v2_Mean_h = new TH1F("v2_Mean_h", ";Centrality (%);#LTv_{2}#GT", nCentBins, centBins);
  TH1F* v2_Raw_h = new TH1F("v2_Raw_h", ";Centrality (%);#LTv_{2}#GT", nCentBins, centBins);
  TH1F* v3_Mean_h = new TH1F("v3_Mean_h", ";Centrality (%);#LTv_{3}#GT", nCentBins, centBins);

  TH1F* v2_HIN16019_h = new TH1F("v2_HIN16019_h", ";#LTv_{2}#GT;Centrality (%)", nCentBins, centBins);

  centerTitles({v2_Mean_h, v2_Raw_h, v3_Mean_h, v2_HIN16019_h});

  Float_t hiHF_;
  Float_t vz_;
  Int_t hiBin_;

  Int_t pprimaryVertexFilter_;
  Int_t HBHENoiseFilterResultRun2Loose_;
  Int_t pclusterCompatibilityFilter_;
  Int_t phfCoincFilter3_;

  goodGlobalSelection sel;
  sel.setIsPbPb(true);

  std::vector<std::vector<double>*> rhoFlowFitParams_p;
  rhoFlowFitParams_p.reserve(nFlowTypes);
  for(Int_t fI = 0; fI < nFlowTypes; ++fI){
    rhoFlowFitParams_p.push_back(NULL);
  }

  const Int_t nMaxEntries = 1000000;
  Int_t totEntries = 0;

  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    std::string tempFileName = mntToXRootdFileString(fileList.at(fI));
    
    std::cout << "File: " << tempFileName << std::endl;
    TFile* inFile_p = TFile::Open(tempFileName.c_str(), "READ");
    TTree* rhoTree_p[nFlowTypes];
    for(Int_t fI = 0; fI < nFlowTypes; ++fI){
      rhoTree_p[fI] = (TTree*)inFile_p->Get(("hiFlowAnalyzer" + flowTypes[fI] + "/t").c_str());
      rhoTree_p[fI]->SetBranchStatus("*", 0);
      rhoTree_p[fI]->SetBranchStatus("rhoFlowFitParams", 1);
      
      rhoTree_p[fI]->SetBranchAddress("rhoFlowFitParams", &(rhoFlowFitParams_p.at(fI)));
    }

    TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
    TTree* skimTree_p = (TTree*)inFile_p->Get("skimanalysis/HltTree");

    hiTree_p->SetBranchStatus("*", 0);
    hiTree_p->SetBranchStatus("vz", 1);
    hiTree_p->SetBranchStatus("hiHF", 1);
    hiTree_p->SetBranchStatus("hiBin", 1);

    hiTree_p->SetBranchAddress("vz", &vz_);
    hiTree_p->SetBranchAddress("hiHF", &hiHF_);
    hiTree_p->SetBranchAddress("hiBin", &hiBin_);

    skimTree_p->SetBranchStatus("*", 0);
    skimTree_p->SetBranchStatus("pprimaryVertexFilter", 1);
    skimTree_p->SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
    skimTree_p->SetBranchStatus("phfCoincFilter3", 1);
    skimTree_p->SetBranchStatus("pclusterCompatibilityFilter", 1);

    skimTree_p->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter_);
    skimTree_p->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose_);
    skimTree_p->SetBranchAddress("phfCoincFilter3", &phfCoincFilter3_);
    skimTree_p->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter_);

    const Int_t nEntries = rhoTree_p[0]->GetEntries();

    for(Int_t entry = 0; entry < nEntries; ++entry){
      if(totEntries%10000 == 0) std::cout << " Entry: " << totEntries << "/" << nMaxEntries << std::endl;

      for(Int_t fI = 0; fI < nFlowTypes; ++fI){
	rhoTree_p[fI]->GetEntry(entry);
      }
      
      hiTree_p->GetEntry(entry);
      skimTree_p->GetEntry(entry);

      sel.setVz(vz_);
      sel.setHiHF(hiHF_);
      sel.setPprimaryVertexFilter(pprimaryVertexFilter_);
      sel.setPhfCoincFilter3(phfCoincFilter3_);
      sel.setHBHENoiseFilterResultRun2Loose(HBHENoiseFilterResultRun2Loose_);
      sel.setPclusterCompatibilityFilter(pclusterCompatibilityFilter_);

      if(!sel.isGood()) continue;

      Int_t centPos = -1;
      for(Int_t cI = 0; cI < nCentBins; ++cI){
	if(hiBin_/2 >= centBinsLow[cI] && hiBin_/2 < centBinsHi[cI]){
	  centPos = cI;
	  break;
	}
      }

      if(centPos < 0) continue;

      for(Int_t fI = 0; fI < nFlowTypes; ++fI){
	if(rhoFlowFitParams_p.at(fI)->size() != 0){
	  double val = ROOT::Math::chisquared_cdf_c(rhoFlowFitParams_p.at(fI)->at(5), rhoFlowFitParams_p.at(fI)->at(6));
	  bool minProb = val > .05;
	  bool maxProb = val < .95;
	  
	  if(minProb && maxProb){
	    v2_h[centPos][fI]->Fill(rhoFlowFitParams_p.at(fI)->at(1));
	    v3_h[centPos][fI]->Fill(rhoFlowFitParams_p.at(fI)->at(3));
	  }
      }
	
	++totEntries;
	if(totEntries >= nMaxEntries) break;
      }
    }

    inFile_p->Close();
    delete inFile_p;

    if(totEntries >= nMaxEntries) break;
  }

  TFile* f_raw = new TFile("/afs/cern.ch/work/c/cmcginn/private/Projects/CMSRooJR/v2v3data/V2_eta1p0_Raw.root", "READ");
  CustomCanvas* c_temp = new CustomCanvas("c_temp","",600,600);
  TF1* f_fit=new TF1("f_gaus","gaus",0,0.3);
  TH1F* h_dummy;
  outFile_p->cd();

  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);

  gSystem->cd("pdfDir");
  
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    for(Int_t fI = 0; fI < nFlowTypes; ++fI){
      v2_h[cI][fI]->Write("", TObject::kOverwrite);
      v3_h[cI][fI]->Write("", TObject::kOverwrite);
    }

    //    v2_h[cI]=(TH1F*)v2_h[cI]->Rebin(100,((std::string)v2_h[cI]->GetName()+"rebin").c_str(),rebinning);
    v2_h[cI][0]->Fit(f_fit);
    TF1* f_v2=new TF1(("f_gaus_v2_"+std::to_string(cI)).c_str(),"gaus",(f_fit->GetParameter(1)>2*f_fit->GetParameter(2))?(f_fit->GetParameter(1)-2*f_fit->GetParameter(2)):0,f_fit->GetParameter(1)+2*f_fit->GetParameter(2));
    v2_h[cI][0]->Fit(f_v2,"r");
    h_dummy=new TH1F(("h_tmp_v2_"+std::to_string(cI)).c_str(),v2_h[cI][0]->GetTitle(),100,(f_fit->GetParameter(1)>2*f_fit->GetParameter(2))?(f_fit->GetParameter(1)-2*f_fit->GetParameter(2)):0,f_fit->GetParameter(1)+2*f_fit->GetParameter(2));
    h_dummy->SetMaximum(1.2*v2_h[cI][0]->GetBinContent(v2_h[cI][0]->GetMaximumBin()));
    h_dummy->Draw();
    v2_h[cI][0]->Draw("e same");
    f_v2->Draw("same");
    c_temp->SaveAs(("v2_Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]) + "_" + dateStr + ".pdf").c_str());
    gPad->SetLogy();
    h_dummy->Draw();
    v2_h[cI][0]->Draw("e same");
    f_v2->Draw("same");
    c_temp->SaveAs(("v2_Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]) + "_" + dateStr + "_log.pdf").c_str());
    gPad->SetLogy(0);
    v2_Mean_h->SetBinContent(cI+1, f_fit->GetParameter(1));
    v2_Mean_h->SetBinError(cI+1, f_fit->GetParError(1));
    //    v2_Mean_h->SetBinContent(cI+1, v2_h[cI][0]->GetMean());
    //    v2_Mean_h->SetBinError(cI+1, v2_h[cI][0]->GetMeanError());

    TH1F* h_tmp=(TH1F*)f_raw->Get(("qwebye/hVnFull_c"+std::to_string(cI+1)).c_str());
    //    h_tmp=(TH1F*)h_tmp->Rebin(100,((std::string)h_tmp->GetName()+"rebin").c_str(),rebinning);
    h_tmp->Fit(f_fit);
    TF1* f_raw=new TF1(("f_gaus_raw_"+std::to_string(cI)).c_str(),"gaus",(f_fit->GetParameter(1)>2*f_fit->GetParameter(2))?(f_fit->GetParameter(1)-2*f_fit->GetParameter(2)):0,f_fit->GetParameter(1)+2*f_fit->GetParameter(2));
    h_tmp->Fit(f_raw,"r");
    h_dummy=new TH1F(("h_tmp_raw_"+std::to_string(cI)).c_str(),h_tmp->GetTitle(),100,(f_fit->GetParameter(1)>2*f_fit->GetParameter(2))?(f_fit->GetParameter(1)-2*f_fit->GetParameter(2)):0,f_fit->GetParameter(1)+2*f_fit->GetParameter(2));
    h_dummy->SetMaximum(1.2*h_tmp->GetBinContent(h_tmp->GetMaximumBin()));
    h_dummy->Draw();
    h_tmp->Draw("e same");
    f_raw->Draw("same");
    c_temp->SaveAs(("raw_Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]) + "_" + dateStr + ".pdf").c_str());
    gPad->SetLogy();
    h_dummy->Draw();
    h_tmp->Draw("e same");
    f_raw->Draw("same");
    c_temp->SaveAs(("raw_Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]) + "_" + dateStr + "_log.pdf").c_str());
    gPad->SetLogy(0);
    v2_Raw_h->SetBinContent(cI+1, f_fit->GetParameter(1));
    v2_Raw_h->SetBinError(cI+1, f_fit->GetParError(1));
    //    v2_Raw_h->SetBinContent(cI+1, h_tmp->GetMean());
    //    v2_Raw_h->SetBinError(cI+1, h_tmp->GetMeanError());

    //    v3_h[cI][0]=(TH1F*)v3_h[cI][0]->Rebin(100,((std::string)v3_h[cI][0]->GetName()+"rebin").c_str(),rebinning);
    v3_h[cI][0]->Fit(f_fit);
    TF1* f_v3=new TF1(("f_gaus_v3_"+std::to_string(cI)).c_str(),"gaus",(f_fit->GetParameter(1)>2*f_fit->GetParameter(2))?(f_fit->GetParameter(1)-2*f_fit->GetParameter(2)):0,f_fit->GetParameter(1)+2*f_fit->GetParameter(2));
    v3_h[cI][0]->Fit(f_v3,"r");
    h_dummy=new TH1F(("h_tmp_v3_"+std::to_string(cI)).c_str(),v3_h[cI][0]->GetTitle(),100,(f_fit->GetParameter(1)>2*f_fit->GetParameter(2))?(f_fit->GetParameter(1)-2*f_fit->GetParameter(2)):0,f_fit->GetParameter(1)+2*f_fit->GetParameter(2));
    h_dummy->SetMaximum(1.2*v3_h[cI][0]->GetBinContent(v3_h[cI][0]->GetMaximumBin()));
    h_dummy->Draw();
    v3_h[cI][0]->Draw("e same");
    f_v3->Draw("same");
    c_temp->SaveAs(("v3_Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]) + "_" + dateStr + ".pdf").c_str());
    gPad->SetLogy();
    h_dummy->Draw();
    v3_h[cI][0]->Draw("e same");
    f_v3->Draw("same");
    c_temp->SaveAs(("v3_Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]) + "_" + dateStr + "_log.pdf").c_str());
    gPad->SetLogy(0);
    v3_Mean_h->SetBinContent(cI+1, f_fit->GetParameter(1));
    v3_Mean_h->SetBinError(cI+1, f_fit->GetParError(1));
    
    v2_HIN16019_h->SetBinContent(cI+1, centYVals16019[cI]);
    v2_HIN16019_h->SetBinError(cI+1, 0.);

    delete v2_h[cI][0];
    delete v3_h[cI][0];
  }

  v2_Mean_h->Write("", TObject::kOverwrite);
  v2_Raw_h->Write("", TObject::kOverwrite);
  v3_Mean_h->Write("", TObject::kOverwrite);

  v2_HIN16019_h->Write("", TObject::kOverwrite);

  TexSlides(new std::vector<std::vector<std::string>*>{c_temp->GetPointer()},("Slides_" + dateStr + ".tex").c_str(),2);
  f_raw->Close();

  gSystem->cd("..");

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

  v2_Raw_h->SetMarkerStyle(20);
  v2_Raw_h->SetMarkerSize(0.8);
  v2_Raw_h->SetMarkerColor(col.getColor(1));
  v2_Raw_h->SetLineColor(col.getColor(1));

  v2_Raw_h->DrawCopy("HIST E1 P SAME");

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
  leg_p->AddEntry(v2_Raw_h, "HIN-16-019 (RAW)", "P L");
  leg_p->AddEntry(v2_HIN16019_h, "HIN-16-019 (UNFOLDED)", "P L");

  leg_p->Draw("SAME");

  std::string saveStr = "pdfDir/v2Overlay" + mcStr + "_" + dateStr + ".pdf";
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
