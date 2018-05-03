#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TMath.h"
#include "TDatime.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TSystem.h"

//#include "TrkCorr_Jun7_Iterative_PbPb_etaLT2p4/getTrkCorr.h""

#include "include/goodGlobalSelection.h"
#include "CustomCanvas.h"
#include "TexSlides.C"

int recreateV2V3(const std::string inFileName)
{
  const Int_t nCentBins = 5;
  const Int_t centBinsLow[nCentBins] = {0, 5, 10, 30, 50};
  const Int_t centBinsHi[nCentBins] = {5, 10, 30, 50, 100};

  CustomCanvas* c1=new CustomCanvas("c1","c1",400,400);

  //  TFile* corrFile_p = new TFile("inputs/EffCorrectionsPixel_NTT_pt_0_10_v2.root", "READ");
  //  TH1F* corrHists_p[nCentBins];
  //  for(Int_t cI = 0; cI < nCentBins; ++cI){
  //    corrHists_p[cI] = (TH1F*)corrFile_p->Get(("Eff_" + std::to_string(centBinsLow[cI]) + "_" + std::to_string(centBinsHi[cI])).c_str());
  //  }

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  TFile* outFile_p = new TFile(("output/v2CrossCheck_" + dateStr + ".root").c_str(), "RECREATE");
  TH1F* v2DeltaPos_h = new TH1F("v2DeltaPos_h", ";#Deltav_{2};Counts", 100, -.2, .2);
  TH1F* v2Pos_h = new TH1F("v2Pos_h", ";v_{2};Counts", 150, 0.0, 0.6);
  //  TH1F* v2EffPos_h = new TH1F("v2EffPos_h", ";v_{2};Counts", 150, 0.0, 0.6);

  TH1F* v2DeltaNeg_h = new TH1F("v2DeltaNeg_h", ";#Deltav_{2};Counts", 100, -.2, .2);
  TH1F* v2Neg_h = new TH1F("v2Neg_h", ";v_{2};Counts", 150, 0.0, 0.6);
  //  TH1F* v2EffNeg_h = new TH1F("v2EffNeg_h", ";v_{2};Counts", 150, 0.0, 0.6);

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* pfTree_p = (TTree*)inFile_p->Get("pfcandAnalyzer/pfTree");
  TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
  TTree* skimTree_p = (TTree*)inFile_p->Get("skimanalysis/HltTree");
  TTree* rhoFlowAnalyzer_p = (TTree*)inFile_p->Get("hiFlowAnalyzer/t");

  std::vector<float>* pfPt_p = NULL;
  std::vector<float>* pfEta_p = NULL;
  std::vector<float>* pfPhi_p = NULL;
  std::vector<int>* pfId_p = NULL;

  pfTree_p->SetBranchStatus("*", 0);
  pfTree_p->SetBranchStatus("pfPt", 1);
  pfTree_p->SetBranchStatus("pfEta", 1);
  pfTree_p->SetBranchStatus("pfPhi", 1);
  pfTree_p->SetBranchStatus("pfId", 1);

  pfTree_p->SetBranchAddress("pfPt", &pfPt_p);
  pfTree_p->SetBranchAddress("pfEta", &pfEta_p);
  pfTree_p->SetBranchAddress("pfPhi", &pfPhi_p);
  pfTree_p->SetBranchAddress("pfId", &pfId_p);

  Float_t hiHF_;
  Float_t vz_;
  Int_t hiBin_;
  Int_t hiNevtPlane_;
  Float_t hiEvtPlanes_[30];

  hiTree_p->SetBranchStatus("*", 0);
  hiTree_p->SetBranchStatus("hiHF", 1);
  hiTree_p->SetBranchStatus("hiBin", 1);
  hiTree_p->SetBranchStatus("vz", 1);
  hiTree_p->SetBranchStatus("hiNevtPlane", 1);
  hiTree_p->SetBranchStatus("hiEvtPlanes_", 1);

  hiTree_p->SetBranchAddress("hiHF", &hiHF_);
  hiTree_p->SetBranchAddress("hiBin", &hiBin_);
  hiTree_p->SetBranchAddress("vz", &vz_);
  hiTree_p->SetBranchAddress("hiNevtPlane", &hiNevtPlane_);
  hiTree_p->SetBranchAddress("hiEvtPlanes_", hiEvtPlanes_);

  Int_t pprimaryVertexFilter_;
  Int_t HBHENoiseFilterResultRun2Loose_;
  Int_t pclusterCompatibilityFilter_;
  Int_t phfCoincFilter3_;

  skimTree_p->SetBranchStatus("*", 0);
  skimTree_p->SetBranchStatus("pprimaryVertexFilter", 1);
  skimTree_p->SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
  skimTree_p->SetBranchStatus("phfCoincFilter3", 1);
  skimTree_p->SetBranchStatus("pclusterCompatibilityFilter", 1);

  skimTree_p->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter_);
  skimTree_p->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose_);
  skimTree_p->SetBranchAddress("phfCoincFilter3", &phfCoincFilter3_);
  skimTree_p->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter_);
  
  std::vector<double>* rhoFlowFitParams_p=NULL;

  rhoFlowAnalyzer_p->SetBranchStatus("*", 0);
  rhoFlowAnalyzer_p->SetBranchStatus("rhoFlowFitParams", 1);

  rhoFlowAnalyzer_p->SetBranchAddress("rhoFlowFitParams", &rhoFlowFitParams_p);

  const Int_t nEntries = TMath::Min((Int_t)pfTree_p->GetEntries(), (Int_t)50000);
  goodGlobalSelection sel;
  sel.setIsPbPb(true);

  //  TrkCorr* trkCorr = new TrkCorr("/export/d00/scratch/cfmcginn/JetRAA/RAACode/TrkCorr_Jun7_Iterative_PbPb_etaLT2p4/");

  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%10000 == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;
    pfTree_p->GetEntry(entry);
    hiTree_p->GetEntry(entry);
    skimTree_p->GetEntry(entry);
    rhoFlowAnalyzer_p->GetEntry(entry);

    sel.setVz(vz_);
    sel.setHiHF(hiHF_);
    sel.setPprimaryVertexFilter(pprimaryVertexFilter_);
    sel.setPhfCoincFilter3(phfCoincFilter3_);
    sel.setHBHENoiseFilterResultRun2Loose(HBHENoiseFilterResultRun2Loose_);
    sel.setPclusterCompatibilityFilter(pclusterCompatibilityFilter_);
    
    if(!sel.isGood()) continue;
    if(hiBin_ > 100) continue;
    if(hiNevtPlane_ < 20) continue;

    double eventPlane2 = hiEvtPlanes_[8];
    double eventPlane3 = hiEvtPlanes_[15];

    std::vector<float> phiPos_;
    std::vector<float> phiNeg_;
    //    std::vector<float> effPos_;
    //    std::vector<float> effNeg_;

    int centPos = -1;
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      if(centBinsLow[cI] <= hiBin_/2 && centBinsHi[cI] > hiBin_/2){
	centPos = cI;
	break;
      }
    }

    for(unsigned pfI = 0; pfI < pfPt_p->size(); pfI++){
      if(pfEta_p->at(pfI) < -1.0) continue;
      if(pfEta_p->at(pfI) > 1.0) continue;
      if(pfPt_p->at(pfI) < .3) continue;
      if(pfPt_p->at(pfI) > 3.) continue;
      if(pfId_p->at(pfI) != 1) continue;

      if (pfEta_p->at(pfI)>0){
	phiPos_.push_back(pfPhi_p->at(pfI));
	//	effPos_.push_back(corrHists_p[centPos]->GetBinContent(corrHists_p[centPos]->FindBin(pfEta_p->at(pfI), pfPt_p->at(pfI))));
      }
      else{
	phiNeg_.push_back(pfPhi_p->at(pfI));
	//        effNeg_.push_back(corrHists_p[centPos]->GetBinContent(corrHists_p[centPos]->FindBin(pfEta_p->at(pfI), pfPt_p->at(pfI))));
      }
    }

    const Int_t nPhiBinsPos = std::fmax(10, phiPos_.size()/30);
    TH1F* phiPos_h = new TH1F("tempPhiPos", ";#phi;Track counts (.3 < p_{T} < 3.,#eta>0)", nPhiBinsPos, -TMath::Pi(), TMath::Pi());
    //    TH1F* phiEffPos_h = new TH1F("tempPhiEffPos", ";#phi;Track counts (.3 < p_{T} < 3.,#eta>0)", nPhiBinsPos, -TMath::Pi(), TMath::Pi());

    const Int_t nPhiBinsNeg = std::fmax(10, phiNeg_.size()/30);
    TH1F* phiNeg_h = new TH1F("tempPhiNeg", ";#phi;Track counts (.3 < p_{T} < 3.,#eta<0)", nPhiBinsNeg, -TMath::Pi(), TMath::Pi());
    //    TH1F* phiEffNeg_h = new TH1F("tempPhiEffNeg", ";#phi;Track counts (.3 < p_{T} < 3.,#eta<0)", nPhiBinsNeg, -TMath::Pi(), TMath::Pi());

    for(unsigned int pfI = 0; pfI < phiPos_.size(); ++pfI){
      phiPos_h->Fill(phiPos_.at(pfI));
      //      phiEffPos_h->Fill(phiPos_.at(pfI), 1./effPos_.at(pfI));
    }

    for(unsigned int pfI = 0; pfI < phiNeg_.size(); ++pfI){
      phiNeg_h->Fill(phiNeg_.at(pfI));
      //      phiEffNeg_h->Fill(phiNeg_.at(pfI), 1./effNeg_.at(pfI));
    }

    std::string flowFitForm = "[0]*(1. + 2.*([1]*TMath::Cos(2.*(x - " + std::to_string(eventPlane2) + ")) + [2]*TMath::Cos(3.*(x - " + std::to_string(eventPlane3) + "))))";
    TF1* flowFitPos_p = new TF1("flowFitPos", flowFitForm.c_str(), -TMath::Pi(), TMath::Pi());
    flowFitPos_p->SetParameter(0, 10);
    flowFitPos_p->SetParameter(1, 0);
    flowFitPos_p->SetParameter(2, 0);

    /*    TF1* flowFitEffPos_p = new TF1("flowFitEffPos", flowFitForm.c_str(), -TMath::Pi(), TMath::Pi());
    flowFitEffPos_p->SetParameter(0, 10);
    flowFitEffPos_p->SetParameter(1, 0);
    flowFitEffPos_p->SetParameter(2, 0);
    */
    phiPos_h->Fit(flowFitPos_p, "Q", "", -TMath::Pi(), TMath::Pi());
    //    phiEffPos_h->Fit(flowFitEffPos_p, "Q", "", -TMath::Pi(), TMath::Pi());

    if(TMath::Abs(flowFitPos_p->GetParameter(1) - rhoFlowFitParams_p->at(1)) > .2) v2DeltaPos_h->Fill(v2DeltaPos_h->GetBinCenter(1));
    else v2DeltaPos_h->Fill(flowFitPos_p->GetParameter(1) - rhoFlowFitParams_p->at(1));

    v2Pos_h->Fill(flowFitPos_p->GetParameter(1));
    //    v2EffPos_h->Fill(flowFitEffPos_p->GetParameter(1));

    TF1* flowFitNeg_p = new TF1("flowFitNeg", flowFitForm.c_str(), -TMath::Pi(), TMath::Pi());
    flowFitNeg_p->SetParameter(0, 10);
    flowFitNeg_p->SetParameter(1, 0);
    flowFitNeg_p->SetParameter(2, 0);

    /*    TF1* flowFitEffNeg_p = new TF1("flowFitEffNeg", flowFitForm.c_str(), -TMath::Pi(), TMath::Pi());
    flowFitEffNeg_p->SetParameter(0, 10);
    flowFitEffNeg_p->SetParameter(1, 0);
    flowFitEffNeg_p->SetParameter(2, 0);
    */
    phiNeg_h->Fit(flowFitNeg_p, "Q", "", -TMath::Pi(), TMath::Pi());
    //    phiEffNeg_h->Fit(flowFitEffNeg_p, "Q", "", -TMath::Pi(), TMath::Pi());

    if(TMath::Abs(flowFitNeg_p->GetParameter(1) - rhoFlowFitParams_p->at(1)) > .2) v2DeltaNeg_h->Fill(v2DeltaNeg_h->GetBinCenter(1));
    else v2DeltaNeg_h->Fill(flowFitNeg_p->GetParameter(1) - rhoFlowFitParams_p->at(1));

    v2Neg_h->Fill(flowFitNeg_p->GetParameter(1));
    //    v2EffNeg_h->Fill(flowFitEffNeg_p->GetParameter(1));

    delete flowFitPos_p;
    //    delete flowFitEffPos_p;
    delete phiPos_h;
    //    delete phiEffPos_h;

    delete flowFitPos_p;
    //    delete flowFitEffPos_p;
    delete phiPos_h;
    //    delete phiEffPos_h;
  }

  inFile_p->Close();
  delete inFile_p;

  gSystem->cd("pdfDir");
  v2Pos_h->Draw();
  c1->SaveAs("v2Pos.jpg");
  v2Neg_h->Draw();
  c1->SaveAs("v2Neg.jpg");
  TexSlides(new std::vector<std::vector<std::string>*> {c1->GetPointer()},"Slides.tex",1);
  gSystem->cd("..");

  outFile_p->cd();
  v2DeltaPos_h->Write("", TObject::kOverwrite);
  v2Pos_h->Write("", TObject::kOverwrite);
  //  v2EffPos_h->Write("", TObject::kOverwrite);
  delete v2DeltaPos_h;
  delete v2Pos_h;
  //  delete v2EffPos_h;
  v2DeltaNeg_h->Write("", TObject::kOverwrite);
  v2Neg_h->Write("", TObject::kOverwrite);
  //  v2EffNeg_h->Write("", TObject::kOverwrite);
  delete v2DeltaNeg_h;
  delete v2Neg_h;
  //  delete v2EffNeg_h;
  outFile_p->Close();
  delete outFile_p;

  //  corrFile_p->Close();
  //  delete corrFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "USAGE: ./recreateV2V3.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += recreateV2V3(argv[1]);
  return retVal;
}
