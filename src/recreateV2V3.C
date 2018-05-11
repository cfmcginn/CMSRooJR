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

//#include "TrkCorr_Jun7_Iterative_PbPb_etaLT2p4/getTrkCorr.h""

#include "include/goodGlobalSelection.h"

int recreateV2V3(const std::string inFileName)
{
  const Int_t nCentBins = 5;
  const Int_t centBinsLow[nCentBins] = {0, 5, 10, 30, 50};
  const Int_t centBinsHi[nCentBins] = {5, 10, 30, 50, 100};

  TFile* corrFile_p = new TFile("input/EffCorrectionsPixel_TT_pt_0_10_v2.root", "READ");
  TH1F* corrHists_p[nCentBins];
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    corrHists_p[cI] = (TH1F*)corrFile_p->Get(("Eff_" + std::to_string(centBinsLow[cI]) + "_" + std::to_string(centBinsHi[cI])).c_str());
  }

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  TFile* outFile_p = new TFile(("output/v2CrossCheck_" + dateStr + ".root").c_str(), "RECREATE");
  TH1F* v2Delta_pf_h = new TH1F("v2Delta_pf_h", ";#Deltav_{2};Counts", 100, -.2, .2);
  TH1F* v2_pf_h = new TH1F("v2_pf_h", ";v_{2};Counts", 150, 0.0, 0.6);
  TH1F* v2Eff_pf_h = new TH1F("v2Eff_pf_h", ";v_{2};Counts", 150, 0.0, 0.6);
  TH1F* v2Delta_trk_h = new TH1F("v2Delta_trk_h", ";#Deltav_{2};Counts", 100, -.2, .2);
  TH1F* v2_trk_h = new TH1F("v2_trk_h", ";v_{2};Counts", 150, 0.0, 0.6);
  TH1F* v2Eff_trk_h = new TH1F("v2Eff_trk_h", ";v_{2};Counts", 150, 0.0, 0.6);

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* pfTree_p = (TTree*)inFile_p->Get("pfcandAnalyzer/pfTree");
  TTree* trackTree_p = (TTree*)inFile_p->Get("anaTrack/trackTree");
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

  const int maxNTrk=10000;
  Int_t nTrk_;
  Float_t trkPt_[maxNTrk];
  Float_t trkPtError_[maxNTrk];
  Float_t trkEta_[maxNTrk];
  Float_t trkPhi_[maxNTrk];
  UChar_t trkNHit_[maxNTrk];
  Float_t trkDxy1_[maxNTrk];
  Float_t trkDxyError1_[maxNTrk];
  Float_t trkDz1_[maxNTrk];
  Float_t trkDzError1_[maxNTrk];
  Float_t trkChi2_[maxNTrk];
  UChar_t trkNdof_[maxNTrk];
  UChar_t trkNlayer_[maxNTrk];
  UChar_t trkAlgo_[maxNTrk];
  Bool_t highPurity_[maxNTrk];
  
  trackTree_p->SetBranchStatus("*", 0);
  trackTree_p->SetBranchStatus("nTrk", 1);
  trackTree_p->SetBranchStatus("trkPt", 1);
  trackTree_p->SetBranchStatus("trkPtError", 1);
  trackTree_p->SetBranchStatus("trkEta", 1);
  trackTree_p->SetBranchStatus("trkPhi", 1);
  trackTree_p->SetBranchStatus("trkNHit", 1);
  trackTree_p->SetBranchStatus("trkDxy1", 1);
  trackTree_p->SetBranchStatus("trkDxyError1", 1);
  trackTree_p->SetBranchStatus("trkDz1", 1);
  trackTree_p->SetBranchStatus("trkDzError1", 1);
  trackTree_p->SetBranchStatus("trkChi2", 1);
  trackTree_p->SetBranchStatus("trkNdof", 1);
  trackTree_p->SetBranchStatus("trkNlayer", 1);
  trackTree_p->SetBranchStatus("trkAlgo", 1);
  trackTree_p->SetBranchStatus("highPurity", 1);

  trackTree_p->SetBranchAddress("nTrk", &nTrk_);
  trackTree_p->SetBranchAddress("trkPt", trkPt_);
  trackTree_p->SetBranchAddress("trkPtError", trkPtError_);
  trackTree_p->SetBranchAddress("trkEta", trkEta_);
  trackTree_p->SetBranchAddress("trkPhi", trkPhi_);
  trackTree_p->SetBranchAddress("trkNHit", trkNHit_);
  trackTree_p->SetBranchAddress("trkDxy1", trkDxy1_);
  trackTree_p->SetBranchAddress("trkDxyError1", trkDxyError1_);
  trackTree_p->SetBranchAddress("trkDz1", trkDz1_);
  trackTree_p->SetBranchAddress("trkDzError1", trkDzError1_);
  trackTree_p->SetBranchAddress("trkChi2", trkChi2_);
  trackTree_p->SetBranchAddress("trkNdof", trkNdof_);
  trackTree_p->SetBranchAddress("trkNlayer", trkNlayer_);
  trackTree_p->SetBranchAddress("trkAlgo", trkAlgo_);
  trackTree_p->SetBranchAddress("highPurity", highPurity_);

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
  hiTree_p->SetBranchStatus("hiEvtPlanes", 1);

  hiTree_p->SetBranchAddress("hiHF", &hiHF_);
  hiTree_p->SetBranchAddress("hiBin", &hiBin_);
  hiTree_p->SetBranchAddress("vz", &vz_);
  hiTree_p->SetBranchAddress("hiNevtPlane", &hiNevtPlane_);
  hiTree_p->SetBranchAddress("hiEvtPlanes", hiEvtPlanes_);

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
    trackTree_p->GetEntry(entry);
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

    std::vector<float> phi_;
    std::vector<float> eff_;
    std::vector<float> phi_trk_;
    std::vector<float> eff_trk_;

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

      phi_.push_back(pfPhi_p->at(pfI));
      eff_.push_back(corrHists_p[centPos]->GetBinContent(corrHists_p[centPos]->FindBin(pfEta_p->at(pfI), pfPt_p->at(pfI))));
    }

    const Int_t nPhiBins = std::fmax(10, phi_.size()/30);
    TH1F* phi_h = new TH1F("tempPhi", ";#phi;Track counts (.3 < p_{T} < 3.)", nPhiBins, -TMath::Pi(), TMath::Pi());
    TH1F* phiEff_h = new TH1F("tempPhiEff", ";#phi;Track counts (.3 < p_{T} < 3.)", nPhiBins, -TMath::Pi(), TMath::Pi());

    for(unsigned int pfI = 0; pfI < phi_.size(); ++pfI){
      phi_h->Fill(phi_.at(pfI));
      phiEff_h->Fill(phi_.at(pfI), 1./eff_.at(pfI));
    }

    std::string flowFitForm = "[0]*(1. + 2.*([1]*TMath::Cos(2.*(x - " + std::to_string(eventPlane2) + ")) + [2]*TMath::Cos(3.*(x - " + std::to_string(eventPlane3) + "))))";
    TF1* flowFit_p = new TF1("flowFit", flowFitForm.c_str(), -TMath::Pi(), TMath::Pi());
    flowFit_p->SetParameter(0, 10);
    flowFit_p->SetParameter(1, 0);
    flowFit_p->SetParameter(2, 0);

    TF1* flowFitEff_p = new TF1("flowFitEff", flowFitForm.c_str(), -TMath::Pi(), TMath::Pi());
    flowFitEff_p->SetParameter(0, 10);
    flowFitEff_p->SetParameter(1, 0);
    flowFitEff_p->SetParameter(2, 0);

    phi_h->Fit(flowFit_p, "Q", "", -TMath::Pi(), TMath::Pi());
    phiEff_h->Fit(flowFitEff_p, "Q", "", -TMath::Pi(), TMath::Pi());

    if(TMath::Abs(flowFit_p->GetParameter(1) - rhoFlowFitParams_p->at(1)) > .2) v2Delta_pf_h->Fill(v2Delta_pf_h->GetBinCenter(1));
    else v2Delta_pf_h->Fill(flowFit_p->GetParameter(1) - rhoFlowFitParams_p->at(1));

    v2_pf_h->Fill(flowFit_p->GetParameter(1));
    v2Eff_pf_h->Fill(flowFitEff_p->GetParameter(1));

    for(Int_t tI = 0; tI < nTrk_; tI++){
      if(trkEta_[tI] < -1.0) continue;
      if(trkEta_[tI] > 1.0) continue;
      if(trkPt_[tI] < 0.3) continue;
      if(trkPt_[tI] > 3.0) continue;
      if(!highPurity_[tI]) continue;

      if(trkNHit_[tI] >= 3 && trkNHit_[tI] <= 6 && trkPt_[tI] < 2.4){
	if(TMath::Abs(trkDz1_[tI]/trkDzError1_[tI]) >= 8.) continue;
	if(trkChi2_[tI]/(Double_t)(trkNdof_[tI]*trkNlayer_[tI]) >= 12.) continue;
      }
      else{
	if(trkPt_[tI] > 2.4 && (trkAlgo_[tI] < 4 || trkAlgo_[tI] > 7)) continue;
	if(TMath::Abs(trkDz1_[tI]/trkDzError1_[tI]) >= 3.) continue;
	if(TMath::Abs(trkDxy1_[tI]/trkDxyError1_[tI]) >= 3.) continue;
	if(TMath::Abs(trkPtError_[tI]/trkPt_[tI]) > .1) continue;
	if(trkNHit_[tI] < 11) continue;
	if(trkChi2_[tI]/(Double_t)(trkNdof_[tI]*trkNlayer_[tI]) >= .15) continue;
      }

      phi_trk_.push_back(trkPhi_[tI]);
      eff_trk_.push_back(corrHists_p[centPos]->GetBinContent(corrHists_p[centPos]->FindBin(trkEta_[tI], trkPt_[tI])));
    }

    const Int_t nTrkBins = std::fmax(10, phi_trk_.size()/30);
    TH1F* phi_trk_h = new TH1F("tempPhiTrk", ";#phi;Track counts (.3 < p_{T} < 3.)", nTrkBins, -TMath::Pi(), TMath::Pi());
    TH1F* phiEff_trk_h = new TH1F("tempPhiEffTrk", ";#phi;Track counts (.3 < p_{T} < 3.)", nTrkBins, -TMath::Pi(), TMath::Pi());

    for(unsigned int tI = 0; tI < phi_trk_.size(); ++tI){
      phi_trk_h->Fill(phi_trk_.at(tI));
      phiEff_trk_h->Fill(phi_trk_.at(tI), 1./eff_trk_.at(tI));
    }

    TF1* flowFit_trk_p = new TF1("flowFitTrk", flowFitForm.c_str(), -TMath::Pi(), TMath::Pi());
    flowFit_trk_p->SetParameter(0, 10);
    flowFit_trk_p->SetParameter(1, 0);
    flowFit_trk_p->SetParameter(2, 0);

    TF1* flowFitEff_trk_p = new TF1("flowFitEffTrk", flowFitForm.c_str(), -TMath::Pi(), TMath::Pi());
    flowFitEff_trk_p->SetParameter(0, 10);
    flowFitEff_trk_p->SetParameter(1, 0);
    flowFitEff_trk_p->SetParameter(2, 0);

    phi_trk_h->Fit(flowFit_trk_p, "Q", "", -TMath::Pi(), TMath::Pi());
    phiEff_trk_h->Fit(flowFitEff_trk_p, "Q", "", -TMath::Pi(), TMath::Pi());

    if(TMath::Abs(flowFit_trk_p->GetParameter(1) - rhoFlowFitParams_p->at(1)) > .2) v2Delta_trk_h->Fill(v2Delta_trk_h->GetBinCenter(1));
    else v2Delta_trk_h->Fill(flowFit_trk_p->GetParameter(1) - rhoFlowFitParams_p->at(1));

    v2_trk_h->Fill(flowFit_trk_p->GetParameter(1));
    v2Eff_trk_h->Fill(flowFitEff_trk_p->GetParameter(1));

    delete flowFit_p;
    delete flowFitEff_p;
    delete flowFit_trk_p;
    delete flowFitEff_trk_p;
    delete phi_h;
    delete phiEff_h;
    delete phi_trk_h;
    delete phiEff_trk_h;
  }

  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();
  v2Delta_pf_h->Write("", TObject::kOverwrite);
  v2_pf_h->Write("", TObject::kOverwrite);
  v2Eff_pf_h->Write("", TObject::kOverwrite);
  v2Delta_trk_h->Write("", TObject::kOverwrite);
  v2_trk_h->Write("", TObject::kOverwrite);
  v2Eff_trk_h->Write("", TObject::kOverwrite);
  delete v2Delta_pf_h;
  delete v2_pf_h;
  delete v2Eff_pf_h;
  delete v2Delta_trk_h;
  delete v2_trk_h;
  delete v2Eff_trk_h;
  outFile_p->Close();
  delete outFile_p;

  corrFile_p->Close();
  delete corrFile_p;

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
