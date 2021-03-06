#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TMath.h"
#include "TDatime.h"
#include "TF1.h"

#include "Math/ProbFuncMathCore.h"

#include "include/doGlobalDebug.h"
#include "include/goodGlobalSelection.h"
#include "include/returnRootFileContentsList.h"
#include "include/vectorStringUtility.h"
#include "include/mntToXRootdFileString.h"

int recreateV2V3Tree(const std::string inFileName, std::string mntFileName)
{
  const Int_t nCentBinsCorr = 5;
  const Int_t centBinsCorrLow[nCentBinsCorr] = {0, 5, 10, 30, 50};
  const Int_t centBinsCorrHi[nCentBinsCorr] = {5, 10, 30, 50, 100};

  //  const std::string fullPath = std::getenv("FULLJRDIR");
  TFile* corrFile_p = new TFile("input/EffCorrectionsPixel_TT_pt_0_10_v2.root", "READ");
  TH1F* corrHists_p[nCentBinsCorr];
  for(Int_t cI = 0; cI < nCentBinsCorr; ++cI){
    corrHists_p[cI] = (TH1F*)corrFile_p->Get(("Eff_" + std::to_string(centBinsCorrLow[cI]) + "_" + std::to_string(centBinsCorrHi[cI])).c_str());
  }

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;


  TFile* outFile_p = new TFile(("/afs/cern.ch/work/s/skanaski/jetstuff/output/v2V3Tree_" + dateStr + ".root").c_str(), "RECREATE");
  TTree* v2V3Tree_p = new TTree("v2V3Tree", "");
  TTree* mntTree_p = new TTree("mntTree", "");

  Int_t hiBinOut_;
  Float_t hiEvt2Plane_;
  Float_t hiEvt3Plane_;
  Float_t v2FromTree_;
  Float_t chi2FromTree_;
  Float_t nDOFFromTree_;
  std::vector<float>* pfPtOut_p=NULL;
  std::vector<float>* pfPhiOut_p=NULL;
  std::vector<float>* pfWeightOut_p=NULL;

  std::vector<float>* trkPtOut_p=NULL;
  std::vector<float>* trkPhiOut_p=NULL;
  std::vector<float>* trkWeightOut_p=NULL;

  int entry_;

  v2V3Tree_p->Branch("hiBin", &hiBinOut_, "hiBin/I");
  v2V3Tree_p->Branch("hiEvt2Plane", &hiEvt2Plane_, "hiEvt2Plane/F");
  v2V3Tree_p->Branch("hiEvt3Plane", &hiEvt3Plane_, "hiEvt3Plane/F");
  v2V3Tree_p->Branch("v2FromTree", &v2FromTree_, "v2FromTree/F");
  v2V3Tree_p->Branch("chi2FromTree", &chi2FromTree_, "chi2FromTree/F");
  v2V3Tree_p->Branch("nDOFFromTree", &nDOFFromTree_, "nDOFFromTree/F");
  v2V3Tree_p->Branch("pfPt", &pfPtOut_p);
  v2V3Tree_p->Branch("pfPhi", &pfPhiOut_p);
  v2V3Tree_p->Branch("pfWeight", &pfWeightOut_p);
  v2V3Tree_p->Branch("entry", &entry_);
 
  TFile* inFile_p =new TFile(inFileName.c_str(), "READ");
  //  TFile* mntFile_p = TFile::Open(mntToXRootdFileString(mntFileName).c_str(),"READ");
  TFile* mntFile_p = new TFile(mntFileName.c_str(), "READ");
  std::vector<std::string> ttreeList = returnRootFileContentsList(inFile_p, "TTree");
  removeVectorDuplicates(&ttreeList);
  //  bool hasTrackTree = vectorContainsString(&ttreeList, "track");
  bool hasTrackTree = true;

  if(hasTrackTree){
    mntTree_p->Branch("hiBin", &hiBinOut_, "hiBin/I");
    mntTree_p->Branch("trkPt", &trkPtOut_p);
    mntTree_p->Branch("trkPhi", &trkPhiOut_p);
    mntTree_p->Branch("trkWeight", &trkWeightOut_p);
    mntTree_p->Branch("entry", &entry_);
    mntTree_p->Branch("hiEvt2Plane", &hiEvt2Plane_, "hiEvt2Plane/F");
  }

  TTree* pfTree_p = (TTree*)inFile_p->Get("pfcandAnalyzer/pfTree");
  TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
  TTree* skimTree_p = (TTree*)inFile_p->Get("skimanalysis/HltTree");
  TTree* rhoFlowAnalyzer_p = (TTree*)inFile_p->Get("hiFlowAnalyzerHFPlane/t");
  TTree* trackTree_p = NULL;
  TTree* mntHiTree_p = NULL;
  TTree* mntSkimTree_p = NULL;
  //  if(hasTrackTree) trackTree_p = (TTree*)mntFile_p->Get("anaTrack/trackTree");
  if (hasTrackTree){
    trackTree_p = (TTree*)mntFile_p->Get("EbyEana/tree");
    mntHiTree_p = (TTree*)mntFile_p->Get("hiEvtAnalyzer/HiTree");
    mntSkimTree_p = (TTree*)mntFile_p->Get("skimanalysis/HltTree");
  }

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

  const Int_t nMaxTrk = 40000;
  Int_t nTrk_;
  Float_t trkPt_[nMaxTrk];
  Float_t trkPhi_[nMaxTrk];
  Float_t trkEta_[nMaxTrk];
  //  UChar_t trkNHit_[nMaxTrk];
  //  Bool_t highPurity_[nMaxTrk];
  //  Float_t trkPtError_[nMaxTrk];
  //  Float_t trkDz1_[nMaxTrk];
  //  Float_t trkDzError1_[nMaxTrk];
  //  Float_t trkDxy1_[nMaxTrk];
  //  Float_t trkDxyError1_[nMaxTrk];
  //  UChar_t trkNlayer_[nMaxTrk];
  //  Float_t trkChi2_[nMaxTrk];
  //  UChar_t trkNdof_[nMaxTrk];
  //  UChar_t trkOriginalAlgo_[nMaxTrk];

  if(hasTrackTree){
    trackTree_p->SetBranchStatus("*", 0);
    trackTree_p->SetBranchStatus("nTrk", 1);
    trackTree_p->SetBranchStatus("trkPt", 1);
    trackTree_p->SetBranchStatus("trkPhi", 1);
    trackTree_p->SetBranchStatus("trkEta", 1);
    //    trackTree_p->SetBranchStatus("trkNHit", 1);
    //    trackTree_p->SetBranchStatus("highPurity", 1);
    //    trackTree_p->SetBranchStatus("trkPtError", 1);
    //    trackTree_p->SetBranchStatus("trkDz1", 1);
    //    trackTree_p->SetBranchStatus("trkDzError1", 1);
    //    trackTree_p->SetBranchStatus("trkDxy1", 1);
    //    trackTree_p->SetBranchStatus("trkDxyError1", 1);
    //    trackTree_p->SetBranchStatus("trkNlayer", 1);
    //    trackTree_p->SetBranchStatus("trkChi2", 1);
    //    trackTree_p->SetBranchStatus("trkNdof", 1);
    //    trackTree_p->SetBranchStatus("trkOriginalAlgo", 1);

    trackTree_p->SetBranchAddress("nTrk", &nTrk_);
    trackTree_p->SetBranchAddress("trkPt", trkPt_);
    trackTree_p->SetBranchAddress("trkPhi", trkPhi_);
    trackTree_p->SetBranchAddress("trkEta", trkEta_);
    //    trackTree_p->SetBranchAddress("trkNHit", trkNHit_);
    //    trackTree_p->SetBranchAddress("highPurity", highPurity_);
    //    trackTree_p->SetBranchAddress("trkPtError", trkPtError_);
    //    trackTree_p->SetBranchAddress("trkDz1", trkDz1_);
    //    trackTree_p->SetBranchAddress("trkDzError1", trkDzError1_);
    //    trackTree_p->SetBranchAddress("trkDxy1", trkDxy1_);
    //    trackTree_p->SetBranchAddress("trkDxyError1", trkDxyError1_);
    //    trackTree_p->SetBranchAddress("trkNlayer", trkNlayer_);
    //    trackTree_p->SetBranchAddress("trkChi2", trkChi2_);
    //    trackTree_p->SetBranchAddress("trkNdof", trkNdof_);
    //    trackTree_p->SetBranchAddress("trkOriginalAlgo", trkOriginalAlgo_);
  }

  mntHiTree_p->SetBranchStatus("*", 0);
  mntHiTree_p->SetBranchStatus("hiBin", 1);
  mntHiTree_p->SetBranchStatus("vz", 1);
  mntHiTree_p->SetBranchStatus("hiHF", 1);
  mntHiTree_p->SetBranchStatus("hiNevtPlane", 1);
  mntHiTree_p->SetBranchStatus("hiEvtPlanes", 1);

  mntHiTree_p->SetBranchAddress("hiBin", &hiBin_);
  mntHiTree_p->SetBranchAddress("vz", &vz_);
  mntHiTree_p->SetBranchAddress("hiHF", &hiHF_);
  mntHiTree_p->SetBranchAddress("hiNevtPlane", &hiNevtPlane_);
  mntHiTree_p->SetBranchAddress("hiEvtPlanes", &hiEvtPlanes_);

  mntSkimTree_p->SetBranchStatus("*", 0);
  mntSkimTree_p->SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
  mntSkimTree_p->SetBranchStatus("pprimaryVertexFilter", 1);
  mntSkimTree_p->SetBranchStatus("phfCoincFilter3", 1);
  mntSkimTree_p->SetBranchStatus("pclusterCompatibilityFilter", 1);

  mntSkimTree_p->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose_);
  mntSkimTree_p->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter_);
  mntSkimTree_p->SetBranchAddress("phfCoincFilter3", &phfCoincFilter3_);
  mntSkimTree_p->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter_);

  const Int_t nEntries = TMath::Min((Int_t)pfTree_p->GetEntries(), (Int_t)10000000000);
  const Int_t nMntEntries = TMath::Min((Int_t)trackTree_p->GetEntries(), (Int_t)10000000000);
  goodGlobalSelection sel;
  sel.setIsPbPb(true);
  int counter = 0;

  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%10000 == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;
    pfTree_p->GetEntry(entry);
    hiTree_p->GetEntry(entry);
    skimTree_p->GetEntry(entry);
    rhoFlowAnalyzer_p->GetEntry(entry);
    //if(hasTrackTree) trackTree_p->GetEntry(entry);

    sel.setVz(vz_);
    sel.setHiHF(hiHF_);
    sel.setPprimaryVertexFilter(pprimaryVertexFilter_);
    sel.setPhfCoincFilter3(phfCoincFilter3_);
    sel.setHBHENoiseFilterResultRun2Loose(HBHENoiseFilterResultRun2Loose_);
    sel.setPclusterCompatibilityFilter(pclusterCompatibilityFilter_);
    
    if(!sel.isGood()) continue;
    if(hiBin_ > 120) continue;
    if(hiNevtPlane_ < 20) continue;

    counter++;

    hiBinOut_ = hiBin_;
    hiEvt2Plane_ = hiEvtPlanes_[8];
    hiEvt3Plane_ = hiEvtPlanes_[15];
    v2FromTree_ = rhoFlowFitParams_p->at(1);
    chi2FromTree_ = rhoFlowFitParams_p->at(5);
    nDOFFromTree_ = rhoFlowFitParams_p->at(6);
    entry_ = entry;

    int centCorrPos = -1;
    for(Int_t cI = 0; cI < nCentBinsCorr; ++cI){
      if(centBinsCorrLow[cI] <= hiBin_/2 && centBinsCorrHi[cI] > hiBin_/2){
	centCorrPos = cI;
	break;
      }
    }
      
    pfPtOut_p->clear();
    pfPhiOut_p->clear();
    pfWeightOut_p->clear();

    for(unsigned int pfI = 0; pfI < pfPt_p->size(); pfI++){
      if(pfEta_p->at(pfI) < -1.0) continue;
      if(pfEta_p->at(pfI) > 1.0) continue;
      if(pfPt_p->at(pfI) < .3) continue;
      if(pfPt_p->at(pfI) > 3.) continue;
      if(pfId_p->at(pfI) != 1) continue;

      double tempWeight = 1./corrHists_p[centCorrPos]->GetBinContent(corrHists_p[centCorrPos]->FindBin(pfEta_p->at(pfI), pfPt_p->at(pfI)));

      pfPtOut_p->push_back(pfPt_p->at(pfI));
      pfPhiOut_p->push_back(pfPhi_p->at(pfI));
      pfWeightOut_p->push_back(tempWeight);
    }

    v2V3Tree_p->Fill();
    //    if (counter==100) break;
  }
  
  counter = 0;
  for(Int_t entry = 0; entry < nMntEntries; ++entry){
    if(entry%10000 == 0) std::cout << " Mnt Entry " << entry << "/" << nMntEntries << std::endl;
    trackTree_p->GetEntry(entry);
    mntHiTree_p->GetEntry(entry);
    mntSkimTree_p->GetEntry(entry);

    sel.setVz(vz_);
    sel.setHiHF(hiHF_);
    sel.setPprimaryVertexFilter(pprimaryVertexFilter_);
    sel.setPhfCoincFilter3(phfCoincFilter3_);
    sel.setHBHENoiseFilterResultRun2Loose(HBHENoiseFilterResultRun2Loose_);
    sel.setPclusterCompatibilityFilter(pclusterCompatibilityFilter_);

    if(!sel.isGood()) continue;
    if(hiBin_ > 120) continue;
    if(hiNevtPlane_ < 20) continue;

    counter++;

    hiBinOut_ = hiBin_;
    entry_ = entry;
    hiEvt2Plane_ = hiEvtPlanes_[8];

    int centCorrPos = -1;
    for(Int_t cI = 0; cI < nCentBinsCorr; ++cI){
      if(centBinsCorrLow[cI] <= hiBin_/2 && centBinsCorrHi[cI] > hiBin_/2){
        centCorrPos = cI;
        break;
      }
    }

    trkPtOut_p->clear();
    trkPhiOut_p->clear();
    trkWeightOut_p->clear();

    for(Int_t tI = 0; tI < nTrk_; tI++){
      if(trkEta_[tI] < -1.0) continue;
      if(trkEta_[tI] > 1.0) continue;
      if(trkPt_[tI] < .3) continue;
      if(trkPt_[tI] > 3.) continue;
      /*      if(!highPurity_[tI]) continue;

      if(trkNHit_[tI] >= 3 && trkNHit_[tI] <= 6 && trkPt_[tI] < 2.4){
	if(TMath::Abs(trkDz1_[tI]/trkDzError1_[tI]) >= 8.) continue;
	if(trkChi2_[tI]/(Double_t)(trkNdof_[tI]*trkNlayer_[tI]) >= 12.) continue;
      }      
      else{
	if(trkPt_[tI] > 2.4 && (trkOriginalAlgo_[tI] < 4 || trkOriginalAlgo_[tI] > 7)) continue;
	//	if(trkOriginalAlgo_[tI] < 4 && trkOriginalAlgo_[tI] > 7) continue;
	if(TMath::Abs(trkDz1_[tI]/trkDzError1_[tI]) >= 3.) continue;
	if(TMath::Abs(trkDxy1_[tI]/trkDxyError1_[tI]) >= 3.) continue;
	if(TMath::Abs(trkPtError_[tI]/trkPt_[tI]) > .1) continue;
	if(trkNHit_[tI] < 11) continue;
	if(trkChi2_[tI]/(Double_t)(trkNdof_[tI]*trkNlayer_[tI]) >= .15) continue;
	}*/

      double tempWeight = 1./corrHists_p[centCorrPos]->GetBinContent(corrHists_p[centCorrPos]->FindBin(trkEta_[tI], trkPt_[tI]));

      trkPtOut_p->push_back(trkPt_[tI]);
      trkPhiOut_p->push_back(trkPhi_[tI]);
      trkWeightOut_p->push_back(tempWeight);
    }

    mntTree_p->Fill();
    //    if (counter==100) break;
  }

  inFile_p->Close();
  delete inFile_p;
  mntFile_p->Close();
  delete mntFile_p;

  outFile_p->cd();
  v2V3Tree_p->Write("", TObject::kOverwrite);
  mntTree_p->Write("", TObject::kOverwrite);
  outFile_p->Close();
  delete outFile_p;

  corrFile_p->Close();
  delete corrFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "USAGE: ./recreateV2V3Tree.exe <inFileName> <mntFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += recreateV2V3Tree(argv[1], argv[2]);
  return retVal;
}
