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

#include "include/doGlobalDebug.h"

bool dR(float eta1, float phi1, float eta2, float phi2, float R){
  return ((eta1-eta2)*(eta1-eta2)+(phi1-phi2)*(phi1-phi2)<=R*R);
}

int recreateV2V3TreeHist(const std::string inFileName)
{
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  TFile* outFile_p = new TFile(("output/v2CrossCheck_TreeHist_" + dateStr + ".root").c_str(), "RECREATE");

  const Int_t nCentBins = 11;
  const Int_t centBinsLow[nCentBins] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55};
  const Int_t centBinsHi[nCentBins] = {10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60};
  const Double_t centBins[nCentBins+1] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60};

  const Int_t nSystR = 3;
  const Double_t systR[nSystR] = {0.2, 0.4, 0.6};
  const Int_t nStrips = 12;

  TH1F* v2Raw_PF_h[nCentBins];
  TH1F* v2RawCorr_PF_h[nCentBins];
  TH1F* v2Fit_PF_h[nCentBins];
  TH1F* v2Diff_PF_h[nCentBins];
  TH1F* v2DiffCorr_PF_h[nCentBins];

  TH1F* systRaw_h[nCentBins][nSystR];
  TH1F* systFit_h[nCentBins][nSystR];

  TH1F* v2Obs_PF_h[nCentBins];
  TH1F* v2ObsCorr_PF_h[nCentBins];

  TH1F* v2Raw_Trk_h[nCentBins];
  TH1F* v2RawCorr_Trk_h[nCentBins];
  
  TH1F* v2Obs_Trk_h[nCentBins];
  TH1F* v2ObsCorr_Trk_h[nCentBins];

  TH1F* compRaw_h[nCentBins];
  TH1F* compRawCorr_h[nCentBins];
  
  TH1F* compObs_h[nCentBins];
  TH1F* compObsCorr_h[nCentBins];

  TH1F* v3Raw_PF_h[nCentBins];
  TH1F* v3RawCorr_PF_h[nCentBins];
  TH1F* v3Fit_PF_h[nCentBins];
  TH1F* v3Diff_PF_h[nCentBins];
  TH1F* v3DiffCorr_PF_h[nCentBins];

  TH1F* v3Obs_PF_h[nCentBins];
  TH1F* v3ObsCorr_PF_h[nCentBins];

  TH1F* v2Raw_Mean_PF_h = new TH1F("v2Raw_Mean_PF_h", ";Centrality (%);#LTv_{2}^{raw}#GT", nCentBins, centBins);
  TH1F* v2Raw_Sigma_PF_h = new TH1F("v2Raw_Sigma_PF_h", ";Centrality (%);#sigma(v_{2}^{raw})", nCentBins, centBins);
  TH1F* v2RawCorr_Mean_PF_h = new TH1F("v2RawCorr_Mean_PF_h", ";Centrality (%);#LTv_{2}^{raw}#GT", nCentBins, centBins);
  TH1F* v2RawCorr_Sigma_PF_h = new TH1F("v2RawCorr_Sigma_PF_h", ";Centrality (%);#sigma(v_{2}^{raw})", nCentBins, centBins);
  TH1F* v2Fit_Mean_PF_h = new TH1F("v2Fit_Mean_PF_h",";Centrality (%);#LTv_{2}^{fit}#GT", nCentBins, centBins);
  TH1F* v2Fit_Sigma_PF_h = new TH1F("v2Fit_Sigma_PF_h",";Centrality (%);#sigma(v_{2}^{fit})", nCentBins, centBins);

  TH1F* v2Obs_Mean_PF_h = new TH1F("v2Obs_Mean_PF_h", ";Centrality (%);#LTv_{2}^{obs}#GT", nCentBins, centBins);
  TH1F* v2Obs_Sigma_PF_h = new TH1F("v2Obs_Sigma_PF_h", ";Centrality (%);#sigma(v_{2}^{obs})", nCentBins, centBins);
  TH1F* v2ObsCorr_Mean_PF_h = new TH1F("v2ObsCorr_Mean_PF_h", ";Centrality (%);#LTv_{2}^{obs}#GT", nCentBins, centBins);
  TH1F* v2ObsCorr_Sigma_PF_h = new TH1F("v2ObsCorr_Sigma_PF_h", ";Centrality (%);#sigma(v_{2}^{obs})", nCentBins, centBins);

  TH1F* v2Raw_Mean_Trk_h = new TH1F("v2Raw_Mean_Trk_h", ";Centrality (%);#LTv_{2}^{raw}#GT", nCentBins, centBins);
  TH1F* v2Raw_Sigma_Trk_h = new TH1F("v2Raw_Sigma_Trk_h", ";Centrality (%);#sigma(v_{2}^{raw})", nCentBins, centBins);
  TH1F* v2RawCorr_Mean_Trk_h = new TH1F("v2RawCorr_Mean_Trk_h", ";Centrality (%);#LTv_{2}^{raw}#GT", nCentBins, centBins);
  TH1F* v2RawCorr_Sigma_Trk_h = new TH1F("v2RawCorr_Sigma_Trk_h", ";Centrality (%);#sigma(v_{2}^{raw})", nCentBins, centBins);

  TH1F* v2Obs_Mean_Trk_h = new TH1F("v2Obs_Mean_Trk_h", ";Centrality (%);#LTv_{2}^{obs}#GT", nCentBins, centBins);
  TH1F* v2Obs_Sigma_Trk_h = new TH1F("v2Obs_Sigma_Trk_h", ";Centrality (%);#sigma(v_{2}^{obs})", nCentBins, centBins);
  TH1F* v2ObsCorr_Mean_Trk_h = new TH1F("v2ObsCorr_Mean_Trk_h", ";Centrality (%);#LTv_{2}^{obs}#GT", nCentBins, centBins);
  TH1F* v2ObsCorr_Sigma_Trk_h = new TH1F("v2ObsCorr_Sigma_Trk_h", ";Centrality (%);#sigma(v_{2}^{obs})", nCentBins, centBins);

  TH1F* compRaw_Mean_h = new TH1F("compRaw_Mean_h", ";Centrality (%);#LTv_{2,trk}^{raw}/v_{2,PF}^{raw}#GT", nCentBins, centBins);
  TH1F* compRaw_Sigma_h = new TH1F("compRaw_Sigma_h", ";Centrality (%);#sigma(v_{2,trk}^{raw}/v_{2,PF}^{raw})", nCentBins, centBins);
  TH1F* compRawCorr_Mean_h = new TH1F("compRawCorr_Mean_h", ";Centrality (%);#LTv_{2,trk}^{raw}/v_{2,PF}^{raw}#GT", nCentBins, centBins);
  TH1F* compRawCorr_Sigma_h = new TH1F("compRawCorr_Sigma_h", ";Centrality (%);#sigma(v_{2,trk}^{raw}/v_{2,PF}^{raw})", nCentBins, centBins);

  TH1F* compObs_Mean_h = new TH1F("compObs_Mean_h", ";Centrality (%);#LTv_{2,trk}^{obs}/v_{2,PF}^{obs}#GT", nCentBins, centBins);
  TH1F* compObs_Sigma_h = new TH1F("compObs_Sigma_h", ";Centrality (%);#sigma(v_{2,trk}^{obs}/v_{2,PF}^{obs})", nCentBins, centBins);
  TH1F* compObsCorr_Mean_h = new TH1F("compObsCorr_Mean_h", ";Centrality (%);#LTv_{2,trk}^{obs}/v_{2,PF}^{obs}#GT", nCentBins, centBins);
  TH1F* compObsCorr_Sigma_h = new TH1F("compObsCorr_Sigma_h", ";Centrality (%);#sigma(v_{2,trk}^{obs}/v_{2,PF}^{obs})", nCentBins, centBins);

  TH1F* v3Raw_Mean_PF_h = new TH1F("v3Raw_Mean_PF_h", ";Centrality (%);#LTv_{3}^{raw}#GT", nCentBins, centBins);
  TH1F* v3Raw_Sigma_PF_h = new TH1F("v3Raw_Sigma_PF_h", ";Centrality (%);#sigma(v_{3}^{raw})", nCentBins, centBins);
  TH1F* v3RawCorr_Mean_PF_h = new TH1F("v3RawCorr_Mean_PF_h", ";Centrality (%);#LTv_{3}^{raw}#GT", nCentBins, centBins);
  TH1F* v3RawCorr_Sigma_PF_h = new TH1F("v3RawCorr_Sigma_PF_h", ";Centrality (%);#sigma(v_{3}^{raw})", nCentBins, centBins);
  TH1F* v3Fit_Mean_PF_h = new TH1F("v3Fit_Mean_PF_h",";Centrality (%);#LTv_{3}^{fit}#GT", nCentBins, centBins);
  TH1F* v3Fit_Sigma_PF_h = new TH1F("v3Fit_Sigma_PF_h",";Centrality (%);#sigma(v_{3}^{fit})", nCentBins, centBins);

  TH1F* v3Obs_Mean_PF_h = new TH1F("v3Obs_Mean_PF_h", ";Centrality (%);#LTv_{3}^{obs}#GT", nCentBins, centBins);
  TH1F* v3Obs_Sigma_PF_h = new TH1F("v3Obs_Sigma_PF_h", ";Centrality (%);#sigma(v_{3}^{obs})", nCentBins, centBins);
  TH1F* v3ObsCorr_Mean_PF_h = new TH1F("v3ObsCorr_Mean_PF_h", ";Centrality (%);#LTv_{3}^{obs}#GT", nCentBins, centBins);
  TH1F* v3ObsCorr_Sigma_PF_h = new TH1F("v3ObsCorr_Sigma_PF_h", ";Centrality (%);#sigma(v_{3}^{obs})", nCentBins, centBins);

  TH1F* systRaw_Mean_h[nSystR];
  TH1F* systRaw_Sigma_h[nSystR];
  TH1F* systFit_Mean_h[nSystR];
  TH1F* systFit_Sigma_h[nSystR];
  TH1F* systematic_h[nSystR];
  for (int i=0;i<nSystR;i++){
    systRaw_Mean_h[i] = new TH1F(("systRaw_Mean_" + std::to_string(10*systR[i]) + "_h").c_str(), ("systRaw Mean, R=" + std::to_string(systR[i])).c_str(), nCentBins, centBins);
    systRaw_Sigma_h[i] = new TH1F(("systRaw_Sigma_" + std::to_string(10*systR[i]) + "_h").c_str(), ("systRaw Sigma, R=" + std::to_string(systR[i])).c_str(), nCentBins, centBins);
    systFit_Mean_h[i] = new TH1F(("systFit_Mean_" + std::to_string(10*systR[i]) + "_h").c_str(), ("systFit Mean, R=" + std::to_string(systR[i])).c_str(), nCentBins, centBins);
    systFit_Sigma_h[i] = new TH1F(("systFit_Sigma_" + std::to_string(10*systR[i]) + "_h").c_str(), ("systFit Sigma, R=" + std::to_string(systR[i])).c_str(), nCentBins, centBins);
    systematic_h[i] = new TH1F(("systematic_" + std::to_string(10*systR[i]) + "_h").c_str(), ("systRaw Sigma-systFit Sigma, R=" + std::to_string(systR[i])).c_str(), nCentBins, centBins);
  }

  TH1F* size_h = new TH1F("size_h", "Vector size, all events", 100, 0, 2000);
  TH1F* size_rawPlus_h = new TH1F("size_rawPlus_h", "Vector size, raw-fit>0.2", 100, 0, 2000);
  TH1F* size_rawMinus_h = new TH1F("size_rawMinus_h", "Vector size, raw-fit<-0.2", 100, 0, 2000);
  TH1F* size_rawCorrPlus_h = new TH1F("size_rawCorrPlus_h", "Vector size, rawCorr-fit>0.2", 100, 0, 2000);
  TH1F* size_rawCorrMinus_h = new TH1F("size_rawCorrMinus_h", "Vector size, rawCorr-fit<-0.2", 100, 0 ,2000);

  TH1F* pt_h = new TH1F("pt_h", "pt, all events", 100, 0, 3.5);
  TH1F* pt_rawPlus_h = new TH1F("pt_rawPlus_h", "pt, raw-fit>0.2", 100, 0, 3.5);
  TH1F* pt_rawMinus_h = new TH1F("pt_rawMinus_h", "pt, raw-fit<-0.2", 100, 0, 3.5);
  TH1F* pt_rawCorrPlus_h = new TH1F("pt_rawCorrPlus_h", "pt, rawCorr-fit>0.2", 100, 0, 3.5);
  TH1F* pt_rawCorrMinus_h = new TH1F("pt_rawCorrMinus_h", "pt, rawCorr-fit<-0.2", 100, 0, 3.5);

  TH1F* phi_h = new TH1F("phi_h", "phi, all events", 100, -TMath::Pi(), TMath::Pi());
  TH1F* phi_rawPlus_h = new TH1F("phi_rawPlus_h", "phi, raw-fit>0.2", 100, -TMath::Pi(), TMath::Pi());
  TH1F* phi_rawMinus_h = new TH1F("phi_rawMinus_h", "phi, raw-fit<-0.2", 100, -TMath::Pi(), TMath::Pi());
  TH1F* phi_rawCorrPlus_h = new TH1F("phi_rawCorrPlus_h", "phi, rawCorr-fit>0.2", 100, -TMath::Pi(), TMath::Pi());
  TH1F* phi_rawCorrMinus_h = new TH1F("phi_rawCorrMinus_h", "phi, rawCorr-fit<-0.2", 100, -TMath::Pi(), TMath::Pi());

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    const std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
    v2Raw_PF_h[cI] = new TH1F(("v2Raw_" + centStr + "_PF_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2RawCorr_PF_h[cI] = new TH1F(("v2RawCorr_" + centStr + "_PF_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2Fit_PF_h[cI] = new TH1F(("v2Fit_" + centStr + "_PF_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2Diff_PF_h[cI] = new TH1F(("v2Diff_" + centStr + "_PF_h").c_str(), ("v2RawEByE-v2FitEByE (" + centStr + ");v_{2};Counts").c_str(), 150, -0.6, 0.6);
    v2DiffCorr_PF_h[cI] = new TH1F(("v2DiffCorr_" + centStr + "_PF_h").c_str(), ("v2RawEByECorr-v2FitEByE (" + centStr + ");v_{2};Counts").c_str(), 150, -0.6, 0.6);

    v2Obs_PF_h[cI] = new TH1F(("v2Obs_" + centStr + "_PF_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2ObsCorr_PF_h[cI] = new TH1F(("v2ObsCorr_" + centStr + "_PF_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);

    v2Raw_Trk_h[cI] = new TH1F(("v2Raw_" + centStr + "_Trk_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2RawCorr_Trk_h[cI] = new TH1F(("v2RawCorr_" + centStr + "_Trk_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2Obs_Trk_h[cI] = new TH1F(("v2Obs_" + centStr + "_Trk_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2ObsCorr_Trk_h[cI] = new TH1F(("v2ObsCorr_" + centStr + "_Trk_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    
    compRaw_h[cI] = new TH1F(("compRaw_" + centStr + "_h").c_str(), ";v_{2,trk}/v_{2,PF};Counts", 150, 0.0, 0.6);
    compRawCorr_h[cI] = new TH1F(("compRawCorr_" + centStr + "_h").c_str(), ";v_{2,trk}/v_{2,PF};Counts", 150, 0.0, 0.6);
    compObs_h[cI] = new TH1F(("compObs_" + centStr + "_h").c_str(), ";v_{2,trk}/v_{2,PF};Counts", 150, 0.0, 0.6);
    compObsCorr_h[cI] = new TH1F(("compObsCorr_" + centStr + "_h").c_str(), ";v_{2,trk}/v_{2,PF};Counts", 150, 0.0, 0.6);

    v3Raw_PF_h[cI] = new TH1F(("v3Raw_" + centStr + "_PF_h").c_str(), ";v_{3};Counts", 150, 0.0, 0.4);
    v3RawCorr_PF_h[cI] = new TH1F(("v3RawCorr_" + centStr + "_PF_h").c_str(), ";v_{3};Counts", 150, 0.0, 0.4);
    v3Fit_PF_h[cI] = new TH1F(("v3Fit_" + centStr + "_PF_h").c_str(), ";v_{3};Counts", 150, 0.0, 0.4);
    v3Diff_PF_h[cI] = new TH1F(("v3Diff_" + centStr + "_PF_h").c_str(), ("v3RawEByE-v3FitEByE (" + centStr + ");v_{3};Counts").c_str(), 150, -0.4, 0.4);
    v3DiffCorr_PF_h[cI] = new TH1F(("v3DiffCorr_" + centStr + "_PF_h").c_str(), ("v3RawEByECorr-v3FitEByE (" + centStr + ");v_{3};Counts").c_str(), 150, -0.4, 0.4);

    v3Obs_PF_h[cI] = new TH1F(("v3Obs_" + centStr + "_PF_h").c_str(), ";v_{3};Counts", 150, 0.0, 0.4);
    v3ObsCorr_PF_h[cI] = new TH1F(("v3ObsCorr_" + centStr + "_PF_h").c_str(), ";v_{3};Counts", 150, 0.0, 0.4);

    for (Int_t i = 0; i < nSystR; i++){
      systRaw_h[cI][i] = new TH1F(("systRaw_" + centStr + "_" + std::to_string(10*systR[i]) + "_h").c_str(), ("systRaw, R=" + std::to_string(systR[i])).c_str(), 150, -10, 10);
      systFit_h[cI][i] = new TH1F(("systFit_" + centStr + "_" + std::to_string(10*systR[i]) + "_h").c_str(), ("systFit, R=" + std::to_string(systR[i])).c_str(), 150, -10, 10);
    }
  }

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* inTree_p = (TTree*)inFile_p->Get("v2V3Tree");
  TObjArray* listOfBranches = (TObjArray*)inTree_p->GetListOfBranches();
  std::vector<std::string> vectorOfBranches;
  for(Int_t oI = 0; oI < listOfBranches->GetEntries(); ++oI){
    std::string branch = listOfBranches->At(oI)->GetName();
    vectorOfBranches.push_back(branch);
  }

  bool hasTrk = false;
  for(unsigned int tI = 0; tI < vectorOfBranches.size(); ++tI){
    if(vectorOfBranches.at(tI).find("trk") != std::string::npos){
      hasTrk = true;
      break;
    }
  }

  Int_t hiBin_;
  Float_t hiEvt2Plane_;
  Float_t hiEvt3Plane_;
  Float_t hiEvt4Plane_;
  Float_t v2FromTree_;  
  std::vector<float>* eByEPt_p=NULL;
  std::vector<float>* eByEEta_p=NULL;
  std::vector<float>* eByEPhi_p=NULL;
  std::vector<float>* eByEWeight_p=NULL;
  std::vector<float>* trkPt_p=NULL;
  std::vector<float>* trkPhi_p=NULL;
  std::vector<float>* trkWeight_p=NULL;  
  
  inTree_p->SetBranchAddress("hiBin", &hiBin_);
  inTree_p->SetBranchAddress("hiEvt2Plane", &hiEvt2Plane_);
  inTree_p->SetBranchAddress("hiEvt3Plane", &hiEvt3Plane_);
  inTree_p->SetBranchAddress("v2FromTree", &v2FromTree_);
  inTree_p->SetBranchAddress("eByEPt", &eByEPt_p);
  inTree_p->SetBranchAddress("eByEEta", &eByEEta_p);
  inTree_p->SetBranchAddress("eByEPhi", &eByEPhi_p);
  inTree_p->SetBranchAddress("eByEWeight", &eByEWeight_p);
  if(hasTrk){
    inTree_p->SetBranchAddress("trkPt", &trkPt_p);
    inTree_p->SetBranchAddress("trkPhi", &trkPhi_p);
    inTree_p->SetBranchAddress("trkWeight", &trkWeight_p);
  }

  const Int_t nEntries = TMath::Min((Int_t)inTree_p->GetEntries(), (Int_t)100000000);

  double globalN[nCentBins];
  double globalV2XRawPF[nCentBins];
  double globalV2YRawPF[nCentBins];
  double globalV2XRawPFCorr[nCentBins];
  double globalV2YRawPFCorr[nCentBins];

  double globalV2XRawTrk[nCentBins];
  double globalV2YRawTrk[nCentBins];
  double globalV2XRawTrkCorr[nCentBins];
  double globalV2YRawTrkCorr[nCentBins];

  double globalV3XRawPF[nCentBins];
  double globalV3YRawPF[nCentBins];
  double globalV3XRawPFCorr[nCentBins];
  double globalV3YRawPFCorr[nCentBins];

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    globalN[cI] = 0.;
    globalV2XRawPF[cI] = 0.;
    globalV2YRawPF[cI] = 0.;
    globalV2XRawPFCorr[cI] = 0.;
    globalV2YRawPFCorr[cI] = 0.;

    globalV2XRawTrk[cI] = 0.;
    globalV2YRawTrk[cI] = 0.;
    globalV2XRawTrkCorr[cI] = 0.;
    globalV2YRawTrkCorr[cI] = 0.;
  }

  std::string flowFitForm = "[0]*(1. + 2.*([1]*TMath::Cos(2.*(x - [2])) + [3]*TMath::Cos(3.*(x - [4])) + [5]*TMath::Cos(4.*(x - [6]))))";
  TF1* flowFit_p = new TF1("flowFit", flowFitForm.c_str(), -TMath::Pi(), TMath::Pi());

  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%10000 == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;
    inTree_p->GetEntry(entry);

    int centPos = -1;
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      if(centBinsLow[cI] <= hiBin_/2 && centBinsHi[cI] > hiBin_/2){
	centPos = cI;
	break;
      }
    }

    if(centPos < 0) continue;

    //    if (eByEPhi_p->size()<150) continue;

    double v2xRawPF = 0.;
    double v2yRawPF = 0.;

    double v2xRawPFCorr = 0.;
    double v2yRawPFCorr = 0.;

    double v2xRawTrk = 0.;
    double v2yRawTrk = 0.;

    double v2xRawTrkCorr = 0.;
    double v2yRawTrkCorr = 0.;

    double v3xRawPF = 0.;
    double v3yRawPF = 0.;

    double v3xRawPFCorr = 0.;
    double v3yRawPFCorr = 0.;

    double weightPF = 0.;
    double weightTrk = 0.;

    const Int_t nPhiBins = std::fmax(10, eByEPhi_p->size()/30);

    TH1F* tmp_phi_h = new TH1F("tempPhi", ";#phi;Track counts (.3 < p_{T} < 3.,#eta>0)", nPhiBins, -TMath::Pi(), TMath::Pi());
  
    double eventPlane4Cos = 0;
    double eventPlane4Sin = 0;

    for(unsigned pfI = 0; pfI < eByEPhi_p->size(); pfI++){
      double tempWeightPF = eByEWeight_p->at(pfI);
      //      double deltaEventPhi = eByEPhi_p->at(pfI) - hiEvt2Plane_;
      double deltaEventPhi = eByEPhi_p->at(pfI);

      v2xRawPF += TMath::Cos(2*(deltaEventPhi));
      v2yRawPF += TMath::Sin(2*(deltaEventPhi));

      v3xRawPF += TMath::Cos(3*deltaEventPhi);
      v3yRawPF += TMath::Sin(3*deltaEventPhi);
	
      v2xRawPFCorr += tempWeightPF*TMath::Cos(2*(deltaEventPhi));
      v2yRawPFCorr += tempWeightPF*TMath::Sin(2*(deltaEventPhi));

      v3xRawPFCorr += tempWeightPF*TMath::Cos(3*deltaEventPhi);
      v3yRawPFCorr += tempWeightPF*TMath::Sin(3*deltaEventPhi);

      weightPF += tempWeightPF;

      tmp_phi_h->Fill(eByEPhi_p->at(pfI));

      eventPlane4Cos += std::cos(4*eByEPhi_p->at(pfI));
      eventPlane4Sin += std::sin(4*eByEPhi_p->at(pfI));
    }

    hiEvt4Plane_ = std::atan2(eventPlane4Sin, eventPlane4Cos)/4;

    //    std::string flowFitForm = "[0]*(1. + 2.*([1]*TMath::Cos(2.*(x - " + std::to_string(hiEvt2Plane_) + ")) + [2]*TMath::Cos(3.*(x - " + std::to_string(hiEvt3Plane_) + "))))";
    //    TF1* flowFit_p = new TF1("flowFit", flowFitForm.c_str(), -TMath::Pi(), TMath::Pi());
    flowFit_p->SetParameter(0, 10);
    flowFit_p->SetParameter(1, 0);
    flowFit_p->FixParameter(2, hiEvt2Plane_);
    flowFit_p->SetParameter(3, 0);
    flowFit_p->FixParameter(4, hiEvt3Plane_);
    flowFit_p->SetParameter(5, 0);
    flowFit_p->FixParameter(6, hiEvt4Plane_);

    tmp_phi_h->Fit(flowFit_p, "Q0", "", -TMath::Pi(), TMath::Pi());

    double v2FitPF = flowFit_p->GetParameter(1);
    double v3FitPF = flowFit_p->GetParameter(3);

    delete tmp_phi_h;
    //    delete flowFit_p;

    Int_t trkCounter = 0;
    if(hasTrk){    
      for(unsigned tI = 0; tI < trkPhi_p->size(); tI++){
	//	if(trkPt_p->at(tI) >= 2.4) continue;
	trkCounter++;
	double tempWeightTrk = trkWeight_p->at(tI);
	//	double deltaEventPhi = trkPhi_p->at(tI) - hiEvt2Plane_;
	double deltaEventPhi = trkPhi_p->at(tI);
	
	v2xRawTrk += TMath::Cos(2*(deltaEventPhi));
	v2yRawTrk += TMath::Sin(2*(deltaEventPhi));

	v2xRawTrkCorr += tempWeightTrk*TMath::Cos(2*(deltaEventPhi));
	v2yRawTrkCorr += tempWeightTrk*TMath::Sin(2*(deltaEventPhi));
	
	weightTrk += tempWeightTrk;
      }
    }

    v2xRawPF /= (double)eByEPhi_p->size();
    v2yRawPF /= (double)eByEPhi_p->size();

    v2xRawPFCorr /= weightPF;
    v2yRawPFCorr /= weightPF;

    v3xRawPF /= (double)eByEPhi_p->size();
    v3yRawPF /= (double)eByEPhi_p->size();

    v3xRawPFCorr /= weightPF;
    v3yRawPFCorr /= weightPF;

    double v2RawPF = TMath::Sqrt(v2xRawPF*v2xRawPF + v2yRawPF*v2yRawPF);
    double v2RawPFCorr = TMath::Sqrt(v2xRawPFCorr*v2xRawPFCorr + v2yRawPFCorr*v2yRawPFCorr);

    double v3RawPF = TMath::Sqrt(v3xRawPF*v3xRawPF + v3yRawPF*v3yRawPF);
    double v3RawPFCorr = TMath::Sqrt(v3xRawPFCorr*v3xRawPFCorr + v3yRawPFCorr*v3yRawPFCorr);

    v2Raw_PF_h[centPos]->Fill(v2RawPF);
    v2RawCorr_PF_h[centPos]->Fill(v2RawPFCorr);
    v2Fit_PF_h[centPos]->Fill(v2FitPF);
    v2Diff_PF_h[centPos]->Fill(v2RawPF-v2FitPF);
    v2DiffCorr_PF_h[centPos]->Fill(v2RawPFCorr-v2FitPF);

    v3Raw_PF_h[centPos]->Fill(v3RawPF);
    v3RawCorr_PF_h[centPos]->Fill(v3RawPFCorr);
    v3Fit_PF_h[centPos]->Fill(v3FitPF);
    v3Diff_PF_h[centPos]->Fill(v3RawPF-v3FitPF);
    v3DiffCorr_PF_h[centPos]->Fill(v3RawPFCorr-v3FitPF);

    size_h->Fill(eByEPhi_p->size());
    for (unsigned pfI = 0; pfI < eByEPhi_p->size(); pfI++){
      pt_h->Fill(eByEPt_p->at(pfI));
      phi_h->Fill(eByEPhi_p->at(pfI));
    }
    if (v2RawPF-v2FitPF>0.2){
      size_rawPlus_h->Fill(eByEPhi_p->size());
      for (unsigned pfI = 0; pfI < eByEPhi_p->size(); pfI++){
	pt_rawPlus_h->Fill(eByEPt_p->at(pfI));
	phi_rawPlus_h->Fill(eByEPhi_p->at(pfI));
      }
    }
    if (v2RawPF-v2FitPF<-0.2){
      size_rawMinus_h->Fill(eByEPhi_p->size());
      for (unsigned pfI = 0; pfI < eByEPhi_p->size(); pfI++){
        pt_rawMinus_h->Fill(eByEPt_p->at(pfI));
        phi_rawMinus_h->Fill(eByEPhi_p->at(pfI));
      }
    }
    if (v2RawPFCorr-v2FitPF>0.2){
      size_rawCorrPlus_h->Fill(eByEPhi_p->size());
      for (unsigned pfI = 0; pfI < eByEPhi_p->size(); pfI++){
        pt_rawCorrPlus_h->Fill(eByEPt_p->at(pfI));
        phi_rawCorrPlus_h->Fill(eByEPhi_p->at(pfI));
      }
    }
    if (v2RawPFCorr-v2FitPF<-0.2){
      size_rawCorrMinus_h->Fill(eByEPhi_p->size());
      for (unsigned pfI = 0; pfI < eByEPhi_p->size(); pfI++){
        pt_rawCorrMinus_h->Fill(eByEPt_p->at(pfI));
        phi_rawCorrMinus_h->Fill(eByEPhi_p->at(pfI));
      }
    }

    globalN[centPos] += 1;
    globalV2XRawPF[centPos] += v2xRawPF;
    globalV2YRawPF[centPos] += v2yRawPF;
    globalV2XRawPFCorr[centPos] += v2xRawPFCorr;
    globalV2YRawPFCorr[centPos] += v2yRawPFCorr;

    globalV3XRawPF[centPos] += v3xRawPF;
    globalV3YRawPF[centPos] += v3yRawPF;
    globalV3XRawPFCorr[centPos] += v3xRawPFCorr;
    globalV3YRawPFCorr[centPos] += v3yRawPFCorr;

    if(hasTrk){
      v2xRawTrk /= (double)trkCounter;
      v2yRawTrk /= (double)trkCounter;
      
      v2xRawTrkCorr /= weightTrk;
      v2yRawTrkCorr /= weightTrk;
      
      double v2RawTrk = TMath::Sqrt(v2xRawTrk*v2xRawTrk + v2yRawTrk*v2yRawTrk);
      double v2RawTrkCorr = TMath::Sqrt(v2xRawTrkCorr*v2xRawTrkCorr + v2yRawTrkCorr*v2yRawTrkCorr);
      
      v2Raw_Trk_h[centPos]->Fill(v2RawTrk);
      v2RawCorr_Trk_h[centPos]->Fill(v2RawTrkCorr);

      compRaw_h[centPos]->Fill(v2RawTrk/v2RawPF);
      compRawCorr_h[centPos]->Fill(v2RawTrkCorr/v2RawPFCorr);
   
      globalN[centPos] += 1;
      globalV2XRawTrk[centPos] += v2xRawTrk;
      globalV2YRawTrk[centPos] += v2yRawTrk;
      globalV2XRawTrkCorr[centPos] += v2xRawTrkCorr;
      globalV2YRawTrkCorr[centPos] += v2yRawTrkCorr;
    }

    for(Int_t sI = 0; sI < nSystR; sI++){
      Double_t R = nSystR[sI];

      // Random eta, phi position
      TRandom3* rand = new TRandom3();

      Double_t eta_c = rand->Rndm()*2*(1.-R)-(1.-R);
      Double_t phi_c = rand->Rndm()*2*TMath::Pi()-TMath::Pi();

      // Calculate sumpt
      Double_t sumpt = 0.;
      for (int pfI = 0; pfI < eByEPt_p->size(); pfI++){
	if (dR(eta_c, phi_c, eByEEta_p->at(pfI), eByEPhi_p->at(pfI), R)) sumpt += eByEPt_p->at(pfI);
      }

      // Fun grid
      Double_t stripSize = R/nStrips;
      Double_t rhoSumRaw = 0.;
      Double_t rhoSumFit = 0.;

      for (Double_t eta = eta_c - R + stripSize/2; eta < eta_c + R; eta += stripSize){
	// Retrieve rho
	Double_t rho_ = 0.;

	// Calculate phi-modulated rho sum
	for (Double_t phi = phi_c - R + stripSize/2; phi < phi_c + R; phi += stripSize){
	  if (dR(eta_c, phi_c, eta, phi, R)){
	    rhoSumRaw += stripSize*stripSize * rho_ * (1+2*(v2RawPF*TMath::Cos(2*(phi-hiEvt2Plane_))+v3RawPF*TMath::Cos(3*(phi-hiEvt3Plane_))));
	    rhoSumFit += stripSize*stripSize * rho_ * (1+2*(v2FitPF*TMath::Cos(2*(phi-hiEvt2Plane_))+v3FitPF*TMath::Cos(3*(phi-hiEvt3Plane_))));
	  }
	}
      }

      // Fill histograms
      systRaw_h[cI][sI]->Fill(sumpt-rhoSumRaw);
      systFit_h[cI][sI]->Fill(sumpt-rhoSumFit);
    }

  }

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    globalV2XRawPF[cI] /= globalN[cI];
    globalV2YRawPF[cI] /= globalN[cI];
    globalV2XRawPFCorr[cI] /= globalN[cI];
    globalV2YRawPFCorr[cI] /= globalN[cI];

    globalV3XRawPF[cI] /= globalN[cI];
    globalV3YRawPF[cI] /= globalN[cI];
    globalV3XRawPFCorr[cI] /= globalN[cI];
    globalV3YRawPFCorr[cI] /= globalN[cI];

    if(hasTrk){
      globalV2XRawTrk[cI] /= globalN[cI];
      globalV2YRawTrk[cI] /= globalN[cI];
      globalV2XRawTrkCorr[cI] /= globalN[cI];
      globalV2YRawTrkCorr[cI] /= globalN[cI];
    }

    std::cout << "Cent, N: " << centBinsLow[cI] << "-" << centBinsHi[cI] << "%, " << globalN[cI] << std::endl;
    std::cout << "  Global v2xRawPF, v2yRawPF, N: " << globalV2XRawPF[cI] << ", " << globalV2YRawPF[cI] << std::endl;
    std::cout << "  Global v2xRawPFCorr, v2yRawPFCorr, N: " << globalV2XRawPFCorr[cI] << ", " << globalV2YRawPFCorr[cI] << std::endl;

    std::cout << "  Global v3xRawPF, v3yRawPF, N: " << globalV3XRawPF[cI] << ", " << globalV3YRawPF[cI] << std::endl;
    std::cout << "  Global v3xRawPFCorr, v3yRawPFCorr, N: " << globalV3XRawPFCorr[cI] << ", " << globalV3YRawPFCorr[cI] << std::endl;

    if(hasTrk){  
      std::cout << "  Global v2xRawTrk, v2yRawTrk, N: " << globalV2XRawTrk[cI] << ", " << globalV2YRawTrk[cI] << std::endl;
      std::cout << "  Global v2xRawTrkCorr, v2yRawTrkCorr, N: " << globalV2XRawTrkCorr[cI] << ", " << globalV2YRawTrkCorr[cI] << std::endl;
    }
  }



  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%10000 == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;
    inTree_p->GetEntry(entry);

    int centPos = -1;
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      if(centBinsLow[cI] <= hiBin_/2 && centBinsHi[cI] > hiBin_/2){
	centPos = cI;
	break;
      }
    }

    if(centPos < 0) continue;

    double v2xObsPF = 0.;
    double v2yObsPF = 0.;

    double v2xObsPFCorr = 0.;
    double v2yObsPFCorr = 0.;

    double v3xObsPF = 0.;
    double v3yObsPF = 0.;

    double v3xObsPFCorr = 0.;
    double v3yObsPFCorr = 0.;

    double weightPF = 0.;
  
    for(unsigned pfI = 0; pfI < eByEPhi_p->size(); pfI++){
      double tempWeightPF = eByEWeight_p->at(pfI);
      //double deltaEventPhi = eByEPhi_p->at(pfI) - hiEvt2Plane_;
      double deltaEventPhi = eByEPhi_p->at(pfI);

      v2xObsPF += TMath::Cos(2*(deltaEventPhi));
      v2yObsPF += TMath::Sin(2*(deltaEventPhi));

      v2xObsPFCorr += tempWeightPF*TMath::Cos(2*(deltaEventPhi));
      v2yObsPFCorr += tempWeightPF*TMath::Sin(2*(deltaEventPhi));

      v3xObsPF += TMath::Cos(3*(deltaEventPhi));
      v3yObsPF += TMath::Sin(3*(deltaEventPhi));

      v3xObsPFCorr += tempWeightPF*TMath::Cos(3*(deltaEventPhi));
      v3yObsPFCorr += tempWeightPF*TMath::Sin(3*(deltaEventPhi));

      weightPF += tempWeightPF;
    }


    double v2xObsTrk = 0.;
    double v2yObsTrk = 0.;

    double v2xObsTrkCorr = 0.;
    double v2yObsTrkCorr = 0.;

    double weightTrk = 0.;

    Int_t trkCounter = 0;
    if(hasTrk){
      for(unsigned tI = 0; tI < trkPhi_p->size(); tI++){
	//	if(trkPt_p->at(tI) >= 2.4) continue;
	trkCounter++;
	double tempWeightTrk = trkWeight_p->at(tI);
	//	double deltaEventPhi = trkPhi_p->at(tI) - hiEvt2Plane_;
	double deltaEventPhi = trkPhi_p->at(tI);
	
	v2xObsTrk += TMath::Cos(2*(deltaEventPhi));
	v2yObsTrk += TMath::Sin(2*(deltaEventPhi));
	
	v2xObsTrkCorr += tempWeightTrk*TMath::Cos(2*(deltaEventPhi));
	v2yObsTrkCorr += tempWeightTrk*TMath::Sin(2*(deltaEventPhi));
	
	weightTrk += tempWeightTrk;
      }
    }

    v2xObsPF /= (double)eByEPhi_p->size();
    v2yObsPF /= (double)eByEPhi_p->size();

    v2xObsPFCorr /= weightPF;
    v2yObsPFCorr /= weightPF;

    v2xObsPF -= globalV2XRawPF[centPos];
    v2yObsPF -= globalV2YRawPF[centPos];

    v2xObsPFCorr -= globalV2XRawPFCorr[centPos];
    v2yObsPFCorr -= globalV2YRawPFCorr[centPos];

    v3xObsPF /= (double)eByEPhi_p->size();
    v3yObsPF /= (double)eByEPhi_p->size();

    v3xObsPFCorr /= weightPF;
    v3yObsPFCorr /= weightPF;

    v3xObsPF -= globalV3XRawPF[centPos];
    v3yObsPF -= globalV3YRawPF[centPos];

    v3xObsPFCorr -= globalV3XRawPFCorr[centPos];
    v3yObsPFCorr -= globalV3YRawPFCorr[centPos];

    double v2ObsPF = TMath::Sqrt(v2xObsPF*v2xObsPF + v2yObsPF*v2yObsPF);
    double v2ObsPFCorr = TMath::Sqrt(v2xObsPFCorr*v2xObsPFCorr + v2yObsPFCorr*v2yObsPFCorr);

    double v3ObsPF = TMath::Sqrt(v3xObsPF*v3xObsPF + v3yObsPF*v3yObsPF);
    double v3ObsPFCorr = TMath::Sqrt(v3xObsPFCorr*v3xObsPFCorr + v3yObsPFCorr*v3yObsPFCorr);

    v2Obs_PF_h[centPos]->Fill(v2ObsPF);
    v2ObsCorr_PF_h[centPos]->Fill(v2ObsPFCorr);

    v3Obs_PF_h[centPos]->Fill(v3ObsPF);
    v3ObsCorr_PF_h[centPos]->Fill(v3ObsPFCorr);


    if(hasTrk){
      v2xObsTrk /= (double)trkCounter;
      v2yObsTrk /= (double)trkCounter;
      
      v2xObsTrkCorr /= weightTrk;
      v2yObsTrkCorr /= weightTrk;
      
      v2xObsTrk -= globalV2XRawTrk[centPos];
      v2yObsTrk -= globalV2YRawTrk[centPos];
      
      v2xObsTrkCorr -= globalV2XRawTrkCorr[centPos];
      v2yObsTrkCorr -= globalV2YRawTrkCorr[centPos];
      
      double v2ObsTrk = TMath::Sqrt(v2xObsTrk*v2xObsTrk + v2yObsTrk*v2yObsTrk);
      double v2ObsTrkCorr = TMath::Sqrt(v2xObsTrkCorr*v2xObsTrkCorr + v2yObsTrkCorr*v2yObsTrkCorr);
      
      v2Obs_Trk_h[centPos]->Fill(v2ObsTrk);
      v2ObsCorr_Trk_h[centPos]->Fill(v2ObsTrkCorr);

      compObs_h[centPos]->Fill(v2ObsTrk/v2ObsPF);
      compObsCorr_h[centPos]->Fill(v2ObsTrkCorr/v2ObsPFCorr);
    }
  }


  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  TF1* tmp_gaus=new TF1("tmp_gaus","gaus");

  for(Int_t cI = 0; cI < nCentBins; ++cI){  
    v2Raw_PF_h[cI]->Write("", TObject::kOverwrite);
    v2RawCorr_PF_h[cI]->Write("", TObject::kOverwrite);
    v2Fit_PF_h[cI]->Write("", TObject::kOverwrite);
    v2Diff_PF_h[cI]->Write("", TObject::kOverwrite);
    v2DiffCorr_PF_h[cI]->Write("", TObject::kOverwrite);

    v2Obs_PF_h[cI]->Write("", TObject::kOverwrite);
    v2ObsCorr_PF_h[cI]->Write("", TObject::kOverwrite);

    v3Raw_PF_h[cI]->Write("", TObject::kOverwrite);
    v3RawCorr_PF_h[cI]->Write("", TObject::kOverwrite);
    v3Fit_PF_h[cI]->Write("", TObject::kOverwrite);
    v3Diff_PF_h[cI]->Write("", TObject::kOverwrite);
    v3DiffCorr_PF_h[cI]->Write("", TObject::kOverwrite);

    v3Obs_PF_h[cI]->Write("", TObject::kOverwrite);
    v3ObsCorr_PF_h[cI]->Write("", TObject::kOverwrite);

    v2Raw_PF_h[cI]->Fit(tmp_gaus,"QN0");
    TF1* v2Raw_PF_f = new TF1("v2Raw_PF_f","gaus",tmp_gaus->GetParameter(1)-2*tmp_gaus->GetParameter(2), tmp_gaus->GetParameter(1)+2*tmp_gaus->GetParameter(2));
    v2Raw_PF_h[cI]->Fit(v2Raw_PF_f,"QNR0");

    v2Raw_Mean_PF_h->SetBinContent(cI+1, v2Raw_PF_f->GetParameter(1));
    v2Raw_Mean_PF_h->SetBinError(cI+1, v2Raw_PF_f->GetParError(1));
    v2Raw_Sigma_PF_h->SetBinContent(cI+1, v2Raw_PF_f->GetParameter(2));
    v2Raw_Sigma_PF_h->SetBinError(cI+1, v2Raw_PF_f->GetParError(2));

    v2RawCorr_PF_h[cI]->Fit(tmp_gaus,"QN0");
    TF1* v2RawCorr_PF_f = new TF1("v2RawCorr_PF_f","gaus",tmp_gaus->GetParameter(1)-2*tmp_gaus->GetParameter(2), tmp_gaus->GetParameter(1)+2*tmp_gaus->GetParameter(2));
    v2RawCorr_PF_h[cI]->Fit(v2RawCorr_PF_f,"QNR0");

    v2RawCorr_Mean_PF_h->SetBinContent(cI+1, v2RawCorr_PF_f->GetParameter(1));
    v2RawCorr_Mean_PF_h->SetBinError(cI+1, v2RawCorr_PF_f->GetParError(1));
    v2RawCorr_Sigma_PF_h->SetBinContent(cI+1, v2RawCorr_PF_f->GetParameter(2));
    v2RawCorr_Sigma_PF_h->SetBinError(cI+1, v2RawCorr_PF_f->GetParError(2));
    
    v2Fit_PF_h[cI]->Fit(tmp_gaus,"QN0");
    TF1* v2Fit_PF_f = new TF1("v2Fit_PF_f","gaus",tmp_gaus->GetParameter(1)-2*tmp_gaus->GetParameter(2), tmp_gaus->GetParameter(1)+2*tmp_gaus->GetParameter(2));
    v2Fit_PF_h[cI]->Fit(v2Fit_PF_f,"QNR0");

    v2Fit_Mean_PF_h->SetBinContent(cI+1, v2Fit_PF_f->GetParameter(1));
    v2Fit_Mean_PF_h->SetBinError(cI+1, v2Fit_PF_f->GetParError(1));
    v2Fit_Sigma_PF_h->SetBinContent(cI+1, v2Fit_PF_f->GetParameter(2));
    v2Fit_Sigma_PF_h->SetBinError(cI+1, v2Fit_PF_f->GetParError(2));

    v2Obs_PF_h[cI]->Fit(tmp_gaus,"QN0");
    TF1* v2Obs_PF_f = new TF1("v2Obs_PF_f","gaus",tmp_gaus->GetParameter(1)-2*tmp_gaus->GetParameter(2), tmp_gaus->GetParameter(1)+2*tmp_gaus->GetParameter(2));
    v2Obs_PF_h[cI]->Fit(v2Obs_PF_f,"QNR0");

    v2Obs_Mean_PF_h->SetBinContent(cI+1, v2Obs_PF_f->GetParameter(1));
    v2Obs_Mean_PF_h->SetBinError(cI+1, v2Obs_PF_f->GetParError(1));
    v2Obs_Sigma_PF_h->SetBinContent(cI+1, v2Obs_PF_f->GetParameter(2));
    v2Obs_Sigma_PF_h->SetBinError(cI+1, v2Obs_PF_f->GetParError(2));

    v2ObsCorr_PF_h[cI]->Fit(tmp_gaus,"QN0");
    TF1* v2ObsCorr_PF_f = new TF1("v2ObsCorr_PF_f","gaus",tmp_gaus->GetParameter(1)-2*tmp_gaus->GetParameter(2), tmp_gaus->GetParameter(1)+2*tmp_gaus->GetParameter(2));
    v2ObsCorr_PF_h[cI]->Fit(v2ObsCorr_PF_f,"QNR0");

    v2ObsCorr_Mean_PF_h->SetBinContent(cI+1, v2ObsCorr_PF_f->GetParameter(1));
    v2ObsCorr_Mean_PF_h->SetBinError(cI+1, v2ObsCorr_PF_f->GetParError(1));
    v2ObsCorr_Sigma_PF_h->SetBinContent(cI+1, v2ObsCorr_PF_f->GetParameter(2));
    v2ObsCorr_Sigma_PF_h->SetBinError(cI+1, v2ObsCorr_PF_f->GetParError(2));

    v3Raw_PF_h[cI]->Fit(tmp_gaus,"QN0");
    TF1* v3Raw_PF_f = new TF1("v3Raw_PF_f","gaus",tmp_gaus->GetParameter(1)-2*tmp_gaus->GetParameter(2), tmp_gaus->GetParameter(1)+2*tmp_gaus->GetParameter(2));
    v3Raw_PF_h[cI]->Fit(v3Raw_PF_f,"QNR0");

    v3Raw_Mean_PF_h->SetBinContent(cI+1, v3Raw_PF_f->GetParameter(1));
    v3Raw_Mean_PF_h->SetBinError(cI+1, v3Raw_PF_f->GetParError(1));
    v3Raw_Sigma_PF_h->SetBinContent(cI+1, v3Raw_PF_f->GetParameter(2));
    v3Raw_Sigma_PF_h->SetBinError(cI+1, v3Raw_PF_f->GetParError(2));

    v3RawCorr_PF_h[cI]->Fit(tmp_gaus,"QN0");
    TF1* v3RawCorr_PF_f = new TF1("v3RawCorr_PF_f","gaus",tmp_gaus->GetParameter(1)-2*tmp_gaus->GetParameter(2), tmp_gaus->GetParameter(1)+2*tmp_gaus->GetParameter(2));
    v3RawCorr_PF_h[cI]->Fit(v3RawCorr_PF_f,"QNR0");

    v3RawCorr_Mean_PF_h->SetBinContent(cI+1, v3RawCorr_PF_f->GetParameter(1));
    v3RawCorr_Mean_PF_h->SetBinError(cI+1, v3RawCorr_PF_f->GetParError(1));
    v3RawCorr_Sigma_PF_h->SetBinContent(cI+1, v3RawCorr_PF_f->GetParameter(2));
    v3RawCorr_Sigma_PF_h->SetBinError(cI+1, v3RawCorr_PF_f->GetParError(2));

    v3Fit_PF_h[cI]->Fit(tmp_gaus,"QN0");
    TF1* v3Fit_PF_f = new TF1("v3Fit_PF_f","gaus",tmp_gaus->GetParameter(1)-2*tmp_gaus->GetParameter(2), tmp_gaus->GetParameter(1)+2*tmp_gaus->GetParameter(2));
    v3Fit_PF_h[cI]->Fit(v3Fit_PF_f,"QNR0");

    v3Fit_Mean_PF_h->SetBinContent(cI+1, v3Fit_PF_f->GetParameter(1));
    v3Fit_Mean_PF_h->SetBinError(cI+1, v3Fit_PF_f->GetParError(1));
    v3Fit_Sigma_PF_h->SetBinContent(cI+1, v3Fit_PF_f->GetParameter(2));
    v3Fit_Sigma_PF_h->SetBinError(cI+1, v3Fit_PF_f->GetParError(2));

    v3Obs_PF_h[cI]->Fit(tmp_gaus,"QN0");
    TF1* v3Obs_PF_f = new TF1("v3Obs_PF_f","gaus",tmp_gaus->GetParameter(1)-2*tmp_gaus->GetParameter(2), tmp_gaus->GetParameter(1)+2*tmp_gaus->GetParameter(2));
    v3Obs_PF_h[cI]->Fit(v3Obs_PF_f,"QNR0");

    v3Obs_Mean_PF_h->SetBinContent(cI+1, v3Obs_PF_f->GetParameter(1));
    v3Obs_Mean_PF_h->SetBinError(cI+1, v3Obs_PF_f->GetParError(1));
    v3Obs_Sigma_PF_h->SetBinContent(cI+1, v3Obs_PF_f->GetParameter(2));
    v3Obs_Sigma_PF_h->SetBinError(cI+1, v3Obs_PF_f->GetParError(2));

    v3ObsCorr_PF_h[cI]->Fit(tmp_gaus,"QN0");
    TF1* v3ObsCorr_PF_f = new TF1("v3ObsCorr_PF_f","gaus",tmp_gaus->GetParameter(1)-2*tmp_gaus->GetParameter(2), tmp_gaus->GetParameter(1)+2*tmp_gaus->GetParameter(2));
    v3ObsCorr_PF_h[cI]->Fit(v3ObsCorr_PF_f,"QNR0");

    v3ObsCorr_Mean_PF_h->SetBinContent(cI+1, v3ObsCorr_PF_f->GetParameter(1));
    v3ObsCorr_Mean_PF_h->SetBinError(cI+1, v3ObsCorr_PF_f->GetParError(1));
    v3ObsCorr_Sigma_PF_h->SetBinContent(cI+1, v3ObsCorr_PF_f->GetParameter(2));
    v3ObsCorr_Sigma_PF_h->SetBinError(cI+1, v3ObsCorr_PF_f->GetParError(2));

    delete v2Raw_PF_h[cI];
    delete v2RawCorr_PF_h[cI];
    delete v2Fit_PF_h[cI];
    delete v2Diff_PF_h[cI];
    delete v2DiffCorr_PF_h[cI];

    delete v2Obs_PF_h[cI];
    delete v2ObsCorr_PF_h[cI];

    delete v2Raw_PF_f;
    delete v2RawCorr_PF_f;
    delete v2Fit_PF_f;
    delete v2Obs_PF_f;
    delete v2ObsCorr_PF_f;

    delete v3Raw_PF_h[cI];
    delete v3RawCorr_PF_h[cI];
    delete v3Fit_PF_h[cI];
    delete v3Diff_PF_h[cI];
    delete v3DiffCorr_PF_h[cI];

    delete v3Obs_PF_h[cI];
    delete v3ObsCorr_PF_h[cI];

    delete v3Raw_PF_f;
    delete v3RawCorr_PF_f;
    delete v3Fit_PF_f;
    delete v3Obs_PF_f;
    delete v3ObsCorr_PF_f;

    for(int i = 0; i < nSystR; i++){
      systRaw_h[cI][i]->Write("", TObject::kOverwrite);
      systFit_h[cI][i]->Write("", TObject::kOverwrite);

      systRaw_h[cI][i]->Fit(tmp_gaus,"QN0");
      TF1* systRaw_f = new TF1("systRaw_f","gaus",tmp_gaus->GetParameter(1)-2*tmp_gaus->GetParameter(2), tmp_gaus->GetParameter(1)+tmp_gaus->GetParameter(2));
      systRaw_h[cI][i]->Fit(systRaw_f,"QNR0");

      systRaw_Mean_h[i]->SetBinContent(cI+1, systRaw_f->GetParameter(1));
      systRaw_Mean_h[i]->SetBinError(cI+1, systRaw_f->GetParError(1));
      systRaw_Sigma_h[i]->SetBinContent(cI+1, systRaw_f->GetParameter(2));
      systRaw_Sigma_h[i]->SetBinError(cI+1, systRaw_f->GetParError(2));

      systFit_h[cI][i]->Fit(tmp_gaus,"QN0");
      TF1* systFit_f = new TF1("systFit_f","gaus",tmp_gaus->GetParameter(1)-2*tmp_gaus->GetParameter(2), tmp_gaus->GetParameter(1)+tmp_gaus->GetParameter(2));
      systFit_h[cI][i]->Fit(systFit_f,"QNR0");

      systFit_Mean_h[i]->SetBinContent(cI+1, systFit_f->GetParameter(1));
      systFit_Mean_h[i]->SetBinError(cI+1, systFit_f->GetParError(1));
      systFit_Sigma_h[i]->SetBinContent(cI+1, systFit_f->GetParameter(2));
      systFit_Sigma_h[i]->SetBinError(cI+1, systFit_f->GetParError(2));

      delete systRaw_h[cI][i];
      delete systFit_h[cI][i];

      delete systRaw_f;
      delete systFit_f;
    }
  }

  v2Raw_Mean_PF_h->Write("", TObject::kOverwrite);
  v2Raw_Sigma_PF_h->Write("", TObject::kOverwrite);

  v2RawCorr_Mean_PF_h->Write("", TObject::kOverwrite);
  v2RawCorr_Sigma_PF_h->Write("", TObject::kOverwrite);

  v2Fit_Mean_PF_h->Write("", TObject::kOverwrite);
  v2Fit_Sigma_PF_h->Write("", TObject::kOverwrite);

  v2Obs_Mean_PF_h->Write("", TObject::kOverwrite);
  v2Obs_Sigma_PF_h->Write("", TObject::kOverwrite);

  v2ObsCorr_Mean_PF_h->Write("", TObject::kOverwrite);
  v2ObsCorr_Sigma_PF_h->Write("", TObject::kOverwrite);

  v3Raw_Mean_PF_h->Write("", TObject::kOverwrite);
  v3Raw_Sigma_PF_h->Write("", TObject::kOverwrite);

  v3RawCorr_Mean_PF_h->Write("", TObject::kOverwrite);
  v3RawCorr_Sigma_PF_h->Write("", TObject::kOverwrite);

  v3Fit_Mean_PF_h->Write("", TObject::kOverwrite);
  v3Fit_Sigma_PF_h->Write("", TObject::kOverwrite);

  v3Obs_Mean_PF_h->Write("", TObject::kOverwrite);
  v3Obs_Sigma_PF_h->Write("", TObject::kOverwrite);

  v3ObsCorr_Mean_PF_h->Write("", TObject::kOverwrite);
  v3ObsCorr_Sigma_PF_h->Write("", TObject::kOverwrite);

  delete v2Raw_Mean_PF_h;
  delete v2Raw_Sigma_PF_h;
  delete v2RawCorr_Mean_PF_h;
  delete v2RawCorr_Sigma_PF_h;
  delete v2Fit_Mean_PF_h;
  delete v2Fit_Sigma_PF_h;

  delete v2Obs_Mean_PF_h;
  delete v2Obs_Sigma_PF_h;
  delete v2ObsCorr_Mean_PF_h;
  delete v2ObsCorr_Sigma_PF_h;

  delete v3Raw_Mean_PF_h;
  delete v3Raw_Sigma_PF_h;
  delete v3RawCorr_Mean_PF_h;
  delete v3RawCorr_Sigma_PF_h;
  delete v3Fit_Mean_PF_h;
  delete v3Fit_Sigma_PF_h;

  delete v3Obs_Mean_PF_h;
  delete v3Obs_Sigma_PF_h;
  delete v3ObsCorr_Mean_PF_h;
  delete v3ObsCorr_Sigma_PF_h;

  for (int i=0; i < nSystR; i++){
    systRaw_Mean_h[i]->Write("", TObject::kOverwrite);
    systRaw_Sigma_h[i]->Write("", TObject::kOverwrite);
    
    systFit_Mean_h[i]->Write("", TObject::kOverwrite);
    systFit_Sigma_h[i]->Write("", TObject::kOverwrite);

    systematic_h[i]->Add(systRaw_Sigma_h[i], systFit_Sigma_h[i], 1, -1);
    systematic_h[i]->Write("", TObject::kOverwrite);

    delete systRaw_Mean_h[i];
    delete systRaw_Sigma_h[i];
    delete systFit_Mean_h[i];
    delete systFit_Sigma_h[i];
    delete systematic_h[i];
  }

  size_h->Write("", TObject::kOverwrite);
  size_rawPlus_h->Write("", TObject::kOverwrite);
  size_rawMinus_h->Write("", TObject::kOverwrite);
  size_rawCorrPlus_h->Write("", TObject::kOverwrite);
  size_rawCorrMinus_h->Write("", TObject::kOverwrite);

  pt_h->Write("", TObject::kOverwrite);
  pt_rawPlus_h->Write("", TObject::kOverwrite);
  pt_rawMinus_h->Write("", TObject::kOverwrite);
  pt_rawCorrPlus_h->Write("", TObject::kOverwrite);
  pt_rawCorrMinus_h->Write("", TObject::kOverwrite);

  phi_h->Write("", TObject::kOverwrite);
  phi_rawPlus_h->Write("", TObject::kOverwrite);
  phi_rawMinus_h->Write("", TObject::kOverwrite);
  phi_rawCorrPlus_h->Write("", TObject::kOverwrite);
  phi_rawCorrMinus_h->Write("", TObject::kOverwrite);

  delete size_h;
  delete size_rawPlus_h;
  delete size_rawMinus_h;
  delete size_rawCorrPlus_h;
  delete size_rawCorrMinus_h;

  delete pt_h;
  delete pt_rawPlus_h;
  delete pt_rawMinus_h;
  delete pt_rawCorrPlus_h;
  delete pt_rawCorrMinus_h;
  
  delete phi_h;
  delete phi_rawPlus_h;
  delete phi_rawMinus_h;
  delete phi_rawCorrPlus_h;
  delete phi_rawCorrMinus_h;

  for(Int_t cI = 0; cI < nCentBins; ++cI){  
    if(hasTrk){
      v2Raw_Trk_h[cI]->Write("", TObject::kOverwrite);
      v2RawCorr_Trk_h[cI]->Write("", TObject::kOverwrite);
      
      v2Obs_Trk_h[cI]->Write("", TObject::kOverwrite);
      v2ObsCorr_Trk_h[cI]->Write("", TObject::kOverwrite);

      v2Raw_Mean_Trk_h->SetBinContent(cI+1, v2Raw_Trk_h[cI]->GetMean());
      v2Raw_Mean_Trk_h->SetBinError(cI+1, v2Raw_Trk_h[cI]->GetMeanError());
      v2Raw_Sigma_Trk_h->SetBinContent(cI+1, v2Raw_Trk_h[cI]->GetStdDev());
      v2Raw_Sigma_Trk_h->SetBinError(cI+1, v2Raw_Trk_h[cI]->GetStdDevError());
      
      v2RawCorr_Mean_Trk_h->SetBinContent(cI+1, v2RawCorr_Trk_h[cI]->GetMean());
      v2RawCorr_Mean_Trk_h->SetBinError(cI+1, v2RawCorr_Trk_h[cI]->GetMeanError());
      v2RawCorr_Sigma_Trk_h->SetBinContent(cI+1, v2RawCorr_Trk_h[cI]->GetStdDev());
      v2RawCorr_Sigma_Trk_h->SetBinError(cI+1, v2RawCorr_Trk_h[cI]->GetStdDevError());
      
      v2Obs_Mean_Trk_h->SetBinContent(cI+1, v2Obs_Trk_h[cI]->GetMean());
      v2Obs_Mean_Trk_h->SetBinError(cI+1, v2Obs_Trk_h[cI]->GetMeanError());
      v2Obs_Sigma_Trk_h->SetBinContent(cI+1, v2Obs_Trk_h[cI]->GetStdDev());
      v2Obs_Sigma_Trk_h->SetBinError(cI+1, v2Obs_Trk_h[cI]->GetStdDevError());
      
      v2ObsCorr_Mean_Trk_h->SetBinContent(cI+1, v2ObsCorr_Trk_h[cI]->GetMean());
      v2ObsCorr_Mean_Trk_h->SetBinError(cI+1, v2ObsCorr_Trk_h[cI]->GetMeanError());
      v2ObsCorr_Sigma_Trk_h->SetBinContent(cI+1, v2ObsCorr_Trk_h[cI]->GetStdDev());
      v2ObsCorr_Sigma_Trk_h->SetBinError(cI+1, v2ObsCorr_Trk_h[cI]->GetStdDevError());

      compRaw_h[cI]->Write("", TObject::kOverwrite);
      compRawCorr_h[cI]->Write("", TObject::kOverwrite);

      compObs_h[cI]->Write("", TObject::kOverwrite);
      compObsCorr_h[cI]->Write("", TObject::kOverwrite);

      compRaw_Mean_h->SetBinContent(cI+1, compRaw_h[cI]->GetMean());
      compRaw_Mean_h->SetBinError(cI+1, compRaw_h[cI]->GetMeanError());
      compRaw_Sigma_h->SetBinContent(cI+1, compRaw_h[cI]->GetStdDev());
      compRaw_Sigma_h->SetBinError(cI+1, compRaw_h[cI]->GetStdDevError());

      compRawCorr_Mean_h->SetBinContent(cI+1, compRawCorr_h[cI]->GetMean());
      compRawCorr_Mean_h->SetBinError(cI+1, compRawCorr_h[cI]->GetMeanError());
      compRawCorr_Sigma_h->SetBinContent(cI+1, compRawCorr_h[cI]->GetStdDev());
      compRawCorr_Sigma_h->SetBinError(cI+1, compRawCorr_h[cI]->GetStdDevError());

      compObs_Mean_h->SetBinContent(cI+1, compObs_h[cI]->GetMean());
      compObs_Mean_h->SetBinError(cI+1, compObs_h[cI]->GetMeanError());
      compObs_Sigma_h->SetBinContent(cI+1, compObs_h[cI]->GetStdDev());
      compObs_Sigma_h->SetBinError(cI+1, compObs_h[cI]->GetStdDevError());

      compObsCorr_Mean_h->SetBinContent(cI+1, compObsCorr_h[cI]->GetMean());
      compObsCorr_Mean_h->SetBinError(cI+1, compObsCorr_h[cI]->GetMeanError());
      compObsCorr_Sigma_h->SetBinContent(cI+1, compObsCorr_h[cI]->GetStdDev());
      compObsCorr_Sigma_h->SetBinError(cI+1, compObsCorr_h[cI]->GetStdDevError());
    }

    delete v2Raw_Trk_h[cI];
    delete v2RawCorr_Trk_h[cI];

    delete v2Obs_Trk_h[cI];
    delete v2ObsCorr_Trk_h[cI];

    delete compRaw_h[cI];
    delete compRawCorr_h[cI];
    
    delete compObs_h[cI];
    delete compObsCorr_h[cI];
  }

  v2Raw_Mean_Trk_h->Write("", TObject::kOverwrite);
  v2Raw_Sigma_Trk_h->Write("", TObject::kOverwrite);

  v2RawCorr_Mean_Trk_h->Write("", TObject::kOverwrite);
  v2RawCorr_Sigma_Trk_h->Write("", TObject::kOverwrite);

  v2Obs_Mean_Trk_h->Write("", TObject::kOverwrite);
  v2Obs_Sigma_Trk_h->Write("", TObject::kOverwrite);

  v2ObsCorr_Mean_Trk_h->Write("", TObject::kOverwrite);
  v2ObsCorr_Sigma_Trk_h->Write("", TObject::kOverwrite);

  compRaw_Mean_h->Write("", TObject::kOverwrite);
  compRaw_Sigma_h->Write("", TObject::kOverwrite);

  compRawCorr_Mean_h->Write("", TObject::kOverwrite);
  compRawCorr_Sigma_h->Write("", TObject::kOverwrite);

  compObs_Mean_h->Write("", TObject::kOverwrite);
  compObs_Sigma_h->Write("", TObject::kOverwrite);

  compObsCorr_Mean_h->Write("", TObject::kOverwrite);
  compObsCorr_Sigma_h->Write("", TObject::kOverwrite);

  if(hasTrk){
    delete v2Raw_Mean_Trk_h;
    delete v2Raw_Sigma_Trk_h;
    delete v2RawCorr_Mean_Trk_h;
    delete v2RawCorr_Sigma_Trk_h;
    
    delete v2Obs_Mean_Trk_h;
    delete v2Obs_Sigma_Trk_h;
    delete v2ObsCorr_Mean_Trk_h;
    delete v2ObsCorr_Sigma_Trk_h;

    delete compRaw_Mean_h;
    delete compRaw_Sigma_h;
    delete compRawCorr_Mean_h;
    delete compRawCorr_Sigma_h;

    delete compObs_Mean_h;
    delete compObs_Sigma_h;
    delete compObsCorr_Mean_h;
    delete compObsCorr_Sigma_h;
  }

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "USAGE: ./recreateV2V3TreeHist.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += recreateV2V3TreeHist(argv[1]);
  return retVal;
}
