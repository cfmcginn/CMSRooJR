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

int recreateV2V3TreeHist(const std::string inFileName)
{
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  TFile* outFile_p = new TFile(("/afs/cern.ch/work/s/skanaski/jetstuff/output/v2CrossCheck_TreeHist_" + dateStr + ".root").c_str(), "RECREATE");

  const Int_t nCentBins = 11;
  const Int_t centBinsLow[nCentBins] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55};
  const Int_t centBinsHi[nCentBins] = {10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60};
  const Double_t centBins[nCentBins+1] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60};

  TH1F* v2Raw_PF_h[nCentBins];
  TH1F* v2RawCorr_PF_h[nCentBins];

  TH1F* v2Obs_PF_h[nCentBins];
  TH1F* v2ObsCorr_PF_h[nCentBins];
  TH1F* v2Comp_h[nCentBins];
  TH1F* v2FT_h[nCentBins];
  TH1F* v2Diff_h[nCentBins];

  TH1F* v2Raw_Trk_h[nCentBins];
  TH1F* v2RawCorr_Trk_h[nCentBins];
  TH1F* v2Fit_Trk_h[nCentBins];
  TH1F* v2Diff_Trk_h[nCentBins];
  TH1F* v2DiffCorr_Trk_h[nCentBins];
  
  TH1F* v2Obs_Trk_h[nCentBins];
  TH1F* v2ObsCorr_Trk_h[nCentBins];

  TH1F* v2Raw_Mean_PF_h = new TH1F("v2Raw_Mean_PF_h", ";Centrality (%);#LTv_{2}^{raw}#GT", nCentBins, centBins);
  TH1F* v2Raw_Sigma_PF_h = new TH1F("v2Raw_Sigma_PF_h", ";Centrality (%);#sigma(v_{2}^{raw})", nCentBins, centBins);
  TH1F* v2RawCorr_Mean_PF_h = new TH1F("v2RawCorr_Mean_PF_h", ";Centrality (%);#LTv_{2}^{raw}#GT", nCentBins, centBins);
  TH1F* v2RawCorr_Sigma_PF_h = new TH1F("v2RawCorr_Sigma_PF_h", ";Centrality (%);#sigma(v_{2}^{raw})", nCentBins, centBins);

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

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    const std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
    v2Raw_PF_h[cI] = new TH1F(("v2Raw_" + centStr + "_PF_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2RawCorr_PF_h[cI] = new TH1F(("v2RawCorr_" + centStr + "_PF_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2Obs_PF_h[cI] = new TH1F(("v2Obs_" + centStr + "_PF_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2ObsCorr_PF_h[cI] = new TH1F(("v2ObsCorr_" + centStr + "_PF_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2Comp_h[cI] = new TH1F(("v2Comp_" + centStr + "_h").c_str(),";v_{2,ObsCorr}/v_{2,Tree};Counts", 150, 0.0, 2.0);
    v2FT_h[cI] = new TH1F(("v2FT_" + centStr + "_h").c_str(),";v_{2,Tree};Counts", 150, 0.0, 0.6);
    v2Diff_h[cI] = new TH1F(("v2Diff_" + centStr + "_h").c_str(),("v2ObsCorr-v2FromTree (" + centStr + ");v_{2};Counts").c_str(), 150, -0.6, 0.6);

    v2Raw_Trk_h[cI] = new TH1F(("v2Raw_" + centStr + "_Trk_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2RawCorr_Trk_h[cI] = new TH1F(("v2RawCorr_" + centStr + "_Trk_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2Fit_Trk_h[cI] = new TH1F(("v2Fit_" + centStr + "_Trk_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2Diff_Trk_h[cI] = new TH1F(("v2Diff_" + centStr + "_Trk_h").c_str(), ("v2RawTrk-v2FitTrk (" + centStr + ");v_{2};Counts").c_str(), 150, -0.6, 0.6);
    v2DiffCorr_Trk_h[cI] = new TH1F(("v2DiffCorr_" + centStr + "_Trk_h").c_str(), ("v2RawTrkCorr-v2FitTrk (" + centStr + ");v_{2};Counts").c_str(), 150, -0.6, 0.6);
    v2Obs_Trk_h[cI] = new TH1F(("v2Obs_" + centStr + "_Trk_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2ObsCorr_Trk_h[cI] = new TH1F(("v2ObsCorr_" + centStr + "_Trk_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
  }

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* inTree_p = (TTree*)inFile_p->Get("v2V3Tree");
  TTree* mntTree_p = (TTree*)inFile_p->Get("mntTree");
  TObjArray* listOfBranches = (TObjArray*)inTree_p->GetListOfBranches();
  std::vector<std::string> vectorOfBranches;
  for(Int_t oI = 0; oI < listOfBranches->GetEntries(); ++oI){
    std::string branch = listOfBranches->At(oI)->GetName();
    vectorOfBranches.push_back(branch);
  }

  /*  bool hasTrk = false;
  for(unsigned int tI = 0; tI < vectorOfBranches.size(); ++tI){
    if(vectorOfBranches.at(tI).find("trk") != std::string::npos){
      hasTrk = true;
      break;
    }
    }*/
  bool hasTrk = true;

  Int_t hiBin_;
  Float_t hiEvt2Plane_;
  Float_t hiEvt3Plane_;
  Float_t v2FromTree_;
  std::vector<float>* pfPhi_p=NULL;
  std::vector<float>* pfWeight_p=NULL;
  std::vector<float>* trkPt_p=NULL;
  std::vector<float>* trkPhi_p=NULL;
  std::vector<float>* trkWeight_p=NULL;  
  int entry_;
  
  inTree_p->SetBranchAddress("hiBin", &hiBin_);
  inTree_p->SetBranchAddress("hiEvt2Plane", &hiEvt2Plane_);
  inTree_p->SetBranchAddress("hiEvt3Plane", &hiEvt3Plane_);
  inTree_p->SetBranchAddress("v2FromTree", &v2FromTree_);
  inTree_p->SetBranchAddress("pfPhi", &pfPhi_p);
  inTree_p->SetBranchAddress("pfWeight", &pfWeight_p);
  inTree_p->SetBranchAddress("entry", &entry_);
  if(hasTrk){
    mntTree_p->SetBranchAddress("hiBin", &hiBin_);
    mntTree_p->SetBranchAddress("trkPt", &trkPt_p);
    mntTree_p->SetBranchAddress("trkPhi", &trkPhi_p);
    mntTree_p->SetBranchAddress("trkWeight", &trkWeight_p);
    mntTree_p->SetBranchAddress("entry", &entry_);
    mntTree_p->SetBranchAddress("hiEvt2Plane", &hiEvt2Plane_);
  }

  const Int_t nEntries = TMath::Min((Int_t)inTree_p->GetEntries(), (Int_t)100000000);
  const Int_t nMntEntries = TMath::Min((Int_t)mntTree_p->GetEntries(), (Int_t)100000000);

  double globalN[nCentBins];
  double globalV2XRawPF[nCentBins];
  double globalV2YRawPF[nCentBins];
  double globalV2XRawPFCorr[nCentBins];
  double globalV2YRawPFCorr[nCentBins];

  double globalV2XRawTrk[nCentBins];
  double globalV2YRawTrk[nCentBins];
  double globalV2XRawTrkCorr[nCentBins];
  double globalV2YRawTrkCorr[nCentBins];

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

    double v2xRawPF = 0.;
    double v2yRawPF = 0.;

    double v2xRawPFCorr = 0.;
    double v2yRawPFCorr = 0.;

    //    double v2xRawTrk = 0.;
    //    double v2yRawTrk = 0.;

    //    double v2xRawTrkCorr = 0.;
    //    double v2yRawTrkCorr = 0.;

    double weightPF = 0.;
    //    double weightTrk = 0.;
  
    for(unsigned pfI = 0; pfI < pfPhi_p->size(); pfI++){
      double tempWeightPF = pfWeight_p->at(pfI);
      //      double deltaEventPhi = pfPhi_p->at(pfI) - hiEvt2Plane_;
      double deltaEventPhi = pfPhi_p->at(pfI);

      v2xRawPF += TMath::Cos(2*(deltaEventPhi));
      v2yRawPF += TMath::Sin(2*(deltaEventPhi));

      v2xRawPFCorr += tempWeightPF*TMath::Cos(2*(deltaEventPhi));
      v2yRawPFCorr += tempWeightPF*TMath::Sin(2*(deltaEventPhi));

      weightPF += tempWeightPF;
    }

    /*    Int_t trkCounter = 0;
    if(hasTrk){    
      for(unsigned tI = 0; tI < trkPhi_p->size(); tI++){
	//	if(trkPt_p->at(tI) >= 2.4) continue;
	trkCounter++;
	double tempWeightTrk = trkWeight_p->at(tI);
	double deltaEventPhi = trkPhi_p->at(tI) - hiEvt2Plane_;
	
4	v2xRawTrk += TMath::Cos(2*(deltaEventPhi));
	v2yRawTrk += TMath::Sin(2*(deltaEventPhi));

	v2xRawTrkCorr += tempWeightTrk*TMath::Cos(2*(deltaEventPhi));
	v2yRawTrkCorr += tempWeightTrk*TMath::Sin(2*(deltaEventPhi));
	
	weightTrk += tempWeightTrk;
      }
      }*/
    if (weightPF==0) continue;

    v2xRawPF /= (double)pfPhi_p->size();
    v2yRawPF /= (double)pfPhi_p->size();

    v2xRawPFCorr /= weightPF;
    v2yRawPFCorr /= weightPF;

    double v2RawPF = TMath::Sqrt(v2xRawPF*v2xRawPF + v2yRawPF*v2yRawPF);
    double v2RawPFCorr = TMath::Sqrt(v2xRawPFCorr*v2xRawPFCorr + v2yRawPFCorr*v2yRawPFCorr);

    v2Raw_PF_h[centPos]->Fill(v2RawPF);
    v2RawCorr_PF_h[centPos]->Fill(v2RawPFCorr);

    //    std::cout<<"PF Event "<<entry<<" v2xRawPF, v2yRawPF, hiBin, hiEvt2Plane: "<<v2xRawPF<<", "<<v2yRawPF<<", "<<hiBin_<<", "<<hiEvt2Plane_<<" (original PF tree entry: "<<entry_<<")"<<std::endl;

    globalN[centPos] += 1;
    globalV2XRawPF[centPos] += v2xRawPF;
    globalV2YRawPF[centPos] += v2yRawPF;
    globalV2XRawPFCorr[centPos] += v2xRawPFCorr;
    globalV2YRawPFCorr[centPos] += v2yRawPFCorr;

    /*
    if(hasTrk){
      v2xRawTrk /= (double)trkCounter;
      v2yRawTrk /= (double)trkCounter;
      
      v2xRawTrkCorr /= weightTrk;
      v2yRawTrkCorr /= weightTrk;
      
      double v2RawTrk = TMath::Sqrt(v2xRawTrk*v2xRawTrk + v2yRawTrk*v2yRawTrk);
      double v2RawTrkCorr = TMath::Sqrt(v2xRawTrkCorr*v2xRawTrkCorr + v2yRawTrkCorr*v2yRawTrkCorr);
      
      v2Raw_Trk_h[centPos]->Fill(v2RawTrk);
      v2RawCorr_Trk_h[centPos]->Fill(v2RawTrkCorr);
      
      globalN[centPos] += 1;
      globalV2XRawTrk[centPos] += v2xRawTrk;
      globalV2YRawTrk[centPos] += v2yRawTrk;
      globalV2XRawTrkCorr[centPos] += v2xRawTrkCorr;
      globalV2YRawTrkCorr[centPos] += v2yRawTrkCorr;
      }*/
  }

  for(Int_t entry = 0; entry < nMntEntries; ++entry){
    if(entry%10000 == 0) std::cout << " Mnt Entry " << entry << "/" << nMntEntries << std::endl;
    mntTree_p->GetEntry(entry);

    int centPos = -1;
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      if(centBinsLow[cI] <= hiBin_/2 && centBinsHi[cI] > hiBin_/2){
        centPos = cI;
        break;
      }
    }

    if(centPos < 0) continue;

    double v2xRawTrk = 0.;
    double v2yRawTrk = 0.;

    double v2xRawTrkCorr = 0.;
    double v2yRawTrkCorr = 0.;

    double weightTrk = 0.;

    Int_t trkCounter = 0;

    const Int_t nPhiBins = std::fmax(10, trkPhi_p->size()/30);

    TH1F* phi_h = new TH1F("tempPhi", ";#phi;Track counts (.3 < p_{T} < 3.,#eta>0)", nPhiBins, -TMath::Pi(), TMath::Pi());

    if(hasTrk){
      for(unsigned tI = 0; tI < trkPhi_p->size(); tI++){
        //      if(trkPt_p->at(tI) >= 2.4) continue;
        trkCounter++;
        double tempWeightTrk = trkWeight_p->at(tI);
	//	double deltaEventPhi = trkPhi_p->at(tI) - hiEvt2Plane_;
 	double deltaEventPhi = trkPhi_p->at(tI);

        v2xRawTrk += TMath::Cos(2*(deltaEventPhi));
        v2yRawTrk += TMath::Sin(2*(deltaEventPhi));

        v2xRawTrkCorr += tempWeightTrk*TMath::Cos(2*(deltaEventPhi));
        v2yRawTrkCorr += tempWeightTrk*TMath::Sin(2*(deltaEventPhi));

        weightTrk += tempWeightTrk;

	phi_h->Fill(trkPhi_p->at(tI));
      }
    }

    if (weightTrk==0) continue;
    /*
    std::string flowFitForm = "[0]*(1. + 2.*([1]*TMath::Cos(2.*(x - " + std::to_string(eventPlane2) + ")) + [2]*TMath::Cos(3.*(x - " + std::to_string(eventPlane3) + "))))";    
    TF1* flowFit_p = new TF1("flowFit", flowFitForm.c_str(), -TMath::Pi(), TMath::Pi());
    flowFit_p->SetParameter(0, 10);
    flowFit_p->SetParameter(1, 0);
    flowFit_p->SetParameter(2, 0);

    phi_h->Fit(flowFit_p, "Q", "", -TMath::Pi(), TMath::Pi());

    double v2FitTrk = flowFit_p->GetParameter(1);
    */
    if(hasTrk){
      v2xRawTrk /= (double)trkCounter;
      v2yRawTrk /= (double)trkCounter;

      v2xRawTrkCorr /= weightTrk;
      v2yRawTrkCorr /= weightTrk;

      //      if (weightTrk==0.) std::cout<<"Divided by 0 "<<entry<<std::endl;
      //      if (v2xRawTrk*v2xRawTrk + v2yRawTrk*v2yRawTrk<0||v2xRawTrkCorr*v2xRawTrkCorr + v2yRawTrkCorr*v2yRawTrkCorr<0) std::cout<<"Sqrt of negative"<<std::endl;

      double v2RawTrk = TMath::Sqrt(v2xRawTrk*v2xRawTrk + v2yRawTrk*v2yRawTrk);
      double v2RawTrkCorr = TMath::Sqrt(v2xRawTrkCorr*v2xRawTrkCorr + v2yRawTrkCorr*v2yRawTrkCorr);

      v2Raw_Trk_h[centPos]->Fill(v2RawTrk);
      v2RawCorr_Trk_h[centPos]->Fill(v2RawTrkCorr);
      //      v2Fit_Trk_h[centPos]->Fill(v2FitTrk);
      //      v2Diff_Trk_h[centPos]->Fill(v2RawTrk-v2FitTrk);
      //      v2DiffCorr_Trk_h[centPos]->Fill(v2RawTrkCorr-v2FitTrk);

      //      std::cout<<"Trk Event "<<entry<<" v2xRawTrk, v2yRawTrk: "<<v2xRawTrk<<", "<<v2yRawTrk<<", "<<hiBin_<<", "<<hiEvt2Plane_<<"  (original trk tree entry: "<<entry_<<")"<<std::endl;

      globalN[centPos] += 1;
      globalV2XRawTrk[centPos] += v2xRawTrk;
      globalV2YRawTrk[centPos] += v2yRawTrk;
      globalV2XRawTrkCorr[centPos] += v2xRawTrkCorr;
      globalV2YRawTrkCorr[centPos] += v2yRawTrkCorr;
    }
  }

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    globalV2XRawPF[cI] /= globalN[cI];
    globalV2YRawPF[cI] /= globalN[cI];
    globalV2XRawPFCorr[cI] /= globalN[cI];
    globalV2YRawPFCorr[cI] /= globalN[cI];

    if(hasTrk){
      globalV2XRawTrk[cI] /= globalN[cI];
      globalV2YRawTrk[cI] /= globalN[cI];
      globalV2XRawTrkCorr[cI] /= globalN[cI];
      globalV2YRawTrkCorr[cI] /= globalN[cI];
    }

    std::cout << "Cent, N: " << centBinsLow[cI] << "-" << centBinsHi[cI] << "%, " << globalN[cI] << std::endl;
    std::cout << "  Global v2xRawPF, v2yRawPF, N: " << globalV2XRawPF[cI] << ", " << globalV2YRawPF[cI] << std::endl;
    std::cout << "  Global v2xRawPFCorr, v2yRawPFCorr, N: " << globalV2XRawPFCorr[cI] << ", " << globalV2YRawPFCorr[cI] << std::endl;

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

    double weightPF = 0.;
  
    for(unsigned pfI = 0; pfI < pfPhi_p->size(); pfI++){
      double tempWeightPF = pfWeight_p->at(pfI);
      //      double deltaEventPhi = pfPhi_p->at(pfI) - hiEvt2Plane_;
      double deltaEventPhi = pfPhi_p->at(pfI);

      v2xObsPF += TMath::Cos(2*(deltaEventPhi));
      v2yObsPF += TMath::Sin(2*(deltaEventPhi));

      v2xObsPFCorr += tempWeightPF*TMath::Cos(2*(deltaEventPhi));
      v2yObsPFCorr += tempWeightPF*TMath::Sin(2*(deltaEventPhi));

      weightPF += tempWeightPF;
    }

    v2xObsPF /= (double)pfPhi_p->size();
    v2yObsPF /= (double)pfPhi_p->size();

    v2xObsPFCorr /= weightPF;
    v2yObsPFCorr /= weightPF;

    v2xObsPF -= globalV2XRawPF[centPos];
    v2yObsPF -= globalV2YRawPF[centPos];

    v2xObsPFCorr -= globalV2XRawPFCorr[centPos];
    v2yObsPFCorr -= globalV2YRawPFCorr[centPos];

    double v2ObsPF = TMath::Sqrt(v2xObsPF*v2xObsPF + v2yObsPF*v2yObsPF);
    double v2ObsPFCorr = TMath::Sqrt(v2xObsPFCorr*v2xObsPFCorr + v2yObsPFCorr*v2yObsPFCorr);

    v2Obs_PF_h[centPos]->Fill(v2ObsPF);
    v2ObsCorr_PF_h[centPos]->Fill(v2ObsPFCorr);
    v2Comp_h[centPos]->Fill(v2ObsPFCorr/v2FromTree_);
    v2FT_h[centPos]->Fill(v2FromTree_);
    v2Diff_h[centPos]->Fill(v2ObsPFCorr-v2FromTree_);

  }

  for(Int_t entry = 0; entry < nMntEntries; ++entry){
    if(entry%10000 == 0) std::cout << " Mnt Entry " << entry << "/" << nMntEntries << std::endl;
    mntTree_p->GetEntry(entry);

    int centPos = -1;
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      if(centBinsLow[cI] <= hiBin_/2 && centBinsHi[cI] > hiBin_/2){
        centPos = cI;
        break;
      }
    }

    if(centPos < 0) continue;

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
    /*
    v2xObsPF /= (double)pfPhi_p->size();
    v2yObsPF /= (double)pfPhi_p->size();

    v2xObsPFCorr /= weightPF;
    v2yObsPFCorr /= weightPF;

    v2xObsPF -= globalV2XRawPF[centPos];
    v2yObsPF -= globalV2YRawPF[centPos];

    v2xObsPFCorr -= globalV2XRawPFCorr[centPos];
    v2yObsPFCorr -= globalV2YRawPFCorr[centPos];

    double v2ObsPF = TMath::Sqrt(v2xObsPF*v2xObsPF + v2yObsPF*v2yObsPF);
    double v2ObsPFCorr = TMath::Sqrt(v2xObsPFCorr*v2xObsPFCorr + v2yObsPFCorr*v2yObsPFCorr);

    v2Obs_PF_h[centPos]->Fill(v2ObsPF);
    v2ObsCorr_PF_h[centPos]->Fill(v2ObsPFCorr);
    */

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
    }
  }


  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  for(Int_t cI = 0; cI < nCentBins; ++cI){  
    v2Raw_PF_h[cI]->Write("", TObject::kOverwrite);
    v2RawCorr_PF_h[cI]->Write("", TObject::kOverwrite);

    v2Obs_PF_h[cI]->Write("", TObject::kOverwrite);
    v2ObsCorr_PF_h[cI]->Write("", TObject::kOverwrite);
    v2Comp_h[cI]->Write("", TObject::kOverwrite);
    v2FT_h[cI]->Write("", TObject::kOverwrite);
    v2Diff_h[cI]->Write("", TObject::kOverwrite);

    v2Raw_Mean_PF_h->SetBinContent(cI+1, v2Raw_PF_h[cI]->GetMean());
    v2Raw_Mean_PF_h->SetBinError(cI+1, v2Raw_PF_h[cI]->GetMeanError());
    v2Raw_Sigma_PF_h->SetBinContent(cI+1, v2Raw_PF_h[cI]->GetStdDev());
    v2Raw_Sigma_PF_h->SetBinError(cI+1, v2Raw_PF_h[cI]->GetStdDevError());

    v2RawCorr_Mean_PF_h->SetBinContent(cI+1, v2RawCorr_PF_h[cI]->GetMean());
    v2RawCorr_Mean_PF_h->SetBinError(cI+1, v2RawCorr_PF_h[cI]->GetMeanError());
    v2RawCorr_Sigma_PF_h->SetBinContent(cI+1, v2RawCorr_PF_h[cI]->GetStdDev());
    v2RawCorr_Sigma_PF_h->SetBinError(cI+1, v2RawCorr_PF_h[cI]->GetStdDevError());

    v2Obs_Mean_PF_h->SetBinContent(cI+1, v2Obs_PF_h[cI]->GetMean());
    v2Obs_Mean_PF_h->SetBinError(cI+1, v2Obs_PF_h[cI]->GetMeanError());
    v2Obs_Sigma_PF_h->SetBinContent(cI+1, v2Obs_PF_h[cI]->GetStdDev());
    v2Obs_Sigma_PF_h->SetBinError(cI+1, v2Obs_PF_h[cI]->GetStdDevError());

    v2ObsCorr_Mean_PF_h->SetBinContent(cI+1, v2ObsCorr_PF_h[cI]->GetMean());
    v2ObsCorr_Mean_PF_h->SetBinError(cI+1, v2ObsCorr_PF_h[cI]->GetMeanError());
    v2ObsCorr_Sigma_PF_h->SetBinContent(cI+1, v2ObsCorr_PF_h[cI]->GetStdDev());
    v2ObsCorr_Sigma_PF_h->SetBinError(cI+1, v2ObsCorr_PF_h[cI]->GetStdDevError());

    delete v2Raw_PF_h[cI];
    delete v2RawCorr_PF_h[cI];
    delete v2Comp_h[cI];
    delete v2FT_h[cI];
    delete v2Diff_h[cI];

    delete v2Obs_PF_h[cI];
    delete v2ObsCorr_PF_h[cI];
  }

  v2Raw_Mean_PF_h->Write("", TObject::kOverwrite);
  v2Raw_Sigma_PF_h->Write("", TObject::kOverwrite);

  v2RawCorr_Mean_PF_h->Write("", TObject::kOverwrite);
  v2RawCorr_Sigma_PF_h->Write("", TObject::kOverwrite);

  v2Obs_Mean_PF_h->Write("", TObject::kOverwrite);
  v2Obs_Sigma_PF_h->Write("", TObject::kOverwrite);

  v2ObsCorr_Mean_PF_h->Write("", TObject::kOverwrite);
  v2ObsCorr_Sigma_PF_h->Write("", TObject::kOverwrite);

  delete v2Raw_Mean_PF_h;
  delete v2Raw_Sigma_PF_h;
  delete v2RawCorr_Mean_PF_h;
  delete v2RawCorr_Sigma_PF_h;

  delete v2Obs_Mean_PF_h;
  delete v2Obs_Sigma_PF_h;
  delete v2ObsCorr_Mean_PF_h;
  delete v2ObsCorr_Sigma_PF_h;

  for(Int_t cI = 0; cI < nCentBins; ++cI){  
    if(hasTrk){
      v2Raw_Trk_h[cI]->Write("", TObject::kOverwrite);
      v2RawCorr_Trk_h[cI]->Write("", TObject::kOverwrite);
      v2Fit_Trk_h[cI]->Write("", TObject::kOverwrite);
      v2Diff_Trk_h[cI]->Write("", TObject::kOverwrite);
      v2DiffCorr_Trk_h[cI]->Write("", TObject::kOverwrite);
      
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
    }

    delete v2Raw_Trk_h[cI];
    delete v2RawCorr_Trk_h[cI];
    delete v2Fit_Trk_h[cI];
    delete v2Diff_Trk_h[cI];
    delete v2DiffCorr_Trk_h[cI];

    delete v2Obs_Trk_h[cI];
    delete v2ObsCorr_Trk_h[cI];
  }

  v2Raw_Mean_Trk_h->Write("", TObject::kOverwrite);
  v2Raw_Sigma_Trk_h->Write("", TObject::kOverwrite);

  v2RawCorr_Mean_Trk_h->Write("", TObject::kOverwrite);
  v2RawCorr_Sigma_Trk_h->Write("", TObject::kOverwrite);

  v2Obs_Mean_Trk_h->Write("", TObject::kOverwrite);
  v2Obs_Sigma_Trk_h->Write("", TObject::kOverwrite);

  v2ObsCorr_Mean_Trk_h->Write("", TObject::kOverwrite);
  v2ObsCorr_Sigma_Trk_h->Write("", TObject::kOverwrite);

  if(hasTrk){
    delete v2Raw_Mean_Trk_h;
    delete v2Raw_Sigma_Trk_h;
    delete v2RawCorr_Mean_Trk_h;
    delete v2RawCorr_Sigma_Trk_h;
    
    delete v2Obs_Mean_Trk_h;
    delete v2Obs_Sigma_Trk_h;
    delete v2ObsCorr_Mean_Trk_h;
    delete v2ObsCorr_Sigma_Trk_h;
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
