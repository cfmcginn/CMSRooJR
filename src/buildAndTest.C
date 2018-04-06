#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLatex.h"

#include "src/RooUnfoldResponse.h"
#include "src/RooUnfoldBayes.h"
#include "src/RooUnfoldSvd.h"

#include "include/histDefUtility.h"
#include "include/plotUtilities.h"
#include "include/kirchnerPalette.h"
#include "include/doGlobalDebug.h"

#include "include/smearingFuncs.h"
#include "include/getLogBins.h"
#include "include/getLinBins.h"

#include "CustomCanvas.h"
#include "TexSlides.C"

int buildAndTest(const std::string inDataName, const std::string inMCName, const bool isPP, const Int_t rVal, const Bool_t isBayes, std::vector<std::vector<std::string>*>* hvv=0)
{
  if (hvv==0) hvv=new std::vector<std::vector<std::string>*>();

  TLatex* label_p = new TLatex();
  label_p->SetTextFont(43);
  label_p->SetTextSize(7);

  TDatime* date = new TDatime();
  kirchnerPalette col;
  TRandom3* randGen_p = new TRandom3(0);
  std::string ppOrPbPbStr = "PbPb";
  if(isPP) ppOrPbPbStr = "PP";
  TFile* outFile_p = new TFile(("output/unfoldTest_" + ppOrPbPbStr + ".root").c_str(), "RECREATE");
  
  std::string mcTreeString = "akCs" + std::to_string(rVal) + "PU3PFFlowJetAnalyzer/akCs" + std::to_string(rVal) + "PU3PFFlowJetAnalyzerOUT";
  if(isPP) mcTreeString = "ak" + std::to_string(rVal) + "PFJetAnalyzer/ak" + std::to_string(rVal) + "PFJetAnalyzerOUT";

  const Float_t jtPtLowPP = 110.;
  const Float_t jtPtLowPbPb = 140.;
  Float_t jtPtLowTemp = jtPtLowPP;
  if(!isPP) jtPtLowTemp = jtPtLowPbPb;
  const Float_t jtPtLow = jtPtLowTemp;

  const Int_t rValReducedCut = 6;

  const Int_t nCentBins = 4;
  const Int_t centBinsLow[nCentBins] = {50, 30, 10, 0};
  const Int_t centBinsHi[nCentBins] = {90, 50, 30, 10};

  const Int_t nAbsEtaBins = 1;
  //  const Float_t absEtaBinsLow[nAbsEtaBins] = {0.0, 0.5, 1.0, 1.6};
  //  const Float_t absEtaBinsHi[nAbsEtaBins] = {0.5, 1.0, 1.6, 2.0};
  const Float_t absEtaBinsLow[nAbsEtaBins] = {0.0};
  const Float_t absEtaBinsHi[nAbsEtaBins] = {2.0};

  const Int_t nTruthBins[nCentBins] = {12, 12, 12, 12};
  const Float_t truthBins5090[nTruthBins[0]+1] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700};
  const Float_t truthBins3050[nTruthBins[1]+1] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700};
  const Float_t truthBins1030[nTruthBins[2]+1] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700};
  const Float_t truthBins010[nTruthBins[3]+1] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700};

  const Int_t nRecoBins[nCentBins] = {10, 10, 10, 10};  
  const Float_t recoBins5090[nRecoBins[0]+1] = {150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650};
  const Float_t recoBins3050[nRecoBins[1]+1] = {150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650};
  const Float_t recoBins1030[nRecoBins[2]+1] = {150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650};
  const Float_t recoBins010[nRecoBins[3]+1] = {150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650};  

  const Int_t nTruthBinsReduced[nCentBins] = {7, 7, 7, 7};  
  const Float_t truthBinsReduced5090[nTruthBinsReduced[0]+1] = {100, 200, 300, 400, 500, 600, 700, 800};
  const Float_t truthBinsReduced3050[nTruthBinsReduced[1]+1] = {100, 200, 300, 400, 500, 600, 700, 800};
  const Float_t truthBinsReduced1030[nTruthBinsReduced[2]+1] = {100, 200, 300, 400, 500, 600, 700, 800};
  const Float_t truthBinsReduced010[nTruthBinsReduced[3]+1] = {100, 200, 300, 400, 500, 600, 700, 800};

  const Int_t nRecoBinsReduced[nCentBins] = {5, 5, 5, 5};  
  const Float_t recoBinsReduced5090[nRecoBinsReduced[0]+1] = {200, 300, 400, 500, 600, 700};
  const Float_t recoBinsReduced3050[nRecoBinsReduced[1]+1] = {200, 300, 400, 500, 600, 700};
  const Float_t recoBinsReduced1030[nRecoBinsReduced[2]+1] = {200, 300, 400, 500, 600, 700};
  const Float_t recoBinsReduced010[nRecoBinsReduced[3]+1] = {200, 300, 400, 500, 600, 700};  

  Float_t truthBins[nCentBins][100];
  Float_t recoBins[nCentBins][100];

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    if(rVal <= rValReducedCut){
      for(Int_t bI = 0; bI < nTruthBins[cI]+1; ++bI){
	if(cI == 0) truthBins[0][bI] = truthBins5090[bI];
	else if(cI == 1) truthBins[1][bI] = truthBins3050[bI];
	else if(cI == 2) truthBins[2][bI] = truthBins1030[bI];
	else if(cI == 3) truthBins[3][bI] = truthBins010[bI];
      }
    }
    else{
      for(Int_t bI = 0; bI < nTruthBinsReduced[cI]+1; ++bI){
	if(cI == 0) truthBins[0][bI] = truthBinsReduced5090[bI];
	else if(cI == 1) truthBins[1][bI] = truthBinsReduced3050[bI];
	else if(cI == 2) truthBins[2][bI] = truthBinsReduced1030[bI];
	else if(cI == 3) truthBins[3][bI] = truthBinsReduced010[bI];
      }
    }

    if(rVal <= rValReducedCut){
      for(Int_t bI = 0; bI < nRecoBins[cI]+1; ++bI){
	if(cI == 0) recoBins[0][bI] = recoBins5090[bI];
	else if(cI == 1) recoBins[1][bI] = recoBins3050[bI];
	else if(cI == 2) recoBins[2][bI] = recoBins1030[bI];
	else if(cI == 3) recoBins[3][bI] = recoBins010[bI];
      }
    }     
    else{
      for(Int_t bI = 0; bI < nRecoBinsReduced[cI]+1; ++bI){
	if(cI == 0) recoBins[0][bI] = recoBinsReduced5090[bI];
	else if(cI == 1) recoBins[1][bI] = recoBinsReduced3050[bI];
	else if(cI == 2) recoBins[2][bI] = recoBinsReduced1030[bI];
	else if(cI == 3) recoBins[3][bI] = recoBinsReduced010[bI];
      }
    }
  }

  TH1D* reco_h[nCentBins][nAbsEtaBins];
  TH1D* recoSymm_h[nCentBins][nAbsEtaBins];
  TH1D* truth_h[nCentBins][nAbsEtaBins];
  TH1D* recoPerp_h[nCentBins][nAbsEtaBins];
  TH1D* recoSymmPerp_h[nCentBins][nAbsEtaBins];
  TH1D* truthPerp_h[nCentBins][nAbsEtaBins];
  TH2D* responseSymm_h[nCentBins][nAbsEtaBins];
  TH2D* responseAsymm_h[nCentBins][nAbsEtaBins];
  RooUnfoldResponse* rooResponse_h[nCentBins][nAbsEtaBins];

  TH1D* recoData_h[nCentBins][nAbsEtaBins];
  TH1D* recoSymmData_h[nCentBins][nAbsEtaBins];

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    const std::string centStr = ppOrPbPbStr + "_Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);

    for(Int_t aI = 0; aI < nAbsEtaBins; ++aI){
      const std::string absEtaStr = "AbsEta" + prettyString(absEtaBinsLow[aI], true, 1) + "to" + prettyString(absEtaBinsHi[aI], true, 1);
      const std::string totStr = centStr + "_" + absEtaStr;

      if(rVal <= rValReducedCut){
	reco_h[cI][aI] = new TH1D(("reco_" + totStr + "_h").c_str(), ";Reco. Jet p_{T};Counts", nRecoBins[cI], recoBins[cI]);
	recoSymm_h[cI][aI] = new TH1D(("recoSymm_" + totStr + "_h").c_str(), ";Reco. Jet p_{T};Counts", nTruthBins[cI], truthBins[cI]);
	truth_h[cI][aI] = new TH1D(("truth_" + totStr + "_h").c_str(), ";Truth. Jet p_{T};Counts", nTruthBins[cI], truthBins[cI]);
	recoPerp_h[cI][aI] = new TH1D(("recoPerp_" + totStr + "_h").c_str(), ";Reco. Jet p_{T};Counts", nRecoBins[cI], recoBins[cI]);
	recoSymmPerp_h[cI][aI] = new TH1D(("recoSymmPerp_" + totStr + "_h").c_str(), ";Reco. Jet p_{T};Counts", nTruthBins[cI], truthBins[cI]);
	truthPerp_h[cI][aI] = new TH1D(("truthPerp_" + totStr + "_h").c_str(), ";Truth. Jet p_{T};Counts", nTruthBins[cI], truthBins[cI]);
	
	responseSymm_h[cI][aI] = new TH2D(("responseSymm_" + totStr + "_h").c_str(), ";Reco. Jet p_{T};Truth Jet p_{T}", nTruthBins[cI], truthBins[cI], nTruthBins[cI], truthBins[cI]);
	responseAsymm_h[cI][aI] = new TH2D(("responseAsymm_" + totStr + "_h").c_str(), ";Reco. Jet p_{T};Truth Jet p_{T}", nRecoBins[cI], recoBins[cI], nTruthBins[cI], truthBins[cI]);
	
	rooResponse_h[cI][aI] = new RooUnfoldResponse(("rooResponse_" + totStr + "_h").c_str(), "");
	rooResponse_h[cI][aI]->Setup(reco_h[cI][aI], truth_h[cI][aI]);
	
	recoData_h[cI][aI] = new TH1D(("recoData_" + totStr + "_h").c_str(), ";Reco. Jet p_{T};Counts", nRecoBins[cI], recoBins[cI]);
	recoSymmData_h[cI][aI] = new TH1D(("recoSymmData_" + totStr + "_h").c_str(), ";Reco. Jet p_{T};Counts", nTruthBins[cI], truthBins[cI]);
      }
      else{
	reco_h[cI][aI] = new TH1D(("reco_" + totStr + "_h").c_str(), ";Reco. Jet p_{T};Counts", nRecoBinsReduced[cI], recoBins[cI]);
	recoSymm_h[cI][aI] = new TH1D(("recoSymm_" + totStr + "_h").c_str(), ";Reco. Jet p_{T};Counts", nTruthBinsReduced[cI], truthBins[cI]);
	truth_h[cI][aI] = new TH1D(("truth_" + totStr + "_h").c_str(), ";Truth. Jet p_{T};Counts", nTruthBinsReduced[cI], truthBins[cI]);
	recoPerp_h[cI][aI] = new TH1D(("recoPerp_" + totStr + "_h").c_str(), ";Reco. Jet p_{T};Counts", nRecoBinsReduced[cI], recoBins[cI]);
	recoSymmPerp_h[cI][aI] = new TH1D(("recoSymmPerp_" + totStr + "_h").c_str(), ";Reco. Jet p_{T};Counts", nTruthBinsReduced[cI], truthBins[cI]);
	truthPerp_h[cI][aI] = new TH1D(("truthPerp_" + totStr + "_h").c_str(), ";Truth. Jet p_{T};Counts", nTruthBinsReduced[cI], truthBins[cI]);
	
	responseSymm_h[cI][aI] = new TH2D(("responseSymm_" + totStr + "_h").c_str(), ";Reco. Jet p_{T};Truth Jet p_{T}", nTruthBinsReduced[cI], truthBins[cI], nTruthBinsReduced[cI], truthBins[cI]);
	responseAsymm_h[cI][aI] = new TH2D(("responseAsymm_" + totStr + "_h").c_str(), ";Reco. Jet p_{T};Truth Jet p_{T}", nRecoBinsReduced[cI], recoBins[cI], nTruthBinsReduced[cI], truthBins[cI]);
	
	rooResponse_h[cI][aI] = new RooUnfoldResponse(("rooResponse_" + totStr + "_h").c_str(), "");
	rooResponse_h[cI][aI]->Setup(reco_h[cI][aI], truth_h[cI][aI]);
	
	recoData_h[cI][aI] = new TH1D(("recoData_" + totStr + "_h").c_str(), ";Reco. Jet p_{T};Counts", nRecoBinsReduced[cI], recoBins[cI]);
	recoSymmData_h[cI][aI] = new TH1D(("recoSymmData_" + totStr + "_h").c_str(), ";Reco. Jet p_{T};Counts", nTruthBinsReduced[cI], truthBins[cI]);
      }

      centerTitles({reco_h[cI][aI], recoSymm_h[cI][aI], truth_h[cI][aI], recoPerp_h[cI][aI], recoSymmPerp_h[cI][aI], truthPerp_h[cI][aI], responseSymm_h[cI][aI], responseAsymm_h[cI][aI], recoData_h[cI][aI], recoSymmData_h[cI][aI]});
      setSumW2({reco_h[cI][aI], recoSymm_h[cI][aI], truth_h[cI][aI], recoPerp_h[cI][aI], recoSymmPerp_h[cI][aI], truthPerp_h[cI][aI], responseSymm_h[cI][aI], responseAsymm_h[cI][aI], recoData_h[cI][aI], recoSymmData_h[cI][aI]});

    }
  }

  TFile* inMCFile_p = new TFile(inMCName.c_str(), "READ");

  Float_t fullWeight_;
  Int_t hiBin_;
  const Int_t nMaxJets = 500;
  Int_t nref_;
  Float_t jtpt_[nMaxJets];
  Float_t jteta_[nMaxJets];
  Float_t refpt_[nMaxJets];

  TTree* inResponseTree_p = (TTree*)inMCFile_p->Get(mcTreeString.c_str());

  inResponseTree_p->SetBranchStatus("*", 0);
  inResponseTree_p->SetBranchStatus("fullWeight", 1);
  inResponseTree_p->SetBranchStatus("hiBin", 1);
  inResponseTree_p->SetBranchStatus("nref", 1);
  inResponseTree_p->SetBranchStatus("jtpt", 1);
  inResponseTree_p->SetBranchStatus("jteta", 1);
  inResponseTree_p->SetBranchStatus("refpt", 1);

  inResponseTree_p->SetBranchAddress("fullWeight", &fullWeight_);
  inResponseTree_p->SetBranchAddress("hiBin", &hiBin_);
  inResponseTree_p->SetBranchAddress("nref", &nref_);
  inResponseTree_p->SetBranchAddress("jtpt", jtpt_);
  inResponseTree_p->SetBranchAddress("jteta", jteta_);
  inResponseTree_p->SetBranchAddress("refpt", refpt_);

  const Int_t nMCEntries = inResponseTree_p->GetEntries();

  std::cout << "Processing MC..." << std::endl;
  for(Int_t entry = 0; entry < nMCEntries; ++entry){
    if(entry%100000 == 0) std::cout << " Entry " << entry << "/" << nMCEntries << std::endl;
    inResponseTree_p->GetEntry(entry);

    std::vector<int> centPos;
    if(isPP){
      for(Int_t cI = 0; cI < nCentBins; ++cI){
	centPos.push_back(cI);
      }
    }
    else{
      for(Int_t cI = 0; cI < nCentBins; ++cI){
	if(centBinsLow[cI]*2 <= hiBin_ && centBinsHi[cI]*2 > hiBin_){
	  centPos.push_back(cI);
	  break;
	}
      }
    }
    
    bool isPerp = randGen_p->Uniform(0.0, 1.0) > 0.9;
    
    for(Int_t jI = 0; jI < nref_; ++jI){
      if(jtpt_[jI] < jtPtLow) continue;

      for(unsigned int cI = 0; cI < centPos.size(); ++cI){
	bool goodReco = false;
	bool goodReco2 = false;
	bool goodTruth = false;

	if(rVal <= rValReducedCut){
	  goodReco = jtpt_[jI] >= recoBins[centPos.at(cI)][0] && jtpt_[jI] < recoBins[centPos.at(cI)][nRecoBins[centPos.at(cI)]];
	  goodReco2 = jtpt_[jI] >= truthBins[centPos.at(cI)][0] && jtpt_[jI] < truthBins[centPos.at(cI)][nTruthBins[centPos.at(cI)]];
	  goodTruth = refpt_[jI] >= truthBins[centPos.at(cI)][0] && refpt_[jI] < truthBins[centPos.at(cI)][nTruthBins[centPos.at(cI)]];
	}
	else{
	  goodReco = jtpt_[jI] >= recoBins[centPos.at(cI)][0] && jtpt_[jI] < recoBins[centPos.at(cI)][nRecoBinsReduced[centPos.at(cI)]];
	  goodReco2 = jtpt_[jI] >= truthBins[centPos.at(cI)][0] && jtpt_[jI] < truthBins[centPos.at(cI)][nTruthBinsReduced[centPos.at(cI)]];
	  goodTruth = refpt_[jI] >= truthBins[centPos.at(cI)][0] && refpt_[jI] < truthBins[centPos.at(cI)][nTruthBinsReduced[centPos.at(cI)]];
	}

	Int_t absEtaPos = -1;
	for(Int_t aI = 0; aI < nAbsEtaBins; ++aI){
	  if(TMath::Abs(jteta_[jI]) >= absEtaBinsLow[aI] && TMath::Abs(jteta_[jI]) < absEtaBinsHi[aI]){
	    absEtaPos = aI;
	    break;
	  }
	}
	if(absEtaPos < 0) continue;

	if(isPerp){
	  if(goodReco && goodTruth){
	    recoPerp_h[centPos.at(cI)][absEtaPos]->Fill(jtpt_[jI], fullWeight_);
	    truthPerp_h[centPos.at(cI)][absEtaPos]->Fill(refpt_[jI], fullWeight_);
	  }
	  else if(goodTruth) truthPerp_h[centPos.at(cI)][absEtaPos]->Fill(refpt_[jI], fullWeight_);

	  if(goodReco2 && goodTruth) recoSymmPerp_h[centPos.at(cI)][absEtaPos]->Fill(jtpt_[jI], fullWeight_);
	}
	else{
	  if(goodReco && goodTruth){
	    reco_h[centPos.at(cI)][absEtaPos]->Fill(jtpt_[jI], fullWeight_);
	    truth_h[centPos.at(cI)][absEtaPos]->Fill(refpt_[jI], fullWeight_);	    
	    responseAsymm_h[centPos.at(cI)][absEtaPos]->Fill(jtpt_[jI], refpt_[jI], fullWeight_);
	    rooResponse_h[centPos.at(cI)][absEtaPos]->Fill(jtpt_[jI], refpt_[jI], fullWeight_);
	  }
	  else if(goodTruth){
	    rooResponse_h[centPos.at(cI)][absEtaPos]->Miss(refpt_[jI], fullWeight_);
	    truth_h[centPos.at(cI)][absEtaPos]->Fill(refpt_[jI], fullWeight_);
	  }

	  if(goodReco2 && goodTruth){
	    recoSymm_h[centPos.at(cI)][absEtaPos]->Fill(jtpt_[jI], fullWeight_);
	    responseSymm_h[centPos.at(cI)][absEtaPos]->Fill(jtpt_[jI], refpt_[jI], fullWeight_);
	  }

	}
      }
    }
  }

  inMCFile_p->Close();
  delete inMCFile_p;

  TFile* inDataFile_p = new TFile(inDataName.c_str(), "READ");
  std::string dataTreeString = "akCs" + std::to_string(rVal) + "PU3PFFlow/akCs" + std::to_string(rVal) + "PU3PFFlowJetAnalyzerOUT";
  if(isPP) dataTreeString = "ak" + std::to_string(rVal) + "PF/ak" + std::to_string(rVal) + "PFJetAnalyzerOUT";
  TTree* inDataTree_p = (TTree*)inDataFile_p->Get(dataTreeString.c_str());

  inDataTree_p->SetBranchStatus("*", 0);
  inDataTree_p->SetBranchStatus("hiBin", 1);
  inDataTree_p->SetBranchStatus("nref", 1);
  inDataTree_p->SetBranchStatus("jtpt", 1);
  inDataTree_p->SetBranchStatus("jteta", 1);

  inDataTree_p->SetBranchAddress("hiBin", &hiBin_);
  inDataTree_p->SetBranchAddress("nref", &nref_);
  inDataTree_p->SetBranchAddress("jtpt", jtpt_);
  inDataTree_p->SetBranchAddress("jteta", jteta_);

  const Int_t nDataEntries = inDataTree_p->GetEntries();

  std::cout << "Processing Data..." << std::endl;
  for(Int_t entry = 0; entry < nDataEntries; ++entry){
    if(entry%100000 == 0) std::cout << " Entry " << entry << "/" << nDataEntries << std::endl;
    inDataTree_p->GetEntry(entry);

    std::vector<int> centPos;
    if(isPP){
      for(Int_t cI = 0; cI < nCentBins; ++cI){
	centPos.push_back(cI);
      }
    }
    else{
      for(Int_t cI = 0; cI < nCentBins; ++cI){
	if(centBinsLow[cI]*2 <= hiBin_ && centBinsHi[cI]*2 > hiBin_){
	  centPos.push_back(cI);
	  break;
	}
      }
    }
    
    for(Int_t jI = 0; jI < nref_; ++jI){
      if(jtpt_[jI] < jtPtLow) continue;

      for(unsigned int cI = 0; cI < centPos.size(); ++cI){
        Int_t absEtaPos = -1;
        for(Int_t aI = 0; aI < nAbsEtaBins; ++aI){
          if(TMath::Abs(jteta_[jI]) >= absEtaBinsLow[aI] && TMath::Abs(jteta_[jI]) < absEtaBinsHi[aI]){
            absEtaPos = aI;
            break;
          }
        }

	if(absEtaPos < 0) std::cout << "Absetapos " << absEtaPos << std::endl;
	if(absEtaPos < 0) continue;

	bool goodReco = false;
	bool goodReco2 = false;

	if(rVal <= rValReducedCut){
	  goodReco = jtpt_[jI] >= recoBins[centPos.at(cI)][0] && jtpt_[jI] < recoBins[centPos.at(cI)][nRecoBins[centPos.at(cI)]];
	  goodReco2 = jtpt_[jI] >= truthBins[centPos.at(cI)][0] && jtpt_[jI] < truthBins[centPos.at(cI)][nTruthBins[centPos.at(cI)]];
	}
	else{
	  goodReco = jtpt_[jI] >= recoBins[centPos.at(cI)][0] && jtpt_[jI] < recoBins[centPos.at(cI)][nRecoBinsReduced[centPos.at(cI)]];
	  goodReco2 = jtpt_[jI] >= truthBins[centPos.at(cI)][0] && jtpt_[jI] < truthBins[centPos.at(cI)][nTruthBinsReduced[centPos.at(cI)]];
	}

	if(goodReco) recoData_h[centPos.at(cI)][absEtaPos]->Fill(jtpt_[jI]);
	if(goodReco2) recoSymmData_h[centPos.at(cI)][absEtaPos]->Fill(jtpt_[jI]);
      }     
    }
  }

  outFile_p->cd();

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    const std::string centStr = ppOrPbPbStr + "_Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);

    for(Int_t aI = 0; aI < nAbsEtaBins; ++aI){
      const std::string absEtaStr = "AbsEta" + prettyString(absEtaBinsLow[aI], true, 1) + "to" + prettyString(absEtaBinsHi[aI], true, 1);
      const std::string totStr = centStr + "_" + absEtaStr;
      
      truth_h[cI][aI]->SetMarkerSize(0.8);
      truth_h[cI][aI]->SetMarkerStyle(24);
      truth_h[cI][aI]->SetMarkerColor(col.getColor(2));
      truth_h[cI][aI]->SetLineColor(col.getColor(2));
      
      truth_h[cI][aI]->GetXaxis()->SetTitleFont(43);
      truth_h[cI][aI]->GetYaxis()->SetTitleFont(43);
      truth_h[cI][aI]->GetXaxis()->SetLabelFont(43);
      truth_h[cI][aI]->GetYaxis()->SetLabelFont(43);
      
      truth_h[cI][aI]->GetXaxis()->SetTitleSize(8);
      truth_h[cI][aI]->GetYaxis()->SetTitleSize(8);
      truth_h[cI][aI]->GetXaxis()->SetLabelSize(8);
      truth_h[cI][aI]->GetYaxis()->SetLabelSize(8);
      
      truthPerp_h[cI][aI]->SetMarkerSize(0.8);
      truthPerp_h[cI][aI]->SetMarkerStyle(24);
      truthPerp_h[cI][aI]->SetMarkerColor(col.getColor(2));
      truthPerp_h[cI][aI]->SetLineColor(col.getColor(2));
      
      truthPerp_h[cI][aI]->GetXaxis()->SetTitleFont(43);
      truthPerp_h[cI][aI]->GetYaxis()->SetTitleFont(43);
      truthPerp_h[cI][aI]->GetXaxis()->SetLabelFont(43);
      truthPerp_h[cI][aI]->GetYaxis()->SetLabelFont(43);
      
      truthPerp_h[cI][aI]->GetXaxis()->SetTitleSize(8);
      truthPerp_h[cI][aI]->GetYaxis()->SetTitleSize(8);
      truthPerp_h[cI][aI]->GetXaxis()->SetLabelSize(8);
      truthPerp_h[cI][aI]->GetYaxis()->SetLabelSize(8);
      
      TH1D* prevHist_p = (TH1D*)(recoData_h[cI][aI]->Clone("prevHist"));
      if(!isBayes) prevHist_p = (TH1D*)(recoSymmData_h[cI][aI]->Clone("prevHist"));

      Int_t bayesStartVal = 1;
      Int_t bayesEndVal = 11;
      if(!isBayes){
	if(rVal <= rValReducedCut) bayesEndVal = nTruthBins[cI];
	else bayesEndVal = nTruthBinsReduced[cI];
      }

      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    
      for(Int_t bI = bayesStartVal; bI <= bayesEndVal; ++bI){
	CustomCanvas* canv_p = new CustomCanvas(("build_" + totStr + "_Bayes" + std::to_string(bI) + "_c").c_str(), "", 4*300, 2*300);
	const Int_t nPads = 12;
	TPad* pads[nPads];
	canv_p->SetTopMargin(0.0);
	canv_p->SetRightMargin(0.0);
	gStyle->SetOptStat(0);
	
	Double_t xValsLow[4] = {0.0, 0.25, 0.50, 0.75};
	Double_t xValsHi[4] = {0.25, 0.50, 0.75, 1.00};
	
	Double_t yValsLow[5] = {0.5, 0.675, 0.0, .175, 0.0};
	Double_t yValsHi[5] = {0.675, 1.0, 0.5, 0.5, .175};
	
	for(Int_t pI = 0; pI < nPads; ++pI){
	  canv_p->cd();
	  if(pI == 0) pads[pI] = new TPad(("pad" + std::to_string(pI)).c_str(), "", xValsLow[0], yValsLow[1], xValsHi[0], yValsHi[1]);
	  else if(pI == 1) pads[pI] = new TPad(("pad" + std::to_string(pI)).c_str(), "", xValsLow[0], yValsLow[0], xValsHi[0], yValsHi[0]);
	  else if(pI == 2) pads[pI] = new TPad(("pad" + std::to_string(pI)).c_str(), "", xValsLow[1], yValsLow[1], xValsHi[1], yValsHi[1]);
	  else if(pI == 3) pads[pI] = new TPad(("pad" + std::to_string(pI)).c_str(), "", xValsLow[1], yValsLow[0], xValsHi[1], yValsHi[0]);
	  else if(pI == 4) pads[pI] = new TPad(("pad" + std::to_string(pI)).c_str(), "", xValsLow[2], yValsLow[1], xValsHi[2], yValsHi[1]);
	  else if(pI == 5) pads[pI] = new TPad(("pad" + std::to_string(pI)).c_str(), "", xValsLow[2], yValsLow[0], xValsHi[2], yValsHi[0]);
	  else if(pI == 6) pads[pI] = new TPad(("pad" + std::to_string(pI)).c_str(), "", xValsLow[3], yValsLow[1], xValsHi[3], yValsHi[1]);
	  else if(pI == 7) pads[pI] = new TPad(("pad" + std::to_string(pI)).c_str(), "", xValsLow[3], yValsLow[0], xValsHi[3], yValsHi[0]);
	  else if(pI == 8) pads[pI] = new TPad(("pad" + std::to_string(pI)).c_str(), "", xValsLow[0], yValsLow[2], xValsHi[0], yValsHi[2]);
	  else if(pI == 9) pads[pI] = new TPad(("pad" + std::to_string(pI)).c_str(), "", xValsLow[1], yValsLow[2], xValsHi[1], yValsHi[2]);
	  else if(pI == 10) pads[pI] = new TPad(("pad" + std::to_string(pI)).c_str(), "", xValsLow[2], yValsLow[2], xValsHi[2], yValsHi[2]);
	  else if(pI == 11) pads[pI] = new TPad(("pad" + std::to_string(pI)).c_str(), "", xValsLow[3], yValsLow[2], xValsHi[3], yValsHi[2]);

	  pads[pI]->Draw("SAME");
	  pads[pI]->cd();
	  if(pI < 8) pads[pI]->SetRightMargin(0.01);
	  pads[pI]->SetTopMargin(0.01);
	  if(pI < 8){
	    if(pI%2 == 0) pads[pI]->SetBottomMargin(0.01);
	    else pads[pI]->SetBottomMargin(pads[pI]->GetLeftMargin()*2.);
	  }
	  else pads[pI]->SetBottomMargin(pads[pI]->GetLeftMargin()*2.);
	}

	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	
	canv_p->cd();
	pads[0]->cd();
	gPad->SetLogy();
	//	gPad->SetLogx();
	
	truth_h[cI][aI]->DrawCopy("HIST E1 P");
	
	canv_p->cd();
	pads[2]->cd();
	gPad->SetLogy();
	//gPad->SetLogx();

	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	
	truthPerp_h[cI][aI]->DrawCopy("HIST E1 P");
	
	RooUnfoldBayes* tempBayes = new RooUnfoldBayes(rooResponse_h[cI][aI], reco_h[cI][aI], bI);
	RooUnfoldSvd* tempSvd = new RooUnfoldSvd(rooResponse_h[cI][aI], reco_h[cI][aI], bI);
	tempBayes->SetVerbose(0);
	tempSvd->SetVerbose(0);

	TH1D* tempReco = NULL;
	if(isBayes) tempReco = (TH1D*)tempBayes->Hreco(RooUnfold::kCovToy);
	else tempReco = (TH1D*)tempSvd->Hreco(RooUnfold::kCovToy);

	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	tempReco->SetTitle("");
	centerTitles(tempReco);

	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	
	tempReco->SetMarkerSize(0.8);
	tempReco->SetMarkerStyle(24);
	tempReco->SetMarkerColor(col.getColor(0));
	tempReco->SetLineColor(col.getColor(0));

	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	
	canv_p->cd();
	pads[0]->cd();
	
	tempReco->DrawCopy("HIST E1 P SAME");
	
	canv_p->cd();
	pads[1]->cd();
	//gPad->SetLogx();
	tempReco->Divide(truth_h[cI][aI]);
	tempReco->SetMaximum(1.2);
	tempReco->SetMinimum(0.8);

	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	
	tempReco->GetYaxis()->SetNdivisions(404);
	
	tempReco->GetYaxis()->SetTitle("Ratio (Orig. MC Unfolded/Truth)");
	
	tempReco->GetXaxis()->SetTitleFont(43);
	tempReco->GetYaxis()->SetTitleFont(43);
	tempReco->GetXaxis()->SetLabelFont(43);
	tempReco->GetYaxis()->SetLabelFont(43);
	
	tempReco->GetXaxis()->SetTitleSize(8);
	tempReco->GetYaxis()->SetTitleSize(8);
	tempReco->GetXaxis()->SetLabelSize(8);
	tempReco->GetYaxis()->SetLabelSize(8);
	
	tempReco->DrawCopy("HIST E1 P");
	

	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	delete tempReco;
	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	delete tempBayes;
	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	delete tempSvd;
	
	tempBayes = new RooUnfoldBayes(rooResponse_h[cI][aI], recoPerp_h[cI][aI], bI);
	tempSvd = new RooUnfoldSvd(rooResponse_h[cI][aI], recoPerp_h[cI][aI], bI);
	tempBayes->SetVerbose(0);
	tempSvd->SetVerbose(0);

	tempReco = NULL;
	if(isBayes) tempReco = (TH1D*)tempBayes->Hreco(RooUnfold::kCovToy);
	else tempReco = (TH1D*)tempSvd->Hreco(RooUnfold::kCovToy);

	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	tempReco->SetTitle("");
	centerTitles(tempReco);
	tempReco->SetMarkerSize(0.8);
	tempReco->SetMarkerStyle(24);
	tempReco->SetMarkerColor(col.getColor(0));
	tempReco->SetLineColor(col.getColor(0));
	
	canv_p->cd();
	pads[2]->cd();
	tempReco->DrawCopy("HIST E1 P SAME");
	
	canv_p->cd();
	pads[3]->cd();
	//gPad->SetLogx();
	tempReco->Divide(truthPerp_h[cI][aI]);
	tempReco->SetMaximum(1.2);
	tempReco->SetMinimum(0.8);
	tempReco->GetYaxis()->SetNdivisions(404);
	
	tempReco->GetYaxis()->SetTitle("Ratio (Perp. MC Unfolded/Truth)");
	
	tempReco->GetXaxis()->SetTitleFont(43);
	tempReco->GetYaxis()->SetTitleFont(43);
	tempReco->GetXaxis()->SetLabelFont(43);
	tempReco->GetYaxis()->SetLabelFont(43);
	
	tempReco->GetXaxis()->SetTitleSize(8);
	tempReco->GetYaxis()->SetTitleSize(8);
	tempReco->GetXaxis()->SetLabelSize(8);
	tempReco->GetYaxis()->SetLabelSize(8);
	
	tempReco->DrawCopy("HIST E1 P");
	
	delete tempReco;
	delete tempBayes;
	delete tempSvd;

	tempBayes = new RooUnfoldBayes(rooResponse_h[cI][aI], recoData_h[cI][aI], bI);
	tempSvd = new RooUnfoldSvd(rooResponse_h[cI][aI], recoData_h[cI][aI], bI);
	tempBayes->SetVerbose(0);
	tempSvd->SetVerbose(0);

	tempReco = NULL;
	if(isBayes) tempReco = (TH1D*)tempBayes->Hreco(RooUnfold::kCovToy);
	else tempReco = (TH1D*)tempSvd->Hreco(RooUnfold::kCovToy);

	tempReco->SetTitle("");
	centerTitles(tempReco);
	tempReco->SetMarkerSize(0.8);
	tempReco->SetMarkerStyle(24);
	tempReco->SetMarkerColor(col.getColor(0));
	tempReco->SetLineColor(col.getColor(0));
	
	canv_p->cd();
	pads[4]->cd();
	gPad->SetLogy();
	//gPad->SetLogx();
	
	tempReco->GetXaxis()->SetTitleFont(43);
	tempReco->GetYaxis()->SetTitleFont(43);
	tempReco->GetXaxis()->SetLabelFont(43);
	tempReco->GetYaxis()->SetLabelFont(43);
	
	tempReco->GetXaxis()->SetTitleSize(8);
	tempReco->GetYaxis()->SetTitleSize(8);
	tempReco->GetXaxis()->SetLabelSize(8);
	tempReco->GetYaxis()->SetLabelSize(8);
	
	tempReco->DrawCopy("HIST E1 P");
	
	if(prevHist_p != NULL){
	  prevHist_p->GetXaxis()->SetTitleFont(43);
	  prevHist_p->GetYaxis()->SetTitleFont(43);
	  prevHist_p->GetXaxis()->SetLabelFont(43);
	  prevHist_p->GetYaxis()->SetLabelFont(43);
	  
	  prevHist_p->GetXaxis()->SetTitleSize(8);
	  prevHist_p->GetYaxis()->SetTitleSize(8);
	  prevHist_p->GetXaxis()->SetLabelSize(8);
	  prevHist_p->GetYaxis()->SetLabelSize(8);
	  
	  prevHist_p->SetMarkerStyle(24);
	  prevHist_p->SetMarkerSize(0.8);
	  prevHist_p->SetMarkerColor(col.getColor(2));
	  prevHist_p->SetLineColor(col.getColor(2));

	  prevHist_p->DrawCopy("HIST E1 P SAME");
	  
	  canv_p->cd();
	  pads[5]->cd();
	  //gPad->SetLogx();
	  
	  prevHist_p->Divide(tempReco);
	  prevHist_p->SetMaximum(1.2);
	  prevHist_p->SetMinimum(0.8);	
	  prevHist_p->GetYaxis()->SetNdivisions(404);
	  
	  prevHist_p->GetYaxis()->SetTitle("Ratio (Prev./Current)");
	  
	  prevHist_p->SetMarkerColor(col.getColor(0));
	  prevHist_p->SetLineColor(col.getColor(0));
	  
	  prevHist_p->DrawCopy("HIST E1 P");
	}
	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	prevHist_p = (TH1D*)tempReco->Clone(("prev" + std::to_string(bI)).c_str());
	
	canv_p->cd();
	pads[6]->cd();
	gPad->SetLogy();
	//gPad->SetLogx();
	
	recoData_h[cI][aI]->GetXaxis()->SetTitleFont(43);
	recoData_h[cI][aI]->GetYaxis()->SetTitleFont(43);
	recoData_h[cI][aI]->GetXaxis()->SetLabelFont(43);
	recoData_h[cI][aI]->GetYaxis()->SetLabelFont(43);
	
	recoData_h[cI][aI]->GetXaxis()->SetTitleSize(8);
	recoData_h[cI][aI]->GetYaxis()->SetTitleSize(8);
	recoData_h[cI][aI]->GetXaxis()->SetLabelSize(8);
	recoData_h[cI][aI]->GetYaxis()->SetLabelSize(8);
	
	recoData_h[cI][aI]->SetMarkerColor(col.getColor(2));
	recoData_h[cI][aI]->SetLineColor(col.getColor(2));
	recoData_h[cI][aI]->SetMarkerStyle(24);
	recoData_h[cI][aI]->SetMarkerSize(0.8);
	
	recoData_h[cI][aI]->DrawCopy("HIST E1 P");
	
	TH1D* refolded_h = (TH1D*)rooResponse_h[cI][aI]->ApplyToTruth(tempReco);
	
	refolded_h->SetMarkerSize(0.8);
	refolded_h->SetMarkerStyle(24);
	refolded_h->SetMarkerColor(col.getColor(0));
	refolded_h->SetLineColor(col.getColor(0));
	
	refolded_h->DrawCopy("HIST E1 P SAME");
	
	canv_p->cd();
	pads[7]->cd();
	//gPad->SetLogx();
	
	refolded_h->SetTitle("");
	refolded_h->Divide(recoData_h[cI][aI]);
	
	refolded_h->GetYaxis()->SetNdivisions(404);
	
	refolded_h->GetXaxis()->SetTitleFont(43);
	refolded_h->GetYaxis()->SetTitleFont(43);
	refolded_h->GetXaxis()->SetLabelFont(43);
	refolded_h->GetYaxis()->SetLabelFont(43);
	
	refolded_h->GetXaxis()->SetTitleSize(8);
	refolded_h->GetYaxis()->SetTitleSize(8);
	refolded_h->GetXaxis()->SetLabelSize(8);
	refolded_h->GetYaxis()->SetLabelSize(8);
	
	centerTitles(refolded_h);
	refolded_h->SetMaximum(1.2);
	refolded_h->SetMinimum(0.8);
	
	refolded_h->DrawCopy("HIST E1 P");
	
	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	
	canv_p->cd();
	pads[8]->cd();
	
	Double_t min = 10000000;
	Double_t max = -1.;
	TH2D* unfolding_h = NULL;
	if(isBayes){
	  unfolding_h = new TH2D(tempBayes->UnfoldingMatrix());
	//		else unfolding_h = new TH2D(tempSvd->UnfoldingMatrix());


	  for(Int_t bIX = 0; bIX < unfolding_h->GetNbinsX(); ++bIX){
	    for(Int_t bIY = 0; bIY < unfolding_h->GetNbinsY(); ++bIY){
	      
	      if(unfolding_h->GetBinContent(bIX+1, bIY+1) < min && unfolding_h->GetBinContent(bIX+1, bIY+1) > 0) min = unfolding_h->GetBinContent(bIX+1, bIY+1);
	      if(unfolding_h->GetBinContent(bIX+1, bIY+1) > max) max = unfolding_h->GetBinContent(bIX+1, bIY+1);
	      
	    }
	  }
	  
	  unfolding_h->SetMaximum(max*5.);
	  unfolding_h->SetMinimum(min/2.);
	  gStyle->SetPaintTextFormat("1.3f");
	  unfolding_h->SetMarkerSize(1.25);
	  unfolding_h->DrawCopy("COLZ TEXT");
	  gPad->SetLogz();
	  delete unfolding_h;
	}
	else{
	  TH1D* temp_h = (TH1D*)(tempSvd->Impl()->GetD());
	  temp_h->DrawCopy("HIST E1");
	  gPad->SetLogy();
	}

	canv_p->cd();
	pads[9]->cd();
	
	for(Int_t bIY = 0; bIY < responseSymm_h[cI][aI]->GetNbinsY(); ++bIY){
	  Double_t total = 0.0;
	  for(Int_t bIX = 0; bIX < responseSymm_h[cI][aI]->GetNbinsX(); ++bIX){
	    total += responseSymm_h[cI][aI]->GetBinContent(bIX+1, bIY+1);
	  }
	  
	  for(Int_t bIX = 0; bIX < responseSymm_h[cI][aI]->GetNbinsX(); ++bIX){
	    responseSymm_h[cI][aI]->SetBinContent(bIX+1, bIY+1, responseSymm_h[cI][aI]->GetBinContent(bIX+1, bIY+1)/total);
	    responseSymm_h[cI][aI]->SetBinError(bIX+1, bIY+1, responseSymm_h[cI][aI]->GetBinError(bIX+1, bIY+1)/total);	  
	  }
	}

	for(Int_t bIY = 0; bIY < responseAsymm_h[cI][aI]->GetNbinsY(); ++bIY){
	  Double_t total = 0.0;
	  for(Int_t bIX = 0; bIX < responseAsymm_h[cI][aI]->GetNbinsX(); ++bIX){
	    total += responseAsymm_h[cI][aI]->GetBinContent(bIX+1, bIY+1);
	  }
	  
	  for(Int_t bIX = 0; bIX < responseAsymm_h[cI][aI]->GetNbinsX(); ++bIX){
	    responseAsymm_h[cI][aI]->SetBinContent(bIX+1, bIY+1, responseAsymm_h[cI][aI]->GetBinContent(bIX+1, bIY+1)/total);
	    responseAsymm_h[cI][aI]->SetBinError(bIX+1, bIY+1, responseAsymm_h[cI][aI]->GetBinError(bIX+1, bIY+1)/total);	  
	  }
	}
	
	min = 100000;
	max = -1;
	for(Int_t bIX = 0; bIX < responseSymm_h[cI][aI]->GetNbinsX(); ++bIX){
	  for(Int_t bIY = 0; bIY < responseSymm_h[cI][aI]->GetNbinsY(); ++bIY){
	    if(min > responseSymm_h[cI][aI]->GetBinContent(bIX+1, bIY+1) && responseSymm_h[cI][aI]->GetBinContent(bIX+1, bIY+1) > 0) min = responseSymm_h[cI][aI]->GetBinContent(bIX+1, bIY+1);
	    if(max < responseSymm_h[cI][aI]->GetBinContent(bIX+1, bIY+1)) max = responseSymm_h[cI][aI]->GetBinContent(bIX+1, bIY+1);
	  }
	}
	
	responseSymm_h[cI][aI]->SetMinimum(min/2.);
	responseSymm_h[cI][aI]->SetMaximum(max*5.);
	gStyle->SetPaintTextFormat("1.3f");
	responseSymm_h[cI][aI]->SetMarkerSize(1.25);
	responseSymm_h[cI][aI]->DrawCopy("COLZ TEXT");
      
	gPad->RedrawAxis();
	gPad->SetLogz();
	//gPad->SetLogx();
	//	gPad->SetLogy();
	
	canv_p->cd();
	pads[10]->cd();

	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	
	TH2D* tempResponseSymm_h = (TH2D*)(rooResponse_h[cI][aI]->Hresponse()->Clone());
	tempResponseSymm_h->SetTitle("");
	for(Int_t bIY = 0; bIY < tempResponseSymm_h->GetNbinsY(); ++bIY){
	  Double_t total = 0.0;
	  for(Int_t bIX = 0; bIX < tempResponseSymm_h->GetNbinsX(); ++bIX){
	    total += tempResponseSymm_h->GetBinContent(bIX+1, bIY+1);
	  }
	  
	  for(Int_t bIX = 0; bIX < tempResponseSymm_h->GetNbinsX(); ++bIX){
	    tempResponseSymm_h->SetBinContent(bIX+1, bIY+1, tempResponseSymm_h->GetBinContent(bIX+1, bIY+1)/total);
	    tempResponseSymm_h->SetBinError(bIX+1, bIY+1, tempResponseSymm_h->GetBinError(bIX+1, bIY+1)/total);	  
	  }
	}
	
	
	tempResponseSymm_h->SetMinimum(min/2.);
	tempResponseSymm_h->SetMaximum(max*5.);
	gStyle->SetPaintTextFormat("1.3f"); 
	tempResponseSymm_h->SetMarkerSize(1.25);
	tempResponseSymm_h->DrawCopy("COLZ TEXT");
	gPad->RedrawAxis();
	gPad->SetLogz();
	//gPad->SetLogx();
	//	gPad->SetLogy();
	//      delete tempResponseSymm_h;
	
	
	canv_p->cd();
	pads[11]->cd();
	
	TMatrixD tempCovBayes = (TMatrixD)tempBayes->Ereco(RooUnfold::kCovToy);
	TMatrixD* pearsonCoefsBayes = (TMatrixD*)tempCovBayes.Clone("pearsonCoefsBayes)");

	TMatrixD tempCovSvd = (TMatrixD)tempSvd->Ereco(RooUnfold::kCovToy);
	TMatrixD* pearsonCoefsSvd = (TMatrixD*)tempCovSvd.Clone("pearsonCoefsSvd)");
	
	max = -1;
	min = 1000000;
	for(Int_t rI = 0; rI < pearsonCoefsBayes->GetNrows(); ++rI){
	  for(Int_t cI = 0; cI < pearsonCoefsBayes->GetNcols(); ++cI){
	    Double_t valBayes = tempCovBayes(rI,cI);
	    
	    if(tempCovBayes(rI,rI) > 0.00000001 && tempCovBayes(cI,cI) > 0.00000001) valBayes /= TMath::Sqrt(tempCovBayes(rI,rI)*tempCovBayes(cI,cI));
	    else if(valBayes < 0.00000001) valBayes = 0;
	    else{
	      std::cout << "UHOH WARNING: PEARSON BAYES DENOM IS 0 WITHOUT NUM BEING 0" << std::endl;
	      valBayes = 0.;
	    }
	    
	    if(valBayes > max) max = valBayes;
	    if(valBayes < min) min = valBayes;
	    
	    (*(pearsonCoefsBayes))(rI,cI) = valBayes;
	  }
	}

	for(Int_t rI = 0; rI < pearsonCoefsSvd->GetNrows(); ++rI){
	  for(Int_t cI = 0; cI < pearsonCoefsSvd->GetNcols(); ++cI){
	    Double_t valSvd = tempCovSvd(rI,cI);
	    
	    if(tempCovSvd(rI,rI) > 0.00000001 && tempCovSvd(cI,cI) > 0.00000001) valSvd /= TMath::Sqrt(tempCovSvd(rI,rI)*tempCovSvd(cI,cI));
	    else if(valSvd < 0.00000001) valSvd = 0;
	    else{
	      std::cout << "UHOH WARNING: PEARSON BAYES DENOM IS 0 WITHOUT NUM BEING 0" << std::endl;
	      valSvd = 0.;
	    }
	    
	    if(valSvd > max) max = valSvd;
	    if(valSvd < min) min = valSvd;
	    
	    (*(pearsonCoefsSvd))(rI,cI) = valSvd;
	  }
	}

	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	
	TH2D* covariancePearson = NULL;
	if(isBayes) covariancePearson = new TH2D(*pearsonCoefsBayes);
	else covariancePearson = new TH2D(*pearsonCoefsSvd);

	covariancePearson->SetMarkerSize(1.25);
	gStyle->SetPaintTextFormat("1.3f");
	covariancePearson->SetMaximum(max + (max - min)/10.);
	covariancePearson->SetMinimum(min - (max - min)/10.);
	covariancePearson->DrawCopy("COLZ");

	for(Int_t bIX = 0; bIX < covariancePearson->GetNbinsX(); ++bIX){
	  for(Int_t bIY = 0; bIY < covariancePearson->GetNbinsY(); ++bIY){
	    Double_t outVal = covariancePearson->GetBinContent(bIX+1, bIY+1);
	    if(TMath::Abs(outVal) > .1) label_p->DrawLatex(bIX+.2, bIY+.5, ("#bf{#color[2]{" + prettyString(outVal, 3, false) + "}}").c_str());
	    else label_p->DrawLatex(bIX+.1, bIY+.5, prettyString(outVal, 3, false).c_str());
	  }
	}
	
	delete covariancePearson;
	
	delete tempReco;
	delete tempBayes;
	delete tempSvd;
        
	std::string bayesStr = "Bayes" + std::to_string(bI);
	if(!isBayes) bayesStr = "Svd" + std::to_string(bI);

	canv_p->SaveAs(("pdfDir/build_R" + std::to_string(rVal) + "_" + totStr + "_" + bayesStr + "_" + std::to_string(date->GetDate()) + ".pdf").c_str());
	for(Int_t pI = 0; pI < nPads; ++pI){
	  delete pads[pI];
	}
      
	hvv->push_back(canv_p->GetPointer());

	delete canv_p;
      }
    }
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    for(Int_t aI = 0; aI < nAbsEtaBins; ++aI){
      reco_h[cI][aI]->Write("", TObject::kOverwrite);
      truth_h[cI][aI]->Write("", TObject::kOverwrite);
      recoPerp_h[cI][aI]->Write("", TObject::kOverwrite);
      truthPerp_h[cI][aI]->Write("", TObject::kOverwrite);
      responseSymm_h[cI][aI]->Write("", TObject::kOverwrite);
      responseAsymm_h[cI][aI]->Write("", TObject::kOverwrite);
      rooResponse_h[cI][aI]->Write("", TObject::kOverwrite);
      recoData_h[cI][aI]->Write("", TObject::kOverwrite);
    }
  }

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    for(Int_t aI = 0; aI < nAbsEtaBins; ++aI){
      delete reco_h[cI][aI];
      delete truth_h[cI][aI];
      delete recoPerp_h[cI][aI];
      delete truthPerp_h[cI][aI];
      delete responseSymm_h[cI][aI];
      delete responseAsymm_h[cI][aI];
      delete rooResponse_h[cI][aI];
      delete recoData_h[cI][aI];
    }
  }

  outFile_p->Close();
  delete outFile_p;
  
  delete randGen_p;
  delete date;

  delete label_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 4){
    std::cout << "Usage: ./buildAndTest.exe <inData> <inMC> <isPP>" << std::endl;
    return 1;
  }
  
  std::vector<std::vector<std::string>*>* hvv=new std::vector<std::vector<std::string>*>();

  const Int_t nRParam = 5;
  const Int_t rParam[nRParam] = {3, 4, 6, 8, 10};
  //  const Int_t rParam[nRParam] = {4};

  int retVal = 0;
  for(Int_t rI = 0; rI < nRParam; ++rI){
    retVal += buildAndTest(argv[1], argv[2], std::stoi(argv[3]), rParam[rI], true, hvv);
    retVal += buildAndTest(argv[1], argv[2], std::stoi(argv[3]), rParam[rI], false, hvv);
  }

  std::cout<<"Generating slides"<<std::endl;
  TexSlides(hvv,"Slides.tex",1);

  return retVal;
}
