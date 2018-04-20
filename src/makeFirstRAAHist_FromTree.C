#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <map>

#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TDatime.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TBox.h"
#include "TMath.h"
#include "TNamed.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMatrixD.h"
#include "TRandom3.h"
#include "TSVDUnfold.h"

#include "src/RooUnfoldResponse.h"
#include "src/RooUnfoldBayes.h"
#include "src/RooUnfoldSvd.h"

#include "include/jetData.h"
#include "include/cppWatch.h"
#include "include/checkMakeDir.h"
#include "include/returnRootFileContentsList.h"
#include "include/doGlobalDebug.h"
#include "include/getLogBins.h"
#include "include/getLinBins.h"
#include "include/kirchnerPalette.h"
#include "include/jtAlgoCentBins.h"
#include "include/histDefUtility.h"
#include "include/plotUtilities.h"
#include "include/scaleErrorTool.h"
#include "include/smearingFuncs.h"

const Float_t mcFraction = 0.9;
const Bool_t doSvd = true;

const Int_t nPadTotValid = 12;

const Int_t nJER = 10;

const Int_t nCentBins = 4;
const Int_t centBinsLow[nCentBins] = {50, 30, 10, 0};
const Int_t centBinsHi[nCentBins] = {90, 50, 30, 10};

const Int_t nAbsEtaBins = 5;
const Float_t absEtaBinsLow[nAbsEtaBins] = {0.0, 0.5, 1.0, 1.5, 0.0};
const Float_t absEtaBinsHi[nAbsEtaBins] = {0.5, 1.0, 1.5, 2.0, 2.0};
const Float_t absEtaBins[nAbsEtaBins] = {0.0, 0.5, 1.0, 1.5, 2.0};

const Float_t jtPtLowGlobal = 75.;
const Float_t jtPtLowPP = 110.;
const Float_t jtPtLowPbPb = 140.;
const Float_t jtAbsEtaMax = 2.;


const Int_t nSyst = 11;
const std::string systStr[nSyst] = {"Nominal", "JECVarUpData", "JECVarDownData", "JECVarUpMC", "JECVarDownMC", "JER1p15", "LumiUp", "LumiDown", "TAAUp", "TAADown", "Fake"};
const Int_t nSystForFill = 6;
std::string systStrForFill[nSystForFill] = {"Nominal", "JECVarUpData", "JECVarDownData", "JECVarUpMC", "JECVarDownMC", "JER1p15"};
const Int_t systPosForFill[nSystForFill] = {0, 1, 2, 3, 4, 5};
const Int_t nSystForCopy = 5;
std::string systStrForCopy[nSystForCopy] = {"LumiUp", "LumiDown", "TAAUp", "TAADown", "Fake"};
const Int_t systPosForCopy[nSystForCopy] = {6, 7, 8, 9, 10};

/*
const Int_t nSyst = 2;
const std::string systStr[nSyst] = {"Nominal", "JER1p15"};
const Int_t nSystForFill = 2;
std::string systStrForFill[nSystForFill] = {"Nominal", "JER1p15"};
const Int_t systPosForFill[nSystForFill] = {0,1};
const Int_t nSystForCopy = 0;
std::string systStrForCopy[nSystForCopy] = {};
const Int_t systPosForCopy[nSystForCopy] = {};
*/

void unfoldBayes(RooUnfoldResponse* rooRes_p, TH1D* jtPt_RAW_p, const Int_t bayesI, TH1D* jtPt_FULL_p, RooUnfoldBayes* inBayes)
{
  inBayes->Setup(rooRes_p, jtPt_RAW_p);
  inBayes->SetIterations(bayesI);
  inBayes->SetSmoothing(false);
  inBayes->SetVerbose(0);

  std::string cloneName = inBayes->GetName();
  cloneName = cloneName + "_TEMP";
  TH1D* hreco = (TH1D*)(inBayes->Hreco(RooUnfold::kCovToy)->Clone(cloneName.c_str()));

  for(Int_t bI = 0; bI < hreco->GetNbinsX(); ++bI){
    jtPt_FULL_p->SetBinContent(bI+1, hreco->GetBinContent(bI+1));
    jtPt_FULL_p->SetBinError(bI+1, hreco->GetBinError(bI+1));
  }
  delete hreco;

  return;
}

void unfoldSvd(RooUnfoldResponse* rooRes_p, TH1D* jtPt_RAW_p, const Int_t kreg, TH1D* jtPt_FULL_p, RooUnfoldSvd* inSvd)
{
  inSvd->Setup(rooRes_p, jtPt_RAW_p);
  inSvd->SetKterm(kreg);
  inSvd->SetVerbose(0);

  std::string cloneName = inSvd->GetName();
  cloneName = cloneName + "_TEMP";
  TH1D* hreco = (TH1D*)(inSvd->Hreco(RooUnfold::kCovToy)->Clone(cloneName.c_str()));

  for(Int_t bI = 0; bI < hreco->GetNbinsX(); ++bI){
    jtPt_FULL_p->SetBinContent(bI+1, hreco->GetBinContent(bI+1));
    jtPt_FULL_p->SetBinError(bI+1, hreco->GetBinError(bI+1));
  }

  delete hreco;

  return;
}


void initTH1(TH1* inHist_p, Int_t style, Double_t size, Int_t colPos)
{
  kirchnerPalette col;

  inHist_p->SetMarkerStyle(style);
  inHist_p->SetMarkerSize(size);
  inHist_p->SetMarkerColor(col.getColor(colPos));
  inHist_p->SetLineColor(col.getColor(colPos));

  inHist_p->GetXaxis()->SetTitleFont(43);
  inHist_p->GetYaxis()->SetTitleFont(43);
  inHist_p->GetXaxis()->SetLabelFont(43);
  inHist_p->GetYaxis()->SetLabelFont(43);

  inHist_p->GetXaxis()->SetTitleSize(10);
  inHist_p->GetYaxis()->SetTitleSize(10);
  inHist_p->GetXaxis()->SetLabelSize(10);
  inHist_p->GetYaxis()->SetLabelSize(10);

  inHist_p->SetTitle("");
 
  return;
}


double getComboBinVal(TH1* inHist1_p, TH1* inHist2_p, Int_t binPos)
{
  double combo = inHist1_p->GetBinContent(binPos) + inHist2_p->GetBinContent(binPos);
  return combo;
}

double getComboBinErr(TH1* inHist1_p, TH1* inHist2_p, Int_t binPos)
{
  double err = TMath::Sqrt(inHist1_p->GetBinError(binPos)*inHist1_p->GetBinError(binPos) + inHist2_p->GetBinError(binPos)*inHist2_p->GetBinError(binPos));
  return err;
}

void doHistCopy(TH1* inHist_p, TH1* outHist_p)
{
  for(Int_t bIX = 0; bIX < inHist_p->GetNbinsX(); ++bIX){
    outHist_p->SetBinContent(bIX+1, inHist_p->GetBinContent(bIX+1));
    outHist_p->SetBinError(bIX+1, inHist_p->GetBinError(bIX+1));
  }

  return;
}

int getRVal(const std::string inStr)
{
  Int_t rVal = -1;
  if(inStr.find("ak1PF") != std::string::npos) rVal = 1;
  else if(inStr.find("ak2PF") != std::string::npos) rVal = 2;
  else if(inStr.find("ak3PF") != std::string::npos) rVal = 3;
  else if(inStr.find("ak4PF") != std::string::npos) rVal = 4;
  else if(inStr.find("ak5PF") != std::string::npos) rVal = 5;
  else if(inStr.find("ak6PF") != std::string::npos) rVal = 6;
  else if(inStr.find("ak8PF") != std::string::npos) rVal = 8;
  else if(inStr.find("ak10PF") != std::string::npos) rVal = 10;
  else if(inStr.find("akCs1P") != std::string::npos) rVal = 1;
  else if(inStr.find("akCs2P") != std::string::npos) rVal = 2;
  else if(inStr.find("akCs3P") != std::string::npos) rVal = 3;
  else if(inStr.find("akCs4P") != std::string::npos) rVal = 4;
  else if(inStr.find("akCs5P") != std::string::npos) rVal = 5;
  else if(inStr.find("akCs6P") != std::string::npos) rVal = 6;
  else if(inStr.find("akCs8P") != std::string::npos) rVal = 8;
  else if(inStr.find("akCs10P") != std::string::npos) rVal = 10;

  return rVal;
}

void doTreeFill(TTree* inTree_p, TTree* outTree_p[nCentBins], jetData* jDataOut[nCentBins][nSystForFill], const Bool_t isMC, const Bool_t isPP, const std::string algoStr, const unsigned int tI, TH1D* jtpt_RAW_h[nCentBins][nAbsEtaBins+1][nSyst], TH1D* refpt_h[nCentBins][nAbsEtaBins+1], TH2D* response_h[nCentBins][nAbsEtaBins+1], TH2D* response_Symm_h[nCentBins][nAbsEtaBins+1], RooUnfoldResponse* rooResponse_h[nCentBins][nAbsEtaBins+1], TH1D* jtpt_RAW_Frac_h[nCentBins][nAbsEtaBins+1][nSyst], TH1D* refpt_Frac_h[nCentBins][nAbsEtaBins+1])
{
  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  TRandom3* randGen_p = new TRandom3(0);

  scaleErrorTool scaleErr("/afs/cern.ch/work/c/cmcginn/private/Projects/CMSRooJR/input/rcDifferences_20180406.txt");
  scaleErr.Init();
  jtAlgoCentBins binner;

  Int_t truthNBins[nCentBins];
  Int_t recoNBins[nCentBins];
  Float_t truthBins[nCentBins][100];
  Float_t recoBins[nCentBins][100];

  Int_t maxTruthBins = -1;

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  const Int_t rVal = getRVal(inTree_p->GetName());
  Float_t jtPtLow = jtPtLowPbPb;
  if(isPP) jtPtLow = jtPtLowPP;
  
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    recoNBins[cI] = binner.getNBins(centBinsLow[cI], centBinsHi[cI], rVal, true);
    truthNBins[cI] = binner.getNBins(centBinsLow[cI], centBinsHi[cI], rVal, false);
    binner.getBins(centBinsLow[cI], centBinsHi[cI], rVal, true, recoNBins[cI], recoBins[cI]);
    binner.getBins(centBinsLow[cI], centBinsHi[cI], rVal, false, truthNBins[cI], truthBins[cI]);
    
    if(truthNBins[cI] > maxTruthBins) maxTruthBins = truthNBins[cI];
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  const Int_t nEntries = inTree_p->GetEntries();
  jetData jData;

  //  inTree_p->SetBranchStatus("*", 0);
  if(!isMC && !isPP) jData.SetStatusAndAddressRead(inTree_p, {"hiBin", "nref", "jtpt", "rawpt", "jteta"});
  else if(isMC && !isPP) jData.SetStatusAndAddressRead(inTree_p, {"hiBin", "fullWeight", "pthatWeight", "nref", "jtpt", "rawpt", "refpt", "jteta"});
  else if(!isMC && isPP) jData.SetStatusAndAddressRead(inTree_p, {"nref", "jtpt", "rawpt", "jteta"});
  else if(isMC && isPP) jData.SetStatusAndAddressRead(inTree_p, {"fullWeight", "pthatWeight", "nref", "jtpt", "rawpt", "refpt", "jteta"});

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%100000 == 0) std::cout << " Entry " << entry << "/" << nEntries << " (" << inTree_p->GetName() << ")" << std::endl;

    inTree_p->GetEntry(entry);

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    std::vector<Int_t> centPos;
    if(!isPP){
      for(Int_t cI = 0; cI < nCentBins; ++cI){
	if(centBinsLow[cI]*2 <= jData.hiBin_ && jData.hiBin_ < centBinsHi[cI]*2){
	  centPos.push_back(cI);
	  break;
	}
      }
      if(centPos.size() == 0) continue;
    }
    else{
      for(Int_t cI = 0; cI < nCentBins; ++cI){
	centPos.push_back(cI);
      }
    }

    Float_t tempWeight = 1.;
    if(isMC) tempWeight = jData.fullWeight_;
    Bool_t isFracEvent = randGen_p->Uniform(0.0, 1.0) > mcFraction;
    Int_t scaleHiBin = 160;
    Int_t jerHiBin = 0;
    if(!isPP){
      scaleHiBin = jData.hiBin_/2;
      jerHiBin = jData.hiBin_/2;
    }

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    for(unsigned int cI = 0; cI < centPos.size(); ++cI){
      Int_t cI2 = centPos.at(cI);
      jDataOut[cI2][0]->nref_ = 0;
    }

    for(Int_t jI = 0; jI < jData.nref_; ++jI){
      if(TMath::Abs(jData.jteta_[jI]) > jtAbsEtaMax) continue;
      if(jData.jtpt_[jI] < jtPtLowGlobal) break;

      std::vector<Int_t> absPos;
      
      for(Int_t aI = 0; aI < nAbsEtaBins; ++aI){
	if(TMath::Abs(jData.jteta_[jI]) >= absEtaBinsLow[aI] && TMath::Abs(jData.jteta_[jI]) < absEtaBinsHi[aI]){
	  absPos.push_back(aI);
	}
      }
      
      Bool_t goodNominalPt = jData.jtpt_[jI] >= jtPtLow;

      for(Int_t sI = 0; sI < nSystForFill; ++sI){
	Float_t jtptTemp = jData.jtpt_[jI];
	Float_t rawptTemp = jData.rawpt_[jI];
	std::vector<Float_t> jerVect;

	if(systStr[sI].find("JECVarUpData") != std::string::npos) jtptTemp *= 1.02;
	else if(systStr[sI].find("JECVarDownData") != std::string::npos) jtptTemp *= 0.98;
	else if(systStr[sI].find("JECVarUpMC") != std::string::npos){
	  Float_t tempScale = jtptTemp/rawptTemp;
	  rawptTemp += TMath::Abs(scaleErr.getMuDataMinusMC(scaleHiBin, jData.jteta_[jI], rVal, "FlowDefaultInRho"));
	  jtptTemp = tempScale*rawptTemp;
	}
	else if(systStr[sI].find("JECVarDownMC") != std::string::npos){
	  Float_t tempScale = jtptTemp/rawptTemp;
	  rawptTemp -= TMath::Abs(scaleErr.getMuDataMinusMC(scaleHiBin, jData.jteta_[jI], rVal, "FlowDefaultInRho"));
	  jtptTemp = tempScale*rawptTemp;
	}
	else if(systStr[sI].find("JER1p15") != std::string::npos){
	  for(Int_t jerI = 0; jerI < nJER; ++jerI){
	    jerVect.push_back(jtptTemp*randGen_p->Gaus(1., getResForPtEtaCentAlgo(jtptTemp, jData.jteta_[jI], jerHiBin, algoStr, true)));
	  }
	}

	Bool_t goodNominal = systStr[sI].find("Nominal") != std::string::npos;

	if(goodNominalPt){
	  for(unsigned int cI = 0; cI < centPos.size(); ++cI){
	    Int_t cI2 = centPos.at(cI);
	    if(goodNominal){
	      jDataOut[cI2][0]->jteta_[jDataOut[cI2][0]->nref_] = jData.jteta_[jI];
	    }

	    jDataOut[cI2][sI]->jtpt_[jDataOut[cI2][0]->nref_] = jtptTemp;
	    if(isMC) jDataOut[cI2][sI]->refpt_[jDataOut[cI2][0]->nref_] = jData.refpt_[jI];
	    jDataOut[cI2][sI]->rawpt_[jDataOut[cI2][0]->nref_] = rawptTemp;
	  }
	}
     
	Float_t fillWeight = tempWeight;

	if(jerVect.size() == 0) jerVect.push_back(jtptTemp);
	else{
	  fillWeight *= (1./(float)nJER);
	  //	  if(!isMC) std::cout << "TEMP WEIGHT: " << fillWeight << ", " << nJER << ", " << std::endl;
	}

	for(unsigned int jerI = 0; jerI < jerVect.size(); ++jerI){
	  if(jerVect.at(jerI) < jtPtLow) continue;	      

	  
	  for(unsigned int cI = 0; cI < centPos.size(); ++cI){	
	    Int_t cI2 = centPos.at(cI);
	    Bool_t goodReco = jerVect.at(jerI) >= recoBins[cI2][0] && jerVect.at(jerI) < recoBins[cI2][recoNBins[cI2]];
	    Bool_t goodTruth = isMC;
	    if(goodTruth) goodTruth = jData.refpt_[jI] >= truthBins[cI2][0] && jData.refpt_[jI] < truthBins[cI2][truthNBins[cI2]];

	    if(goodNominal){
	      jDataOut[cI2][0]->jteta_[jDataOut[cI2][0]->nref_] = jData.jteta_[jI];
	    }

	    if(jerI == 0){
	      jDataOut[cI2][sI]->jtpt_[jDataOut[cI2][0]->nref_] = jerVect.at(jerI);
	      jDataOut[cI2][sI]->rawpt_[jDataOut[cI2][0]->nref_] = rawptTemp;
	    }
	    
	    for(unsigned int aI = 0; aI < absPos.size(); ++aI){
	      if(goodReco && (goodTruth || !isMC)){
		if(isFracEvent && isMC){
		  jtpt_RAW_Frac_h[cI2][absPos.at(aI)][systPosForFill[sI]]->Fill(jerVect.at(jerI), fillWeight);
		  if(goodNominal) refpt_Frac_h[cI2][absPos.at(aI)]->Fill(jData.refpt_[jI], fillWeight);
		}
		else{
		  jtpt_RAW_h[cI2][absPos.at(aI)][systPosForFill[sI]]->Fill(jerVect.at(jerI), fillWeight);
		  if(goodNominal && isMC){
		    refpt_h[cI2][absPos.at(aI)]->Fill(jData.refpt_[jI], fillWeight);
		    response_h[cI2][absPos.at(aI)]->Fill(jerVect.at(jerI), jData.refpt_[jI], fillWeight);
		    response_Symm_h[cI2][absPos.at(aI)]->Fill(jerVect.at(jerI), jData.refpt_[jI], fillWeight);
		    rooResponse_h[cI2][absPos.at(aI)]->Fill(jerVect.at(jerI), jData.refpt_[jI], fillWeight);	      
		  }
		}
	      }
	      else if(goodTruth && isMC && goodNominal){
		if(isFracEvent){
		  refpt_Frac_h[cI2][absPos.at(aI)]->Fill(jData.refpt_[jI], fillWeight);
		}
		else{
		  response_Symm_h[cI2][absPos.at(aI)]->Fill(jerVect.at(jerI), jData.refpt_[jI], fillWeight);
		  rooResponse_h[cI2][absPos.at(aI)]->Miss(jData.refpt_[jI], fillWeight);
		  refpt_h[cI2][absPos.at(aI)]->Fill(jData.refpt_[jI], fillWeight);	    
		}
	      }
	    }
	  }		
	}

	jerVect.clear();
      }

      if(goodNominalPt){
	for(unsigned int cI = 0; cI < centPos.size(); ++cI){
	  Int_t cI2 = centPos.at(cI);
	  jDataOut[cI2][0]->nref_ += 1;
	}
      }
    }
  
    for(unsigned int cI = 0; cI < centPos.size(); ++cI){
      Int_t cI2 = centPos.at(cI);
      if(!isPP) jDataOut[cI2][0]->hiBin_ = jData.hiBin_;
      if(isMC) jDataOut[cI2][0]->pthatWeight_ = jData.pthatWeight_;
      if(isMC) jDataOut[cI2][0]->fullWeight_ = jData.fullWeight_;

      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << ", " << cI2 << std::endl;
      if(doGlobalDebug){
	if(!isPP) std::cout << jDataOut[cI2][0]->hiBin_ << std::endl;
	if(isMC) std::cout << jDataOut[cI2][0]->pthatWeight_ << std::endl;
	if(isMC) std::cout << jDataOut[cI2][0]->fullWeight_ << std::endl;
	std::cout << jDataOut[cI2][0]->nref_ << std::endl;
      }

      if(randGen_p->Uniform(0,10) <= 1) outTree_p[cI2]->Fill();
      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    }
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  delete randGen_p;

  return;
}

void makeEtaHist(TH1D* jteta_h, std::vector<TH1D*> jtpt_h)
{
  for(Int_t aI = 0; aI < jteta_h->GetNbinsX(); ++aI){
    Double_t val = 0;
    Double_t valErr = 0;
    
    for(Int_t bIX = 0; bIX < jtpt_h.at(aI)->GetNbinsX(); ++bIX){
      if(jtpt_h.at(aI)->GetBinCenter(bIX+1) > 400.) break;
      if(jtpt_h.at(aI)->GetBinCenter(bIX+1) > 200.){
	val += jtpt_h.at(aI)->GetBinContent(bIX+1);	
	valErr += jtpt_h.at(aI)->GetBinError(bIX+1)*jtpt_h.at(aI)->GetBinError(bIX+1);
      }
    }
    
    jteta_h->SetBinContent(aI+1, val);
    jteta_h->SetBinError(aI+1, TMath::Sqrt(valErr));
  }

  return;
}


double getAbsEtaDelta(const std::string inStr)
{
  double retVal = -1.;
  if(inStr.find("AbsEta0p0to2p0") != std::string::npos) retVal = 2.;
  else if(inStr.find("AbsEtaTot") != std::string::npos) retVal = 2.;
  else if(inStr.find("AbsEta0p0to0p5") != std::string::npos) retVal = 0.5;
  else if(inStr.find("AbsEta0p5to1p0") != std::string::npos) retVal = 0.5;
  else if(inStr.find("AbsEta1p0to1p5") != std::string::npos) retVal = 0.5;
  else if(inStr.find("AbsEta1p5to2p0") != std::string::npos) retVal = 0.5;
  else{
    std::cout << "WARNING: \'" << inStr << "\' is not found, return -1" << std::endl;
  }
  return retVal;
}

double getCentDelta(const std::string inStr)
{
  double retVal = -1.;
  if(inStr.find("Cent0to10") != std::string::npos) retVal = 10.;
  else if(inStr.find("Cent10to30") != std::string::npos) retVal = 20.;
  else if(inStr.find("Cent30to50") != std::string::npos) retVal = 20.;
  else if(inStr.find("Cent50to90") != std::string::npos) retVal = 40.;
  else{
    std::cout << "WARNING: \'" << inStr << "\' is not found, return -1" << std::endl;
  }
  return retVal;
}

double getTAAScaleFactor(const std::string inStr)
{
  double scaleFactor = -1.;
  if(inStr.find("Cent0to10") != std::string::npos) scaleFactor = 23.22;
  else if(inStr.find("Cent10to30") != std::string::npos) scaleFactor = 11.51;
  else if(inStr.find("Cent30to50") != std::string::npos) scaleFactor = 3.819;
  else if(inStr.find("Cent50to90") != std::string::npos) scaleFactor = 0.543 ;
  else{
    std::cout << "WARNING: \'" << inStr << "\' is not found, return -1" << std::endl;
  }
  return scaleFactor;
}

double getTAAScaleFactorUp(const std::string inStr)
{
  double scaleFactor = -1.;
  if(inStr.find("Cent0to10") != std::string::npos) scaleFactor = .019;
  else if(inStr.find("Cent10to30") != std::string::npos) scaleFactor = .026;
  else if(inStr.find("Cent30to50") != std::string::npos) scaleFactor = .054;
  else if(inStr.find("Cent50to90") != std::string::npos) scaleFactor = .112;
  else{
    std::cout << "WARNING: \'" << inStr << "\' is not found, return -1" << std::endl;
  }
  return scaleFactor;
}

double getTAAScaleFactorDown(const std::string inStr)
{
  double scaleFactor = -1.;
  if(inStr.find("Cent0to10") != std::string::npos) scaleFactor = .03;
  else if(inStr.find("Cent10to30") != std::string::npos) scaleFactor = .034;
  else if(inStr.find("Cent30to50") != std::string::npos) scaleFactor = .054;
  else if(inStr.find("Cent50to90") != std::string::npos) scaleFactor = .073;
  else{
    std::cout << "WARNING: \'" << inStr << "\' is not found, return -1" << std::endl;
  }
  return scaleFactor;
}

void divBinWidth(TH1* inHist_p)
{
  for(Int_t bIX = 0; bIX < inHist_p->GetNbinsX(); ++bIX){
    Double_t val = inHist_p->GetBinContent(bIX+1);
    Double_t err = inHist_p->GetBinError(bIX+1);
    val /= inHist_p->GetBinWidth(bIX+1);
    err /= inHist_p->GetBinWidth(bIX+1);
    inHist_p->SetBinContent(bIX+1, val);
    inHist_p->SetBinError(bIX+1, err);
  }

  return;
}

void setRAABins(TH1* inRAA_p, TH1* inPbPb_p, TH1* inPP_p)
{
  for(Int_t bIX = 0; bIX < inRAA_p->GetNbinsX(); ++bIX){
    Double_t val = inPbPb_p->GetBinContent(bIX+1)/inPP_p->GetBinContent(bIX+1);
    Double_t errPbPb = inPbPb_p->GetBinError(bIX+1)/inPbPb_p->GetBinContent(bIX+1);
    Double_t errPP = inPP_p->GetBinError(bIX+1)/inPP_p->GetBinContent(bIX+1);
    Double_t err = val*TMath::Sqrt(errPbPb*errPbPb + errPP*errPP);

    inRAA_p->SetBinContent(bIX+1, val);
    inRAA_p->SetBinError(bIX+1, err);
  }

  return;
}


void dumpTiming(std::vector<cppWatch> watches, std::vector<std::string> names)
{
  std::cout << "DUMP TIMING" << std::endl;
  Double_t total = 0.;
  for(unsigned int tI = 0; tI < watches.size(); ++tI){
    total += watches.at(tI).total();
  }

  for(unsigned int tI = 0; tI < watches.size(); ++tI){
    std::cout << " Timer " << names.at(tI) << ": " << watches.at(tI).total() << " (" << watches.at(tI).total()/total << ")" << std::endl;
  }

  return;
}


int makeFirstRAAHist_FromTree(const std::string inFileNamePbPb, const std::string inFileNamePP, const std::string inFileNamePbPbMC, const std::string inFileNamePPMC, const std::string inAlgoName)
{
  cppWatch totalWatch;
  if(!checkFile(inFileNamePbPb) || inFileNamePbPb.find(".root") == std::string::npos){
    std::cout << "inFileNamePbPb \'" << inFileNamePbPb << "\' is invalid return 1" << std::endl;
    return 1;
  }
  if(!checkFile(inFileNamePP) || inFileNamePP.find(".root") == std::string::npos){
    std::cout << "inFileNamePP \'" << inFileNamePP << "\' is invalid return 1" << std::endl;
    return 1;
  }

  cppWatch startWatch;
  cppWatch histWatch;
  cppWatch pbpbDataWatch;
  cppWatch pbpbMCWatch;
  cppWatch ppDataWatch;
  cppWatch ppMCWatch;
  cppWatch endWatch;
  cppWatch countWatch;
  cppWatch unfoldLoopWatch;
  cppWatch plotValidWatch;
  cppWatch copyWatch;
  cppWatch finalEtaWatch;
  cppWatch makeEtaHistWatch;
  cppWatch rescaleWatch;
  cppWatch raaWatch;
  cppWatch writeWatch;
  cppWatch deleteWatch;
  cppWatch closeWatch;
  cppWatch bayesWatch;
  cppWatch svdWatch;

  startWatch.start();

  TDatime* date = new TDatime();
  std::string dirName1 = "pdfDir/";
  std::string dirName2 = dirName1 + std::to_string(date->GetDate()) + "/";
  checkMakeDir(dirName1);
  checkMakeDir(dirName2);

  TRandom3* randGen_p = new TRandom3(0);
  kirchnerPalette colors;
  jtAlgoCentBins binner;
  scaleErrorTool scaleErr("/afs/cern.ch/work/c/cmcginn/private/Projects/CMSRooJR/input/rcDifferences_20180322.txt");
  scaleErr.Init();

  std::map<std::string, std::string> jetMapPbPbToPP;
  jetMapPbPbToPP["akCs3P"] = "ak3PF";
  jetMapPbPbToPP["akCs4P"] = "ak4PF";
  jetMapPbPbToPP["akCs6P"] = "ak6PF";
  jetMapPbPbToPP["akCs8P"] = "ak8PF";
  jetMapPbPbToPP["akCs10P"] = "ak10PF";

  TFile* inFilePbPb_p = new TFile(inFileNamePbPb.c_str(), "READ");
  std::vector<std::string> jetTreeListPbPb = returnRootFileContentsList(inFilePbPb_p, "TTree", "JetA");

  TFile* inFilePP_p = new TFile(inFileNamePP.c_str(), "READ");
  std::vector<std::string> jetTreeListPPInit = returnRootFileContentsList(inFilePP_p, "TTree", "JetA");

  TFile* inFilePbPbMC_p = new TFile(inFileNamePbPbMC.c_str(), "READ");
  std::vector<std::string> jetTreeListPbPbMCInit = returnRootFileContentsList(inFilePbPbMC_p, "TTree", "JetA");

  TFile* inFilePPMC_p = new TFile(inFileNamePPMC.c_str(), "READ");
  std::vector<std::string> jetTreeListPPMCInit = returnRootFileContentsList(inFilePPMC_p, "TTree", "JetA");

  unsigned int pos = 0;
  while(jetTreeListPbPb.size() > pos){
    if(jetTreeListPbPb.at(pos).find(inAlgoName) != std::string::npos) ++pos;
    else jetTreeListPbPb.erase(jetTreeListPbPb.begin()+pos);
  }

  std::vector<std::string> jetTreeListPP;
  std::vector<std::string> jetTreeListPbPbMC;
  std::vector<std::string> jetTreeListPPMC;

  pos = 0;
  while(pos < jetTreeListPbPb.size()){
    int isFoundPP = -1;
    int isFoundPbPbMC = -1;
    int isFoundPPMC = -1;

    std::string jetName = jetTreeListPbPb.at(pos).substr(0, jetTreeListPbPb.at(pos).find("/"));

    for(auto it = jetMapPbPbToPP.begin(); it != jetMapPbPbToPP.end(); ++it){
      if(jetTreeListPbPb.at(pos).find(it->first) != std::string::npos){
	
	for(unsigned int i = 0; i < jetTreeListPPInit.size(); ++i){
	  if(jetTreeListPPInit.at(i).find(it->second) != std::string::npos){
	    isFoundPP = i;
	    break;
	  }
	}	
	break;
      }
    }
	
    for(unsigned int i = 0; i < jetTreeListPbPbMCInit.size(); ++i){
      if(jetTreeListPbPbMCInit.at(i).find(jetName) != std::string::npos){
	isFoundPbPbMC = i;
	break;
      }
    }	

    for(auto it = jetMapPbPbToPP.begin(); it != jetMapPbPbToPP.end(); ++it){
      if(jetTreeListPbPb.at(pos).find(it->first) != std::string::npos){
	
	for(unsigned int i = 0; i < jetTreeListPPMCInit.size(); ++i){
	  if(jetTreeListPPMCInit.at(i).find(it->second) != std::string::npos){
	    isFoundPPMC = i;
	    break;
	  }
	}	
	break;
      }
    }

    if(isFoundPP >= 0 && isFoundPbPbMC >= 0 && isFoundPPMC >= 0){
      jetTreeListPP.push_back(jetTreeListPPInit.at(isFoundPP));
      jetTreeListPbPbMC.push_back(jetTreeListPbPbMCInit.at(isFoundPbPbMC));
      jetTreeListPPMC.push_back(jetTreeListPPMCInit.at(isFoundPPMC));
      ++pos;
    }
    else{
      if(!isFoundPP){
	std::cout << "WARNING: Cannot find partner for PbPb Jet Algo \'" << jetTreeListPbPb.at(pos) << "\' from possible PP set: ";
	for(unsigned int j = 0; j < jetTreeListPP.size(); ++j){
	  std::cout << jetTreeListPP.at(j) << ",";
	}
	std::cout << std::endl;
      }
      if(!isFoundPbPbMC){
	std::cout << "WARNING: Cannot find partner for PbPb Jet Algo \'" << jetTreeListPbPb.at(pos) << "\' from possible PbPb MC set: ";
	for(unsigned int j = 0; j < jetTreeListPbPbMC.size(); ++j){
	  std::cout << jetTreeListPbPbMC.at(j) << ",";
	}
	std::cout << std::endl;
      }
      if(!isFoundPPMC){
	std::cout << "WARNING: Cannot find partner for PbPb Jet Algo \'" << jetTreeListPbPb.at(pos) << "\' from possible PP MC set: ";
	for(unsigned int j = 0; j < jetTreeListPPMC.size(); ++j){
	  std::cout << jetTreeListPPMC.at(j) << ",";
	}
	std::cout << std::endl;
      }
      
      std::cout << " Removing algo" << std::endl;      
      jetTreeListPbPb.erase(jetTreeListPbPb.begin()+pos);
    }
  }

  std::cout << "Processing PbPb TTrees..." << std::endl;
  for(unsigned int i = 0; i < jetTreeListPbPb.size(); ++i){
    std::cout << " " << i << "/" << jetTreeListPbPb.size() << ": " << jetTreeListPbPb.at(i) << std::endl;
  }
  std::cout << "Processing PP TTrees..." << std::endl;
  for(unsigned int i = 0; i < jetTreeListPP.size(); ++i){
    std::cout << " " << i << "/" << jetTreeListPP.size() << ": " << jetTreeListPP.at(i) << std::endl;
  }
  std::cout << "Processing PbPb MC TTrees..." << std::endl;
  for(unsigned int i = 0; i < jetTreeListPbPbMC.size(); ++i){
    std::cout << " " << i << "/" << jetTreeListPbPbMC.size() << ": " << jetTreeListPbPbMC.at(i) << std::endl;
  }
  std::cout << "Processing PP MC TTrees..." << std::endl;
  for(unsigned int i = 0; i < jetTreeListPPMC.size(); ++i){
    std::cout << " " << i << "/" << jetTreeListPPMC.size() << ": " << jetTreeListPPMC.at(i) << std::endl;
  }
  
  if(jetTreeListPbPb.size() == 0){
    std::cout << "No valid algos, return 1" << std::endl;
    return 1;
  }

  Int_t rVals = getRVal(jetTreeListPP.at(0));

  Int_t truthNBins[nCentBins];
  Int_t recoNBins[nCentBins];
  Float_t truthBins[nCentBins][100];
  Float_t recoBins[nCentBins][100];

  Int_t maxTruthBins = -1;
  
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    recoNBins[cI] = binner.getNBins(centBinsLow[cI], centBinsHi[cI], rVals, true);
    truthNBins[cI] = binner.getNBins(centBinsLow[cI], centBinsHi[cI], rVals, false);
    binner.getBins(centBinsLow[cI], centBinsHi[cI], rVals, true, recoNBins[cI], recoBins[cI]);
    binner.getBins(centBinsLow[cI], centBinsHi[cI], rVals, false, truthNBins[cI], truthBins[cI]);
    
    if(truthNBins[cI] > maxTruthBins) maxTruthBins = truthNBins[cI];
  }

  startWatch.stop();
  dumpTiming({startWatch}, {"start"});
  histWatch.start();

  const Int_t nBayesIter = 9;
  const Int_t nSvd = maxTruthBins-1;

  const std::string dateTimeStr = std::to_string(date->GetDate()) + "_" + std::to_string(date->GetHour());
  TFile* outFile_p = new TFile(("output/firstRAAHist_nBayes" + std::to_string(nBayesIter) + "_DoSvd" + std::to_string(doSvd) + "_nSvd" + std::to_string(nSvd) + "_" + inAlgoName + "_" + dateTimeStr + ".root").c_str(), "RECREATE");

  TDirectory* jetPbPbDirs_p;
  TDirectory* jetPPDirs_p;

  TDirectory* centPbPbDirs_p[nCentBins];
  TDirectory* centPPDirs_p[nCentBins];

  TTree* treeCentSystCheckPbPbData_p[nCentBins];
  TTree* treeCentSystCheckPPData_p[nCentBins];
  TTree* treeCentSystCheckPbPbMC_p[nCentBins];
  TTree* treeCentSystCheckPPMC_p[nCentBins];

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  jetData* jDataOutPbPbData[nCentBins][nSystForFill];
  jetData* jDataOutPPData[nCentBins][nSystForFill];
  jetData* jDataOutPbPbMC[nCentBins][nSystForFill];
  jetData* jDataOutPPMC[nCentBins][nSystForFill];


  TH1D* jtptPbPbData_RAW_h[nCentBins][nAbsEtaBins+1][nSyst];
  TH1D* jtptPbPbData_FULLBayes_h[nCentBins][nAbsEtaBins+1][nSyst][nBayesIter];
  TH1D* jtptPbPbData_FULLSvd_h[nCentBins][nAbsEtaBins+1][nSyst][nSvd];
  TH1D* jtetaPbPbData_FULLBayes_h[nCentBins][nSyst][nBayesIter];
  RooUnfoldBayes* rooUnfoldBayesPbPbData_FULLBayes_h[nCentBins][nAbsEtaBins+1][nSyst][nBayesIter];
  RooUnfoldSvd* rooUnfoldSvdPbPbData_FULLSvd_h[nCentBins][nAbsEtaBins+1][nSyst][nSvd];
  TH1D* refptPbPbData_h[nCentBins][nAbsEtaBins+1];
  TH2D* responsePbPbData_h[nCentBins][nAbsEtaBins+1];
  TH2D* responsePbPbData_Symm_h[nCentBins][nAbsEtaBins+1];
  RooUnfoldResponse* rooResponsePbPbData_h[nCentBins][nAbsEtaBins+1];
  TH1D* jtptPbPbData_RAW_Frac_h[nCentBins][nAbsEtaBins+1][nSyst];
  TH1D* refptPbPbData_Frac_h[nCentBins][nAbsEtaBins+1];


  TH1D* jtptPPData_RAW_h[nCentBins][nAbsEtaBins+1][nSyst];
  TH1D* jtptPPData_FULLBayes_h[nCentBins][nAbsEtaBins+1][nSyst][nBayesIter];
  TH1D* jtptPPData_FULLSvd_h[nCentBins][nAbsEtaBins+1][nSyst][nSvd];
  TH1D* jtetaPPData_FULLBayes_h[nCentBins][nSyst][nBayesIter];
  RooUnfoldBayes* rooUnfoldBayesPPData_FULLBayes_h[nCentBins][nAbsEtaBins+1][nSyst][nBayesIter];
  RooUnfoldSvd* rooUnfoldSvdPPData_FULLSvd_h[nCentBins][nAbsEtaBins+1][nSyst][nSvd];
  TH1D* refptPPData_h[nCentBins][nAbsEtaBins+1];
  TH2D* responsePPData_h[nCentBins][nAbsEtaBins+1];
  TH2D* responsePPData_Symm_h[nCentBins][nAbsEtaBins+1];
  RooUnfoldResponse* rooResponsePPData_h[nCentBins][nAbsEtaBins+1];
  TH1D* jtptPPData_RAW_Frac_h[nCentBins][nAbsEtaBins+1][nSyst];
  TH1D* refptPPData_Frac_h[nCentBins][nAbsEtaBins+1];

  TH1D* jtptPbPbMC_RAW_h[nCentBins][nAbsEtaBins+1][nSyst];
  TH1D* jtptPbPbMC_FULLBayes_h[nCentBins][nAbsEtaBins+1][nSyst][nBayesIter];
  RooUnfoldBayes* rooUnfoldBayesPbPbMC_FULLBayes_h[nCentBins][nAbsEtaBins+1][nSyst][nBayesIter];
  RooUnfoldSvd* rooUnfoldSvdPbPbMC_FULLSvd_h[nCentBins][nAbsEtaBins+1][nSyst][nSvd];
  TH1D* jtptPbPbMC_FULLSvd_h[nCentBins][nAbsEtaBins+1][nSyst][nSvd];
  TH1D* refptPbPbMC_h[nCentBins][nAbsEtaBins+1];
  TH1D* jtptPbPbMC_RAW_Frac_h[nCentBins][nAbsEtaBins+1][nSyst];
  TH1D* jtptPbPbMC_FULLBayes_Frac_h[nCentBins][nAbsEtaBins+1][nSyst][nBayesIter];
  RooUnfoldBayes* rooUnfoldBayesPbPbMC_FULLBayes_Frac_h[nCentBins][nAbsEtaBins+1][nSyst][nBayesIter];
  RooUnfoldSvd* rooUnfoldSvdPbPbMC_FULLSvd_Frac_h[nCentBins][nAbsEtaBins+1][nSyst][nSvd];
  TH1D* jtptPbPbMC_FULLSvd_Frac_h[nCentBins][nAbsEtaBins+1][nSyst][nSvd];
  TH1D* refptPbPbMC_Frac_h[nCentBins][nAbsEtaBins+1];
  TH2D* responsePbPbMC_h[nCentBins][nAbsEtaBins+1];
  TH2D* responsePbPbMC_Symm_h[nCentBins][nAbsEtaBins+1];
  RooUnfoldResponse* rooResponsePbPbMC_h[nCentBins][nAbsEtaBins+1];
  TH1D* jtetaPbPbMC_FULLBayes_h[nCentBins][nSyst][nBayesIter];
  TH1D* jtetaPbPbMC_FULLBayes_Frac_h[nCentBins][nSyst][nBayesIter];

  TH1D* jtptPPMC_RAW_h[nCentBins][nAbsEtaBins+1][nSyst]; 
  TH1D* jtptPPMC_FULLBayes_h[nCentBins][nAbsEtaBins+1][nSyst][nBayesIter];
  RooUnfoldBayes* rooUnfoldBayesPPMC_FULLBayes_h[nCentBins][nAbsEtaBins+1][nSyst][nBayesIter];
  RooUnfoldSvd* rooUnfoldSvdPPMC_FULLSvd_h[nCentBins][nAbsEtaBins+1][nSyst][nSvd];
  TH1D* jtptPPMC_FULLSvd_h[nCentBins][nAbsEtaBins+1][nSyst][nSvd];
  TH1D* refptPPMC_h[nCentBins][nAbsEtaBins+1];
  TH1D* jtptPPMC_RAW_Frac_h[nCentBins][nAbsEtaBins+1][nSyst]; 
  TH1D* jtptPPMC_FULLBayes_Frac_h[nCentBins][nAbsEtaBins+1][nSyst][nBayesIter];
  RooUnfoldBayes* rooUnfoldBayesPPMC_FULLBayes_Frac_h[nCentBins][nAbsEtaBins+1][nSyst][nBayesIter];
  RooUnfoldSvd* rooUnfoldSvdPPMC_FULLSvd_Frac_h[nCentBins][nAbsEtaBins+1][nSyst][nSvd];
  TH1D* jtptPPMC_FULLSvd_Frac_h[nCentBins][nAbsEtaBins+1][nSyst][nSvd];
  TH1D* refptPPMC_Frac_h[nCentBins][nAbsEtaBins+1];
  TH2D* responsePPMC_h[nCentBins][nAbsEtaBins+1];
  TH2D* responsePPMC_Symm_h[nCentBins][nAbsEtaBins+1];
  RooUnfoldResponse* rooResponsePPMC_h[nCentBins][nAbsEtaBins+1];
  TH1D* jtetaPPMC_FULLBayes_h[nCentBins][nSyst][nBayesIter];
  TH1D* jtetaPPMC_FULLBayes_Frac_h[nCentBins][nSyst][nBayesIter];

  const std::string algoStrPbPb = jetTreeListPbPb.at(0).substr(0, jetTreeListPbPb.at(0).find("/"));
  const std::string algoStrPP = jetTreeListPP.at(0).substr(0, jetTreeListPP.at(0).find("/"));
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    const std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
      
    const std::string dirStrPbPb = algoStrPbPb + "/" + centStr;
    const std::string dirStrPP = algoStrPP + "/" + centStr;
    
    outFile_p->cd();
    jetPbPbDirs_p = outFile_p->GetDirectory(algoStrPbPb.c_str());
    if(jetPbPbDirs_p) jetPbPbDirs_p->cd();
    else{
      jetPbPbDirs_p = outFile_p->mkdir(algoStrPbPb.c_str());
      jetPbPbDirs_p->cd();
    }
    
    outFile_p->cd();
    centPbPbDirs_p[cI] = outFile_p->GetDirectory(dirStrPbPb.c_str());
    if(centPbPbDirs_p[cI]) centPbPbDirs_p[cI]->cd();
    else{
      centPbPbDirs_p[cI] = outFile_p->mkdir(dirStrPbPb.c_str());
      centPbPbDirs_p[cI]->cd();
    }
    
    std::string treeNamePbPbData = algoStrPbPb + "_Data_" + centStr;
    std::string treeNamePPData = algoStrPP + "_Data_" + centStr;

    std::string treeNamePbPbMC = algoStrPbPb + "_MC_" + centStr;
    std::string treeNamePPMC = algoStrPP + "_MC_" + centStr;
    
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    treeCentSystCheckPbPbData_p[cI] = new TTree(treeNamePbPbData.c_str(), "");

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    for(Int_t sI = 0; sI < nSystForFill; ++sI){
      jDataOutPbPbData[cI][sI] = new jetData();
    }
    
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    
    treeCentSystCheckPbPbData_p[cI]->Branch("hiBin", &(jDataOutPbPbData[cI][0]->hiBin_), "hiBin/I");
    treeCentSystCheckPbPbData_p[cI]->Branch("nref", &(jDataOutPbPbData[cI][0]->nref_), "nref/I");
    treeCentSystCheckPbPbData_p[cI]->Branch("jteta", (jDataOutPbPbData[cI][0]->jteta_), "jteta[nref]/F");
    
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    
    for(Int_t sI = 0; sI < nSystForFill; ++sI){
      std::string jtBranchName1 = "jtpt_" + systStrForFill[sI];
      std::string jtBranchName2 = "jtpt_" + systStrForFill[sI] + "[nref]/F";
      
      std::string rawBranchName1 = "rawpt_" + systStrForFill[sI];
      std::string rawBranchName2 = "rawpt_" + systStrForFill[sI] + "[nref]/F";
      
      treeCentSystCheckPbPbData_p[cI]->Branch(jtBranchName1.c_str(), (jDataOutPbPbData[cI][sI]->jtpt_), jtBranchName2.c_str());
      treeCentSystCheckPbPbData_p[cI]->Branch(rawBranchName1.c_str(), (jDataOutPbPbData[cI][sI]->rawpt_), rawBranchName2.c_str());
    }
    

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    treeCentSystCheckPbPbMC_p[cI] = new TTree(treeNamePbPbMC.c_str(), "");
    for(Int_t sI = 0; sI < nSystForFill; ++sI){
      jDataOutPbPbMC[cI][sI] = new jetData();
    }
    
    treeCentSystCheckPbPbMC_p[cI]->Branch("hiBin", &(jDataOutPbPbMC[cI][0]->hiBin_), "hiBin/I");
    treeCentSystCheckPbPbMC_p[cI]->Branch("pthatWeight", &(jDataOutPbPbMC[cI][0]->pthatWeight_), "pthatWeight/F");
    treeCentSystCheckPbPbMC_p[cI]->Branch("fullWeight", &(jDataOutPbPbMC[cI][0]->fullWeight_), "fullWeight/F");
    treeCentSystCheckPbPbMC_p[cI]->Branch("nref", &(jDataOutPbPbMC[cI][0]->nref_), "nref/I");
    treeCentSystCheckPbPbMC_p[cI]->Branch("jteta", (jDataOutPbPbMC[cI][0]->jteta_), "jteta[nref]/F");
    treeCentSystCheckPbPbMC_p[cI]->Branch("refpt", (jDataOutPbPbMC[cI][0]->refpt_), "refpt[nref]/F");
    
    for(Int_t sI = 0; sI < nSystForFill; ++sI){
      std::string jtBranchName1 = "jtpt_" + systStrForFill[sI];
      std::string jtBranchName2 = "jtpt_" + systStrForFill[sI] + "[nref]/F";
      
      std::string rawBranchName1 = "rawpt_" + systStrForFill[sI];
      std::string rawBranchName2 = "rawpt_" + systStrForFill[sI] + "[nref]/F";
      
      treeCentSystCheckPbPbMC_p[cI]->Branch(jtBranchName1.c_str(), (jDataOutPbPbMC[cI][sI]->jtpt_), jtBranchName2.c_str());
      treeCentSystCheckPbPbMC_p[cI]->Branch(rawBranchName1.c_str(), (jDataOutPbPbMC[cI][sI]->rawpt_), rawBranchName2.c_str());
    }

    treeCentSystCheckPPData_p[cI] = new TTree(treeNamePPData.c_str(), "");
    for(Int_t sI = 0; sI < nSystForFill; ++sI){
      jDataOutPPData[cI][sI] = new jetData();
    }
    
    treeCentSystCheckPPData_p[cI]->Branch("nref", &(jDataOutPPData[cI][0]->nref_), "nref/I");
    treeCentSystCheckPPData_p[cI]->Branch("jteta", (jDataOutPPData[cI][0]->jteta_), "jteta[nref]/F");
    
    for(Int_t sI = 0; sI < nSystForFill; ++sI){
      std::string jtBranchName1 = "jtpt_" + systStrForFill[sI];
      std::string jtBranchName2 = "jtpt_" + systStrForFill[sI] + "[nref]/F";
      
      std::string rawBranchName1 = "rawpt_" + systStrForFill[sI];
      std::string rawBranchName2 = "rawpt_" + systStrForFill[sI] + "[nref]/F";
      
      treeCentSystCheckPPData_p[cI]->Branch(jtBranchName1.c_str(), (jDataOutPPData[cI][sI]->jtpt_), jtBranchName2.c_str());
      treeCentSystCheckPPData_p[cI]->Branch(rawBranchName1.c_str(), (jDataOutPPData[cI][sI]->rawpt_), rawBranchName2.c_str());
    }

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    treeCentSystCheckPPMC_p[cI] = new TTree(treeNamePPMC.c_str(), "");
    for(Int_t sI = 0; sI < nSystForFill; ++sI){
      jDataOutPPMC[cI][sI] = new jetData();
    }
    
    treeCentSystCheckPPMC_p[cI]->Branch("pthatWeight", &(jDataOutPPMC[cI][0]->pthatWeight_), "pthatWeight/F");
    treeCentSystCheckPPMC_p[cI]->Branch("fullWeight", &(jDataOutPPMC[cI][0]->fullWeight_), "fullWeight/F");
    treeCentSystCheckPPMC_p[cI]->Branch("nref", &(jDataOutPPMC[cI][0]->nref_), "nref/I");
    treeCentSystCheckPPMC_p[cI]->Branch("jteta", (jDataOutPPMC[cI][0]->jteta_), "jteta[nref]/F");
    treeCentSystCheckPPMC_p[cI]->Branch("refpt", (jDataOutPPMC[cI][0]->refpt_), "refpt[nref]/F");
    
    for(Int_t sI = 0; sI < nSystForFill; ++sI){
      std::string jtBranchName1 = "jtpt_" + systStrForFill[sI];
      std::string jtBranchName2 = "jtpt_" + systStrForFill[sI] + "[nref]/F";
      
      std::string rawBranchName1 = "rawpt_" + systStrForFill[sI];
      std::string rawBranchName2 = "rawpt_" + systStrForFill[sI] + "[nref]/F";
      
      treeCentSystCheckPPMC_p[cI]->Branch(jtBranchName1.c_str(), (jDataOutPPMC[cI][sI]->jtpt_), jtBranchName2.c_str());
      treeCentSystCheckPPMC_p[cI]->Branch(rawBranchName1.c_str(), (jDataOutPPMC[cI][sI]->rawpt_), rawBranchName2.c_str());
    }
    

    for(Int_t aI = 0; aI < nAbsEtaBins+1; ++aI){
      std::string absEtaStr = "AbsEtaTot";
      if(aI < nAbsEtaBins) absEtaStr = "AbsEta" + prettyString(absEtaBinsLow[aI], 1, true) + "to" + prettyString(absEtaBinsHi[aI], 1, true);
      const std::string binStrPbPb = algoStrPbPb + "_" + centStr + "_" + absEtaStr;
      
      for(Int_t sI = 0; sI < nSyst; ++sI){
	jtptPbPbData_RAW_h[cI][aI][sI] = new TH1D(("jtptPbPbData_" + binStrPbPb + "_RAW_" + systStr[sI] + "_h").c_str(), ";Reco. Jet p_{T};Counts", recoNBins[cI], recoBins[cI]);
	jtptPbPbMC_RAW_h[cI][aI][sI] = new TH1D(("jtptPbPbMC_" + binStrPbPb + "_RAW_" + systStr[sI] + "_h").c_str(), ";Reco. Jet p_{T};Counts", recoNBins[cI], recoBins[cI]);
	jtptPbPbData_RAW_Frac_h[cI][aI][sI] = NULL;
	jtptPbPbMC_RAW_Frac_h[cI][aI][sI] = new TH1D(("jtptPbPbMC_" + binStrPbPb + "_RAW_" + systStr[sI] + "_Frac_h").c_str(), ";Reco. Jet p_{T};Counts", recoNBins[cI], recoBins[cI]);
	
	for(Int_t bI = 0; bI < nBayesIter; ++bI){
	  jtptPbPbData_FULLBayes_h[cI][aI][sI][bI] = new TH1D(("jtptPbPbData_" + binStrPbPb + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[sI] + "_h").c_str(), ";Reco. Jet p_{T};Counts", truthNBins[cI], truthBins[cI]);
	  jtptPbPbMC_FULLBayes_h[cI][aI][sI][bI] = new TH1D(("jtptPbPbMC_" + binStrPbPb + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[sI] + "_h").c_str(), ";Reco. Jet p_{T};Counts", truthNBins[cI], truthBins[cI]);
	  jtptPbPbMC_FULLBayes_Frac_h[cI][aI][sI][bI] = new TH1D(("jtptPbPbMC_" + binStrPbPb + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[sI] + "_Frac_h").c_str(), ";Reco. Jet p_{T};Counts", truthNBins[cI], truthBins[cI]);
	  
	  rooUnfoldBayesPbPbData_FULLBayes_h[cI][aI][sI][bI] = new RooUnfoldBayes("rooUnfoldBayesPbPbData_" + binStrPbPb + "_FullBayes" + std::to_string(bI+1) + "_" + systStr[sI] + "_h", "");
	  rooUnfoldBayesPbPbMC_FULLBayes_h[cI][aI][sI][bI] = new RooUnfoldBayes("rooUnfoldBayesPbPbMC_" + binStrPbPb + "_FullBayes" + std::to_string(bI+1) + "_" + systStr[sI] + "_h", "");
	  rooUnfoldBayesPbPbMC_FULLBayes_Frac_h[cI][aI][sI][bI] = new RooUnfoldBayes("rooUnfoldBayesPbPbMC_" + binStrPbPb + "_FullBayes" + std::to_string(bI+1) + "_" + systStr[sI] + "_Frac_h", "");
	}
	
	
	for(Int_t bI = 0; bI < nSvd; ++bI){
	  if(bI+1 > truthNBins[cI]) break;
	  
	  jtptPbPbData_FULLSvd_h[cI][aI][sI][bI] = new TH1D(("jtptPbPbData_" + binStrPbPb + "_FULLSvd" + std::to_string(bI+1) + "_" + systStr[sI] + "_h").c_str(), ";Reco. Jet p_{T};Counts", truthNBins[cI], truthBins[cI]);
	  jtptPbPbMC_FULLSvd_h[cI][aI][sI][bI] = new TH1D(("jtptPbPbMC_" + binStrPbPb + "_FULLSvd" + std::to_string(bI+1) + "_" + systStr[sI] + "_h").c_str(), ";Reco. Jet p_{T};Counts", truthNBins[cI], truthBins[cI]);
	  jtptPbPbMC_FULLSvd_Frac_h[cI][aI][sI][bI] = new TH1D(("jtptPbPbMC_" + binStrPbPb + "_FULLSvd" + std::to_string(bI+1) + "_" + systStr[sI] + "_Frac_h").c_str(), ";Reco. Jet p_{T};Counts", truthNBins[cI], truthBins[cI]);
	  
	  rooUnfoldSvdPbPbData_FULLSvd_h[cI][aI][sI][bI] = new RooUnfoldSvd(("rooUnfoldSvdPbPbData_" + binStrPbPb + "_FullSvd" + std::to_string(bI+1) + "_" + systStr[sI] + "_h").c_str(), "");
	  rooUnfoldSvdPbPbMC_FULLSvd_h[cI][aI][sI][bI] = new RooUnfoldSvd(("rooUnfoldSvdPbPbMC_" + binStrPbPb + "_FullSvd" + std::to_string(bI+1) + "_" + systStr[sI] + "_h").c_str(), "");
	  rooUnfoldSvdPbPbMC_FULLSvd_Frac_h[cI][aI][sI][bI] = new RooUnfoldSvd(("rooUnfoldSvdPbPbMC_" + binStrPbPb + "_FullSvd" + std::to_string(bI+1) + "_" + systStr[sI] + "_Frac_h").c_str(), "");
	}
      }
      
      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
      
      refptPbPbMC_h[cI][aI] = new TH1D(("refptPbPbMC_" + binStrPbPb + "_h").c_str(), ";Truth Jet p_{T};Counts", truthNBins[cI], truthBins[cI]);
      refptPbPbMC_Frac_h[cI][aI] = new TH1D(("refptPbPbMC_" + binStrPbPb + "_Frac_h").c_str(), ";Truth Jet p_{T};Counts", truthNBins[cI], truthBins[cI]);
      
      refptPbPbData_h[cI][aI] = NULL;
      refptPbPbData_Frac_h[cI][aI] = NULL;
      responsePbPbData_h[cI][aI] = NULL;
      responsePbPbData_Symm_h[cI][aI] = NULL;
      rooResponsePbPbData_h[cI][aI] = NULL;
      
      responsePbPbMC_h[cI][aI] = new TH2D(("responsePbPbMC_" + binStrPbPb + "_h").c_str(), ";Reco. Jet p_{T};Truth Jet p_{T}", recoNBins[cI], recoBins[cI], truthNBins[cI], truthBins[cI]);
      responsePbPbMC_Symm_h[cI][aI] = new TH2D(("responsePbPbMC_" + binStrPbPb + "_Symm_h").c_str(), ";Reco. Jet p_{T};Truth Jet p_{T}", truthNBins[cI], truthBins[cI], truthNBins[cI], truthBins[cI]);
      rooResponsePbPbMC_h[cI][aI] = new RooUnfoldResponse(("rooResponsePbPbMC_" + binStrPbPb + "_h").c_str(), "");
      rooResponsePbPbMC_h[cI][aI]->Setup(jtptPbPbMC_RAW_h[cI][aI][0], refptPbPbMC_h[cI][aI]);	
    }
    
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    
    for(Int_t sI = 0; sI < nSyst; ++sI){
      for(Int_t bI = 0; bI < nBayesIter; ++bI){
	jtetaPbPbData_FULLBayes_h[cI][sI][bI] = new TH1D(("jtetaPbPbData_" + algoStrPbPb + "_" + centStr + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[sI] + "_h").c_str(), ";|#eta_{Jet}|;Counts", nAbsEtaBins-1, absEtaBins);
	
	jtetaPbPbMC_FULLBayes_h[cI][sI][bI] = new TH1D(("jtetaPbPbMC_" + algoStrPbPb + "_" + centStr + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[sI] + "_h").c_str(), ";|#eta_{Jet}|;Counts", nAbsEtaBins-1, absEtaBins);
	
	jtetaPbPbMC_FULLBayes_Frac_h[cI][sI][bI] = new TH1D(("jtetaPbPbMC_" + algoStrPbPb + "_" + centStr + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[sI] + "_Frac_h").c_str(), ";|#eta_{Jet}|;Counts", nAbsEtaBins-1, absEtaBins);
      }
    }
    
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
      
    outFile_p->cd();
    jetPPDirs_p = outFile_p->GetDirectory(algoStrPP.c_str());
    if(jetPPDirs_p) jetPPDirs_p->cd();
    else{
      jetPPDirs_p = outFile_p->mkdir(algoStrPP.c_str());
      jetPPDirs_p->cd();
    }
    
    outFile_p->cd();
    centPPDirs_p[cI] = outFile_p->GetDirectory(dirStrPP.c_str());
    if(centPPDirs_p[cI]) centPPDirs_p[cI]->cd();
    else{
      centPPDirs_p[cI] = outFile_p->mkdir(dirStrPP.c_str());
      centPPDirs_p[cI]->cd();
    }
    
    
    for(Int_t aI = 0; aI < nAbsEtaBins+1; ++aI){	
      std::string absEtaStr = "AbsEtaTot";
      if(aI < nAbsEtaBins) absEtaStr = "AbsEta" + prettyString(absEtaBinsLow[aI], 1, true) + "to" + prettyString(absEtaBinsHi[aI], 1, true);
      const std::string binStrPP = algoStrPP + "_" + centStr + "_" + absEtaStr;
      
      for(Int_t sI = 0; sI < nSyst; ++sI){
	jtptPPData_RAW_h[cI][aI][sI] = new TH1D(("jtptPPData_" + binStrPP + "_RAW_" + systStr[sI] + "_h").c_str(), ";Reco. Jet p_{T};Counts", recoNBins[cI], recoBins[cI]);
	jtptPPMC_RAW_h[cI][aI][sI] = new TH1D(("jtptPPMC_" + binStrPP + "_RAW_" + systStr[sI] + "_h").c_str(), ";Reco. Jet p_{T};Counts", recoNBins[cI], recoBins[cI]);
	jtptPPMC_RAW_Frac_h[cI][aI][sI] = new TH1D(("jtptPPMC_" + binStrPP + "_RAW_" + systStr[sI] + "_Frac_h").c_str(), ";Reco. Jet p_{T};Counts", recoNBins[cI], recoBins[cI]);
	jtptPPData_RAW_Frac_h[cI][aI][sI] = NULL;
	
	for(Int_t bI = 0; bI < nBayesIter; ++bI){
	  jtptPPData_FULLBayes_h[cI][aI][sI][bI] = new TH1D(("jtptPPData_" + binStrPP + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[sI] + "_h").c_str(), ";Reco. Jet p_{T};Counts", truthNBins[cI], truthBins[cI]);
	  jtptPPMC_FULLBayes_h[cI][aI][sI][bI] = new TH1D(("jtptPPMC_" + binStrPP + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[sI] + "_h").c_str(), ";Reco. Jet p_{T};Counts", truthNBins[cI], truthBins[cI]);
	  jtptPPMC_FULLBayes_Frac_h[cI][aI][sI][bI] = new TH1D(("jtptPPMC_" + binStrPP + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[sI] + "_Frac_h").c_str(), ";Reco. Jet p_{T};Counts", truthNBins[cI], truthBins[cI]);
	  
	  rooUnfoldBayesPPData_FULLBayes_h[cI][aI][sI][bI] = new RooUnfoldBayes("rooUnfoldBayesPPData_" + binStrPP + "_FullBayes" + std::to_string(bI+1) + "_" + systStr[sI] + "_h", "");
	  rooUnfoldBayesPPMC_FULLBayes_h[cI][aI][sI][bI] = new RooUnfoldBayes("rooUnfoldBayesPPMC_" + binStrPP + "_FullBayes" + std::to_string(bI+1) + "_" + systStr[sI] + "_h", "");
	  rooUnfoldBayesPPMC_FULLBayes_Frac_h[cI][aI][sI][bI] = new RooUnfoldBayes("rooUnfoldBayesPPMC_" + binStrPP + "_FullBayes" + std::to_string(bI+1) + "_" + systStr[sI] + "_Frac_h", "");
	}
	
	for(Int_t bI = 0; bI < nSvd; ++bI){
	  if(bI+1 > truthNBins[cI]) break;
	  
	  jtptPPData_FULLSvd_h[cI][aI][sI][bI] = new TH1D(("jtptPPData_" + binStrPP + "_FULLSvd" + std::to_string(bI+1) + "_" + systStr[sI] + "_h").c_str(), ";Reco. Jet p_{T};Counts", truthNBins[cI], truthBins[cI]);
	  jtptPPMC_FULLSvd_h[cI][aI][sI][bI] = new TH1D(("jtptPPMC_" + binStrPP + "_FULLSvd" + std::to_string(bI+1) + "_" + systStr[sI] + "_h").c_str(), ";Reco. Jet p_{T};Counts", truthNBins[cI], truthBins[cI]);
	  jtptPPMC_FULLSvd_Frac_h[cI][aI][sI][bI] = new TH1D(("jtptPPMC_" + binStrPP + "_FULLSvd" + std::to_string(bI+1) + "_" + systStr[sI] + "_Frac_h").c_str(), ";Reco. Jet p_{T};Counts", truthNBins[cI], truthBins[cI]);
	  
	  rooUnfoldSvdPPData_FULLSvd_h[cI][aI][sI][bI] = new RooUnfoldSvd(("rooUnfoldSvdPPData_" + binStrPP + "_FullSvd" + std::to_string(bI+1) + "_" + systStr[sI] + "_h").c_str(), "");
	  rooUnfoldSvdPPMC_FULLSvd_h[cI][aI][sI][bI] = new RooUnfoldSvd(("rooUnfoldSvdPPMC_" + binStrPP + "_FullSvd" + std::to_string(bI+1) + "_" + systStr[sI] + "_h").c_str(), "");
	  rooUnfoldSvdPPMC_FULLSvd_Frac_h[cI][aI][sI][bI] = new RooUnfoldSvd(("rooUnfoldSvdPPMC_" + binStrPP + "_FullSvd" + std::to_string(bI+1) + "_" + systStr[sI] + "_Frac_h").c_str(), "");
	}
      }
      
      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
      
      refptPPMC_h[cI][aI] = new TH1D(("refptPPMC_" + binStrPP + "_h").c_str(), ";Truth Jet p_{T};Counts", truthNBins[cI], truthBins[cI]);
      refptPPMC_Frac_h[cI][aI] = new TH1D(("refptPPMC_" + binStrPP + "_Frac_h").c_str(), ";Truth Jet p_{T};Counts", truthNBins[cI], truthBins[cI]);
      responsePPMC_h[cI][aI] = new TH2D(("responsePPMC_" + binStrPP + "_h").c_str(), ";Reco. Jet p_{T};Truth Jet p_{T}", recoNBins[cI], recoBins[cI], truthNBins[cI], truthBins[cI]);
      responsePPMC_Symm_h[cI][aI] = new TH2D(("responsePPMC_" + binStrPP + "_Symm_h").c_str(), ";Reco. Jet p_{T};Truth Jet p_{T}", truthNBins[cI], truthBins[cI], truthNBins[cI], truthBins[cI]);
      rooResponsePPMC_h[cI][aI] = new RooUnfoldResponse(("rooResponsePPMC_" + binStrPP + "_h").c_str(), "");
      rooResponsePPMC_h[cI][aI]->Setup(jtptPPMC_RAW_h[cI][aI][0], refptPPMC_h[cI][aI]);
      
      refptPPData_h[cI][aI] = NULL;
      refptPPData_Frac_h[cI][aI] = NULL;
      rooResponsePPData_h[cI][aI] = NULL;
      responsePPData_h[cI][aI] = NULL;
      responsePPData_Symm_h[cI][aI] = NULL;
      
      centerTitles({refptPbPbMC_h[cI][aI], refptPbPbMC_Frac_h[cI][aI], responsePbPbMC_h[cI][aI], responsePbPbMC_Symm_h[cI][aI], refptPPMC_h[cI][aI], refptPPMC_Frac_h[cI][aI], responsePPMC_h[cI][aI], responsePPMC_Symm_h[cI][aI]});
      setSumW2({refptPbPbMC_h[cI][aI], refptPbPbMC_Frac_h[cI][aI], responsePbPbMC_h[cI][aI], responsePbPbMC_Symm_h[cI][aI], refptPPMC_h[cI][aI], refptPPMC_Frac_h[cI][aI], responsePPMC_h[cI][aI], responsePPMC_Symm_h[cI][aI]});
      
      for(Int_t sI = 0; sI < nSyst; ++sI){
	centerTitles({jtptPbPbData_RAW_h[cI][aI][sI], jtptPbPbMC_RAW_h[cI][aI][sI], jtptPbPbMC_RAW_Frac_h[cI][aI][sI], jtptPPData_RAW_h[cI][aI][sI], jtptPPMC_RAW_h[cI][aI][sI], jtptPPMC_RAW_Frac_h[cI][aI][sI]});
	setSumW2({jtptPbPbData_RAW_h[cI][aI][sI], jtptPbPbMC_RAW_h[cI][aI][sI], jtptPbPbMC_RAW_Frac_h[cI][aI][sI], jtptPPData_RAW_h[cI][aI][sI], jtptPPMC_RAW_h[cI][aI][sI], jtptPPMC_RAW_Frac_h[cI][aI][sI]});
	
	for(Int_t bI = 0; bI < nBayesIter; ++bI){
	  centerTitles({jtptPbPbData_FULLBayes_h[cI][aI][sI][bI], jtptPbPbMC_FULLBayes_h[cI][aI][sI][bI], jtptPbPbMC_FULLBayes_Frac_h[cI][aI][sI][bI], jtptPPData_FULLBayes_h[cI][aI][sI][bI], jtptPPMC_FULLBayes_h[cI][aI][sI][bI], jtptPPMC_FULLBayes_Frac_h[cI][aI][sI][bI]});
	  setSumW2({jtptPbPbData_FULLBayes_h[cI][aI][sI][bI], jtptPbPbMC_FULLBayes_h[cI][aI][sI][bI], jtptPbPbMC_FULLBayes_Frac_h[cI][aI][sI][bI], jtptPPData_FULLBayes_h[cI][aI][sI][bI], jtptPPMC_FULLBayes_h[cI][aI][sI][bI], jtptPPMC_FULLBayes_Frac_h[cI][aI][sI][bI]});
	}
	
	for(Int_t bI = 0; bI < nSvd; ++bI){
	  if(bI+1 > truthNBins[cI]) break;
	  
	  centerTitles({jtptPbPbData_FULLSvd_h[cI][aI][sI][bI], jtptPbPbMC_FULLSvd_h[cI][aI][sI][bI], jtptPbPbMC_FULLSvd_Frac_h[cI][aI][sI][bI], jtptPPData_FULLSvd_h[cI][aI][sI][bI], jtptPPMC_FULLSvd_h[cI][aI][sI][bI], jtptPPMC_FULLSvd_Frac_h[cI][aI][sI][bI]});
	  setSumW2({jtptPbPbData_FULLSvd_h[cI][aI][sI][bI], jtptPbPbMC_FULLSvd_h[cI][aI][sI][bI], jtptPbPbMC_FULLSvd_Frac_h[cI][aI][sI][bI], jtptPPData_FULLSvd_h[cI][aI][sI][bI], jtptPPMC_FULLSvd_h[cI][aI][sI][bI], jtptPPMC_FULLSvd_Frac_h[cI][aI][sI][bI]});
	}
      }
    }
    
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    
    for(Int_t sI = 0; sI < nSyst; ++sI){
      for(Int_t bI = 0; bI < nBayesIter; ++bI){
	jtetaPPData_FULLBayes_h[cI][sI][bI] = new TH1D(("jtetaPPData_" + algoStrPP + "_" + centStr + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[sI] + "_h").c_str(), ";|#eta_{Jet}|;Counts", nAbsEtaBins-1, absEtaBins);

	jtetaPPMC_FULLBayes_h[cI][sI][bI] = new TH1D(("jtetaPPMC_" + algoStrPP + "_" + centStr + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[sI] + "_h").c_str(), ";|#eta_{Jet}|;Counts", nAbsEtaBins-1, absEtaBins);
	
	jtetaPPMC_FULLBayes_Frac_h[cI][sI][bI] = new TH1D(("jtetaPPMC_" + algoStrPP + "_" + centStr + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[sI] + "_Frac_h").c_str(), ";|#eta_{Jet}|;Counts", nAbsEtaBins-1, absEtaBins);
      }
      
      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    }
  }

  histWatch.stop();
  dumpTiming({startWatch, histWatch}, {"start", "hist"});
  pbpbDataWatch.start();

  TTree* pbpbTrees_p = (TTree*)inFilePbPb_p->Get(jetTreeListPbPb.at(0).c_str());
  TTree* pbpbMCTrees_p = (TTree*)inFilePbPbMC_p->Get(jetTreeListPbPbMC.at(0).c_str());
  TTree* ppTrees_p = (TTree*)inFilePP_p->Get(jetTreeListPP.at(0).c_str());
  TTree* ppMCTrees_p = (TTree*)inFilePPMC_p->Get(jetTreeListPPMC.at(0).c_str());

  std::cout << "Processing PbPbData" << std::endl;
  inFilePbPb_p->cd();
  std::string algoStr = jetTreeListPbPb.at(0).substr(0, jetTreeListPbPb.at(0).find("/"));
  doTreeFill(pbpbTrees_p, treeCentSystCheckPbPbData_p, (jDataOutPbPbData), false, false, algoStr, 0, jtptPbPbData_RAW_h, refptPbPbData_h, responsePbPbData_h, responsePbPbData_Symm_h, rooResponsePbPbData_h, jtptPbPbData_RAW_Frac_h, refptPbPbData_Frac_h);
  inFilePbPb_p->Close();
  delete inFilePbPb_p;

  pbpbDataWatch.stop();
  dumpTiming({startWatch, histWatch, pbpbDataWatch}, {"start", "hist", "pbpbData"});
  ppDataWatch.start();

  std::cout << "Processing PP" << std::endl;
  inFilePP_p->cd();
  doTreeFill(ppTrees_p, treeCentSystCheckPPData_p, (jDataOutPPData), false, true, algoStr, 0, jtptPPData_RAW_h, refptPPData_h, responsePPData_h, responsePPData_Symm_h,  rooResponsePPData_h, jtptPPData_RAW_Frac_h, refptPPData_Frac_h);
  inFilePP_p->Close();
  delete inFilePP_p;

  ppDataWatch.stop();
  dumpTiming({startWatch, histWatch, pbpbDataWatch, ppDataWatch}, {"start", "hist", "pbpbData", "ppData"});
  ppMCWatch.start();

  std::cout << "Processing PPMC" << std::endl;
  inFilePPMC_p->cd();
  doTreeFill(ppMCTrees_p, treeCentSystCheckPPMC_p, (jDataOutPPMC), true, true, algoStr, 0, jtptPPMC_RAW_h, refptPPMC_h, responsePPMC_h, responsePPMC_Symm_h,  rooResponsePPMC_h, jtptPPMC_RAW_Frac_h, refptPPMC_Frac_h);
  inFilePPMC_p->Close();
  delete inFilePPMC_p;

  ppMCWatch.stop();
  dumpTiming({startWatch, histWatch, pbpbDataWatch, ppDataWatch, ppMCWatch}, {"start", "hist", "pbpbData", "ppData", "ppMC"});
  pbpbMCWatch.start();

  std::cout << "Processing PbPbMC" << std::endl;
  inFilePbPbMC_p->cd();
  doTreeFill(pbpbMCTrees_p, treeCentSystCheckPbPbMC_p, (jDataOutPbPbMC), true, false, algoStr, 0, jtptPbPbMC_RAW_h, refptPbPbMC_h, responsePbPbMC_h, responsePbPbMC_Symm_h,  rooResponsePbPbMC_h, jtptPbPbMC_RAW_Frac_h, refptPbPbMC_Frac_h);
  inFilePbPbMC_p->Close();
  delete inFilePbPbMC_p;

  pbpbMCWatch.stop();
  dumpTiming({startWatch, histWatch, pbpbDataWatch, ppDataWatch, ppMCWatch, pbpbMCWatch}, {"start", "hist", "pbpbData", "ppData", "ppMC", "pbpbMC"});
  endWatch.start();
  countWatch.start();

  countWatch.stop();
  dumpTiming({startWatch, histWatch, pbpbDataWatch, ppDataWatch, ppMCWatch, pbpbMCWatch, countWatch}, {"start", "hist", "pbpbData", "ppData", "ppMC", "pbpbMC", "count"});
  unfoldLoopWatch.start();


  std::cout << "START UNFOLDING" << std::endl;
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    for(Int_t aI = 0; aI < nAbsEtaBins; ++aI){
      std::cout << "ON cent, absEta: " << cI << "/" << nCentBins << ", " << aI << "/" << nAbsEtaBins << std::endl;

      for(Int_t sI = 0; sI < nSystForFill; ++sI){
	
	if(doGlobalDebug) std::cout << "FIRST THREADS (" << cI << "/" << nCentBins << ", " << aI << "/" << nAbsEtaBins << ", " << sI << "/" << nSystForFill << ")" << std::endl;

	for(Int_t bI = 0; bI < nBayesIter; ++bI){
	  unfoldBayes(rooResponsePbPbMC_h[cI][aI], jtptPbPbData_RAW_h[cI][aI][systPosForFill[sI]], bI+1, jtptPbPbData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI], rooUnfoldBayesPbPbData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI]);
	}

	if(doSvd){
	  for(Int_t bI = 1; bI < nSvd; ++bI){
	    unfoldSvd(rooResponsePbPbMC_h[cI][aI], jtptPbPbData_RAW_h[cI][aI][systPosForFill[sI]], bI+1, jtptPbPbData_FULLSvd_h[cI][aI][systPosForFill[sI]][bI], rooUnfoldSvdPbPbData_FULLSvd_h[cI][aI][systPosForFill[sI]][bI]);
	  }
	}

	if(doGlobalDebug) std::cout << "SECOND THREADS (" << cI << "/" << nCentBins << ", " << aI << "/" << nAbsEtaBins << ", " << sI << "/" << nSystForFill << ")" << std::endl;

	for(Int_t bI = 0; bI < nBayesIter; ++bI){
	  unfoldBayes(rooResponsePbPbMC_h[cI][aI], jtptPbPbMC_RAW_h[cI][aI][systPosForFill[sI]], bI+1, jtptPbPbMC_FULLBayes_h[cI][aI][systPosForFill[sI]][bI], rooUnfoldBayesPbPbMC_FULLBayes_h[cI][aI][systPosForFill[sI]][bI]);
	}

	if(doSvd){
	  for(Int_t bI = 1; bI < nSvd; ++bI){
	    unfoldSvd(rooResponsePbPbMC_h[cI][aI], jtptPbPbMC_RAW_h[cI][aI][systPosForFill[sI]], bI+1, jtptPbPbMC_FULLSvd_h[cI][aI][systPosForFill[sI]][bI], rooUnfoldSvdPbPbMC_FULLSvd_h[cI][aI][systPosForFill[sI]][bI]);
	  }
	}

	if(doGlobalDebug) std::cout << "THIRD THREADS (" << cI << "/" << nCentBins << ", " << aI << "/" << nAbsEtaBins << ", " << sI << "/" << nSystForFill << ")" << std::endl;

	for(Int_t bI = 0; bI < nBayesIter; ++bI){
	  unfoldBayes(rooResponsePbPbMC_h[cI][aI], jtptPbPbMC_RAW_Frac_h[cI][aI][systPosForFill[sI]], bI+1, jtptPbPbMC_FULLBayes_Frac_h[cI][aI][systPosForFill[sI]][bI], rooUnfoldBayesPbPbMC_FULLBayes_Frac_h[cI][aI][systPosForFill[sI]][bI]);
	}

	if(doSvd){
	  for(Int_t bI = 1; bI < nSvd; ++bI){
	    unfoldSvd(rooResponsePbPbMC_h[cI][aI], jtptPbPbMC_RAW_Frac_h[cI][aI][systPosForFill[sI]], bI+1, jtptPbPbMC_FULLSvd_Frac_h[cI][aI][systPosForFill[sI]][bI], rooUnfoldSvdPbPbMC_FULLSvd_Frac_h[cI][aI][systPosForFill[sI]][bI]);
	  }
	}

	if(doGlobalDebug) std::cout << "FOURTH THREADS (" << cI << "/" << nCentBins << ", " << aI << "/" << nAbsEtaBins << ", " << sI << "/" << nSystForFill << ")" << std::endl;

	for(Int_t bI = 0; bI < nBayesIter; ++bI){
	  unfoldBayes(rooResponsePPMC_h[cI][aI], jtptPPData_RAW_h[cI][aI][systPosForFill[sI]], bI+1, jtptPPData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI], rooUnfoldBayesPPData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI]);
	}

	if(doSvd){
	  for(Int_t bI = 1; bI < nSvd; ++bI){
	    unfoldSvd(rooResponsePPMC_h[cI][aI], jtptPPData_RAW_h[cI][aI][systPosForFill[sI]], bI+1, jtptPPData_FULLSvd_h[cI][aI][systPosForFill[sI]][bI], rooUnfoldSvdPPData_FULLSvd_h[cI][aI][systPosForFill[sI]][bI]);
	  }
	}


	if(doGlobalDebug) std::cout << "FIFTH THREADS (" << cI << "/" << nCentBins << ", " << aI << "/" << nAbsEtaBins << ", " << sI << "/" << nSystForFill << ")" << std::endl;

	for(Int_t bI = 0; bI < nBayesIter; ++bI){
	  unfoldBayes(rooResponsePPMC_h[cI][aI], jtptPPMC_RAW_h[cI][aI][systPosForFill[sI]], bI+1, jtptPPMC_FULLBayes_h[cI][aI][systPosForFill[sI]][bI], rooUnfoldBayesPPMC_FULLBayes_h[cI][aI][systPosForFill[sI]][bI]);
	}

	if(doSvd){
	  for(Int_t bI = 1; bI < nSvd; ++bI){
	    unfoldSvd(rooResponsePPMC_h[cI][aI], jtptPPMC_RAW_h[cI][aI][systPosForFill[sI]], bI+1, jtptPPMC_FULLSvd_h[cI][aI][systPosForFill[sI]][bI], rooUnfoldSvdPPMC_FULLSvd_h[cI][aI][systPosForFill[sI]][bI]);
	  }
	}

	if(doGlobalDebug) std::cout << "SIXTH THREADS (" << cI << "/" << nCentBins << ", " << aI << "/" << nAbsEtaBins << ", " << sI << "/" << nSystForFill << ")" << std::endl;

	for(Int_t bI = 0; bI < nBayesIter; ++bI){
	  unfoldBayes(rooResponsePPMC_h[cI][aI], jtptPPMC_RAW_Frac_h[cI][aI][systPosForFill[sI]], bI+1, jtptPPMC_FULLBayes_Frac_h[cI][aI][systPosForFill[sI]][bI], rooUnfoldBayesPPMC_FULLBayes_Frac_h[cI][aI][systPosForFill[sI]][bI]);
	}    

	if(doSvd){
	  for(Int_t bI = 1; bI < nSvd; ++bI){
	    unfoldSvd(rooResponsePPMC_h[cI][aI], jtptPPMC_RAW_Frac_h[cI][aI][systPosForFill[sI]], bI+1, jtptPPMC_FULLSvd_Frac_h[cI][aI][systPosForFill[sI]][bI], rooUnfoldSvdPPMC_FULLSvd_Frac_h[cI][aI][systPosForFill[sI]][bI]);
	  }
	}

      }
    }
  }

  unfoldLoopWatch.stop();
  dumpTiming({startWatch, histWatch, pbpbDataWatch, ppDataWatch, ppMCWatch, pbpbMCWatch, countWatch, unfoldLoopWatch}, {"start", "hist", "pbpbData", "ppData", "ppMC", "pbpbMC", "count", "unfold"});
  plotValidWatch.start();
  
  std::cout << "THREADS COMPLETE" << std::endl;

  plotValidWatch.stop();
  dumpTiming({startWatch, histWatch, pbpbDataWatch, ppDataWatch, ppMCWatch, pbpbMCWatch, countWatch, unfoldLoopWatch, plotValidWatch}, {"start", "hist", "pbpbData", "ppData", "ppMC", "pbpbMC", "count", "unfold", "plotValid"});
  copyWatch.start();

  
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    const std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
    const std::string dirStrPbPb = algoStrPbPb + "/" + centStr;
    const std::string dirStrPP = algoStrPP + "/" + centStr;
    
    for(Int_t aI = 0; aI < nAbsEtaBins; ++aI){
      std::string absEtaStr = "AbsEta" + prettyString(absEtaBinsLow[aI],1,true) + "to" + prettyString(absEtaBinsHi[aI],1,true);
      if(absEtaStr.find("AbsEta0p0to2p0") != std::string::npos) continue;

      for(Int_t sI = 0; sI < nSystForCopy; ++sI){
	doHistCopy(jtptPbPbData_RAW_h[cI][aI][0], jtptPbPbData_RAW_h[cI][aI][systPosForCopy[sI]]);
	doHistCopy(jtptPPData_RAW_h[cI][aI][0], jtptPPData_RAW_h[cI][aI][systPosForCopy[sI]]);
	doHistCopy(jtptPbPbMC_RAW_h[cI][aI][0], jtptPbPbMC_RAW_h[cI][aI][systPosForCopy[sI]]);
	doHistCopy(jtptPPMC_RAW_h[cI][aI][0], jtptPPMC_RAW_h[cI][aI][systPosForCopy[sI]]);
	doHistCopy(jtptPbPbMC_RAW_Frac_h[cI][aI][0], jtptPbPbMC_RAW_Frac_h[cI][aI][systPosForCopy[sI]]);
	doHistCopy(jtptPPMC_RAW_Frac_h[cI][aI][0], jtptPPMC_RAW_Frac_h[cI][aI][systPosForCopy[sI]]);
	
	for(Int_t bI = 0; bI < nBayesIter; ++bI){
	  doHistCopy(jtptPbPbData_FULLBayes_h[cI][aI][0][bI], jtptPbPbData_FULLBayes_h[cI][aI][systPosForCopy[sI]][bI]);
	  doHistCopy(jtptPPData_FULLBayes_h[cI][aI][0][bI], jtptPPData_FULLBayes_h[cI][aI][systPosForCopy[sI]][bI]);
	  doHistCopy(jtptPbPbMC_FULLBayes_h[cI][aI][0][bI], jtptPbPbMC_FULLBayes_h[cI][aI][systPosForCopy[sI]][bI]);
	  doHistCopy(jtptPPMC_FULLBayes_h[cI][aI][0][bI], jtptPPMC_FULLBayes_h[cI][aI][systPosForCopy[sI]][bI]);
	  doHistCopy(jtptPbPbMC_FULLBayes_Frac_h[cI][aI][0][bI], jtptPbPbMC_FULLBayes_Frac_h[cI][aI][systPosForCopy[sI]][bI]);
	  doHistCopy(jtptPPMC_FULLBayes_Frac_h[cI][aI][0][bI], jtptPPMC_FULLBayes_Frac_h[cI][aI][systPosForCopy[sI]][bI]);
	}
	if(doSvd){
	  for(Int_t bI = 1; bI < nSvd; ++bI){
	    doHistCopy(jtptPbPbData_FULLSvd_h[cI][aI][0][bI], jtptPbPbData_FULLSvd_h[cI][aI][systPosForCopy[sI]][bI]);
	    doHistCopy(jtptPPData_FULLSvd_h[cI][aI][0][bI], jtptPPData_FULLSvd_h[cI][aI][systPosForCopy[sI]][bI]);
	    doHistCopy(jtptPbPbMC_FULLSvd_h[cI][aI][0][bI], jtptPbPbMC_FULLSvd_h[cI][aI][systPosForCopy[sI]][bI]);
	    doHistCopy(jtptPPMC_FULLSvd_h[cI][aI][0][bI], jtptPPMC_FULLSvd_h[cI][aI][systPosForCopy[sI]][bI]);
	    doHistCopy(jtptPbPbMC_FULLSvd_Frac_h[cI][aI][0][bI], jtptPbPbMC_FULLSvd_Frac_h[cI][aI][systPosForCopy[sI]][bI]);
	    doHistCopy(jtptPPMC_FULLSvd_Frac_h[cI][aI][0][bI], jtptPPMC_FULLSvd_Frac_h[cI][aI][systPosForCopy[sI]][bI]);
	  }
	}
      }	
    }
  }

  copyWatch.stop();
  dumpTiming({startWatch, histWatch, pbpbDataWatch, ppDataWatch, ppMCWatch, pbpbMCWatch, countWatch, unfoldLoopWatch, plotValidWatch, copyWatch}, {"start", "hist", "pbpbData", "ppData", "ppMC", "pbpbMC", "count", "unfold", "plotValid", "copy"});
  finalEtaWatch.start();
  std::cout << "WRITE COMPLETE" << std::endl;

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    const std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
    const std::string dirStrPbPb = algoStrPbPb + "/" + centStr;
    const std::string dirStrPP = algoStrPP + "/" + centStr;
    
    for(Int_t bIX = 0; bIX < refptPbPbMC_h[cI][nAbsEtaBins]->GetNbinsX(); ++bIX){
      refptPbPbMC_h[cI][nAbsEtaBins]->SetBinContent(bIX+1, 0.0);
      refptPbPbMC_Frac_h[cI][nAbsEtaBins]->SetBinContent(bIX+1, 0.0);
      
      refptPbPbMC_h[cI][nAbsEtaBins]->SetBinError(bIX+1, 0.0);
      refptPbPbMC_Frac_h[cI][nAbsEtaBins]->SetBinError(bIX+1, 0.0);
    }
    
    for(Int_t bIX = 0; bIX < refptPPMC_h[cI][nAbsEtaBins]->GetNbinsX(); ++bIX){
      refptPPMC_h[cI][nAbsEtaBins]->SetBinContent(bIX+1, 0.0);
      refptPPMC_Frac_h[cI][nAbsEtaBins]->SetBinContent(bIX+1, 0.0);
      
      refptPPMC_h[cI][nAbsEtaBins]->SetBinError(bIX+1, 0.0);
      refptPPMC_Frac_h[cI][nAbsEtaBins]->SetBinError(bIX+1, 0.0);
    }
    
    for(Int_t sI = 0; sI < nSyst; ++sI){
      for(Int_t bIX = 0; bIX < jtptPbPbData_RAW_h[cI][nAbsEtaBins][sI]->GetNbinsX(); ++bIX){
	jtptPbPbData_RAW_h[cI][nAbsEtaBins][sI]->SetBinContent(bIX+1, 0.0);
	jtptPbPbMC_RAW_h[cI][nAbsEtaBins][sI]->SetBinContent(bIX+1, 0.0);
	jtptPbPbMC_RAW_Frac_h[cI][nAbsEtaBins][sI]->SetBinContent(bIX+1, 0.0);
	
	jtptPbPbData_RAW_h[cI][nAbsEtaBins][sI]->SetBinError(bIX+1, 0.0);
	jtptPbPbMC_RAW_h[cI][nAbsEtaBins][sI]->SetBinError(bIX+1, 0.0);
	jtptPbPbMC_RAW_Frac_h[cI][nAbsEtaBins][sI]->SetBinError(bIX+1, 0.0);
      }
      

      for(Int_t bI = 0; bI < nBayesIter; ++bI){
	for(Int_t bIX = 0; bIX < jtptPbPbData_FULLBayes_h[cI][nAbsEtaBins][sI][bI]->GetNbinsX(); ++bIX){
	  jtptPbPbData_FULLBayes_h[cI][nAbsEtaBins][sI][bI]->SetBinContent(bIX+1, 0.0);
	  jtptPbPbMC_FULLBayes_h[cI][nAbsEtaBins][sI][bI]->SetBinContent(bIX+1, 0.0);
	  jtptPbPbMC_FULLBayes_Frac_h[cI][nAbsEtaBins][sI][bI]->SetBinContent(bIX+1, 0.0);
	  
	  jtptPbPbData_FULLBayes_h[cI][nAbsEtaBins][sI][bI]->SetBinError(bIX+1, 0.0);
	  jtptPbPbMC_FULLBayes_h[cI][nAbsEtaBins][sI][bI]->SetBinError(bIX+1, 0.0);
	  jtptPbPbMC_FULLBayes_Frac_h[cI][nAbsEtaBins][sI][bI]->SetBinError(bIX+1, 0.0);
	}	  
      }
      
      if(doSvd){
	for(Int_t bI = 1; bI < nSvd; ++bI){
	  for(Int_t bIX = 0; bIX < jtptPbPbData_FULLSvd_h[cI][nAbsEtaBins][sI][bI]->GetNbinsX(); ++bIX){
	    jtptPbPbData_FULLSvd_h[cI][nAbsEtaBins][sI][bI]->SetBinContent(bIX+1, 0.0);
	    jtptPbPbMC_FULLSvd_h[cI][nAbsEtaBins][sI][bI]->SetBinContent(bIX+1, 0.0);
	    jtptPbPbMC_FULLSvd_Frac_h[cI][nAbsEtaBins][sI][bI]->SetBinContent(bIX+1, 0.0);
	    
	    jtptPbPbData_FULLSvd_h[cI][nAbsEtaBins][sI][bI]->SetBinError(bIX+1, 0.0);
	    jtptPbPbMC_FULLSvd_h[cI][nAbsEtaBins][sI][bI]->SetBinError(bIX+1, 0.0);
	    jtptPbPbMC_FULLSvd_Frac_h[cI][nAbsEtaBins][sI][bI]->SetBinError(bIX+1, 0.0);
	  }	  
	}
      }
      
      for(Int_t bIX = 0; bIX < jtptPPData_RAW_h[cI][nAbsEtaBins][sI]->GetNbinsX(); ++bIX){
	jtptPPData_RAW_h[cI][nAbsEtaBins][sI]->SetBinContent(bIX+1, 0.0);
	jtptPPMC_RAW_h[cI][nAbsEtaBins][sI]->SetBinContent(bIX+1, 0.0);
	jtptPPMC_RAW_Frac_h[cI][nAbsEtaBins][sI]->SetBinContent(bIX+1, 0.0);
	
	jtptPPData_RAW_h[cI][nAbsEtaBins][sI]->SetBinError(bIX+1, 0.0);
	jtptPPMC_RAW_h[cI][nAbsEtaBins][sI]->SetBinError(bIX+1, 0.0);
	jtptPPMC_RAW_Frac_h[cI][nAbsEtaBins][sI]->SetBinError(bIX+1, 0.0);
      }
      
      for(Int_t bI = 0; bI < nBayesIter; ++bI){
	for(Int_t bIX = 0; bIX < jtptPPData_FULLBayes_h[cI][nAbsEtaBins][sI][bI]->GetNbinsX(); ++bIX){
	  jtptPPData_FULLBayes_h[cI][nAbsEtaBins][sI][bI]->SetBinContent(bIX+1, 0.0);
	  jtptPPMC_FULLBayes_h[cI][nAbsEtaBins][sI][bI]->SetBinContent(bIX+1, 0.0);
	  jtptPPMC_FULLBayes_Frac_h[cI][nAbsEtaBins][sI][bI]->SetBinContent(bIX+1, 0.0);
	  
	  jtptPPData_FULLBayes_h[cI][nAbsEtaBins][sI][bI]->SetBinError(bIX+1, 0.0);
	  jtptPPMC_FULLBayes_h[cI][nAbsEtaBins][sI][bI]->SetBinError(bIX+1, 0.0);
	  jtptPPMC_FULLBayes_Frac_h[cI][nAbsEtaBins][sI][bI]->SetBinError(bIX+1, 0.0);
	}	  
      }
      
      if(doSvd){
	for(Int_t bI = 1; bI < nSvd; ++bI){
	  for(Int_t bIX = 0; bIX < jtptPPData_FULLSvd_h[cI][nAbsEtaBins][sI][bI]->GetNbinsX(); ++bIX){
	    jtptPPData_FULLSvd_h[cI][nAbsEtaBins][sI][bI]->SetBinContent(bIX+1, 0.0);
	    jtptPPMC_FULLSvd_h[cI][nAbsEtaBins][sI][bI]->SetBinContent(bIX+1, 0.0);
	    jtptPPMC_FULLSvd_Frac_h[cI][nAbsEtaBins][sI][bI]->SetBinContent(bIX+1, 0.0);
	    
	    jtptPPData_FULLSvd_h[cI][nAbsEtaBins][sI][bI]->SetBinError(bIX+1, 0.0);
	    jtptPPMC_FULLSvd_h[cI][nAbsEtaBins][sI][bI]->SetBinError(bIX+1, 0.0);
	    jtptPPMC_FULLSvd_Frac_h[cI][nAbsEtaBins][sI][bI]->SetBinError(bIX+1, 0.0);
	  }	  
	}
      }
    }
    
    for(Int_t aI = 0; aI < nAbsEtaBins-1; ++aI){
      for(Int_t bIX = 0; bIX < refptPbPbMC_h[cI][nAbsEtaBins]->GetNbinsX(); ++bIX){
	Double_t mcVal = getComboBinVal(refptPbPbMC_h[cI][nAbsEtaBins], refptPbPbMC_h[cI][aI], bIX+1);
	Double_t mcFracVal = getComboBinVal(refptPbPbMC_Frac_h[cI][nAbsEtaBins], refptPbPbMC_Frac_h[cI][aI], bIX+1);
	Double_t mcErr = getComboBinErr(refptPbPbMC_h[cI][nAbsEtaBins], refptPbPbMC_h[cI][aI], bIX+1);
	Double_t mcFracErr = getComboBinErr(refptPbPbMC_Frac_h[cI][nAbsEtaBins], refptPbPbMC_Frac_h[cI][aI], bIX+1);
	
	refptPbPbMC_h[cI][nAbsEtaBins]->SetBinContent(bIX+1, mcVal);
	refptPbPbMC_Frac_h[cI][nAbsEtaBins]->SetBinContent(bIX+1, mcFracVal);
	
	refptPbPbMC_h[cI][nAbsEtaBins]->SetBinError(bIX+1, mcErr);
	refptPbPbMC_Frac_h[cI][nAbsEtaBins]->SetBinError(bIX+1, mcFracErr);
      }
     
      for(Int_t sI = 0; sI < nSyst; ++sI){
	for(Int_t bIX = 0; bIX < jtptPbPbData_RAW_h[cI][nAbsEtaBins][sI]->GetNbinsX(); ++bIX){
	  Double_t dataVal = getComboBinVal(jtptPbPbData_RAW_h[cI][nAbsEtaBins][sI], jtptPbPbData_RAW_h[cI][aI][sI], bIX+1);
	  Double_t mcVal = getComboBinVal(jtptPbPbMC_RAW_h[cI][nAbsEtaBins][sI], jtptPbPbMC_RAW_h[cI][aI][sI], bIX+1);
	  Double_t mcFracVal = getComboBinVal(jtptPbPbMC_RAW_Frac_h[cI][nAbsEtaBins][sI], jtptPbPbMC_RAW_Frac_h[cI][aI][sI], bIX+1);
	  Double_t dataErr = getComboBinErr(jtptPbPbData_RAW_h[cI][nAbsEtaBins][sI], jtptPbPbData_RAW_h[cI][aI][sI], bIX+1);
	  Double_t mcErr = getComboBinErr(jtptPbPbMC_RAW_h[cI][nAbsEtaBins][sI], jtptPbPbMC_RAW_h[cI][aI][sI], bIX+1);
	  Double_t mcFracErr = getComboBinErr(jtptPbPbMC_RAW_Frac_h[cI][nAbsEtaBins][sI], jtptPbPbMC_RAW_Frac_h[cI][aI][sI], bIX+1);
	  
	  
	  jtptPbPbData_RAW_h[cI][nAbsEtaBins][sI]->SetBinContent(bIX+1, dataVal);
	  jtptPbPbMC_RAW_h[cI][nAbsEtaBins][sI]->SetBinContent(bIX+1, mcVal);
	  jtptPbPbMC_RAW_Frac_h[cI][nAbsEtaBins][sI]->SetBinContent(bIX+1, mcFracVal);
	  
	  jtptPbPbData_RAW_h[cI][nAbsEtaBins][sI]->SetBinError(bIX+1, dataErr);
	  jtptPbPbMC_RAW_h[cI][nAbsEtaBins][sI]->SetBinError(bIX+1, mcErr);
	  jtptPbPbMC_RAW_Frac_h[cI][nAbsEtaBins][sI]->SetBinError(bIX+1, mcFracErr);
	}
	
	for(Int_t bI = 0; bI < nBayesIter; ++bI){
	  for(Int_t bIX = 0; bIX < jtptPbPbData_FULLBayes_h[cI][nAbsEtaBins][sI][bI]->GetNbinsX(); ++bIX){
	    Double_t dataVal = getComboBinVal(jtptPbPbData_FULLBayes_h[cI][nAbsEtaBins][sI][bI], jtptPbPbData_FULLBayes_h[cI][aI][sI][bI], bIX+1);
	    Double_t mcVal = getComboBinVal(jtptPbPbMC_FULLBayes_h[cI][nAbsEtaBins][sI][bI], jtptPbPbMC_FULLBayes_h[cI][aI][sI][bI], bIX+1);
	    Double_t mcFracVal = getComboBinVal(jtptPbPbMC_FULLBayes_Frac_h[cI][nAbsEtaBins][sI][bI], jtptPbPbMC_FULLBayes_Frac_h[cI][aI][sI][bI], bIX+1);
	    Double_t dataErr = getComboBinErr(jtptPbPbData_FULLBayes_h[cI][nAbsEtaBins][sI][bI], jtptPbPbData_FULLBayes_h[cI][aI][sI][bI], bIX+1);
	    Double_t mcErr = getComboBinErr(jtptPbPbMC_FULLBayes_h[cI][nAbsEtaBins][sI][bI], jtptPbPbMC_FULLBayes_h[cI][aI][sI][bI], bIX+1);
	    Double_t mcFracErr = getComboBinErr(jtptPbPbMC_FULLBayes_Frac_h[cI][nAbsEtaBins][sI][bI], jtptPbPbMC_FULLBayes_Frac_h[cI][aI][sI][bI], bIX+1);
	    
	    jtptPbPbData_FULLBayes_h[cI][nAbsEtaBins][sI][bI]->SetBinContent(bIX+1, dataVal);
	    jtptPbPbMC_FULLBayes_h[cI][nAbsEtaBins][sI][bI]->SetBinContent(bIX+1, mcVal);
	    jtptPbPbMC_FULLBayes_Frac_h[cI][nAbsEtaBins][sI][bI]->SetBinContent(bIX+1, mcFracVal);
	    
	    jtptPbPbData_FULLBayes_h[cI][nAbsEtaBins][sI][bI]->SetBinError(bIX+1, dataErr);
	    jtptPbPbMC_FULLBayes_h[cI][nAbsEtaBins][sI][bI]->SetBinError(bIX+1, mcErr);
	    jtptPbPbMC_FULLBayes_Frac_h[cI][nAbsEtaBins][sI][bI]->SetBinError(bIX+1, mcFracErr);
	  }
	}
	
	if(doSvd){
	  for(Int_t bI = 1; bI < nSvd; ++bI){
	    for(Int_t bIX = 0; bIX < jtptPbPbData_FULLSvd_h[cI][nAbsEtaBins][sI][bI]->GetNbinsX(); ++bIX){
	      Double_t dataVal = getComboBinVal(jtptPbPbData_FULLSvd_h[cI][nAbsEtaBins][sI][bI], jtptPbPbData_FULLSvd_h[cI][aI][sI][bI], bIX+1);
	      Double_t mcVal = getComboBinVal(jtptPbPbMC_FULLSvd_h[cI][nAbsEtaBins][sI][bI], jtptPbPbMC_FULLSvd_h[cI][aI][sI][bI], bIX+1);
	      Double_t mcFracVal = getComboBinVal(jtptPbPbMC_FULLSvd_Frac_h[cI][nAbsEtaBins][sI][bI], jtptPbPbMC_FULLSvd_Frac_h[cI][aI][sI][bI], bIX+1);
	      Double_t dataErr = getComboBinErr(jtptPbPbData_FULLSvd_h[cI][nAbsEtaBins][sI][bI], jtptPbPbData_FULLSvd_h[cI][aI][sI][bI], bIX+1);
	      Double_t mcErr = getComboBinErr(jtptPbPbMC_FULLSvd_h[cI][nAbsEtaBins][sI][bI], jtptPbPbMC_FULLSvd_h[cI][aI][sI][bI], bIX+1);
	      Double_t mcFracErr = getComboBinErr(jtptPbPbMC_FULLSvd_Frac_h[cI][nAbsEtaBins][sI][bI], jtptPbPbMC_FULLSvd_Frac_h[cI][aI][sI][bI], bIX+1);
	      
	      jtptPbPbData_FULLSvd_h[cI][nAbsEtaBins][sI][bI]->SetBinContent(bIX+1, dataVal);
	      jtptPbPbMC_FULLSvd_h[cI][nAbsEtaBins][sI][bI]->SetBinContent(bIX+1, mcVal);
	      jtptPbPbMC_FULLSvd_Frac_h[cI][nAbsEtaBins][sI][bI]->SetBinContent(bIX+1, mcFracVal);
	      
	      jtptPbPbData_FULLSvd_h[cI][nAbsEtaBins][sI][bI]->SetBinError(bIX+1, dataErr);
	      jtptPbPbMC_FULLSvd_h[cI][nAbsEtaBins][sI][bI]->SetBinError(bIX+1, mcErr);
	      jtptPbPbMC_FULLSvd_Frac_h[cI][nAbsEtaBins][sI][bI]->SetBinError(bIX+1, mcFracErr);
	    }
	  }
	}
      }
      
      
      for(Int_t bIX = 0; bIX < refptPPMC_h[cI][nAbsEtaBins]->GetNbinsX(); ++bIX){
	Double_t mcVal = getComboBinVal(refptPPMC_h[cI][nAbsEtaBins], refptPPMC_h[cI][aI], bIX+1);
	Double_t mcFracVal = getComboBinVal(refptPPMC_Frac_h[cI][nAbsEtaBins], refptPPMC_Frac_h[cI][aI], bIX+1);
	Double_t mcErr = getComboBinErr(refptPPMC_h[cI][nAbsEtaBins], refptPPMC_h[cI][aI], bIX+1);
	Double_t mcFracErr = getComboBinErr(refptPPMC_Frac_h[cI][nAbsEtaBins], refptPPMC_Frac_h[cI][aI], bIX+1);
	
	refptPPMC_h[cI][nAbsEtaBins]->SetBinContent(bIX+1, mcVal);
	refptPPMC_Frac_h[cI][nAbsEtaBins]->SetBinContent(bIX+1, mcFracVal);
	
	refptPPMC_h[cI][nAbsEtaBins]->SetBinError(bIX+1, mcErr);
	refptPPMC_Frac_h[cI][nAbsEtaBins]->SetBinError(bIX+1, mcFracErr);
      }
      
      
      for(Int_t sI = 0; sI < nSyst; ++sI){
	for(Int_t bIX = 0; bIX < jtptPPData_RAW_h[cI][nAbsEtaBins][sI]->GetNbinsX(); ++bIX){
	  Double_t dataVal = getComboBinVal(jtptPPData_RAW_h[cI][nAbsEtaBins][sI], jtptPPData_RAW_h[cI][aI][sI], bIX+1);
	  Double_t mcVal = getComboBinVal(jtptPPMC_RAW_h[cI][nAbsEtaBins][sI], jtptPPMC_RAW_h[cI][aI][sI], bIX+1);
	  Double_t mcFracVal = getComboBinVal(jtptPPMC_RAW_Frac_h[cI][nAbsEtaBins][sI], jtptPPMC_RAW_Frac_h[cI][aI][sI], bIX+1);
	  Double_t dataErr = getComboBinErr(jtptPPData_RAW_h[cI][nAbsEtaBins][sI], jtptPPData_RAW_h[cI][aI][sI], bIX+1);
	  Double_t mcErr = getComboBinErr(jtptPPMC_RAW_h[cI][nAbsEtaBins][sI], jtptPPMC_RAW_h[cI][aI][sI], bIX+1);
	  Double_t mcFracErr = getComboBinErr(jtptPPMC_RAW_Frac_h[cI][nAbsEtaBins][sI], jtptPPMC_RAW_Frac_h[cI][aI][sI], bIX+1);
	    
	  
	  jtptPPData_RAW_h[cI][nAbsEtaBins][sI]->SetBinContent(bIX+1, dataVal);
	  jtptPPMC_RAW_h[cI][nAbsEtaBins][sI]->SetBinContent(bIX+1, mcVal);
	  jtptPPMC_RAW_Frac_h[cI][nAbsEtaBins][sI]->SetBinContent(bIX+1, mcFracVal);

	  jtptPPData_RAW_h[cI][nAbsEtaBins][sI]->SetBinError(bIX+1, dataErr);
	  jtptPPMC_RAW_h[cI][nAbsEtaBins][sI]->SetBinError(bIX+1, mcErr);
	  jtptPPMC_RAW_Frac_h[cI][nAbsEtaBins][sI]->SetBinError(bIX+1, mcFracErr);
	}
	
	for(Int_t bI = 0; bI < nBayesIter; ++bI){
	  for(Int_t bIX = 0; bIX < jtptPPData_FULLBayes_h[cI][nAbsEtaBins][sI][bI]->GetNbinsX(); ++bIX){
	    Double_t dataVal = getComboBinVal(jtptPPData_FULLBayes_h[cI][nAbsEtaBins][sI][bI], jtptPPData_FULLBayes_h[cI][aI][sI][bI], bIX+1);
	    Double_t mcVal = getComboBinVal(jtptPPMC_FULLBayes_h[cI][nAbsEtaBins][sI][bI], jtptPPMC_FULLBayes_h[cI][aI][sI][bI], bIX+1);
	    Double_t mcFracVal = getComboBinVal(jtptPPMC_FULLBayes_Frac_h[cI][nAbsEtaBins][sI][bI], jtptPPMC_FULLBayes_Frac_h[cI][aI][sI][bI], bIX+1);
	    Double_t dataErr = getComboBinErr(jtptPPData_FULLBayes_h[cI][nAbsEtaBins][sI][bI], jtptPPData_FULLBayes_h[cI][aI][sI][bI], bIX+1);
	    Double_t mcErr = getComboBinErr(jtptPPMC_FULLBayes_h[cI][nAbsEtaBins][sI][bI], jtptPPMC_FULLBayes_h[cI][aI][sI][bI], bIX+1);
	    Double_t mcFracErr = getComboBinErr(jtptPPMC_FULLBayes_Frac_h[cI][nAbsEtaBins][sI][bI], jtptPPMC_FULLBayes_Frac_h[cI][aI][sI][bI], bIX+1);
	    
	    jtptPPData_FULLBayes_h[cI][nAbsEtaBins][sI][bI]->SetBinContent(bIX+1, dataVal);
	    jtptPPMC_FULLBayes_h[cI][nAbsEtaBins][sI][bI]->SetBinContent(bIX+1, mcVal);
	    jtptPPMC_FULLBayes_Frac_h[cI][nAbsEtaBins][sI][bI]->SetBinContent(bIX+1, mcFracVal);
	    
	    jtptPPData_FULLBayes_h[cI][nAbsEtaBins][sI][bI]->SetBinError(bIX+1, dataErr);
	    jtptPPMC_FULLBayes_h[cI][nAbsEtaBins][sI][bI]->SetBinError(bIX+1, mcErr);
	    jtptPPMC_FULLBayes_Frac_h[cI][nAbsEtaBins][sI][bI]->SetBinError(bIX+1, mcFracErr);
	  }
	}
	
	if(doSvd){
	  for(Int_t bI = 1; bI < nSvd; ++bI){
	    for(Int_t bIX = 0; bIX < jtptPPData_FULLSvd_h[cI][nAbsEtaBins][sI][bI]->GetNbinsX(); ++bIX){
	      Double_t dataVal = getComboBinVal(jtptPPData_FULLSvd_h[cI][nAbsEtaBins][sI][bI], jtptPPData_FULLSvd_h[cI][aI][sI][bI], bIX+1);
	      Double_t mcVal = getComboBinVal(jtptPPMC_FULLSvd_h[cI][nAbsEtaBins][sI][bI], jtptPPMC_FULLSvd_h[cI][aI][sI][bI], bIX+1);
	      Double_t mcFracVal = getComboBinVal(jtptPPMC_FULLSvd_Frac_h[cI][nAbsEtaBins][sI][bI], jtptPPMC_FULLSvd_Frac_h[cI][aI][sI][bI], bIX+1);
	      Double_t dataErr = getComboBinErr(jtptPPData_FULLSvd_h[cI][nAbsEtaBins][sI][bI], jtptPPData_FULLSvd_h[cI][aI][sI][bI], bIX+1);
	      Double_t mcErr = getComboBinErr(jtptPPMC_FULLSvd_h[cI][nAbsEtaBins][sI][bI], jtptPPMC_FULLSvd_h[cI][aI][sI][bI], bIX+1);
	      Double_t mcFracErr = getComboBinErr(jtptPPMC_FULLSvd_Frac_h[cI][nAbsEtaBins][sI][bI], jtptPPMC_FULLSvd_Frac_h[cI][aI][sI][bI], bIX+1);
	      
	      jtptPPData_FULLSvd_h[cI][nAbsEtaBins][sI][bI]->SetBinContent(bIX+1, dataVal);
	      jtptPPMC_FULLSvd_h[cI][nAbsEtaBins][sI][bI]->SetBinContent(bIX+1, mcVal);
	      jtptPPMC_FULLSvd_Frac_h[cI][nAbsEtaBins][sI][bI]->SetBinContent(bIX+1, mcFracVal);
	      
	      jtptPPData_FULLSvd_h[cI][nAbsEtaBins][sI][bI]->SetBinError(bIX+1, dataErr);
	      jtptPPMC_FULLSvd_h[cI][nAbsEtaBins][sI][bI]->SetBinError(bIX+1, mcErr);
	      jtptPPMC_FULLSvd_Frac_h[cI][nAbsEtaBins][sI][bI]->SetBinError(bIX+1, mcFracErr);
	    }
	  }
	}
      }      
    }      
  }

  finalEtaWatch.stop();
  dumpTiming({startWatch, histWatch, pbpbDataWatch, ppDataWatch, ppMCWatch, pbpbMCWatch, countWatch, unfoldLoopWatch, plotValidWatch, copyWatch, finalEtaWatch}, {"start", "hist", "pbpbData", "ppData", "ppMC", "pbpbMC", "count", "unfold", "plotValid", "copy", "finalEta"});
  makeEtaHistWatch.start();
  
  std::cout << "ETA HISTS" << std::endl;
  //Setting eta hists
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    for(Int_t sI = 0; sI < nSyst; ++sI){
      for(Int_t bI = 0; bI < nBayesIter; ++bI){
	
	std::vector<TH1D*> jtptPbPbData_Vect;
	std::vector<TH1D*> jtptPPData_Vect;
	std::vector<TH1D*> jtptPbPbMC_Vect;
	std::vector<TH1D*> jtptPPMC_Vect;
	std::vector<TH1D*> jtptPbPbMC_Frac_Vect;
	std::vector<TH1D*> jtptPPMC_Frac_Vect;
	for(Int_t aI = 0; aI < jtetaPbPbData_FULLBayes_h[cI][sI][bI]->GetNbinsX(); ++aI){
	  jtptPbPbData_Vect.push_back(jtptPbPbData_FULLBayes_h[cI][aI][sI][bI]);
	  jtptPPData_Vect.push_back(jtptPPData_FULLBayes_h[cI][aI][sI][bI]);
	  jtptPbPbMC_Vect.push_back(jtptPbPbMC_FULLBayes_h[cI][aI][sI][bI]);
	  jtptPPMC_Vect.push_back(jtptPPMC_FULLBayes_h[cI][aI][sI][bI]);
	  jtptPbPbMC_Frac_Vect.push_back(jtptPbPbMC_FULLBayes_Frac_h[cI][aI][sI][bI]);
	  jtptPPMC_Frac_Vect.push_back(jtptPPMC_FULLBayes_Frac_h[cI][aI][sI][bI]);
	}
	
	makeEtaHist(jtetaPbPbData_FULLBayes_h[cI][sI][bI], jtptPbPbData_Vect);
	makeEtaHist(jtetaPPData_FULLBayes_h[cI][sI][bI], jtptPPData_Vect);
	makeEtaHist(jtetaPbPbMC_FULLBayes_h[cI][sI][bI], jtptPbPbMC_Vect);
	makeEtaHist(jtetaPPMC_FULLBayes_h[cI][sI][bI], jtptPPMC_Vect);
	makeEtaHist(jtetaPbPbMC_FULLBayes_Frac_h[cI][sI][bI], jtptPbPbMC_Frac_Vect);
	makeEtaHist(jtetaPPMC_FULLBayes_Frac_h[cI][sI][bI], jtptPPMC_Frac_Vect);
		
	jtptPbPbData_Vect.clear();
	jtptPPData_Vect.clear();
	jtptPbPbMC_Vect.clear();
	jtptPPMC_Vect.clear();
	jtptPbPbMC_Frac_Vect.clear();
	jtptPPMC_Frac_Vect.clear();
      }
    }    
  }

  makeEtaHistWatch.stop();
  dumpTiming({startWatch, histWatch, pbpbDataWatch, ppDataWatch, ppMCWatch, pbpbMCWatch, countWatch, unfoldLoopWatch, plotValidWatch, copyWatch, finalEtaWatch, makeEtaHistWatch}, {"start", "hist", "pbpbData", "ppData", "ppMC", "pbpbMC", "count", "unfold", "plotValid", "copy", "finalEta", "makeEtaHist"});
  writeWatch.start();
  std::cout << "WRITING" << std::endl;

  outFile_p->cd();
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    const std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
    const std::string dirStrPbPb = algoStrPbPb + "/" + centStr;
    const std::string dirStrPP = algoStrPP + "/" + centStr;
    
    if(doGlobalDebug) std::cout << cI << "/" << nCentBins << ", " << __FILE__ << ", " << __LINE__ << std::endl;
    
    outFile_p->cd();
    jetPbPbDirs_p = outFile_p->GetDirectory(algoStrPbPb.c_str());
    if(jetPbPbDirs_p) jetPbPbDirs_p->cd();
    else{
      jetPbPbDirs_p = outFile_p->mkdir(algoStrPbPb.c_str());
      jetPbPbDirs_p->cd();
    }
    
    outFile_p->cd();
    centPbPbDirs_p[cI] = outFile_p->GetDirectory(dirStrPbPb.c_str());
    if(centPbPbDirs_p[cI]) centPbPbDirs_p[cI]->cd();
    else{
      centPbPbDirs_p[cI] = outFile_p->mkdir(dirStrPbPb.c_str());
      centPbPbDirs_p[cI]->cd();
    }
    
    for(Int_t aI = 0; aI < nAbsEtaBins+1; ++aI){      
      if(doGlobalDebug) std::cout << cI << "/" << nCentBins << ", " << aI << "/" << nAbsEtaBins << ", " << __FILE__ << ", " << __LINE__ << std::endl;
      
      for(Int_t sI = 0; sI < nSyst; ++sI){
	if(doGlobalDebug) std::cout << cI << "/" << nCentBins << ", " << aI << "/" << nAbsEtaBins << ", " << sI << "/" << nSyst << ", " << __FILE__ << ", " << __LINE__ << std::endl;
	
	
	jtptPbPbData_RAW_h[cI][aI][sI]->Write("", TObject::kOverwrite);
	jtptPbPbMC_RAW_h[cI][aI][sI]->Write("", TObject::kOverwrite);
	jtptPbPbMC_RAW_Frac_h[cI][aI][sI]->Write("", TObject::kOverwrite);
	
	for(Int_t bI = 0; bI < nBayesIter; ++bI){
	  jtptPbPbData_FULLBayes_h[cI][aI][sI][bI]->Write("", TObject::kOverwrite);
	  jtptPbPbMC_FULLBayes_h[cI][aI][sI][bI]->Write("", TObject::kOverwrite);
	  jtptPbPbMC_FULLBayes_Frac_h[cI][aI][sI][bI]->Write("", TObject::kOverwrite);
	  
	  if(sI < nSystForFill && aI < nAbsEtaBins){
	    rooUnfoldBayesPbPbData_FULLBayes_h[cI][aI][sI][bI]->Write("", TObject::kOverwrite);
	    rooUnfoldBayesPbPbMC_FULLBayes_h[cI][aI][sI][bI]->Write("", TObject::kOverwrite);
	    rooUnfoldBayesPbPbMC_FULLBayes_Frac_h[cI][aI][sI][bI]->Write("", TObject::kOverwrite);
	  }
	}
		
	if(doSvd){
	  for(Int_t bI = 1; bI < nSvd; ++bI){
	    if(bI+1 > truthNBins[cI]) break;
	    jtptPbPbData_FULLSvd_h[cI][aI][sI][bI]->Write("", TObject::kOverwrite);
	    jtptPbPbMC_FULLSvd_h[cI][aI][sI][bI]->Write("", TObject::kOverwrite);
	    jtptPbPbMC_FULLSvd_Frac_h[cI][aI][sI][bI]->Write("", TObject::kOverwrite);

	    if(sI < nSystForFill && aI < nAbsEtaBins){
	      rooUnfoldSvdPbPbData_FULLSvd_h[cI][aI][sI][bI]->Write("", TObject::kOverwrite);
	      rooUnfoldSvdPbPbMC_FULLSvd_h[cI][aI][sI][bI]->Write("", TObject::kOverwrite);
	      rooUnfoldSvdPbPbMC_FULLSvd_Frac_h[cI][aI][sI][bI]->Write("", TObject::kOverwrite);

	      std::string svdDVPbPbData = rooUnfoldSvdPbPbData_FULLSvd_h[cI][aI][sI][bI]->GetName();
	      std::string svdDVPbPbMC = rooUnfoldSvdPbPbMC_FULLSvd_h[cI][aI][sI][bI]->GetName();
	      std::string svdDVPbPbMCFrac= rooUnfoldSvdPbPbMC_FULLSvd_Frac_h[cI][aI][sI][bI]->GetName();

	      svdDVPbPbData = svdDVPbPbData + "_DV";
	      svdDVPbPbMC = svdDVPbPbMC + "_DV";
	      svdDVPbPbMCFrac = svdDVPbPbMCFrac + "_DV";

	      TH1D* rooUnfoldSvdPbPbData_FULLSvd_DV_h = (TH1D*)rooUnfoldSvdPbPbData_FULLSvd_h[cI][aI][sI][bI]->Impl()->GetD()->Clone(svdDVPbPbData.c_str());
	      TH1D* rooUnfoldSvdPbPbMC_FULLSvd_DV_h = (TH1D*)rooUnfoldSvdPbPbMC_FULLSvd_h[cI][aI][sI][bI]->Impl()->GetD()->Clone(svdDVPbPbMC.c_str());
	      TH1D* rooUnfoldSvdPbPbMCFrac_FULLSvd_DV_h = (TH1D*)rooUnfoldSvdPbPbMC_FULLSvd_Frac_h[cI][aI][sI][bI]->Impl()->GetD()->Clone(svdDVPbPbMCFrac.c_str());

	      rooUnfoldSvdPbPbData_FULLSvd_DV_h->Write("", TObject::kOverwrite);
	      rooUnfoldSvdPbPbMC_FULLSvd_DV_h->Write("", TObject::kOverwrite);
	      rooUnfoldSvdPbPbMCFrac_FULLSvd_DV_h->Write("", TObject::kOverwrite);

	      delete rooUnfoldSvdPbPbData_FULLSvd_DV_h;
	      delete rooUnfoldSvdPbPbMC_FULLSvd_DV_h;
	      delete rooUnfoldSvdPbPbMCFrac_FULLSvd_DV_h;


	      std::string svdSVPbPbData = rooUnfoldSvdPbPbData_FULLSvd_h[cI][aI][sI][bI]->GetName();
	      std::string svdSVPbPbMC = rooUnfoldSvdPbPbMC_FULLSvd_h[cI][aI][sI][bI]->GetName();
	      std::string svdSVPbPbMCFrac= rooUnfoldSvdPbPbMC_FULLSvd_Frac_h[cI][aI][sI][bI]->GetName();

	      svdSVPbPbData = svdSVPbPbData + "_SV";
	      svdSVPbPbMC = svdSVPbPbMC + "_SV";
	      svdSVPbPbMCFrac = svdSVPbPbMCFrac + "_SV";

	      TH1D* rooUnfoldSvdPbPbData_FULLSvd_SV_h = (TH1D*)rooUnfoldSvdPbPbData_FULLSvd_h[cI][aI][sI][bI]->Impl()->GetSV()->Clone(svdSVPbPbData.c_str());
	      TH1D* rooUnfoldSvdPbPbMC_FULLSvd_SV_h = (TH1D*)rooUnfoldSvdPbPbMC_FULLSvd_h[cI][aI][sI][bI]->Impl()->GetSV()->Clone(svdSVPbPbMC.c_str());
	      TH1D* rooUnfoldSvdPbPbMCFrac_FULLSvd_SV_h = (TH1D*)rooUnfoldSvdPbPbMC_FULLSvd_Frac_h[cI][aI][sI][bI]->Impl()->GetSV()->Clone(svdSVPbPbMCFrac.c_str());

	      rooUnfoldSvdPbPbData_FULLSvd_SV_h->Write("", TObject::kOverwrite);
	      rooUnfoldSvdPbPbMC_FULLSvd_SV_h->Write("", TObject::kOverwrite);
	      rooUnfoldSvdPbPbMCFrac_FULLSvd_SV_h->Write("", TObject::kOverwrite);

	      delete rooUnfoldSvdPbPbData_FULLSvd_SV_h;
	      delete rooUnfoldSvdPbPbMC_FULLSvd_SV_h;
	      delete rooUnfoldSvdPbPbMCFrac_FULLSvd_SV_h;
	    }
	  }	   
	}
      }   
      
      if(doGlobalDebug) std::cout << __FILE__ << ", " << __FILE__ << ", " << cI << ", " << aI << std::endl;
      
      refptPbPbMC_h[cI][aI]->Write("", TObject::kOverwrite);
      refptPbPbMC_Frac_h[cI][aI]->Write("", TObject::kOverwrite);
      responsePbPbMC_h[cI][aI]->Write("", TObject::kOverwrite);
      responsePbPbMC_Symm_h[cI][aI]->Write("", TObject::kOverwrite);
      rooResponsePbPbMC_h[cI][aI]->Write("", TObject::kOverwrite);
      
      if(doGlobalDebug) std::cout << __FILE__ << ", " << __FILE__ << ", " << cI << ", " << aI << std::endl;
    }
    
    for(Int_t sI = 0; sI < nSyst; ++sI){
      for(Int_t bI = 0; bI < nBayesIter; ++bI){
	jtetaPbPbData_FULLBayes_h[cI][sI][bI]->Write("", TObject::kOverwrite);
	jtetaPbPbMC_FULLBayes_h[cI][sI][bI]->Write("", TObject::kOverwrite);
	jtetaPbPbMC_FULLBayes_Frac_h[cI][sI][bI]->Write("", TObject::kOverwrite);
      }
    }
    
    if(doGlobalDebug) std::cout << cI << "/" << nCentBins << ", " << __FILE__ << ", " << __LINE__ << std::endl;
    
    outFile_p->cd();
    jetPPDirs_p = outFile_p->GetDirectory(algoStrPP.c_str());
    if(jetPPDirs_p) jetPPDirs_p->cd();
    else{
      jetPPDirs_p = outFile_p->mkdir(algoStrPP.c_str());
      jetPPDirs_p->cd();
    }
    
    outFile_p->cd();
    centPPDirs_p[cI] = outFile_p->GetDirectory(dirStrPP.c_str());
    if(centPPDirs_p[cI]) centPPDirs_p[cI]->cd();
    else{
      centPPDirs_p[cI] = outFile_p->mkdir(dirStrPP.c_str());
      centPPDirs_p[cI]->cd();
    }
    
    
    for(Int_t aI = 0; aI < nAbsEtaBins+1; ++aI){
      for(Int_t sI = 0; sI < nSyst; ++sI){
	jtptPPData_RAW_h[cI][aI][sI]->Write("", TObject::kOverwrite);
	jtptPPMC_RAW_h[cI][aI][sI]->Write("", TObject::kOverwrite);
	jtptPPMC_RAW_Frac_h[cI][aI][sI]->Write("", TObject::kOverwrite);
	
	for(Int_t bI = 0; bI < nBayesIter; ++bI){
	  jtptPPData_FULLBayes_h[cI][aI][sI][bI]->Write("", TObject::kOverwrite);
	  jtptPPMC_FULLBayes_h[cI][aI][sI][bI]->Write("", TObject::kOverwrite);
	  jtptPPMC_FULLBayes_Frac_h[cI][aI][sI][bI]->Write("", TObject::kOverwrite);
	  
	  if(sI < nSystForFill && aI < nAbsEtaBins){
	    rooUnfoldBayesPPData_FULLBayes_h[cI][aI][sI][bI]->Write("", TObject::kOverwrite);
	    rooUnfoldBayesPPMC_FULLBayes_h[cI][aI][sI][bI]->Write("", TObject::kOverwrite);
	    rooUnfoldBayesPPMC_FULLBayes_Frac_h[cI][aI][sI][bI]->Write("", TObject::kOverwrite);
	  }
	}
	
	if(doSvd){
	  for(Int_t bI = 1; bI < nSvd; ++bI){
	    if(bI+1 > truthNBins[cI]) break;

	    jtptPPData_FULLSvd_h[cI][aI][sI][bI]->Write("", TObject::kOverwrite);
	    jtptPPMC_FULLSvd_h[cI][aI][sI][bI]->Write("", TObject::kOverwrite);
	    jtptPPMC_FULLSvd_Frac_h[cI][aI][sI][bI]->Write("", TObject::kOverwrite);

	    if(sI < nSystForFill && aI < nAbsEtaBins){
	      rooUnfoldSvdPPData_FULLSvd_h[cI][aI][sI][bI]->Write("", TObject::kOverwrite);
	      rooUnfoldSvdPPMC_FULLSvd_h[cI][aI][sI][bI]->Write("", TObject::kOverwrite);
	      rooUnfoldSvdPPMC_FULLSvd_Frac_h[cI][aI][sI][bI]->Write("", TObject::kOverwrite);

	      std::string svdDVPPData = rooUnfoldSvdPPData_FULLSvd_h[cI][aI][sI][bI]->GetName();
	      std::string svdDVPPMC = rooUnfoldSvdPPMC_FULLSvd_h[cI][aI][sI][bI]->GetName();
	      std::string svdDVPPMCFrac= rooUnfoldSvdPPMC_FULLSvd_Frac_h[cI][aI][sI][bI]->GetName();

	      svdDVPPData = svdDVPPData + "_DV";
	      svdDVPPMC = svdDVPPMC + "_DV";
	      svdDVPPMCFrac = svdDVPPMCFrac + "_DV";

	      TH1D* rooUnfoldSvdPPData_FULLSvd_DV_h = (TH1D*)rooUnfoldSvdPPData_FULLSvd_h[cI][aI][sI][bI]->Impl()->GetD()->Clone(svdDVPPData.c_str());
	      TH1D* rooUnfoldSvdPPMC_FULLSvd_DV_h = (TH1D*)rooUnfoldSvdPPMC_FULLSvd_h[cI][aI][sI][bI]->Impl()->GetD()->Clone(svdDVPPMC.c_str());
	      TH1D* rooUnfoldSvdPPMCFrac_FULLSvd_DV_h = (TH1D*)rooUnfoldSvdPPMC_FULLSvd_Frac_h[cI][aI][sI][bI]->Impl()->GetD()->Clone(svdDVPPMCFrac.c_str());

	      rooUnfoldSvdPPData_FULLSvd_DV_h->Write("", TObject::kOverwrite);
	      rooUnfoldSvdPPMC_FULLSvd_DV_h->Write("", TObject::kOverwrite);
	      rooUnfoldSvdPPMCFrac_FULLSvd_DV_h->Write("", TObject::kOverwrite);

	      delete rooUnfoldSvdPPData_FULLSvd_DV_h;
	      delete rooUnfoldSvdPPMC_FULLSvd_DV_h;
	      delete rooUnfoldSvdPPMCFrac_FULLSvd_DV_h;


	      std::string svdSVPPData = rooUnfoldSvdPPData_FULLSvd_h[cI][aI][sI][bI]->GetName();
	      std::string svdSVPPMC = rooUnfoldSvdPPMC_FULLSvd_h[cI][aI][sI][bI]->GetName();
	      std::string svdSVPPMCFrac= rooUnfoldSvdPPMC_FULLSvd_Frac_h[cI][aI][sI][bI]->GetName();

	      svdSVPPData = svdSVPPData + "_SV";
	      svdSVPPMC = svdSVPPMC + "_SV";
	      svdSVPPMCFrac = svdSVPPMCFrac + "_SV";

	      TH1D* rooUnfoldSvdPPData_FULLSvd_SV_h = (TH1D*)rooUnfoldSvdPPData_FULLSvd_h[cI][aI][sI][bI]->Impl()->GetSV()->Clone(svdSVPPData.c_str());
	      TH1D* rooUnfoldSvdPPMC_FULLSvd_SV_h = (TH1D*)rooUnfoldSvdPPMC_FULLSvd_h[cI][aI][sI][bI]->Impl()->GetSV()->Clone(svdSVPPMC.c_str());
	      TH1D* rooUnfoldSvdPPMCFrac_FULLSvd_SV_h = (TH1D*)rooUnfoldSvdPPMC_FULLSvd_Frac_h[cI][aI][sI][bI]->Impl()->GetSV()->Clone(svdSVPPMCFrac.c_str());

	      rooUnfoldSvdPPData_FULLSvd_SV_h->Write("", TObject::kOverwrite);
	      rooUnfoldSvdPPMC_FULLSvd_SV_h->Write("", TObject::kOverwrite);
	      rooUnfoldSvdPPMCFrac_FULLSvd_SV_h->Write("", TObject::kOverwrite);

	      delete rooUnfoldSvdPPData_FULLSvd_SV_h;
	      delete rooUnfoldSvdPPMC_FULLSvd_SV_h;
	      delete rooUnfoldSvdPPMCFrac_FULLSvd_SV_h;
	    }
	  }	   
	}
      }
      
      refptPPMC_h[cI][aI]->Write("", TObject::kOverwrite);
      refptPPMC_Frac_h[cI][aI]->Write("", TObject::kOverwrite);
      responsePPMC_h[cI][aI]->Write("", TObject::kOverwrite);
      responsePPMC_Symm_h[cI][aI]->Write("", TObject::kOverwrite);
      rooResponsePPMC_h[cI][aI]->Write("", TObject::kOverwrite);
    }
    
    for(Int_t sI = 0; sI < nSyst; ++sI){
      for(Int_t bI = 0; bI < nBayesIter; ++bI){
	jtetaPPData_FULLBayes_h[cI][sI][bI]->Write("", TObject::kOverwrite);
	jtetaPPMC_FULLBayes_h[cI][sI][bI]->Write("", TObject::kOverwrite);
	jtetaPPMC_FULLBayes_Frac_h[cI][sI][bI]->Write("", TObject::kOverwrite);
      }
    }

    treeCentSystCheckPbPbData_p[cI]->Write("", TObject::kOverwrite);
    treeCentSystCheckPPData_p[cI]->Write("", TObject::kOverwrite);
    treeCentSystCheckPbPbMC_p[cI]->Write("", TObject::kOverwrite);
    treeCentSystCheckPPMC_p[cI]->Write("", TObject::kOverwrite);
  }

  std::cout << "FIRST ROUND DONE" << std::endl;  

  writeWatch.stop();
  dumpTiming({startWatch, histWatch, pbpbDataWatch, ppDataWatch, ppMCWatch, pbpbMCWatch, countWatch, unfoldLoopWatch, plotValidWatch, copyWatch, finalEtaWatch, makeEtaHistWatch, rescaleWatch, raaWatch, writeWatch}, {"start", "hist", "pbpbData", "ppData", "ppMC", "pbpbMC", "count", "unfold", "plotValid", "copy", "finalEta", "makeEtaHist", "rescale", "raa", "write"});
  deleteWatch.start();

  std::cout << "DELETING" << std::endl;

  dumpTiming({startWatch, histWatch, pbpbDataWatch, ppDataWatch, ppMCWatch, pbpbMCWatch, countWatch, unfoldLoopWatch, plotValidWatch, copyWatch, finalEtaWatch, makeEtaHistWatch, rescaleWatch, raaWatch, writeWatch, deleteWatch}, {"start", "hist", "pbpbData", "ppData", "ppMC", "pbpbMC", "count", "unfold", "plotValid", "copy", "finalEta", "makeEtaHist", "rescale", "raa", "write", "delete"});
  deleteWatch.stop();
  closeWatch.start();

  std::cout << "DONE DELETING" << std::endl;

  outFile_p->cd();
  TNamed nCentBinsName("nCentBins", std::to_string(nCentBins));
  nCentBinsName.Write("", TObject::kOverwrite);

  std::string centBinsLowStr = "";
  std::string centBinsHiStr = "";
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    centBinsLowStr = centBinsLowStr + std::to_string(centBinsLow[cI]) + ",";
    centBinsHiStr = centBinsHiStr + std::to_string(centBinsHi[cI]) + ",";
  }
  TNamed centBinsLowName("centBinsLow", centBinsLowStr.c_str());
  TNamed centBinsHiName("centBinsHi", centBinsHiStr.c_str());
  centBinsLowName.Write("", TObject::kOverwrite);
  centBinsHiName.Write("", TObject::kOverwrite);

  TNamed nAbsEtaBinsName("nAbsEtaBins", std::to_string(nAbsEtaBins));
  nAbsEtaBinsName.Write("", TObject::kOverwrite);

  std::string absEtaBinsLowStr = "";
  std::string absEtaBinsHiStr = "";
  for(Int_t aI = 0; aI < nAbsEtaBins; ++aI){
    absEtaBinsLowStr = absEtaBinsLowStr + prettyString(absEtaBinsLow[aI], 1, true) + ",";
    absEtaBinsHiStr = absEtaBinsHiStr + prettyString(absEtaBinsHi[aI], 1, true) + ",";
  }
  TNamed absEtaBinsLowName("absEtaBinsLow", absEtaBinsLowStr.c_str());
  TNamed absEtaBinsHiName("absEtaBinsHi", absEtaBinsHiStr.c_str());
  absEtaBinsLowName.Write("", TObject::kOverwrite);
  absEtaBinsHiName.Write("", TObject::kOverwrite);

  TNamed nSystName("nSyst", std::to_string(nSyst));
  nSystName.Write("", TObject::kOverwrite);

  std::string systOutStr = "";
  for(Int_t sI = 0; sI < nSyst; ++sI){
    systOutStr = systOutStr + systStr[sI] + ",";
  }

  TNamed systName("systStr", systOutStr.c_str());
  systName.Write("", TObject::kOverwrite);


  TNamed nSystForFillName("nSystForFill", std::to_string(nSystForFill));
  nSystForFillName.Write("", TObject::kOverwrite);

  std::string systStrForFillOut = "";
  for(Int_t sI = 0; sI < nSystForFill; ++sI){
    systStrForFillOut = systStrForFillOut + systStrForFill[sI] + ",";
  }

  TNamed systForFillName("systStrForFill", systStrForFillOut.c_str());
  systForFillName.Write("", TObject::kOverwrite);

  TNamed nSystForCopyName("nSystForCopy", std::to_string(nSystForCopy));
  nSystForCopyName.Write("", TObject::kOverwrite);

  std::string systStrForCopyOut = "";
  for(Int_t sI = 0; sI < nSystForCopy; ++sI){
    systStrForCopyOut = systStrForCopyOut + systStrForCopy[sI] + ",";
  }

  TNamed systForCopyName("systStrForCopy", systStrForCopyOut.c_str());
  systForCopyName.Write("", TObject::kOverwrite);

  TNamed nBayesIterName("nBayesIter", std::to_string(nBayesIter));
  nBayesIterName.Write("", TObject::kOverwrite);

  TNamed doSvdName("doSvd", std::to_string(doSvd));
  doSvdName.Write("", TObject::kOverwrite);

  TNamed nSvdName("nSvd", std::to_string(nSvd));
  nSvdName.Write("", TObject::kOverwrite);

  outFile_p->Close();
  std::cout << "FILE CLOSED" << std::endl;
//  delete outFile_p;
  std::cout << "FILE DELETED" << std::endl;

  closeWatch.stop();
  dumpTiming({startWatch, histWatch, pbpbDataWatch, ppDataWatch, ppMCWatch, pbpbMCWatch, countWatch, unfoldLoopWatch, plotValidWatch, copyWatch, finalEtaWatch, makeEtaHistWatch, rescaleWatch, raaWatch, writeWatch, deleteWatch, closeWatch}, {"start", "hist", "pbpbData", "ppData", "ppMC", "pbpbMC", "count", "unfold", "plotValid", "copy", "finalEta", "makeEtaHist", "rescale", "raa", "write", "delete", "close"});

  delete randGen_p;
  delete date;
  std::cout << "DATE DELETED" << std::endl;

  endWatch.stop();

  Double_t total = startWatch.total() + histWatch.total() + pbpbDataWatch.total() + pbpbMCWatch.total() + ppDataWatch.total() + ppMCWatch.total() + endWatch.total();

  std::cout << "Timer start: " << startWatch.total() << " (" << startWatch.total()/total << ")" << std::endl; 
  std::cout << "Timer hist: " << histWatch.total() << " (" << histWatch.total()/total << ")" << std::endl;
  std::cout << "Timer pbpbData: " << pbpbDataWatch.total() << " (" << pbpbDataWatch.total()/total << ")" << std::endl;
  std::cout << "Timer pbpbMC: " << pbpbMCWatch.total() << " (" << pbpbMCWatch.total()/total << ")" << std::endl;
  std::cout << "Timer ppData: " << ppDataWatch.total() << " (" << ppDataWatch.total()/total << ")" << std::endl;
  std::cout << "Timer ppMC: " << ppMCWatch.total() << " (" << ppMCWatch.total()/total << ")" << std::endl;
  std::cout << "Timer count: " << countWatch.total() << " (" << countWatch.total()/total << ")" << std::endl;
  std::cout << "Timer unfoldLoop: " << unfoldLoopWatch.total() << " (" << unfoldLoopWatch.total()/total << ")" << std::endl;
  std::cout << "Timer copy: " << copyWatch.total() << " (" << copyWatch.total()/total << ")" << std::endl;
  std::cout << "Timer finalEta: " << finalEtaWatch.total() << " (" << finalEtaWatch.total()/total << ")" << std::endl;
  std::cout << "Timer makeEtaHist: " << makeEtaHistWatch.total() << " (" << makeEtaHistWatch.total()/total << ")" << std::endl;
  std::cout << "Timer rescale: " << rescaleWatch.total() << " (" << rescaleWatch.total()/total << ")" << std::endl;
  std::cout << "Timer raa: " << raaWatch.total() << " (" << raaWatch.total()/total << ")" << std::endl;
  std::cout << "Timer write: " << writeWatch.total() << " (" << writeWatch.total()/total << ")" << std::endl;
  std::cout << "Timer delete: " << deleteWatch.total() << " (" << deleteWatch.total()/total << ")" << std::endl;

  totalWatch.total();
  std::cout << "Timer total: " << totalWatch.total() << " (" << totalWatch.total()/total << ")" << std::endl;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 6){
    std::cout << "Usage: ./bin/makeFirstRAAHist_FromTree.exe <inFileNamePbPb> <inFileNamePP> <inFileNameMCPbPb> <inFileNamePPMC> <algoName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += makeFirstRAAHist_FromTree(argv[1], argv[2], argv[3], argv[4], argv[5]);
  return retVal;
}
