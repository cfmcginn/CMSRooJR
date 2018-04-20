#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLatex.h"

#include "Math/ProbFuncMathCore.h"

#include "include/inToOutFileString.h"
#include "include/mntToXRootdFileString.h"
#include "include/plotUtilities.h"
#include "include/histDefUtility.h"
#include "include/etaPhiFunc.h"
#include "include/returnRootFileContentsList.h"
#include "include/pseudoTowerGeometry.h"
#include "include/getLinBins.h"
#include "include/findBinPos.h"
#include "include/doGlobalDebug.h"
#include "include/kirchnerPalette.h"

std::vector<std::string> getUniqueStrings(std::vector< std::string > inVect)
{
  unsigned int pos = 0;
  while(inVect.size() > pos){
    bool isUnique = true;
    for(unsigned int pI = pos+1; pI < inVect.size(); ++pI){
      if(inVect.at(pos).size() == inVect.at(pI).size() && inVect.at(pos).find(inVect.at(pI)) != std::string::npos){
	inVect.erase(inVect.begin()+pI);
	isUnique = false;
	break;
      }
    }
    if(isUnique) ++pos;
  }

  return inVect;
}

double getModFactor(const double phi, const double eventPlane2, const double eventPlane3, const double par1, const double par2)
{
  double mod = 1. + 2.*(par1*cos(2.*(phi - eventPlane2))) + par2*cos(3.*(phi - eventPlane3));
  return mod;
}

int rc(const std::string inFileName, bool isMC = false, bool doFlowExamples = false)
{
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

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  bool isPP = true;
  bool hasPF = true;
  bool hasPFCs = true;
  bool hasRho = true;
  bool hasHLT = true;
  bool hasSkim = true;

  std::string mcStr = "_DATA";
  if(isMC) mcStr = "_MC";

  bool isMB = false;
  std::string mbStr = "";
  if(inFileName.find("_MB_") != std::string::npos || inFileName.find("HIMinimumBias") != std::string::npos){
    mbStr = "_MB";
    isMB = true;
  }

  TFile* inFile_p = TFile::Open(mntToXRootdFileString(fileList.at(0)).c_str(), "READ");
  TTree* hiTree_p = NULL;

  std::vector<std::string> allTrees = getUniqueStrings(returnRootFileContentsList(inFile_p, "TTree"));
  std::cout << "All trees in file: " << std::endl;
  for(unsigned int aI = 0; aI < allTrees.size(); ++aI){
    std::cout << " " << aI << "/" << allTrees.size() << ": " << allTrees.at(aI) << std::endl;
  }

  std::vector<std::string> hiTreesGlobal = getUniqueStrings(returnRootFileContentsList(inFile_p, "TTree", "hiEvtAnalyzer"));
  std::vector<std::string> pfTreesGlobal = getUniqueStrings(returnRootFileContentsList(inFile_p, "TTree", "pfcandAnalyzer"));
  unsigned int pos = 0;
  while(pfTreesGlobal.size() > pos){
    if(pfTreesGlobal.at(pos).find("pfcandAnalyzer/") == std::string::npos) pfTreesGlobal.erase(pfTreesGlobal.begin()+pos);
    else ++pos;
  }

  std::vector<std::string> flowGlobal = getUniqueStrings(returnRootFileContentsList(inFile_p, "TTree", "FlowAnalyzer"));
  std::vector<std::string> pfCsTreesGlobal = getUniqueStrings(returnRootFileContentsList(inFile_p, "TTree", "pfcandAnalyzerCs"));
  std::vector<std::string> rhoTreesGlobal = getUniqueStrings(returnRootFileContentsList(inFile_p, "TTree", "PuRho"));
  std::vector<std::string> hltTreesGlobal = getUniqueStrings(returnRootFileContentsList(inFile_p, "TTree", "hltanalysis"));
  std::vector<std::string> skimTreesGlobal = getUniqueStrings(returnRootFileContentsList(inFile_p, "TTree", "skimanalysis"));
  if(pfTreesGlobal.size() == 0) hasPF = false;
  if(pfCsTreesGlobal.size() == 0) hasPFCs = false;
  if(rhoTreesGlobal.size() == 0) hasRho = false;
  if(hltTreesGlobal.size() == 0) hasHLT = false;
  if(skimTreesGlobal.size() == 0) hasSkim = false;

  hasPF == true ? std::cout << "Has PF Tree" << std::endl : std::cout << "Doesnt have PF Tree" << std::endl;

  if(hiTreesGlobal.size() == 0) std::cout << "WARNING: no hiEvtAnalyzer in file \'" << fileList.at(0) << "\', will run as pp mode" << std::endl;
  else{
    hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
    if(hiTree_p->GetMinimum("hiBin") != hiTree_p->GetMaximum("hiBin")) isPP = false;
  }
  inFile_p->Close();
  delete inFile_p;
  
  const Int_t nEvtPlaneBins = 4;
  const std::string evtPlaneBins[nEvtPlaneBins] = {"AllPlane", "InPlane", "OutPlane", "MidPlane"};

  const Int_t nCentBinsPerma = 4;
  Int_t nCentBinsTemp = 1;
  if(!isPP) nCentBinsTemp = nCentBinsPerma;
  const Int_t nCentBins = nCentBinsTemp;
  const Int_t centBinsLow[nCentBinsPerma] = {50, 30, 10, 0};
  const Int_t centBinsHi[nCentBinsPerma] = {100, 50, 30, 10};


  const Int_t nRParam = 5;
  const Float_t rParam[nRParam] = {0.3, 0.4, 0.6, 0.8, 1.0};
  const Int_t nRCSlices[nRParam] = {100, 100, 100, 200, 200};

  const Int_t nRCBins[nRParam][nCentBinsPerma] = { {101, 101, 51, 51}, {101, 101, 51, 51}, {101, 101, 51, 51}, {101, 101, 51, 51}, {101, 101, 51, 51} };
  const Float_t rcBinsLow[nRParam][nCentBinsPerma] = { {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0} };
  const Float_t rcBinsHi[nRParam][nCentBinsPerma] = { {20, 40, 75, 75}, {20, 40, 75, 75}, {20, 50, 75, 100}, {20, 50, 100, 150}, {20, 50, 125, 200} };

  const Int_t nFlow = 1 + (Int_t)flowGlobal.size();
  TTree* flowTree_p[nFlow];
  for(Int_t i = 0; i < nFlow; ++i){
    flowTree_p[i] = NULL;
  }

  std::cout << "Is PP: " << isPP << std::endl;
  
  kirchnerPalette col;
  TDatime* date = new TDatime();
  TRandom3* randGen_p = new TRandom3(0);
  pseudoTowGeo tow;

  std::vector<double> etaTowBounds = tow.getEtaTowBounds();
  std::vector<int> etaNTowInPhi = tow.getNTowInPhi();
  std::vector< std::vector<double> > etaPhiTowSumMap;
  std::vector<float> etaTowSumMap;
  std::vector< std::vector<int> > etaTowNMapExcludeR;
  std::vector< std::vector<float> > etaTowSumMapExcludeR;

  const Int_t nEtaTowBoundsHist = etaTowBounds.size() - 1;
  Double_t etaTowBoundsHist[nEtaTowBoundsHist+1];
  for(unsigned int bI = 0; bI < etaTowBounds.size(); ++bI){
    etaTowBoundsHist[bI] = etaTowBounds.at(bI);
  }

  etaPhiTowSumMap.reserve(etaNTowInPhi.size());
  etaTowSumMap.reserve(etaNTowInPhi.size());
  etaTowSumMapExcludeR.reserve(nRParam);

  for(unsigned int i = 0; i < etaNTowInPhi.size(); ++i){
    etaPhiTowSumMap.push_back({});
    etaTowSumMap.push_back(0);

    for(Int_t j = 0; j < etaNTowInPhi.at(i); ++j){
      etaPhiTowSumMap.at(i).push_back(0.0);
    }
  }

  for(Int_t j = 0; j < nRParam; ++j){
    etaTowSumMapExcludeR.push_back({});
    etaTowNMapExcludeR.push_back({});

    for(unsigned int i = 0; i < etaNTowInPhi.size(); ++i){
      etaTowNMapExcludeR.at(j).push_back(0);
      etaTowSumMapExcludeR.at(j).push_back(0.0);
    }
  }

  Double_t phiTowBounds72[72+1];
  Double_t phiTowBounds36[36+1];
  Double_t phiTowBounds18[18+1];
  getLinBins(-TMath::Pi(), TMath::Pi(), 72, phiTowBounds72);
  getLinBins(-TMath::Pi(), TMath::Pi(), 36, phiTowBounds36);
  getLinBins(-TMath::Pi(), TMath::Pi(), 18, phiTowBounds18);

  UInt_t run_, lumi_;
  ULong64_t evt_;

  const Int_t nAbsEtaBins = 3;
  const Float_t absEtaBinsLow[nAbsEtaBins] = {0.0, 1.0, 0.0};
  const Float_t absEtaBinsHi[nAbsEtaBins] = {1.0, 2.0, 2.0};

  const Int_t nRCGen = 1;
  
  const Int_t nMaxJets = 500;
  const std::string outFileName = "output/" + inToOutFileString(inFileName, "RC" + mcStr);
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");

  //  const Int_t nRCSlices = 100;
  std::vector< std::vector<float> > areaSlices;
  for(Int_t i = 0; i < nRParam; ++i){
    areaSlices.push_back({});

   
    Double_t width = 2.*rParam[i]/(Double_t)nRCSlices[i];
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    std::cout << "Width for R=" << prettyString(rParam[i], 1, false) << ": " << width << std::endl;

    for(Int_t j = 0; j < nRCSlices[i]/2; ++j){
      Double_t x1= rParam[i] - width*j;
      Double_t x2= rParam[i] - width*(j+1);

      Double_t area = 0;
      if(3 == (int)(rParam[i]*10)){
	if(.3 - x1 < .0001) x1 = .2999;
	if(x2  < .0001) x2 = .0001;

	area = TMath::Abs(.5*TMath::Sqrt(.09 - x1*x1)*x1 + .045*TMath::ASin(10.*x1/3.) - (.5*TMath::Sqrt(.09 - x2*x2)*x2 + .045*TMath::ASin(10.*x2/3.)));
      }
      else if(4 == (int)(rParam[i]*10)){
	if(.4 - x1 < .0001) x1 = .3999;
	if(x2  < .0001) x2 = .0001;

	area = TMath::Abs(.1*TMath::Sqrt(4. - 25.*x1*x1)*x1 + .08*TMath::ASin(2.5*x1) - (.1*TMath::Sqrt(4. - 25.*x2*x2)*x2 + .08*TMath::ASin(2.5*x2)));
      }
      else if(6 == (int)(rParam[i]*10)){
	if(.6 - x1 < .0001) x1 = .5999;
	if(x2  < .0001) x2 = .0001;

	area = TMath::Abs(.1*TMath::Sqrt(9. - 25.*x1*x1)*x1 + .18*TMath::ASin(5.*x1/3.) - (.1*TMath::Sqrt(9. - 25.*x2*x2)*x2 + .18*TMath::ASin(5.*x2/3.)));
      }
      else if(8 == (int)(rParam[i]*10)){
	if(.8 - x1 < .0001) x1 = .7999;
	if(x2  < .0001) x2 = .0001;

	area = TMath::Abs(.1*TMath::Sqrt(16. - 25.*x1*x1)*x1 + .32*TMath::ASin(1.25*x1) - (.1*TMath::Sqrt(16. - 25.*x2*x2)*x2 + .32*TMath::ASin(1.25*x2)));
      }
      else if(10 == (int)(rParam[i]*10)){
	if(1. - x1 < .0001) x1 = .9999;
	if(x2  < .0001) x2 = .0001;

	area = TMath::Abs(.5*TMath::Sqrt(1. - x1*x1)*x1 + .5*TMath::ASin(x1) - (.5*TMath::Sqrt(1. - x2*x2)*x2 + .5*TMath::ASin(x2)));
      }
      area *= 2.;

      areaSlices.at(i).push_back(area);
    }
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    for(Int_t j = 0; j < nRCSlices[i]/2; ++j){
      areaSlices.at(i).push_back(areaSlices.at(i).at(nRCSlices[i]/2 - j - 1));
    }
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::cout << "AREA CHECK" << std::endl;
  for(Int_t i = 0; i < nRParam; ++i){
    Double_t totArea = 0.;
    for(Int_t j = 0; j < nRCSlices[i]; ++j){
      //      std::cout << " j " << areaSlices.at(i).at(j) << std::endl;
      totArea += areaSlices.at(i).at(j);
    }

    std::cout << " R=" << prettyString(rParam[i],1,false) << "(True/False): " << TMath::Pi()*rParam[i]*rParam[i] << "/" << totArea << "=" << TMath::Pi()*rParam[i]*rParam[i]/totArea << std::endl;
  }


  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  TH2F* rcSum_VHiBin_h[nAbsEtaBins][nRParam];
  TH2F* pu_VHiBin_h[nAbsEtaBins][nRParam];
  TH2F* puFlow_VHiBin_h[nAbsEtaBins][nRParam];
  TH2F* rcSumMinPU_VHiBin_h[nAbsEtaBins][nRParam];
  TH2F* rcSumMinPUFlow_VHiBin_h[nAbsEtaBins][nRParam];

  TH2F* rcSum_VHiBin_GenExclude_h[nAbsEtaBins][nRParam];
  TH2F* pu_VHiBin_GenExclude_h[nAbsEtaBins][nRParam];
  TH2F* puFlow_VHiBin_GenExclude_h[nAbsEtaBins][nRParam];
  TH2F* rcSumMinPU_VHiBin_GenExclude_h[nAbsEtaBins][nRParam];
  TH2F* rcSumMinPUFlow_VHiBin_GenExclude_h[nAbsEtaBins][nRParam];
  
  TH1F* rcSum_VHiBin_Mean_h[nAbsEtaBins][nRParam];
  TH1F* pu_VHiBin_Mean_h[nAbsEtaBins][nRParam];
  TH1F* puFlow_VHiBin_Mean_h[nAbsEtaBins][nRParam];
  TH1F* rcSumMinPU_VHiBin_Mean_h[nAbsEtaBins][nRParam];
  TH1F* rcSumMinPUFlow_VHiBin_Mean_h[nAbsEtaBins][nRParam];

  TH1F* rcSum_VHiBin_Mean_GenExclude_h[nAbsEtaBins][nRParam];
  TH1F* pu_VHiBin_Mean_GenExclude_h[nAbsEtaBins][nRParam];
  TH1F* puFlow_VHiBin_Mean_GenExclude_h[nAbsEtaBins][nRParam];
  TH1F* rcSumMinPU_VHiBin_Mean_GenExclude_h[nAbsEtaBins][nRParam];
  TH1F* rcSumMinPUFlow_VHiBin_Mean_GenExclude_h[nAbsEtaBins][nRParam];

  TH1F* rcSum_VHiBin_Sigma_h[nAbsEtaBins][nRParam];
  TH1F* pu_VHiBin_Sigma_h[nAbsEtaBins][nRParam];
  TH1F* puFlow_VHiBin_Sigma_h[nAbsEtaBins][nRParam];
  TH1F* rcSumMinPU_VHiBin_Sigma_h[nAbsEtaBins][nRParam];
  TH1F* rcSumMinPUFlow_VHiBin_Sigma_h[nAbsEtaBins][nRParam];

  TH1F* rcSum_VHiBin_Sigma_GenExclude_h[nAbsEtaBins][nRParam];
  TH1F* pu_VHiBin_Sigma_GenExclude_h[nAbsEtaBins][nRParam];
  TH1F* puFlow_VHiBin_Sigma_GenExclude_h[nAbsEtaBins][nRParam];
  TH1F* rcSumMinPU_VHiBin_Sigma_GenExclude_h[nAbsEtaBins][nRParam];
  TH1F* rcSumMinPUFlow_VHiBin_Sigma_GenExclude_h[nAbsEtaBins][nRParam];

  TH1F* rcSum_VHiBin_MeanRed_h[nAbsEtaBins][nRParam];
  TH1F* pu_VHiBin_MeanRed_h[nAbsEtaBins][nRParam];
  TH1F* puFlow_VHiBin_MeanRed_h[nAbsEtaBins][nRParam];
  TH1F* rcSumMinPU_VHiBin_MeanRed_h[nAbsEtaBins][nRParam];
  TH1F* rcSumMinPUFlow_VHiBin_MeanRed_h[nAbsEtaBins][nRParam];

  TH1F* rcSum_VHiBin_MeanRed_GenExclude_h[nAbsEtaBins][nRParam];
  TH1F* pu_VHiBin_MeanRed_GenExclude_h[nAbsEtaBins][nRParam];
  TH1F* puFlow_VHiBin_MeanRed_GenExclude_h[nAbsEtaBins][nRParam];
  TH1F* rcSumMinPU_VHiBin_MeanRed_GenExclude_h[nAbsEtaBins][nRParam];
  TH1F* rcSumMinPUFlow_VHiBin_MeanRed_GenExclude_h[nAbsEtaBins][nRParam];

  TH1F* rcSum_VHiBin_SigmaRed_h[nAbsEtaBins][nRParam];
  TH1F* pu_VHiBin_SigmaRed_h[nAbsEtaBins][nRParam];
  TH1F* puFlow_VHiBin_SigmaRed_h[nAbsEtaBins][nRParam];
  TH1F* rcSumMinPU_VHiBin_SigmaRed_h[nAbsEtaBins][nRParam];
  TH1F* rcSumMinPUFlow_VHiBin_SigmaRed_h[nAbsEtaBins][nRParam];

  TH1F* rcSum_VHiBin_SigmaRed_GenExclude_h[nAbsEtaBins][nRParam];
  TH1F* pu_VHiBin_SigmaRed_GenExclude_h[nAbsEtaBins][nRParam];
  TH1F* puFlow_VHiBin_SigmaRed_GenExclude_h[nAbsEtaBins][nRParam];
  TH1F* rcSumMinPU_VHiBin_SigmaRed_GenExclude_h[nAbsEtaBins][nRParam];
  TH1F* rcSumMinPUFlow_VHiBin_SigmaRed_GenExclude_h[nAbsEtaBins][nRParam];

  TH1F* rcSum_VHiBin_SigmaOverMean_h[nAbsEtaBins][nRParam];
  TH1F* pu_VHiBin_SigmaOverMean_h[nAbsEtaBins][nRParam];
  TH1F* puFlow_VHiBin_SigmaOverMean_h[nAbsEtaBins][nRParam];

  
  TH1F* rcSum_VHiBin_Points_h[nAbsEtaBins][nRParam][200];
  TH1F* pu_VHiBin_Points_h[nAbsEtaBins][nRParam][200];
  TH1F* puFlow_VHiBin_Points_h[nAbsEtaBins][nRParam][200];
  TH1F* rcSumMinPU_VHiBin_Points_h[nAbsEtaBins][nRParam][200];
  TH1F* rcSumMinPUFlow_VHiBin_Points_h[nAbsEtaBins][nRParam][200];

  TH1F* rcSum_VHiBin_Points_GenExclude_h[nAbsEtaBins][nRParam][200];
  TH1F* pu_VHiBin_Points_GenExclude_h[nAbsEtaBins][nRParam][200];
  TH1F* puFlow_VHiBin_Points_GenExclude_h[nAbsEtaBins][nRParam][200];
  TH1F* rcSumMinPU_VHiBin_Points_GenExclude_h[nAbsEtaBins][nRParam][200];
  TH1F* rcSumMinPUFlow_VHiBin_Points_GenExclude_h[nAbsEtaBins][nRParam][200];
  
  const Int_t nRedCentBins = 40;
  const Double_t redCentBinsMeanR0p3[nRedCentBins] = {58.9507,50.5478,45.352,40.6612,34.8815,30.9559,29.0054,26.1559,23.1831,20.6756,17.4449,15.9349,14.4703,12.9198,11.4067,10.2323,9.12026,7.48967,6.59092,5.73161,5.19153,4.36589,3.66128,2.86646,2.70771,2.20059,1.9296,1.38559,1.12771,0.918846,0.739115,0.648358,0.536478,0.256043,0.174638,0.164799,0.176926,0.0966774,0.0395185};
  const Double_t redCentBinsSigmaR0p3[nRedCentBins] = {13.7254,12.1122,11.4963,12.0154,10.7128,9.17153,9.14818,8.62232,7.46348,7.14227,6.88655,6.76017,6.39801,5.8351,5.31991,4.79249,5.07358,4.27094,3.54909,3.60555,3.01206,3.1044,2.65711,2.13227,2.10842,1.87699,1.86255,1.37924,1.74245,1.32724,1.02438,0.91799,0.906352,0.564485,0.41157,0.385843,0.523353,0.346894,0.154523};

  TH1F* rcSum_VHiBin_PointsRed_h[nAbsEtaBins][nRParam][nRedCentBins];
  TH1F* pu_VHiBin_PointsRed_h[nAbsEtaBins][nRParam][nRedCentBins];
  TH1F* puFlow_VHiBin_PointsRed_h[nAbsEtaBins][nRParam][nRedCentBins];
  TH1F* rcSumMinPU_VHiBin_PointsRed_h[nAbsEtaBins][nRParam][nRedCentBins];
  TH1F* rcSumMinPUFlow_VHiBin_PointsRed_h[nAbsEtaBins][nRParam][nRedCentBins];

  TH1F* rcSum_VHiBin_PointsRed_GenExclude_h[nAbsEtaBins][nRParam][nRedCentBins];
  TH1F* pu_VHiBin_PointsRed_GenExclude_h[nAbsEtaBins][nRParam][nRedCentBins];
  TH1F* puFlow_VHiBin_PointsRed_GenExclude_h[nAbsEtaBins][nRParam][nRedCentBins];
  TH1F* rcSumMinPU_VHiBin_PointsRed_GenExclude_h[nAbsEtaBins][nRParam][nRedCentBins];
  TH1F* rcSumMinPUFlow_VHiBin_PointsRed_GenExclude_h[nAbsEtaBins][nRParam][nRedCentBins];
 
  TH1F* rcSum_h[nCentBins][nEvtPlaneBins][nAbsEtaBins][nRParam];
  TH1F* rcSum_GenExclude_h[nCentBins][nEvtPlaneBins][nAbsEtaBins][nRParam];

  TH1F* rcSumMinPU_h[nCentBins][nEvtPlaneBins][nAbsEtaBins][nRParam];
  TH1F* rcSumMinPU_GenExclude_h[nCentBins][nEvtPlaneBins][nAbsEtaBins][nRParam];

  TH1F* rcSumMinPUFlow_h[nCentBins][nEvtPlaneBins][nAbsEtaBins][nRParam][nFlow];
  TH1F* rcSumMinPUFlow_GenExclude_h[nCentBins][nEvtPlaneBins][nAbsEtaBins][nRParam][nFlow];

  const Int_t nPFID = 7;
  TH1F* pfCandMap_Mean_h[nCentBins][nPFID+1];
  TH1F* pfCandMap_Sigma_h[nCentBins][nPFID+1];
  TH1F* pfCandMap_Points_h[nCentBins][nPFID+1][nEtaTowBoundsHist];

  TH1F* pfCsCandMap_Mean_h[nCentBins][nPFID+1];
  TH1F* pfCsCandMap_Sigma_h[nCentBins][nPFID+1];
  TH1F* pfCsCandMap_Points_h[nCentBins][nPFID+1][nEtaTowBoundsHist];

  for(Int_t i = 0; i < nAbsEtaBins; ++i){
    std::string absEtaStr = "AbsEta" + prettyString(absEtaBinsLow[i], 1, true) + "to" + prettyString(absEtaBinsHi[i], 1, true);
    std::string absEtaStr2 = prettyString(absEtaBinsLow[i], 1, false) + " < |#eta| < " + prettyString(absEtaBinsHi[i], 1, false);
    
    for(Int_t j = 0; j < nRParam; ++j){
      std::string rStr = "R" + prettyString(rParam[j], 1, true);
      std::string rStr2 = "R=" + prettyString(rParam[j], 1, false);

      Float_t rHi = 200;
      if(int(rParam[j]*10)==6) rHi = 400;
      else if(int(rParam[j]*10)==8) rHi = 600;
      else if(int(rParam[j]*10)==10) rHi = 800;

      Float_t rMinHi = 400;
      
      rcSum_VHiBin_h[i][j] = new TH2F(("rcSum_VHiBin_" + absEtaStr + "_" + rStr + "_h").c_str(), (";hiBin;RC #Sigmap_{T} (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5, 200, 0, rHi);
      pu_VHiBin_h[i][j] = new TH2F(("pu_VHiBin_" + absEtaStr + "_" + rStr + "_h").c_str(), (";hiBin;RC #rhoA (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5, 200, 0, rHi);
      puFlow_VHiBin_h[i][j] = new TH2F(("puFlow_VHiBin_" + absEtaStr + "_" + rStr + "_h").c_str(), (";hiBin;RC #rhoA (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5, 200, 0, rHi);
      rcSumMinPU_VHiBin_h[i][j] = new TH2F(("rcSumMinPU_VHiBin_" + absEtaStr + "_" + rStr + "_h").c_str(), (";hiBin;RC #Sigmap_{T} - #rhoA (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5, 400, -rMinHi, rMinHi);
      rcSumMinPUFlow_VHiBin_h[i][j] = new TH2F(("rcSumMinPUFlow_VHiBin_" + absEtaStr + "_" + rStr + "_h").c_str(), (";hiBin;RC #Sigmap_{T} - #rhoA (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5, 400, -rMinHi, rMinHi);

      if(isMC && !isMB){
	rcSum_VHiBin_GenExclude_h[i][j] = new TH2F(("rcSum_VHiBin_" + absEtaStr + "_" + rStr + "_GenExclude_h").c_str(), (";hiBin;RC #Sigmap_{T} (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5, 200, 0, rHi);
	pu_VHiBin_GenExclude_h[i][j] = new TH2F(("pu_VHiBin_" + absEtaStr + "_" + rStr + "_GenExclude_h").c_str(), (";hiBin;RC #rhoA (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5, 200, 0, rHi);
	puFlow_VHiBin_GenExclude_h[i][j] = new TH2F(("puFlow_VHiBin_" + absEtaStr + "_" + rStr + "_GenExclude_h").c_str(), (";hiBin;RC #rhoA (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5, 200, 0, rHi);
	rcSumMinPU_VHiBin_GenExclude_h[i][j] = new TH2F(("rcSumMinPU_VHiBin_" + absEtaStr + "_" + rStr + "_GenExclude_h").c_str(), (";hiBin;RC #Sigmap_{T} - #rhoA (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5, 400, -rMinHi, rMinHi);
	rcSumMinPUFlow_VHiBin_GenExclude_h[i][j] = new TH2F(("rcSumMinPUFlow_VHiBin_" + absEtaStr + "_" + rStr + "_GenExclude_h").c_str(), (";hiBin;RC #Sigmap_{T} - #rhoA (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5, 400, -rMinHi, rMinHi);
      }

      rcSum_VHiBin_Mean_h[i][j] = new TH1F(("rcSum_VHiBin_" + absEtaStr + "_" + rStr + "_Mean_h").c_str(), (";hiBin;#LTRC #Sigmap_{T}#GT (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5);
      pu_VHiBin_Mean_h[i][j] = new TH1F(("pu_VHiBin_" + absEtaStr + "_" + rStr + "_Mean_h").c_str(), (";hiBin;#LTRC #rhoA#GT (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5);
      puFlow_VHiBin_Mean_h[i][j] = new TH1F(("puFlow_VHiBin_" + absEtaStr + "_" + rStr + "_Mean_h").c_str(), (";hiBin;#LTRC #rhoA#GT (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5);
      rcSumMinPU_VHiBin_Mean_h[i][j] = new TH1F(("rcSumMinPU_VHiBin_" + absEtaStr + "_" + rStr + "_Mean_h").c_str(), (";hiBin;#LTRC #rhoA#GT (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5);
      rcSumMinPUFlow_VHiBin_Mean_h[i][j] = new TH1F(("rcSumMinPUFlow_VHiBin_" + absEtaStr + "_" + rStr + "_Mean_h").c_str(), (";hiBin;#LTRC #rhoA#GT (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5);

      if(isMC && !isMB){
	rcSum_VHiBin_Mean_GenExclude_h[i][j] = new TH1F(("rcSum_VHiBin_" + absEtaStr + "_" + rStr + "_Mean_GenExclude_h").c_str(), (";hiBin;#LTRC #Sigmap_{T}#GT (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5);
	pu_VHiBin_Mean_GenExclude_h[i][j] = new TH1F(("pu_VHiBin_" + absEtaStr + "_" + rStr + "_Mean_GenExclude_h").c_str(), (";hiBin;#LTRC #rhoA#GT (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5);
	puFlow_VHiBin_Mean_GenExclude_h[i][j] = new TH1F(("puFlow_VHiBin_" + absEtaStr + "_" + rStr + "_Mean_GenExclude_h").c_str(), (";hiBin;#LTRC #rhoA#GT (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5);
	rcSumMinPU_VHiBin_Mean_GenExclude_h[i][j] = new TH1F(("rcSumMinPU_VHiBin_" + absEtaStr + "_" + rStr + "_Mean_GenExclude_h").c_str(), (";hiBin;#LTRC #rhoA#GT (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5);
	rcSumMinPUFlow_VHiBin_Mean_GenExclude_h[i][j] = new TH1F(("rcSumMinPUFlow_VHiBin_" + absEtaStr + "_" + rStr + "_Mean_GenExclude_h").c_str(), (";hiBin;#LTRC #rhoA#GT (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5);
      }


      rcSum_VHiBin_MeanRed_h[i][j] = new TH1F(("rcSum_VHiBin_" + absEtaStr + "_" + rStr + "_MeanRed_h").c_str(), (";hiBin;#LTRC #Sigmap_{T}#GT (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), nRedCentBins, -0.5, 199.5);
      pu_VHiBin_MeanRed_h[i][j] = new TH1F(("pu_VHiBin_" + absEtaStr + "_" + rStr + "_MeanRed_h").c_str(), (";hiBin;#LTRC #rhoA#GT (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), nRedCentBins, -0.5, 199.5);
      puFlow_VHiBin_MeanRed_h[i][j] = new TH1F(("puFlow_VHiBin_" + absEtaStr + "_" + rStr + "_MeanRed_h").c_str(), (";hiBin;#LTRC #rhoA#GT (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), nRedCentBins, -0.5, 199.5);
      rcSumMinPU_VHiBin_MeanRed_h[i][j] = new TH1F(("rcSumMinPU_VHiBin_" + absEtaStr + "_" + rStr + "_MeanRed_h").c_str(), (";hiBin;#LTRC #rhoA#GT (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), nRedCentBins, -0.5, 199.5);
      rcSumMinPUFlow_VHiBin_MeanRed_h[i][j] = new TH1F(("rcSumMinPUFlow_VHiBin_" + absEtaStr + "_" + rStr + "_MeanRed_h").c_str(), (";hiBin;#LTRC #rhoA#GT (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), nRedCentBins, -0.5, 199.5);

      if(isMC && !isMB){
	rcSum_VHiBin_MeanRed_GenExclude_h[i][j] = new TH1F(("rcSum_VHiBin_" + absEtaStr + "_" + rStr + "_MeanRed_GenExclude_h").c_str(), (";hiBin;#LTRC #Sigmap_{T}#GT (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), nRedCentBins, -0.5, 199.5);
	pu_VHiBin_MeanRed_GenExclude_h[i][j] = new TH1F(("pu_VHiBin_" + absEtaStr + "_" + rStr + "_MeanRed_GenExclude_h").c_str(), (";hiBin;#LTRC #rhoA#GT (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), nRedCentBins, -0.5, 199.5);
	puFlow_VHiBin_MeanRed_GenExclude_h[i][j] = new TH1F(("puFlow_VHiBin_" + absEtaStr + "_" + rStr + "_MeanRed_GenExclude_h").c_str(), (";hiBin;#LTRC #rhoA#GT (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), nRedCentBins, -0.5, 199.5);
	rcSumMinPU_VHiBin_MeanRed_GenExclude_h[i][j] = new TH1F(("rcSumMinPU_VHiBin_" + absEtaStr + "_" + rStr + "_MeanRed_GenExclude_h").c_str(), (";hiBin;#LTRC #rhoA#GT (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), nRedCentBins, -0.5, 199.5);
	rcSumMinPUFlow_VHiBin_MeanRed_GenExclude_h[i][j] = new TH1F(("rcSumMinPUFlow_VHiBin_" + absEtaStr + "_" + rStr + "_MeanRed_GenExclude_h").c_str(), (";hiBin;#LTRC #rhoA#GT (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), nRedCentBins, -0.5, 199.5);
      }

      rcSum_VHiBin_Sigma_h[i][j] = new TH1F(("rcSum_VHiBin_" + absEtaStr + "_" + rStr + "_Sigma_h").c_str(), (";hiBin;#sigma(RC #Sigmap_{T}) (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5);
      pu_VHiBin_Sigma_h[i][j] = new TH1F(("pu_VHiBin_" + absEtaStr + "_" + rStr + "_Sigma_h").c_str(), (";hiBin;#sigma(RC #rhoA) (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5);
      puFlow_VHiBin_Sigma_h[i][j] = new TH1F(("puFlow_VHiBin_" + absEtaStr + "_" + rStr + "_Sigma_h").c_str(), (";hiBin;#sigma(RC #rhoA) (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5);
      rcSumMinPU_VHiBin_Sigma_h[i][j] = new TH1F(("rcSumMinPU_VHiBin_" + absEtaStr + "_" + rStr + "_Sigma_h").c_str(), (";hiBin;#sigma(RC #rhoA) (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5);
      rcSumMinPUFlow_VHiBin_Sigma_h[i][j] = new TH1F(("rcSumMinPUFlow_VHiBin_" + absEtaStr + "_" + rStr + "_Sigma_h").c_str(), (";hiBin;#sigma(RC #rhoA) (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5);

      if(isMC && !isMB){
	rcSum_VHiBin_Sigma_GenExclude_h[i][j] = new TH1F(("rcSum_VHiBin_" + absEtaStr + "_" + rStr + "_Sigma_GenExclude_h").c_str(), (";hiBin;#sigma(RC #Sigmap_{T}) (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5);
	pu_VHiBin_Sigma_GenExclude_h[i][j] = new TH1F(("pu_VHiBin_" + absEtaStr + "_" + rStr + "_Sigma_GenExclude_h").c_str(), (";hiBin;#sigma(RC #rhoA) (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5);
	puFlow_VHiBin_Sigma_GenExclude_h[i][j] = new TH1F(("puFlow_VHiBin_" + absEtaStr + "_" + rStr + "_Sigma_GenExclude_h").c_str(), (";hiBin;#sigma(RC #rhoA) (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5);
	rcSumMinPU_VHiBin_Sigma_GenExclude_h[i][j] = new TH1F(("rcSumMinPU_VHiBin_" + absEtaStr + "_" + rStr + "_Sigma_GenExclude_h").c_str(), (";hiBin;#sigma(RC #rhoA) (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5);
	rcSumMinPUFlow_VHiBin_Sigma_GenExclude_h[i][j] = new TH1F(("rcSumMinPUFlow_VHiBin_" + absEtaStr + "_" + rStr + "_Sigma_GenExclude_h").c_str(), (";hiBin;#sigma(RC #rhoA) (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5);
      }


      rcSum_VHiBin_SigmaRed_h[i][j] = new TH1F(("rcSum_VHiBin_" + absEtaStr + "_" + rStr + "_SigmaRed_h").c_str(), (";hiBin;#sigma(RC #SigmaRedp_{T}) (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), nRedCentBins, -0.5, 199.5);
      pu_VHiBin_SigmaRed_h[i][j] = new TH1F(("pu_VHiBin_" + absEtaStr + "_" + rStr + "_SigmaRed_h").c_str(), (";hiBin;#sigma(RC #rhoA) (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), nRedCentBins, -0.5, 199.5);
      puFlow_VHiBin_SigmaRed_h[i][j] = new TH1F(("puFlow_VHiBin_" + absEtaStr + "_" + rStr + "_SigmaRed_h").c_str(), (";hiBin;#sigma(RC #rhoA) (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), nRedCentBins, -0.5, 199.5);
      rcSumMinPU_VHiBin_SigmaRed_h[i][j] = new TH1F(("rcSumMinPU_VHiBin_" + absEtaStr + "_" + rStr + "_SigmaRed_h").c_str(), (";hiBin;#sigma(RC #rhoA) (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), nRedCentBins, -0.5, 199.5);
      rcSumMinPUFlow_VHiBin_SigmaRed_h[i][j] = new TH1F(("rcSumMinPUFlow_VHiBin_" + absEtaStr + "_" + rStr + "_SigmaRed_h").c_str(), (";hiBin;#sigma(RC #rhoA) (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), nRedCentBins, -0.5, 199.5);

      if(isMC && !isMB){
	rcSum_VHiBin_SigmaRed_GenExclude_h[i][j] = new TH1F(("rcSum_VHiBin_" + absEtaStr + "_" + rStr + "_SigmaRed_GenExclude_h").c_str(), (";hiBin;#sigma(RC #SigmaRedp_{T}) (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), nRedCentBins, -0.5, 199.5);
	pu_VHiBin_SigmaRed_GenExclude_h[i][j] = new TH1F(("pu_VHiBin_" + absEtaStr + "_" + rStr + "_SigmaRed_GenExclude_h").c_str(), (";hiBin;#sigma(RC #rhoA) (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), nRedCentBins, -0.5, 199.5);
	puFlow_VHiBin_SigmaRed_GenExclude_h[i][j] = new TH1F(("puFlow_VHiBin_" + absEtaStr + "_" + rStr + "_SigmaRed_GenExclude_h").c_str(), (";hiBin;#sigma(RC #rhoA) (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), nRedCentBins, -0.5, 199.5);
	rcSumMinPU_VHiBin_SigmaRed_GenExclude_h[i][j] = new TH1F(("rcSumMinPU_VHiBin_" + absEtaStr + "_" + rStr + "_SigmaRed_GenExclude_h").c_str(), (";hiBin;#sigma(RC #rhoA) (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), nRedCentBins, -0.5, 199.5);
	rcSumMinPUFlow_VHiBin_SigmaRed_GenExclude_h[i][j] = new TH1F(("rcSumMinPUFlow_VHiBin_" + absEtaStr + "_" + rStr + "_SigmaRed_GenExclude_h").c_str(), (";hiBin;#sigma(RC #rhoA) (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), nRedCentBins, -0.5, 199.5);
      }

      rcSum_VHiBin_SigmaOverMean_h[i][j] = new TH1F(("rcSum_VHiBin_" + absEtaStr + "_" + rStr + "_SigmaOverMean_h").c_str(), (";hiBin;#sigma(RC #Sigmap_{T})/#LTRC #rhoA#GT (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5);
      pu_VHiBin_SigmaOverMean_h[i][j] = new TH1F(("pu_VHiBin_" + absEtaStr + "_" + rStr + "_SigmaOverMean_h").c_str(), (";hiBin;#sigma(RC #rhoA)/#LTRC #rhoA#GT (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5);
      puFlow_VHiBin_SigmaOverMean_h[i][j] = new TH1F(("puFlow_VHiBin_" + absEtaStr + "_" + rStr + "_SigmaOverMean_h").c_str(), (";hiBin;#sigma(RC #rhoA)/#LTRC #rhoA#GT (" + rStr2 + ", " + absEtaStr2 + ")").c_str(), 200, -0.5, 199.5);

      centerTitles({rcSum_VHiBin_h[i][j], pu_VHiBin_h[i][j], puFlow_VHiBin_h[i][j], rcSumMinPU_VHiBin_h[i][j], rcSumMinPUFlow_VHiBin_h[i][j], rcSum_VHiBin_Mean_h[i][j], pu_VHiBin_Mean_h[i][j], puFlow_VHiBin_Mean_h[i][j], rcSumMinPU_VHiBin_Mean_h[i][j], rcSumMinPUFlow_VHiBin_Mean_h[i][j], rcSum_VHiBin_Sigma_h[i][j], pu_VHiBin_Sigma_h[i][j], puFlow_VHiBin_Sigma_h[i][j], rcSumMinPU_VHiBin_Sigma_h[i][j], rcSumMinPUFlow_VHiBin_Sigma_h[i][j], rcSum_VHiBin_MeanRed_h[i][j], pu_VHiBin_MeanRed_h[i][j], puFlow_VHiBin_MeanRed_h[i][j], rcSumMinPU_VHiBin_MeanRed_h[i][j], rcSumMinPUFlow_VHiBin_MeanRed_h[i][j], rcSum_VHiBin_SigmaRed_h[i][j], pu_VHiBin_SigmaRed_h[i][j], puFlow_VHiBin_SigmaRed_h[i][j], rcSumMinPU_VHiBin_SigmaRed_h[i][j], rcSumMinPUFlow_VHiBin_SigmaRed_h[i][j], rcSum_VHiBin_SigmaOverMean_h[i][j], pu_VHiBin_SigmaOverMean_h[i][j], puFlow_VHiBin_SigmaOverMean_h[i][j]});

      if(isMC && !isMB) centerTitles({rcSum_VHiBin_GenExclude_h[i][j], pu_VHiBin_GenExclude_h[i][j], puFlow_VHiBin_GenExclude_h[i][j], rcSumMinPU_VHiBin_GenExclude_h[i][j], rcSumMinPUFlow_VHiBin_GenExclude_h[i][j], rcSum_VHiBin_Mean_GenExclude_h[i][j], pu_VHiBin_Mean_GenExclude_h[i][j], puFlow_VHiBin_Mean_GenExclude_h[i][j], rcSumMinPU_VHiBin_Mean_GenExclude_h[i][j], rcSumMinPUFlow_VHiBin_Mean_GenExclude_h[i][j], rcSum_VHiBin_Sigma_GenExclude_h[i][j], pu_VHiBin_Sigma_GenExclude_h[i][j], puFlow_VHiBin_Sigma_GenExclude_h[i][j], rcSumMinPU_VHiBin_Sigma_GenExclude_h[i][j], rcSumMinPUFlow_VHiBin_Sigma_GenExclude_h[i][j], rcSum_VHiBin_MeanRed_GenExclude_h[i][j], pu_VHiBin_MeanRed_GenExclude_h[i][j], puFlow_VHiBin_MeanRed_GenExclude_h[i][j], rcSumMinPU_VHiBin_MeanRed_GenExclude_h[i][j], rcSumMinPUFlow_VHiBin_MeanRed_GenExclude_h[i][j], rcSum_VHiBin_SigmaRed_GenExclude_h[i][j], pu_VHiBin_SigmaRed_GenExclude_h[i][j], puFlow_VHiBin_SigmaRed_GenExclude_h[i][j], rcSumMinPU_VHiBin_SigmaRed_GenExclude_h[i][j], rcSumMinPUFlow_VHiBin_SigmaRed_GenExclude_h[i][j]});

           
      for(Int_t pI = 0; pI < 200; ++pI){
	rcSum_VHiBin_Points_h[i][j][pI] = new TH1F(("rcSum_VHiBin_" + absEtaStr + "_" + rStr + "_hiBin" + std::to_string(pI) + "_Points_h").c_str(), (";RC #Sigmap_{T} (" + rStr2 + ", " + absEtaStr2 + ");Counts").c_str(), 200, 0, 1000);
	pu_VHiBin_Points_h[i][j][pI] = new TH1F(("pu_VHiBin_" + absEtaStr + "_" + rStr + "_hiBin" + std::to_string(pI) + "_Points_h").c_str(), (";RC #rhoA (" + rStr2 + ", " + absEtaStr2 + ");Counts").c_str(), 200, 0, 1000);
	puFlow_VHiBin_Points_h[i][j][pI] = new TH1F(("puFlow_VHiBin_" + absEtaStr + "_" + rStr + "_hiBin" + std::to_string(pI) + "_Points_h").c_str(), (";RC #rhoA (" + rStr2 + ", " + absEtaStr2 + ");Counts").c_str(), 200, 0, 1000);
	rcSumMinPU_VHiBin_Points_h[i][j][pI] = new TH1F(("rcSumMinPU_VHiBin_" + absEtaStr + "_" + rStr + "_hiBin" + std::to_string(pI) + "_Points_h").c_str(), (";RC #rhoA (" + rStr2 + ", " + absEtaStr2 + ");Counts").c_str(), 200, -1000, 1000);
	rcSumMinPUFlow_VHiBin_Points_h[i][j][pI] = new TH1F(("rcSumMinPUFlow_VHiBin_" + absEtaStr + "_" + rStr + "_hiBin" + std::to_string(pI) + "_Points_h").c_str(), (";RC #rhoA (" + rStr2 + ", " + absEtaStr2 + ");Counts").c_str(), 200, -1000, 1000);

	if(isMC && !isMB){
	rcSum_VHiBin_Points_GenExclude_h[i][j][pI] = new TH1F(("rcSum_VHiBin_" + absEtaStr + "_" + rStr + "_hiBin" + std::to_string(pI) + "_Points_GenExclude_h").c_str(), (";RC #Sigmap_{T} (" + rStr2 + ", " + absEtaStr2 + ");Counts").c_str(), 200, 0, 1000);
	pu_VHiBin_Points_GenExclude_h[i][j][pI] = new TH1F(("pu_VHiBin_" + absEtaStr + "_" + rStr + "_hiBin" + std::to_string(pI) + "_Points_GenExclude_h").c_str(), (";RC #rhoA (" + rStr2 + ", " + absEtaStr2 + ");Counts").c_str(), 200, 0, 1000);
	puFlow_VHiBin_Points_GenExclude_h[i][j][pI] = new TH1F(("puFlow_VHiBin_" + absEtaStr + "_" + rStr + "_hiBin" + std::to_string(pI) + "_Points_GenExclude_h").c_str(), (";RC #rhoA (" + rStr2 + ", " + absEtaStr2 + ");Counts").c_str(), 200, 0, 1000);
	rcSumMinPU_VHiBin_Points_GenExclude_h[i][j][pI] = new TH1F(("rcSumMinPU_VHiBin_" + absEtaStr + "_" + rStr + "_hiBin" + std::to_string(pI) + "_Points_GenExclude_h").c_str(), (";RC #rhoA (" + rStr2 + ", " + absEtaStr2 + ");Counts").c_str(), 200, -1000, 1000);
	rcSumMinPUFlow_VHiBin_Points_GenExclude_h[i][j][pI] = new TH1F(("rcSumMinPUFlow_VHiBin_" + absEtaStr + "_" + rStr + "_hiBin" + std::to_string(pI) + "_Points_GenExclude_h").c_str(), (";RC #rhoA (" + rStr2 + ", " + absEtaStr2 + ");Counts").c_str(), 200, -1000, 1000);
	}

	centerTitles({rcSum_VHiBin_Points_h[i][j][pI], pu_VHiBin_Points_h[i][j][pI], puFlow_VHiBin_Points_h[i][j][pI], rcSumMinPU_VHiBin_Points_h[i][j][pI], rcSumMinPUFlow_VHiBin_Points_h[i][j][pI]});

	if(isMC && !isMB) centerTitles({rcSum_VHiBin_Points_GenExclude_h[i][j][pI], pu_VHiBin_Points_GenExclude_h[i][j][pI], puFlow_VHiBin_Points_GenExclude_h[i][j][pI], rcSumMinPU_VHiBin_Points_GenExclude_h[i][j][pI], rcSumMinPUFlow_VHiBin_Points_GenExclude_h[i][j][pI]});
      }
         
    
      for(Int_t pI = 0; pI < nRedCentBins; ++pI){
	Double_t tempMeanVal = redCentBinsMeanR0p3[pI]*rParam[j]*rParam[j]/(.3*.3);
	Double_t tempSigmaVal = redCentBinsSigmaR0p3[pI]*rParam[j]*rParam[j]/(.3*.3);

	rcSum_VHiBin_PointsRed_h[i][j][pI] = new TH1F(("rcSum_VHiBin_" + absEtaStr + "_" + rStr + "_hiBin" + std::to_string(pI) + "_PointsRed_h").c_str(), (";RC #Sigmap_{T} (" + rStr2 + ", " + absEtaStr2 + ");Counts").c_str(), 100, TMath::Max(0., tempMeanVal - tempSigmaVal*4.), tempMeanVal + tempSigmaVal*4.);
	pu_VHiBin_PointsRed_h[i][j][pI] = new TH1F(("pu_VHiBin_" + absEtaStr + "_" + rStr + "_hiBin" + std::to_string(pI) + "_PointsRed_h").c_str(), (";RC #rhoA (" + rStr2 + ", " + absEtaStr2 + ");Counts").c_str(), 100, TMath::Max(0., tempMeanVal - tempSigmaVal*4.), tempMeanVal + tempSigmaVal*4.);
	puFlow_VHiBin_PointsRed_h[i][j][pI] = new TH1F(("puFlow_VHiBin_" + absEtaStr + "_" + rStr + "_hiBin" + std::to_string(pI) + "_PointsRed_h").c_str(), (";RC #rhoA (" + rStr2 + ", " + absEtaStr2 + ");Counts").c_str(), 100, TMath::Max(0., tempMeanVal - tempSigmaVal*4.), tempMeanVal + tempSigmaVal*4.);
	rcSumMinPU_VHiBin_PointsRed_h[i][j][pI] = new TH1F(("rcSumMinPU_VHiBin_" + absEtaStr + "_" + rStr + "_hiBin" + std::to_string(pI) + "_PointsRed_h").c_str(), (";RC #rhoA (" + rStr2 + ", " + absEtaStr2 + ");Counts").c_str(), 100, -4.*tempSigmaVal, 4.*tempSigmaVal);
	rcSumMinPUFlow_VHiBin_PointsRed_h[i][j][pI] = new TH1F(("rcSumMinPUFlow_VHiBin_" + absEtaStr + "_" + rStr + "_hiBin" + std::to_string(pI) + "_PointsRed_h").c_str(), (";RC #rhoA (" + rStr2 + ", " + absEtaStr2 + ");Counts").c_str(), 100, -4.*tempSigmaVal, 4.*tempSigmaVal);

	if(isMC && !isMB){
	rcSum_VHiBin_PointsRed_GenExclude_h[i][j][pI] = new TH1F(("rcSum_VHiBin_" + absEtaStr + "_" + rStr + "_hiBin" + std::to_string(pI) + "_PointsRed_GenExclude_h").c_str(), (";RC #Sigmap_{T} (" + rStr2 + ", " + absEtaStr2 + ");Counts").c_str(), 100, TMath::Max(0., tempMeanVal - tempSigmaVal*4.), tempMeanVal + tempSigmaVal*4.);
	pu_VHiBin_PointsRed_GenExclude_h[i][j][pI] = new TH1F(("pu_VHiBin_" + absEtaStr + "_" + rStr + "_hiBin" + std::to_string(pI) + "_PointsRed_GenExclude_h").c_str(), (";RC #rhoA (" + rStr2 + ", " + absEtaStr2 + ");Counts").c_str(), 100, TMath::Max(0., tempMeanVal - tempSigmaVal*4.), tempMeanVal + tempSigmaVal*4.);
	puFlow_VHiBin_PointsRed_GenExclude_h[i][j][pI] = new TH1F(("puFlow_VHiBin_" + absEtaStr + "_" + rStr + "_hiBin" + std::to_string(pI) + "_PointsRed_GenExclude_h").c_str(), (";RC #rhoA (" + rStr2 + ", " + absEtaStr2 + ");Counts").c_str(), 100, TMath::Max(0., tempMeanVal - tempSigmaVal*4.), tempMeanVal + tempSigmaVal*4.);
	rcSumMinPU_VHiBin_PointsRed_GenExclude_h[i][j][pI] = new TH1F(("rcSumMinPU_VHiBin_" + absEtaStr + "_" + rStr + "_hiBin" + std::to_string(pI) + "_PointsRed_GenExclude_h").c_str(), (";RC #rhoA (" + rStr2 + ", " + absEtaStr2 + ");Counts").c_str(), 100, -4.*tempSigmaVal, 4.*tempSigmaVal);
	rcSumMinPUFlow_VHiBin_PointsRed_GenExclude_h[i][j][pI] = new TH1F(("rcSumMinPUFlow_VHiBin_" + absEtaStr + "_" + rStr + "_hiBin" + std::to_string(pI) + "_PointsRed_GenExclude_h").c_str(), (";RC #rhoA (" + rStr2 + ", " + absEtaStr2 + ");Counts").c_str(), 100, -4.*tempSigmaVal, 4.*tempSigmaVal);
	}

	centerTitles({rcSum_VHiBin_PointsRed_h[i][j][pI], pu_VHiBin_PointsRed_h[i][j][pI], puFlow_VHiBin_PointsRed_h[i][j][pI], rcSumMinPU_VHiBin_PointsRed_h[i][j][pI], rcSumMinPUFlow_VHiBin_PointsRed_h[i][j][pI]});

	if(isMC && !isMB) centerTitles({rcSum_VHiBin_PointsRed_GenExclude_h[i][j][pI], pu_VHiBin_PointsRed_GenExclude_h[i][j][pI], puFlow_VHiBin_PointsRed_GenExclude_h[i][j][pI], rcSumMinPU_VHiBin_PointsRed_GenExclude_h[i][j][pI], rcSumMinPUFlow_VHiBin_PointsRed_GenExclude_h[i][j][pI]});
      }

      for(Int_t cI = 0; cI < nCentBins; ++cI){
	std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
	if(nCentBins == 1) centStr = "PP";

	for(Int_t eI = 0; eI < nEvtPlaneBins; ++eI){	
	  rcSum_h[cI][eI][i][j] = new TH1F(("rcSum_" + absEtaStr + "_" + rStr + "_" + centStr + "_" + evtPlaneBins[eI] + "_h").c_str(), (";RC #Sigmap_{T} (" + rStr2 + ", " + absEtaStr2 + ");Counts").c_str(), nRCBins[j][cI], rcBinsLow[j][cI], rcBinsHi[j][cI]);
	  
	  if(isMC && !isMB) rcSum_GenExclude_h[cI][eI][i][j] = new TH1F(("rcSum_" + absEtaStr + "_" + rStr + "_" + centStr + "_" + evtPlaneBins[eI] + "_GenExclude_h").c_str(), (";RC #Sigmap_{T} (" + rStr2 + ", " + absEtaStr2 + ", Exclude Gen.);Counts").c_str(), nRCBins[j][cI], rcBinsLow[j][cI], rcBinsHi[j][cI]);
	  
	  rcSumMinPU_h[cI][eI][i][j] = new TH1F(("rcSumMinPU_" + absEtaStr + "_" + rStr + "_" + centStr + "_" + evtPlaneBins[eI] + "_h").c_str(), (";RC #Sigmap_{T} - #rhoA (" + rStr2 + ", " + absEtaStr2 + ");Counts").c_str(), 2*nRCBins[j][cI], -rcBinsHi[j][cI], rcBinsHi[j][cI]);
	  
	  if(isMC && !isMB) rcSumMinPU_GenExclude_h[cI][eI][i][j] = new TH1F(("rcSumMinPU_" + absEtaStr + "_" + rStr + "_" + centStr + "_" + evtPlaneBins[eI] + "_GenExclude_h").c_str(), (";RC #Sigmap_{T} - #rhoA (" + rStr2 + ", " + absEtaStr2 + ", Exclude Gen.);Counts").c_str(), 2*nRCBins[j][cI], -rcBinsHi[j][cI], rcBinsHi[j][cI]);
	  
	  centerTitles({rcSum_h[cI][eI][i][j], rcSumMinPU_h[cI][eI][i][j]});
	  if(isMC && !isMB)centerTitles({rcSum_GenExclude_h[cI][eI][i][j], rcSumMinPU_GenExclude_h[cI][eI][i][j]});

	  for(Int_t fI = 0; fI < nFlow; ++fI){
	    std::string flowStr = "FlowDefaultInRho";
	    if(fI != 0){
	      std::string flowAnaStr = flowGlobal.at(fI-1);

	      flowStr = "Flow";
	      if(flowAnaStr.find("Free") != std::string::npos) flowStr = flowStr + "FreePlane";
	      else if(flowAnaStr.find("HF") != std::string::npos) flowStr = flowStr + "HFPlane";
	      if(flowAnaStr.find("Jetty") != std::string::npos) flowStr = flowStr + "JettyExclude";
	    }
	    
	    
	    rcSumMinPUFlow_h[cI][eI][i][j][fI] = new TH1F(("rcSumMinPUFlow_" + absEtaStr + "_" + rStr + "_" + centStr + "_" + evtPlaneBins[eI] + "_" + flowStr + "_h").c_str(), (";RC #Sigmap_{T} - #rhoA (" + rStr2 + ", " + absEtaStr2 + ");Counts").c_str(), 2*nRCBins[j][cI], -rcBinsHi[j][cI], rcBinsHi[j][cI]);
	    
	    if(isMC && !isMB) rcSumMinPUFlow_GenExclude_h[cI][eI][i][j][fI] = new TH1F(("rcSumMinPUFlow_" + absEtaStr + "_" + rStr + "_" + centStr + "_" + evtPlaneBins[eI] + "_GenExclude_" + flowStr + "_h").c_str(), (";RC #Sigmap_{T} - #rhoA (" + rStr2 + ", " + absEtaStr2 + ", Exclude Gen.);Counts").c_str(), 2*nRCBins[j][cI], -rcBinsHi[j][cI], rcBinsHi[j][cI]);
	    
	    centerTitles({rcSumMinPUFlow_h[cI][eI][i][j][fI]});
	    if(isMC && !isMB)centerTitles({rcSumMinPUFlow_GenExclude_h[cI][eI][i][j][fI]});
	  }
	}
      }
    }
  }

  for(Int_t cI = 0; cI < nCentBins; ++cI){ 
    std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
    if(nCentBins == 1) centStr = "PP";
    
    for(Int_t pI = 0; pI < nPFID+1; ++pI){
      std::string pfStr = "AllPFID";
      if(pI != 0) pfStr = "PFID" + std::to_string(pI);
      
      pfCandMap_Mean_h[cI][pI] = new TH1F(("pfCandMap_Mean_" + centStr + "_" + pfStr + "_h").c_str(), ";Particle Flow #eta;.087xd#LTp_{T}#GT/d#eta", nEtaTowBoundsHist, etaTowBoundsHist);
      pfCandMap_Sigma_h[cI][pI] = new TH1F(("pfCandMap_Sigma_" + centStr + "_" + pfStr + "_h").c_str(), ";Particle Flow #eta;#sigma(.087xdp_{T}/d#eta)", nEtaTowBoundsHist, etaTowBoundsHist);
      
      pfCsCandMap_Mean_h[cI][pI] = new TH1F(("pfCsCandMap_Mean_" + centStr + "_" + pfStr + "_h").c_str(), ";Particle Flow #eta;.087xd#LTp_{T}#GT/d#eta", nEtaTowBoundsHist, etaTowBoundsHist);
      pfCsCandMap_Sigma_h[cI][pI] = new TH1F(("pfCsCandMap_Sigma_" + centStr + "_" + pfStr + "_h").c_str(), ";Particle Flow #eta;#sigma(.087xdp_{T}/d#eta)", nEtaTowBoundsHist, etaTowBoundsHist);
      
      centerTitles({pfCandMap_Mean_h[cI][pI], pfCandMap_Sigma_h[cI][pI], pfCsCandMap_Mean_h[cI][pI], pfCsCandMap_Sigma_h[cI][pI]});
      
      
      for(Int_t bI = 0; bI < nEtaTowBoundsHist; ++bI){
	pfCandMap_Points_h[cI][pI][bI] = new TH1F(("pfCandMap_Points_" + centStr + "_" + pfStr + "_Bin" + std::to_string(bI) + "_h").c_str(), ";PF p_{T} Sum / #Delta#eta;Counts", 500, 0, 500);
	pfCsCandMap_Points_h[cI][pI][bI] = new TH1F(("pfCsCandMap_Points_" + centStr + "_" + pfStr + "_Bin" + std::to_string(bI) + "_h").c_str(), ";PF p_{T} Sum / #Delta#eta;Counts", 500, 0, 500);
	centerTitles({pfCandMap_Points_h[cI][pI][bI], pfCsCandMap_Points_h[cI][pI][bI]});
      }
    }
  }
  
  inFile_p = TFile::Open(mntToXRootdFileString(fileList.at(0)).c_str(), "READ");
  std::vector<std::string> jetTreesGlobal = getUniqueStrings(returnRootFileContentsList(inFile_p, "TTree", "JetAnalyzer"));
  if(isPP){
    unsigned int jetPos = 0;
    while(jetPos < jetTreesGlobal.size()){
      if(jetTreesGlobal.at(jetPos).find("Cs") != std::string::npos) jetTreesGlobal.erase(jetTreesGlobal.begin()+jetPos);
      else if(jetTreesGlobal.at(jetPos).find("Pu") != std::string::npos) jetTreesGlobal.erase(jetTreesGlobal.begin()+jetPos);
      else if(jetTreesGlobal.at(jetPos).find("CS") != std::string::npos) jetTreesGlobal.erase(jetTreesGlobal.begin()+jetPos);
      else if(jetTreesGlobal.at(jetPos).find("PU") != std::string::npos) jetTreesGlobal.erase(jetTreesGlobal.begin()+jetPos);
      else jetPos++;
    }
  }
  else{
    unsigned int jetPos = 0;
    while(jetPos < jetTreesGlobal.size()){
      if(jetTreesGlobal.at(jetPos).find("ak1") != std::string::npos) jetTreesGlobal.erase(jetTreesGlobal.begin()+jetPos);
      else if(jetTreesGlobal.at(jetPos).find("ak2") != std::string::npos) jetTreesGlobal.erase(jetTreesGlobal.begin()+jetPos);
      else if(jetTreesGlobal.at(jetPos).find("ak3") != std::string::npos) jetTreesGlobal.erase(jetTreesGlobal.begin()+jetPos);
      else if(jetTreesGlobal.at(jetPos).find("ak4") != std::string::npos) jetTreesGlobal.erase(jetTreesGlobal.begin()+jetPos);
      else if(jetTreesGlobal.at(jetPos).find("ak5") != std::string::npos) jetTreesGlobal.erase(jetTreesGlobal.begin()+jetPos);
      else if(jetTreesGlobal.at(jetPos).find("ak6") != std::string::npos) jetTreesGlobal.erase(jetTreesGlobal.begin()+jetPos);
      else if(jetTreesGlobal.at(jetPos).find("ak8") != std::string::npos) jetTreesGlobal.erase(jetTreesGlobal.begin()+jetPos);
      else if(jetTreesGlobal.at(jetPos).find("ak10") != std::string::npos) jetTreesGlobal.erase(jetTreesGlobal.begin()+jetPos);
      else jetPos++;
    }
  }

  const Int_t nJetAlgo = jetTreesGlobal.size();
  Int_t rParamToJetAlgo[nRParam];
  for(Int_t rI = 0; rI < nRParam; ++rI){
    
    for(Int_t jI = 0; jI < nJetAlgo; ++jI){
      if(jetTreesGlobal.at(jI).find(std::to_string(((int)(rParam[rI]*10)))) != std::string::npos){
	std::cout << "Mapping rparam \'" << rParam[rI] << "\' onto jetalgo \'" << jetTreesGlobal.at(jI) << "\'." << std::endl;
	rParamToJetAlgo[rI] = jI;
	break;
      }
    }
  }

  inFile_p->Close();
  delete inFile_p;

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  Int_t nFiles = TMath::Min((Int_t)fileList.size(), 1000);
  for(Int_t fI = 0; fI < nFiles; ++fI){
    std::cout << "Processing file " << fI << "/" << nFiles << " \'" << fileList.at(fI) << "\'..." << std::endl;

    inFile_p = TFile::Open(mntToXRootdFileString(fileList.at(fI)).c_str(), "READ");
    std::vector<std::string> jetTreesLocal = getUniqueStrings(returnRootFileContentsList(inFile_p, "TTree", "JetAnalyzer"));
    Bool_t isGood = true;
    for(Int_t i = 0; i < nJetAlgo; ++i){
      bool isGood2 = false;
      for(unsigned int j = 0; j < jetTreesLocal.size(); ++j){
	if(jetTreesGlobal.at(i).size() == jetTreesLocal.at(j).size() && jetTreesLocal.at(j).find(jetTreesGlobal.at(i)) != std::string::npos){
	  isGood2 = true;
	  break;
	}
      }

      isGood = isGood && isGood2;
      if(!isGood) break;
    }

    if(!isGood){
      std::cout << "Warning, file \'" << fileList.at(fI) << "\' is missing jet collections, continue" << std::endl;
      inFile_p->Close();
      delete inFile_p;
      continue;
    }

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
    TTree* pfTree_p = NULL; 
    TTree* pfCsTree_p = NULL; 
    TTree* rhoTree_p = NULL;
    TTree* hltTree_p = NULL;
    TTree* skimTree_p = NULL;
    std::vector<std::vector<double>*> rhoFlowFitParamsMore_p;
    rhoFlowFitParamsMore_p.reserve(nFlow);
    for(int fI = 0; fI < nFlow; ++fI){
      flowTree_p[fI] = NULL;
      rhoFlowFitParamsMore_p.push_back(NULL);
    }

    if(hasPF) pfTree_p = (TTree*)inFile_p->Get("pfcandAnalyzer/pfTree");
    if(hasPFCs) pfCsTree_p = (TTree*)inFile_p->Get(pfCsTreesGlobal.at(0).c_str());
    if(hasRho) rhoTree_p = (TTree*)inFile_p->Get("hiPuRhoR3Analyzer/t");
    if(hasHLT) hltTree_p = (TTree*)inFile_p->Get("hltanalysis/HltTree");
    if(hasSkim) skimTree_p = (TTree*)inFile_p->Get("skimanalysis/HltTree");
    if(flowGlobal.size() > 0){
      for(unsigned int fI = 0; fI < flowGlobal.size(); ++fI){
	std::cout << " Flow " << fI << ": " << flowGlobal.at(fI) << std::endl;
	flowTree_p[fI] = (TTree*)inFile_p->Get(flowGlobal.at(fI).c_str());
      }
    }

    TTree* jetTree_p[nJetAlgo];

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    for(Int_t i = 0; i < nRParam; ++i){
      jetTree_p[i] = NULL;
      if(isMC && !isMB) jetTree_p[i] = (TTree*)inFile_p->Get(jetTreesGlobal.at(rParamToJetAlgo[i]).c_str());
    }

    Float_t vz_;
    Int_t hiBin_;
    Int_t hiNevtPlane_;
    Float_t hiEvtPlanes_[29];

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    
    hiTree_p->SetBranchStatus("*", 0);
    hiTree_p->SetBranchStatus("run", 1);
    hiTree_p->SetBranchStatus("lumi", 1);
    hiTree_p->SetBranchStatus("evt", 1);
    hiTree_p->SetBranchStatus("vz", 1);
    hiTree_p->SetBranchStatus("hiBin", 1);
    hiTree_p->SetBranchStatus("hiNevtPlane", 1);
    hiTree_p->SetBranchStatus("hiEvtPlanes", 1);

    hiTree_p->SetBranchAddress("run", &run_);
    hiTree_p->SetBranchAddress("lumi", &lumi_);
    hiTree_p->SetBranchAddress("evt", &evt_);
    hiTree_p->SetBranchAddress("vz", &vz_);
    hiTree_p->SetBranchAddress("hiBin", &hiBin_);
    hiTree_p->SetBranchAddress("hiNevtPlane", &hiNevtPlane_);
    hiTree_p->SetBranchAddress("hiEvtPlanes", hiEvtPlanes_);

    std::vector<float>* pfPt_p = NULL;
    std::vector<float>* pfPhi_p = NULL;
    std::vector<float>* pfEta_p = NULL;
    std::vector<int>* pfId_p = NULL;

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    if(hasPF){
      pfTree_p->SetBranchStatus("*", 0);
      pfTree_p->SetBranchStatus("pfPt", 1);
      pfTree_p->SetBranchStatus("pfPhi", 1);
      pfTree_p->SetBranchStatus("pfEta", 1);
      pfTree_p->SetBranchStatus("pfId", 1);
      
      pfTree_p->SetBranchAddress("pfPt", &pfPt_p);
      pfTree_p->SetBranchAddress("pfPhi", &pfPhi_p);
      pfTree_p->SetBranchAddress("pfEta", &pfEta_p);
      pfTree_p->SetBranchAddress("pfId", &pfId_p);
    }

    std::vector<float>* pfPtCs_p = NULL;
    std::vector<float>* pfPhiCs_p = NULL;
    std::vector<float>* pfEtaCs_p = NULL;
    std::vector<int>* pfIdCs_p = NULL;
    
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    if(hasPFCs){
      pfCsTree_p->SetBranchStatus("*", 0);
      pfCsTree_p->SetBranchStatus("pfPt", 1);
      pfCsTree_p->SetBranchStatus("pfPhi", 1);
      pfCsTree_p->SetBranchStatus("pfEta", 1);
      pfCsTree_p->SetBranchStatus("pfId", 1);
      
      pfCsTree_p->SetBranchAddress("pfPt", &pfPtCs_p);
      pfCsTree_p->SetBranchAddress("pfPhi", &pfPhiCs_p);
      pfCsTree_p->SetBranchAddress("pfEta", &pfEtaCs_p);
      pfCsTree_p->SetBranchAddress("pfId", &pfIdCs_p);
    }

    std::vector<double>* rho_p = NULL;
    std::vector<double>* etaMin_p = NULL;
    std::vector<double>* etaMax_p = NULL;
    std::vector<double>* rhoFlowFitParams_p = NULL;

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    if(hasRho){
      rhoTree_p->SetBranchStatus("*", 0);
      rhoTree_p->SetBranchStatus("rho", 1);
      rhoTree_p->SetBranchStatus("etaMin", 1);
      rhoTree_p->SetBranchStatus("etaMax", 1);
      rhoTree_p->SetBranchStatus("rhoFlowFitParams", 1);

      rhoTree_p->SetBranchAddress("rho", &rho_p);
      rhoTree_p->SetBranchAddress("etaMin", &etaMin_p);
      rhoTree_p->SetBranchAddress("etaMax", &etaMax_p);
      rhoTree_p->SetBranchAddress("rhoFlowFitParams", &rhoFlowFitParams_p);
    }

    if(flowGlobal.size() > 0){

      for(unsigned int fI = 0; fI < flowGlobal.size(); ++fI){
	flowTree_p[fI]->SetBranchAddress("rhoFlowFitParams", &(rhoFlowFitParamsMore_p.at(fI)) );
      }
    }

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;


    Int_t HLT_HIPuAK4CaloJet80_Eta5p1_v;
    Int_t HLT_HIPuAK4CaloJet100_Eta5p1_v;

    if(hasHLT){
      hltTree_p->SetBranchStatus("*", 0);
      if(hltTree_p->GetListOfBranches()->FindObject("HLT_HIPuAK4CaloJet80_Eta5p1_v1")) hltTree_p->SetBranchStatus("HLT_HIPuAK4CaloJet80_Eta5p1_v1", 1);
      else if(hltTree_p->GetListOfBranches()->FindObject("HLT_HIPuAK4CaloJet80_Eta5p1_v2")) hltTree_p->SetBranchStatus("HLT_HIPuAK4CaloJet80_Eta5p1_v2", 1);
      if(hltTree_p->GetListOfBranches()->FindObject("HLT_HIPuAK4CaloJet100_Eta5p1_v1")) hltTree_p->SetBranchStatus("HLT_HIPuAK4CaloJet100_Eta5p1_v1", 1);
      else if(hltTree_p->GetListOfBranches()->FindObject("HLT_HIPuAK4CaloJet100_Eta5p1_v2")) hltTree_p->SetBranchStatus("HLT_HIPuAK4CaloJet100_Eta5p1_v2", 1);

      if(hltTree_p->GetListOfBranches()->FindObject("HLT_HIPuAK4CaloJet80_Eta5p1_v1")) hltTree_p->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_v1", &HLT_HIPuAK4CaloJet80_Eta5p1_v);
      else if(hltTree_p->GetListOfBranches()->FindObject("HLT_HIPuAK4CaloJet80_Eta5p1_v2")) hltTree_p->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_v2", &HLT_HIPuAK4CaloJet80_Eta5p1_v);
      if(hltTree_p->GetListOfBranches()->FindObject("HLT_HIPuAK4CaloJet100_Eta5p1_v1")) hltTree_p->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_v1", &HLT_HIPuAK4CaloJet100_Eta5p1_v);
      else if(hltTree_p->GetListOfBranches()->FindObject("HLT_HIPuAK4CaloJet100_Eta5p1_v2")) hltTree_p->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_v2", &HLT_HIPuAK4CaloJet100_Eta5p1_v);
    }


    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    Int_t pprimaryVertexFilter_;
    Int_t pcollisionEventSelection_;
    Int_t HBHENoiseFilterResultRun2Loose_;
    Int_t pclusterCompatibilityFilter_;
    Int_t phfCoincFilter3_;

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    if(hasSkim){
      skimTree_p->SetBranchStatus("*", 0);
      skimTree_p->SetBranchStatus("pcollisionEventSelection", 1);
      skimTree_p->SetBranchStatus("pprimaryVertexFilter", 1);
      skimTree_p->SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
      skimTree_p->SetBranchStatus("phfCoincFilter3", 1);
      skimTree_p->SetBranchStatus("pclusterCompatibilityFilter", 1);

      skimTree_p->SetBranchAddress("pcollisionEventSelection", &pcollisionEventSelection_);
      skimTree_p->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter_);
      skimTree_p->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose_);
      skimTree_p->SetBranchAddress("phfCoincFilter3", &phfCoincFilter3_);
      skimTree_p->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter_);
    }

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    Int_t ngen_[nJetAlgo];
    Float_t genpt_[nJetAlgo][nMaxJets];
    Float_t genphi_[nJetAlgo][nMaxJets];
    Float_t geneta_[nJetAlgo][nMaxJets];
    Int_t gensubid_[nJetAlgo][nMaxJets];

    if(isMC && !isMB){
      for(Int_t i = 0; i < nRParam; ++i){
	jetTree_p[i]->SetBranchStatus("*", 0);
	jetTree_p[i]->SetBranchStatus("ngen", 1);
	jetTree_p[i]->SetBranchStatus("genpt", 1);
	jetTree_p[i]->SetBranchStatus("genphi", 1);
	jetTree_p[i]->SetBranchStatus("geneta", 1);
	jetTree_p[i]->SetBranchStatus("gensubid", 1);
	
	jetTree_p[i]->SetBranchAddress("ngen", &(ngen_[i]));
	jetTree_p[i]->SetBranchAddress("genpt", genpt_[i]);
	jetTree_p[i]->SetBranchAddress("genphi", genphi_[i]);
	jetTree_p[i]->SetBranchAddress("geneta", geneta_[i]);
	jetTree_p[i]->SetBranchAddress("gensubid", gensubid_[i]);
      }
    }

    std::cout << "Processing..." << std::endl;
    
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    const ULong64_t maxNEntries = 100000000;
    const ULong64_t nEntries = TMath::Min((ULong64_t)hiTree_p->GetEntries(), maxNEntries);
    for(ULong64_t entry = 0; entry < nEntries; ++entry){
      if(entry%1000 == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;
      hiTree_p->GetEntry(entry);

      if(run_ != 262620 || lumi_ != 204 || evt_ != 10438474) continue;

      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      if(TMath::Abs(vz_) > 15.) continue;

      if(hasSkim){
	skimTree_p->GetEntry(entry);
	if(!pcollisionEventSelection_) continue;
	if(!HBHENoiseFilterResultRun2Loose_) continue;
	if(!pprimaryVertexFilter_) continue;
	if(!pclusterCompatibilityFilter_) continue;
	if(!phfCoincFilter3_) continue;
      }

      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
      
      if(hasHLT){
	hltTree_p->GetEntry(entry);
	if(isMB && !isMC){
	  if(HLT_HIPuAK4CaloJet80_Eta5p1_v) continue;
	  else if(HLT_HIPuAK4CaloJet100_Eta5p1_v) continue;
	}
      }

      if(hasPF) pfTree_p->GetEntry(entry);
      if(hasPFCs) pfCsTree_p->GetEntry(entry);
      if(hasRho) rhoTree_p->GetEntry(entry);

      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      if(flowGlobal.size() > 0){
	for(unsigned int fI = 0; fI < flowGlobal.size(); ++fI){
	  flowTree_p[fI]->GetEntry(entry);
	}
      }
      
      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      Int_t centPos = 0;
      if(!isPP){
	for(Int_t cI = 0; cI < nCentBins; ++cI){
	  if(hiBin_/2 >= centBinsLow[cI] && hiBin_/2 < centBinsHi[cI]){
	    centPos = cI;
	  }
	}
      }

      const Double_t nHistProb = 100.;
      const Double_t histNum = randGen_p->Uniform(0, nEntries);

      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      if(hasRho){
	double val = ROOT::Math::chisquared_cdf_c(rhoFlowFitParams_p->at(5), rhoFlowFitParams_p->at(6));
	bool minProb = val > .05;
	bool maxProb = val < .95;

	if(doFlowExamples && ((hasPF && minProb && maxProb && histNum < nHistProb) || (run_ == 262620 && lumi_ == 204 && evt_ == 10438474))){
	  std::vector<float> finalPhis;
	  std::vector<float> finalPts;
	  
	  std::vector<float> finalPhis_Side;
	  std::vector<float> finalPts_Side;
	  
	  for(unsigned int pfI = 0; pfI < pfPt_p->size(); ++pfI){
	    if(pfPt_p->at(pfI) < .3) continue;
	    if(pfPt_p->at(pfI) >= 3.) continue;
	    if(pfId_p->at(pfI) != 1) continue;
	    
	    if(TMath::Abs(pfEta_p->at(pfI)) > 1.){
	      finalPhis_Side.push_back(pfPhi_p->at(pfI));
	      finalPts_Side.push_back(pfPt_p->at(pfI));
	    }
	    else{
	      finalPhis.push_back(pfPhi_p->at(pfI));
	      finalPts.push_back(pfPt_p->at(pfI));
	    }
	  }
	  
	  const Int_t nBins = TMath::Min(21, TMath::Max(10, (Int_t)finalPhis.size()/10));
	  const Int_t nBins_Side = TMath::Min(21, TMath::Max(10, (Int_t)finalPhis_Side.size()/10));
	  
	  TCanvas* tempCanv_p = new TCanvas("temp", "temp", 500, 500);
	  tempCanv_p->SetTopMargin(0.075);
	  tempCanv_p->SetRightMargin(0.01);
	  tempCanv_p->SetBottomMargin(0.10);
	  tempCanv_p->SetLeftMargin(0.125);
	  gStyle->SetOptStat(0);
	  TCanvas* tempCanvWeight_p = new TCanvas("tempWeight", "tempWeight", 500, 500);
	  tempCanvWeight_p->SetTopMargin(0.075);
	  tempCanvWeight_p->SetRightMargin(0.01);
	  tempCanvWeight_p->SetBottomMargin(0.10);
	  tempCanvWeight_p->SetLeftMargin(0.125);
	  gStyle->SetOptStat(0);
	  TCanvas* tempCanvSide_p = new TCanvas("tempSide", "tempSide", 500, 500);
	  tempCanvSide_p->SetTopMargin(0.075);
	  tempCanvSide_p->SetRightMargin(0.01);
	  tempCanvSide_p->SetBottomMargin(0.10);
	  tempCanvSide_p->SetLeftMargin(0.125);
	  gStyle->SetOptStat(0);
	  TH1F* tempHist_p = new TH1F("tempHist_h", (";#phi;Counts (.3 < p_{T} < 3., #color[" + std::to_string(kRed) + "]{|#eta| < 1.}, ID=1)").c_str(), nBins, -TMath::Pi(), TMath::Pi());
	  TH1F* tempHistWeight_p = new TH1F("tempHistWeight_h", ";#phi;Counts(PtWeight)", nBins, -TMath::Pi(), TMath::Pi());
	  TH1F* tempHistSide_p = new TH1F("tempHistSide_h", (";#phi;Counts (.3 < p_{T} < 3., #color[" + std::to_string(kRed) + "]{1. < |#eta| < 2.}, ID=1)").c_str(), nBins_Side, -TMath::Pi(), TMath::Pi());
	  tempHistWeight_p->Sumw2();

	  std::string flowFitForm = "[0]*(1 + 2*(" + std::to_string(rhoFlowFitParams_p->at(1)) + "*TMath::Cos(2*(x - " + std::to_string(rhoFlowFitParams_p->at(2)) + ")) + " + std::to_string(rhoFlowFitParams_p->at(3)) + "*TMath::Cos(3*(x - " + std::to_string(rhoFlowFitParams_p->at(4)) + "))))";
	  TF1* overlay_p = new TF1("overlay_p", flowFitForm.c_str(), -TMath::Pi(), TMath::Pi());
	  overlay_p->SetParameter(0, rhoFlowFitParams_p->at(0));
	  
	  TF1* overlayWeight_p = new TF1("overlayWeight_p", flowFitForm.c_str(), -TMath::Pi(), TMath::Pi());
	  overlayWeight_p->SetParameter(0, rhoFlowFitParams_p->at(0));
	  
	  TF1* overlaySide_p = new TF1("overlaySide_p", flowFitForm.c_str(), -TMath::Pi(), TMath::Pi());
	  overlaySide_p->SetParameter(0, rhoFlowFitParams_p->at(0));
	  
	  for(unsigned int pfI = 0; pfI < finalPhis.size(); ++pfI){
	    tempHist_p->Fill(finalPhis.at(pfI));
	    tempHistWeight_p->Fill(finalPhis.at(pfI), pfPt_p->at(pfI));
	  }
	  
	  for(unsigned int pfI = 0; pfI < finalPhis_Side.size(); ++pfI){
	    tempHistSide_p->Fill(finalPhis_Side.at(pfI));
	  }
	  
	  tempHist_p->Fit("overlay_p", "", "Q N", -TMath::Pi(), TMath::Pi());
	  tempHistWeight_p->Fit("overlayWeight_p", "", "Q N", -TMath::Pi(), TMath::Pi());
	  tempHistSide_p->Fit("overlaySide_p", "", "Q N", -TMath::Pi(), TMath::Pi());
	  finalPhis.clear();
	  finalPts.clear();
	  finalPhis_Side.clear();
	  finalPts_Side.clear();
	  
	  std::string canvSaveName = "pdfDir/phiFit_hiBin" + std::to_string(hiBin_) + "_Run" + std::to_string(run_) + "_Lumi" + std::to_string(lumi_) + "_Evt" + std::to_string(evt_) + "_Entry" + std::to_string(entry) + mbStr + mcStr + "_" + std::to_string(date->GetDate()) + ".pdf";
	  tempCanv_p->cd();
	  tempCanv_p->SetTicks(0, 1);
	  tempHist_p->SetMinimum(0.0);
	  tempHist_p->SetMarkerStyle(20);
	  tempHist_p->SetMarkerSize(0.6);
	  tempHist_p->SetMarkerColor(col.getColor(2));
	  tempHist_p->SetLineColor(col.getColor(2));
	  
	  tempHist_p->GetXaxis()->CenterTitle();
	  tempHist_p->GetYaxis()->CenterTitle();
	  tempHist_p->DrawCopy("E1 P");
	  overlay_p->SetMarkerColor(1);
	  overlay_p->SetLineColor(1);
	  //	overlay_p->SetLineStyle(2);
	  overlay_p->DrawCopy("SAME");
	  
	  flowFitForm = std::to_string(overlay_p->GetParameter(0)) + "*(1 + 2*(" + std::to_string(rhoFlowFitParams_p->at(1)) + "*TMath::Cos(2*(x - " + std::to_string(rhoFlowFitParams_p->at(2)) + "))))";
	  TF1* overlay2_p = new TF1("overlay2_p", flowFitForm.c_str(), -TMath::Pi(), TMath::Pi());
	  overlay2_p->SetLineStyle(2);
	  overlay2_p->SetLineColor(kRed);
	  overlay2_p->SetMarkerColor(kRed);
	  overlay2_p->DrawCopy("SAME");
	  
	  flowFitForm = std::to_string(overlay_p->GetParameter(0)) + "*(1 + 2*(" + std::to_string(rhoFlowFitParams_p->at(3)) + "*TMath::Cos(3*(x - " + std::to_string(rhoFlowFitParams_p->at(4)) + "))))";
	  TF1* overlay3_p = new TF1("overlay3_p", flowFitForm.c_str(), -TMath::Pi(), TMath::Pi());
	  overlay3_p->SetLineStyle(2);
	  overlay3_p->SetLineColor(kBlue);
	  overlay3_p->SetMarkerColor(kBlue);
	  overlay3_p->DrawCopy("SAME");
	  
	  TLatex* label_p = new TLatex();
	  label_p->SetTextFont(43);
	  label_p->SetTextSize(12);
	  label_p->SetNDC();
	  
	  label_p->DrawLatex(.13, .95, ("#rho_{0}(1 + 2#color[" + std::to_string(kRed) + "]{v_{2}}cos(2[#phi - #color[" + std::to_string(kRed) + "]{#phi_{2}}]) + 2#color[" + std::to_string(kBlue) + "]{v_{3}}cos(3[#phi - #color[" + std::to_string(kBlue) + "]{#phi_{3}}]))").c_str());
	  label_p->DrawLatex(.58, .9, ("Run=" + std::to_string(run_) + ";Lumi=" + std::to_string(lumi_) + ";Evt=" + std::to_string(evt_)).c_str());
	  label_p->DrawLatex(.13, .9, ("#color[" + std::to_string(kRed) + "]{v_{2}=" + prettyString(rhoFlowFitParams_p->at(1), 3, false) + "}").c_str());
	  label_p->DrawLatex(.23, .9, ("#color[" + std::to_string(kRed) + "]{#phi_{2}=" + prettyString(rhoFlowFitParams_p->at(2), 3, false) + "}").c_str());
	  label_p->DrawLatex(.33, .9, ("#color[" + std::to_string(kBlue) + "]{v_{3}=" + prettyString(rhoFlowFitParams_p->at(3), 3, false) + "}").c_str());
	  label_p->DrawLatex(.43, .9, ("#color[" + std::to_string(kBlue) + "]{#phi_{3}=" + prettyString(rhoFlowFitParams_p->at(4), 3, false) + "}").c_str());
	  label_p->DrawLatex(.6, .95, ("hiBin=" + std::to_string(hiBin_)).c_str());
	  
	  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	  
	  tempCanv_p->SaveAs(canvSaveName.c_str());
	  canvSaveName.replace(canvSaveName.find(".pdf"), 4, ".C");
	  tempCanv_p->SaveAs(canvSaveName.c_str());
	  
	  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	  
	  canvSaveName = "pdfDir/phiFitWeight_hiBin" + std::to_string(hiBin_) + "_Run" + std::to_string(run_) + "_Lumi" + std::to_string(lumi_) + "_Evt" + std::to_string(evt_) + mbStr + mcStr + "_" + std::to_string(date->GetDate()) + ".pdf";
	  tempCanvWeight_p->cd();
	  tempHistWeight_p->SetMinimum(0.0);
	  tempHistWeight_p->SetMarkerStyle(20);
	  tempHistWeight_p->SetMarkerSize(0.6);
	  tempHistWeight_p->SetMarkerColor(col.getColor(2));
	  tempHistWeight_p->SetLineColor(col.getColor(2));
	  tempHistWeight_p->GetXaxis()->CenterTitle();
	  tempHistWeight_p->GetYaxis()->CenterTitle();
	  tempHistWeight_p->DrawCopy("E1 P");
	  overlayWeight_p->SetLineStyle(2);
	  overlayWeight_p->DrawCopy("SAME");
	  // tempCanvWeight_p->SaveAs(canvSaveName.c_str());
	  
	  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	  
	  canvSaveName = "pdfDir/phiFitSide_hiBin" + std::to_string(hiBin_) + "_Run" + std::to_string(run_) + "_Lumi" + std::to_string(lumi_) + "_Evt" + std::to_string(evt_) + "_Entry" + std::to_string(entry) + mbStr + mcStr + "_" + std::to_string(date->GetDate()) + ".pdf";
	  tempCanvSide_p->cd();
	  tempHistSide_p->SetMinimum(0.0);
	  tempHistSide_p->SetMarkerStyle(20);
	  tempHistSide_p->SetMarkerSize(0.6);
	  tempHistSide_p->SetMarkerColor(col.getColor(2));
	  tempHistSide_p->SetLineColor(col.getColor(2));
	  tempHistSide_p->GetXaxis()->CenterTitle();
	  tempHistSide_p->GetYaxis()->CenterTitle();
	  tempHistSide_p->DrawCopy("E1 P");
	  overlaySide_p->SetLineStyle(2);
	  overlaySide_p->DrawCopy("SAME");
	  
	  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	  
	  label_p->DrawLatex(.13, .95, ("#rho_{0}(1 + 2#color[" + std::to_string(kRed) + "]{v_{2}}cos(2[#phi - #color[" + std::to_string(kRed) + "]{#phi_{2}}]) + 2#color[" + std::to_string(kBlue) + "]{v_{3}}cos(3[#phi - #color[" + std::to_string(kBlue) + "]{#phi_{3}}]))").c_str());
	  label_p->DrawLatex(.58, .9, ("Run=" + std::to_string(run_) + ";Lumi=" + std::to_string(lumi_) + ";Evt=" + std::to_string(evt_)).c_str());
	  label_p->DrawLatex(.13, .9, ("#color[" + std::to_string(kRed) + "]{v_{2}=" + prettyString(rhoFlowFitParams_p->at(1), 3, false) + "}").c_str());
	  label_p->DrawLatex(.23, .9, ("#color[" + std::to_string(kRed) + "]{#phi_{2}=" + prettyString(rhoFlowFitParams_p->at(2), 3, false) + "}").c_str());
	  label_p->DrawLatex(.33, .9, ("#color[" + std::to_string(kBlue) + "]{v_{3}=" + prettyString(rhoFlowFitParams_p->at(3), 3, false) + "}").c_str());
	  label_p->DrawLatex(.43, .9, ("#color[" + std::to_string(kBlue) + "]{#phi_{3}=" + prettyString(rhoFlowFitParams_p->at(4), 3, false) + "}").c_str());
	  label_p->DrawLatex(.6, .95, ("hiBin=" + std::to_string(hiBin_)).c_str());
	  
	  tempCanvSide_p->SaveAs(canvSaveName.c_str());                                                              
	  canvSaveName.replace(canvSaveName.find(".pdf"), 4, ".C");
          tempCanvSide_p->SaveAs(canvSaveName.c_str());

	  delete overlayWeight_p;
	  delete tempHistWeight_p;
	  delete tempCanvWeight_p;
	  
	  delete overlaySide_p;
	  delete tempHistSide_p;
	  delete tempCanvSide_p;
	  
	  delete overlay_p;
	  delete overlay2_p;
	  delete overlay3_p;
	  delete tempHist_p;
	  delete tempCanv_p;
	  
	  delete label_p;
	}
      }
    
      if(hasPF){
	Float_t pfEtaSum[nPFID+1][nEtaTowBoundsHist];
	for(Int_t idI = 0; idI < nPFID+1; ++idI){
	  for(Int_t pfI = 0; pfI < nEtaTowBoundsHist; ++pfI){pfEtaSum[idI][pfI] = 0.0;}
	}

	for(unsigned int pI = 0; pI < pfPt_p->size(); ++pI){
	  Int_t etaBinPos = findBinPos(pfEta_p->at(pI), nEtaTowBoundsHist, etaTowBoundsHist);
	  pfEtaSum[0][etaBinPos] += pfPt_p->at(pI);
	  pfEtaSum[pfId_p->at(pI)][etaBinPos] += pfPt_p->at(pI);
	}

	for(Int_t idI = 0; idI < nPFID+1; ++idI){
	  for(Int_t pfI = 0; pfI < nEtaTowBoundsHist; ++pfI){
	    Float_t etaBinWidth = (etaTowBoundsHist[pfI+1] - etaTowBoundsHist[pfI])/.087;
	    pfCandMap_Points_h[centPos][idI][pfI]->Fill(pfEtaSum[idI][pfI]/etaBinWidth);
	  }
	}
      }

      if(hasPFCs){
	Float_t pfEtaSum[nPFID+1][nEtaTowBoundsHist];
	for(Int_t idI = 0; idI < nPFID+1; ++idI){
	  for(Int_t pfI = 0; pfI < nEtaTowBoundsHist; ++pfI){pfEtaSum[idI][pfI] = 0.0;}
	}

	for(unsigned int pI = 0; pI < pfPtCs_p->size(); ++pI){
	  Int_t etaBinPos = findBinPos(pfEtaCs_p->at(pI), nEtaTowBoundsHist, etaTowBoundsHist);
	  pfEtaSum[0][etaBinPos] += pfPtCs_p->at(pI);
	  pfEtaSum[pfIdCs_p->at(pI)][etaBinPos] += pfPtCs_p->at(pI);
	}

	for(Int_t idI = 0; idI < nPFID+1; ++idI){
	  for(Int_t pfI = 0; pfI < nEtaTowBoundsHist; ++pfI){
	    Float_t etaBinWidth = (etaTowBoundsHist[pfI+1] - etaTowBoundsHist[pfI])/.087;
	    pfCandMap_Points_h[centPos][idI][pfI]->Fill(pfEtaSum[idI][pfI]/etaBinWidth);
	  }
	}
      }

      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
      
      if(isMC && !isMB) for(Int_t i = 0; i < nRParam; ++i){jetTree_p[i]->GetEntry(entry);} 
     
      for(unsigned int i = 0; i < etaNTowInPhi.size(); ++i){
	etaTowSumMap.at(i) = 0;

	for(Int_t j = 0; j < etaNTowInPhi.at(i); ++j){
	  etaPhiTowSumMap.at(i).at(j) = 0.0;
	}
      }

      for(Int_t i = 0; i < nRParam; ++i){
	for(unsigned int j = 0; j < etaNTowInPhi.size(); ++j){
	  etaTowNMapExcludeR.at(i).at(j) = 0;
	  etaTowSumMapExcludeR.at(i).at(j) = 0.0;
	}
      }


      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      //Temp for cone check
      float globalEta = -100;
      float globalPhi = -100;

      for(Int_t i = 0; i < nRParam; ++i){
	std::vector<float> etaVal;
	std::vector<float> phiVal;
	std::vector<float> sumVal;
	std::vector< std::vector<int> > etaPoses;
	std::vector< std::vector<int> > evtPlanePoses;

	std::vector<float> etaVal_GenExclude;
	std::vector<float> phiVal_GenExclude;
	std::vector<float> sumVal_GenExclude;
	std::vector< std::vector<int> > etaPoses_GenExclude;
	std::vector< std::vector<int> > evtPlanePoses_GenExclude;
      
	while(etaVal.size() < (unsigned int)nRCGen){

	  float tempEta = -100;
	  float tempPhi = -100;
	  if(i == 0){
	    tempEta = randGen_p->Uniform(-absEtaBinsHi[nAbsEtaBins-1], absEtaBinsHi[nAbsEtaBins-1]);
	    tempPhi = randGen_p->Uniform(-TMath::Pi(), TMath::Pi());	  
	    bool isGood = true;
	    for(unsigned int l = 0; l < etaVal.size(); ++l){
	      if(getDR(tempEta, tempPhi, etaVal[l], phiVal[l]) < rParam[i]){
		isGood = false;
		break;
	      }
	    }

	    if(!isGood) continue;

	    globalEta = tempEta;
	    globalPhi = tempPhi;
	  }
	  else{
	    tempEta = globalEta;
	    tempPhi = globalPhi;
	  }
	  	
	  etaVal.push_back(tempEta);
	  phiVal.push_back(tempPhi);
	  sumVal.push_back(0.0);
	  etaPoses.push_back({});
	  evtPlanePoses.push_back({0});
	  
	  for(Int_t aI = 0; aI < nAbsEtaBins; ++aI){
	    if(TMath::Abs(tempEta) >= absEtaBinsLow[aI] && TMath::Abs(tempEta) < absEtaBinsHi[aI]){
	      etaPoses[etaPoses.size()-1].push_back(aI);
	    }
	  }

	  if(hiNevtPlane_ > 0){
	    Double_t deltaPhi = TMath::Abs(getDPHI(tempPhi, hiEvtPlanes_[8]));
	    if(deltaPhi < TMath::Pi()/6. || deltaPhi > 5.*TMath::Pi()/6.) evtPlanePoses[evtPlanePoses.size()-1].push_back(1);
	    else if(deltaPhi < 2.*TMath::Pi()/6. || deltaPhi > 4.*TMath::Pi()/6.) evtPlanePoses[evtPlanePoses.size()-1].push_back(3);
	    else evtPlanePoses[evtPlanePoses.size()-1].push_back(2);
	  }
	}

	while(etaVal_GenExclude.size() < (unsigned int)nRCGen && isMC && !isMB){
	  float tempEta = randGen_p->Uniform(-absEtaBinsHi[nAbsEtaBins-1], absEtaBinsHi[nAbsEtaBins-1]);
	  float tempPhi = randGen_p->Uniform(-TMath::Pi(), TMath::Pi());
	  
	  bool isGood = true;
	  for(unsigned int l = 0; l < etaVal_GenExclude.size(); ++l){
	    if(getDR(tempEta, tempPhi, etaVal_GenExclude[l], phiVal_GenExclude[l]) < rParam[i]){
	      isGood = false;
	      break;
	    }
	  }

	  //	  if(i == nRParam - 1) std::cout << tempEta << ", " << tempPhi << std::endl;
	
	  
	  for(Int_t gI = 0; gI < ngen_[i]; ++gI){
	    if(genpt_[i][gI] < 15.) continue;
	    if(gensubid_[i][gI] != 0) continue;
	    
	    if(getDR(tempEta, tempPhi, geneta_[i][gI], genphi_[i][gI]) >= rParam[i]) continue;
		
	    isGood = false;
	    break;
	  }

	  if(!isGood) continue;
	  
	  etaVal_GenExclude.push_back(tempEta);
	  phiVal_GenExclude.push_back(tempPhi);
	  sumVal_GenExclude.push_back(0.0);
	  etaPoses_GenExclude.push_back({});
	  evtPlanePoses_GenExclude.push_back({0});

	  for(Int_t aI = 0; aI < nAbsEtaBins; ++aI){
	    if(TMath::Abs(tempEta) >= absEtaBinsLow[aI] && TMath::Abs(tempEta) < absEtaBinsHi[aI]){
	      etaPoses_GenExclude[etaPoses_GenExclude.size()-1].push_back(aI);
	    }
	  }

	  if(hiNevtPlane_ > 0){
            Double_t deltaPhi = TMath::Abs(getDPHI(tempPhi, hiEvtPlanes_[8]));
            if(deltaPhi < TMath::Pi()/4. || deltaPhi > 3.*TMath::Pi()/4.) evtPlanePoses_GenExclude[evtPlanePoses_GenExclude.size()-1].push_back(1);
            else evtPlanePoses_GenExclude[evtPlanePoses_GenExclude.size()-1].push_back(2);
          }
	}

	if(hasPF){
	  for(unsigned int pfI = 0; pfI < pfPt_p->size(); ++pfI){
	    for(unsigned int aI = 0; aI < etaVal.size(); ++aI){
	      if(getDR(pfEta_p->at(pfI), pfPhi_p->at(pfI), etaVal[aI], phiVal[aI]) < rParam[i]){
		sumVal[aI] += pfPt_p->at(pfI);
	      }
	    }

	    for(unsigned int aI = 0; aI < etaVal_GenExclude.size(); ++aI){
	      if(getDR(pfEta_p->at(pfI), pfPhi_p->at(pfI), etaVal_GenExclude[aI], phiVal_GenExclude[aI]) < rParam[i]){
		sumVal_GenExclude[aI] += pfPt_p->at(pfI);
	      }
	    }
	  }
	}

	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	Double_t width = 2.*rParam[i]/(Double_t)nRCSlices[i];
      	
	for(unsigned int aI = 0; aI < etaVal.size(); ++aI){
	  Double_t rcEstimate = 0.0;
	  Double_t rcEstimateFlow[nFlow];
	  for(int fI = 0; fI < nFlow; ++fI){
	    rcEstimateFlow[fI] = 0.0;
	  }

	  if(hasRho){
	    for(Int_t rcI = 0; rcI < nRCSlices[i]; ++rcI){
	      Float_t tempEtaVal1 = etaVal.at(aI) - rParam[i] + rcI*width;
	      Float_t tempEtaVal2 = etaVal.at(aI) - rParam[i] + (rcI+1)*width;
	      Float_t tempEtaVal = (tempEtaVal1 + tempEtaVal2)/2.;
	      Int_t etaPos = findBinPos(tempEtaVal, (*etaMin_p), (*etaMax_p));
	      
	      Float_t flowMod[nFlow];
	      for(int fI = 0; fI < nFlow; ++fI){
		flowMod[fI] = 1.0;
	      }

	      double val = ROOT::Math::chisquared_cdf_c(rhoFlowFitParams_p->at(5), rhoFlowFitParams_p->at(6));
	      bool minProb = val > .05;
	      bool maxProb = val < .95;

	      if(minProb && maxProb){
		Float_t tempFlowSum = 0;
		Float_t tempStartPhi = TMath::Sqrt(rParam[i]*rParam[i] - (etaVal.at(aI) - tempEtaVal)*(etaVal.at(aI) - tempEtaVal));
		
		Double_t nSamplesPhi = 25;
		Float_t phiWidth = 2.*tempStartPhi/nSamplesPhi;
		for(Int_t phiI = 0; phiI < nSamplesPhi; ++phiI){
		  Float_t currentPhi = phiVal.at(aI) - tempStartPhi + phiI*phiWidth;
		  
		  if(currentPhi < -2.*TMath::Pi()) currentPhi += 2.*TMath::Pi();
		  else if(currentPhi > 2.*TMath::Pi()) currentPhi -= 2.*TMath::Pi();
		  Float_t currModFactor = getModFactor(currentPhi, rhoFlowFitParams_p->at(2),  rhoFlowFitParams_p->at(4),  rhoFlowFitParams_p->at(1),  rhoFlowFitParams_p->at(3));
		  tempFlowSum += currModFactor;
		  
		}
		flowMod[0] = tempFlowSum/nSamplesPhi;	      
	      }

	      if(flowGlobal.size() > 0){
		for(unsigned int fI = 0; fI < flowGlobal.size(); ++fI){
		  double val = ROOT::Math::chisquared_cdf_c(rhoFlowFitParamsMore_p.at(fI)->at(5), rhoFlowFitParamsMore_p.at(fI)->at(6));
		  bool minProb = val > .05;
		  bool maxProb = val < .95;

		  if(minProb && maxProb){
		    Float_t tempFlowSum = 0;
		    Float_t tempStartPhi = TMath::Sqrt(rParam[i]*rParam[i] - (etaVal.at(aI) - tempEtaVal)*(etaVal.at(aI) - tempEtaVal));

		    Double_t nSamplesPhi = 25;
		    Float_t phiWidth = 2.*tempStartPhi/nSamplesPhi;
		    for(Int_t phiI = 0; phiI < nSamplesPhi; ++phiI){
		      Float_t currentPhi = phiVal.at(aI) - tempStartPhi + phiI*phiWidth;

		      if(currentPhi < -2.*TMath::Pi()) currentPhi += 2.*TMath::Pi();
		      else if(currentPhi > 2.*TMath::Pi()) currentPhi -= 2.*TMath::Pi();
		      Float_t currModFactor = getModFactor(currentPhi, rhoFlowFitParamsMore_p.at(fI)->at(2),  rhoFlowFitParamsMore_p.at(fI)->at(4), rhoFlowFitParamsMore_p.at(fI)->at(1),  rhoFlowFitParamsMore_p.at(fI)->at(3));
		      tempFlowSum += currModFactor;

		    }
		    flowMod[fI+1] = tempFlowSum/nSamplesPhi;
		  }

		}
	      }
	      
	      rcEstimate += rho_p->at(etaPos)*areaSlices.at(i).at(rcI);
	      for(int fI = 0; fI < nFlow; ++fI){
		rcEstimateFlow[fI] += rho_p->at(etaPos)*areaSlices.at(i).at(rcI)*flowMod[fI];
	      }
	    }
	  }
	  
	  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;		  
	  for(unsigned int eI = 0; eI < etaPoses[aI].size(); ++eI){
	    rcSum_VHiBin_h[(etaPoses[aI])[eI]][i]->Fill(hiBin_, sumVal[aI]);
	    pu_VHiBin_h[(etaPoses[aI])[eI]][i]->Fill(hiBin_, rcEstimate);
	    puFlow_VHiBin_h[(etaPoses[aI])[eI]][i]->Fill(hiBin_, rcEstimateFlow[0]);
	    rcSumMinPU_VHiBin_h[(etaPoses[aI])[eI]][i]->Fill(hiBin_, sumVal[aI] - rcEstimate);
	    rcSumMinPUFlow_VHiBin_h[(etaPoses[aI])[eI]][i]->Fill(hiBin_, sumVal[aI] - rcEstimateFlow[0]);

	    
	    rcSum_VHiBin_Points_h[(etaPoses[aI])[eI]][i][hiBin_]->Fill(sumVal[aI]);
	    pu_VHiBin_Points_h[(etaPoses[aI])[eI]][i][hiBin_]->Fill(rcEstimate);
	    puFlow_VHiBin_Points_h[(etaPoses[aI])[eI]][i][hiBin_]->Fill(rcEstimateFlow[0]);
	    rcSumMinPU_VHiBin_Points_h[(etaPoses[aI])[eI]][i][hiBin_]->Fill(sumVal[aI] - rcEstimate);
	    rcSumMinPUFlow_VHiBin_Points_h[(etaPoses[aI])[eI]][i][hiBin_]->Fill(sumVal[aI] - rcEstimateFlow[0]);
	    

	    rcSum_VHiBin_PointsRed_h[(etaPoses[aI])[eI]][i][hiBin_/5]->Fill(sumVal[aI]);
	    pu_VHiBin_PointsRed_h[(etaPoses[aI])[eI]][i][hiBin_/5]->Fill(rcEstimate);
	    puFlow_VHiBin_PointsRed_h[(etaPoses[aI])[eI]][i][hiBin_/5]->Fill(rcEstimateFlow[0]);
	    rcSumMinPU_VHiBin_PointsRed_h[(etaPoses[aI])[eI]][i][hiBin_/5]->Fill(sumVal[aI] - rcEstimate);
	    rcSumMinPUFlow_VHiBin_PointsRed_h[(etaPoses[aI])[eI]][i][hiBin_/5]->Fill(sumVal[aI] - rcEstimateFlow[0]);
	    
	    for(unsigned int pI = 0; pI < evtPlanePoses[aI].size(); ++pI){
	      rcSum_h[centPos][evtPlanePoses[aI].at(pI)][(etaPoses[aI])[eI]][i]->Fill(sumVal[aI]);
	      rcSumMinPU_h[centPos][evtPlanePoses[aI].at(pI)][(etaPoses[aI])[eI]][i]->Fill(sumVal[aI] - rcEstimate);
	      for(int fI = 0; fI < nFlow; ++fI){
		rcSumMinPUFlow_h[centPos][evtPlanePoses[aI].at(pI)][(etaPoses[aI])[eI]][i][fI]->Fill(sumVal[aI] - rcEstimateFlow[fI]);
	      }
	    }
	  }
	  etaPoses[aI].clear();
	  evtPlanePoses[aI].clear();
	}

	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	
	for(unsigned int aI = 0; aI < etaVal_GenExclude.size(); ++aI){
	  Double_t rcEstimate = 0.0;
          Double_t rcEstimateFlow[nFlow];
          for(int fI = 0; fI < nFlow; ++fI){
            rcEstimateFlow[fI] = 0.0;
          }

	  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	  if(hasRho){
	    for(Int_t rcI = 0; rcI < nRCSlices[i]; ++rcI){
	      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	      
	      Float_t tempEtaVal1 = etaVal_GenExclude.at(aI) - rParam[i] + rcI*width;
	      Float_t tempEtaVal2 = etaVal_GenExclude.at(aI) - rParam[i] + (rcI+1)*width;
	      Float_t tempEtaVal = (tempEtaVal1 + tempEtaVal2)/2.;
	      
	      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	      Int_t etaPos = findBinPos(tempEtaVal, (*etaMin_p), (*etaMax_p));
	      
	      Float_t flowMod[nFlow];
	      for(int fI = 0; fI < nFlow; ++fI){
		flowMod[fI] = 1.0;
	      }

	      if(rhoFlowFitParams_p->at(5) < 2.){
		Float_t tempFlowSum = 0;
		Float_t tempStartPhi = TMath::Sqrt(rParam[i]*rParam[i] - (etaVal_GenExclude.at(aI) - tempEtaVal)*(etaVal_GenExclude.at(aI) - tempEtaVal));
		Double_t nSamplesPhi = 25;
		Float_t phiWidth = 2.*tempStartPhi/nSamplesPhi;
		for(Int_t phiI = 0; phiI < nSamplesPhi; ++phiI){
		  Float_t currentPhi = phiVal_GenExclude.at(aI) - tempStartPhi + phiI*phiWidth;
		  if(currentPhi < -2.*TMath::Pi()) currentPhi += 2.*TMath::Pi();
		  else if(currentPhi > 2.*TMath::Pi()) currentPhi -= 2.*TMath::Pi();
		  tempFlowSum += getModFactor(currentPhi, rhoFlowFitParams_p->at(2),  rhoFlowFitParams_p->at(4),  rhoFlowFitParams_p->at(1),  rhoFlowFitParams_p->at(3));
		}
		flowMod[0] = tempFlowSum/nSamplesPhi;	      
	      }
	      

	      if(flowGlobal.size() > 0){
		for(unsigned int fI = 0; fI < flowGlobal.size(); ++fI){
		  double val = ROOT::Math::chisquared_cdf_c(rhoFlowFitParamsMore_p.at(fI)->at(5), rhoFlowFitParamsMore_p.at(fI)->at(6));
		  bool minProb = val > .05;
		  bool maxProb = val < .95;

		  if(minProb && maxProb){
		    Float_t tempFlowSum = 0;
		    Float_t tempStartPhi = TMath::Sqrt(rParam[i]*rParam[i] - (etaVal.at(aI) - tempEtaVal)*(etaVal.at(aI) - tempEtaVal));

		    Double_t nSamplesPhi = 25;
		    Float_t phiWidth = 2.*tempStartPhi/nSamplesPhi;
		    for(Int_t phiI = 0; phiI < nSamplesPhi; ++phiI){
		      Float_t currentPhi = phiVal.at(aI) - tempStartPhi + phiI*phiWidth;

		      if(currentPhi < -2.*TMath::Pi()) currentPhi += 2.*TMath::Pi();
		      else if(currentPhi > 2.*TMath::Pi()) currentPhi -= 2.*TMath::Pi();
		      Float_t currModFactor = getModFactor(currentPhi, rhoFlowFitParamsMore_p.at(fI)->at(2),  rhoFlowFitParamsMore_p.at(fI)->at(4), rhoFlowFitParamsMore_p.at(fI)->at(1),  rhoFlowFitParamsMore_p.at(fI)->at(3));
		      tempFlowSum += currModFactor;

		    }
		    flowMod[fI+1] = tempFlowSum/nSamplesPhi;
		  }

		}
	      }

	      rcEstimate += rho_p->at(etaPos)*areaSlices.at(i).at(rcI);
              for(int fI = 0; fI < nFlow; ++fI){
                rcEstimateFlow[fI] += rho_p->at(etaPos)*areaSlices.at(i).at(rcI)*flowMod[fI];
              }
	    }
	  }
	
	  //	  if(i == 0 && hiBin_ < 60) std::cout << "Entry: " << entry << std::endl;

	  for(unsigned int eI = 0; eI < etaPoses_GenExclude[aI].size(); ++eI){
	    rcSum_VHiBin_GenExclude_h[(etaPoses_GenExclude[aI])[eI]][i]->Fill(hiBin_, sumVal_GenExclude[aI]);
	    pu_VHiBin_GenExclude_h[(etaPoses_GenExclude[aI])[eI]][i]->Fill(hiBin_, rcEstimate);
	    puFlow_VHiBin_GenExclude_h[(etaPoses_GenExclude[aI])[eI]][i]->Fill(hiBin_, rcEstimateFlow[0]);
	    rcSumMinPU_VHiBin_GenExclude_h[(etaPoses_GenExclude[aI])[eI]][i]->Fill(hiBin_, sumVal_GenExclude[aI] - rcEstimate);
	    rcSumMinPUFlow_VHiBin_GenExclude_h[(etaPoses_GenExclude[aI])[eI]][i]->Fill(hiBin_, sumVal_GenExclude[aI] - rcEstimateFlow[0]);

	    rcSum_VHiBin_Points_GenExclude_h[(etaPoses[aI])[eI]][i][hiBin_]->Fill(sumVal_GenExclude[aI]);
	    pu_VHiBin_Points_GenExclude_h[(etaPoses[aI])[eI]][i][hiBin_]->Fill(rcEstimate);
	    puFlow_VHiBin_Points_GenExclude_h[(etaPoses[aI])[eI]][i][hiBin_]->Fill(rcEstimateFlow[0]);
	    rcSumMinPU_VHiBin_Points_GenExclude_h[(etaPoses[aI])[eI]][i][hiBin_]->Fill(sumVal_GenExclude[aI] - rcEstimate);
	    rcSumMinPUFlow_VHiBin_Points_GenExclude_h[(etaPoses[aI])[eI]][i][hiBin_]->Fill(sumVal_GenExclude[aI] - rcEstimateFlow[0]);

	    rcSum_VHiBin_PointsRed_GenExclude_h[(etaPoses[aI])[eI]][i][hiBin_]->Fill(sumVal_GenExclude[aI]);
	    pu_VHiBin_PointsRed_GenExclude_h[(etaPoses[aI])[eI]][i][hiBin_]->Fill(rcEstimate);
	    puFlow_VHiBin_PointsRed_GenExclude_h[(etaPoses[aI])[eI]][i][hiBin_]->Fill(rcEstimateFlow[0]);
	    rcSumMinPU_VHiBin_PointsRed_GenExclude_h[(etaPoses[aI])[eI]][i][hiBin_]->Fill(sumVal_GenExclude[aI] - rcEstimate);
	    rcSumMinPUFlow_VHiBin_PointsRed_GenExclude_h[(etaPoses[aI])[eI]][i][hiBin_]->Fill(sumVal_GenExclude[aI] - rcEstimateFlow[0]);
	    

	    for(unsigned int pI = 0; pI < evtPlanePoses_GenExclude[aI].size(); ++pI){
	      rcSum_GenExclude_h[centPos][evtPlanePoses_GenExclude[aI].at(pI)][(etaPoses_GenExclude[aI])[eI]][i]->Fill(sumVal_GenExclude[aI]);
	      rcSumMinPU_GenExclude_h[centPos][evtPlanePoses_GenExclude[aI].at(pI)][(etaPoses_GenExclude[aI])[eI]][i]->Fill(sumVal_GenExclude[aI] - rcEstimate);
	      for(int fI = 0; fI < nFlow; ++fI){
		rcSumMinPUFlow_GenExclude_h[centPos][evtPlanePoses_GenExclude[aI].at(pI)][(etaPoses_GenExclude[aI])[eI]][i][fI]->Fill(sumVal_GenExclude[aI] - rcEstimateFlow[fI]);
	      }
	    }
	  }
	  etaPoses_GenExclude[aI].clear();
	}

	etaVal.clear();
	phiVal.clear();
	sumVal.clear();
	etaPoses.clear();

	etaVal_GenExclude.clear();
	phiVal_GenExclude.clear();
	sumVal_GenExclude.clear();
	etaPoses_GenExclude.clear();
      }
    }
  
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    inFile_p->Close();
    delete inFile_p;
  }

  outFile_p->cd();

  for(Int_t i = 0; i < nAbsEtaBins; ++i){
    for(Int_t j = 0; j < nRParam; ++j){

      rcSum_VHiBin_h[i][j]->Write("", TObject::kOverwrite);
      pu_VHiBin_h[i][j]->Write("", TObject::kOverwrite);
      puFlow_VHiBin_h[i][j]->Write("", TObject::kOverwrite);
      rcSumMinPU_VHiBin_h[i][j]->Write("", TObject::kOverwrite);
      rcSumMinPUFlow_VHiBin_h[i][j]->Write("", TObject::kOverwrite);

      if(isMC && !isMB){
	rcSum_VHiBin_GenExclude_h[i][j]->Write("", TObject::kOverwrite);
	pu_VHiBin_GenExclude_h[i][j]->Write("", TObject::kOverwrite);
	puFlow_VHiBin_GenExclude_h[i][j]->Write("", TObject::kOverwrite);
	rcSumMinPU_VHiBin_GenExclude_h[i][j]->Write("", TObject::kOverwrite);
	rcSumMinPUFlow_VHiBin_GenExclude_h[i][j]->Write("", TObject::kOverwrite);
      }

      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      
      for(Int_t pI = 0; pI < 200; ++pI){
	rcSum_VHiBin_Mean_h[i][j]->SetBinContent(pI+1, rcSum_VHiBin_Points_h[i][j][pI]->GetMean());
	rcSum_VHiBin_Mean_h[i][j]->SetBinError(pI+1, rcSum_VHiBin_Points_h[i][j][pI]->GetMeanError());

	pu_VHiBin_Mean_h[i][j]->SetBinContent(pI+1, pu_VHiBin_Points_h[i][j][pI]->GetMean());
	pu_VHiBin_Mean_h[i][j]->SetBinError(pI+1, pu_VHiBin_Points_h[i][j][pI]->GetMeanError());

	puFlow_VHiBin_Mean_h[i][j]->SetBinContent(pI+1, puFlow_VHiBin_Points_h[i][j][pI]->GetMean());
	puFlow_VHiBin_Mean_h[i][j]->SetBinError(pI+1, puFlow_VHiBin_Points_h[i][j][pI]->GetMeanError());

	rcSumMinPU_VHiBin_Mean_h[i][j]->SetBinContent(pI+1, rcSumMinPU_VHiBin_Points_h[i][j][pI]->GetMean());
	rcSumMinPU_VHiBin_Mean_h[i][j]->SetBinError(pI+1, rcSumMinPU_VHiBin_Points_h[i][j][pI]->GetMeanError());

	rcSumMinPUFlow_VHiBin_Mean_h[i][j]->SetBinContent(pI+1, rcSumMinPUFlow_VHiBin_Points_h[i][j][pI]->GetMean());
	rcSumMinPUFlow_VHiBin_Mean_h[i][j]->SetBinError(pI+1, rcSumMinPUFlow_VHiBin_Points_h[i][j][pI]->GetMeanError());

	if(isMC && !isMB){
	  rcSum_VHiBin_Mean_GenExclude_h[i][j]->SetBinContent(pI+1, rcSum_VHiBin_Points_GenExclude_h[i][j][pI]->GetMean());
	  rcSum_VHiBin_Mean_GenExclude_h[i][j]->SetBinError(pI+1, rcSum_VHiBin_Points_GenExclude_h[i][j][pI]->GetMeanError());
	  
	  pu_VHiBin_Mean_GenExclude_h[i][j]->SetBinContent(pI+1, pu_VHiBin_Points_GenExclude_h[i][j][pI]->GetMean());
	  pu_VHiBin_Mean_GenExclude_h[i][j]->SetBinError(pI+1, pu_VHiBin_Points_GenExclude_h[i][j][pI]->GetMeanError());
	  
	  puFlow_VHiBin_Mean_GenExclude_h[i][j]->SetBinContent(pI+1, puFlow_VHiBin_Points_GenExclude_h[i][j][pI]->GetMean());
	  puFlow_VHiBin_Mean_GenExclude_h[i][j]->SetBinError(pI+1, puFlow_VHiBin_Points_GenExclude_h[i][j][pI]->GetMeanError());
	  
	  rcSumMinPU_VHiBin_Mean_GenExclude_h[i][j]->SetBinContent(pI+1, rcSumMinPU_VHiBin_Points_GenExclude_h[i][j][pI]->GetMean());
	  rcSumMinPU_VHiBin_Mean_GenExclude_h[i][j]->SetBinError(pI+1, rcSumMinPU_VHiBin_Points_GenExclude_h[i][j][pI]->GetMeanError());
	  
	  rcSumMinPUFlow_VHiBin_Mean_GenExclude_h[i][j]->SetBinContent(pI+1, rcSumMinPUFlow_VHiBin_Points_GenExclude_h[i][j][pI]->GetMean());
	  rcSumMinPUFlow_VHiBin_Mean_GenExclude_h[i][j]->SetBinError(pI+1, rcSumMinPUFlow_VHiBin_Points_GenExclude_h[i][j][pI]->GetMeanError());
	}
      
	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	rcSum_VHiBin_Sigma_h[i][j]->SetBinContent(pI+1, rcSum_VHiBin_Points_h[i][j][pI]->GetStdDev());
	rcSum_VHiBin_Sigma_h[i][j]->SetBinError(pI+1, rcSum_VHiBin_Points_h[i][j][pI]->GetStdDevError());

	pu_VHiBin_Sigma_h[i][j]->SetBinContent(pI+1, pu_VHiBin_Points_h[i][j][pI]->GetStdDev());
	pu_VHiBin_Sigma_h[i][j]->SetBinError(pI+1, pu_VHiBin_Points_h[i][j][pI]->GetStdDevError());

	puFlow_VHiBin_Sigma_h[i][j]->SetBinContent(pI+1, puFlow_VHiBin_Points_h[i][j][pI]->GetStdDev());
	puFlow_VHiBin_Sigma_h[i][j]->SetBinError(pI+1, puFlow_VHiBin_Points_h[i][j][pI]->GetStdDevError());

	rcSumMinPU_VHiBin_Sigma_h[i][j]->SetBinContent(pI+1, rcSumMinPU_VHiBin_Points_h[i][j][pI]->GetStdDev());
	rcSumMinPU_VHiBin_Sigma_h[i][j]->SetBinError(pI+1, rcSumMinPU_VHiBin_Points_h[i][j][pI]->GetStdDevError());

	rcSumMinPUFlow_VHiBin_Sigma_h[i][j]->SetBinContent(pI+1, rcSumMinPUFlow_VHiBin_Points_h[i][j][pI]->GetStdDev());
	rcSumMinPUFlow_VHiBin_Sigma_h[i][j]->SetBinError(pI+1, rcSumMinPUFlow_VHiBin_Points_h[i][j][pI]->GetStdDevError());


	if(isMC && !isMB){
	  rcSum_VHiBin_Sigma_GenExclude_h[i][j]->SetBinContent(pI+1, rcSum_VHiBin_Points_GenExclude_h[i][j][pI]->GetStdDev());
	  rcSum_VHiBin_Sigma_GenExclude_h[i][j]->SetBinError(pI+1, rcSum_VHiBin_Points_GenExclude_h[i][j][pI]->GetStdDevError());
	  
	  pu_VHiBin_Sigma_GenExclude_h[i][j]->SetBinContent(pI+1, pu_VHiBin_Points_GenExclude_h[i][j][pI]->GetStdDev());
	  pu_VHiBin_Sigma_GenExclude_h[i][j]->SetBinError(pI+1, pu_VHiBin_Points_GenExclude_h[i][j][pI]->GetStdDevError());

	  puFlow_VHiBin_Sigma_GenExclude_h[i][j]->SetBinContent(pI+1, puFlow_VHiBin_Points_GenExclude_h[i][j][pI]->GetStdDev());
	  puFlow_VHiBin_Sigma_GenExclude_h[i][j]->SetBinError(pI+1, puFlow_VHiBin_Points_GenExclude_h[i][j][pI]->GetStdDevError());
	  
	  rcSumMinPU_VHiBin_Sigma_GenExclude_h[i][j]->SetBinContent(pI+1, rcSumMinPU_VHiBin_Points_GenExclude_h[i][j][pI]->GetStdDev());
	  rcSumMinPU_VHiBin_Sigma_GenExclude_h[i][j]->SetBinError(pI+1, rcSumMinPU_VHiBin_Points_GenExclude_h[i][j][pI]->GetStdDevError());
	  
	  rcSumMinPUFlow_VHiBin_Sigma_GenExclude_h[i][j]->SetBinContent(pI+1, rcSumMinPUFlow_VHiBin_Points_GenExclude_h[i][j][pI]->GetStdDev());
	  rcSumMinPUFlow_VHiBin_Sigma_GenExclude_h[i][j]->SetBinError(pI+1, rcSumMinPUFlow_VHiBin_Points_GenExclude_h[i][j][pI]->GetStdDevError());
	}

	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	if(rcSum_VHiBin_Points_h[i][j][pI]->GetMean() > 0){
	  rcSum_VHiBin_SigmaOverMean_h[i][j]->SetBinContent(pI+1, rcSum_VHiBin_Points_h[i][j][pI]->GetStdDev()/rcSum_VHiBin_Points_h[i][j][pI]->GetMean());
	  rcSum_VHiBin_SigmaOverMean_h[i][j]->SetBinError(pI+1, rcSum_VHiBin_Points_h[i][j][pI]->GetStdDevError()/rcSum_VHiBin_Points_h[i][j][pI]->GetMean());
	}
	else{
	  rcSum_VHiBin_SigmaOverMean_h[i][j]->SetBinContent(pI+1, 0);
	  rcSum_VHiBin_SigmaOverMean_h[i][j]->SetBinError(pI+1, 0);
	}

	if(pu_VHiBin_Points_h[i][j][pI]->GetMean() > 0){
	  pu_VHiBin_SigmaOverMean_h[i][j]->SetBinContent(pI+1, pu_VHiBin_Points_h[i][j][pI]->GetStdDev()/pu_VHiBin_Points_h[i][j][pI]->GetMean());
	  pu_VHiBin_SigmaOverMean_h[i][j]->SetBinError(pI+1, pu_VHiBin_Points_h[i][j][pI]->GetStdDevError()/pu_VHiBin_Points_h[i][j][pI]->GetMean());
	}
	else{
	  pu_VHiBin_SigmaOverMean_h[i][j]->SetBinContent(pI+1, 0);
	  pu_VHiBin_SigmaOverMean_h[i][j]->SetBinError(pI+1, 0);
	}

	if(puFlow_VHiBin_Points_h[i][j][pI]->GetMean() > 0){
	  puFlow_VHiBin_SigmaOverMean_h[i][j]->SetBinContent(pI+1, puFlow_VHiBin_Points_h[i][j][pI]->GetStdDev()/puFlow_VHiBin_Points_h[i][j][pI]->GetMean());
	  puFlow_VHiBin_SigmaOverMean_h[i][j]->SetBinError(pI+1, puFlow_VHiBin_Points_h[i][j][pI]->GetStdDevError()/puFlow_VHiBin_Points_h[i][j][pI]->GetMean());
	}
	else{
	  puFlow_VHiBin_SigmaOverMean_h[i][j]->SetBinContent(pI+1, 0);
	  puFlow_VHiBin_SigmaOverMean_h[i][j]->SetBinError(pI+1, 0);	  
	}

	//	rcSumMinPU_VHiBin_Points_h[i][j][pI]->Write("", TObject::kOverwrite);
	//	rcSumMinPUFlow_VHiBin_Points_h[i][j][pI]->Write("", TObject::kOverwrite);
      }
      

      for(Int_t pI = 0; pI < nRedCentBins; ++pI){
	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	rcSum_VHiBin_MeanRed_h[i][j]->SetBinContent(pI+1, rcSum_VHiBin_PointsRed_h[i][j][pI]->GetMean());
	rcSum_VHiBin_MeanRed_h[i][j]->SetBinError(pI+1, rcSum_VHiBin_PointsRed_h[i][j][pI]->GetMeanError());

	pu_VHiBin_MeanRed_h[i][j]->SetBinContent(pI+1, pu_VHiBin_PointsRed_h[i][j][pI]->GetMean());
	pu_VHiBin_MeanRed_h[i][j]->SetBinError(pI+1, pu_VHiBin_PointsRed_h[i][j][pI]->GetMeanError());

	puFlow_VHiBin_MeanRed_h[i][j]->SetBinContent(pI+1, puFlow_VHiBin_PointsRed_h[i][j][pI]->GetMean());
	puFlow_VHiBin_MeanRed_h[i][j]->SetBinError(pI+1, puFlow_VHiBin_PointsRed_h[i][j][pI]->GetMeanError());

	rcSumMinPU_VHiBin_MeanRed_h[i][j]->SetBinContent(pI+1, rcSumMinPU_VHiBin_PointsRed_h[i][j][pI]->GetMean());
	rcSumMinPU_VHiBin_MeanRed_h[i][j]->SetBinError(pI+1, rcSumMinPU_VHiBin_PointsRed_h[i][j][pI]->GetMeanError());

	rcSumMinPUFlow_VHiBin_MeanRed_h[i][j]->SetBinContent(pI+1, rcSumMinPUFlow_VHiBin_PointsRed_h[i][j][pI]->GetMean());
	rcSumMinPUFlow_VHiBin_MeanRed_h[i][j]->SetBinError(pI+1, rcSumMinPUFlow_VHiBin_PointsRed_h[i][j][pI]->GetMeanError());

	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	if(isMC && !isMB){
	  rcSum_VHiBin_MeanRed_GenExclude_h[i][j]->SetBinContent(pI+1, rcSum_VHiBin_PointsRed_GenExclude_h[i][j][pI]->GetMean());
	  rcSum_VHiBin_MeanRed_GenExclude_h[i][j]->SetBinError(pI+1, rcSum_VHiBin_PointsRed_GenExclude_h[i][j][pI]->GetMeanError());
	  
	  pu_VHiBin_MeanRed_GenExclude_h[i][j]->SetBinContent(pI+1, pu_VHiBin_PointsRed_GenExclude_h[i][j][pI]->GetMean());
	  pu_VHiBin_MeanRed_GenExclude_h[i][j]->SetBinError(pI+1, pu_VHiBin_PointsRed_GenExclude_h[i][j][pI]->GetMeanError());
	  
	  puFlow_VHiBin_MeanRed_GenExclude_h[i][j]->SetBinContent(pI+1, puFlow_VHiBin_PointsRed_GenExclude_h[i][j][pI]->GetMean());
	  puFlow_VHiBin_MeanRed_GenExclude_h[i][j]->SetBinError(pI+1, puFlow_VHiBin_PointsRed_GenExclude_h[i][j][pI]->GetMeanError());
	  
	  rcSumMinPU_VHiBin_MeanRed_GenExclude_h[i][j]->SetBinContent(pI+1, rcSumMinPU_VHiBin_PointsRed_GenExclude_h[i][j][pI]->GetMean());
	  rcSumMinPU_VHiBin_MeanRed_GenExclude_h[i][j]->SetBinError(pI+1, rcSumMinPU_VHiBin_PointsRed_GenExclude_h[i][j][pI]->GetMeanError());
	  
	  rcSumMinPUFlow_VHiBin_MeanRed_GenExclude_h[i][j]->SetBinContent(pI+1, rcSumMinPUFlow_VHiBin_PointsRed_GenExclude_h[i][j][pI]->GetMean());
	  rcSumMinPUFlow_VHiBin_MeanRed_GenExclude_h[i][j]->SetBinError(pI+1, rcSumMinPUFlow_VHiBin_PointsRed_GenExclude_h[i][j][pI]->GetMeanError());
	}

	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	rcSum_VHiBin_SigmaRed_h[i][j]->SetBinContent(pI+1, rcSum_VHiBin_PointsRed_h[i][j][pI]->GetStdDev());
	rcSum_VHiBin_SigmaRed_h[i][j]->SetBinError(pI+1, rcSum_VHiBin_PointsRed_h[i][j][pI]->GetStdDevError());

	pu_VHiBin_SigmaRed_h[i][j]->SetBinContent(pI+1, pu_VHiBin_PointsRed_h[i][j][pI]->GetStdDev());
	pu_VHiBin_SigmaRed_h[i][j]->SetBinError(pI+1, pu_VHiBin_PointsRed_h[i][j][pI]->GetStdDevError());

	puFlow_VHiBin_SigmaRed_h[i][j]->SetBinContent(pI+1, puFlow_VHiBin_PointsRed_h[i][j][pI]->GetStdDev());
	puFlow_VHiBin_SigmaRed_h[i][j]->SetBinError(pI+1, puFlow_VHiBin_PointsRed_h[i][j][pI]->GetStdDevError());

	rcSumMinPU_VHiBin_SigmaRed_h[i][j]->SetBinContent(pI+1, rcSumMinPU_VHiBin_PointsRed_h[i][j][pI]->GetStdDev());
	rcSumMinPU_VHiBin_SigmaRed_h[i][j]->SetBinError(pI+1, rcSumMinPU_VHiBin_PointsRed_h[i][j][pI]->GetStdDevError());

	rcSumMinPUFlow_VHiBin_SigmaRed_h[i][j]->SetBinContent(pI+1, rcSumMinPUFlow_VHiBin_PointsRed_h[i][j][pI]->GetStdDev());
	rcSumMinPUFlow_VHiBin_SigmaRed_h[i][j]->SetBinError(pI+1, rcSumMinPUFlow_VHiBin_PointsRed_h[i][j][pI]->GetStdDevError());


	if(isMC && !isMB){
	  rcSum_VHiBin_SigmaRed_GenExclude_h[i][j]->SetBinContent(pI+1, rcSum_VHiBin_PointsRed_GenExclude_h[i][j][pI]->GetStdDev());
	  rcSum_VHiBin_SigmaRed_GenExclude_h[i][j]->SetBinError(pI+1, rcSum_VHiBin_PointsRed_GenExclude_h[i][j][pI]->GetStdDevError());
	  
	  pu_VHiBin_SigmaRed_GenExclude_h[i][j]->SetBinContent(pI+1, pu_VHiBin_PointsRed_GenExclude_h[i][j][pI]->GetStdDev());
	  pu_VHiBin_SigmaRed_GenExclude_h[i][j]->SetBinError(pI+1, pu_VHiBin_PointsRed_GenExclude_h[i][j][pI]->GetStdDevError());

	  puFlow_VHiBin_SigmaRed_GenExclude_h[i][j]->SetBinContent(pI+1, puFlow_VHiBin_PointsRed_GenExclude_h[i][j][pI]->GetStdDev());
	  puFlow_VHiBin_SigmaRed_GenExclude_h[i][j]->SetBinError(pI+1, puFlow_VHiBin_PointsRed_GenExclude_h[i][j][pI]->GetStdDevError());
	  
	  rcSumMinPU_VHiBin_SigmaRed_GenExclude_h[i][j]->SetBinContent(pI+1, rcSumMinPU_VHiBin_PointsRed_GenExclude_h[i][j][pI]->GetStdDev());
	  rcSumMinPU_VHiBin_SigmaRed_GenExclude_h[i][j]->SetBinError(pI+1, rcSumMinPU_VHiBin_PointsRed_GenExclude_h[i][j][pI]->GetStdDevError());
	  
	  rcSumMinPUFlow_VHiBin_SigmaRed_GenExclude_h[i][j]->SetBinContent(pI+1, rcSumMinPUFlow_VHiBin_PointsRed_GenExclude_h[i][j][pI]->GetStdDev());
	  rcSumMinPUFlow_VHiBin_SigmaRed_GenExclude_h[i][j]->SetBinError(pI+1, rcSumMinPUFlow_VHiBin_PointsRed_GenExclude_h[i][j][pI]->GetStdDevError());
	}

	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
      
	rcSum_VHiBin_PointsRed_h[i][j][pI]->Write("", TObject::kOverwrite);;
	pu_VHiBin_PointsRed_h[i][j][pI]->Write("", TObject::kOverwrite);;
	puFlow_VHiBin_PointsRed_h[i][j][pI]->Write("", TObject::kOverwrite);;
	rcSumMinPU_VHiBin_PointsRed_h[i][j][pI]->Write("", TObject::kOverwrite);;
	rcSumMinPUFlow_VHiBin_PointsRed_h[i][j][pI]->Write("", TObject::kOverwrite);;

	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	if(isMC && !isMB){
	  rcSum_VHiBin_PointsRed_GenExclude_h[i][j][pI]->Write("", TObject::kOverwrite);;
	  pu_VHiBin_PointsRed_GenExclude_h[i][j][pI]->Write("", TObject::kOverwrite);;
	  puFlow_VHiBin_PointsRed_GenExclude_h[i][j][pI]->Write("", TObject::kOverwrite);;
	  rcSumMinPU_VHiBin_PointsRed_GenExclude_h[i][j][pI]->Write("", TObject::kOverwrite);;
	  rcSumMinPUFlow_VHiBin_PointsRed_GenExclude_h[i][j][pI]->Write("", TObject::kOverwrite);
	}
	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
      }

      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      rcSum_VHiBin_Mean_h[i][j]->Write("", TObject::kOverwrite);
      rcSum_VHiBin_Sigma_h[i][j]->Write("", TObject::kOverwrite);
      rcSum_VHiBin_SigmaOverMean_h[i][j]->Write("", TObject::kOverwrite);

      pu_VHiBin_Mean_h[i][j]->Write("", TObject::kOverwrite);
      pu_VHiBin_Sigma_h[i][j]->Write("", TObject::kOverwrite);
      pu_VHiBin_SigmaOverMean_h[i][j]->Write("", TObject::kOverwrite);

      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      puFlow_VHiBin_Mean_h[i][j]->Write("", TObject::kOverwrite);
      puFlow_VHiBin_Sigma_h[i][j]->Write("", TObject::kOverwrite);
      puFlow_VHiBin_SigmaOverMean_h[i][j]->Write("", TObject::kOverwrite);

      rcSumMinPU_VHiBin_Mean_h[i][j]->Write("", TObject::kOverwrite);
      rcSumMinPU_VHiBin_Sigma_h[i][j]->Write("", TObject::kOverwrite);

      rcSumMinPUFlow_VHiBin_Mean_h[i][j]->Write("", TObject::kOverwrite);
      rcSumMinPUFlow_VHiBin_Sigma_h[i][j]->Write("", TObject::kOverwrite);


      rcSum_VHiBin_MeanRed_h[i][j]->Write("", TObject::kOverwrite);
      rcSum_VHiBin_SigmaRed_h[i][j]->Write("", TObject::kOverwrite);

      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      pu_VHiBin_MeanRed_h[i][j]->Write("", TObject::kOverwrite);
      pu_VHiBin_SigmaRed_h[i][j]->Write("", TObject::kOverwrite);

      puFlow_VHiBin_MeanRed_h[i][j]->Write("", TObject::kOverwrite);
      puFlow_VHiBin_SigmaRed_h[i][j]->Write("", TObject::kOverwrite);

      rcSumMinPU_VHiBin_MeanRed_h[i][j]->Write("", TObject::kOverwrite);
      rcSumMinPU_VHiBin_SigmaRed_h[i][j]->Write("", TObject::kOverwrite);

      rcSumMinPUFlow_VHiBin_MeanRed_h[i][j]->Write("", TObject::kOverwrite);
      rcSumMinPUFlow_VHiBin_SigmaRed_h[i][j]->Write("", TObject::kOverwrite);

      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      if(isMC && !isMB){
	rcSum_VHiBin_Mean_GenExclude_h[i][j]->Write("", TObject::kOverwrite);
	rcSum_VHiBin_Sigma_GenExclude_h[i][j]->Write("", TObject::kOverwrite);
	
	pu_VHiBin_Mean_GenExclude_h[i][j]->Write("", TObject::kOverwrite);
	pu_VHiBin_Sigma_GenExclude_h[i][j]->Write("", TObject::kOverwrite);
	
	puFlow_VHiBin_Mean_GenExclude_h[i][j]->Write("", TObject::kOverwrite);
	puFlow_VHiBin_Sigma_GenExclude_h[i][j]->Write("", TObject::kOverwrite);
	
	rcSumMinPU_VHiBin_Mean_GenExclude_h[i][j]->Write("", TObject::kOverwrite);
	rcSumMinPU_VHiBin_Sigma_GenExclude_h[i][j]->Write("", TObject::kOverwrite);
	
	rcSumMinPUFlow_VHiBin_Mean_GenExclude_h[i][j]->Write("", TObject::kOverwrite);
	rcSumMinPUFlow_VHiBin_Sigma_GenExclude_h[i][j]->Write("", TObject::kOverwrite);


	rcSum_VHiBin_MeanRed_GenExclude_h[i][j]->Write("", TObject::kOverwrite);
	rcSum_VHiBin_SigmaRed_GenExclude_h[i][j]->Write("", TObject::kOverwrite);
	
	pu_VHiBin_MeanRed_GenExclude_h[i][j]->Write("", TObject::kOverwrite);
	pu_VHiBin_SigmaRed_GenExclude_h[i][j]->Write("", TObject::kOverwrite);
	
	puFlow_VHiBin_MeanRed_GenExclude_h[i][j]->Write("", TObject::kOverwrite);
	puFlow_VHiBin_SigmaRed_GenExclude_h[i][j]->Write("", TObject::kOverwrite);
	
	rcSumMinPU_VHiBin_MeanRed_GenExclude_h[i][j]->Write("", TObject::kOverwrite);
	rcSumMinPU_VHiBin_SigmaRed_GenExclude_h[i][j]->Write("", TObject::kOverwrite);
	
	rcSumMinPUFlow_VHiBin_MeanRed_GenExclude_h[i][j]->Write("", TObject::kOverwrite);
	rcSumMinPUFlow_VHiBin_SigmaRed_GenExclude_h[i][j]->Write("", TObject::kOverwrite);
      }

      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      for(Int_t cI = 0; cI < nCentBins; ++cI){

	for(Int_t eI = 0; eI < nEvtPlaneBins; ++eI){
	  rcSum_h[cI][eI][i][j]->Write("", TObject::kOverwrite);
	  if(isMC && !isMB) rcSum_GenExclude_h[cI][eI][i][j]->Write("", TObject::kOverwrite);
	  
	  rcSumMinPU_h[cI][eI][i][j]->Write("", TObject::kOverwrite);
	  if(isMC && !isMB) rcSumMinPU_GenExclude_h[cI][eI][i][j]->Write("", TObject::kOverwrite);

	  for(Int_t fI = 0; fI < nFlow; ++fI){
	    rcSumMinPUFlow_h[cI][eI][i][j][fI]->Write("", TObject::kOverwrite);
	    if(isMC && !isMB) rcSumMinPUFlow_GenExclude_h[cI][eI][i][j][fI]->Write("", TObject::kOverwrite);
	  }
	}
      }
    }
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    for(Int_t pI = 0; pI < nPFID+1; ++pI){
      for(Int_t bI = 0; bI < nEtaTowBoundsHist; ++bI){

	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << ", " << cI << ", " << pI << ", " << bI << std::endl;

	pfCandMap_Mean_h[cI][pI]->SetBinContent(bI+1, pfCandMap_Points_h[cI][pI][bI]->GetMean());
	pfCandMap_Mean_h[cI][pI]->SetBinError(bI+1, pfCandMap_Points_h[cI][pI][bI]->GetMeanError());

	pfCandMap_Sigma_h[cI][pI]->SetBinContent(bI+1, pfCandMap_Points_h[cI][pI][bI]->GetStdDev());
	pfCandMap_Sigma_h[cI][pI]->SetBinError(bI+1, pfCandMap_Points_h[cI][pI][bI]->GetStdDevError());

	pfCandMap_Points_h[cI][pI][bI]->Write("", TObject::kOverwrite);
	delete pfCandMap_Points_h[cI][pI][bI];

	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << ", " << cI << ", " << pI << ", " << bI << std::endl;

	if(hasPFCs){
	  pfCsCandMap_Mean_h[cI][pI]->SetBinContent(bI+1, pfCsCandMap_Points_h[cI][pI][bI]->GetMean());
	  pfCsCandMap_Mean_h[cI][pI]->SetBinError(bI+1, pfCsCandMap_Points_h[cI][pI][bI]->GetMeanError());
	  
	  
	  pfCsCandMap_Sigma_h[cI][pI]->SetBinContent(bI+1, pfCsCandMap_Points_h[cI][pI][bI]->GetStdDev());
	  pfCsCandMap_Sigma_h[cI][pI]->SetBinError(bI+1, pfCsCandMap_Points_h[cI][pI][bI]->GetStdDevError());
	  
	  pfCsCandMap_Points_h[cI][pI][bI]->Write("", TObject::kOverwrite);
	}

	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << ", " << cI << ", " << pI << ", " << bI << std::endl;

	delete pfCsCandMap_Points_h[cI][pI][bI];
      }


      pfCandMap_Mean_h[cI][pI]->Write("", TObject::kOverwrite);
      delete pfCandMap_Mean_h[cI][pI];

      pfCandMap_Sigma_h[cI][pI]->Write("", TObject::kOverwrite);
      delete pfCandMap_Sigma_h[cI][pI];

      pfCsCandMap_Mean_h[cI][pI]->Write("", TObject::kOverwrite);
      delete pfCsCandMap_Mean_h[cI][pI];

      pfCsCandMap_Sigma_h[cI][pI]->Write("", TObject::kOverwrite);
      delete pfCsCandMap_Sigma_h[cI][pI];
    }
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;


  for(Int_t i = 0; i < nAbsEtaBins; ++i){
    for(Int_t j = 0; j < nRParam; ++j){

      delete rcSum_VHiBin_h[i][j];
      delete pu_VHiBin_h[i][j];
      delete puFlow_VHiBin_h[i][j];
      delete rcSumMinPU_VHiBin_h[i][j];
      delete rcSumMinPUFlow_VHiBin_h[i][j];

      if(isMC && !isMB){
	delete rcSum_VHiBin_GenExclude_h[i][j];
	delete pu_VHiBin_GenExclude_h[i][j];
	delete puFlow_VHiBin_GenExclude_h[i][j];
	delete rcSumMinPU_VHiBin_GenExclude_h[i][j];
	delete rcSumMinPUFlow_VHiBin_GenExclude_h[i][j];
      }

      delete rcSum_VHiBin_Mean_h[i][j];
      delete pu_VHiBin_Mean_h[i][j];
      delete puFlow_VHiBin_Mean_h[i][j];
      delete rcSumMinPU_VHiBin_Mean_h[i][j];
      delete rcSumMinPUFlow_VHiBin_Mean_h[i][j];

      delete rcSum_VHiBin_Sigma_h[i][j];
      delete pu_VHiBin_Sigma_h[i][j];
      delete puFlow_VHiBin_Sigma_h[i][j];
      delete rcSumMinPU_VHiBin_Sigma_h[i][j];
      delete rcSumMinPUFlow_VHiBin_Sigma_h[i][j];

      if(isMC && !isMB){
	delete rcSum_VHiBin_Mean_GenExclude_h[i][j];
	delete pu_VHiBin_Mean_GenExclude_h[i][j];
	delete puFlow_VHiBin_Mean_GenExclude_h[i][j];
	delete rcSumMinPU_VHiBin_Mean_GenExclude_h[i][j];
	delete rcSumMinPUFlow_VHiBin_Mean_GenExclude_h[i][j];
      
	delete rcSum_VHiBin_Sigma_GenExclude_h[i][j];
	delete pu_VHiBin_Sigma_GenExclude_h[i][j];
	delete puFlow_VHiBin_Sigma_GenExclude_h[i][j];
	delete rcSumMinPU_VHiBin_Sigma_GenExclude_h[i][j];
	delete rcSumMinPUFlow_VHiBin_Sigma_GenExclude_h[i][j];
      }

      delete rcSum_VHiBin_MeanRed_h[i][j];
      delete pu_VHiBin_MeanRed_h[i][j];
      delete puFlow_VHiBin_MeanRed_h[i][j];
      delete rcSumMinPU_VHiBin_MeanRed_h[i][j];
      delete rcSumMinPUFlow_VHiBin_MeanRed_h[i][j];

      delete rcSum_VHiBin_SigmaRed_h[i][j];
      delete pu_VHiBin_SigmaRed_h[i][j];
      delete puFlow_VHiBin_SigmaRed_h[i][j];
      delete rcSumMinPU_VHiBin_SigmaRed_h[i][j];
      delete rcSumMinPUFlow_VHiBin_SigmaRed_h[i][j];

      if(isMC && !isMB){
	delete rcSum_VHiBin_MeanRed_GenExclude_h[i][j];
	delete pu_VHiBin_MeanRed_GenExclude_h[i][j];
	delete puFlow_VHiBin_MeanRed_GenExclude_h[i][j];
	delete rcSumMinPU_VHiBin_MeanRed_GenExclude_h[i][j];
	delete rcSumMinPUFlow_VHiBin_MeanRed_GenExclude_h[i][j];
      
	delete rcSum_VHiBin_SigmaRed_GenExclude_h[i][j];
	delete pu_VHiBin_SigmaRed_GenExclude_h[i][j];
	delete puFlow_VHiBin_SigmaRed_GenExclude_h[i][j];
	delete rcSumMinPU_VHiBin_SigmaRed_GenExclude_h[i][j];
	delete rcSumMinPUFlow_VHiBin_SigmaRed_GenExclude_h[i][j];
      }

      delete rcSum_VHiBin_SigmaOverMean_h[i][j];
      delete pu_VHiBin_SigmaOverMean_h[i][j];
      delete puFlow_VHiBin_SigmaOverMean_h[i][j];

      
      for(Int_t pI = 0; pI < 200; ++pI){
	delete rcSum_VHiBin_Points_h[i][j][pI];
	delete pu_VHiBin_Points_h[i][j][pI];
	delete puFlow_VHiBin_Points_h[i][j][pI];
	delete rcSumMinPU_VHiBin_Points_h[i][j][pI];
	delete rcSumMinPUFlow_VHiBin_Points_h[i][j][pI];

	if(isMC && !isMB){
	  delete rcSum_VHiBin_Points_GenExclude_h[i][j][pI];
	  delete pu_VHiBin_Points_GenExclude_h[i][j][pI];
	  delete puFlow_VHiBin_Points_GenExclude_h[i][j][pI];
	  delete rcSumMinPU_VHiBin_Points_GenExclude_h[i][j][pI];
	  delete rcSumMinPUFlow_VHiBin_Points_GenExclude_h[i][j][pI];
	}
      }
     

      for(Int_t pI = 0; pI < nRedCentBins; ++pI){
	delete rcSum_VHiBin_PointsRed_h[i][j][pI];
	delete pu_VHiBin_PointsRed_h[i][j][pI];
	delete puFlow_VHiBin_PointsRed_h[i][j][pI];
	delete rcSumMinPU_VHiBin_PointsRed_h[i][j][pI];
	delete rcSumMinPUFlow_VHiBin_PointsRed_h[i][j][pI];

	if(isMC && !isMB){
	  delete rcSum_VHiBin_PointsRed_GenExclude_h[i][j][pI];
	  delete pu_VHiBin_PointsRed_GenExclude_h[i][j][pI];
	  delete puFlow_VHiBin_PointsRed_GenExclude_h[i][j][pI];
	  delete rcSumMinPU_VHiBin_PointsRed_GenExclude_h[i][j][pI];
	  delete rcSumMinPUFlow_VHiBin_PointsRed_GenExclude_h[i][j][pI];
	}
      }

      for(Int_t cI = 0; cI < nCentBins; ++cI){
	for(Int_t eI = 0; eI < nEvtPlaneBins; ++eI){
	  delete rcSum_h[cI][eI][i][j];
	  if(isMC && !isMB) delete rcSum_GenExclude_h[cI][eI][i][j];
	  delete rcSumMinPU_h[cI][eI][i][j];
	  if(isMC && !isMB) delete rcSumMinPU_GenExclude_h[cI][eI][i][j];

	  for(Int_t fI = 0; fI < nFlow; ++fI){
	    delete rcSumMinPUFlow_h[cI][eI][i][j][fI];
	    if(isMC && !isMB) delete rcSumMinPUFlow_GenExclude_h[cI][eI][i][j][fI];
	  }
	}
      }
    }
  }

      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  outFile_p->Close();
  delete outFile_p;

  delete date;
  delete randGen_p;

  std::cout << "Job complete. Return 0" << std::endl;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2 && argc != 3 && argc != 4){
    std::cout << "Usage: ./rc.exe <inFileName> <isMC-optional> <doFlowExamples-optional>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 2) retVal += rc(argv[1]);
  else if(argc == 3) retVal += rc(argv[1], std::stoi(argv[2]));
  else if(argc == 4) retVal += rc(argv[1], std::stoi(argv[2]), std::stoi(argv[3]));
  return retVal;
}
