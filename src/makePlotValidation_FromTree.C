#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <thread>

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

//const Bool_t doSvd = false;

const Int_t nPadTotValid = 13;

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


void comboValidation(TH1* histTot_p, TH1* comboTot_p)
{
  TDatime* date = new TDatime();
  const Int_t nPadTot = 4;
  const Int_t nPadX = 2;
  TPad* pads[nPadTot];
  TCanvas* validation_p = new TCanvas("temp", "temp", nPadX*400, 400);
  validation_p->SetTopMargin(0.01);
  validation_p->SetBottomMargin(0.01);
  validation_p->SetLeftMargin(0.01);
  validation_p->SetRightMargin(0.01);

  for(Int_t i = 0; i < nPadTot; ++i){
    validation_p->cd();

    if(i == 0) pads[i] = new TPad(("pad" + std::to_string(i)).c_str(), "", 0.0, 0.35, 0.5, 1.0);
    else if(i == 1) pads[i] = new TPad(("pad" + std::to_string(i)).c_str(), "", 0.0, 0.0, 0.5, 0.35);
    else if(i == 2) pads[i] = new TPad(("pad" + std::to_string(i)).c_str(), "", 0.5, 0.35, 1.0, 1.0);
    else if(i == 3) pads[i] = new TPad(("pad" + std::to_string(i)).c_str(), "", 0.5, 0.0, 1.0, 0.35);

    pads[i]->Draw("SAME");
    pads[i]->cd();
  }

  validation_p->cd();
  pads[0]->cd();
  pads[0]->SetBottomMargin(0.01);
  gStyle->SetOptStat(0);


  histTot_p->GetYaxis()->SetTitleOffset(histTot_p->GetYaxis()->GetTitleOffset()*1.5);
  histTot_p->DrawCopy("HIST E1 P");
  comboTot_p->DrawCopy("SAME HIST E1 P");

  validation_p->cd();
  pads[1]->cd();
  pads[1]->SetTopMargin(0.01);
  gStyle->SetOptStat(0);

  TH1D* div = (TH1D*)histTot_p->Clone("div");
  div->Divide(comboTot_p);
  div->GetYaxis()->SetTitle("Inclusive/Combo");
  //  div->GetYaxis()->SetTitleOffset(div->GetYaxis()->GetTitleOffset()*1.5);
  div->GetYaxis()->SetNdivisions(404);
  div->DrawCopy("HIST E1 P");
  delete div;
  
  validation_p->cd();
  pads[2]->cd();
  pads[2]->SetBottomMargin(0.01);
  gStyle->SetOptStat(0);

  TH1D* errHistTot_p = (TH1D*)histTot_p->Clone("errHistTot_p");
  TH1D* errComboTot_p = (TH1D*)comboTot_p->Clone("errComboTot_p");

  for(Int_t bIX = 0; bIX < errHistTot_p->GetNbinsX(); ++bIX){
    errHistTot_p->SetBinContent(bIX+1, errHistTot_p->GetBinError(bIX+1));
    errHistTot_p->SetBinError(bIX+1, 0.0);

    errComboTot_p->SetBinContent(bIX+1, errComboTot_p->GetBinError(bIX+1));
    errComboTot_p->SetBinError(bIX+1, 0.0);
  }

  errHistTot_p->GetYaxis()->SetTitle("Error");
  errHistTot_p->DrawCopy("HIST E1 P");
  errComboTot_p->DrawCopy("HIST E1 P SAME");

  validation_p->cd();
  pads[3]->cd();
  pads[3]->SetTopMargin(0.01);
  gStyle->SetOptStat(0);

  errHistTot_p->Divide(errComboTot_p);
  errHistTot_p->GetYaxis()->SetTitle("Inclusive/Combo");
  //  errHistTot_p->GetYaxis()->SetTitleOffset(errHistTot_p->GetYaxis()->GetTitleOffset()*1.5);
  errHistTot_p->GetYaxis()->SetNdivisions(404);
  errHistTot_p->SetMaximum(1.1);
  errHistTot_p->SetMinimum(0.9);
  errHistTot_p->DrawCopy("HIST E1 P");

  delete errHistTot_p;
  delete errComboTot_p;

  std::string dirName1 = "pdfDir/";
  std::string dirName2 = dirName1 + std::to_string(date->GetDate()) + "/";
  checkMakeDir(dirName1);
  checkMakeDir(dirName2);
  std::string saveName = dirName2 + "comboCheck_" + std::string(histTot_p->GetName()) + "_" + std::to_string(date->GetDate()) + ".pdf";

  validation_p->SaveAs(saveName.c_str());

  for(Int_t i = 0; i < nPadTot; ++i){
    delete pads[i];
  }
  delete validation_p;
  delete date;

  return;
}

void plotValidation(TCanvas* validation_p, TPad* pads[nPadTotValid], std::vector<TH1D*> inTH1_p, RooUnfoldResponse* rooRes, RooUnfoldBayes* rooUnfBayes, RooUnfoldSvd* rooUnfSvd, const bool isBayes, const unsigned int iterPos, const std::string ppPbPbStr)
{
  const Int_t nTotHist = 7;

  if(inTH1_p.size() != nTotHist){
    std::cout << "WARNING IN PLOT VALIDATION: INSUFFICIENT HIST PASSED, GIVEN \'" << inTH1_p.size() << "\', require \'" << nTotHist << "\'. return" << std::endl;
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  TLatex label_p;
  label_p.SetTextFont(43);
  label_p.SetTextSize(14);
  label_p.SetNDC();

  TDatime date;
  kirchnerPalette col;

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::string dirName1 = "pdfDir/";
  std::string dirName2 = dirName1 + std::to_string(date.GetDate()) + "/";
  std::string saveName = dirName2 + std::string(inTH1_p.at(0)->GetName()) + "_" + std::to_string(date.GetDate()) + ".pdf";

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  const std::string uniqueStr = std::to_string(iterPos) + "_" + ppPbPbStr;

  validation_p->SetTopMargin(0.01);
  validation_p->SetBottomMargin(0.01);
  validation_p->SetRightMargin(0.01);
  validation_p->SetLeftMargin(0.01);
  
  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  validation_p->cd();
  pads[0]->cd();

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << ", " << validation_p->GetName() << std::endl;

  initTH1(inTH1_p.at(0), 24, .8, 0);
  initTH1(inTH1_p.at(1), 21, .8, 2);
  initTH1(inTH1_p.at(2), 24, .8, 0);
  initTH1(inTH1_p.at(3), 21, .8, 2);
  initTH1(inTH1_p.at(4), 24, .8, 0);
  if(inTH1_p.at(5) != NULL) initTH1(inTH1_p.at(5), 21, .8, 2);
  initTH1(inTH1_p.at(6), 21, .8, 2);

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
 

  inTH1_p.at(0)->DrawCopy("HIST E1 P");
  inTH1_p.at(1)->DrawCopy("HIST E1 P SAME");
  label_p.DrawLatex(.3, .9, "MC (Self-consistent)");
  pads[0]->SetLogy();


  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  validation_p->cd();
  pads[1]->cd();

  TH1D* tempDiv_p = (TH1D*)inTH1_p.at(0)->Clone(("tempDiv_" + uniqueStr + "_p").c_str());
  tempDiv_p->Divide(inTH1_p.at(1));

  tempDiv_p->SetMaximum(1.2);
  tempDiv_p->SetMinimum(0.8);

  tempDiv_p->DrawCopy("HIST E1 P");
  delete tempDiv_p;

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  validation_p->cd();
  pads[2]->cd();
  inTH1_p.at(2)->DrawCopy("HIST E1 P");
  inTH1_p.at(3)->DrawCopy("HIST E1 P SAME");
  label_p.DrawLatex(.3, .9, "MC (Orthogonal)");
  pads[2]->SetLogy();

  validation_p->cd();
  pads[3]->cd();

  tempDiv_p = (TH1D*)inTH1_p.at(2)->Clone(("tempDiv_" + uniqueStr + "_p").c_str());
  tempDiv_p->Divide(inTH1_p.at(3));

  tempDiv_p->SetMaximum(1.2);
  tempDiv_p->SetMinimum(0.8);

  tempDiv_p->DrawCopy("HIST E1 P");
  delete tempDiv_p;

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  validation_p->cd();
  pads[4]->cd();
  inTH1_p.at(4)->DrawCopy("HIST E1 P");
  if(inTH1_p.at(5) != NULL) inTH1_p.at(5)->DrawCopy("HIST E1 P SAME");

  label_p.DrawLatex(.3, .9, "Data (Prev. Iter)");
  pads[4]->SetLogy();

  if(inTH1_p.at(5) != NULL){
    validation_p->cd();
    pads[5]->cd();    
    tempDiv_p = (TH1D*)inTH1_p.at(4)->Clone(("tempDiv_" + uniqueStr + "_p").c_str());
    tempDiv_p->Divide(inTH1_p.at(5));
    tempDiv_p->SetMaximum(1.2);
    tempDiv_p->SetMinimum(0.8);
    tempDiv_p->DrawCopy("HIST E1 P");    
    delete tempDiv_p;
  }

  validation_p->cd();
  pads[6]->cd();

  TH1D* refolded_h = (TH1D*)rooRes->ApplyToTruth(inTH1_p.at(4));
  initTH1(refolded_h, 24, .8, 0);
  centerTitles({refolded_h, inTH1_p.at(6)});

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  refolded_h->DrawCopy("HIST E1");
  inTH1_p.at(6)->DrawCopy("HIST E1 P SAME");
  pads[6]->SetLogy();

  label_p.DrawLatex(.3, .9, "Data (Refolded)");

  validation_p->cd();
  pads[7]->cd();

  refolded_h->Divide(inTH1_p.at(6));
  refolded_h->SetMaximum(1.2);
  refolded_h->SetMinimum(0.8);
  refolded_h->DrawCopy("HIST E1 P");
  delete refolded_h;

  validation_p->cd();
  pads[8]->cd();
  pads[8]->SetTopMargin(pads[8]->GetBottomMargin());

  if(!isBayes && rooUnfSvd != NULL){
    std::string dVectName = rooUnfSvd->GetName();
    dVectName = dVectName + "_D_TEMP";
    TH1D* dVect_p = (TH1D*)rooUnfSvd->Impl()->GetD()->Clone(dVectName.c_str());
    dVect_p->DrawCopy("HIST E1 P");
    label_p.DrawLatex(.3, .93, "D-Vector");

    delete dVect_p;
  }

  validation_p->cd();
  pads[9]->cd();
  pads[9]->SetTopMargin(pads[9]->GetBottomMargin());


  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << ", " << inTH1_p.at(6)->GetName() << std::endl;

  TH2D* res = (TH2D*)rooRes->Hresponse()->Clone(("tempRes_" + uniqueStr).c_str());
  res->SetTitle("");
  centerTitles(res);
  for(Int_t bIY = 0; bIY < res->GetNbinsY(); ++bIY){
    Double_t total = 0.0;
    for(Int_t bIX = 0; bIX < res->GetNbinsX(); ++bIX){      
      total += res->GetBinContent(bIX+1, bIY+1);
    }

    for(Int_t bIX = 0; bIX < res->GetNbinsX(); ++bIX){      
      res->SetBinContent(bIX+1, bIY+1, res->GetBinContent(bIX+1, bIY+1)/total);
    }
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  res->SetMarkerSize(1.25);
  gStyle->SetPaintTextFormat("1.3f");
  res->DrawCopy("COLZ TEXT");
  label_p.DrawLatex(.3, .93, "Response");

  pads[9]->SetLogz();
  delete res;

  validation_p->cd();
  pads[10]->cd();
  pads[10]->SetTopMargin(pads[10]->GetBottomMargin());

  if(!isBayes && rooUnfSvd != NULL){
    std::string sVectName = rooUnfSvd->GetName();
    sVectName = sVectName + "_SV_TEMP";
    TH1D* sVect_p = (TH1D*)rooUnfSvd->Impl()->GetSV()->Clone(sVectName.c_str());


    sVect_p->DrawCopy("HIST E1 P");
    label_p.DrawLatex(.3, .93, "Singluar Values");

    delete sVect_p;
  }

  validation_p->cd();
  pads[11]->cd();
  pads[11]->SetTopMargin(pads[11]->GetBottomMargin());

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  TMatrixD* tempCov = NULL;
  if(isBayes){
    if(rooUnfBayes != NULL) tempCov = (TMatrixD*)(rooUnfBayes->Ereco(RooUnfold::kCovToy).Clone("tempCov"));
  }
  else{
    if(rooUnfSvd != NULL) tempCov = (TMatrixD*)(rooUnfSvd->Ereco(RooUnfold::kCovToy).Clone("tempCov"));
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  if(tempCov != NULL){
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    TMatrixD* tempCovPears = (TMatrixD*)tempCov->Clone("tempCovPears");
    
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    
    for(Int_t rI = 0; rI < tempCovPears->GetNrows(); ++rI){
      for(Int_t cI = 0; cI < tempCovPears->GetNcols(); ++cI){
	Double_t val = (*(tempCov))(rI,cI);
	if(TMath::Abs((*(tempCov))(rI,rI)) > 0.00000001 && TMath::Abs((*(tempCov))(cI,cI)) > 0.00000001) val /= TMath::Sqrt((*(tempCov))(cI,cI)*(*(tempCov))(rI,rI));
	else if(TMath::Abs(val) < 0.00000001) val = 0;
	else{
	  //	  std::cout << "UHOH WARNING: PEARSON DENOM IS 0 WITHOUT NUM BEING 0" << std::endl;
	  val = 0.;
	}
	
	(*(tempCovPears))(rI,cI) = val;
      }
    }
  
    TH2D* cov = new TH2D(*tempCovPears);
    if(doGlobalDebug) std::cout << "PRINT cov" << std::endl;
    if(doGlobalDebug) cov->Print("ALL");
    cov->DrawCopy("COLZ");

    label_p.DrawLatex(.3, .93, "Pearson (#bf{#color[2]{Red>=10%}}");

    TLatex label2_p;
    label2_p.SetTextFont(43);
    label2_p.SetTextSize(9);
  
    for(Int_t bIX = 0; bIX < cov->GetNbinsX(); ++bIX){
      for(Int_t bIY = 0; bIY < cov->GetNbinsY(); ++bIY){
	Double_t val = cov->GetBinContent(bIX+1, bIY+1);
	
	if(val < .1) label2_p.DrawLatex(bIX+.2, bIY+.5, prettyString(val, 3, false).c_str());
	else label2_p.DrawLatex(bIX+.2, bIY+.5, ("#bf{#color[2]{" + prettyString(val, 3, false) + "}}").c_str());
      }
    }
    
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    delete cov;
  }

  return;
}

double getComboBinVal(TH1* inHist1_p, TH1* inHist2_p, Int_t binPos)
{
  double combo = inHist1_p->GetBinContent(binPos) + inHist2_p->GetBinContent(binPos);
  return combo;
}

double getComboBinErr(TH1* inHist1_p, TH1* inHist2_p, Int_t binPos)
{
  double err = TMath::Sqrt(inHist1_p->GetBinContent(binPos)*inHist1_p->GetBinContent(binPos) + inHist2_p->GetBinContent(binPos)*inHist2_p->GetBinContent(binPos));
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


void drawSyst(TCanvas* canv_p, std::vector<TH1*> inSystHist)
{
  if(inSystHist.size() <= 1) return;

  const Int_t nBins = inSystHist.at(0)->GetNbinsX();
  const Int_t nSubSyst = 8;
  const std::string subSystStr[nSubSyst] = {"JECVarData", "JECVarMC", "JER1p15", "LumiUp", "LumiDown", "TAAUp", "TAADown", "Fake"};
  const int subSystCol[nSubSyst] = {1, 0, 2, 4, 4, 3, 3, 6};
  const int subSystStyle[nSubSyst] = {1, 1, 1, 2, 2, 2, 2, 1};

  Double_t systVal[nSubSyst][nBins];
  for(Int_t sI = 0; sI < nSubSyst; ++sI){
    for(Int_t bIX = 0; bIX < nBins; ++bIX){
      systVal[sI][bIX] = 0;
    }
  }

  TDatime* date = new TDatime();
  std::string dirName1 = "pdfDir/";
  std::string dirName2 = dirName1 + std::to_string(date->GetDate()) + "/";
  std::string saveName = canv_p->GetName();
  saveName = dirName2 + "syst_" + saveName + "_" + std::to_string(date->GetDate()) + ".pdf";
  delete date;
  TCanvas* systCanv_p = new TCanvas("temp", "temp", 400, 400);

  std::vector<double> systVect;
  for(unsigned int i = 1; i < inSystHist.size(); ++i){
    std::string algoStr = inSystHist.at(i)->GetName();
    int pos = -1;
    for(Int_t sI = 0; sI < nSubSyst; ++sI){
      if(algoStr.find(subSystStr[sI]) != std::string::npos) pos = sI;
    }

    if(pos < 0){
      if(algoStr.find("Data") != std::string::npos) pos = 0;
      else if(algoStr.find("MC") != std::string::npos) pos = 1;
    }

    for(Int_t bIX = 0; bIX < inSystHist.at(0)->GetNbinsX(); ++bIX){
      Double_t val = TMath::Abs(inSystHist.at(0)->GetBinContent(bIX+1) - inSystHist.at(i)->GetBinContent(bIX+1));
      if(val > systVal[pos][bIX]) systVal[pos][bIX] = val;
    }
  }

  std::string algoStr = inSystHist.at(0)->GetName();

  Int_t colPos = 0;
  if(algoStr.find("akCs4P") != std::string::npos) colPos = 4;
  else if(algoStr.find("ak4PF") != std::string::npos) colPos = 4;
  else if(algoStr.find("akCs6P") != std::string::npos) colPos = 6;
  else if(algoStr.find("ak6PF") != std::string::npos) colPos = 6;
  else if(algoStr.find("akCs8P") != std::string::npos) colPos = 5;
  else if(algoStr.find("ak8PF") != std::string::npos) colPos = 5;
  else if(algoStr.find("akCs10P") != std::string::npos) colPos = 2;
  else if(algoStr.find("ak10PF") != std::string::npos) colPos = 2;

  canv_p->cd();

  kirchnerPalette col;

  TBox* tempBox_p = new TBox();
  tempBox_p->SetFillColorAlpha(col.getColor(colPos), .25);

  for(Int_t bIX = 0; bIX < inSystHist.at(0)->GetNbinsX(); ++bIX){
    Double_t totVal = 0.;
    for(Int_t sI = 0; sI < nSubSyst; ++sI){
      if(subSystStr[sI].find("TAA") != std::string::npos) continue;
      else if(subSystStr[sI].find("Lumi") != std::string::npos) continue;

      totVal = TMath::Sqrt(totVal*totVal + systVal[sI][bIX]*systVal[sI][bIX]);
    }
    if(inSystHist.at(0)->GetBinLowEdge(bIX+1)>=199){
      tempBox_p->DrawBox(inSystHist.at(0)->GetBinLowEdge(bIX+1), inSystHist.at(0)->GetBinContent(bIX+1) - totVal, inSystHist.at(0)->GetBinLowEdge(bIX+2), inSystHist.at(0)->GetBinContent(bIX+1) + totVal);
    }
  }

  Double_t lumiUpErr = TMath::Abs(1 - systVal[3][2]/inSystHist.at(0)->GetBinContent(3));
  Double_t lumiDownErr = TMath::Abs(1 - systVal[4][2]/inSystHist.at(0)->GetBinContent(3));
  Double_t taaUpErr = TMath::Abs(1 - systVal[5][2]/inSystHist.at(0)->GetBinContent(3));
  Double_t taaDownErr = TMath::Abs(1 - systVal[6][2]/inSystHist.at(0)->GetBinContent(3));

  Double_t histXBins[20+1];
  getLogBins(inSystHist.at(0)->GetBinLowEdge(1), inSystHist.at(0)->GetBinLowEdge(inSystHist.at(0)->GetNbinsX()+1), 20, histXBins);

  tempBox_p->SetFillColor(col.getColor(1));
  tempBox_p->DrawBox(histXBins[15], 1 - lumiDownErr, histXBins[17], 1 + lumiUpErr);

  tempBox_p->SetFillColor(col.getColor(2));
  tempBox_p->DrawBox(histXBins[18], 1 - taaDownErr, histXBins[20], 1 + taaUpErr);
  
  gPad->RedrawAxis();

  delete tempBox_p;

  systCanv_p->cd();

  TH1D* tempTotalSyst_p = (TH1D*)inSystHist.at(0)->Clone("tempTotalSyst");
  for(Int_t bIX = 0; bIX < tempTotalSyst_p->GetNbinsX(); ++bIX){
    tempTotalSyst_p->SetBinContent(bIX+1, 0.0);
  }

  for(Int_t sI = 0; sI < nSubSyst; ++sI){
    TH1D* tempSyst_p = (TH1D*)inSystHist.at(0)->Clone("tempSyst");

    for(Int_t bIX = 0; bIX < tempSyst_p->GetNbinsX(); ++bIX){
      tempSyst_p->SetBinContent(bIX+1, systVal[sI][bIX]);

      Double_t tot = TMath::Sqrt(tempTotalSyst_p->GetBinContent(bIX+1)*tempTotalSyst_p->GetBinContent(bIX+1) + systVal[sI][bIX]*systVal[sI][bIX]);
      tempTotalSyst_p->SetBinContent(bIX+1, tot);
    }

    tempSyst_p->SetLineWidth(2);
    tempSyst_p->SetLineStyle(subSystStyle[sI]);
    tempSyst_p->SetMarkerColor(col.getColor(subSystCol[sI]));
    tempSyst_p->SetLineColor(col.getColor(subSystCol[sI]));

    if(sI) tempSyst_p->DrawCopy("HIST L");
    else tempSyst_p->DrawCopy("HIST L SAME");
    delete tempSyst_p;
  }

  tempTotalSyst_p->SetLineWidth(2);
  tempTotalSyst_p->SetMarkerColor(1);
  tempTotalSyst_p->SetLineColor(1);
  tempTotalSyst_p->DrawCopy("HIST L SAME");

  systCanv_p->SaveAs(saveName.c_str());
  delete systCanv_p;

  return;
}


int makePlotValidation_FromTree(std::string inFileName)
{
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

  std::string pbpbTreeStr;
  std::string ppTreeStr;
  Int_t rVals;

  for(std::map<std::string, std::string>::iterator it = jetMapPbPbToPP.begin(); it != jetMapPbPbToPP.end(); ++it){
    if(inFileName.find(it->first) != std::string::npos){
      rVals = getRVal(it->first);
      break;
    }
  }


  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  std::vector<std::string> retDir = returnRootFileContentsList(inFile_p, "TDirectoryFile");

  if(retDir.size() < 2) std::cout << "WARNING LESS THAN TWO DIR" << std::endl;

  for(unsigned int rI = 0; rI < retDir.size(); ++rI){
    if(retDir.at(rI).find("Cs") != std::string::npos) pbpbTreeStr = retDir.at(rI).substr(0, retDir.at(rI).find("/"));
    else ppTreeStr = retDir.at(rI).substr(0, retDir.at(rI).find("/"));
  }

  inFile_p->cd();

  const Int_t nCentBins = std::stoi(std::string(inFile_p->Get("nCentBins")->GetTitle()));
  const Int_t nAbsEtaBins = std::stoi(std::string(inFile_p->Get("nAbsEtaBins")->GetTitle()));
  const Int_t nSyst = std::stoi(std::string(inFile_p->Get("nSyst")->GetTitle()));
  const Int_t nSystForFill = std::stoi(std::string(inFile_p->Get("nSystForFill")->GetTitle()));
  const Int_t nSystForCopy = std::stoi(std::string(inFile_p->Get("nSystForCopy")->GetTitle()));
  const Int_t nBayesIter = std::stoi(std::string(inFile_p->Get("nBayesIter")->GetTitle())); 
  const Bool_t doSvd = false;//std::stoi(std::string(inFile_p->Get("doSvd")->GetTitle())); 
  const Int_t nSvd = std::stoi(std::string(inFile_p->Get("nSvd")->GetTitle())); 
 
  std::string centBinsLowStr = std::string(inFile_p->Get("centBinsLow")->GetTitle());
  std::string centBinsHiStr = std::string(inFile_p->Get("centBinsHi")->GetTitle());
  Int_t centBinsLow[nCentBins];
  Int_t centBinsHi[nCentBins];

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    centBinsLow[cI] = std::stoi(centBinsLowStr.substr(0,centBinsLowStr.find(",")));
    centBinsHi[cI] = std::stoi(centBinsHiStr.substr(0,centBinsHiStr.find(",")));
    centBinsLowStr.replace(0, centBinsLowStr.find(",")+1, "");
    centBinsHiStr.replace(0, centBinsHiStr.find(",")+1, "");
  }

  std::string absEtaBinsLowStr = std::string(inFile_p->Get("absEtaBinsLow")->GetTitle());
  std::string absEtaBinsHiStr = std::string(inFile_p->Get("absEtaBinsHi")->GetTitle());
  std::string absEtaBinsLow[nAbsEtaBins];
  std::string absEtaBinsHi[nAbsEtaBins];

  for(Int_t aI = 0; aI < nAbsEtaBins; ++aI){
    absEtaBinsLow[aI] = absEtaBinsLowStr.substr(0,absEtaBinsLowStr.find(","));
    absEtaBinsHi[aI] = absEtaBinsHiStr.substr(0,absEtaBinsHiStr.find(","));
    absEtaBinsLowStr.replace(0, absEtaBinsLowStr.find(",")+1, "");
    absEtaBinsHiStr.replace(0, absEtaBinsHiStr.find(",")+1, "");
  }

  std::string systStrTemp = std::string(inFile_p->Get("systStr")->GetTitle());
  std::string systStrForFillTemp = std::string(inFile_p->Get("systStrForFill")->GetTitle());
  std::string systStrForCopyTemp = std::string(inFile_p->Get("systStrForCopy")->GetTitle());

  std::string systStr[nSyst];
  std::string systStrForFill[nSystForFill];
  std::string systStrForCopy[nSystForCopy];
  Int_t systPosForFill[nSystForFill];
  Int_t systPosForCopy[nSystForCopy];
 

  for(Int_t sI = 0; sI < nSyst; ++sI){
    systStr[sI] = systStrTemp.substr(0, systStrTemp.find(","));
    systStrTemp.replace(0,systStrTemp.find(",")+1,"");
  }

  for(Int_t sI = 0; sI < nSystForFill; ++sI){
    systStrForFill[sI] = systStrForFillTemp.substr(0, systStrForFillTemp.find(","));
    systStrForFillTemp.replace(0,systStrForFillTemp.find(",")+1,"");

    for(Int_t sI2 = 0; sI2 < nSyst; ++sI2){
      if(systStrForFill[sI].size() == systStr[sI2].size() && systStrForFill[sI].find(systStr[sI2]) != std::string::npos){
	systPosForFill[sI] = sI2;
	break;
      }
    }
  }

  for(Int_t sI = 0; sI < nSystForCopy; ++sI){
    systStrForCopy[sI] = systStrForCopyTemp.substr(0, systStrForCopyTemp.find(","));
    systStrForCopyTemp.replace(0,systStrForCopyTemp.find(",")+1,"");

    for(Int_t sI2 = 0; sI2 < nSyst; ++sI2){
      if(systStrForCopy[sI].size() == systStr[sI2].size() && systStrForCopy[sI].find(systStr[sI2]) != std::string::npos){
	systPosForCopy[sI] = sI2;
	break;
      }
    }
  }

  Int_t truthNBins[nCentBins];
  Int_t recoNBins[nCentBins];
  Float_t truthBins[nCentBins][100];
  Float_t recoBins[nCentBins][100];

  Int_t maxTruthBins = -1;

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    const Int_t rVal = getRVal(pbpbTreeStr);
    recoNBins[cI] = binner.getNBins(centBinsLow[cI], centBinsHi[cI], rVal, true);
    truthNBins[cI] = binner.getNBins(centBinsLow[cI], centBinsHi[cI], rVal, false);
    binner.getBins(centBinsLow[cI], centBinsHi[cI], rVal, true, recoNBins[cI], recoBins[cI]);
    binner.getBins(centBinsLow[cI], centBinsHi[cI], rVal, false, truthNBins[cI], truthBins[cI]);
    
    if(truthNBins[cI] > maxTruthBins) maxTruthBins = truthNBins[cI];
  }

  std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  TCanvas* canvPbPbData_FULLBayes_h[nCentBins][nAbsEtaBins+1][nSyst][nBayesIter];
  TCanvas* canvPbPbData_FULLSvd_h[nCentBins][nAbsEtaBins+1][nSyst][nSvd];
  TPad* padsPbPbData_FULLBayes_h[nCentBins][nAbsEtaBins+1][nSyst][nBayesIter][nPadTotValid];
  TPad* padsPbPbData_FULLSvd_h[nCentBins][nAbsEtaBins+1][nSyst][nSvd][nPadTotValid];

  TCanvas* canvPPData_FULLBayes_h[nCentBins][nAbsEtaBins+1][nSyst][nBayesIter];
  TCanvas* canvPPData_FULLSvd_h[nCentBins][nAbsEtaBins+1][nSyst][nSvd];
  TPad* padsPPData_FULLBayes_h[nCentBins][nAbsEtaBins+1][nSyst][nBayesIter][nPadTotValid];
  TPad* padsPPData_FULLSvd_h[nCentBins][nAbsEtaBins+1][nSyst][nSvd][nPadTotValid];
  
  std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  const Int_t nRow = 2;
  const Int_t nColumn = 4;

  const Int_t nPadX = 4;
  const Int_t nPadY = 4;
  Double_t xLow[nPadX] = {0.0, 0.25, 0.5, .75};
  Double_t xHi[nPadX] = {0.25, 0.5, .75, 1.0};
  Double_t yLow[nPadY] = {0.675, 0.5, 0.0, 0.0};
  Double_t yHi[nPadY] = {1.0, 0.675, 0.5, 0.25};

  const Double_t xLowTot[nPadTotValid] = {xLow[0], xLow[0], xLow[1], xLow[1], xLow[2], xLow[2], xLow[3], xLow[3], xLow[0], xLow[1], xLow[2], xLow[3], xLow[0]};
  const Double_t xHiTot[nPadTotValid] = {xHi[0], xHi[0], xHi[1], xHi[1], xHi[2], xHi[2], xHi[3], xHi[3], xHi[0], xHi[1], xHi[2], xHi[3], xHi[0]};
  const Double_t yLowTot[nPadTotValid] = {yLow[0], yLow[1], yLow[0], yLow[1], yLow[0], yLow[1], yLow[0], yLow[1], yLow[2], yLow[2], yLow[2], yLow[2], yLow[3]};
  const Double_t yHiTot[nPadTotValid] = {yHi[0], yHi[1], yHi[0], yHi[1], yHi[0], yHi[1], yHi[0], yHi[1], yHi[2], yHi[2], yHi[2], yHi[2], yHi[3]};

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  
  const std::string algoStrPbPb = pbpbTreeStr;
  const std::string algoStrPP = ppTreeStr;

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    const std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
    
    for(Int_t aI = 0; aI < nAbsEtaBins+1; ++aI){
      std::string absEtaStr = "AbsEtaTot";
      if(aI < nAbsEtaBins) absEtaStr = "AbsEta" + absEtaBinsLow[aI] + "to" + absEtaBinsHi[aI];
      const std::string binStrPbPb = algoStrPbPb + "_" + centStr + "_" + absEtaStr;
      
      for(Int_t sI = 0; sI < nSystForFill; ++sI){	  
	for(Int_t bI = 0; bI < nBayesIter; ++bI){
	  canvPbPbData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI] = new TCanvas(("canvPbPbData_" + binStrPbPb + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_h").c_str(), "", nColumn*400, nRow*250);
	  
	  for(Int_t pI = 0; pI < nPadTotValid; ++pI){
	    canvPbPbData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI]->cd();	      
	    padsPbPbData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI][pI] = new TPad(("padsPbPbData_" + binStrPbPb + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_" + std::to_string(pI) + "_h").c_str(), "", xLowTot[pI], yLowTot[pI], xHiTot[pI], yHiTot[pI]);

	    if(pI != 12) padsPbPbData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI][pI]->Draw("SAME");
	    padsPbPbData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI][pI]->SetTopMargin(0.01);
	    if(pI < 8) padsPbPbData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI][pI]->SetRightMargin(0.01);
	    if(pI < 8 && pI%2 == 0) padsPbPbData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI][pI]->SetBottomMargin(0.01);
	  }
	}
	
	for(Int_t bI = 0; bI < nSvd; ++bI){
	  if(bI+1 > truthNBins[cI]) break;
	  canvPbPbData_FULLSvd_h[cI][aI][systPosForFill[sI]][bI] = new TCanvas(("canvPbPbData_" + binStrPbPb + "_FULLSvd" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_h").c_str(), "", nColumn*400, nRow*250);
	  
	  for(Int_t pI = 0; pI < nPadTotValid; ++pI){
	    canvPbPbData_FULLSvd_h[cI][aI][systPosForFill[sI]][bI]->cd();
	    padsPbPbData_FULLSvd_h[cI][aI][systPosForFill[sI]][bI][pI] = new TPad(("padsPbPbData_" + binStrPbPb + "_FULLSvd" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_" + std::to_string(pI) + "_h").c_str(), "", xLowTot[pI], yLowTot[pI], xHiTot[pI], yHiTot[pI]);
	    padsPbPbData_FULLSvd_h[cI][aI][systPosForFill[sI]][bI][pI]->Draw("SAME");
	    padsPbPbData_FULLSvd_h[cI][aI][systPosForFill[sI]][bI][pI]->SetTopMargin(0.01);
	    if(pI < 8) padsPbPbData_FULLSvd_h[cI][aI][systPosForFill[sI]][bI][pI]->SetRightMargin(0.01);
	    if(pI < 8 && pI%2 == 0) padsPbPbData_FULLSvd_h[cI][aI][systPosForFill[sI]][bI][pI]->SetBottomMargin(0.01);
	  }
	}
      }
      
      const std::string binStrPP = algoStrPP + "_" + centStr + "_" + absEtaStr;
      
      for(Int_t sI = 0; sI < nSystForFill; ++sI){	  
	for(Int_t bI = 0; bI < nBayesIter; ++bI){
	  canvPPData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI] = new TCanvas(("canvPPData_" + binStrPP + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_h").c_str(), "", nColumn*400, nRow*250);
	  
	  for(Int_t pI = 0; pI < nPadTotValid; ++pI){
	    canvPPData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI]->cd();
	    padsPPData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI][pI] = new TPad(("padsPPData_" + binStrPP + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_" + std::to_string(pI) + "_h").c_str(), "", xLowTot[pI], yLowTot[pI], xHiTot[pI], yHiTot[pI]);
	    padsPPData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI][pI]->Draw("SAME");
	    padsPPData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI][pI]->SetTopMargin(0.01);
	    if(pI < 8) padsPPData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI][pI]->SetRightMargin(0.01);
	    if(pI < 8 && pI%2 == 0) padsPPData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI][pI]->SetBottomMargin(0.01);
	  }
	}
	
	for(Int_t bI = 0; bI < nSvd; ++bI){
	  if(bI+1 > truthNBins[cI]) break;
	  
	  canvPPData_FULLSvd_h[cI][aI][systPosForFill[sI]][bI] = new TCanvas(("canvPPData_" + binStrPP + "_FULLSvd" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_h").c_str(), "", nColumn*400, nRow*250);
	  
	  for(Int_t pI = 0; pI < nPadTotValid; ++pI){
	    canvPPData_FULLSvd_h[cI][aI][systPosForFill[sI]][bI]->cd();
	    padsPPData_FULLSvd_h[cI][aI][systPosForFill[sI]][bI][pI] = new TPad(("padsPPData_" + binStrPP + "_FULLSvd" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_" + std::to_string(pI) + "_h").c_str(), "", xLowTot[pI], yLowTot[pI], xHiTot[pI], yHiTot[pI]);
	    padsPPData_FULLSvd_h[cI][aI][systPosForFill[sI]][bI][pI]->Draw("SAME");
	    padsPPData_FULLSvd_h[cI][aI][systPosForFill[sI]][bI][pI]->SetTopMargin(0.01);
	    if(pI < 8) padsPPData_FULLSvd_h[cI][aI][systPosForFill[sI]][bI][pI]->SetRightMargin(0.01);
	    if(pI < 8 && pI%2 == 0) padsPPData_FULLSvd_h[cI][aI][systPosForFill[sI]][bI][pI]->SetBottomMargin(0.01);
	  }
	}
      }
    }
  }

  std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  TH1D* jtptPbPbData_RAW_h[nCentBins][nAbsEtaBins+1][nSystForFill];
  TH1D* jtptPbPbData_FULLBayes_h[nCentBins][nAbsEtaBins+1][nSystForFill][nBayesIter];
  TH1D* jtptPbPbData_FULLSvd_h[nCentBins][nAbsEtaBins+1][nSystForFill][nSvd];
  TH1D* jtetaPbPbData_FULLBayes_h[nCentBins][nSystForFill][nBayesIter];
  RooUnfoldBayes* rooUnfoldBayesPbPbData_FULLBayes_h[nCentBins][nAbsEtaBins+1][nSystForFill][nBayesIter];
  RooUnfoldSvd* rooUnfoldSvdPbPbData_FULLSvd_h[nCentBins][nAbsEtaBins+1][nSystForFill][nSvd];

  TH1D* jtptPPData_RAW_h[nCentBins][nAbsEtaBins+1][nSystForFill];
  TH1D* jtptPPData_FULLBayes_h[nCentBins][nAbsEtaBins+1][nSystForFill][nBayesIter];
  TH1D* jtptPPData_FULLSvd_h[nCentBins][nAbsEtaBins+1][nSystForFill][nSvd];
  TH1D* jtetaPPData_FULLBayes_h[nCentBins][nSystForFill][nBayesIter];
  RooUnfoldBayes* rooUnfoldBayesPPData_FULLBayes_h[nCentBins][nAbsEtaBins+1][nSystForFill][nBayesIter];
  RooUnfoldSvd* rooUnfoldSvdPPData_FULLSvd_h[nCentBins][nAbsEtaBins+1][nSystForFill][nSvd];

  TH1D* jtptPbPbMC_RAW_h[nCentBins][nAbsEtaBins+1][nSystForFill];
  TH1D* jtptPbPbMC_FULLBayes_h[nCentBins][nAbsEtaBins+1][nSystForFill][nBayesIter];
  RooUnfoldBayes* rooUnfoldBayesPbPbMC_FULLBayes_h[nCentBins][nAbsEtaBins+1][nSystForFill][nBayesIter];
  TH1D* jtptPbPbMC_FULLSvd_h[nCentBins][nAbsEtaBins+1][nSystForFill][nSvd];
  TH1D* refptPbPbMC_h[nCentBins][nAbsEtaBins+1];
  TH1D* jtptPbPbMC_RAW_Frac_h[nCentBins][nAbsEtaBins+1][nSystForFill];
  TH1D* jtptPbPbMC_FULLBayes_Frac_h[nCentBins][nAbsEtaBins+1][nSystForFill][nBayesIter];
  RooUnfoldBayes* rooUnfoldBayesPbPbMC_FULLBayes_Frac_h[nCentBins][nAbsEtaBins+1][nSystForFill][nBayesIter];
  TH1D* jtptPbPbMC_FULLSvd_Frac_h[nCentBins][nAbsEtaBins+1][nSystForFill][nSvd];
  TH1D* refptPbPbMC_Frac_h[nCentBins][nAbsEtaBins+1];
  TH2D* responsePbPbMC_h[nCentBins][nAbsEtaBins+1];
  RooUnfoldResponse* rooResponsePbPbMC_h[nCentBins][nAbsEtaBins+1];
  TH1D* jtetaPbPbMC_FULLBayes_h[nCentBins][nSystForFill][nBayesIter];
  TH1D* jtetaPbPbMC_FULLBayes_Frac_h[nCentBins][nSystForFill][nBayesIter];

  TH1D* jtptPPMC_RAW_h[nCentBins][nAbsEtaBins+1][nSystForFill]; 
  TH1D* jtptPPMC_FULLBayes_h[nCentBins][nAbsEtaBins+1][nSystForFill][nBayesIter];
  RooUnfoldBayes* rooUnfoldBayesPPMC_FULLBayes_h[nCentBins][nAbsEtaBins+1][nSystForFill][nBayesIter];
  TH1D* jtptPPMC_FULLSvd_h[nCentBins][nAbsEtaBins+1][nSystForFill][nSvd];
  TH1D* refptPPMC_h[nCentBins][nAbsEtaBins+1];
  TH1D* jtptPPMC_RAW_Frac_h[nCentBins][nAbsEtaBins+1][nSystForFill]; 
  TH1D* jtptPPMC_FULLBayes_Frac_h[nCentBins][nAbsEtaBins+1][nSystForFill][nBayesIter];
  RooUnfoldBayes* rooUnfoldBayesPPMC_FULLBayes_Frac_h[nCentBins][nAbsEtaBins+1][nSystForFill][nBayesIter];
  TH1D* jtptPPMC_FULLSvd_Frac_h[nCentBins][nAbsEtaBins+1][nSystForFill][nSvd];
  TH1D* refptPPMC_Frac_h[nCentBins][nAbsEtaBins+1];
  TH2D* responsePPMC_h[nCentBins][nAbsEtaBins+1];
  RooUnfoldResponse* rooResponsePPMC_h[nCentBins][nAbsEtaBins+1];
  TH1D* jtetaPPMC_FULLBayes_h[nCentBins][nSystForFill][nBayesIter];
  TH1D* jtetaPPMC_FULLBayes_Frac_h[nCentBins][nSystForFill][nBayesIter];
  
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    const std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
    
    for(Int_t aI = 0; aI < nAbsEtaBins+1; ++aI){
      std::string absEtaStr = "AbsEtaTot";
      if(aI < nAbsEtaBins) absEtaStr = "AbsEta" +absEtaBinsLow[aI] + "to" + absEtaBinsHi[aI];
      const std::string binStrPbPb = algoStrPbPb + "_" + centStr + "_" + absEtaStr;
      
      for(Int_t sI = 0; sI < nSystForFill; ++sI){
	jtptPbPbData_RAW_h[cI][aI][systPosForFill[sI]] = (TH1D*)inFile_p->Get((algoStrPbPb + "/" + centStr + "/jtptPbPbData_" + binStrPbPb + "_RAW_" + systStr[systPosForFill[sI]] + "_h").c_str());
	jtptPbPbMC_RAW_h[cI][aI][systPosForFill[sI]] = (TH1D*)inFile_p->Get((algoStrPbPb + "/" + centStr + "/jtptPbPbMC_" + binStrPbPb + "_RAW_" + systStr[systPosForFill[sI]] + "_h").c_str());
	jtptPbPbMC_RAW_Frac_h[cI][aI][systPosForFill[sI]] = (TH1D*)inFile_p->Get((algoStrPbPb + "/" + centStr + "/jtptPbPbMC_" + binStrPbPb + "_RAW_" + systStr[systPosForFill[sI]] + "_Frac_h").c_str());
	
	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	
	for(Int_t bI = 0; bI < nBayesIter; ++bI){
	  jtptPbPbData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI] = (TH1D*)inFile_p->Get((algoStrPbPb + "/" + centStr + "/jtptPbPbData_" + binStrPbPb + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_h").c_str());
	  
	  jtptPbPbMC_FULLBayes_h[cI][aI][systPosForFill[sI]][bI] = (TH1D*)inFile_p->Get((algoStrPbPb + "/" + centStr + "/jtptPbPbMC_" + binStrPbPb + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_h").c_str());
	  jtptPbPbMC_FULLBayes_Frac_h[cI][aI][systPosForFill[sI]][bI] = (TH1D*)inFile_p->Get((algoStrPbPb + "/" + centStr + "/jtptPbPbMC_" + binStrPbPb + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_Frac_h").c_str());
	  
	  rooUnfoldBayesPbPbData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI] = (RooUnfoldBayes*)inFile_p->Get((algoStrPbPb + "/" + centStr +"/rooUnfoldBayesPbPbData_" + binStrPbPb + "_FullBayes" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_h").c_str());
	  
	}
	
	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	
	for(Int_t bI = 0; bI < nSvd; ++bI){
	  if(bI+1 > truthNBins[cI]) break;
	  
	  jtptPbPbData_FULLSvd_h[cI][aI][systPosForFill[sI]][bI] = (TH1D*)inFile_p->Get((algoStrPbPb + "/" + centStr + "/jtptPbPbData_" + binStrPbPb + "_FULLSvd" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_h").c_str());
	  jtptPbPbMC_FULLSvd_h[cI][aI][systPosForFill[sI]][bI] = (TH1D*)inFile_p->Get((algoStrPbPb + "/" + centStr + "/jtptPbPbMC_" + binStrPbPb + "_FULLSvd" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_h").c_str());
	  jtptPbPbMC_FULLSvd_Frac_h[cI][aI][systPosForFill[sI]][bI] = (TH1D*)inFile_p->Get((algoStrPbPb + "/" + centStr + "/jtptPbPbMC_" + binStrPbPb + "_FULLSvd" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_Frac_h").c_str());
	}
      }
      
      refptPbPbMC_h[cI][aI] = (TH1D*)inFile_p->Get((algoStrPbPb + "/" + centStr + "/refptPbPbMC_" + binStrPbPb + "_h").c_str());
      refptPbPbMC_Frac_h[cI][aI] = (TH1D*)inFile_p->Get((algoStrPbPb + "/" + centStr + "/refptPbPbMC_" + binStrPbPb + "_Frac_h").c_str());
      
      responsePbPbMC_h[cI][aI] = (TH2D*)inFile_p->Get((algoStrPbPb + "/" + centStr + "/responsePbPbMC_" + binStrPbPb + "_h").c_str());
      rooResponsePbPbMC_h[cI][aI] = (RooUnfoldResponse*)inFile_p->Get((algoStrPbPb + "/" + centStr + "/rooResponsePbPbMC_" + binStrPbPb + "_h").c_str());
    }
    
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    
    for(Int_t sI = 0; sI < nSystForFill; ++sI){
      for(Int_t bI = 0; bI < nBayesIter; ++bI){
	jtetaPbPbData_FULLBayes_h[cI][systPosForFill[sI]][bI] = (TH1D*)inFile_p->Get((algoStrPbPb + "/" + centStr + "/jtetaPbPbData_" + algoStrPbPb + "_" + centStr + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_h").c_str());
	
	jtetaPbPbMC_FULLBayes_h[cI][systPosForFill[sI]][bI] = (TH1D*)inFile_p->Get((algoStrPbPb + "/" + centStr + "/jtetaPbPbMC_" + algoStrPbPb + "_" + centStr + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_h").c_str());
	
	jtetaPbPbMC_FULLBayes_Frac_h[cI][systPosForFill[sI]][bI] = (TH1D*)inFile_p->Get((algoStrPbPb + "/" + centStr + "/jtetaPbPbMC_" + algoStrPbPb + "_" + centStr + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_Frac_h").c_str());
      }
    }
    
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    
    for(Int_t aI = 0; aI < nAbsEtaBins+1; ++aI){	
      std::string absEtaStr = "AbsEtaTot";
      if(aI < nAbsEtaBins) absEtaStr = "AbsEta" +absEtaBinsLow[aI] + "to" + absEtaBinsHi[aI];
      const std::string binStrPP = algoStrPP + "_" + centStr + "_" + absEtaStr;
      
      for(Int_t sI = 0; sI < nSystForFill; ++sI){
	jtptPPData_RAW_h[cI][aI][systPosForFill[sI]] = (TH1D*)inFile_p->Get((algoStrPP + "/" + centStr + "/jtptPPData_" + binStrPP + "_RAW_" + systStr[systPosForFill[sI]] + "_h").c_str());
	jtptPPMC_RAW_h[cI][aI][systPosForFill[sI]] = (TH1D*)inFile_p->Get((algoStrPP + "/" + centStr + "/jtptPPMC_" + binStrPP + "_RAW_" + systStr[systPosForFill[sI]] + "_h").c_str());
	jtptPPMC_RAW_Frac_h[cI][aI][systPosForFill[sI]] = (TH1D*)inFile_p->Get((algoStrPP + "/" + centStr + "/jtptPPMC_" + binStrPP + "_RAW_" + systStr[systPosForFill[sI]] + "_Frac_h").c_str());
	
	for(Int_t bI = 0; bI < nBayesIter; ++bI){
	  jtptPPData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI] = (TH1D*)inFile_p->Get((algoStrPP + "/" + centStr + "/jtptPPData_" + binStrPP + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_h").c_str());

	  jtptPPMC_FULLBayes_h[cI][aI][systPosForFill[sI]][bI] = (TH1D*)inFile_p->Get((algoStrPP + "/" + centStr + "/jtptPPMC_" + binStrPP + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_h").c_str());
	  jtptPPMC_FULLBayes_Frac_h[cI][aI][systPosForFill[sI]][bI] = (TH1D*)inFile_p->Get((algoStrPP + "/" + centStr + "/jtptPPMC_" + binStrPP + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_Frac_h").c_str());
	  
	  rooUnfoldBayesPPData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI] = (RooUnfoldBayes*)inFile_p->Get((algoStrPP + "/" + centStr +"/rooUnfoldBayesPPData_" + binStrPP + "_FullBayes" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_h").c_str());
	}
	
	for(Int_t bI = 0; bI < nSvd; ++bI){
	  if(bI+1 > truthNBins[cI]) break;
	  
	  jtptPPData_FULLSvd_h[cI][aI][systPosForFill[sI]][bI] = (TH1D*)inFile_p->Get((algoStrPP + "/" + centStr + "/jtptPPData_" + binStrPP + "_FULLSvd" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_h").c_str());
	  jtptPPMC_FULLSvd_h[cI][aI][systPosForFill[sI]][bI] = (TH1D*)inFile_p->Get((algoStrPP + "/" + centStr + "/jtptPPMC_" + binStrPP + "_FULLSvd" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_h").c_str());
	  jtptPPMC_FULLSvd_Frac_h[cI][aI][systPosForFill[sI]][bI] = (TH1D*)inFile_p->Get((algoStrPP + "/" + centStr + "/jtptPPMC_" + binStrPP + "_FULLSvd" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_Frac_h").c_str());
	}
      }
      
      refptPPMC_h[cI][aI] = (TH1D*)inFile_p->Get((algoStrPP + "/" + centStr + "/refptPPMC_" + binStrPP + "_h").c_str());
      refptPPMC_Frac_h[cI][aI] = (TH1D*)inFile_p->Get((algoStrPP + "/" + centStr + "/refptPPMC_" + binStrPP + "_Frac_h").c_str());
      responsePPMC_h[cI][aI] = (TH2D*)inFile_p->Get((algoStrPP + "/" + centStr + "/responsePPMC_" + binStrPP + "_h").c_str());
      rooResponsePPMC_h[cI][aI] = (RooUnfoldResponse*)inFile_p->Get((algoStrPP + "/" + centStr + "/rooResponsePPMC_" + binStrPP + "_h").c_str());	
    }
    
    for(Int_t sI = 0; sI < nSystForFill; ++sI){
      for(Int_t bI = 0; bI < nBayesIter; ++bI){
	jtetaPPData_FULLBayes_h[cI][systPosForFill[sI]][bI] = (TH1D*)inFile_p->Get((algoStrPP + "/" + centStr + "/jtetaPPData_" + algoStrPP + "_" + centStr + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_h").c_str());
	
	jtetaPPMC_FULLBayes_h[cI][systPosForFill[sI]][bI] = (TH1D*)inFile_p->Get((algoStrPP + "/" + centStr + "/jtetaPPMC_" + algoStrPP + "_" + centStr + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_h").c_str());
	
	jtetaPPMC_FULLBayes_Frac_h[cI][systPosForFill[sI]][bI] = (TH1D*)inFile_p->Get((algoStrPP + "/" + centStr + "/jtetaPPMC_" + algoStrPP + "_" + centStr + "_FULLBayes" + std::to_string(bI+1) + "_" + systStr[systPosForFill[sI]] + "_Frac_h").c_str());
      }
    }
  }
    
  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  /*  
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    const std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
    //      const std::string dirStrPbPb = algoStrPbPb + "/" + centStr;
    //      const std::string dirStrPP = algoStrPP + "/" + centStr;
    
    for(Int_t aI = 0; aI < nAbsEtaBins; ++aI){
      std::string absEtaStr = "AbsEta" +absEtaBinsLow[aI] + "to" + absEtaBinsHi[aI];
      if(absEtaStr.find("AbsEta0p0to2p0") != std::string::npos) continue;

      for(Int_t sI = 0; sI < nSystForFillForCopy; ++sI){
	for(Int_t bI = 0; bI < nBayesIter; ++bI){
	  doHistCopy(jtptPbPbData_FULLBayes_h[cI][aI][0][bI], jtptPbPbData_FULLBayes_h[cI][aI][systPosForCopy[sI]][bI]);
	  doHistCopy(jtptPPData_FULLBayes_h[cI][aI][0][bI], jtptPPData_FULLBayes_h[cI][aI][systPosForCopy[sI]][bI]);
	  doHistCopy(jtptPbPbMC_FULLBayes_h[cI][aI][0][bI], jtptPbPbMC_FULLBayes_h[cI][aI][systPosForCopy[sI]][bI]);
	  doHistCopy(jtptPPMC_FULLBayes_h[cI][aI][0][bI], jtptPPMC_FULLBayes_h[cI][aI][systPosForCopy[sI]][bI]);
	  doHistCopy(jtptPbPbMC_FULLBayes_Frac_h[cI][aI][0][bI], jtptPbPbMC_FULLBayes_Frac_h[cI][aI][systPosForCopy[sI]][bI]);
	  doHistCopy(jtptPPMC_FULLBayes_Frac_h[cI][aI][0][bI], jtptPPMC_FULLBayes_Frac_h[cI][aI][systPosForCopy[sI]][bI]);
	}
	if(doSvd){
	  for(Int_t bI = 0; bI < nSvd; ++bI){
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
  */

  std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    std::cout << "Validating R=" << rVals << ", Cent=" << centBinsLow[cI] << "-" << centBinsHi[cI] << std::endl;
    
    for(Int_t aI = 0; aI < nAbsEtaBins; ++aI){
      for(Int_t sI = 0; sI < nSystForFill; ++sI){
	for(Int_t bI = 0; bI < nBayesIter; ++bI){
	  
	  TH1D* prevHist_p = NULL;
	  if(bI != 0) prevHist_p = jtptPbPbData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI-1];

	  RooUnfoldBayes* tempBayes = NULL;
	  if(aI < nAbsEtaBins) tempBayes = rooUnfoldBayesPbPbData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI];

	  plotValidation(canvPbPbData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI], padsPbPbData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI], {jtptPbPbMC_FULLBayes_h[cI][aI][systPosForFill[sI]][bI], refptPbPbMC_h[cI][aI], jtptPbPbMC_FULLBayes_Frac_h[cI][aI][systPosForFill[sI]][bI], refptPbPbMC_Frac_h[cI][aI], jtptPbPbData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI], prevHist_p, jtptPbPbData_RAW_h[cI][aI][systPosForFill[sI]]}, rooResponsePbPbMC_h[cI][aI], tempBayes, NULL, true, bI, "PbPb");

	  prevHist_p = NULL;
	  if(bI != 0) prevHist_p = jtptPPData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI-1];
	  plotValidation(canvPPData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI], padsPPData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI], {jtptPPMC_FULLBayes_h[cI][aI][systPosForFill[sI]][bI], refptPPMC_h[cI][aI], jtptPPMC_FULLBayes_Frac_h[cI][aI][systPosForFill[sI]][bI], refptPPMC_Frac_h[cI][aI], jtptPPData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI], prevHist_p, jtptPPData_RAW_h[cI][aI][systPosForFill[sI]]}, rooResponsePPMC_h[cI][aI], tempBayes, NULL, true, bI, "PP");
	}
      }
    }
  }

  

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    comboValidation(refptPbPbMC_h[cI][nAbsEtaBins-1], refptPbPbMC_h[cI][nAbsEtaBins]);
    comboValidation(refptPbPbMC_Frac_h[cI][nAbsEtaBins-1], refptPbPbMC_Frac_h[cI][nAbsEtaBins]);
    
    comboValidation(refptPPMC_h[cI][nAbsEtaBins-1], refptPPMC_h[cI][nAbsEtaBins]);
    comboValidation(refptPPMC_Frac_h[cI][nAbsEtaBins-1], refptPPMC_Frac_h[cI][nAbsEtaBins]);

    
    for(Int_t bI = 0; bI < nBayesIter; ++bI){
      comboValidation(jtptPbPbData_FULLBayes_h[cI][nAbsEtaBins-1][0][bI], jtptPbPbData_FULLBayes_h[cI][nAbsEtaBins][0][bI]);      
      comboValidation(jtptPPData_FULLBayes_h[cI][nAbsEtaBins-1][0][bI], jtptPPData_FULLBayes_h[cI][nAbsEtaBins][0][bI]);      

      comboValidation(jtptPbPbMC_FULLBayes_h[cI][nAbsEtaBins-1][0][bI], jtptPbPbMC_FULLBayes_h[cI][nAbsEtaBins][0][bI]);      
      comboValidation(jtptPPMC_FULLBayes_h[cI][nAbsEtaBins-1][0][bI], jtptPPMC_FULLBayes_h[cI][nAbsEtaBins][0][bI]);      

      comboValidation(jtptPbPbMC_FULLBayes_Frac_h[cI][nAbsEtaBins-1][0][bI], jtptPbPbMC_FULLBayes_Frac_h[cI][nAbsEtaBins][0][bI]);      
      comboValidation(jtptPPMC_FULLBayes_Frac_h[cI][nAbsEtaBins-1][0][bI], jtptPPMC_FULLBayes_Frac_h[cI][nAbsEtaBins][0][bI]);      
    }    
  }
     

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    for(Int_t aI = 0; aI < nAbsEtaBins+1; ++aI){
      for(Int_t sI = 0; sI < nSystForFill; ++sI){
	for(Int_t bI = 0; bI < nBayesIter; ++bI){
	  std::string saveName = canvPbPbData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI]->GetName();
	  saveName = dirName2 + saveName + std::to_string(date->GetDate()) + ".pdf";
	  
	  for(Int_t pI = 0; pI < nPadTotValid; ++pI){
	    canvPbPbData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI]->cd();
	    padsPbPbData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI][pI]->cd();
	    
	    gStyle->SetOptStat(0);
	  }
	  
	  canvPbPbData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI]->SaveAs(saveName.c_str());

	  saveName = canvPPData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI]->GetName();
	  saveName = dirName2 + saveName + std::to_string(date->GetDate()) + ".pdf";
	  
	  for(Int_t pI = 0; pI < nPadTotValid; ++pI){
	    canvPPData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI]->cd();
	    padsPPData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI][pI]->cd();
	    
	    gStyle->SetOptStat(0);
	  }
	  canvPPData_FULLBayes_h[cI][aI][systPosForFill[sI]][bI]->SaveAs(saveName.c_str());
	}
      }       
    }
  }

  std::cout << "FIRST ROUND DONE" << std::endl;

  //  inFile_p->Close();
  //  delete inFile_p;

  delete randGen_p;
  delete date;
  std::cout << "DATE DELETED" << std::endl;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./makePlotValidation_FromTree.exe <inFileName>" << std::endl;
    return 1;
  }


  int retVal = 0;
  retVal += makePlotValidation_FromTree(argv[1]);
  std::cout << "LAST VALUE RETURNED" << std::endl;
  return retVal;
}
