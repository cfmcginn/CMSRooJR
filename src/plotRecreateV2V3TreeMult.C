#include <iostream>
#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TDatime.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"

#include "include/histDefUtility.h"
#include "include/kirchnerPalette.h"
#include "include/vanGoghPalette.h"
#include "CustomCanvas.h"
#include "TexSlides.C"


void setColorStyleLabelTitle(TH1F* inHist_p, const Int_t color, const Int_t style, const Int_t labelSize, const Int_t titleSize)
{
  inHist_p->SetLineColor(color);
  inHist_p->SetMarkerColor(color);
  inHist_p->SetMarkerStyle(style);

  inHist_p->GetXaxis()->SetTitleFont(43);
  inHist_p->GetYaxis()->SetTitleFont(43);
  inHist_p->GetXaxis()->SetLabelFont(43);
  inHist_p->GetYaxis()->SetLabelFont(43);

  inHist_p->GetXaxis()->SetTitleSize(titleSize);
  inHist_p->GetYaxis()->SetTitleSize(titleSize);
  inHist_p->GetXaxis()->SetLabelSize(labelSize);
  inHist_p->GetYaxis()->SetLabelSize(labelSize);

  return;
}


int plotRecreateV2V3(const std::string inFileName, const std::string inJamesFileName)
{
  vanGoghPalette vg;
  const Int_t flowStyleSet[4] = {33, 34, 21, 22};
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  TLatex* label_PF_p = new TLatex();
  label_PF_p->SetNDC();
  label_PF_p->SetTextFont(43);
  label_PF_p->SetTextSize(18);
  /*
  TLatex* label_Trk_p = new TLatex();
  label_Trk_p->SetNDC();
  label_Trk_p->SetTextFont(43);
  label_Trk_p->SetTextFont(18);
  */
  TLegend* leg_PF_p = new TLegend(0.6, 0.55, 0.9, 0.85);
  leg_PF_p->SetTextFont(43);
  leg_PF_p->SetTextSize(18);
  leg_PF_p->SetBorderSize(0);
  leg_PF_p->SetFillColor(0);
  leg_PF_p->SetFillStyle(0);
  /*
  TLegend* leg_Trk_p = new TLegend(0.6, 0.55, 0.9, 0.85);
  leg_Trk_p->SetTextFont(43);
  leg_Trk_p->SetTextSize(18);
  leg_Trk_p->SetBorderSize(0);
  leg_Trk_p->SetFillColor(0);
  leg_Trk_p->SetFillStyle(0);
  */
  const Int_t nLowCuts = 7;
  Int_t lowCuts[nLowCuts] = {0,100,150,200,250,300,400};
  const Int_t nHighCuts = 6;
  Int_t highCuts[nHighCuts] = {100000000,1600,1500,1400,1300,1200};

  const Int_t nCentBins = 11;
  const Int_t centBinsLow[nCentBins] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55};
  const Int_t centBinsHi[nCentBins] = {10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60};
  //  const Double_t centBins[nCentBins+1] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60};

  TFile* jamesFile_p = new TFile(inJamesFileName.c_str(), "READ");
  TH1F* jamesHist_p[nCentBins];
  //  TH1F* jamesMean_p = new TH1F("jamesMean_p", ";Centrality (%);#LTv_{2}^{obs}#GT", nCentBins, centBins);
  //  TH1F* jamesSigma_p = new TH1F("jamesSigma_p", ";Centrality (%);#sigma(v_{2}^{obs})", nCentBins, centBins);

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    jamesHist_p[cI] = (TH1F*)jamesFile_p->Get(("qwebye/hVnFull_c" + std::to_string(cI+1)).c_str());

    //    jamesMean_p->SetBinContent(cI+1, jamesHist_p[cI]->GetMean());
    //    jamesMean_p->SetBinError(cI+1, jamesHist_p[cI]->GetMeanError());

    //    jamesSigma_p->SetBinContent(cI+1, jamesHist_p[cI]->GetStdDev());
    //    jamesSigma_p->SetBinError(cI+1, jamesHist_p[cI]->GetStdDevError());

    centerTitles(jamesHist_p[cI]);
    setSumW2(jamesHist_p[cI]);


    jamesHist_p[cI]->Scale(1./jamesHist_p[cI]->Integral());
    jamesHist_p[cI]->SetMaximum(jamesHist_p[cI]->GetMaximum()*1.5);

    setColorStyleLabelTitle(jamesHist_p[cI], 1, 20, 12, 14);
    jamesHist_p[cI]->GetYaxis()->SetTitleOffset(jamesHist_p[cI]->GetYaxis()->GetTitleOffset()*1.5);
    jamesHist_p[cI]->SetTitle("");

    if(cI == 0){
      leg_PF_p->AddEntry(jamesHist_p[cI], "HIN-16-019", "P L");
      //      leg_Trk_p->AddEntry(jamesHist_p[cI], "HIN-16-019", "P L");
    }
  }

  //  centerTitles({jamesMean_p, jamesSigma_p});
  //  setSumW2({jamesMean_p, jamesSigma_p});
  //  jamesMean_p->SetMaximum(jamesMean_p->GetMaximum()*1.5);
  //  jamesMean_p->SetMinimum(0.0);

  //  setColorStyleLabelTitle(jamesMean_p, 1, 20, 12, 14);
  //  jamesMean_p->GetYaxis()->SetTitleOffset(jamesHist_p[0]->GetYaxis()->GetTitleOffset());
  
  //  jamesSigma_p->SetMaximum(jamesSigma_p->GetMaximum()*1.5);
  //  jamesSigma_p->SetMinimum(0.0);
  //  setColorStyleLabelTitle(jamesSigma_p, 1, 20, 12, 14);
  //  jamesSigma_p->GetYaxis()->SetTitleOffset(jamesHist_p[0]->GetYaxis()->GetTitleOffset());

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TH1F* v2Raw_PF_h[nLowCuts][nHighCuts][nCentBins];
  TH1F* v2RawCorr_PF_h[nLowCuts][nHighCuts][nCentBins];
  //  TH1F* v2ObsCorr_PF_h[nLowCut][nHighCuts][nCentBins];
  TH1F* v2Fit_PF_h[nLowCuts][nHighCuts][nCentBins];
  //  TH1F* v2Diff_PF_h[nCentBins];
  //  TH1F* v2DiffCorr_PF_h[nCentBins];
  //  TH1F* v2Comp_h[nCentBins];
  //  TH1F* v2FT_h[nCentBins];
  //  TH1F* v2Diff_h[nCentBins];

  //  TH1F* v2Raw_Trk_h[nCentBins];
  //  TH1F* v2RawCorr_Trk_h[nCentBins];
  //  TH1F* v2ObsCorr_Trk_h[nCentBins];

  //  TH1F* v2Raw_Mean_PF_h = (TH1F*)inFile_p->Get("v2Raw_Mean_PF_h");
  //  TH1F* v2Raw_Sigma_PF_h = (TH1F*)inFile_p->Get("v2Raw_Sigma_PF_h");

  //  TH1F* v2RawCorr_Mean_PF_h = (TH1F*)inFile_p->Get("v2RawCorr_Mean_PF_h");
  //  TH1F* v2RawCorr_Sigma_PF_h = (TH1F*)inFile_p->Get("v2RawCorr_Sigma_PF_h");

  //  TH1F* v2Fit_Mean_PF_h = (TH1F*)inFile_p->Get("v2Fit_Mean_PF_h");
  //  TH1F* v2Fit_Sigma_PF_h = (TH1F*)inFile_p->Get("v2Fit_Sigma_PF_h");

  //  TH1F* v2ObsCorr_Mean_PF_h = (TH1F*)inFile_p->Get("v2ObsCorr_Mean_PF_h");
  //  TH1F* v2ObsCorr_Sigma_PF_h = (TH1F*)inFile_p->Get("v2ObsCorr_Sigma_PF_h");

  ///  TH1F* v2Raw_Mean_Trk_h = (TH1F*)inFile_p->Get("v2Raw_Mean_Trk_h");
  //  TH1F* v2Raw_Sigma_Trk_h = (TH1F*)inFile_p->Get("v2Raw_Sigma_Trk_h");

  //  TH1F* v2RawCorr_Mean_Trk_h = (TH1F*)inFile_p->Get("v2RawCorr_Mean_Trk_h");
  //  TH1F* v2RawCorr_Sigma_Trk_h = (TH1F*)inFile_p->Get("v2RawCorr_Sigma_Trk_h");

  //  TH1F* v2ObsCorr_Mean_Trk_h = (TH1F*)inFile_p->Get("v2ObsCorr_Mean_Trk_h");
  //  TH1F* v2ObsCorr_Sigma_Trk_h = (TH1F*)inFile_p->Get("v2ObsCorr_Sigma_Trk_h");

  //  TH1F* size_h = (TH1F*)inFile_p->Get("size_h");
  //  TH1F* size_rawPlus_h = (TH1F*)inFile_p->Get("size_rawPlus_h");
  //  TH1F* size_rawMinus_h = (TH1F*)inFile_p->Get("size_rawMinus_h");
  //  TH1F* size_rawCorrPlus_h = (TH1F*)inFile_p->Get("size_rawCorrPlus_h");
  //  TH1F* size_rawCorrMinus_h = (TH1F*)inFile_p->Get("size_rawCorrMinus_h");

  for(Int_t nL = 0; nL < nLowCuts; nL++){
    const std::string lowString = (nL==0)?"":std::to_string(lowCuts[nL]);
    for (Int_t nH = 0; nH < nHighCuts; nH++){
      const std::string highString = (nH==0)?"":std::to_string(highCuts[nH]);
      for(Int_t cI = 0; cI < nCentBins; ++cI){
	const std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
	v2Raw_PF_h[nL][nH][cI] = (TH1F*)inFile_p->Get(("v2Raw_" + lowString + "_" + highString + "_" + centStr + "_PF_h").c_str());
	v2RawCorr_PF_h[nL][nH][cI] = (TH1F*)inFile_p->Get(("v2RawCorr_" + lowString + "_" + highString + "_" + centStr + "_PF_h").c_str());
	//    v2ObsCorr_PF_h[cI] = (TH1F*)inFile_p->Get(("v2ObsCorr_" + centStr + "_PF_h").c_str());
	v2Fit_PF_h[nL][nH][cI] = (TH1F*)inFile_p->Get(("v2Fit_" + lowString + "_" + highString + "_" + centStr + "_PF_h").c_str());
	//    v2Diff_PF_h[cI] = (TH1F*)inFile_p->Get(("v2Diff_" + centStr + "_PF_h").c_str());
	//    v2DiffCorr_PF_h[cI] = (TH1F*)inFile_p->Get(("v2DiffCorr_" + centStr + "_PF_h").c_str());
    //    v2Comp_h[cI] = (TH1F*)inFile_p->Get(("v2Comp_" + centStr + "_h").c_str());
    //    v2FT_h[cI] = (TH1F*)inFile_p->Get(("v2FT_" + centStr + "_h").c_str());
    //    v2Diff_h[cI] = (TH1F*)inFile_p->Get(("v2Diff_" + centStr + "_h").c_str());
    
    //    v2Raw_Trk_h[cI] = (TH1F*)inFile_p->Get(("v2Raw_" + centStr + "_Trk_h").c_str());
    //    v2RawCorr_Trk_h[cI] = (TH1F*)inFile_p->Get(("v2RawCorr_" + centStr + "_Trk_h").c_str());
    //    v2ObsCorr_Trk_h[cI] = (TH1F*)inFile_p->Get(("v2ObsCorr_" + centStr + "_Trk_h").c_str());

    //    centerTitles({v2Raw_PF_h[cI], v2RawCorr_PF_h[cI], v2ObsCorr_PF_h[cI], v2Comp_h[cI], v2FT_h[cI], v2Diff_h[cI], v2Raw_Trk_h[cI], v2RawCorr_Trk_h[cI], v2ObsCorr_Trk_h[cI]});
	centerTitles({v2Raw_PF_h[nL][nH][cI], v2RawCorr_PF_h[nL][nH][cI], v2Fit_PF_h[nL][nH][cI]});
    //    setSumW2({v2Raw_PF_h[cI], v2RawCorr_PF_h[cI], v2ObsCorr_PF_h[cI], v2Comp_h[cI], v2FT_h[cI], v2Diff_h[cI], v2Raw_Trk_h[cI], v2RawCorr_Trk_h[cI], v2ObsCorr_Trk_h[cI]});
	setSumW2({v2Raw_PF_h[nL][nH][cI], v2RawCorr_PF_h[nL][nH][cI], v2Fit_PF_h[nL][nH][cI]});

	setColorStyleLabelTitle(v2Raw_PF_h[nL][nH][cI], vg.getColor(0), flowStyleSet[0], 12, 14);
	setColorStyleLabelTitle(v2RawCorr_PF_h[nL][nH][cI], vg.getColor(1), flowStyleSet[1], 12, 14);
    //    setColorStyleLabelTitle(v2ObsCorr_PF_h[cI], vg.getColor(2), flowStyleSet[2], 12, 14);
	setColorStyleLabelTitle(v2Fit_PF_h[nL][nH][cI], vg.getColor(2), flowStyleSet[2], 12, 14);
    //    setColorStyleLabelTitle(v2Diff_PF_h[cI], vg.getColor(0), flowStyleSet[0], 12, 14);
    //    setColorStyleLabelTitle(v2DiffCorr_PF_h[cI], vg.getColor(0), flowStyleSet[0], 12, 14);
    //    setColorStyleLabelTitle(v2Comp_h[cI], vg.getColor(0), flowStyleSet[0], 12, 14);
    //    setColorStyleLabelTitle(v2FT_h[cI], vg.getColor(0), flowStyleSet[0], 12, 14);
    //    setColorStyleLabelTitle(v2Diff_h[cI], vg.getColor(0), flowStyleSet[0], 12, 14);

    //    setColorStyleLabelTitle(v2Raw_Trk_h[cI], vg.getColor(0), flowStyleSet[0], 12, 14);
    //    setColorStyleLabelTitle(v2RawCorr_Trk_h[cI], vg.getColor(1), flowStyleSet[1], 12, 14);
    //    setColorStyleLabelTitle(v2ObsCorr_Trk_h[cI], vg.getColor(2), flowStyleSet[2], 12, 14);
	
	v2Raw_PF_h[nL][nH][cI]->Scale(1./v2Raw_PF_h[nL][nH][cI]->Integral());
	v2RawCorr_PF_h[nL][nH][cI]->Scale(1./v2RawCorr_PF_h[nL][nH][cI]->Integral());
    //    v2ObsCorr_PF_h[cI]->Scale(1./v2ObsCorr_PF_h[cI]->Integral());
	v2Fit_PF_h[nL][nH][cI]->Scale(1./v2Fit_PF_h[nL][nH][cI]->Integral());

    //    v2Raw_Trk_h[cI]->Scale(1./v2Raw_Trk_h[cI]->Integral());
    //    v2RawCorr_Trk_h[cI]->Scale(1./v2RawCorr_Trk_h[cI]->Integral());
    //    v2ObsCorr_Trk_h[cI]->Scale(1./v2ObsCorr_Trk_h[cI]->Integral());

	v2Raw_PF_h[nL][nH][cI]->GetYaxis()->SetTitleOffset(jamesHist_p[cI]->GetYaxis()->GetTitleOffset());
	v2Raw_PF_h[nL][nH][cI]->GetXaxis()->SetTitleOffset(jamesHist_p[cI]->GetYaxis()->GetTitleOffset()*1.5);       

    //    v2Raw_Trk_h[cI]->GetYaxis()->SetTitleOffset(jamesHist_p[cI]->GetYaxis()->GetTitleOffset());
    //    v2Raw_Trk_h[cI]->GetXaxis()->SetTitleOffset(jamesHist_p[cI]->GetYaxis()->GetTitleOffset()*1.5);

	if(cI+nL+nH == 0){
	  leg_PF_p->AddEntry(v2Raw_PF_h[nL][nH][cI], "Jet Extraction", "P L");
	  leg_PF_p->AddEntry(v2RawCorr_PF_h[nL][nH][cI], "Eff. Corr.", "P L");
	  leg_PF_p->AddEntry(v2Fit_PF_h[nL][nH][cI], "Fit", "P L");
	  
      //      leg_Trk_p->AddEntry(v2Raw_Trk_h[cI], "Jet Extraction", "P L");
      //      leg_Trk_p->AddEntry(v2RawCorr_Trk_h[cI], "Eff. Corr.", "P L");
      //      leg_Trk_p->AddEntry(v2ObsCorr_Trk_h[cI], "Eff.,Det.,Trk.", "P L");
	}
      }
    }
  }

      //  centerTitles({v2Raw_Mean_PF_h, v2Raw_Sigma_PF_h, v2RawCorr_Mean_PF_h, v2RawCorr_Sigma_PF_h, v2Fit_Mean_PF_h, v2Fit_Sigma_PF_h, v2ObsCorr_Mean_PF_h, v2ObsCorr_Sigma_PF_h});//, v2Raw_Mean_Trk_h, v2Raw_Sigma_Trk_h, v2RawCorr_Mean_Trk_h, v2RawCorr_Sigma_Trk_h, v2ObsCorr_Mean_Trk_h, v2ObsCorr_Sigma_Trk_h});
      //  setSumW2({v2Raw_Mean_PF_h, v2Raw_Sigma_PF_h, v2RawCorr_Mean_PF_h, v2RawCorr_Sigma_PF_h, v2Fit_Mean_PF_h, v2Fit_Sigma_PF_h, v2ObsCorr_Mean_PF_h, v2ObsCorr_Sigma_PF_h});//, v2Raw_Mean_Trk_h, v2Raw_Sigma_Trk_h, v2RawCorr_Mean_Trk_h, v2RawCorr_Sigma_Trk_h, v2ObsCorr_Mean_Trk_h, v2ObsCorr_Sigma_Trk_h});

  //  setColorStyleLabelTitle(v2Raw_Mean_PF_h, vg.getColor(0), flowStyleSet[0], 12, 14);
  //  setColorStyleLabelTitle(v2RawCorr_Mean_PF_h, vg.getColor(1), flowStyleSet[1], 12, 14);
  //  setColorStyleLabelTitle(v2Fit_Mean_PF_h, vg.getColor(2), flowStyleSet[2], 12, 14);
  //  setColorStyleLabelTitle(v2ObsCorr_Mean_PF_h, vg.getColor(2), flowStyleSet[2], 12, 14);

  //  setColorStyleLabelTitle(v2Raw_Mean_Trk_h, vg.getColor(0), flowStyleSet[0], 12, 14);
  //  setColorStyleLabelTitle(v2RawCorr_Mean_Trk_h, vg.getColor(1), flowStyleSet[1], 12, 14);
  //  setColorStyleLabelTitle(v2ObsCorr_Mean_Trk_h, vg.getColor(2), flowStyleSet[2], 12, 14);

  //  setColorStyleLabelTitle(v2Raw_Sigma_PF_h, vg.getColor(0), flowStyleSet[0], 12, 14);
  //  setColorStyleLabelTitle(v2RawCorr_Sigma_PF_h, vg.getColor(1), flowStyleSet[1], 12, 14);
  //  setColorStyleLabelTitle(v2Fit_Sigma_PF_h, vg.getColor(2), flowStyleSet[2], 12, 14);
  //  setColorStyleLabelTitle(v2ObsCorr_Sigma_PF_h, vg.getColor(2), flowStyleSet[2], 12, 14);

  //  setColorStyleLabelTitle(v2Raw_Sigma_Trk_h, vg.getColor(0), flowStyleSet[0], 12, 14);
  //  setColorStyleLabelTitle(v2RawCorr_Sigma_Trk_h, vg.getColor(1), flowStyleSet[1], 12, 14);
  //  setColorStyleLabelTitle(v2ObsCorr_Sigma_Trk_h, vg.getColor(2), flowStyleSet[2], 12, 14);

  gSystem->cd("pdfDir");

  CustomCanvas* canv_p = new CustomCanvas("canv_p", "", 450, 450);
  /*
  size_h->Draw();
  canv_p->SaveAs("size.pdf");
  size_rawPlus_h->Draw();
  canv_p->SaveAs("size_rawPlus.pdf");
  size_rawMinus_h->Draw();
  canv_p->SaveAs("size_rawMinus.pdf");
  size_rawCorrPlus_h->Draw();
  canv_p->SaveAs("size_rawCorrPlus.pdf");
  size_rawCorrMinus_h->Draw();
  canv_p->SaveAs("size_rawCorrMinus.pdf");
  */
  gStyle->SetOptStat(0);
  /*
  Float_t max=jamesMean_p->GetMaximum();
  if (v2Raw_Mean_PF_h->GetMaximum()>max) max=v2Raw_Mean_PF_h->GetMaximum();
  if (v2RawCorr_Mean_PF_h->GetMaximum()>max) max=v2RawCorr_Mean_PF_h->GetMaximum();
  if (v2ObsCorr_Mean_PF_h->GetMaximum()>max) max=v2ObsCorr_Mean_PF_h->GetMaximum();
  jamesMean_p->SetMaximum(1.5*max);
  jamesMean_p->SetTitle("Mean (PF)");
  jamesMean_p->Draw("HIST E1 P");
  v2Raw_Mean_PF_h->Draw("SAME HIST E1 P");
  v2RawCorr_Mean_PF_h->Draw("SAME HIST E1 P");
  v2Fit_Mean_PF_h->Draw("SAME HIST E1 P");
  leg_PF_p->Draw("SAME");
  canv_p->SaveAs("Mean_PF.pdf");
  
  max=jamesMean_p->GetMaximum();
  if (v2Raw_Mean_Trk_h->GetMaximum()>max) max=v2Raw_Mean_Trk_h->GetMaximum();
  if (v2RawCorr_Mean_Trk_h->GetMaximum()>max) max=v2RawCorr_Mean_Trk_h->GetMaximum();
  if (v2ObsCorr_Mean_Trk_h->GetMaximum()>max) max=v2ObsCorr_Mean_Trk_h->GetMaximum();
  jamesMean_p->SetMaximum(1.5*max);
  jamesMean_p->SetTitle("Mean (Trk)");
  jamesMean_p->Draw("HIST E1 P");
  v2Raw_Mean_Trk_h->Draw("SAME HIST E1 P");
  v2RawCorr_Mean_Trk_h->Draw("SAME HIST E1 P");
  v2ObsCorr_Mean_Trk_h->Draw("SAME HIST E1 P");
  leg_Trk_p->Draw("SAME");
  canv_p->SaveAs("Mean_Trk.pdf");
  
  max=jamesSigma_p->GetMaximum();
  if (v2Raw_Sigma_PF_h->GetMaximum()>max) max=v2Raw_Sigma_PF_h->GetMaximum();
  if (v2RawCorr_Sigma_PF_h->GetMaximum()>max) max=v2RawCorr_Sigma_PF_h->GetMaximum();
  if (v2ObsCorr_Sigma_PF_h->GetMaximum()>max) max=v2ObsCorr_Sigma_PF_h->GetMaximum();
  jamesSigma_p->SetMaximum(1.5*max);
  jamesSigma_p->SetTitle("Sigma (PF)");
  jamesSigma_p->Draw("HIST P");
  v2Raw_Sigma_PF_h->Draw("SAME HIST P");
  v2RawCorr_Sigma_PF_h->Draw("SAME HIST P");
  v2Fit_Sigma_PF_h->Draw("SAME HIST P");
  leg_PF_p->Draw("SAME");
  canv_p->SaveAs("Sigma_PF.pdf");
  
  max=jamesSigma_p->GetMaximum();
  if (v2Raw_Sigma_Trk_h->GetMaximum()>max) max=v2Raw_Sigma_Trk_h->GetMaximum();
  if (v2RawCorr_Sigma_Trk_h->GetMaximum()>max) max=v2RawCorr_Sigma_Trk_h->GetMaximum();
  if (v2ObsCorr_Sigma_Trk_h->GetMaximum()>max) max=v2ObsCorr_Sigma_Trk_h->GetMaximum();
  jamesSigma_p->SetMaximum(1.5*max);
  jamesSigma_p->SetTitle("Sigma (Trk)");
  jamesSigma_p->Draw("HIST P");
  v2Raw_Sigma_Trk_h->Draw("SAME HIST P");
  v2RawCorr_Sigma_Trk_h->Draw("SAME HIST P");
  v2ObsCorr_Sigma_Trk_h->Draw("SAME HIST P");
  leg_Trk_p->Draw("SAME");
  canv_p->SaveAs("Sigma_Trk.pdf");
  
  gStyle->SetOptStat(1111);
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    const std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
    const std::string saveName = "v2Diff_" + centStr + "_" + dateStr + "_PF.pdf";

    canv_p->SetTopMargin(0.1);
    canv_p->SetBottomMargin(0.1);
    canv_p->SetRightMargin(0.1);
    canv_p->SetLeftMargin(0.1);

    //    TH1F* tmp = new TH1F(("tmp" + centStr).c_str(), ("v2ObsCorr_PF/v2FromTree (" + centStr + ");v_{2};").c_str(), 150, 0.0, 0.6);
    //    v2FT_h[cI]->Divide(v2ObsCorr_PF_h[cI],v2FT_h[cI]);
    //    v2FT_h[cI]->SetTitle(("v2ObsCorr_PF/v2FromTree (" + centStr + ");v_{2};").c_str());
    //    v2FT_h[cI]->SetMaximum(3.0);
    v2Diff_PF_h[cI]->DrawCopy("HIST P");

    canv_p->SaveAs(saveName.c_str());
  }
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    const std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
    const std::string saveName = "v2DiffCorr_" + centStr + "_" + dateStr + "_PF.pdf";

    canv_p->SetTopMargin(0.1);
    canv_p->SetBottomMargin(0.1);
    canv_p->SetRightMargin(0.1);
    canv_p->SetLeftMargin(0.1);

    //    TH1F* tmp = new TH1F(("tmp" + centStr).c_str(), ("v2ObsCorr_PF/v2FromTree (" + centStr + ");v_{2};").c_str(), 150, 0.0, 0.6);
    //    v2FT_h[cI]->Divide(v2ObsCorr_PF_h[cI],v2FT_h[cI]);
    //    v2FT_h[cI]->SetTitle(("v2ObsCorr_PF/v2FromTree (" + centStr + ");v_{2};").c_str());
    //    v2FT_h[cI]->SetMaximum(3.0);
    v2DiffCorr_PF_h[cI]->DrawCopy("HIST P");

    canv_p->SaveAs(saveName.c_str());
  }
  */
  gStyle->SetOptStat(0);
  Float_t yPadBottomFrac = 0.35;
  
  for(Int_t nL = 0; nL < nLowCuts; nL++){
    for(Int_t nH = 0; nH < nHighCuts; nH++){
      for(Int_t cI = 0; cI < nCentBins; ++cI){
	//    TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
	canv_p->SetTopMargin(0.01);
	canv_p->SetBottomMargin(0.01);
	canv_p->SetRightMargin(0.01);
	canv_p->SetLeftMargin(0.01);

	TPad* pads[2];
	pads[0] = new TPad("pad0", "", 0.0, yPadBottomFrac, 1.0, 1.0);
	canv_p->cd();
	pads[0]->SetRightMargin(0.01);
	pads[0]->SetBottomMargin(0.001);
	pads[0]->SetLeftMargin(pads[0]->GetLeftMargin()*1.3);
	pads[0]->SetTopMargin(pads[0]->GetLeftMargin()/(1.0-yPadBottomFrac)*3./6.);
	pads[0]->Draw("SAME");
	pads[0]->cd();
    
	jamesHist_p[cI]->DrawCopy("HIST P");
	v2Raw_PF_h[nL][nH][cI]->DrawCopy("HIST P SAME");
	v2RawCorr_PF_h[nL][nH][cI]->DrawCopy("HIST P SAME");
	v2Fit_PF_h[nL][nH][cI]->DrawCopy("HIST P SAME");

	const std::string centStr2 = std::to_string(centBinsLow[cI]) + "-" + std::to_string(centBinsHi[cI]) + "%";
	std::string cutStr = ", Mult";
	if (nL==0){
	  if (nH==0) cutStr += " All";
	  else cutStr += "<" + std::to_string(highCuts[nH]);
	}
	else{
	  if (nH==0) cutStr += ">" + std::to_string(lowCuts[nL]);
	  else cutStr += " " + std::to_string(lowCuts[nL]) + "-" + std::to_string(highCuts[nH]);
	}

	label_PF_p->DrawLatex(0.15, 0.94, "#bf{CMS Preliminary}");
	label_PF_p->DrawLatex(0.20, 0.82, ("#bf{" + centStr2 + cutStr + "}").c_str());
	leg_PF_p->Draw("SAME");

	pads[1] = new TPad("pad1", "", 0.0, 0.0, 1.0, yPadBottomFrac);
	canv_p->cd();
	pads[1]->SetTopMargin(0.001);
	pads[1]->SetRightMargin(0.01);
	pads[1]->SetBottomMargin(pads[0]->GetLeftMargin()/yPadBottomFrac);
	pads[1]->SetLeftMargin(pads[0]->GetLeftMargin());
	pads[1]->Draw("SAME");
	pads[1]->cd();

	v2Raw_PF_h[nL][nH][cI]->Divide(jamesHist_p[cI]);
	v2RawCorr_PF_h[nL][nH][cI]->Divide(jamesHist_p[cI]);
	v2Fit_PF_h[nL][nH][cI]->Divide(jamesHist_p[cI]);
    
	v2Raw_PF_h[nL][nH][cI]->SetMinimum(0.55);
	v2Raw_PF_h[nL][nH][cI]->SetMaximum(1.45);
    
	v2Raw_PF_h[nL][nH][cI]->GetYaxis()->SetNdivisions(505);

	v2Raw_PF_h[nL][nH][cI]->DrawCopy("HIST P");
	v2RawCorr_PF_h[nL][nH][cI]->DrawCopy("SAME HIST P");
	v2Fit_PF_h[nL][nH][cI]->DrawCopy("SAME HIST P");
    
	gStyle->SetOptStat(0);

	const std::string lowString = (nL==0)?"":std::to_string(lowCuts[nL]);
	const std::string highString = (nH==0)?"":std::to_string(highCuts[nH]);
	const std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
	const std::string saveName = "plotRecreateV2James_" + lowString + "_" + highString + "_" + centStr + "_" + dateStr + "_PF.pdf";

	canv_p->SaveAs(saveName.c_str());

	delete pads[0];
	delete pads[1];
      }
    }
  }

  /*
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    const std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
    const std::string saveName = "v2FT_" + centStr + "_" + dateStr + "_PF.pdf";

    canv_p->SetTopMargin(0.1);
    canv_p->SetBottomMargin(0.1);
    canv_p->SetRightMargin(0.1);
    canv_p->SetLeftMargin(0.1);

    //    TH1F* tmp = new TH1F(("tmp" + centStr).c_str(), ("v2ObsCorr_PF/v2FromTree (" + centStr + ");v_{2};").c_str(), 150, 0.0, 0.6);
    v2FT_h[cI]->Divide(v2ObsCorr_PF_h[cI],v2FT_h[cI]);
    v2FT_h[cI]->SetTitle(("v2ObsCorr_PF/v2FromTree (" + centStr + ");v_{2};").c_str());
    v2FT_h[cI]->DrawCopy("HIST P");

    canv_p->SaveAs(saveName.c_str());
    }*/
  /*
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    //    TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
    canv_p->SetTopMargin(0.01);
    canv_p->SetBottomMargin(0.01);
    canv_p->SetRightMargin(0.01);
    canv_p->SetLeftMargin(0.01);

    TPad* pads[2];
    pads[0] = new TPad("pad0", "", 0.0, yPadBottomFrac, 1.0, 1.0);
    canv_p->cd();
    pads[0]->SetRightMargin(0.01);
    pads[0]->SetBottomMargin(0.001);
    pads[0]->SetLeftMargin(pads[0]->GetLeftMargin()*1.3);
    pads[0]->SetTopMargin(pads[0]->GetLeftMargin()/(1.0-yPadBottomFrac)*3./6.);
    pads[0]->Draw("SAME");
    pads[0]->cd();

    jamesHist_p[cI]->DrawCopy("HIST P");
    v2Raw_Trk_h[cI]->DrawCopy("HIST P SAME");
    v2RawCorr_Trk_h[cI]->DrawCopy("HIST P SAME");
    v2ObsCorr_Trk_h[cI]->DrawCopy("HIST P SAME");

    const std::string centStr2 = std::to_string(centBinsLow[cI]) + "-" + std::to_string(centBinsHi[cI]) + "%";

    label_Trk_p->DrawLatex(0.15, 0.94, "#bf{CMS Preliminary}");
    label_Trk_p->DrawLatex(0.25, 0.82, ("#bf{" + centStr2 + "}").c_str());
    label_Trk_p->DrawLatex(0.35, 0.82, "With Trk");
    leg_Trk_p->Draw("SAME");

    pads[1] = new TPad("pad1", "", 0.0, 0.0, 1.0, yPadBottomFrac);
    canv_p->cd();
    pads[1]->SetTopMargin(0.001);
    pads[1]->SetRightMargin(0.01);
    pads[1]->SetBottomMargin(pads[0]->GetLeftMargin()/yPadBottomFrac);
    pads[1]->SetLeftMargin(pads[0]->GetLeftMargin());
    pads[1]->Draw("SAME");
    pads[1]->cd();

    v2Raw_Trk_h[cI]->Divide(jamesHist_p[cI]);
    v2RawCorr_Trk_h[cI]->Divide(jamesHist_p[cI]);
    v2ObsCorr_Trk_h[cI]->Divide(jamesHist_p[cI]);

    v2Raw_Trk_h[cI]->SetMinimum(0.55);
    v2Raw_Trk_h[cI]->SetMaximum(1.45);

    v2Raw_Trk_h[cI]->GetYaxis()->SetNdivisions(505);

    v2Raw_Trk_h[cI]->DrawCopy("HIST P");
    v2RawCorr_Trk_h[cI]->DrawCopy("SAME HIST P");
    v2ObsCorr_Trk_h[cI]->DrawCopy("SAME HIST P");

    gStyle->SetOptStat(0);

    const std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
    const std::string saveName = "plotRecreateV2James_" + centStr + "_" + dateStr + "_Trk.pdf";

    canv_p->SaveAs(saveName.c_str());

    delete pads[0];
    delete pads[1];
  }
  */

  TexSlides(new std::vector<std::vector<std::string>*>{canv_p->GetPointer()},"Slides.tex",1);
  delete canv_p;

  gSystem->cd("..");

  inFile_p->Close();
  delete inFile_p;

  jamesFile_p->Close();
  delete jamesFile_p;


  delete leg_PF_p;
  delete label_PF_p;
  //delete leg_Trk_p;
  //  delete label_Trk_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./plotRecreateV2V3.exe <inFileName> <inJamesFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += plotRecreateV2V3(argv[1], argv[2]);
  return retVal;
}
