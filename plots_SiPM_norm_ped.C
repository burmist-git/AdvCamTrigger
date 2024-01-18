//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include <time.h>

using namespace std;

Int_t plots_SiPM_norm_ped(){

  TString fileN;
  fileN = "./hist_SiPM.root";

  TFile *f1 = new TFile(fileN.Data());

  TGraph *gr = (TGraph*)f1->Get("gr_aml_gain_pedestal");

  TCanvas *c1 = new TCanvas("c1",fileN.Data(),10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
 
  gr->SetLineColor(kBlack);
  gr->SetLineWidth(3.0);

  gr->SetMaximum(0.0018);
  gr->SetMinimum(0.0);
  gr->SetTitle("");
  
  gr->Draw("APL");

  return 0;
}
