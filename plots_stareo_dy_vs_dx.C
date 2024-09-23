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

Int_t plots_stareo_dy_vs_dx(){

  TString fileN01;
  //
  fileN01 = "./hist_proton_st.root";
  //
  TFile *f01 = new TFile(fileN01.Data());
  //
  TH2D *h2_01 = (TH2D*)f01->Get("h2_y_vs_x_LST3_LST4");
  //h1_01->Rebin(4);
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,1000);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE); 
  //
  h2_01->SetTitle("");
  //
  //h1_01->SetMaximum(1.13);
  h2_01->Draw("ZCOLOR");
  h2_01->GetXaxis()->SetTitle("Cent. of gravity x_{LST3} - Cent. of gravity x_{LST4}");
  h2_01->GetYaxis()->SetTitle("Cent. of gravity y_{LST3} - Cent. of gravity y_{LST4}");
  //
  gPad->SetGridx();
  gPad->SetGridy();
  //  
  return 0;
}
