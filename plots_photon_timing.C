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

Int_t plots_photon_timing(){
  //
  TString fileN01;
  //fileN01 = "./hist_superfast_proton_nsb_1x_0pecut.root";
  //fileN01 = "./hist_superfast_proton_nsb_1x_50pecut.root";
  //fileN01 = "./hist_superfast_gamma_on_nsb_1x_50pecut.root";
  fileN01 = "./hist_superfast_gamma_diffuse_nsb_1x_50pecut.root";
  TFile *f01 = new TFile(fileN01.Data());
  //
  TH1D *h1_01 = (TH1D*)f01->Get("h1_ev_time");
  TH1D *h1_02 = (TH1D*)f01->Get("h1_pe_time");
  TH1D *h1_03 = (TH1D*)f01->Get("h1_pe_time_shift");
  //  
  TCanvas *c1 = new TCanvas("c1","c1",10,10,1200,1000);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
  //
  h1_01->SetLineColor(kBlack);
  h1_01->SetLineWidth(3.0);
  h1_02->SetLineColor(kBlue+2);
  h1_02->SetLineWidth(3.0);
  h1_03->SetLineColor(kRed+2);
  h1_03->SetLineWidth(3.0);
  //  
  //gr->SetMaximum(0.0018);
  //gr->SetMinimum(0.0);  
  //gr->Draw("APL");
  //
  h1_03->SetTitle("");
  h1_03->Draw();
  h1_01->Draw("sames");
  h1_02->Draw("sames");  
  //  
  return 0;
}
