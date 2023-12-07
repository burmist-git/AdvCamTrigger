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

Int_t plot_digital_sum_NSB(){

  TString fileN01;
  
  fileN01 = "./hist_trg_proton_nsb_1x_corsika_0000ID.root";

  TFile *f01 = new TFile(fileN01.Data());


  TH1D *h1_01 = (TH1D*)f01->Get("h1_digital_sum");
  TH1D *h1_02 = (TH1D*)f01->Get("h1_digital_sum_3ns");
  TH1D *h1_03 = (TH1D*)f01->Get("h1_digital_sum_5ns");
  
  //h1_1->Rebin(4);
  //h1_1->SetTitle("");

  //TCanvas *c1 = new TCanvas("c1",fileN.Data(),10,10,800,800);
  TCanvas *c1 = new TCanvas("c1",fileN01.Data(),10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
 
  h1_01->SetLineColor(kBlack);
  h1_01->SetLineWidth(3.0);
  h1_02->SetLineColor(kRed);
  h1_02->SetLineWidth(3.0);
  h1_03->SetLineColor(kBlue+2);
  h1_03->SetLineWidth(3.0);
  //
  h1_03->Draw();
  h1_02->Draw("sames");
  h1_01->Draw("sames");
  //
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  //
  leg->AddEntry(h1_01, "1ns", "apl");
  leg->AddEntry(h1_02, "3ns", "apl");
  leg->AddEntry(h1_03, "5ns", "apl");
  leg->Draw();
  return 0;
}
