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

void norm_hist(TH1D *h1, TH1D *h1_norm);

Int_t plots_v(){
  //
  TString fileN01;
  TString fileN02;
  //
  fileN01 = "./hist_corsika_run307.root";
  fileN02 = "../terzina_wfSim/wfSim_simtelarray_sipm_Rate_386MHz_el_noise.root";
  //fileN02 = "./hist_corsika_run307_no_nsb_cut.root";
  //
  TFile *f1 = new TFile(fileN01.Data());
  TFile *f2 = new TFile(fileN02.Data());
  //
  //TH1D *h1_01 = (TH1D*)f1->Get("h1_v");
  //TH1D *h1_02 = (TH1D*)f2->Get("h1_v");
  TH1D *h1_01 = (TH1D*)f1->Get("h1_d_v");
  TH1D *h1_02 = (TH1D*)f2->Get("h1_d_v");
  //
  TH1D *h1_01_norm = new TH1D();
  TH1D *h1_02_norm = new TH1D();
  //
  norm_hist(h1_01,h1_01_norm);
  norm_hist(h1_02,h1_02_norm);
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE);
  // 
  h1_01_norm->SetLineColor(kBlack);
  h1_01_norm->SetLineWidth(3.0);
  //
  h1_02_norm->SetLineColor(kRed);
  h1_02_norm->SetLineWidth(3.0);
  //
  h1_02_norm->SetTitle("");
  //
  h1_01_norm->Draw();
  h1_02_norm->Draw("sames");
  //
  h1_01_norm->GetXaxis()->SetTitle("FADC counts");
  //  
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(h1_01_norm, "386 MHz (simtelarray)", "apl");
  leg->AddEntry(h1_02_norm, "386 MHz (wfsim)", "apl");
  leg->Draw();
  //
  return 0;
}

void norm_hist(TH1D *h1, TH1D *h1_norm){
  h1_norm->SetBins(h1->GetNbinsX(),
		   h1->GetBinLowEdge(1),
		   h1->GetBinLowEdge(h1->GetNbinsX())+h1->GetBinWidth(h1->GetNbinsX()));
  for(Int_t i = 1;i<=h1->GetNbinsX();i++){
    h1_norm->SetBinContent(i,h1->GetBinContent(i)/h1->GetEntries());
  }
}
