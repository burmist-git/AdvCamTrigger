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

Int_t plots_SuperFast_eff(){
  //
  TString fileN01;
  fileN01 = "./hist_superfast_gamma_on_nsb_1x.root";
  //
  TFile *f1 = new TFile(fileN01.Data());
  //TH1D *h1_all = (TH1D*)f1->Get("_h1_E_evH_sim_more_then_100pe_all_soft");
  //TH1D *h1_more = (TH1D*)f1->Get("_h1_E_evH_sim_more_then_100pe_100peminperch_soft");
  //TH1D *h1_all = (TH1D*)f1->Get("_h1_E_evH_sim_more_then_100pe_all");
  //TH1D *h1_more = (TH1D*)f1->Get("_h1_E_evH_sim_more_then_100pe_100peminperch");
  TH1D *h1_all = (TH1D*)f1->Get("_h1_E_evH_sim_more_then_100pe_100peminperch_eff");
  TH1D *h1_more = (TH1D*)f1->Get("_h1_E_evH_sim_more_then_100pe_100peminperch_soft_eff");
  //
  TCanvas *c1 = new TCanvas("c1",fileN01.Data(),10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE);
  //
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  //
  h1_all->SetLineWidth(2.0);
  h1_all->SetLineColor(kBlack);
  h1_more->SetLineWidth(2.0);
  h1_more->SetLineColor(kRed+2);
  h1_all->SetMarkerColor(kBlack);
  h1_more->SetMarkerColor(kRed+2);
  //  
  h1_all->SetTitle("");
  h1_all->Draw();
  h1_more->Draw("sames");
  h1_all->SetMinimum(1.0e-5);
  h1_all->SetMaximum(1.0);
  
  //
  //h1_all->GetYaxis()->SetTitle("Number of on-axis gammas");
  h1_all->GetYaxis()->SetTitle("Ratio");
  h1_all->GetXaxis()->SetTitle("Energy, GeV");
  //
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  //leg->AddEntry(h1_all,  "All events  (>100 p.e. in total) (1/E^2.7)", "apl");
  //leg->AddEntry(h1_more, "Events with (>100 p.e./ channel) (1/E^2.7)", "apl");
  leg->AddEntry(h1_all,  "Ratio (1/E^2.7)", "apl");
  leg->AddEntry(h1_more, "Ratio (1/E^2.2)", "apl");
  leg->Draw();  

  //
  return 0;
}
