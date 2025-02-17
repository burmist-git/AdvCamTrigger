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

void copyHistogram( TH1D *h1, TH1D *h1_copy, TString h1_name_title, bool ifBinsOnly = false, Double_t norm = 1.0);
void get_cumulative(TH1D *h1, TH1D *h1_cumulative, Double_t nevsim = 1.0);

Int_t plots_SuperFast_eff(){
  //
  TString fileN01;
  fileN01 = "./hist_superfast_gamma_on_nsb_1x.root";
  //
  TFile *f1 = new TFile(fileN01.Data());
  //TH1D *h1_all = (TH1D*)f1->Get("_h1_E_evH_sim_more_then_100pe_all_soft");
  //TH1D *h1_more = (TH1D*)f1->Get("_h1_E_evH_sim_more_then_100pe_100peminperch_soft");
  //TH1D *h1_all = (TH1D*)f1->Get("_h1_E_evH_sim_more_then_100pe_all_soft");
  TH1D *h1_more_hard = (TH1D*)f1->Get("_h1_E_evH_sim_more_then_100pe_100peminperch_soft");
  TH1D *h1_energy = (TH1D*)f1->Get("h1_energy");
  TH1D *h1_all = (TH1D*)f1->Get("_h1_E_evH_sim_more_then_100pe_all");
  TH1D *h1_more = (TH1D*)f1->Get("_h1_E_evH_sim_more_then_100pe_100peminperch");
  //TH1D *h1_all = (TH1D*)f1->Get("_h1_E_evH_sim_more_then_100pe_100peminperch_eff");
  //TH1D *h1_more = (TH1D*)f1->Get("_h1_E_evH_sim_more_then_100pe_100peminperch_soft_eff");
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
  h1_more->SetTitle("");
  h1_all->Draw();
  h1_more->Draw("sames");
  //h1_more->Draw();
  //h1_all->SetMinimum(1.0e+2);
  //h1_all->SetMaximum(1.0e+7);
  
  //
  h1_all->GetYaxis()->SetTitle("Number of on-axis gammas");
  //h1_all->GetYaxis()->SetTitle("Ratio");
  h1_all->GetXaxis()->SetTitle("Energy, GeV");
  //
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(h1_all,  "All events  (>100 p.e. in total) (1/E^2.2)", "apl");
  leg->AddEntry(h1_more, "Events with (>100 p.e./ channel in 10 channels) (1/E^2.2)", "apl");
  //leg->AddEntry(h1_all,  "Ratio (1/E^2.7)", "apl");
  //leg->AddEntry(h1_more, "Ratio (1/E^2.2)", "apl");
  //leg->Draw();  


  //
  //

  TCanvas *c2 = new TCanvas("c2",fileN01.Data(),10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE);
  //gStyle->SetOptStat(kTRUE);
  //
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  //
  //
  TH1D *h1_more_copy = new TH1D();
  get_cumulative(h1_more, h1_more_copy, h1_energy->Integral());
  TH1D *h1_more_hard_copy = new TH1D();
  get_cumulative(h1_more_hard, h1_more_hard_copy, h1_energy->Integral());
  //
  h1_more_copy->Draw();
  h1_more_hard_copy->Draw("same");
  //
  h1_more_copy->SetLineWidth(2.0);
  h1_more_copy->SetLineColor(kBlack);
  h1_more_hard_copy->SetLineWidth(2.0);
  h1_more_hard_copy->SetLineColor(kRed+2);
  //
  TLegend *leg02 = new TLegend( 0.6, 0.6, 0.9, 0.9,"","brNDC");
  leg02->AddEntry(h1_more_copy,      "1/E^2.7 (soft)", "apl");
  leg02->AddEntry(h1_more_hard_copy, "1/E^2.2 (hard)", "apl");
  leg02->Draw();
  //
  //
  return 0;
}

//void get_cumulative( TH1D *h1, TH1D *h1_cumulative);
void copyHistogram( TH1D *h1, TH1D *h1_copy, TString h1_name_title, bool ifBinsOnly, Double_t norm){
  h1_copy->SetNameTitle( h1_name_title.Data(), h1_name_title.Data());
  //
  Int_t n_bins = h1->GetNbinsX()+1;
  Double_t *bins_low_edge = new Double_t[n_bins];
  for(int i = 0;i<h1->GetNbinsX();i++)
    bins_low_edge[i] = h1->GetBinLowEdge(i+1);
  bins_low_edge[h1->GetNbinsX()] = h1->GetBinLowEdge(h1->GetNbinsX()) + h1->GetBinWidth(h1->GetNbinsX());
  //
  h1_copy->SetBins(h1->GetNbinsX(),bins_low_edge);
  if(!ifBinsOnly && norm>0.0)
    for(int i = 1;i<=h1->GetNbinsX();i++)
      h1_copy->SetBinContent(i,h1->GetBinContent(i)/norm);
}

void get_cumulative(TH1D *h1, TH1D *h1_cumulative, Double_t nevsim){
  TString nametitle=h1->GetName();
  nametitle += "_copy";
  copyHistogram( h1, h1_cumulative, nametitle, true, 1.0);
  for(int i = 1;i<h1->GetNbinsX();i++)
    h1_cumulative->SetBinContent(i,h1->Integral(i,h1->GetNbinsX())/nevsim);
}
