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

Int_t plots_k_dist(){

  TString fileN01;
  TString fileN02;
  //
  fileN01 = "./hist_trg_counter_1ev.root";
  fileN02 = "./hist_50pe_trg_counter_1ev.root"; 

  //
  TFile *f01 = new TFile(fileN01.Data());
  TFile *f02 = new TFile(fileN02.Data());
  //
  TGraph *gr01 = (TGraph*)f01->Get("gr_k_dist_graph");
  TGraph *gr02 = (TGraph*)f02->Get("gr_k_dist_graph");
  //
  TCanvas *c1 = new TCanvas("c1",fileN01.Data(),10,10,800,800);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
  //
  gr01->SetMarkerColor(kBlack);
  gr01->SetMarkerStyle(24);
  //
  gr01->SetMarkerColor(kBlack);
  gr01->SetMarkerSize(1);
  gr01->SetLineWidth(2);
  gr01->SetMarkerStyle(20);
  //
  gr02->SetMarkerColor(kRed+1);
  gr02->SetMarkerSize(1);
  gr02->SetLineWidth(2);
  gr02->SetMarkerStyle(20);
  //
  TMultiGraph *mg = new TMultiGraph();  
  mg->Add(gr01);
  mg->Add(gr02);
  mg->Draw("ap");
  //
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(gr01, "14-dist", "ap");
  leg->AddEntry(gr02, "4-dist", "ap");
  leg->Draw();  
  //
  //TString file_out = fileN01;
  //file_out += ".pdf";
  //c1->SaveAs(file_out.Data());
  //
  return 0;
}
