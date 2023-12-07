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

Int_t plots_corey_vs_corex(){
  //
  TString fileN01;
  TString fileN02;
  TString fileN03;
  TString fileN04;
  fileN01 = "./hist_fast_gamma_on_nsb_1x_42987507ev.root";
  //fileN02 = fileN01;
  //fileN03 = fileN01;
  //fileN04 = fileN01;
  fileN02 = "./hist_fast_gamma_diffuse_nsb_1x_138771820ev.root";
  fileN03 = "./hist_fast_electron_nsb_1x_48347162ev.root";
  fileN04 = "./hist_fast_proton_nsb_1x_112643300ev.root";
  //
  TFile *f01 = new TFile(fileN01.Data());
  TFile *f02 = new TFile(fileN02.Data());
  TFile *f03 = new TFile(fileN03.Data());
  TFile *f04 = new TFile(fileN04.Data());
  //
  TH2D *h2_g = (TH2D*)f01->Get("h2_ycore_vs_xcore");
  TH2D *h2_gr = (TH2D*)f01->Get("h2_ycore_vs_xcore_ring");
  //
  TH2D *h2_gd = (TH2D*)f02->Get("h2_ycore_vs_xcore");
  TH2D *h2_gdr = (TH2D*)f02->Get("h2_ycore_vs_xcore_ring");
  //
  TH2D *h2_e = (TH2D*)f03->Get("h2_ycore_vs_xcore");
  TH2D *h2_er = (TH2D*)f03->Get("h2_ycore_vs_xcore_ring");
  //
  TH2D *h2_p = (TH2D*)f04->Get("h2_ycore_vs_xcore");
  TH2D *h2_pr = (TH2D*)f04->Get("h2_ycore_vs_xcore_ring");
  //
  h2_gr->SetMarkerColor(kRed);
  h2_gdr->SetMarkerColor(kRed);
  h2_er->SetMarkerColor(kRed);
  h2_pr->SetMarkerColor(kRed);
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,1000);
  c1->Divide(2,2);
  c1->cd(1);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE); 
  //
  gPad->SetGridx();
  gPad->SetGridy();
  //
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.1);
  gPad->SetTopMargin(0.15);
  gPad->SetBottomMargin(0.1);
  //
  h2_g->SetTitle("");
  h2_g->Draw("ZCOLOR");
  h2_gr->SetMarkerColor(kRed);
  h2_gr->Draw("same");
  //
  h2_g->GetXaxis()->SetTitle("x, m");
  h2_g->GetYaxis()->SetTitle("y, m");
  //
  //

  c1->cd(2);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE); 
  //
  gPad->SetGridx();
  gPad->SetGridy();
  //
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.1);
  gPad->SetTopMargin(0.15);
  gPad->SetBottomMargin(0.1);
  //
  h2_gd->SetTitle("");
  h2_gd->Draw("ZCOLOR");
  h2_gdr->SetMarkerColor(kRed);
  h2_gdr->Draw("same");
  //
  h2_gd->GetXaxis()->SetTitle("x, m");
  h2_gd->GetYaxis()->SetTitle("y, m");
  //
  //
  c1->cd(3);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE); 
  //
  gPad->SetGridx();
  gPad->SetGridy();
  //
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.1);
  gPad->SetTopMargin(0.15);
  gPad->SetBottomMargin(0.1);
  //
  h2_e->SetTitle("");
  h2_e->Draw("ZCOLOR");
  h2_er->SetMarkerColor(kRed);
  h2_er->Draw("same");
  //
  h2_e->GetXaxis()->SetTitle("x, m");
  h2_e->GetYaxis()->SetTitle("y, m");
  //
  //
  c1->cd(4);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE); 
  //
  gPad->SetGridx();
  gPad->SetGridy();
  //
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.1);
  gPad->SetTopMargin(0.15);
  gPad->SetBottomMargin(0.1);
  //
  h2_p->SetTitle("");
  h2_p->Draw("ZCOLOR");
  h2_pr->SetMarkerColor(kRed);
  h2_pr->Draw("same");
  //
  h2_p->GetXaxis()->SetTitle("x, m");
  h2_p->GetYaxis()->SetTitle("y, m");
  //
  //
  TCanvas *c01 = new TCanvas("c01","c01",10,10,800,800);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE); 
  gPad->SetGridx();
  gPad->SetGridy();
  c01->SetRightMargin(0.15);
  c01->SetLeftMargin(0.15);
  c01->SetTopMargin(0.15);
  c01->SetBottomMargin(0.1);
  //
  h2_g->SetTitle("");
  h2_g->Draw("ZCOLOR");
  h2_gr->SetMarkerColor(kRed);
  h2_gr->Draw("same");
  h2_g->GetXaxis()->SetTitle("x, m");
  h2_g->GetYaxis()->SetTitle("y, m");
  c01->SaveAs("c01_h2_g.pdf");
  //
  //
  TCanvas *c02 = new TCanvas("c02","c02",10,10,800,800);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE); 
  gPad->SetGridx();
  gPad->SetGridy();
  c02->SetRightMargin(0.15);
  c02->SetLeftMargin(0.15);
  c02->SetTopMargin(0.15);
  c02->SetBottomMargin(0.1);
  //
  h2_gd->SetTitle("");
  h2_gd->Draw("ZCOLOR");
  h2_gdr->SetMarkerColor(kRed);
  h2_gdr->Draw("same");
  h2_gd->GetXaxis()->SetTitle("x, m");
  h2_gd->GetYaxis()->SetTitle("y, m");
  c02->SaveAs("c01_h2_gd.pdf");
  //
  //
  TCanvas *c03 = new TCanvas("c03","c03",10,10,800,800);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE); 
  gPad->SetGridx();
  gPad->SetGridy();
  c03->SetRightMargin(0.15);
  c03->SetLeftMargin(0.15);
  c03->SetTopMargin(0.15);
  c03->SetBottomMargin(0.1);
  //
  h2_e->SetTitle("");
  h2_e->Draw("ZCOLOR");
  h2_er->SetMarkerColor(kRed);
  h2_er->Draw("same");
  h2_e->GetXaxis()->SetTitle("x, m");
  h2_e->GetYaxis()->SetTitle("y, m");
  c03->SaveAs("c01_h2_e.pdf");
  //
  //
  TCanvas *c04 = new TCanvas("c04","c04",10,10,800,800);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE); 
  gPad->SetGridx();
  gPad->SetGridy();
  c04->SetRightMargin(0.15);
  c04->SetLeftMargin(0.15);
  c04->SetTopMargin(0.15);
  c04->SetBottomMargin(0.1);
  //
  h2_p->SetTitle("");
  h2_p->Draw("ZCOLOR");
  h2_pr->SetMarkerColor(kRed);
  h2_pr->Draw("same");
  h2_p->GetXaxis()->SetTitle("x, m");
  h2_p->GetYaxis()->SetTitle("y, m");
  c04->SaveAs("c01_h2_p.pdf");

  
  
  //
  //
  //
  //
  /*
  TCanvas *c2 = new TCanvas("c2","c2",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE); 
  //
  gPad->SetGridx();
  gPad->SetGridy();
  //
  c2->SetRightMargin(0.15);
  c2->SetLeftMargin(0.1);
  c2->SetTopMargin(0.15);
  c2->SetBottomMargin(0.1);
  //
  h2_gd->SetTitle("");
  h2_gd->Draw("ZCOLOR");
  h2_gdr->SetMarkerColor(kRed);
  h2_gdr->Draw("same");
  //
  h2_gd->GetXaxis()->SetTitle("x, m");
  h2_gd->GetYaxis()->SetTitle("y, m");
  //
  */
  return 0;
}
