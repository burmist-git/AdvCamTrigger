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

Int_t plots_proton_diff_DBSCAN(){

  TString fileN01;
  //
  fileN01 = "./hist_rate_proton_diff.root";
  //
  TFile *f01 = new TFile(fileN01.Data());
  //
  TH1D *h1_tot = (TH1D*)f01->Get("_h1_rate_tot_vs_th");
  TH1D *h1_protons = (TH1D*)f01->Get("_h1_rate_vs_th");
  TH1D *h1_NSB = (TH1D*)f01->Get("_h1_rate_NSB_vs_th");
  //
  TCanvas *c1 = new TCanvas("c1",fileN01.Data(),10,10,800,800);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE);
  //
  h1_tot->SetLineColor(kBlack);
  h1_tot->SetLineWidth(3);
  //
  h1_protons->SetLineColor(kBlue+2);
  h1_protons->SetLineWidth(3);
  //
  h1_NSB->SetLineColor(kRed+2);
  h1_NSB->SetLineWidth(3);
  //
  h1_tot->SetTitle("");
  h1_tot->Draw();
  h1_protons->Draw("same");
  h1_NSB->Draw("same");
  //
  h1_tot->SetMinimum(10);
  h1_tot->SetMaximum(20.0*1.0e+6);
  h1_tot->SetMinimum(1.0e+3);
  //
  h1_tot->GetXaxis()->SetTitle("Number of points in the cluster");
  h1_tot->GetYaxis()->SetTitle("Rate, Hz");
  //
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(h1_tot, "NSB(268 MHz) + protons + NSB(~50 KHz)", "apl");
  leg->AddEntry(h1_protons, "protons + NSB(~50 KHz)", "apl");
  leg->AddEntry(h1_NSB, "NSB(268 MHz)", "apl");
  leg->Draw();  
  //
  //TString file_out = fileN01;
  //file_out += ".pdf";
  //c1->SaveAs(file_out.Data());
  //
  return 0;
}
