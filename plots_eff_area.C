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

Int_t plots_eff_area(){

  TString fileN01;
  TString fileN02;
  //
  fileN01 = "./hist_rate_gamma_DBscan.root";
  fileN02 = "./hist_rate_gamma_superflower.root";
  //fileN01 = "../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_on_nsb_1x/trgA_test_DBscan/hist_trgA_corsika_0binE_0binTheta_0binDist_0001ID.root";
  //fileN02 = "../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_on_nsb_1x/trgA_test_superflower/hist_trgA_corsika_0binE_0binTheta_0binDist_0001ID.root";
  //
  TFile *f01 = new TFile(fileN01.Data());
  TFile *f02 = new TFile(fileN02.Data());
  //
  TH1D *h1_eff_area_01 = (TH1D*)f01->Get("_h1_energy_eff_r");
  TH1D *h1_eff_area_02 = (TH1D*)f02->Get("_h1_energy_eff_r");
  //
  TCanvas *c1 = new TCanvas("c1",fileN01.Data(),10,10,800,800);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE);
  //
  h1_eff_area_01->SetLineColor(kBlack);
  h1_eff_area_01->SetLineWidth(3);
  //
  h1_eff_area_02->SetLineColor(kRed);
  h1_eff_area_02->SetLineWidth(3);
  //
  h1_eff_area_01->SetTitle("");
  h1_eff_area_01->Draw();
  h1_eff_area_02->Draw("same");
  h1_eff_area_01->SetMinimum(1.0e+4);
  h1_eff_area_01->SetMinimum(1.0e+7);

  /*
  //
  //h1_protons->SetTitle("");
  //h1_protons->Draw();
  h1_protons02->Draw("same");
  //
  //h1_protons->Draw("same");
  h1_NSB->Draw("same");
  //
  //h1_tot->SetMinimum(1.0e+3);
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
  */
  return 0;
}
