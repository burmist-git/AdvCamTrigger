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

Int_t plots_stareo_dangle(){

  TString fileN01;
  //
  fileN01 = "./hist_proton_st.root";
  //
  TFile *f01 = new TFile(fileN01.Data());
  //
  TH1D *h1_01 = (TH1D*)f01->Get("h1_dalpha_LST1_LST2");
  TH1D *h1_02 = (TH1D*)f01->Get("h1_dalpha_LST1_LST3");
  TH1D *h1_03 = (TH1D*)f01->Get("h1_dalpha_LST1_LST4");
  //
  h1_01->Rebin(8);
  h1_02->Rebin(8);
  h1_03->Rebin(8);
  //
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE); 
  //
  h1_01->SetLineColor(kBlack);
  h1_01->SetLineWidth(3.0);
  h1_02->SetLineColor(kRed+2);
  h1_02->SetLineWidth(3.0);
  h1_03->SetLineColor(kBlue+2);
  h1_03->SetLineWidth(3.0);
  //
  h1_01->SetTitle("");
  //
  h1_01->SetMaximum(30000);
  h1_01->Draw("HISTE");
  h1_02->Draw("HISTEsames");
  h1_03->Draw("HISTEsames");
  //
  h1_01->GetXaxis()->SetTitle("Angle between cent. of gravity, deg");
  //
  gPad->SetGridx();
  gPad->SetGridy();
  //
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(h1_01, "LST1 - LST2", "apl");
  leg->AddEntry(h1_02, "LST1 - LST3", "apl");
  leg->AddEntry(h1_03, "LST1 - LST4", "apl");
  leg->Draw();
  //
  return 0;
}
