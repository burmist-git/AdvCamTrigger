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

TString getObjectName(TString namein, Int_t i);

Int_t plots_wf(){

  TString fileN;
  TString objectName;
  fileN = "./hist_corsika_run307_gamma_no_nsb_cut.root";

  TFile *f1 = new TFile(fileN.Data());

  TGraph *gr = (TGraph*)f1->Get("gr_ch_3336");
  TGraph *gr_sim = (TGraph*)f1->Get("gr_sim_ch_3336");  

  
  TCanvas *c1 = new TCanvas("c1",fileN.Data(),10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
  // 
  gr->SetLineColor(kBlack);
  gr->SetLineWidth(3.0);
  gr->SetMarkerColor(kBlack);
  gr->SetMarkerStyle(20);
  //
  gr_sim->SetLineColor(kRed);
  gr_sim->SetLineWidth(3.0);
  gr_sim->SetMarkerColor(kRed);
  gr_sim->SetMarkerStyle(20);
  //
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr);
  mg->Add(gr_sim);
  //
  mg->Draw("apl");
  //
  mg->GetXaxis()->SetTitle("Time, a.u.");
  mg->GetYaxis()->SetTitle("FADC counts");
  //  
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(gr, "simtelarr", "apl");
  leg->AddEntry(gr_sim, "sim", "apl");
  leg->Draw();  

  return 0;
}
