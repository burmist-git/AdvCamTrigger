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
  fileN = "./hist_corsika_run307.root";

  TFile *f1 = new TFile(fileN.Data());

  //
  //TGraph *gr_01 = (TGraph*)f1->Get("gr0pe_980");
  //TGraph *gr_02 = (TGraph*)f1->Get("gr1pe_980");
  //TGraph *gr_03 = (TGraph*)f1->Get("gr2pe_980");
  //TGraph *gr_04 = (TGraph*)f1->Get("gr3pe_980");
  //

  //TGraph *gr_01 = (TGraph*)f1->Get("gr1pe_980");
  //TGraph *gr_02 = (TGraph*)f1->Get("gr_0pe_av");
  TGraph *gr_03 = (TGraph*)f1->Get("gr_1pe_av");
  TGraph *gr_04 = (TGraph*)f1->Get("gr_2pe_av");
  TGraph *gr_05 = (TGraph*)f1->Get("gr_3pe_av");  

  TGraph *gr_01 = (TGraph*)f1->Get("gr_1pe_av");
  TGraph *gr_02 = (TGraph*)f1->Get("gr_template");

  
    
  //h1_1->Rebin(4);
  //h1_1->SetTitle("");

  //TCanvas *c1 = new TCanvas("c1",fileN.Data(),10,10,800,800);
  TCanvas *c1 = new TCanvas("c1",fileN.Data(),10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
  // 
  gr_01->SetLineColor(kBlack);
  gr_01->SetLineWidth(3.0);
  gr_01->SetMarkerColor(kBlack);
  gr_01->SetMarkerStyle(20);
  //
  gr_02->SetLineColor(kRed);
  gr_02->SetLineWidth(3.0);
  gr_02->SetMarkerColor(kRed);
  gr_02->SetMarkerStyle(20);
  //
  gr_03->SetLineColor(kBlue+2);
  gr_03->SetLineWidth(3.0);
  gr_03->SetMarkerColor(kBlue+2);
  gr_03->SetMarkerStyle(20);
  //
  gr_04->SetLineColor(kMagenta+2);
  gr_04->SetLineWidth(3.0);
  gr_04->SetMarkerColor(kMagenta+2);
  gr_04->SetMarkerStyle(20);
  //
  gr_05->SetLineColor(kGreen+2);
  gr_05->SetLineWidth(3.0);
  gr_05->SetMarkerColor(kGreen+2);
  gr_05->SetMarkerStyle(20);
  //
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr_01);
  mg->Add(gr_02);
  //mg->Add(gr_03);
  //mg->Add(gr_04);
  //mg->Add(gr_05);
  //
  mg->Draw("apl");
  //
  mg->GetXaxis()->SetTitle("Time, ns");
  mg->GetYaxis()->SetTitle("FADC counts");
  //  
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(gr_01, "1 p.e. (av.)", "apl");
  leg->AddEntry(gr_02, "1 p.e. (Template)", "apl");
  //leg->AddEntry(gr_01, "1 p.e.", "apl");
  //leg->AddEntry(gr_02, "0 p.e. (av.)", "apl");
  //leg->AddEntry(gr_03, "1 p.e. (av.)", "apl");
  //leg->AddEntry(gr_04, "2 p.e. (av.)", "apl");
  //leg->AddEntry(gr_05, "3 p.e. (av.)", "apl");
  leg->Draw();
  

  return 0;
}
