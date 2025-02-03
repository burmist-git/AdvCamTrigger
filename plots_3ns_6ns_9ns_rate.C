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

void save_hist_to_csv(TString outfilename, TH1D *h1);

Int_t plots_3ns_6ns_9ns_rate(){

  TString fileN01;
  TString fileN02;
  TString fileN03;
  fileN01="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/nsb_1x_268MHz/trgA_test/hist_trgA_corsika_0000ID_3ns.root";
  fileN02="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/nsb_1x_268MHz/trgA_test/hist_trgA_corsika_0000ID_6ns.root";
  fileN03="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/nsb_1x_268MHz/trgA_test/hist_trgA_corsika_0000ID_9ns.root";
  //  
  TFile *f1 = new TFile(fileN01.Data());
  TFile *f2 = new TFile(fileN02.Data());
  TFile *f3 = new TFile(fileN03.Data());
  //
  //TGraph *gr01 = (TGraph*)f1->Get("_gr_wf_tmpl");
  //TGraph *gr02 = (TGraph*)f2->Get("_gr_wf_tmpl");
  //TGraph *gr03 = (TGraph*)f3->Get("_gr_wf_tmpl");
  //
  //TH1D *gr01 = (TH1D*)f1->Get("h1_rate_digital_sum_pe");
  //TH1D *gr02 = (TH1D*)f2->Get("h1_rate_digital_sum_pe");
  //TH1D *gr03 = (TH1D*)f3->Get("h1_rate_digital_sum_pe");
  //
  //
  //TH1D *gr01 = (TH1D*)f1->Get("h1_digital_sum");
  //TH1D *gr02 = (TH1D*)f2->Get("h1_digital_sum");
  //TH1D *gr03 = (TH1D*)f3->Get("h1_digital_sum");
  //
  TH1D *gr01 = (TH1D*)f1->Get("h1_fadc_val");
  TH1D *gr02 = (TH1D*)f2->Get("h1_fadc_val");
  TH1D *gr03 = (TH1D*)f3->Get("h1_fadc_val");
  //
  //save_hist_to_csv("rate_digital_sum_pe_3ns.csv",gr01);
  //save_hist_to_csv("rate_digital_sum_pe_6ns.csv",gr02);
  //save_hist_to_csv("rate_digital_sum_pe_9ns.csv",gr03);
  //  
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
  // 
  gr01->SetLineColor(kBlack);
  gr01->SetLineWidth(3.0);
  gr01->SetMarkerColor(kBlack);
  gr01->SetMarkerStyle(20);
  //
  gr02->SetLineColor(kRed+2);
  gr02->SetLineWidth(3.0);
  gr02->SetMarkerColor(kRed+2);
  gr02->SetMarkerStyle(20);
  //
  gr03->SetLineColor(kBlue+2);
  gr03->SetLineWidth(3.0);
  gr03->SetMarkerColor(kBlue+2);
  gr03->SetMarkerStyle(20);
  //
  //
  //
  //TMultiGraph *mg = new TMultiGraph();
  //mg->Add(gr01);
  //mg->Add(gr02);
  //mg->Add(gr03);
  //mg->Draw("apl");
  //
  //
  //
  gr03->SetMaximum(6.0e+5);
  gr03->Draw();
  gr02->Draw("sames");
  gr01->Draw("sames");  
  //
  //mg->GetXaxis()->SetTitle("Time, a.u.");
  //mg->GetYaxis()->SetTitle("FADC counts");
  //mg->GetXaxis()->SetTitle("Time, ns");
  //mg->GetYaxis()->SetTitle("Amplitude, a.u.");
  //  
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(gr01, "FWHM 3 ns", "apl");
  leg->AddEntry(gr02, "FWHM 6 ns", "apl");
  leg->AddEntry(gr03, "FWHM 9 ns", "apl");
  leg->Draw();  

  return 0;
}

void save_hist_to_csv(TString outfilename, TH1D *h1){
  ofstream outfile;
  outfile.open(outfilename.Data());
  //
  outfile<<"x,y"<<endl;
  for(Int_t i = 1;i<=h1->GetNbinsX();i++){
    Double_t val = h1->GetBinContent(i);
    Double_t valx = h1->GetBinCenter(i);
    if(valx>0.0)
      outfile<<valx<<","<<val<<endl;
  }
  outfile.close();
}
