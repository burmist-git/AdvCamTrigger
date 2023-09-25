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

Int_t plot_digital_sum(){

  TString fileN01;
  TString fileN02;
  TString fileN03;
  TString fileN04;
  TString fileN05;
  
  TString objectName;
  fileN01 = "./hist_gamma_on_nsb_1x_NSB.root";
  fileN02 = "./hist_gamma_on_nsb_1x_20p.e.root";
  fileN03 = "./hist_gamma_on_nsb_1x_40p.e.root";
  fileN04 = "./hist_gamma_on_nsb_1x_60p.e.root";
  fileN05 = "./hist_gamma_on_nsb_1x_100p.e.root";

  TFile *f01 = new TFile(fileN01.Data());
  TFile *f02 = new TFile(fileN02.Data());
  TFile *f03 = new TFile(fileN03.Data());  
  TFile *f04 = new TFile(fileN04.Data());
  TFile *f05 = new TFile(fileN05.Data());

  TH1D *h1_01 = (TH1D*)f01->Get("h1_digital_sum");
  TH1D *h1_02 = (TH1D*)f02->Get("h1_digital_sum");
  TH1D *h1_03 = (TH1D*)f03->Get("h1_digital_sum");
  TH1D *h1_04 = (TH1D*)f04->Get("h1_digital_sum");
  TH1D *h1_05 = (TH1D*)f05->Get("h1_digital_sum");
  
  //h1_1->Rebin(4);
  //h1_1->SetTitle("");

  //TCanvas *c1 = new TCanvas("c1",fileN.Data(),10,10,800,800);
  TCanvas *c1 = new TCanvas("c1",fileN01.Data(),10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
 
  h1_01->SetLineColor(kBlack);
  h1_01->SetLineWidth(3.0);
  h1_02->SetLineColor(kRed);
  h1_02->SetLineWidth(3.0);
  h1_03->SetLineColor(kBlue);
  h1_03->SetLineWidth(3.0);
  h1_04->SetLineColor(kGreen+2);
  h1_04->SetLineWidth(3.0);
  h1_05->SetLineColor(kMagenta+2);
  h1_05->SetLineWidth(3.0);
  //
  h1_01->Draw();
  h1_02->Draw("same");
  h1_03->Draw("same");
  h1_04->Draw("same");
  h1_05->Draw("same");
  //
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  //
  leg->AddEntry(h1_01, "NSB", "apl");
  leg->AddEntry(h1_02, "20p.e.", "apl");
  leg->AddEntry(h1_03, "40p.e.", "apl");
  leg->AddEntry(h1_04, "60p.e.", "apl");
  leg->AddEntry(h1_05, "100p.e.", "apl");
  //  
  leg->Draw();
  //h1_1->GetXaxis()->SetTitle("Value, Unit");
  //h1_1->GetXaxis()->SetRangeUser(-0.13,-0.12);

  /*
  for(i = 0;i<nChannels;i++){
    h1_Arr[i]->SetLineColor(colorArr[i]);
    h1_Arr[i]->SetLineWidth(3.0);
    if(i == 0){
      h1_Arr[i]->Draw();
      h1_Arr[i]->GetXaxis()->SetTitle("Value, Unit");
      //h1_1->GetXaxis()->SetRangeUser(-0.13,-0.12);
    }
    else
      h1_Arr[i]->Draw("same");
  }
  
  //h1_1->Fit("gaus");
  TString legInfo;
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  for(i = 0;i<nChannels;i++){
    legInfo = "ch ";legInfo += i;
    leg->AddEntry(h1_Arr[i], legInfo.Data(), "l");
  }
  leg->Draw();
  */

  /*
  TMultiGraph *mg = new TMultiGraph();
  
  for(i = 0;i<nChannels;i++){
    gr_Arr[i]->SetLineColor(colorArr[i]);
    gr_Arr[i]->SetLineWidth(3.0);
    gr_Arr[i]->SetMarkerColor(colorArr[i]);
    gr_Arr[i]->SetMarkerStyle(markerArr[i]);
    mg->Add(gr_Arr[i]);
  }

  mg->Draw("apl");
  
  mg->GetXaxis()->SetTitle("ValueX, Unit");
  mg->GetYaxis()->SetTitle("ValueY, Unit");
  
  TString legInfo;

  for(i = 0;i<nChannels;i++){
    legInfo = "ch ";legInfo += i;
    leg->AddEntry(gr_Arr[i], legInfo.Data(), "apl");
  }
  leg->Draw();
  
  */
  return 0;
}
