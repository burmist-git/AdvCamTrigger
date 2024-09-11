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

void coincident_rate( TGraph *gr, Double_t dt, Double_t Rmin, Double_t Rmax, Int_t n);
void coincident_ratio_rate( TGraph *gr, Double_t dt, Double_t Rmin, Double_t Rmax, Int_t n);
void lin_rate( TGraph *gr, Double_t Rmin, Double_t Rmax, Int_t n);

Int_t plots_coincident_rate(){
  //
  TRandom3 *rnd = new TRandom3(1231231);
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE); 
  //
  //
  //
  TGraph *gr_200ns = new TGraph();
  gr_200ns->SetNameTitle("gr_200ns","gr_200ns");
  TGraph *gr_100ns = new TGraph();
  gr_100ns->SetNameTitle("gr_100ns","gr_100ns");
  TGraph *gr_50ns = new TGraph();
  gr_50ns->SetNameTitle("gr_50ns","gr_50ns");
  TGraph *gr_lin = new TGraph();
  gr_lin->SetNameTitle("gr_lin","gr_lin");
  //
  //
  TGraph *gr_200ns_R = new TGraph();
  gr_200ns_R->SetNameTitle("gr_200ns_R","gr_200ns_R");
  TGraph *gr_100ns_R = new TGraph();
  gr_100ns_R->SetNameTitle("gr_100ns_R","gr_100ns_R");
  TGraph *gr_50ns_R = new TGraph();
  gr_50ns_R->SetNameTitle("gr_50ns_R","gr_50ns_R");  
  //
  //
  coincident_rate( gr_200ns, 200*1.0e-9, 80000, 50.0*1.0e+6, 1000);
  coincident_rate( gr_100ns, 100*1.0e-9, 80000, 50.0*1.0e+6, 1000);
  coincident_rate( gr_50ns,   50*1.0e-9, 80000, 50.0*1.0e+6, 1000);
  //
  coincident_ratio_rate( gr_200ns_R, 200*1.0e-9, 80000, 2.0*1.0e+6, 1000);
  coincident_ratio_rate( gr_100ns_R, 100*1.0e-9, 80000, 2.0*1.0e+6, 1000);
  coincident_ratio_rate( gr_50ns_R,   50*1.0e-9, 80000, 2.0*1.0e+6, 1000);
  //
  //
  lin_rate( gr_lin, 80000, 50.0*1.0e+6, 1000);
  //
  gr_200ns->SetLineColor(kMagenta+2);
  gr_200ns->SetLineWidth(3.0);
  //
  gr_100ns->SetLineColor(kBlue+2);
  gr_100ns->SetLineWidth(3.0);
  //
  gr_50ns->SetLineColor(kRed+2);
  gr_50ns->SetLineWidth(3.0);
  //
  gr_lin->SetLineColor(kBlack);
  gr_lin->SetLineWidth(3.0);
  //
  gr_200ns_R->SetLineColor(kMagenta+2);
  gr_200ns_R->SetLineWidth(3.0);
  //
  gr_100ns_R->SetLineColor(kBlue+2);
  gr_100ns_R->SetLineWidth(3.0);
  //
  gr_50ns_R->SetLineColor(kRed+2);
  gr_50ns_R->SetLineWidth(3.0);
  //
  //gr->Draw("APL");
  //
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr_200ns);
  mg->Add(gr_100ns);
  mg->Add(gr_50ns);
  mg->Add(gr_lin);
  //
  mg->SetMinimum(1.1e+9);
  mg->SetMinimum(500.0);
  //
  mg->Draw("al");
  //
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  //
  mg->GetXaxis()->SetTitle("Rate, Hz");
  mg->GetYaxis()->SetTitle("Coincidence rate, Hz");
  //
  //
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(gr_200ns, "#tau = 200 ns", "al");
  leg->AddEntry(gr_100ns, "#tau = 100 ns", "al");
  leg->AddEntry(gr_50ns,  "#tau = 50 ns", "al");
  leg->AddEntry(gr_lin,   "Coinc. Rate = Rate", "al");
  leg->Draw();
  //
  //
  //
  TCanvas *c2 = new TCanvas("c2","c2",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE); 
  //
  TMultiGraph *mg02 = new TMultiGraph();
  mg02->Add(gr_200ns_R);
  mg02->Add(gr_100ns_R);
  mg02->Add(gr_50ns_R);
  //
  mg02->SetMinimum(1.100);
  mg02->SetMinimum(0.005);
  //
  mg02->Draw("al");
  //
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  //
  mg02->GetXaxis()->SetTitle("Rate, Hz");
  mg02->GetYaxis()->SetTitle("Coincidence rate/Rate");  
  //
  //
  TLegend *leg02 = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg02->AddEntry(gr_200ns, "#tau = 200 ns", "al");
  leg02->AddEntry(gr_100ns, "#tau = 100 ns", "al");
  leg02->AddEntry(gr_50ns,  "#tau = 50 ns", "al");
  leg02->Draw();
  //
  //
  Double_t p_rate_hz = 30000.0; //Hz
  Double_t nSec = 100.0;

  int nn = nSec*p_rate_hz;
  Double_t *evt_time = new Double_t[nn];

  Double_t n_ns_Sec = nSec/1.0e-9;
  cout<<"nSec "<<nSec<<endl;
  for (int k = 0; k < nn; ++k){
    evt_time[k] = rnd->Uniform(0.0,nSec);
  }
  long long size_nn = (long long)nn;
  long long ind[nn];
  TMath::Sort(size_nn,evt_time,ind,kFALSE);
  //
  TH1D *h1_dt = new TH1D("h1_dt","h1_dt", (Int_t)500000.0/4000.0, 0.0, 500000.0);
  TH1D *h1_t = new TH1D("h1_t","h1_t",(Int_t)((n_ns_Sec+2.0/1.0e-9)/1.0e+9),-1.0/1.0e-9, n_ns_Sec+1.0/1.0e-9);
  //
  TCanvas *c3 = new TCanvas("c3","c3",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE); 
  //
  Double_t evt_time_ns;
  Double_t evt_time_ns_t0;
  Double_t evt_time_ns_t1;
  //
  for (int k = 0; k < nn; ++k){
    //
    evt_time_ns = evt_time[k]/1.0e-9;
    h1_t->Fill(evt_time_ns);
    //
    if(k<(nn-1)){
      evt_time_ns_t0 = evt_time[ind[k]]/1.0e-9;
      evt_time_ns_t1 = evt_time[ind[k+1]]/1.0e-9;
      h1_dt->Fill(evt_time_ns_t1-evt_time_ns_t0);
    }
  }
  //
  h1_t->Draw();
  


  TCanvas *c4 = new TCanvas("c4","c4",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  h1_dt->Draw();
  
  /*
  //rnd
  float a[8]={2, 1, 5, 4,7,15, 12, 9};
  //for (int k = 0; k < size; ++k){
  //std::cout <<a[k]<< " INDEX   "<< ind[k] << std::endl;
  //}
  */  
  //
  //
  return 0;
}

void coincident_rate( TGraph *gr, Double_t dt, Double_t Rmin, Double_t Rmax, Int_t n){
  Double_t R;
  Double_t rc;
  Double_t R1;
  Double_t R2;
  for( Int_t i = 0; i < n; i++){
    //
    R = (Rmax - Rmin)/(n-1)*i + Rmin;
    R1 = R;
    R2 = R;
    //
    rc = 2*R1*R2*dt;
    gr->SetPoint(i,R,rc);
    //if(rc<R)
    //gr->SetPoint(i,R,rc/R);
    //else
    //gr->SetPoint(i,R,1);
  }
}

void coincident_ratio_rate( TGraph *gr, Double_t dt, Double_t Rmin, Double_t Rmax, Int_t n){
  Double_t R;
  Double_t rc;
  Double_t R1;
  Double_t R2;
  for( Int_t i = 0; i < n; i++){
    //
    R = (Rmax - Rmin)/(n-1)*i + Rmin;
    R1 = R;
    R2 = R;
    //
    rc = 2*R1*R2*dt;
    if(rc<R)
      gr->SetPoint(i,R,rc/R);
    else
      gr->SetPoint(i,R,1);
  }
}

void lin_rate( TGraph *gr, Double_t Rmin, Double_t Rmax, Int_t n){
  Double_t R;
  for( Int_t i = 0; i < n; i++){
    R = (Rmax - Rmin)/(n-1)*i + Rmin;
    gr->SetPoint(i,R,R);
  }
}

