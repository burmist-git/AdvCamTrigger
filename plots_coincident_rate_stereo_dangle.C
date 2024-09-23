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

void coincident_rate( TGraph *gr, Double_t dt, Double_t dalpha, Double_t Rmin, Double_t Rmax, Int_t n);
void coincident_three_rate( TGraph *gr, Double_t dt, Double_t dalpha, Double_t Rmin, Double_t Rmax, Int_t n);

Int_t plots_coincident_rate_stereo_dangle(){
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
  TGraph *gr_80ns_1deg = new TGraph();
  gr_80ns_1deg->SetNameTitle("gr_80ns_1deg","gr_80ns_1deg");
  TGraph *gr_70ns_07deg = new TGraph();
  gr_70ns_07deg->SetNameTitle("gr_70ns_07deg","gr_70ns_07deg");
  TGraph *gr_80ns_1deg_three = new TGraph();
  gr_80ns_1deg_three->SetNameTitle("gr_80ns_1deg_three","gr_80ns_1deg_three");
  TGraph *gr_70ns_07deg_three = new TGraph();
  gr_70ns_07deg_three->SetNameTitle("gr_70ns_07deg_three","gr_70ns_07deg_three");
  //
  coincident_rate( gr_80ns_1deg,  80*1.0e-9, 1.0, 100000, 3.0*1.0e+6, 1000);
  coincident_rate( gr_70ns_07deg,  70*1.0e-9, 0.7, 100000, 3.0*1.0e+6, 1000);
  coincident_three_rate(gr_80ns_1deg_three, 80*1.0e-9, 1.0, 100000, 10.0*1.0e+6, 1000);
  coincident_three_rate(gr_70ns_07deg_three, 70*1.0e-9, 0.7, 100000, 10.0*1.0e+6, 1000);
  //
  gr_80ns_1deg->SetLineColor(kBlue);
  gr_80ns_1deg->SetLineWidth(3.0);
  //
  gr_70ns_07deg->SetLineColor(kRed);
  gr_70ns_07deg->SetLineWidth(3.0);
  //
  gr_80ns_1deg_three->SetLineColor(kBlue+2);
  gr_80ns_1deg_three->SetLineWidth(3.0);
  //
  gr_70ns_07deg_three->SetLineColor(kRed+2);
  gr_70ns_07deg_three->SetLineWidth(3.0);
  //
  //gr->Draw("APL");
  //
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr_80ns_1deg);
  mg->Add(gr_70ns_07deg);
  mg->Add(gr_80ns_1deg_three);
  mg->Add(gr_70ns_07deg_three);
  //
  //mg->SetMinimum(1.1e+9);
  mg->SetMaximum(1.0e+8);
  //
  mg->Draw("al");
  //
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  //
  mg->GetXaxis()->SetTitle("L0 rate per LST, Hz");
  mg->GetYaxis()->SetTitle("Total fake rate of the 4 LSTs, Hz");
  //
  //
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(gr_80ns_1deg,  "#tau = 80 ns and 1.0 deg topological trigger", "al");
  leg->AddEntry(gr_70ns_07deg, "#tau = 70 ns and 0.7 deg topological trigger", "al");
  leg->AddEntry(gr_80ns_1deg_three,  "#tau = 80 ns and 1.0 deg topological trigger (triple coinc.)", "al");
  leg->AddEntry(gr_70ns_07deg_three, "#tau = 70 ns and 0.7 deg topological trigger (triple coinc.)", "al");
  leg->Draw();
  //
  //
  return 0;
}

void coincident_rate( TGraph *gr, Double_t dt, Double_t dalpha, Double_t Rmin, Double_t Rmax, Int_t n){
  Double_t R;
  Double_t rc;
  Double_t R1;
  Double_t R2;
  Double_t dalpha_reduction = 4.0/dalpha/dalpha;  
  for( Int_t i = 0; i < n; i++){
    //
    R = (Rmax - Rmin)/(n-1)*i + Rmin;
    R1 = R;
    R2 = R;
    //
    rc = 6*2*R1*R2*dt/dalpha_reduction;
    gr->SetPoint(i,R,rc);
  }
}

void coincident_three_rate( TGraph *gr, Double_t dt, Double_t dalpha, Double_t Rmin, Double_t Rmax, Int_t n){
  Double_t R;
  Double_t rc;
  Double_t R1;
  Double_t R2;
  Double_t R3;
  Double_t dalpha_reduction = 4.0/dalpha/dalpha;  
  for( Int_t i = 0; i < n; i++){
    //
    R = (Rmax - Rmin)/(n-1)*i + Rmin;
    R1 = R;
    R2 = R;
    R3 = R;
    //
    rc = 4*3*R1*R2*R2*dt*dt/dalpha_reduction/dalpha_reduction;
    gr->SetPoint(i,R,rc);
  }
}
