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
void coincident_rate_exp( TGraph *gr, Double_t dt, Double_t Rmin, Double_t Rmax, Int_t n);
Double_t calculate_coincidence_rate(Double_t dt, Double_t R, bool def = true);

Int_t plots_coincident_rate_stereo_2Dplot(){
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
  TGraph *gr_cr_100ns = new TGraph();
  gr_cr_100ns->SetNameTitle("gr_rc_100ns","gr_rc_100ns");
  TGraph *gr_crexp_100ns = new TGraph();
  gr_crexp_100ns->SetNameTitle("gr_rcexp_100ns","gr_rcexp_100ns");
  //
  //
  coincident_rate( gr_cr_100ns, 100*1.0e-9, 0.1*1.0e+5, 5.0*1.0e+5, 1000);
  coincident_rate_exp( gr_crexp_100ns, 100*1.0e-9, 0.1*1.0e+5, 5.0*1.0e+5, 1000);
  //
  //
  gr_cr_100ns->SetLineColor(kBlack);
  gr_cr_100ns->SetLineWidth(2.0);
  //
  gr_crexp_100ns->SetLineColor(kRed+2);
  gr_crexp_100ns->SetLineWidth(2.0);
  //
  //
  //
  //gr->Draw("APL");
  //
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr_cr_100ns);
  mg->Add(gr_crexp_100ns);
  //mg->Add(gr_70ns_07deg);
  //mg->Add(gr_80ns_1deg_three);
  //mg->Add(gr_70ns_07deg_three);
  //
  //mg->SetMinimum(1.1e+9);
  //mg->SetMaximum(1.0e+8);
  //
  mg->Draw("al");
  //
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  //
  mg->GetXaxis()->SetTitle("Rate, Hz");
  mg->GetYaxis()->SetTitle("Random coincidence, Hz");
  //mg->GetXaxis()->SetTitle("L0 rate per LST, Hz");
  //mg->GetYaxis()->SetTitle("Total fake rate of the 4 LSTs, Hz");
  //
  //
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(gr_cr_100ns,    "#tau = 100 ns Formula: 2 #times R^{2} #times #tau", "al");
  leg->AddEntry(gr_crexp_100ns, "#tau = 100 ns Formula: R #times (1 - exp(-2 #times R #times #tau)) ", "al");
  //leg->AddEntry(gr_70ns_07deg, "#tau = 70 ns and 0.7 deg topological trigger", "al");
  //leg->AddEntry(gr_80ns_1deg_three,  "#tau = 80 ns and 1.0 deg topological trigger (triple coinc.)", "al");
  //leg->AddEntry(gr_70ns_07deg_three, "#tau = 70 ns and 0.7 deg topological trigger (triple coinc.)", "al");
  leg->Draw();
  //
  //
  return 0;
}

void coincident_rate( TGraph *gr, Double_t dt, Double_t Rmin, Double_t Rmax, Int_t n){
  Double_t R;
  Double_t rc;
  for( Int_t i = 0; i < n; i++){
    //
    R = (Rmax - Rmin)/(n-1)*i + Rmin;
    //
    rc = calculate_coincidence_rate(dt,  R);
    gr->SetPoint(i,R,rc);
  }
}

void coincident_rate_exp( TGraph *gr, Double_t dt, Double_t Rmin, Double_t Rmax, Int_t n){
  Double_t R;
  Double_t rc;
  for( Int_t i = 0; i < n; i++){
    //
    R = (Rmax - Rmin)/(n-1)*i + Rmin;
    //
    rc = calculate_coincidence_rate(dt,  R, false);
    gr->SetPoint(i,R,rc);
  }
}

Double_t calculate_coincidence_rate(Double_t dt, Double_t R, bool def = true){
  if(def)
    return 2*R*R*dt;
  return R*(1.0 - TMath::Exp(-R*2.0*dt));
}
