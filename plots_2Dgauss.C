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

Int_t plots_2Dgauss(){
  //
  TRandom3 *rnd = new TRandom3(1231231);
  //
  Int_t nn = 10000000;
  //
  TH1D *h1_x = new TH1D("h1_x","h1_x",200, -0.0005, 0.0005);
  TH1D *h1_y = new TH1D("h1_y","h1_y",200, -0.0005, 0.0005);
  TH1D *h1_ax = new TH1D("h1_ax","h1_ax",200, -0.02, 0.02);
  TH1D *h1_ay = new TH1D("h1_ay","h1_ay",200, -0.02, 0.02);
  TH2D *h2_ay_vs_ax = new TH2D("h2_ay_vs_ax","h2_ay_vs_ax",200, -0.02, 0.02, 200, -0.02, 0.02);
  TH1D *h1_r = new TH1D("h1_r","h1_r",800,  0.0, 0.001);
  TH1D *h1_sum = new TH1D("h1_sum","h1_sum",1600,  -0.001, 0.001);
  TH1D *h1_theta = new TH1D("h1_theta","h1_theta",1600,  0.0, 0.02);
  TH2D *h2_y_vs_x = new TH2D("h2_y_vs_x","h2_y_vs_x",200, -0.0005, 0.0005, 200, -0.0005, 0.0005);
  //
  Double_t x;
  Double_t y;
  Double_t ax;
  Double_t ay;
  Double_t r;
  Double_t sum;
  Double_t theta;
  for(Int_t i = 0;i<nn;i++){
    ax  = rnd->Gaus(0.0,0.0046);
    ay  = rnd->Gaus(0.0,0.0046);
    x = TMath::Sin(ax/180.0*TMath::Pi());
    y = TMath::Sin(ay/180.0*TMath::Pi());
    //
    h1_ax->Fill(ax);
    h1_ay->Fill(ay);
    //
    r = TMath::Sqrt(x*x + y*y);
    theta = TMath::ASin(r)*180.0/TMath::Pi();
    sum = x + y;
    //
    h1_r->Fill(r);
    h1_sum->Fill(sum);
    h1_x->Fill(x);
    h1_y->Fill(y);
    h2_y_vs_x->Fill(x,y);
    h2_ay_vs_ax->Fill(ax,ay);
    h1_theta->Fill(theta);
  }
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,900,900);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE); 
  //
  c1->Divide(3,3);
  c1->cd(1);
  h1_x->Draw();
  c1->cd(2);
  h1_y->Draw();
  c1->cd(3);
  h2_y_vs_x->Draw("ZCOLOR");
  c1->cd(4);
  h1_ax->Draw();
  c1->cd(5);
  h1_ay->Draw();
  c1->cd(6);
  h2_ay_vs_ax->Draw("ZCOLOR");
  c1->cd(7);
  h1_sum->Draw();
  c1->cd(8);
  h1_r->Draw();
  //
  c1->cd(9);
  h1_theta->Draw();
  //
  return 0;
}
