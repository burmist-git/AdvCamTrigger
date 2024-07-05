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

Double_t power_law(Double_t *x, Double_t *par){
  Double_t A  = par[0];
  Double_t B  = par[1];
  Double_t val = A/(TMath::Power(x[0],B));
  return val;
}

Int_t plots_power_law_dist(){

  Double_t eMin = 10.0;
  Double_t eMax = 100000.0;
  
  TH1D *h1 = new TH1D("h1","h1",10000,eMin,eMax);
  //
  TRandom3 *rnd = new TRandom3(123123);
  //for(Int_t i = 0;i<100000000;i++)
  //h1->Fill(TMath::Sqrt(1.0/(rnd->Uniform(1.0/eMax/eMax,1.0/eMin/eMin))));
  for(Int_t i = 0;i<1000000000;i++)
    h1->Fill(1.0/(rnd->Uniform(1.0/eMax,1.0/eMin)));
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //
  //h1->SetMinimum(0);
  //h1->Draw("errors");
  h1->Draw();
  //
  const Int_t npar = 2;
  Double_t inParameters[npar];
  inParameters[0] = 1.0e+10;
  inParameters[1] = 2.0;
  //
  TF1 *f_power_law = new TF1( "f_power_law", power_law, eMin, eMax, npar);
  f_power_law->SetParameters(inParameters);
  f_power_law->SetParName(0, "A");
  f_power_law->SetParName(1, "B");
  //f_power_law->FixParameter(0,inParameters[0]);
  //f_power_law->FixParameter(1,inParameters[1]);
  h1->Fit("f_power_law","","",eMin+20, eMax-1000);


  //
  return 0;
}
