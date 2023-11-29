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

Int_t plots_unirandom(){

  TH1D *h1 = new TH1D("h1","h1",3,-1.5,1.5);
  //
  TRandom3 *rnd = new TRandom3(123123);
  for(Int_t i = 0;i<100000;i++)
    h1->Fill(((Int_t)rnd->Uniform(0,3)-1));
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //
  h1->SetMinimum(0);
  h1->Draw("errors");
  //
  return 0;
}
