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

Int_t plots_shift_fit(){

  TString fileN;
  fileN = "./hist_gamma_diffuse_nsb_1x_PCAp.root";
  //fileN = "./hist_gamma_on_nsb_1x_PCAp.root";

  TFile *f1 = new TFile(fileN.Data());

  TH2Poly *sc_shift = (TH2Poly*)f1->Get("sipm_cam_shift");
  TH2Poly *sc_norot = (TH2Poly*)f1->Get("sipm_cam_norot");

  TCanvas *c1 = new TCanvas("c1",fileN.Data(),10,10,1200,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);


  c1->Divide(2,1);
  c1->cd(1);
  sc_norot->Draw("ZCOLOR");
  c1->cd(2);
  sc_shift->Draw("ZCOLOR");
  
  return 0;
}
