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

Int_t plots_shift_scan(){
  //
  TString fileN01;
  TString fileN02;
  TString fileN03;
  TString fileN04;
  TString fileN05;
  TString fileN06;
  TString fileN07;
  TString fileN08;
  fileN01 = "./hist_gamma_diffuse_nsb_1x_PCAp_70.0phizero.root";
  fileN02 = "./hist_gamma_diffuse_nsb_1x_PCAp_80.0phizero.root";
  fileN03 = "./hist_gamma_diffuse_nsb_1x_PCAp_90.0phizero.root";
  fileN04 = "./hist_gamma_diffuse_nsb_1x_PCAp_100.0phizero.root";
  fileN05 = "./hist_gamma_diffuse_nsb_1x_PCAp_110.0phizero.root";
  fileN06 = "./hist_gamma_diffuse_nsb_1x_PCAp_120.0phizero.root";
  fileN07 = "./hist_gamma_diffuse_nsb_1x_PCAp_130.0phizero.root";
  fileN08 = "./hist_gamma_diffuse_nsb_1x_PCAp_140.0phizero.root";
  TFile *f01 = new TFile(fileN01.Data());
  TFile *f02 = new TFile(fileN02.Data());
  TFile *f03 = new TFile(fileN03.Data());
  TFile *f04 = new TFile(fileN04.Data());
  TFile *f05 = new TFile(fileN05.Data());
  TFile *f06 = new TFile(fileN06.Data());
  TFile *f07 = new TFile(fileN07.Data());
  TFile *f08 = new TFile(fileN08.Data());

  TH2Poly *sc_shift_01 = (TH2Poly*)f01->Get("sipm_cam_shift");
  TH2Poly *sc_norot_01 = (TH2Poly*)f01->Get("sipm_cam_norot");
  TH2Poly *sc_shift_02 = (TH2Poly*)f02->Get("sipm_cam_shift");
  TH2Poly *sc_norot_02 = (TH2Poly*)f02->Get("sipm_cam_norot");
  TH2Poly *sc_shift_03 = (TH2Poly*)f03->Get("sipm_cam_shift");
  TH2Poly *sc_norot_03 = (TH2Poly*)f03->Get("sipm_cam_norot");
  TH2Poly *sc_shift_04 = (TH2Poly*)f04->Get("sipm_cam_shift");
  TH2Poly *sc_norot_04 = (TH2Poly*)f04->Get("sipm_cam_norot");
  TH2Poly *sc_shift_05 = (TH2Poly*)f05->Get("sipm_cam_shift");
  TH2Poly *sc_norot_05 = (TH2Poly*)f05->Get("sipm_cam_norot");
  TH2Poly *sc_shift_06 = (TH2Poly*)f06->Get("sipm_cam_shift");
  TH2Poly *sc_norot_06 = (TH2Poly*)f06->Get("sipm_cam_norot");
  TH2Poly *sc_shift_07 = (TH2Poly*)f07->Get("sipm_cam_shift");
  TH2Poly *sc_norot_07 = (TH2Poly*)f07->Get("sipm_cam_norot");
  TH2Poly *sc_shift_08 = (TH2Poly*)f08->Get("sipm_cam_shift");
  TH2Poly *sc_norot_08 = (TH2Poly*)f08->Get("sipm_cam_norot");
  
  TCanvas *c1 = new TCanvas("c1","c1",10,10,200,1000);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);


  c1->Divide(2,8);
  c1->cd(1);
  sc_norot_01->Draw("ZCOLOR");
  c1->cd(2);
  sc_shift_01->Draw("ZCOLOR");
  c1->cd(3);
  sc_norot_02->Draw("ZCOLOR");
  c1->cd(4);
  sc_shift_02->Draw("ZCOLOR");
  c1->cd(5);
  sc_norot_03->Draw("ZCOLOR");
  c1->cd(6);
  sc_shift_03->Draw("ZCOLOR");
  c1->cd(7);
  sc_norot_04->Draw("ZCOLOR");
  c1->cd(8);
  sc_shift_04->Draw("ZCOLOR");
  c1->cd(9);
  sc_norot_05->Draw("ZCOLOR");
  c1->cd(10);
  sc_shift_05->Draw("ZCOLOR");
  c1->cd(11);
  sc_norot_06->Draw("ZCOLOR");
  c1->cd(12);
  sc_shift_06->Draw("ZCOLOR");
  c1->cd(13);
  sc_norot_07->Draw("ZCOLOR");
  c1->cd(14);
  sc_shift_07->Draw("ZCOLOR");
  c1->cd(15);
  sc_norot_08->Draw("ZCOLOR");
  c1->cd(16);
  sc_shift_08->Draw("ZCOLOR");
  
  return 0;
}
