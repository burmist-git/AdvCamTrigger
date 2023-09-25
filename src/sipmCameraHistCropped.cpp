//my
#include "sipmCameraHistCropped.hh"
#include "sipmCameraHist.hh"

//c, c++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <time.h>
#include <math.h>
#include <vector>

//root
#include <TVector2.h>
#include <TPolyLine.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TText.h>
#include <TMath.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TCrown.h>
#include <TArc.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TPad.h>
#include <TString.h>
#include <TFile.h>
#include <TAxis.h>
#include <TVector2.h>
#include <TImage.h>
#include <TColor.h>

using namespace std;

sipmCameraHistCropped::sipmCameraHistCropped(const char* name, const char* title, const sipmCameraHist *sipmHist) : TH2Poly()
{
  //
  SetName(name);
  SetTitle(title);
  //
  _name = name;
  _title = title;
  //
  //Double_t phi0 = 4.900;
  //Double_t phi0 = 0;
  //Double_t phi_max = phi0 + 0.175/3.0;
  //Double_t phi_min = phi0 - 0.175/3.0;
  Double_t x_mean = 0.0;
  Double_t y_mean = 0.0;
  //
  _n_pixels = 0;
  //
  //for(unsigned int i = 0;i<sipmHist->get_pixel_vec().size();i++){
  //
  //if(sipmHist->get_pixel_vec().at(i).pix_phi>phi_min &&
  //sipmHist->get_pixel_vec().at(i).pix_phi<phi_max){
  //
  //AddBin(sipmHist->get_pixel_vec().at(0).n,
  //	   sipmHist->get_pixel_vec().at(i).xp,
  //	   sipmHist->get_pixel_vec().at(i).yp);
  //_n_pixels++; 
  //}
  //
  vector<unsigned int> pixel_seed_v;
  pixel_seed_v.push_back(5439);
  pixel_seed_v.push_back(5467);
  //pixel_seed_v.push_back(5483);
  //pixel_seed_v.push_back(5473);
  //
  unsigned int pix_id;
  //
  for(unsigned int i = 0;i<pixel_seed_v.size();i++){
    pix_id = pixel_seed_v.at(i);
    if(check_map(pix_id))
      _pixel_map.push_back(pix_id);
    for(unsigned int j = 0;j<sipmHist->get_pixel_vec().at(pixel_seed_v.at(i)).v_pixel_super_flower.size();j++){
      pix_id = (unsigned int)sipmHist->get_pixel_vec().at(pixel_seed_v.at(i)).v_pixel_super_flower.at(j).pixel_id;
      if(check_map(pix_id))
	_pixel_map.push_back(pix_id);
    }
  }  
  //
  for( unsigned int i = 0; i < _pixel_map.size(); i++){
    x_mean += sipmHist->get_pixel_vec().at(_pixel_map.at(i)).x;
    y_mean += sipmHist->get_pixel_vec().at(_pixel_map.at(i)).y;
  }
  //
  x_mean /= _pixel_map.size();
  y_mean /= _pixel_map.size();
  //
  for( unsigned int i = 0; i < _pixel_map.size(); i++){
    Double_t *xp = new Double_t[sipmHist->get_pixel_vec().at(_pixel_map.at(i)).n];
    Double_t *yp = new Double_t[sipmHist->get_pixel_vec().at(_pixel_map.at(i)).n];
    for(unsigned int j = 0; j<(unsigned int)sipmHist->get_pixel_vec().at(_pixel_map.at(i)).n; j++){
      xp[j] = sipmHist->get_pixel_vec().at(_pixel_map.at(i)).xp[j] - x_mean;
      yp[j] = sipmHist->get_pixel_vec().at(_pixel_map.at(i)).yp[j] - y_mean;
    }
    AddBin(sipmHist->get_pixel_vec().at(_pixel_map.at(i)).n,xp,yp);
  }
  _n_pixels = _pixel_map.size();
  cout<<"_n_pixels = "<<_n_pixels<<endl;
}

bool sipmCameraHistCropped::check_map(unsigned int pix_id){
  if(_pixel_map.size() == 0)
    return true;
  for(unsigned int i = 0; i < _pixel_map.size(); i++){
    if(_pixel_map.at(i) == pix_id)
      return false;
  }
  return true;
}

sipmCameraHistCropped::~sipmCameraHistCropped(){
}

void sipmCameraHistCropped::Clean(){
  for(Int_t i = 0;i<=GetNcells();i++){
    SetBinContent(i,0);
  }
}

void sipmCameraHistCropped::test(TString pdf_out_name){
  TRandom3 *rnd = new TRandom3(123123); 
  //cout<<"GetN() "<<GetNcells()<<endl;
  for(Int_t i = 0; i < GetNcells(); i++)
    SetBinContent(i,(Int_t)rnd->Uniform(1,10));
  //
  Draw_cam("ZCOLOR",pdf_out_name.Data());
}

void sipmCameraHistCropped::test01(const sipmCameraHist *sipmHist, TString pdf_out_name){
  //
  Double_t phi0 = 4.900;
  //Double_t phi0 = 0;
  Double_t phi_max = phi0 + 0.175/3.0;
  Double_t phi_min = phi0 - 0.175/3.0;
  //
  Clean();
  //
  for(unsigned int i = 0;i<sipmHist->get_pixel_vec().size();i++){
    if(sipmHist->get_pixel_vec().at(i).pix_phi>phi_min &&
       sipmHist->get_pixel_vec().at(i).pix_phi<phi_max){
      Fill(sipmHist->get_pixel_vec().at(i).x,sipmHist->get_pixel_vec().at(i).y,i);
    }
  }
  //
  SetMaximum(10000.0);
  SetMarkerSize(0.1);
  //SetLineWidth(0);
  Draw_cam("ZCOLOR TEXT",pdf_out_name.Data());
}

void sipmCameraHistCropped::test02(TString pdf_out_name){
  Clean();
  for(Int_t i = 0;i<=GetNcells();i++)
    SetBinContent(i+1,i);
  SetMarkerSize(0.8);
  Draw_cam("ZCOLOR TEXT",pdf_out_name.Data());
}

void sipmCameraHistCropped::Draw_cam(TString settings,
				     TString pdf_out_file){
  //
  //Double_t lx_camera = 2.5;
  //Double_t ly_camera = 2.5;
  Double_t lx_camera = 0.25;
  Double_t ly_camera = 0.25;
  Double_t d_frame = 0.1;
  //
  //gStyle->SetPalette(kRainBow);
  //gStyle->SetPalette(kCool);
  //gStyle->SetPalette(kIsland);
  //gStyle->SetPalette(kCherry);
  //TColor::InvertPalette();
  //
  gStyle->SetPalette(kInvertedDarkBodyRadiator);
  //
  gStyle->SetOptStat(kFALSE);
  SetTitle("");
  SetName("");
  //
  TCanvas *c1 = new TCanvas("c1","c1",700,700);
  //
  c1->SetRightMargin(0.12);
  c1->SetLeftMargin(0.12);
  c1->SetTopMargin(0.1);
  c1->SetBottomMargin(0.15);
  //
  //gPad->SetGridx();
  //gPad->SetGridy();
  //gPad->SetLogz();
  //
  //SetMaximum(500.0);
  //SetMinimum(300.0);
  //SetMinimum(280.0);
  //SetMinimum(299.0);
  //SetMaximum(308.0);
  //
  TH2F *frame = new TH2F( "h2", "h2", 40, -lx_camera/2.0-d_frame,lx_camera/2.0+d_frame,40, -ly_camera/2.0-d_frame,ly_camera/2.0+d_frame);
  frame->SetTitle("");
  frame->GetXaxis()->SetTitle("x, m");
  frame->GetYaxis()->SetTitle("y, m");
  frame->GetXaxis()->CenterTitle();
  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetTitleOffset(1.5);
  frame->SetStats(kFALSE);
  frame->Draw();
  //
  //settings += " same TEXT";
  settings += " same";
  //
  Draw(settings.Data());
  //
  if(pdf_out_file != "")
    c1->SaveAs(pdf_out_file.Data());
}
