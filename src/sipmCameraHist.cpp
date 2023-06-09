//my
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

using namespace std;

sipmCameraHist::sipmCameraHist(const char* name, const char* title, const char* mapping_csv_file, Double_t rot_alpha_deg) : TH2Poly(), _rot_alpha_deg(rot_alpha_deg)
{
  //
  _pixel_size = 0.02329999953508377;
  //
  load_mapping( mapping_csv_file);
  //
  SetName(name);
  SetTitle(title);
  //
   _name = name;
   _title = title;
  //
  for(unsigned int i = 0;i<_pixel_vec.size();i++)
    AddBin(_pixel_vec.at(0).n,_pixel_vec.at(i).xp,_pixel_vec.at(i).yp);
}

void sipmCameraHist::load_mapping(const char* mapping_csv_file){
  //
  ifstream fFile(mapping_csv_file);
  cout<<mapping_csv_file<<std::endl;
  //
  Float_t x, y, drawer_id;
  Int_t pixel_id = 0;  
  //
  if(fFile.is_open()){
    while(fFile>>x>>y>>drawer_id){
      pixel_info pix_i;      
      pix_i.pixel_id = pixel_id;
      pix_i.x = x;
      pix_i.y = y;
      pix_i.drawer_id = (Int_t)drawer_id;
      pix_i.rotatePix(_rot_alpha_deg);
      //
      pix_i.build_Cell(0, _pixel_size);
      //
      pixel_id++;
      _pixel_vec.push_back(pix_i);
    }
    fFile.close();
  }
  //
}

void sipmCameraHist::dump_mapping_info(){
  pixel_info::print_info_header();
  for(unsigned int i = 0;i<_pixel_vec.size();i++)
    _pixel_vec.at(i).print_info();
  //
  cout<<"_n_pixels       "<<_n_pixels<<endl
      <<"_n_drawers      "<<_n_drawers<<endl
      <<"_name           "<<_name<<endl
      <<"_title          "<<_title<<endl
      <<"_rot_alpha_deg  "<<_rot_alpha_deg<<endl;
  //
}
 
sipmCameraHist::~sipmCameraHist(){
}

void sipmCameraHist::count_signal(Double_t th_val, Int_t &nch, Int_t &npe){
  /*
  nch = 0;
  npe = 0;
  for(Int_t i = 0;i<GetNcells();i++){
    //cout<<GetNcells()<<endl;
    if(GetBinContent(i)>=th_val){
      npe = npe + GetBinContent(i);
      nch++;
    }
  }
  */
}

void sipmCameraHist::Clean(){
  for(Int_t i = 0;i<GetNcells();i++){
    SetBinContent(i,0);
  }
}

void sipmCameraHist::Draw_cam( TString settings, TString pdf_out_file){
  //
  Double_t lx_camera = 2.5;
  Double_t ly_camera = 2.5;
  Double_t d_frame = 0.1;
  //
  gStyle->SetPalette(kRainBow);
  gStyle->SetOptStat(kFALSE);
  SetTitle("");
  SetName("");
  //
  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  c1->SetRightMargin(0.12);
  c1->SetLeftMargin(0.12);
  c1->SetTopMargin(0.1);
  c1->SetBottomMargin(0.15);
  //
  c1->SetGridx();
  c1->SetGridy();
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
  settings += " same";
  Draw(settings.Data());
  if(pdf_out_file != "")
    c1->SaveAs(pdf_out_file.Data());
}

void sipmCameraHist::test(){
  TRandom3 *rnd = new TRandom3(123123); 
  //cout<<"GetN() "<<GetNcells()<<endl;
  for(Int_t i = 0;i<GetNcells();i++)
    SetBinContent(i,(Int_t)rnd->Uniform(1,10));
  //
  Draw_cam("ZCOLOR","sipmCameraHist_test.pdf");
}
