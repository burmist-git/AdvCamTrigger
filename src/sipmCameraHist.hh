#pragma once

//root
#include <TObject.h>
#include <TH2Poly.h>
#include <TGraph.h>
#include <TVector2.h>
#include <TCanvas.h>
#include <TMath.h>

//c, c++
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>

struct pixel_info {
  //
  Int_t pixel_id;
  Float_t x;
  Float_t y;
  Int_t drawer_id;
  //
  const Int_t n = 7;
  Double_t *xp = new Double_t[n];
  Double_t *yp = new Double_t[n];
  //
  pixel_info(){
    //
    pixel_id = -999;
    x = -999.0;
    y = -999.0;
    drawer_id = -999;
    //
    for(Int_t i = 0;i<n;i++){
      xp[i] = -999.0;
      yp[i] = -999.0;
    }
    //
  };
  void print_info(){
    std::cout<<std::setw(10)<<pixel_id
	     <<std::setw(10)<<x
      	     <<std::setw(10)<<y
	     <<std::setw(10)<<drawer_id
	     <<std::endl;
  };
  static void print_info_header(){
    std::cout<<std::setw(10)<<"pixel_id"
	     <<std::setw(10)<<"x"
      	     <<std::setw(10)<<"y"
	     <<std::setw(10)<<"drawer_id"
	     <<std::endl;
  };
  void rotatePix(Double_t rot_alpha_deg){
    if(rot_alpha_deg != 0.0){
      TVector2 v(x,y);
      x = v.Rotate(rot_alpha_deg/180.0*TMath::Pi()).X();
      y = v.Rotate(rot_alpha_deg/180.0*TMath::Pi()).Y();
    }
  };
  void build_Cell(Int_t cell_type_id, Double_t l){
    if(cell_type_id == 0){
      Double_t alpha   = 2.0*TMath::Pi()/6.0;
      Double_t alpha_2 = alpha/2.0;
      Double_t alpha0  = alpha_2;
      Double_t l_2 = l/2.0;
      Double_t r = l_2/TMath::Cos(alpha_2);
      Double_t theta = 0.0;
      for(Int_t i = 0;i<n;i++){
	theta = alpha0 + alpha*i;
	xp[i] = r*TMath::Cos(theta) + x;
	yp[i] = r*TMath::Sin(theta) + y;
      }
    }
    else{
      std::cout<<"  ---> ERROR : cell_type_id = "<<cell_type_id<<std::endl;
      assert(0);
    }
  }
};

class sipmCameraHist: public TH2Poly {
 public:
  
  sipmCameraHist(const char* name, const char* title, const char* mapping_csv_file, Double_t rot_alpha_deg);
  ~sipmCameraHist();
  void dump_mapping_info();
  void test();
  void test02();
  void test03();
  void test_drawer_id();
  void Clean();
  void count_signal(Double_t th_val, Int_t &nch, Int_t &npe);
  void Draw_cam( TString settings, TString pdf_out_file);
  void Draw_cam( TString settings,
		 TString pdf_out_file,
		 TString particle_type,
		 Int_t wf_time_id,
		 Int_t   event_id,
		 Float_t energy,
		 Float_t xcore,
		 Float_t ycore,
		 Float_t ev_time,
		 Int_t nphotons,
		 Int_t n_pe,
		 Int_t n_pixels);
  //  
  TString _name;
  TString _title;
  
 private:

  static const unsigned int _n_pixels = 7987;
  static const unsigned int _n_drawers = 1141;
  double _pixel_size;
  std::vector<pixel_info> _pixel_vec;
  void load_mapping(const char* mapping_csv_file);
  Double_t _rot_alpha_deg;
  
};
