#pragma once

//root
#include <TObject.h>
#include <TH2Poly.h>
#include <TGraph.h>
#include <TVector2.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLine.h>

//c, c++
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>

class TVector3;

struct mirror_info {
  //
  Int_t mirror_id;
  Float_t x_down;
  Float_t y_right_cam;
  Float_t x;
  Float_t y;
  Float_t z;
  Float_t flat_to_flat;
  Float_t focal_length;
  Float_t phi;
  Float_t r;
  Float_t n_x;
  Float_t n_y;
  Float_t n_z;
  Float_t dist_of_optical_axis_to_mirror_center_cm;
  //
  const Int_t n = 7;
  Double_t *xp = new Double_t[n];
  Double_t *yp = new Double_t[n];
  //
  mirror_info(){
    //
    mirror_id = -999;
    x_down = -999.0;
    y_right_cam = -999.0;
    x = -999.0;
    y = -999.0;
    z = -999.0;
    flat_to_flat = -999.0;
    focal_length = -999.0;
    phi = -999.0;
    r = -999.0;
    //
    n_x = -999.0;
    n_y = -999.0;
    n_z = -999.0;
    dist_of_optical_axis_to_mirror_center_cm = -999.0;
    //
    for(Int_t i = 0;i<n;i++){
      xp[i] = -999.0;
      yp[i] = -999.0;
    }
  }
  void test(){
    //
    mirror_id = 0;
    x_down = 0.0;
    y_right_cam = 0.0;
    x = 0.0;
    y = 0.0;
    z = 0.0;
    flat_to_flat = 151.0; //cm
    focal_length = 2920;  //cm
    phi = 0.0;
    r = 0.0;
    //
    for(Int_t i = 0;i<n;i++){
      xp[i] = -999.0;
      yp[i] = -999.0;
    }
  }
  void print_info(){
    std::cout<<std::setw(15)<<mirror_id
	     <<std::setw(15)<<x_down
      	     <<std::setw(15)<<y_right_cam
	     <<std::setw(15)<<x
      	     <<std::setw(15)<<y
	     <<std::setw(15)<<flat_to_flat
      	     <<std::setw(15)<<focal_length
	     <<std::setw(15)<<phi
      	     <<std::setw(15)<<r
	     <<std::setw(15)<<n_x
	     <<std::endl;
  }
  static void print_info_header(){
    std::cout<<std::setw(15)<<"mirror_id"
	     <<std::setw(15)<<"x_down"
      	     <<std::setw(15)<<"y_right_cam"
	     <<std::setw(15)<<"x"
      	     <<std::setw(15)<<"y"
	     <<std::setw(15)<<"flat_to_flat"
      	     <<std::setw(15)<<"focal_length"
	     <<std::setw(15)<<"phi"
      	     <<std::setw(15)<<"r"
	     <<std::setw(15)<<"n_x"
	     <<std::endl;
  }
  //
  //
  //void rotatePix(Double_t rot_alpha_deg){
  //  if(rot_alpha_deg != 0.0){
  //   TVector2 v(x,y);
  //   x = v.Rotate(rot_alpha_deg/180.0*TMath::Pi()).X();
  //   y = v.Rotate(rot_alpha_deg/180.0*TMath::Pi()).Y();
  //}
  //}
  void build_Cell(){
    Double_t alpha   = 2.0*TMath::Pi()/6.0;
    Double_t alpha_2 = alpha/2.0;
    Double_t alpha0  = alpha_2;
    Double_t l_2 = flat_to_flat/2.0;
    Double_t rr = l_2/TMath::Cos(alpha_2);
    Double_t theta = 0.0;
    for(Int_t i = 0;i<n;i++){
      theta = alpha0 + alpha*i;
      xp[i] = rr*TMath::Cos(theta) + x;
      yp[i] = rr*TMath::Sin(theta) + y;
    }
  }
  void get_mirror_x0_y0_from_sim_telarray_cfg(){
    Double_t alpha   = 2.0*TMath::Pi()/6.0;
    Double_t alpha_2 = alpha/2.0;
    Double_t alpha0  = alpha_2;
    Double_t l_2 = flat_to_flat/2.0;
    Double_t rr = l_2/TMath::Cos(alpha_2);
    x = y_right_cam;
    y = x_down - rr/2.0;
  }
};

class lstMirrorHist: public TH2Poly {

public:
  
  lstMirrorHist(Int_t dummy);
  lstMirrorHist(const char* name, const char* title, lstMirrorHist *mirrHist);
  lstMirrorHist(const char* name = "lstMirrorHist",
		const char* title = "lstMirrorHist",
		const char* mapping_csv_file = "mirror_CTA-N-LST1_v2019-03-31_format.dat",
		bool if_mapping_ideal = false,
		Double_t rot_alpha_deg = 0.0);  
  ~lstMirrorHist();
  void dump_mapping_info();
  void Clean();
  void test(TString pdf_out_name = "lstMirrorHist_test.pdf", TString hist_out_name = "hist_lstMirrorHist_test.root");
  void test_ideal(TString pdf_out_name = "lstMirrorHist_test_ideal.pdf", TString hist_out_name = "hist_lstMirrorHist_test_ideal.root");
  TCanvas *Draw_Mirrors(TString settings, TString pdf_out_file);
  //
  inline const std::vector<mirror_info> &get_mirror_vec() const {return _mirror_vec;}
  
  void setHistNameTitle(TString namestr = " ", TString titlestr = " ");

  void save_to_csv_mirror_vec(TString out_file_name = "mirror_CTA-LST-v20141201-198_short.csv");

  void fill_mirror_vec_with_ideal_val(Double_t focus_val_cm);
  Double_t get_ideal_z( Double_t focus_val_cm, Double_t x, Double_t y);
  void get_ideal_nx_ny_nz( Double_t focus_val_cm,
			   Double_t x, Double_t y, Double_t z,
			   Double_t &nx, Double_t &ny, Double_t &nz);
  //
  Double_t get_ideal_curvature( Double_t focus_val_cm, Double_t x, Double_t y);
  Double_t get_spherical_mirror_focus_form_curvature( Double_t curvature);
  Double_t get_ideal_angle_deg_mirror_oa( Double_t focus_val_cm, Double_t x, Double_t y);
  
private:
  
  std::vector<mirror_info> _mirror_vec;
  void load_mapping(const char* mapping_csv_file);
  void load_mapping_ideal(const char* mapping_csv_file);
  //  
  TString _name;
  TString _title;
  Double_t _rot_alpha_deg;

};
