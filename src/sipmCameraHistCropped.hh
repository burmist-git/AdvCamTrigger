#pragma once

//my
#include "sipmCameraHist.hh"

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

class sipmCameraHist;
class TPrincipal;

class sipmCameraHistCropped: public TH2Poly {
  
public:
  
  sipmCameraHistCropped( const char* name, const char* title, const sipmCameraHist *sipmHist, bool if_centrate = false);
  sipmCameraHistCropped( const char* name, const char* title, const sipmCameraHist *sipmHist, TString in_map_file_name);
  sipmCameraHistCropped( const char* name, const char* title, const sipmCameraHist *sipmHist, const std::vector<unsigned int> &pixel_map_v);
  
  ~sipmCameraHistCropped();
  //
  static const Int_t _dd_im = 887;
  //
  void Clean();
  void test(TString pdf_out_name = "sipmCameraHistCropped_test.pdf");
  void test01(const sipmCameraHist *sipmHist, TString pdf_out_name = "sipmCameraHistCropped_test01.pdf");
  void test02(TString pdf_out_name = "sipmCameraHistCropped_test02.pdf");
  void Draw_cam(TString settings, TString pdf_out_file);
  //  
  TString _name;
  TString _title;
  //  
  void Fill_pe(const Int_t npixels_n, const Int_t *pix_id, const Float_t *pix_pe_time,
	       const Double_t ev_time, const Double_t time_offset, const Double_t alpha, TPrincipal *principal = NULL, bool if_centrate = false, Double_t x_shift = 0.0, Double_t y_shift = 0.0);
  //  
  inline const unsigned int get_n_pixels() const {return _n_pixels;}
  inline const std::vector<unsigned int> &get_pixel_map() const {return _pixel_map;}
  //
  void Fill_principal( std::vector<sipmCameraHistCropped*> &sipm_cam_hist_v, const TPrincipal *principal);
  void Fill_principal( std::vector<sipmCameraHistCropped*> &sipm_cam_principal_hist_v, Double_t data_Vh[_dd_im][_dd_im]);

  static void save_pixel_map( const std::vector<unsigned int> &pixel_map_v, TString out_map_file_name="sipmCameraHistCropped_pix.map");
  static void load_pixel_map( std::vector<unsigned int> &pixel_map_v, TString in_map_file_name="sipmCameraHistCropped_pix.map");

  void Save_to_csv(TString csvname, const std::vector<sipmCameraHistCropped*> simp_hist_crop_v);

  void draw_crop_vector( Int_t nx, Int_t ny, const std::vector<sipmCameraHistCropped*> simp_hist_crop_v, TCanvas *c1);
  
private:

  bool check_map(unsigned int pix_id);
  
  unsigned int _n_pixels;
  std::vector<unsigned int> _pixel_map;

  const sipmCameraHist *_sipm_cam;

  void build_pixel_seed(std::vector<unsigned int> &pixel_seed_v);

};
