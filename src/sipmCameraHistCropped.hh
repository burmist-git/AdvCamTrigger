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

class sipmCameraHistCropped: public TH2Poly {
  
public:
  
  sipmCameraHistCropped( const char* name, const char* title, const sipmCameraHist *sipmHist);
  ~sipmCameraHistCropped();
  //
  void Clean();
  void test(TString pdf_out_name = "sipmCameraHistCropped_test.pdf");
  void test01(const sipmCameraHist *sipmHist, TString pdf_out_name = "sipmCameraHistCropped_test01.pdf");
  void test02(TString pdf_out_name = "sipmCameraHistCropped_test02.pdf");
  void Draw_cam(TString settings, TString pdf_out_file);
  //  
  TString _name;
  TString _title;
  
 private:

  bool check_map(unsigned int pix_id);
  
  unsigned int _n_pixels;
  std::vector<unsigned int> _pixel_map;
  
};
