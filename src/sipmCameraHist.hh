#pragma once

//my 
#include "anabase.hh"

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

namespace line_info{
  struct line_flower_info {
    Float_t x;
    Float_t y;
    Float_t phi;
    line_flower_info(){
      x = -999.0;
      y = -999.0;
      phi = -999.0;
    }  
    static void print_info_header(){
      std::cout<<std::setw(15)<<"x"
	       <<std::setw(15)<<"y"
	       <<std::setw(15)<<"phi"
	       <<std::endl;
    }
    void print_info(){
      std::cout<<std::setw(15)<<x
	       <<std::setw(15)<<y
	       <<std::setw(15)<<phi
	       <<std::endl;
    }
    static void print_info(const line_flower_info a){
      std::cout<<std::setw(15)<<a.x
	       <<std::setw(15)<<a.y
	       <<std::setw(15)<<a.phi
	       <<std::endl;
    }
    static void print_info(const std::vector<line_flower_info> &a){
      print_info_header();
      for( unsigned int i = 0; i < a.size(); i++)
	print_info(a.at(i));
    }
    static void bubbleSort(std::vector<line_flower_info> &a){
      bool swapp = true;
      while(swapp){
	swapp = false;
	for( unsigned int i = 0; i < a.size()-1; i++){
	  if (a.at(i).phi>a.at(i+1).phi){
	    //
	    line_flower_info tmp;
	    tmp = a.at(i);
	    //
	    a.at(i) = a.at(i+1);
	    a.at(i+1) = tmp;
	    swapp = true;
	  }
	}
      }
    }
  };
}

struct pixel_neighbors_info {
  Int_t pixel_id;
  Float_t distance_between_pixels;
  Float_t x;
  Float_t y;
  Float_t x0;
  Float_t y0;
  Float_t phi;
  Float_t phi_deg;
  pixel_neighbors_info(){
    pixel_id = -999;
    distance_between_pixels = -999.0;
    x = -999.0;
    y = -999.0;
    x0 = -999.0;
    y0 = -999.0;
    phi = -999.0;
    phi_deg = -999.0;
  }
  static void print_info_header(){
    std::cout<<std::setw(15)<<"pixel_id"
	     <<std::setw(15)<<"x"
      	     <<std::setw(15)<<"y"
	     <<std::setw(15)<<"dist_btw_pix"
      	     <<std::setw(15)<<"phi_deg"
	     <<std::endl;
  }
  void print_info(){
    std::cout<<std::setw(15)<<pixel_id
	     <<std::setw(15)<<x
      	     <<std::setw(15)<<y
	     <<std::setw(15)<<distance_between_pixels
      	     <<std::setw(15)<<phi_deg
	     <<std::endl;
  }
};

struct pixel_info {
  //
  Int_t pixel_id;
  Float_t x;
  Float_t y;
  Float_t pix_phi;
  Float_t pix_r;
  Int_t drawer_id;
  //
  const Int_t n = 7;
  Double_t *xp = new Double_t[n];
  Double_t *yp = new Double_t[n];
  std::vector<pixel_neighbors_info> v_pixel_neighbors;
  std::vector<pixel_neighbors_info> v_pixel_neighbors_second;
  std::vector<pixel_neighbors_info> v_pixel_neighbors_third;
  std::vector<pixel_neighbors_info> v_pixel_flower;
  std::vector<pixel_neighbors_info> v_pixel_super_flower;
  std::vector<TLine> v_line_flower;
  pixel_neighbors_info self_neighbors_info;
  //
  pixel_info(){
    //
    pixel_id = -999;
    x = -999.0;
    y = -999.0;
    drawer_id = -999;
    //
    pix_phi = -999.0;
    pix_r = -999.0;
    //
    for(Int_t i = 0;i<n;i++){
      xp[i] = -999.0;
      yp[i] = -999.0;
    }
    //
  }
  void print_info(){
    std::cout<<std::setw(10)<<pixel_id
	     <<std::setw(10)<<x
      	     <<std::setw(10)<<y
	     <<std::setw(10)<<drawer_id
	     <<std::endl;
  }
  static void print_info_header(){
    std::cout<<std::setw(10)<<"pixel_id"
	     <<std::setw(10)<<"x"
      	     <<std::setw(10)<<"y"
	     <<std::setw(10)<<"drawer_id"
	     <<std::endl;
  }
  void rotatePix(Double_t rot_alpha_deg){
    if(rot_alpha_deg != 0.0){
      TVector2 v(x,y);
      x = v.Rotate(rot_alpha_deg/180.0*TMath::Pi()).X();
      y = v.Rotate(rot_alpha_deg/180.0*TMath::Pi()).Y();
    }
  }
  void build_Cell(Int_t cell_type_id, Double_t l, Float_t x_v, Float_t y_v){
    if(cell_type_id == 0){
      Double_t alpha   = 2.0*TMath::Pi()/6.0;
      Double_t alpha_2 = alpha/2.0;
      Double_t alpha0  = alpha_2;
      Double_t l_2 = l/2.0;
      Double_t r = l_2/TMath::Cos(alpha_2);
      Double_t theta = 0.0;
      for(Int_t i = 0;i<n;i++){
	theta = alpha0 + alpha*i;
	xp[i] = r*TMath::Cos(theta) + x_v;
	yp[i] = r*TMath::Sin(theta) + y_v;
      }
    }
    else{
      std::cout<<"  ---> ERROR : cell_type_id = "<<cell_type_id<<std::endl;
      assert(0);
    }
  }
  void build_Cell(Int_t cell_type_id, Double_t l){
    build_Cell(cell_type_id, l, x, y);
  }
  void find_pixel_neighbors( const std::vector<pixel_info> &pixel_vec, double pixel_pitch){
    find_pixel_neighbors( pixel_vec, pixel_pitch, NULL);
  }
  void find_pixel_neighbors( const std::vector<pixel_info> &pixel_vec, double pixel_pitch, TH1D *h1_distance_between_pixels){
    //
    self_neighbors_info.pixel_id = pixel_id;
    self_neighbors_info.distance_between_pixels = 0.0;
    self_neighbors_info.x = x;
    self_neighbors_info.y = y;
    self_neighbors_info.x0 = 0.0;
    self_neighbors_info.y0 = 0.0;
    self_neighbors_info.phi = 0.0;
    self_neighbors_info.phi_deg = 0.0;
    //    
    //std::vector<pixel_neighbors_info> v_pixel_neighbors;
    Float_t tollerance = 0.3;
    Float_t dl = pixel_pitch*tollerance;
    Float_t distance_between_pixels = 0.0;
    for( unsigned int i = 0; i < pixel_vec.size(); i++){
      //for( unsigned int i = 0; i < 10; i++){
      distance_between_pixels = get_dist_pixel(pixel_vec.at(i).x, pixel_vec.at(i).y);
      if(h1_distance_between_pixels != NULL)
	h1_distance_between_pixels->Fill(distance_between_pixels);
      if(distance_between_pixels>(pixel_pitch-dl) && distance_between_pixels<(pixel_pitch+dl)){
	//std::cout<<"distance_between_pixels  = "<<distance_between_pixels<<std::endl
	//	 <<"pixel_vec.at(i).pixel_id = "<<pixel_vec.at(i).pixel_id<<std::endl
	//	 <<"pixel_vec.at(i).x        = "<<pixel_vec.at(i).x<<std::endl
	//	 <<"pixel_vec.at(i).y        = "<<pixel_vec.at(i).y<<std::endl
	//	 <<"pixel_id                 = "<<pixel_id<<std::endl
	//	 <<"x                        = "<<x<<std::endl
	//	 <<"y                        = "<<y<<std::endl;	
	pixel_neighbors_info pix_neighb_inf;
	pix_neighb_inf.pixel_id = pixel_vec.at(i).pixel_id;
	pix_neighb_inf.x = pixel_vec.at(i).x;
	pix_neighb_inf.y = pixel_vec.at(i).y;
	pix_neighb_inf.distance_between_pixels = distance_between_pixels;
	pix_neighb_inf.x0 = pixel_vec.at(i).x - x;
	pix_neighb_inf.y0 = pixel_vec.at(i).y - y;
	TVector2 nei_v(pix_neighb_inf.x0,pix_neighb_inf.y0);
	pix_neighb_inf.phi = nei_v.Phi();
	pix_neighb_inf.phi_deg = pix_neighb_inf.phi*180.0/TMath::Pi();
	v_pixel_neighbors.push_back(pix_neighb_inf);
	v_pixel_flower.push_back(pix_neighb_inf);
      }      
      if(distance_between_pixels>(2*pixel_pitch-dl) && distance_between_pixels<(2*pixel_pitch+dl)){
	pixel_neighbors_info pix_neighb_inf;
	pix_neighb_inf.pixel_id = pixel_vec.at(i).pixel_id;
	pix_neighb_inf.x = pixel_vec.at(i).x;
	pix_neighb_inf.y = pixel_vec.at(i).y;
	pix_neighb_inf.distance_between_pixels = distance_between_pixels;
	pix_neighb_inf.x0 = pixel_vec.at(i).x - x;
	pix_neighb_inf.y0 = pixel_vec.at(i).y - y;
	TVector2 nei_v(pix_neighb_inf.x0,pix_neighb_inf.y0);
	pix_neighb_inf.phi = nei_v.Phi();
	pix_neighb_inf.phi_deg = pix_neighb_inf.phi*180.0/TMath::Pi();
	v_pixel_neighbors_second.push_back(pix_neighb_inf);
      }
      if(distance_between_pixels>0.062 && distance_between_pixels<0.074){
	pixel_neighbors_info pix_neighb_inf;
	pix_neighb_inf.pixel_id = pixel_vec.at(i).pixel_id;
	pix_neighb_inf.x = pixel_vec.at(i).x;
	pix_neighb_inf.y = pixel_vec.at(i).y;
	pix_neighb_inf.distance_between_pixels = distance_between_pixels;
	pix_neighb_inf.x0 = pixel_vec.at(i).x - x;
	pix_neighb_inf.y0 = pixel_vec.at(i).y - y;
	TVector2 nei_v(pix_neighb_inf.x0,pix_neighb_inf.y0);
	pix_neighb_inf.phi = nei_v.Phi();
	pix_neighb_inf.phi_deg = pix_neighb_inf.phi*180.0/TMath::Pi();
	v_pixel_neighbors_third.push_back(pix_neighb_inf);
      }
    }
    bubbleSort(v_pixel_neighbors_third);
  }
  void build_pixel_super_flower( const std::vector<pixel_info> &pixel_vec){
    for( unsigned int i = 0; i < v_pixel_flower.size(); i++){
      pixel_neighbors_info pix_neighb_inf = v_pixel_flower.at(i);
      v_pixel_super_flower.push_back(pix_neighb_inf);
    }
    std::vector<unsigned int> super_flower_pixel_v;
    //get_flower_pixel_v(super_flower_pixel_v);
    //get_flower_pixel_v(super_flower_pixel_v,true);
    get_flower_pixel_v(super_flower_pixel_v, false);
    for( unsigned int i = 0; i < super_flower_pixel_v.size(); i++){      
      pixel_neighbors_info pix_neighb_inf = pixel_vec.at(super_flower_pixel_v.at(i)).self_neighbors_info;
      v_pixel_super_flower.push_back(pix_neighb_inf);
      for( unsigned int j = 0; j < pixel_vec.at(super_flower_pixel_v.at(i)).v_pixel_flower.size(); j++){
	pixel_neighbors_info pix_neighb_inf = pixel_vec.at(super_flower_pixel_v.at(i)).v_pixel_flower.at(j);
	v_pixel_super_flower.push_back(pix_neighb_inf);
      }
    }
  }
  Float_t get_dist_pixel(Float_t px, Float_t py){
    return TMath::Sqrt((px - x)*(px - x) + (py - y)*(py - y));
  }
  static void bubbleSort(std::vector<pixel_neighbors_info> &a){
    bool swapp = true;
    while(swapp){
      swapp = false;
      for( unsigned int i = 0; i < a.size()-1; i++){
	if (a.at(i).phi_deg>a.at(i+1).phi_deg){
	  //
	  pixel_neighbors_info tmp;
	  tmp = a.at(i);
	  //
	  a.at(i) = a.at(i+1);
	  a.at(i+1) = tmp;
	  swapp = true;
	}
      }
    }
  }
  static Float_t find_min_dist(const std::vector<pixel_neighbors_info> a){
    Float_t min_dist = 999.0;
    for( unsigned int i = 0; i < a.size(); i++)
      if(a.at(i).distance_between_pixels<min_dist)
	min_dist = a.at(i).distance_between_pixels;
    return min_dist;
  }
  //void get_flower_pixel_v(std::vector<unsigned int> &super_flower_pixel_v){
  //get_flower_pixel_v(super_flower_pixel_v, true);
  //}
  void get_flower_pixel_v(std::vector<unsigned int> &super_flower_pixel_v, Bool_t first_yes = true){
    //super_flower_pixel_v.push_back(13);
    //super_flower_pixel_v.push_back(15);
    //super_flower_pixel_v.push_back(23);
    //super_flower_pixel_v.push_back(31);
    //super_flower_pixel_v.push_back(39);
    //super_flower_pixel_v.push_back(47);
    Float_t min_dist = find_min_dist(v_pixel_neighbors_third);
    Float_t phi_0_deg = -999.9;
    Int_t counter = 0;
    for( unsigned int i = 0; i < v_pixel_neighbors_third.size(); i++){
      if(v_pixel_neighbors_third.at(i).distance_between_pixels>=min_dist*0.99 &&
	 v_pixel_neighbors_third.at(i).distance_between_pixels<=min_dist*1.01){
	if(first_yes){
	  if(counter == 0)
	    phi_0_deg = v_pixel_neighbors_third.at(i).phi_deg;
	}
	else{
	  if(counter == 1)
	    phi_0_deg = v_pixel_neighbors_third.at(i).phi_deg;
	}
	counter++;
      }
    }
    counter = 0;
    for( unsigned int i = 0; i < v_pixel_neighbors_third.size(); i++){
      Float_t phi_tmp = v_pixel_neighbors_third.at(i).phi_deg - phi_0_deg;
      if((phi_tmp>= 60.0*counter-1.0) && (phi_tmp<= 60.0*counter+1.0)){
	super_flower_pixel_v.push_back(v_pixel_neighbors_third.at(i).pixel_id);
	counter++;
      }
    }
  }
  //static void print_info(const std::vector<line_info::line_flower_info> &a){
  //line_info::print_info_header();
  //for( unsigned int i = 0; i < a.size(); i++)
  //line_info::print_info(a.at(i))
  //}
  //static void bubbleSort(std::vector<line_info::line_flower_info> &a){
  void get_flower_contour_lines(double pixel_pitch){
    //
    Double_t alpha   = 2.0*TMath::Pi()/6.0;
    Double_t alpha_2 = alpha/2.0;
    Double_t alpha0  = alpha_2;
    Double_t pixel_pitch_2 = pixel_pitch/2.0;
    Double_t r = pixel_pitch_2/TMath::Cos(alpha_2);
    Double_t theta = 0.0;
    std::vector<Double_t> xv_v;
    std::vector<Double_t> yv_v;
    for(Int_t i = 0;i<n;i++){
      theta = alpha0 + alpha*i;
      xv_v.push_back(r*TMath::Cos(theta));
      yv_v.push_back(r*TMath::Sin(theta));
    }
    std::vector<Double_t> xc_v;
    std::vector<Double_t> yc_v;
    for(Int_t i = 0;i<n;i++){
      theta = alpha*i;
      xc_v.push_back(pixel_pitch*TMath::Cos(theta));
      yc_v.push_back(pixel_pitch*TMath::Sin(theta));
    }
    //
    TLine ln01(x+xc_v.at(0)+xv_v.at(0),
	       y+yc_v.at(0)+yv_v.at(0),
	       x+xc_v.at(0)+xv_v.at(1),
	       y+yc_v.at(0)+yv_v.at(1));
    v_line_flower.push_back(ln01);
    TLine ln02(x+xc_v.at(0)+xv_v.at(1),
	       y+yc_v.at(0)+yv_v.at(1),
	       x+xc_v.at(1)+xv_v.at(0),
	       y+yc_v.at(1)+yv_v.at(0));
    v_line_flower.push_back(ln02);
    TLine ln03(x+xc_v.at(1)+xv_v.at(0),
	       y+yc_v.at(1)+yv_v.at(0),
	       x+xc_v.at(1)+xv_v.at(1),
	       y+yc_v.at(1)+yv_v.at(1));
    v_line_flower.push_back(ln03);
    TLine ln04(x+xc_v.at(1)+xv_v.at(1),
	       y+yc_v.at(1)+yv_v.at(1),
	       x+xc_v.at(1)+xv_v.at(2),
	       y+yc_v.at(1)+yv_v.at(2));
    v_line_flower.push_back(ln04);
    TLine ln05(x+xc_v.at(1)+xv_v.at(2),
	       y+yc_v.at(1)+yv_v.at(2),
	       x+xc_v.at(2)+xv_v.at(1),
	       y+yc_v.at(2)+yv_v.at(1));
    v_line_flower.push_back(ln05);
    TLine ln06(x+xc_v.at(2)+xv_v.at(1),
	       y+yc_v.at(2)+yv_v.at(1),
	       x+xc_v.at(2)+xv_v.at(2),
	       y+yc_v.at(2)+yv_v.at(2));
    v_line_flower.push_back(ln06);
    TLine ln07(x+xc_v.at(2)+xv_v.at(2),
	       y+yc_v.at(2)+yv_v.at(2),
	       x+xc_v.at(2)+xv_v.at(3),
	       y+yc_v.at(2)+yv_v.at(3));
    v_line_flower.push_back(ln07);
    TLine ln08(x+xc_v.at(2)+xv_v.at(3),
	       y+yc_v.at(2)+yv_v.at(3),
	       x+xc_v.at(3)+xv_v.at(2),
	       y+yc_v.at(3)+yv_v.at(2));
    v_line_flower.push_back(ln08);
    TLine ln09(x+xc_v.at(3)+xv_v.at(2),
	       y+yc_v.at(3)+yv_v.at(2),
	       x+xc_v.at(3)+xv_v.at(3),
	       y+yc_v.at(3)+yv_v.at(3));
    v_line_flower.push_back(ln09);
    TLine ln10(x+xc_v.at(3)+xv_v.at(3),
	       y+yc_v.at(3)+yv_v.at(3),
	       x+xc_v.at(3)+xv_v.at(4),
	       y+yc_v.at(3)+yv_v.at(4));
    v_line_flower.push_back(ln10);
    TLine ln11(x+xc_v.at(3)+xv_v.at(4),
	       y+yc_v.at(3)+yv_v.at(4),
	       x+xc_v.at(4)+xv_v.at(3),
	       y+yc_v.at(4)+yv_v.at(3));
    v_line_flower.push_back(ln11);
    TLine ln12(x+xc_v.at(4)+xv_v.at(3),
	       y+yc_v.at(4)+yv_v.at(3),
	       x+xc_v.at(4)+xv_v.at(4),
	       y+yc_v.at(4)+yv_v.at(4));
    v_line_flower.push_back(ln12);
    TLine ln13(x+xc_v.at(4)+xv_v.at(4),
	       y+yc_v.at(4)+yv_v.at(4),
	       x+xc_v.at(4)+xv_v.at(5),
	       y+yc_v.at(4)+yv_v.at(5));
    v_line_flower.push_back(ln13);
    TLine ln14(x+xc_v.at(4)+xv_v.at(5),
	       y+yc_v.at(4)+yv_v.at(5),
	       x+xc_v.at(5)+xv_v.at(4),
	       y+yc_v.at(5)+yv_v.at(4));
    v_line_flower.push_back(ln14);
    TLine ln15(x+xc_v.at(5)+xv_v.at(4),
	       y+yc_v.at(5)+yv_v.at(4),
	       x+xc_v.at(5)+xv_v.at(5),
	       y+yc_v.at(5)+yv_v.at(5));
    v_line_flower.push_back(ln15);
    TLine ln16(x+xc_v.at(5)+xv_v.at(5),
	       y+yc_v.at(5)+yv_v.at(5),
	       x+xc_v.at(5)+xv_v.at(6),
	       y+yc_v.at(5)+yv_v.at(6));
    v_line_flower.push_back(ln16);
    TLine ln17(x+xc_v.at(5)+xv_v.at(6),
    	       y+yc_v.at(5)+yv_v.at(6),
    	       x+xc_v.at(6)+xv_v.at(5),
    	       y+yc_v.at(6)+yv_v.at(5));
    v_line_flower.push_back(ln17);
    TLine ln18(x+xc_v.at(6)+xv_v.at(5),
    	       y+yc_v.at(6)+yv_v.at(5),
    	       x+xc_v.at(6)+xv_v.at(6),
    	       y+yc_v.at(6)+yv_v.at(6));
    v_line_flower.push_back(ln18);
    //
    //Double_t rp_tmp;
    //std::vector<line_info::line_flower_info> v_line_flower_points_info;
    //for( unsigned int i = 0; i < pixel_vec.at(0).v_pixel_flower.size(); i++){
    //for(Int_t j = 0;j<n;j++){
    //rp_tmp = TMath::Sqrt(pixel_vec.at(0).v_pixel_flower.at(i).xp[j]*pixel_vec.at(0).v_pixel_flower.at(i).xp[j] +
    //pixel_vec.at(0).v_pixel_flower.at(i).yp[j]*pixel_vec.at(0).v_pixel_flower.at(i).yp[j]);
    //if(rp_tmp>pixel_pitch){
    //line_info::line_flower_info line_str;
    //line_str.x = pixel_vec.at(0).v_pixel_flower.at(i).xp[i]*pixel_vec.at(0).v_pixel_flower.at(i).xp[j];
    //line_str.y = pixel_vec.at(0).v_pixel_flower.at(i).xp[i]*pixel_vec.at(0).v_pixel_flower.at(i).yp[j];
    //TVector2 v(line_str.x,line_str.y);
    //line_str.phi = v.Phi();
    //v_line_flower_points_info.push_back(line_str);
    //}
    //}
    //}
    //std::vector<TLine> v_line_flower;    
    //assert(0);
  }
};

class sipmCameraHist: public TH2Poly {

public:
  
  sipmCameraHist(const char* name, const char* title, sipmCameraHist *sipmHist);
  sipmCameraHist(const char* name, const char* title, const char* mapping_csv_file, Double_t rot_alpha_deg);
  sipmCameraHist(const char* name, const char* title, const char* mapping_csv_file, Double_t rot_alpha_deg, TH1D *h1_distance_between_pixels);
  ~sipmCameraHist();
  void dump_mapping_info();
  void test(TString pdf_out_name = "sipmCameraHist_test.pdf");
  void test02();
  void test03();
  void test04();
  void test05();
  void test055();
  void save_pixel_neighbors_to_csv(TString outfilename = "pixel_mapping_neighbors.csv", Int_t npix_neighbors = 6);
  void save_trigger_channel_mask_isolated_flower(TString file_out_name = "trigger_channel_mask_isolated_flower.list");
  void save_trigger_channel_mask_all_pixels(TString file_out_name="trigger_channel_mask_all_pixels.list");
  void save_isolated_flower_seed_flower(TString file_out_name = "isolated_flower_seed_flower.list");
  void save_all_seed_flower(TString file_out_name = "all_seed_flower.list");
  void save_isolated_flower_seed_super_flower(TString file_out_name = "isolated_flower_seed_super_flower.list");
  void test_drawer_id();
  void test_pixel_neighbors_id();
  void test_pixel_neighbors_id(Int_t pix_id);
  void test_pixel_neighbors_id(Int_t npixels_n, Int_t *pix_id);
  void test_pixel_neighbors_second_id();
  void test_pixel_neighbors_second_id(Int_t pix_id);
  void test_pixel_neighbors_second_id(Int_t npixels_n, Int_t *pix_id);
  void test_pixel_neighbors_third_id();
  void test_pixel_neighbors_third_id(Int_t pix_id);
  void test_pixel_neighbors_third_id(Int_t npixels_n, Int_t *pix_id);
  void test_pixel_super_flower();
  void test_pixel_super_flower(Int_t pix_id);
  void test_pixel_super_flower(Int_t npixels_n, Int_t *pix_id);
  void test_pixel_neighbors_bubbleSort(Int_t pix_id);
  void test_trigger_channel_mask_isolated_flower(TString pdf_out_name = "test_trigger_channel_mask_isolated_flower.pdf");
  void test_trigger_channel_mask_isolated_flower_plus_super_flower(TString pdf_out_name = "test_trigger_channel_mask_isolated_flower_plus_super_flower.pdf", unsigned int seedID = 0);
  void test_of_inefficient_regions_isolated_flower_plus_super_flower(TString file_out_name = "test_of_inefficient_regions_isolated_flower_plus_super_flower.pdf");
  void Clean();
  void count_signal(Double_t th_val, Int_t &nch, Int_t &npe);
  //
  std::vector<Int_t> get_trigger_channel_mask_isolated_flower();
  //
  TCanvas *Draw_cam_pixID();
  void Draw_cam(TString settings, TString pdf_out_file);
  void Draw_cam(TString settings, TString pdf_out_file, sipmCameraHist *simp_ref_hist, const anabase *ab = NULL);
  void Draw_cam(TString settings, TString pdf_out_file, sipmCameraHist *simp_ref_hist, const std::vector<unsigned int> &pixel_line_flower_vec, const anabase *ab = NULL);
  void Draw_cam(TString settings, TString pdf_out_file, const std::vector<unsigned int> &pixel_line_flower_vec);
  void Draw_cam(TString settings,
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
  void Draw_cam(TString settings,
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
		Int_t n_pixels,
		const std::vector<unsigned int> &pixel_line_flower_vec,
		sipmCameraHist *simp_ref_hist);
  void Draw_cam(TString settings,
		TString pdf_out_file,
		TString particle_type,
		Int_t wf_time_id,
		Int_t event_id,
		Float_t energy,
		Float_t xcore,
		Float_t ycore,
		Float_t ev_time,
		Int_t nphotons,
		Int_t n_pe,
		Int_t n_pixels,
		const std::vector<unsigned int> &pixel_line_flower_vec);
  //
  inline void set_wf_time_id(Int_t wf_time_id){_wf_time_id = wf_time_id;}
  inline void set_anabase(const anabase *ab){_ab = ab;}
  //  
  TString _name;
  TString _title;
  //
  static const unsigned int _n_pixels = 7987;
  static const unsigned int _n_drawers = 1141;
  //
  inline const std::vector<pixel_info> &get_pixel_vec() const {return _pixel_vec;}
  void Fill_wf(const std::vector<std::vector<Int_t>> &wf);
  void Fill_wf(const std::vector<Int_t> &wf);
  void Fill_pe(const Int_t npixels_n, const Int_t *pix_id);
  void Fill_pe(const Int_t npixels_n, const Int_t *pix_id, const Double_t alpha, const Double_t x_shift, const Double_t y_shift);
  void Fill_pe(const Int_t npixels_n, const Int_t *pix_id, const Double_t alpha, TH1D *h1_theta = NULL, TH1D *h1_theta_deg = NULL, TH1D *h1_r = NULL);
  void Fill_pe_center(const Int_t npixels_n, const Int_t *pix_id);
  void Fill_pix_x_y_hist( const Int_t npixels_n, const Int_t *pix_id, TH1D *h1_x, TH1D *h1_y);
  void Fill_pix_hist2D_y_vs_x( const Int_t npixels_n, const Int_t *pix_id, TH2D *h2_y_vs_x);
  void get_pix_mean( const Int_t npixels_n, const Int_t *pix_id, Double_t &x_mean, Double_t &y_mean);
  //
  void get_pix_density_info(const Int_t npixels_n, const Int_t *pix_id,
			    Double_t &x_mean, Double_t &y_mean,
			    Double_t &x_min, Double_t &x_max,
			    Double_t &y_min, Double_t &y_max,
			    Double_t &dx, Double_t &dy,
			    Double_t &x_std, Double_t &y_std, Int_t verbosity = 0);
  //
  void get_pix_time_info(const Int_t npixels_n, const Float_t *pe_time,
			 const Float_t ev_time,
			 const Float_t time_offset,
			 Double_t &t_min, Double_t &t_max,
			 Double_t &t_mean, Double_t &t_std,
			 Int_t &dt, Int_t verbosity = 0);
  //
  void simulateFlover_ideal_resp(Int_t pixelID, Int_t npixels_n,
				 Int_t *pix_id, Float_t *pe_time);
  //  
  static void rotatePix(Double_t alpha, const Double_t xo, const Double_t yo, Double_t &xn, Double_t &yn);
  const bool check_ch_ID(const Int_t chIDval) const;
  const bool check_ch_ID(const unsigned int chIDval) const;
  const bool check_ch_ID() const;
  
  static void get_part_coordinates_in_tel_frame( const TVector3 &vx_tel, const TVector3 &vy_tel, const TVector3 &vz_tel,
						 Double_t theta, Double_t phi, Double_t &theta_in_tel, Double_t &phi_in_tel);
  static void get_tel_frame( Double_t tel_theta, Double_t tel_phi, TVector3 &vx_tel, TVector3 &vy_tel, TVector3 &vz_tel);
  static void get_x_y_shift(const TVector3 &vx_tel, const TVector3 &vy_tel, const TVector3 &vz_tel,
			    Double_t azimuth, Double_t altitude, Double_t &x_shift, Double_t &y_shift, Double_t phi0_shift = 90.0);
  static Double_t angle_between_optical_axis_and_particle( Double_t tel_theta, Double_t tel_phi, Double_t azimuth, Double_t altitude);
  static Double_t get_theta_p_t_anaFast(Double_t azimuth, Double_t altitude);
  
  private:

  double _pixel_size;
  double _pixel_pitch;
  std::vector<pixel_info> _pixel_vec;
  void load_mapping(const char* mapping_csv_file);
  Double_t _rot_alpha_deg;
  Int_t _wf_time_id;
  
  const anabase *_ab;
};
