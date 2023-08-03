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
#include <TColor.h>

using namespace std;

sipmCameraHist::sipmCameraHist(const char* name, const char* title, const char* mapping_csv_file, Double_t rot_alpha_deg, TH1D *h1_distance_between_pixels) : TH2Poly(), _rot_alpha_deg(rot_alpha_deg)
{
  //
  //
  _pixel_size = 0.02329999953508377;
  _pixel_pitch = 0.024300;
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
  //
  //for(unsigned int i = 0;i<1;i++)
  for(unsigned int i = 0;i<_pixel_vec.size();i++) 
    _pixel_vec.at(i).find_pixel_neighbors(_pixel_vec,_pixel_pitch, h1_distance_between_pixels);
  for(unsigned int i = 0;i<_pixel_vec.size();i++)
    _pixel_vec.at(i).build_pixel_super_flower(_pixel_vec);
  //
  for(unsigned int i = 0;i<_pixel_vec.size();i++)
    _pixel_vec.at(i).get_flower_contour_lines(_pixel_pitch);
}

sipmCameraHist::sipmCameraHist(const char* name, const char* title, const char* mapping_csv_file, Double_t rot_alpha_deg) : sipmCameraHist(name, title, mapping_csv_file, rot_alpha_deg, NULL)
{
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
      TVector2 vv(pix_i.x,pix_i.y);
      pix_i.pix_phi = vv.Phi();
      pix_i.pix_r = vv.Mod();
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

void sipmCameraHist::Draw_cam( TString settings,
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
			       Int_t n_pixels){
  std::vector<unsigned int> pixel_line_flower_vec;
  Draw_cam(settings, pdf_out_file, particle_type, wf_time_id, event_id, energy, xcore, ycore, ev_time, nphotons, n_pe, n_pixels, pixel_line_flower_vec);
}

void sipmCameraHist::Draw_cam( TString settings,
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
			       const std::vector<unsigned int> &pixel_line_flower_vec){
  Draw_cam( settings, pdf_out_file,
	    particle_type, wf_time_id, event_id,
	    energy, xcore, ycore, ev_time, nphotons,
	    n_pe, n_pixels, pixel_line_flower_vec, NULL);
}

void sipmCameraHist::Draw_cam( TString settings,
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
			       const std::vector<unsigned int> &pixel_line_flower_vec,
			       sipmCameraHist *simp_ref_hist = NULL){
  //
  Double_t lx_camera = 2.5;
  Double_t ly_camera = 2.5;
  //Double_t lx_camera = 0.5;
  //Double_t ly_camera = 0.5;
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
  TCanvas *c1;
  if(simp_ref_hist != NULL){
    c1 = new TCanvas("c1","c1",1800,600);
    c1->Divide(3,1);
    //std::cout<<"simp_ref_hist != NULL"<<std::endl;
  }
  else{
    c1 = new TCanvas("c1","c1",1400,700);
    //std::cout<<"simp_ref_hist == NULL"<<std::endl;
    c1->Divide(2,1);
  }
  c1->cd(1);
  gPad->SetRightMargin(0.12);
  gPad->SetLeftMargin(0.12);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.15);
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
  if(simp_ref_hist != NULL)
    simp_ref_hist->Draw(settings.Data());
  else
    Draw(settings.Data());
  //
  //cout<<"pixel_line_flower_vec.size() --> "<<pixel_line_flower_vec.size()<<endl;
  for( unsigned int i = 0; i < pixel_line_flower_vec.size(); i++){
    //cout<<_pixel_vec.at(pixel_line_flower_vec.at(i)).v_line_flower.size()<<endl;
    for( unsigned int j = 0; j < _pixel_vec.at(pixel_line_flower_vec.at(i)).v_line_flower.size(); j++){
      //cout<<_pixel_vec.at(pixel_line_flower_vec.at(i)).v_line_flower.size()<<endl;
      _pixel_vec.at(pixel_line_flower_vec.at(i)).v_line_flower.at(j).SetLineColor(kRed);
      _pixel_vec.at(pixel_line_flower_vec.at(i)).v_line_flower.at(j).SetLineWidth(1.0);
      _pixel_vec.at(pixel_line_flower_vec.at(i)).v_line_flower.at(j).Draw("same");
    }
  }
  //
  if(simp_ref_hist != NULL){
    c1->cd(2);
    frame->Draw();
    Draw(settings.Data());
    for( unsigned int i = 0; i < pixel_line_flower_vec.size(); i++){
      for( unsigned int j = 0; j < _pixel_vec.at(pixel_line_flower_vec.at(i)).v_line_flower.size(); j++){
	//cout<<_pixel_vec.at(pixel_line_flower_vec.at(i)).v_line_flower.size()<<endl;
	_pixel_vec.at(pixel_line_flower_vec.at(i)).v_line_flower.at(j).SetLineColor(kRed);
	_pixel_vec.at(pixel_line_flower_vec.at(i)).v_line_flower.at(j).SetLineWidth(1.0);
	//_pixel_vec.at(pixel_line_flower_vec.at(i)).v_line_flower.at(j).Draw("same");
      }
    }
  }
  //
  if(simp_ref_hist != NULL)
    c1->cd(3);
  else
    c1->cd(2);
  //
  TString wf_time_str  = "wf_time  : "; wf_time_str += wf_time_id; wf_time_str += " ns";
  TString event_id_str = "event_id : "; event_id_str += event_id;
  TString energy_str   = "energy   : "; energy_str += (Int_t)(energy*1000.0); energy_str += " GeV";
  TString xcore_str    = "xcore    : "; xcore_str += (Int_t)xcore; xcore_str += " m";
  TString ycore_str    = "ycore    : "; ycore_str += (Int_t)ycore; ycore_str += " m";
  TString ev_time_str  = "ev_time  : "; ev_time_str += (Int_t)ev_time; ev_time_str += " ns";
  TString nphotons_str = "nphotons : "; nphotons_str += nphotons;
  TString n_pe_str     = "n_pe     : "; n_pe_str += n_pe;
  TString n_pixels_str = "n_pixels : "; n_pixels_str += n_pixels;
  //
  TText *t0 = new TText(0.5,0.95,wf_time_str.Data());
  t0->SetTextAlign(22);
  t0->SetTextFont(43);
  t0->SetTextSize(40);
  t0->Draw();
  TText *t = new TText(0.5,0.9,particle_type.Data());
  t->SetTextAlign(22);
  t->SetTextFont(43);
  t->SetTextSize(40);
  t->Draw();
  TText *t2 = new TText(0.5,0.85,event_id_str.Data());
  t2->SetTextAlign(22);
  t2->SetTextFont(43);
  t2->SetTextSize(40);
  t2->Draw("same");
  TText *t3 = new TText(0.5,0.80,energy_str.Data());
  t3->SetTextAlign(22);
  t3->SetTextFont(43);
  t3->SetTextSize(40);
  t3->Draw("same");
  TText *t4 = new TText(0.5,0.75,xcore_str.Data());
  t4->SetTextAlign(22);
  t4->SetTextFont(43);
  t4->SetTextSize(40);
  t4->Draw("same");
  TText *t5 = new TText(0.5,0.70,ycore_str.Data());
  t5->SetTextAlign(22);
  t5->SetTextFont(43);
  t5->SetTextSize(40);
  t5->Draw("same");
  TText *t6 = new TText(0.5,0.65,ev_time_str.Data());
  t6->SetTextAlign(22);
  t6->SetTextFont(43);
  t6->SetTextSize(40);
  t6->Draw("same");
  TText *t7 = new TText(0.5,0.60,nphotons_str.Data());
  t7->SetTextAlign(22);
  t7->SetTextFont(43);
  t7->SetTextSize(40);
  t7->Draw("same");
  TText *t8 = new TText(0.5,0.55,n_pe_str.Data());
  t8->SetTextAlign(22);
  t8->SetTextFont(43);
  t8->SetTextSize(40);
  t8->Draw("same");
  TText *t9 = new TText(0.5,0.50,n_pixels_str.Data());
  t9->SetTextAlign(22);
  t9->SetTextFont(43);
  t9->SetTextSize(40);
  t9->Draw("same");
  //
  //
  if(pdf_out_file != "")
    c1->SaveAs(pdf_out_file.Data());
}

void sipmCameraHist::Draw_cam( TString settings,
			       TString pdf_out_file){
  Draw_cam( settings, pdf_out_file, "NONE", -999, -999, -999.0, -999.0, -999.0, -999.0, -999, -999, -999);
}

void sipmCameraHist::Draw_cam( TString settings,
			       TString pdf_out_file,
			       sipmCameraHist *simp_ref_hist){
  std::vector<unsigned int> pixel_line_flower_vec;
  Draw_cam( settings, pdf_out_file, "NONE", -999, -999, -999.0, -999.0, -999.0, -999.0, -999, -999, -999, pixel_line_flower_vec, simp_ref_hist);
}

void sipmCameraHist::Draw_cam(TString settings, TString pdf_out_file, sipmCameraHist *simp_ref_hist, const std::vector<unsigned int> &pixel_line_flower_vec){
  Draw_cam( settings, pdf_out_file, "NONE", -999, -999, -999.0, -999.0, -999.0, -999.0, -999, -999, -999, pixel_line_flower_vec, simp_ref_hist);
}

void sipmCameraHist::Draw_cam( TString settings,
			       TString pdf_out_file,
			       const std::vector<unsigned int> &pixel_line_flower_vec){
  Draw_cam( settings, pdf_out_file, "NONE", -999, -999, -999.0, -999.0, -999.0, -999.0, -999, -999, -999, pixel_line_flower_vec, NULL);
}

void sipmCameraHist::test(){
  TRandom3 *rnd = new TRandom3(123123); 
  //cout<<"GetN() "<<GetNcells()<<endl;
  for(Int_t i = 0;i<GetNcells();i++)
    SetBinContent(i,(Int_t)rnd->Uniform(1,10));
  //
  Draw_cam("ZCOLOR","sipmCameraHist_test.pdf");
}

void sipmCameraHist::test02(){
  for(unsigned int i = 0;i<_pixel_vec.size();i++){
    if(i<98)
      SetBinContent(i+1,i+1);
    else
      SetBinContent(i+1,0);
  }
  Draw_cam("ZCOLOR","sipmCameraHist_test02.pdf");
}

void sipmCameraHist::test03(){
  for(unsigned int i = 0;i<_pixel_vec.size();i++){
      SetBinContent(i+1,0);
  }
  SetMinimum(1.0);
  Draw_cam("","sipmCameraHist_test03.pdf");
}

void sipmCameraHist::test04(){
  for(unsigned int i = 0;i<10;i++){
      SetBinContent(i+1,i+1);
  }
  SetMinimum(1.0);
  Draw_cam("text","sipmCameraHist_test04.pdf");
}

void sipmCameraHist::test05(){
  for(unsigned int i = 0;i<10;i++){
    SetBinContent(i+1,_pixel_vec.at(i).pixel_id);
  }
  SetMinimum(1.0);
  Draw_cam("text","sipmCameraHist_test05.pdf");
}

void sipmCameraHist::test_drawer_id(){
  for(unsigned int i = 0;i<(unsigned int)GetNcells();i++){
    if(i<_pixel_vec.size()){
      //if(_pixel_vec.at(i).drawer_id < 1000)
      SetBinContent(i+1,(_pixel_vec.at(i).drawer_id+1)%10+1);
      //else
      //SetBinContent(i+1,0);
    }
  }
  SetMaximum(20.0);
  Draw_cam("ZCOLOR","sipmCameraHist_test_drawer_id.pdf");
}

void sipmCameraHist::test_pixel_neighbors_id(){
  test_pixel_neighbors_id(0);
}

//pix_id : [0 , _n_pixels)
void sipmCameraHist::test_pixel_neighbors_id(Int_t pix_id){
  if( (pix_id < 0) || ((unsigned int)pix_id > _n_pixels)){
    std::cout<<" ERROR --> (pix_id < 0) || (pix_id > _n_pixels)"<<std::endl
	     <<"                            pix_id = "<<pix_id<<std::endl;
    assert(0);
  }
  //
  SetBinContent(pix_id+1,20);
  for(unsigned int i = 0;i<_pixel_vec.size();i++){
    if(_pixel_vec.at(i).pixel_id == pix_id){
      for(unsigned int j = 0;j<_pixel_vec.at(i).v_pixel_neighbors.size();j++)
	SetBinContent(_pixel_vec.at(i).v_pixel_neighbors.at(j).pixel_id+1,10);
    }
    //
  }
  //    
  SetMaximum(20.0);
  Draw_cam("ZCOLOR","sipmCameraHist_test_pixel_neighbors_id.pdf");
}

void sipmCameraHist::test_pixel_neighbors_second_id(){
  test_pixel_neighbors_second_id(0);
}
  
//pix_id : [0 , _n_pixels)
void sipmCameraHist::test_pixel_neighbors_second_id(Int_t pix_id){
  if( (pix_id < 0) || ((unsigned int)pix_id > _n_pixels)){
    std::cout<<" ERROR --> (pix_id < 0) || (pix_id > _n_pixels)"<<std::endl
	     <<"                            pix_id = "<<pix_id<<std::endl;
    assert(0);
  }
  //
  SetBinContent(pix_id+1,20);
  for(unsigned int i = 0;i<_pixel_vec.size();i++){
    if(_pixel_vec.at(i).pixel_id == pix_id){
      for(unsigned int j = 0;j<_pixel_vec.at(i).v_pixel_neighbors.size();j++)
	SetBinContent(_pixel_vec.at(i).v_pixel_neighbors.at(j).pixel_id+1,10);
      for(unsigned int j = 0;j<_pixel_vec.at(i).v_pixel_neighbors_second.size();j++)
	SetBinContent(_pixel_vec.at(i).v_pixel_neighbors_second.at(j).pixel_id+1,5);
    }
    //
  }
  //    
  SetMaximum(20.0);
  Draw_cam("ZCOLOR","sipmCameraHist_test_pixel_neighbors_second_id.pdf");
}

void sipmCameraHist::test_pixel_neighbors_id(Int_t npixels_n,Int_t *pix_id){
  for(Int_t i = 0;i<npixels_n;i++)
    test_pixel_neighbors_id(pix_id[i]);
}

void sipmCameraHist::test_pixel_neighbors_second_id(Int_t npixels_n,Int_t *pix_id){
  for(Int_t i = 0;i<npixels_n;i++)
    test_pixel_neighbors_second_id(pix_id[i]);
}

void sipmCameraHist::test_pixel_neighbors_third_id(){
  test_pixel_neighbors_third_id(0);
}
  
//pix_id : [0 , _n_pixels)
void sipmCameraHist::test_pixel_neighbors_third_id(Int_t pix_id){
  if( (pix_id < 0) || ((unsigned int)pix_id > _n_pixels)){
    std::cout<<" ERROR --> (pix_id < 0) || (pix_id > _n_pixels)"<<std::endl
	     <<"                            pix_id = "<<pix_id<<std::endl;
    assert(0);
  }
  //
  SetBinContent(pix_id+1,30);
  for(unsigned int i = 0;i<_pixel_vec.size();i++){
    if(_pixel_vec.at(i).pixel_id == pix_id){
      for(unsigned int j = 0;j<_pixel_vec.at(i).v_pixel_neighbors.size();j++)
	SetBinContent(_pixel_vec.at(i).v_pixel_neighbors.at(j).pixel_id+1,20);
      for(unsigned int j = 0;j<_pixel_vec.at(i).v_pixel_neighbors_second.size();j++)
	SetBinContent(_pixel_vec.at(i).v_pixel_neighbors_second.at(j).pixel_id+1,10);
      for(unsigned int j = 0;j<_pixel_vec.at(i).v_pixel_neighbors_third.size();j++)
	SetBinContent(_pixel_vec.at(i).v_pixel_neighbors_third.at(j).pixel_id+1,5);
    }
    //
  }
  //    
  SetMaximum(30.0);
  Draw_cam("ZCOLOR","sipmCameraHist_test_pixel_neighbors_third_id.pdf");
}

void sipmCameraHist::test_pixel_neighbors_third_id(Int_t npixels_n,Int_t *pix_id){
  for(Int_t i = 0;i<npixels_n;i++)
    test_pixel_neighbors_third_id(pix_id[i]);
}

void sipmCameraHist::test_pixel_super_flower(){
  test_pixel_super_flower(0);
}
  
//pix_id : [0 , _n_pixels)
void sipmCameraHist::test_pixel_super_flower(Int_t pix_id){
  if( (pix_id < 0) || ((unsigned int)pix_id > _n_pixels)){
    std::cout<<" ERROR --> (pix_id < 0) || (pix_id > _n_pixels)"<<std::endl
	     <<"                            pix_id = "<<pix_id<<std::endl;
    assert(0);
  }
  //
  SetBinContent(pix_id+1,20);
  for(unsigned int i = 0;i<_pixel_vec.size();i++){
    if(_pixel_vec.at(i).pixel_id == pix_id){
      for(unsigned int j = 0;j<_pixel_vec.at(i).v_pixel_super_flower.size();j++){
	SetBinContent(_pixel_vec.at(i).v_pixel_super_flower.at(j).pixel_id+1,
		      (GetBinContent(_pixel_vec.at(i).v_pixel_super_flower.at(j).pixel_id+1)+10));
      }
    }
  }
  //    
  SetMaximum(0.0);
  SetMaximum(20.0);
  std::vector<unsigned int> pixel_line_flower_vec;
  pixel_line_flower_vec.push_back(pix_id);
  Draw_cam("ZCOLOR","sipmCameraHist_test_pixel_super_flower_id.pdf",pixel_line_flower_vec);
}

void sipmCameraHist::test_pixel_super_flower(Int_t npixels_n,Int_t *pix_id){
  for(Int_t i = 0;i<npixels_n;i++)
    test_pixel_super_flower(pix_id[i]);
}

void sipmCameraHist::test_pixel_neighbors_bubbleSort(Int_t pix_id){
  for(unsigned int i = 0;i<_pixel_vec.size();i++){
    if(_pixel_vec.at(i).pixel_id == pix_id){
      if(_pixel_vec.at(i).v_pixel_neighbors_third.size()>0){
	_pixel_vec.at(i).v_pixel_neighbors_third.at(0).print_info_header();
	for(unsigned int j = 0; j < _pixel_vec.at(i).v_pixel_neighbors_third.size(); j++)
	  _pixel_vec.at(i).v_pixel_neighbors_third.at(j).print_info();
	//
	pixel_info::bubbleSort(_pixel_vec.at(i).v_pixel_neighbors_third);
	//
	_pixel_vec.at(i).v_pixel_neighbors_third.at(0).print_info_header();
	for(unsigned int j = 0; j < _pixel_vec.at(i).v_pixel_neighbors_third.size(); j++)
	  _pixel_vec.at(i).v_pixel_neighbors_third.at(j).print_info();
      }
    }
  }
}

void sipmCameraHist::Fill_wf(const std::vector<std::vector<Int_t>> &wf){
  for( unsigned int i = 0; i < wf.size(); i++ ){
    for( unsigned int j = 0; j < wf.at(i).size(); j++){
      SetBinContent(i+1,GetBinContent(i+1) + wf.at(i).at(j));
    }
  }
}

void sipmCameraHist::Fill_wf(const std::vector<Int_t> &wf){
  for( unsigned int i = 0; i < wf.size(); i++ ){
    SetBinContent(i+1,GetBinContent(i+1) + wf.at(i));
  }
}

void sipmCameraHist::Fill_pe(const Int_t npixels_n, const Int_t *pix_id){
  for(Int_t i = 0;i<npixels_n;i++)
    Fill((Double_t)_pixel_vec.at((unsigned int)pix_id[i]).x,
	 (Double_t)_pixel_vec.at((unsigned int)pix_id[i]).y);
}

void sipmCameraHist::Fill_pe(const Int_t npixels_n, const Int_t *pix_id, const Double_t alpha){
  for(Int_t i = 0;i<npixels_n;i++){
    Double_t xn;
    Double_t yn;
    rotatePix(alpha,
	      (Double_t)_pixel_vec.at((unsigned int)pix_id[i]).x,
	      (Double_t)_pixel_vec.at((unsigned int)pix_id[i]).y,
	      xn, yn);    
    Fill(xn,yn);
  }
}

void sipmCameraHist::get_pix_mean( const Int_t npixels_n, const Int_t *pix_id, Double_t &x_mean, Double_t &y_mean){
  x_mean = 0.0;
  y_mean = 0.0;
  if(npixels_n<=0.0)
    return;
  for(Int_t i = 0;i<npixels_n;i++){
    x_mean += _pixel_vec.at((unsigned int)pix_id[i]).x;
    y_mean += _pixel_vec.at((unsigned int)pix_id[i]).y;
  }
  x_mean /= npixels_n;
  y_mean /= npixels_n;
}

void sipmCameraHist::Fill_pe_center(const Int_t npixels_n, const Int_t *pix_id){
  Double_t x_mean = 0.0;
  Double_t y_mean = 0.0;
  get_pix_mean( npixels_n, pix_id, x_mean, y_mean);
  for(Int_t i = 0;i<npixels_n;i++)
    Fill(((Double_t)_pixel_vec.at((unsigned int)pix_id[i]).x - x_mean),
	 ((Double_t)_pixel_vec.at((unsigned int)pix_id[i]).y - y_mean));
}

void sipmCameraHist::rotatePix(Double_t alpha, const Double_t xo, const Double_t yo, Double_t &xn, Double_t &yn){
  if(alpha != 0.0){
    TVector2 v( xo, yo);
    xn = v.Rotate(alpha).X();
    yn = v.Rotate(alpha).Y();
  }
}
