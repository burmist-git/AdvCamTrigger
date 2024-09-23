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
#include <TVector3.h>
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

sipmCameraHist::sipmCameraHist(const char* name, const char* title, sipmCameraHist *sipmHist) : TH2Poly(), _ab(NULL)
{
  //
  
  //
  SetName(name);
  SetTitle(title);
  //
  _name = name;
  _title = title;
  //
  for(unsigned int i = 0;i<sipmHist->get_pixel_vec().size();i++)
    AddBin(sipmHist->get_pixel_vec().at(0).n,
	   sipmHist->get_pixel_vec().at(i).xp,
	   sipmHist->get_pixel_vec().at(i).yp);
}

sipmCameraHist::sipmCameraHist(const char* name, const char* title, const char* mapping_csv_file, Double_t rot_alpha_deg, TH1D *h1_distance_between_pixels) : TH2Poly(), _rot_alpha_deg(rot_alpha_deg), _ab(NULL)
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

const bool sipmCameraHist::check_ch_ID(const unsigned int chIDval) const {
  if(chIDval>=0 && chIDval<_n_pixels)
    return true;
  return false;
}

const bool sipmCameraHist::check_ch_ID(const Int_t chIDval) const {
  return check_ch_ID((unsigned int )chIDval);
}

const bool sipmCameraHist::check_ch_ID() const {
  return true;
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

std::vector<Int_t> sipmCameraHist::get_trigger_channel_mask_isolated_flower(){
  //
  std::vector<Int_t> isolated_flower_chID;
  std::vector<Int_t> black_list_chID;
  //
  for(unsigned int i = 0;i<_pixel_vec.size();i++){
    //if(_pixel_vec.at(i).pixel_id == 0){
    if(isolated_flower_chID.size() == 0){
      isolated_flower_chID.push_back(_pixel_vec.at(i).pixel_id);
      for(unsigned int kk = 0;kk<_pixel_vec.at(i).v_pixel_flower.size();kk++){
	black_list_chID.push_back(_pixel_vec.at(i).v_pixel_flower.at(kk).pixel_id);
      }
    }
    else{
      bool isOk = true;
      for(unsigned int jj = 0;jj<black_list_chID.size();jj++){
	if(black_list_chID.at(jj) == _pixel_vec.at(i).pixel_id){
	  isOk = false;	  
	  break;
	}
      }
      if(isOk){
	isolated_flower_chID.push_back(_pixel_vec.at(i).pixel_id);
	for(unsigned int kk = 0;kk<_pixel_vec.at(i).v_pixel_flower.size();kk++){
	  black_list_chID.push_back(_pixel_vec.at(i).v_pixel_flower.at(kk).pixel_id);
	}
      }
    }
  }
  //
  return isolated_flower_chID;
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
  TString azimuth_str;
  TString altitude_str;
  TString h_first_int_str;
  TString hmax_str;
  //
  if(_ab != NULL){
    azimuth_str     = "azimuth     : "; azimuth_str += (Int_t)(_ab->azimuth*180.0/TMath::Pi()*10); azimuth_str += "/10 deg";
    //azimuth_str     = "azimuth     : "; azimuth_str += _ab->event_id;
    altitude_str    = "altitude    : "; altitude_str += (Int_t)(_ab->altitude*180.0/TMath::Pi()*10); altitude_str += "/10 deg";
    h_first_int_str = "h_first_int : "; h_first_int_str += (Int_t)_ab->h_first_int; h_first_int_str += " km";
    hmax_str        = "hmax        : "; hmax_str += (Int_t)_ab->hmax; hmax_str += " km";
  }
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
  if(_ab != NULL){
    TText *t10 = new TText(0.5,0.45,azimuth_str.Data());
    t10->SetTextAlign(22);
    t10->SetTextFont(43);
    t10->SetTextSize(40);
    t10->Draw("same");
    TText *t11 = new TText(0.5,0.40,altitude_str.Data());
    t11->SetTextAlign(22);
    t11->SetTextFont(43);
    t11->SetTextSize(40);
    t11->Draw("same");
    TText *t12 = new TText(0.5,0.35,h_first_int_str.Data());
    t12->SetTextAlign(22);
    t12->SetTextFont(43);
    t12->SetTextSize(40);
    t12->Draw("same");
    TText *t13 = new TText(0.5,0.30,hmax_str.Data());
    t13->SetTextAlign(22);
    t13->SetTextFont(43);
    t13->SetTextSize(40);
    t13->Draw("same");
  }
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
			       sipmCameraHist *simp_ref_hist,
			       const anabase *ab){
  std::vector<unsigned int> pixel_line_flower_vec;
  if(ab == NULL)
    Draw_cam( settings, pdf_out_file, "NONE", -999, -999, -999.0, -999.0, -999.0, -999.0, -999, -999, -999, pixel_line_flower_vec, simp_ref_hist);
  else
    Draw_cam( settings, pdf_out_file, ab->_particle_type_name.Data(), _wf_time_id, ab->event_id, ab->energy, ab->xcore, ab->ycore, ab->ev_time, ab->nphotons, ab->n_pe, ab->n_pixels, pixel_line_flower_vec, simp_ref_hist);
}

void sipmCameraHist::Draw_cam(TString settings, TString pdf_out_file, sipmCameraHist *simp_ref_hist, const std::vector<unsigned int> &pixel_line_flower_vec, const anabase *ab){
  if(ab == NULL)
    Draw_cam( settings, pdf_out_file, "NONE", -999, -999, -999.0, -999.0, -999.0, -999.0, -999, -999, -999, pixel_line_flower_vec, simp_ref_hist);
  else
    Draw_cam( settings, pdf_out_file, ab->_particle_type_name.Data(), _wf_time_id, ab->event_id, ab->energy, ab->xcore, ab->ycore, ab->ev_time, ab->nphotons, ab->n_pe, ab->n_pixels, pixel_line_flower_vec, simp_ref_hist);
}

void sipmCameraHist::Draw_cam( TString settings,
			       TString pdf_out_file,
			       const std::vector<unsigned int> &pixel_line_flower_vec){
  Draw_cam( settings, pdf_out_file, "NONE", -999, -999, -999.0, -999.0, -999.0, -999.0, -999, -999, -999, pixel_line_flower_vec, NULL);
}

void sipmCameraHist::save_trigger_channel_mask_isolated_flower(TString file_out_name){
  std::vector<Int_t> isolated_flower_seeds = get_trigger_channel_mask_isolated_flower(); 
  std::cout<<"isolated_flower_seeds.size() = "<<isolated_flower_seeds.size()<<std::endl;  
  //
  ofstream outfile;
  outfile.open(file_out_name.Data());
  //
  for(unsigned int i = 0;i<isolated_flower_seeds.size();i++)
    outfile<<std::setw(10)<<isolated_flower_seeds.at(i)<<std::endl;
  //
  outfile.close();
}

void sipmCameraHist::save_isolated_flower_seed_flower(TString file_out_name){
  std::vector<Int_t> isolated_flower_seeds = get_trigger_channel_mask_isolated_flower(); 
  std::cout<<"save_isolated_flower_seed_flower"<<std::endl;  
  unsigned int pix_seed;
  //
  //
  ofstream outfile;
  outfile.open(file_out_name.Data());
  //
  //
  //std::vector<pixel_neighbors_info> v_pixel_neighbors;
  //std::vector<pixel_neighbors_info> v_pixel_neighbors_second;
  //std::vector<pixel_neighbors_info> v_pixel_neighbors_third;
  //std::vector<pixel_neighbors_info> v_pixel_flower;
  //std::vector<pixel_neighbors_info> v_pixel_super_flower;
  //
  //
  for(unsigned int ii = 0;ii<isolated_flower_seeds.size();ii++){
    pix_seed = isolated_flower_seeds.at(ii);
    outfile<<std::setw(10)<<_pixel_vec.at(pix_seed).v_pixel_neighbors.size()+1;
    outfile<<std::setw(10)<<pix_seed;
    for(unsigned int j = 0;j<_pixel_vec.at(pix_seed).v_pixel_neighbors.size();j++)
      outfile<<std::setw(10)<<_pixel_vec.at(pix_seed).v_pixel_neighbors.at(j).pixel_id;
    outfile<<std::endl;
  }
  //
  outfile.close();
}

void sipmCameraHist::save_isolated_flower_seed_super_flower(TString file_out_name){
  std::vector<Int_t> isolated_flower_seeds = get_trigger_channel_mask_isolated_flower(); 
  std::cout<<"save_isolated_flower_seed_super_flower"<<std::endl;  
  unsigned int pix_seed;
  //
  //
  ofstream outfile;
  outfile.open(file_out_name.Data());
  //
  //
  //std::vector<pixel_neighbors_info> v_pixel_neighbors;
  //std::vector<pixel_neighbors_info> v_pixel_neighbors_second;
  //std::vector<pixel_neighbors_info> v_pixel_neighbors_third;
  //std::vector<pixel_neighbors_info> v_pixel_flower;
  //std::vector<pixel_neighbors_info> v_pixel_super_flower;
  //
  //
  Int_t npointsCounter = 0;
  Int_t npointsCounterMax = 49; // 49 pixels for super flower
  for(unsigned int ii = 0;ii<isolated_flower_seeds.size();ii++){
    pix_seed = isolated_flower_seeds.at(ii);
    outfile<<std::setw(10)<<_pixel_vec.at(pix_seed).v_pixel_super_flower.size()+1;
    outfile<<std::setw(10)<<pix_seed;
    npointsCounter = 1;
    for(unsigned int j = 0;j<_pixel_vec.at(pix_seed).v_pixel_super_flower.size();j++){
      outfile<<std::setw(10)<<_pixel_vec.at(pix_seed).v_pixel_super_flower.at(j).pixel_id;
      npointsCounter++;
    }
    if(npointsCounter == npointsCounterMax){
      outfile<<std::endl;
    }
    else if(npointsCounter<npointsCounterMax){
      for(int ii = npointsCounter; ii<npointsCounterMax;ii++)
	outfile<<std::setw(10)<<-999;
      outfile<<std::endl;
    }
    else if(npointsCounter>npointsCounterMax){
      cout<<"ERROR --> npointsCounter > npointsCounterMax"<<endl
	  <<"          npointsCounter = "<<npointsCounter<<endl
	  <<"       npointsCounterMax = "<<npointsCounterMax<<endl;
      assert(0);
    }
  }
  //
  outfile.close();
}

void sipmCameraHist::save_trigger_channel_mask_all_pixels(TString file_out_name){
  //
  ofstream outfile;
  outfile.open(file_out_name.Data());
  //
  for(unsigned int i = 0;i<_pixel_vec.size();i++)
    outfile<<std::setw(10)<<_pixel_vec.at(i).pixel_id<<std::endl;
  //
  outfile.close();
}

void sipmCameraHist::test_trigger_channel_mask_isolated_flower(TString pdf_out_name){
  std::vector<Int_t> isolated_flower_seeds = get_trigger_channel_mask_isolated_flower(); 
  //
  std::cout<<"isolated_flower_seeds.size() = "<<isolated_flower_seeds.size()<<std::endl;  
  //
  for(unsigned int i = 0;i<isolated_flower_seeds.size();i++)
      SetBinContent(isolated_flower_seeds.at(i)+1,1);
  //
  Draw_cam("ZCOLOR",pdf_out_name.Data());
}

void sipmCameraHist::test_trigger_channel_mask_isolated_flower_plus_super_flower(TString pdf_out_name, unsigned int seedID){
  std::vector<Int_t> isolated_flower_seeds = get_trigger_channel_mask_isolated_flower(); 
  std::cout<<"isolated_flower_seeds.size() = "<<isolated_flower_seeds.size()<<std::endl;  
  //
  for(unsigned int i = 0;i<isolated_flower_seeds.size();i++)
      SetBinContent(isolated_flower_seeds.at(i)+1,10);
  //
  //std::vector<unsigned int> pix_seed_v;
  //pix_seed_v.push_back(0);
  //pix_seed_v.push_back(77);
  //pix_seed_v.push_back(119);
  //
  unsigned int pix_seed = 0;
  //
  //
  //for(unsigned int ii = 0;ii<pix_seed_v.size();ii++){
  //pix_seed = pix_seed_v.at(ii);
  pix_seed=seedID;
  for(unsigned int j = 0;j<_pixel_vec.at(pix_seed).v_pixel_super_flower.size();j++){
    if(GetBinContent(_pixel_vec.at(pix_seed).v_pixel_super_flower.at(j).pixel_id+1)<10)
      SetBinContent(_pixel_vec.at(pix_seed).v_pixel_super_flower.at(j).pixel_id+1,5);
  }
  //}
  //
  //    
  SetMaximum(0.0);
  SetMaximum(10.0);
  Draw_cam("ZCOLOR",pdf_out_name.Data());
}

void sipmCameraHist::test(TString pdf_out_name){
  TRandom3 *rnd = new TRandom3(123123); 
  //cout<<"GetN() "<<GetNcells()<<endl;
  for(Int_t i = 0;i<GetNcells();i++)
    SetBinContent(i,(Int_t)rnd->Uniform(1,10));
  //
  Draw_cam("ZCOLOR",pdf_out_name.Data());
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

void sipmCameraHist::test055(){
  std::cout<<"_pixel_vec.size() = "<<_pixel_vec.size()<<std::endl;
  for(unsigned int i = 0;i<_pixel_vec.size();i++){
    SetBinContent(i+1,_pixel_vec.at(i).pixel_id);
  }
  //Draw_cam("text","sipmCameraHist_test_pix_ID.pdf");
  TCanvas *c1 = Draw_cam_pixID();
  c1->SaveAs("sipmCameraHist_test_pix_ID.pdf");
}


/*
  for(unsigned int i = 0;i<_pixel_vec.size();i++){
    if(_pixel_vec.at(i).pixel_id == pix_id){
      for(unsigned int j = 0;j<_pixel_vec.at(i).v_pixel_neighbors.size();j++)
	SetBinContent(_pixel_vec.at(i).v_pixel_neighbors.at(j).pixel_id+1,10);
      for(unsigned int j = 0;j<_pixel_vec.at(i).v_pixel_neighbors_second.size();j++)
	SetBinContent(_pixel_vec.at(i).v_pixel_neighbors_second.at(j).pixel_id+1,5);
    }
    //
  }
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
*/

void sipmCameraHist::save_pixel_neighbors_to_csv(TString outfilename, Int_t npix_neighbors){
  ofstream outfile;
  outfile.open(outfilename.Data());
  //
  Int_t neighbors_counter = 0;
  //
  for(unsigned int i = 0;i<_pixel_vec.size();i++){
    neighbors_counter = 0;
    for(unsigned int j = 0;j<_pixel_vec.at(i).v_pixel_neighbors.size();j++){
      outfile<<std::setw(10)<<_pixel_vec.at(i).v_pixel_neighbors.at(j).pixel_id;
      neighbors_counter++;
    }
    for(unsigned int i = neighbors_counter ; i < npix_neighbors ; i++){
      outfile<<std::setw(10)<<"NAN";
    }
    outfile<<std::endl;
  }
  outfile.close();
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
  Draw_cam("ZCOLOR","sipmCameraHist_test_pixel_super_flower_id.pdf", pixel_line_flower_vec);
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
    if(pix_id[i]>=0 && pix_id[i]<(Int_t)_n_pixels)
      Fill((Double_t)_pixel_vec.at((unsigned int)pix_id[i]).x,
	   (Double_t)_pixel_vec.at((unsigned int)pix_id[i]).y);
}

Double_t sipmCameraHist::angle_between_optical_axis_and_particle( Double_t tel_theta, Double_t tel_phi,
								  Double_t azimuth, Double_t altitude){
  //
  Double_t part_theta = TMath::Pi()-altitude;
  Double_t part_phi = azimuth;
  //
  TVector3 tel;
  TVector3 part;
  tel.SetMagThetaPhi(1.0,tel_theta,tel_phi);
  part.SetMagThetaPhi(1.0,part_theta,part_phi);
  return TMath::ACos(tel.Dot(part));
}

void sipmCameraHist::get_tel_frame( Double_t tel_theta, Double_t tel_phi, TVector3 &vx_tel, TVector3 &vy_tel, TVector3 &vz_tel){
  vx_tel.RotateZ(tel_phi);
  vx_tel.RotateY(-tel_theta);
  vy_tel.RotateZ(tel_phi);
  vy_tel.RotateY(-tel_theta);
  vz_tel.RotateZ(tel_phi);
  vz_tel.RotateY(-tel_theta);
}

void sipmCameraHist::get_part_coordinates_in_tel_frame( const TVector3 &vx_tel, const TVector3 &vy_tel, const TVector3 &vz_tel,
							Double_t theta, Double_t phi, Double_t &theta_in_tel, Double_t &phi_in_tel){
  TVector3 part;
  part.SetMagThetaPhi(1.0,theta,phi);
  TVector3 v_in_tel(vx_tel.Dot(part),vy_tel.Dot(part),vz_tel.Dot(part));
  theta_in_tel = v_in_tel.Theta();
  phi_in_tel = v_in_tel.Phi();
}

Double_t sipmCameraHist::get_theta_p_t_anaFast(Double_t azimuth, Double_t altitude){
  TVector3 v_det(1.0*TMath::Sin(20.0/180.0*TMath::Pi()),0,1.0*TMath::Cos(20.0/180.0*TMath::Pi()));
  TVector3 v_prot;
  v_prot.SetMagThetaPhi(1.0,TMath::Pi()/2.0-altitude,TMath::Pi() - azimuth);
  TVector3 v_prot_inv(v_prot.x(),v_prot.y(),v_prot.z());
  return TMath::ACos(v_prot_inv.Dot(v_det)/v_prot_inv.Mag()/v_det.Mag());
}

void sipmCameraHist::get_x_y_shift(const TVector3 &vx_tel, const TVector3 &vy_tel, const TVector3 &vz_tel,
				   Double_t azimuth, Double_t altitude, Double_t &x_shift, Double_t &y_shift,
				   Double_t phi0_shift){
  Double_t R_m_LST = 56.4;        //m
  Double_t F_m_LST = R_m_LST/2.0; //m
  Double_t theta = TMath::Pi()/2.0-altitude;
  Double_t phi = azimuth;
  Double_t theta_in_tel;
  Double_t phi_in_tel;
  sipmCameraHist::get_part_coordinates_in_tel_frame(vx_tel, vy_tel, vz_tel,
						    theta, phi,
						    theta_in_tel, phi_in_tel);
  Double_t phi0 = phi0_shift/180.0*TMath::Pi();
  Double_t r_shift = TMath::Tan(theta_in_tel)*F_m_LST;
  TVector2 v_shift;
  v_shift.SetMagPhi( r_shift, phi_in_tel + phi0);
  x_shift = -v_shift.X();
  y_shift = v_shift.Y();
  //cout<<"theta    "<<theta*180.0/TMath::Pi()<<endl
  //<<"phi          "<<phi*180.0/TMath::Pi()<<endl
  //<<"theta_in_tel "<<theta_in_tel*180.0/TMath::Pi()<<endl
  //<<"phi_in_tel   "<<phi_in_tel*180.0/TMath::Pi()<<endl;
}

void sipmCameraHist::Fill_pe(const Int_t npixels_n, const Int_t *pix_id, const Double_t alpha, const Double_t x_shift, const Double_t y_shift){
  for(Int_t i = 0;i<npixels_n;i++){
    if(pix_id[i]>=0 && pix_id[i]<(Int_t)_n_pixels){
      Double_t xn;
      Double_t yn;
      if(alpha != 0.0){
	rotatePix(alpha,
		  ((Double_t)_pixel_vec.at((unsigned int)pix_id[i]).x + x_shift),
		  ((Double_t)_pixel_vec.at((unsigned int)pix_id[i]).y + y_shift),
		  xn, yn);
      }
      else{
	xn = (Double_t)_pixel_vec.at((unsigned int)pix_id[i]).x + x_shift;
	yn = (Double_t)_pixel_vec.at((unsigned int)pix_id[i]).y + y_shift;
      }
      Fill(xn,yn);
    }
  }
}

void sipmCameraHist::Fill_pe(const Int_t npixels_n, const Int_t *pix_id, const Double_t alpha, TH1D *h1_theta, TH1D *h1_theta_deg, TH1D *h1_r){
  for(Int_t i = 0;i<npixels_n;i++){
    if(pix_id[i]>=0 && pix_id[i]<(Int_t)_n_pixels){
      Double_t xn;
      Double_t yn;
      rotatePix(alpha,
		(Double_t)_pixel_vec.at((unsigned int)pix_id[i]).x,
		(Double_t)_pixel_vec.at((unsigned int)pix_id[i]).y,
		xn, yn);    
      Fill(xn,yn);
      if(h1_theta != NULL){
	TVector2 vn( xn, yn);
	h1_theta->Fill(vn.Phi());
	h1_theta_deg->Fill(vn.Phi()*180.0/TMath::Pi());
	h1_r->Fill(vn.Mod());
      }
    } 
  }
}

void sipmCameraHist::Fill_pix_hist2D_y_vs_x( const Int_t npixels_n, const Int_t *pix_id, TH2D *h2_y_vs_x){
  for(Int_t i = 0;i<npixels_n;i++)
    if(pix_id[i]>=0 && pix_id[i]<(Int_t)_n_pixels)
      h2_y_vs_x->Fill((Double_t)_pixel_vec.at((unsigned int)pix_id[i]).x, (Double_t)_pixel_vec.at((unsigned int)pix_id[i]).y);
}

void sipmCameraHist::Fill_pix_x_y_hist( const Int_t npixels_n, const Int_t *pix_id, TH1D *h1_x, TH1D *h1_y){
  for(Int_t i = 0;i<npixels_n;i++){
    if(pix_id[i]>=0 && pix_id[i]<(Int_t)_n_pixels){
      h1_x->Fill((Double_t)_pixel_vec.at((unsigned int)pix_id[i]).x);
      h1_y->Fill((Double_t)_pixel_vec.at((unsigned int)pix_id[i]).y);
    }
  }
}

void sipmCameraHist::Fill_pe_center(const Int_t npixels_n, const Int_t *pix_id){
  Double_t x_mean = 0.0;
  Double_t y_mean = 0.0;
  get_pix_mean( npixels_n, pix_id, x_mean, y_mean);
  for(Int_t i = 0;i<npixels_n;i++)
    if(pix_id[i]>=0 && pix_id[i]<(Int_t)_n_pixels)
      Fill(((Double_t)_pixel_vec.at((unsigned int)pix_id[i]).x - x_mean),
	   ((Double_t)_pixel_vec.at((unsigned int)pix_id[i]).y - y_mean));
}

void sipmCameraHist::get_pix_density_info( const Int_t npixels_n, const Int_t *pix_id,
					   Double_t &x_mean, Double_t &y_mean,
					   Double_t &x_min, Double_t &x_max,
					   Double_t &y_min, Double_t &y_max,
					   Double_t &dx, Double_t &dy,
					   Double_t &x_std, Double_t &y_std, Int_t verbosity){
  x_mean = 0.0;
  y_mean = 0.0;
  x_min = 0.0;
  x_max = 0.0;
  y_min = 0.0;
  y_max = 0.0;
  dx = 0.0;
  dy = 0.0;
  x_std = 0.0;
  y_std = 0.0;
  //
  Double_t x_mean_sq = 0.0;
  Double_t y_mean_sq = 0.0;
  //  
  if(npixels_n<=0)
    return;
  //
  Double_t x_val = _pixel_vec.at((unsigned int)pix_id[0]).x;
  Double_t y_val = _pixel_vec.at((unsigned int)pix_id[0]).y;
  x_min = x_val;
  x_max = x_val;
  y_min = y_val;
  y_max = y_val;
  //  
  for(Int_t i = 0;i<npixels_n;i++){
    //
    x_val = _pixel_vec.at((unsigned int)pix_id[i]).x;
    y_val = _pixel_vec.at((unsigned int)pix_id[i]).y;
    //
    x_mean += x_val;
    y_mean += y_val;
    x_mean_sq += x_val*x_val;
    y_mean_sq += y_val*y_val;
    //
    if(x_min>x_val)
      x_min = x_val;
    //
    if(y_min>y_val)
      y_min = y_val;
    //
    if(x_max<x_val)
      x_max = x_val;
    //
    if(y_max<y_val)
      y_max = y_val;
  }
  //
  x_mean /= npixels_n;
  y_mean /= npixels_n;
  //
  x_mean_sq /= npixels_n;
  y_mean_sq /= npixels_n;
  //
  dx = TMath::Abs((x_max - x_min));
  dy = TMath::Abs((y_max - y_min));
  //
  //for(Int_t i = 0;i<npixels_n;i++){
  //x_val = _pixel_vec.at((unsigned int)pix_id[i]).x;
  //y_val = _pixel_vec.at((unsigned int)pix_id[i]).y;
  //x_std += (x_val-x_mean)*(x_val-x_mean);
  //y_std += (y_val-y_mean)*(y_val-y_mean);
  //}  
  //
  //x_std /= npixels_n;
  //y_std /= npixels_n;
  //
  //x_std = TMath::Sqrt(x_std);
  //y_std = TMath::Sqrt(y_std);
  //  
  x_std = TMath::Sqrt(x_mean_sq - x_mean*x_mean);
  y_std = TMath::Sqrt(y_mean_sq - y_mean*y_mean);
  //
  if(verbosity>0){
    cout<<"x_min  "<<x_min<<endl
	<<"x_max  "<<x_max<<endl
	<<"y_min  "<<y_min<<endl
	<<"y_max  "<<y_max<<endl
	<<"x_mean "<<x_mean<<endl
	<<"y_mean "<<y_mean<<endl
	<<"dx     "<<dx<<endl
	<<"dy     "<<dy<<endl
	<<"x_std  "<<x_std<<endl
	<<"y_std  "<<y_std<<endl;  
  }
}

void sipmCameraHist::get_pix_time_info( const Int_t npixels_n, const Float_t *pe_time,
					const Float_t ev_time,
					const Float_t time_offset,
					Double_t &t_min, Double_t &t_max,
					Double_t &t_mean, Double_t &t_std,
					Int_t &dt, Int_t verbosity){
  t_min = 0.0;
  t_max = 0.0;
  t_mean = 0.0;
  t_std = 0.0;
  dt = 0;
  //  
  if(npixels_n<=0)
    return;
  //
  t_min = pe_time[0] - ev_time + time_offset;
  t_max = pe_time[0] - ev_time + time_offset;
  //
  Double_t t_val;
  Double_t t_mean_sq = 0.0;
  //
  for(Int_t i = 0;i<npixels_n;i++){
    //
    t_val = pe_time[i] - ev_time + time_offset;
    //
    t_mean += t_val;
    t_mean_sq += t_val*t_val;
    //
    if(t_min>t_val)
      t_min = t_val;
    //
    if(t_max<t_val)
      t_max = t_val;
  }
  //
  t_mean /= npixels_n;
  //
  t_mean_sq /= npixels_n;
  //
  dt = floor(TMath::Abs((t_max - t_min))+1);
  //
  t_std = TMath::Sqrt(t_mean_sq - t_mean*t_mean);
  //
  if(verbosity>0){
    cout<<"t_min  "<<t_min<<endl
	<<"t_max  "<<t_max<<endl
	<<"t_mean "<<t_mean<<endl
	<<"t_std  "<<t_std<<endl
	<<"dt     "<<dt<<endl;
  }
}

void sipmCameraHist::get_pix_mean( const Int_t npixels_n, const Int_t *pix_id, Double_t &x_mean, Double_t &y_mean){
  x_mean = 0.0;
  y_mean = 0.0;
  if(npixels_n<=0)
    return;
  for(Int_t i = 0;i<npixels_n;i++){
    x_mean += _pixel_vec.at((unsigned int)pix_id[i]).x;
    y_mean += _pixel_vec.at((unsigned int)pix_id[i]).y;
  }
  x_mean /= npixels_n;
  y_mean /= npixels_n;
}

void sipmCameraHist::rotatePix(Double_t alpha, const Double_t xo, const Double_t yo, Double_t &xn, Double_t &yn){
  if(alpha != 0.0){
    TVector2 v( xo, yo);
    xn = v.Rotate(alpha).X();
    yn = v.Rotate(alpha).Y();
  }
  else{
    xn = xo;
    yn = yo;
  }
}

void sipmCameraHist::simulateFlover_ideal_resp(Int_t pixelID, Int_t npixels_n,
					       Int_t *pix_id, Float_t *pe_time){
  if(npixels_n<=0)
    return;
  Int_t i_counter = 0;
  //cout<<"v_pixel_flower_size : "<<_pixel_vec.at((unsigned int)pixelID).v_pixel_flower.size()<<endl;
  //assert(0);
  unsigned int j_max = (unsigned int)npixels_n/(_pixel_vec.at((unsigned int)pixelID).v_pixel_flower.size()+1);
  for(unsigned int j = 0;j<j_max;j++){
    for(unsigned int i = 0;i<(_pixel_vec.at((unsigned int)pixelID).v_pixel_flower.size()+1);i++){
      if(i == 0){
	pix_id[i_counter] = pixelID;
	pe_time[i_counter] = 30.0;
	i_counter++;
      }
      if(i_counter<npixels_n){
	pix_id[i_counter] = _pixel_vec.at((unsigned int)pixelID).v_pixel_flower.at(i).pixel_id;
	pe_time[i_counter] = 30.0;
	i_counter++;
      }
    }
  }
}

TCanvas *sipmCameraHist::Draw_cam_pixID(){
  //
  Double_t lx_camera = 2.5;
  Double_t ly_camera = 2.5;
  //Double_t lx_camera = 0.5;
  //Double_t ly_camera = 0.5;
  Double_t d_frame = 0.1;
  //
  gStyle->SetPalette(kRainBow);
  //gStyle->SetPalette(kCool);
  //gStyle->SetPalette(kIsland);
  //gStyle->SetPalette(kCherry);
  //TColor::InvertPalette();
  //
  //gStyle->SetPalette(kInvertedDarkBodyRadiator);
  //
  gStyle->SetOptStat(kFALSE);
  SetTitle("");
  SetName("");
  //
  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  gPad->SetRightMargin(0.12);
  gPad->SetLeftMargin(0.12);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.15);
  //
  //gStyle->SetBarWidth(0.05);
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
  /*
  TH2F *frame = new TH2F( "h2", "h2", 40, -lx_camera/2.0-d_frame,lx_camera/2.0+d_frame,40, -ly_camera/2.0-d_frame,ly_camera/2.0+d_frame);
  frame->SetTitle("");
  frame->GetXaxis()->SetTitle("x, m");
  frame->GetYaxis()->SetTitle("y, m");
  frame->GetXaxis()->CenterTitle();
  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetTitleOffset(1.5);
  frame->SetStats(kFALSE);
  frame->Draw();
  */
  //
  SetMinimum(0);
  SetMaximum(8000);
  SetMarkerSize(0.1);
  GetYaxis()->SetTickLength(0);
  //SetLineColorAlpha(kBlack,0.1);
  //SetLineWidth(0.01);
  //
  Draw("TEXT ZCOLOR");
  //if(pdf_out_file != "")
  //c1->SaveAs("sipmCameraHist_test_pix_ID_test.pdf");
  return c1;
}
