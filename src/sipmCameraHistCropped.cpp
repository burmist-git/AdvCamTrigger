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
#include <TPrincipal.h>

using namespace std;

sipmCameraHistCropped::sipmCameraHistCropped( const char* name, const char* title, const sipmCameraHist *sipmHist, const std::vector<unsigned int> &pixel_map_v) : TH2Poly()
{
  //
  SetName(name);
  SetTitle(title);
  //
  _name = name;
  _title = title;
  _sipm_cam = sipmHist;
  //
  _n_pixels = pixel_map_v.size();
  //
  for( unsigned int i = 0; i < pixel_map_v.size(); i++)
    AddBin(sipmHist->get_pixel_vec().at(pixel_map_v.at(i)).n,
	   sipmHist->get_pixel_vec().at(pixel_map_v.at(i)).xp,
	   sipmHist->get_pixel_vec().at(pixel_map_v.at(i)).yp);
}

sipmCameraHistCropped::sipmCameraHistCropped(const char* name, const char* title, const sipmCameraHist *sipmHist, TString in_map_file_name) : TH2Poly()
{
  //
  SetName(name);
  SetTitle(title);
  //
  _name = name;
  _title = title;
  _sipm_cam = sipmHist;
  //
  _n_pixels = 0;
  load_pixel_map( _pixel_map, in_map_file_name);
  _n_pixels = _pixel_map.size();
  //
  for( unsigned int i = 0; i < _pixel_map.size(); i++)
    AddBin(sipmHist->get_pixel_vec().at(_pixel_map.at(i)).n,
	   sipmHist->get_pixel_vec().at(_pixel_map.at(i)).xp,
	   sipmHist->get_pixel_vec().at(_pixel_map.at(i)).yp);
}

sipmCameraHistCropped::sipmCameraHistCropped(const char* name, const char* title, const sipmCameraHist *sipmHist, bool if_centrate) : TH2Poly()
{
  //
  SetName(name);
  SetTitle(title);
  //
  _name = name;
  _title = title;
  _sipm_cam = sipmHist;
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
  build_pixel_seed(pixel_seed_v);
  //pixel_seed_v.push_back(5439);
  //pixel_seed_v.push_back(5467);
  //
  ////pixel_seed_v.push_back(5483);
  ////pixel_seed_v.push_back(5473);
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
      if(if_centrate){
	xp[j] = sipmHist->get_pixel_vec().at(_pixel_map.at(i)).xp[j] - x_mean;
	yp[j] = sipmHist->get_pixel_vec().at(_pixel_map.at(i)).yp[j] - y_mean;
      }
      else {
	xp[j] = sipmHist->get_pixel_vec().at(_pixel_map.at(i)).xp[j];
	yp[j] = sipmHist->get_pixel_vec().at(_pixel_map.at(i)).yp[j];
      }
    }
    AddBin(sipmHist->get_pixel_vec().at(_pixel_map.at(i)).n,xp,yp);
  }
  _n_pixels = _pixel_map.size();
  //cout<<"_n_pixels = "<<_n_pixels<<endl;
  save_pixel_map(_pixel_map);
}

void sipmCameraHistCropped::save_pixel_map( const vector<unsigned int> &pixel_map_v, TString out_map_file_name){
  ofstream fileout;
  fileout.open(out_map_file_name.Data());
  for(unsigned int i = 0;i<pixel_map_v.size();i++)
    fileout<<pixel_map_v.at(i)<<endl;
  fileout.close();
}

void sipmCameraHistCropped::load_pixel_map( vector<unsigned int> &pixel_map_v, TString in_map_file_name){
  ifstream filein(in_map_file_name.Data());
  unsigned int pixIDval;
  if(filein.is_open()){
    while(filein>>pixIDval)
      pixel_map_v.push_back(pixIDval);
  }
  else{
    cout<<" ERROR --> Unable to open file: "<<in_map_file_name<<endl; 
    assert(0);
  }
  filein.close();
}

void sipmCameraHistCropped::build_pixel_seed(vector<unsigned int> &pixel_seed_v){
  //
  Double_t phi0 = 4.900;
  Double_t phi_max = phi0 + 0.175;
  Double_t phi_min = phi0 - 0.175;
  //
  pixel_seed_v.push_back(0);
  for(unsigned int i = 1;i<_sipm_cam->get_pixel_vec().size();i++){
    if(_sipm_cam->get_pixel_vec().at(i).pix_phi>phi_min &&
       _sipm_cam->get_pixel_vec().at(i).pix_phi<phi_max){
      pixel_seed_v.push_back(i);
    }
  }
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
  SetMarkerSize(0.1);
  Draw_cam("ZCOLOR TEXT",pdf_out_name.Data());
}

void sipmCameraHistCropped::Draw_cam(TString settings,
				     TString pdf_out_file){
  //
  Double_t lx_camera = 2.5;
  Double_t ly_camera = 2.5;
  //Double_t lx_camera = 0.25;
  //Double_t ly_camera = 0.25;
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

void sipmCameraHistCropped::Fill_pe(const Int_t npixels_n, const Int_t *pix_id, const Float_t *pix_pe_time,
				    const Double_t ev_time, const Double_t time_offset, const Double_t alpha,
				    TPrincipal *principal, bool if_centrate,
				    const Double_t x_shift, const Double_t y_shift){
  
  Double_t *xn = new Double_t[npixels_n];
  Double_t *yn = new Double_t[npixels_n];
  //
  Double_t x_mean = 0.0;
  Double_t y_mean = 0.0;
  Int_t n_pe_ok = 0;
  //
  for(Int_t i = 0;i<npixels_n;i++){
    if(_sipm_cam->check_ch_ID(pix_id[i])){
      _sipm_cam->rotatePix(alpha,
			   ((Double_t)_sipm_cam->get_pixel_vec().at((unsigned int)pix_id[i]).x + x_shift),
			   ((Double_t)_sipm_cam->get_pixel_vec().at((unsigned int)pix_id[i]).y + y_shift),
			   xn[i], yn[i]);
      x_mean += xn[i];
      y_mean += yn[i];
      n_pe_ok++;
    }
  }
  if(n_pe_ok == 0)
    return;
  x_mean /= n_pe_ok;
  y_mean /= n_pe_ok;
  //
  Clean();
  //
  for(Int_t i = 0;i<npixels_n;i++){
    if(_sipm_cam->check_ch_ID(pix_id[i])){
      if(if_centrate){
	Fill((xn[i] - x_mean),(yn[i] - y_mean));
      }
      else{
	Fill(xn[i], yn[i]);
	//for(Int_t ii = 0; ii < 100; ii++){
	//Fill( rnd_x, rnd_y);
	//}
      }
    }
  }
  //
  if(principal != NULL){
    Double_t *data = new Double_t[_n_pixels];
    for(unsigned int i = 0;i<_n_pixels;i++)
      data[i]=GetBinContent(i+1);
    principal->AddRow(data);
  }
  //for(Int_t i = 0;i<npixels_n;i++)
  //if(_sipm_cam->check_ch_ID())
  //Fill((xn[i] - x_mean),(yn[i] - y_mean));
}

void sipmCameraHistCropped::Fill_principal( std::vector<sipmCameraHistCropped*> &simp_hist_crop_eigenVectors_v, const TPrincipal *principal){
  //
  if(principal == NULL)
    return;
  const double *data_v = principal->GetEigenVectors()->GetMatrixArray();
  TString name;
  //
  for(Int_t i = 0;i<principal->GetEigenVectors()->GetNcols();i++){
    name = "principal_eigenVectors";
    name += "_v";
    name += i;
    sipmCameraHistCropped* simp_hist_eigenVectors = new sipmCameraHistCropped(name.Data(),name.Data(), _sipm_cam, get_pixel_map());
    for(Int_t j = 0;j<principal->GetEigenVectors()->GetNrows();j++)
      simp_hist_eigenVectors->SetBinContent(j+1,data_v[i*+principal->GetEigenVectors()->GetNrows()+j]);
    simp_hist_crop_eigenVectors_v.push_back(simp_hist_eigenVectors);
  }
}

void sipmCameraHistCropped::Fill_principal( std::vector<sipmCameraHistCropped*> &sipm_cam_principal_hist_v, Double_t data_Vh[_dd_im][_dd_im]){
  TString name;
  //
  for(Int_t i = 0;i<_dd_im;i++){
    name = "principal";
    name += "_v";
    name += i;
    sipmCameraHistCropped* simp_hist_principal = new sipmCameraHistCropped(name.Data(),name.Data(), _sipm_cam, get_pixel_map());
    for(Int_t j = 0;j<_dd_im;j++)
      simp_hist_principal->SetBinContent(j+1,data_Vh[i][j]);
    sipm_cam_principal_hist_v.push_back(simp_hist_principal);
  }  
}

void sipmCameraHistCropped::Fill_reco( std::vector<sipmCameraHistCropped*> &sipm_cam_reco_hist_v, std::vector<std::vector<Double_t>> &reco_v){
  TString name;
  //
  for(unsigned int i = 0;i<reco_v.size();i++){
    name = "reco_";
    name += i;
    sipmCameraHistCropped* simp_hist_reco = new sipmCameraHistCropped(name.Data(),name.Data(), _sipm_cam, get_pixel_map());
    for(Int_t j = 0;j<_dd_im;j++)
      simp_hist_reco->SetBinContent(j+1,reco_v.at(i).at((unsigned int)j));
    sipm_cam_reco_hist_v.push_back(simp_hist_reco);
  }  
}

void sipmCameraHistCropped::Save_to_csv(TString csvname, const std::vector<sipmCameraHistCropped*> simp_hist_crop_v){
  std::ofstream csvfile;
  //
  cout<<"simp_hist_crop_v.size() = "<<simp_hist_crop_v.size()<<endl
      <<"_n_pixels               = "<<_n_pixels<<endl;
  //
  csvfile.open(csvname.Data(), std::ios_base::app);
  for(unsigned int i = 0; i <simp_hist_crop_v.size() ; i++){
    for(unsigned int j = 0; j < _n_pixels; j++){
      if(j == (_n_pixels-1))
	csvfile<<simp_hist_crop_v.at(i)->GetBinContent(j+1);
      else
	csvfile<<simp_hist_crop_v.at(i)->GetBinContent(j+1)<<" ";
    }
    csvfile<<endl;
  }
  csvfile.close();  
}

void sipmCameraHistCropped::draw_crop_vector( Int_t nx, Int_t ny, const std::vector<sipmCameraHistCropped*> simp_hist_crop_v, TCanvas *c1, Int_t faceID_shift){
  c1->cd();
  //gStyle->SetPalette(kGreyScale);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE);
  c1->Divide(nx,ny,0.0001,0.0001,0);
  Int_t padID = 1;
  Int_t faceID = faceID_shift;
  for( Int_t i = 0; i<nx; i++){
    for( Int_t j = 0; j<ny; j++){
      c1->cd(padID);
      gStyle->SetOptStat(kFALSE);
      gPad->SetLeftMargin(0.0);
      gPad->SetTopMargin(0.0);
      gPad->SetBottomMargin(0.0);
      gPad->SetRightMargin(0.0);
      if(faceID<(Int_t)simp_hist_crop_v.size()){
	simp_hist_crop_v.at(faceID)->SetStats(0);
	simp_hist_crop_v.at(faceID)->SetTitle("");
	simp_hist_crop_v.at(faceID)->GetXaxis()->SetLabelSize(0);
	simp_hist_crop_v.at(faceID)->GetYaxis()->SetLabelSize(0);
	simp_hist_crop_v.at(faceID)->GetXaxis()->SetLabelOffset(999);
	simp_hist_crop_v.at(faceID)->GetYaxis()->SetLabelOffset(999);
	simp_hist_crop_v.at(faceID)->GetXaxis()->SetTickLength(0);
	simp_hist_crop_v.at(faceID)->GetYaxis()->SetTickLength(0);
	simp_hist_crop_v.at(faceID)->GetZaxis()->SetTickLength(0);
	simp_hist_crop_v.at(faceID)->Draw("ZCOLOR");
      }
      padID++;
      faceID++;
    }
  }
}
