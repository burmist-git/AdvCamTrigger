//my
#include "lstMirrorHist.hh"

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
#include <TGraphErrors.h>

using namespace std;

lstMirrorHist::lstMirrorHist(Int_t dummy) : TH2Poly()
{
  //
  setHistNameTitle();
  //
  mirror_info mirr;
  mirr.test();
  mirr.build_Cell();
  _mirror_vec.push_back(mirr);
  //
  for(unsigned int i = 0;i<_mirror_vec.size();i++)
    AddBin(_mirror_vec.at(0).n,_mirror_vec.at(i).xp,_mirror_vec.at(i).yp);
  //
}

lstMirrorHist::lstMirrorHist(const char* name, const char* title, lstMirrorHist *mirrHist) : TH2Poly()
{
  //
  setHistNameTitle(name,title);
  //
  for(unsigned int i = 0;i<mirrHist->get_mirror_vec().size();i++)
    AddBin(mirrHist->get_mirror_vec().at(0).n,
	   mirrHist->get_mirror_vec().at(i).xp,
	   mirrHist->get_mirror_vec().at(i).yp);
}

lstMirrorHist::lstMirrorHist(const char* name, const char* title, const char* mapping_csv_file, bool if_mapping_ideal, Double_t rot_alpha_deg) : TH2Poly()
{
  _rot_alpha_deg = rot_alpha_deg;
  _name = name;
  _title = title;
  //
  if(!if_mapping_ideal)
    load_mapping(mapping_csv_file);
  else
    load_mapping_ideal(mapping_csv_file);
  //
  SetName(name);
  SetTitle(title);
  for(unsigned int i = 0;i<_mirror_vec.size();i++)
    AddBin(_mirror_vec.at(0).n,_mirror_vec.at(i).xp,_mirror_vec.at(i).yp);
}

//../LST/lst-sim-config/mirror_CTA-N-LST1_v2019-03-31.dat
//mirror_CTA-N-LST1_v2019-03-31_format.dat
//# Mirror positions and focal lengths for the LST1 prototype 
//# on La Palma, as based on the file newMirrorList.txt provided by LST team,
//# reformated to required format by 
//#   awk '{ printf "  %8.2f\t  %8.2f\t  %6.2f\t  %8.2f\t    %d\t   0.0\t    #%%  id=%d\n",$1,$2,$3,100.*$4,$5,$6 }'
//# 
//# Updated on 31/03/19. Two mirrors (IDs 162 and 172) were exchanged for ones with similar focal lengths.
//#
//# Size and area covered as reported by 'draw_mirrors':
//# 198 mirrors are inside a radius of 12.393 m
//# Total surface area: 390.98 m**2
//# Total projected area: 388.38 m**2
//# Center of gravity at x=0.269 m, y=-0.000 m
//#
//# Columns:
//#   1: x pos. [cm] (note: x -> down)
//#   2: y pos. [cm] (note: y -> right when looking from camera)
//#   3: diameter [cm], flat-to-flat
//#   4: focal length [cm] or 0. 
//#      Note: positive, no RANDOM_FOCAL_LENGTH based values added in sim_telarray.
//#   5: shape code (3)
//#   6: z pos. = 0.0 -> will be calculated automatically based on paraboloid
//#   #% separator keeping further data as comments.
//#   7: mirror number in original list (not used)
//#
//  -1022.49	   -462.00	  151.00	   2920.00	    3	   0.0	    #%  198
//  -1022.49	   -308.00	  151.00	   2910.00	    3	   0.0	    #%  197
//
//
void lstMirrorHist::load_mapping(const char* mapping_csv_file){
  //
  ifstream fFile(mapping_csv_file);
  cout<<mapping_csv_file<<std::endl;
  //
  Int_t mirror_id;
  //
  Float_t x_down;
  Float_t y_right_cam;
  //
  //Float_t x;
  //Float_t y;
  Float_t flat_to_flat;
  Float_t focal_length;
  //
  TString mot;
  //
  if(fFile.is_open()){

    while(fFile>>mot)
      if(mot=="7:")
	break;
    fFile>>mot; //mirror
    fFile>>mot; //number
    fFile>>mot; //in
    fFile>>mot; //original
    fFile>>mot; //list
    fFile>>mot; //(not
    fFile>>mot; //used)
    if(mot != "used)")
      assert(0);
    fFile>>mot;
    while(fFile>>x_down>>y_right_cam>>flat_to_flat>>focal_length>>mot>>mot>>mot>>mirror_id){
      mirror_info mirror_i;
      //
      mirror_i.mirror_id = mirror_id;
      mirror_i.x_down = x_down;
      mirror_i.y_right_cam = y_right_cam;
      mirror_i.flat_to_flat = flat_to_flat;
      mirror_i.focal_length = focal_length;
      mirror_i.get_mirror_x0_y0_from_sim_telarray_cfg();
      TVector2 vv(mirror_i.x,mirror_i.y);
      mirror_i.phi = vv.Phi();
      mirror_i.r = vv.Mod();
      //
      mirror_i.build_Cell();
      //
      _mirror_vec.push_back(mirror_i);
    }
    fFile.close();
  }
  //
}

//../LST/lst-sim-config/mirror_CTA-LST-v20141201-198.dat
//# File provided differs in 4 positions from earlier mirror_CTA-LST-1.51-198-0.03.dat as used in prod-2.
//# The earlier file was derived from the spreadsheat provided by the CTA Project Office. The same mirrors could be selected with:
//#   Parameters for plot_mirrors: -h -p --xy --dslen 56.000000 -s -0.460000,0.000000 1.5 22.6 0 0 28 1.51 1.54 12.4 12.3 12.3
//#   Total projected area: 386.9 m**2 for f = 28.0 m
//#   Total mirror surface area: 391.0 m**2
//#   Rmax = 12.12 m (Dmax = 24.25 m for symmetric dish)
//#   Average distance of mirrors to focus: 28.58 m (parabolic, R=56.00 m)
//#   Using 198 hexagonal mirrors inside a radius of 12.12 m (mirror centers inside 11.26 m)
//#
//# Columns:
//#   1: x pos. [cm] (note: x->down)
//#   2: y pos. [cm]
//#   3: diameter [cm], flat-to-flat
//#   4: focal length [cm] or 0
//#   5: shape code (3)
//#   #% separator to avoid having the given z position to be used by sim_telarray
//#   6: z pos. + 92.550000 [cm] (not used if preceding '#' is preserved)
//#   7: mirror number in original list (not used)
//#   8: n_x (not used)
//#   9: n_y (not used)
//#  10: n_z (not used)
//#  11: distance of optical axis to mirror center [cm] (not used)
//#
//-1022.49         -462.00        151.00  0       3       #%      19.86   127     0.00811         -0.01795        0.09804 1122.01795
//-1022.49         -308.00        151.00  0       3       #%      9.27    123     0.00542         -0.01799        0.09822 1067.86904
//-1022.49         -154.00        151.00  0       3       #%      2.92    134     0.00271         -0.01801        0.09833 1034.01948
//
//
void lstMirrorHist::load_mapping_ideal(const char* mapping_csv_file){
  //
  ifstream fFile(mapping_csv_file);
  cout<<mapping_csv_file<<std::endl;
  //
  Int_t mirror_id;
  //
  Float_t x_down;
  Float_t y_right_cam;
  //
  //Float_t x;
  //Float_t y;
  Float_t z;
  Float_t flat_to_flat;
  Float_t focal_length;
  Float_t n_x;
  Float_t n_y;
  Float_t n_z;
  Float_t dist_of_optical_axis_to_mirror_center_cm;
  //
  TString mot;
  //
  if(fFile.is_open()){
    while(fFile>>mot)
      if(mot=="11:")
	break;
    fFile>>mot; //distance
    fFile>>mot; //of
    fFile>>mot; //optical
    fFile>>mot; //axis
    fFile>>mot; //to
    fFile>>mot; //mirror
    fFile>>mot; //center
    fFile>>mot; //[cm]
    fFile>>mot; //(not
    fFile>>mot; //used)
    if(mot != "used)")
      assert(0);
    fFile>>mot;
    while(fFile>>x_down>>y_right_cam>>flat_to_flat>>focal_length>>mot>>mot>>z>>mirror_id>>n_x>>n_y>>n_z>>dist_of_optical_axis_to_mirror_center_cm){
      mirror_info mirror_i;
      //
      mirror_i.mirror_id = mirror_id;
      mirror_i.x_down = x_down;
      mirror_i.y_right_cam = y_right_cam;
      mirror_i.z = (z + 92.550000);
      mirror_i.flat_to_flat = flat_to_flat;
      mirror_i.focal_length = focal_length;
      mirror_i.n_x = n_y;
      mirror_i.n_y = n_x;
      mirror_i.n_z = n_z;
      mirror_i.dist_of_optical_axis_to_mirror_center_cm = dist_of_optical_axis_to_mirror_center_cm;
      mirror_i.get_mirror_x0_y0_from_sim_telarray_cfg();
      TVector2 vv(mirror_i.x,mirror_i.y);
      mirror_i.phi = vv.Phi();
      mirror_i.r = vv.Mod();
      //
      mirror_i.build_Cell();
      //
      _mirror_vec.push_back(mirror_i);
    }
    fFile.close();
  }
  //
}

void lstMirrorHist::dump_mapping_info(){
  mirror_info::print_info_header();
  for(unsigned int i = 0;i<_mirror_vec.size();i++)
    _mirror_vec.at(i).print_info();
}
 
lstMirrorHist::~lstMirrorHist(){
}

void lstMirrorHist::Clean(){
  for(Int_t i = 0;i<GetNcells();i++)
    SetBinContent(i,0);
}

void lstMirrorHist::test(TString pdf_out_name, TString hist_out_name){
  TFile* rootFile = new TFile(hist_out_name.Data(), "RECREATE", " Histograms", 1);
  rootFile->cd();
  if (rootFile->IsZombie()){
    cout<<"  ERROR ---> file "<<hist_out_name.Data()<<" is zombi"<<endl;
    assert(0);
  }
  else
    cout<<"  Output Histos file ---> "<<hist_out_name.Data()<<endl;
  //
  TString pdf_out_name_full;
  //
  for(Int_t i = 0;i<GetNcells();i++)
    SetBinContent(i,i);
  pdf_out_name_full = pdf_out_name;
  pdf_out_name_full += "_local_numbering.pdf";
  setHistNameTitle("local_numbering","local_numbering");
  Draw_Mirrors("TEXT", pdf_out_name_full)->Write();
  //
  Clean();
  pdf_out_name_full = pdf_out_name;
  pdf_out_name_full += "_plane.pdf";
  setHistNameTitle("plane","plane");
  Draw_Mirrors("TEXT", pdf_out_name_full)->Write();  
  //
  Clean();
  pdf_out_name_full = pdf_out_name;
  pdf_out_name_full += "_lst_numbering.pdf";
  setHistNameTitle("lst_numbering","lst_numbering");
  for(unsigned int i = 0;i<_mirror_vec.size();i++)
    SetBinContent(i+1,_mirror_vec.at(i).mirror_id);
  Draw_Mirrors("TEXT", pdf_out_name_full)->Write();  
  //
  Clean();
  pdf_out_name_full = pdf_out_name;
  pdf_out_name_full += "_lst_focal_length.pdf";
  setHistNameTitle("lst_focal_length","lst_focal_length");
  for(unsigned int i = 0;i<_mirror_vec.size();i++)
      SetBinContent(i+1,_mirror_vec.at(i).focal_length);
  Draw_Mirrors("ZCOLOR", pdf_out_name_full)->Write();  
  //
  rootFile->Close();
}

void lstMirrorHist::test_ideal(TString pdf_out_name, TString hist_out_name){
  TFile* rootFile = new TFile(hist_out_name.Data(), "RECREATE", " Histograms", 1);
  rootFile->cd();
  if (rootFile->IsZombie()){
    cout<<"  ERROR ---> file "<<hist_out_name.Data()<<" is zombi"<<endl;
    assert(0);
  }
  else
    cout<<"  Output Histos file ---> "<<hist_out_name.Data()<<endl;
  //
  TString pdf_out_name_full;
  //
  fill_mirror_vec_with_ideal_val(2800.0);
  //
  Clean();
  pdf_out_name_full = pdf_out_name;
  pdf_out_name_full += "_plane_ideal.pdf";
  setHistNameTitle("plane","plane");
  Draw_Mirrors("TEXT", pdf_out_name_full)->Write();  
  //
  Clean();
  pdf_out_name_full = pdf_out_name;
  pdf_out_name_full += "_lst_numbering_ideal.pdf";
  setHistNameTitle("lst_numbering","lst_numbering");
  for(unsigned int i = 0;i<_mirror_vec.size();i++)
    SetBinContent(i+1,_mirror_vec.at(i).mirror_id);
  Draw_Mirrors("TEXT", pdf_out_name_full)->Write();  
  //
  Clean();
  pdf_out_name_full = pdf_out_name;
  pdf_out_name_full += "_lst_z_ideal.pdf";
  setHistNameTitle("lst_z","lst_z");
  TGraphErrors *gr_z_vs_r = new TGraphErrors();
  gr_z_vs_r->SetNameTitle("gr_z_vs_r","gr_z_vs_r");
  for(unsigned int i = 0;i<_mirror_vec.size();i++){
    SetBinContent(i+1,_mirror_vec.at(i).z);
    gr_z_vs_r->SetPoint(i,TMath::Sqrt(_mirror_vec.at(i).x*_mirror_vec.at(i).x + _mirror_vec.at(i).y*_mirror_vec.at(i).y),
			_mirror_vec.at(i).z);
    gr_z_vs_r->SetPointError(i, 40, 4 + TMath::Abs(_mirror_vec.at(i).z)*0.03);
  }
  Draw_Mirrors("ZCOLOR", pdf_out_name_full)->Write();
  gr_z_vs_r->Write();
  //
  Clean();
  pdf_out_name_full = pdf_out_name;
  pdf_out_name_full += "_lst_dist_oa_ideal.pdf";
  setHistNameTitle("lst_dist_oa","lst_dist_oa");
  for(unsigned int i = 0;i<_mirror_vec.size();i++)
    SetBinContent(i+1,_mirror_vec.at(i).dist_of_optical_axis_to_mirror_center_cm);
  Draw_Mirrors("ZCOLOR", pdf_out_name_full)->Write();
  //
  Clean();
  pdf_out_name_full = pdf_out_name;
  pdf_out_name_full += "_lst_nz_ideal.pdf";
  setHistNameTitle("lst_nz","lst_nz");
  //
  //TH1D *h1_nx = new TH1D("h1_nx","h1_nx",1000, -0.03, 0.03);
  //TH1D *h1_ny = new TH1D("h1_ny","h1_ny",1000, -0.03, 0.03);
  //TH1D *h1_nz = new TH1D("h1_nz","h1_nz",1000, 0.095, 0.105);
  //TH1D *h1_nmod = new TH1D("h1_nmod","h1_nmod",1000, 0.0, 0.2);
  //
  TH1D *h1_nx = new TH1D("h1_nx","h1_nx",1000, -1.1, 1.1);
  TH1D *h1_ny = new TH1D("h1_ny","h1_ny",1000, -1.1, 1.1);
  TH1D *h1_nz = new TH1D("h1_nz","h1_nz",1000, -1.1, 1.1);
  TH1D *h1_nmod = new TH1D("h1_nmod","h1_nmod",1000, 0.0, 1.1);
  //
  for(unsigned int i = 0;i<_mirror_vec.size();i++){
    SetBinContent(i+1,_mirror_vec.at(i).n_z);
    ///////////////////////
    h1_nx->Fill(_mirror_vec.at(i).n_x);
    h1_ny->Fill(_mirror_vec.at(i).n_y);
    h1_nz->Fill(_mirror_vec.at(i).n_z);
    h1_nmod->Fill(TMath::Sqrt(_mirror_vec.at(i).n_x*_mirror_vec.at(i).n_x +
			      _mirror_vec.at(i).n_y*_mirror_vec.at(i).n_y +
			      _mirror_vec.at(i).n_z*_mirror_vec.at(i).n_z));
    ///////////////////////
  }
  SetMinimum(0.980);
  SetMaximum(1.020);
  Draw_Mirrors("ZCOLOR", pdf_out_name_full)->Write();
  h1_nx->Write();
  h1_ny->Write();
  h1_nz->Write();
  h1_nmod->Write();
  //
  Clean();
  pdf_out_name_full = pdf_out_name;
  pdf_out_name_full += "_lst_nx_ideal.pdf";
  setHistNameTitle("lst_nx","lst_nx");
  for(unsigned int i = 0;i<_mirror_vec.size();i++)
    SetBinContent(i+1,_mirror_vec.at(i).n_x);
  SetMinimum(-0.2);
  SetMaximum( 0.2);
  Draw_Mirrors("ZCOLOR", pdf_out_name_full)->Write();  
  //
  Clean();
  pdf_out_name_full = pdf_out_name;
  pdf_out_name_full += "_lst_ny_ideal.pdf";
  setHistNameTitle("lst_ny","lst_ny");
  for(unsigned int i = 0;i<_mirror_vec.size();i++)
    SetBinContent(i+1,_mirror_vec.at(i).n_y);
  SetMinimum(-0.2);
  SetMaximum( 0.2);
  Draw_Mirrors("ZCOLOR", pdf_out_name_full)->Write();  
  //
  Clean();
  pdf_out_name_full = pdf_out_name;
  pdf_out_name_full += "_lst_spherical_focus.pdf";
  setHistNameTitle("lst_spherical_focus_ny","lst_spherical_focus");
  for(unsigned int i = 0;i<_mirror_vec.size();i++){
    Double_t focus = 1.0/(1.0/2800.0 + 1.0/12.0/100000.0);
    Double_t curvature =  get_ideal_curvature( focus,  _mirror_vec.at(i).x,  _mirror_vec.at(i).y);
    SetBinContent(i+1,get_spherical_mirror_focus_form_curvature(curvature));
  }
  SetMinimum(2790);
  SetMaximum(2990);
  Draw_Mirrors("ZCOLOR", pdf_out_name_full)->Write();  
  //
  TH1D *h1_angle_deg_mirror_oa = new TH1D("h1_angle_deg_mirror_oa","h1_angle_deg_mirror_oa",1000,0.0,45.0);
  TH1D *h1_angle_deg_mirror_oa_delta = new TH1D("h1_angle_deg_mirror_oa_delta","h1_angle_deg_mirror_oa_delta",1000,-0.1,0.1);
  for(unsigned int i = 0;i<_mirror_vec.size();i++){
    Double_t focus = 1.0/(1.0/2800.0);
    h1_angle_deg_mirror_oa->Fill(get_ideal_angle_deg_mirror_oa( focus, _mirror_vec.at(i).x,  _mirror_vec.at(i).y));
    Double_t focus01 = 1.0/(1.0/2800.0);
    Double_t focus02 = 1.0/(1.0/2800.0 + 1.0/12.0/100000.0);
    h1_angle_deg_mirror_oa_delta->Fill(get_ideal_angle_deg_mirror_oa( focus01, _mirror_vec.at(i).x,  _mirror_vec.at(i).y) -
				       get_ideal_angle_deg_mirror_oa( focus02, _mirror_vec.at(i).x,  _mirror_vec.at(i).y));

  }
  h1_angle_deg_mirror_oa->Write();
  h1_angle_deg_mirror_oa_delta->Write();
  //
  save_to_csv_mirror_vec();
  //
  rootFile->Close();
}

TCanvas *lstMirrorHist::Draw_Mirrors(TString settings, TString pdf_out_file){
  //
  //Double_t frame_lx = 200;//cm
  //Double_t frame_ly = 200;//cm
  //
  Double_t frame_lx = 3000;//cm
  Double_t frame_ly = 3000;//cm
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
  //SetTitle("");
  //SetName("");
  //
  TString c1_name = "c1_";
  c1_name += _name;
  TString c1_title = "c1_";
  c1_title += _title;
  TCanvas *c1 = new TCanvas(c1_name.Data(),c1_title.Data(),600,600);
  //
  gPad->SetRightMargin(0.12);
  gPad->SetLeftMargin(0.12);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.15);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  //gPad->SetLogz();
  //
  //SetMaximum(500.0);
  //SetMinimum(300.0);
  //SetMinimum(280.0);
  //SetMinimum(299.0);
  //SetMaximum(308.0);
  //
  TH2F *frame = new TH2F( "h2", "h2",
			  10, -frame_lx/2.0, frame_lx/2.0,
			  10, -frame_ly/2.0, frame_ly/2.0);
  //frame->SetTitle("");
  frame->SetTitle(_title.Data());
  frame->GetXaxis()->SetTitle("x, cm");
  frame->GetYaxis()->SetTitle("y, cm");
  frame->GetXaxis()->CenterTitle();
  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetTitleOffset(1.5);
  frame->SetStats(kFALSE);
  frame->Draw();
  //
  //settings += " same TEXT";
  settings += " same";
  //
  //if(simp_ref_hist != NULL)
  //simp_ref_hist->Draw(settings.Data());
  //else
  Draw(settings.Data());
  c1->SaveAs(pdf_out_file.Data());
  //
  return c1;
}

void lstMirrorHist::setHistNameTitle(TString namestr, TString titlestr){
  SetName(namestr.Data());
  SetTitle(titlestr.Data());
  _name = namestr.Data();
  _title = titlestr.Data();
}

void lstMirrorHist::save_to_csv_mirror_vec(TString out_file_name){
  std::ofstream outfile;
  outfile.open(out_file_name.Data());
  outfile<<"mirror_id,x_down,y_right_cam,x,y,z,flat_to_flat,focal_length,phi,r,n_x,n_y,n_z,dist_oa"<<endl;
  for(unsigned int i = 0;i<_mirror_vec.size();i++){
    outfile<<_mirror_vec.at(i).mirror_id<<","
	   <<_mirror_vec.at(i).x_down<<","
	   <<_mirror_vec.at(i).y_right_cam<<","
	   <<_mirror_vec.at(i).x<<","
	   <<_mirror_vec.at(i).y<<","
	   <<_mirror_vec.at(i).z<<","
	   <<_mirror_vec.at(i).flat_to_flat<<","
	   <<_mirror_vec.at(i).focal_length<<","
	   <<_mirror_vec.at(i).phi<<","
	   <<_mirror_vec.at(i).r<<","
	   <<_mirror_vec.at(i).n_x<<","
	   <<_mirror_vec.at(i).n_y<<","
	   <<_mirror_vec.at(i).n_z<<","
	   <<_mirror_vec.at(i).dist_of_optical_axis_to_mirror_center_cm<<"\n";
  }
  outfile.close();
}

void lstMirrorHist::fill_mirror_vec_with_ideal_val(Double_t focus_val_cm){
  Double_t n_x, n_y, n_z;
  for(unsigned int i = 0;i<_mirror_vec.size();i++){
    _mirror_vec.at(i).z = get_ideal_z( focus_val_cm, _mirror_vec.at(i).x, _mirror_vec.at(i).y);    
    get_ideal_nx_ny_nz( focus_val_cm,
			_mirror_vec.at(i).x, _mirror_vec.at(i).y, _mirror_vec.at(i).z, n_x, n_y, n_z);
    _mirror_vec.at(i).n_x = (Float_t)n_x;
    _mirror_vec.at(i).n_y = (Float_t)n_y;
    _mirror_vec.at(i).n_z = (Float_t)n_z;
    //<<_mirror_vec.at(i).focal_length<<","
    //<<_mirror_vec.at(i).dist_of_optical_axis_to_mirror_center_cm<<"\n";
  }
}

Double_t lstMirrorHist::get_ideal_z( Double_t focus_val_cm, Double_t x, Double_t y){
  Double_t k = 1.0/4.0/focus_val_cm;
  return k*(x*x + y*y);
}

void lstMirrorHist::get_ideal_nx_ny_nz( Double_t focus_val_cm,
					Double_t x, Double_t y, Double_t z,
					Double_t &nx, Double_t &ny, Double_t &nz){
  Double_t dd = 0.1;   // cm
  TVector3 r(x,y,z);
  TVector3 rdx( x+dd, y   , get_ideal_z( focus_val_cm, x+dd, y   ));
  TVector3 rdy( x   , y+dd, get_ideal_z( focus_val_cm, x   , y+dd));
  TVector3 drx = rdx - r;
  TVector3 dry = rdy - r;
  TVector3 n = drx.Cross(dry);
  n.SetMag(1.0);
  nx = n.X();
  ny = n.Y();
  nz = n.Z();
}

Double_t lstMirrorHist::get_ideal_angle_deg_mirror_oa( Double_t focus_val_cm, Double_t x, Double_t y){
  Double_t z = get_ideal_z( focus_val_cm, x, y);
  Double_t nx, ny, nz;
  get_ideal_nx_ny_nz( focus_val_cm, x,  y,  z, nx, ny, nz);
  TVector3 nm(nx, ny, nz);
  TVector3 oa(0.0, 0.0, 1.0);  
  return TMath::ACos(oa.Dot(nm))*180.0/TMath::Pi();
}

Double_t lstMirrorHist::get_ideal_curvature( Double_t focus_val_cm, Double_t x, Double_t y){
  Double_t k = 1.0/4.0/focus_val_cm;
  Double_t r = TMath::Sqrt(x*x+y*y);
  return 2*k/(TMath::Power((1.0 + 4*k*k*r*r),3.0/2.0));
}

Double_t lstMirrorHist::get_spherical_mirror_focus_form_curvature( Double_t curvature){
  if(curvature>0.0)
    return 1.0/curvature/2.0;  
  return -999.0;
}
