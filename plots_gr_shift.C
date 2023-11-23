//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include <time.h>

using namespace std;

void get_x_y_shift(const TVector3 &vx_tel, const TVector3 &vy_tel, const TVector3 &vz_tel,
		   Double_t azimuth, Double_t altitude, Double_t &x_shift, Double_t &y_shift,
		   Double_t phi0_shift);

void get_tel_frame( Double_t tel_theta, Double_t tel_phi, TVector3 &vx_tel, TVector3 &vy_tel, TVector3 &vz_tel);

void get_part_coordinates_in_tel_frame( const TVector3 &vx_tel, const TVector3 &vy_tel, const TVector3 &vz_tel,
					Double_t theta, Double_t phi, Double_t &theta_in_tel, Double_t &phi_in_tel);

Double_t plots_gr_shift_loop(Double_t phi0_shift=77.4);

Int_t plots_gr_shift(){

  /*
  TGraph *gr_l2 = new TGraph();
  
  Double_t phi0_shift;
  Double_t phi0_shift_min = 70.0;
  Double_t phi0_shift_max = 90.0;
  Int_t nn = 1000;
  for(Int_t i = 0; i<nn;i++){
    phi0_shift = phi0_shift_min + (phi0_shift_max - phi0_shift_min)/(nn-1)*i;
    gr_l2->SetPoint(i,phi0_shift,plots_gr_shift_loop(phi0_shift));
  }  
  gr_l2->Draw("APL");
*/

  plots_gr_shift_loop(77.4);
  return 0;
}

Double_t plots_gr_shift_loop(Double_t phi0_shift){

  Double_t tel_theta = 20.0/180.0*TMath::Pi();
  Double_t tel_phi = 180.0/180.0*TMath::Pi();
  
  TVector3 vx_tel(1.0,0.0,0.0);
  TVector3 vy_tel(0.0,1.0,0.0);
  TVector3 vz_tel(0.0,0.0,1.0);
  get_tel_frame( tel_theta, tel_phi, vx_tel, vy_tel, vz_tel);
  
  TGraph *gr_alt_x = new TGraph();
  gr_alt_x->SetNameTitle("gr_alt","gr_alt");
  gr_alt_x->SetPoint(0,68,-0.2);
  gr_alt_x->SetPoint(1,69,-0.1);
  gr_alt_x->SetPoint(2,70,0.0);
  gr_alt_x->SetPoint(3,71,0.1);
  gr_alt_x->SetPoint(4,72,0.2);

  TGraph *gr_alt_y = new TGraph();
  gr_alt_y->SetNameTitle("gr_alt","gr_alt");
  gr_alt_y->SetPoint(0,68, 0.9);
  gr_alt_y->SetPoint(1,69, 0.45);
  gr_alt_y->SetPoint(2,70, 0.0);
  gr_alt_y->SetPoint(3,71,-0.5);
  gr_alt_y->SetPoint(4,72, -1.0);

  //
  TGraph *gr_az_x = new TGraph();
  gr_az_x->SetNameTitle("gr_az","gr_az");
  gr_az_x->SetPoint(0,177.5, -0.40);
  gr_az_x->SetPoint(1,  180, 0);
  gr_az_x->SetPoint(2,182.5, 0.35);

  TGraph *gr_az_y = new TGraph();
  gr_az_y->SetNameTitle("gr_az","gr_az");
  gr_az_y->SetPoint(0,177.5, -0.1);
  gr_az_y->SetPoint(1,  180, 0.0);
  gr_az_y->SetPoint(2,182.5, 0.1);
  

  TGraph *gr_alt_x_teo = new TGraph();
  gr_alt_x_teo->SetNameTitle("gr_alt_x_teo","gr_alt_x_teo");
  TGraph *gr_alt_y_teo = new TGraph();
  gr_alt_y_teo->SetNameTitle("gr_alt_y_teo","gr_alt_y_teo");
  //
  TGraph *gr_az_x_teo = new TGraph();
  gr_az_x_teo->SetNameTitle("gr_az_x_teo","gr_az_x_teo");
  TGraph *gr_az_y_teo = new TGraph();
  gr_az_y_teo->SetNameTitle("gr_az_y_teo","gr_az_y_teo");
  //
  Double_t azimuth = 180.0/180.0*TMath::Pi();
  Double_t altitude_min = 67.0/180.0*TMath::Pi();
  Double_t altitude_max = 73.0/180.0*TMath::Pi();
  Double_t altitude;
  Double_t x_shift;
  Double_t y_shift;
  //Double_t phi0_shift = 77.0;
  //Double_t phi0_shift = -101.0;
  Int_t nn = 100;
  for(Int_t i = 0; i<nn;i++){
    altitude = altitude_min + (altitude_max - altitude_min)/(nn-1)*i;
    get_x_y_shift(vx_tel, vy_tel, vz_tel,
		  azimuth,  altitude, x_shift, y_shift,
		  phi0_shift);
    gr_alt_x_teo->SetPoint(i,altitude*180.0/TMath::Pi(),x_shift);
    gr_alt_y_teo->SetPoint(i,altitude*180.0/TMath::Pi(),y_shift);
  }
  azimuth = 180.0/180.0*TMath::Pi();
  Double_t azimuth_min = 177.0/180.0*TMath::Pi();
  Double_t azimuth_max = 183.0/180.0*TMath::Pi();
  altitude = 70.0/180.0*TMath::Pi();
  for(Int_t i = 0; i<nn;i++){
    azimuth = azimuth_min + (azimuth_max - azimuth_min)/(nn-1)*i;
    get_x_y_shift(vx_tel, vy_tel, vz_tel,
		  azimuth,  altitude, x_shift, y_shift,
		  phi0_shift);
    gr_az_x_teo->SetPoint(i,azimuth*180.0/TMath::Pi(),x_shift);
    gr_az_y_teo->SetPoint(i,azimuth*180.0/TMath::Pi(),y_shift);
  }

  Double_t l2 = 0.0;
  //
  l2+=(gr_alt_x->Eval(68.0) - gr_alt_x_teo->Eval(68.0))*(gr_alt_x->Eval(68.0) - gr_alt_x_teo->Eval(68.0));
  l2+=(gr_alt_x->Eval(69.0) - gr_alt_x_teo->Eval(69.0))*(gr_alt_x->Eval(69.0) - gr_alt_x_teo->Eval(69.0));
  l2+=(gr_alt_x->Eval(70.0) - gr_alt_x_teo->Eval(70.0))*(gr_alt_x->Eval(70.0) - gr_alt_x_teo->Eval(70.0))*1000;
  l2+=(gr_alt_x->Eval(71.0) - gr_alt_x_teo->Eval(71.0))*(gr_alt_x->Eval(71.0) - gr_alt_x_teo->Eval(71.0));
  l2+=(gr_alt_x->Eval(71.0) - gr_alt_x_teo->Eval(71.0))*(gr_alt_x->Eval(71.0) - gr_alt_x_teo->Eval(71.0));
  //
  l2+=(gr_alt_y->Eval(68.0) - gr_alt_y_teo->Eval(68.0))*(gr_alt_y->Eval(68.0) - gr_alt_y_teo->Eval(68.0));
  l2+=(gr_alt_y->Eval(69.0) - gr_alt_y_teo->Eval(69.0))*(gr_alt_y->Eval(69.0) - gr_alt_y_teo->Eval(69.0));
  l2+=(gr_alt_y->Eval(70.0) - gr_alt_y_teo->Eval(70.0))*(gr_alt_y->Eval(70.0) - gr_alt_y_teo->Eval(70.0))*1000;
  l2+=(gr_alt_y->Eval(71.0) - gr_alt_y_teo->Eval(71.0))*(gr_alt_y->Eval(71.0) - gr_alt_y_teo->Eval(71.0));
  l2+=(gr_alt_y->Eval(71.0) - gr_alt_y_teo->Eval(71.0))*(gr_alt_y->Eval(71.0) - gr_alt_y_teo->Eval(71.0));
  //
  l2+=(gr_az_x->Eval(177.5) - gr_az_x_teo->Eval(177.5))*(gr_az_x->Eval(177.5) - gr_az_x_teo->Eval(177.5));
  l2+=(gr_az_x->Eval(180.0) - gr_az_x_teo->Eval(180.0))*(gr_az_x->Eval(180.0) - gr_az_x_teo->Eval(180.0))*1000;
  l2+=(gr_az_x->Eval(182.5) - gr_az_x_teo->Eval(182.5))*(gr_az_x->Eval(182.5) - gr_az_x_teo->Eval(182.5));
  //
  l2+=(gr_az_y->Eval(177.5) - gr_az_y_teo->Eval(177.5))*(gr_az_y->Eval(177.5) - gr_az_y_teo->Eval(177.5));
  l2+=(gr_az_y->Eval(180.0) - gr_az_y_teo->Eval(180.0))*(gr_az_y->Eval(180.0) - gr_az_y_teo->Eval(180.0))*1000;
  l2+=(gr_az_y->Eval(182.5) - gr_az_y_teo->Eval(182.5))*(gr_az_y->Eval(182.5) - gr_az_y_teo->Eval(182.5));

  
  cout<<phi0_shift<<" "<<l2<<endl;
  

  TCanvas *c1 = new TCanvas("c1","c1",10,10,1200,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
  // 
  c1->Divide(2,1);
  //
  c1->cd(1);
  gr_alt_x->SetLineColor(kBlack);
  gr_alt_x->SetLineWidth(3.0);
  gr_alt_x->SetMarkerColor(kBlack);
  gr_alt_x->SetMarkerStyle(20);
  //
  TMultiGraph *mg01 = new TMultiGraph();
  mg01->Add(gr_alt_x);
  mg01->Add(gr_alt_x_teo);
  mg01->Draw("apl");
  //mg->GetXaxis()->SetTitle("Time, a.u.");
  //mg->GetYaxis()->SetTitle("FADC counts");

  //
  c1->cd(2);
  gr_alt_y->SetLineColor(kBlack);
  gr_alt_y->SetLineWidth(3.0);
  gr_alt_y->SetMarkerColor(kBlack);
  gr_alt_y->SetMarkerStyle(20);
  TMultiGraph *mg02 = new TMultiGraph();
  mg02->Add(gr_alt_y);
  mg02->Add(gr_alt_y_teo);
  mg02->Draw("apl");


  TCanvas *c2 = new TCanvas("c2","c2",10,10,1200,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
  // 
  c2->Divide(2,1);
  //
  c2->cd(1);
  gr_az_x->SetLineColor(kBlack);
  gr_az_x->SetLineWidth(3.0);
  gr_az_x->SetMarkerColor(kBlack);
  gr_az_x->SetMarkerStyle(20);
  //
  TMultiGraph *mg03 = new TMultiGraph();
  mg03->Add(gr_az_x);
  mg03->Add(gr_az_x_teo);
  mg03->Draw("apl");
  //mg->GetXaxis()->SetTitle("Time, a.u.");
  //mg->GetYaxis()->SetTitle("FADC counts");

  //
  c2->cd(2);
  gr_az_y->SetLineColor(kBlack);
  gr_az_y->SetLineWidth(3.0);
  gr_az_y->SetMarkerColor(kBlack);
  gr_az_y->SetMarkerStyle(20);
  TMultiGraph *mg04 = new TMultiGraph();
  mg04->Add(gr_az_y);
  mg04->Add(gr_az_y_teo);
  mg04->Draw("apl");


  
  /*
  //
  gr_sim->SetLineColor(kRed);
  gr_sim->SetLineWidth(3.0);
  gr_sim->SetMarkerColor(kRed);
  gr_sim->SetMarkerStyle(20);
  */
  //
  //  

  //TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  //leg->AddEntry(gr, "simtelarr", "apl");
  //leg->AddEntry(gr_sim, "sim", "apl");
  //leg->Draw();  

  return l2;
}

void get_x_y_shift(const TVector3 &vx_tel, const TVector3 &vy_tel, const TVector3 &vz_tel,
		   Double_t azimuth, Double_t altitude, Double_t &x_shift, Double_t &y_shift,
		   Double_t phi0_shift){
  Double_t R_m_LST = 56.4;        //m
  Double_t F_m_LST = R_m_LST/2.0; //m
  Double_t theta = TMath::Pi()/2.0-altitude;
  Double_t phi = azimuth;
  Double_t theta_in_tel;
  Double_t phi_in_tel;
  get_part_coordinates_in_tel_frame(vx_tel, vy_tel, vz_tel,
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

void get_tel_frame( Double_t tel_theta, Double_t tel_phi, TVector3 &vx_tel, TVector3 &vy_tel, TVector3 &vz_tel){
  vx_tel.RotateZ(tel_phi);
  vx_tel.RotateY(-tel_theta);
  vy_tel.RotateZ(tel_phi);
  vy_tel.RotateY(-tel_theta);
  vz_tel.RotateZ(tel_phi);
  vz_tel.RotateY(-tel_theta);
}

void get_part_coordinates_in_tel_frame( const TVector3 &vx_tel, const TVector3 &vy_tel, const TVector3 &vz_tel,
					Double_t theta, Double_t phi, Double_t &theta_in_tel, Double_t &phi_in_tel){
  TVector3 part;
  part.SetMagThetaPhi(1.0,theta,phi);
  TVector3 v_in_tel(vx_tel.Dot(part),vy_tel.Dot(part),vz_tel.Dot(part));
  theta_in_tel = v_in_tel.Theta();
  phi_in_tel = v_in_tel.Phi();
}
