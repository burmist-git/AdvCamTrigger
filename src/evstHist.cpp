//my
#include "evstHist.hh"

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

using namespace std;

evstHist::evstHist(): evstHist("","")
{
}

evstHist::evstHist(const char* name, const char* title,
		   Double_t val_Emin, Double_t val_Emax, Int_t val_N_bins_E,
		   Double_t val_Thetamin, Double_t val_Thetamax, Int_t val_N_bins_t) : TH2Poly()
{    
  //
  _hist_name = name;
  _hist_title = title;
  _title="";
  //
  SetName(name);
  SetTitle(title);
  //
  //_Emin = 1.0;     //1   GeV
  //_Emax = 100000;  //100 TeV
  //_N_bins_E  = 25; //
  _Emin = val_Emin;
  _Emax = val_Emax;
  _N_bins_E  = val_N_bins_E;
  //
  //
  _Eminarr = new Double_t[_N_bins_E];
  _Emaxarr = new Double_t[_N_bins_E];
  Double_t *E_bins = new Double_t[_N_bins_E+1];
  //
  //
  Double_t bin_log;
  Double_t E_log_l;
  Double_t E_log_r;
  bin_log = (TMath::Log10(_Emax) - TMath::Log10(_Emin))/_N_bins_E;
  for(Int_t i = 0;i<_N_bins_E; i++){
    E_log_l = TMath::Log10(_Emin) + bin_log*i;
    E_log_r = TMath::Log10(_Emin) + bin_log*(i+1);
    _Eminarr[i] = TMath::Power(10,E_log_l);
    _Emaxarr[i] = TMath::Power(10,E_log_r); 
    E_bins[i] = _Eminarr[i];
    //cout<<E_bins[i]<<endl;
  }
  E_bins[_N_bins_E] = _Emax;
  //
  TString namehist=name;
  TString nametitle=title;
  //  
  TString namehist_E="_h1_E_";
  namehist_E += namehist;
  TString nametitle_E="_h1_E_";
  nametitle_E +=nametitle;
  _h1_E = new TH1D(namehist_E.Data(),nametitle_E.Data(), _N_bins_E, E_bins);
  //  
  //_Thetamin = 0.0;  //deg
  //_Thetamax = 10.0; //deg
  //_N_bins_t = 10;   //
  //
  _Thetamin = val_Thetamin;
  _Thetamax = val_Thetamax;
  _N_bins_t = val_N_bins_t;
  //
  TString namehist_theta="_h1_theta_";
  namehist_theta += namehist;
  TString nametitle_theta="_h1_theta_";
  nametitle_theta +=nametitle;
  _h1_theta = new TH1D(namehist_theta.Data(),nametitle_theta.Data(), _N_bins_t, _Thetamin, _Thetamax);
  //
  _Thetaminarr = new Double_t[_N_bins_t];
  _Thetamaxarr = new Double_t[_N_bins_t];
  for(Int_t i = 0;i<_N_bins_t; i++){
    _Thetaminarr[i] = (_Thetamax - _Thetamin)/_N_bins_t*i;
    _Thetamaxarr[i] = (_Thetamax - _Thetamin)/_N_bins_t*(i+1);
  }
  //
  Int_t n = 5;  
  Double_t* x;
  Double_t* y;
  x = new Double_t [n];
  y = new Double_t [n];
  for(Int_t i_E = 0;i_E<_N_bins_E;i_E++){
    for(Int_t i_t = 0;i_t<_N_bins_t;i_t++){
      //
      //we go clockwise
      x[0] = _Thetaminarr[i_t];
      y[0] = _Eminarr[i_E];
      //
      x[1] = x[0];
      y[1] = _Emaxarr[i_E];
      //
      x[2] = _Thetamaxarr[i_t];
      y[2] = y[1];
      //
      x[3] = x[2];
      y[3] = y[0];
      //
      x[4] = x[0];
      y[4] = y[0];
      //
      AddBin(n,x,y);
      //
      _v_cellID_th_bin_ID.push_back(i_t);
      _v_cellID_e_bin_ID.push_back(i_E);
    }
  }
  //
  TString h1_r_name_title;
  for(Int_t i_E = 0;i_E<_N_bins_E;i_E++){
    for(Int_t i_t = 0;i_t<_N_bins_t;i_t++){
      TH1D *h1 = new TH1D();
      h1_r_name_title = "h1_r_";
      h1_r_name_title += namehist;
      h1_r_name_title += "_i_E_"; h1_r_name_title += i_E;
      h1_r_name_title += "_i_t_"; h1_r_name_title += i_t;
      h1->SetNameTitle(h1_r_name_title.Data(),h1_r_name_title.Data());
      set_r_core_bins(h1);
      init_core_hist(h1);
      _v_r.push_back(h1);
    }
  }
  //
  _r_core_max_val = Get_r_core_max_val();
  _N_bins_r_core = get_v_r().at(0)->GetNbinsX();
}

evstHist::~evstHist(){
}

Double_t evstHist::Get_r_core_max_val(){
  Double_t bin_l = 0.0;
  Double_t bin_r = 0.0;
  get_Bin_Edge(get_v_r().at(0), get_v_r().at(0)->GetNbinsX(), bin_l, bin_r);
  return bin_r;
}

bool evstHist::Get_th_bin_ID_and_e_bin_ID( Int_t cellID, Int_t &th_bin_ID, Int_t &e_bin_ID){
  th_bin_ID = -999;
  e_bin_ID = -999;
  if( (cellID>=1) && ((unsigned int)cellID<=_v_cellID_th_bin_ID.size()) ){
    th_bin_ID = _v_cellID_th_bin_ID.at(((unsigned int)cellID-1));
    if( (unsigned int)cellID<=_v_cellID_e_bin_ID.size() ){
      e_bin_ID = _v_cellID_e_bin_ID.at(((unsigned int)cellID)-1);
      return true;
    }
  }
  return false;
}

void evstHist::Fill_rcore( Double_t th, Double_t E, Double_t r_core, Double_t weight_val){
  _v_r.at(get_bin_ID(E,th)-1)->Fill( r_core, weight_val);
}

void evstHist::set_r_core_bins(TH1D *h1, Double_t r_core_max){
  const Int_t nn = 11;
  Double_t* rBins = new Double_t[nn];
  //
  rBins[0] = 0;
  rBins[1] = 150;
  rBins[2] = 250;
  rBins[3] = 350;
  rBins[4] = 500;
  rBins[5] = 650;
  rBins[6] = 800;
  rBins[7] = 1000;
  rBins[8] = 1300;
  rBins[9] = 1600;
  rBins[10] = r_core_max;
  h1->SetBins(nn-1,rBins);
}

void evstHist::init_core_hist(TH1D *h1){
  for(Int_t i = 1;i<=h1->GetNbinsX();i++)
    h1->SetBinContent(i,0.0);
}

void evstHist::test(){
  cout<<"GetNcells() = "<<GetNcells()<<endl;
  for(Int_t i = 0;i<GetNcells();i++)
    SetBinContent(i,i);
  //SetBinContent(0,-999);
  //for(Int_t i = 0;i<GetNcells();i++)
  //cout<<i<<" "<<GetBinContent(i)<<endl;
  //cout<<259<<" "<<GetBinContent(259)<<endl;
  //cout<<1000<<" "<<GetBinContent(1000)<<endl;
  //cout<<-1<<" "<<GetBinContent(-1)<<endl;
  Draw_hist("evH_test.pdf");
}

void evstHist::test_Get_th_bin_ID_and_e_bin_ID(){
  Int_t th_bin_ID;
  Int_t e_bin_ID;
  for(Int_t i = 0;i<=GetNcells();i++){
    Get_th_bin_ID_and_e_bin_ID(i, th_bin_ID, e_bin_ID);
    SetBinContent(i,e_bin_ID+1);
  }
  Draw_hist("evH_test_Get_th_bin_ID_and_e_bin_ID_E.pdf");
  for(Int_t i = 0;i<=GetNcells();i++){
    Get_th_bin_ID_and_e_bin_ID(i, th_bin_ID, e_bin_ID);
    SetBinContent(i,th_bin_ID+1);
  }
  Draw_hist("evH_test_Get_th_bin_ID_and_e_bin_ID_th.pdf");
}

void evstHist::test_get_bin(Double_t E, Double_t th, Double_t val){
  SetBinContent(get_bin_ID(E,th),val);
}

Int_t evstHist::get_bin_ID( Double_t E, Double_t th){
  Int_t E_bin = _h1_E->FindBin(E);
  Int_t th_bin = _h1_theta->FindBin(th);
  return _N_bins_t*(E_bin-1) + th_bin;
}

void evstHist::init_h1_arb_v(const char* name, const char* title, Int_t nBins, Double_t valmin, Double_t valmax){
  Double_t E;
  Double_t th;
  Double_t r_core;
  unsigned int arbitrary_hist_ID;
  //
  for(Int_t i_E = 0;i_E<_N_bins_E;i_E++){
    for(Int_t i_t = 0;i_t<_N_bins_t;i_t++){
      for(Int_t i_r = 0;i_r<_N_bins_r_core;i_r++){      
	TH1D *h1 = new TH1D();
	_h1_arb_v.push_back(h1);
      }
    }
  }
  //  
  Int_t iBin;
  Int_t jBin;
  Int_t kBin;
  //
  TString h1_arb_name = name;
  TString h1_arb_title = title;
  //
  for(Int_t i = 1; i<=get_E_hist()->GetNbinsX(); i++){
    E = get_E_hist()->GetBinCenter(i);
    for(Int_t j = 1; j<=get_theta_hist()->GetNbinsX(); j++){
      th = get_theta_hist()->GetBinCenter(j);
      for(Int_t k = 1; k<=get_v_r().at(0)->GetNbinsX(); k++){
	r_core = get_v_r().at(0)->GetBinCenter(k);
	if(get_arbitrary_hist_ID( E, th, r_core, arbitrary_hist_ID)){
	  iBin = i-1;
	  jBin = j-1;
	  kBin = k-1;
	  //
	  h1_arb_name = name;
	  h1_arb_name += "_i_E_"; h1_arb_name += iBin;
	  h1_arb_name += "_i_t_"; h1_arb_name += jBin;
	  h1_arb_name += "_i_r_"; h1_arb_name += kBin;
	  //
	  h1_arb_title = title;
	  h1_arb_title += "_i_E_"; h1_arb_title += iBin;
	  h1_arb_title += "_i_t_"; h1_arb_title += jBin;
	  h1_arb_title += "_i_r_"; h1_arb_title += kBin;
	  //
	  get_h1_arb_v().at(arbitrary_hist_ID)->SetNameTitle(h1_arb_name.Data(),h1_arb_title.Data());
	  get_h1_arb_v().at(arbitrary_hist_ID)->SetBins(nBins, valmin, valmax);
	}
      }
    }
  }
}

void evstHist::fill_h1_arb_v(Double_t E, Double_t th, Double_t r_core, Double_t val, Double_t valweight){
  unsigned int arbitrary_hist_ID;
  if(get_arbitrary_hist_ID( E, th, r_core, arbitrary_hist_ID)){
    if(get_h1_arb_v().size()>arbitrary_hist_ID){
      if(valweight!=-999.0)
	get_h1_arb_v().at(arbitrary_hist_ID)->Fill(val,valweight);
      else
	get_h1_arb_v().at(arbitrary_hist_ID)->Fill(val);
    }
  }
}

void evstHist::test_get_arbitrary_hist_ID(){
  Double_t E;
  Double_t th;
  Double_t r_core = 10.0;
  unsigned int arbitrary_hist_ID;
  for(Int_t i = 1; i<=get_theta_hist()->GetNbinsX(); i++){
    th = get_theta_hist()->GetBinCenter(i);
    for(Int_t j = 1; j<=get_E_hist()->GetNbinsX(); j++){
      E = get_E_hist()->GetBinCenter(j);
      if(get_arbitrary_hist_ID( E, th, r_core, arbitrary_hist_ID))
	SetBinContent(get_bin_ID(E,th),arbitrary_hist_ID);
    }
  }
  Draw_hist("evH_test_get_arbitrary_hist_ID_r_core_10m.pdf");
  //
  r_core = 200.0;
  for(Int_t i = 1; i<=get_theta_hist()->GetNbinsX(); i++){
    th = get_theta_hist()->GetBinCenter(i);
    for(Int_t j = 1; j<=get_E_hist()->GetNbinsX(); j++){
      E = get_E_hist()->GetBinCenter(j);
      if(get_arbitrary_hist_ID( E, th, r_core, arbitrary_hist_ID))
	SetBinContent(get_bin_ID(E,th),arbitrary_hist_ID);
    }
  }
  Draw_hist("evH_test_get_arbitrary_hist_ID_r_core_200m.pdf");
  //
  r_core = 1800.0;
  for(Int_t i = 1; i<=get_theta_hist()->GetNbinsX(); i++){
    th = get_theta_hist()->GetBinCenter(i);
    for(Int_t j = 1; j<=get_E_hist()->GetNbinsX(); j++){
      E = get_E_hist()->GetBinCenter(j);
      if(get_arbitrary_hist_ID( E, th, r_core, arbitrary_hist_ID))
	SetBinContent(get_bin_ID(E,th),arbitrary_hist_ID);
    }
  }
  Draw_hist("evH_test_get_arbitrary_hist_ID_r_core_1800m.pdf");
}

bool evstHist::get_arbitrary_hist_ID( Double_t E, Double_t th, Double_t r_core, unsigned int &arbitrary_hist_ID){
  Int_t E_bin  = _h1_E->FindBin(E);
  Int_t th_bin = _h1_theta->FindBin(th);
  Int_t r_bin  = 0;
  arbitrary_hist_ID = -1;  
  if(get_v_r().size()>0)
    r_bin = get_v_r().at(0)->FindBin(r_core);
  if( E_bin < 1 || th_bin < 1 || r_bin < 1)
    return false;
  if( E_bin > _N_bins_E || th_bin > _N_bins_t || r_bin > _N_bins_r_core)
    return false;
  //
  arbitrary_hist_ID = (unsigned int)(r_bin-1)*(_N_bins_E*_N_bins_t) + (unsigned int)(_N_bins_t*(E_bin-1) + th_bin) - 1;
  //
  //cout<<setw(15)<<"E"
  //    <<setw(15)<<"th"
  //    <<setw(15)<<"r_core"
  //    <<setw(15)<<"hist_ID"<<endl;
  //cout<<setw(15)<<E
  //    <<setw(15)<<th
  //    <<setw(15)<<r_core<<endl;
  //cout<<setw(15)<<E_bin
  //    <<setw(15)<<th_bin
  //    <<setw(15)<<r_bin
  //    <<setw(15)<<arbitrary_hist_ID<<endl;
  //
  return true;
}

TCanvas* evstHist::Draw_hist(TString fileName, TString frame_title){

  //kDeepSea=51,          kGreyScale=52,    kDarkBodyRadiator=53,
  //kBlueYellow= 54,      kRainBow=55,      kInvertedDarkBodyRadiator=56,
  //kBird=57,             kCubehelix=58,    kGreenRedViolet=59,
  //kBlueRedYellow=60,    kOcean=61,        kColorPrintableOnGrey=62,
  //kAlpine=63,           kAquamarine=64,   kArmy=65,
  //kAtlantic=66,         kAurora=67,       kAvocado=68,
  //kBeach=69,            kBlackBody=70,    kBlueGreenYellow=71,
  //kBrownCyan=72,        kCMYK=73,         kCandy=74,
  //kCherry=75,           kCoffee=76,       kDarkRainBow=77,
  //kDarkTerrain=78,      kFall=79,         kFruitPunch=80,
  //kFuchsia=81,          kGreyYellow=82,   kGreenBrownTerrain=83,
  //kGreenPink=84,        kIsland=85,       kLake=86,
  //kLightTemperature=87, kLightTerrain=88, kMint=89,
  //kNeon=90,             kPastel=91,       kPearl=92,
  //kPigeon=93,           kPlum=94,         kRedBlue=95,
  //kRose=96,             kRust=97,         kSandyTerrain=98,
  //kSienna=99,           kSolar=100,       kSouthWest=101,
  //kStarryNight=102,     kSunset=103,      kTemperatureMap=104,
  //kThermometer=105,     kValentine=106,   kVisibleSpectrum=107,
  //kWaterMelon=108,      kCool=109,        kCopper=110,
  //kGistEarth=111,       kViridis=112,     kCividis=113

  //gStyle->SetPalette(kInvertedDarkBodyRadiator);

  //
  TString c1_name = "c1_";
  c1_name += _hist_name;
  TString c1_title = "c1_";
  c1_title += _hist_title;
  //
  TCanvas *c1 = new TCanvas(c1_name.Data(),c1_title.Data(),1500,1000);
  c1->SetRightMargin(0.12);
  c1->SetLeftMargin(0.12);
  c1->SetTopMargin(0.1);
  c1->SetBottomMargin(0.15);
  c1->SetLogy();
  c1->SetLogz();
  
  TH2F *frame = new TH2F("h2","h2", 40, _Thetamin, _Thetamax, 40, _Emin, _Emax);
  //frame->SetTitle("Proton rate, Hz");
  //frame->SetTitle(_title.Data());
  frame->SetTitle(frame_title.Data());
  frame->GetXaxis()->SetTitle("Theta, deg");
  frame->GetYaxis()->SetTitle("Energy, GeV");
  //frame->GetXaxis()->CenterTitle();
  //frame->GetYaxis()->CenterTitle();
  //frame->GetYaxis()->SetTitleOffset(1.5);
  frame->SetStats(kFALSE);
  frame->Draw();
  //
  //SetMaximum(50);
  //SetMinimum(0);
  //
  SetMarkerSize(0.8);
  Draw("same ZCOLOR TEXT");
  //Draw("same TEXT");
  //Draw("same ZCOLOR");
  //c1->SaveAs("Draw_map_test.pdf");
  if(fileName != "")
  c1->SaveAs(fileName.Data());
  return c1;
}

TCanvas* evstHist::Draw_hist_core(Int_t e_bin_i, TString fileName, TString frame_title){
  //
  TString c1_name = "c1_";
  c1_name += _hist_name;
  c1_name += "_";
  c1_name += e_bin_i;
  TString c1_title = "c1_";
  c1_title += _hist_title;
  c1_title += "_";
  c1_title += e_bin_i;
  //
  TCanvas *c1 = new TCanvas(c1_name.Data(),c1_title.Data(),1500,1000);
  c1->SetRightMargin(0.12);
  c1->SetLeftMargin(0.12);
  c1->SetTopMargin(0.1);
  c1->SetBottomMargin(0.15);
  //
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE);
  //
  Int_t th_bin_ID;
  Int_t e_bin_ID;
  //
  Int_t n_plots = 0;
  //
  vector<Int_t> color_id;
  color_id.push_back(kBlack);     //0
  color_id.push_back(kBlue);      //1
  color_id.push_back(kRed);       //2
  color_id.push_back(kGreen);     //3
  color_id.push_back(kMagenta);   //4
  color_id.push_back(kYellow+2);  //5
  color_id.push_back(kBlue+2);    //6
  color_id.push_back(kRed+2);     //7
  color_id.push_back(kGreen+2);   //8
  color_id.push_back(kMagenta+2); //9
  //
  TLegend *leg = new TLegend(0.7,0.7,0.99,0.99,"","brNDC");
  //
  TString n_plots_str = " ";
  //
  for(Int_t i = 0;i<_N_bins_E*_N_bins_t;i++){
    Get_th_bin_ID_and_e_bin_ID( i, th_bin_ID, e_bin_ID);
    if(e_bin_i == e_bin_ID){
      if(n_plots == 0){
	_v_r.at(i)->SetMaximum(1.1);
	_v_r.at(i)->Draw();
      }
      else{
	_v_r.at(i)->Draw("same");
      }
      if(n_plots<color_id.size())
	_v_r.at(i)->SetLineColor(color_id.at((unsigned int)n_plots));
      _v_r.at(i)->SetLineWidth(2.0);
      n_plots_str = " ";
      n_plots_str += n_plots;
      leg->AddEntry(_v_r.at(i), n_plots_str.Data(), "apl");
      n_plots++;
    }
  }
  leg->Draw("same");
  //SetMarkerSize(0.8);
  //Draw("same ZCOLOR TEXT");
  //Draw("same TEXT");
  //Draw("same ZCOLOR");
  //c1->SaveAs("Draw_map_test.pdf");
  if(fileName != "")
  c1->SaveAs(fileName.Data());
  return c1;
}

void evstHist::Divideh1(evstHist *evH_cut, evstHist *evH_all, Double_t normval){
  //
  Double_t nev;
  Double_t nev_norm;
  //
  for(Int_t i = 1;i<=_h1_E->GetNbinsX();i++){
    nev=evH_cut->get_E_hist()->GetBinContent(i);
    nev_norm=evH_all->get_E_hist()->GetBinContent(i);
    //    
    if(nev_norm != 0.0)
      _h1_E->SetBinContent(i,nev/nev_norm*normval);
    else
      _h1_E->SetBinContent(i,0.0);
  }
  //
  for(Int_t i = 1;i<=_h1_theta->GetNbinsX();i++){
    nev=evH_cut->get_theta_hist()->GetBinContent(i);
    nev_norm=evH_all->get_theta_hist()->GetBinContent(i);
    //    
    if(nev_norm != 0.0)
      _h1_theta->SetBinContent(i,nev/nev_norm*normval);
    else
      _h1_theta->SetBinContent(i,0.0);
  }
}

void evstHist::Divide(evstHist *evH_cut, evstHist *evH_all, bool with_r_core){
  Double_t nev;
  Double_t nev_norm;
  for(Int_t i = 1;i<=GetNcells();i++){
    nev=evH_cut->GetBinContent(i);
    nev_norm=evH_all->GetBinContent(i);
    if(nev_norm != 0.0)
      SetBinContent(i,nev/nev_norm);
    else
      SetBinContent(i,0.0);
  }
  if(with_r_core){
    if(evH_cut->get_v_r().size() == evH_all->get_v_r().size()){
      if(check_bin_compatibility(evH_cut, with_r_core)){
	if(check_bin_compatibility(evH_all, with_r_core)){
	  for(unsigned int ii = 0;ii<get_v_r().size();ii++){
	    for(Int_t i = 1;i<=evH_cut->get_v_r().at(ii)->GetNbinsX();i++){
	      nev=evH_cut->get_v_r().at(ii)->GetBinContent(i);
	      nev_norm=evH_all->get_v_r().at(ii)->GetBinContent(i);
	      if(nev_norm != 0.0)
		get_v_r().at(ii)->SetBinContent(i,nev/nev_norm);
	      else
		get_v_r().at(ii)->SetBinContent(i,0.0);	      
	    }
	  }
	}
      }
    }
    else{
      cout<<"evH_cut->get_v_r().size() == evH_all->get_v_r().size()"<<endl
	  <<"evH_cut->get_v_r().size() == "<<evH_cut->get_v_r().size()<<endl
	  <<"evH_all->get_v_r().size() == "<<evH_all->get_v_r().size()<<endl;
    }
  }
}

bool evstHist::check_bin_compatibility(const evstHist *evH, bool with_r_core){
  if(evH->get_Emin() != get_Emin())
    return false;
  if(evH->get_Emax() != get_Emax())
    return false;
  if(evH->get_N_bins_E() !=get_N_bins_E())
    return false;
  if(evH->get_Thetamin() != get_Thetamin())
    return false;
  if(evH->get_Thetamax() != get_Thetamax())
    return false;
  if(evH->get_N_bins_t() != get_N_bins_t())
    return false;
  if(with_r_core)
    if(!check_r_core_bin_compatibility(evH))
      return false;
  return true;
}

bool evstHist::check_r_core_bin_compatibility(const evstHist *evH){
  return true;
}

void evstHist::Multiply_core(evstHist *evH){
  Double_t nev;
  Double_t nev_norm;
  for(unsigned int ii = 0;ii<get_v_r().size();ii++){
    nev_norm=evH->GetBinContent(ii+1);
    SetBinContent(ii+1,nev_norm);
    for(Int_t i = 1;i<=get_v_r().at(ii)->GetNbinsX();i++){
      nev=get_v_r().at(ii)->GetBinContent(i);
      get_v_r().at(ii)->SetBinContent(i,nev*nev_norm);
    }
  }
}

void evstHist::Multiply(evstHist *evH_eff, evstHist *evH_flux, bool with_r_core){
  Double_t nev;
  Double_t nev_norm;
  for(Int_t i = 1;i<=GetNcells();i++){
    nev=evH_eff->GetBinContent(i);
    nev_norm=evH_flux->GetBinContent(i);
    SetBinContent(i,nev*nev_norm);
  }
  if(with_r_core){
    if(evH_eff->get_v_r().size() == evH_flux->get_v_r().size()){
      if(check_bin_compatibility(evH_eff, with_r_core)){
	if(check_bin_compatibility(evH_flux, with_r_core)){
	  for(unsigned int ii = 0;ii<get_v_r().size();ii++){
	    for(Int_t i = 1;i<=evH_eff->get_v_r().at(ii)->GetNbinsX();i++){
	      nev=evH_eff->get_v_r().at(ii)->GetBinContent(i);
	      nev_norm=evH_flux->get_v_r().at(ii)->GetBinContent(i);
	      get_v_r().at(ii)->SetBinContent(i,nev*nev_norm);
	    }
	  }
	}
      }
    }
    else{
      cout<<"evH_eff->get_v_r().size()  == evH_flux->get_v_r().size()"<<endl
	  <<"evH_eff->get_v_r().size()  == "<<evH_eff->get_v_r().size()<<endl
	  <<"evH_flux->get_v_r().size() == "<<evH_flux->get_v_r().size()<<endl;
    }
  }
}

Double_t evstHist::get_Weight_ETeV(Double_t ETev){
  if(ETev>0.0)
    return 1.0/TMath::Power(ETev,0.68);
  return 0.0;
}

Double_t evstHist::get_Weight_ETeV_soft(Double_t ETev){
  if(ETev>0.0)
    return 1.0/TMath::Power(ETev,0.22);
  return 0.0;
}

Double_t evstHist::Weight_integral_GeV( Double_t e_r_GeV, Double_t e_l_GeV){
  if(e_r_GeV<0.0 || e_l_GeV<0.0){
    cout<<" ERROR --> e_r_GeV<0.0 || e_l_GeV<0.0"<<endl
	<<"           e_r_GeV = "<<e_r_GeV<<endl
      	<<"           e_l_GeV = "<<e_l_GeV<<endl;
    assert(0);
  }
  if(e_r_GeV <= e_l_GeV){
    cout<<" ERROR --> e_r_GeV <= e_l_GeV"<<endl
	<<"           e_r_GeV = "<<e_r_GeV<<endl
      	<<"           e_l_GeV = "<<e_l_GeV<<endl;
    assert(0);
  }
  return 3.125*(TMath::Power(e_r_GeV/1000.0,0.32)-TMath::Power(e_l_GeV/1000.0,0.32));
}

void evstHist::Weight_TeV(){
  //
  Double_t nev;
  Double_t nev_norm;
  Double_t e_l;
  Double_t e_r;
  Int_t th_bin_ID;
  Int_t e_bin_ID;
  for(Int_t i = 0;i<=GetNcells();i++){
    nev=GetBinContent(i);
    if(Get_th_bin_ID_and_e_bin_ID(i, th_bin_ID, e_bin_ID)){
      e_l = get_E_hist()->GetBinCenter(e_bin_ID+1) - get_E_hist()->GetBinWidth(e_bin_ID+1)/2.0;
      e_r = get_E_hist()->GetBinCenter(e_bin_ID+1) + get_E_hist()->GetBinWidth(e_bin_ID+1)/2.0;    
      nev_norm=Weight_integral_GeV(e_r,e_l)/((e_r - e_l)/1000.0);
      SetBinContent(i,nev*nev_norm);
      unsigned int ii = i-1;
      for(Int_t j = 1;j<=get_v_r().at(ii)->GetNbinsX();j++){
	nev=get_v_r().at(ii)->GetBinContent(j);
	get_v_r().at(ii)->SetBinContent(j,nev*nev_norm);
      }
    }
  }
}

void evstHist::selfNorm(){
  Double_t nev;
  Double_t nev_norm;
  for(unsigned int ii = 0;ii<get_v_r().size();ii++){
    nev_norm=get_v_r().at(ii)->Integral();
    for(Int_t i = 1;i<=get_v_r().at(ii)->GetNbinsX();i++){
      nev=get_v_r().at(ii)->GetBinContent(i);
      if(nev_norm != 0.0)
	get_v_r().at(ii)->SetBinContent(i,nev/nev_norm);
      else
	get_v_r().at(ii)->SetBinContent(i,0.0);	      
    }
  }
  for(Int_t i = 1;i<=GetNcells();i++){
    nev=GetBinContent(i);
    nev_norm=GetBinContent(i);
    if(nev_norm != 0.0)
      SetBinContent(i,nev/nev_norm);
    else
      SetBinContent(i,0.0);
  }
}

void evstHist::DumpBinContent(TString data_out, bool with_r_core){
  ofstream myfile;
  myfile.open (data_out.Data());
  if(!with_r_core){
    for(Int_t i = 1;i<=GetNcells();i++){
      myfile<<setprecision(20)<<GetBinContent(i)<<endl;
    }
  }
  //
  if(with_r_core){
    for(Int_t i = 1;i<=GetNcells();i++){
      myfile<<setprecision(20)<<GetBinContent(i)<<" ";
      if(i<=_N_bins_E*_N_bins_t){
	for(Int_t j = 1;j<=_v_r.at((unsigned int)(i-1))->GetNbinsX();j++)
	  myfile<<setprecision(20)<<_v_r.at((unsigned int)(i-1))->GetBinContent(j)<<" ";
      }
      myfile<<endl;
    }
  }
  myfile.close();
}
  
void evstHist::LoadBinContent(TString data_in, bool with_r_core){
  ifstream myfile(data_in.Data());
  Double_t val;
  Int_t i = 1;
  if (myfile.is_open()) {
    if(!with_r_core){
      while (myfile>>val){
	SetBinContent(i,val);
	i++;
      }
    }
    if(with_r_core){
      while (myfile>>val){
	SetBinContent(i,val);
	if(i<=_N_bins_E*_N_bins_t){
	  for(Int_t j = 1;j<=_v_r.at((unsigned int)(i-1))->GetNbinsX();j++){
	    myfile>>val;
	    _v_r.at((unsigned int)(i-1))->SetBinContent(j,val);
	  }
	  i++;
	}
      }
    }
    myfile.close();
  }
}

const void evstHist::get_Bin_Edge(const TH1D *h1, Int_t binID, Double_t &bin_l, Double_t &bin_r){
  bin_l = h1->GetBinLowEdge(binID+1);
  bin_r = h1->GetBinLowEdge(binID+1) + h1->GetBinWidth(binID+1);
}
  
const void evstHist::Print_hist_BinsInfo(Int_t binE, Int_t binTheta, Int_t binDist){
  PrintBinsInfo(get_E_hist(),binE);
  PrintBinsInfo(get_theta_hist(), binTheta);
  PrintBinsInfo(get_v_r().at(0), binDist);
}

const void evstHist::PrintBinsInfo(const TH1D *h1, Int_t binID){
  Double_t bin_l;
  Double_t bin_r;
  if(binID<0)
    std::cout<<"Name  : "<<h1->GetName()<<std::endl
	     <<"Title : "<<h1->GetTitle()<<std::endl;
  for(Int_t i = 1; i<=h1->GetNbinsX();i++){
    bin_l = h1->GetBinLowEdge(i);
    bin_r = h1->GetBinLowEdge(i) + h1->GetBinWidth(i);
    if(binID<0)
      std::cout<<setw(5)<<(i-1)<<" "<<setw(20)<<bin_l<<" "<<setw(20)<<bin_r<<std::endl;
    else
      if(binID == (i-1))
	std::cout<<setw(5)<<(i-1)<<" "<<setw(20)<<bin_l<<" "<<setw(20)<<bin_r<<std::endl;
  }
}

Double_t evstHist::GetTotIntegral(){
  Double_t integral_v = 0.0;
  for(Int_t ii = 1;ii<=(_N_bins_E*_N_bins_t);ii++)
    integral_v+=GetBinContent(ii);
  return integral_v;
}

Double_t evstHist::GetIntegral(Double_t e_min, Double_t e_max, Double_t theta_min, Double_t theta_max) const {
  Double_t integral_v = 0.0;
  Int_t i_cell;
  //
  Double_t e_min_bin;
  Double_t e_max_bin;
  Double_t theta_min_bin;
  Double_t theta_max_bin;
  //
  for(Int_t i_theta = 0;i_theta<_N_bins_t;i_theta++){
    for(Int_t i_E = 0;i_E<_N_bins_E;i_E++){
      i_cell = (i_E)*_N_bins_t+(i_theta+1);
      //
      e_min_bin = _h1_E->GetBinLowEdge(i_E+1);
      e_max_bin = _h1_E->GetBinLowEdge(i_E+1) + _h1_E->GetBinWidth(i_E+1);
      //
      theta_min_bin = _h1_theta->GetBinLowEdge(i_theta+1);
      theta_max_bin = _h1_theta->GetBinLowEdge(i_theta+1) + _h1_theta->GetBinWidth(i_theta+1);
      //
      if(e_min_bin>=e_min && e_max_bin<=e_max){
	if(theta_min_bin>=theta_min && theta_max_bin<=theta_max){
	  //cout<<"i_cell = "<<i_cell<<endl;
	  integral_v+=GetBinContent(i_cell);
	}
      }
    }
  }
  return integral_v;
}
