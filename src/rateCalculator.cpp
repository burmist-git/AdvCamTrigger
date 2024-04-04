//my
#include "rateCalculator.hh"
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
#include <vector>

//root
#include <TString.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TSystemFile.h>
#include <TFile.h>
#include <TMath.h>

using namespace std;

rateCalculator::rateCalculator():
  _name(""),
  _title(""),
  _hist_file_prefix(""),
  _output_hist_file(""),
  _evH_flux(NULL),
  _h1_rate(NULL),
  _disable_energy_theta_rcore_binwise_cuts(false),
  _h1_N_dbc(NULL),
  _h1_ev_per_job(NULL),
  _h1_N_dbc_mean(NULL),
  _h1_dbc_number_of_points_tot_w(NULL),
  _h1_rate_vs_th(NULL),
  _h1_rate_NSB_vs_th(NULL),
  _h1_dbc_number_of_points_NSB(NULL),
  _h1_N_dbc_mean_NSB(NULL),
  _h1_rate_tot_vs_th(NULL),
  _h1_dbc_mean_time_ii(NULL),
  _h1_dbc_mean_time_ii_NSB(NULL),
  __h1_energy_eff_r(NULL),
  _h1_energy_all(NULL),
  _h1_energy_trg(NULL),
  _h1_energy_eff_r(NULL),
  _n_ev_tot(0.0),
  _n_jobs(-999),
  _rsimulation(801)
{
}

rateCalculator::rateCalculator( const char* name, const char* title, TString hist_file_prefix, TString output_hist_file, evstHist* evH_flux,
				Int_t n_jobs, bool disable_energy_theta_rcore_binwise_cuts, Int_t n_jobs_NSB, TString hist_file_dir_NSB, Double_t rsimulation):
  _name(name),
  _title(title),
  _hist_file_prefix(hist_file_prefix),
  _output_hist_file(output_hist_file),
  _evH_flux(evH_flux),
  _h1_rate(new TH1D),
  _disable_energy_theta_rcore_binwise_cuts(disable_energy_theta_rcore_binwise_cuts),
  _h1_N_dbc(new TH1D),
  _h1_dbc_number_of_points_tot_w(new TH1D),
  _h1_rate_vs_th(new TH1D),
  _h1_rate_NSB_vs_th(new TH1D),
  _h1_dbc_number_of_points_NSB(new TH1D),
  _h1_rate_tot_vs_th(new TH1D),
  _h1_dbc_mean_time_ii(new TH1D),
  _h1_dbc_mean_time_ii_NSB(new TH1D),
  __h1_energy_eff_r(new TH1D),
  _h1_energy_all(new TH1D),
  _h1_energy_trg(new TH1D),
  _h1_energy_eff_r(new TH1D),
  _n_ev_tot(0.0),
  _n_jobs(n_jobs),
  _rsimulation(801)
{
  TString h1_rate_name = "_h1_rate";
  h1_rate_name += "_";
  h1_rate_name += _name;
  TString h1_rate_title = "_h1_rate";
  h1_rate_title += "_";
  h1_rate_title += _title;
  _h1_rate->SetNameTitle(h1_rate_name.Data(),h1_rate_title.Data());
  //
  _h1_ev_per_job = new TH1D("_h1_ev_per_job","_h1_ev_per_job",1000,0.0,20000.0);  
  //
  _h1_N_dbc->SetNameTitle("_h1_N_dbc","_h1_N_dbc");
  _h1_N_dbc_mean = new TH1D("_h1_N_dbc_mean","_h1_N_dbc_mean",5000,0.0,4.0);
  _h1_dbc_number_of_points_tot_w->SetNameTitle("_h1_dbc_number_of_points_tot_w","_h1_dbc_number_of_points_tot_w");
  _h1_rate_vs_th->SetNameTitle("_h1_rate_vs_th","_h1_rate_vs_th");
  _h1_rate_NSB_vs_th->SetNameTitle("_h1_rate_NSB_vs_th","_h1_rate_NSB_vs_th");
  _h1_dbc_number_of_points_NSB->SetNameTitle("_h1_dbc_number_of_points_NSB","_h1_dbc_number_of_points_NSB");
  _h1_N_dbc_mean_NSB = new TH1D("_h1_N_dbc_mean_NSB","_h1_N_dbc_mean_NSB",5000,0.0,2.0);
  _h1_rate_tot_vs_th->SetNameTitle("_h1_rate_tot_vs_th","_h1_rate_tot_vs_th");
  //
  _h1_dbc_mean_time_ii->SetNameTitle("_h1_dbc_mean_time_ii","_h1_dbc_mean_time_ii");
  _h1_dbc_mean_time_ii_NSB->SetNameTitle("_h1_dbc_mean_time_ii_NSB","_h1_dbc_mean_time_ii_NSB");      
  //
  __h1_energy_eff_r->SetNameTitle("__h1_energy_eff_r","__h1_energy_eff_r");
  _h1_energy_all->SetNameTitle("_h1_energy_all","_h1_energy_all");
  _h1_energy_trg->SetNameTitle("_h1_energy_trg","_h1_energy_trg");
  _h1_energy_eff_r->SetNameTitle("_h1_energy_eff_r","_h1_energy_eff_r");  
  //
  if(!_disable_energy_theta_rcore_binwise_cuts){
    calculate_rate_binwise();
  }
  else{
    calculate_rate();
    calculate_rate_NSB( n_jobs_NSB, hist_file_dir_NSB);
    calculate_rate_tot();
  }
  save_output();
}

rateCalculator::~rateCalculator(){
}

void rateCalculator::save_output(){
  TFile* rootFile = new TFile(_output_hist_file.Data(), "RECREATE", " Histograms", 1);
  rootFile->cd();
  if (rootFile->IsZombie()){
    cout<<"  ERROR ---> file "<<_output_hist_file.Data()<<" is zombi"<<endl;
    assert(0);
  }
  else
    cout<<"  Output Histos file ---> "<<_output_hist_file.Data()<<endl;
  //
  for(unsigned int i = 0;i<_h1_dbc_number_of_points_w_v.size();i++)
    _h1_dbc_number_of_points_w_v.at(i)->Write();
  //
  _h1_ev_per_job->Write();
  _h1_N_dbc->Write();
  _h1_rate->Write();
  _h1_N_dbc_mean->Write();
  _h1_dbc_number_of_points_tot_w->Write();
  //
  _h1_rate_vs_th->Write();
  _h1_dbc_number_of_points_NSB->Write();
  _h1_N_dbc_mean_NSB->Write();
  _h1_rate_NSB_vs_th->Write();
  _h1_rate_tot_vs_th->Write();
  //
  _h1_dbc_mean_time_ii->Write();
  _h1_dbc_mean_time_ii_NSB->Write();
  //
  __h1_energy_eff_r->Write();
  _h1_energy_all->Write();
  _h1_energy_trg->Write();
  _h1_energy_eff_r->Write();
  //
  rootFile->Close();  
  //
  //cout<<"_n_ev_tot = "<<_n_ev_tot/112643298.0<<endl;
}

void rateCalculator::copyHist(TH1D *h1, TH1D *h1_to_cp){
  int nBins = h1_to_cp->GetNbinsX();
  double *bins_low_edge= new double[nBins+1];
  for(int i = 1;i<=nBins;i++)
    bins_low_edge[i-1] = h1_to_cp->GetBinLowEdge(i);
  bins_low_edge[nBins] = h1_to_cp->GetBinLowEdge(nBins) + h1_to_cp->GetBinWidth(nBins);
  h1->SetBins(nBins,bins_low_edge);
  //  
  for(Int_t i = 1;i<=h1_to_cp->GetNbinsX();i++)
    h1->SetBinContent(i,h1_to_cp->GetBinContent(i));
  h1->SetEntries(h1_to_cp->GetEntries());
}

void rateCalculator::addHist(TH1D *h1, TH1D *h1_to_add){
  for(Int_t i = 1;i<=h1_to_add->GetNbinsX();i++)
    h1->SetBinContent(i,h1_to_add->GetBinContent(i)+h1->GetBinContent(i));
}

void rateCalculator::get_histogram( TString hist_file_name, TH1D *h1, TString hist_name_to_get){
  TFile *f1 = new TFile(hist_file_name.Data());
  TH1D *h1_tmp = (TH1D*)f1->Get(hist_name_to_get.Data());
  copyHist(h1,h1_tmp);
  f1->Close();
}

void rateCalculator::append_histogram(TH1D *h1, TH1D *h1_to_add){
  if(h1->GetNbinsX()==1)
    copyHist(h1,h1_to_add);
  else
    addHist(h1,h1_to_add);
}

void rateCalculator::calculate_rate_binwise(){
  //
  TString hist_file_name;
  Int_t i_binE;
  Int_t i_binTheta;
  Int_t i_binDist;
  //
  for(i_binE = 0;i_binE<_evH_flux->get_E_hist()->GetNbinsX();i_binE++){
    for(i_binTheta = 0;i_binTheta<_evH_flux->get_theta_hist()->GetNbinsX();i_binTheta++){
      for(i_binDist = 0;i_binDist<_evH_flux->get_v_r().at(0)->GetNbinsX();i_binDist++){
	if(get_histogram_file_name(_hist_file_prefix, i_binE, i_binTheta, i_binDist, hist_file_name)){
	  TH1D *h1 = new TH1D();
	  TString h1_name_title = "h1_";
	  h1_name_title += i_binE; h1_name_title += "binE_";
	  h1_name_title += i_binTheta; h1_name_title += "binTheta_";
	  h1_name_title += i_binDist; h1_name_title += "binDist";
	  h1->SetNameTitle(h1_name_title.Data(),h1_name_title.Data());
	  get_histogram(hist_file_name, h1, "h1_dbc_number_of_points_w");
	  _h1_dbc_number_of_points_w_v.push_back(h1);
	  append_histogram( _h1_rate, h1);
	}
      }
    }
  } 
  //if(get_histogram_file_name(_hist_file_prefix, 5, 0, 0, hist_file_name))
  //cout<<"hist_file_name = "<<hist_file_name<<endl;
  //if(get_histogram_file_name(_hist_file_prefix, 15, 3, 9, hist_file_name))
  //cout<<"hist_file_name = "<<hist_file_name<<endl;
}

void rateCalculator::calculate_rate_tot(){
  copyHist(_h1_rate_tot_vs_th,_h1_rate_vs_th);
  addHist(_h1_rate_tot_vs_th,_h1_rate_NSB_vs_th);
}

void rateCalculator::calculate_rate(){
  //
  TString hist_file_name;
  Int_t i_binE = 0;
  Int_t i_binTheta = 0;
  Int_t i_binDist = 0;
  //
  for(Int_t i = 0;i<_n_jobs;i++){
    if(get_histogram_file_name(_hist_file_prefix, i_binE, i_binTheta, i_binDist, hist_file_name, i)){
      cout<<hist_file_name<<endl;
      TH1D *h1 = new TH1D();
      TString h1_name_title = "h1_dbc_number_of_points_w_";
      h1_name_title += i; h1_name_title += "ID";
      h1->SetNameTitle(h1_name_title.Data(),h1_name_title.Data());
      get_histogram(hist_file_name, h1, "h1_dbc_number_of_points_w");
      _h1_dbc_number_of_points_w_v.push_back(h1);
      append_histogram(_h1_rate,h1);
      //
      TH1D *htmp = new TH1D();
      get_histogram(hist_file_name, htmp, "h1_npe");
      _n_ev_tot += htmp->GetEntries();
      _h1_ev_per_job->Fill(htmp->GetEntries());
      delete htmp;
      //
      htmp = new TH1D();
      get_histogram(hist_file_name, htmp, "h1_N_dbc");
      append_histogram(_h1_N_dbc,htmp);
      _h1_N_dbc_mean->Fill(htmp->GetMean());
      delete htmp;      
      //
      htmp = new TH1D();
      get_histogram(hist_file_name, htmp, "h1_dbc_number_of_points_tot_w");
      append_histogram(_h1_dbc_number_of_points_tot_w,htmp);
      delete htmp;      
      //
      htmp = new TH1D();
      get_histogram(hist_file_name, htmp, "h1_dbc_mean_time_ii");
      append_histogram(_h1_dbc_mean_time_ii,htmp);
      delete htmp;      
      //
      htmp = new TH1D();
      get_histogram(hist_file_name, htmp, "_h1_energy_eff_r");
      append_histogram(__h1_energy_eff_r,htmp);
      delete htmp;      
      //
      htmp = new TH1D();
      get_histogram(hist_file_name, htmp, "_h1_energy_all");
      append_histogram(_h1_energy_all,htmp);
      delete htmp;      
      //
      htmp = new TH1D();
      get_histogram(hist_file_name, htmp, "_h1_energy_trg");
      append_histogram(_h1_energy_trg,htmp);
      delete htmp;
    }
  } 
  //////////////////////
  Double_t val_old;
  Double_t val_new;
  for(Int_t i = 1;i<=_h1_rate->GetNbinsX();i++){
    val_old = _h1_rate->GetBinContent(i);
    if(_h1_N_dbc_mean->GetMean()>0.0)
      val_new = val_old/_h1_N_dbc_mean->GetMean();
    else
      val_new = 0.0;
    _h1_rate->SetBinContent(i,val_new);
  }
  _h1_rate->SetBinContent(1,_h1_dbc_number_of_points_tot_w->GetBinContent(1));
  //copyHist( _h1_rate_vs_th, _h1_rate);
  copyHist( _h1_rate_vs_th, _h1_dbc_number_of_points_tot_w);
  for(Int_t i = 1;i<=_h1_rate->GetNbinsX();i++){
    //val_new = _h1_rate->Integral(i,_h1_rate->GetNbinsX());
    val_new = _h1_dbc_number_of_points_tot_w->Integral(i,_h1_rate->GetNbinsX());
    _h1_rate_vs_th->SetBinContent(i,val_new);
  }
  Double_t rate_20ep = 48976.454;
  Double_t norm = _h1_rate_vs_th->GetBinContent(1);
  for(Int_t i = 1;i<=_h1_rate_vs_th->GetNbinsX();i++){
    val_new = _h1_rate_vs_th->GetBinContent(i)/norm*rate_20ep;
    _h1_rate_vs_th->SetBinContent(i,val_new);
  }  
  //////////////////////
  //////////////////////
  if(_n_jobs > 0){
    for(Int_t i = 1;i<=__h1_energy_eff_r->GetNbinsX();i++){
      val_new = __h1_energy_eff_r->GetBinContent(i)/_n_jobs;
      __h1_energy_eff_r->SetBinContent(i,val_new);
    }
  }
  //////////////////////
  //////////////////////
  copyHist( _h1_energy_eff_r, __h1_energy_eff_r);
  Double_t n_all_ev;
  Double_t n_trg_ev;
  Double_t trf_eff_norm;
  Double_t trf_eff_norm_err;
  Double_t generation_area = TMath::Pi()*_rsimulation*_rsimulation;
  for(Int_t i = 1;i<=_h1_energy_all->GetNbinsX();i++){
    n_all_ev = _h1_energy_all->GetBinContent(i);
    n_trg_ev = _h1_energy_trg->GetBinContent(i);
    trf_eff_norm = 0.0;
    if(n_all_ev>0){
      trf_eff_norm = n_trg_ev/n_all_ev;
      trf_eff_norm_err = TMath::Sqrt(trf_eff_norm*(1.0-trf_eff_norm)/n_trg_ev);
      _h1_energy_eff_r->SetBinContent(i,trf_eff_norm*generation_area);
      _h1_energy_eff_r->SetBinError(i,trf_eff_norm_err*generation_area);
    }
  }
  TString out_file_name = _output_hist_file;
  out_file_name += ".txt";
  print_hist_to_csv(_h1_energy_eff_r, out_file_name);
  //////////////////////
}

void rateCalculator::calculate_rate_NSB(Int_t n_jobs_NSB, TString hist_file_dir){
  TString hist_file_name;
  for(Int_t i = 0;i<n_jobs_NSB;i++){
    if(get_histogram_file_name(hist_file_dir, "hist_trgA_corsika_", hist_file_name, i)){
      cout<<hist_file_name<<endl;
      //
      TH1D *htmp = new TH1D();
      get_histogram(hist_file_name, htmp, "h1_dbc_number_of_points");
      append_histogram(_h1_dbc_number_of_points_NSB,htmp);
      delete htmp;
      //
      htmp = new TH1D();
      get_histogram(hist_file_name, htmp, "h1_N_dbc");
      _h1_N_dbc_mean_NSB->Fill(htmp->GetMean());
      delete htmp;
      //
      htmp = new TH1D();
      get_histogram(hist_file_name, htmp, "h1_dbc_mean_time_ii");
      append_histogram(_h1_dbc_mean_time_ii_NSB,htmp);
      delete htmp;      
    }
  }
  Int_t fadc_MHz = 1024;
  Float_t fadc_sample_in_ns = 1000.0/fadc_MHz;
  Int_t nn_fadc_point = 75.0;
  _NSB_rate_lowth = 1.0/((nn_fadc_point - 2*2)*fadc_sample_in_ns/_h1_N_dbc_mean_NSB->GetMean()*1.0e-9);
  //
  copyHist( _h1_rate_NSB_vs_th, _h1_dbc_number_of_points_NSB);
  Double_t val_new;
  for(Int_t i = 1;i<=_h1_dbc_number_of_points_NSB->GetNbinsX();i++){
    //val_new = _h1_rate->Integral(i,_h1_rate->GetNbinsX());
    val_new = _h1_dbc_number_of_points_NSB->Integral(i,_h1_dbc_number_of_points_NSB->GetNbinsX());
    _h1_rate_NSB_vs_th->SetBinContent(i,val_new);
  }
  Double_t norm = _h1_dbc_number_of_points_NSB->Integral();
  for(Int_t i = 1;i<=_h1_rate_NSB_vs_th->GetNbinsX();i++){
    val_new = _h1_rate_NSB_vs_th->GetBinContent(i)/norm*_NSB_rate_lowth;
    _h1_rate_NSB_vs_th->SetBinContent(i,val_new);
  }  
}

//hist_trgA_corsika_0038ID.root for NSB
bool rateCalculator::get_histogram_file_name(TString  hist_file_dir, TString file_prefix, TString &hist_file_name, Int_t jobID){
  hist_file_name="";
  TSystemDirectory dire(hist_file_dir.Data(), hist_file_dir.Data()); 
  TList *files = dire.GetListOfFiles();
  TString file_suff;
  if(jobID<10){
    file_suff = "000";
    file_suff += jobID;
    file_suff += "ID.root";
  }
  else if(jobID>=10 && jobID<100){
    file_suff = "00";
    file_suff += jobID;
    file_suff += "ID.root";
  }
  if(files){
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if(fname.BeginsWith(file_prefix.Data()) &&
	 fname.EndsWith(file_suff.Data())){
	hist_file_name=hist_file_dir;
	hist_file_name+="/";
	hist_file_name+=fname.Data();
	//cout<<file_suff<<endl;
	return true;
      }
    }
  }
  return false;
}

bool rateCalculator::get_histogram_file_name(TString hist_file_prefix, Int_t i_binE, Int_t i_binTheta, Int_t i_binDist, TString &hist_file_name, Int_t jobID){
  hist_file_name="";
  TSystemDirectory dire(hist_file_prefix.Data(), hist_file_prefix.Data()); 
  TList *files = dire.GetListOfFiles(); 
  //hist_trgA_corsika_18binE_8binTheta_4binDist_1384ID.root
  //hist_trgA_corsika_0binE_0binTheta_0binDist_0656ID.root
  TString file_pref = "hist_trgA_corsika_";
  file_pref += i_binE;     file_pref += "binE_";
  file_pref += i_binTheta; file_pref += "binTheta_";
  file_pref += i_binDist;  file_pref += "binDist_";
  TString file_suff;
  if(jobID>=0){
    file_suff = "";
    file_suff += jobID;
    file_suff += "ID.root";
    //cout<<"file_suff = "<<file_suff<<endl;
  }
  else{
    file_suff = "ID.root";
  }
  if(files){
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if(fname.BeginsWith(file_pref.Data()) &&
	 fname.EndsWith(file_suff.Data())){
	hist_file_name=hist_file_prefix;
	hist_file_name+="/";
	hist_file_name+=fname.Data();
	//cout<<fname.Data()<<endl;
	return true;
      }
    }
  }
  //files->Print();
  //TIter iter(&files);
  //cout<<"files->GetSize()"<<files->GetSize()<<endl;
  //cout<<"files->GetSize()"<<files->GetSize()<<endl;
  //gSystem->AccessPathName("../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/trgA/hist_trgA_corsika_17binE_3binTheta_6binDist_1236ID.root");
  //TSystemFile file;
  //TIter next(files);
  //while ( (file = (TSystemFile)next()) ) {
  //TString fname = file->GetName();
  //cout<<"fname "<<fname<<endl;
  //if ( !(file->IsDirectory()) &&
  //(!prefix || !(*prefix) || fname.BeginsWith(prefix)) &&
  //(!suffix || !(*suffix) || ) ) {
  //cout << fname.Data() << endl;
  //}
  return false;
}

void rateCalculator::print_hist_to_csv(TH1D *h1, TString out_file_name){
  std::ofstream outfile;
  outfile.open(out_file_name.Data());
  //////////////////////
  Double_t trf_eff_norm;
  Double_t trf_eff_norm_err;
  Double_t energy_min;
  Double_t energy;
  Double_t energy_max;
  outfile<<setw(20)<<"E_min_GeV"
	 <<setw(20)<<"E_GeV"
	 <<setw(20)<<"E_max_GeV"
	 <<setw(20)<<"Area_m2"
	 <<setw(20)<<"Area_error_m2"<<endl;
  for(Int_t i = 1;i<=h1->GetNbinsX();i++){
    trf_eff_norm = h1->GetBinContent(i);
    trf_eff_norm_err = h1->GetBinError(i);
    energy_min = h1->GetBinLowEdge(i);
    energy = h1->GetBinCenter(i);
    energy_max = h1->GetBinLowEdge(i) + h1->GetBinWidth(i);
    outfile<<setw(20)<<energy_min
	   <<setw(20)<<energy
      	   <<setw(20)<<energy_max
	   <<setw(20)<<trf_eff_norm
      	   <<setw(20)<<trf_eff_norm_err<<endl;
  }
  outfile.close();
}
