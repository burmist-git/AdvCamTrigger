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

using namespace std;

rateCalculator::rateCalculator():
  _name(""),
  _title(""),
  _hist_file_prefix(""),
  _output_hist_file(""),
  _evH_flux(NULL),
  _h1_rate(NULL)
{
}

rateCalculator::rateCalculator( const char* name, const char* title, TString hist_file_prefix, TString output_hist_file, evstHist* evH_flux):
  _name(name),
  _title(title),
  _hist_file_prefix(hist_file_prefix),
  _output_hist_file(output_hist_file),
  _evH_flux(evH_flux),
  _h1_rate(new TH1D)
{
  TString h1_rate_name = "_h1_rate";
  h1_rate_name += "_";
  h1_rate_name += _name;
  TString h1_rate_title = "_h1_rate";
  h1_rate_title += "_";
  h1_rate_title += _title;
  _h1_rate->SetNameTitle(h1_rate_name.Data(),h1_rate_title.Data());
  calculate_rate();
  save_output();
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
  for(unsigned int i = 0;i<_h1_rates_per_bin_v.size();i++)
    _h1_rates_per_bin_v.at(i)->Write();
  //
  _h1_rate->Write();  
  //
  rootFile->Close();  
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
}

void rateCalculator::addHist(TH1D *h1, TH1D *h1_to_add){
  for(Int_t i = 1;i<=h1_to_add->GetNbinsX();i++)
    h1->SetBinContent(i,h1_to_add->GetBinContent(i)+h1->GetBinContent(i));
}

void rateCalculator::get_histogram(TString hist_file_name, TH1D *h1){
  TFile *f1 = new TFile(hist_file_name.Data());
  TH1D *h1_tmp = (TH1D*)f1->Get("h1_dbc_number_of_points_w");
  copyHist(h1,h1_tmp);
  _h1_rates_per_bin_v.push_back(h1);
  f1->Close();
}

void rateCalculator::append_histogram(TH1D *h1){
  if(_h1_rate->GetNbinsX()==1)
    copyHist(_h1_rate,h1);
  else
    addHist(_h1_rate,h1);
}

void rateCalculator::calculate_rate(){
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
	  get_histogram(hist_file_name, h1);
	  append_histogram(h1);
	}
      }
    }
  } 
  //if(get_histogram_file_name(_hist_file_prefix, 5, 0, 0, hist_file_name))
  //cout<<"hist_file_name = "<<hist_file_name<<endl;
  //if(get_histogram_file_name(_hist_file_prefix, 15, 3, 9, hist_file_name))
  //cout<<"hist_file_name = "<<hist_file_name<<endl;
}

bool rateCalculator::get_histogram_file_name(TString hist_file_prefix, Int_t i_binE, Int_t i_binTheta, Int_t i_binDist, TString &hist_file_name){
  hist_file_name="";
  TSystemDirectory dire(hist_file_prefix.Data(), hist_file_prefix.Data()); 
  TList *files = dire.GetListOfFiles(); 
  //hist_trgA_corsika_18binE_8binTheta_4binDist_1384ID.root
  TString file_pref = "hist_trgA_corsika_";
  file_pref += i_binE;     file_pref += "binE_";
  file_pref += i_binTheta; file_pref += "binTheta_";
  file_pref += i_binDist;  file_pref += "binDist_";
  TString file_suff = "ID.root";
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
	return true;
	cout<<fname.Data()<<endl;
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

rateCalculator::~rateCalculator(){
}
