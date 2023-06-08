//root
#include "TROOT.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <TStyle.h>
#include <TString.h>
#include <TCanvas.h>
#include <TTree.h>

//C, C++
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <fstream>
#include <iomanip>
#include <time.h>

const Int_t nn_max = 1000000;
const Int_t nn_fadc_point = 75;
const Int_t nn_PMT_channels = 7987;
//const Int_t nn_PMT_channels = 3000;

struct wfheaderStr {
  Int_t event_id;
  Float_t energy;
  Float_t xcore;
  Float_t ycore;
  Float_t ev_time;
  Int_t nphotons;
  Int_t n_pe;
  Int_t n_pixels;
};
std::vector<wfheaderStr> header_vec;

struct pe_info_str {
  Int_t event_id;
  std::vector<Int_t> chID;
  std::vector<Float_t> time;
};
std::vector<pe_info_str> pe_info_vec;

struct wf_str {
  Int_t wf[nn_PMT_channels][nn_fadc_point];
};
std::vector<wf_str> wf_str_vec;
  
void read_header(TString file_name);
void read_pe_info(TString file_name);
void read_wf(TString wf_file_list);

int main(int argc, char *argv[]){
  //
  clock_t start, finish;
  start = clock();
  //
  if(argc == 6 && atoi(argv[1]) == 0){
    //
    TString header_file = argv[2];
    TString pe_info_file = argv[3];
    TString wf_file_list = argv[4];
    TString outputRootFile = argv[5];
    std::cout<<std::endl
	     <<"header_file    "<<header_file<<std::endl
	     <<"pe_info_file   "<<pe_info_file<<std::endl
	     <<"wf_file_list   "<<wf_file_list<<std::endl
	     <<"outputRootFile "<<outputRootFile<<std::endl;
    ////////////////
    read_header(header_file);
    read_pe_info(pe_info_file);
    read_wf(wf_file_list);
    ////////////////
    //
    ///////////////////Root file with data/////////////////
    TFile *hfile = new TFile(outputRootFile, "RECREATE", "Simtel log data", 1);
    if (hfile->IsZombie()) {
      std::cout<<" ---> ERROR : PROBLEM with the initialization of the output ROOT file : "<<std::endl
	       <<outputRootFile
	       <<std::endl;
      assert(0);
    }
    TTree *tree = new TTree("T", "Simtel data");
    hfile->SetCompressionLevel(2);
    tree->SetAutoSave(1000000);  
    // Create new event
    TTree::SetBranchStyle(0);
    ///////////////////////////////////////////////////////
    //
    Int_t event_id;
    Float_t energy;
    Float_t xcore;
    Float_t ycore;
    Float_t ev_time;
    //
    Int_t nphotons;
    Int_t n_pe;
    Int_t n_pixels;
    //
    Int_t pe_chID[nn_max];
    Float_t pe_time[nn_max];
    //
    Int_t wf[nn_PMT_channels][nn_fadc_point];
    //
    //Event////////////////////////////////////////////////
    tree->Branch("event_id",&event_id, "event_id/I");
    tree->Branch("energy",&energy, "energy/F");
    tree->Branch("xcore", &xcore, "xcore/F");
    tree->Branch("ycore", &ycore, "ycore/F");
    tree->Branch("ev_time",&ev_time, "ev_time/F");
    tree->Branch("nphotons",&nphotons, "nphotons/I");
    tree->Branch("n_pe",&n_pe, "n_pe/I");
    tree->Branch("n_pixels",&n_pixels, "n_pixels/I");
    //
    tree->Branch("pe_chID",pe_chID, "pe_chID[n_pe]/I");
    tree->Branch("pe_time",pe_time, "pe_time[n_pe]/F");
    //
    tree->Branch("wf", wf, "wf[7987][75]/I");
    //
    ///////////////////////////////////////////////////////
    //
    std::cout<<"header_vec.size()  "<<header_vec.size()<<std::endl
	     <<"pe_info_vec.size() "<<pe_info_vec.size()<<std::endl
    	     <<"wf_str_vec.size()  "<<wf_str_vec.size()<<std::endl;
    if(header_vec.size() != pe_info_vec.size()){
      std::cout<<" --> ERROR : header_vec.size() != pe_info_vec.size() "<<std::endl
	       <<"             header_vec.size()  = "<<header_vec.size()<<std::endl
	       <<"            pe_info_vec.size()  = "<<pe_info_vec.size()<<std::endl;
      assert(0);
    }     
    //
    for(unsigned int i = 0; i < header_vec.size(); i++){
      if(i%10000 == 0)
        std::cout<<i<<std::endl;
      //
      event_id = header_vec.at(i).event_id;
      energy = header_vec.at(i).energy;
      xcore = header_vec.at(i).xcore;
      ycore = header_vec.at(i).ycore;
      ev_time = header_vec.at(i).ev_time;
      nphotons = header_vec.at(i).nphotons;
      n_pe = header_vec.at(i).n_pe;
      n_pixels = header_vec.at(i).n_pixels;
      //
      if(pe_info_vec.at(i).event_id != event_id){
	std::cout<<" --> ERROR : pe_info_vec.at(i).event_id != event_id "<<std::endl
		 <<"             pe_info_vec.at(i).event_id  = "<<pe_info_vec.at(i).event_id<<std::endl
		 <<"                               event_id  = "<<event_id<<std::endl;
	assert(0);
      }
      //
      if(pe_info_vec.at(i).chID.size() != (unsigned int)n_pe ||
	 pe_info_vec.at(i).time.size() != (unsigned int)n_pe){
	std::cout<<" --> ERROR : pe_info_vec.at(i).chID.size() != (unsigned int)n_pe || "<<std::endl
		 <<"             pe_info_vec.at(i).time.size() != (unsigned int)n_pe"<<std::endl
		 <<"             pe_info_vec.at(i).chID.size()  = "<<pe_info_vec.at(i).chID.size()<<std::endl
		 <<"             pe_info_vec.at(i).time.size()  = "<<pe_info_vec.at(i).time.size()<<std::endl
		 <<"                                      n_pe  = "<<n_pe<<std::endl;
	assert(0);
      }
      for(int j = 0; j < nn_max; j++){
	pe_chID[j] = 0;
	pe_time[j] = 0.0;
      }
      ////////////////////
      for(Int_t ii = 0; ii<nn_PMT_channels;ii++)
	for(Int_t jj = 0; jj<nn_fadc_point;jj++)
	  wf[ii][jj] = wf_str_vec.at(i).wf[ii][jj];
      ////////////////////
      if(n_pe<nn_max){
	for(unsigned int j = 0; j < pe_info_vec.at(i).chID.size(); j++){
	  pe_chID[j] = pe_info_vec.at(i).chID.at(j);
	  pe_time[j] = pe_info_vec.at(i).time.at(j);
	}
	//
	tree->Fill();
      }
    }
    hfile = tree->GetCurrentFile();
    hfile->Write();
    hfile->Close();
  }
  else if(argc == 6 && atoi(argv[1]) == 1){
  }
  else{
    std::cout<<"  runID [1] = 0        "<<std::endl
	     <<"        [2] - header_file "<<std::endl
	     <<"        [2] - pe_info_file "<<std::endl
	     <<"        [3] - outputRootFile "<<std::endl;
    std::cout<<"  runID [1] = 1        "<<std::endl
	     <<"        [2] - path     "<<std::endl
	     <<"        [3] - prefix   "<<std::endl
    	     <<"        [4] - nfiles   "<<std::endl;
  }  //
  finish = clock();
  std::cout<<"-------------------------"<<std::endl
  	   <<"Working time : "<<((finish - start)/CLOCKS_PER_SEC)<<" (sec)"<<std::endl
  	   <<"-------------------------"<<std::endl;
  //
  return 0;
}

void read_header(TString file_name){
  //
  std::ifstream fFile(file_name.Data());
  std::cout<<file_name<<std::endl;
  //
  std::string mot;
  Float_t event_id;
  Float_t energy;
  Float_t xcore;
  Float_t ycore;
  Float_t ev_time;
  Float_t nphotons;
  Float_t n_pe;
  Float_t n_pixels;
  //
  if(fFile.is_open()){
    while(fFile>>mot>>
	  event_id>>
	  energy>>
	  xcore>>
	  ycore>>
	  ev_time>>
	  nphotons>>
	  n_pe>>
	  n_pixels){
      wfheaderStr tmp;      
      tmp.event_id=event_id;
      tmp.energy=energy;
      tmp.xcore=xcore;
      tmp.ycore=ycore;
      tmp.ev_time=ev_time;
      tmp.nphotons=nphotons;
      tmp.n_pe=n_pe;
      tmp.n_pixels=n_pixels;
      //std::cout<<tmp.event_id<<std::endl;
      header_vec.push_back(tmp);
    }
    fFile.close();
  }  
}

void read_pe_info(TString file_name){
  //
  std::ifstream fFile(file_name.Data());
  std::cout<<file_name<<std::endl;
  //
  Int_t verbosity = 0;
  //
  Int_t id;
  Float_t event_id;
  Float_t ch_id;
  Float_t pe_time;
  //
  if(fFile.is_open()){
    for(unsigned int i = 0; i < header_vec.size(); i++){
      pe_info_str pe_info;
      pe_info.event_id = header_vec.at(i).event_id;
      for(unsigned int j = 0; j < (unsigned int)header_vec.at(i).n_pe; j++){
	fFile>>id>>event_id>>ch_id>>pe_time;
	if(verbosity>1){
	  if(j%1000000==0)
	    std::cout<<std::setw(10)<<"id       "<<std::setw(10)<<id
		     <<std::setw(10)<<"event_id "<<std::setw(10)<<event_id
		     <<std::setw(10)<<"ch_id    "<<std::setw(10)<<ch_id
		     <<std::setw(10)<<"pe_time  "<<std::setw(10)<<pe_time<<std::endl;
	}
	if((Int_t)event_id != header_vec.at(i).event_id){
	  std::cout<<" ERROR ---> (Int_t)event_id != header_vec.at(i).event_id "<<std::endl
		   <<"            (Int_t)event_id "<<(Int_t)event_id<<std::endl
		   <<"  header_vec.at(i).event_id "<<header_vec.at(i).event_id<<std::endl;
	  assert(0);
	}
	pe_info.chID.push_back((Int_t)ch_id);
	pe_info.time.push_back(pe_time);
      }
      pe_info_vec.push_back(pe_info);
    }
    fFile.close();
  }
}

void read_wf(TString wf_file_list){
  std::ifstream indataList;
  std::ifstream indataFile;
  indataList.open(wf_file_list.Data());
  std::string mot;
  TString inputDataFile;
  Int_t nEv = 0;
  Int_t vVal;
  //
  while (indataList >> mot >> nEv ){
    inputDataFile = mot;
    std::cout<<" ---> Conversion of "<<inputDataFile<<std::endl;
    indataFile.open(inputDataFile.Data());
    for( Int_t i_ev = 0; i_ev < nEv; i_ev++){
      wf_str wf;
      for(Int_t i_ch = 0;i_ch<nn_PMT_channels; i_ch++){
	//indataFile>>vVal;
	for(Int_t i_adc = 0;i_adc<nn_fadc_point; i_adc++){
	  indataFile>>vVal;
	  wf.wf[i_ch][i_adc] = vVal;
	}
      }
      wf_str_vec.push_back(wf);
    }
    indataFile.close();
  }
  indataList.close();
}
