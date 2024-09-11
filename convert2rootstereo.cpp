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
#include <vector>

const Int_t nn_max = 1000000;
const Int_t nn_fadc_point = 75;
const Int_t nn_PMT_channels = 7987;

struct wfheaderStr {
  Int_t event_id;
  Float_t energy;
  Float_t azimuth;
  Float_t altitude;
  Float_t h_first_int;
  Float_t xmax;
  Float_t hmax;
  Float_t emax;
  Float_t cmax;
  Float_t xcore;
  Float_t ycore;
  Float_t ev_time_LST1;
  Float_t ev_time_LST2;
  Float_t ev_time_LST3;
  Float_t ev_time_LST4;  
  Int_t nphotons_LST1;
  Int_t nphotons_LST2;
  Int_t nphotons_LST3;
  Int_t nphotons_LST4;  
  Int_t n_pe_LST1;
  Int_t n_pe_LST2;
  Int_t n_pe_LST3;
  Int_t n_pe_LST4;
  Int_t n_pixels_LST1;
  Int_t n_pixels_LST2;
  Int_t n_pixels_LST3;
  Int_t n_pixels_LST4;
};
std::vector<wfheaderStr> header_vec;

struct pe_info_str {
  Int_t event_id;
  std::vector<Int_t> chID;
  std::vector<Float_t> time;
};
std::vector<pe_info_str> pe_info_vec_LST1;
std::vector<pe_info_str> pe_info_vec_LST2;
std::vector<pe_info_str> pe_info_vec_LST3;
std::vector<pe_info_str> pe_info_vec_LST4;
  
void read_header(TString file_name, std::vector<wfheaderStr> &header_v);
void read_pe_info(TString file_name, const std::vector<wfheaderStr> &header_v, Int_t LST_ID, std::vector<pe_info_str> &pe_info_vec);
void convert2root(const std::vector<TString> &header_file_v,
		  const std::vector<TString> &pe_info_file_v_LST1,
		  const std::vector<TString> &pe_info_file_v_LST2,
		  const std::vector<TString> &pe_info_file_v_LST3,
		  const std::vector<TString> &pe_info_file_v_LST4,
		  TString outputRootFile);

int main(int argc, char *argv[]){
  //
  clock_t start, finish;
  start = clock();
  //
  if(argc == 8 && atoi(argv[1]) == 0){
    TString header_file = argv[2];
    TString pe_info_file_LST1 = argv[3];
    TString pe_info_file_LST2 = argv[4];
    TString pe_info_file_LST3 = argv[5];
    TString pe_info_file_LST4 = argv[6];
    TString outputRootFile = argv[7];
    std::cout<<std::endl
	     <<"header_file       "<<header_file<<std::endl
	     <<"pe_info_file_LST1 "<<pe_info_file_LST1<<std::endl
      	     <<"pe_info_file_LST2 "<<pe_info_file_LST2<<std::endl
      	     <<"pe_info_file_LST3 "<<pe_info_file_LST3<<std::endl
      	     <<"pe_info_file_LST4 "<<pe_info_file_LST4<<std::endl
	     <<"outputRootFile    "<<outputRootFile<<std::endl;
    //
    std::vector<TString> header_file_v;
    std::vector<TString> pe_info_file_v_LST1;
    std::vector<TString> pe_info_file_v_LST2;
    std::vector<TString> pe_info_file_v_LST3;
    std::vector<TString> pe_info_file_v_LST4;
    header_file_v.push_back(header_file);
    pe_info_file_v_LST1.push_back(pe_info_file_LST1);
    pe_info_file_v_LST2.push_back(pe_info_file_LST2);
    pe_info_file_v_LST3.push_back(pe_info_file_LST3);
    pe_info_file_v_LST4.push_back(pe_info_file_LST4);
    convert2root(header_file_v,
		 pe_info_file_v_LST1,
		 pe_info_file_v_LST2,
		 pe_info_file_v_LST3,
		 pe_info_file_v_LST4,		 
		 outputRootFile);
  }
  else{
    std::cout<<"  runID [1] = 0        "<<std::endl
	     <<"        [2] - header_file       "<<std::endl
	     <<"        [3] - pe_info_file_LST1 "<<std::endl
	     <<"        [3] - pe_info_file_LST2 "<<std::endl
      	     <<"        [3] - pe_info_file_LST3 "<<std::endl
      	     <<"        [3] - pe_info_file_LST4 "<<std::endl
	     <<"        [5] - outputRootFile    "<<std::endl;
  }  //
  finish = clock();
  std::cout<<"-------------------------"<<std::endl
  	   <<"Working time : "<<((finish - start)/CLOCKS_PER_SEC)<<" (sec)"<<std::endl
  	   <<"-------------------------"<<std::endl;
  //
  return 0;
}

void read_header( TString file_name, std::vector<wfheaderStr> &header_v){
  //
  std::ifstream fFile(file_name.Data());
  std::cout<<file_name<<std::endl;
  //
  std::string mot;       //0
  Float_t event_id;      //1
  Float_t energy;        //2
  Float_t azimuth;       //3
  Float_t altitude;      //4
  Float_t h_first_int;   //5
  Float_t xmax;          //6
  Float_t hmax;          //7
  Float_t emax;          //8
  Float_t cmax;          //9
  Float_t xcore;         //10
  Float_t ycore;         //11
  //
  Float_t ev_time_LST1;  //12
  Float_t ev_time_LST2;  //13
  Float_t ev_time_LST3;  //14
  Float_t ev_time_LST4;  //15
  //
  Float_t nphotons_LST1; //16
  Float_t nphotons_LST2; //17
  Float_t nphotons_LST3; //18
  Float_t nphotons_LST4; //19
  //  
  Float_t n_pe_LST1;     //20
  Float_t n_pe_LST2;     //21
  Float_t n_pe_LST3;     //22
  Float_t n_pe_LST4;     //23
  //
  Float_t n_pixels_LST1; //24
  Float_t n_pixels_LST2; //25
  Float_t n_pixels_LST3; //26
  Float_t n_pixels_LST4; //27
  //
  //0  0
  //1  2707.0
  //2  0.1709720939397812
  //3  3.392895221710205
  //4  1.3400671482086182
  //5  14761.0888671875
  //6  542.0
  //7  5353.32568359375
  //8  547.1428833007812
  //9  564.1544189453125
  //10 -504.9872741699219
  //11 -129.00323486328125
  //12 0.0
  //13 0.0
  //14 603.0427856445312
  //15 0.0
  //16 0.0
  //17 0.0
  //18 47.0
  //19 0.0
  //20 0.0
  //21 0.0
  //22 20.0
  //23 0.0
  //24 0.0
  //25 0.0
  //26 15.0
  //27 0.0
  //
  //
  if(fFile.is_open()){
    while(fFile>>mot>>
	  event_id>>
	  energy>>
	  azimuth>>
	  altitude>>
	  h_first_int>>
	  xmax>>
	  hmax>>
	  emax>>
	  cmax>>
	  xcore>>
	  ycore>>
	  ev_time_LST1>>
	  ev_time_LST2>>
	  ev_time_LST3>>
	  ev_time_LST4>>	  
	  nphotons_LST1>>
	  nphotons_LST2>>
	  nphotons_LST3>>
	  nphotons_LST4>>	  
	  n_pe_LST1>>
	  n_pe_LST2>>
	  n_pe_LST3>>
	  n_pe_LST4>>	  
	  n_pixels_LST1>>
	  n_pixels_LST2>>
	  n_pixels_LST3>>
	  n_pixels_LST4){
      wfheaderStr tmp;      
      tmp.event_id=event_id;
      tmp.energy=energy;
      tmp.azimuth=azimuth;
      tmp.altitude=altitude;
      tmp.h_first_int=h_first_int;
      tmp.xmax=xmax;
      tmp.hmax=hmax;
      tmp.emax=emax;
      tmp.cmax=cmax;
      tmp.xcore=xcore;
      tmp.ycore=ycore;
      tmp.ev_time_LST1=ev_time_LST1;
      tmp.ev_time_LST2=ev_time_LST2;
      tmp.ev_time_LST3=ev_time_LST3;
      tmp.ev_time_LST4=ev_time_LST4;
      tmp.nphotons_LST1=nphotons_LST1;
      tmp.nphotons_LST2=nphotons_LST2;
      tmp.nphotons_LST3=nphotons_LST3;
      tmp.nphotons_LST4=nphotons_LST4;
      tmp.n_pe_LST1=n_pe_LST1;
      tmp.n_pe_LST2=n_pe_LST2;
      tmp.n_pe_LST3=n_pe_LST3;
      tmp.n_pe_LST4=n_pe_LST4;
      tmp.n_pixels_LST1=n_pixels_LST1;
      tmp.n_pixels_LST2=n_pixels_LST2;
      tmp.n_pixels_LST3=n_pixels_LST3;
      tmp.n_pixels_LST4=n_pixels_LST4;
      //std::cout<<tmp.event_id<<std::endl;
      //assert(0);
      header_v.push_back(tmp);
    }
    fFile.close();
  }  
}

void read_pe_info( TString file_name, const std::vector<wfheaderStr> &header_v, Int_t LST_ID, std::vector<pe_info_str> &pe_info_vec){
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
  unsigned int npe_lst = 0;
  //
  if(fFile.is_open()){
    for(unsigned int i = 0; i < header_v.size(); i++){
      pe_info_str pe_info;
      pe_info.event_id = header_v.at(i).event_id;
      //
      //
      if(LST_ID == 1)
	npe_lst = (unsigned int)header_v.at(i).n_pe_LST1;
      else if (LST_ID == 2)
	npe_lst = (unsigned int)header_v.at(i).n_pe_LST2;
      else if (LST_ID == 3)
	npe_lst = (unsigned int)header_v.at(i).n_pe_LST3;
      else if (LST_ID == 4)
	npe_lst = (unsigned int)header_v.at(i).n_pe_LST4;
      else
	assert(0);
      //
      //
      for(unsigned int j = 0; j < npe_lst; j++){
	fFile>>id>>event_id>>ch_id>>pe_time;
	//std::cout<<std::setw(10)<<"id       "<<std::setw(10)<<id
	//	 <<std::setw(10)<<"event_id "<<std::setw(10)<<event_id
	//	 <<std::setw(10)<<"ch_id    "<<std::setw(10)<<ch_id
	//	 <<std::setw(10)<<"pe_time  "<<std::setw(10)<<pe_time<<std::endl;
	if(verbosity>1){
	  if(j%1000000==0)
	    std::cout<<std::setw(10)<<"id       "<<std::setw(10)<<id
		     <<std::setw(10)<<"event_id "<<std::setw(10)<<event_id
		     <<std::setw(10)<<"ch_id    "<<std::setw(10)<<ch_id
		     <<std::setw(10)<<"pe_time  "<<std::setw(10)<<pe_time<<std::endl;
	}
	if((Int_t)event_id != header_v.at(i).event_id){
	  std::cout<<" ERROR ---> (Int_t)event_id != header_v.at(i).event_id "<<std::endl
		   <<"            (Int_t)event_id "<<(Int_t)event_id<<std::endl
		   <<"  header_vec.at(i).event_id "<<header_v.at(i).event_id<<std::endl
		   <<"  header_vec.at(0).event_id "<<header_v.at(0).event_id<<std::endl;
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

void convert2root(const std::vector<TString> &header_file_v,
		  const std::vector<TString> &pe_info_file_v_LST1,
		  const std::vector<TString> &pe_info_file_v_LST2,
		  const std::vector<TString> &pe_info_file_v_LST3,
		  const std::vector<TString> &pe_info_file_v_LST4,
		  TString outputRootFile) {
  //
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
  unsigned int npe_current = 0;
  //
  Int_t event_id;
  Float_t energy;
  Float_t azimuth;
  Float_t altitude;
  Float_t h_first_int;
  Float_t xmax;
  Float_t hmax;
  Float_t emax;
  Float_t cmax;
  Float_t xcore;
  Float_t ycore;
  Float_t ev_time;
  //
  Int_t nphotons_LST1;
  Int_t nphotons_LST2;
  Int_t nphotons_LST3;
  Int_t nphotons_LST4;
  //
  Int_t n_pe_LST1;
  Int_t n_pe_LST2;
  Int_t n_pe_LST3;
  Int_t n_pe_LST4;
  //
  Int_t n_pixels_LST1;
  Int_t n_pixels_LST2;
  Int_t n_pixels_LST3;
  Int_t n_pixels_LST4;
  //
  Int_t pe_chID_LST1[nn_max];
  Int_t pe_chID_LST2[nn_max];
  Int_t pe_chID_LST3[nn_max];
  Int_t pe_chID_LST4[nn_max];  
  //
  Float_t pe_time_LST1[nn_max];
  Float_t pe_time_LST2[nn_max];
  Float_t pe_time_LST3[nn_max];
  Float_t pe_time_LST4[nn_max];
  //
  //Event////////////////////////////////////////////////
  tree->Branch("event_id",&event_id, "event_id/I");
  tree->Branch("energy",&energy, "energy/F");
  tree->Branch("azimuth",&azimuth, "azimuth/F");
  tree->Branch("altitude",&altitude, "altitude/F");
  tree->Branch("h_first_int",&h_first_int, "h_first_int/F");
  tree->Branch("xmax",&xmax, "xmax/F");
  tree->Branch("hmax",&hmax, "hmax/F");
  tree->Branch("emax",&emax, "emax/F");
  tree->Branch("cmax",&cmax, "cmax/F");
  tree->Branch("xcore", &xcore, "xcore/F");
  tree->Branch("ycore", &ycore, "ycore/F");
  //
  tree->Branch("ev_time_LST1",&ev_time_LST1,"ev_time_LST1/F");
  tree->Branch("ev_time_LST2",&ev_time_LST2,"ev_time_LST2/F");
  tree->Branch("ev_time_LST3",&ev_time_LST3,"ev_time_LST3/F");
  tree->Branch("ev_time_LST4",&ev_time_LST4,"ev_time_LST4/F");  
  //
  tree->Branch("nphotons_LST1",&nphotons_LST1, "nphotons_LST1/I");
  tree->Branch("nphotons_LST2",&nphotons_LST2, "nphotons_LST2/I");
  tree->Branch("nphotons_LST3",&nphotons_LST3, "nphotons_LST3/I");
  tree->Branch("nphotons_LST4",&nphotons_LST4, "nphotons_LST4/I");  
  //
  tree->Branch("n_pe_LST1",&n_pe_LST1, "n_pe_LST1/I");
  tree->Branch("n_pe_LST2",&n_pe_LST2, "n_pe_LST2/I");
  tree->Branch("n_pe_LST3",&n_pe_LST3, "n_pe_LST3/I");
  tree->Branch("n_pe_LST4",&n_pe_LST4, "n_pe_LST4/I");
  //
  tree->Branch("n_pixels_LST1",&n_pixels_LST1, "n_pixels_LST1/I");
  tree->Branch("n_pixels_LST2",&n_pixels_LST2, "n_pixels_LST2/I");
  tree->Branch("n_pixels_LST3",&n_pixels_LST3, "n_pixels_LST3/I");
  tree->Branch("n_pixels_LST4",&n_pixels_LST4, "n_pixels_LST4/I");
  //
  tree->Branch("pe_chID_LST1",pe_chID_LST1, "pe_chID_LST1[n_pe_LST1]/I");
  tree->Branch("pe_chID_LST2",pe_chID_LST2, "pe_chID_LST2[n_pe_LST2]/I");
  tree->Branch("pe_chID_LST3",pe_chID_LST3, "pe_chID_LST3[n_pe_LST3]/I");
  tree->Branch("pe_chID_LST4",pe_chID_LST4, "pe_chID_LST4[n_pe_LST4]/I");
  //
  tree->Branch("pe_time_LST1",pe_time_LST1, "pe_time_LST1[n_pe_LST1]/F");
  tree->Branch("pe_time_LST2",pe_time_LST2, "pe_time_LST2[n_pe_LST2]/F");
  tree->Branch("pe_time_LST3",pe_time_LST3, "pe_time_LST3[n_pe_LST3]/F");
  tree->Branch("pe_time_LST4",pe_time_LST4, "pe_time_LST4[n_pe_LST4]/F");
  //
  ///////////////////////////////////////////////////////
  //
  ////////////////
  for(unsigned int j = 0; j<header_file_v.size();j++){
    //
    header_vec.clear();
    pe_info_vec_LST1.clear();
    pe_info_vec_LST2.clear();
    pe_info_vec_LST3.clear();
    pe_info_vec_LST4.clear();
    //
    std::cout<<std::endl
	     <<"header_file       "<<header_file_v.at(j)<<std::endl
	     <<"pe_info_file_LST1 "<<pe_info_file_v_LST1.at(j)<<std::endl
      	     <<"pe_info_file_LST2 "<<pe_info_file_v_LST2.at(j)<<std::endl
	     <<"pe_info_file_LST3 "<<pe_info_file_v_LST3.at(j)<<std::endl
      	     <<"pe_info_file_LST4 "<<pe_info_file_v_LST4.at(j)<<std::endl;
    //
    read_header(header_file_v.at(j), header_vec);
    read_pe_info(pe_info_file_v_LST1.at(j), header_vec, 1, pe_info_vec_LST1);
    read_pe_info(pe_info_file_v_LST2.at(j), header_vec, 2, pe_info_vec_LST2);
    read_pe_info(pe_info_file_v_LST3.at(j), header_vec, 3, pe_info_vec_LST3);
    read_pe_info(pe_info_file_v_LST4.at(j), header_vec, 4, pe_info_vec_LST4);
    //    
    std::cout<<"header_vec.size()  "<<header_vec.size()<<std::endl
	     <<"pe_info_vec_LST1.size() "<<pe_info_vec_LST1.size()<<std::endl
      	     <<"pe_info_vec_LST2.size() "<<pe_info_vec_LST2.size()<<std::endl
      	     <<"pe_info_vec_LST3.size() "<<pe_info_vec_LST3.size()<<std::endl
	     <<"pe_info_vec_LST4.size() "<<pe_info_vec_LST4.size()<<std::endl;
    //
    for(unsigned int i = 0; i < header_vec.size(); i++){
      if(i%10000 == 0)
	std::cout<<i<<std::endl;
      //
      event_id = header_vec.at(i).event_id;
      energy = header_vec.at(i).energy;
      azimuth = header_vec.at(i).azimuth;
      altitude = header_vec.at(i).altitude;
      h_first_int = header_vec.at(i).h_first_int;
      xmax = header_vec.at(i).xmax;
      hmax = header_vec.at(i).hmax;
      emax = header_vec.at(i).emax;
      cmax = header_vec.at(i).cmax;
      xcore = header_vec.at(i).xcore;
      ycore = header_vec.at(i).ycore;
      //
      ev_time_LST1 = header_vec.at(i).ev_time_LST1;
      ev_time_LST2 = header_vec.at(i).ev_time_LST2;
      ev_time_LST3 = header_vec.at(i).ev_time_LST3;
      ev_time_LST4 = header_vec.at(i).ev_time_LST4;
      //
      nphotons_LST1 = header_vec.at(i).nphotons_LST1;
      nphotons_LST2 = header_vec.at(i).nphotons_LST2;
      nphotons_LST3 = header_vec.at(i).nphotons_LST3;
      nphotons_LST4 = header_vec.at(i).nphotons_LST4;
      //
      n_pe_LST1 = header_vec.at(i).n_pe_LST1;
      n_pe_LST2 = header_vec.at(i).n_pe_LST2;
      n_pe_LST3 = header_vec.at(i).n_pe_LST3;
      n_pe_LST4 = header_vec.at(i).n_pe_LST4;
      //
      n_pixels_LST1 = header_vec.at(i).n_pixels_LST1;
      n_pixels_LST2 = header_vec.at(i).n_pixels_LST2;
      n_pixels_LST3 = header_vec.at(i).n_pixels_LST3;
      n_pixels_LST4 = header_vec.at(i).n_pixels_LST4;
      //
      for(int j = 0; j < nn_max; j++){
	//
	pe_chID_LST1[j] = 0;
	pe_time_LST1[j] = 0.0;
	//
	pe_chID_LST2[j] = 0;
	pe_time_LST2[j] = 0.0;
	//
	pe_chID_LST3[j] = 0;
	pe_time_LST3[j] = 0.0;
	//
	pe_chID_LST4[j] = 0;
	pe_time_LST4[j] = 0.0;
	//
      }
      ////////////////////
      /*
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
      */
      ////////////////////
      //LST1
      if(n_pe_LST1>0){
	npe_current = pe_info_vec_LST1.at(i).chID.size();
	//
	if(pe_info_vec_LST1.at(i).event_id != event_id){
	  std::cout<<" --> ERROR : pe_info_vec_LST1.at(i).event_id != event_id "<<std::endl
		   <<"             pe_info_vec_LST1.at(i).event_id  = "<<pe_info_vec_LST1.at(i).event_id<<std::endl
		   <<"                                    event_id  = "<<event_id<<std::endl;
	  assert(0);
	}
	//
	if(pe_info_vec_LST1.at(i).chID.size() != (unsigned int)n_pe_LST1 ||
	   pe_info_vec_LST1.at(i).time.size() != (unsigned int)n_pe_LST1){
	  std::cout<<" --> ERROR : pe_info_vec_LST1.at(i).chID.size() != (unsigned int)n_pe_LST1 || "<<std::endl
		   <<"             pe_info_vec_LST1.at(i).time.size() != (unsigned int)n_pe_LST1"<<std::endl
		   <<"             pe_info_vec_LST1.at(i).chID.size()  = "<<pe_info_vec_LST1.at(i).chID.size()<<std::endl
		   <<"             pe_info_vec_LST1.at(i).time.size()  = "<<pe_info_vec_LST1.at(i).time.size()<<std::endl
		   <<"                                           n_pe  = "<<n_pe_LST1<<std::endl;
	  assert(0);
	}
	if(npe_current<nn_max)
	  npe_current = nn_max;
	for(unsigned int j = 0; j < npe_current; j++){
	  pe_chID_LST1[j] = pe_info_vec_LST1.at(i).chID.at(j);
	  pe_time_LST1[j] = pe_info_vec_LST1.at(i).time.at(j);
	}
      }
      //LST2
      if(n_pe_LST2>0){
	npe_current = pe_info_vec_LST2.at(i).chID.size();
	//
	if(pe_info_vec_LST2.at(i).event_id != event_id){
	  std::cout<<" --> ERROR : pe_info_vec_LST2.at(i).event_id != event_id "<<std::endl
		   <<"             pe_info_vec_LST2.at(i).event_id  = "<<pe_info_vec_LST2.at(i).event_id<<std::endl
		   <<"                                    event_id  = "<<event_id<<std::endl;
	  assert(0);
	}
	//
	if(pe_info_vec_LST2.at(i).chID.size() != (unsigned int)n_pe_LST2 ||
	   pe_info_vec_LST2.at(i).time.size() != (unsigned int)n_pe_LST2){
	  std::cout<<" --> ERROR : pe_info_vec_LST2.at(i).chID.size() != (unsigned int)n_pe_LST2 || "<<std::endl
		   <<"             pe_info_vec_LST2.at(i).time.size() != (unsigned int)n_pe_LST2"<<std::endl
		   <<"             pe_info_vec_LST2.at(i).chID.size()  = "<<pe_info_vec_LST2.at(i).chID.size()<<std::endl
		   <<"             pe_info_vec_LST2.at(i).time.size()  = "<<pe_info_vec_LST2.at(i).time.size()<<std::endl
		   <<"                                           n_pe  = "<<n_pe_LST2<<std::endl;
	  assert(0);
	}
	if(npe_current<nn_max)
	  npe_current = nn_max;
	for(unsigned int j = 0; j < npe_current; j++){
	  pe_chID_LST2[j] = pe_info_vec_LST2.at(i).chID.at(j);
	  pe_time_LST2[j] = pe_info_vec_LST2.at(i).time.at(j);
	}
      }
      //LST3
      if(n_pe_LST3>0){
	npe_current = pe_info_vec_LST3.at(i).chID.size();
	//
	if(pe_info_vec_LST3.at(i).event_id != event_id){
	  std::cout<<" --> ERROR : pe_info_vec_LST3.at(i).event_id != event_id "<<std::endl
		   <<"             pe_info_vec_LST3.at(i).event_id  = "<<pe_info_vec_LST3.at(i).event_id<<std::endl
		   <<"                                    event_id  = "<<event_id<<std::endl;
	  assert(0);
	}
	//
	if(pe_info_vec_LST3.at(i).chID.size() != (unsigned int)n_pe_LST3 ||
	   pe_info_vec_LST3.at(i).time.size() != (unsigned int)n_pe_LST3){
	  std::cout<<" --> ERROR : pe_info_vec_LST3.at(i).chID.size() != (unsigned int)n_pe_LST3 || "<<std::endl
		   <<"             pe_info_vec_LST3.at(i).time.size() != (unsigned int)n_pe_LST3"<<std::endl
		   <<"             pe_info_vec_LST3.at(i).chID.size()  = "<<pe_info_vec_LST3.at(i).chID.size()<<std::endl
		   <<"             pe_info_vec_LST3.at(i).time.size()  = "<<pe_info_vec_LST3.at(i).time.size()<<std::endl
		   <<"                                           n_pe  = "<<n_pe_LST3<<std::endl;
	  assert(0);
	}
	if(npe_current<nn_max)
	  npe_current = nn_max;
	for(unsigned int j = 0; j < npe_current; j++){
	  pe_chID_LST3[j] = pe_info_vec_LST3.at(i).chID.at(j);
	  pe_time_LST3[j] = pe_info_vec_LST3.at(i).time.at(j);
	}
      }
      //LST4
      if(n_pe_LST4>0){
	npe_current = pe_info_vec_LST4.at(i).chID.size();
	//
	if(pe_info_vec_LST4.at(i).event_id != event_id){
	  std::cout<<" --> ERROR : pe_info_vec_LST4.at(i).event_id != event_id "<<std::endl
		   <<"             pe_info_vec_LST4.at(i).event_id  = "<<pe_info_vec_LST4.at(i).event_id<<std::endl
		   <<"                                    event_id  = "<<event_id<<std::endl;
	  assert(0);
	}
	//
	if(pe_info_vec_LST4.at(i).chID.size() != (unsigned int)n_pe_LST4 ||
	   pe_info_vec_LST4.at(i).time.size() != (unsigned int)n_pe_LST4){
	  std::cout<<" --> ERROR : pe_info_vec_LST4.at(i).chID.size() != (unsigned int)n_pe_LST4 || "<<std::endl
		   <<"             pe_info_vec_LST4.at(i).time.size() != (unsigned int)n_pe_LST4"<<std::endl
		   <<"             pe_info_vec_LST4.at(i).chID.size()  = "<<pe_info_vec_LST4.at(i).chID.size()<<std::endl
		   <<"             pe_info_vec_LST4.at(i).time.size()  = "<<pe_info_vec_LST4.at(i).time.size()<<std::endl
		   <<"                                           n_pe  = "<<n_pe_LST4<<std::endl;
	  assert(0);
	}
	if(npe_current<nn_max)
	  npe_current = nn_max;
	for(unsigned int j = 0; j < npe_current; j++){
	  pe_chID_LST4[j] = pe_info_vec_LST4.at(i).chID.at(j);
	  pe_time_LST4[j] = pe_info_vec_LST4.at(i).time.at(j);
	}
      }
      //
      tree->Fill();      
    }
  }
  hfile = tree->GetCurrentFile();
  hfile->Write();
  hfile->Close();
}
