//my
#include "anaEvPerEv.hh"
#include "ana.hh"
#include "anabase.hh"
#include "sipmCameraHist.hh"
#include "wfCamSim.hh"
#include "triggerSim.hh"
#include "evstHist.hh"
#include "dbscan.hh"

//root
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TVector3.h>

//C, C++
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <bits/stdc++.h>
#include <sys/stat.h>

using namespace std;

void anaEvPerEv::get_events_map(TString txtFileOut){
  //
  clock_t start, finish;
  clock_t start_sim, finish_sim;
  start = clock();
  //
  finish = clock();
  cout<<"initialization time : "<<((finish - start)/CLOCKS_PER_SEC)<<" (sec)"<<endl;	
  //  
  start_sim = clock();
  //
  ofstream fileOut;
  fileOut.open(txtFileOut.Data());
  fileOut<<setw(15)<<"jentry_int"
	 <<setw(15)<<"energy"
	 <<setw(15)<<"n_pe"
	 <<setw(15)<<"n_pixels"
	 <<setw(15)<<endl;
  //
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"nentries = "<<nentries<<endl;
  Long64_t nbytes = 0, nb = 0;
  Int_t jentry_int;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry); nbytes += nb;
    jentry_int = (int)jentry;
    fileOut<<setw(15)<<jentry_int
	   <<setw(15)<<energy
	   <<setw(15)<<n_pe
      	   <<setw(15)<<n_pixels
	   <<setw(15)<<endl;
  }
  finish_sim = clock();
  cout<<" "<<((finish_sim - start_sim)/(CLOCKS_PER_SEC))<<" (sec)"<<endl;    
  fileOut.close();
}
 
void anaEvPerEv::save_event_to_bin_file(TString binFileOut, Int_t event_ID, Int_t rndseed){
  //
  clock_t start, finish;
  clock_t start_sim, finish_sim;
  start = clock();
  //
  const unsigned int nn_PMT_channels = nChannels;
  Int_t fadc_sum_offset = 15;
  Int_t fadc_MHz = 1024;
  Int_t fadc_offset = 300;
  Float_t fadc_sample_in_ns = 1000.0/fadc_MHz;
  Float_t time_offset = fadc_sum_offset*fadc_sample_in_ns;
  Float_t NGB_rate_in_MHz = 268.0;
  //Float_t NGB_rate_in_MHz = 0.1;
  Float_t fadc_electronic_noise_RMS = 3.8082498; //takes into account 3.5/sqrt(12)
  //Float_t fadc_electronic_noise_RMS = 0.01;
  //
  TRandom3 *rnd = new TRandom3(rndseed);
  //
  vector<vector<Int_t>> wfcam(nn_PMT_channels, vector<Int_t>(nn_fadc_point));
  wfCamSim *wfc = new wfCamSim(rnd, "Template_CTA_SiPM.txt", "spe.dat",
			       nn_fadc_point, nn_PMT_channels, fadc_offset, fadc_sample_in_ns, NGB_rate_in_MHz, fadc_electronic_noise_RMS);
  wfc->print_wfCamSim_configure();  
  //
  sipmCameraHist *sipm_cam = new sipmCameraHist("sipm_cam","sipm_cam","pixel_mapping.csv",0);
  finish = clock();
  cout<<"initialization time : "<<((finish - start)/CLOCKS_PER_SEC)<<" (sec)"<<endl;	
  //  
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"nentries = "<<nentries<<endl;
  Long64_t nbytes = 0, nb = 0;
  Long64_t jentry = (Long64_t)event_ID;
  Long64_t ientry = LoadTree(jentry);
  nb = fChain->GetEntry(jentry); nbytes += nb;
  start_sim = clock();
  wfc->simulate_cam_event(nn_fadc_point,
			  nn_PMT_channels,
			  wfcam,
			  ev_time,
			  time_offset,
			  n_pe,
			  pe_chID,
			  pe_time);
  finish_sim = clock();
  cout<<" "<<((finish_sim - start_sim)/(CLOCKS_PER_SEC/1000))<<" (msec)"<<endl;    
  //
  int val;
  //
  int event_ID_int = (int)event_ID;
  float energy_float = (float)energy;
  int n_pe_int = (int)n_pe;
  int n_pixels_int = (int)n_pixels;
  //
  int nn_PMT_channels_int = (int)nn_PMT_channels;
  int nn_fadc_point_int = (int)nn_fadc_point;
  float fadc_sample_in_ns_float = (float)fadc_sample_in_ns;
  float NGB_rate_in_MHz_float = (float)NGB_rate_in_MHz;
  float fadc_electronic_noise_RMS_float = (float)fadc_electronic_noise_RMS;  
  //  
  //
  //int(4 bytes)    event_ID_int
  //float(4 bytes)  energy_float
  //int(4 bytes)    n_pe_int
  //int(4 bytes)    n_pixels_int
  //
  //int(4 bytes)   nn_PMT_channels_int
  //int(4 bytes)   nn_fadc_point_int
  //float(4 bytes) fadc_sample_in_ns_float
  //float(4 bytes) NGB_rate_in_MHz_float
  //float(4 bytes) fadc_electronic_noise_RMS_float
  //
  ofstream outBinfile;
  outBinfile.open(binFileOut.Data(), ios::out | ios::app | ios::binary);
  //
  outBinfile.write((char*)&event_ID_int, sizeof(event_ID_int));
  outBinfile.write((char*)&energy_float, sizeof(energy_float));
  outBinfile.write((char*)&n_pe_int, sizeof(n_pe_int));
  outBinfile.write((char*)&n_pixels_int, sizeof(n_pixels_int));    
  //
  outBinfile.write((char*)&nn_PMT_channels_int, sizeof(nn_PMT_channels_int));
  outBinfile.write((char*)&nn_fadc_point_int, sizeof(nn_fadc_point_int));
  outBinfile.write((char*)&fadc_sample_in_ns_float, sizeof(fadc_sample_in_ns_float));
  outBinfile.write((char*)&NGB_rate_in_MHz_float, sizeof(NGB_rate_in_MHz_float));
  outBinfile.write((char*)&fadc_electronic_noise_RMS_float, sizeof(fadc_electronic_noise_RMS_float));
  //
  int ch_id_int;
  float ch_t_float;
  //
  //vector<TGraph*> v_gr;
  //
  for(unsigned int i = 0;i<nn_PMT_channels;i++){
    //TGraph *gr = new TGraph();
    //TString grnameTitle = "gr_";
    //grnameTitle += v_gr.size();
    //gr->SetNameTitle(grnameTitle.Data(),grnameTitle.Data());
    for(unsigned int j = 0;j<(unsigned int)nn_fadc_point;j++){
      val=wfcam.at(i).at(j);
      //
      //gr->SetPoint(j,j*fadc_sample_in_ns_float,val);
      //
      outBinfile.write((char*)&val, sizeof(val));
    }
    //v_gr.push_back(gr);
  }
  //
  for(unsigned int i = 0;i<n_pe;i++){
    ch_id_int = (int)pe_chID[i];
    ch_t_float = (float)(pe_time[i] - ev_time + time_offset);
    outBinfile.write((char*)&ch_id_int, sizeof(ch_id_int));
    outBinfile.write((char*)&ch_t_float, sizeof(ch_t_float));
    //cout<<ch_id_int<<" "<<ch_t_float<<endl;
  }  
  //
  outBinfile.close();
  //
  //TFile* rootFile = new TFile("gr_test.root", "RECREATE", " Histograms", 1);
  //rootFile->cd();
  //cout<<"  Output Histos file ---> "<<"gr_test.root"<<endl;
  //for(unsigned int i = 0;i<nn_PMT_channels;i++)
  //  v_gr.at(i)->Write();
  //rootFile->Close();
  //
}

