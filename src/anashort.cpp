//my
#include "anashort.hh"
#include "sipmCameraHist.hh"
#include "wfCamSim.hh"

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

void anashort::Loop(TString histOut){
  //
  const unsigned int nn_fadc_point = 75;
  const unsigned int nn_PMT_channels = 7987;
  //
  Int_t fadc_sum_offset = 15;
  //Float_t fadc_amplitude = 8.25;
  Int_t fadc_MHz = 1024;
  Float_t fadc_sample_in_ns = 1000.0/fadc_MHz;
  Float_t time_offset = fadc_sum_offset*fadc_sample_in_ns;
  Float_t NGB_rate_in_MHz = 386.0;
  //
  //auto wfcam = new Int_t [nn_PMT_channels][nn_fadc_point];
  vector<vector<Int_t>> wfcam(nn_PMT_channels, vector<Int_t>(nn_fadc_point));
  //
  TRandom3 *rnd = new TRandom3(123123);
  wfCamSim *wf = new wfCamSim( rnd, "Template_CTA_SiPM.txt", "spe.dat");
  //
  TH1D *h1_wf_ampl_ADC_test = new TH1D(); 
  h1_wf_ampl_ADC_test->SetNameTitle("h1_wf_ampl_ADC_test","h1_wf_ampl_ADC_test");
  wf->test_single_pe_amplitude_generator(h1_wf_ampl_ADC_test, 1000000);
  //wf->test_single_pe_amplitude_from_hist_generator(h1_wf_ampl_ADC_test, 10000000);
  // 
  TH1D *h1_nphotons = new TH1D("h1_nphotons","h1_nphotons",1000,0.0,10000);
  TH1D *h1_n_pe = new TH1D("h1_n_pe","h1_n_pe",1000,0.0,1000);
  TH1D *h1_n_pixels = new TH1D("h1_n_pixels","h1_n_pixels",1000,0.0,10000);
  //
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"nentries = "<<nentries<<endl;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if(jentry%100 == 0)
      cout<<jentry<<endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    //if (ientry > 10) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //
    h1_nphotons->Fill(nphotons);
    h1_n_pe->Fill(n_pe);
    h1_n_pixels->Fill(n_pixels);
    //
    //5363 -25.0553
    //5363 -25.0553
    //
    //5363 -25.0034
    //5363 -25.0034
    //
    //5380 -25.0254
    //5380 -25.0254
    //  
    //6685 -25.1979
    //6685 -25.1979
    //
    //6685 -25.1899
    //6685 -25.1899
    //
    if(n_pe<101){
      cout<<"n_pe "<<n_pe<<endl;
      for(Int_t i = 0;i<n_pe;i++){
	cout<<pe_chID[i]<<endl
	    <<pe_time[i]<<endl;
      }
      //
      wf->simulate_cam_event(nn_fadc_point,
			     nn_PMT_channels,
			     wfcam,
			     NGB_rate_in_MHz,
			     ev_time,
			     time_offset,
			     n_pe,
			     pe_chID,
			     pe_time);
      assert(0);
    }
    //
  }
  //
  //----------------------
  TFile* rootFile = new TFile(histOut.Data(), "RECREATE", " Histograms", 1);
  rootFile->cd();
  if (rootFile->IsZombie()){
    cout<<"  ERROR ---> file "<<histOut.Data()<<" is zombi"<<endl;
    assert(0);
  }
  else
    cout<<"  Output Histos file ---> "<<histOut.Data()<<endl;
  //
  h1_nphotons->Write();
  h1_n_pe->Write();
  h1_n_pixels->Write();
  //
  wf->getTemplate()->Write();
  wf->get_gr_wf_ampl()->Write();
  wf->get_h1_wf_ampl()->Write();
  wf->get_h1_wf_ampl_ADC()->Write();
  //
  h1_wf_ampl_ADC_test->Write();
  //
  rootFile->Close();
}
