//my
#include "src/ana.hh"
#include "src/anaTrg.hh"
#include "src/anaTrgA.hh"
#include "src/anashort.hh"
#include "src/anaPCA.hh"
#include "src/anaPCAp.hh"
#include "src/anaFast.hh"
#include "src/sipmCameraHist.hh"
#include "src/sipmCameraHistCropped.hh"
#include "src/evstHist.hh"
#include "src/rateCalculator.hh"


//root
#include "TROOT.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TGraph.h"

//C, C++
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <fstream>
#include <iomanip>
#include <time.h>

using namespace std;

int main(int argc, char *argv[]){
  cout<<"--> main "<<endl;
  if(argc == 4 && atoi(argv[1])==0){
    TString rootFilesList = argv[2];
    TString outRootFileF = argv[3];
    cout<<"--> Parameters <--"<<endl
	<<"rootFilesList : "<<rootFilesList<<endl
	<<"outRootFileF  : "<<outRootFileF<<endl;
    ana a(rootFilesList);
    a.Loop(outRootFileF);
  }
  else if(argc == 4 && atoi(argv[1])==100){
    TString inRootFiles = argv[2];
    TString outRootFileF = argv[3];
    cout<<"--> Parameters <--"<<endl
	<<"inRootFiles   : "<<inRootFiles<<endl
	<<"outRootFileF  : "<<outRootFileF<<endl;
    anaTrgA a( inRootFiles, 1);
    a.SiPM_dist(outRootFileF);
  }
  else if(argc == 4 && atoi(argv[1])==1){
    TString inRootFiles = argv[2];
    TString outRootFileF = argv[3];
    cout<<"--> Parameters <--"<<endl
	<<"inRootFiles   : "<<inRootFiles<<endl
	<<"outRootFileF  : "<<outRootFileF<<endl;
    ana a( inRootFiles, atoi(argv[1]));
    a.Loop(outRootFileF);
  }
  else if(argc == 4 && atoi(argv[1])==11){
    TString inRootFiles = argv[2];
    TString outRootFileF = argv[3];
    cout<<"--> Parameters <--"<<endl
	<<"inRootFiles   : "<<inRootFiles<<endl
	<<"outRootFileF  : "<<outRootFileF<<endl;
    anaTrg a( inRootFiles, 1);
    a.Loop(outRootFileF);
  }
  else if(argc == 12 && atoi(argv[1])==111){
    TString inRootFiles = argv[2];
    TString outRootFileF = argv[3];
    Int_t binE = atoi(argv[4]);
    Int_t binTheta = atoi(argv[5]);
    Int_t binDist = atoi(argv[6]);
    Int_t npe_min = atoi(argv[7]);
    Int_t npe_max = atoi(argv[8]);
    Int_t nEv_max = atoi(argv[9]);
    Int_t rndseed = atoi(argv[10]);
    Int_t data_chunk_ID = atoi(argv[11]);
    cout<<"--> Parameters <--"<<endl
	<<"inRootFiles   : "<<inRootFiles<<endl
	<<"outRootFileF  : "<<outRootFileF<<endl
	<<"binE          : "<<binE<<endl
	<<"binTheta      : "<<binTheta<<endl
	<<"binDist       : "<<binDist<<endl
	<<"npe_min       : "<<npe_min<<endl
	<<"npe_max       : "<<npe_max<<endl
	<<"nEv_max       : "<<nEv_max<<endl
    	<<"rndseed       : "<<rndseed<<endl
	<<"data_chunk_ID : "<<data_chunk_ID<<endl;
    //
    anaTrgA a( inRootFiles, 1);
    a.set_disable_energy_theta_rcore_binwise_cuts(true);
    a.Loop(outRootFileF, binE, binTheta, binDist, npe_min, npe_max, nEv_max, rndseed, data_chunk_ID);
  }
  else if(argc == 6 && atoi(argv[1])==112){
    TString inRootFiles = argv[2];
    TString outRootFileF = argv[3];
    Int_t nEv_max = atoi(argv[4]);
    Int_t rndseed = atoi(argv[5]);
    Int_t binE = 9;
    Int_t binTheta = 4;
    Int_t binDist = 7;
    Int_t npe_min = 1;
    Int_t npe_max = 100000;    
    Int_t data_chunk_ID = -999;
    cout<<"--> Parameters <--"<<endl
	<<"inRootFiles   : "<<inRootFiles<<endl
	<<"outRootFileF  : "<<outRootFileF<<endl
	<<"binE          : "<<binE<<endl
	<<"binTheta      : "<<binTheta<<endl
	<<"binDist       : "<<binDist<<endl
	<<"npe_min       : "<<npe_min<<endl
	<<"npe_max       : "<<npe_max<<endl
	<<"nEv_max       : "<<nEv_max<<endl
    	<<"rndseed       : "<<rndseed<<endl
	<<"data_chunk_ID : "<<data_chunk_ID<<endl;
    //
    anaTrgA a( inRootFiles, 1);
    a.Loop(outRootFileF, binE, binTheta, binDist, npe_min, npe_max, nEv_max, rndseed, data_chunk_ID, true);
  }
  else if(argc == 6 && atoi(argv[1])==113){
    TString inRootFiles = argv[2];
    TString outRootFileF = argv[3];
    Int_t nEv_max = atoi(argv[4]);
    Int_t rndseed = atoi(argv[5]);
    Int_t binE = 9;
    Int_t binTheta = 4;
    Int_t binDist = 7;
    Int_t npe_min = 1;
    Int_t npe_max = 100000;    
    Int_t data_chunk_ID = -999;
    bool k_dist_graph_flag = true;
    cout<<"--> Parameters <--"<<endl
	<<"inRootFiles       : "<<inRootFiles<<endl
	<<"outRootFileF      : "<<outRootFileF<<endl
	<<"binE              : "<<binE<<endl
	<<"binTheta          : "<<binTheta<<endl
	<<"binDist           : "<<binDist<<endl
	<<"npe_min           : "<<npe_min<<endl
	<<"npe_max           : "<<npe_max<<endl
	<<"nEv_max           : "<<nEv_max<<endl
    	<<"rndseed           : "<<rndseed<<endl
	<<"data_chunk_ID     : "<<data_chunk_ID<<endl
      	<<"k_dist_graph_flag : "<<k_dist_graph_flag<<endl;
    //
    anaTrgA a( inRootFiles, 1);
    a.set_k_dist_graph_flag(k_dist_graph_flag);
    a.Loop(outRootFileF, binE, binTheta, binDist, npe_min, npe_max, nEv_max, rndseed, data_chunk_ID, true);
  }
  else if(argc == 4 && atoi(argv[1])==1111){
    TString inRootFiles = argv[2];
    TString outRootFileF = argv[3];
    cout<<"--> Parameters <--"<<endl
	<<"inRootFiles   : "<<inRootFiles<<endl
	<<"outRootFileF  : "<<outRootFileF<<endl;
    //
    anaTrgA a( inRootFiles, 1);
    a.test_single_pe_amplitude_generator(outRootFileF);
  }
  else if(argc == 4 && atoi(argv[1])==2){
    TString rootFilesList = argv[2];
    TString outRootFileF = argv[3];
    cout<<"--> Parameters <--"<<endl
	<<"rootFilesList : "<<rootFilesList<<endl
	<<"outRootFileF  : "<<outRootFileF<<endl;
    anashort a(rootFilesList);
    a.Loop(outRootFileF);
  }
  else if(argc == 4 && atoi(argv[1])==3){
    TString inRootFiles = argv[2];
    TString outRootFileF = argv[3];
    cout<<"--> Parameters <--"<<endl
	<<"inRootFiles   : "<<inRootFiles<<endl
	<<"outRootFileF  : "<<outRootFileF<<endl;
    anashort a( inRootFiles, atoi(argv[1]));
    a.Loop(outRootFileF);
  }
  else if(argc == 6 && atoi(argv[1])==4){
    TString inRootFiles = argv[2];
    TString outRootFileF = argv[3];
    Long64_t evID = (Long64_t)atoi(argv[4]);
    TString particle_type_name = argv[5];
    cout<<"--> Parameters <--"<<endl
	<<"inRootFiles        : "<<inRootFiles<<endl
	<<"outRootFileF       : "<<outRootFileF<<endl
      	<<"evID               : "<<evID<<endl
	<<"particle_type_name : "<<particle_type_name<<endl;
    //ana a(inRootFiles, 1);
    ana a(inRootFiles);
    a.set_particle_type_name(particle_type_name);
    a.save_wf_for_event(outRootFileF, evID);
  }
  else if(argc == 2 && atoi(argv[1])==5){
    cout<<"--> Parameters <--"<<endl
	<<"argv[1] : "<<atoi(argv[1])<<endl;
    //sipmCameraHist *sipm_cam = new sipmCameraHist("sipm_cam","sipm_cam","pixel_mapping.csv",0);
    //sipm_cam->test();
    //
    //sipmCameraHist *sipm_cam = new sipmCameraHist("sipm_cam","sipm_cam","pixel_mapping.csv",0);
    //sipm_cam->test();
    //sipmCameraHist *sipm_cam_form_simpHist = new sipmCameraHist("sipm_cam_form_simpHist","sipm_cam_form_simpHist",sipm_cam);
    //sipm_cam_form_simpHist->test("sipmCameraHist_formSimpHist_test.pdf");
    //
    //sipmCameraHist *sipm_cam_test02 = new sipmCameraHist("sipm_cam_test02","sipm_cam_test02","pixel_mapping.csv",0);
    //sipm_cam_test02->test02();
    //
    //sipmCameraHist *sipm_cam_test03 = new sipmCameraHist("sipm_cam_test03","sipm_cam_test03","pixel_mapping.csv",0);
    //sipm_cam_test03->test03();
    //
    //sipmCameraHist *sipm_cam_test04 = new sipmCameraHist("sipm_cam_test04","sipm_cam_test04","pixel_mapping.csv",0);
    //sipm_cam_test04->test04();
    //
    //sipmCameraHist *sipm_cam_test05 = new sipmCameraHist("sipm_cam_test05","sipm_cam_test05","pixel_mapping.csv",0);
    //sipm_cam_test05->test05();
    //
    //sipmCameraHist *sipm_cam_test06 = new sipmCameraHist("sipm_cam_test06","sipm_cam_test06","pixel_mapping.csv",0);
    //sipm_cam_test06->test_drawer_id();
    //
    //sipmCameraHist *sipm_cam_test07 = new sipmCameraHist("sipm_cam_test07","sipm_cam_test07","pixel_mapping.csv",0);
    //sipm_cam_test07->test_pixel_neighbors_id();
    //sipmCameraHist *sipm_cam_test08 = new sipmCameraHist("sipm_cam_test08","sipm_cam_test08","pixel_mapping.csv",0);
    //sipm_cam_test08->test_pixel_neighbors_id(1000);
    //
    Int_t npix = 10;
    int *pix_id = new int[10];
    pix_id[0] = 0;
    pix_id[1] = 50;
    pix_id[2] = 100;
    pix_id[3] = 300;
    pix_id[4] = 400;
    pix_id[5] = 500;
    pix_id[6] = 800;
    pix_id[7] = 1000;
    pix_id[8] = 2000;
    pix_id[9] = 5000;
    //sipmCameraHist *sipm_cam_test09 = new sipmCameraHist("sipm_cam_test09","sipm_cam_test09","pixel_mapping.csv",0);
    //sipm_cam_test09->test_pixel_neighbors_id(npix, pix_id);
    //
    //sipmCameraHist *sipm_cam_test10 = new sipmCameraHist("sipm_cam_test10","sipm_cam_test10","pixel_mapping.csv",0);
    //sipm_cam_test10->test_pixel_neighbors_second_id(npix, pix_id);
    //
    //sipmCameraHist *sipm_cam_test11 = new sipmCameraHist("sipm_cam_test11","sipm_cam_test11","pixel_mapping.csv",0);
    //sipm_cam_test11->test_pixel_neighbors_third_id(npix, pix_id);
    //
    //TH1D *h1_distance_between_pixels = new TH1D("h1_distance_between_pixels","h1_distance_between_pixels",100000,-0.001,3.0);
    //sipmCameraHist *sipm_cam_test10 = new sipmCameraHist("sipm_cam_test10","sipm_cam_test10","pixel_mapping.csv",0,h1_distance_between_pixels);
    //h1_distance_between_pixels->SaveAs("h1_distance_between_pixels.C");
    //
    sipmCameraHist *sipm_cam_test12 = new sipmCameraHist("sipm_cam_test12","sipm_cam_test12","pixel_mapping.csv",0);
    sipm_cam_test12->test_pixel_super_flower(npix, pix_id);
    //sipm_cam_test12->test_pixel_neighbors_bubbleSort(0);
    //
    sipmCameraHistCropped *sipm_cam_crop = new sipmCameraHistCropped("sipm_cam_crop","sipm_cam_crop",sipm_cam_test12);
    sipm_cam_crop->test();
    sipm_cam_crop->test01(sipm_cam_test12);
    sipm_cam_crop->test02();
  }
  else if(argc == 4 && atoi(argv[1])==6){
    TString rootFilesList = argv[2];
    TString outRootFileF = argv[3];
    cout<<"--> Parameters <--"<<endl
	<<"rootFilesList : "<<rootFilesList<<endl
	<<"outRootFileF  : "<<outRootFileF<<endl;
    anaPCA a(rootFilesList,"anaPCA.conf");
    a.Loop(outRootFileF);
  }
  else if(argc == 4 && atoi(argv[1])==61){
    TString rootFilesList = argv[2];
    TString outRootFileF = argv[3];
    cout<<"--> Parameters <--"<<endl
	<<"rootFilesList : "<<rootFilesList<<endl
	<<"outRootFileF  : "<<outRootFileF<<endl;
      //<<"phi0_shift    : "<<atof(argv[4])<<endl;
    anaPCAp a(rootFilesList,"anaPCA.conf");
    //a._phi0_shift=atof(argv[4]);
    a.Loop(outRootFileF);
  }
  else if(argc == 4 && atoi(argv[1])==62){
    TString rootFilesList = argv[2];
    TString outRootFileF = argv[3];
    cout<<"--> Parameters <--"<<endl
	<<"rootFilesList : "<<rootFilesList<<endl
	<<"outRootFileF  : "<<outRootFileF<<endl;
    anaPCAp a(rootFilesList,"anaPCA.conf");
    a.draw_principal(outRootFileF);
  }
  else if(argc == 5 && atoi(argv[1])==7){
    TString rootFilesList = argv[2];
    TString outRootFileF = argv[3];
    TString conf_file = argv[4];
    cout<<"--> Parameters <--"<<endl
	<<"rootFilesList : "<<rootFilesList<<endl
	<<"outRootFileF  : "<<outRootFileF<<endl
    	<<"conf_file     : "<<conf_file<<endl;
    anaFast a(rootFilesList, conf_file.Data());
    a.Loop(outRootFileF);
  }
  else if(argc == 7 && atoi(argv[1])==8){
    TString rootFilesList = argv[2];
    TString outRootFileF = argv[3];
    TString conf_file = argv[4];
    TString simtel_all_dat = argv[5];
    TString flux_dat = argv[6];
    cout<<"--> Parameters <--"<<endl
	<<"rootFilesList  : "<<rootFilesList<<endl
	<<"outRootFileF   : "<<outRootFileF<<endl
    	<<"conf_file      : "<<conf_file<<endl
	<<"simtel_all_dat : "<<simtel_all_dat<<endl
	<<"flux_dat       : "<<flux_dat<<endl;
    anaFast a(rootFilesList, conf_file.Data());
    a.Loop_scan(outRootFileF, simtel_all_dat, flux_dat);
  }
  else if(argc == 2 && atoi(argv[1])==888){
    Double_t val_Emin = 1.0;      // GeV
    Double_t val_Emax = 100000;   // GeV
    Int_t val_N_bins_E = 25;
    Double_t val_Thetamin = 0.0;  //deg
    Double_t val_Thetamax = 10.0; //deg
    Int_t val_N_bins_t = 10;
    //
    evstHist *evH = new evstHist("evH","evH",
				 val_Emin, val_Emax, val_N_bins_E,
				 val_Thetamin, val_Thetamax, val_N_bins_t);
    evH->test();
    evH->test_Get_th_bin_ID_and_e_bin_ID();
    //
    evstHist *evH02 = new evstHist("evH02","evH02",
				   val_Emin, val_Emax, val_N_bins_E,
				   val_Thetamin, val_Thetamax, val_N_bins_t);
    evH02->test_get_bin(1.1, 0.1, 1);
    evH02->test_get_bin(10.1, 0.1, 2);
    evH02->test_get_bin(10.1, 9.1, 3);
    evH02->test_get_bin(100.1, 9.1, 4);
    evH02->test_get_bin(1000.1, 9.1, 5);
    evH02->test_get_bin(10000.1, 8.1, 6);
    evH02->test_get_bin(99999.1, 7.1, 7);
    evH02->Draw_hist("evH_test_get_bin.pdf");    
    //
    evstHist *evH03 = new evstHist("evH03","evH03",
				   val_Emin, val_Emax, val_N_bins_E,
				   val_Thetamin, val_Thetamax, val_N_bins_t);
    evH03->test_get_arbitrary_hist_ID();
  }
  else if(argc == 5 && atoi(argv[1])==999){
    //
    TString particle_type = argv[2];
    TString hist_file_prefix = argv[3];
    TString outRootFileF = argv[4];
    TString rateCalc_name_title = "rateCalc";
    //
    Double_t val_Emin = 1.0;      // GeV
    Double_t val_Emax = 100000;   // GeV
    Int_t val_N_bins_E = 25;      //
    Double_t val_Thetamin = 0.0;  //deg
    Double_t val_Thetamax = 10.0; //deg
    Int_t val_N_bins_t = 10;
    evstHist *evH_flux= new evstHist("evH_flux","evH_flux",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
    //
    cout<<"--> Parameters <--"<<endl
	<<"particle_type    : "<<particle_type<<endl
	<<"hist_file_prefix : "<<hist_file_prefix<<endl;
    //
    if(particle_type = "p"){
      rateCalc_name_title += "_proton";
      evH_flux->LoadBinContent("../cosmique_gamma_hadron_generator/flux_diff_protons.dat", true);
    }
    else{
      assert(0);
    }
    //
    rateCalculator r(rateCalc_name_title.Data(), rateCalc_name_title.Data(), hist_file_prefix, outRootFileF, evH_flux);
  }
  else if(argc == 6 && atoi(argv[1])==9999){
    //
    TString particle_type = argv[2];
    TString hist_file_prefix = argv[3];
    TString outRootFileF = argv[4];
    Int_t n_jobs = atoi(argv[5]);
    TString rateCalc_name_title = "rateCalc";
    //
    Double_t val_Emin = 1.0;      // GeV
    Double_t val_Emax = 100000;   // GeV
    Int_t val_N_bins_E = 25;      //
    Double_t val_Thetamin = 0.0;  //deg
    Double_t val_Thetamax = 10.0; //deg
    Int_t val_N_bins_t = 10;
    evstHist *evH_flux= new evstHist("evH_flux","evH_flux",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
    //
    cout<<"--> Parameters: <--"<<endl
	<<"particle_type    : "<<particle_type<<endl
	<<"hist_file_prefix : "<<hist_file_prefix<<endl
	<<"outRootFileF     : "<<outRootFileF<<endl
	<<"n_jobs           : "<<n_jobs<<endl;
    //
    if(particle_type = "p"){
      rateCalc_name_title += "_proton";
      evH_flux->LoadBinContent("../cosmique_gamma_hadron_generator/flux_diff_protons.dat", true);
    }
    else{
      assert(0);
    }
    //
    bool disable_energy_theta_rcore_binwise_cuts = true;
    rateCalculator r(rateCalc_name_title.Data(), rateCalc_name_title.Data(), hist_file_prefix, outRootFileF, evH_flux, n_jobs, disable_energy_theta_rcore_binwise_cuts);
  }
  else{
    cout<<" --> ERROR in input arguments "<<endl
	<<" runID [1] = 0 (execution ID number)"<<endl
      	<<"       [2] - file with list of the root files"<<endl
	<<"       [3] - name of root file with histograms"<<endl;
    cout<<" runID [1] = 1 (execution ID number)"<<endl
      	<<"       [2] - in root files"<<endl
	<<"       [3] - name of root file with histograms"<<endl;
    cout<<" runID [1] = 2 (execution ID number) short file format"<<endl
      	<<"       [2] - file with list of the root files"<<endl
	<<"       [3] - name of root file with histograms"<<endl;
    cout<<" runID [1] = 3 (execution ID number) short file format"<<endl
      	<<"       [2] - in root files"<<endl
	<<"       [3] - name of root file with histograms"<<endl;
    cout<<" runID [1] = 4 (execution ID number) short file format"<<endl
      	<<"       [2] - in root files"<<endl
	<<"       [3] - name of root file with histograms"<<endl
    	<<"       [4] - evID"<<endl;
    cout<<" runID [1] = 5 (execution ID number) test of camera hist."<<endl;
    cout<<" runID [1] = 6 (execution ID number) short file format"<<endl
      	<<"       [2] - file with list of the root files"<<endl
	<<"       [3] - name of root file with histograms"<<endl;
    cout<<" runID [1] = 61 (execution ID number) short file format"<<endl
      	<<"       [2] - file with list of the root files"<<endl
	<<"       [3] - name of root file with histograms"<<endl;
    cout<<" runID [1] = 62 (execution ID number) draw principal"<<endl
      	<<"       [2] - file with list of the root files"<<endl
	<<"       [3] - name of root file with histograms"<<endl;
    cout<<" runID [1] = 7 (execution ID number) short file format: fast"<<endl
      	<<"       [2] - file with list of the root files"<<endl
	<<"       [3] - name of root file with histograms"<<endl;
    cout<<" runID [1] = 8 (execution ID number) short file format: fast npe scan"<<endl
      	<<"       [2] - file with list of the root files"<<endl
	<<"       [3] - name of root file with histograms"<<endl;
    cout<<" runID [1] = 11 (execution ID number) Trg "<<endl
      	<<"       [2] - in root file"<<endl
	<<"       [3] - name of root file with histograms"<<endl;
    cout<<" runID [1] = 111 (execution ID number) TrgA "<<endl
      	<<"       [2] - in root file"<<endl
	<<"       [3] - name of root file with histograms"<<endl
      	<<"       [4] - binE"<<endl
	<<"       [5] - binTheta"<<endl
      	<<"       [6] - binDist"<<endl
	<<"       [7] - npe_min"<<endl
      	<<"       [8] - npe_max"<<endl
	<<"       [9] - nEv_max"<<endl
    	<<"       [10]- rndseed"<<endl
	<<"       [11]- data_chunk_ID"<<endl;
    cout<<" runID [1] = 112 (execution ID number) TrgA NGB"<<endl
      	<<"       [2] - in root file"<<endl
	<<"       [3] - name of root file with histograms"<<endl
	<<"       [4] - nEv_max"<<endl
	<<"       [5]- rndseed"<<endl;
    cout<<" runID [1] = 113 (execution ID number) TrgA NGB k-dist plots"<<endl
      	<<"       [2] - in root file"<<endl
	<<"       [3] - name of root file with histograms"<<endl
	<<"       [4] - nEv_max"<<endl
	<<"       [5]- rndseed"<<endl;
    cout<<" runID [1] = 1111 (execution ID number) test single pe amplitude generator "<<endl
      	<<"       [2] - in root file"<<endl
	<<"       [3] - name of root file with histograms"<<endl;
    cout<<" runID [1] = 100 (execution ID number) SiPM distributions "<<endl
	<<"       [2] - in root file"<<endl
	<<"       [3] - name of root file with histograms"<<endl;
    cout<<" runID [1] = 888 (execution ID number) test of evstHist"<<endl;
    cout<<" runID [1] = 999 (execution ID number) rate calculator binwise"<<endl
	<<"       [2] - particle type (g,gd,e,p)"<<endl
      	<<"       [3] - histogram file prefix"<<endl
    	<<"       [4] - name of root file with histograms"<<endl;
    cout<<" runID [1] = 9999 (execution ID number) rate calculator"<<endl
	<<"       [2] - particle type (g,gd,e,p)"<<endl
      	<<"       [3] - histogram file prefix"<<endl
    	<<"       [4] - name of root file with histograms"<<endl
	<<"       [5] - n_jons"<<endl;
  }
  return 0;
}
