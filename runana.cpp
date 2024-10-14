//my
#include "src/ana.hh"
#include "src/anaTrg.hh"
#include "src/anaTrgA.hh"
#include "src/anaTrgB.hh"
#include "src/anaEvPerEv.hh"
#include "src/anashort.hh"
#include "src/anaPCA.hh"
#include "src/anaPCAp.hh"
#include "src/anaFast.hh"
#include "src/anaSuperFast.hh"
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
  else if(argc == 14 && atoi(argv[1])==111){
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
    Double_t rsimulation = atof(argv[12]);
    TString trgSetup = argv[13];
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
	<<"data_chunk_ID : "<<data_chunk_ID<<endl
      	<<"rsimulation   : "<<rsimulation<<endl
	<<"trgSetup      : "<<trgSetup<<endl;
    //
    anaTrgA a( inRootFiles, 1);
    a.set_disable_energy_theta_rcore_binwise_cuts(true);
    a.set_rsimulation(rsimulation);
    a.set_trg_conf_file(trgSetup);
    a.Loop(outRootFileF, binE, binTheta, binDist, npe_min, npe_max, nEv_max, rndseed, data_chunk_ID);
  }
  else if(argc == 7 && atoi(argv[1])==112){
    //
    cout<<" runID = "<<atoi(argv[1])<<endl
	<<" runID = 112 --> TrgA NGB "<<endl;
    //
    TString inRootFiles = argv[2];
    TString outRootFileF = argv[3];
    Int_t nEv_max = atoi(argv[4]);
    Int_t rndseed = atoi(argv[5]);
    TString trgSetup = argv[6];
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
	<<"data_chunk_ID : "<<data_chunk_ID<<endl
      	<<"trgSetup      : "<<trgSetup<<endl;    
    //
    anaTrgA a( inRootFiles, 1);
    a.set_trg_conf_file(trgSetup);
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
    sipmCameraHist *sipm_cam = new sipmCameraHist("sipm_cam","sipm_cam","pixel_mapping.csv",0);
    sipm_cam->test();
    sipmCameraHist *sipm_cam_form_simpHist = new sipmCameraHist("sipm_cam_form_simpHist","sipm_cam_form_simpHist",sipm_cam);
    sipm_cam_form_simpHist->test("sipmCameraHist_formSimpHist_test.pdf");
    //
    sipmCameraHist *sipm_cam_test02 = new sipmCameraHist("sipm_cam_test02","sipm_cam_test02","pixel_mapping.csv",0);
    sipm_cam_test02->test02();
    //
    sipmCameraHist *sipm_cam_test03 = new sipmCameraHist("sipm_cam_test03","sipm_cam_test03","pixel_mapping.csv",0);
    sipm_cam_test03->test03();
    //
    sipmCameraHist *sipm_cam_test04 = new sipmCameraHist("sipm_cam_test04","sipm_cam_test04","pixel_mapping.csv",0);
    sipm_cam_test04->test04();
    //
    sipmCameraHist *sipm_cam_test05 = new sipmCameraHist("sipm_cam_test05","sipm_cam_test05","pixel_mapping.csv",0);
    sipm_cam_test05->test05();
    //
    sipmCameraHist *sipm_cam_test06 = new sipmCameraHist("sipm_cam_test06","sipm_cam_test06","pixel_mapping.csv",0);
    sipm_cam_test06->test_drawer_id();
    //
    sipmCameraHist *sipm_cam_test07 = new sipmCameraHist("sipm_cam_test07","sipm_cam_test07","pixel_mapping.csv",0);
    sipm_cam_test07->test_pixel_neighbors_id();
    sipmCameraHist *sipm_cam_test08 = new sipmCameraHist("sipm_cam_test08","sipm_cam_test08","pixel_mapping.csv",0);
    sipm_cam_test08->test_pixel_neighbors_id(1000);
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
    sipmCameraHist *sipm_cam_test09 = new sipmCameraHist("sipm_cam_test09","sipm_cam_test09","pixel_mapping.csv",0);
    sipm_cam_test09->test_pixel_neighbors_id(npix, pix_id);
    //
    sipmCameraHist *sipm_cam_test10 = new sipmCameraHist("sipm_cam_test10","sipm_cam_test10","pixel_mapping.csv",0);
    sipm_cam_test10->test_pixel_neighbors_second_id(npix, pix_id);
    //
    sipmCameraHist *sipm_cam_test11 = new sipmCameraHist("sipm_cam_test11","sipm_cam_test11","pixel_mapping.csv",0);
    sipm_cam_test11->test_pixel_neighbors_third_id(npix, pix_id);
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
  else if(argc == 2 && atoi(argv[1])==55){
    cout<<"--> Parameters <--"<<endl
	<<"argv[1] : "<<atoi(argv[1])<<endl;
    //
    sipmCameraHist *sipm_cam = new sipmCameraHist("sipm_cam","sipm_cam","pixel_mapping.csv",0);
    sipm_cam->save_pixel_neighbors_to_csv();
    sipm_cam->test055();
    //
    Int_t npix = 19;
    int *pix_id = new int[19];
    pix_id[0] = 0;
    pix_id[1] = 287;
    pix_id[2] = 6902;
    pix_id[3] = 5579;
    pix_id[4] = 4256;
    pix_id[5] = 2933;
    pix_id[6] = 1610;
    pix_id[7] = 1353;
    pix_id[8] = 7982;
    pix_id[9] = 6649;
    pix_id[10] = 5234;
    pix_id[11] = 4008;
    pix_id[12] = 2692;
    pix_id[13] = 153;
    pix_id[14] = 5449;
    pix_id[15] = 5444;
    pix_id[16] = 2820;
    pix_id[17] = 1497;
    pix_id[18] = 187;
    //
    sipmCameraHist *sipm_cam_test09 = new sipmCameraHist("sipm_cam_test09","sipm_cam_test09","pixel_mapping.csv",0);
    sipm_cam_test09->test_pixel_neighbors_id(npix, pix_id);
    //
    sipmCameraHist *sipm_cam_test10 = new sipmCameraHist("sipm_cam_test10","sipm_cam_test10","pixel_mapping.csv",0);
    sipm_cam_test10->test_pixel_neighbors_second_id(npix, pix_id);
    //
    sipmCameraHist *sipm_cam_test11 = new sipmCameraHist("sipm_cam_test11","sipm_cam_test11","pixel_mapping.csv",0);
    sipm_cam_test11->test_pixel_neighbors_third_id(npix, pix_id);
    //
    sipmCameraHist *sipm_cam_test12 = new sipmCameraHist("sipm_cam_test12","sipm_cam_test12","pixel_mapping.csv",0);
    sipm_cam_test12->test_pixel_super_flower(npix, pix_id);
  }
  else if(argc == 2 && atoi(argv[1])==555){
    cout<<"--> Parameters <--"<<endl
	<<"argv[1] : "<<atoi(argv[1])<<endl;
    //
    sipmCameraHist *sipm_cam = new sipmCameraHist("sipm_cam","sipm_cam","pixel_mapping.csv",0);
    sipm_cam->test_trigger_channel_mask_isolated_flower();
  }
  else if(argc == 2 && atoi(argv[1])==556){
    cout<<"--> Parameters <--"<<endl
	<<"argv[1] : "<<atoi(argv[1])<<endl;
    //
    sipmCameraHist *sipm_cam = new sipmCameraHist("sipm_cam","sipm_cam","pixel_mapping.csv",0);
    std::vector<Int_t> isolated_flower_seeds = sipm_cam->get_trigger_channel_mask_isolated_flower();
    TString gif_out_name_preff = "test_trigger_channel_mask_isolated_flower_plus_super_flower";
    TString gif_out_name;
    unsigned int seedID;
    for(unsigned int i = 0;i<2;i++){
      seedID = isolated_flower_seeds.at(i);
      gif_out_name = gif_out_name_preff;
      gif_out_name += seedID;
      gif_out_name += ".gif";
      sipm_cam->Clean();
      sipm_cam->test_trigger_channel_mask_isolated_flower_plus_super_flower(gif_out_name,seedID);
    }
    sipm_cam->Clean();
    sipm_cam->test_of_inefficient_regions_isolated_flower_plus_super_flower();    
  }
  else if(argc == 2 && atoi(argv[1])==5550){
    cout<<"--> Parameters <--"<<endl
	<<"argv[1] : "<<atoi(argv[1])<<endl;
    //
    sipmCameraHist *sipm_cam = new sipmCameraHist("sipm_cam","sipm_cam","pixel_mapping.csv",0);
    sipm_cam->save_trigger_channel_mask_isolated_flower();
    sipm_cam->save_trigger_channel_mask_all_pixels();
    sipm_cam->save_isolated_flower_seed_flower();
    sipm_cam->save_isolated_flower_seed_super_flower();
    sipm_cam->save_all_seed_flower();
    //
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
  else if(argc == 6 && atoi(argv[1])==63){
    TString rootFilesList = argv[2];
    TString outRootFileF = argv[3];
    Int_t nMaxEv = atoi(argv[4]);
    TString recoFile = argv[5];
    cout<<"--> Parameters <--"<<endl
	<<"rootFilesList : "<<rootFilesList<<endl
	<<"outRootFileF  : "<<outRootFileF<<endl
    	<<"nMaxEv        : "<<nMaxEv<<endl
	<<"recoFile      : "<<recoFile<<endl;
    anaPCAp a(rootFilesList,"anaPCA.conf");
    a.draw_reco(outRootFileF, nMaxEv, recoFile);
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
  else if(argc == 5 && atoi(argv[1])==77){
    TString rootFilesList = argv[2];
    TString outRootFileF = argv[3];
    TString conf_file = argv[4];
    cout<<"--> Parameters <--"<<endl
	<<"rootFilesList : "<<rootFilesList<<endl
	<<"outRootFileF  : "<<outRootFileF<<endl
    	<<"conf_file     : "<<conf_file<<endl;
    anaSuperFast a(rootFilesList, conf_file.Data());
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
  else if(argc == 9 && atoi(argv[1])==9999){
    //
    TString particle_type = argv[2];
    TString hist_file_prefix = argv[3];
    TString outRootFileF = argv[4];
    Int_t n_jobs = atoi(argv[5]);
    //
    Int_t n_jobs_NSB = atoi(argv[6]);
    TString hist_file_dir_NSB = argv[7];
    //
    Double_t rsimulation = atof(argv[8]);
    //
    TString rateCalc_name_title = "rateCalc";
    //
    Double_t val_Emin = 1.0;      // GeV
    Double_t val_Emax = 100000;   // GeV
    Int_t val_N_bins_E = 25;      //
    Double_t val_Thetamin = 0.0;  //deg
    Double_t val_Thetamax = 10.0; //deg
    Int_t val_N_bins_t = 10;
    //
    //
    evstHist *evH_flux= new evstHist("evH_flux","evH_flux",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
    //
    cout<<"--> Parameters: <--"<<endl
	<<"particle_type     : "<<particle_type<<endl
	<<"hist_file_prefix  : "<<hist_file_prefix<<endl
	<<"outRootFileF      : "<<outRootFileF<<endl
	<<"n_jobs            : "<<n_jobs<<endl
	<<"n_jobs_NSB        : "<<n_jobs_NSB<<endl
	<<"hist_file_dir_NSB : "<<hist_file_dir_NSB<<endl
	<<"rsimulation       : "<<rsimulation<<endl;
    //
    if(particle_type == "p"){
      rateCalc_name_title += "_proton";
      evH_flux->LoadBinContent("../cosmique_gamma_hadron_generator/flux_diff_protons.dat", true);
    }
    else{
      cout<<"--> Warning this is for gamma diffuse <--"<<endl
	  <<"but we consider diffused proton flux "<<endl
      	  <<"it is done for trigger effective area culculation only !!!"<<endl;
      rateCalc_name_title += "_gamma";
      evH_flux->LoadBinContent("../cosmique_gamma_hadron_generator/flux_diff_protons.dat", true);
      //assert(0);
    }
    //
    bool disable_energy_theta_rcore_binwise_cuts = true;
    rateCalculator r(rateCalc_name_title.Data(), rateCalc_name_title.Data(), hist_file_prefix, outRootFileF, evH_flux,
		     n_jobs, disable_energy_theta_rcore_binwise_cuts,
		     n_jobs_NSB, hist_file_dir_NSB, rsimulation);
  }
  else if(argc == 6 && atoi(argv[1])==333){
    TString inRootFiles = argv[2];
    Int_t event_ID = atoi(argv[3]);
    TString binFileOut = argv[4];
    Int_t rndseed = atoi(argv[5]);
    //
    cout<<"--> Parameters <--"<<endl
	<<"inRootFiles : "<<inRootFiles<<endl
      	<<"event_ID    : "<<event_ID<<endl
	<<"binFileOut  : "<<binFileOut<<endl
    	<<"rndseed     : "<<rndseed<<endl;
    //
    anaEvPerEv a( inRootFiles, 1);
    a.save_event_to_bin_file(binFileOut, event_ID, rndseed);
  }
  else if(argc == 4 && atoi(argv[1])==3333){
    TString inRootFiles = argv[2];
    TString txtFileOut = argv[3];
    //
    cout<<"--> Parameters <--"<<endl
	<<"inRootFiles : "<<inRootFiles<<endl
	<<"txtFileOut  : "<<txtFileOut<<endl;
    //
    anaEvPerEv a( inRootFiles, 1);
    a.get_events_map(txtFileOut);
  }
  else if(argc == 14 && atoi(argv[1])==222){
    TString inRootFile = argv[2];
    TString outRootFileHist = argv[3];
    Float_t NGB_rate_in_MHz = (Float_t)atof(argv[4]);
    Float_t fadc_electronic_noise_RMS = (Float_t)atof(argv[5]);
    Int_t npe_min = atoi(argv[6]);
    Int_t npe_max = atoi(argv[7]);
    Int_t nEvSim_max = atoi(argv[8]);
    Int_t nEv_max = atoi(argv[9]);
    Int_t rndseed = atoi(argv[10]);
    bool NGBsim;
    Int_t NGBsim_val = atoi(argv[11]);
    if(NGBsim_val == 1)
      NGBsim = true;
    else
      NGBsim = false;
    TString trgSetup = argv[12];
    TString name_ana_conf_file = argv[13];
    //
    cout<<"--> Parameters <--"<<endl
	<<"inRootFile                "<<inRootFile<<endl
	<<"outRootFileHist           "<<outRootFileHist<<endl
	<<"NGB_rate_in_MHz           "<<NGB_rate_in_MHz<<endl
	<<"fadc_electronic_noise_RMS "<<fadc_electronic_noise_RMS<<endl
	<<"npe_min                   "<<npe_min<<endl
	<<"npe_max                   "<<npe_max<<endl
	<<"nEvSim_max                "<<nEvSim_max<<endl
	<<"nEv_max                   "<<nEv_max<<endl
	<<"rndseed                   "<<rndseed<<endl
	<<"NGBsim                    "<<NGBsim<<endl
	<<"trgSetup                  "<<trgSetup<<endl
    	<<"name_ana_conf_file        "<<name_ana_conf_file<<endl;    
    //
    anaTrgB a( inRootFile, 1);
    a.set_trg_conf_file(trgSetup);
    a.set_ana_conf_file(name_ana_conf_file);
    a.Loop(outRootFileHist,
	   NGB_rate_in_MHz, fadc_electronic_noise_RMS,
	   npe_min, npe_max, nEvSim_max, nEv_max,
	   rndseed, NGBsim);
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
    cout<<" runID [1] = 5    (execution ID number) test of camera hist."<<endl;
    cout<<" runID [1] = 55   (execution ID number) test of camera hist: flower pixID (out csv file with mapping out csv file with mapping)."<<endl;
    cout<<" runID [1] = 555  (execution ID number) test of camera hist: test_trigger_channel_mask_isolated_flower."<<endl;
    cout<<" runID [1] = 556  (execution ID number) test of camera hist: test_trigger_channel_mask_isolated_flower_plus_super_flower."<<endl;
    cout<<" runID [1] = 5550 (execution ID number) save trigger channel mask isolated flower, all pixes, isolated flower seed flower and isolated flower seed super flower "<<endl;
    cout<<" runID [1] = 6 (execution ID number) short file format"<<endl
      	<<"       [2] - file with list of the root files"<<endl
	<<"       [3] - name of root file with histograms"<<endl;
    cout<<" runID [1] = 61 (execution ID number) short file format"<<endl
      	<<"       [2] - file with list of the root files"<<endl
	<<"       [3] - name of root file with histograms"<<endl;
    cout<<" runID [1] = 62 (execution ID number) draw principal"<<endl
      	<<"       [2] - file with list of the root files"<<endl
	<<"       [3] - name of root file with histograms"<<endl;
    cout<<" runID [1] = 63 (execution ID number) draw reco. shower (from SVD - principal)"<<endl
      	<<"       [2] - file with list of the root files"<<endl
	<<"       [3] - name of root file with histograms"<<endl
	<<"       [4] - number of showers (max)"<<endl
    	<<"       [5] - recoFile (ex: reco.csv)"<<endl;
    cout<<" runID [1] = 7 (execution ID number) short file format: fast"<<endl
      	<<"       [2] - file with list of the root files"<<endl
	<<"       [3] - name of root file with histograms"<<endl;
    cout<<" runID [1] = 77 (execution ID number) short file format: super fast"<<endl
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
	<<"       [11]- data_chunk_ID"<<endl
    	<<"       [12]- rsimulation"<<endl
    	<<"       [13]- trgSetup"<<endl;
    cout<<" runID [1] = 112 (execution ID number) TrgA NGB"<<endl
      	<<"       [2] - in root file"<<endl
	<<"       [3] - name of root file with histograms"<<endl
	<<"       [4] - nEv_max"<<endl
	<<"       [5]- rndseed"<<endl
    	<<"       [6]- trgSetup"<<endl;
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
	<<"       [5] - n_jons"<<endl
      	<<"       [6] - n_jobs_NSB"<<endl
	<<"       [7] - hist_file_dir_NSB"<<endl
	<<"       [8] - rsimulation"<<endl;
    cout<<" runID [1] = 333 (execution ID number) anaEvPerEv"<<endl
      	<<"       [2] - in root file"<<endl
      	<<"       [3] - event ID"<<endl
	<<"       [4] - name of the binary file with event"<<endl
    	<<"       [5]- rndseed"<<endl;
    cout<<" runID [1] = 3333 (execution ID number) anaEvPerEv to build get_events_map"<<endl
      	<<"       [2] - in root file"<<endl
	<<"       [3] - name of txt file with event map"<<endl;
    cout<<" runID [1] = 222 (execution ID number) TrgB "<<endl
      	<<"       [2] - in root file"<<endl
	<<"       [3] - name of root file with histograms"<<endl
      	<<"       [4] - NGB_rate_in_MHz"<<endl
      	<<"       [5] - fadc_electronic_noise_RMS"<<endl
	<<"       [6] - npe_min"<<endl
      	<<"       [7] - npe_max"<<endl
	<<"       [8] - nEvSim_max"<<endl
	<<"       [9] - nEv_max"<<endl
    	<<"       [10]- rndseed"<<endl
      	<<"       [11]- NGBsim"<<endl
    	<<"       [12]- trgSetup"<<endl;
  }
  return 0;
}
