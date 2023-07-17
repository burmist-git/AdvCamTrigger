//my
#include "src/ana.hh"
#include "src/anashort.hh"
#include "src/sipmCameraHist.hh"

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
  if(argc == 4 && atoi(argv[1])==0){
    TString rootFilesList = argv[2];
    TString outRootFileF = argv[3];
    cout<<"--> Parameter calculation from the WF <--"<<endl
	<<"rootFilesList : "<<rootFilesList<<endl
	<<"outRootFileF  : "<<outRootFileF<<endl;
    ana a(rootFilesList);
    a.Loop(outRootFileF);
  }
  else if(argc == 4 && atoi(argv[1])==1){
    TString inRootFiles = argv[2];
    TString outRootFileF = argv[3];
    cout<<"--> Parameter calculation from the WF <--"<<endl
	<<"inRootFiles   : "<<inRootFiles<<endl
	<<"outRootFileF  : "<<outRootFileF<<endl;
    ana a( inRootFiles, atoi(argv[1]));
    a.Loop(outRootFileF);
  }
  else if(argc == 4 && atoi(argv[1])==2){
    TString rootFilesList = argv[2];
    TString outRootFileF = argv[3];
    cout<<"--> Parameter calculation from the WF <--"<<endl
	<<"rootFilesList : "<<rootFilesList<<endl
	<<"outRootFileF  : "<<outRootFileF<<endl;
    anashort a(rootFilesList);
    a.Loop(outRootFileF);
  }
  else if(argc == 4 && atoi(argv[1])==3){
    TString inRootFiles = argv[2];
    TString outRootFileF = argv[3];
    cout<<"--> Parameter calculation from the WF <--"<<endl
	<<"inRootFiles   : "<<inRootFiles<<endl
	<<"outRootFileF  : "<<outRootFileF<<endl;
    anashort a( inRootFiles, atoi(argv[1]));
    a.Loop(outRootFileF);
  }
  else if(argc == 5 && atoi(argv[1])==4){
    TString inRootFiles = argv[2];
    TString outRootFileF = argv[3];
    Long64_t evID = (Long64_t)atoi(argv[4]);
    cout<<"--> Parameter calculation from the WF <--"<<endl
	<<"inRootFiles   : "<<inRootFiles<<endl
	<<"outRootFileF  : "<<outRootFileF<<endl
      	<<"evID          : "<<evID<<endl;
    ana a(inRootFiles, 1);
    a.save_wf_for_event(outRootFileF, evID);
  }
  else if(argc == 2 && atoi(argv[1])==5){
    cout<<"--> Parameter calculation from the WF <--"<<endl
	<<"argv[1] : "<<atoi(argv[1])<<endl;
    /*
    sipmCameraHist *sipm_cam = new sipmCameraHist("sipm_cam","sipm_cam","pixel_mapping.csv",0);
    sipm_cam->test();
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
    */
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
  }
  return 0;
}
