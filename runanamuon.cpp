//my
#include "src/anamuon.hh"

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
    //ana a(rootFilesList);
    //a.Loop(outRootFileF);
  }
  else if(argc == 4 && atoi(argv[1])==1){
    TString inRootFile = argv[2];
    TString outRootFileHist = argv[3];
    cout<<"--> Parameters <--"<<endl
	<<"inRootFile      : "<<inRootFile<<endl
	<<"outRootFileHist  : "<<outRootFileHist<<endl;
    anamuon a( inRootFile, atoi(argv[1]));
    a.Loop(outRootFileHist);
  }
  else if(argc == 6 && atoi(argv[1])==2){
    TString inRootFile = argv[2];
    TString outRootFileHist = argv[3];
    Long64_t evID = (Long64_t)atoi(argv[4]);
    TString particle_type_name = argv[5];
    cout<<"--> Parameters <--"<<endl
	<<"inRootFile         : "<<inRootFile<<endl
	<<"outRootFileHist    : "<<outRootFileHist<<endl
      	<<"evID               : "<<evID<<endl
	<<"particle_type_name : "<<particle_type_name<<endl;
    anamuon a(inRootFile, 1);
    a.set_particle_type_name(particle_type_name);
    a.save_wf_for_event(outRootFileHist, evID);
  }
  else{
    cout<<" --> ERROR in input arguments "<<endl
	<<" runID [1] = 0 (execution ID number)"<<endl
      	<<"       [2] - file with list of the root files"<<endl
	<<"       [3] - name of root file with histograms"<<endl;
    cout<<" runID [1] = 1 (execution ID number)"<<endl
      	<<"       [2] - in root file"<<endl
	<<"       [3] - name of root file with histograms"<<endl;
    cout<<" runID [1] = 2 (execution ID number) short file format"<<endl
      	<<"       [2] - in root file"<<endl
	<<"       [3] - name of root file with histograms"<<endl
    	<<"       [4] - evID"<<endl
	<<"       [5] - particle_type_name"<<endl;
  }
  return 0;
}
