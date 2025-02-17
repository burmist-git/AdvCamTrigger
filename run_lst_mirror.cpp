//my
#include "src/lstMirrorHist.hh"

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
  if(argc == 2 && atoi(argv[1])==0){
    cout<<"--> lstMirrorHist <--"<<endl;
    //lstMirrorHist *lstmirr = new lstMirrorHist(1);
    //
    lstMirrorHist *lstmirr = new lstMirrorHist();
    lstmirr->dump_mapping_info();
    lstmirr->test();
  }
  else if (argc == 2 && atoi(argv[1])==1){
    lstMirrorHist *lstmirr = new lstMirrorHist( "lstMirrorHist_ideal", "lstMirrorHist_ideal", "mirror_CTA-LST-v20141201-198.dat", true, 0.0);
    lstmirr->dump_mapping_info();
    lstmirr->test_ideal();
  }
  else{
    cout<<" --> ERROR in input arguments "<<endl
	<<" runID [1] = 0 (lstMirrorHist)"<<endl
    	<<" runID [1] = 1 (lstMirrorHist ideal)"<<endl;
  }
  return 0;
}
