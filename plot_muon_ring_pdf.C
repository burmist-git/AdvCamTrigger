//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>

#include <time.h>

using namespace std;

Int_t plot_muon_ring_pdf(){
  //
  TString fileN01;
  fileN01 = "../scratch/simtel_data/muon/hist/hist_run1_muon.root";
  //
  TFile *f01 = new TFile(fileN01.Data());
  //
  Int_t n_pdf = 10000;
  //gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=082_Reading_NTC_SiPM_Tile_PCB_all.pdf 082_Reading_NTC_SiPM_Tile_PCB.pdf TDK_NTCG_series.pdf
  ofstream outfile;
  outfile.open("merg_muon_pdf.sh");
  outfile<<"gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=c1_all.pdf ";
  //
  for(Int_t i = 0;i<n_pdf;i++){
    TCanvas *c1 = NULL;
    TString c1name = "trueRing/c1_";
    c1name += i;
    c1 = (TCanvas*)f01->Get(c1name.Data());
    if(c1 == NULL){
      //cout<<"c1 == NULL for i = "<<i<<endl;
    }
    else{
      //cout<<"c1 != NULL"<<endl;
      c1name += ".pdf";
      c1->SaveAs(c1name.Data());
      outfile<<c1name<<" ";
    }
  }
  outfile.close();
  //
  return 0;
}
