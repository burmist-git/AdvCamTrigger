//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include <time.h>

using namespace std;

Int_t plots_stareo_dt_fit(){

  TString fileN01;
  //
  fileN01 = "./hist_proton_st.root";
  //
  TFile *f01 = new TFile(fileN01.Data());
  //
  TH1D *h1_01 = (TH1D*)f01->Get("h1_dtime_LST1_m_LST2");
  //h1_01->Rebin(4);
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE); 
  //
  h1_01->SetLineColor(kBlack);
  h1_01->SetLineWidth(3.0);
  //
  h1_01->SetTitle("");
  //
  //h1_01->SetMaximum(1.13);
  h1_01->Draw();
  //
  Double_t xmean = h1_01->GetMean();
  Double_t threerms = h1_01->GetRMS()*3.0;
  cout<<"h1_01->GetMean() = "<<h1_01->GetMean()<<endl;
  //
  // Three TF1 objects are created, one for each subrange.
  TF1 *g1 = new TF1("g1", "gaus", xmean - 1.5*threerms, xmean + 1.5*threerms);
  TF1 *g2 = new TF1("g2", "gaus", xmean - 1.5*threerms, xmean + 1.5*threerms);

  TF1 *total = new TF1("total", "gaus(0)+gaus(3)", xmean - 1.5*threerms, xmean + 1.5*threerms);

  total->SetParName(0,"A1");
  total->SetParName(1,"mean1");
  total->SetParName(2,"sigma1");
  total->SetParName(3,"A2");
  total->SetParName(4,"mean2");
  total->SetParName(5,"sigma2");
  
  
  double par[6];
  
  par[0] = 2.96607e+03;
  par[1] = -7.51993e+01;
  par[2] = 1.34863e+01;

  par[3] = 2.96607e+03/3.0;
  par[4] = -7.51993e+01;
  par[5] = 1.34863e+01*2.0;
    
  
  //g1->GetParameters(&par[0]);
  //g2->GetParameters(&par[3]);
  
  // Use the parameters on the sum.
  total->SetParameters(par);
  h1_01->Fit(total, "R+","xmean - 1.5*threerms, xmean + 1.5*threerms");

  cout<<total->GetParameter(0)<<endl;

  g1->SetParameter(0,  total->GetParameter(0));
  g1->SetParameter(1,  total->GetParameter(1));
  g1->SetParameter(2,  total->GetParameter(2));
  g2->SetParameter(0,  total->GetParameter(3));
  g2->SetParameter(1,  total->GetParameter(4));
  g2->SetParameter(2,  total->GetParameter(5));

  g1->SetLineColor(kBlue+2);
  g1->Draw("sames");

  g2->SetLineColor(kMagenta+2);
  g2->Draw("sames");
  
  /*

  
g2->SetParameters(par[3]);

  */  
   
  
  //
  //h1_01->GetXaxis()->SetTitle("Number of photons on the telescope sphere");
  //gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  //
  /*
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(h1_01, "LST1 - LST2", "apl");
  leg->AddEntry(h1_02, "LST1 - LST3", "apl");
  leg->AddEntry(h1_03, "LST1 - LST4", "apl");
  leg->AddEntry(h1_04, "LST2 - LST3", "apl");
  leg->AddEntry(h1_05, "LST2 - LST4", "apl");
  leg->AddEntry(h1_06, "LST3 - LST4", "apl");
  leg->Draw();
  //
  //LST1 -> LST2 -75.1
  //LST1 -> LST3 -212.2
  //LST1 -> LST4 -150.7
  //LST2 -> LST3 -137.1
  //LST2 -> LST4 -75.6
  //LST3 -> LST4  61.5
  */
  //  
  return 0;
}
