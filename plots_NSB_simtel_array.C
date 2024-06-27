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

void read_file(TString name,TH1D *h1);
void copyHist(TH1D *h1, TH1D *h1_to_cp);
void norm_hist(TH1D *h1, Double_t norm);

Int_t plots_NSB_simtel_array(){
  ////////////////////
  TString fileN0101;
  TString fileN0102;
  TString fileN0103;
  TString fileN0201;
  TString fileN0202;
  TString fileN0203;
  ////////////////////
  //
  fileN0101 = "../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/nsb_1x_386MHz/trgA_test/hist_trgA_corsika_0000ID.root";
  fileN0102 = "../DBscan_on_simtel_data/hist_wf_val_NSB386MHz.csv";
  fileN0103 = "../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/trgA_test_386MHz/hist_trgA_corsika_0binE_0binTheta_0binDist_0000ID.root";
  //
  fileN0201 = "../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/nsb_1x_268MHz/trgA_test/hist_trgA_corsika_0000ID.root";
  fileN0202 = "../DBscan_on_simtel_data/hist_wf_val_NSB268MHz.csv";
  fileN0203 = "../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/trgA_test_268MHz/hist_trgA_corsika_0binE_0binTheta_0binDist_0000ID.root";
  ////////////////////
  TFile *f0101 = new TFile(fileN0101.Data());
  TFile *f0103 = new TFile(fileN0103.Data());
  TFile *f0201 = new TFile(fileN0201.Data());
  TFile *f0203 = new TFile(fileN0203.Data());
  ////////////////////
  TH1D *h1_fadc_val_386MHz = (TH1D*)f0101->Get("h1_fadc_val");
  TH1D *h1_fadc_val_386MHz_proton = (TH1D*)f0103->Get("h1_fadc_val");
  TH1D *h1_fadc_val_simtel_386MHz = new TH1D( "h1_fadc_val_simtel_386MHz", "h1_fadc_val_simtel_386MHz", (Int_t)(400.5-249.5), 249.5, 400.5);
  //
  TH1D *h1_fadc_val_268MHz = (TH1D*)f0201->Get("h1_fadc_val");
  TH1D *h1_fadc_val_268MHz_proton = (TH1D*)f0203->Get("h1_fadc_val");
  TH1D *h1_fadc_val_simtel_268MHz = new TH1D( "h1_fadc_val_simtel_268MHz", "h1_fadc_val_simtel_268MHz", (Int_t)(400.5-249.5), 249.5, 400.5);
  //
  TH1D *h1_fadc_val_386MHz_norm = new TH1D();
  h1_fadc_val_386MHz_norm->SetNameTitle("h1_fadc_val_386MHz_norm","h1_fadc_val_386MHz_norm");
  TH1D *h1_fadc_val_386MHz_proton_norm = new TH1D();
  h1_fadc_val_386MHz_proton_norm->SetNameTitle("h1_fadc_val_386MHz_proton_norm","h1_fadc_val_386MHz_proton_norm");
  TH1D *h1_fadc_val_simtel_386MHz_norm = new TH1D();
  h1_fadc_val_simtel_386MHz_norm->SetNameTitle("h1_fadc_val_simtel_386MHz_norm","h1_fadc_val_simtel_386MHz_norm");
  TH1D *h1_fadc_val_268MHz_norm = new TH1D();
  h1_fadc_val_268MHz_norm->SetNameTitle("h1_fadc_val_268MHz_norm","h1_fadc_val_268MHz_norm");
  TH1D *h1_fadc_val_268MHz_proton_norm = new TH1D();
  h1_fadc_val_268MHz_proton_norm->SetNameTitle("h1_fadc_val_268MHz_proton_norm","h1_fadc_val_268MHz_proton_norm");
  TH1D *h1_fadc_val_simtel_268MHz_norm = new TH1D();
  h1_fadc_val_simtel_268MHz_norm->SetNameTitle("h1_fadc_val_simtel_268MHz_norm","h1_fadc_val_simtel_268MHz_norm");
  //
  ////////////////////
  TCanvas *c1 = new TCanvas("c1",fileN0101.Data(), 10, 10, 800, 800);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
  ////////////////////
  ////////////////////
  read_file(fileN0102, h1_fadc_val_simtel_386MHz);
  read_file(fileN0202, h1_fadc_val_simtel_268MHz);
  ////////////////////
  ////////////////////
  copyHist(h1_fadc_val_386MHz_norm, h1_fadc_val_386MHz);  
  copyHist(h1_fadc_val_386MHz_proton_norm, h1_fadc_val_386MHz_proton);
  copyHist(h1_fadc_val_simtel_386MHz_norm, h1_fadc_val_simtel_386MHz);
  copyHist(h1_fadc_val_268MHz_norm, h1_fadc_val_268MHz);
  copyHist(h1_fadc_val_268MHz_proton_norm, h1_fadc_val_268MHz_proton);
  copyHist(h1_fadc_val_simtel_268MHz_norm, h1_fadc_val_simtel_268MHz);
  norm_hist(h1_fadc_val_386MHz_norm, h1_fadc_val_386MHz->Integral());
  norm_hist(h1_fadc_val_386MHz_proton_norm, h1_fadc_val_386MHz_proton->Integral());
  norm_hist(h1_fadc_val_simtel_386MHz_norm, h1_fadc_val_simtel_386MHz->Integral());
  norm_hist(h1_fadc_val_268MHz_norm, h1_fadc_val_268MHz->Integral());
  norm_hist(h1_fadc_val_268MHz_proton_norm, h1_fadc_val_268MHz_proton->Integral());
  norm_hist(h1_fadc_val_simtel_268MHz_norm, h1_fadc_val_simtel_268MHz->Integral());  
  //
  ////////////////////
  h1_fadc_val_386MHz->SetLineColor(kBlack);
  h1_fadc_val_386MHz->SetLineWidth(3);
  h1_fadc_val_simtel_386MHz->SetLineColor(kRed);
  h1_fadc_val_simtel_386MHz->SetLineWidth(3);
  h1_fadc_val_268MHz->SetLineColor(kBlack);
  h1_fadc_val_268MHz->SetLineWidth(3);
  h1_fadc_val_simtel_268MHz->SetLineColor(kRed);
  h1_fadc_val_simtel_268MHz->SetLineWidth(3);
  //
  h1_fadc_val_386MHz_norm->SetLineColor(kBlack);
  h1_fadc_val_386MHz_norm->SetLineWidth(3);
  h1_fadc_val_simtel_386MHz_norm->SetLineColor(kRed);
  h1_fadc_val_simtel_386MHz_norm->SetLineWidth(3);
  h1_fadc_val_268MHz_norm->SetLineColor(kBlack);
  h1_fadc_val_268MHz_norm->SetLineWidth(3);
  h1_fadc_val_simtel_268MHz_norm->SetLineColor(kRed);
  h1_fadc_val_simtel_268MHz_norm->SetLineWidth(3);
  //
  h1_fadc_val_386MHz_proton_norm->SetLineColor(kBlue);
  h1_fadc_val_386MHz_proton_norm->SetLineWidth(3);
  //
  h1_fadc_val_268MHz_proton_norm->SetLineColor(kBlue);
  h1_fadc_val_268MHz_proton_norm->SetLineWidth(3);
  ////////////////////
  //
  //h1_fadc_val_386MHz->Draw();
  //h1_fadc_val_simtel_386MHz->Draw("same");
  //
  //h1_fadc_val_386MHz_norm->Draw();
  //h1_fadc_val_simtel_386MHz_norm->Draw("sames");
  //h1_fadc_val_386MHz_proton_norm->Draw("sames");
  //
  h1_fadc_val_268MHz_norm->Draw();
  h1_fadc_val_simtel_268MHz_norm->Draw("sames");
  //h1_fadc_val_268MHz_proton_norm->Draw("sames");
  //
  //h1_fadc_val_268MHz->Draw();
  //h1_fadc_val_386MHz->Draw("same");
  //
  //h1_fadc_val_simtel_268MHz->Draw();
  //h1_fadc_val_simtel_386MHz->Draw("same");
  //
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(h1_fadc_val_386MHz, "Stand alone simulation", "apl");
  leg->AddEntry(h1_fadc_val_simtel_386MHz, "sim_telarray ", "apl");
  leg->Draw();  
  //
  //
  /*
  //
  //
  h1_eff_area_01->SetTitle("");
  h1_eff_area_01->Draw();
  h1_eff_area_02->Draw("same");
  h1_eff_area_01->SetMinimum(1.0e+4);
  h1_eff_area_01->SetMinimum(1.0e+7);
  //
  //h1_protons->SetTitle("");
  //h1_protons->Draw();
  h1_protons02->Draw("same");
  //
  //h1_protons->Draw("same");
  h1_NSB->Draw("same");
  //
  //h1_tot->SetMinimum(1.0e+3);
  //
  h1_tot->GetXaxis()->SetTitle("Number of points in the cluster");
  h1_tot->GetYaxis()->SetTitle("Rate, Hz");
  //
  //TString file_out = fileN01;
  //file_out += ".pdf";
  //c1->SaveAs(file_out.Data());
  //
  */
  return 0;
}

void read_file( TString name, TH1D *h1){
  std::ifstream fFile(name.Data());
  std::string mot;
  Double_t xbin, val;
  Int_t i = 1;
  if(fFile.is_open()){
    fFile>>mot>>mot;
    //cout<<"mot "<<mot<<endl;
    while(fFile>>xbin>>val){
      //cout<<val<<endl;
      h1->SetBinContent(i,val);
      i++;
    }
    fFile.close();
  }
}

void copyHist(TH1D *h1, TH1D *h1_to_cp){
  int nBins = h1_to_cp->GetNbinsX();
  double *bins_low_edge = new double[nBins+1];
  for(int i = 1;i<=nBins;i++)
    bins_low_edge[i-1] = h1_to_cp->GetBinLowEdge(i);
  bins_low_edge[nBins] = h1_to_cp->GetBinLowEdge(nBins) + h1_to_cp->GetBinWidth(nBins);
  h1->SetBins(nBins,bins_low_edge);
  //  
  for(Int_t i = 1;i<=h1_to_cp->GetNbinsX();i++)
    h1->SetBinContent(i,h1_to_cp->GetBinContent(i));
}

void norm_hist(TH1D *h1, Double_t norm = 1){
  for(Int_t i = 1;i<=h1->GetNbinsX();i++){
    h1->SetBinContent(i,h1->GetBinContent(i)/norm);
  }
}
