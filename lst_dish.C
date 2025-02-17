//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include <time.h>

TGraphErrors *_gr_z_vs_r;

Double_t dish_parabola(Double_t x, Double_t *par);
void fcn(int &npar, double *gin, double &f, double *par, int iflag);

Int_t lst_dish(){
  //
  TString fileN01;
  fileN01 = "./hist_lstMirrorHist_test_ideal.root";
  //
  TFile *f01 = new TFile(fileN01.Data());
  _gr_z_vs_r = (TGraphErrors*)f01->Get("gr_z_vs_r");
  //
  TGraph *gr_fit = new TGraph();
  gr_fit->SetNameTitle("gr_fit","gr_fit");
  ///////////////////////////
  //Fit
  const Int_t npar = 3;
  Double_t x_min = -100.0;
  Double_t x_max =  1200.0;
  Double_t inParameters[npar];
  Double_t outParameters[npar];
  Double_t outParametersError[npar];
  inParameters[0] = 0.0;
  inParameters[1] = 0.0;
  inParameters[2] = 90.0/1000.0/1000.0;
  //
  TMinuit *minuit_m = new TMinuit(npar);
  //gMinuit->SetPrintLevel(-1.0);
  minuit_m->SetFCN(fcn); 
  double arglist[10];
  int ierflg = 0;
  arglist[0] = 1;
  minuit_m->mnexcm("SET ERR", arglist ,1,ierflg);
  // 
  // Set starting values and step sizes for parameters
  //   3  k  8.55621e-05   3.98315e-05  -3.27795e-06  -4.18921e-04
  //   3  k  8.91666e-05   3.81152e-06   3.49399e-09   1.27101e-06
  //   3  k  8.91714e-05   9.39948e-08   3.17300e-09   7.16589e-04
  minuit_m->mnparm(0, "r0", inParameters[0], 0.0, 0,0,ierflg);
  minuit_m->mnparm(1, "l",  inParameters[1], 0.0, 0,0,ierflg);
  minuit_m->mnparm(2, "k",  inParameters[2], 0.0001, 0,0,ierflg);
  //
  
  // Now ready for minimization step
  arglist[0] = 500000;
  arglist[1] = 1.;
  minuit_m->mnexcm("MIGRAD", arglist ,2,ierflg);
  
  // Print results
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  minuit_m->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  //gMinuit->mnprin(3,amin);
  //
  for(Int_t i = 0;i<npar;i++)  
    minuit_m->GetParameter(i, outParameters[i], outParametersError[i]);

  double x_tmp;
  for(Int_t i = 0;i<1000;i++){  
    x_tmp = x_min + (x_max - x_min)/(1000-1)*i;
    gr_fit->SetPoint(i,x_tmp,dish_parabola(x_tmp,outParameters));
  }
  //
  //cout<<x0out<<endl
  //  <<y0out<<endl
  //  <<Rout<<endl;
  //  


  //
  //TF1 *f_dish = new TF1( "dish_parabola", dish_parabola, x_min, x_max, npar);
  //TF1 *f3 = gre3->GetFunction("pol3");
  //f_dish->SetParameters(inParameters);
  //f_dish->SetParName(0, "r0");
  //f_dish->SetParName(1, "k");
  //f_dish->FixParameter(0,inParameters[0]);
  //gr_z_vs_r->Fit("dish_parabola","","",x_min, x_max);
  //_gr_z_vs_r->Fit("pol1");
  ///////////////////////////  
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  //gPad->SetLogx();
  // gPad->SetLogy();
  //
  //_gr_z_vs_r->Fit("pol2", "F");
  _gr_z_vs_r->SetTitle("");
  //_gr_z_vs_r->Draw("AP");
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(_gr_z_vs_r);
  mg->Add( gr_fit);
  //
  mg->Draw("ap");
  //
  //mg->GetXaxis()->SetTitle("Time, a.u.");
  //mg->GetYaxis()->SetTitle("FADC counts");
  //  
  //TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  //leg->AddEntry(gr, "simtelarr", "apl");
  //leg->AddEntry(gr_sim, "sim", "apl");
  //leg->Draw();  
  //_gr_z_vs_r->Draw("AP");
  //f_dish  h1_E_proton->GetXaxis()->SetTitle("E, TeV");
  cout<<" f = "<<1.0/4.0/outParameters[2]<<endl;
  return 0;
}

Double_t dish_parabola(Double_t x, Double_t *par){
  return x*x*par[2] + x*par[1] + par[0];
}

void fcn(int &npar, double *gin, double &f, double *par, int iflag){
  double chisq = 0.0;
  double x, y;
  double delta;
  double xerr, yerr;
  for (int i = 0; i<_gr_z_vs_r->GetN(); i++){
    _gr_z_vs_r->GetPoint( i, x, y);
    xerr = _gr_z_vs_r->GetErrorX(i);
    yerr = _gr_z_vs_r->GetErrorY(i);    
    //delta = (y - dish_parabola( x, par))/TMath::Sqrt(xerr*xerr + yerr*yerr);
    delta = (y - dish_parabola( x, par));
    chisq += delta*delta;
  }
  f = chisq;
}
