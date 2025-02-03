//my
#include "muonRingFitter.hh"

//root
#include <TMinuit.h>
#include <TGraph.h>

//C, C++
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <fstream>
#include <iomanip>
#include <vector>
#include <time.h>
#include <bits/stdc++.h>
#include <sys/stat.h>

using namespace std;

TGraph* muonRingFitter::_gr_to_fit = new TGraph();
Double_t muonRingFitter::amplifySmallerOrBiggerRadius = 0.5;

void muonRingFitter::set_gr(TGraph *gr){
  double x, y;
  for (int i = 0; i<_gr_to_fit->GetN(); i++)
    _gr_to_fit->SetPoint( i, -999.0, -999.0);
  for (int i = 0; i<gr->GetN(); i++){
     gr->GetPoint( i, x, y);
    _gr_to_fit->SetPoint( i, x, y);
  }
}

muonRingFitter::muonRingFitter(Int_t npar) : TMinuit(npar)
{
  SetPrintLevel(-1);  
}

muonRingFitter::~muonRingFitter(){
}

double muonRingFitter::equation_of_circle(double x, double y, double *par){
  return par[2]*par[2] - (x-par[0])*(x-par[0]) - (y-par[1])*(y-par[1]);
}

void muonRingFitter::fcn_ring(int &npar, double *gin, double &f, double *par, int iflag){
  double chisq = 0.0;
  double x, y;
  double delta;
  for (int i = 0; i<_gr_to_fit->GetN(); i++){
    _gr_to_fit->GetPoint( i, x, y);
    if(x>-999.0 && y>-999.0){
      delta = equation_of_circle(x, y, par);
      if(delta>0.0)
	delta = delta*amplifySmallerOrBiggerRadius;
      chisq += TMath::Abs(delta);
    }
  }
  f = chisq;
}

void muonRingFitter::fit_ring(Double_t x0in, Double_t y0in, Double_t Rin,
			      Double_t &x0out, Double_t &y0out, Double_t &Rout,
			      Double_t &x0outerr, Double_t &y0outerr, Double_t &Routerr){
  //  
  SetFCN(fcn_ring); 
  //
  double arglist[10];
  int ierflg = 0;
  arglist[0] = 1;
  mnexcm("SET ERR", arglist ,1,ierflg);
  // 
  //Set starting values and step sizes for parameters
  mnparm(0, "x0", x0in, 0.001, 0,0,ierflg);
  mnparm(1, "y0", y0in, 0.001, 0,0,ierflg);
  mnparm(2, "R", Rin, 0.001, 0,0,ierflg);
  //
  //Now ready for minimization step
  arglist[0] = 100000;
  arglist[1] = 1.0;
  mnexcm("MIGRAD", arglist ,2,ierflg);
  //
  GetParameter(0, x0out, x0outerr);
  GetParameter(1, y0out, y0outerr);
  GetParameter(2, Rout, Routerr);
  //
}

