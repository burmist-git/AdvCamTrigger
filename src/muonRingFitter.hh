#ifndef muonRingFitter_hh
#define muonRingFitter_hh

//root
#include <TROOT.h>
#include <TMinuit.h>
#include <TGraph.h>

//C, C++
#include <vector>

class TString;
class TGraph;

class muonRingFitter: public TMinuit {
public:

  muonRingFitter(Int_t npar = 3);
  ~muonRingFitter();
  
  void fit_ring(Double_t x0in, Double_t y0in, Double_t Rin,
  		Double_t &x0out, Double_t &y0out, Double_t &Rout,
  		Double_t &x0outerr, Double_t &y0outerr, Double_t &Routerr);

  static void fcn_ring(int &npar, double *gin, double &f, double *par, int iflag);
  static double equation_of_circle(double x, double y, double *par);

  static TGraph* _gr_to_fit;
  static void set_gr(TGraph *gr);

  static Double_t amplifySmallerOrBiggerRadius;
  
};

#endif
