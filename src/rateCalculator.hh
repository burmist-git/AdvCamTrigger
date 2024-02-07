#pragma once

//root
#include <TString.h>

//c, c++
#include <string>
#include <vector>

class evstHist;
class TH1D;

class rateCalculator {
public:
  
  rateCalculator();
  rateCalculator( const char* name, const char* title, TString hist_file_prefix, TString output_hist_file, evstHist* evH_flux);
  ~rateCalculator();
  
private:
  
  TString _name;
  TString _title;
  TString _hist_file_prefix;
  TString _output_hist_file;
  evstHist *_evH_flux;
  
  bool get_histogram_file_name(TString hist_file_prefix, Int_t i_binE, Int_t i_binTheta, Int_t i_binDist, TString &hist_file_name);
  void calculate_rate();
  void append_histogram(TH1D *h1);
  void get_histogram(TString hist_file_name, TH1D *h1);

  void copyHist(TH1D *h1, TH1D *h1_to_cp);
  void addHist(TH1D *h1, TH1D *h1_to_add);
  
  TH1D *_h1_rate;
  std::vector<TH1D*> _h1_rates_per_bin_v;

  void save_output();
  
};
