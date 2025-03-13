#ifndef anabase_hh
#define anabase_hh

#include <TROOT.h>

#include <vector>

class TChain;
class TFile;
class TTree;
class TString;
class TBranch;
class TGraph;
class TH1D;
class TH2D;
class TProfile;

class anabase {

public :
  anabase(TString fileList, Bool_t short_format_flag);
  anabase(TString inFileName, Int_t keyID, Bool_t short_format_flag);
  anabase(TString fileListm);
  anabase(TString inFileName, Int_t keyID);
  ~anabase();
  TString _particle_type_name;
  inline void set_particle_type_name(TString particle_type_name){_particle_type_name = particle_type_name;};
  void printEv();  
  void getCore_rel_R_theta(const Double_t x0, const Double_t y0, const Double_t xx, const Double_t yy, Double_t &rr, Double_t &theta);
  Int_t get_current_data_chunk_ID(Long64_t nentries, Long64_t jentry);
  
public :
  void Loop(TString histOut);
  Long64_t LoadTree(Long64_t entry);
  
public :
  Bool_t _short_format_flag;
  Int_t GetEntry(Long64_t entry);
  void Init(TTree *tree);
  Bool_t Notify();
  void Show(Long64_t entry = -1);
  Int_t Cut(Long64_t entry);
  
public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  //Int_t           evt;
  //Int_t           run;
  //Float_t         pValue;
  //...
  //...
  //
  //---------------------------------------------------
  // ADD HERE :
  //Tree name
  //const TString treeName = "arich";
  const TString treeName = "T";
  //
  ////////////////////////////////////
  Int_t   event_id;
  Float_t energy;
  //
  Float_t azimuth;
  Float_t altitude;
  Float_t h_first_int;
  Float_t xmax;
  Float_t hmax;
  Float_t emax;
  Float_t cmax;
  //
  Float_t xcore;
  Float_t ycore;
  Float_t ev_time;
  Int_t   nphotons;
  Int_t   n_pe;
  Int_t   n_pixels;
  Int_t   pe_chID[946466];
  Float_t pe_time[946466];
  Int_t   wf[7987][75];
  ////////////////////////////////////  
  static const Int_t nChannels = 7987;
  static const Int_t nn_fadc_point = 75;

  Int_t _n_data_chunks;
  Int_t _data_chunk_ID;
  
  //
  //---------------------------------------------------
  
  // List of branches
  //TBranch        *b_evt;
  //TBranch        *b_run;
  //TBranch        *b_pValue;
  //...
  //...
  //
  //---------------------------------------------------
  // ADD HERE :
  TBranch *b_event_id;   //!
  TBranch *b_energy;   //!
  //
  TBranch *b_azimuth;   //!
  TBranch *b_altitude;   //!
  TBranch *b_h_first_int;   //!
  TBranch *b_xmax;   //!
  TBranch *b_hmax;   //!
  TBranch *b_emax;   //!
  TBranch *b_cmax;   //!
  //
  TBranch *b_xcore;   //!
  TBranch *b_ycore;   //!
  TBranch *b_ev_time;   //!
  TBranch *b_nphotons;   //!
  TBranch *b_n_pe;   //!
  TBranch *b_n_pixels;   //!
  TBranch *b_pe_chID;   //!
  TBranch *b_pe_time;   //!
  TBranch *b_wf;   //!
  //---------------------------------------------------
  void tGraphInit(TGraph *gr[nChannels], TString grName, TString grTitle);
  void h1D1Init(TH1D *h1D1[nChannels],TString h1name, TString h1Title,
		Int_t Nbin, Float_t Vmin, Float_t Vmax);
  void h1D1Init(std::vector<TH1D*> &h1D1_v, unsigned int n_len, TString h1name, TString h1Title,
		Int_t Nbin, Float_t Vmin, Float_t Vmax);
  void h2D2Init(TH2D *h2D1[nChannels], TString h2name, TString h2Title,
                Int_t Nbin1, Float_t Vmin1, Float_t Vmax1,
                Int_t Nbin2, Float_t Vmin2, Float_t Vmax2);
  void h2D2Init(std::vector<TH2D*> &h2D1_v, unsigned int n_len, TString h2name, TString h2Title,
                Int_t Nbin1, Float_t Vmin1, Float_t Vmax1,
                Int_t Nbin2, Float_t Vmin2, Float_t Vmax2);
  void tProfInit(TProfile *tprof[nChannels],TString prname, TString prTitle,
                 Int_t Nbin, Float_t Vmin, Float_t Vmax);
  void tProfInit(std::vector<TProfile*> &tprof_v, unsigned int n_len,
		 TString prname, TString prTitle, Int_t Nbin, Float_t Vmin, Float_t Vmax);
  double getUnixTimeFromTime(double d_year, double d_month, double d_day, double d_hour, double d_min, double d_sec);  
  //
  static void TH2D_divide( TH2D *h2_w, TH2D *h2,TH2D *h2_norm);
  static void TH1D_divide( TH1D *h1_w, TH1D *h1,TH1D *h1_norm);
  static void TH1D_divide( TH1D *h1, TH1D *h1_norm, Double_t norm);
  //
};

#endif
