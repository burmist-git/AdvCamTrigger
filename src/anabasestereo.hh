#ifndef anabasestereo_hh
#define anabasestereo_hh

#include <TROOT.h>
#include <TVector3.h>

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
class TVector3;

class anabasestereo {

public :
  anabasestereo(TString fileList);
  anabasestereo(TString inFileName, Int_t keyID);
  ~anabasestereo();
  TString _particle_type_name;
  inline void set_particle_type_name(TString particle_type_name){_particle_type_name = particle_type_name;};
  static Double_t getExpectedTimeDelayBetweenTwoLST(const TVector3 &v1, const TVector3 &v2,
						    Double_t tel_arr_pointing_azimuth, Double_t tel_arr_pointing_altitude,
						    Double_t hmax_approximate_m = 10000.0);
  Int_t _n_file_counter;
  //
  Bool_t if_single_LST(Int_t &coincidenceID);
  Bool_t if_double_LST(Int_t &coincidenceID);
  Bool_t if_triple_LST(Int_t &coincidenceID);
  Bool_t if_four_LST(Int_t &coincidenceID);
  Bool_t if_stereo_trigger_pne(Int_t npecut);
  
public :
  void Loop(TString histOut);
  Long64_t LoadTree(Long64_t entry);
  
public :
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
  static const Int_t _npeMax = 1000000;
  static const Int_t nChannels = 7987;
  static const Int_t nn_fadc_point = 75;
  //
  const TVector3 LST1_r0;
  const TVector3 LST2_r0;
  const TVector3 LST3_r0;
  const TVector3 LST4_r0;
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
  Float_t xcore;
  Float_t ycore;
  //
  Float_t ev_time_LST1;
  Float_t ev_time_LST2;
  Float_t ev_time_LST3;
  Float_t ev_time_LST4;
  Int_t   nphotons_LST1;
  Int_t   nphotons_LST2;
  Int_t   nphotons_LST3;
  Int_t   nphotons_LST4;
  Int_t   n_pe_LST1;
  Int_t   n_pe_LST2;
  Int_t   n_pe_LST3;
  Int_t   n_pe_LST4;
  Int_t   n_pixels_LST1;
  Int_t   n_pixels_LST2;
  Int_t   n_pixels_LST3;
  Int_t   n_pixels_LST4;
  Int_t   pe_chID_LST1[_npeMax]; //[n_pe_LST1]
  Int_t   pe_chID_LST2[_npeMax]; //[n_pe_LST2]
  Int_t   pe_chID_LST3[_npeMax]; //[n_pe_LST3]
  Int_t   pe_chID_LST4[_npeMax]; //[n_pe_LST4]
  Float_t pe_time_LST1[_npeMax]; //[n_pe_LST1]
  Float_t pe_time_LST2[_npeMax]; //[n_pe_LST2]
  Float_t pe_time_LST3[_npeMax]; //[n_pe_LST3]
  Float_t pe_time_LST4[_npeMax]; //[n_pe_LST4]
  ////////////////////////////////////    
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
  TBranch *b_event_id;      //!
  TBranch *b_energy;        //!
  //
  TBranch *b_azimuth;       //!
  TBranch *b_altitude;      //!
  TBranch *b_h_first_int;   //!
  TBranch *b_xmax;          //!
  TBranch *b_hmax;          //!
  TBranch *b_emax;          //!
  TBranch *b_cmax;          //!
  TBranch *b_xcore;         //!
  TBranch *b_ycore;         //!
  //
  TBranch *b_ev_time_LST1;  //!
  TBranch *b_ev_time_LST2;  //!
  TBranch *b_ev_time_LST3;  //!
  TBranch *b_ev_time_LST4;  //!
  TBranch *b_nphotons_LST1; //!
  TBranch *b_nphotons_LST2; //!
  TBranch *b_nphotons_LST3; //!
  TBranch *b_nphotons_LST4; //!
  TBranch *b_n_pe_LST1;     //!
  TBranch *b_n_pe_LST2;     //!
  TBranch *b_n_pe_LST3;     //!
  TBranch *b_n_pe_LST4;     //!
  TBranch *b_n_pixels_LST1; //!
  TBranch *b_n_pixels_LST2; //!
  TBranch *b_n_pixels_LST3; //!
  TBranch *b_n_pixels_LST4; //!
  TBranch *b_pe_chID_LST1;  //!
  TBranch *b_pe_chID_LST2;  //!
  TBranch *b_pe_chID_LST3;  //!
  TBranch *b_pe_chID_LST4;  //!
  TBranch *b_pe_time_LST1;  //!
  TBranch *b_pe_time_LST2;  //!
  TBranch *b_pe_time_LST3;  //!
  TBranch *b_pe_time_LST4;  //!
  //---------------------------------------------------
};

#endif
