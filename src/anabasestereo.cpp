//my
#include "anabasestereo.hh"

//root
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TVector2.h>
#include <TVector3.h>

//C, C++
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <bits/stdc++.h>
#include <vector>

using namespace std;

anabasestereo::anabasestereo(TString fileList) : _particle_type_name("NONE"),
						 LST1_r0(-70.93, -52.07, 43.0),
						 LST2_r0(-35.27,  66.14, 32.0),
						 LST3_r0( 75.28,  50.49, 28.7),
						 LST4_r0( 30.91, -64.54, 32.0),
						 _n_file_counter(-1)
{
  ifstream indata;
  TString rootFileName;
  TChain *theChain = new TChain(treeName.Data());
  indata.open(fileList.Data()); 
  assert(indata.is_open());  
  while (indata  >> rootFileName ){
    if(indata.eof()){
      std::cout<<"EOF"<<std::endl;
      break;
    }
    cout<<"        adding "<<rootFileName<<endl;
    theChain->Add(rootFileName.Data(),-1);
  }
  indata.close();
  Init(theChain);
}

anabasestereo::anabasestereo(TString inFileName, Int_t keyID) : _particle_type_name("NONE"),
								fChain(0),
								LST1_r0(-70.93, -52.07, 43.0),
								LST2_r0(-35.27,  66.14, 32.0),
								LST3_r0( 75.28,  50.49, 28.7),
								LST4_r0( 30.91, -64.54, 32.0),
								_n_file_counter(-1)
{
  if(keyID == 1){
    ifstream indata;
    TChain *theChain = new TChain(treeName.Data());
    cout<<"        adding "<<inFileName<<endl;
    theChain->Add(inFileName.Data(),-1);
    Init(theChain);
  }
  else
    assert(0);
}

anabasestereo::~anabasestereo(){
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

void anabasestereo::Loop(TString histOut){
}

Int_t anabasestereo::GetEntry(Long64_t entry){
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t anabasestereo::LoadTree(Long64_t entry){
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void anabasestereo::Init(TTree *tree){
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  //fChain->SetBranchAddress("evt", &evt, &b_evt);
  //fChain->SetBranchAddress("run", &run, &b_run);
  //fChain->SetBranchAddress("pValue", &pValue, &b_pValue);
  //...
  //...
  //
  //---------------------------------------------------
  // ADD HERE :
  fChain->SetBranchAddress("event_id", &event_id, &b_event_id);
  fChain->SetBranchAddress("energy", &energy, &b_energy);
  fChain->SetBranchAddress("azimuth", &azimuth, &b_azimuth);
  fChain->SetBranchAddress("altitude", &altitude, &b_altitude);
  fChain->SetBranchAddress("h_first_int", &h_first_int, &b_h_first_int);
  fChain->SetBranchAddress("xmax", &xmax, &b_xmax);
  fChain->SetBranchAddress("hmax", &hmax, &b_hmax);
  fChain->SetBranchAddress("emax", &emax, &b_emax);
  fChain->SetBranchAddress("cmax", &cmax, &b_cmax);
  fChain->SetBranchAddress("xcore", &xcore, &b_xcore);
  fChain->SetBranchAddress("ycore", &ycore, &b_ycore);
  fChain->SetBranchAddress("ev_time_LST1", &ev_time_LST1, &b_ev_time_LST1);
  fChain->SetBranchAddress("ev_time_LST2", &ev_time_LST2, &b_ev_time_LST2);
  fChain->SetBranchAddress("ev_time_LST3", &ev_time_LST3, &b_ev_time_LST3);
  fChain->SetBranchAddress("ev_time_LST4", &ev_time_LST4, &b_ev_time_LST4);
  fChain->SetBranchAddress("nphotons_LST1", &nphotons_LST1, &b_nphotons_LST1);
  fChain->SetBranchAddress("nphotons_LST2", &nphotons_LST2, &b_nphotons_LST2);
  fChain->SetBranchAddress("nphotons_LST3", &nphotons_LST3, &b_nphotons_LST3);
  fChain->SetBranchAddress("nphotons_LST4", &nphotons_LST4, &b_nphotons_LST4);
  fChain->SetBranchAddress("n_pe_LST1", &n_pe_LST1, &b_n_pe_LST1);
  fChain->SetBranchAddress("n_pe_LST2", &n_pe_LST2, &b_n_pe_LST2);
  fChain->SetBranchAddress("n_pe_LST3", &n_pe_LST3, &b_n_pe_LST3);
  fChain->SetBranchAddress("n_pe_LST4", &n_pe_LST4, &b_n_pe_LST4);
  fChain->SetBranchAddress("n_pixels_LST1", &n_pixels_LST1, &b_n_pixels_LST1);
  fChain->SetBranchAddress("n_pixels_LST2", &n_pixels_LST2, &b_n_pixels_LST2);
  fChain->SetBranchAddress("n_pixels_LST3", &n_pixels_LST3, &b_n_pixels_LST3);
  fChain->SetBranchAddress("n_pixels_LST4", &n_pixels_LST4, &b_n_pixels_LST4);
  fChain->SetBranchAddress("pe_chID_LST1", pe_chID_LST1, &b_pe_chID_LST1);
  fChain->SetBranchAddress("pe_chID_LST2", pe_chID_LST2, &b_pe_chID_LST2);
  fChain->SetBranchAddress("pe_chID_LST3", pe_chID_LST3, &b_pe_chID_LST3);
  fChain->SetBranchAddress("pe_chID_LST4", pe_chID_LST4, &b_pe_chID_LST4);
  fChain->SetBranchAddress("pe_time_LST1", pe_time_LST1, &b_pe_time_LST1);
  fChain->SetBranchAddress("pe_time_LST2", pe_time_LST2, &b_pe_time_LST2);
  fChain->SetBranchAddress("pe_time_LST3", pe_time_LST3, &b_pe_time_LST3);
  fChain->SetBranchAddress("pe_time_LST4", pe_time_LST4, &b_pe_time_LST4);
  //---------------------------------------------------
  Notify();
}

Bool_t anabasestereo::Notify(){
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  _n_file_counter++;
  return kTRUE;
}

void anabasestereo::Show(Long64_t entry){
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t anabasestereo::Cut(Long64_t entry){
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

Double_t anabasestereo::getExpectedTimeDelayBetweenTwoLST(const TVector3 &v1, const TVector3 &v2,
							  Double_t tel_arr_pointing_azimuth, Double_t tel_arr_pointing_altitude,
							  Double_t hmax_approximate_m){
  TVector3 R;   // approximate hmax position
  TVector3 dR1;
  TVector3 dR2;
  Double_t Rmag;
  if(TMath::Sin(tel_arr_pointing_altitude) > 0)
    Rmag = hmax_approximate_m/TMath::Sin(tel_arr_pointing_altitude); //m
  else
    Rmag = hmax_approximate_m*10.0; //m
  R.SetMagThetaPhi(Rmag,TMath::Pi()/2.0 - tel_arr_pointing_altitude, tel_arr_pointing_azimuth);
  dR1 = R - v1;
  dR2 = R - v2;
  cout<<"dR1.Mag() "<<dR1.Mag()<<endl
      <<"dR2.Mag() "<<dR2.Mag()<<endl;
  return (dR1.Mag() - dR2.Mag())/(TMath::C()*1.0e-9);
}

//coincidenceID     telID
//      0        -    1
//      1        -    2
//      2        -    3
//      3        -    4
Bool_t anabasestereo::if_single_LST(Int_t &coincidenceID){
  coincidenceID = -1;
  const Int_t nLST_tel = 4;
  Int_t tel_with_signals[nLST_tel] = {0,0,0,0};
  Int_t tel_npe[nLST_tel] = {n_pe_LST1,n_pe_LST2,n_pe_LST3,n_pe_LST4};
  //
  Int_t n_sig_tel = 0;
  for(Int_t i = 0;i<nLST_tel;i++){
    if(tel_npe[i]>0){
      tel_with_signals[i] = 1;
      n_sig_tel++;
    }
  }
  if(n_sig_tel == 1){
    for(Int_t i = 0;i<nLST_tel;i++){
      if(tel_with_signals[i]>0){
	coincidenceID = i;
	return true;
      }
    }
  }
  //
  return false;
}

//coincidenceID      telID
//      0        -    1 2
//      1        -    1 3
//      2        -    1 4
//      3        -    2 3
//      4        -    2 4
//      5        -    3 4
Bool_t anabasestereo::if_double_LST(Int_t &coincidenceID){
  coincidenceID = -1;
  const Int_t nLST_tel = 4;
  Int_t tel_with_signals[nLST_tel] = {0,0,0,0};
  Int_t tel_npe[nLST_tel] = {n_pe_LST1,n_pe_LST2,n_pe_LST3,n_pe_LST4};
  //
  Int_t n_sig_tel = 0;
  for(Int_t i = 0;i<nLST_tel;i++){
    if(tel_npe[i]>0){
      tel_with_signals[i] = 1;
      n_sig_tel++;
    }
  }
  if(n_sig_tel == 2){
    // tel 1 and tel 2
    if( tel_with_signals[0] == 1 &&
        tel_with_signals[1] == 1 ){
	coincidenceID = 0;
      return true;      
    } // tel 1 and tel 3
    else if(tel_with_signals[0] == 1 &&
	    tel_with_signals[2] == 1){
	coincidenceID = 1;
      return true;      
    } // tel 1 and tel 4
    else if(tel_with_signals[0] == 1 &&
	    tel_with_signals[3] == 1){
	coincidenceID = 2;
      return true;      
    } // tel 2 and tel 3
    else if(tel_with_signals[1] == 1 &&
	    tel_with_signals[2] == 1){
	coincidenceID = 3;
      return true;      
    } // tel 2 and tel 4
    else if(tel_with_signals[1] == 1 &&
	    tel_with_signals[3] == 1){
	coincidenceID = 4;
      return true;      
    } // tel 3 and tel 4
    else if(tel_with_signals[2] == 1 &&
	    tel_with_signals[3] == 1){
	coincidenceID = 5;
      return true;      
    }
    else{
      cout<<"tel_with_signals = "
	  <<setw(10)<<tel_with_signals[0]
	  <<setw(10)<<tel_with_signals[1]
	  <<setw(10)<<tel_with_signals[2]
	  <<setw(10)<<tel_with_signals[3]<<endl;
      cout<<"tel_npe          = "
	  <<setw(10)<<tel_npe[0]
	  <<setw(10)<<tel_npe[1]
	  <<setw(10)<<tel_npe[2]
	  <<setw(10)<<tel_npe[3]<<endl;
      assert(0);
    }
  }
  //
  return false;
}

//coincidenceID       telID
//      0        -    1 2 3
//      1        -    2 3 4
//      2        -    3 4 1
//      3        -    1 2 4
Bool_t anabasestereo::if_triple_LST(Int_t &coincidenceID){
  coincidenceID = -1;
  const Int_t nLST_tel = 4;
  Int_t tel_with_signals[nLST_tel] = {0,0,0,0};
  Int_t tel_npe[nLST_tel] = {n_pe_LST1,n_pe_LST2,n_pe_LST3,n_pe_LST4};
  //
  Int_t n_sig_tel = 0;
  for(Int_t i = 0;i<nLST_tel;i++){
    if(tel_npe[i]>0){
      tel_with_signals[i] = 1;
      n_sig_tel++;
    }
  }
  if(n_sig_tel == 3){
    // tel 1 and tel 2 and tel 3
    if( tel_with_signals[0] == 1 &&
        tel_with_signals[1] == 1 &&
	tel_with_signals[2] == 1 ){
	coincidenceID = 0;
      return true;      
    }  // tel 2 and tel 3 and tel 4
    else if(tel_with_signals[1] == 1 &&
	    tel_with_signals[2] == 1 &&
	    tel_with_signals[3] == 1){
	coincidenceID = 1;
      return true;      
    } // tel 3 and tel 4 and tel 1
    else if(tel_with_signals[2] == 1 &&
	    tel_with_signals[3] == 1 &&
	    tel_with_signals[0] == 1){
	coincidenceID = 2;
      return true;      
    } // tel 1 and tel 2 and tel 4
    else if(tel_with_signals[0] == 1 &&
	    tel_with_signals[1] == 1 &&
	    tel_with_signals[3] == 1){
	coincidenceID = 3;
      return true;      
    }
    else{
      cout<<"tel_with_signals = "
	  <<setw(10)<<tel_with_signals[0]
	  <<setw(10)<<tel_with_signals[1]
	  <<setw(10)<<tel_with_signals[2]
	  <<setw(10)<<tel_with_signals[3]<<endl;
      cout<<"tel_npe          = "
	  <<setw(10)<<tel_npe[0]
	  <<setw(10)<<tel_npe[1]
	  <<setw(10)<<tel_npe[2]
	  <<setw(10)<<tel_npe[3]<<endl;
      assert(0);
    }
  }
  //
  return false;
}

//coincidenceID         telID
//      0         -    1 2 3 4
Bool_t anabasestereo::if_four_LST(Int_t &coincidenceID){
  coincidenceID = -1;
  const Int_t nLST_tel = 4;
  Int_t tel_npe[nLST_tel] = {n_pe_LST1,n_pe_LST2,n_pe_LST3,n_pe_LST4};
  //
  Int_t n_sig_tel = 0;
  for(Int_t i = 0;i<nLST_tel;i++){
    if(tel_npe[i]>0)
      n_sig_tel++;
  }
  if(n_sig_tel == 4){
    coincidenceID = 0;
    return true;
  }
  //
  return false;
}

Bool_t anabasestereo::if_stereo_trigger_pne(Int_t npecut){
  const Int_t nLST_tel = 4;
  Int_t tel_npe[nLST_tel] = {n_pe_LST1,n_pe_LST2,n_pe_LST3,n_pe_LST4};
  //
  Int_t n_sig_tel = 0;
  for(Int_t i = 0;i<nLST_tel;i++){
    if(tel_npe[i]>=npecut)
      n_sig_tel++;
  }
  if(n_sig_tel > 1)
    return true;
  //
  return false;
}
