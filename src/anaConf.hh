#ifndef anaConf_hh
#define anaConf_hh

//root
#include "TROOT.h"

//C, C++
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <assert.h>

struct anaConf {
  TString confFileName;
  Long64_t nentries_max;
  Long64_t jentry_modulo;
  Int_t n_ev_cuts_max;
  Bool_t wf_trg_sim;
  Bool_t disable_all_cuts;
  Bool_t cuts_set_to_false;
  Bool_t if_cout_jentry;
  //
  Double_t ellipse_A;
  Double_t ellipse_B;
  anaConf(){
    confFileName="NONE";
    nentries_max = -1;
    jentry_modulo = 100000;
    n_ev_cuts_max = -1;
    wf_trg_sim = false;
    disable_all_cuts = true;
    cuts_set_to_false = false;
    if_cout_jentry = false;
    ellipse_A = 851.3;
    ellipse_B = 800.0;
  }
  void printInfo(){
    std::cout<<"confFileName      "<<confFileName<<std::endl
	     <<"nentries_max      "<<nentries_max<<std::endl
      	     <<"jentry_modulo     "<<jentry_modulo<<std::endl
	     <<"n_ev_cuts_max     "<<n_ev_cuts_max<<std::endl
	     <<"wf_trg_sim        "<<wf_trg_sim<<std::endl
    	     <<"disable_all_cuts  "<<disable_all_cuts<<std::endl
	     <<"cuts_set_to_false "<<cuts_set_to_false<<std::endl
      	     <<"if_cout_jentry    "<<if_cout_jentry<<std::endl
    	     <<"ellipse_A         "<<ellipse_A<<std::endl
	     <<"ellipse_B         "<<ellipse_B<<std::endl;
  }  
  void readFromFile(TString name){
    std::ifstream confFile(name.Data());
    confFileName = name;
    if (confFile.is_open()) {
      std::string mot;
      while(confFile>>mot){
	if(mot == "nentries_max:")
	  confFile>>nentries_max;
	if(mot == "jentry_modulo:")
	  confFile>>jentry_modulo;
	if(mot == "n_ev_cuts_max:")
	  confFile>>n_ev_cuts_max;
	if(mot == "wf_trg_sim:")
	  confFile>>wf_trg_sim;
	if(mot == "disable_all_cuts:")
	  confFile>>disable_all_cuts;
	if(mot == "cuts_set_to_false:")
	  confFile>>cuts_set_to_false;
	if(mot == "if_cout_jentry:")
	  confFile>>if_cout_jentry;
	if(mot == "ellipse_A:")
	  confFile>>ellipse_A;
	if(mot == "ellipse_B:")
	  confFile>>ellipse_B;
      }
      confFile.close();
    }
    else {
      std::cout << "Unable to open file"<<std::endl; 
      assert(0);
    }
  }
};

#endif
