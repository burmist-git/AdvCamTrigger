#ifndef triggerSim_hh
#define triggerSim_hh

//my
#include "dbscan.hh"

//C, C++
#include <vector>
#include <array>
#include <iostream>
#include <fstream>

//root
#include <TROOT.h>

class sipmCameraHist;
class dbscan;
class TH1D;
class TH2D;
class TString;

struct trg_setup {
  // digital filter per channel
  static const Int_t nFilterPerChannelMax = 10; 
  Int_t nFilterPerChannel;
  Double_t filterPerChannel_val[nFilterPerChannelMax];
  //
  // digital sum
  Int_t sum_type;
  Bool_t norm_yes;
  Int_t digital_sum_cut;
  //
  // digital filter per digital sum
  static const Int_t nFilterPerDigitalSumMax = 10; 
  Int_t nFilterPerDigitalSum;
  Double_t filterPerDigitalSum_val[nFilterPerDigitalSumMax];
  // dbscan
  Bool_t if_dbscan;
  Double_t dbscan_time_todist = 0.05;
  unsigned int dbscan_minPts = 22;
  float dbscan_eps = 0.1;
  //
  TString trigger_channel_mask_file_list;
  //
  trg_setup(){
    //
    nFilterPerChannel = 0;
    filterPerChannel_val[0] = 1.0;
    for (Int_t ii = 1;ii<nFilterPerChannelMax;ii++)
      filterPerChannel_val[ii] = 0.0;    
    //
    // digital sum
    sum_type = 0;
    norm_yes = true;
    digital_sum_cut = 307;
    //
    //Digital filter per digital sum
    nFilterPerDigitalSum = 0;
    filterPerDigitalSum_val[0] = 1.0;
    for (Int_t ii = 1;ii<nFilterPerDigitalSumMax;ii++)
      filterPerDigitalSum_val[ii] = 0.0;    
    //
    // dbscan
    if_dbscan = 1;
    dbscan_time_todist = 0.05;
    dbscan_minPts = 22;
    dbscan_eps = 0.1;
    //
    trigger_channel_mask_file_list="";
  }
  //
  void trg_setup_info(){
    std::cout<<"digital filter per channel"<<std::endl;
    std::cout<<std::setw(25)<<"nFilterPerChannel: "<<nFilterPerChannel<<std::endl;
    for (Int_t ii = 0;ii<nFilterPerChannel;ii++)
      std::cout<<" "<<filterPerChannel_val[ii];    
    std::cout<<std::endl;    
    std::cout<<"digital sum"<<std::endl;
    std::cout<<std::setw(10)<<"sum_type"
	     <<std::setw(10)<<"norm_yes"
	     <<std::setw(20)<<"digital_sum_cut"
	     <<std::endl;
    std::cout<<std::setw(10)<<sum_type
	     <<std::setw(10)<<norm_yes
	     <<std::setw(20)<<digital_sum_cut
	     <<std::endl;
    std::cout<<"digital filter per digital sum"<<std::endl;
    std::cout<<std::setw(25)<<"nFilterPerDigitalSum: "<<nFilterPerDigitalSum<<std::endl;
    for (Int_t ii = 0;ii<nFilterPerDigitalSum;ii++)
      std::cout<<" "<<filterPerDigitalSum_val[ii];        
    std::cout<<"dbscan"<<std::endl;
    std::cout<<std::setw(20)<<"if_dbscan"
	     <<std::setw(25)<<"dbscan_time_todist"
	     <<std::setw(20)<<"dbscan_minPts"
	     <<std::setw(20)<<"dbscan_eps"
	     <<std::endl;
    std::cout<<std::setw(20)<<if_dbscan
	     <<std::setw(25)<<dbscan_time_todist
	     <<std::setw(20)<<dbscan_minPts
	     <<std::setw(20)<<dbscan_eps
	     <<std::endl;
    std::cout<<"trigger_channel_mask_file_list"<<trigger_channel_mask_file_list<<std::endl;
  }
  //
  void load_trg_setup(const char* in_file){
    //
    ifstream fFile(in_file);
    //std::cout<<" Trigger setup file : "<<in_file<<std::endl;
    //
    Float_t x, y, drawer_id;
    Int_t pixel_id = 0;
    string mot;
    Double_t val;
    Double_t val_tmp;
    Double_t val_int;
    //
    if(fFile.is_open()){
      while(fFile>>mot){
	//cout<<"motB "<<mot<<endl;
	if(mot == "nFilterPerChannel:"){
	  fFile>>nFilterPerChannel;
	  cout<<"---> nFilterPerChannel "<<nFilterPerChannel<<endl;
	}
	else if(mot == "filterPerChannel_val:"){
	  for(Int_t ii = 0;ii<nFilterPerChannelMax;ii++){
	    fFile>>val;
	    filterPerChannel_val[ii] = val;
	  }
	  cout<<"---> filterPerChannel_val ";
	  for(Int_t ii = 0;ii<nFilterPerChannelMax;ii++)
	    cout<<filterPerChannel_val[ii]<<" ";
	  cout<<endl;
	}	
	else if(mot == "sum_type:"){
	  fFile>>sum_type;	  
	  cout<<"---> sum_type "<<sum_type<<endl;
	}
	else if(mot == "norm_yes:"){
	  fFile>>val_int;
	  if(val_int == 1)
	    norm_yes = true;
	  else
	    norm_yes = false;
	  cout<<"---> norm_yes "<<norm_yes<<endl;
	}
	else if(mot == "digital_sum_cut:"){
	  fFile>>digital_sum_cut;
	  cout<<"---> digital_sum_cut "<<digital_sum_cut<<endl;
	}
	else if(mot == "nFilterPerDigitalSum:"){
	  fFile>>nFilterPerDigitalSum;
	  cout<<"---> nFilterPerDigitalSum "<<nFilterPerDigitalSum<<endl;
	}
	else if(mot == "filterPerDigitalSum_val:"){
	  for(Int_t ii = 0;ii<nFilterPerDigitalSumMax;ii++){	   	    
	    fFile>>val;
	    filterPerDigitalSum_val[ii] = val;
	  }
	  cout<<"---> filterPerDigitalSum_val ";
	  for(Int_t ii = 0;ii<nFilterPerDigitalSumMax;ii++)
	    cout<<filterPerDigitalSum_val[ii]<<" ";
	  cout<<endl;
	}
	else if(mot == "if_dbscan:"){
	  fFile>>val_int;
	  if(val_int == 1)
	    if_dbscan = true;
	  else
	    if_dbscan = false;
	  cout<<"---> if_dbscan "<<if_dbscan<<endl;
	}
	else if(mot == "dbscan_time_todist:"){
	  fFile>>dbscan_time_todist;
	  cout<<"---> dbscan_time_todist "<<dbscan_time_todist<<endl;
	}
	else if(mot == "dbscan_minPts:"){
	  fFile>>val_int;
	  dbscan_minPts = (unsigned int)val_int;
	  cout<<"---> dbscan_minPts "<<dbscan_minPts<<endl;
	}
	else if(mot == "dbscan_eps:"){
	  fFile>>val_tmp;
	  dbscan_eps = (float)val_tmp;
	  cout<<"---> dbscan_eps "<<dbscan_eps<<endl;
	}
	else if(mot == "trigger_channel_mask_file_list:"){
	  fFile>>trigger_channel_mask_file_list;
	  cout<<"---> trigger_channel_mask_file_list "<<trigger_channel_mask_file_list<<endl;
	}	
      }
      fFile.close();
    }
  }
};
  
class triggerSim {
    
  public :
  //
  triggerSim(const sipmCameraHist* simphist);
  ~triggerSim();

  //std::vector<std::array<int, 2>> get_trigger(const std::vector<std::vector<int>> &wf);
  std::vector<std::vector<int>> get_trigger_test(const std::vector<std::vector<int>> &wf);
  std::vector<std::vector<int>> get_trigger_test();
  std::vector<std::vector<unsigned int>> get_trigger(const std::vector<std::vector<int>> &wf,
						     TH1D *h1_digital_sum = NULL,
						     TH1D *h1_fadc_val = NULL);
  static void print_trigger_vec(const std::vector<std::array<int, 2>> &trg_vector);
  static void print_trigger_vec(const std::vector<std::vector<unsigned int>> &trg_vector);
  static void print_trigger_vec_to_csv(const std::vector<std::vector<unsigned int>> &trg_vector, const sipmCameraHist *sipm_cam, TString out_file_name, bool if_short_format = true);

  void plot_and_save_to_hist_root(TString outrootFile, vector<Double_t> &k_dist_graph);

  inline const vector<cluster_info>& get_dbclusters() const {return _dbclusters_v;};
  
  void const fill_fadc_val_vs_time( const std::vector<std::vector<int>> &wf, TH2D *h2) const;
  
  unsigned int const get_n_skip_edge_points(unsigned int val) const {return _n_skip_edge_points;};
  inline void set_n_skip_edge_points(unsigned int val){_n_skip_edge_points = val;};

  inline const Int_t get_dbscan_run_time_musec() const {return _dbscan_run_time_musec;};
  inline const Int_t get_dbscan_N_points() const {return _dbscan_N_points;};
  inline const Int_t get_n_digital_sum_micro_clusters() const {return _n_digital_sum_micro_clusters;};
  
  inline void set_k_dist_graph_flag(bool val){_k_dist_graph_flag = val;};

  inline const Int_t get_digital_sum_max() const {return _digital_sum_max;};

  void autofill_trg_channel_mask(unsigned int nch_max = 8000, unsigned int ival = 1);
  void fill_trg_channel_mask_from_file(TString mask_file);
  void print_trg_channel_mask(int verbosity = 0);
  
  trg_setup _trg_setup;
  
private:

  int get_flower_digital_sum(const unsigned int ch_i, const unsigned int wf_j, const std::vector<std::vector<int>> &wf, Bool_t norm_yes);
  int get_digital_sum(const unsigned int ch_i, const unsigned int wf_j, const std::vector<std::vector<int>> &wf, Bool_t norm_yes, Int_t sum_type);

  std::vector<std::vector<unsigned int>> build_spatial_cluster(const std::vector<std::vector<unsigned int>> &trg_vector);
  std::vector<std::vector<unsigned int>> build_spatial_time_cluster_dbscan(const std::vector<std::vector<unsigned int>> &trg_vector);
  
  const sipmCameraHist* _simphist;

  dbscan *_dbs;
  vector<cluster_info> _dbclusters_v;
  
  Int_t _trg_counter;

  Int_t _dbscan_run_time_musec;
  Int_t _dbscan_N_points;
  Int_t _n_digital_sum_micro_clusters;
  
  unsigned int _n_skip_edge_points;

  bool _k_dist_graph_flag;

  bool _digital_sum_max_only;

  Int_t _digital_sum_max;

  vector<unsigned int> _trg_channel_mask;

};

#endif
