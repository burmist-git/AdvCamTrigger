//my
#include "triggerSim.hh"
#include "sipmCameraHist.hh"
#include "dbscan.hh"

//root
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TFile.h"

//C, C++
#include <iostream>
#include <vector>
#include <array>
#include <assert.h>
#include <iomanip>
#include <stdlib.h>

triggerSim::triggerSim(const sipmCameraHist* simphist) : _simphist(simphist), _dbs(new dbscan()), _trg_counter(0), _n_skip_edge_points(0), _k_dist_graph_flag(false), _digital_sum_max_only(false)
{
  //_dbs->print_cluster_stats();
  autofill_trg_channel_mask();
}

triggerSim::~triggerSim(){
}

void triggerSim::autofill_trg_channel_mask( unsigned int nch_max, unsigned int ival){
  if(_trg_channel_mask.size()>0)
    _trg_channel_mask.clear();
  for( unsigned int ii = 0; ii < nch_max; ii++)
    _trg_channel_mask.push_back(ival);
}

void triggerSim::fill_trg_channel_mask_from_file(TString mask_file){
  int val_chID;
  if(_trg_channel_mask.size()>0)
    _trg_channel_mask.clear();
  autofill_trg_channel_mask(8000,0);
  //
  ifstream fFile(mask_file.Data());
  //
  cout<<"mask_file = "<<mask_file<<endl;
  //
  if(fFile.is_open()){
    while(fFile>>val_chID){
      _trg_channel_mask.at((unsigned int)val_chID) = 1;
    }
    fFile.close();
  }
  else{
    std::cout<<" ERROR --> file : "<<mask_file<<std::endl
             <<" does not exist."<<std::endl;
    assert(0);
  }
}

void triggerSim::print_trg_channel_mask(int verbosity){
  std::cout<<"Trigger channel mask file : "<<_trg_setup.trigger_channel_mask_file_list<<std::endl
	   <<"_trg_channel_mask.size()  : "<<_trg_channel_mask.size()<<std::endl;
  if(_trg_setup.trigger_channel_mask_file_list.Sizeof()<=1)
    std::cout<<"The trigger channel mask file is not specified - the autofill function is executed."<<std::endl;
  if(verbosity>1){
    std::cout<<"channel list "<<std::endl;
    for( unsigned int ii = 0; ii < _trg_channel_mask.size(); ii++)
      std::cout<<_trg_channel_mask.at(ii)<<std::endl;
  }
}

int triggerSim::get_flower_digital_sum( const unsigned int ch_i, const unsigned int wf_j, const std::vector<std::vector<int>> &wf, Bool_t norm_yes = true){
  return get_digital_sum( ch_i, wf_j, wf, norm_yes, 0);
}

std::vector<std::vector<int>> triggerSim::get_trigger_test(){
  std::vector<std::vector<int>> trg_vector;
  return trg_vector;
}

void triggerSim::print_trigger_vec(const std::vector<std::array<int, 2>> &trg_vector){
  for(unsigned int i = 0;i<trg_vector.size();i++)
    std::cout<<trg_vector.at(i).at(0)<<" "<<trg_vector.at(i).at(1)<<std::endl;
}

void triggerSim::print_trigger_vec(const std::vector<std::vector<unsigned int>> &trg_vector){
  std::cout<<"print_trigger_vec"<<std::endl;
  for(unsigned int i = 0;i<trg_vector.size();i++){
    std::cout<<i;
    for(unsigned int j = 0;j<trg_vector.at(i).size();j++)
      std::cout<<" "<<trg_vector.at(i).at(j);
    std::cout<<std::endl;
  }
}

void triggerSim::print_trigger_vec_to_csv(const std::vector<std::vector<unsigned int>> &trg_vector, const sipmCameraHist *sipm_cam, TString out_file_name, bool if_short_format){
  std::ofstream outfile;
  outfile.open(out_file_name.Data());
  if(!if_short_format){
    for(unsigned int i = 0;i<trg_vector.size();i++){
      outfile<<std::setw(4)<<i;
      outfile<<std::setw(10)<<trg_vector.at(i).size()<<std::endl;
      if(trg_vector.at(i).size()>0){
	for(unsigned int j = 0;j<trg_vector.at(i).size();j++)
	  outfile<<" "<<sipm_cam->get_pixel_vec().at(trg_vector.at(i).at(j)).x;
	outfile<<std::endl;
	for(unsigned int j = 0;j<trg_vector.at(i).size();j++)
	  outfile<<" "<<sipm_cam->get_pixel_vec().at(trg_vector.at(i).at(j)).y;    
	outfile<<std::endl;
      }
    }
  }
  else{
    for(unsigned int i = 0;i<trg_vector.size();i++){
      if(trg_vector.at(i).size()>0){
	for(unsigned int j = 0;j<trg_vector.at(i).size();j++){
	  outfile<<sipm_cam->get_pixel_vec().at(trg_vector.at(i).at(j)).x
		 <<" "<<sipm_cam->get_pixel_vec().at(trg_vector.at(i).at(j)).y
		 <<" "<<i
		 <<" "<<j
		 <<" "<<trg_vector.at(i).at(j)<<std::endl;
	}
      }
    }
  }
  outfile.close();
}

void const triggerSim::fill_fadc_val_vs_time( const std::vector<std::vector<int>> &wf, TH2D *h2) const {
  for(unsigned int ch_i = 0;ch_i<wf.size();ch_i++)
    for(unsigned int wf_j = 0;wf_j<wf.at(ch_i).size();wf_j++)
      h2->Fill(wf_j,wf.at(ch_i).at(wf_j));
}

std::vector<std::vector<unsigned int>> triggerSim::get_trigger(const std::vector<std::vector<int>> &wf,
							       TH1D *h1_digital_sum,
							       TH1D *h1_fadc_val,
							       TH1D *h1_digital_sum_pe,
							       TH1D *h1_fadc_val_pe,
							       Double_t adc_per_pe,
							       Double_t pedestal_digital_sum_pe,
							       Double_t pedestal_fadc_pe){
  //
  std::vector<std::vector<unsigned int>> trg_vector;
  std::vector<std::vector<unsigned int>> trg_vector_empty;
  std::vector<unsigned int> trg_chID;
  trg_vector.clear();
  trg_chID.clear();
  trg_vector_empty.clear();
  int digital_sum = 0;
  int fadc_val;
  //
  for(unsigned int wf_j = (0 + _n_skip_edge_points);wf_j<(wf.at(0).size() - _n_skip_edge_points);wf_j++){
    trg_chID.clear();
    for(unsigned int ch_i = 0;ch_i<wf.size();ch_i++){
      if(_trg_channel_mask.at(ch_i) == 1){
	//
	digital_sum = get_digital_sum( ch_i, wf_j, wf, _trg_setup.norm_yes, _trg_setup.sum_type);
	fadc_val = wf.at(ch_i).at(wf_j);
	if(h1_digital_sum != NULL)
	  h1_digital_sum->Fill(digital_sum);
	if(h1_fadc_val != NULL)
	  h1_fadc_val->Fill(fadc_val);
	//
	//if(h1_digital_sum_pe != NULL)
	//h1_digital_sum_pe->Fill(digital_sum/adc_per_pe);
	//if(h1_fadc_val_pe != NULL)
	//h1_fadc_val_pe->Fill(fadc_val/adc_per_pe);
	//
	if(h1_digital_sum_pe != NULL)
	  h1_digital_sum_pe->Fill((digital_sum/adc_per_pe - pedestal_digital_sum_pe));
	if(h1_fadc_val_pe != NULL)
	  h1_fadc_val_pe->Fill((fadc_val/adc_per_pe - pedestal_fadc_pe));
	//
	if(digital_sum>_trg_setup.digital_sum_cut)
	  trg_chID.push_back(ch_i);
      }
    }
    trg_vector.push_back(trg_chID);
  }
  //
  //print_trigger_vec(trg_vector);
  //return build_spatial_cluster(trg_vector);
  //return build_spatial_time_cluster(trg_vector);
  //std::cout<<"trg_vector.size() = "<<trg_vector.size()<<std::endl;
  _n_digital_sum_micro_clusters = 0;
  for(unsigned int jj = 0; jj<trg_vector.size();jj++)
    _n_digital_sum_micro_clusters += (Int_t)trg_vector.at(jj).size();
  if(_trg_setup.if_dbscan)
    return build_spatial_time_cluster_dbscan(trg_vector);
  else
    return trg_vector;
  //
  return trg_vector_empty;
}

std::vector<std::vector<unsigned int>> triggerSim::build_spatial_time_cluster_dbscan(const std::vector<std::vector<unsigned int>> &trg_vector){
  vector<point> points;
  points.clear();
  //
  Double_t i_time_todist = _trg_setup.dbscan_time_todist;
  unsigned int minPts = _trg_setup.dbscan_minPts;
  float eps = _trg_setup.dbscan_eps;
  //
  for(unsigned int i = 0;i<trg_vector.size();i++){
    for(unsigned int j = 0;j<trg_vector.at(i).size();j++){
      point singlepoint;
      singlepoint.x = _simphist->get_pixel_vec().at(trg_vector.at(i).at(j)).x;
      singlepoint.y = _simphist->get_pixel_vec().at(trg_vector.at(i).at(j)).y;
      singlepoint.z = i*i_time_todist;
      singlepoint.time_ii = i;
      singlepoint.pixel_id = trg_vector.at(i).at(j);
      singlepoint.point_id = (Int_t)points.size();
      points.push_back(singlepoint);      
    }
  }
  _dbs->run( minPts, eps, points);
  //print_trigger_vec(trg_vector);
  //_dbs->print_points_info();
  _dbs->get_cluster_stats();
  _dbs->print_cluster_stats();
  //
  if(_k_dist_graph_flag){
    TString outrootFile;
    outrootFile = "hist_k_dist_graph_";
    outrootFile += _trg_counter;
    outrootFile += "ev";
    outrootFile += ".root";
    vector<Double_t> k_dist_graph;
    k_dist_graph = _dbs->build_k_dist_graph(3);
    plot_and_save_to_hist_root(outrootFile,k_dist_graph);
    _trg_counter++;
  }
  //
  std::vector<std::vector<unsigned int>> cam_trg_vector;
  std::vector<unsigned int> spatial_cluster;
  cam_trg_vector.clear();
  spatial_cluster.clear();
  //
  for(unsigned int i = 0;i<trg_vector.size();i++){
    spatial_cluster.clear();
    for(unsigned int j = 0;j<trg_vector.at(i).size();j++){
      for(unsigned int k = 0;k<_dbs->get_points_v().size();k++){
	if(_dbs->get_points_v().at(k).pixel_id == trg_vector.at(i).at(j) && _dbs->get_points_v().at(k).time_ii == i){
	  if(_dbs->get_points_v().at(k).clusterID > -1)
	    spatial_cluster.push_back(trg_vector.at(i).at(j));
	}
      }
    }
    cam_trg_vector.push_back(spatial_cluster);
  }
  _dbclusters_v.clear();
  _dbclusters_v = _dbs->get_clusters_v();
  _dbscan_run_time_musec = _dbs->get_dbscan_run_time_musec();
  _dbscan_N_points = _dbs->get_points_v().size();
  //
  _dbs->clear();
  //
  //std::cout<<"_dbscan_N_points "<<_dbscan_N_points<<std::endl
  //	   <<"points.size()    "<<points.size()<<std::endl;
  //
  //dbscan::print_cluster_stats(_dbclusters_v);
  //print_trigger_vec(cam_trg_vector);
  //
  return cam_trg_vector;
}

void triggerSim::plot_and_save_to_hist_root(TString outrootFile, vector<Double_t> &k_dist_graph){
  if(k_dist_graph.size()>0){
    TGraph *gr_k_dist_graph = new TGraph();
    gr_k_dist_graph->SetNameTitle("gr_k_dist_graph","gr_k_dist_graph");
    //
    for(unsigned int k = 0; k < k_dist_graph.size(); k++)
      gr_k_dist_graph->SetPoint( k, k, k_dist_graph.at(k_dist_graph.size()-1-k));
    //  
    TFile* rootFile = new TFile(outrootFile.Data(), "RECREATE", " Histograms", 1);
    rootFile->cd();
    //
    gr_k_dist_graph->Write();
    //
    rootFile->Close();
  }
}

std::vector<std::vector<unsigned int>> triggerSim::build_spatial_cluster(const std::vector<std::vector<unsigned int>> &trg_vector){
  std::vector<std::vector<unsigned int>> cam_trg_vector;
  std::vector<unsigned int> spatial_cluster;
  Double_t x1;
  Double_t y1;
  Double_t x2;
  Double_t y2;
  Double_t dist;
  Int_t n = 0;
  for(unsigned int i = 0;i<trg_vector.size();i++){
    spatial_cluster.clear();
    for(unsigned int j = 0;j<trg_vector.at(i).size();j++){
      //std::cout<<_simphist->get_pixel_vec().at(trg_vector.at(i).at(j)).pixel_id<<" "<<trg_vector.at(i).at(j)<<" ";
      n=0;
      x1=_simphist->get_pixel_vec().at(trg_vector.at(i).at(j)).x;
      y1=_simphist->get_pixel_vec().at(trg_vector.at(i).at(j)).y;
      for(unsigned int k = 0;k<trg_vector.at(i).size();k++){
	x2=_simphist->get_pixel_vec().at(trg_vector.at(i).at(k)).x;
	y2=_simphist->get_pixel_vec().at(trg_vector.at(i).at(k)).y;
	dist = TMath::Sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
	if(dist<0.1)
	  n++;
      }
      if(n>=2)
	spatial_cluster.push_back(trg_vector.at(i).at(j));
    }
    cam_trg_vector.push_back(spatial_cluster);
  }
  return cam_trg_vector;
}

//
// sum_type = 0 pixel flower                         (first  neighbors)
// sum_type = 1 pixel flower + neighbors             (second neighbors)
// sum_type = 2 pixel flower + neighbors + neighbors (third  neighbors) 
// sum_type = 3 flower of flowers                    (flower of flowers)
// sum_type = 4 one pixel                            (one pixel)
// norm_yes = true | false                           (normalise by number of operations)
//
int triggerSim::get_digital_sum( const unsigned int ch_i, const unsigned int wf_j,
				 const std::vector<std::vector<int>> &wf,
				 Bool_t norm_yes = true, Int_t sum_type = 0){
  int digital_sum = 0;
  int norm = 0;
  int j_start_int;
  int j_stop_int;
  int w_l = 0;
  int w_r = _trg_setup.nFilterPerChannel;
  unsigned int j_start;
  unsigned int j_stop;
  //
  int filter_val_iter = 0;
  //
  j_start_int = wf_j + w_l;
  j_stop_int  = wf_j + w_r;
  j_start_int = (j_start_int >= 0 ? j_start_int : 0);
  j_stop_int  = (j_stop_int  < (Int_t)wf.at(ch_i).size() ? j_stop_int : (Int_t)(wf.at(ch_i).size()-1));
  //
  j_start = (unsigned int)j_start_int;
  j_stop = (unsigned int)j_stop_int;
  //
  //cout<<"j_start "<<j_start<<endl
  //    <<"j_stop  "<<j_stop<<endl;
  //
  filter_val_iter = 0;
  for(unsigned int j = j_start; j<=j_stop; j++){
    if(_trg_setup.nFilterPerChannel>0){
      if(filter_val_iter<_trg_setup.nFilterPerChannel)
	digital_sum += wf.at(ch_i).at(j)*_trg_setup.filterPerChannel_val[filter_val_iter];
      else
	assert(0);
      filter_val_iter++;
    }
    else{
      digital_sum += wf.at(ch_i).at(j);
    }
    norm++;
  }
  if(sum_type == 0){
    for(unsigned int i = 0 ;i<_simphist->get_pixel_vec().at(ch_i).v_pixel_flower.size();i++){
      if(_trg_setup.nFilterPerChannel>0){
	filter_val_iter = 0;
	for(unsigned int j = j_start; j<=j_stop; j++){
	  if(filter_val_iter<_trg_setup.nFilterPerChannel)
	    digital_sum += wf.at(_simphist->get_pixel_vec().at(ch_i).v_pixel_flower.at(i).pixel_id).at(j)*_trg_setup.filterPerChannel_val[filter_val_iter];
	  else
	    assert(0);	
	  filter_val_iter++;
	  norm++;
	}
      }
      else{
	for(unsigned int j = j_start; j<=j_stop; j++){
	  digital_sum += wf.at(_simphist->get_pixel_vec().at(ch_i).v_pixel_flower.at(i).pixel_id).at(j);
	  norm++;
	}
      }
    }
  }
  else if(sum_type == 1){
    for(unsigned int i = 0 ;i<_simphist->get_pixel_vec().at(ch_i).v_pixel_neighbors_second.size();i++){
      for(unsigned int j = j_start; j<=j_stop; j++){
	digital_sum += wf.at(_simphist->get_pixel_vec().at(ch_i).v_pixel_neighbors_second.at(i).pixel_id).at(j);
	norm++;
      }
    }
  }
  else if(sum_type == 2){
    for(unsigned int i = 0 ;i<_simphist->get_pixel_vec().at(ch_i).v_pixel_neighbors_third.size();i++){
      for(unsigned int j = j_start; j<=j_stop; j++){
	digital_sum += wf.at(_simphist->get_pixel_vec().at(ch_i).v_pixel_neighbors_third.at(i).pixel_id).at(j);
	norm++;
      }
    }
  }
  else if(sum_type == 3){
    for(unsigned int i = 0 ;i<_simphist->get_pixel_vec().at(ch_i).v_pixel_super_flower.size();i++){
      for(unsigned int j = j_start; j<=j_stop; j++){
	digital_sum += wf.at(_simphist->get_pixel_vec().at(ch_i).v_pixel_super_flower.at(i).pixel_id).at(j);
	norm++;
      }
    }
  }
  else if(sum_type == 4){
    norm++;
  }
  else{
    assert(0);
  }
  if(norm_yes && norm>0)
    return digital_sum/norm;
  return digital_sum;
}

std::vector<std::vector<int>> triggerSim::get_trigger_test(const std::vector<std::vector<int>> &wf){
  //
  std::vector<std::vector<int>> trg_vector;
  std::vector<int> trg_chID;
  std::vector<int> trg_chID_empty;
  //
  trg_chID.clear();
  trg_chID_empty.clear();
  for(unsigned int i = 0;i<wf.at(0).size();i++){
    if(i == 36){
      trg_chID.clear();
      trg_chID.push_back(100);
      trg_chID.push_back(10);
      trg_chID.push_back(0);
      trg_vector.push_back(trg_chID);
    }
    else if(i == 37){
      trg_chID.clear();
      trg_chID.push_back(100);
      trg_chID.push_back(10);
      trg_chID.push_back(0);
      trg_vector.push_back(trg_chID);
    }
    else if(i == 38){
      trg_chID.clear();
      trg_chID.push_back(101);
      trg_chID.push_back(11);
      trg_chID.push_back(1);
      trg_vector.push_back(trg_chID);
    }
    else if(i == 39){
      trg_chID.clear();
      trg_chID.push_back(101);
      trg_chID.push_back(11);
      trg_chID.push_back(1);
      trg_vector.push_back(trg_chID);
    }
    else{
      trg_vector.push_back(trg_chID_empty);
    }
  }
  //
  return trg_vector;
}
