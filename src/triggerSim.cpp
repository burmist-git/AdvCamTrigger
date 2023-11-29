//my
#include "triggerSim.hh"
#include "sipmCameraHist.hh"

//root
#include "TH1D.h"

//C, C++
#include <iostream>
#include <vector>
#include <array>
#include <assert.h>
#include <iomanip>
#include <stdlib.h>

triggerSim::triggerSim(const sipmCameraHist* simphist) : _simphist(simphist)
{
}

triggerSim::~triggerSim(){
}

//std::vector<std::array<int, 2>> triggerSim::get_trigger(const std::vector<std::vector<int>> &wf){
//
//std::vector<std::array<int, 2>> trg_vector;
//std::array<int, 2> chID_and_timeBin;
//
//chID_and_timeBin[0] = 100;
//chID_and_timeBin[1] = 35;
//trg_vector.push_back(chID_and_timeBin);  
//
//chID_and_timeBin[0] = 0;
//chID_and_timeBin[1] = 36;
//trg_vector.push_back(chID_and_timeBin);  
//
//chID_and_timeBin[0] = 10;
//chID_and_timeBin[1] = 37;
//trg_vector.push_back(chID_and_timeBin);  
//
//return trg_vector;
//}

std::vector<std::vector<int>> triggerSim::get_trigger_test(){
  std::vector<std::vector<int>> trg_vector;
  return trg_vector;
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
  
std::vector<std::vector<unsigned int>> triggerSim::get_trigger(const std::vector<std::vector<int>> &wf,
							       TH1D *h1_digital_sum,
							       TH1D *h1_digital_sum_3ns,
							       TH1D *h1_digital_sum_5ns,
							       TH1D *h1_fadc_val){
  std::vector<std::vector<unsigned int>> trg_vector;
  std::vector<unsigned int> trg_chID;
  int digital_sum;
  int digital_sum_3ns;
  int digital_sum_5ns;
  int fadc_val;
  for(unsigned int wf_j = 0;wf_j<wf.at(0).size();wf_j++){
    trg_chID.clear();
    for(unsigned int ch_i = 0;ch_i<wf.size();ch_i++){
      //digital_sum = get_flower_digital_sum(ch_i,wf_j,wf,0,0,true);
      //digital_sum_3ns = get_flower_digital_sum(ch_i,wf_j,wf,-1,1,true);
      //digital_sum_5ns = get_flower_digital_sum(ch_i,wf_j,wf,-2,2,true);
      digital_sum     = get_digital_sum( ch_i, wf_j, wf,  0, 0, true, 0);
      digital_sum_3ns = get_digital_sum( ch_i, wf_j, wf, -1, 1, true, 0);
      digital_sum_5ns = get_digital_sum( ch_i, wf_j, wf, -2, 2, true, 0);
      fadc_val = wf.at(ch_i).at(wf_j);
      if(h1_digital_sum != NULL)
	h1_digital_sum->Fill(digital_sum);
      if(h1_digital_sum_3ns != NULL)
	h1_digital_sum_3ns->Fill(digital_sum_3ns);
      if(h1_digital_sum_5ns != NULL)
	h1_digital_sum_5ns->Fill(digital_sum_5ns);
      if(h1_fadc_val != NULL)
	h1_fadc_val ->Fill(fadc_val);
      if(digital_sum_5ns>308){
	trg_chID.push_back(ch_i);
	//std::cout<<ch_i<<std::endl;
      }
      //std::cout<<digital_sum<<std::endl;
    }
    trg_vector.push_back(trg_chID);
  }
  return trg_vector;
}

int triggerSim::get_digital_sum( const unsigned int ch_i, const unsigned int wf_j, const std::vector<std::vector<int>> &wf, Int_t w_l=-1, Int_t w_r=1,
				 Bool_t norm_yes = true, Int_t sum_type = 0){
  int digital_sum = 0;
  int norm = 0;
  //std::cout<<"_simphist->get_pixel_vec().size() "<<_simphist->get_pixel_vec().size()<<std::endl;
  int j_start_int;
  int j_stop_int;
  unsigned int j_start;
  unsigned int j_stop;
  //
  j_start_int = wf_j + w_l;
  j_stop_int  = wf_j + w_r;
  j_start_int = (j_start_int >= 0 ? j_start_int : 0);
  j_stop_int  = (j_stop_int  < (Int_t)wf.at(ch_i).size() ? j_stop_int : (Int_t)(wf.at(ch_i).size()-1));
  //
  j_start = (unsigned int)j_start_int;
  j_stop = (unsigned int)j_stop_int;
  //
  for(unsigned int j = j_start; j<=j_stop; j++){
    digital_sum += wf.at(ch_i).at(j);
    norm++;
  }
  if(sum_type == 0){
    for(unsigned int i = 0 ;i<_simphist->get_pixel_vec().at(ch_i).v_pixel_flower.size();i++){
      for(unsigned int j = j_start; j<=j_stop; j++){
	digital_sum += wf.at(_simphist->get_pixel_vec().at(ch_i).v_pixel_flower.at(i).pixel_id).at(j);
	norm++;
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
  else{
    assert(0);
  }
  //if(digital_sum == 0)
  //std::cout<<ch_i<<" "<<wf_j<<" "<<j_start<<" "<<j_stop<<std::endl;
  if(norm_yes && norm>0)
    return digital_sum/norm;
  return digital_sum;
}
  
int triggerSim::get_flower_digital_sum( const unsigned int ch_i, const unsigned int wf_j, const std::vector<std::vector<int>> &wf,
					Int_t w_l=-1, Int_t w_r=1, Bool_t norm_yes = true){
  return get_digital_sum( ch_i, wf_j, wf, w_l, w_r, norm_yes, 0);
}
