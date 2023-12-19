//my
#include "triggerSim.hh"
#include "sipmCameraHist.hh"
#include "dbscan.hh"

//root
#include "TH1D.h"
#include "TGraph2D.h"

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

void triggerSim::print_trigger_vec_to_csv(const std::vector<std::vector<unsigned int>> &trg_vector, const sipmCameraHist *sipm_cam, TString out_file_name, bool if_short_format){
  std::ofstream outfile;
  outfile.open(out_file_name.Data());
  //TGraph2D *gr2D = new TGraph2D();
  //TString out_gr2D_file_name = "gr2D_";
  //out_gr2D_file_name += out_file_name;
  //gr2D->SetNameTitle(out_gr2D_file_name.Data(),out_gr2D_file_name.Data());
  //TString out_h1_file_name = "h1_";
  //out_h1_file_name += out_file_name;
  //TH1D *h1_time_cluster = new TH1D(out_h1_file_name.Data(),out_h1_file_name.Data(),75,0,75);
  if(!if_short_format){
    for(unsigned int i = 0;i<trg_vector.size();i++){
      outfile<<std::setw(4)<<i;
      outfile<<std::setw(10)<<trg_vector.at(i).size()<<std::endl;
      //h1_time_cluster->SetBinContent(i+1,trg_vector.at(i).size());
      if(trg_vector.at(i).size()>0){
	//for(unsigned int j = 0;j<trg_vector.at(i).size();j++)
	//gr2D->SetPoint(gr2D->GetN(),sipm_cam->get_pixel_vec().at(trg_vector.at(i).at(j)).x,sipm_cam->get_pixel_vec().at(trg_vector.at(i).at(j)).y,i);
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
  //out_gr2D_file_name += ".C";
  //gr2D->SetDrawOption("P");
  //gr2D->SaveAs(out_gr2D_file_name.Data());
  //out_h1_file_name += ".C";
  //h1_time_cluster->SaveAs(out_h1_file_name.Data());
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
      //if(digital_sum_5ns>308){
      //trg_chID.push_back(ch_i);
      //std::cout<<ch_i<<std::endl;
      //}
      if(digital_sum_3ns>307){
	trg_chID.push_back(ch_i);
	//std::cout<<ch_i<<std::endl;
      }
      //std::cout<<digital_sum<<std::endl;
    }
    trg_vector.push_back(trg_chID);
  }
  //return trg_vector;
  //print_trigger_vec(trg_vector);
  //return build_spatial_time_cluster(trg_vector);
  return build_spatial_time_cluster_dbscan(trg_vector);
}

std::vector<std::vector<unsigned int>> triggerSim::build_spatial_time_cluster_dbscan(const std::vector<std::vector<unsigned int>> &trg_vector){
  vector<Point> points;
  //
  Double_t i_time_todist=0.05;
  //
  for(unsigned int i = 0;i<trg_vector.size();i++){
    for(unsigned int j = 0;j<trg_vector.at(i).size();j++){
      Point singlepoint;
      singlepoint.x = _simphist->get_pixel_vec().at(trg_vector.at(i).at(j)).x;
      singlepoint.y = _simphist->get_pixel_vec().at(trg_vector.at(i).at(j)).y;
      singlepoint.z = i*i_time_todist;
      singlepoint.clusterID = UNCLASSIFIED;
      singlepoint.pixel_id = trg_vector.at(i).at(j);
      points.push_back(singlepoint);      
    }
  }
  unsigned int minPts = 15;
  float eps = 0.1;
  DBSCAN ds(minPts, eps, points);
  ds.run();
  //
  ds.print_points_info();
  //
  std::vector<std::vector<unsigned int>> cam_trg_vector;
  std::vector<unsigned int> spatial_cluster;
  //
  for(unsigned int i = 0;i<trg_vector.size();i++){
    spatial_cluster.clear();
    for(unsigned int j = 0;j<trg_vector.at(i).size();j++){
      for(unsigned int k = 0;k<ds.m_points.size();k++){
	if(ds.m_points.at(k).pixel_id == trg_vector.at(i).at(j)){
	  //std::cout<<"ds.m_points.at(k).pixel_id = "<<ds.m_points.at(k).pixel_id<<std::endl;
	  if(ds.m_points.at(k).clusterID > 0)
	    spatial_cluster.push_back(trg_vector.at(i).at(j));
	}
      }
    }
    cam_trg_vector.push_back(spatial_cluster);
  }
  //
  //print_trigger_vec(cam_trg_vector);
  //
  return cam_trg_vector;
  //return trg_vector;
}

std::vector<std::vector<unsigned int>> triggerSim::build_spatial_time_cluster(const std::vector<std::vector<unsigned int>> &trg_vector){
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
	//std::cout<<"dist "<<dist<<std::endl;
	if(dist<0.1)
	  n++;
      }
      //if(n>=3 && i>=36 && i<=39)
      if(n>=3)
	spatial_cluster.push_back(trg_vector.at(i).at(j));
    }
    //std::cout<<std::endl;
    cam_trg_vector.push_back(spatial_cluster);
  }
  return cam_trg_vector;
  //return build_time_cluster(cam_trg_vector);
}

std::vector<std::vector<unsigned int>> triggerSim::build_time_cluster(const std::vector<std::vector<unsigned int>> &trg_vector){
  std::vector<std::vector<unsigned int>> cam_trg_vector;
  std::vector<unsigned int> spatial_cluster;
  std::vector<unsigned int> tmp_time_l;
  Double_t x1;
  Double_t y1;
  Double_t x2;
  Double_t y2;
  Double_t x3;
  Double_t y3;
  Double_t dist;
  for(unsigned int i = 0;i<trg_vector.size();i++){
    tmp_time_l.push_back(0);
    spatial_cluster.clear();
    for(unsigned int j = i;j<trg_vector.size();j++){
      if(trg_vector.at(j).size()>0)
	tmp_time_l.at(i)++;
      else
	break;
    }
    cam_trg_vector.push_back(spatial_cluster);
  }
  for(unsigned int i = 0;i<tmp_time_l.size();i++){
    if(tmp_time_l.at(i)>=3){
      for(unsigned int j = i;j<(i+tmp_time_l.at(i));j++){
	for(unsigned int k = 0;k<trg_vector.at(j).size();k++){
	  cam_trg_vector.at(j).push_back(trg_vector.at(j).at(k));
	}
      }
    }
  }
  /*
    for(unsigned int j = 0;j<;j++){
      for(unsigned int k = 0;k<trg_vector.at(i).size();k++){
	x2=_simphist->get_pixel_vec().at(trg_vector.at(i).at(k)).x;
	y2=_simphist->get_pixel_vec().at(trg_vector.at(i).at(k)).y;
	dist = TMath::Sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
	//std::cout<<"dist "<<dist<<std::endl;
	if(dist<0.1)
	  n++;
      }
      if(n>=3)
	spatial_cluster.push_back(trg_vector.at(i).at(j));
    }
    //std::cout<<std::endl;
    cam_trg_vector.push_back(spatial_cluster);
  }
  //return cam_trg_vector;
  */
  return cam_trg_vector;
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
