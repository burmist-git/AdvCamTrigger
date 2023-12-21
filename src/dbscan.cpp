//my
#include "dbscan.hh"

//root
#include "TROOT.h"
#include "TMath.h"

//C, C++
#include <iostream>
#include <assert.h>
#include <iomanip>
#include <stdlib.h>

DBSCAN::DBSCAN(){
}

DBSCAN::DBSCAN(unsigned int minPts, float eps, vector<Point> points){
  m_minPoints = minPts;
  m_epsilon = eps;
  m_points = points;
  m_pointSize = points.size();
}

DBSCAN::~DBSCAN(){;}

const Int_t DBSCAN::get_nclusters() const {
  Int_t nclus = 0;
  for(unsigned int k = 0; k<m_points.size(); k++){
    if(m_points.at(k).clusterID > 0)
      nclus++;
  }
  return nclus;
}

int DBSCAN::run(unsigned int minPts, float eps, vector<Point> points){
  m_minPoints = minPts;
  m_epsilon = eps;
  m_points = points;
  m_pointSize = points.size();
  return run();
}

int DBSCAN::run(){
  int clusterID = 1;
  vector<Point>::iterator iter;
  for(iter = m_points.begin(); iter != m_points.end(); ++iter){
    if ( iter->clusterID == UNCLASSIFIED ){
      if ( expandCluster(*iter, clusterID) != FAILURE ){
	clusterID += 1;
      }
    }
  }
  return 0;
}

int DBSCAN::expandCluster(Point point, int clusterID) {    
  vector<int> clusterSeeds = calculateCluster(point);
  if ( clusterSeeds.size() < m_minPoints ){
    point.clusterID = NOISE;
    return FAILURE;
  }
  else{
    int index = 0, indexCorePoint = 0;
    vector<int>::iterator iterSeeds;
    for( iterSeeds = clusterSeeds.begin(); iterSeeds != clusterSeeds.end(); ++iterSeeds){
      m_points.at(*iterSeeds).clusterID = clusterID;
      if (m_points.at(*iterSeeds).x == point.x && m_points.at(*iterSeeds).y == point.y && m_points.at(*iterSeeds).z == point.z ){
	indexCorePoint = index;
      }
      ++index;
    }
    clusterSeeds.erase(clusterSeeds.begin()+indexCorePoint);
    for( vector<int>::size_type i = 0, n = clusterSeeds.size(); i < n; ++i ){
      vector<int> clusterNeighors = calculateCluster(m_points.at(clusterSeeds[i]));
      if ( clusterNeighors.size() >= m_minPoints ){
	vector<int>::iterator iterNeighors;
	for ( iterNeighors = clusterNeighors.begin(); iterNeighors != clusterNeighors.end(); ++iterNeighors ){
	  if ( m_points.at(*iterNeighors).clusterID == UNCLASSIFIED || m_points.at(*iterNeighors).clusterID == NOISE ){
	    if ( m_points.at(*iterNeighors).clusterID == UNCLASSIFIED ){
	      clusterSeeds.push_back(*iterNeighors);
	      n = clusterSeeds.size();
	    }
	    m_points.at(*iterNeighors).clusterID = clusterID;
	  }
	}
      }
    }
    return SUCCESS;
  }
}

vector<int> DBSCAN::calculateCluster(Point point){
  int index = 0;
  vector<Point>::iterator iter;
  vector<int> clusterIndex;
  for( iter = m_points.begin(); iter != m_points.end(); ++iter){
    if ( calculateDistance(point, *iter) <= m_epsilon ){
      clusterIndex.push_back(index);
    }
    index++;
  }
  return clusterIndex;
}

inline double DBSCAN::calculateDistance(const Point& pointCore, const Point& pointTarget){
  return TMath::Sqrt(pow(pointCore.x - pointTarget.x,2)+pow(pointCore.y - pointTarget.y,2)+pow(pointCore.z - pointTarget.z,2));
}

const void DBSCAN::print_points_info() const {
  std::cout<<" nclusters : "<<get_nclusters()<<std::endl;
  for(unsigned int k = 0;k<m_points.size();k++)
    std::cout<<setw(10)<<m_points.at(k).x
	     <<setw(10)<<m_points.at(k).y
	     <<setw(10)<<m_points.at(k).z
      	     <<setw(10)<<m_points.at(k).ii
      	     <<setw(10)<<m_points.at(k).pixel_id
	     <<setw(10)<<m_points.at(k).clusterID<<std::endl;
}
