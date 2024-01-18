//https://cdn.aaai.org/KDD/1996/KDD96-037.pdf
#ifndef dbscan_hh
#define dbscan_hh

//root
#include "TROOT.h"

//C, C++
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

#define CORE_POINT 1
#define BORDER_POINT 2
#define UNCLASSIFIED -1
#define NOISE -2

using namespace std;

struct cluster_info {
  Int_t clusterID;
  Int_t number_of_points;
  Int_t number_of_CORE_POINT;
  Int_t number_of_BORDER_POINT;
  cluster_info(){
    clusterID = -999;
    number_of_points = 0;
    number_of_CORE_POINT = 0;
    number_of_BORDER_POINT = 0;
  }
  void print_cluster_info_header() const {
    std::cout<<setw(30)<<"clusterID"
	     <<setw(30)<<"number_of_points"
	     <<setw(30)<<"number_of_CORE_POINT"
	     <<setw(30)<<"number_of_BORDER_POINT"
	     <<std::endl;
  }
  void print_cluster_info() const {
    std::cout<<setw(30)<<clusterID
	     <<setw(30)<<number_of_points
	     <<setw(30)<<number_of_CORE_POINT
	     <<setw(30)<<number_of_BORDER_POINT
	     <<std::endl;    
  }
};

struct point {
  Int_t pixel_id;    //
  Double_t x, y, z;  // x, y, z position
  Int_t time_ii;     // time bin
  Int_t clusterID;   // clustered ID
  Int_t point_type;  // type (UNCLASSIFIED|CORE_POINT|BORDER_POINT|define NOISE)
  Int_t point_id;    //
  point(){
    pixel_id = -999;
    x = -999.0;
    y = -999.0;
    z = -999.0;
    time_ii=-999;
    clusterID = -999;
    point_type = UNCLASSIFIED;
    point_id = -999;
  }
};

class dbscan {
public:    

  dbscan();
  ~dbscan();
  
  Int_t run(unsigned int minPts, Double_t eps, vector<point> point_v);

  vector<unsigned int> regionQuery(point p);
  bool expandCluster(point &p, int clusterID);

  inline Double_t calculateDistance(const point& pointCore, const point& pointTarget);  
  const void print_points_info() const;
  void print_points_info(vector<point> point_v);
  static void print_points_info_st(const vector<point> point_v);
  Int_t get_nclusters() const;

  void get_cluster_stats();
  void print_cluster_stats() const;
  Int_t get_number_of_NOISE() const;

  inline const vector<point> get_points_v() const {return _points_v;};

  void clear();
  
private:

  vector<point> _points_v;
  vector<cluster_info> _clusters_v;
  
  unsigned int _minPts;
  Double_t _eps;
  
};

#endif // dbscan_hh
