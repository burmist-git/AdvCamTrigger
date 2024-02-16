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
  Int_t number_of_pixels;
  Int_t n_t_ii;
  Int_t number_of_CORE_POINT;
  Int_t number_of_BORDER_POINT;
  Double_t mean_x;
  Double_t mean_y;
  Double_t mean_time_ii;
  vector<Int_t> pixel_id_v;
  vector<Int_t> time_ii_v;
  cluster_info(){
    clusterID = -999;
    number_of_points = 0;
    number_of_pixels = 0;
    n_t_ii = 0;
    number_of_CORE_POINT = 0;
    number_of_BORDER_POINT = 0;
    mean_x       = 0.0;
    mean_y       = 0.0;
    mean_time_ii = 0.0;
  }
  void print_cluster_info_header() const {
    std::cout<<setw(13)<<"clusterID"
	     <<setw(19)<<"number_of_points"
      	     <<setw(19)<<"number_of_pixels"
	     <<setw(8)<<"n_t_ii"
	     <<setw(24)<<"number_of_CORE_POINT"
	     <<setw(24)<<"number_of_BORDER_POINT"
	     <<setw(14)<<"mean_x"
	     <<setw(14)<<"mean_y"
	     <<setw(14)<<"mean_time_ii"
	     <<std::endl;
  }
  void print_cluster_info() const {
    std::cout<<setw(13)<<clusterID
	     <<setw(19)<<number_of_points
      	     <<setw(19)<<number_of_pixels
	     <<setw(8)<<n_t_ii
	     <<setw(24)<<number_of_CORE_POINT
	     <<setw(24)<<number_of_BORDER_POINT
	     <<setw(14)<<mean_x
	     <<setw(14)<<mean_y
	     <<setw(14)<<mean_time_ii
	     <<std::endl;    
  }
  static void bubbleSort(vector<Int_t> &a){
    bool swapp = true;
    while(swapp){
      swapp = false;
      for( unsigned int i = 0; i < a.size()-1; i++){
	if( a.at(i)>a.at(i+1)){
	  //
	  Int_t tmp;
	  tmp = a.at(i);
	  //
	  a.at(i) = a.at(i+1);
	  a.at(i+1) = tmp;
	  swapp = true;
	}
      }
    }
  }
  void get_n_t_ii(){
    cluster_info::bubbleSort(time_ii_v);
    n_t_ii = time_ii_v.back() - time_ii_v.front();
  }
  void get_number_of_pixels(){
    cluster_info::bubbleSort(pixel_id_v);
    Int_t npx = (Int_t)pixel_id_v.back() + 1;
    if(npx>0){
      Int_t *tmpint = new Int_t[npx];
      //
      for(Int_t i = 0;i<npx;i++)
	tmpint[i] = 0;
      //
      for( unsigned int k = 0; k < pixel_id_v.size(); k++)
	tmpint[pixel_id_v.at(k)] = 1;
      for( Int_t i = 0; i < npx; i++)
	number_of_pixels += tmpint[i];
    }
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

struct point_kdis {
  point p;
  vector<Double_t> dists;
  static void bubbleSort(vector<Double_t> &a){
    bool swapp = true;
    while(swapp){
      swapp = false;
      for( unsigned int i = 0; i < a.size()-1; i++){
	if( a.at(i)>a.at(i+1)){
	  //
	  Double_t tmp;
	  tmp = a.at(i);
	  //
	  a.at(i) = a.at(i+1);
	  a.at(i+1) = tmp;
	  swapp = true;
	}
      }
    }
  }
};
  
class dbscan {
public:    

  dbscan();
  ~dbscan();
  
  Int_t run(unsigned int minPts, Double_t eps, vector<point> point_v);
  void set_points(vector<point> point_v);

  vector<unsigned int> regionQuery(const point &p);
  bool expandCluster(point &p, int clusterID);

  inline Double_t calculateDistance(const point& pointCore, const point& pointTarget);  
  const void print_points_info() const;
  void print_points_info(vector<point> point_v);
  static void print_points_info_st(const vector<point> point_v);
  Int_t get_nclusters() const;

  void get_cluster_stats();
  static void print_cluster_stats(vector<cluster_info> clusters_v);
  void print_cluster_stats() const;
  Int_t get_number_of_NOISE() const;
  void print_single_point_info( const point& p) const;
  
  inline const vector<point>& get_points_v() const {return _points_v;};
  inline const vector<cluster_info> get_clusters_v() const {return _clusters_v;};

  inline const Int_t get_dbscan_run_time_musec() const {return _dbscan_run_time_musec;};
  
  void clear();

  vector<Double_t> build_k_dist_graph(Int_t kk);

  Int_t _dbscan_run_time_musec;
  
private:

  vector<point> _points_v;
  vector<cluster_info> _clusters_v;
  vector<point_kdis> _points_kdis_v;
  
  void erase_point_from_seeds(const point &p, vector<unsigned int> &seeds_v);

  unsigned int _minPts;
  Double_t _eps;
  
};

#endif // dbscan_hh
