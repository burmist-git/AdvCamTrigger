// initial code : https://github.com/james-yoo/DBSCAN
#ifndef DBSCAN_HH
#define DBSCAN_HH

//root
#include "TROOT.h"

//C, C++
#include <vector>
#include <cmath>

#define UNCLASSIFIED -1
#define CORE_POINT 1
#define BORDER_POINT 2
#define NOISE -2
#define SUCCESS 0
#define FAILURE -3

using namespace std;

typedef struct Point_
{
  Int_t pixel_id;
  float x, y, z;   // X, Y, Z position
  int ii;          // time bin
  int clusterID;   // clustered ID
}Point;

class DBSCAN {
public:    

  DBSCAN();
  DBSCAN(unsigned int minPts, float eps, vector<Point> points);
  ~DBSCAN();
  
  int run();
  int run(unsigned int minPts, float eps, vector<Point> points);
  vector<int> calculateCluster(Point point);
  int expandCluster(Point point, int clusterID);
  inline double calculateDistance(const Point& pointCore, const Point& pointTarget);
  
  int getTotalPointSize() {return m_pointSize;}
  int getMinimumClusterSize() {return m_minPoints;}
  int getEpsilonSize() {return m_epsilon;}
  
public:
  vector<Point> m_points;
  const void print_points_info() const;
  const Int_t get_nclusters() const;
  
private:    
  unsigned int m_pointSize;
  unsigned int m_minPoints;
  float m_epsilon;
};

#endif // DBSCAN_HH
