#ifndef CSCUCLA_CSCDIGITUPLES_ANALYZETUPLES_GEMINFO_H
#define CSCUCLA_CSCDIGITUPLES_ANALYZETUPLES_GEMINFO_H

#include "GEMGeometry.h"
#include "GEMChannelMapping.h"
#include "GEMConfigInfo.h"
#include "GEMConfigInfo.h"
#include <vector>

class Event;

namespace CSCGEMTuples {
class VFAT {
public:
  VFAT(): iEta(-1), iPhi(-1) {}
  VFAT(int iEta, int iPhi, std::vector<int>& onStrips);
  int nStrips() const {return  strips.size();}
  int iEta;
  int iPhi;
  std::vector<int> strips;
};

class GEMCluster {
public:
  GEMCluster(int iEta, int idx,  int nStrips) : iEta(iEta), idx(idx),  nStrips(nStrips) {}
  int getNStrips() const {return nStrips;}
  int getFirstStrip() const {return idx;}
  int getLastStrip() const {return idx + nStrips - 1;}

  //Local coordinates accessors
  ErrorPoint2D& localCoords() { return lc;}
  const ErrorPoint2D& localCoords() const { return lc;}
  double x() const {return lc.x();}
  double y() const {return lc.y();}
  double error_x() const {return lc.error_x();}
  double error_y() const {return lc.error_y();}

  int iEta= -1;
  int idx=-1; //In partition strip number (1 -> 384)
  int nStrips = 0;
  ErrorPoint2D lc; //Need to set manually
};
class GEMInfo {
public:
  GEMInfo(const GEMConfigInfo& gemInfo);
  void build(Event * event);
  std::vector<VFAT> vFats;
  std::vector<GEMCluster> clusters;
  int BX = -1;
  int evtN = -1;
  GEMGeometry geo;
  GEMChannelMapping channelMapping;

private:
  void getClusters();

};
}

#endif
