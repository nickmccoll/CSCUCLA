#ifndef CSCUCLA_GEMCSCTUPLEANALYSIS_GEMGEOMETRY_H
#define CSCUCLA_GEMCSCTUPLEANALYSIS_GEMGEOMETRY_H
#include "GEMPartitionGeometry.h"
#include<vector>
#include<string>
namespace CSCGEMTuples {

class GEMGeometry {
public:
  GEMGeometry() : distBtwPartitions(0), nPhis(0), nPartitions(0), nStrips(0) { build();};
  GEMGeometry(float gemXScale) : distBtwPartitions(0), nPhis(0), nPartitions(0), nStrips(0)  { build(gemXScale);};
  unsigned int getNPartitions() const {return nPartitions;}
  unsigned int getNPhis() const {return nPhis;}
  unsigned int getNStrips() const {return nStrips;}
  void getVFATIEtaIPhi(int vFATChannelNumber, int & iEta, int & iPhi ) const;
  int geVFATChannelNumber(int iEta, int iPhi) const;
  unsigned int getNStripsPerVFAT() const {return nStrips/nPhis;}
  unsigned int getPartitionStripNumber(int iPhi, int vFatStripNumber) const; // go from 0-127 on each vfat to 1-384 for the entire eta partition

  const GEMPartitionGeometry& getPartition(int iEta) const {return partitions[iEta-1];}
  const Point2D& getPartitionCenter(int iEta) const {return partitionCenters[iEta-1];}

  Point2D getClusterPosition(int iEta,int firstStrip, int nStrips) const;
  Error2D getClusterError(int iEta,int firstStrip, int nStrips) const;
  ErrorPoint2D getClusterInfo(int iEta,int firstStrip, int nStrips) const;

  //General information
  Point2D getStripCenter(int iEta, int stripNumber) const;
  float getStripAngle(int iEta, float strip) const;
  float getStripPitch(int iEta, float yHeight) const;
  float getPartitionHeight(int iEta) const;
  float getPartitionTop(int iEta) const;
  float getPartitionBottom(int iEta) const;
  float getPartitionBottomEdge(int iEta) const;
  float getPartitionTopEdge(int iEta) const;

  //Finders
  //return eta partion local y belongs to
  // -1 means that it is above the widest vfat
  // nRows means it is below the narrowest one
  // (I know....)
  int findEtaPartition(float yValue) const;


  //Numbering information








private:
  void build(float gemXScale = -1);
  std::vector<GEMPartitionGeometry> partitions;
  std::vector<Point2D> partitionCenters;
  float distBtwPartitions;
  unsigned int nPhis;
  unsigned int nPartitions;
  unsigned int nStrips;
};


};
#endif
