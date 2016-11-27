
#include "../include/GEMGeometry.h"
#include <iostream>
#include <TString.h>
namespace CSCGEMTuples {

void GEMGeometry::getVFATIEtaIPhi(int vFATChannelNumber, int & iEta, int & iPhi ) const  {
  iPhi = vFATChannelNumber/nPartitions +1;
  iEta = nPartitions - (vFATChannelNumber % nPartitions);
}

int GEMGeometry::geVFATChannelNumber(int iEta, int iPhi) const{
  return nPartitions*(iPhi - 1) + nPartitions - iEta;
}
unsigned int GEMGeometry::getPartitionStripNumber(int iPhi, int vFatStripNumber) const {
  return vFatStripNumber + 1 + getNStripsPerVFAT()*(nPhis-iPhi);
}
Point2D GEMGeometry::getClusterPosition(int iEta,int firstStrip, int nStrips) const {
  return getPartition(iEta).clusterPosition(firstStrip,nStrips) + getPartitionCenter(iEta);
}
Error2D GEMGeometry::getClusterError(int iEta,int firstStrip, int nStrips) const {
  return getPartition(iEta).clusterError(firstStrip,nStrips);
}
ErrorPoint2D GEMGeometry::getClusterInfo(int iEta,int firstStrip, int nStrips) const {
  return ErrorPoint2D( getPartition(iEta).clusterPosition(firstStrip,nStrips) + getPartitionCenter(iEta),
      getPartition(iEta).clusterError(firstStrip,nStrips)
  );

}

Point2D GEMGeometry::getStripCenter(int iEta, int stripCenter) const {
  return getPartition(iEta).centerOfStrip(stripCenter) +   getPartitionCenter(iEta);
}

float GEMGeometry::getStripAngle(int iEta, float strip) const {
  return getPartition(iEta).stripAngle(strip);
}
float GEMGeometry::getPartitionHeight(int iEta) const {
  return getPartition(iEta).height();
}
float GEMGeometry::getPartitionTop(int iEta) const {
  return getPartitionCenter(iEta).y() + getPartition(iEta).height()/2;
}
float GEMGeometry::getPartitionBottom(int iEta) const {
  return getPartitionCenter(iEta).y() - getPartition(iEta).height()/2;
}
float GEMGeometry::getPartitionBottomEdge(int iEta) const {
  return getPartition(iEta).bottomEdgeSize();
}
float GEMGeometry::getPartitionTopEdge(int iEta) const {
  return getPartition(iEta).topEdgeSize();
}

int GEMGeometry::findEtaPartition(float yValue) const {
  if(yValue >=  getPartitionTop(1) ) return 0;
  for(unsigned int iR = 1; iR <= getNPartitions(); ++iR){
    if(yValue >=  getPartitionBottom(iR) ) return iR;
  }
  return getNPartitions() + 1;
}


void GEMGeometry::build(float gemXScale) {
  //Hard code paramters because generality is not needed!
  nPartitions = 8;
  distBtwPartitions = 0.05;
  nPhis = 3;
  nStrips = 384;
  std::vector<float> botEdge(nPartitions);
  std::vector<float> topEdge(nPartitions);
  std::vector<float> height(nPartitions);

  //Numbering in CMSSW goes from wide to narrow
  botEdge[0] = 20.5701; topEdge[0] = 22.2929; height[0] = 9.7025;
  botEdge[1] = 18.8429; topEdge[1] = 20.5657; height[1] = 9.7025;
  botEdge[2] = 17.4084; topEdge[2] = 18.8384; height[2] = 8.0535;
  botEdge[3] = 15.974 ; topEdge[3] = 17.404 ; height[3] = 8.0535;
  botEdge[4] = 14.7749; topEdge[4] = 15.9696; height[4] = 6.728 ;
  botEdge[5] = 13.5759; topEdge[5] = 14.7705; height[5] = 6.728 ;
  botEdge[6] = 12.5676; topEdge[6] = 13.5714; height[6] = 5.6535;
  botEdge[7] = 11.5593; topEdge[7] = 12.5631; height[7] = 5.6535;


//  botEdge[0] = 19.5003	; topEdge[0] = 20.97895; height[0] = 8.3275;
//  botEdge[1] = 18.0172	; topEdge[1] = 19.49585; height[1] = 8.3275;
//  botEdge[2] = 16.7648	; topEdge[2] = 18.01275; height[2] = 7.0285;
//  botEdge[3] = 15.51235 ; topEdge[3] = 16.76035; height[3] = 7.0285;
//  botEdge[4] = 14.44635 ; topEdge[4] = 15.5079;  height[4] = 5.9785;
//  botEdge[5] = 13.38035 ; topEdge[5] = 14.4419;  height[5] = 5.9785;
//  botEdge[6] = 12.4698	; topEdge[6] = 13.3759;  height[6] = 5.1030;
//  botEdge[7] = 11.5593	; topEdge[7] = 12.4654;  height[7] = 5.1030;

//
//  botEdge[0] = 19.5003  ; topEdge[0] = 20.97895; height[0] = 9.7025;
//  botEdge[1] = 18.0172  ; topEdge[1] = 19.49585; height[1] = 9.7025;
//  botEdge[2] = 16.7648  ; topEdge[2] = 18.01275; height[2] = 8.0535;
//  botEdge[3] = 15.51235 ; topEdge[3] = 16.76035; height[3] = 8.0535;
//  botEdge[4] = 14.44635 ; topEdge[4] = 15.5079;  height[4] = 6.728 ;
//  botEdge[5] = 13.38035 ; topEdge[5] = 14.4419;  height[5] = 6.728 ;
//  botEdge[6] = 12.4698  ; topEdge[6] = 13.3759;  height[6] = 5.6535;
//  botEdge[7] = 11.5593  ; topEdge[7] = 12.4654;  height[7] = 5.6535;

//
//  botEdge[0] = 18.8429; topEdge[0] = 20.5657; height[0] = 9.7025;
//  botEdge[1] = 17.4084; topEdge[1] = 18.8384; height[1] = 9.7025;
//  botEdge[2] = 15.974 ; topEdge[2] = 17.404 ; height[2] = 8.0535;
//  botEdge[3] = 14.7749; topEdge[3] = 15.9696; height[3] = 8.0535;
//  botEdge[4] = 13.5759; topEdge[4] = 14.7705; height[4] = 6.728 ;
//  botEdge[5] = 12.5676; topEdge[5] = 13.5714; height[5] = 6.728 ;
//  botEdge[6] = 11.5593; topEdge[6] = 12.5631; height[6] = 5.6535;
//  botEdge[7] = 11.5593; topEdge[7] = 12.5631; height[7] = 5.6535;
  if(gemXScale > 0)
  for(unsigned int iE = 0; iE < nPartitions; ++iE){
    topEdge[iE] *= gemXScale;
    botEdge[iE] *= gemXScale;
  }

  for(unsigned int iL = 0; iL < nPartitions; ++iL){
    partitions.emplace_back(botEdge[iL],topEdge[iL],height[iL],nStrips);
  }

  float totalHeight = 0;
  partitionCenters.resize(nPartitions);
  //First build it bottom up!
  for(int iL = nPartitions -1; iL >= 0; --iL){
    partitionCenters[iL].set(0,totalHeight + height[iL]);
    totalHeight += height[iL]*2;
    if(iL) totalHeight += distBtwPartitions;
  }

  //Print in tabel as shown in the twiki
//  for(int iL = nLay -1; iL >= 0; --iL){
//    std::cout << 10.0 *  partitions[iL].bottomEdgeSize() <<","<<10.0 *  partitions[iL].topEdgeSize() <<","
//        <<10.0 * partionCenters[iL].y() -   10.0 *partitions[iL].height()/2<<","<<10.0 * partionCenters[iL].y() +   10.0 *partitions[iL].height()/2<<std::endl;
//  }

  //Now adjust it so that it is centered at 0
  Point2D displacement(0.,totalHeight/2);
  for(unsigned int iL = 0; iL < nPartitions; ++iL){
    partitionCenters[iL] -= displacement;
  }
//  std::cout << std::endl;
  //Print in tabel as we want it!
//  for(int iL = nLay -1; iL >= 0; --iL){
//    std::cout << 10.0 *  partitions[iL].bottomEdgeSize() <<","<<10.0 *  partitions[iL].topEdgeSize() <<","
//        <<10.0 * partionCenters[iL].y() -  10.0 * partitions[iL].height()/2<<","<<10.0 * partionCenters[iL].y() +   10.0 *partitions[iL].height()/2<<std::endl;
//  }


}



}
