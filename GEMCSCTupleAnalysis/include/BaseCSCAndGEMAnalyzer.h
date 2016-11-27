#ifndef CSCUCLA_CSCGEMTUPLEANALYZER_BASECSCANDGEMANALYZER_H
#define CSCUCLA_CSCGEMTUPLEANALYZER_BASECSCANDGEMANALYZER_H
#include "BaseTupleAnalyzer.h"
#include "CSCInfo.h"
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <assert.h>
#include <iostream>

#include "GEMInfo.h"
#include "Segment.h"

namespace CSCGEMTuples {
class AnalyzeCSC : public BaseTupleAnalyzer{
public:
  AnalyzeCSC(std::string fileName, std::string treeName);
  virtual void analyze(int reportFrequency = 1000000);
  virtual void runAEvent() {buildSegments();}
  void buildSegments();

  std::vector<Segment> segments;

  EventInfo eventInfo;
  RecHitInfo recHitInfo;
  StripInfo  stripInfo;
  CompInfo  compInfo;
  WireInfo  wireInfo;
  LCTInfo  lctInfo;
  SegmentInfo  segmentInfo;
  CLCTInfo  clctInfo;
};



class AnalyzeGEM : public BaseTupleAnalyzer{
public:
  AnalyzeGEM(std::string fileName, std::string treeName,const GEMConfigInfo* gemInfo);
  virtual void runAEvent();
  Event * event;
  GEMInfo gemInfo;
};

class AnalyzeTMB : public BaseTupleAnalyzer{
public:
  AnalyzeTMB(std::string fileName, std::string treeName);
  virtual void runAEvent();

  std::vector<int> * gem_chamber  ;
  std::vector<int> * gem_partition;
  std::vector<int> * gem_BX       ;
  std::vector<int> * gem_strip    ;
  std::vector<int> * gem_nStrips  ;
  std::vector<int> * lct_chamber  ;
  std::vector<int> * lct_BX       ;
  std::vector<int> * lct_strip    ;
  std::vector<int> * lct_wg       ;
  std::vector<int> * lct_valid    ;


};


class AnalyzeBoth {
public:

  AnalyzeBoth(std::string cscFile, std::string gemFile, const GEMConfigInfo* gemInfo = new GEMConfigInfo);
  virtual ~AnalyzeBoth(){}

  virtual void runAEvent() {};
  void analyze(int reportFrequency = 1000000);

  AnalyzeCSC csc;
  AnalyzeGEM gem;

  //Easy accessors
  std::vector<Segment> * cscSegments;
  GEMGeometry * gemGeo;
  std::vector<GEMCluster>* gemClusters;



};
}
#endif
