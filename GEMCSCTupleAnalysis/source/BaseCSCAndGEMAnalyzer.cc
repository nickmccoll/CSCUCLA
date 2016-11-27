#include "../include/BaseCSCAndGEMAnalyzer.h"
#include "../../gem-light-dqm/gemtreewriter/include/Event.h"

namespace CSCGEMTuples {
    AnalyzeCSC::AnalyzeCSC(std::string fileName, std::string treeName) : BaseTupleAnalyzer(fileName,treeName){
        eventInfo.load(this);
        recHitInfo.load(this);
        stripInfo.load(this);
        compInfo.load(this);
        wireInfo.load(this);
        lctInfo.load(this);
        segmentInfo.load(this);
        clctInfo.load(this);
    }

    void AnalyzeCSC::analyze(int reportFrequency) {
        while(nextEvent(reportFrequency)){runAEvent();eventNumber++;}
    }

    void AnalyzeCSC::buildSegments() {
      segments.clear();
      segments.resize(segmentInfo.segment_cov_x->size());
      for(unsigned int iS = 0; iS < segments.size(); ++iS){
        segments[iS].setCoords(segmentInfo.segment_pos_x->at(iS),segmentInfo.segment_pos_y->at(iS),segmentInfo.segment_dxdz->at(iS),segmentInfo.segment_dydz->at(iS));
        segments[iS].setCov(segmentInfo.segment_cov_x->at(iS),segmentInfo.segment_cov_x_y->at(iS),segmentInfo.segment_cov_dxdz_x->at(iS),segmentInfo.segment_cov_dydz_x->at(iS),
            segmentInfo.segment_cov_y->at(iS),segmentInfo.segment_cov_dxdz_y->at(iS),segmentInfo.segment_cov_dydz_y->at(iS),
            segmentInfo.segment_cov_dxdz->at(iS),segmentInfo.segment_cov_dxdz_dydz->at(iS),segmentInfo.segment_cov_dydz->at(iS));
      }
    }

    AnalyzeGEM::AnalyzeGEM(std::string fileName, std::string treeName,const GEMConfigInfo * gemInfo)
        : BaseTupleAnalyzer(fileName,treeName),
        gemInfo(*gemInfo){
            event = new Event();
            setBranchAddress("GEMEvents",&event,true);
        }
    void AnalyzeGEM::runAEvent() {
        gemInfo.build(event);
    }

    AnalyzeBoth::AnalyzeBoth(std::string cscFile, std::string gemFile,const GEMConfigInfo* gemInfo) :
        csc(cscFile,"CSCDigiTree"),
        gem(gemFile,"GEMtree", gemInfo),
        cscSegments(&csc.segments),
        gemGeo(&gem.gemInfo.geo),
        gemClusters(&gem.gemInfo.clusters)
    {
    }

    void AnalyzeBoth::analyze(int reportFrequency) {
        gem.eventNumber +=1;
        while(csc.nextEvent(reportFrequency) && gem.nextEvent(reportFrequency)){
            if(csc.eventInfo.Event_EventNumber != gem.event->GetEventNumber()){
                gem.eventNumber +=1;
                continue;
            }
            csc.runAEvent();
            gem.runAEvent();
            runAEvent();
            csc.eventNumber++;
            gem.eventNumber++;
        }
    }




    AnalyzeTMB::AnalyzeTMB(std::string fileName, std::string treeName)
        : BaseTupleAnalyzer(fileName,treeName) {

       gem_chamber   = new std::vector<int> ;
       gem_partition = new std::vector<int> ;
       gem_BX        = new std::vector<int> ;
       gem_strip     = new std::vector<int> ;
       gem_nStrips   = new std::vector<int> ;
       lct_chamber   = new std::vector<int> ;
       lct_BX        = new std::vector<int> ;
       lct_strip     = new std::vector<int> ;
       lct_wg        = new std::vector<int> ;
       lct_valid     = new std::vector<int> ;
      setBranchAddress("gem_chamber"   ,&gem_chamber  , true);
      setBranchAddress("gem_partition" ,&gem_partition, true);
      setBranchAddress("gem_BX"        ,&gem_BX       , true);
      setBranchAddress("gem_strip"     ,&gem_strip    , true);
      setBranchAddress("gem_nStrips"   ,&gem_nStrips  , true);
      setBranchAddress("lct_chamber"   ,&lct_chamber  , true);
      setBranchAddress("lct_BX"        ,&lct_BX       , true);
      setBranchAddress("lct_strip"     ,&lct_strip    , true);
      setBranchAddress("lct_wg"        ,&lct_wg       , true);
      setBranchAddress("lct_valid"     ,&lct_valid    , true);

        }
    void AnalyzeTMB::runAEvent() {
    }



}
