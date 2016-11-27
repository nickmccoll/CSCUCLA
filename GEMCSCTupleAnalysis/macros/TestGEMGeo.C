
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../include/BaseCSCAndGEMAnalyzer.h"
#include "include/GEMPlottingInfo.h"
#include "../include/GEMGeometry.h"
#include "../include/Segment.h"

#include<iostream>

using namespace std;
using namespace CSCGEMTuples;

class Analyze : public AnalyzeBoth {
public:
  int nEvt = 0;
  Analyze(std::string cscFileName, std::string gemFileName,const GEMConfigInfo* info) : AnalyzeBoth(cscFileName,gemFileName,info)
  {

    //Test projection
    plotter.book2D("new_proj_SegmentMap"            ,";x[cm];y[cm]"     ,90,-90,90,90,-90,90);
    plotter.book2D("move_proj_SegmentMap"            ,";x[cm];y[cm]"     ,90,-90,90,90,-90,90);
    plotter.book2D("rot_proj_SegmentMap"            ,";x[cm];y[cm]"     ,90,-90,90,90,-90,90);
    plotter.book2D("tiltx_proj_SegmentMap"            ,";x[cm];y[cm]"     ,90,-90,90,90,-90,90);
    plotter.book2D("tilty_proj_SegmentMap"            ,";x[cm];y[cm]"     ,90,-90,90,90,-90,90);
    //end projection

    plotter.book1D("clus_x"            ,";X [cm]"     ,100,-90,90);
    plotter.book1D("clus_xError"            ,";error [cm]"     ,500,0,10);
    plotter.book1D("clus_yError"            ,";error [cm]"     ,500,0,10);
    plotter.book2D("proj_SegmentMap"            ,";x[cm];y[cm]"     ,90,-90,90,90,-90,90);
    plotter.book1D("proj_xError"            ,";error [cm]"     ,500,0,10);
    plotter.book1D("proj_yError"            ,";error [cm]"     ,500,0,10);

    plotter.book1D("delta_xError"            ,";error [cm]"     ,500,-5,5);
    plotter.book1D("delta_yError"            ,";error [cm]"     ,500,-20,20);

    plotter.book1D("deltaX"            ,";#Delta X [cm]"     ,500,-90,90);
    plotter.book1D("deltaY"            ,";#Delta Y [cm]"     ,500,-90,90);

    plotter.book1D("deltaXoSig"            ,";#Delta X / #sigma"     ,500,-10,10);
    plotter.book1D("deltaYoSig"            ,";#Delta Y / #sigma"     ,500,-10,10);

    for(unsigned int iV =  0 ; iV < gemGeo->getNPartitions() * gemGeo->getNPhis(); ++iV){
      TString name = TString::Format("deltaXoSig_%u",iV);
      plotter.book1D(name            ,";#Delta X / #sigma"     ,1000,-10,10);
    }

    for(unsigned int iE =  1 ; iE <= gemGeo->getNPartitions(); ++iE){
      TString name1 = TString::Format("deltaYoSig_eta_%u",iE);
      TString name2 = TString::Format("deltaY_eta_%u",iE);
      plotter.book1D(name1              ,";(y_{g} - y_{p})/#sigma(y_{g} - y_{p});a.u."     ,250,-5,5);
      plotter.book1D(name2              ,";y_{g} - y_{p} [cm];a.u."     ,100,-30,30);
    }

    if (true){
      double xBins[] = {-22,-16.5,-11,-5.5,0,5.5,11,16.5,22};
      for(unsigned int iX =  0 ; iX < 8; ++iX){

        TString name1 = TString::Format("deltaXoSig_x_%.0f_%.0f",xBins[iX]*10,xBins[iX+1]*10);
        TString name2 = TString::Format("deltaX_x_%.0f_%.0f",xBins[iX]*10,xBins[iX+1]*10);
        TString name3 = TString::Format("deltaXoSig_const_x_%.0f_%.0f",xBins[iX]*10,xBins[iX+1]*10);
        TString name4 = TString::Format("deltaX_const_x_%.0f_%.0f",xBins[iX]*10,xBins[iX+1]*10);

        TString name5 = TString::Format("deltaXoSig_const_oneStrip_x_%.0f_%.0f",xBins[iX]*10,xBins[iX+1]*10);
        TString name6 = TString::Format("deltaX_const_oneStrip_x_%.0f_%.0f",xBins[iX]*10,xBins[iX+1]*10);

        plotter.book1D(name1              ,";(x_{g} - x_{p})/#sigma(x_{g} - x_{p});a.u."     ,100,-15,15);
        plotter.book1D(name2              ,";x_{g} - x_{p} [cm];a.u."     ,60,-5,5);
        plotter.book1D(name3              ,";(x_{g} - x_{p})/#sigma(x_{g} - x_{p});a.u."     ,100,-15,15);
        plotter.book1D(name4              ,";x_{g} - x_{p} [cm];a.u."     ,60,-5,5);
        plotter.book1D(name5              ,";(x_{g} - x_{p})/#sigma(x_{g} - x_{p});a.u."     ,100,-15,15);
        plotter.book1D(name6              ,";x_{g} - x_{p} [cm];a.u."     ,60,-5,5);
      }
    }

        plotter.book2D("proj_x_r1_1"            ,";x [cm];y [cm]"     ,90,-90,90,90,-90,90);
        plotter.book2D("proj_x_r1_384"            ,";x [cm];y [cm]"     ,90,-90,90,90,-90,90);
        plotter.book2D("proj_x_r4_20"            ,";x [cm];y [cm]"     ,90,-90,90,90,-90,90);
  }
  virtual  ~Analyze() {};

  void write(TString outFileName){cout << nEvt << endl; plotter.write(outFileName);}


  void testProjection(const Segment& seg) {
    //test projection
    auto p = seg.project(0,0,30,0);
    plotter.get2D("new_proj_SegmentMap")->Fill(p.x(),p.y());

    p = seg.project(10,10,30,0);
    plotter.get2D("move_proj_SegmentMap")->Fill(p.x(),p.y());

    p = seg.project(0,0,30,3.14/4.0);
    plotter.get2D("rot_proj_SegmentMap")->Fill(p.x(),p.y());

    p = seg.project2(0,0,30,3.14/4.0,0,0);
    plotter.get2D("tiltx_proj_SegmentMap")->Fill(p.x(),p.y());

    p = seg.project2(0,0,30,0,3.14/4.0,0);
    plotter.get2D("tilty_proj_SegmentMap")->Fill(p.x(),p.y());
  }

  virtual void runAEvent() {
    if(!pureSample(csc)) return;
    nEvt++;
    testProjection(cscSegments->at(0));

    if(gemClusters->size() != 1) return;

    auto gemClus = gemClusters->at(0).localCoords();
    const double clusX = gemClus.x();
    const double clusY = gemClus.y();

//    cout <<transRow <<" "<< gemClus.x() <<" "<<" "<<transStrip <<" "<< gem.gemInfo.clusters[0].getNStrips()<<" " << gemClus.y()<<  endl;

    plotter.get1D("clus_x")->Fill(gemClus.x());
    plotter.get1D("clus_xError")->Fill(gemClus.error_x());
    plotter.get1D("clus_yError")->Fill(gemClus.error_y());
    auto cscSeg = cscSegments->at(0).project(-0.66767,1.1331,33.2863,-0.0274707);
    double origX = gemClus.x();

    plotter.get2D("proj_SegmentMap")->Fill(cscSeg.x(),cscSeg.y());

    plotter.get1D("proj_xError")->Fill(cscSeg.error_x());
    plotter.get1D("proj_yError")->Fill(cscSeg.error_y());

    gemClus -= cscSeg;

    plotter.get1D("delta_xError")->Fill(gemClus.error_x());
    plotter.get1D("delta_yError")->Fill(gemClus.error_y());

    plotter.get1D("deltaX")->Fill(gemClus.x());
    plotter.get1D("deltaY")->Fill(gemClus.y());

    plotter.get1D("deltaXoSig")->Fill(gemClus.x()/gemClus.error_x());
    plotter.get1D("deltaYoSig")->Fill(gemClus.y()/gemClus.error_y());


    TString name = TString::Format("deltaXoSig_%u",gemGeo->geVFATChannelNumber(gem.gemInfo.vFats[0].iEta,gem.gemInfo.vFats[0].iPhi));
    plotter.get1D(name)->Fill(gemClus.x()/gemClus.error_x());

    if( true) {
      unsigned int iE = gemClusters->at(0).iEta;
      TString name1 = TString::Format("deltaYoSig_eta_%u",iE);
      TString name2 = TString::Format("deltaY_eta_%u",iE);
      plotter.get1D(name2)->Fill(gemClus.y());
      plotter.get1D(name1)->Fill(gemClus.y()/gemClus.error_y());
    }



    if (true){
      double xBins[] = {-22,-16.5,-11,-5.5,0,5.5,11,16.5,22};
      for(unsigned int iX =  0 ; iX < 8; ++iX){
        if(clusX >= xBins[iX+1] ) continue;
        TString name1 = TString::Format("deltaXoSig_x_%.0f_%.0f",xBins[iX]*10,xBins[iX+1]*10);
        TString name2 = TString::Format("deltaX_x_%.0f_%.0f",xBins[iX]*10,xBins[iX+1]*10);
        TString name3 = TString::Format("deltaXoSig_const_x_%.0f_%.0f",xBins[iX]*10,xBins[iX+1]*10);
        TString name4 = TString::Format("deltaX_const_x_%.0f_%.0f",xBins[iX]*10,xBins[iX+1]*10);
        TString name5 = TString::Format("deltaXoSig_const_oneStrip_x_%.0f_%.0f",xBins[iX]*10,xBins[iX+1]*10);
        TString name6 = TString::Format("deltaX_const_oneStrip_x_%.0f_%.0f",xBins[iX]*10,xBins[iX+1]*10);

        plotter.get1D(name1)->Fill(gemClus.x()/gemClus.error_x());
        plotter.get1D(name2)->Fill(gemClus.x());

            if(TMath::Abs(csc.segmentInfo.segment_dxdz->at(0)) < 0.008777){
              plotter.get1D(name3)->Fill(gemClus.x()/gemClus.error_x());
              plotter.get1D(name4)->Fill(gemClus.x());
              if(gemClusters->at(0).getNStrips() == 1){
                plotter.get1D(name5)->Fill(gemClus.x()/gemClus.error_x());
                plotter.get1D(name6)->Fill(gemClus.x());
              }
            }


        break;

      }
    }

      if(gemClusters->at(0).iEta==1){
        if(gemClusters->at(0).getFirstStrip()==1) plotter.get2D("proj_x_r1_1")->Fill(cscSeg.x(),cscSeg.y());
        if(gemClusters->at(0).getFirstStrip()==384) plotter.get2D("proj_x_r1_384")->Fill(cscSeg.x(),cscSeg.y());
      }
      if(gemClusters->at(0).iEta==4){
        if(gemClusters->at(0).getFirstStrip()==20) plotter.get2D("proj_x_r4_20")->Fill(cscSeg.x(),cscSeg.y());
      }


  }

  HistGetter plotter;
};

#endif

void TestGEMGeo(std::string cscfileName="csc_forsync.root",std::string gemfilename = "gem_forsync.root",std::string outFileName = "plotCSCAndGEM_out.root"){
  GEMConfigInfo info;
  info.vFATChanMapName       = "slot_table_904_june09.csv";
  info.gemXScale = -1;
  Analyze a(cscfileName,gemfilename,&info);
  a.analyze();
  a.write(outFileName);
//  a.draw();
}
