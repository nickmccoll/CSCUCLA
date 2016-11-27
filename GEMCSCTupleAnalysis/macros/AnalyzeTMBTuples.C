
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../include/BaseCSCAndGEMAnalyzer.h"
#include "include/HistGetter.h"

using namespace std;
using namespace CSCGEMTuples;

class Analyzer : public AnalyzeTMB{
public:
  Analyzer(std::string fileName, std::string treeName) : AnalyzeTMB(fileName,treeName){
  }





  virtual void runAEvent() {
    static int nEntries = tree->GetEntries();
    if(eventNumber > 650000 )return;
    if(lct_chamber->size() == 0 ) return;
    if(lct_chamber->size() > 1 && lct_valid->at(1)) return;
    plotter.getOrMake1D("nEvents", ";# of events", 1,0,2)->Fill(1);


    auto convP = [] (int in) -> int {
//      return in;
        return -1 * (in -7);
    };
    auto toPad = [] (int in) -> int {
      in /= 8;
        return -1 * (in -23);
//      return in;
    };


    vector<int> nClusters(20);
    int nN = 0;
    int nGT8 = 0;
    int nNoiseVFAT = 0;
    int nGood = 0;
    for(unsigned int iG =0; iG < gem_chamber->size(); ++iG){
      if(gem_chamber->at(iG) == 1){continue;}
      int pad = toPad(gem_strip->at(iG));
      int part = convP(gem_partition->at(iG));
      if(pad == 16 && part == 2 &&  gem_BX->at(iG) >= 5) nN++;
      if(pad >= 16 && part == 2 &&  gem_BX->at(iG) >= 2 && gem_BX->at(iG) <= 4 ) nNoiseVFAT++;
      if( gem_BX->at(iG) >= 2 && gem_BX->at(iG) <= 4 ) nGood++;
      if(gem_BX->at(iG) >= 8) nGT8++;
      nClusters[gem_BX->at(iG)]++;
    }
    int nF = 0;
    int b1 = -1;
    int b2 = -1;
    for(unsigned int iB = 0; iB < nClusters.size(); ++iB){
      if(!nClusters[iB]) continue;
      nF++;
      if(b1 < 0) b1 = iB;
      else if (b2 < 0) b2 =iB;
    }


    plotter.getOrMake2D("incl_HS_WG", ";HS # ; WG #", 200,-0.5,199.5,50,-0.5,49.5)->Fill(lct_strip->at(0), lct_wg->at(0));
    if(nGT8) plotter.getOrMake2D("withBXGT8_HS_WG", ";HS # ; WG #", 50,-0.5,199.5,25,-0.5,49.5)->Fill(lct_strip->at(0), lct_wg->at(0));;
    int hs = lct_strip->at(0); int wg = lct_wg->at(0);
    int n0 = 0;
    int n1 = 0;
    for(unsigned int iG =0; iG < gem_chamber->size(); ++iG){
      if(gem_chamber->at(iG) == 1){
        n1++;
        continue;
      }
      n0++;
      int pad = toPad(gem_strip->at(iG));
      int part = convP(gem_partition->at(iG));



      if(nF > 1) plotter.getOrMake1D("doubleBX_gemBX", ";# of BX relative to ALCT data; # of GEM clusters", 20,-.5,19.5)->Fill(gem_BX->at(iG));
      plotter.getOrMake1D("gemBX", ";# of BX relative to ALCT data; # of GEM clusters", 20,-.5,19.5)->Fill(gem_BX->at(iG));
      plotter.getOrMake1D("gemClusterSize", ";# strips per cluster", 400,-.5,399.5)->Fill(gem_nStrips->at(iG));
      plotter.getOrMake2D("incl_pad_x_partition", ";pad # ; partition #", 24,-0.5,23.5,8,-0.5,7.5)->Fill(pad,part);

      plotter.getOrMake2D(TString::Format("incl_BX%i_pad_x_partition",gem_BX->at(iG)), ";pad # ; partition #", 24,-0.5,23.5,8,-0.5,7.5)->Fill(pad,part);

      if(gem_BX->at(iG) >= 8) plotter.getOrMake2D("incl_BX8to19_pad_x_partition", ";pad # ; partition #", 24,-0.5,23.5,8,-0.5,7.5)->Fill(pad,part);
      if(gem_BX->at(iG) >= 8 && gem_BX->at(iG) <= 15) plotter.getOrMake2D("incl_BX8to15_pad_x_partition", ";pad # ; partition #", 24,-0.5,23.5,8,-0.5,7.5)->Fill(pad,part);
      if(gem_BX->at(iG) >= 2 && gem_BX->at(iG) <= 4) plotter.getOrMake2D("incl_BX2to4_pad_x_partition", ";pad # ; partition #", 24,-0.5,23.5,8,-0.5,7.5)->Fill(pad,part);
      if(hs >= 20 && hs <= 40 && wg >= 35 && wg <= 40 ) if(gem_BX->at(iG) >= 2 && gem_BX->at(iG) <= 4) plotter.getOrMake2D("hs20to40wg35to40_BX2to4_pad_x_partition", ";pad # ; partition #", 24,-0.5,23.5,8,-0.5,7.5)->Fill(pad,part);
      if(hs >= 190 && hs <= 200 && wg >= 1 && wg <= 5 ) if(gem_BX->at(iG) >= 2 && gem_BX->at(iG) <= 4) plotter.getOrMake2D("hs190to200wg1to5_BX2to4_pad_x_partition", ";pad # ; partition #", 24,-0.5,23.5,8,-0.5,7.5)->Fill(pad,part);

      if(nF > 1 ) if(gem_BX->at(iG) >= 2 && gem_BX->at(iG) <= 4) plotter.getOrMake2D("doubleBX_BX2to4_pad_x_partition", ";pad # ; partition #", 24,-0.5,23.5,8,-0.5,7.5)->Fill(pad,part);
      if(nN)       if(gem_BX->at(iG) >= 2 && gem_BX->at(iG) <= 4) plotter.getOrMake2D("noiseStrip_BX2to4_pad_x_partition", ";pad # ; partition #", 24,-0.5,23.5,8,-0.5,7.5)->Fill(pad,part);
      if(!(pad ==16 && part == 2) ) plotter.getOrMake1D("maskNoiseStrip_gemBX", ";# of BX relative to ALCT data; # of GEM clusters", 20,-.5,19.5)->Fill(gem_BX->at(iG));
      if((pad ==16 && part == 2) ) plotter.getOrMake1D("onlyNoiseStrip_gemBX", ";# of BX relative to ALCT data; # of GEM clusters", 20,-.5,19.5)->Fill(gem_BX->at(iG));
      if(nGT8) if(gem_BX->at(iG) >= 2 && gem_BX->at(iG) <= 4) plotter.getOrMake2D("wBXGT8_BX2to4_pad_x_partition", ";pad # ; partition #", 24,-0.5,23.5,8,-0.5,7.5)->Fill(pad,part);
      if(!nNoiseVFAT) plotter.getOrMake1D("maskNoiseVFat_gemBX", ";# of BX relative to ALCT data; # of GEM clusters", 20,-.5,19.5)->Fill(gem_BX->at(iG));
      plotter.getOrMake2D("incl_event_bx", ";event # ; BX #", 1000,0,nEntries,20,-.5,19.5)->Fill(eventNumber,gem_BX->at(iG));
      plotter.getOrMake2D("gem_BX_v_nGood", "; BX # ;N Good", 20,-.5,19.5,20,-.5,19.5)->Fill(gem_BX->at(iG),nGood);

      bool shouldKill =false;
      bool shouldKill2 =false;
      auto getVF = [](int s)-> int {
        if(s <=7 ) return 0;
        if(s <=15 ) return 1;
         return 2;
      };
      for(unsigned int iG2 =0; iG2 < iG; ++iG2){
        if(gem_chamber->at(iG2) == 1) continue;
                int pad2 = toPad(gem_strip->at(iG2));
                int part2 = convP(gem_partition->at(iG2));
                if(pad2 == pad && part2 == part ) shouldKill = true;
                if(getVF(pad2) == getVF(pad) && part2 == part && gem_BX->at(iG2) != gem_BX->at(iG) ) shouldKill2 = true;
      }
      if(!shouldKill)plotter.getOrMake1D("deadtime_gemBX", ";# of BX relative to ALCT data; # of GEM clusters", 20,-.5,19.5)->Fill(gem_BX->at(iG));
      if(!shouldKill2)plotter.getOrMake1D("deadtime2_gemBX", ";# of BX relative to ALCT data; # of GEM clusters", 20,-.5,19.5)->Fill(gem_BX->at(iG));

    }

//    for(unsigned int iG =0; iG < gem_chamber->size(); ++iG){
//      if(gem_chamber->at(iG) == 1){continue;}
//      int pad = toPad(gem_strip->at(iG));
//      int part = convP(gem_partition->at(iG));
//      int BX =gem_BX->at(iG);
//      for(unsigned int iG2 =iG+1; iG2 < gem_chamber->size(); ++iG2){
//        if(gem_chamber->at(iG2) == 1){continue;}
//        int pad2 = toPad(gem_strip->at(iG2));
//        int part2 = convP(gem_partition->at(iG2));
//        int BX2 =gem_BX->at(iG);
//        if(BX != BX2) continue;
//        if(part != part2) continue;
//        std::cout <<BX <<" "<< part <<" "<< pad <<" "<< pad2 <<std::endl;
//      }
//    }

    if(nClusters[11]){
      cout << endl << "<<<<<<<<<<<<<<<<<<";
      cout << lct_strip->at(0) <<" "<< lct_wg->at(0) <<" "<<endl;
      for(unsigned int iG = 0; iG < gem_chamber->size(); ++iG){
        if(gem_chamber->at(iG) == 0)
        cout <<"("<<gem_BX->at(iG)<<","<<gem_partition->at(iG)<<","<<gem_strip->at(iG)<<") ";
      }
    }

    plotter.getOrMake1D("gemChamberDiff", ";Number of strips per chamber", 39,-19.5,19.5)->Fill(n0-n1);

    for(unsigned int iB = 0; iB < nClusters.size(); ++iB)
      plotter.getOrMake2D("nBX_nClus",";# of BX relative to ALCT data; # of GEM clusters", 20,-.5,19.5, 20,-.5,19.5)->Fill(iB,nClusters[iB]);



    plotter.getOrMake1D("nBXwS", ";# of BX with gem data", 20,-.5,19.5)->Fill(nF);
    if(nF > 1)plotter.getOrMake2D("gemBXvBX", ";BX1 ;BX2",  20,-.5,19.5, 20,-.5,19.5)->Fill(b1,b2);





//    std::vector<int> gem_chamber  ;
//    std::vector<int> gem_partition;
//    std::vector<int> gem_BX       ;
//    std::vector<int> gem_strip    ;
//    std::vector<int> gem_nStrips  ;
//
//    std::vector<int> lct_chamber  ;
//    std::vector<int> lct_BX       ;
//    std::vector<int> lct_strip    ;
//    std::vector<int> lct_wg       ;
//    std::vector<int> lct_valid    ;

  }

  void write(TString fileName){ plotter.write(fileName);}
  HistGetter plotter;
};

#endif

void AnalyzeTMBTuples(std::string fileName,std::string outFileName = "plots.root"){
  Analyzer a(fileName,"Events");
  a.analyze();
  a.write(outFileName);
}
