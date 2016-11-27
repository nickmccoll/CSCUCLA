
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../include/BaseCSCAndGEMAnalyzer.h"
#include "include/GEMPlottingInfo.h"
#include "../include/GEMGeometry.h"
#include "../include/Segment.h"

#include<iostream>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"

using namespace std;
using namespace CSCGEMTuples;

int nPairs = 0;
int nMaxPairs =5000;
bool first = true;
double setZ = 50;
Segment * segments= new Segment[nMaxPairs];
ErrorPoint2D * clusters= new ErrorPoint2D[nMaxPairs];

double totalDist(const double *xx ){
  double dist = 0;
  for(unsigned int iP = 0; iP < nPairs; ++iP){
    if(first) cout << segments[iP].project(0,0,33,0).x() <<" "<<segments[iP].project(0,0,33,0).y()<<" "<< segments[iP].project(0,0,33,0).error_x() <<" "<<segments[iP].project(0,0,33,0).error_y() <<" -> ";
    auto p = segments[iP].project(xx[0],xx[1],xx[2],xx[3]);
    if(first) cout << p.x() <<" "<<p.y()<<" "<< p.error_x() <<" "<<p.error_y() <<" -> ";
    if(first) cout << clusters[iP].x() <<" "<<clusters[iP].y()<<" "<< clusters[iP].error_x() <<" "<<clusters[iP].error_y() <<" -> ";
    p -= clusters[iP];
    if(first) cout << p.x() <<" "<<p.y()<<" "<< p.error_x() <<" "<<p.error_y() <<" -> ";
    dist += p.x()*p.x() + p.y()*p.y() ;
    if(first) cout << dist <<endl;
  }
//  cout << xx[0]<<" "<<xx[1]<<" "<<xx[2]<<" "<< xx[3] << " " << dist<< endl;
  first = false;
  return dist;
}


double totalDist2(const double *xx ){
  double dist = 0;
  for(unsigned int iP = 0; iP < nPairs; ++iP){
    auto p = segments[iP].project2(xx[0],xx[1],xx[2],xx[3],xx[4],xx[5]);
    if(first) cout << p.x() <<" "<<p.y()<<" "<< p.error_x() <<" "<<p.error_y() <<" -> ";
    if(first) cout << clusters[iP].x() <<" "<<clusters[iP].y()<<" "<< clusters[iP].error_x() <<" "<<clusters[iP].error_y() <<" -> ";
    p -= clusters[iP];
    if(first) cout << p.x() <<" "<<p.y()<<" "<< p.error_x() <<" "<<p.error_y() <<" -> ";
//    dist += p.x()*p.x()/p.cov_xx() + p.y()*p.y()/p.cov_yy() ;
    dist += p.x()*p.x() + p.y()*p.y() ;
    if(first) cout <<  p.x()*p.x()/p.cov_xx() <<" "<<p.y()*p.y()/p.cov_yy()<<" "<< dist <<endl;
  }
//  cout << xx[0]<<" "<<xx[1]<<" "<<xx[2]<<" "<< xx[3] << " " << dist<< endl;
  first = false;
  return dist;
}




double totalDistWithoutZ(const double *xx ){
  double dist = 0;
  for(unsigned int iP = 0; iP < nPairs; ++iP){
    auto p = segments[iP].project(xx[0],xx[1],setZ,xx[2]);
    if(first) cout << p.x() <<" "<<p.y()<<" "<< p.error_x() <<" "<<p.error_y() <<" -> ";
    if(first) cout << clusters[iP].x() <<" "<<clusters[iP].y()<<" "<< clusters[iP].error_x() <<" "<<clusters[iP].error_y() <<" -> ";
    p -= clusters[iP];
    if(first) cout << p.x() <<" "<<p.y()<<" "<< p.error_x() <<" "<<p.error_y() <<" -> ";
//    dist += p.x()*p.x()/p.cov_xx() + p.y()*p.y()/p.cov_yy() ;
//    dist += p.x()*p.x()/p.cov_xx() + p.y()*p.y()/p.cov_yy() ;
    dist += (p.x()*p.x() + p.y()*p.y() ) *(p.x()*p.x() + p.y()*p.y() ) /( p.x()*p.x()*p.cov_xx()  + p.y()*p.y()*p.cov_yy() + 2*p.x()*p.y()*p.cov_xy() )  ;
    if(first) cout <<  p.x()*p.x()/p.cov_xx() <<" "<<p.y()*p.y()/p.cov_yy()<<" "<< dist <<endl;
  }
//  cout << xx[0]<<" "<<xx[1]<<" "<<xx[2]<<" "<< xx[3] << " " << dist<< endl;
  first = false;
  return dist;
}



class Analyze : public AnalyzeBoth {
public:
  Analyze(std::string cscFileName, std::string gemFileName,const GEMConfigInfo* info) : AnalyzeBoth(cscFileName,gemFileName,info)
  {
    cout << nPairs <<" "<< nMaxPairs <<endl;
  }
  virtual  ~Analyze() {};

  virtual void runAEvent() {
    if(nPairs >= nMaxPairs) return;
    if(!pureSample(csc)) return;
    if(gemClusters->size() != 1) return;


    auto p = cscSegments->at(0).project(0,0,33,0);
    if(TMath::Abs(p.x()) > 25) return;
    if(TMath::Abs(p.y()) > 65) return;


    //Clean up some bad strips
    if(gemClusters->at(0).iEta == 4  && gemClusters->at(0).getFirstStrip() <= 20) return;
    if(gemClusters->at(0).getFirstStrip() >= 65 && gemClusters->at(0).getFirstStrip() <= 80) return;
    if(gemClusters->at(0).getFirstStrip() >= 112 && gemClusters->at(0).getFirstStrip() <= 128) return;

    clusters[nPairs] = gemClusters->at(0).localCoords();
    segments[nPairs] = cscSegments->at(0);

    nPairs++;
  }

  void solve(TString outFileName){


    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "MIGRAD");
    min->SetMaxFunctionCalls(1000000);
    min->SetTolerance(0.001);
    min->SetPrintLevel(10);

        ROOT::Math::Functor f(&totalDist,4);
    double step[4] = {0.01,0.01,.01,.01};
    double variable[4] = { 0,0,50,0};
    min->SetFunction(f);
    min->SetVariable(0,"x"  ,variable[0], step[0]);
    min->SetVariable(1,"y"  ,variable[1], step[1]);
    min->SetVariable(2,"z"  ,variable[2], step[2]);
    min->SetVariable(3,"phi",variable[3], step[3]);

       // do the minimization
       min->Minimize();

       const double *xs = min->X();
       std::cout << "Minimum: f(" << xs[0] << "," << xs[1]<< "," << xs[2]<< "," << xs[3] << "): "
                 << min->MinValue()  << std::endl;


//           setZ=xs[2];
//
//           ROOT::Math::Minimizer* min2 = ROOT::Math::Factory::CreateMinimizer("Minuit2", "MIGRAD");
//           min2->SetMaxFunctionCalls(1000000);
//           min2->SetTolerance(0.001);
//           min2->SetPrintLevel(10);
//
//               ROOT::Math::Functor f2(&totalDistWithoutZ,3);
//           double step2[3] = {0.01,0.01,.01};
//           double variable2[3] = { xs[0],xs[1],xs[3]};
//           min2->SetFunction(f2);
//           min2->SetVariable(0,"x"  ,variable2[0], step2[0]);
//           min2->SetVariable(1,"y"  ,variable2[1], step2[1]);
//           min2->SetVariable(2,"phi",variable2[2], step2[2]);
//
//              // do the minimization
//              min2->Minimize();
//
//              const double *xs2 = min2->X();
//              std::cout << "Minimum: f(" << xs2[0] << "," << xs2[1]<< "," << xs2[2] << "): "
//                        << min2->MinValue()  << std::endl;

  }

};

#endif

void GEMAlignment(std::string cscfileName="csc_forsync.root",std::string gemfilename = "gem_forsync.root",std::string outFileName = "plotCSCAndGEM_out.root"){
  GEMConfigInfo info;
  info.vFATChanMapName       = "slot_table_904_june09.csv";
  Analyze a(cscfileName,gemfilename,&info);
  a.analyze();
  a.solve(outFileName);
//  a.draw();
}
