
#include "../include/GEMInfo.h"
#include "../../gem-light-dqm/gemtreewriter/include/Event.h"

namespace CSCGEMTuples {
VFAT::VFAT(int iEta, int iPhi, std::vector<int>& onStrips) : iEta(iEta), iPhi(iPhi),  strips(onStrips){
  std::sort(strips.begin(),strips.end());
}

GEMInfo::GEMInfo(const GEMConfigInfo& gemInfo) :geo(gemInfo.gemXScale), channelMapping(gemInfo,geo.getNPartitions(), geo.getNPhis()) {}
void GEMInfo::build(Event * event) {
  vFats.clear();
  std::vector<AMC13Event> v_amc13 = event->amc13s();
  std::vector<AMCdata> v_amc = v_amc13[0].amcs();
  BX = v_amc[0].BX();
  evtN = event->GetEventNumber();
  std::vector<GEBdata> v_geb;
      v_geb = v_amc[0].gebs();
      for (unsigned int j = 0; j < v_geb.size(); j++) {
        // create vector of VFATdata. For data format details look at Event.h
        std::vector<VFATdata> v_vfat;
        v_vfat = v_geb.at(j).vfats();

        // loop over vfats
        for (unsigned int k = 0; k < v_vfat.size(); k++) {
          auto * m_vfat = &v_vfat[k];
          std::vector<int> firedChannels;
          uint16_t chan0xfFiredchip = 0;
          for (int chan = 0; chan < 128; ++chan) {
            if (chan < 64) {
              chan0xfFiredchip = ((m_vfat->lsData() >> chan) & 0x1);
              if (chan0xfFiredchip) {
                firedChannels.push_back(chan);
              }
            } else {
              chan0xfFiredchip = ((m_vfat->msData() >> (chan - 64)) & 0x1);
              if (chan0xfFiredchip) {
                firedChannels.push_back(chan);
              }
            }
          }

          int sn_ = channelMapping.getGEBSlotIndex(v_vfat[k].ChipID());
//          std::cout << v_vfat[k].ChipID()<< sn_ << std::endl;
          int iEta,iPhi;
          geo.getVFATIEtaIPhi(sn_,iEta,iPhi);
          std::vector<int> strips;
          if (firedChannels.size()) {
//
            for(unsigned int iC = 0; iC < firedChannels.size(); ++iC){
              int stripNumberOnVF = channelMapping.getStripNumber(sn_,firedChannels[iC]);
              strips.push_back(stripNumberOnVF);
            }
            VFAT fat(iEta,iPhi,strips);
            vFats.push_back(fat);
          }
        }
      }
      getClusters();

}

void GEMInfo::getClusters() {
clusters.clear();
std::vector<std::vector<int> > stripsPerPartition (geo.getNPartitions(), std::vector<int>(0));
//  std::cout << std::endl << "Start event!"<<std::endl;
for(const auto& vf : vFats){
//    std:: cout <<vf.idx <<" "<< vf.nRow <<" ";
  for(unsigned int iS = 0; iS < vf.strips.size(); ++iS){
    int globalStripNumber = geo.getPartitionStripNumber(vf.iPhi,vf.strips[iS]);
    stripsPerPartition[vf.iEta -1].push_back(globalStripNumber); //iEtas go from 1-8
//      std::cout <<"("<<vf.strips[iS]<<","<<globalStripNumber<<") ";
  }
//    std::cout << std::endl;
}

for(unsigned int iE = 0; iE < stripsPerPartition.size(); ++iE){
  std::sort(stripsPerPartition[iE].begin(),stripsPerPartition[iE].end(), std::less<int>());
//    std::cout << iR <<" -> ";
  for(unsigned int iC = 0; iC < stripsPerPartition[iE].size(); ++iC){
//      std::cout << " " << stripsPerRow[iR][iC];
    if(iC == 0 || stripsPerPartition[iE][iC] - stripsPerPartition[iE][iC -1]  > 1 ){
      clusters.emplace_back(iE + 1,stripsPerPartition[iE][iC],1 );
    } else{
      clusters.back().nStrips++;
    }
  }
//    std::cout << std::endl;
}

//  for(unsigned int iC = 0; iC < clusters.size(); ++iC){
//    std::cout << "(" << clusters[iC].nRow <<","<<clusters[iC].getNStrips()<<","<<clusters[iC].getFirstStrip()<<","<<clusters[iC].getLastStrip() <<") ";
//  }
//  std::cout << std::endl;

//now set the local coordinates;
for(auto& c : clusters){ c.localCoords() = geo.getClusterInfo(c.iEta,c.getFirstStrip(),c.getNStrips());}

}

}
