#include "SimMuon/GEMDigitizer/interface/GEMDigiPreRecoModel.h"
#include "SimDataFormats/EncodedEventId/interface/EncodedEventId.h"

void 
GEMDigiPreRecoModel::fillDigis(int rollDetId, GEMDigiPreRecoCollection& digis)
{
  for (auto d: digi_)
  {
    digis.insertDigi(GEMDetId(rollDetId), d);
  }
  digi_.clear();
}

