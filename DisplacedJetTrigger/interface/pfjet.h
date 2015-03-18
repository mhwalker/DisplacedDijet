#ifndef DisplacedDijet_DisplacedJetTrigger_pfjet_h
#define DisplacedDijet_DisplacedJetTrigger_pfjet_h

#include <vector>
#include "FWCore/Utilities/interface/typedefs.h"
#include "DisplacedDijet/DisplacedJetTrigger/interface/track.h"

struct pfjet {

   double energy,pt,eta,phi;
   double phFrac,neuHadFrac,chgHadFrac,eleFrac,muFrac;
   int phN,neuHadN,chgHadN,eleN,muN;
   std::vector<track> tracks;

};

#endif
