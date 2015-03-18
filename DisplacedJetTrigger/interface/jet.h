#ifndef DisplacedDijet_DisplacedJetTrigger_jet_h
#define DisplacedDijet_DisplacedJetTrigger_jet_h

#include <vector>
#include "FWCore/Utilities/interface/typedefs.h"
#include "DisplacedDijet/DisplacedJetTrigger/interface/track.h"

struct jet {

   double energy,pt,eta,phi;
   std::vector<track> tracks;

};

#endif
