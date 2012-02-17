#include "MyAnalysis/DisplacedJetAnlzr/interface/genjet.h"
#include "MyAnalysis/DisplacedJetAnlzr/interface/pfjet.h"
#include "MyAnalysis/DisplacedJetAnlzr/interface/exotic.h"
#include "DataFormats/Common/interface/Wrapper.h"

namespace {

  struct dictionary {
    edm::Wrapper<exotic> exo;
    edm::Wrapper<std::vector<exotic> > exos;
    edm::Wrapper<genjet> genj;
    edm::Wrapper<std::vector<genjet> > genjs;
    edm::Wrapper<track> trk;
    edm::Wrapper<std::vector<track> > trks;
    edm::Wrapper<pfjet> pfj;
    edm::Wrapper<std::vector<pfjet> > pfjets;
  };

}
