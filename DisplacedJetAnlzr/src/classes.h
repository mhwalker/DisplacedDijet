#ifndef DJ_CLASSES_H
#define DJ_CLASSES_H

#include "DisplacedDijet/DisplacedJetAnlzr/interface/DisplacedDijet.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include <vector>
#include <string>
#include <map>

namespace {

  struct dictionary {
    std::map<std::string,bool> dummy0;
    edm::Wrapper<std::map<std::string,bool> > dummy1;

    std::map<std::string,int> dummi0;
    edm::Wrapper<std::map<std::string,int> > dummi1;

    std::map<std::string,std::string> dummee0;
    edm::Wrapper<std::map<std::string,std::string> > dummee1;

    edm::Wrapper<DisplacedDijet> dd;
    edm::Wrapper<std::vector<DisplacedDijet> >dds;

  };

}

#endif
