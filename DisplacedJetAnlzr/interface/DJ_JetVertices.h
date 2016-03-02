#ifndef DJ_JETVERTICEsS
#define DJ_JETVERTICESS

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//Data Formats and Tools 
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "RecoVertex/ConfigurableVertexReco/interface/ConfigurableVertexFitter.h"
#include "PhysicsTools/RecoUtils/interface/CheckHitPattern.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "DisplacedDijet/DisplacedJetAnlzr/interface/helpers.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"

class DJ_JetVertices : public edm::EDProducer {
   public:
      explicit DJ_JetVertices(const edm::ParameterSet&);

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&);

      void GetMothers(const TrackingParticle *p, std::vector<std::pair<int,double> > &moms);
      void GetEventInfo(const edm::Event&, const edm::EventSetup&);


      // configurables
      const edm::InputTag patJetCollectionTag_;
      bool useTrackingParticles_;
      double PromptTrackDxyCut_,TrackPtCut_,TrackingEfficiencyFactor_,vtxWeight_;
      unsigned int PV_;
      const edm::ParameterSet vtxconfig_;
      ConfigurableVertexFitter vtxfitter_;

      // stuff to use over the different functions
      reco::RecoToSimCollection RecoToSimColl;

      // global edm objects
      reco::Vertex pv;
      edm::ESHandle<TransientTrackBuilder> theB;
      CheckHitPattern checkHitPattern_;

      std::string associatorName_;

};

#endif
