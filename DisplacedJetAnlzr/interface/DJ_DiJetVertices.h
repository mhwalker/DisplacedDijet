#ifndef DJ_DIJETVERTICEsS
#define DJ_DIJETVERTICESS

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//Data Formats and Tools 
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "RecoVertex/ConfigurableVertexReco/interface/ConfigurableVertexFitter.h"
#include "PhysicsTools/RecoUtils/interface/CheckHitPattern.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/GeomPropagators/interface/StateOnTrackerBound.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "DisplacedDijet/DisplacedJetAnlzr/interface/helpers.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"

#include <TVector3.h>

#include <algorithm>
#include <math.h>

class DJ_DiJetVertices : public edm::EDProducer {
   public:
      explicit DJ_DiJetVertices(const edm::ParameterSet&);

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&);

      void GetMothers(const TrackingParticle* p, std::vector<std::pair<int,double> > &moms);
      void GetEventInfo(const edm::Event&, const edm::EventSetup&);
      reco::Candidate::LorentzVector detectorP4(pat::Jet &jet, reco::Vertex &vtx, TransientVertex &tvtx, int CorrectTracks);

      void deltaVertex2D(GlobalPoint secVert, std::vector<reco::TransientTrack> tracks, double& dPhi, double& medianDPhi);
      double calculateAlphaMax(std::vector<reco::TransientTrack>tracks,std::vector<int> whichVertex);

      // configurables
      const edm::InputTag patJetCollectionTag_;
      bool useTrackingParticles_;
      double PromptTrackDxyCut_,TrackPtCut_,TrackingEfficiencyFactor_,vtxWeight_,maxTrackToJetDeltaR_;
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
      edm::InputTag vertexCollectionTag_;
      edm::InputTag trackCollectionTag_;
      edm::InputTag truthTrackCollectionTag_;

      edm::Handle<reco::BeamSpot> beamspotHandle_;
      edm::Handle<edm::View<reco::Track> > trackHandle_;
      edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
      edm::EDGetTokenT<edm::View<reco::Track> >trackToken_;
      edm::EDGetTokenT<edm::View<reco::Vertex> > vertexToken_;
      edm::Handle<edm::View<reco::Vertex> > recVtxs;
      edm::InputTag beamspotLabel_;

      std::vector<int> exoticMotherIDs_;
};

#endif
