#include "DisplacedDijet/DisplacedJetAnlzr/interface/DJ_DiJetVertices.h"
#include "DisplacedDijet/DisplacedJetAnlzr/interface/DisplacedDijet.h"

using namespace std;

DJ_DiJetVertices::DJ_DiJetVertices(const edm::ParameterSet& iConfig) : 
patJetCollectionTag_(iConfig.getParameter<edm::InputTag>("patJetCollectionTag")),
useTrackingParticles_(iConfig.getParameter<bool>("useTrackingParticles")),
PromptTrackDxyCut_(iConfig.getParameter<double>("PromptTrackDxyCut")),
TrackPtCut_(iConfig.getParameter<double>("TrackPtCut")),
TrackingEfficiencyFactor_(iConfig.getParameter<double>("TrackingEfficiencyFactor")),
vtxWeight_(iConfig.getParameter<double>("vtxWeight")),
PV_(iConfig.getParameter<unsigned int>("PV")),
vtxconfig_(iConfig.getParameter<edm::ParameterSet>("vertexfitter")),
vtxfitter_(vtxconfig_) {

  maxTrackToJetDeltaR_ = 0.4;

  if(iConfig.exists("maxTrackToJetDeltaR"))maxTrackToJetDeltaR_ = iConfig.getParameter<double>("maxTrackToJetDeltaR");


  associatorName_ = "TrackAssociatorByHits";
  if(iConfig.exists("associator"))associatorName_ = iConfig.getParameter<string>("associator");

  vertexCollectionTag_ = edm::InputTag("offlinePrimaryVertices");
  trackCollectionTag_ = edm::InputTag("generalTracks");
  truthTrackCollectionTag_ = edm::InputTag("mergedtruth","MergedTrackTruth","HLT");
  exoticMotherIDs_ = {6000111,6000112};
  if(iConfig.exists("motherIDs"))exoticMotherIDs_ = iConfig.getParameter<vector<int> >("motherIDs");

  beamspotLabel_ = edm::InputTag("offlineBeamSpot");
  if(iConfig.exists("beamSpotInputTag"))beamspotLabel_ = iConfig.getParameter<edm::InputTag>("beamSpotInputTag");
  beamspotToken_ = consumes<reco::BeamSpot>(beamspotLabel_);

  produces<std::vector<DisplacedDijet> >("DisplacedDijets");
  consumes<reco::TrackToTrackingParticleAssociator>(edm::InputTag(associatorName_));
  consumes<std::vector<reco::CaloJet> >(patJetCollectionTag_);
  //consumes<std::vector<pat::Jet> >(patJetCollectionTag_);
  vertexToken_ = consumes<edm::View<reco::Vertex> >(vertexCollectionTag_);
  //vertexToken_ = consumes<reco::VertexCollection>(vertexCollectionTag_);
  //consumes<reco::TrackCollection>(trackCollectionTag_);
  trackToken_ = consumes<edm::View<reco::Track> >(trackCollectionTag_);
  consumes<std::vector<TrackingParticle> >(truthTrackCollectionTag_);

}

void
DJ_DiJetVertices::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::auto_ptr<std::vector<DisplacedDijet> > output_DisplacedDijets(new std::vector<DisplacedDijet>);

  bool hasTracks = iEvent.getByToken(trackToken_,trackHandle_);

  if(!hasTracks){
   iEvent.put(output_DisplacedDijets,"DisplacedDijets");
   return;
  }

   GetEventInfo(iEvent,iSetup);

   vector<int> whichVertex(trackHandle_->size(),-1);
   for(int i = 0; i < (int)trackHandle_->size(); i++){
     double maxWeight = 0;
     int jj = -1;
     reco::TrackBaseRef tref(trackHandle_,i);
     for(int j = 0; j < (int)recVtxs->size();j++){
       if(recVtxs->at(j).trackWeight(tref) > maxWeight){
	 maxWeight = recVtxs->at(j).trackWeight(tref);
	 jj = j;
       }
     }
     whichVertex[i] = jj;
   }

   edm::ESHandle<MagneticField> magneticField;
   iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
   const MagneticField* magneticField_ = &*magneticField;

   edm::ESHandle<Propagator> thePropagator_;
   iSetup.get<TrackingComponentsRecord>().get("PropagatorWithMaterial",thePropagator_);
   StateOnTrackerBound stateOnTracker(thePropagator_.product());

   edm::Handle<std::vector<reco::CaloJet> > patJetsHandle;
   iEvent.getByLabel(patJetCollectionTag_,patJetsHandle);
   iEvent.getByToken(beamspotToken_, beamspotHandle_);

   //iEvent.getByLabel(trackCollectionTag_,trackHandle_);

   std::vector<std::vector<reco::TransientTrack> > associatedTrackVector;
   std::vector<std::vector<int> > whichVertexVector;

   for(int i = 0; i < int(patJetsHandle->size()); i++){
     const reco::CaloJet& calo_jet = patJetsHandle->at(i);
     TVector3 jetVec;
     jetVec.SetPtEtaPhi(calo_jet.detectorP4().Pt(),calo_jet.detectorP4().Eta(),calo_jet.detectorP4().phi());
     std::vector<reco::TransientTrack> associatedTracks;
     std::vector<int> vertexVector;
     for(int j = 0; j < (int)trackHandle_->size(); j++){
       reco::TrackBaseRef tref(trackHandle_,j);
       if(tref->pt() < TrackPtCut_)continue;
       if (!tref->quality(reco::TrackBase::highPurity)) continue;
       FreeTrajectoryState fts = trajectoryStateTransform::initialFreeState(trackHandle_->at(j),magneticField_);
       TrajectoryStateOnSurface outer = stateOnTracker(fts);
       if(!outer.isValid())continue;
       GlobalPoint outerPos = outer.globalPosition();
       TVector3 trackPos(outerPos.x(),outerPos.y(),outerPos.z());
       if(trackPos.DeltaR(jetVec) > maxTrackToJetDeltaR_)continue;
       reco::TransientTrack tt(trackHandle_->at(j),magneticField_);
       associatedTracks.push_back(tt);
       vertexVector.push_back(whichVertex[j]);
     }
     associatedTrackVector.push_back(associatedTracks);
     whichVertexVector.push_back(vertexVector);
   }

   for (int i=0;i<int(patJetsHandle->size()-1);i++){
    for (size_t k=i+1;k<patJetsHandle->size();k++){

     reco::CaloJet jet1 = patJetsHandle->at(i);
     reco::CaloJet jet2 = patJetsHandle->at(k);
     reco::Candidate::LorentzVector p4 = jet1.p4() + jet2.p4();

     //DEBUG
     //std::cout << "jet1: pt,eta,phi" << jet1.pt() << " " <<jet1.eta() << " " <<jet1.phi() << std::endl; 
     //std::cout << "jet2: pt,eta,phi" << jet2.pt() << " " <<jet2.eta() << " " <<jet2.phi() << std::endl; 
     //ENDDEBUG

     GlobalVector direction(p4.px(), p4.py(), p4.pz());
     direction = direction.unit();

     //tracks selection
     //reco::TrackRefVector dijettrks = jet1.associatedTracks();
     //reco::TrackRefVector dijettrks2 = jet2.associatedTracks();
     std::vector<reco::TransientTrack> dijettrks = associatedTrackVector[i];
     std::vector<reco::TransientTrack> dijettrks2 = associatedTrackVector[k];

     ////TO REMOVE
     std::vector<int> indices1(dijettrks.size());
     std::vector<int> indices2(dijettrks2.size());
     std::fill(indices1.begin(),indices1.end(),1);
     std::fill(indices2.begin(),indices2.end(),2);

     for(size_t j=0;j<dijettrks2.size();j++) dijettrks.push_back(dijettrks2[j]);
     //indices for jet1 and jet 2
     std::vector<int> indices(dijettrks.size());
     std::fill(indices.begin(),indices.end()-dijettrks2.size(),1);
     std::fill(indices.end()-dijettrks2.size(),indices.end(),2);

     std::vector<int> vertexVector = whichVertexVector[i];
     std::vector<int> vertexVector2 = whichVertexVector[k];
     for(size_t j = 0; j < vertexVector2.size(); j++)vertexVector.push_back(vertexVector2[j]);

     std::vector<reco::TransientTrack> trksToVertex;     
     std::vector<int> indicesToVertex;
     std::vector<float> glxysToVertex;
     std::vector<float> ip2dsToVertex;
     
     int nPromptTracks=0;
     int nPromptTracks1=0;
     int nPromptTracks2=0;
     float PromptEnergy=0;
     float PromptEnergy1=0;
     float PromptEnergy2=0;
     float trkAvgPt=0;

     std::vector<reco::TransientTrack> matchedTracks;
     std::vector<int> matchedVertexVector;

     for (size_t j=0;j<dijettrks.size();j++){
        
       const reco::TransientTrack trk = dijettrks[j];

       reco::TransientTrack t_trk = trk;
       matchedTracks.push_back(t_trk);
       matchedVertexVector.push_back(vertexVector[j]);


       //reco::TransientTrack t_trk = theB->build(trk);
       Measurement1D ip2d = IPTools::signedTransverseImpactParameter(t_trk,direction,pv).second;
       Measurement1D ip3d = IPTools::signedImpactParameter3D(t_trk,direction,pv).second;
       if (fabs(ip2d.value())<PromptTrackDxyCut_){ 
	 PromptEnergy += sqrt(0.1396*0.1396 + trk.track().p()*trk.track().p());
	 nPromptTracks+=1;
	 if(indices[j] == 1){
	   nPromptTracks1 += 1;
	   PromptEnergy1 += sqrt(0.1396*0.1396 + trk.track().p()*trk.track().p());
	 } else if(indices[j] == 2){
	   nPromptTracks2 += 1;
	   PromptEnergy2 += sqrt(0.1396*0.1396 + trk.track().p()*trk.track().p());
	 }
	 continue;
       }
	// tracking inefficiency factor
	if (rand()/float(RAND_MAX) > TrackingEfficiencyFactor_) continue;

	trkAvgPt+=trk.track().pt();
        float r = 100*3.3*trk.track().pt()/3.8;
        float guesslxy = ip2d.value()/sin(trk.track().phi()-direction.phi())*(1-2.5*fabs(ip2d.value())/r);

        ip2dsToVertex.push_back(ip2d.value());
        glxysToVertex.push_back(fabs(guesslxy));
        trksToVertex.push_back(t_trk);
        indicesToVertex.push_back(indices[j]);

     }

     if(trksToVertex.size() < 2)continue;
     vector<float> tracksIPLogSig;
     vector<float> tracksIPLog10Sig;

     DisplacedDijet dijetHolder;
     dijetHolder.TrkAvgPt = trksToVertex.size() > 0 ? trkAvgPt/trksToVertex.size() : -1;
     dijetHolder.NPromptTracks = nPromptTracks;
     dijetHolder.PromptEnergyFrac = PromptEnergy/(jet1.energy()+jet2.energy());
     dijetHolder.NDispTracks = trksToVertex.size();
     dijetHolder.NPromptTracks1 = nPromptTracks1;
     dijetHolder.NPromptTracks2 = nPromptTracks2;
     dijetHolder.PromptEnergyFrac1 = PromptEnergy1/jet1.energy();
     dijetHolder.PromptEnergyFrac2 = PromptEnergy2/jet2.energy();


     bool goodVtx = false;
     TransientVertex jvtx = vtxfitter_.vertex(trksToVertex);
     reco::Vertex vtx(jvtx);
     if (jvtx.isValid() && 
          !vtx.isFake() &&
          (vtx.nTracks(vtxWeight_)>1) &&
          (vtx.normalizedChi2()>0) &&
          (vtx.normalizedChi2()<10)) goodVtx = true;

     ROOT::Math::SVector<double,3> vector(vtx.position().x() - pv.x(),vtx.position().y()-pv.y(),0);
     float lxy = ROOT::Math::Mag(vector);
     reco::Candidate::CovarianceMatrix matrix = vtx.covariance() + pv.covariance();
     float err = sqrt(ROOT::Math::Similarity(matrix,vector))/lxy;
     float sig = lxy/err;
     
     int n=vtx.nTracks(vtxWeight_);
     int charge=0;
     int nposip2d=0;
     int hitsInFrontOfVert=0;
     int missHitsAfterVert=0;
     int FromExo=0;
     int VtxN1=0;
     int VtxN2=0;
     std::vector<float> glxysVertex;
     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vtxP4;
     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vtxP4_1;
     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vtxP4_2;

     //DEBUG
     //std::cout << "vtx Lxy: " << lxy << std::endl; 
     //ENDDEBUG
     for (size_t j=0;j<trksToVertex.size();j++){

       reco::TransientTrack t_trk = trksToVertex[j];

       if (jvtx.trackWeight(t_trk)>vtxWeight_){
         GlobalVector p3 = t_trk.trajectoryStateClosestToPoint(jvtx.position()).momentum();
 	 vtxP4 += ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double> >(p3.x(),p3.y(),p3.z(),0.13957018);
         charge+=t_trk.track().charge();
	 tracksIPLogSig.push_back(log(t_trk.stateAtBeamLine().transverseImpactParameter().significance()));
	 tracksIPLog10Sig.push_back(log10(t_trk.stateAtBeamLine().transverseImpactParameter().significance()));
	 // DEBUG
         //std::cout << "vtx_trk: pt,eta,phi,N,weight " << t_trk.track().pt() << " " << t_trk.track().eta() << " " << t_trk.track().phi() << " " << indicesToVertex[j]
         //<< " " << jvtx.trackWeight(t_trk) << std::endl;
	 // END DEBUG
         if (indicesToVertex[j] == 1){
	   VtxN1+=1;
	   vtxP4_1 += ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double> >(p3.x(),p3.y(),p3.z(),0.13957018);

	 }
         if (indicesToVertex[j] == 2){
	   VtxN2+=1;
	   vtxP4_2 += ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double> >(p3.x(),p3.y(),p3.z(),0.13957018);

	 }
         if (ip2dsToVertex[j]>0) nposip2d+=1;
         // hitPattern
         CheckHitPattern::Result res = checkHitPattern_.analyze(iSetup,t_trk.track(),jvtx.vertexState(),false);
         hitsInFrontOfVert += res.hitsInFrontOfVert;
         missHitsAfterVert += res.missHitsAfterVert;
         //glxys
         glxysVertex.push_back(glxysToVertex[j]);
         // tracks from Exotic
         edm::RefToBase<reco::Track> ref_trk = t_trk.trackBaseRef();
         if(RecoToSimColl.find(ref_trk) != RecoToSimColl.end()){
           TrackingParticleRef tp = RecoToSimColl[ref_trk].begin()->first;

           if (tp->genParticles().size()>0){
             std::vector<std::pair<int,double> > moms;
             const reco::GenParticle *gp = tp->genParticles().at(0).get();
             moms.push_back(std::pair<int,double> (gp->pdgId(),gp->vertex().rho()));
             GetMothers(tp.get(),moms);
	     if(find(exoticMotherIDs_.begin(),exoticMotherIDs_.end(),moms.back().first) != exoticMotherIDs_.end())FromExo += 1;
	     //if (moms.back().first == 6000111 || moms.back().first==6000112) FromExo+=1;
           } // genParticle found
         } // TrackAssociation  
       } // vtx tracks
     } // tracks to Vertex

     float dR = deltaR(direction.eta(),direction.phi(),vtxP4.eta(),vtxP4.phi());

     reco::Candidate::LorentzVector physicsP4;
     //reco::Candidate::LorentzVector physicsP41,physicsP42,physicsP43;
     if (goodVtx){
       // corrected P4
       physicsP4 = jet1.physicsP4(vtx.position())
	 +jet2.physicsP4(vtx.position());

       //physicsP41 = detectorP4(jet1,vtx,jvtx,0) + detectorP4(jet2,vtx,jvtx,0);
       //physicsP42 = detectorP4(jet1,vtx,jvtx,1) + detectorP4(jet2,vtx,jvtx,1);
       //physicsP43 = detectorP4(jet1,vtx,jvtx,2) + detectorP4(jet2,vtx,jvtx,2);
     }

     double deltaPhi = 0;
     double medianDeltaPhi = 0;
     double alphaMax = calculateAlphaMax(matchedTracks,matchedVertexVector);
     if(goodVtx)deltaVertex2D(jvtx.position(),matchedTracks,deltaPhi,medianDeltaPhi);

     dijetHolder.CorrPt = (goodVtx ? physicsP4.Pt() : -1);
     dijetHolder.CorrEta = (goodVtx ? physicsP4.Eta() : -999);
     dijetHolder.CorrPhi = (goodVtx ? physicsP4.Phi() : -999);
     dijetHolder.CorrMass = (goodVtx ? physicsP4.M() : -1);

     dijetHolder.Lxy = goodVtx ? lxy : -1;
     dijetHolder.Lxysig = goodVtx ? sig : -1;
     dijetHolder.VtxX = goodVtx ? (vtx.position().x()-pv.x()) : 999;
     dijetHolder.VtxY = goodVtx ? (vtx.position().y()-pv.y()) : 999;
     dijetHolder.VtxZ = goodVtx ? (vtx.position().z()-pv.z()) : 999;
     dijetHolder.VtxChi2 = goodVtx ? vtx.normalizedChi2() : -1;
     dijetHolder.Vtxmass = goodVtx ? vtxP4.M() : -1;
     dijetHolder.Vtxpt = goodVtx ? vtxP4.Pt() : -1;
     dijetHolder.VtxdR = goodVtx ? dR : -1;
     dijetHolder.VtxN = goodVtx ? n : -1;
     dijetHolder.VtxCharge = goodVtx ? charge : 999;
     dijetHolder.NAvgHitsInFrontOfVert = goodVtx ? hitsInFrontOfVert/float(n) : -1;
     dijetHolder.NAvgMissHitsAfterVert = goodVtx ? missHitsAfterVert/float(n) : -1;
     dijetHolder.ExoVtxFrac = goodVtx ? FromExo/float(n) : -1;
     dijetHolder.Posip2dFrac = goodVtx ? nposip2d/float(n) : -1;
     dijetHolder.VtxN1 = goodVtx ? VtxN1 : -1;
     dijetHolder.VtxN2 = goodVtx ? VtxN2 : -1;

     dijetHolder.trackMass1 = goodVtx ? vtxP4_1.M() : -1;
     dijetHolder.trackEnergy1 = goodVtx ? vtxP4_1.E() : -1;
     dijetHolder.trackMass2 = goodVtx ? vtxP4_2.M() : -1;
     dijetHolder.trackEnergy2 = goodVtx ? vtxP4_2.E() : -1;

     dijetHolder.openingAngle = goodVtx ? ROOT::Math::VectorUtil::Angle(vtxP4_1,vtxP4_2) : -999;

     sort(tracksIPLogSig.begin(), tracksIPLogSig.end());
     sort(tracksIPLog10Sig.begin(), tracksIPLog10Sig.end());

     if(tracksIPLog10Sig.size() == 0){
       dijetHolder.medianIPLog10Sig = -999;
     }else if((tracksIPLog10Sig.size()%2 == 0)){
       dijetHolder.medianIPLog10Sig = (tracksIPLog10Sig.at(tracksIPLog10Sig.size()/2-1)+tracksIPLog10Sig.at((tracksIPLog10Sig.size()/2)))/2;
     }else{
       dijetHolder.medianIPLog10Sig = tracksIPLog10Sig.at((tracksIPLog10Sig.size()-1)/2);
     }
     dijetHolder.alphaMax = alphaMax;     
     dijetHolder.totalTrackAngle = goodVtx ? deltaPhi : -999;
     dijetHolder.medianTrackAngle = goodVtx ? medianDeltaPhi : -999;
     dijetHolder.medianTrackAngleLog10 = goodVtx ? log10(medianDeltaPhi) : -999;
     // do glxy stuff here
     helpers help;

     dijetHolder.glxydistall = goodVtx ? help.AvgDistance(glxysToVertex,lxy) : -1;
     dijetHolder.glxydistvtx = goodVtx ? help.AvgDistance(glxysVertex,lxy) : -1;
     dijetHolder.glxyrmsall = goodVtx ? help.RMS(glxysToVertex,lxy) : -1;
     dijetHolder.glxyrmsvtx = goodVtx ? help.RMS(glxysVertex,lxy) : -1;

     // clusters
     std::vector<std::vector<float> > clusters = help.clusters(glxysToVertex,0.15*lxy);
     std::vector<float> bestcluster = help.bestcluster(clusters,lxy);

     int clrN1=0;
     int clrN2=0;
     for (size_t l=0;l<bestcluster.size();l++){
       float glxy = bestcluster[l];
       for (size_t m=0;m<glxysToVertex.size();m++)
         if(fabs(glxy-glxysToVertex[m])<1e-5){
           if (indicesToVertex[m]==1) clrN1++;
           else clrN2++;
           break;
       }
     }

     dijetHolder.glxyrmsclr = goodVtx ? help.RMS(bestcluster,lxy) : -1;
     dijetHolder.glxydistclr = goodVtx ? help.AvgDistance(bestcluster,lxy) : -1;
     dijetHolder.Nclusters = goodVtx ? help.Nclusters(clusters) : -1;
     dijetHolder.bestclusterlxy = goodVtx ? ( (bestcluster.size()>0 ) ? help.Avg(bestcluster)/lxy : -1) : -1;
     dijetHolder.bestclusterN = goodVtx ? bestcluster.size() : -1;
     dijetHolder.bestclusterN1 = goodVtx ? clrN1 : -1;
     dijetHolder.bestclusterN2 = goodVtx ? clrN2 : -1;

     output_DisplacedDijets->push_back(dijetHolder);
   } 
  }// dijet loop

   iEvent.put(output_DisplacedDijets,"DisplacedDijets");

}


reco::Candidate::LorentzVector DJ_DiJetVertices::detectorP4(pat::Jet &jet, reco::Vertex &vtx, TransientVertex &tvtx, int CorrectTracks=0){

  reco::Candidate::LorentzVector P4;
  // correct track momentas with pca to the vertex
  const std::vector<reco::PFCandidatePtr> & parts = jet.getPFConstituents();
  for(std::vector<reco::PFCandidatePtr>::const_iterator part = parts.begin(); part!= parts.end(); ++part){
    const reco::PFCandidate *pfCand = part->get();
    if (pfCand->charge() == 0){
      P4 += reco::Jet::physicsP4(vtx.position(),*pfCand,pfCand->vertex());
    } else {
      if (CorrectTracks == 0) P4 += pfCand->p4();
      else if (CorrectTracks == 1) {
        P4 += reco::Jet::physicsP4(vtx.position(),*pfCand,pfCand->vertex());
      }
      else if (CorrectTracks == 2){
        if (pfCand->trackRef().isNull()) continue;
        reco::TransientTrack t_trk = theB->build(pfCand->trackRef());
        //const reco::Track *trk = pfCand->trackRef().get();
        GlobalVector p3 = t_trk.trajectoryStateClosestToPoint(tvtx.position()).momentum();
        P4 += ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double> >(p3.x(),p3.y(),p3.z(),pfCand->mass());
      }
    }
  }

  return P4;
}

void DJ_DiJetVertices::GetEventInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup){

// Vertices and tracks

  //edm::Handle<reco::VertexCollection> recVtxs;
   iEvent.getByToken(vertexToken_, recVtxs);
   if (recVtxs->size()>PV_) 
     pv=recVtxs->at(PV_); 
   else
     pv = recVtxs->front();

// TransientTrack Builder
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

// Reco to Sim Track association
   if(!iEvent.isRealData() && useTrackingParticles_){
     edm::Handle<std::vector<TrackingParticle> > TPCollectionH ;
     try{
       if(iEvent.getByLabel(truthTrackCollectionTag_,TPCollectionH)){
	 const reco::TrackToTrackingParticleAssociator* m_associator;
	 edm::Handle<reco::TrackToTrackingParticleAssociator> assocHandle;
	 if(iEvent.getByLabel(associatorName_,assocHandle)){
	   m_associator = assocHandle.product();
	   RecoToSimColl = m_associator->associateRecoToSim(trackHandle_,TPCollectionH);
	 }
       }
     } catch (...) {;}
   }

}

void DJ_DiJetVertices::GetMothers(const TrackingParticle* gp, std::vector<std::pair<int,double> > &moms){

  const TrackingVertex* gv = gp->parentVertex().get();
  if(gv != 0 ){
    for(TrackingVertex::tp_iterator mom = gv->daughterTracks_begin(); mom != gv->daughterTracks_end(); mom++){
      moms.push_back(std::pair<int,double> ( (*mom).get()->genParticles().at(0).get()->pdgId(), gv->position().rho() ));
      if (find(exoticMotherIDs_.begin(),exoticMotherIDs_.end(),moms.back().first) != exoticMotherIDs_.end())
	return;
      GetMothers((*mom).get(),moms);
      break;
    }
  }      
  return ;
}

void DJ_DiJetVertices::deltaVertex2D(GlobalPoint secVert, std::vector<reco::TransientTrack> tracks, double& dPhi, double& medianDPhi)
{
  const reco::BeamSpot& pat_beamspot = (*beamspotHandle_);
  TVector2 bmspot(pat_beamspot.x0(),pat_beamspot.y0());
  TVector2 sv(secVert.x(),secVert.y());
  TVector2 diff = (sv-bmspot);
  TVector2 trackPt(0,0);
  vector<double>angles;
  for(int i = 0; i < (int)tracks.size(); i++){
    TVector2 tt;
    tt.SetMagPhi(tracks[i].trajectoryStateClosestToPoint(secVert).momentum().transverse(),tracks[i].trajectoryStateClosestToPoint(secVert).momentum().phi());
    angles.push_back(fabs(diff.DeltaPhi(tt)));
    trackPt += tt;
  }
  sort(angles.begin(),angles.end());
  dPhi = diff.DeltaPhi(trackPt);
  if(angles.size()==0){
    medianDPhi = 1e-10;
  }else if(angles.size() %2 == 0){
    medianDPhi = angles.at(angles.size()/2-1);
  }else{
    medianDPhi = angles.at((angles.size()-1)/2);
  }
  //pt = (trackPt - ((trackPt * diff)/(diff * diff) * diff)).Mod();
}

double DJ_DiJetVertices::calculateAlphaMax(vector<reco::TransientTrack>tracks,vector<int> whichVertex)
{
  double total = 0;
  vector<double> alphas(recVtxs->size(),0);
  for(int i = 0; i < (int)tracks.size(); i++){
    double pt = tracks[i].initialFreeState().momentum().transverse();
    total += pt;
    if(whichVertex[i] < 0)continue;
    alphas[whichVertex[i]] += pt;
  }
  double alphaMax = 0;
  for(int i = 0; i < (int)alphas.size(); i++){
    if(alphas[i] > alphaMax)alphaMax = alphas[i];
  }
  return alphaMax / total;
}
