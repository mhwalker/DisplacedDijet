#include "DisplacedDijet/DisplacedJetAnlzr/interface/DJ_DiJetVertices.h"
#include "DisplacedDijet/DisplacedJetAnlzr/interface/DisplacedDijet.h"

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

   produces<std::vector<DisplacedDijet> >("DisplacedDijets");
}

void
DJ_DiJetVertices::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::auto_ptr<std::vector<DisplacedDijet> > output_DisplacedDijets(new std::vector<DisplacedDijet>);

   GetEventInfo(iEvent,iSetup);

   edm::Handle<std::vector<pat::Jet> > patJetsHandle;
   iEvent.getByLabel(patJetCollectionTag_,patJetsHandle);

   for (int i=0;i<int(patJetsHandle->size()-1);i++){
    for (size_t k=i+1;k<patJetsHandle->size();k++){
     pat::Jet jet1 = patJetsHandle->at(i);
     pat::Jet jet2 = patJetsHandle->at(k);
     reco::Candidate::LorentzVector p4 = jet1.p4() + jet2.p4();

     //DEBUG
     //std::cout << "jet1: pt,eta,phi" << jet1.pt() << " " <<jet1.eta() << " " <<jet1.phi() << std::endl; 
     //std::cout << "jet2: pt,eta,phi" << jet2.pt() << " " <<jet2.eta() << " " <<jet2.phi() << std::endl; 
     //ENDDEBUG

     GlobalVector direction(p4.px(), p4.py(), p4.pz());
     direction = direction.unit();

     //tracks selection
     reco::TrackRefVector dijettrks = jet1.associatedTracks();
     reco::TrackRefVector dijettrks2 = jet2.associatedTracks();

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

     std::vector<reco::TransientTrack> trksToVertex;     
     std::vector<int> indicesToVertex;
     std::vector<float> glxysToVertex;
     std::vector<float> ip2dsToVertex;
     
     int nPromptTracks=0;
     float PromptEnergy=0;
     float trkAvgPt=0;
     for (size_t j=0;j<dijettrks.size();j++){
        
       const reco::TrackRef trk = dijettrks[j];

       if (!trk->quality(reco::TrackBase::highPurity)) continue;
       if (trk->pt() < TrackPtCut_) continue;

       reco::TransientTrack t_trk = theB->build(trk);
       Measurement1D ip2d = IPTools::signedTransverseImpactParameter(t_trk,direction,pv).second;
       Measurement1D ip3d = IPTools::signedImpactParameter3D(t_trk,direction,pv).second;
        if (fabs(ip3d.value())<0.03) nPromptTracks+=1;
        if (fabs(ip2d.value())<PromptTrackDxyCut_){ 
          PromptEnergy += sqrt(0.1396*0.1396 + trk->p()*trk->p());
          continue;
        }
	// tracking inefficiency factor
	if (rand()/float(RAND_MAX) > TrackingEfficiencyFactor_) continue;

	trkAvgPt+=trk->pt();
        float r = 100*3.3*trk->pt()/3.8;
        float guesslxy = ip2d.value()/sin(trk->phi()-direction.phi())*(1-2.5*fabs(ip2d.value())/r);

        ip2dsToVertex.push_back(ip2d.value());
        glxysToVertex.push_back(fabs(guesslxy));
        trksToVertex.push_back(t_trk);
        indicesToVertex.push_back(indices[j]);

     }

     if(trksToVertex.size() < 2)continue;

     DisplacedDijet dijetHolder;
     dijetHolder.TrkAvgPt = trksToVertex.size() > 0 ? trkAvgPt/trksToVertex.size() : -1;
     dijetHolder.NPromptTracks = nPromptTracks;
     dijetHolder.PromptEnergyFrac = PromptEnergy/(jet1.energy()+jet2.energy());
     dijetHolder.NDispTracks = trksToVertex.size();

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

     //DEBUG
     //std::cout << "vtx Lxy: " << lxy << std::endl; 
     //ENDDEBUG
     for (size_t j=0;j<trksToVertex.size();j++){
       reco::TransientTrack t_trk = trksToVertex[j];
       if (jvtx.trackWeight(t_trk)>vtxWeight_){
         GlobalVector p3 = t_trk.trajectoryStateClosestToPoint(jvtx.position()).momentum();
 	 vtxP4 += ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double> >(p3.x(),p3.y(),p3.z(),0.13957018);
         charge+=t_trk.track().charge();
	 // DEBUG
         //std::cout << "vtx_trk: pt,eta,phi,N,weight " << t_trk.track().pt() << " " << t_trk.track().eta() << " " << t_trk.track().phi() << " " << indicesToVertex[j]
         //<< " " << jvtx.trackWeight(t_trk) << std::endl;
	 // END DEBUG
         if (indicesToVertex[j] == 1) VtxN1+=1;
         if (indicesToVertex[j] == 2) VtxN2+=1;
         if (ip2dsToVertex[j]>0) nposip2d+=1;
         // hitPattern
         CheckHitPattern::Result res = checkHitPattern_.analyze(iSetup,t_trk.track(),jvtx.vertexState());
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
            if (moms.back().first == 6000111 || moms.back().first==6000112) FromExo+=1;
           } // genParticle found
         } // TrackAssociation  
       } // vtx tracks
     } // tracks to Vertex

     float dR = deltaR(direction.eta(),direction.phi(),vtxP4.eta(),vtxP4.phi());

     reco::Candidate::LorentzVector physicsP4;
     //reco::Candidate::LorentzVector physicsP41,physicsP42,physicsP43;
     if (goodVtx){
       // corrected P4
       physicsP4 = jet1.physicsP4(vtx.position(),jet1,jet1.vertex())
                                                 +jet2.physicsP4(vtx.position(),jet2,jet2.vertex());

       //physicsP41 = detectorP4(jet1,vtx,jvtx,0) + detectorP4(jet2,vtx,jvtx,0);
       //physicsP42 = detectorP4(jet1,vtx,jvtx,1) + detectorP4(jet2,vtx,jvtx,1);
       //physicsP43 = detectorP4(jet1,vtx,jvtx,2) + detectorP4(jet2,vtx,jvtx,2);
     }

     dijetHolder.CorrPt = (goodVtx ? physicsP4.Pt() : -1);
     dijetHolder.CorrEta = (goodVtx ? physicsP4.Eta() : -1);
     dijetHolder.CorrPhi = (goodVtx ? physicsP4.Phi() : -1);
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

   edm::Handle<reco::VertexCollection> recVtxs;
   iEvent.getByLabel("offlinePrimaryVertices", recVtxs);
   if (recVtxs->size()>PV_) 
     pv=recVtxs->at(PV_); 
   else
     pv = recVtxs->front();

   edm::Handle<edm::View<reco::Track> > generalTracks;
   iEvent.getByLabel("generalTracks",generalTracks);
// TransientTrack Builder
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

// Reco to Sim Track association
   if(!iEvent.isRealData() && useTrackingParticles_){
     edm::Handle<std::vector<TrackingParticle> > TPCollectionH ;
     try{
       iEvent.getByLabel(edm::InputTag("mergedtruth","MergedTrackTruth","HLT"),TPCollectionH);
       edm::ESHandle<TrackAssociatorBase> myAssociator;
       iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits", myAssociator);
       RecoToSimColl = myAssociator->associateRecoToSim(generalTracks,TPCollectionH,&iEvent,&iSetup );
     } catch (...) {;}
   }

}

void DJ_DiJetVertices::GetMothers(const TrackingParticle* gp, std::vector<std::pair<int,double> > &moms){

  const TrackingVertex* gv = gp->parentVertex().get();
  if(gv != 0 ){
    for(TrackingVertex::tp_iterator mom = gv->daughterTracks_begin(); mom != gv->daughterTracks_end(); mom++){
      moms.push_back(std::pair<int,double> ( (*mom).get()->genParticles().at(0).get()->pdgId(), gv->position().rho() ));
      if (moms.back().first == 6000111 || moms.back().first == 6000112)
	return;
      GetMothers((*mom).get(),moms);
      break;
    }
  }      
  return ;
}

