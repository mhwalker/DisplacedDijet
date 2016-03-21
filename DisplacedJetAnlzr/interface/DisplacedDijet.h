#ifndef DisplacedDijet_DisplacedJetAnlzr_DisplacedDijet_hh
#define DisplacedDijet_DisplacedJetAnlzr_DisplacedDijet_hh


struct DisplacedDijet {
  float CorrPt,CorrEta,CorrPhi,CorrMass;
  int NPromptTracks,NDispTracks,NPromptTracks1,NPromptTracks2;
  float PromptEnergyFrac,Lxy,Lxysig,VtxX,VtxY,VtxZ,VtxChi2,Vtxmass,Vtxpt;
  float PromptEnergyFrac1,PromptEnergyFrac2;
  int VtxN,VtxN1,VtxN2;
  float VtxdR,VtxCharge,TrkAvgPt,Posip2dFrac,NAvgMissHitsAfterVert;
  float NAvgHitsInFrontOfVert,ExoVtxFrac;
  float glxydistall,glxydistvtx,glxydistclr,glxyrmsall,glxyrmsvtx,glxyrmsclr;
  int Nclusters,bestclusterN,bestclusterN1,bestclusterN2;
  float bestclusterlxy;
  float alphaMax;
  float openingAngle;
  float medianIPLog10Sig;
  float medianTrackAngle;
  float medianTrackAngleLog10;
  float totalTrackAngle;
  float trackMass1,trackMass2;
  float trackEnergy1,trackEnergy2;

};

#endif
