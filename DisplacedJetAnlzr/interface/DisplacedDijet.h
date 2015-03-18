#ifndef DisplacedDijet_DisplacedJetAnlzr_DisplacedDijet_hh
#define DisplacedDijet_DisplacedJetAnlzr_DisplacedDijet_hh


struct DisplacedDijet {
  float CorrPt,CorrEta,CorrPhi,CorrMass;
  int NPromptTracks,NDispTracks;
  float PromptEnergyFrac,Lxy,Lxysig,VtxX,VtxY,VtxZ,VtxChi2,Vtxmass,Vtxpt;
  int VtxN,VtxN1,VtxN2;
  float VtxdR,VtxCharge,TrkAvgPt,Posip2dFrac,NAvgMissHitsAfterVert;
  float NAvgHitsInFrontOfVert,ExoVtxFrac;
  float glxydistall,glxydistvtx,glxydistclr,glxyrmsall,glxyrmsvtx,glxyrmsclr;
  int Nclusters,bestclusterN,bestclusterN1,bestclusterN2;
  float bestclusterlxy;

};

#endif
