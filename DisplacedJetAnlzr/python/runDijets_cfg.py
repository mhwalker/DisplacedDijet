import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("RutgersAOD")
options = VarParsing.VarParsing ('analysis')

#set default arguments
options.inputFiles='/store/user/mwalker/SIGNAL/Displaced/LL_sbottom_500_100.0mm/step2_LL_sbottom_500_100.0mm_9.root'
options.outputFile='test.root'
#options.inputFiles= '/store/relval/CMSSW_7_0_0/RelValProdTTbar_13/AODSIM/POSTLS170_V3-v2/00000/40D11F5C-EA98-E311-BE17-02163E00E964.root'
#options.inputFiles= 'file:/cms/thomassen/2012/Signal/StopRPV/store/aodsim/LLE122/StopRPV_8TeV_chunk3_stop950_bino800_LLE122_aodsim.root'
#options.maxEvents = 100 # -1 means all events
#options.maxEvents = 100

# get and parse the command line arguments
options.parseArguments()

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
)

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("TrackingTools.MaterialEffects.MaterialPropagator_cfi")
process.load('PhysicsTools.PatAlgos.patSequences_cff')

process.load("DisplacedDijet.DisplacedJetAnlzr.DJ_DiJetVertices_cfi")

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
process.out = cms.OutputModule(
  "PoolOutputModule"
, fileName = cms.untracked.string(options.outputFile )
, outputCommands = cms.untracked.vstring(
    *patEventContentNoCleaning
  )
)
process.out.outputCommands += [ 'drop recoGenJets_*_*_*', 'keep *_*_DisplacedDijets_*' ]

print process.out.outputCommands

process.p = cms.Path(
    process.particleFlowPtrs *
    process.patCandidates *
    process.selectedPatCandidates *
    process.djdijetvertices
    )
process.outpath = cms.EndPath(process.out)
