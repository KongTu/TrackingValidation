import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.Timing = cms.Service("Timing")

process.ana = cms.EDAnalyzer('validation',
                      
                		  vertexSrc = cms.string('hiSelectedVertex'),
                		  trackSrc = cms.InputTag('hiGeneralTracks'),             
				  pfCandSrc = cms.untracked.InputTag('particleFlowTmp'),
				  centralitySrc = cms.InputTag("hiCentrality"),
          						centralityBinSrc = cms.InputTag("centralityBin","HFtowers"),
 
				  
				  offlineDCA = cms.untracked.double(3.0),#3.0
				  offlineChi2 = cms.untracked.double(0.15),#0.15
				  offlineptErr = cms.untracked.double(0.1),#0.05
				  offlinenhits = cms.untracked.double(11),#10

				  doCentrality = cms.untracked.bool(True)

					
)

### standard includes
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

### conditions
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '75X_mcRun2_HeavyIon_v1','')

process.options   = cms.untracked.PSet( wantSummary =
cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1 ) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'root://xrootd3.cmsaf.mit.edu//store/user/qwang/HIHardProbes/HIHardProbes_FullTrackSkim2015_v3/151216_192437/0000/FullTrack_1.root'
#'root://xrootd3.cmsaf.mit.edu//store/user/qwang/HIHardProbes/HIHardProbes_FullTrackSkim2015_v3/151216_192437/0000/FullTrack_10.root'
'/store/user/qwang/HIHardProbes/HIHardProbes_FullTrackSkim2015_v3/151216_192437/0000/FullTrack_100.root'
))

process.p1 = cms.Path(process.ana)


process.TFileService = cms.Service("TFileService",fileName = cms.string("validation.root"))

