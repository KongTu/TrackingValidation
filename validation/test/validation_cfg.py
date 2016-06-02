import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.Timing = cms.Service("Timing")

process.ana = cms.EDAnalyzer('validation',
                      
                		  vertexName = cms.InputTag('offlinePrimaryVertices'),
                		  trackName = cms.InputTag('generalTracks'),             
				  pfCandName = cms.untracked.InputTag('particleFlow'),
				  centralitySrc = cms.InputTag("hiCentrality"),
          			  centralityBinSrc = cms.InputTag("centralityBin","HFtowers"),
				  genSrc = cms.InputTag("genParticles"), 
   				  beamSpotSrc = cms.InputTag("offlineBeamSpot"),				 
				  towerName = cms.InputTag("towerMaker"),
 
				  offlineDCA = cms.untracked.double(3.0),#3.0
				  offlineChi2 = cms.untracked.double(999.0),#0.15
				  offlineptErr = cms.untracked.double(0.1),#0.05
				  offlinenhits = cms.untracked.double(-1),#10

				  doCentrality = cms.untracked.bool(False),
				  doGenParticle = cms.untracked.bool(False)

					
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
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
### conditions
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '75X_dataRun2_v6','')

process.options   = cms.untracked.PSet( wantSummary =
cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1 ) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'/store/user/gsfs/Hydjet_Quenched_MinBias_Drum5Nv75/Hydjet_Quenched_MinBias_Drum5Nv75_RECO_20160318/160318_231321/0000/step3_MinBias_PbPb_RAW2DIGI_L1Reco_RECO_1.root'
#'/store/user/qwang/HIMinimumBias3/HIMinBias_v2/160128_201904/0000/HIMinBias_105.root'

#'root://xrootd-cms.infn.it//store/user/davidlw/Hydjet_Quenched_MinBias_5020GeV_750/ppRECO_std_v3/160216_232421/0000/step2pp_RAW2DIGI_L1Reco_RECO_1.root',
#'/store/user/davidlw/Hydjet_Quenched_MinBias_5020GeV_750/ppRECO_std_v3/160216_232421/0000/step2pp_RAW2DIGI_L1Reco_RECO_102.root',
#'/store/user/davidlw/Hydjet_Quenched_MinBias_5020GeV_750/ppRECO_std_v3/160216_232421/0000/step2pp_RAW2DIGI_L1Reco_RECO_103.root',
#'/store/user/davidlw/Hydjet_Quenched_MinBias_5020GeV_750/ppRECO_std_v3/160216_232421/0000/step2pp_RAW2DIGI_L1Reco_RECO_104.root',
#'/store/user/davidlw/Hydjet_Quenched_MinBias_5020GeV_750/ppRECO_std_v3/160216_232421/0000/step2pp_RAW2DIGI_L1Reco_RECO_105.root'

'/store/hidata/HIRun2015/HIMinimumBias5/AOD/02May2016-v1/10000/006477CE-3326-E611-8C08-003048F317E6.root'
))

process.p1 = cms.Path(process.ana)


process.TFileService = cms.Service("TFileService",fileName = cms.string("test.root"))

