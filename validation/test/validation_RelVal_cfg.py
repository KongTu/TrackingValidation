import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.Timing = cms.Service("Timing")

process.ana = cms.EDAnalyzer('validation',
                      
                		  vertexName = cms.InputTag('hiSelectedVertex'),
                		  trackName = cms.InputTag('hiGeneralTracks'),             
				  pfCandName = cms.untracked.InputTag('particleFlowTmp'),
				  centralitySrc = cms.InputTag("hiCentrality"),
          			  centralityBinSrc = cms.InputTag("centralityBin","HFtowers"),
				  towerName = cms.InputTag("towerMaker"), 
				  
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
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
### conditions
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_v8','')

process.options   = cms.untracked.PSet( wantSummary =
cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1 ) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'/store/group/phys_heavyions/velicanu/trackingvalidation/CMSSW_8_0_0_pre1/76X_dataRun1_v10/step2_1.root',
#'/store/group/phys_heavyions/velicanu/trackingvalidation/CMSSW_8_0_0_pre1/76X_dataRun1_v10/step2_2.root',
#'/store/group/phys_heavyions/velicanu/trackingvalidation/CMSSW_8_0_0_pre1/76X_dataRun1_v10/step2_3.root',
#'/store/group/phys_heavyions/velicanu/trackingvalidation/CMSSW_8_0_0_pre1/76X_dataRun1_v10/step2_4.root',
#'/store/group/phys_heavyions/velicanu/trackingvalidation/CMSSW_8_0_0_pre1/76X_dataRun1_v10/step2_5.root',
#'/store/group/phys_heavyions/velicanu/trackingvalidation/CMSSW_8_0_0_pre1/76X_dataRun1_v10/step2_6.root',
#'/store/group/phys_heavyions/velicanu/trackingvalidation/CMSSW_8_0_0_pre1/76X_dataRun1_v10/step2_7.root',
#'/store/group/phys_heavyions/velicanu/trackingvalidation/CMSSW_8_0_0_pre1/76X_dataRun1_v10/step2_8.root',
#'/store/group/phys_heavyions/velicanu/trackingvalidation/CMSSW_8_0_0_pre1/76X_dataRun1_v10/step2_9.root',
#'/store/group/phys_heavyions/velicanu/trackingvalidation/CMSSW_8_0_0_pre1/76X_dataRun1_v10/step2_10.root'

'/store/group/phys_heavyions/velicanu/trackingvalidation/CMSSW_8_0_0_pre6/80X_dataRun2_v8/step2_1.root',
'/store/group/phys_heavyions/velicanu/trackingvalidation/CMSSW_8_0_0_pre6/80X_dataRun2_v8/step2_2.root',
'/store/group/phys_heavyions/velicanu/trackingvalidation/CMSSW_8_0_0_pre6/80X_dataRun2_v8/step2_3.root',
'/store/group/phys_heavyions/velicanu/trackingvalidation/CMSSW_8_0_0_pre6/80X_dataRun2_v8/step2_4.root',
'/store/group/phys_heavyions/velicanu/trackingvalidation/CMSSW_8_0_0_pre6/80X_dataRun2_v8/step2_5.root',
'/store/group/phys_heavyions/velicanu/trackingvalidation/CMSSW_8_0_0_pre6/80X_dataRun2_v8/step2_6.root',
'/store/group/phys_heavyions/velicanu/trackingvalidation/CMSSW_8_0_0_pre6/80X_dataRun2_v8/step2_7.root',
'/store/group/phys_heavyions/velicanu/trackingvalidation/CMSSW_8_0_0_pre6/80X_dataRun2_v8/step2_8.root',
'/store/group/phys_heavyions/velicanu/trackingvalidation/CMSSW_8_0_0_pre6/80X_dataRun2_v8/step2_9.root',
'/store/group/phys_heavyions/velicanu/trackingvalidation/CMSSW_8_0_0_pre6/80X_dataRun2_v8/step2_10.root'
))

process.p1 = cms.Path(process.ana)


process.TFileService = cms.Service("TFileService",fileName = cms.string("test.root"))

