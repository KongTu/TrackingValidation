### this is an example for running on RECO
### the name must be changed crab.cfg for actual running

from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()
config.General.requestName = 'HydjetTune_Drum5N_v2'
config.General.workArea = 'HydjetTune_Drum5N_v2'
config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.allowUndistributedCMSSW = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'validation_cfg.py'

config.Data.inputDBS = 'phys03'
config.Data.inputDataset = '/Hydjet_Quenched_MinBias_Drum5Nv75/gsfs-Hydjet_Quenched_MinBias_Drum5Nv75_RECO_20160318-a78b3072e46c33a0fe1fffcdcc303b7d/USER'
config.Data.splitting = 'FileBased'
config.Data.ignoreLocality = False
config.Data.unitsPerJob = 40
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Site.storageSite = 'T2_US_MIT'
config.Site.whitelist = ['T2_US_MIT']
