from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config

config = config()

#datetime object
import datetime
timestamp = datetime.datetime.now().strftime("_%Y%m%d_%H%M%S")

#user inputs
datasetName="DYTo2L_M50"
datasetStr = "/DYJetsToLL_M-50_TuneCH3_13TeV-madgraphMLM-herwig7/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM"

config.General.requestName = 'NanoAnaCrabJob_'+datasetName
config.General.workArea = 'CrabJobs_'+datasetName+'_'+timestamp
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'crab_script.sh'
config.JobType.inputFiles = ['VLLAna_C.so','VLLAna.C','VLLAna.h','ana_crab.C','runana_crab.C','shell_instructions.sh','FrameworkJobReport.xml']
config.JobType.outputFiles = ['skimFile.root']

config.Data.inputDataset = datasetStr
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1  # Number of files per job
config.Data.outLFNDirBase = '/store/user/alaha/NanoAnaCrabJobs'
config.Data.publication = False
config.Data.outputDatasetTag = 'NanoAnaCrabJob_'+datasetName+'_'+timestamp

config.Site.storageSite = 'T3_CH_CERNBOX'
