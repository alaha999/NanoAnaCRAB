import FWCore.ParameterSet.Config as cms

process = cms.Process("NANO")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

# Define input source (managed by process.source)
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring())

# Locally, you can specify the input files for debugging purposes
process.source.fileNames = [
    # Add files here for local testing:
    # Example file path for local debugging
    #"file:/afs/cern.ch/user/p/phazarik/work/VLLSearch/AnaScript/Skimmer/test_inputs/DYJetsToLL_M-50_2018UL.root"
    "/store/mc/RunIISummer20UL18NanoAODv9/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2540000/10995124-9BEC-5A48-8D83-30542D27979D.root"
]

# Define output file (using TFileService for ROOT output)
#process.TFileService = cms.Service("TFileService", fileName = cms.string("skimFile.root"))

# Define custom parameters inside a cms.PSet
#process.custom_params = cms.PSet(
#    codedir_ = cms.string("/afs/cern.ch/user/p/phazarik/work/VLLSearch/AnaScript/Skimmer"),
#    data_flag_ = cms.string("0"),  # Changed to cms.string to be consistent
#    campaign_ = cms.string("2018_UL"),
#    flag_ = cms.string("dy"),
#    samplename_ = cms.string("DYJetsToLL_M50")
#)
