import FWCore.ParameterSet.Config as cms

processName = "MuonSuite"
process = cms.Process(processName)

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( (
       '/store/user/calabria/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/calabria_DYToMuMu_AODSIM_CMSSW_6_2_0_SLHC12_2023Scenario_Case1_4/314265c4c851b22933fa9c86eb7294b1/step3_1000_3_8Gu.root',
       '/store/user/calabria/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/calabria_DYToMuMu_AODSIM_CMSSW_6_2_0_SLHC12_2023Scenario_Case1_4/314265c4c851b22933fa9c86eb7294b1/step3_1001_3_Emn.root',
       '/store/user/calabria/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/calabria_DYToMuMu_AODSIM_CMSSW_6_2_0_SLHC12_2023Scenario_Case1_4/314265c4c851b22933fa9c86eb7294b1/step3_1002_2_d8r.root',
       '/store/user/calabria/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/calabria_DYToMuMu_AODSIM_CMSSW_6_2_0_SLHC12_2023Scenario_Case1_4/314265c4c851b22933fa9c86eb7294b1/step3_1003_1_6sG.root',
       '/store/user/calabria/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/calabria_DYToMuMu_AODSIM_CMSSW_6_2_0_SLHC12_2023Scenario_Case1_4/314265c4c851b22933fa9c86eb7294b1/step3_1005_2_ODk.root',
       '/store/user/calabria/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/calabria_DYToMuMu_AODSIM_CMSSW_6_2_0_SLHC12_2023Scenario_Case1_4/314265c4c851b22933fa9c86eb7294b1/step3_1006_2_BLd.root',
       '/store/user/calabria/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/calabria_DYToMuMu_AODSIM_CMSSW_6_2_0_SLHC12_2023Scenario_Case1_4/314265c4c851b22933fa9c86eb7294b1/step3_1007_4_HGF.root',
       '/store/user/calabria/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/calabria_DYToMuMu_AODSIM_CMSSW_6_2_0_SLHC12_2023Scenario_Case1_4/314265c4c851b22933fa9c86eb7294b1/step3_1008_1_Xds.root',
       '/store/user/calabria/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/calabria_DYToMuMu_AODSIM_CMSSW_6_2_0_SLHC12_2023Scenario_Case1_4/314265c4c851b22933fa9c86eb7294b1/step3_1009_3_FKR.root',
    ))
secFiles.extend((

    ))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('drop *', "keep *_MEtoEDMConverter_*_"+processName),
    fileName = cms.untracked.string('validationEDM.root')
)
process.outpath = cms.EndPath(process.out)

process.load('Configuration/StandardSequences/RawToDigi_cff')
process.raw2digi_step = cms.Path(process.RawToDigi)

process.load("Configuration/StandardSequences/SimulationRandomNumberGeneratorSeeds_cff")

process.load("DQMServices.Components.MEtoEDMConverter_cfi")
process.MEtoEDMConverter_step = cms.Path(process.MEtoEDMConverter)

process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.Geometry.GeometryExtended2023MuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023Muon_cff')
#process.load('Configuration.Geometry.GeometryExtended2019Reco_cff')
#process.load('Configuration.Geometry.GeometryExtended2019_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgrade2019', '')

#---- Validation stuffs ----#
## Default validation modules
process.load("Configuration.StandardSequences.Validation_cff")
process.validation_step = cms.Path(process.validation)
## Load muon validation modules
#process.recoMuonVMuAssoc.outputFileName = 'validationME.root'
process.muonValidation_step = cms.Path(process.recoMuonValidation)

process.schedule = cms.Schedule(
    #process.raw2digi_step,
    #process.validation_step,
    process.muonValidation_step,
    process.MEtoEDMConverter_step,process.outpath)

