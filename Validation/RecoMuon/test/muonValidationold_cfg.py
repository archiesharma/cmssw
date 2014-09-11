import FWCore.ParameterSet.Config as cms

processName = "MuonSuite"
process = cms.Process(processName)

#readFiles = cms.untracked.vstring()
#secFiles = cms.untracked.vstring()
#process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
#readFiles.extend( (
#      '/store/user/calabria/calabria_SingleMuPt10_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case5/calabria_SingleMuPt10_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case5/b6da911bb96a528ccff7d001f1ad0373/out_reco_102_2_w1A.root'
#    '/store/relval/CMSSW_2_2_0/RelValSingleMuPt10/GEN-SIM-RECO/IDEAL_V9_v1/0000/10C12A24-74B9-DD11-85B2-001617DBCF6A.root',
#    '/store/relval/CMSSW_2_2_0/RelValSingleMuPt10/GEN-SIM-RECO/IDEAL_V9_v1/0000/3A14ADED-B4B9-DD11-8F0B-001617E30D40.root' ,
#    ))
#secFiles.extend((
#     '/store/user/calabria/calabria_SingleMuPt10_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case5/calabria_SingleMuPt10_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case5/b6da911bb96a528ccff7d001f1ad0373/out_reco_103_2_tsF.root'
#    '/store/relval/CMSSW_2_2_0/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/IDEAL_V9_v1/0000/526D7CD6-68B9-DD11-886D-001617DBD224.root',
#    '/store/relval/CMSSW_2_2_0/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/IDEAL_V9_v1/0000/A2C70EEE-B4B9-DD11-8170-001617DBD316.root',
#    '/store/relval/CMSSW_2_2_0/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/IDEAL_V9_v1/0000/D437AC21-6FB9-DD11-BEA1-001617E30CC8.root' 
#    ))

process.source = cms.Source("PoolSource",

    fileNames = cms.untracked.vstring(
#     '/store/relval/CMSSW_2_2_0/RelValSingleMuPt10/GEN-SIM-RECO/IDEAL_V9_v1/0000/10C12A24-74B9-DD11-85B2-001617DBCF6A.root'

#    '/store/user/calabria/calabria_SingleMuPt10_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case5/calabria_SingleMuPt10_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case5/b6da911bb96a528ccff7d001f1ad0373/out_reco_103_2_tsF.root',
#'/store/user/calabria/calabria_SingleMuPt10_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case5/calabria_SingleMuPt10_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case5/b6da911bb96a528ccff7d001f1ad0373/out_reco_102_2_w1A.root',
#'/store/user/calabria/calabria_SingleMuPt10_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case5/calabria_SingleMuPt10_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case5/b6da911bb96a528ccff7d001f1ad0373/out_reco_104_1_Urh.root',
#'/store/user/calabria/calabria_SingleMuPt10_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case5/calabria_SingleMuPt10_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case5/b6da911bb96a528ccff7d001f1ad0373/out_reco_105_1_Y0x.root',
#'/store/user/calabria/calabria_SingleMuPt10_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case5/calabria_SingleMuPt10_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case5/b6da911bb96a528ccff7d001f1ad0373/out_reco_106_1_eAg.root',
#'/store/user/calabria/calabria_SingleMuPt10_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case5/calabria_SingleMuPt10_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case5/b6da911bb96a528ccff7d001f1ad0373/out_reco_107_1_90G.root',
#'/store/user/calabria/calabria_SingleMuPt10_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case5/calabria_SingleMuPt10_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case5/b6da911bb96a528ccff7d001f1ad0373/out_reco_108_1_uEx.root',
#'/store/user/calabria/calabria_SingleMuPt10_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case5/calabria_SingleMuPt10_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case5/b6da911bb96a528ccff7d001f1ad0373/out_reco_109_1_Uqq.root',
#' /store/user/calabria/calabria_SingleMuPt10_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case5/calabria_SingleMuPt10_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case5/b6da911bb96a528ccff7d001f1ad0373/out_reco_10_2_tr8.root',
#'/store/user/calabria/calabria_SingleMuPt10_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case5/calabria_SingleMuPt10_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case5/b6da911bb96a528ccff7d001f1ad0373/out_reco_110_1_SOe.root'
#'/store/user/calabria/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step_NoUpg/751194ca5d3aeb41ed7baca383591b5f/step3_1000_1_PE3.root',
#'/store/user/calabria/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step_NoUpg/751194ca5d3aeb41ed7baca383591b5f/step3_100_1_Jfn.root'
#       '/store/user/calabria/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step_NoUpg/751194ca5d3aeb41ed7baca383591b5f/step3_1000_1_PE3.root',
#       '/store/user/calabria/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step_NoUpg/751194ca5d3aeb41ed7baca383591b5f/step3_100_1_Jfn.root'

#'/store/user/calabria/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step_NoUpg/751194ca5d3aeb41ed7baca383591b5f/step3_1000_1_PE3.root',
       '/store/user/calabria/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step_NoUpg/751194ca5d3aeb41ed7baca383591b5f/step3_100_1_Jfn.root'

#'file:step3_1000_1_PE3.root'

)
)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('drop *', "keep *_MEtoEDMConverter_*_"+processName),
    fileName = cms.untracked.string('validationEDM.root')
)
process.outpath = cms.EndPath(process.out)

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.MessageLogger.categories = cms.untracked.vstring('MuonAssociatorEDProducer', 'MuonTrackProducer', 'MuonAssociatorByHits', 'DTHitAssociator', 'RPCHitAssociator', 'GEMHitAssociator', 'MuonTruth', 'MuonTrackValidator', 'FwkJob', 'FwkReport', 'FwkSummary', 'Root_NoDictionary')

process.MessageLogger.cerr = cms.untracked.PSet(
    noTimeStamps = cms.untracked.bool(True),

    threshold = cms.untracked.string('WARNING'),
 
    MuonAssociatorEDProducer = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ), 
    MuonTrackProducer = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
    MuonAssociatorByHits = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ), 
    DTHitAssociator = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
    RPCHitAssociator = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
    GEMHitAssociator = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
   MuonTruth = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
    MuonTrackValidator = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    )
)

process.MessageLogger.cout = cms.untracked.PSet(
    noTimeStamps = cms.untracked.bool(True),

#    threshold = cms.untracked.string('DEBUG'),
    threshold = cms.untracked.string('INFO'),

    default = cms.untracked.PSet(
       limit = cms.untracked.int32(0)
#     limit = cms.untracked.int32(10000000)
    ),
    MuonAssociatorEDProducer = cms.untracked.PSet(
        limit = cms.untracked.int32(10000000)
    ),
    MuonTrackProducer = cms.untracked.PSet(
        limit = cms.untracked.int32(10000000)
    ),
    MuonAssociatorByHits = cms.untracked.PSet(
        limit = cms.untracked.int32(10000000)
    ),
    DTHitAssociator = cms.untracked.PSet(
        limit = cms.untracked.int32(10000000)
    ),
    RPCHitAssociator = cms.untracked.PSet(
        limit = cms.untracked.int32(10000000)
    ),
   GEMHitAssociator = cms.untracked.PSet(
        limit = cms.untracked.int32(10000000)
    ),
   MuonTruth = cms.untracked.PSet(
        limit = cms.untracked.int32(10000000)
    ),
    MuonTrackValidator = cms.untracked.PSet(
        limit = cms.untracked.int32(10000000)
    ),
    FwkReport = cms.untracked.PSet(
        reportEvery = cms.untracked.int32(1),
        limit = cms.untracked.int32(10000000)
    ),
    FwkSummary = cms.untracked.PSet(
        reportEvery = cms.untracked.int32(1),
        limit = cms.untracked.int32(10000000)
    ),
    FwkJob = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
    Root_NoDictionary = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    )
)

process.MessageLogger.statistics = cms.untracked.vstring('cout')

process.load('Configuration/StandardSequences/RawToDigi_cff')
process.raw2digi_step = cms.Path(process.RawToDigi)

process.load("Configuration/StandardSequences/SimulationRandomNumberGeneratorSeeds_cff")

process.load("DQMServices.Components.MEtoEDMConverter_cfi")
process.MEtoEDMConverter_step = cms.Path(process.MEtoEDMConverter)

process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.Geometry.GeometryExtended2023MuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023Muon_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
#process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')


#process.GlobalTag.globaltag = "IDEAL_V9::All"

#---- Validation stuffs ----#
## Default validation modules
process.load("Configuration.StandardSequences.Validation_cff")
process.validation_step = cms.Path(process.validation)
## Load muon validation modules
#process.recoMuonVMuAssoc.outputFileName = 'validationME.root'
#process.muonValidation_step = cms.Path(cms.SequencePlaceholder("mix")+process.recoMuonValidation)
process.muonValidation_step = cms.Path(process.recoMuonValidation)

process.schedule = cms.Schedule(
    process.raw2digi_step,
#    process.validation_step,
    process.muonValidation_step,
    process.MEtoEDMConverter_step,process.outpath)

