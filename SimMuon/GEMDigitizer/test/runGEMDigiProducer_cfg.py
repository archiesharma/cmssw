import FWCore.ParameterSet.Config as cms

process = cms.Process("NewGEMDIGI")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023MuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023Muon_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

# customization of the process.pdigi sequence to add the GEM digitizer
from SimMuon.GEMDigitizer.customizeGEMDigi import customize_digi_addGEM_muon_only
process = customize_digi_addGEM_muon_only(process)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'file:out_sim.root'
         '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_1.root',
          '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_10.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_101.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_102.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_103.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_104.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_105.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_107.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_108.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_11.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_110.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_111.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_112.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_113.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_114.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_115.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_116.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_117.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_118.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_119.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_12.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_120.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_121.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_122.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_123.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_124.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_125.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_126.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_127.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_128.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_129.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_13.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_130.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_131.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_132.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_133.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_134.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_135.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_136.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_137.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_138.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_139.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_14.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_140.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_141.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_142.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_143.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_144.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_145.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_146.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_147.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_148.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_149.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_15.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_150.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_151.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_152.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_153.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_154.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_155.root',
       '/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_50ps_1ns_v2_RECO/151021_081213/0000/out_reco_156.root'

    )
)

process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(
        'file:/tmp/archie/out_digi.root'

    ),

    outputCommands = cms.untracked.vstring(
        'keep  *_*_*_*',
        #'drop CastorDataFramesSorted_simCastorDigis_*_GEMDIGI'
        # drop all CF stuff
        ##'drop *_mix_*_*',
        # drop tracker simhits
        ##'drop PSimHits_*_Tracker*_*',
        # drop calorimetry stuff
        ##'drop PCaloHits_*_*_*',
        # clean up simhits from other detectors
        ##'drop PSimHits_*_Totem*_*',
        ##'drop PSimHits_*_FP420*_*',
        ##'drop PSimHits_*_BSC*_*',
        # drop some not useful muon digis and links
        ##'drop *_*_MuonCSCStripDigi_*',
        ##'drop *_*_MuonCSCStripDigiSimLinks_*',
        #'drop *SimLink*_*_*_*',
        ##'drop *RandomEngineStates_*_*_*',
        ##'drop *_randomEngineStateProducer_*_*'
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('digi_step')
    )
)

process.digi_step    = cms.Path(process.pdigi)
process.endjob_step  = cms.Path(process.endOfProcess)
process.out_step     = cms.EndPath(process.output)

process.schedule = cms.Schedule(
    process.digi_step,
    process.endjob_step,
    process.out_step
)
