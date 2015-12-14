import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuon_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PH2_1K_FB_V6::All', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#  'file:/afs/cern.ch/work/a/archie/public/out_digi.root'                
    'file:/tmp/archie/out_digi.root'
#        'root://xrootd.unl.edu//store/user/piet/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_RECO/a6c1ab73bd1959e4a7fbbca874362562/out_reco_1000_1_vCL.root' 
#   '/store/group/upgrade/muon/ME0GlobalReco/ME0MuonReRun_DY_SLHC23patch1_SegmentReRunFullRun_ForPublish/M-20_TuneZ2star_14TeV_6_2_0_SLHC23patch1_2023/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_2023SHCalNoTaper_PU140_Selectors_RECO/b52ce42d5986c94dc336f39e015d825e/output_896_2_3B7.root',
#   '/store/group/upgrade/muon/ME0GlobalReco/ME0MuonReRun_DY_SLHC23patch1_SegmentReRunFullRun_ForPublish/M-20_TuneZ2star_14TeV_6_2_0_SLHC23patch1_2023/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_2023SHCalNoTaper_PU140_Selectors_RECO/b52ce42d5986c94dc336f39e015d825e/output_897_2_1fp.root',
#  '/store/group/upgrade/muon/ME0GlobalReco/ME0MuonReRun_DY_SLHC23patch1_SegmentReRunFullRun_ForPublish/M-20_TuneZ2star_14TeV_6_2_0_SLHC23patch1_2023/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_2023SHCalNoTaper_PU140_Selectors_RECO/b52ce42d5986c94dc336f39e015d825e/output_898_2_zko.root'               

   ),
)

process.digiSeq = cms.EDAnalyzer('GE21DigiAnalyzer',
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("OutDigiTree.root")
)

process.p = cms.Path(process.digiSeq)
