import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuon_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoMuon.MuonIdentification.me0MuonReco_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PH2_1K_FB_V6::All', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
                 

        'root://xrootd.unl.edu//store/user/piet/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_RECO/a6c1ab73bd1959e4a7fbbca874362562/out_reco_1000_1_vCL.root' 
#   '/store/group/upgrade/muon/ME0GlobalReco/ME0MuonReRun_DY_SLHC23patch1_SegmentReRunFullRun_ForPublish/M-20_TuneZ2star_14TeV_6_2_0_SLHC23patch1_2023/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_2023SHCalNoTaper_PU140_Selectors_RECO/b52ce42d5986c94dc336f39e015d825e/output_896_2_3B7.root',
#   '/store/group/upgrade/muon/ME0GlobalReco/ME0MuonReRun_DY_SLHC23patch1_SegmentReRunFullRun_ForPublish/M-20_TuneZ2star_14TeV_6_2_0_SLHC23patch1_2023/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_2023SHCalNoTaper_PU140_Selectors_RECO/b52ce42d5986c94dc336f39e015d825e/output_897_2_1fp.root',
#  '/store/group/upgrade/muon/ME0GlobalReco/ME0MuonReRun_DY_SLHC23patch1_SegmentReRunFullRun_ForPublish/M-20_TuneZ2star_14TeV_6_2_0_SLHC23patch1_2023/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_2023SHCalNoTaper_PU140_Selectors_RECO/b52ce42d5986c94dc336f39e015d825e/output_898_2_zko.root'               

   ),
                            
#    eventsToProcess = cms.untracked.VEventRange('1:115618'),
                            
)

process.looseSeq = cms.EDAnalyzer('ME0TimingAnalysis',
    wp = cms.string("loose"),
    fD = cms.string("ID"), 
    TD = cms.string("false"),
    etaMin = cms.double(2.0),
    etaMax = cms.double(2.8),
    dr = cms.double(0.25),
    ptMin = cms.double(5.0),
   ptMinGen = cms.double(0.0), 
   timeMin = cms.double(5.5),
    timeMax = cms.double(30.5)
)

process.tightSeq = cms.EDAnalyzer('ME0TimingAnalysis',
    wp = cms.string("tight"),
    fD = cms.string("ID"),
    TD = cms.string("false"),
    etaMin = cms.double(2.0),
    etaMax = cms.double(2.8),
    dr = cms.double(0.25),
    ptMin = cms.double(5.0),
    ptMinGen = cms.double(0.0),
    timeMin = cms.double(5.5),
    timeMax = cms.double(30.5)
)
 
process.noIdSeq = cms.EDAnalyzer('ME0TimingAnalysis',
    wp = cms.string("tight"),
    fD = cms.string("noID"), 
    TD = cms.string("false"),
    etaMin = cms.double(2.0),
    etaMax = cms.double(2.8),
    dr = cms.double(0.25),
    ptMin = cms.double(5.0),
    ptMinGen = cms.double(0.0),
    timeMin = cms.double(5.5),
    timeMax = cms.double(30.5)
)

process.assinTimeSeq = cms.EDAnalyzer('ME0TimingAnalysis',
    wp = cms.string("loose"),
    fD = cms.string("ID"),
    TD = cms.string("true"),
    etaMin = cms.double(2.0),
    etaMax = cms.double(2.8),
    dr = cms.double(0.25),
    ptMin = cms.double(5.0),
    ptMinGen = cms.double(0.0),
    timeMin = cms.double(5.5),
    timeMax = cms.double(30.5)
)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string("OutTree.root")
)

process.p = cms.Path(process.looseSeq + process.tightSeq + process.noIdSeq + process.assinTimeSeq)
#process.p = cms.Path(process.looseSeq)
