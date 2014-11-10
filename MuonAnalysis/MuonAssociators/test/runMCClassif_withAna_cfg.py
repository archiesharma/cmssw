import FWCore.ParameterSet.Config as cms
process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.Geometry.GeometryExtended2023Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.GlobalTag.globaltag = 'DES23_62_V1::All'
process.maxEvents = cms.untracked.PSet(
     input = cms.untracked.int32(-1) 
)

#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
#process.Tracer = cms.Service("Tracer")

readFiles = cms.untracked.vstring()
process.source = cms.Source ("PoolSource",fileNames = readFiles)
readFiles.extend( [

#'file:C04E0DBB-8FDB-E311-8C30-003048FFD720.root'

#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_100_1_YJi.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_101_1_VPF.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_102_1_900.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_103_1_GSS.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_104_1_Prc.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_105_1_Pre.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_106_1_sIy.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_107_1_9RE.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_108_1_LQf.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_109_1_1f7.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_10_1_dtP.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_110_1_8ry.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_111_1_2tG.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_112_1_pSe.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_113_1_yKy.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_114_1_Cg5.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_115_1_GLp.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_116_1_xzu.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_117_1_3QP.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_118_1_u4y.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_119_1_Khn.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_11_1_vVC.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_120_1_0rv.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_121_1_dFn.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_122_1_hLE.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_123_1_v28.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_124_1_7kV.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_125_1_rNb.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_126_1_b99.root',
#'root://xrootd.ba.infn.it//store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_127_1_JY7.root'



#'file:/afs/cern.ch/work/m/mileva/punchThrough/CMSSW_6_2_0_SLHC12/src/files_bkg/step3_1000_1_m8t.root',
#'file:/afs/cern.ch/work/m/mileva/punchThrough/CMSSW_6_2_0_SLHC12/src/files_bkg/step3_1001_1_KZe.root',
#'file:/afs/cern.ch/work/m/mileva/punchThrough/CMSSW_6_2_0_SLHC12/src/files_bkg/step3_100_1_5UV.root',
#'file:/afs/cern.ch/work/m/mileva/punchThrough/CMSSW_6_2_0_SLHC12/src/files_bkg/step3_101_1_E8S.root',
#'file:/afs/cern.ch/work/m/mileva/punchThrough/CMSSW_6_2_0_SLHC12/src/files_bkg/step3_102_1_jUo.root',
#'file:/afs/cern.ch/work/m/mileva/punchThrough/CMSSW_6_2_0_SLHC12/src/files_bkg/step3_103_1_MNx.root'


##read files from eos:
#'root://eoscms//eos/cms/store/user/mileva/gemPunchThrough/step3_1000_1_m8t.root',
#'root://eoscms//eos/cms/store/user/mileva/gemPunchThrough/step3_1001_1_KZe.root',
#'root://eoscms//eos/cms/store/user/mileva/gemPunchThrough/step3_100_1_5UV.root'


##without background:
#'file:/afs/cern.ch/work/m/mileva/punchThrough/CMSSW_6_2_0_SLHC12/src/files_Nobkg/step3_1000_1_o1X.root',
#'file:/afs/cern.ch/work/m/mileva/punchThrough/CMSSW_6_2_0_SLHC12/src/files_Nobkg/step3_1001_1_At2.root',
#'file:/afs/cern.ch/work/m/mileva/punchThrough/CMSSW_6_2_0_SLHC12/src/files_Nobkg/step3_100_1_j9X.root',
#'file:/afs/cern.ch/work/m/mileva/punchThrough/CMSSW_6_2_0_SLHC12/src/files_Nobkg/step3_101_1_ENu.root',
#'file:/afs/cern.ch/work/m/mileva/punchThrough/CMSSW_6_2_0_SLHC12/src/files_Nobkg/step3_102_1_UIO.root'

#'/store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_100_1_YJi.root'

#      '/store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_100_1_YJi.root',
#       '/store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_101_1_VPF.root',
#       '/store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_102_1_900.root',
#       '/store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_103_1_GSS.root',
#       '/store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_104_1_Prc.root',
#       '/store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_105_1_Pre.root',
#       '/store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_106_1_sIy.root',
#       '/store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_107_1_9RE.root',
#       '/store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_108_1_LQf.root',
#       '/store/group/upgrade/muon/RecoFolder/MuMinusPt100_2023_3Step/Muminus_Pt100-gun/calabria_MuMinusPt100_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/33104faeb9d0af795150c166b6ba827d/step3_109_1_1f7.root'
 
#       '/store/user/calabria/MinBias_TuneZ2star_14TeV-pythia6/calabria_MinBias2_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/751194ca5d3aeb41ed7baca383591b5f/step3_1000_1_m8t.root',
#       '/store/user/calabria/MinBias_TuneZ2star_14TeV-pythia6/calabria_MinBias2_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/751194ca5d3aeb41ed7baca383591b5f/step3_1001_1_KZe.root',
#       '/store/user/calabria/MinBias_TuneZ2star_14TeV-pythia6/calabria_MinBias2_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/751194ca5d3aeb41ed7baca383591b5f/step3_100_1_5UV.root',
#       '/store/user/calabria/MinBias_TuneZ2star_14TeV-pythia6/calabria_MinBias2_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/751194ca5d3aeb41ed7baca383591b5f/step3_101_1_E8S.root',
#       '/store/user/calabria/MinBias_TuneZ2star_14TeV-pythia6/calabria_MinBias2_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/751194ca5d3aeb41ed7baca383591b5f/step3_102_1_jUo.root',
#       '/store/user/calabria/MinBias_TuneZ2star_14TeV-pythia6/calabria_MinBias2_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/751194ca5d3aeb41ed7baca383591b5f/step3_103_1_MNx.root',
#       '/store/user/calabria/MinBias_TuneZ2star_14TeV-pythia6/calabria_MinBias2_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/751194ca5d3aeb41ed7baca383591b5f/step3_104_1_57O.root',
#       '/store/user/calabria/MinBias_TuneZ2star_14TeV-pythia6/calabria_MinBias2_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/751194ca5d3aeb41ed7baca383591b5f/step3_105_1_SeP.root',
#       '/store/user/calabria/MinBias_TuneZ2star_14TeV-pythia6/calabria_MinBias2_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/751194ca5d3aeb41ed7baca383591b5f/step3_106_1_rBt.root',
#       '/store/user/calabria/MinBias_TuneZ2star_14TeV-pythia6/calabria_MinBias2_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/751194ca5d3aeb41ed7baca383591b5f/step3_107_1_cuW.root'
       '/store/user/calabria/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/calabria_DYToMuMu_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/7e82c0353afc1625567b3844e4c947b8/step3_1000_1_HIo.root',
       '/store/user/calabria/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/calabria_DYToMuMu_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/7e82c0353afc1625567b3844e4c947b8/step3_1001_1_PG2.root',
       '/store/user/calabria/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/calabria_DYToMuMu_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/7e82c0353afc1625567b3844e4c947b8/step3_1002_1_tVX.root',
       '/store/user/calabria/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/calabria_DYToMuMu_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/7e82c0353afc1625567b3844e4c947b8/step3_1003_1_0rr.root',
       '/store/user/calabria/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/calabria_DYToMuMu_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_3Step/7e82c0353afc1625567b3844e4c947b8/step3_1004_1_uAL.root'


]);

process.checkMuSize = cms.EDFilter("MuonsizeFilter"
                                            )
## require scraping filter
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                     applyfilter = cms.untracked.bool(True),
                                     debugOn = cms.untracked.bool(False),
                                     numtrack = cms.untracked.uint32(10),
                                     thresh = cms.untracked.double(0.2)
                                     )

process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                            vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                            minimumNDOF = cms.uint32(4) ,
                                            maxAbsZ = cms.double(24),
                                            maxd0 = cms.double(2)
                                            )


#from SimGeneral.MixingModule.mixNoPU_cfi                          import *
#from SimGeneral.TrackingAnalysis.trackingParticlesNoSimHits_cfi   import * 
#from SimMuon.MCTruth.MuonAssociatorByHitsESProducer_NoSimHits_cfi import *

process.load('MuonAnalysis.MuonAssociators.muonClassifwithAna_cfi')

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('Test_DY_withUpg.root')
)

process.p = cms.Path(
                     process.checkMuSize*
                     process.scrapingVeto*
        	     process.primaryVertexFilter*
                     process.classByHits
		 ) 
