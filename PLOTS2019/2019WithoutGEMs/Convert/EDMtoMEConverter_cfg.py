import FWCore.ParameterSet.Config as cms

process = cms.Process("EDMtoMEConvert")

from FWCore.ParameterSet.VarParsing import VarParsing

import sys

#---------------------------------------------------------------------------------------------------------

def readFileList(fileList, inputFileName, fileNamePrefix):
    inputFile = open(inputFileName, 'r')
    for name in inputFile:
	name1 = name.replace("',", "")
	name2 = name1.replace("'", "")
        fileList.extend([ fileNamePrefix + name2 ])
    inputFile.close()

def parseAndApplyOptions(process) :
    options.register ('globalTag', 'START53_V7A::All', VarParsing.multiplicity.singleton,
                      VarParsing.varType.string, "Global Tag to use.")
    options.register ('fileList', 'fileList.txt', VarParsing.multiplicity.singleton,
                      VarParsing.varType.string, "List of root files to process.")
    options.register ('fileNamePrefix', '', VarParsing.multiplicity.singleton,
                      VarParsing.varType.string, "Prefix to add to input file names.")
    options.register ('output', 'Tree.root', VarParsing.multiplicity.singleton,
                      VarParsing.varType.string, "Tree root file.")

    if len(sys.argv) > 0:
        last = sys.argv.pop()
        sys.argv.extend(last.split(","))

    options.parseArguments()

    #process.GlobalTag.globaltag = options.globalTag

    readFileList(process.source.fileNames, options.fileList, options.fileNamePrefix)
    #process.out.fileName = options.outputFile
    process.dqmSaver.workflow = "/GlobalValidation/Test/RECO" + options.output
    #process.maxEvents.input = options.maxEvents

    return

#---------------------------------------------------------------------------------------------------------

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("DQMServices.Components.EDMtoMEConverter_cff")
process.load("DQMServices.Components.DQMEnvironment_cfi")
#process.load("Configuration.StandardSequences.FakeConditions_cff")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
process.load("Validation.Configuration.postValidation_cff")
process.load("Validation.RecoMuon.PostProcessorHLT_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(



    ),
    #processingMode = cms.untracked.string('RunsAndLumis')
)

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    fileMode = cms.untracked.string('FULLMERGE')
)

process.DQMStore.referenceFileName = ""
process.DQMStore.collateHistograms = True

process.dqmSaver.convention = "Offline"
#Settings equivalent to 'RelVal' convention:
process.dqmSaver.saveByRun = cms.untracked.int32(-1)
process.dqmSaver.saveAtJobEnd = cms.untracked.bool(True)
process.dqmSaver.forceRunNumber = cms.untracked.int32(1)
#End of 'RelVal convention settings
process.dqmSaver.workflow = "/GlobalValidation/Test/RECO"

options = VarParsing ('python')
parseAndApplyOptions(process)

from Validation.RecoMuon.PostProcessor_cff import *
process.p1 = cms.Path(process.EDMtoMEConverter*
                      #process.postValidation*
		      process.recoMuonPostProcessors*
                      #process.recoMuonPostProcessorsHLT*
                      process.dqmSaver)
