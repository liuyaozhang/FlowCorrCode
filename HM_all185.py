import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM.HLTPaths = ['HLT_PAPixelTracks_Multiplicity100_v*','HLT_PAPixelTracks_Multiplicity130_v*','HLT_PAPixelTracks_Multiplicity160_v*']
process.hltHM.andOr = cms.bool(True)
process.hltHM.throw = cms.bool(False)

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
        )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(5000)
process.options   = cms.untracked.PSet( wantSummary = 
cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) 
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
                )
                            )
#import sys
#import os
#inputfilename = open("HighMultCorrelation.txt", "r")
#ffrom=int(sys.argv[2])
#fto=int(sys.argv[3])
#for j in range(1,ffrom+1):
#    inputfilename.readline()
#for i in range(ffrom,fto):
#    process.source.fileNames.append(inputfilename.readline())

process.demo = cms.EDAnalyzer('DemoAnalyzerAll',
multMax = cms.untracked.double(220),
                multMin = cms.untracked.double(185),
ptcut_ks = cms.untracked.vdouble(0.2,0.6,1.2,1.8,2.6,3.6,4.6,6.0,9.0,12.0),
ptcut_la = cms.untracked.vdouble(0.2,0.6,1.2,1.8,2.6,3.6,4.6,6.0,9.0,12.0),
ptbin_n = cms.untracked.int32(9),
sigma_ks = 
cms.untracked.vdouble(0.00923285,0.00762145,0.00671416,0.0065142,0.00664755,0.00685524,0.00711999,0.00759312,0.00735706),
mean_ks = cms.untracked.vdouble(0.497098,0.497614,0.497662,0.497585,0.497546,0.497512,0.49757,0.497572,0.497687),
sigma_la = 
cms.untracked.vdouble(0.00361575,0.00361575,0.00322582,0.00305862,0.00289019,0.00285282,0.00294214,0.00291627,0.00461972),
mean_la = cms.untracked.vdouble(1.11617,1.11617,1.11592,1.11593,1.11588,1.11583,1.11582,1.11579,1.11581)
                              )

process.TFileService = cms.Service("TFileService",
                                       fileName = 
cms.string('HighMultCorrelation185_220_all.root')
                                   )


process.p = cms.Path(process.hltHM * process.demo)
