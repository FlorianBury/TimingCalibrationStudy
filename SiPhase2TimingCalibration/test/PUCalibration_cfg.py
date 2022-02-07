from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')

options.register('N',
                 100,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Number of events to be processed")
options.register('threshold',
                 5800,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.float,
                 "Value of the threshold")
options.register('thresholdsmearing',
                 0.,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.float,
                 "Value of the threshold smearing")
options.register('tofsmearing',
                 0.,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.float,
                 "Value of the tof smearing")
options.register('mode',
                 'scan',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Mode : scan (scan delay values) or emulate (use a random value of delay)") 
options.register('offset',
                 -1.,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.float,
                 "Specific offset value only used in emulate mode, if left to default (-1) will take random value")
options.register('verbose',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Verbose level : 0 (nothing) | 1 (track info) | 2 (BX scan info) | 3 (Firing of the hit detect) | 4 (full detail on algo)")



options.parseArguments()

import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
process = cms.Process('MIX',eras.Phase2)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D21Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.N),
)


# Input source
input_filename = 'file:/nfs/scratch/fynu/fbury/UpgradeCalibrations/Files/CMSSW_106X/SingleMuPt_100_GEN_SIM_1.root',
process.source = cms.Source("PoolSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    fileNames = cms.untracked.vstring(input_filename),
    inputCommands = cms.untracked.vstring(
        'keep *', 
        'drop *_genParticles_*_*', 
        'drop *_genParticlesForJets_*_*', 
        'drop *_kt4GenJets_*_*', 
        'drop *_kt6GenJets_*_*', 
        'drop *_iterativeCone5GenJets_*_*', 
        'drop *_ak4GenJets_*_*', 
        'drop *_ak7GenJets_*_*', 
        'drop *_ak8GenJets_*_*', 
        'drop *_ak4GenJetsNoNu_*_*', 
        'drop *_ak8GenJetsNoNu_*_*', 
        'drop *_genCandidatesForMET_*_*', 
        'drop *_genParticlesForMETAllVisible_*_*', 
        'drop *_genMetCalo_*_*', 
        'drop *_genMetCaloAndNonPrompt_*_*', 
        'drop *_genMetTrue_*_*', 
        'drop *_genMetIC5GenJs_*_*'
    ),
    secondaryFileNames = cms.untracked.vstring()
)


process.options = cms.untracked.PSet()

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
    oncePerEventMode = cms.untracked.bool(True),
    ignoreTotal = cms.untracked.int32(1)
)

# Production Info 
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step2 nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)


# Customization of Mixing Module
process.mix.input.nbPileupEvents.averageNumber = cms.double(200.0)
process.mix.bunchspace = cms.int32(25)
process.mix.minBunch = cms.int32(0)
process.mix.maxBunch = cms.int32(0)
process.mix.input.fileNames = cms.untracked.vstring([
   'file:/nfs/scratch/fynu/fbury/UpgradeCalibrations/Files/CMSSW_106X/SingleMuPt_100_GEN_SIM_1.root',
   'file:/nfs/scratch/fynu/fbury/UpgradeCalibrations/Files/CMSSW_106X/SingleMuPt_100_GEN_SIM_2.root',
   'file:/nfs/scratch/fynu/fbury/UpgradeCalibrations/Files/CMSSW_106X/SingleMuPt_100_GEN_SIM_3.root',
   'file:/nfs/scratch/fynu/fbury/UpgradeCalibrations/Files/CMSSW_106X/SingleMuPt_100_GEN_SIM_4.root',
   'file:/nfs/scratch/fynu/fbury/UpgradeCalibrations/Files/CMSSW_106X/SingleMuPt_100_GEN_SIM_5.root',
   'file:/nfs/scratch/fynu/fbury/UpgradeCalibrations/Files/CMSSW_106X/SingleMuPt_100_GEN_SIM_6.root',
   'file:/nfs/scratch/fynu/fbury/UpgradeCalibrations/Files/CMSSW_106X/SingleMuPt_100_GEN_SIM_7.root',
   'file:/nfs/scratch/fynu/fbury/UpgradeCalibrations/Files/CMSSW_106X/SingleMuPt_100_GEN_SIM_8.root',
   'file:/nfs/scratch/fynu/fbury/UpgradeCalibrations/Files/CMSSW_106X/SingleMuPt_100_GEN_SIM_9.root',
   'file:/nfs/scratch/fynu/fbury/UpgradeCalibrations/Files/CMSSW_106X/SingleMuPt_100_GEN_SIM_10.root',
   'file:/nfs/scratch/fynu/fbury/UpgradeCalibrations/Files/CMSSW_106X/SingleMuPt_100_GEN_SIM_11.root',
   'file:/nfs/scratch/fynu/fbury/UpgradeCalibrations/Files/CMSSW_106X/SingleMuPt_100_GEN_SIM_12.root',
   'file:/nfs/scratch/fynu/fbury/UpgradeCalibrations/Files/CMSSW_106X/SingleMuPt_100_GEN_SIM_13.root',
   'file:/nfs/scratch/fynu/fbury/UpgradeCalibrations/Files/CMSSW_106X/SingleMuPt_100_GEN_SIM_14.root',
   'file:/nfs/scratch/fynu/fbury/UpgradeCalibrations/Files/CMSSW_106X/SingleMuPt_100_GEN_SIM_15.root',
   'file:/nfs/scratch/fynu/fbury/UpgradeCalibrations/Files/CMSSW_106X/SingleMuPt_100_GEN_SIM_16.root',
   'file:/nfs/scratch/fynu/fbury/UpgradeCalibrations/Files/CMSSW_106X/SingleMuPt_100_GEN_SIM_17.root',
   'file:/nfs/scratch/fynu/fbury/UpgradeCalibrations/Files/CMSSW_106X/SingleMuPt_100_GEN_SIM_18.root',
   'file:/nfs/scratch/fynu/fbury/UpgradeCalibrations/Files/CMSSW_106X/SingleMuPt_100_GEN_SIM_19.root',
   'file:/nfs/scratch/fynu/fbury/UpgradeCalibrations/Files/CMSSW_106X/SingleMuPt_100_GEN_SIM_20.root'])

process.mix.mixObjects.mixSH.crossingFrames.extend([
        'TrackerHitsPixelBarrelHighTof', 
        'TrackerHitsPixelBarrelLowTof', 
        'TrackerHitsPixelEndcapHighTof', 
        'TrackerHitsPixelEndcapLowTof', 
        'TrackerHitsTECHighTof', 
        'TrackerHitsTECLowTof', 
        'TrackerHitsTIBHighTof', 
        'TrackerHitsTIBLowTof', 
        'TrackerHitsTIDHighTof', 
        'TrackerHitsTIDLowTof', 
        'TrackerHitsTOBHighTof',
        'TrackerHitsTOBLowTof',
])
process.mix.digitizers.mergedtruth.select.signalOnlyTP = cms.bool(False)

process.dump = cms.EDAnalyzer("EventContentAnalyzer")


# Output definition
import random
if options.mode == 'scan':
    filename = 'BXHistScan_N_{:d}_threshold_{:d}_thresholdsmearing_{:0.1f}_tofsmearing_{:0.1f}_raw'.format(
                    options.N,
                    int(options.threshold),
                    options.thresholdsmearing,
                    options.tofsmearing)
elif options.mode == 'emulate':
    if options.offset == -1.:
        offset_emulate = round(random.random()*50,2)
    else:
        offset_emulate = round(options.offset,2)
        
    filename = 'BXHistEmulateDelay_{:0.2f}_N_{:d}_threshold{:d}_thresholdsmearing_{:0.1f}_tofsmearing_{:0.1f}_raw'.format(
                    offset_emulate,
                    options.N,
                    int(options.threshold),
                    options.thresholdsmearing,
                    options.tofsmearing)
else:
    raise RuntimeError("mode {} argument not understood".format(options.mode))

filename = filename.replace('.','p')+'.root'
print ("Producing file %s"%filename)

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),  
    fileName = cms.untracked.string(filename),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)


process.load('SimTracker.SiPhase2TimingCalibration.Phase2TrackerBXHistogram_cfi')
process.digiana_seq = cms.Sequence(process.timeCalib)

### Modify Hit parameters 
process.timeCalib.ThresholdInElectrons_Barrel = cms.double(options.threshold)
process.timeCalib.ThresholdInElectrons_Endcap = cms.double(options.threshold)
process.timeCalib.ThresholdSmearing_Barrel = cms.double(options.thresholdsmearing)
process.timeCalib.ThresholdSmearing_Endcap = cms.double(options.thresholdsmearing)
process.timeCalib.TOFSmearing = cms.double(options.tofsmearing)
process.timeCalib.UseMixing = cms.bool(True)
process.timeCalib.Mode = cms.string(options.mode)
process.timeCalib.VerbosityLevel = cms.int32(options.verbose)
if options.mode == 'emulate':
    process.timeCalib.OffsetEmulate = cms.double(offset_emulate)

process.load('IOMC.RandomEngine.IOMC_cff')
process.RandomNumberGeneratorService.generator.initialSeed  = random.randrange(1,10e07)
process.RandomNumberGeneratorService.VtxSmeared.initialSeed = random.randrange(1,10e07)
process.RandomNumberGeneratorService.g4SimHits.initialSeed  = random.randrange(1,10e07)
setattr(process.RandomNumberGeneratorService,'timeCalib',cms.PSet(
                    initialSeed = cms.untracked.uint32(random.randrange(1,10e07)),
                    engineName  = cms.untracked.string('TRandom3'))  
)


process.load('DQMServices.Components.DQMEventInfo_cfi')
process.dqmEnv.subSystemFolder = cms.untracked.string('Ph2TkTB')

process.dqm_comm = cms.Sequence(process.dqmEnv)

# Path and EndPath definitions
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)


process.dqm_step =  cms.Path(process.digiana_seq * process.dqm_comm )

process.mix_step = cms.Path(process.mix)
# Schedule definition
process.schedule = cms.Schedule(
    process.mix_step,
    process.dqm_step,
    process.FEVTDEBUGoutput_step,
    process.endjob_step
)

