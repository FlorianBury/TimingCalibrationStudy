from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')

options.register('N',
                 100,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Number of events to be processed")
options.register('pt',
                 2.,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.float,
                 "Value of the pt cut")
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
options.register('subdet',
                 'ALL',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Part of the OT to use for the hits : ALL | BLL | BHL | BLH | BHH | ELL | EHL | ELH | EHH")
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
#process = cms.Process('MIX',eras.Phase2)
process = cms.Process('MIX',eras.Phase2C17I13M9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D88Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21','')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.N),
)


# Input source
process.source = cms.Source("PoolSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    fileNames = cms.untracked.vstring(
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/00a6ee79-da73-4cf3-beb6-4915933d2baf.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/00d31a1e-da64-44ca-bf3d-edbfb23be061.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/0151b925-28d5-4bbd-a05f-2d163d8dfbf4.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/04488096-cf82-4b5a-ba75-95c0dd7024bb.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/0fc60f0a-09c4-46bf-8a67-55d1111ecba7.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/12d7f9ab-ff27-4b7c-b5a9-efc03cc594c6.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/14d10810-db85-45bd-944c-9a516c82269d.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/1c2bb522-d885-406b-847a-07b2ef427847.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/1f0056ce-5ccd-43ee-b547-3fcbefb6799c.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/2882e082-ceea-4e03-af19-433f71ed860b.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/29b73b0f-6bf9-49d3-bf78-0b9c8a1f7985.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/2c6ece2c-93e7-44d8-9e16-a69a7c94e069.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/2f1e9846-1a61-4c7b-8165-9e3291f66f50.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/324feadf-f877-4d49-b119-2ebdcb6fd2ee.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/35c7859e-0919-4aeb-9726-ea4616d6159e.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/3722c427-89ac-4f39-9265-95a303e477c8.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/47ab1cd6-2f52-4948-aeac-4256fbfa8e21.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/518e0259-fb8e-4c1d-b894-a379e0447db4.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/528014ba-f803-4067-a831-88f784ee3979.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/53816c6f-ee39-4967-895d-8ecd074f80b7.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/587b7d74-c940-406c-b25a-382d0e2c915e.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/5c63b089-a807-42f2-8842-dd7c30567d83.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/5cfcf255-adbb-44e3-97f4-86a19d91c407.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/60427347-d3e3-40eb-984b-9e70c48ebdb8.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/6836a877-a2a9-4ea4-a7af-540ac616eda2.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/7c921071-f927-446b-b3d1-dbc2743ce0ee.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/8231f1c7-f3ad-4f8a-82fa-acfbacccce86.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/89fa055a-16b2-460a-b76c-f0250b3ab12d.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/8b5649d0-d0d0-4860-985b-8235b0b5ac53.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/91c85fbb-f61f-4bb4-9446-e5b02471360e.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/937ef1c9-8088-4623-992d-8d8665d9ec51.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/97557031-d017-4f8c-bceb-697dffb97d65.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/99b7cdf2-1568-4197-9cec-c51b2311ed8d.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/9c9d076a-62a4-4f27-b4a7-0208d27d6e74.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/9ca19975-e2f7-4649-9525-1e47382e4e97.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/a0fbaf4d-d7ae-4172-997c-347055fc7107.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/aa8679f3-39c8-45d7-a0b4-e7e4315fad60.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/abac8e4e-5985-46e4-b3fd-10fac267c568.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/adda46d6-0b92-45cc-9731-7a163cfef1f2.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/bf427181-f567-4a11-a6df-65dabd65332b.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/c3c53074-c86a-4cd0-8d2d-d3d1c7287b54.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/d61c911f-cf14-41ad-909f-58a4e0f25c7c.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/e7d8ce40-ac11-418a-964c-46b491c022d0.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/e9672321-3a6e-4e39-831a-c7c60f587509.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/eb5c3b98-ffe3-4e9c-a1b3-2fa7a308d289.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/f0fb66dd-0504-4f4b-bb66-4eb35c093c72.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/f5f3d4a5-a40a-42e8-abf2-9f9fa1627aad.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/f6aa2409-b476-4788-b1c4-4eea6287a743.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/f90f4a1e-8ff7-4d55-9875-31e3d3ff0397.root',
        'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValTTbar_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/fd2243be-880a-426d-9e78-41fd9b5c1f40.root',
    ),
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
        'drop *_genMetIC5GenJs_*_*',
    ),
    secondaryFileNames = cms.untracked.vstring()
)


process.options = cms.untracked.PSet(
)

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
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/02d8d7a3-95da-49d8-8444-68de41aab483.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/03daf578-0b8a-4f69-b248-0959fd2c42af.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/102c44a4-918a-4f03-adfa-87878b830d88.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/112fb487-7f18-4ccc-acb2-9be6c9e5250f.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/129d181d-40bb-455a-868b-eef467d891d3.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/1767ef65-e1cc-48b3-a183-d3cd99686efe.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/18c2ef6f-594a-493b-b640-4a2895dd5975.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/19055773-0e9d-4ec1-a843-cda0bb9b334a.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/1ddd6611-d838-4665-b197-5813e45bcf96.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/1f499606-55ea-497b-82b4-f55106db0cca.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/21682149-6e58-4291-9c93-9338d04bd20e.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/227afd35-465b-4817-9499-a45a2aee6efd.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/260834a3-ff85-4d34-bdf7-84e3358be8ca.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/286ea1be-2010-4003-94b2-f5b17ff0b02e.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/2a2714fc-aec4-4d31-8fdd-8b290b21779a.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/2a8e35ca-4430-4f51-ad64-25229c44572b.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/2abf9f47-4527-4186-af59-0519231f1d5b.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/2b69d9f8-71b0-45e1-bcd5-249e54046e93.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/2c075cff-b6bf-45bb-bcab-52e8a6b73000.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/2e2d49b7-7255-4ec2-8bb2-5affafa030a3.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/2f721d7c-6bf5-4679-af92-b730603276e2.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/32cb2fe5-cf5a-4493-9b0d-c4eb7a483698.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/35202f14-9c44-43f7-82ae-96e0e25b2472.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/3e5852c6-4423-4eb9-837f-ca0d04b6dddb.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/42db1128-4ff7-4f13-be1d-c531c89b09c4.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/44cc1705-0090-4a0f-906e-261f982109ec.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/4ade1f74-809e-4e12-8376-8c180a1da0c2.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/4ed6441b-4866-4912-a1a4-ac33b700f17b.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/50004d3c-6d56-4b38-8ad4-7917c9629f9c.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/53fcb3ff-299a-4910-8eee-d97357d5b617.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/58ac0b23-51b7-44ce-9d90-c8f46c5e1cf4.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/602ba932-5041-4675-ae9c-ca7f78bcc50a.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/60f6423c-2ab4-4c37-841d-7d8dee78a408.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/62fd357c-d323-435a-863f-f923a3208730.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/641711e6-956b-4ceb-9d54-6beb1a15660b.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/6d3af75b-3d19-49a8-8173-3d2f8b9ba079.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/6fcfab21-c996-4503-b331-53161da2b701.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/72004f26-4adc-4447-8390-e9606f7acf73.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/750d1ef5-d177-4c34-8c89-7e1801a3004c.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/77b5f091-6819-4d36-9035-c53b6ffb57a4.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/79dd6761-326a-45dc-9955-fe595d9778ed.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/7ab02707-3c5b-41a4-ab21-2eb7515e1e09.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/7c436661-ad26-486c-bfae-6a72785cb16c.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/80346df5-d76e-437f-a8dc-40657b2ffe52.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/8069244c-95a5-49d6-ac4c-fc2162fa55ec.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/828b04ce-2aac-4af3-9a1e-513d10f5181e.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/83089dc5-b242-4ab3-95ae-c2ac7869c27d.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/84a18f03-911c-49df-ba8e-03e21a2785d8.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/85c7cdc2-5397-445f-8bb1-a477973a4dde.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/86c68b6a-3b68-4607-a6c7-9ef100408bf7.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/88d9f1c1-77c8-4179-bf7b-d9d367c31fba.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/89fbca33-962a-4801-86cd-590cfb11f94a.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/8ded18b1-707c-411a-b6a0-46b0535f8bff.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/8e5b07a8-3ba1-46e3-a8cb-51174405cf35.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/921c0e2b-5b55-4e54-bb6f-db2e707fa90b.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/959fdc30-494a-4881-8394-bd7dfac7d425.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/a0d522c9-d9da-4482-8648-bc570c436cbf.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/a119aed5-3d10-4bef-a666-548169251df4.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/a1fe1c74-9ca5-4364-9eee-d6cc7d95c456.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/b2a1d533-b477-47ec-93a2-9201fc85f591.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/b925dad2-32b6-4601-b15a-cada8db7d3a8.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/bffdd57d-b0c6-4828-a99d-909c003f2fe6.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/c5d23120-b272-4b57-a2e8-fce826419cbf.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/c780fda8-4259-46e7-9a0c-784148f7e5fc.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/c84b4133-cd4a-40df-930a-d36158d4fd20.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/d46684bf-e91b-4d81-84a3-7b1f1b19a8be.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/dcf164ec-6d33-400f-9627-c9e1f2a11bfa.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/de159e7b-c110-4bca-90b5-f4123be118fe.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/e2736381-5f99-4b34-ac69-5fd92fab685d.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/e67e8559-41fc-43c9-b11a-7d2cee804748.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/e88bec4c-5837-4988-92b5-ea491403e0c7.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/e931e4b8-9b8a-4fd8-84bb-1936c371f247.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/ed80409d-e80d-47a2-9233-caf0d6fb0842.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/eedf3fd2-afa4-49c8-8494-28ab6689a67d.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/f14655ac-2db4-44cd-8fbc-708d347b0aa7.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/f1ffcb9b-c19c-4fc9-847f-ed898f3b199b.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/f4103ee3-fe7b-4f92-bd11-4affb42011be.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/f87e4295-7d9d-49fa-ba84-b285c11afd67.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/f9db95d8-cd2d-4375-a76a-c151754b96d6.root',
   'file:/home/ucl/cp3/fbury/scratch/UpgradeCalibrations/Files/CMSSW_12_5_0/RelValMinBias_14TeV/GEN-SIM/125X_mcRun4_realistic_v2_2026D88noPU-v1/ff5397bd-28bc-471b-b728-ed59bed93e20.root',
])



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
    filename = 'BXHistScan_subdet_{:s}_N_{:d}_pt_{:.01f}_threshold_{:d}_thresholdsmearing_{:0.1f}_tofsmearing_{:0.1f}_raw'.format(
                    options.subdet,
                    options.N,
                    options.pt,
                    int(options.threshold),
                    options.thresholdsmearing,
                    options.tofsmearing)
elif options.mode == 'emulate':
    if options.offset == -1.:
        offset_emulate = round(random.random()*50,2)
    else:
        offset_emulate = round(options.offset,2)
    filename = 'BXHistEmulateDelay_{:0.2f}_subdet_{:s}_N_{:d}_pt_{:.01f}_threshold{:d}_thresholdsmearing_{:0.1f}_tofsmearing_{:0.1f}_raw'.format(
                    offset_emulate,
                    options.subdet,
                    options.N,
                    options.pt,
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
    fileName = cms.untracked.string('file:'+filename),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)


process.load('SimTracker.SiPhase2TimingCalibration.Phase2TrackerBXHistogram_cfi')
process.digiana_seq = cms.Sequence(process.timeCalib)

### Modify Hit parameters
process.timeCalib.PTCut = cms.double(options.pt)
process.timeCalib.ThresholdInElectrons_Barrel = cms.double(options.threshold)
process.timeCalib.ThresholdInElectrons_Endcap = cms.double(options.threshold)
process.timeCalib.ThresholdSmearing_Barrel = cms.double(options.thresholdsmearing)
process.timeCalib.ThresholdSmearing_Endcap = cms.double(options.thresholdsmearing)
process.timeCalib.TOFSmearing = cms.double(options.tofsmearing)
process.timeCalib.UseMixing = cms.bool(True)
process.timeCalib.Mode = cms.string(options.mode)
process.timeCalib.Subdetector = cms.string(options.subdet)
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

