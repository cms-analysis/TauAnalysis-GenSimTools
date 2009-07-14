import FWCore.ParameterSet.Config as cms

process = cms.Process("PROD")

# The number of events to be processed.
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )
    
# For valgrind studies
# process.ProfilerService = cms.Service("ProfilerService",
#    lastEvent = cms.untracked.int32(13),
#    firstEvent = cms.untracked.int32(3),
#    paths = cms.untracked.vstring('p1')
#)

process.load("SimGeneral/HepPDTESSource/pdt_cfi")

# Include the RandomNumberGeneratorService definition
process.load("FastSimulation/Configuration/RandomServiceInitialization_cff")


# Generate multijet events with different ptHAT bins
process.load("FastSimulation.Configuration.QCDpt15_cfi")

# Famos sequences (Frontier conditions)
process.load("FastSimulation/Configuration/CommonInputs_cff")
process.load("FastSimulation/Configuration/FamosSequences_cff")
process.GlobalTag.globaltag = "IDEAL_V12::All"

# Parametrized magnetic field (new mapping, 4.0 and 3.8T)
#process.load("Configuration.StandardSequences.MagneticField_40T_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.VolumeBasedMagneticFieldESProducer.useParametrizedTrackerField = True

# If you want to turn on/off pile-up
# process.famosPileUp.PileUpSimulator.averageNumber = 5.0    
# You may not want to simulate everything for your study
process.famosSimHits.SimulateCalorimetry = True
process.famosSimHits.SimulateTracking = True
process.famosSimHits.SimulateMuons = True

# Generator level filter
process.genLeadTrackFilter = cms.EDFilter('GenLeadTrackFilter',
  HepMCProduct             = cms.InputTag("source"),
  GenLeadTrackPt           = cms.double(12),
  GenEta                   = cms.double(2.5)
)
# GenLeadTrackFilter

# Produce Tracks and Clusters
process.p1 = cms.Path(process.genLeadTrackFilter*process.famosWithTracksAndEcalClusters)

# To write out events (not need: FastSimulation _is_ fast!)
process.o1 = cms.OutputModule("PoolOutputModule",
    cms.PSet(
        outputCommands = cms.untracked.vstring("keep *", "drop *_mix_*_*")
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p1')
    ),
    fileName = cms.untracked.string("/home/sdas/FastSimProduction/QCD_PtTrack15_FASTSIM_5000.root")
    #fileName = cms.untracked.string("QCD_PtTrack15_FASTSIM.root")
)

process.outpath = cms.EndPath(process.o1)

# Keep the logging output to a nice level #
# process.Timing =  cms.Service("Timing")
# process.load("FWCore/MessageService/MessageLogger_cfi")
# process.MessageLogger.destinations = cms.untracked.vstring("pyDetailedInfo.txt")

# Make the job crash in case of missing product
process.options = cms.untracked.PSet( Rethrow = cms.untracked.vstring('ProductNotFound') )
