import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames = cms.vstring(),
    skipEvents = cms.uint32(0),
    maxEvents = cms.int32(-1),
    outputEvery = cms.uint32(1000)
)

process.fwliteOutput = cms.PSet(
    fileName = cms.string('')
)

process.addMEM_3l_1tau = cms.PSet(
    treeName = cms.string('Events'),

    era = cms.string(''),
    isMC = cms.bool(True),

    leptonSelection = cms.string(''),
    hadTauSelection = cms.string(''),

    branchName_electrons = cms.string('Electron'),
    branchName_muons = cms.string('Muon'),
    branchName_hadTaus = cms.string('Tau'),
    branchName_jets = cms.string('Jet'),
    branchName_met = cms.string('MET'),
    branchName_vertex = cms.string('PV'),

    copy_all_branches = cms.bool(True),
    copy_histograms = cms.vstring('Count.*'),

    selEventsFileName_input = cms.string(''),
    isDEBUG = cms.bool(False),
    readGenObjects = cms.bool(True),
    jetCleaningByIndex = cms.bool(True),

    central_or_shift = cms.vstring(),
    useNonNominal = cms.bool(False),
    dryRun = cms.bool(False),
)
