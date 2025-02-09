import FWCore.ParameterSet.Config as cms
import os

from tthAnalysis.NanoAOD.LeptonFakeRate_trigger_cfi import *
from tthAnalysis.HiggsToTauTau.configs.recommendedMEtFilters_cfi import *
from tthAnalysis.HiggsToTauTau.configs.EvtYieldHistManager_cfi import *
from tthAnalysis.HiggsToTauTau.analysisSettings import *

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames = cms.vstring(),
    maxEvents = cms.int32(-1),
    outputEvery = cms.uint32(100000)
)

process.fwliteOutput = cms.PSet(
    fileName = cms.string('')
)

process.analyze_LeptonFakeRate = cms.PSet(
    treeName = cms.string('Events'),

    process = cms.string(''),
    #histogramDir = cms.string(''), ## Not Relevant to the LeptonFakeRate code
    era = cms.string(''),

    triggers_1e = cms.vstring(),
    use_triggers_1e = cms.bool(True),
    triggers_1mu = cms.vstring(),
    use_triggers_1mu = cms.bool(True),
    triggers_2e = cms.vstring(),
    use_triggers_2e = cms.bool(True),
    triggers_2mu = cms.vstring(),
    use_triggers_2mu = cms.bool(True),
    ## ---- Not Relevant to the LeptonFakeRate code-------------------------------------------##
    ##--- (since not present in tthAnalysis/NanoAOD/python/LeptonFakeRate_trigger_cfi.py) ----## 
    #triggers_1e1mu = cms.vstring(),      
    #use_triggers_1e1mu = cms.bool(True), ## Not Relevant to the LeptonFakeRate code 
    ## ---------------------------------------------------------------------------------------##
    triggers_mu_cfg = cms.VPSet(),
    triggers_e_cfg = cms.VPSet(),

    absEtaBins_e = cms.vdouble(0., 1.479, 9.9),
    ptBins_e = cms.vdouble(15., 20., 30., 45., 65., 100000.),
    absEtaBins_mu = cms.vdouble(0., 1.479, 9.9),
    ptBins_mu = cms.vdouble(10., 15., 20., 30., 45., 65., 100000.),

    minConePt_global_e = cms.double(10),
    minRecoPt_global_e = cms.double(7),
    minConePt_global_mu = cms.double(10),
    minRecoPt_global_mu = cms.double(5),

    lep_mva_cut_mu = cms.double(1.),
    lep_mva_cut_e  = cms.double(1.),
    lep_mva_wp = cms.string(''),
    METScaleSyst   = cms.double(0.10), ## MET Syst set to 10%

    isMC = cms.bool(True),
    central_or_shift = cms.string(''),
    lumiScale = cms.VPSet(),
    ref_genWeight = cms.double(0.),

    apply_genWeight = cms.bool(True),
    apply_DYMCReweighting = cms.bool(False),
    apply_DYMCNormScaleFactors = cms.bool(False),
    apply_topPtReweighting = cms.string(''),
    apply_l1PreFireWeight = cms.bool(True),
    apply_met_filters = cms.bool(True),
    min_PV_ndof = cms.double(100.),
    cfgMEtFilter = cms.PSet(),
    fillGenEvtHistograms = cms.bool(True),
    cfgEvtYieldHistManager = cms.PSet(),
    triggerWhiteList = cms.PSet(),

    branchName_electrons = cms.string('Electron'),
    branchName_muons = cms.string('Muon'),
    branchName_hadTaus = cms.string('Tau'),
    branchName_jets = cms.string('Jet'),
    branchName_met = cms.string('MET'),
    branchName_genmet = cms.string('GenMET'),
    branchName_vertex = cms.string('PV'),

    branchName_genTauLeptons = cms.string('GenTau'),
    branchName_genLeptons = cms.string('GenLep'),
    branchName_genHadTaus = cms.string('GenVisTau'),
    branchName_genPhotons = cms.string('GenPhoton'),
    branchName_genJets = cms.string('GenJet'),
    branchName_muonGenMatch     = cms.string('MuonGenMatch'),
    branchName_electronGenMatch = cms.string('ElectronGenMatch'),
    branchName_jetGenMatch      = cms.string('JetGenMatch'),
    redoGenMatching = cms.bool(True),
    genMatchingByIndex = cms.bool(True),
    jetCleaningByIndex = cms.bool(True),

    selEventsFileName_input = cms.string(''),
    selEventsFileName_output = cms.string(''),
    useNonNominal = cms.bool(False), ## Added from HH 3l channel
    isDEBUG = cms.bool(False),
    applyMETFilters = cms.bool(True),
    hasLHE = cms.bool(True),
    useObjectMultiplicity = cms.bool(False),

    evtWeight = cms.PSet(
        apply = cms.bool(False),
        histogramFile = cms.string(''),
        histogramName = cms.string(''),
        branchNameXaxis = cms.string(''),
        branchNameYaxis = cms.string(''),
        branchTypeXaxis = cms.string(''),
        branchTypeYaxis = cms.string(''),
    ),
    tHweights = cms.VPSet(),
    fillNtuple = cms.bool(False),  ## Boolean handle to include Ntuples for optmizing Lepton I.D. cuts for MC Closure sidebands (2lss, TTHadronic)
    enable_MC_Closure_sidebands = cms.bool(False) ## Boolean handle to control inclusion of MC Closure sidebands (2lss, TTHadronic)
)
