import logging
import re

from tthAnalysis.HiggsToTauTau.configs.analyzeConfig import *
from tthAnalysis.HiggsToTauTau.jobTools import create_if_not_exists
from tthAnalysis.HiggsToTauTau.analysisTools import initDict, getKey, create_cfg, createFile, generateInputFileList, is_dymc_reweighting

def get_lepton_and_hadTau_selection_and_frWeight(lepton_and_hadTau_selection, lepton_and_hadTau_frWeight):
  lepton_and_hadTau_selection_and_frWeight = lepton_and_hadTau_selection
  if lepton_and_hadTau_selection.startswith("Fakeable"):
    if lepton_and_hadTau_frWeight == "enabled":
      lepton_and_hadTau_selection_and_frWeight += "_wFakeRateWeights"
    elif lepton_and_hadTau_frWeight == "disabled":
      lepton_and_hadTau_selection_and_frWeight += "_woFakeRateWeights"
  lepton_and_hadTau_selection_and_frWeight = lepton_and_hadTau_selection_and_frWeight.replace("|", "_")
  return lepton_and_hadTau_selection_and_frWeight

def getHistogramDirList(lepton_selection, hadTau_selection, lepton_and_hadTau_frWeight, lepton_charge_selection, hadTau_charge_selection, chargeSumSelection, subcategories):
  histogramDirListL = []
  for cc, cat in enumerate(subcategories) :
      histogramDir = "%s_%s" % (cat, hadTau_selection)
      histogramDir += "%s_%s" % (cat, lepton_selection)
      if lepton_charge_selection != "disabled":
        histogramDir += "_lep%s" % lepton_charge_selection
      if hadTau_charge_selection != "disabled":
        histogramDir += "_hadTau%s" % hadTau_charge_selection
      if lepton_selection.find("Fakeable") != -1 or hadTau_selection.find("Fakeable") != -1:
        if lepton_and_hadTau_frWeight == "enabled":
          histogramDir += "_wFakeRateWeights"
        elif lepton_and_hadTau_frWeight == "disabled":
          histogramDir += "_woFakeRateWeights"
      histogramDirListL+=[histogramDir]
  return histogramDirListL

class analyzeConfig_inclusive_toNN(analyzeConfig):
  """Configuration metadata needed to run analysis in a single go.

  Sets up a folder structure by defining full path names; no directory creation is delegated here.

  Args specific to analyzeConfig_inclusive_toNN:
    hadTau_selection: either `Tight`, `Loose` or `Fakeable`
    hadTau_charge_selection: either `OS` or `SS` (opposite-sign or same-sign)

  See $CMSSW_BASE/src/tthAnalysis/HiggsToTauTau/python/analyzeConfig.py
  for documentation of further Args.

  """
  def __init__(self,
        configDir,
        outputDir,
        executable_analyze,
        cfgFile_analyze,
        samples,
        lep_mva_wp,
        lepton_charge_selections,
        hadTau_selection,
        hadTau_charge_selections,
        applyFakeRateWeights,
        chargeSumSelections,
        central_or_shifts,
        max_files_per_job,
        era,
        use_lumi,
        lumi,
        check_output_files,
        running_method,
        num_parallel_jobs,
        executable_addBackgrounds,
        executable_addBackgroundJetToTauFakes,
        histograms_to_fit,
        select_rle_output         = False,
        executable_prep_dcard     = "prepareDatacards",
        executable_add_syst_dcard = "addSystDatacards",
        verbose                   = False,
        dry_run                   = False,
        do_sync                   = False,
        isDebug                   = False,
        rle_select                = '',
        use_nonnominal            = False,
        hlt_filter                = False,
        use_home                  = True,
      ):
    analyzeConfig.__init__(self,
      configDir                 = configDir,
      outputDir                 = outputDir,
      executable_analyze        = executable_analyze,
      channel                   = "inclusive_toNN",
      subcategories             = ["all"],
      samples                   = samples,
      lep_mva_wp                = lep_mva_wp,
      central_or_shifts         = central_or_shifts,
      max_files_per_job         = max_files_per_job,
      era                       = era,
      use_lumi                  = use_lumi,
      lumi                      = lumi,
      check_output_files        = check_output_files,
      running_method            = running_method,
      num_parallel_jobs         = num_parallel_jobs,
      histograms_to_fit         = histograms_to_fit,
      triggers                  = [ '1e', '1mu', '2e', '2mu', '1e1mu', '2tau', '3e', '3mu', '1e2mu', '2e1mu', '1e1tau', '1mu1tau' ],
      executable_prep_dcard     = executable_prep_dcard,
      executable_add_syst_dcard = executable_add_syst_dcard,
      verbose                   = verbose,
      dry_run                   = dry_run,
      do_sync                   = do_sync,
      isDebug                   = isDebug,
      use_home                  = use_home,
    )

    self.lepton_and_hadTau_selections = [ "Tight", "Fakeable" ]
    self.lepton_and_hadTau_frWeights = [ "enabled", "disabled" ]
    self.hadTau_selection_part2 = hadTau_selection
    self.applyFakeRateWeights = applyFakeRateWeights
    self.lepton_charge_selections = lepton_charge_selections
    self.hadTau_charge_selections = hadTau_charge_selections
    run_mcClosure = 'central' not in self.central_or_shifts or len(central_or_shifts) > 1 or self.do_sync
    if self.era != '2017':
      logging.warning('mcClosure for lepton FR not possible for era %s' % self.era)
      run_mcClosure = False
    if run_mcClosure:
      # Run MC closure jobs only if the analysis is run w/ (at least some) systematic uncertainties
      #self.lepton_and_hadTau_selections.extend([ "Fakeable_mcClosure_all" ]) #TODO
      pass

    self.lepton_genMatches = [ "2l0g0j", "1l1g0j", "1l0g1j", "0l2g0j", "0l1g1j", "0l0g2j" ]
    self.hadTau_genMatches = [ "2t0e0m0j", "1t1e0m0j", "1t0e1m0j", "1t0e0m1j", "0t2e0m0j", "0t1e1m0j", "0t1e0m1j", "0t0e2m0j", "0t0e1m1j", "0t0e0m2j" ]

    self.apply_leptonGenMatching = None
    self.apply_hadTauGenMatching = None
    self.lepton_and_hadTau_genMatches_nonfakes = []
    self.lepton_and_hadTau_genMatches_conversions = []
    self.lepton_and_hadTau_genMatches_fakes = []
    self.lepton_and_hadTau_genMatches_gentau = []
    self.lepton_and_hadTau_genMatches_faketau = []
    if self.applyFakeRateWeights == "4L":
      self.apply_leptonGenMatching = True
      self.apply_hadTauGenMatching = True
      for lepton_genMatch in self.lepton_genMatches:
        for hadTau_genMatch in self.hadTau_genMatches:
          lepton_and_hadTau_genMatch = "&".join([ lepton_genMatch, hadTau_genMatch ])
          if lepton_genMatch.endswith("0g0j") and hadTau_genMatch.endswith("0j"):
            self.lepton_and_hadTau_genMatches_nonfakes.append(lepton_and_hadTau_genMatch)
          elif lepton_genMatch.endswith("0j") and hadTau_genMatch.endswith("0j"):
            self.lepton_and_hadTau_genMatches_conversions.append(lepton_and_hadTau_genMatch)
          else:
            self.lepton_and_hadTau_genMatches_fakes.append(lepton_and_hadTau_genMatch)
      if run_mcClosure:
        self.lepton_and_hadTau_selections.extend([ "Fakeable_mcClosure_e", "Fakeable_mcClosure_m", "Fakeable_mcClosure_t" ])
    elif applyFakeRateWeights == "2lepton":
      self.apply_leptonGenMatching = True
      self.apply_hadTauGenMatching = True
      for lepton_genMatch in self.lepton_genMatches:
        for hadTau_genMatch in self.hadTau_genMatches:
          lepton_and_hadTau_genMatch = "&".join([ lepton_genMatch, hadTau_genMatch ])
          if lepton_genMatch.endswith("0g0j"):
            self.lepton_and_hadTau_genMatches_nonfakes.append(lepton_and_hadTau_genMatch)
            if hadTau_genMatch.endswith("0j"):
              self.lepton_and_hadTau_genMatches_gentau.append(lepton_and_hadTau_genMatch)
            else:
              self.lepton_and_hadTau_genMatches_faketau.append(lepton_and_hadTau_genMatch)
          elif lepton_genMatch.endswith("0j"):
            self.lepton_and_hadTau_genMatches_conversions.append(lepton_and_hadTau_genMatch)
          else:
            self.lepton_and_hadTau_genMatches_fakes.append(lepton_and_hadTau_genMatch)
      if run_mcClosure:
        self.lepton_and_hadTau_selections.extend([ "Fakeable_mcClosure_e", "Fakeable_mcClosure_m" ])
    elif applyFakeRateWeights == "2tau":
      self.apply_leptonGenMatching = True
      self.apply_hadTauGenMatching = True
      for lepton_genMatch in self.lepton_genMatches:
        for hadTau_genMatch in self.hadTau_genMatches:
          if lepton_genMatch.find("0g") != -1 and hadTau_genMatch.endswith("0j"):
            self.lepton_and_hadTau_genMatches_nonfakes.append(hadTau_genMatch)
          elif hadTau_genMatch.endswith("0j"):
            self.lepton_and_hadTau_genMatches_conversions.append(hadTau_genMatch)
          else:
            self.lepton_and_hadTau_genMatches_fakes.append(hadTau_genMatch)
      if run_mcClosure:
        self.lepton_and_hadTau_selections.extend([ "Fakeable_mcClosure_t" ])
    else:
      raise ValueError("Invalid Configuration parameter 'applyFakeRateWeights' = %s !!" % applyFakeRateWeights)

    self.chargeSumSelections = chargeSumSelections

    self.executable_addBackgrounds = executable_addBackgrounds
    self.executable_addFakes = executable_addBackgroundJetToTauFakes

    self.nonfake_backgrounds = [ "TT", "TTW", "TTZ", "TTWW", "EWK", "Rares", "tHq", "tHW", "VH" ]

    self.cfgFile_analyze = os.path.join(self.template_dir, cfgFile_analyze)
    self.prep_dcard_processesToCopy = [ "data_obs" ] + self.nonfake_backgrounds + [ "conversions", "fakes_data", "fakes_mc" ]
    histogramDir_prep_dcard_local = []
    for cc, cat in enumerate(self.subcategories) :
        histogramDir_prep_dcard_local+=[self.subcategories[cc]+"_Tight"]
    self.histogramDir_prep_dcard = histogramDir_prep_dcard_local
    self.make_plots_backgrounds = [ "TTW", "TTZ", "TTWW", "EWK", "Rares", "tHq", "tHW" ] + [ "conversions", "fakes_data" ]
    self.cfgFile_make_plots = os.path.join(self.template_dir, "makePlots_inclusive_toNN_cfg.py")
    self.cfgFile_make_plots_mcClosure = os.path.join(self.template_dir, "makePlots_mcClosure_inclusive_toNN_cfg.py") #TODO

    self.select_rle_output = select_rle_output
    self.rle_select = rle_select
    self.use_nonnominal = use_nonnominal
    self.hlt_filter = hlt_filter

    self.isBDTtraining = False

  def set_BDT_training(self, hadTau_selection_relaxed):
    """Run analysis with loose selection criteria for leptons and hadronic taus,
       for the purpose of preparing event list files for BDT training.
    """
    self.lepton_and_hadTau_selections = [ "forBDTtraining" ]
    self.lepton_and_hadTau_frWeights  = [ "disabled" ]
    super(analyzeConfig_inclusive_toNN, self).set_BDT_training(hadTau_selection_relaxed)

  def createCfg_analyze(self, jobOptions, sample_info, lepton_and_hadTau_selection):
    """Create python configuration file for the analyze_inclusive_toNN executable (analysis code)

    Args:
      inputFiles: list of input files (Ntuples)
      outputFile: output file of the job -- a ROOT file containing histogram
      process: either `TT`, `TTW`, `TTZ`, `EWK`, `Rares`, `data_obs`, `ttH_hww`, 'ttH_hzg', 'ttH_hmm', `ttH_hzz` or `ttH_htt`
      is_mc: flag indicating whether job runs on MC (True) or data (False)
      lumi_scale: event weight (= xsection * luminosity / number of events)
      central_or_shift: either 'central' or one of the systematic uncertainties defined in $CMSSW_BASE/src/tthAnalysis/HiggsToTauTau/bin/analyze_inclusive_toNN.cc
    """
    lepton_and_hadTau_frWeight = "disabled" if jobOptions['applyFakeRateWeights'] == "disabled" else "enabled"

    jobOptions['histogramDir'] = getHistogramDirList(
      lepton_and_hadTau_selection, jobOptions['hadTauSelection'], lepton_and_hadTau_frWeight,
      jobOptions['leptonChargeSelection'], jobOptions['hadTauChargeSelection'], jobOptions['chargeSumSelection'], [self.channel]
    )[0]
    if 'mcClosure' in lepton_and_hadTau_selection:
      self.mcClosure_dir['%s_%s_%s' % (lepton_and_hadTau_selection, jobOptions['chargeSumSelection'], jobOptions['hadTauChargeSelection'])] = jobOptions['histogramDir']

    self.set_leptonFakeRateWeightHistogramNames(jobOptions['central_or_shift'], lepton_and_hadTau_selection)
    jobOptions['leptonFakeRateWeight.inputFileName'] = self.leptonFakeRateWeight_inputFile
    jobOptions['leptonFakeRateWeight.histogramName_e'] = self.leptonFakeRateWeight_histogramName_e
    jobOptions['leptonFakeRateWeight.histogramName_mu'] = self.leptonFakeRateWeight_histogramName_mu

    jobOptions['hadTauFakeRateWeight.inputFileName'] = self.hadTauFakeRateWeight_inputFile
    graphName = 'jetToTauFakeRate/%s/$etaBin/jetToTauFakeRate_mc_hadTaus_pt' % self.hadTau_selection_part2
    jobOptions['hadTauFakeRateWeight.lead.graphName'] = graphName
    jobOptions['hadTauFakeRateWeight.sublead.graphName'] = graphName
    fitFunctionName = 'jetToTauFakeRate/%s/$etaBin/fitFunction_data_div_mc_hadTaus_pt' % self.hadTau_selection_part2
    jobOptions['hadTauFakeRateWeight.lead.fitFunctionName'] = fitFunctionName
    jobOptions['hadTauFakeRateWeight.sublead.fitFunctionName'] = fitFunctionName
    if "mcClosure" in lepton_and_hadTau_selection:
      jobOptions['hadTauFakeRateWeight.applyGraph_lead'] = True
      jobOptions['hadTauFakeRateWeight.applyGraph_sublead'] = True
      jobOptions['hadTauFakeRateWeight.applyFitFunction_lead'] = False
      jobOptions['hadTauFakeRateWeight.applyFitFunction_sublead'] = False
      if self.applyFakeRateWeights not in [ "4L", "2tau" ] and not self.isBDTtraining:
        # We want to preserve the same logic as running in SR and applying the FF method only to leptons [*]
        jobOptions['hadTauFakeRateWeight.applyFitFunction_lead'] = True
        jobOptions['hadTauFakeRateWeight.applyFitFunction_sublead'] = True
    if jobOptions['hadTauSelection'].find("Tight") != -1 and self.applyFakeRateWeights not in [ "4L", "2tau" ] and not self.isBDTtraining:
      # [*] SR and applying the FF method only to leptons
      jobOptions['hadTauFakeRateWeight.applyGraph_lead'] = False # FR in MC for the leading tau
      jobOptions['hadTauFakeRateWeight.applyGraph_sublead'] = False
      jobOptions['hadTauFakeRateWeight.applyFitFunction_lead'] = True # data-to-MC SF for the leading tau
      jobOptions['hadTauFakeRateWeight.applyFitFunction_sublead'] = True
      jobOptions['apply_hadTauFakeRateSF'] = True

    lines = super(analyzeConfig_inclusive_toNN, self).createCfg_analyze(jobOptions, sample_info)
    create_cfg(self.cfgFile_analyze, jobOptions['cfgFile_modified'], lines)

  def create(self):
    """Creates all necessary config files and runs the complete analysis workfow -- either locally or on the batch system
    """

    for sample_name, sample_info in self.samples.items():
      if not sample_info["use_it"] or sample_info["sample_category"] in [ "additional_signal_overlap", "background_data_estimate" ]:
        continue
      process_name = sample_info["process_name_specific"]
      for lepton_charge_selection in self.lepton_charge_selections:
        for hadTau_charge_selection in self.hadTau_charge_selections:
          for lepton_and_hadTau_selection in self.lepton_and_hadTau_selections:
            for lepton_and_hadTau_frWeight in self.lepton_and_hadTau_frWeights:
              if lepton_and_hadTau_frWeight == "enabled" and not lepton_and_hadTau_selection.startswith("Fakeable"):
                continue
              lepton_and_hadTau_selection_and_frWeight = get_lepton_and_hadTau_selection_and_frWeight(lepton_and_hadTau_selection, lepton_and_hadTau_frWeight)
              for chargeSumSelection in self.chargeSumSelections:
                key_dir = getKey(process_name, lepton_charge_selection, hadTau_charge_selection, lepton_and_hadTau_selection_and_frWeight, chargeSumSelection)
                lepton_and_hadTau_charge_selection = ""
                if lepton_charge_selection != "disabled":
                  lepton_and_hadTau_charge_selection += "_lep%s" % lepton_charge_selection
                if hadTau_charge_selection != "disabled":
                  lepton_and_hadTau_charge_selection += "_hadTau%s" % hadTau_charge_selection
                for dir_type in [ DKEY_CFGS, DKEY_HIST, DKEY_LOGS, DKEY_ROOT, DKEY_RLES, DKEY_SYNC ]:
                  initDict(self.dirs, [ key_dir, dir_type ])
                  if dir_type in [ DKEY_CFGS, DKEY_LOGS ]:
                    self.dirs[key_dir][dir_type] = os.path.join(self.configDir, dir_type, self.channel,
                      "_".join([ lepton_and_hadTau_selection_and_frWeight + lepton_and_hadTau_charge_selection ]), process_name)
                  else:
                    self.dirs[key_dir][dir_type] = os.path.join(self.outputDir, dir_type, self.channel,
                      "_".join([ lepton_and_hadTau_selection_and_frWeight + lepton_and_hadTau_charge_selection ]), process_name)
    for dir_type in [ DKEY_CFGS, DKEY_SCRIPTS, DKEY_HIST, DKEY_LOGS, DKEY_DCRD, DKEY_PLOT, DKEY_HADD_RT, DKEY_SYNC ]:
      initDict(self.dirs, [ dir_type ])
      if dir_type in [ DKEY_CFGS, DKEY_SCRIPTS, DKEY_LOGS, DKEY_DCRD, DKEY_PLOT, DKEY_HADD_RT ]:
        self.dirs[dir_type] = os.path.join(self.configDir, dir_type, self.channel)
      else:
        self.dirs[dir_type] = os.path.join(self.outputDir, dir_type, self.channel)

    for key in self.dirs.keys():
      if type(self.dirs[key]) == dict:
        for dir_type in self.dirs[key].keys():
          create_if_not_exists(self.dirs[key][dir_type])
      else:
        create_if_not_exists(self.dirs[key])

    inputFileLists = {}
    for sample_name, sample_info in self.samples.items():
      if not sample_info["use_it"] or sample_info["sample_category"] in [ "additional_signal_overlap", "background_data_estimate" ]:
        continue
      logging.info("Checking input files for sample %s" % sample_info["process_name_specific"])
      inputFileLists[sample_name] = generateInputFileList(sample_info, self.max_files_per_job)

    mcClosure_regex = re.compile('Fakeable_mcClosure_(?P<type>m|e|t)_wFakeRateWeights')
    for lepton_charge_selection in self.lepton_charge_selections:
      for hadTau_charge_selection in self.hadTau_charge_selections:
        for lepton_and_hadTau_selection in self.lepton_and_hadTau_selections:
          lepton_selection = lepton_and_hadTau_selection
          hadTau_selection = lepton_and_hadTau_selection
          electron_selection = lepton_selection
          muon_selection = lepton_selection

          if self.applyFakeRateWeights == "2tau":
            lepton_selection = "Tight"
            electron_selection = "Tight"
            muon_selection = "Tight"
          elif self.applyFakeRateWeights == "2lepton":
            hadTau_selection = "Tight"
          hadTau_selection = "|".join([ hadTau_selection, self.hadTau_selection_part2 ])

          if "forBDTtraining" in lepton_and_hadTau_selection :
            electron_selection = "Loose" #"Tight" #
            muon_selection = "Loose" #"Tight" #"Loose"
            hadTau_selection = "Tight|%s" % self.hadTau_selection_relaxed
          elif lepton_and_hadTau_selection == "Fakeable_mcClosure_e":
            electron_selection = "Fakeable"
            muon_selection = "Tight"
            hadTau_selection = "Tight"
            hadTau_selection = "|".join([hadTau_selection, self.hadTau_selection_part2])
          elif lepton_and_hadTau_selection == "Fakeable_mcClosure_m":
            electron_selection = "Tight"
            muon_selection = "Fakeable"
            hadTau_selection = "Tight"
            hadTau_selection = "|".join([hadTau_selection, self.hadTau_selection_part2])
          elif lepton_and_hadTau_selection == "Fakeable_mcClosure_t":
            electron_selection = "Tight"
            muon_selection = "Tight"
            hadTau_selection = "Fakeable"
            hadTau_selection = "|".join([hadTau_selection, self.hadTau_selection_part2])

          for lepton_and_hadTau_frWeight in self.lepton_and_hadTau_frWeights:
            if lepton_and_hadTau_frWeight == "enabled" and not lepton_and_hadTau_selection.startswith("Fakeable"):
              continue
            if lepton_and_hadTau_frWeight == "disabled" and not lepton_and_hadTau_selection in [ "Tight", "forBDTtraining", "forBDTtraining_VHbb" ]:
              continue
            lepton_and_hadTau_selection_and_frWeight = get_lepton_and_hadTau_selection_and_frWeight(lepton_and_hadTau_selection, lepton_and_hadTau_frWeight)

            for chargeSumSelection in self.chargeSumSelections:

              for sample_name, sample_info in self.samples.items():
                if not sample_info["use_it"] or sample_info["sample_category"] in [ "additional_signal_overlap", "background_data_estimate" ]:
                  continue
                process_name = sample_info["process_name_specific"]
                logging.info("Creating configuration files to run '%s' for sample %s" % (self.executable_analyze, process_name))

                sample_category = sample_info["sample_category"]
                is_mc = (sample_info["type"] == "mc")
                is_signal = (sample_category == "signal")

                for central_or_shift in self.central_or_shifts:

                  inputFileList = inputFileLists[sample_name]
                  for jobId in inputFileList.keys():
                    if central_or_shift != "central":
                      isFR_shape_shift = (central_or_shift in systematics.FR_all)
                      if not ((lepton_and_hadTau_selection == "Fakeable" and chargeSumSelection == "OS" and isFR_shape_shift) or
                              (lepton_and_hadTau_selection == "Tight"    and chargeSumSelection == "OS")):
                        continue
                      if not is_mc and not isFR_shape_shift:
                        continue

                    if central_or_shift in systematics.LHE().ttH and sample_category != "signal":
                      continue
                    if central_or_shift in systematics.LHE().ttW and sample_category != "TTW":
                      continue
                    if central_or_shift in systematics.LHE().ttZ and sample_category != "TTZ":
                      continue
                    if central_or_shift in systematics.DYMCReweighting and not is_dymc_reweighting(sample_name):
                      continue

                    logging.info(" ... for '%s' and systematic uncertainty option '%s'" % (lepton_and_hadTau_selection_and_frWeight, central_or_shift))

                    # build config files for executing analysis code
                    key_dir = getKey(process_name, lepton_charge_selection, hadTau_charge_selection, lepton_and_hadTau_selection_and_frWeight, chargeSumSelection)
                    key_analyze_job = getKey(process_name, lepton_charge_selection, hadTau_charge_selection,
                      lepton_and_hadTau_selection_and_frWeight, chargeSumSelection, central_or_shift, jobId)
                    ntupleFiles = inputFileList[jobId]
                    if len(ntupleFiles) == 0:
                      logging.warning("No input ntuples for %s --> skipping job !!" % (key_analyze_job))
                      continue

                    syncOutput = ''
                    syncTree = ''
                    syncRequireGenMatching = True
                    mcClosure_match = mcClosure_regex.match(lepton_and_hadTau_selection_and_frWeight)
                    if self.do_sync:
                      if chargeSumSelection != 'OS':
                        continue

                    if syncTree and central_or_shift != "central":
                      syncTree = os.path.join(central_or_shift, syncTree)

                    syncRLE = ''
                    if self.do_sync and self.rle_select:
                      syncRLE = self.rle_select % syncTree
                      if not os.path.isfile(syncRLE):
                        logging.warning("Input RLE file for the sync is missing: %s; skipping the job" % syncRLE)
                        continue

                    if syncOutput:
                      self.inputFiles_sync['sync'].append(syncOutput)

                    cfg_key = getKey(
                       self.channel, process_name, lepton_charge_selection, hadTau_charge_selection,
                       lepton_and_hadTau_selection_and_frWeight, chargeSumSelection, central_or_shift,
                       jobId
                    )
                    cfgFile_modified_path = os.path.join(self.dirs[key_dir][DKEY_CFGS], "analyze_%s_cfg.py" % cfg_key)
                    histogramFile_path = os.path.join(self.dirs[key_dir][DKEY_HIST], "%s.root" % key_analyze_job)
                    logFile_path = os.path.join(self.dirs[key_dir][DKEY_LOGS], "analyze_%s.log" % cfg_key)
                    rleOutputFile_path = os.path.join(self.dirs[key_dir][DKEY_RLES], "rle_%s.txt" % cfg_key) \
                                         if self.select_rle_output else ""
                    applyFakeRateWeights = self.applyFakeRateWeights  \
                      if self.isBDTtraining or not (lepton_selection == "Tight" and hadTau_selection.find("Tight") != -1) \
                      else "disabled"

                    self.jobOptions_analyze[key_analyze_job] = {
                      'ntupleFiles'              : ntupleFiles,
                      'cfgFile_modified'         : cfgFile_modified_path,
                      'histogramFile'            : histogramFile_path,
                      'logFile'                  : logFile_path,
                      'selEventsFileName_output' : rleOutputFile_path,
                      'leptonChargeSelection'    : lepton_charge_selection,
                      'electronSelection'        : electron_selection,
                      'muonSelection'            : muon_selection,
                      'lep_mva_cut'              : self.lep_mva_cut,
                      'apply_leptonGenMatching'  : self.apply_leptonGenMatching,
                      'hadTauChargeSelection'    : hadTau_charge_selection,
                      'hadTauSelection'          : hadTau_selection,
                      'apply_hadTauGenMatching'  : self.apply_hadTauGenMatching,
                      'chargeSumSelection'       : chargeSumSelection,
                      'applyFakeRateWeights'     : applyFakeRateWeights,
                      'central_or_shift'         : central_or_shift,
                      'selectBDT'                : self.isBDTtraining,
                      'syncOutput'               : syncOutput,
                      'syncTree'                 : syncTree,
                      'syncRLE'                  : syncRLE,
                      'syncRequireGenMatching'   : syncRequireGenMatching,
                      'apply_hlt_filter'         : self.hlt_filter,
                      'useNonNominal'            : self.use_nonnominal,
                      'fillGenEvtHistograms'     : True,
                    }
                    self.createCfg_analyze(self.jobOptions_analyze[key_analyze_job], sample_info, lepton_and_hadTau_selection)

                    # initialize input and output file names for hadd_stage1
                    key_hadd_stage1 = getKey(process_name, lepton_charge_selection, hadTau_charge_selection, lepton_and_hadTau_selection_and_frWeight, chargeSumSelection)
                    if not key_hadd_stage1 in self.inputFiles_hadd_stage1:
                      self.inputFiles_hadd_stage1[key_hadd_stage1] = []
                    self.inputFiles_hadd_stage1[key_hadd_stage1].append(self.jobOptions_analyze[key_analyze_job]['histogramFile'])
                    self.outputFile_hadd_stage1[key_hadd_stage1] = os.path.join(self.dirs[DKEY_HIST], "histograms_harvested_stage1_%s_%s_%s_%s_%s_%s.root" % \
                      (self.channel, process_name, lepton_charge_selection, hadTau_charge_selection, lepton_and_hadTau_selection_and_frWeight, chargeSumSelection))

                    if self.isBDTtraining:
                      self.targets.append(self.outputFile_hadd_stage1[key_hadd_stage1])

                if self.isBDTtraining or self.do_sync:
                  continue

              if self.isBDTtraining or self.do_sync:
                continue

    if self.isBDTtraining or self.do_sync:
      if self.is_sbatch:
        logging.info("Creating script for submitting '%s' jobs to batch system" % self.executable_analyze)
        self.sbatchFile_analyze = os.path.join(self.dirs[DKEY_SCRIPTS], "sbatch_analyze_%s.py" % self.channel)
        if self.isBDTtraining:
          self.createScript_sbatch_analyze(self.executable_analyze, self.sbatchFile_analyze, self.jobOptions_analyze)
        elif self.do_sync:
          self.createScript_sbatch_syncNtuple(self.executable_analyze, self.sbatchFile_analyze, self.jobOptions_analyze)
      logging.info("Creating Makefile")
      lines_makefile = []
      if self.isBDTtraining:
        self.addToMakefile_analyze(lines_makefile)
        self.addToMakefile_hadd_stage1(lines_makefile)
      elif self.do_sync:
        self.addToMakefile_syncNtuple(lines_makefile)
        outputFile_sync_path = os.path.join(self.outputDir, DKEY_SYNC, '%s.root' % self.channel)
        self.outputFile_sync['sync'] = outputFile_sync_path
        self.targets.append(outputFile_sync_path)
        self.addToMakefile_hadd_sync(lines_makefile)
      self.createMakefile(lines_makefile)
      logging.info("Done")
      return self.num_jobs

    if self.is_sbatch:
      logging.info("Creating script for submitting '%s' jobs to batch system" % self.executable_analyze)
      self.sbatchFile_analyze = os.path.join(self.dirs[DKEY_SCRIPTS], "sbatch_analyze_%s.py" % self.channel)
      self.createScript_sbatch_analyze(self.executable_analyze, self.sbatchFile_analyze, self.jobOptions_analyze)
      logging.info("Creating script for submitting '%s' jobs to batch system" % self.executable_addBackgrounds)

    logging.info("Creating Makefile")
    lines_makefile = []
    self.addToMakefile_analyze(lines_makefile)
    self.addToMakefile_hadd_stage1(lines_makefile)
    self.addToMakefile_hadd_stage2(lines_makefile)
    self.addToMakefile_add_syst_fakerate(lines_makefile)
    self.createMakefile(lines_makefile)

    logging.info("Done")

    return self.num_jobs
