#!/usr/bin/env python

from tthAnalysis.HiggsToTauTau.configs.analyzeConfig_2lss_1tau import analyzeConfig_2lss_1tau
from tthAnalysis.HiggsToTauTau.jobTools import query_yes_no
from tthAnalysis.HiggsToTauTau.analysisSettings import systematics, get_lumi
from tthAnalysis.HiggsToTauTau.runConfig import tthAnalyzeParser, filter_samples
from tthAnalysis.HiggsToTauTau.common import logging, load_samples

import os
import sys
import getpass
import re

# E.g. to run: ./test/tthAnalyzeRun_2lss_1tau.py -v 2017Dec13 -m default -e 2017

mode_choices         = [
  'default', 'addMEM', 'forBDTtraining_beforeAddMEM', 'forBDTtraining_afterAddMEM', 'sync', 'sync_wMEM',
]
sys_choices      = [ 'full', 'internal' ] + systematics.an_tth_opts + [ 'MEM_2lss_1tau' ]
systematics.full = systematics.an_tth
systematics.internal = systematics.an_internal_no_mem

parser = tthAnalyzeParser()
parser.add_modes(mode_choices)
parser.add_sys(sys_choices)
parser.add_preselect()
parser.add_rle_select()
parser.add_nonnominal()
parser.add_tau_id_wp()
parser.add_tau_id()
parser.add_hlt_filter()
parser.add_files_per_job()
parser.add_use_home()
parser.add_jet_cleaning()
parser.add_gen_matching()
parser.add_sideband()
parser.enable_regrouped_jerc()
parser.add_split_trigger_sys()
args = parser.parse_args()

# Common arguments
era                = args.era
version            = args.version
dry_run            = args.dry_run
no_exec            = args.no_exec
auto_exec          = args.auto_exec
check_output_files = not args.not_check_input_files
debug              = args.debug
sample_filter      = args.filter
num_parallel_jobs  = args.num_parallel_jobs
running_method     = args.running_method
keep_logs          = args.keep_logs

# Additional arguments
mode              = args.mode
systematics_label = args.systematics
use_preselected   = args.use_preselected
rle_select        = os.path.expanduser(args.rle_select)
use_nonnominal    = args.original_central
hlt_filter        = args.hlt_filter
files_per_job     = args.files_per_job
use_home          = args.use_home
jet_cleaning      = args.jet_cleaning
gen_matching      = args.gen_matching
sideband          = args.sideband
tau_id            = args.tau_id
regroup_jerc      = args.enable_regrouped_jerc
split_trigger_sys = args.split_trigger_sys

if regroup_jerc:
  if 'full' not in systematics_label:
    raise RuntimeError("Regrouped JEC or split JER was enabled but not running with full systematics")
  if regroup_jerc == 'both':
    systematics.full.extend(systematics.JEC_regrouped + systematics.JER_split)
  elif regroup_jerc == 'jes':
    systematics.full.extend(systematics.JEC_regrouped)
  elif regroup_jerc == 'jer':
    systematics.full.extend(systematics.JER_split)
  else:
    raise RuntimeError("Invalid choice: %s" % regroup_jerc)
if split_trigger_sys == 'yes':
  for trigger_sys in systematics.triggerSF:
    del systematics.internal[systematics.internal.index(trigger_sys)]
    del systematics.full[systematics.full.index(trigger_sys)]
if split_trigger_sys in [ 'yes', 'both' ]:
  systematics.internal.extend(systematics.triggerSF_2lss)
  systematics.full.extend(systematics.triggerSF_2lss)

# Use the arguments
if "MEM" in mode:
  systematics.full += systematics.MEM_2lss_1tau
  systematics.internal += systematics.MEM_2lss_1tau
central_or_shifts = []
for systematic_label in systematics_label:
  for central_or_shift in getattr(systematics, systematic_label):
    if central_or_shift not in central_or_shifts:
      if central_or_shift in systematics.MEM_2lss_1tau and "MEM" not in mode:
        continue
      central_or_shifts.append(central_or_shift)
do_sync = mode.startswith('sync')
lumi = get_lumi(era)
jet_cleaning_by_index = (jet_cleaning == 'by_index')
gen_matching_by_index = (gen_matching == 'by_index')

MEMbranch                = ''
lepton_charge_selections = [ "SS" ] if mode.find("forBDTtraining") != -1 else [ "OS", "SS" ]

hadTauWP_map = {
  'dR03mva' : 'Loose',
  'deepVSj' : 'VLoose',
}
hadTau_selection = tau_id + hadTauWP_map[tau_id]

hadTauWP_veto_map = {
  'dR03mva' : 'Medium',
  'deepVSj' : 'Medium',
}
hadTau_selection_veto = tau_id + hadTauWP_veto_map[tau_id]

if sideband == 'disabled':
  chargeSumSelections = [ "OS" ]
elif sideband == 'enabled':
  chargeSumSelections = [ "OS", "SS" ]
elif sideband == 'only':
  chargeSumSelections = [ "SS" ]
else:
  raise ValueError("Invalid choice for the sideband: %s" % sideband)

MEMsample_base = "addMEM_2lss1tau_{}".format(hadTau_selection)

if mode == "default":
  samples = load_samples(era, suffix = "preselected" if use_preselected else "")
elif mode == "addMEM":
  if not use_preselected:
    raise ValueError("MEM branches can be read only from preselected Ntuples")
  samples = load_samples(era, suffix = MEMsample_base)
  MEMbranch = 'memObjects_2lss_1tau_lepFakeable_tauTight_{}'.format(hadTau_selection)
elif mode == "forBDTtraining_beforeAddMEM":
  if use_preselected:
    raise ValueError("Makes no sense to use preselected samples w/ BDT training mode")
  samples = load_samples(era, suffix = "BDT")
  if args.tau_id_wp:
    tau_id = args.tau_id[:7]
  hadTau_selection_relaxed = tau_id + hadTauWP_map[tau_id]
elif mode == "forBDTtraining_afterAddMEM":
  if use_preselected:
    raise ValueError("Makes no sense to use preselected samples w/ BDT training mode")
  samples = load_samples(era, suffix = "{}_BDT".format(MEMsample_base))
  if args.tau_id_wp:
    tau_id = args.tau_id[:7]
  hadTau_selection_relaxed = tau_id + hadTauWP_map[tau_id]
  MEMbranch                = 'memObjects_2lss_1tau_lepLoose_tauTight_{}'.format(hadTau_selection_relaxed)
elif mode == "sync_wMEM":
  if not use_preselected:
    raise ValueError("MEM branches can be read only from preselected Ntuples")
  samples = load_samples(era, suffix = "{}_sync{}".format(MEMsample_base, '' if use_nonnominal else "_nom"))
  MEMbranch = 'memObjects_2lss_1tau_lepFakeable_tauTight_{}'.format(hadTau_selection)
elif mode == "sync":
  sample_suffix = "sync" if use_nonnominal else "sync_nom"
  if use_preselected:
    sample_suffix = "preselected_{}".format(sample_suffix)
  samples = load_samples(era, suffix = sample_suffix)
else:
  raise ValueError("Invalid mode: %s" % mode)

if not mode.startswith("sync"):
  for sample_name, sample_info in samples.items():
    if sample_name == 'sum_events':
      continue
    if re.match("(^WZTo3LNu$|^WZTo3LNu_ext(\d)?$)", sample_info["process_name_specific"]):
      sample_info["use_it"] = True
      sample_info["skipEvery"] = 2
    if era == "2018" and sample_info["sample_category"] == "tHq" and sample_info["use_it"]:
      sample_info["skipEvery"] = 3

if __name__ == '__main__':
  logging.info(
    "Running the jobs with the following systematic uncertainties enabled: %s" % \
    ', '.join(central_or_shifts)
  )
  if not use_preselected:
    logging.warning('Running the analysis on fully inclusive samples!')

  if sample_filter:
    samples = filter_samples(samples, sample_filter)

  if args.tau_id_wp:
    logging.info("Changing tau ID working point: %s -> %s" % (hadTau_selection, args.tau_id_wp))
    hadTau_selection = args.tau_id_wp
    if MEMbranch:
      MEMbranch = 'memObjects_2lss_1tau_lepFakeable_tauTight_{}'.format(hadTau_selection)

  configDir = os.path.join("/scratch-persistent", getpass.getuser(), "ttHAnalysis", era, version),
  localDir  = os.path.join("/home",               getpass.getuser(), "ttHAnalysis", era, version)
  outputDir = os.path.join("/local",              getpass.getuser(), "ttHAnalysis", era, version)

  analysis = analyzeConfig_2lss_1tau(
    configDir                 = configDir,
    localDir                  = localDir,
    outputDir                 = outputDir,
    executable_analyze        = "analyze_2lss_1tau",
    cfgFile_analyze           = "analyze_2lss_1tau_cfg.py",
    samples                   = samples,
    MEMbranch                 = MEMbranch,
    lepton_charge_selections  = lepton_charge_selections,
    hadTau_selection          = hadTau_selection,
    hadTau_selection_veto     = hadTau_selection_veto,
    # CV: apply "fake" background estimation to leptons only and not to hadronic taus, as discussed on slide 10 of
    #     https://indico.cern.ch/event/597028/contributions/2413742/attachments/1391684/2120220/16.12.22_ttH_Htautau_-_Review_of_systematics.pdf
    applyFakeRateWeights      = "2lepton",
    chargeSumSelections       = chargeSumSelections,
    jet_cleaning_by_index     = jet_cleaning_by_index,
    gen_matching_by_index     = gen_matching_by_index,
    central_or_shifts         = central_or_shifts,
    max_files_per_job         = files_per_job,
    era                       = era,
    use_lumi                  = True,
    lumi                      = lumi,
    check_output_files        = check_output_files,
    running_method            = running_method,
    num_parallel_jobs         = num_parallel_jobs,
    executable_addBackgrounds = "addBackgrounds",
    # CV: use common executable for estimating jet->lepton and jet->tau_h fake background
    executable_addFakes       = "addBackgroundLeptonFakes",
    executable_addFlips       = "addBackgroundLeptonFlips",
    histograms_to_fit         = {
      "EventCounter"    : {},
      "mTauTauVis1"     : {},
      "mTauTauVis2"     : {},
      "mTauTauVis"      : {},
      "memOutput_LR"    : {},
      "mvaOutput_final" : {},
      "output_NN_ttH" : {},
      "output_NN_tH" : {},
      "output_NN_rest" : {}
    },
    select_rle_output         = True,
    dry_run                   = dry_run,
    do_sync                   = do_sync,
    isDebug                   = debug,
    rle_select                = rle_select,
    use_nonnominal            = use_nonnominal,
    hlt_filter                = hlt_filter,
    use_home                  = use_home,
    keep_logs                 = keep_logs,
    submission_cmd            = sys.argv,
  )

  if mode.find("forBDTtraining") != -1:
    analysis.set_BDT_training(hadTau_selection_relaxed)

  job_statistics = analysis.create()
  for job_type, num_jobs in job_statistics.items():
    logging.info(" #jobs of type '%s' = %i" % (job_type, num_jobs))

  if auto_exec:
    run_analysis = True
  elif no_exec:
    run_analysis = False
  else:
    run_analysis = query_yes_no("Start jobs ?")
  if run_analysis:
    analysis.run()
  else:
    sys.exit(0)
