#!/usr/bin/env python

from tthAnalysis.HiggsToTauTau.jobTools import run_cmd, human_size, getmtime
from tthAnalysis.HiggsToTauTau.analysisSettings import Triggers, HTXS_BINS
from tthAnalysis.HiggsToTauTau.safe_root import ROOT
from tthAnalysis.HiggsToTauTau.common import logging, SmartFormatter

from tthAnalysis.NanoAODTools.tHweights_cfi import tHweights

import argparse
import os.path
import sys
import imp
import jinja2
import re
import copy
import itertools
import time
import shutil
import datetime
import math
import collections
import array

HISTOGRAM_COUNT                                         = 'Count'
HISTOGRAM_COUNTWEIGHTED                                 = 'CountWeighted'
HISTOGRAM_COUNTWEIGHTED_FULL                            = "{}Full".format(HISTOGRAM_COUNTWEIGHTED)
HISTOGRAM_COUNTWEIGHTED_L1PREFIRE_NOM                   = 'CountWeightedL1PrefireNom'
HISTOGRAM_COUNTWEIGHTED_L1PREFIRE                       = 'CountWeightedL1Prefire'
HISTOGRAM_COUNTWEIGHTED_LHESCALE                        = 'CountWeightedLHEWeightScale'
HISTOGRAM_COUNTWEIGHTED_LHESCALE_FULL                   = HISTOGRAM_COUNTWEIGHTED_LHESCALE.replace(HISTOGRAM_COUNTWEIGHTED, HISTOGRAM_COUNTWEIGHTED_FULL)
HISTOGRAM_COUNTWEIGHTED_LHESCALE_L1PREFIRE_NOM          = 'CountWeightedLHEWeightScaleL1PrefireNom'
HISTOGRAM_COUNTWEIGHTED_LHESCALE_L1PREFIRE_NOM_FULL     = HISTOGRAM_COUNTWEIGHTED_LHESCALE_L1PREFIRE_NOM.replace(HISTOGRAM_COUNTWEIGHTED, HISTOGRAM_COUNTWEIGHTED_FULL)
HISTOGRAM_COUNTWEIGHTED_LHEENVELOPE                     = 'CountWeightedLHEEnvelope'
HISTOGRAM_COUNTWEIGHTED_LHEENVELOPE_L1PREFIRE_NOM       = 'CountWeightedLHEEnvelopeL1PrefireNom'
HISTOGRAM_COUNTWEIGHTED_PSWEIGHT                        = 'CountWeightedPSWeight'
HISTOGRAM_COUNTWEIGHTED_PSWEIGHT_FULL                   = HISTOGRAM_COUNTWEIGHTED_PSWEIGHT.replace(HISTOGRAM_COUNTWEIGHTED, HISTOGRAM_COUNTWEIGHTED_FULL)
HISTOGRAM_COUNTWEIGHTED_PSWEIGHT_L1PREFIRE_NOM          = 'CountWeightedPSWeightL1PrefireNom'
HISTOGRAM_COUNTWEIGHTED_PSWEIGHT_ORIGINAL               = 'CountWeightedPSWeightOriginalXWGTUP'
HISTOGRAM_COUNTWEIGHTED_PSWEIGHT_ORIGINAL_L1PREFIRE_NOM = 'CountWeightedPSWeightOriginalXWGTUPL1PrefireNom'

EVENTS_TREE = 'Events'
BRANCH_LHEPDFWEIGHT = 'LHEPdfWeight'
BRANCH_NLHEREWEIGHTINGWEIGHT = 'nLHEReweightingWeight'
BRANCH_NPSWEIGHTS = 'nPSWeight'

HISTOGRAM_COUNT_KEY = 'histogram_count'
TREE_COUNT_KEY      = 'tree_count'
FSIZE_KEY           = 'fsize'
BRANCH_NAMES_KEY    = 'branch_names'

LHE_REGEX     = re.compile('(n|)LHE(Scale|Pdf)Weight')
LHE_DOC_REGEX = re.compile('LHE pdf variation weights \(w_var \/ w\_nominal\) for LHA IDs (?P<lha_start>[0-9]+) - (?P<lha_end>[0-9]+)')

HISTOGRAM_COUNT_COMMON_MC_TMP = [
  HISTOGRAM_COUNTWEIGHTED,
  HISTOGRAM_COUNTWEIGHTED_LHESCALE,
  HISTOGRAM_COUNTWEIGHTED_LHEENVELOPE,
  HISTOGRAM_COUNTWEIGHTED_PSWEIGHT,
  HISTOGRAM_COUNTWEIGHTED_PSWEIGHT_ORIGINAL,
]
HISTOGRAM_COUNT_COMMON_MC = copy.deepcopy(HISTOGRAM_COUNT_COMMON_MC_TMP)
HISTOGRAM_COUNT_COMMON_MC.extend(
  histogram_name.replace(HISTOGRAM_COUNTWEIGHTED, HISTOGRAM_COUNTWEIGHTED_FULL) for histogram_name in HISTOGRAM_COUNT_COMMON_MC_TMP
)
HISTOGRAM_COUNT_EXTENDED_MC_TMP = [
  HISTOGRAM_COUNTWEIGHTED_L1PREFIRE_NOM,
  HISTOGRAM_COUNTWEIGHTED_L1PREFIRE,
  HISTOGRAM_COUNTWEIGHTED_LHESCALE_L1PREFIRE_NOM,
  HISTOGRAM_COUNTWEIGHTED_LHEENVELOPE_L1PREFIRE_NOM,
  HISTOGRAM_COUNTWEIGHTED_PSWEIGHT_L1PREFIRE_NOM,
  HISTOGRAM_COUNTWEIGHTED_PSWEIGHT_ORIGINAL_L1PREFIRE_NOM,
]
HISTOGRAM_COUNT_EXTENDED_MC = copy.deepcopy(HISTOGRAM_COUNT_EXTENDED_MC_TMP)
HISTOGRAM_COUNT_EXTENDED_MC.extend(
  histogram_name.replace(HISTOGRAM_COUNTWEIGHTED, HISTOGRAM_COUNTWEIGHTED_FULL) for histogram_name in HISTOGRAM_COUNT_EXTENDED_MC_TMP
)
LHESCALEARR = [
  HISTOGRAM_COUNTWEIGHTED_LHESCALE,
  HISTOGRAM_COUNTWEIGHTED_LHESCALE_FULL,
  HISTOGRAM_COUNTWEIGHTED_LHESCALE_L1PREFIRE_NOM,
  HISTOGRAM_COUNTWEIGHTED_LHESCALE_L1PREFIRE_NOM_FULL,
]
TH_INDICES = [ coupling.idx.value() for coupling in tHweights if coupling.idx.value() >= 0 ]

# see https://github.com/cms-nanoAOD/cmssw/blob/9a2728ac9f44fc45ba1aa56389e28c594207c0fe/PhysicsTools/NanoAOD/python/nano_cff.py#L99-L104
LHE_DOC = {
  91400  : { 'name' : 'PDF4LHC15_nnlo_30_pdfas',    'count' : 33  },
  306000 : { 'name' : 'NNPDF31_nnlo_hessian_pdfas', 'count' : 103 },
  260000 : { 'name' : 'NNPDF30_nlo_as_0118',        'count' : 101 },
  262000 : { 'name' : 'NNPDF30_lo_as_0130',         'count' : 101 },
  292000 : { 'name' : 'NNPDF30_nlo_nf_4_pdfas',     'count' : 103 },
  292200 : { 'name' : 'NNPDF30_nlo_nf_5_pdfas',     'count' : 103 },
}

try:
    from urllib.parse import urlparse
except ImportError:
     from urlparse import urlparse

class path_info:
  def __init__(self, name):
    self.name          = name
    self.kind          = 'F' if os.path.isfile(self.name) else 'D'
    self.size          = os.path.getsize(name)
    self.basename      = os.path.basename(self.name)
    self.depth         = len(self.name.split(os.path.sep)) - 1
    self.sparent_depth = -1

  def __str__(self):
    return self.name

  def __repr__(self):
    return self.name

  def isfile(self):
    return self.kind == 'F'

  def isdir(self):
    return self.kind == 'D'


class FileTracker:
  def __init__(self):
    self.zero_file_size  = []
    self.zombie_files    = []
    self.corrupted_files = []

def load_dict(path, name):
  if not os.path.isfile(path):
    logging.error("No such dictionary file: {dict_path}".format(dict_path = path))
    sys.exit(1)
  imp_dict = imp.load_source('', path)
  if not hasattr(imp_dict, name):
    logging.error("No such dictionary in the file '{dict_path}': {dict_name}".format(
      dict_path = path, dict_name = name,
    ))
    sys.exit(1)
  samples = getattr(imp_dict, name)
  return samples

header_str = """from collections import OrderedDict as OD

# file generated at {{ date }} with the following command:
# {{ command }}

{{ dict_name }} = OD()

"""

dictionary_entry_str = """{{ dict_name }}["{{ dbs_name }}"] = OD([
  ("type",                            "{{ sample_type }}"),
  ("sample_category",                 "{{ sample_category }}"),
  ("process_name_specific",           "{{ process_name_specific }}"),
  ("nof_files",                       {{ nof_files }}),
  ("nof_db_files",                    {{ nof_db_files }}),
  ("nof_events",                      { {%- for histogram_name, event_counts in nof_events.items() %}
    {{ "%-80s"|format("'%s'"|format(histogram_name)) }} : [ {% for event_count in event_counts -%}{{ event_count }}, {% endfor %}],
  {%- endfor %}
  }),
  ("nof_tree_events",                 {{ nof_tree_events }}),
  ("nof_db_events",                   {{ nof_db_events }}),
  ("fsize_local",                     {{ fsize_local }}), # {{ fsize_local_human }}, avg file size {{ avg_fsize_local_human }}
  ("fsize_db",                        {{ fsize_db }}), # {{ fsize_db_human }}, avg file size {{ avg_fsize_db_human }}
  ("use_it",                          {{ use_it }}),{% if sample_type == "mc" %}
  ("xsection",                        {{ xsection }}),
  ("genWeight",                       {{ genWeight }}),{% endif %}
  ("triggers",                        {{ triggers }}),
  ("has_LHE",                         {{ has_LHE }}),
  ("nof_PSweights",                   {{ nof_PSweights }}),
  ("LHE_set",                         "{{ LHE_set }}"),
  ("nof_reweighting",                 {{ nof_reweighting }}),
  ("local_paths",
    [
{{ paths }}
    ]
  ),
  ("missing_completely",           [
{{ missing_completely }}
  ]),
  ("missing_from_superset",        [
{{ missing_from_superset }}
  ]),
  ("missing_hlt_paths",            [
{{ missing_hlt_paths }}
  ]),
  ("hlt_paths",                    [
{{ hlt_paths }}
  ]),
])
"""

dictionary_sum_events_str = """{{ dict_name }}["sum_events"] = [{%- for sample_list in sample_lists %}
  [ {% for sample in sample_list %}{{ "%-50s"|format("'%s',"|format(sample)) }} {% endfor %} ],
{%- endfor %}
]
"""

path_entry_str = """      OD([
        ("path",      "{{ path }}"),
        ("selection", "{{ selection }}"),
        ("blacklist", {{ blacklist }}),
      ]),
"""

missing_branches_str = """{%- if is_available -%}
  {%- for missing_branch in missing_branches %}
    "{{ missing_branch }}",
  {%- endfor -%}
{%- else %}
    # not computed
{%- endif -%}
"""

sh_str = """#!/bin/bash

{{ cmd }}
"""

class PathEntry:
  def __init__(self, path, indices, histogram_names, blacklist_ext):
    self.path            = path
    self.indices         = indices
    self.is_presel       = 'wPresel' in self.path

    nof_events_transposed = {
      histogram_name : [ [] for _ in range(nBins) ] \
      for histogram_name, nBins in histogram_names.items() if nBins > 0
    }
    for index_entry in self.indices.values():
      for histogram_name in index_entry[HISTOGRAM_COUNT_KEY]:
        for idxBin, nof_events in enumerate(index_entry[HISTOGRAM_COUNT_KEY][histogram_name]):
          nof_events_transposed[histogram_name][idxBin].append(nof_events)
    self.nof_events = {
      histogram_name : [ math.fsum(entry[idxBin]) for idxBin in range(len(entry)) ] \
      for histogram_name, entry in nof_events_transposed.items()
    }

    self.nof_tree_events = sum(index_entry[TREE_COUNT_KEY]      for index_entry in self.indices.values())
    self.fsize           = sum(index_entry[FSIZE_KEY]           for index_entry in self.indices.values())
    self.nof_files       = max(self.indices.keys())
    self.blacklist_ext   = blacklist_ext
    self.blacklist       = []
    self.selection       = [] # if empty, select all

  def __repr__(self):
    return self.path

def get_triggers(process_name_specific, is_data, era):
  if 'SingleElec' in process_name_specific:
    return ['1e', '1e1tau']
  if 'SingleMuon' in process_name_specific:
    return ['1mu', '1mu1tau']
  if 'DoubleEG' in process_name_specific:
    return ['2e', '3e']
  if 'DoubleMuon' in process_name_specific:
    return ['2mu', '3mu']
  if 'MuonEG' in process_name_specific:
    return ['1e1mu', '2e1mu', '1e2mu']
  if 'Tau' in process_name_specific:
    return ['1e1tau', '1mu1tau', '2tau']
  if 'EGamma' in process_name_specific: # merge of SingleElectron and DoubleEG PDs
    return ['1e', '1e1tau', '2e', '3e']
  if is_data:
    raise ValueError("Expected MC!")
  return [
    '1e', '1mu', '2e', '2mu', '1e1mu', '3e', '3mu', '2e1mu', '1e2mu', '1e1tau', '1mu1tau', '2tau'
  ]

def get_array_type(tree, branch_name, array_multiplier = 1):
  branch = tree.GetBranch(branch_name)
  leaf = branch.GetLeaf(branch_name)
  leaf_type = leaf.GetTypeName()

  if leaf_type == 'UInt_t':
    arr_type = 'I'
  elif leaf_type == 'Int_t':
    arr_type = 'i'
  elif leaf_type == 'ULong64_t':
    arr_type = 'L'
  elif leaf_type == 'Long64_t':
    arr_type = 'l'
  elif leaf_type == 'Float_t':
    arr_type = 'f'
  elif leaf_type == 'Double_t':
    arr_type = 'd'
  elif leaf_type == 'UChar_t':
    arr_type = 'B'
  elif leaf_type == 'Char_t':
    arr_type = 'b'
  elif leaf_type == 'UShort_t':
    arr_type = 'H'
  elif leaf_type == 'Short_t':
    arr_type = 'h'
  elif leaf_type == 'Bool_t':
    arr_type = 'b' # use char; bool allocates 8 bits anyways
  else:
    raise ValueError("Invalid leaf type: %s" % leaf_type)
  if arr_type in ['d', 'f']:
    return (arr_type, [0.] * array_multiplier)
  else:
    return (arr_type, [0] * array_multiplier)

def get_nof_events_sum_str(nof_events_sum, histogram_name):
  if histogram_name == HISTOGRAM_COUNT or histogram_name.startswith('{}_'.format(HISTOGRAM_COUNT)):
    nof_events_sum_str = str(int(nof_events_sum))
  else:
    nof_events_sum_str = '{:.8e}'.format(nof_events_sum)
  if nof_events_sum == 0:
    nof_events_sum_str = "0."
  assert(nof_events_sum_str)
  return nof_events_sum_str

def process_paths(meta_dict, key, count_histograms):
  local_paths = meta_dict[key]['paths']
  nof_events = max(path_entry.nof_tree_events for path_entry in local_paths)
  is_presel_list = [ path_entry.is_presel for path_entry in local_paths ]
  assert(all(is_presel_list) or all(not flag for flag in is_presel_list))
  is_presel = all(is_presel_list)

  # build the blacklists for all the paths
  for local_path in local_paths:
    local_path.blacklist = list(sorted(
      list(set(range(1, local_path.nof_files + 1)) - set(local_path.indices.keys()))
    ))

  if len(local_paths) > 1:
    # sort the paths by the largest coverage only if the Ntuples are unskimmed
    local_paths_sorted = list(sorted(
      local_paths,
      key     = lambda local_path: local_path.nof_tree_events,
      reverse = True,
    ))

    local_path_cand_idxs = []
    if not is_presel:
      for local_path_cand_idx, local_path_sorted in enumerate(local_paths_sorted):
        if nof_events == local_paths_sorted[local_path_cand_idx].nof_tree_events:
          # the path with the largest coverage already spans all possible files
          local_path_cand_idxs.append(local_path_cand_idx)
    else:
      local_path_cand_idxs = list(range(len(local_paths_sorted)))
    assert(local_path_cand_idxs)
    # if there are multiple paths with the same coverage, pick the later one
    local_paths_sorted_by_date = sorted(
      local_path_cand_idxs,
      key = lambda local_path_cand_idx: getmtime(local_paths_sorted[local_path_cand_idx].path),
      reverse = True
    )
    local_paths_sorted = [ local_paths_sorted[local_paths_sorted_by_date[0]] ]
  elif len(local_paths) == 1:
    local_paths_sorted = local_paths
  else:
    assert(False)

  assert(len(local_paths_sorted) == 1)
  local_path_choice = local_paths_sorted[0]
  # let's compute the number of files, events and the list of blacklisted files
  nof_events = collections.OrderedDict()
  histogram_names = meta_dict[key]['histogram_names']
  for histogram_name, nBins in histogram_names.items():
    if nBins < 0:
      continue
    nof_events[histogram_name] = []
    for idxBin in range(nBins):
      nof_events_sum = math.fsum(
        index_entry[HISTOGRAM_COUNT_KEY][histogram_name][idxBin] for index_entry in local_path_choice.indices.values() \
        if histogram_name in index_entry[HISTOGRAM_COUNT_KEY]
      )
      nof_events[histogram_name].append(get_nof_events_sum_str(nof_events_sum, histogram_name))

  nof_tree_events = sum(index_entry[TREE_COUNT_KEY] for index_entry in local_path_choice.indices.values())
  fsize           = sum(index_entry[FSIZE_KEY]      for index_entry in local_path_choice.indices.values())

  process_name = meta_dict[key]['process_name_specific']
  if count_histograms and process_name in count_histograms:
    count_histograms_process = count_histograms[process_name]
    # the assumption is that all necessary event counts are already stored in the auxiliary file
    for count_histogram_name in count_histograms_process:
      if count_histogram_name in nof_events and HISTOGRAM_COUNTWEIGHTED_FULL not in count_histogram_name:
        event_counts_int = [ float(event_count) for event_count in nof_events[count_histogram_name] ]
        event_counts_ext = count_histograms_process[count_histogram_name]
        event_counts_int_len = len(event_counts_int)
        event_counts_ext_len = len(event_counts_ext)
        if event_counts_int_len != event_counts_ext_len:
          raise RuntimeError(
            "Expected %d event counts but got %d instead in %s of process %s" % \
            (event_counts_int_len, event_counts_ext_len, count_histogram_name, process_name)
          )
        for count_idx in range(event_counts_int_len):
          event_count_int = event_counts_int[count_idx]
          event_count_ext = event_counts_ext[count_idx]
          if event_count_int != 0.:
            event_count_diff = abs(event_count_int - event_count_ext) / event_count_int
            if event_count_diff > 1.e-2:
              error_msg = "Observed too large difference of %.3f%% in %s at index %d of process %s" % \
                (event_count_diff * 100., count_histogram_name, count_idx, process_name)
              if local_path_choice.blacklist_ext:
                logging.error(error_msg)
              else:
                raise RuntimeError(error_msg)
          else:
            if event_count_int != event_count_ext:
              raise RuntimeError(
                "Expected 0 events but observed %.3f in %s at index %d of process %s" % \
                (event_count_ext, count_histogram_name, count_idx, process_name)
              )
        #continue
      if 'Pdf' in count_histogram_name:
        continue
      if count_histogram_name in LHESCALEARR and len(count_histograms_process[count_histogram_name]) != 9:
        # Ignore LHE scale weights if we don't have the correct number of them
        continue
      nof_events[count_histogram_name] = [
        get_nof_events_sum_str(nof_events_sum, count_histogram_name) \
        for nof_events_sum in count_histograms_process[count_histogram_name]
      ]

  nof_events_zeroes = []
  for count_histogram_name in nof_events:
    if all(nof_events_sum == '0.' for nof_events_sum in nof_events[count_histogram_name]):
      nof_events_zeroes.append(count_histogram_name)
  for count_histogram_name in nof_events_zeroes:
    del nof_events[count_histogram_name]

  meta_dict[key]['nof_events']      = nof_events
  meta_dict[key]['nof_tree_events'] = nof_tree_events
  meta_dict[key]['fsize_local']     = fsize
  meta_dict[key]['local_paths'] = [{
    'path'      : local_path_choice.path,
    'selection' : '*',
    'blacklist' : list(sorted(set(local_path_choice.blacklist_ext) | set(local_path_choice.blacklist))),
  }]
  meta_dict[key]['nof_files'] = local_path_choice.nof_files

def get_lhe_set(tree):
  lhe_branch = tree.GetBranch(BRANCH_LHEPDFWEIGHT)
  if lhe_branch:
    lhe_doc = lhe_branch.GetTitle()
    lhe_match = LHE_DOC_REGEX.match(lhe_doc)
    if lhe_match:
      lhe_start = int(lhe_match.group('lha_start'))
      lhe_end = int(lhe_match.group('lha_end'))
      lhe_doc = 'LHA IDs {} - {}'.format(lhe_start, lhe_end)
      lhe_count = lhe_end - lhe_start + 1
      lhe_val = {}
      if lhe_start in LHE_DOC:
        lhe_val = LHE_DOC[lhe_start]
      elif (lhe_start - 1) in LHE_DOC:
        lhe_val = LHE_DOC[lhe_start - 1]
      if lhe_val:
        lhe_doc += ' -> {name} PDF set, expecting {count} weights'.format(**lhe_val)
      else:
        lhe_doc += ' -> unrecognizable PDF set'
      lhe_doc += ' (counted {} weights)'.format(lhe_count)
    return lhe_doc
  return ''

def has_LHE(indices):
  branch_names = set()
  for index_entry in indices.values():
    branch_names.update(index_entry[BRANCH_NAMES_KEY])
  return any(map(lambda branch_name: bool(LHE_REGEX.match(branch_name)), branch_names))

def get_missing_from_superset(indices):
  branch_names_superset = set()
  for index_entry in indices.values():
    branch_names_superset.update(index_entry[BRANCH_NAMES_KEY])
  missing_branches_superset = set()
  for file_idx, index_entry in indices.items():
    missing_branches = branch_names_superset - set(index_entry[BRANCH_NAMES_KEY])
    missing_branches_superset.update(missing_branches)
  return list(sorted(list(missing_branches_superset)))

def get_hlt_paths(indices):
  branch_names_union = set.union(*[
    set(index_entry[BRANCH_NAMES_KEY]) for index_entry in indices.values()
  ])
  hlt_union = [ branch_name for branch_name in branch_names_union if branch_name.startswith('HLT_') ]
  return hlt_union

def get_missing_hlt_paths(required_triggers, indices, all_paths):
  branch_names_intersection = set.intersection(*[
    set(index_entry[BRANCH_NAMES_KEY]) for index_entry in indices.values()
  ])
  required_paths = set.union(*[ all_paths[trigger_name] for trigger_name in required_triggers ])
  missing_paths = list(sorted(list(required_paths - branch_names_intersection)))
  return missing_paths

def get_dir_entries(path):
  return list(map(lambda dir_entry: path_info(os.path.join(path.name, dir_entry)), os.listdir(path.name)))

def get_is_njet(process_name):
  return process_name.startswith(
    tuple('DYToLL_{}J'.format(i) for i in range(3)) + \
    ('DYJetsToLL_M-50_amcatnloFXFX', 'WJetsToLNu_HT', 'DYJetsToLL_M50_HT', 'DYJetsToLL_M-10to50')
  )

def traverse_single(meta_dict, path_obj, key, check_every_event, missing_branches,
                    filetracker, file_idx, era, triggerTable, count_histograms, lost_ntuples):
  ''' Assume that the following subdirectories are of the form: 0000, 0001, 0002, ...
      In these directories we expect root files of the form: tree_1.root, tree_2.root, ...
      If either of those assumptions doesn't hold, we bail out; no clever event count needed
  :param meta_dict:         Meta-dictionary
  :param path_obj:          Contains meta-information about a path
  :param key:               Key to the meta-dictionary the entry of which will be updated
  :param check_every_event: Loop over all events for error checking purposes
  :param missing_branches:  Find missing branches from the superset of branches in a sample
  :param filetracker:       An instance of FileTracker() for logging broken files
  :param file_idx:          Index of the corrupted file
  :param triggerTable:      Trigger instance containing a list of required triggers for the era
  :param count_histograms   Externally provided event sums
  :return: None
  '''
  if 'paths' not in meta_dict[key]:
    meta_dict[key]['paths'] = []
  if path_obj.name in meta_dict[key]['paths']:
    logging.warning("Path {path} has already been traversed".format(path = path_obj.name))
    return

  logging.info("Single-traversing {path}".format(path = path_obj.name))
  entries = get_dir_entries(path_obj)
  entries_valid = []
  for entry in entries:
    if not entry.isdir():
      continue
    if len(entry.basename) != 4:
      continue
    try:
      int(entry.basename)
    except:
      continue
    entries_valid.append(entry)

  digit_regex = re.compile(r"tree_(?P<i>\d+)\.root$")
  is_data = meta_dict[key]['sample_category'] == 'data_obs'
  is_rwgt = meta_dict[key]['sample_category'] in [ "tHq", "tHW", "ttH_ctcvcp" ]
  is_htxs = meta_dict[key]['sample_category'].startswith('ttH')
  process_name = meta_dict[key]['process_name_specific']
  is_lo = 'amcatnlo' not in key
  is_njet = get_is_njet(process_name)
  is_ht = process_name.startswith(
    tuple('W{}JetsToLNu'.format(i) for i in range(1, 5)) + tuple('DY{}JetsToLL_M-50'.format(i) for i in range(1, 5))
  )
  is_njet_ht = process_name.startswith('WJetsToLNu_madgraphMLM') or (process_name.startswith('DYJetsToLL_M-50') and is_lo)
  assert(not (is_htxs and is_njet))

  lheScaleArr = copy.deepcopy(LHESCALEARR)
  th_arr = [ -1 ]
  if is_rwgt:
    th_arr.extend(TH_INDICES)
  aux_arr = [ "" ]

  ht_arr = [ 0, 70, 100, 200, 400, 600, 800, 1200, 2500 ]
  aux_njet_arr = [ "LHENjet{}".format(njet) for njet in range(5) ]
  aux_ht_arr = [ "LHEHT{}to{}".format(ht_arr[ht_idx], ht_arr[ht_idx + 1]) for ht_idx in range(len(ht_arr) - 1) ] + \
               [ "LHEHT{}toInf".format(ht_arr[-1]) ]

  if is_htxs:
    aux_arr.extend(HTXS_BINS)
  elif is_njet:
    aux_arr.extend(aux_njet_arr)
  elif is_ht:
    aux_arr.extend(aux_ht_arr)
  elif is_njet_ht:
    for aux_njet_bin in aux_njet_arr:
      for aux_ht_bin in aux_ht_arr:
        aux_arr.append("{}_{}".format(aux_njet_bin, aux_ht_bin))

  histogram_names = collections.OrderedDict([ ( HISTOGRAM_COUNT, -1 ) ])
  if not is_data:
    for tH_idx in th_arr:
      for aux_bin in aux_arr:
        for histogram_name_tmp in HISTOGRAM_COUNT_COMMON_MC:
          histogram_name = histogram_name_tmp
          if tH_idx >= 0:
            histogram_name += "_rwgt{}".format(tH_idx)
          if aux_bin:
            histogram_name += "_{}".format(aux_bin)
            histogram_name_count = "{}_{}".format(HISTOGRAM_COUNT, aux_bin)
            if histogram_name_count not in histogram_names:
              histogram_names[histogram_name_count] = -1
          histogram_names[histogram_name] = -1
        if era in [ 2016, 2017 ]:
          for histogram_name_tmp in HISTOGRAM_COUNT_EXTENDED_MC:
            histogram_name = histogram_name_tmp
            if tH_idx >= 0:
              histogram_name += "_rwgt{}".format(tH_idx)
            if aux_bin:
              histogram_name += "_{}".format(aux_bin)
            histogram_names[histogram_name] = -1

  if key.startswith('/TTTo'):
    histogram_names_extend = []
    for histogram_name in histogram_names:
      if histogram_name == HISTOGRAM_COUNT:
        continue
      for topPtPrefix in [ "TOP16011", "Linear", "Quadratic", "HighPt" ]:
        for topPtSuffix in [ "TopPtRwgtSF", "TopPtRwgtSFSquared" ]:
          histogram_names_extend.append('{}{}{}'.format(histogram_name, topPtPrefix, topPtSuffix))
    for histogram_name in histogram_names_extend:
        histogram_names[histogram_name] = -1

  count_histograms_process = count_histograms[process_name] if process_name in count_histograms else {}

  indices = {}
  lhe_set = ''
  lhe_set_tried = False
  lhe_correct_binning = True
  nof_reweighting_weights = 0
  reweighting_tried = not is_rwgt
  nof_PSweights = 0
  PS_tried = False
  blacklist = []
  for entry in entries_valid:

    subentries = get_dir_entries(entry)
    subentry_files = filter(lambda path: path.isfile(), subentries)
    for subentry_file in subentry_files:
      index_entry = {
        HISTOGRAM_COUNT_KEY : {},
        TREE_COUNT_KEY      : -1,
        FSIZE_KEY           : -1,
        BRANCH_NAMES_KEY    : [],
      }

      digit_match = digit_regex.search(subentry_file.basename)
      if not digit_match:
        continue
      matched_idx = int(digit_match.group('i'))
      if file_idx > 0 and matched_idx != file_idx:
        logging.debug("Skipping file {path}".format(path = subentry_file.name))
        continue

      if subentry_file.name in lost_ntuples:
        logging.info("Excluding file {} as it has been found in the list of lost Ntuples".format(subentry_file.name))
        blacklist.append(matched_idx)
        continue

      if subentry_file.size == 0:
        logging.debug("File {path} has a file size of 0".format(path = subentry_file.name))
        filetracker.zero_file_size.append(subentry_file.name)
        continue
      index_entry[FSIZE_KEY] = subentry_file.size

      logging.debug("Opening file {path}".format(path = subentry_file.name))
      root_file = ROOT.TFile.Open(subentry_file.name, "read")
      if not root_file:
        logging.warning("Could not open {path}".format(path = subentry_file.name))
        filetracker.corrupted_files.append((subentry_file.name))
        continue
      if root_file.IsZombie():
        logging.warning("File {path} is a zombie".format(path = subentry_file.name))
        root_file.Close()
        del root_file
        filetracker.zombie_files.append(subentry_file.name)
        continue

      if EVENTS_TREE not in root_file.GetListOfKeys():
        raise ValueError("Tree of the name {tree} is not in file {path}".format(
          tree = EVENTS_TREE,
          path = subentry_file.name,
        ))
      tree = root_file.Get(EVENTS_TREE)
      if not tree:
        raise ValueError("Could not find tree of the name {tree} in file {path}".format(
          tree = check_every_event,
          path = subentry_file.name,
        ))
      index_entry[TREE_COUNT_KEY] = tree.GetEntries()
      index_entry[BRANCH_NAMES_KEY] = [ branch.GetName() for branch in tree.GetListOfBranches() ]

      if not PS_tried and not is_data:
        if BRANCH_NPSWEIGHTS in index_entry[BRANCH_NAMES_KEY]:
          nof_psweights_br = array.array('I', [0])
          tree.SetBranchAddress(BRANCH_NPSWEIGHTS, nof_psweights_br)
          tree.GetEntry(0)
          nof_PSweights = nof_psweights_br[0]
        PS_tried = True

      if not reweighting_tried:
        if BRANCH_NLHEREWEIGHTINGWEIGHT in index_entry[BRANCH_NAMES_KEY]:
          nof_reweighting_weights_br = array.array('I', [0])
          tree.SetBranchAddress(BRANCH_NLHEREWEIGHTINGWEIGHT, nof_reweighting_weights_br)
          tree.GetEntry(0)
          nof_reweighting_weights = nof_reweighting_weights_br[0]
        reweighting_tried = True

      if check_every_event:
        logging.info("Inspecting file {path} for corruption".format(path = subentry_file.name))
        no_errors_midfile = True
        for i in range(0, index_entry[TREE_COUNT_KEY]):
          nof_bytes_read = tree.GetEntry(i)
          if nof_bytes_read < 0:
            filetracker.corrupted_files.append(subentry_file.name)
            logging.debug("File {path} seems to be corrupted starting from event {idx}".format(
              path = subentry_file.name,
              idx  = i,
            ))
            no_errors_midfile = False
            break
        if not no_errors_midfile:
          # closing the ttree and file
          root_file.Close()
          del tree
          del root_file
          continue

      root_file_keys = [ root_file_key.GetName() for root_file_key in root_file.GetListOfKeys() ]
      has_any_histograms = any(HISTOGRAM_COUNT in histogram_name for histogram_name in histogram_names)
      if has_any_histograms:
        for histogram_name in histogram_names:
          if histogram_name not in root_file_keys:
            if HISTOGRAM_COUNT not in histogram_name:
              continue
            if ((HISTOGRAM_COUNTWEIGHTED_PSWEIGHT in histogram_name or
                 HISTOGRAM_COUNTWEIGHTED_PSWEIGHT_FULL in histogram_name) and nof_PSweights != 4):
              continue
            if histogram_name not in count_histograms_process:
              logging.warning("Histogram of the name {histogram_name} is not in file {path}".format(
                histogram_name = histogram_name,
                path           = subentry_file.name,
              ))
          else:
            histogram = root_file.Get(histogram_name)
            if not histogram:
              raise ValueError("Could not find histogram of the name {histogram_name} in file {path}".format(
                histogram_name = histogram_name,
                path           = subentry_file.name,
              ))
            nBins = histogram.GetNbinsX()
            if nBins not in [ 8, 9 ] and histogram_name in lheScaleArr:
              logging.warning("Expected 8 or 9 bins but found {nBins} bins in {histogram_name}".format(
                nBins          = nBins,
                histogram_name = histogram_name,
              ))
              lhe_correct_binning = False
              continue
            index_entry[HISTOGRAM_COUNT_KEY][histogram_name] = [
              histogram.GetBinContent(idxBin) for idxBin in range(1, nBins + 1)
            ]
            if histogram_names[histogram_name] < 0:
              histogram_names[histogram_name] = nBins
            else:
              if histogram_names[histogram_name] != nBins:
                raise RuntimeError(
                  "Expected to find {nBins_expected} bins in histogram {histogram_name} from file {path} "
                  "but got {nBins_actual} bins instead".format(
                    nBins_expected = histogram_names[histogram_name],
                    histogram_name = histogram_name,
                    path           = subentry_file.name,
                    nBins_actual   = nBins,
                  )
                )
            del histogram

        for histogram_name in index_entry[HISTOGRAM_COUNT_KEY]:
          if histogram_names[histogram_name] == 8 and (histogram_name in lheScaleArr or is_njet):
            histogram_name_nolhe = histogram_name.replace('LHEWeightScale', '')
            assert(histogram_name_nolhe != histogram_name)
            assert(histogram_name_nolhe in index_entry[HISTOGRAM_COUNT_KEY])
            # use nominal weight as the 5th LHE scale weight
            index_entry[HISTOGRAM_COUNT_KEY][histogram_name].insert(4, index_entry[HISTOGRAM_COUNT_KEY][histogram_name_nolhe][0])

      # this was probably a success: record the results
      indices[matched_idx] = copy.deepcopy(index_entry)
      logging.debug(
        "Found {nof_tree_events} tree events in file {filename}".format(
          nof_tree_events = index_entry[TREE_COUNT_KEY],
          filename        = subentry_file,
        )
      )

      if not is_data and not lhe_set_tried:
        lhe_set = get_lhe_set(tree)
        lhe_set_tried = True

      root_file.Close()
      del tree
      del root_file

  if not indices:
    logging.debug("Path {path} contains no ROOT files".format(path = path_obj.name))
    return

  for matched_idx in indices:
    for histogram_name in indices[matched_idx][HISTOGRAM_COUNT_KEY]:
      if histogram_names[histogram_name] == 8 and (histogram_name in lheScaleArr or is_njet):
        histogram_names[histogram_name] += 1

  logging.debug("Found total {nof_tree_events} tree events in {nof_files} files in "
                "{path} for entry {key}".format(
    nof_tree_events = sum([index_entry[TREE_COUNT_KEY]      for index_entry in indices.values()]),
    nof_files       = len(indices.keys()),
    path            = path_obj.name,
    key             = key,
  ))

  overlap_with_triggers = []
  if not meta_dict[key]['located']:
    missing_from_superset = [] if not missing_branches else get_missing_from_superset(indices)
    overlap_with_triggers = list(triggerTable.triggers_flat & set(missing_from_superset))
    if overlap_with_triggers:
      logging.error(
        "Found an overlap b/w the list of required triggers and the list of missing branches in "
        "sample {}: {}".format(meta_dict[key]['process_name_specific'], ', '.join(overlap_with_triggers))
       )
    meta_dict[key]['triggers']                        = get_triggers(
      meta_dict[key]['process_name_specific'], is_data, era
    )
    meta_dict[key]['missing_hlt_paths']               = get_missing_hlt_paths(
      get_triggers('', False, era), indices, triggerTable.triggers_all
    )
    meta_dict[key]['hlt_paths']                       = get_hlt_paths(indices) if is_data else []
    meta_dict[key]['genWeight']                       = not is_data
    meta_dict[key]['type']                            = 'data' if is_data else 'mc'
    meta_dict[key]['reHLT']                           = True
    meta_dict[key]['located']                         = True
    meta_dict[key]['has_LHE']                         = False if is_data else (lhe_correct_binning and has_LHE(indices))
    meta_dict[key]['nof_PSweights']                   = nof_PSweights
    meta_dict[key]['missing_from_superset']           = missing_from_superset
    meta_dict[key]['missing_completely']              = overlap_with_triggers
    meta_dict[key]['histogram_names']                 = histogram_names
    meta_dict[key]['LHE_set']                         = lhe_set
    meta_dict[key]['nof_reweighting']                 = nof_reweighting_weights
  meta_dict[key]['paths'].append(
    PathEntry(path_obj.name, indices, histogram_names, blacklist)
  )

  return

def traverse_double(meta_dict, path_obj, key, check_every_event, missing_branches,
                    filetracker, file_idx, era, triggerTable, count_histograms, lost_ntuples):
  ''' Assume that the name of the following subdirectories are the CRAB job IDs
      The tree structure inside those directories should be the same as described in
      traverse_single()
      Therefore, we loop over the CRAB job IDs and pass each subfolder to traverse_single()
  :param meta_dict:         Meta-dictionary
  :param path_obj:          Contains meta-information about a path
  :param key:               Key to the meta-dictionary the entry of which will be updated
  :param check_every_event: Loop over all events for error checking purposes
  :param missing_branches:  Find missing branches from the superset of branches in a sample
  :param filetracker:       An instance of FileTracker() for logging broken files
  :param file_idx:          Index of the corrupted file
  :param triggerTable:      Trigger instance containing a list of required triggers for the era
  :param count_histograms   Externally provided event sums
  :return: None
  '''
  logging.info("Double-traversing {path}".format(path = path_obj.name))
  entries = get_dir_entries(path_obj)
  for entry in entries:
    traverse_single(
      meta_dict, entry, key, check_every_event, missing_branches,
      filetracker, file_idx, era, triggerTable, count_histograms, lost_ntuples
    )
  return

def get_path_info(path):
  return path_info(path)

def obtain_paths(input_path):
  paths = []
  if input_path:
    # check if the input path is a path or a file
    path = input_path[0]
    if os.path.isfile(path):
      with open(path, 'r') as f:
        for line in f:
          line_stripped = line.rstrip('\n').rstrip(os.path.sep)
          if line_stripped:
            paths.append(line_stripped.split()[0])
    else:
      paths = input_path
  else:
    paths = input_path
  return paths

def round_sign(x, sign_digits = 6):
  return round(x, max(int(abs(math.floor(math.log10(x)))) + sign_digits, 0))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(
    formatter_class = lambda prog: SmartFormatter(prog, max_help_position = 35))
  parser.add_argument('-m', '--meta-dictionary', dest = 'meta_dictionary', metavar = 'file',
                      required = True, type = str,
                      help = 'R|Path to the meta dictionary')
  parser.add_argument('-N', '--output-dict-name', dest = 'output_dict_name', metavar = 'name',
                      type = str, default = 'sample',
                      help = 'R|Name of the output dictionary')
  parser.add_argument('-p', '--path', dest = 'path', metavar = 'directory', required = True,
                      type = str, nargs = '+',
                      help = 'R|List of full paths to the directories to scan')
  parser.add_argument('-e', '--exclude-path', dest = 'exclude_path', metavar = 'directory',
                      required = False, type = str, nargs = '+', default = [])
  parser.add_argument('-d', '--depth', dest = 'depth', metavar = 'number', required = False,
                      type = int, default = -1,
                      help = 'R|How many directory levels to traverse (default: all)')
  parser.add_argument('-f', '--filter', dest = 'filter', metavar = 'name', required = False,
                      type = str, default = '.*',
                      help = 'R|Regular expression for selecting only specific samples')
  parser.add_argument('-o', '--output-directory', dest = 'output_directory', metavar = 'path',
                      type = str, default = '.',
                      help = 'R|Output directory')
  parser.add_argument('-g', '--generate-python', dest = 'generate_python', metavar = 'name',
                      type = str, default = 'dict.py',
                      help = 'R|File name of the new python dictionary')
  parser.add_argument('-c', '--check-every-event', dest = 'check_every_event', metavar = 'name',
                      type = str, default = "", required = False,
                      help = 'R|Supply TTree name to check every event (NB! Extremely slow!)')
  parser.add_argument('-z', '--save-zombies', dest = 'save_zombies', metavar = 'save_zombies',
                      type = str, default = '',
                      help = 'R|Save the list of zombie files')
  parser.add_argument('-Z', '--save-zeroes', dest = 'save_zeroes', metavar = 'save_zeros',
                      type = str, default = '',
                      help = 'R|Save the list of files with zero file size')
  parser.add_argument('-C', '--save-corrupted', dest = 'save_corrupted', metavar = 'save_corrupted',
                      type = str, default = '',
                      help = 'R|Save the list of corrupted files')
  parser.add_argument('-j', '--file-idx', dest = 'file_idx', metavar = 'number', type = int,
                      default = -1,
                      help = 'R|Check files at specified index (default: all files)')
  parser.add_argument('-s', '--skip-header', dest = 'skip_header', action = 'store_true',
                      default = False,
                      help = 'R|Skip dictionary definitions in the output')
  parser.add_argument('-J', '--generate-jobs', dest = 'generate_jobs', metavar = 'generate_jobs',
                      type = str, default = '', required = False,
                      help = 'R|Generate SLURM jobs instead of running locally')
  parser.add_argument('-E', '--era', dest = 'era', metavar = 'era', type = int, default = -1,
                      required = True, choices = (2016,2017,2018),
                      help = 'R|Era of the samples')
  parser.add_argument('-M', '--find-missing-branches', dest = 'missing_branches', action = 'store_true',
                      default = False,
                      help = 'R|Find missing branches from the superset of branches in a sample')
  parser.add_argument('-x', '--clean', dest = 'clean', action = 'store_true', default = False,
                      help = 'R|Clean the temporary SLURM directory specified by -J')
  parser.add_argument('-F', '--force', dest = 'force', action = 'store_true', default = False,
                      help = 'R|Force the creation of missing directories')
  parser.add_argument('-q', '--count-histograms', dest = 'count_histograms', metavar = 'file',
                      type = str, default = '', required = False,
                      help = 'R|A ROOT file with histograms storing the event sums')
  parser.add_argument('-l', '--lost', dest = 'lost', metavar = 'file', type = str, default = '', required = False,
                      help = 'R|File containing list of lost Ntuples')
  parser.add_argument('-v', '--verbose', dest = 'verbose', action = 'store_true', default = False,
                      help = 'R|Enable verbose printout')
  args = parser.parse_args()

  if args.verbose:
    logging.getLogger().setLevel(logging.DEBUG)

  if not os.path.isdir(args.output_directory):
    if not args.force:
      raise parser.error("Directory %s does not exist (use -F/--force to create it)" % args.output_directory)
    else:
      os.makedirs(args.output_directory)

  if args.generate_jobs and not os.path.isdir(args.generate_jobs):
    if not args.force:
      raise parser.error("Directory %s does not exist" % args.generate_jobs)
    else:
      os.makedirs(args.generate_jobs)

  if (args.file_idx < 0 or not args.filter) and args.check_every_event:
    raise parser.error("Checking all files for data corruption is extremely slow! "
                       "Use -j/--file-idx and -f/--filter options on specific files if "
                       "you really need to check such files!")

  if args.save_corrupted and not args.check_every_event:
    logging.warning("The flag -C/--save-corrupted has no effect w/o -c/--check-every-event option")

  lost_ntuples = []
  if args.lost:
    with open(args.lost, 'r') as lost_file:
      for line in lost_file:
        lost_ntuples.append(line.strip())

  count_histograms = {}
  count_histogram_input = args.count_histograms
  if count_histogram_input:
    if not os.path.isfile(count_histogram_input):
      raise RuntimeError("No such file: %s" % count_histogram_input)
    count_histogram_file = ROOT.TFile.Open(count_histogram_input)
    count_process_names = [ key.GetName() for key in count_histogram_file.GetListOfKeys() ]
    nof_count_processes = len(count_process_names)
    logging.info("Reading event counts of {} processes from {}".format(nof_count_processes, count_histogram_input))
    for count_process_idx, count_process_name in enumerate(count_process_names):
      logging.info("Reading events counts of process {} ({}/{})".format(
        count_process_name, count_process_idx + 1, nof_count_processes
      ))
      count_histogram_dir = count_histogram_file.Get(count_process_name)
      count_histogram_names = [ key.GetName() for key in count_histogram_dir.GetListOfKeys() ]
      process_content = collections.OrderedDict()
      for count_histogram_name in count_histogram_names:
        count_histogram = count_histogram_dir.Get(count_histogram_name)
        count_histogram.SetDirectory(0)
        process_content[count_histogram_name] = [
          count_histogram.GetBinContent(bin_idx) for bin_idx in range(1, count_histogram.GetNbinsX() + 1)
        ]
      for count_histogram_name in process_content:
        if len(process_content[count_histogram_name]) == 8 and \
           (count_histogram_name in LHESCALEARR or get_is_njet(count_process_name)):
          histogram_name_nolhe = count_histogram_name.replace('LHEWeightScale', '')
          assert(histogram_name_nolhe != count_histogram_name)
          assert(histogram_name_nolhe in process_content)
          # use nominal weight as the 5th LHE scale weight
          process_content[count_histogram_name].insert(4, process_content[histogram_name_nolhe][0])
          logging.info("Got 8 weights in {} -> added the nominal weight to make it 9".format(count_histogram_name))
      count_histograms[count_process_name] = process_content
    count_histogram_file.Close()

  filetracker = FileTracker()

  paths_unchecked = obtain_paths(args.path)
  excluded_paths  = obtain_paths(args.exclude_path)

  # check if the given paths actually exist
  paths = [get_path_info(path) for path in paths_unchecked]
  invalid_paths = filter(lambda path: not path.isdir(), paths)
  if invalid_paths:
    raise parser.error('The following paths do not exist: %s' % ', '.join(invalid_paths))
  for path in paths:
    path.sparent_depth = path.depth

  # load the dictionaries
  meta_dict = load_dict(args.meta_dictionary, "meta_dictionary")
  sum_events = load_dict(args.meta_dictionary, "sum_events")
  process_names = { entry['process_name_specific']         : dbs_name for dbs_name, entry in meta_dict.items() }
  crab_strings  = { os.path.basename(entry['crab_string']) : dbs_name for dbs_name, entry in meta_dict.items() if entry['crab_string'] != "" }
  crab_strings_single = [ os.path.basename(entry['crab_string']) for entry in meta_dict.values() if "/" in entry['crab_string'] ]
  for key, entry in meta_dict.items():
    entry['located'] = False
  for key_arr in sum_events:
    for key in key_arr:
      if key not in process_names:
        raise ValueError("No such key in meta_dictionary: %s" % key)

  # set up the regex object
  name_regex = re.compile(args.filter)
  triggerTable = Triggers(str(args.era))

  # process the directory structure of each path
  paths_to_traverse = {}
  while paths:
    path = paths.pop(0)
    if path in excluded_paths:
      logging.info("Skipping path {path} since it is in the exclusion list".format(path = path.name))
      continue
    logging.debug("Considering path {path}".format(path = path.name))
    if args.depth > 0 and (path.depth - path.sparent_depth) >= args.depth:
      continue
    if path.basename in process_names:
      expected_key = meta_dict[process_names[path.basename]]['process_name_specific']
      is_match = name_regex.match(expected_key)
      if is_match:
        if args.generate_jobs:
          if expected_key not in paths_to_traverse:
            paths_to_traverse[expected_key] = []
          paths_to_traverse[expected_key].append(path.name)
        else:
          traverse_single(
            meta_dict, path, process_names[path.basename],
            args.check_every_event, args.missing_branches, filetracker, args.file_idx, args.era,
            triggerTable, count_histograms, lost_ntuples
          )
    elif path.basename in crab_strings:
      expected_key = meta_dict[crab_strings[path.basename]]['process_name_specific']
      is_match = name_regex.match(expected_key)
      if is_match:
        if args.generate_jobs:
          if expected_key not in paths_to_traverse:
            paths_to_traverse[expected_key] = []
          paths_to_traverse[expected_key].append(path.name)
        else:
          if path.basename in crab_strings_single:
            traverse_double(
              meta_dict, path, crab_strings[path.basename],
              args.check_every_event, args.missing_branches, filetracker, args.file_idx, args.era,
              triggerTable, count_histograms, lost_ntuples
            )
            traverse_single(
              meta_dict, path, crab_strings[path.basename],
              args.check_every_event, args.missing_branches, filetracker, args.file_idx, args.era,
              triggerTable, count_histograms, lost_ntuples
            )
          else:
            traverse_double(
              meta_dict, path, crab_strings[path.basename],
              args.check_every_event, args.missing_branches, filetracker, args.file_idx, args.era,
              triggerTable, count_histograms, lost_ntuples
            )
    else:
      entries = get_dir_entries(path)
      entries_dirs = filter(
        lambda entry: entry.isdir() and os.path.basename(entry.name) not in ["failed", "log"] and \
                      not any(map(
                        lambda path_to_traverse: entry.name.startswith(path_to_traverse),
                        list(itertools.chain.from_iterable(paths_to_traverse.values()))
                      )),
        entries
      )
      for entry in entries_dirs:
        if entry not in paths:
          if entry.sparent_depth < 0:
            entry.sparent_depth = path.sparent_depth
          logging.debug(
            "Adding entry {entry} ({sparent_depth}/{depth})".format(
              entry         = entry.name,
              sparent_depth = entry.sparent_depth,
              depth         = entry.depth,
            )
          )
          paths.append(entry)

  output = jinja2.Template(header_str).render(
    command   = ' '.join([os.path.basename(__file__)] + sys.argv[1:]),
    date      = '{date:%Y-%m-%d %H:%M:%S}'.format(date = datetime.datetime.now()),
    dict_name = args.output_dict_name,
  ) if not args.skip_header else ''

  if args.generate_jobs:
    # divide the paths according to their phase space
    # first, let's check if we've complete phase space together
    path_arrs = []
    for key_arr in sum_events:
      key_intersection = set(key_arr) & set(paths_to_traverse.keys())
      if len(key_intersection) == len(key_arr):
        # we've complete phase space together
        path_arr = []
        for key in key_arr:
          path_arr.extend(paths_to_traverse[key])
          del paths_to_traverse[key]
        path_arrs.append(path_arr)
      elif len(key_intersection) == 0:
        # the phase space is completely missing
        pass
      else:
        raise ValueError("Incomplete phase space: %s (should be: %s)" % (
          ', '.join(key_intersection),
          ', '.join(key_arr)
        ))
    # set the remaining paths
    path_arrs.extend(paths_to_traverse.values())

    commands = {}
    to_cat = { 'dicts' : [], }
    for arg in vars(args):
      args_attr = getattr(args, arg)
      if args_attr:
        option_key = ''
        option_default = None
        for option in parser._optionals._actions:
          if option.dest == arg:
            option_key = option.option_strings[0]
            option_default = option.default
            break
        if not option_key:
          raise ValueError("Internal error: inconsistencies in ArgumentParser!")
        if args_attr != option_default:
          if type(args_attr) is not bool:
            if type(args_attr) is list:
              commands[option_key] = ' '.join(map(str, args_attr))
            else:
              commands[option_key] = str(args_attr)

    # copy the supplied CLI parameters over, modify them such that the intermediate files
    # would be stored in a different directory specified by args.generate_jobs and construct
    # the list of shell commands that will be submitted to SLURM system
    job_params = []
    for path_idx, path_arr in enumerate(path_arrs):
      commands_cp = copy.deepcopy(commands)
      commands_cp['-p'] = os.path.join(os.path.realpath(args.generate_jobs), ' '.join(path_arr))
      commands_cp['-g'] = os.path.join(os.path.realpath(args.generate_jobs), 'dict.py.%i' % path_idx)
      to_cat['dicts'].append(commands_cp['-g'])


      for key in ['-z', '-Z', '-C']:
        if key in commands_cp:
          commands_cp[key] = os.path.join(os.path.realpath(args.generate_jobs), '%s.%i' % (commands_cp[key], path_idx))
          if key not in to_cat:
            to_cat[key] = []
          to_cat[key].append(commands_cp[key])
      commands_cp['-m'] = os.path.join(os.getcwd(), commands_cp['-m'])
      commands_cp['-s'] = ''
      del commands_cp['-J']

      cmd = ' '.join(['python', sys.argv[0]] + [k + ' ' + v for k, v in commands_cp.items()])
      sh = jinja2.Template(sh_str).render(cmd = cmd)
      sh_file = os.path.join(args.generate_jobs, 'job_%i.sh' % path_idx)
      with open(sh_file, 'w') as f:
        f.write(sh)
      log_file = os.path.join(args.generate_jobs, 'log_%i.txt' % path_idx)
      job_params.append((log_file, sh_file))

    # submit the jobs
    submit_cmds = list(map(
      lambda job_param: 'sbatch --mem=1800M --partition=small --output=%s %s' % job_param,
      job_params
    ))
    squeue_codes = []
    for submit_cmd in submit_cmds:
      squeue_code = run_cmd(submit_cmd).split()[-1]
      squeue_codes.append(squeue_code)
      logging.info("Submitted sbatch job {jobId}".format(jobId = squeue_code))

    has_completed = not bool(squeue_codes)
    while not has_completed:
      squeue = run_cmd("squeue -j {jobIds} -h | wc -l".format(jobIds = ','.join(squeue_codes))).rstrip('\n')
      if squeue == '0':
        has_completed = True
      logging.debug("{nofJobs} job(s) still running...".format(nofJobs = squeue))
      time.sleep(5)
    logging.info("All jobs have been finished")

    # cat the dictionary
    for dict_file in to_cat['dicts']:
      if not os.path.exists(dict_file):
        raise ValueError("Missing temporary dictionary: %s" % dict_file)
      with open(dict_file, 'r') as f:
        output += '\n'.join(map(lambda line: line.rstrip('\n'), f.readlines()))

    # cat the faulty files
    for key in ['-z', '-Z', '-C']:
      if key not in to_cat:
        continue
      for faulty_list_file in to_cat[key]:
        if os.path.exists(faulty_list_file):
          with open(faulty_list_file, 'r') as f:
            lines = map(lambda line: line.rstrip('\n'), f.readlines())
            if key == '-z':
              filetracker.zombie_files.extend(lines)
            if key == '-Z':
              filetracker.zero_file_size.extend(lines)
            if key == '-C':
              filetracker.corrupted_files.extend(lines)

    if args.clean:
      shutil.rmtree(args.generate_jobs)
  else:
    # we need to post-process the meta dictionary
    for key, entry in meta_dict.items():
      if not name_regex.match(entry['process_name_specific']):
        continue
      if entry['located']:
        process_paths(meta_dict, key, count_histograms)

    for key, entry in meta_dict.items():
      if not name_regex.match(entry['process_name_specific']):
        continue
      if entry['located']:
        path_entries_arr = []
        for path_entry in meta_dict[key]['local_paths']:
          path_entries_arr.append(jinja2.Template(path_entry_str).render(
            path      = path_entry['path'],
            selection = path_entry['selection'],
            blacklist = path_entry['blacklist'], #TODO: format properly
          ))
        is_mc = meta_dict[key]['type'] == 'mc'
        missing_branches_template_filled = jinja2.Template(missing_branches_str).render(
          is_available     = args.missing_branches and not is_mc,
          missing_branches = sorted(meta_dict[key]['missing_from_superset'], key = lambda s: s.lower()),
        ).lstrip('\n')
        completely_missing_branches_template_filled = jinja2.Template(missing_branches_str).render(
          is_available     = args.missing_branches and not is_mc,
          missing_branches = sorted(meta_dict[key]['missing_completely'], key = lambda s: s.lower()),
        ).lstrip('\n')
        missing_hlt_paths_filled = jinja2.Template(missing_branches_str).render(
          is_available     = True,
          missing_branches = sorted(meta_dict[key]['missing_hlt_paths'], key = lambda s: s.lower()),
        ).lstrip('\n')
        hlt_paths_filled = jinja2.Template(missing_branches_str).render(
          is_available = not is_mc,
          missing_branches = sorted(meta_dict[key]['hlt_paths'], key = lambda  s: s.lower()),
        ).lstrip('\n')
        output += jinja2.Template(dictionary_entry_str).render(
          dict_name                       = args.output_dict_name,
          dbs_name                        = key,
          sample_type                     = meta_dict[key]['type'],
          sample_category                 = meta_dict[key]['sample_category'],
          process_name_specific           = meta_dict[key]['process_name_specific'],
          nof_files                       = meta_dict[key]['nof_files'],
          nof_events                      = meta_dict[key]['nof_events'] if is_mc else {},
          nof_tree_events                 = meta_dict[key]['nof_tree_events'],
          nof_db_events                   = meta_dict[key]['nof_db_events'],
          nof_db_files                    = meta_dict[key]['nof_db_files'],
          fsize_db                        = meta_dict[key]['fsize_db'],
          fsize_db_human                  = human_size(meta_dict[key]['fsize_db']),
          avg_fsize_db_human              = human_size(float(meta_dict[key]['fsize_db']) / meta_dict[key]['nof_db_files']),
          fsize_local                     = meta_dict[key]['fsize_local'],
          fsize_local_human               = human_size(meta_dict[key]['fsize_local']),
          avg_fsize_local_human           = human_size(float(meta_dict[key]['fsize_local']) / meta_dict[key]['nof_files']),
          use_it                          = meta_dict[key]['use_it'],
          xsection                        = round_sign(meta_dict[key]['xsection'], 6) if is_mc else None,
          genWeight                       = meta_dict[key]['genWeight'],
          triggers                        = meta_dict[key]['triggers'],
          has_LHE                         = meta_dict[key]['has_LHE'],
          nof_PSweights                   = meta_dict[key]['nof_PSweights'],
          LHE_set                         = meta_dict[key]['LHE_set'],
          nof_reweighting                 = meta_dict[key]['nof_reweighting'],
          missing_from_superset           = missing_branches_template_filled,
          missing_completely              = completely_missing_branches_template_filled,
          missing_hlt_paths               = missing_hlt_paths_filled,
          hlt_paths                       = hlt_paths_filled,
          paths                           = '\n'.join(path_entries_arr),
        ) + '\n\n'
      else:
        logging.warning("Could not locate paths for {key}".format(key = key))

  output += jinja2.Template(dictionary_sum_events_str).render(
    dict_name    = args.output_dict_name,
    sample_lists = sum_events,
  ) + '\n\n'

  dictionary_path = os.path.join(args.output_directory, args.generate_python)
  with open(dictionary_path, 'w') as f:
    f.write(output)
  logging.info("Wrote the dictionary to {path}".format(path = dictionary_path))

  if filetracker.zero_file_size:
    logging.info("The following files had file size of zero:\n{zero_fs}".format(
      zero_fs = '\n'.join(filetracker.zero_file_size),
    ))
    if args.save_zeroes:
      zeroes_path = os.path.join(args.output_directory, args.save_zeroes)
      with open(zeroes_path, 'w') as f:
        f.write('\n'.join(filetracker.zero_file_size) + '\n')
      logging.info("Saved the list of files with zero file size to {path}".format(
        path = zeroes_path,
      ))
  if filetracker.zombie_files:
    logging.info("The following files were zombies:\n{zombies}".format(
      zombies = '\n'.join(filetracker.zombie_files),
    ))
    if args.save_zombies:
      zombies_path = os.path.join(args.output_directory, args.save_zombies)
      with open(zombies_path, 'w') as f:
        f.write('\n'.join(filetracker.zombie_files) + '\n')
      logging.info("Saved the list of zombie files to {path}".format(path = zombies_path))
  if filetracker.corrupted_files:
    logging.info("The following files were corrupted:\n{corrupted}".format(
      corrupted = '\n'.join(filetracker.corrupted_files),
    ))
    if args.save_corrupted:
      corrupted_path = os.path.join(args.output_directory, args.save_corrupted)
      with open(corrupted_path, 'r') as f:
        f.write('\n'.join(filetracker.corrupted_files) + '\n')
      logging.info("Saved the list of corrupted files to {path}".format(path = corrupted_path))
