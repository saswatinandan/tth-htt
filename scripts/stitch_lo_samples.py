#!/usr/bin/env python

from tthAnalysis.HiggsToTauTau.safe_root import ROOT
from tthAnalysis.HiggsToTauTau.common import SmartFormatter, load_samples

import array
import itertools
import copy
import os
import argparse

def is_valid_event_type(nof_key):
  return 'PSWeight' not in nof_key and '_' not in nof_key

def copy_nof_events(sample_entry):
  return {
      nof_key : sample_entry['nof_events'][nof_key] \
        for nof_key in sample_entry['nof_events'] \
        if is_valid_event_type(nof_key)
    }

def comp_weights_1(f, samples, samples_to_stitch, split_var, apply_sf = True):
  inclusive_samples  = samples_to_stitch['inclusive']['samples']
  inclusive_binning  = samples_to_stitch['inclusive'][split_var]

  split_dict       = samples_to_stitch['exclusive'][split_var]
  split_binning    = [ sample['value'] for sample in split_dict ]
  complete_binning = list(sorted(list(
    map(float, set(inclusive_binning) | set(list(itertools.chain.from_iterable(split_binning))))
  )))

  inclusive_xs = -1
  for sample_key, sample_entry in samples.items():
    if sample_key == 'sum_events': continue
    if sample_entry['process_name_specific'] == inclusive_samples[0]:
      inclusive_xs = sample_entry['xsection']
  assert(inclusive_xs > 0)

  # sum the inclusive nof events
  inclusive_nof_events = {}
  for sample_key, sample_entry in samples.items():
    if sample_key == 'sum_events': continue
    if sample_entry['process_name_specific'] in inclusive_samples:
      if not inclusive_nof_events:
        inclusive_nof_events = copy_nof_events(sample_entry)
      else:
        nof_events_keys = set(nof_key for nof_key in sample_entry['nof_events'] if is_valid_event_type(nof_key))
        assert(nof_events_keys == set(inclusive_nof_events))
        for nof_key, nof_arr in sample_entry['nof_events'].items():
          if not is_valid_event_type(nof_key):
            continue
          assert(len(nof_arr) == len(inclusive_nof_events[nof_key]))
          for idx, nof in enumerate(nof_arr):
            assert(nof > 0)
            inclusive_nof_events[nof_key][idx] += nof

  # sum the binned nof events
  for binned_sample in split_dict:
    nof_events = {}
    xs = -1
    for sample_key, sample_entry in samples.items():
      if sample_key == 'sum_events': continue
      if sample_entry['process_name_specific'] in binned_sample['samples']:
        if not nof_events:
          nof_events = copy_nof_events(sample_entry)
          inclusive_nof_events_type = set(event_type for event_type in inclusive_nof_events.keys() if is_valid_event_type(event_type))
          nof_events_type = set(event_type for event_type in nof_events.keys() if is_valid_event_type(event_type))
          assert(inclusive_nof_events_type == nof_events_type)
        else:
          nof_events_keys = set(nof_key for nof_key in sample_entry['nof_events'] if is_valid_event_type(nof_key))
          assert(nof_events_keys == set(nof_events.keys()))
          for nof_key, nof_arr in sample_entry['nof_events'].items():
            if not is_valid_event_type(nof_key):
              continue
            assert(len(nof_arr) == len(nof_events[nof_key]))
            for idx, nof in enumerate(nof_arr):
              assert(nof > 0)
              nof_events[nof_key][idx] += nof
        if xs < 0:
          xs = sample_entry['xsection']
    assert(xs > 0)
    binned_sample['xsection'] = xs
    binned_sample['nof_events'] = nof_events

    lumis = {}
    for nof_key, nof_arr in binned_sample['nof_events'].items():
      lumis[nof_key] = list(map(lambda nof: nof / binned_sample['xsection'], nof_arr))
    binned_sample['lumis'] = lumis

  # compute integrated luminosities for the inclusive sample
  inclusive_lumis = {}
  for nof_key, nof_arr in inclusive_nof_events.items():
    inclusive_lumis[nof_key] = list(map(lambda nof: nof / inclusive_xs, nof_arr))

  # decide on the bin indices
  idxs_split_sample = []
  for binned_sample in split_dict:
    binned_idx = complete_binning.index(binned_sample['value'][0]) + 1
    binned_sample['idx'] = binned_idx
    idxs_split_sample.append(binned_idx)

  for inclusive_sample in inclusive_samples:
    if inclusive_sample not in [ key.GetName() for key in f.GetListOfKeys() ]:
      histogram_dir_root = f.mkdir(inclusive_sample)
    else:
      histogram_dir_root = f.Get(inclusive_sample)
    if split_var not in [ key.GetName() for key in histogram_dir_root.GetListOfKeys() ]:
      histogram_dir = histogram_dir_root.mkdir(split_var)
    else:
      histogram_dir = histogram_dir_root.Get(split_var)
    histogram_dir.cd()

    for nof_key, lumi_arr in inclusive_lumis.items():
      if not is_valid_event_type(nof_key):
        continue
      for idx, lumi_incl in enumerate(lumi_arr):
        histogram_name = '%s_%d' % (nof_key, idx)
        binning = array.array('f', complete_binning)

        histogram = ROOT.TH1D(histogram_name, histogram_name, len(binning) - 1, binning)
        histogram.SetDirectory(histogram_dir)
        histogram.SetXTitle(split_var)

        for split_idx in range(1, len(binning)):
          if split_idx in idxs_split_sample:
            for binned_sample in split_dict:
              if split_idx == binned_sample['idx']:
                lumi_split = binned_sample['lumis'][nof_key][idx]
                if binning[split_idx] > inclusive_binning[-1] or \
                   binning[split_idx] < inclusive_binning[0]:
                  lumi_incl_calc = 0.
                else:
                  lumi_incl_calc = lumi_incl
                weight = lumi_incl_calc / (lumi_incl_calc + lumi_split)
                histogram.SetBinContent(split_idx, weight)
          else:
            histogram.SetBinContent(split_idx, 1.)
          histogram.GetXaxis().SetBinLabel(
            split_idx, '%.0f <= %s < %.0f' % (
              complete_binning[split_idx - 1], split_var, complete_binning[split_idx]
            )
          )

        histogram.Write()

  for binned_sample in split_dict:
    for sample_name in binned_sample['samples']:
      if sample_name not in [ key.GetName() for key in f.GetListOfKeys() ]:
        histogram_dir_root = f.mkdir(sample_name)
      else:
        histogram_dir_root = f.Get(sample_name)
      if split_var not in [ key.GetName() for key in histogram_dir_root.GetListOfKeys() ]:
        histogram_dir = histogram_dir_root.mkdir(split_var)
      else:
        histogram_dir = histogram_dir_root.Get(split_var)
      histogram_dir.cd()

      for nof_key, lumi_arr in inclusive_lumis.items():
        if not is_valid_event_type(nof_key):
          continue
        for idx, lumi_incl in enumerate(lumi_arr):

          histogram_name = '%s_%d' % (nof_key, idx)
          binning = array.array('f', complete_binning)

          histogram = ROOT.TH1D(histogram_name, histogram_name, len(binning) - 1, binning)
          histogram.SetDirectory(histogram_dir)
          histogram.SetXTitle(split_var)

          for split_idx in range(1, len(binning)):
            if split_idx == binned_sample['idx']:
              lumi_split = binned_sample['lumis'][nof_key][idx]
              if binning[split_idx] > inclusive_binning[-1] or \
                 binning[split_idx] < inclusive_binning[0]:
                weight = 1.
              else:
                if apply_sf:
                  weight = lumi_split / (lumi_incl + lumi_split)
                else:
                  weight = lumi_incl / (lumi_incl + lumi_split)
              assert(weight >= 0.)
              histogram.SetBinContent(split_idx, weight)
            else:
              histogram.SetBinContent(split_idx, 0.)
            histogram.GetXaxis().SetBinLabel(
              split_idx, '%.0f <= %s < %.0f' % (
                complete_binning[split_idx - 1], split_var, complete_binning[split_idx]
              )
            )

          histogram.Write()

def comp_weights_2(f, samples, samples_to_stitch, split_var_1, split_var_2, apply_sf = True):
  inclusive_samples   = samples_to_stitch['inclusive']['samples']
  inclusive_binning_1 = samples_to_stitch['inclusive'][split_var_1]
  inclusive_binning_2 = samples_to_stitch['inclusive'][split_var_2]

  split_dict_1 = samples_to_stitch['exclusive'][split_var_1]
  split_dict_2 = samples_to_stitch['exclusive'][split_var_2]
  split_binning_1 = [ sample['value'] for sample in split_dict_1 ]
  split_binning_2 = [ sample['value'] for sample in split_dict_2 ]
  complete_binning_1 = list(sorted(list(
    map(float, set(inclusive_binning_1) | set(list(itertools.chain.from_iterable(split_binning_1))))
  )))
  complete_binning_2 = list(sorted(list(
    map(float, set(inclusive_binning_2) | set(list(itertools.chain.from_iterable(split_binning_2))))
  )))

  inclusive_xs = -1
  for sample_key, sample_entry in samples.items():
    if sample_key == 'sum_events': continue
    if sample_entry['process_name_specific'] == inclusive_samples[0]:
      inclusive_xs = sample_entry['xsection']
  assert(inclusive_xs > 0)

  # sum the inclusive nof events
  inclusive_nof_events = {}
  for sample_key, sample_entry in samples.items():
    if sample_key == 'sum_events': continue
    if sample_entry['process_name_specific'] in inclusive_samples:
      if not inclusive_nof_events:
        inclusive_nof_events = copy_nof_events(sample_entry)
      else:
        nof_events_keys = set(nof_key for nof_key in sample_entry['nof_events'] if is_valid_event_type(nof_key))
        assert(nof_events_keys == set(inclusive_nof_events))
        for nof_key, nof_arr in sample_entry['nof_events'].items():
          if not is_valid_event_type(nof_key):
            continue
          assert(len(nof_arr) == len(inclusive_nof_events[nof_key]))
          for idx, nof in enumerate(nof_arr):
            assert(nof > 0)
            inclusive_nof_events[nof_key][idx] += nof

  # sum the binned nof events
  for binned_sample in split_dict_1:
    nof_events = {}
    xs = -1
    for sample_key, sample_entry in samples.items():
      if sample_key == 'sum_events': continue
      if sample_entry['process_name_specific'] in binned_sample['samples']:
        if not nof_events:
          nof_events = copy_nof_events(sample_entry)
          inclusive_nof_events_type = set(event_type for event_type in inclusive_nof_events.keys() if is_valid_event_type(event_type))
          nof_events_type = set(event_type for event_type in nof_events.keys() if is_valid_event_type(event_type))
          assert(inclusive_nof_events_type == nof_events_type)
        else:
          nof_events_keys = set(nof_key for nof_key in sample_entry['nof_events'] if is_valid_event_type(nof_key))
          assert(nof_events_keys == set(nof_events.keys()))
          for nof_key, nof_arr in sample_entry['nof_events'].items():
            if not is_valid_event_type(nof_key):
              continue
            assert(len(nof_arr) == len(nof_events[nof_key]))
            for idx, nof in enumerate(nof_arr):
              assert(nof > 0)
              nof_events[nof_key][idx] += nof
        if xs < 0:
          xs = sample_entry['xsection']
    assert(xs > 0)
    binned_sample['xsection'] = xs
    binned_sample['nof_events'] = nof_events

    lumis = {}
    for nof_key, nof_arr in binned_sample['nof_events'].items():
      lumis[nof_key] = list(map(lambda nof: nof / binned_sample['xsection'], nof_arr))
    binned_sample['lumis'] = lumis

  for binned_sample in split_dict_2:
    nof_events = {}
    xs = -1
    for sample_key, sample_entry in samples.items():
      if sample_key == 'sum_events': continue
      if sample_entry['process_name_specific'] in binned_sample['samples']:
        if not nof_events:
          nof_events = copy_nof_events(sample_entry)
          inclusive_nof_events_type = set(event_type for event_type in inclusive_nof_events.keys() if is_valid_event_type(event_type))
          nof_events_type = set(event_type for event_type in nof_events.keys() if is_valid_event_type(event_type))
          assert(inclusive_nof_events_type == nof_events_type)
        else:
          nof_events_keys = set(nof_key for nof_key in sample_entry['nof_events'] if is_valid_event_type(nof_key))
          assert(nof_events_keys == set(nof_events.keys()))
          for nof_key, nof_arr in sample_entry['nof_events'].items():
            if not is_valid_event_type(nof_key):
              continue
            assert(len(nof_arr) == len(nof_events[nof_key]))
            for idx, nof in enumerate(nof_arr):
              assert(nof > 0)
              nof_events[nof_key][idx] += nof
        if xs < 0:
          xs = sample_entry['xsection']
    assert(xs > 0)
    binned_sample['xsection'] = xs
    binned_sample['nof_events'] = nof_events

    lumis = {}
    for nof_key, nof_arr in binned_sample['nof_events'].items():
      lumis[nof_key] = list(map(lambda nof: nof / binned_sample['xsection'], nof_arr))
    binned_sample['lumis'] = lumis

  # compute integrated luminosities for the inclusive sample
  inclusive_lumis = {}
  for nof_key, nof_arr in inclusive_nof_events.items():
    inclusive_lumis[nof_key] = list(map(lambda nof: nof / inclusive_xs, nof_arr))

  # decide on the bin indices
  idxs_split_sample_1 = []
  idxs_split_sample_2 = []
  for binned_sample in split_dict_1:
    binned_idx = complete_binning_1.index(binned_sample['value'][0]) + 1
    binned_sample['idx'] = binned_idx
    idxs_split_sample_1.append(binned_idx)
  for binned_sample in split_dict_2:
    binned_idx = complete_binning_2.index(binned_sample['value'][0]) + 1
    binned_sample['idx'] = binned_idx
    idxs_split_sample_2.append(binned_idx)

  for inclusive_sample in inclusive_samples:
    if inclusive_sample not in [ key.GetName() for key in f.GetListOfKeys() ]:
      histogram_dir_root = f.mkdir(inclusive_sample)
    else:
      histogram_dir_root = f.Get(inclusive_sample)
    subdir_name = '%s_v_%s' % (split_var_1, split_var_2)
    if subdir_name not in [ key.GetName() for key in histogram_dir_root.GetListOfKeys() ]:
      histogram_dir = histogram_dir_root.mkdir(subdir_name)
    else:
      histogram_dir = histogram_dir_root.Get(subdir_name)
    histogram_dir.cd()

    for nof_key, lumi_arr in inclusive_lumis.items():
      if not is_valid_event_type(nof_key):
        continue
      for idx, lumi_incl in enumerate(lumi_arr):

        histogram_name = '%s_%d' % (nof_key, idx)
        binning_1 = array.array('f', complete_binning_1)
        binning_2 = array.array('f', complete_binning_2)

        histogram = ROOT.TH2D(
          histogram_name, histogram_name, len(binning_1) - 1, binning_1, len(binning_2) - 1, binning_2
        )
        histogram.SetDirectory(histogram_dir)
        histogram.SetXTitle(split_var_1)
        histogram.SetYTitle(split_var_2)

        for split_idx_1 in range(1, len(binning_1)):
          lumi_split_1 = 0.
          if split_idx_1 in idxs_split_sample_1:
            for binned_sample in split_dict_1:
              if split_idx_1 == binned_sample['idx']:
                lumi_split_1 = binned_sample['lumis'][nof_key][idx]
                break
          for split_idx_2 in range(1, len(binning_2)):
            lumi_split_2 = 0.
            if split_idx_2 in idxs_split_sample_2:
              for binned_sample in split_dict_2:
                if split_idx_2 == binned_sample['idx']:
                  lumi_split_2 = binned_sample['lumis'][nof_key][idx]
                  break

            if binning_1[split_idx_1] > inclusive_binning_1[-1] or \
               binning_1[split_idx_1] < inclusive_binning_1[0] or \
               binning_2[split_idx_2] > inclusive_binning_2[-1] or \
               binning_2[split_idx_2] < inclusive_binning_2[0]:
              weight = 0.
            else:
              weight = lumi_incl / (lumi_incl + lumi_split_1 + lumi_split_2)
            assert(weight >= 0.)
            histogram.SetBinContent(split_idx_1, split_idx_2, weight)
            histogram.GetXaxis().SetBinLabel(
              split_idx_1, '%.0f <= %s < %.0f' % (
                complete_binning_1[split_idx_1 - 1], split_var_1, complete_binning_1[split_idx_1]
              )
            )
            histogram.GetYaxis().SetBinLabel(
              split_idx_2, '%.0f <= %s < %.0f' % (
                complete_binning_2[split_idx_2 - 1], split_var_2, complete_binning_2[split_idx_2]
              )
            )

        histogram.Write()

  for binned_sample_1 in split_dict_1:
    for sample_name in binned_sample_1['samples']:
      if sample_name not in [ key.GetName() for key in f.GetListOfKeys() ]:
        histogram_dir_root = f.mkdir(sample_name)
      else:
        histogram_dir_root = f.Get(sample_name)
      subdir_name = '%s_v_%s' % (split_var_1, split_var_2)
      if subdir_name not in [ key.GetName() for key in histogram_dir_root.GetListOfKeys() ]:
        histogram_dir = histogram_dir_root.mkdir(subdir_name)
      else:
        histogram_dir = histogram_dir_root.Get(subdir_name)
      histogram_dir.cd()

      for nof_key, lumi_arr in inclusive_lumis.items():
        if not is_valid_event_type(nof_key):
          continue
        for idx, lumi_incl in enumerate(lumi_arr):

          histogram_name = '%s_%d' % (nof_key, idx)
          binning_1 = array.array('f', complete_binning_1)
          binning_2 = array.array('f', complete_binning_2)

          histogram = ROOT.TH2D(
            histogram_name, histogram_name, len(binning_1) - 1, binning_1, len(binning_2) - 1, binning_2
          )
          histogram.SetDirectory(histogram_dir)
          histogram.SetXTitle(split_var_1)
          histogram.SetYTitle(split_var_2)

          lumi_split_1 = -1
          for split_idx_1 in range(1, len(binning_1)):
            if split_idx_1 == binned_sample_1['idx']:
              lumi_split_1 = binned_sample_1['lumis'][nof_key][idx]
              for split_idx_2 in range(1, len(binning_2)):
                lumi_split_2 = 0.
                for binned_sample_2 in split_dict_2:
                  if split_idx_2 == binned_sample_2['idx']:
                    lumi_split_2 = binned_sample_2['lumis'][nof_key][idx]
                    break
                if binning_1[split_idx_1] > inclusive_binning_1[-1] or \
                   binning_1[split_idx_1] < inclusive_binning_1[0] or \
                   binning_2[split_idx_2] > inclusive_binning_2[-1] or \
                   binning_2[split_idx_2] < inclusive_binning_2[0]:
                  if lumi_split_2 == 0.:
                    weight = 1.
                  else:
                    if apply_sf:
                      weight = lumi_split_1 / (lumi_split_1 + lumi_split_2)
                    else:
                      weight = lumi_incl / (lumi_split_1 + lumi_split_2)
                else:
                  if apply_sf:
                    weight = lumi_split_1 / (lumi_incl + lumi_split_1 + lumi_split_2)
                  else:
                    weight = lumi_incl / (lumi_incl + lumi_split_1 + lumi_split_2)
                assert(weight >= 0.)
                histogram.SetBinContent(split_idx_1, split_idx_2, weight)
            else:
              for split_idx_2 in range(1, len(binning_2)):
                histogram.SetBinContent(split_idx_1, split_idx_2, 0.)

            histogram.GetXaxis().SetBinLabel(
              split_idx_1, '%.0f <= %s < %.0f' % (
                complete_binning_1[split_idx_1 - 1], split_var_1, complete_binning_1[split_idx_1]
              )
            )
          for split_idx_2 in range(1, len(binning_2)):
            histogram.GetYaxis().SetBinLabel(
              split_idx_2, '%.0f <= %s < %.0f' % (
                complete_binning_2[split_idx_2 - 1], split_var_2, complete_binning_2[split_idx_2]
              )
            )

          histogram.Write()

  for binned_sample_2 in split_dict_2:
    for sample_name in binned_sample_2['samples']:
      if sample_name not in [ key.GetName() for key in f.GetListOfKeys() ]:
        histogram_dir_root = f.mkdir(sample_name)
      else:
        histogram_dir_root = f.Get(sample_name)
      subdir_name = '%s_v_%s' % (split_var_1, split_var_2)
      if subdir_name not in [ key.GetName() for key in histogram_dir_root.GetListOfKeys() ]:
        histogram_dir = histogram_dir_root.mkdir(subdir_name)
      else:
        histogram_dir = histogram_dir_root.Get(subdir_name)
      histogram_dir.cd()

      for nof_key, lumi_arr in inclusive_lumis.items():
        if not is_valid_event_type(nof_key):
          continue
        for idx, lumi_incl in enumerate(lumi_arr):

          histogram_name = '%s_%d' % (nof_key, idx)
          binning_1 = array.array('f', complete_binning_1)
          binning_2 = array.array('f', complete_binning_2)

          histogram = ROOT.TH2D(
            histogram_name, histogram_name, len(binning_1) - 1, binning_1, len(binning_2) - 1, binning_2
          )
          histogram.SetDirectory(histogram_dir)
          histogram.SetXTitle(split_var_1)
          histogram.SetYTitle(split_var_2)

          lumi_split_2 = -1
          for split_idx_2 in range(1, len(binning_2)):
            if split_idx_2 == binned_sample_2['idx']:
              lumi_split_2 = binned_sample_2['lumis'][nof_key][idx]
              for split_idx_1 in range(1, len(binning_1)):
                lumi_split_1 = 0.
                for binned_sample_1 in split_dict_1:
                  if split_idx_1 == binned_sample_1['idx']:
                    lumi_split_1 = binned_sample_1['lumis'][nof_key][idx]
                    break
                if binning_1[split_idx_1] > inclusive_binning_1[-1] or \
                   binning_1[split_idx_1] < inclusive_binning_1[0] or \
                   binning_2[split_idx_2] > inclusive_binning_2[-1] or \
                   binning_2[split_idx_2] < inclusive_binning_2[0]:
                  if lumi_split_1 == 0.:
                    weight = 1.
                  else:
                    if apply_sf:
                      weight = lumi_split_2 / (lumi_split_1 + lumi_split_2)
                    else:
                      weight = lumi_incl / (lumi_split_1 + lumi_split_2)
                else:
                  if apply_sf:
                    weight = lumi_split_2 / (lumi_incl + lumi_split_1 + lumi_split_2)
                  else:
                    weight = lumi_incl / (lumi_incl + lumi_split_1 + lumi_split_2)
                histogram.SetBinContent(split_idx_1, split_idx_2, weight)
            else:
              for split_idx_1 in range(1, len(binning_1)):
                histogram.SetBinContent(split_idx_1, split_idx_2, 0.)

            histogram.GetYaxis().SetBinLabel(
              split_idx_2, '%.0f <= %s < %.0f' % (
                complete_binning_2[split_idx_2 - 1], split_var_2, complete_binning_2[split_idx_2]
              )
            )
          for split_idx_1 in range(1, len(binning_1)):
            histogram.GetXaxis().SetBinLabel(
              split_idx_1, '%.0f <= %s < %.0f' % (
                complete_binning_1[split_idx_1 - 1], split_var_1, complete_binning_1[split_idx_1]
              )
            )

          histogram.Write()

def comp_weights_2_wo_inclusive(f, samples, samples_to_stitch, split_var_1, split_var_2):
  inclusive_samples   = samples_to_stitch['inclusive']['samples']
  inclusive_binning_1 = samples_to_stitch['inclusive'][split_var_1]
  inclusive_binning_2 = samples_to_stitch['inclusive'][split_var_2]

  split_dict_1 = samples_to_stitch['exclusive'][split_var_1]
  split_dict_2 = samples_to_stitch['exclusive'][split_var_2]
  split_binning_1 = [ sample['value'] for sample in split_dict_1 ]
  split_binning_2 = [ sample['value'] for sample in split_dict_2 ]
  complete_binning_1 = list(sorted(list(
    map(float, set(inclusive_binning_1) | set(list(itertools.chain.from_iterable(split_binning_1))))
  )))
  complete_binning_2 = list(sorted(list(
    map(float, set(inclusive_binning_2) | set(list(itertools.chain.from_iterable(split_binning_2))))
  )))

  inclusive_xs = -1
  for sample_key, sample_entry in samples.items():
    if sample_key == 'sum_events': continue
    if sample_entry['process_name_specific'] == inclusive_samples[0]:
      inclusive_xs = sample_entry['xsection']
  assert(inclusive_xs > 0)

  # sum the inclusive nof events
  inclusive_nof_events = {}
  for sample_key, sample_entry in samples.items():
    if sample_key == 'sum_events': continue
    if sample_entry['process_name_specific'] in inclusive_samples:
      if not inclusive_nof_events:
        inclusive_nof_events = copy_nof_events(sample_entry)
      else:
        nof_events_keys = set(nof_key for nof_key in sample_entry['nof_events'] if is_valid_event_type(nof_key))
        assert(nof_events_keys == set(inclusive_nof_events))
        for nof_key, nof_arr in sample_entry['nof_events'].items():
          if not is_valid_event_type(nof_key):
            continue
          assert(len(nof_arr) == len(inclusive_nof_events[nof_key]))
          for idx, nof in enumerate(nof_arr):
            assert(nof > 0)
            inclusive_nof_events[nof_key][idx] += nof

  # sum the binned nof events
  for binned_sample in split_dict_1:
    nof_events = {}
    xs = -1
    for sample_key, sample_entry in samples.items():
      if sample_key == 'sum_events': continue
      if sample_entry['process_name_specific'] in binned_sample['samples']:
        if not nof_events:
          nof_events = copy_nof_events(sample_entry)
          inclusive_nof_events_type = set(event_type for event_type in inclusive_nof_events.keys() if is_valid_event_type(event_type))
          nof_events_type = set(event_type for event_type in nof_events.keys() if is_valid_event_type(event_type))
          assert(inclusive_nof_events_type == nof_events_type)
        else:
          nof_events_keys = set(nof_key for nof_key in sample_entry['nof_events'] if is_valid_event_type(nof_key))
          assert(nof_events_keys == set(nof_events.keys()))
          for nof_key, nof_arr in sample_entry['nof_events'].items():
            if not is_valid_event_type(nof_key):
              continue
            assert(len(nof_arr) == len(nof_events[nof_key]))
            for idx, nof in enumerate(nof_arr):
              assert(nof > 0)
              nof_events[nof_key][idx] += nof
        if xs < 0:
          xs = sample_entry['xsection']
    assert(xs > 0)
    binned_sample['xsection'] = xs
    binned_sample['nof_events'] = nof_events

    lumis = {}
    for nof_key, nof_arr in binned_sample['nof_events'].items():
      lumis[nof_key] = list(map(lambda nof: nof / binned_sample['xsection'], nof_arr))
    binned_sample['lumis'] = lumis

  for binned_sample in split_dict_2:
    nof_events = {}
    xs = -1
    for sample_key, sample_entry in samples.items():
      if sample_key == 'sum_events': continue
      if sample_entry['process_name_specific'] in binned_sample['samples']:
        if not nof_events:
          nof_events = copy_nof_events(sample_entry)
          inclusive_nof_events_type = set(event_type for event_type in inclusive_nof_events.keys() if is_valid_event_type(event_type))
          nof_events_type = set(event_type for event_type in nof_events.keys() if is_valid_event_type(event_type))
          assert(inclusive_nof_events_type == nof_events_type)
        else:
          nof_events_keys = set(nof_key for nof_key in sample_entry['nof_events'] if is_valid_event_type(nof_key))
          assert(nof_events_keys == set(nof_events.keys()))
          for nof_key, nof_arr in sample_entry['nof_events'].items():
            if not is_valid_event_type(nof_key):
              continue
            assert(len(nof_arr) == len(nof_events[nof_key]))
            for idx, nof in enumerate(nof_arr):
              assert(nof > 0)
              nof_events[nof_key][idx] += nof
        if xs < 0:
          xs = sample_entry['xsection']
    assert(xs > 0)
    binned_sample['xsection'] = xs
    binned_sample['nof_events'] = nof_events

    lumis = {}
    for nof_key, nof_arr in binned_sample['nof_events'].items():
      lumis[nof_key] = list(map(lambda nof: nof / binned_sample['xsection'], nof_arr))
    binned_sample['lumis'] = lumis

  # compute integrated luminosities for the inclusive sample
  inclusive_lumis = {}
  for nof_key, nof_arr in inclusive_nof_events.items():
    inclusive_lumis[nof_key] = list(map(lambda nof: nof / inclusive_xs, nof_arr))

  # decide on the bin indices
  idxs_split_sample_1 = []
  idxs_split_sample_2 = []
  for binned_sample in split_dict_1:
    binned_idx = complete_binning_1.index(binned_sample['value'][0]) + 1
    binned_sample['idx'] = binned_idx
    idxs_split_sample_1.append(binned_idx)
  for binned_sample in split_dict_2:
    binned_idx = complete_binning_2.index(binned_sample['value'][0]) + 1
    binned_sample['idx'] = binned_idx
    idxs_split_sample_2.append(binned_idx)

  for binned_sample_1 in split_dict_1:
    for sample_name in binned_sample_1['samples']:
      if sample_name not in [ key.GetName() for key in f.GetListOfKeys() ]:
        histogram_dir_root = f.mkdir(sample_name)
      else:
        histogram_dir_root = f.Get(sample_name)
      subdir_name = '%s_v_%s_wo_inclusive' % (split_var_1, split_var_2)
      if subdir_name not in [ key.GetName() for key in histogram_dir_root.GetListOfKeys() ]:
        histogram_dir = histogram_dir_root.mkdir(subdir_name)
      else:
        histogram_dir = histogram_dir_root.Get(subdir_name)
      histogram_dir.cd()

      for nof_key, lumi_arr in inclusive_lumis.items():
        if not is_valid_event_type(nof_key):
          continue
        for idx, lumi_incl in enumerate(lumi_arr):

          histogram_name = '%s_%d' % (nof_key, idx)
          binning_1 = array.array('f', complete_binning_1)
          binning_2 = array.array('f', complete_binning_2)

          histogram = ROOT.TH2D(
            histogram_name, histogram_name, len(binning_1) - 1, binning_1, len(binning_2) - 1, binning_2
          )
          histogram.SetDirectory(histogram_dir)
          histogram.SetXTitle(split_var_1)
          histogram.SetYTitle(split_var_2)

          lumi_split_1 = -1
          for split_idx_1 in range(1, len(binning_1)):
            if split_idx_1 == binned_sample_1['idx']:
              lumi_split_1 = binned_sample_1['lumis'][nof_key][idx]
              for split_idx_2 in range(1, len(binning_2)):
                lumi_split_2 = 0.
                for binned_sample_2 in split_dict_2:
                  if split_idx_2 == binned_sample_2['idx']:
                    lumi_split_2 = binned_sample_2['lumis'][nof_key][idx]
                    break
                if binning_1[split_idx_1] > inclusive_binning_1[-1] or \
                   binning_1[split_idx_1] < inclusive_binning_1[0] or \
                   binning_2[split_idx_2] > inclusive_binning_2[-1] or \
                   binning_2[split_idx_2] < inclusive_binning_2[0]:
                  if lumi_split_2 == 0.:
                    weight = 1.
                  else:
                    weight = lumi_split_1 / (lumi_split_1 + lumi_split_2)
                else:
                  weight = lumi_split_1 / (lumi_split_1 + lumi_split_2)
                assert(weight >= 0.)
                histogram.SetBinContent(split_idx_1, split_idx_2, weight)
            else:
              for split_idx_2 in range(1, len(binning_2)):
                histogram.SetBinContent(split_idx_1, split_idx_2, 0.)

            histogram.GetXaxis().SetBinLabel(
              split_idx_1, '%.0f <= %s < %.0f' % (
                complete_binning_1[split_idx_1 - 1], split_var_1, complete_binning_1[split_idx_1]
              )
            )
          for split_idx_2 in range(1, len(binning_2)):
            histogram.GetYaxis().SetBinLabel(
              split_idx_2, '%.0f <= %s < %.0f' % (
                complete_binning_2[split_idx_2 - 1], split_var_2, complete_binning_2[split_idx_2]
              )
            )

          histogram.Write()

  for binned_sample_2 in split_dict_2:
    for sample_name in binned_sample_2['samples']:
      if sample_name not in [ key.GetName() for key in f.GetListOfKeys() ]:
        histogram_dir_root = f.mkdir(sample_name)
      else:
        histogram_dir_root = f.Get(sample_name)
      subdir_name = '%s_v_%s_wo_inclusive' % (split_var_1, split_var_2)
      if subdir_name not in [ key.GetName() for key in histogram_dir_root.GetListOfKeys() ]:
        histogram_dir = histogram_dir_root.mkdir(subdir_name)
      else:
        histogram_dir = histogram_dir_root.Get(subdir_name)
      histogram_dir.cd()

      for nof_key, lumi_arr in inclusive_lumis.items():
        if not is_valid_event_type(nof_key):
          continue
        for idx, lumi_incl in enumerate(lumi_arr):

          histogram_name = '%s_%d' % (nof_key, idx)
          binning_1 = array.array('f', complete_binning_1)
          binning_2 = array.array('f', complete_binning_2)

          histogram = ROOT.TH2D(
            histogram_name, histogram_name, len(binning_1) - 1, binning_1, len(binning_2) - 1, binning_2
          )
          histogram.SetDirectory(histogram_dir)
          histogram.SetXTitle(split_var_1)
          histogram.SetYTitle(split_var_2)

          lumi_split_2 = -1
          for split_idx_2 in range(1, len(binning_2)):
            if split_idx_2 == binned_sample_2['idx']:
              lumi_split_2 = binned_sample_2['lumis'][nof_key][idx]
              for split_idx_1 in range(1, len(binning_1)):
                lumi_split_1 = 0.
                for binned_sample_1 in split_dict_1:
                  if split_idx_1 == binned_sample_1['idx']:
                    lumi_split_1 = binned_sample_1['lumis'][nof_key][idx]
                    break
                if binning_1[split_idx_1] > inclusive_binning_1[-1] or \
                   binning_1[split_idx_1] < inclusive_binning_1[0] or \
                   binning_2[split_idx_2] > inclusive_binning_2[-1] or \
                   binning_2[split_idx_2] < inclusive_binning_2[0]:
                  if lumi_split_1 == 0.:
                    weight = 1.
                  else:
                    weight = lumi_split_2 / (lumi_split_1 + lumi_split_2)
                else:
                  weight = lumi_split_2 / (lumi_split_1 + lumi_split_2)
                histogram.SetBinContent(split_idx_1, split_idx_2, weight)
            else:
              for split_idx_1 in range(1, len(binning_1)):
                histogram.SetBinContent(split_idx_1, split_idx_2, 0.)

            histogram.GetYaxis().SetBinLabel(
              split_idx_2, '%.0f <= %s < %.0f' % (
                complete_binning_2[split_idx_2 - 1], split_var_2, complete_binning_2[split_idx_2]
              )
            )
          for split_idx_1 in range(1, len(binning_1)):
            histogram.GetXaxis().SetBinLabel(
              split_idx_1, '%.0f <= %s < %.0f' % (
                complete_binning_1[split_idx_1 - 1], split_var_1, complete_binning_1[split_idx_1]
              )
            )

          histogram.Write()

parser = argparse.ArgumentParser(
    formatter_class = lambda prog: SmartFormatter(prog, max_help_position = 35)
)
parser.add_argument('-e', '--era',
  type = str, dest = 'era', metavar = 'era', required = True, choices = [ '2016', '2017', '2018' ],
  help = 'R|Era',
)
parser.add_argument('-o', '--output',
  type = str, dest = 'output', metavar = 'path', required = False,
  help = 'R|Output file name',
)
args = parser.parse_args()

era = args.era
output = args.output

if not output:
  output = os.path.join(
    os.environ['CMSSW_BASE'], 'src', 'tthAnalysis', 'HiggsToTauTau', 'data', 'stitched_weights_lo_{}.root'.format(era)
  )

if era == '2016':
  from tthAnalysis.HiggsToTauTau.samples.stitch import samples_to_stitch_lo_2016 as samples_to_stitch
elif era == '2017':
  from tthAnalysis.HiggsToTauTau.samples.stitch import samples_to_stitch_lo_2017 as samples_to_stitch
elif era == '2018':
  from tthAnalysis.HiggsToTauTau.samples.stitch import samples_to_stitch_lo_2018 as samples_to_stitch
else:
  raise RuntimeError("Invalid era: %s" % era)

samples = load_samples(era)

apply_sf = True
fp = ROOT.TFile.Open(output, 'recreate')
for sample_to_stitch in samples_to_stitch:
  binning_vars = sorted([ var for var in sample_to_stitch['exclusive'] ], reverse = True)
  if len(binning_vars) == 1:
    comp_weights_1(fp, samples, sample_to_stitch, binning_vars[0], apply_sf)
  elif len(binning_vars) == 2:
    for binning_var in binning_vars:
      comp_weights_1(fp, samples, sample_to_stitch, binning_var, apply_sf)
    comp_weights_2(fp, samples, sample_to_stitch, binning_vars[0], binning_vars[1], apply_sf)
    comp_weights_2_wo_inclusive(fp, samples, sample_to_stitch, binning_vars[0], binning_vars[1])
  else:
    raise ValueError('More than 2 variables used for the binning: %s' % ', '.join(binning_vars))
fp.Close()

