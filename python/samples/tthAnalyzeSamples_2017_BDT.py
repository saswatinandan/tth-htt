from tthAnalysis.HiggsToTauTau.samples.tthAnalyzeSamples_2017 import samples_2017

bdt_samples = [
  "ttHToNonbb_M125_powheg",
  "ttHToNonbb_M125_powheg_ext1",
  "TTZJets_LO",
  "TTZJets_LO_ext1",
  "TTWJets_LO",
  "TTWJets_LO_ext1",
  "TTTo2L2Nu",
  "TTTo2L2Nu_PSweights",
  "TTToSemiLeptonic",
  "TTToSemiLeptonic_PSweights",
  "TTToHadronic",
  "TTToHadronic_PSweights",
  'signal_ggf_nonresonant_node_sm_hh_2b2t',
  'signal_ggf_nonresonant_node_sm_hh_2b2v',
  'signal_ggf_nonresonant_node_sm_hh_2v2t',
  'signal_ggf_nonresonant_node_sm_hh_4t',
  'signal_ggf_nonresonant_node_sm_hh_4v' ,
  "GluGluHToTauTau",
  "GluGluHToTauTau_ext1",
  "GluGluHToZZTo4L",
  "GluGluHToZZTo2L2Q_M125",
  "GluGluHToWWToLNuQQ_M125_PSweights",
  "GluGluHToWWTo2L2Nu_M125",
  "VBFHToTauTau",
  "VBF_HToZZTo4L",
  "VBFHToWWToLNuQQ_M125_PSweights",
  "VBFHToWWTo2L2Nu_M125",
  "VHToNonbb_M125",
  "ZH_HToBB_ZToLL",
  "ZHToTauTau"
  "TTWH",
  "TTZH"
]

for sample_name, sample_info in samples_2017.items():
  if sample_name == 'sum_events': continue
  sample_info["use_it"] = sample_info["process_name_specific"] in bdt_samples
