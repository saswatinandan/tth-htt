from tthAnalysis.HiggsToTauTau.samples.tthAnalyzeSamples_2018_base import samples_2018

bdt_samples = [
  "ttHToNonbb_M125_powheg",
  "THQ_ctcvcp",
  "THW_ctcvcp",
  "TTZJets_LO_ext1",
  "TTWJets_LO_ext1",
  "TTTo2L2Nu",
  "TTToSemiLeptonic",
  "TTToHadronic",
  "GluGluHToTauTau",
  "GluGluHToZZTo4L",
  "GluGluHToZZTo2L2Q_M125",
  "GluGluHToWWToLNuQQ_M125",
  "GluGluHToWWTo2L2Nu_M125",
  "VBFHToTauTau_ext1",
  "VBF_HToZZTo4L",
  "VBFHToWWToLNuQQ_M125",
  "VBFHToWWTo2L2Nu_M125",
  "VHToNonbb_M125",
  "ZH_HToBB_ZToLL",
  "ZHToTauTau"
  "TTWH_ext1",
  "TTZH_ext1"
]

for sample_name, sample_info in samples_2018.items():
  if sample_name == 'sum_events': continue
  sample_info["use_it"] = sample_info["process_name_specific"] in bdt_samples
