from tthAnalysis.HiggsToTauTau.samples.tthAnalyzeSamples_2017 import samples_2017

bdt_samples = [
  "ttHToNonbb_M125_powheg",
  "ttHToNonbb_M125_powheg_ext1",
  "TTZJets_LO",
  "TTZJets_LO_ext1",
  "TTWJets_LO",
  "TTWJets_LO_ext1",
  "TTJets_DiLept",
  "TTJets_SingleLeptFromT",
  "TTJets_SingleLeptFromTbar",
]

for sample_name, sample_info in samples_2017.items():
  if sample_name == 'sum_events': continue
  sample_info["use_it"] = sample_info["process_name_specific"] in bdt_samples
