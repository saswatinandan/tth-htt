from tthAnalysis.HiggsToTauTau.samples.tthAnalyzeSamples_2016_base import samples_2016 as samples_2016_base
from hhAnalysis.multilepton.samples.hhAnalyzeSamples_2016_hh import samples_2016 as samples_2016_hh_multilepton
from hhAnalysis.bbww.samples.hhAnalyzeSamples_2016_hh import samples_2016 as samples_2016_hh_bbww

from tthAnalysis.HiggsToTauTau.samples.reclassifySamples import reclassifySamples
samples_2016     = reclassifySamples(samples_2016_base, samples_2016_hh_multilepton, samples_2016_hh_bbww)
samples_2016_aux = reclassifySamples(samples_2016_base, samples_2016_hh_multilepton, samples_2016_hh_bbww, analysis_type = 'aux')
