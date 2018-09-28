import sys , time
from collections import OrderedDict
import ROOT

tree_base = "0l_2tau_OS_Tight/sel/evt/"
treetoread=[
    "EWK2t0e0m0j",
    "EWK1t1e0m0j",
    "EWK1t0e1m0j",
    "EWK1t0e0m1j",
    "EWK0t2e0m0j",
    "EWK0t1e1m0j",
    "EWK0t1e0m1j",
    "EWK0t0e2m0j",
    "EWK0t0e1m1j",
    "EWK0t0e0m2j"
    ]

mom = "/hdfs/local/acaan/ttHAnalysis/2017/0l_2tau_datacards_categories_2018August25/histograms/0l_2tau/"
sources = OrderedDict()

sources["DY"] = [
    [
    #"histograms_harvested_stage1_0l_2tau_DYJetsToLL_M-10to50_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_DYJetsToLL_M-4to50_HT-100to200_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_DYJetsToLL_M-4to50_HT-100to200_ext1_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_DYJetsToLL_M-4to50_HT-200to400_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_DYJetsToLL_M-4to50_HT-200to400_ext1_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_DYJetsToLL_M-4to50_HT-400to600_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_DYJetsToLL_M-4to50_HT-400to600_ext1_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_DYJetsToLL_M-4to50_HT-600toInf_Tight_OS.root",
    ]
    ]

sources["DY_M50"] = [
    [
    "histograms_harvested_stage1_0l_2tau_DY1JetsToLL_M-50_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_DY1JetsToLL_M-50_ext1_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_DY2JetsToLL_M-50_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_DY2JetsToLL_M-50_ext1_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_DY3JetsToLL_M-50_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_DY3JetsToLL_M-50_ext1_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_DY4JetsToLL_M-50_Tight_OS.root",
    #"histograms_harvested_stage1_0l_2tau_DYJetsToLL_M-10to50_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_DYJetsToLL_M-50_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_DYJetsToLL_M-50_ext1_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_DYJetsToLL_M50_HT100to200_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_DYJetsToLL_M50_HT100to200_ext1_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_DYJetsToLL_M50_HT1200to2500_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_DYJetsToLL_M50_HT200to400_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_DYJetsToLL_M50_HT200to400_ext1_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_DYJetsToLL_M50_HT2500toInf_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_DYJetsToLL_M50_HT400to600_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_DYJetsToLL_M50_HT400to600_ext1_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_DYJetsToLL_M50_HT600to800_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_DYJetsToLL_M50_HT800to1200_Tight_OS.root",
    ]
    ]

sources["WJets"] = [
    [
    "histograms_harvested_stage1_0l_2tau_W1JetsToLNu_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_W2JetsToLNu_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_W3JetsToLNu_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_W4JetsToLNu_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_WJetsToLNu_HT100To200_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_WJetsToLNu_HT1200To2500_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_WJetsToLNu_HT200To400_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_WJetsToLNu_HT2500ToInf_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_WJetsToLNu_HT400To600_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_WJetsToLNu_HT600To800_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_WJetsToLNu_HT800To1200_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_WJetsToLNu_Tight_OS.root",
    ]
]

sources["WW"] = [
    [
    "histograms_harvested_stage1_0l_2tau_WWTo2L2Nu_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_WWToLNuQQ_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_WWToLNuQQ_ext1_Tight_OS.root",
    ]
]

sources["WZ"] = [
    [
    "histograms_harvested_stage1_0l_2tau_WZTo3LNu_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_WZTo1L1Nu2Q_Tight_OS.root",
    ]
]

sources["ZZ"] = [
    [
    "histograms_harvested_stage1_0l_2tau_ZZTo2L2Nu_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_ZZTo4L_Tight_OS.root",
    "histograms_harvested_stage1_0l_2tau_ZZTo4L_ext1_Tight_OS.root",
    ]
]

countall = 0
for key in  sources.keys() :
    counter = 0
    nev = 0
    for file in sources[key][0] :
      tfile = ROOT.TFile(mom + file)
      for genbase in treetoread :
        #print tree_base+genbase+"/EventCounter"
        histo = tfile.Get(tree_base+genbase+"/EventCounter")
        counter += histo.Integral()
        nev += histo.GetEntries()
    print (key, len(sources[key][0]), counter, nev)
    countall += counter
print ("total", countall)
