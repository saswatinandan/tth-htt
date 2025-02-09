from collections import OrderedDict as OD

# file generated at 2020-10-26 12:42:25 with the following command:
# create_dictionary.py -m python/samples/metaDict_2018_sync.py -p /local/karl/ttHNtupleProduction/2018/2020Oct26_wPresel_nonNom_tth_sync/ntuples -N samples_2018 -E 2018 -o python/samples -g tthAnalyzeSamples_2018_preselected_sync.py -M

samples_2018 = OD()
samples_2018["/ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM"] = OD([
  ("type",                            "mc"),
  ("sample_category",                 "ttH"),
  ("process_name_specific",           "ttHToNonbb_M125_powheg"),
  ("nof_files",                       1),
  ("nof_db_files",                    224),
  ("nof_events",                      {
    'Count'                                                                          : [ 63000, ],
    'CountWeighted'                                                                  : [ 6.16550078e+04, 6.16721758e+04, 6.16409766e+04, ],
    'CountWeightedLHEWeightScale'                                                    : [ 8.97077031e+04, 8.70789922e+04, 8.50663281e+04, 6.31800820e+04, 6.16550078e+04, 6.05346953e+04, 5.17520078e+04, 5.03899258e+04, 4.93522227e+04, ],
    'CountWeightedLHEEnvelope'                                                       : [ 9.08448516e+04, 4.86382656e+04, ],
    'CountWeightedPSWeight'                                                          : [ 6.16897617e+04, 6.18149453e+04, 9.19098047e+04, 6.16137773e+04, 6.11918203e+04, 3.42826406e+04, ],
    'CountWeightedPSWeightOriginalXWGTUP'                                            : [ 3.39102539e+04, 3.40079375e+04, 5.05702305e+04, 3.38740039e+04, 3.36386172e+04, 1.88316230e+04, ],
    'CountWeightedFull'                                                              : [ 3.31697422e+04, 3.31749844e+04, 3.31615977e+04, ],
    'CountWeightedFullLHEWeightScale'                                                : [ 4.82600938e+04, 4.68459570e+04, 4.57633672e+04, 3.39892148e+04, 3.31697422e+04, 3.25662012e+04, 2.78409805e+04, 2.71084746e+04, 2.65502188e+04, ],
    'CountWeightedFullLHEEnvelope'                                                   : [ 4.88720117e+04, 2.61659824e+04, ],
    'CountWeightedFullPSWeight'                                                      : [ 3.31873281e+04, 3.32545625e+04, 4.94446953e+04, 3.31465586e+04, 3.29198047e+04, 1.84431895e+04, ],
    'CountWeightedFullPSWeightOriginalXWGTUP'                                        : [ 1.82427949e+04, 1.82951855e+04, 2.72055410e+04, 1.82231113e+04, 1.80966641e+04, 1.01310186e+04, ],
    'Count_fwd'                                                                      : [ 952, ],
    'CountWeighted_fwd'                                                              : [ 9.38807007e+02, 9.34249634e+02, 9.43504089e+02, ],
    'CountWeightedLHEWeightScale_fwd'                                                : [ 1.33891028e+03, 1.28863757e+03, 1.25051086e+03, 9.71085815e+02, 9.38807007e+02, 9.14465454e+02, 7.97404358e+02, 7.69096252e+02, 7.47273987e+02, ],
    'CountWeightedLHEEnvelope_fwd'                                                   : [ 1.34829053e+03, 7.41973022e+02, ],
    'CountWeightedPSWeight_fwd'                                                      : [ 9.42432617e+02, 9.30160217e+02, 1.34558240e+03, 9.33966187e+02, 9.29800537e+02, 5.54618530e+02, ],
    'CountWeightedPSWeightOriginalXWGTUP_fwd'                                        : [ 5.17883179e+02, 5.13209717e+02, 7.40759094e+02, 5.14576416e+02, 5.09347137e+02, 3.04778625e+02, ],
    'CountWeightedFull_fwd'                                                          : [ 5.05048950e+02, 5.02603851e+02, 5.07569397e+02, ],
    'CountWeightedFullLHEWeightScale_fwd'                                            : [ 7.20296143e+02, 6.93249573e+02, 6.72736877e+02, 5.22414795e+02, 5.05048950e+02, 4.91953583e+02, 4.28979523e+02, 4.13749786e+02, 4.02010162e+02, ],
    'CountWeightedFullLHEEnvelope_fwd'                                               : [ 7.25342712e+02, 3.99158051e+02, ],
    'CountWeightedFullPSWeight_fwd'                                                  : [ 5.06997101e+02, 5.00390594e+02, 7.23908936e+02, 5.02447327e+02, 5.00244019e+02, 2.98364136e+02, ],
    'CountWeightedFullPSWeightOriginalXWGTUP_fwd'                                    : [ 2.78603912e+02, 2.76087585e+02, 3.98520752e+02, 2.76828308e+02, 2.74035492e+02, 1.63960052e+02, ],
    'Count_pt0to60'                                                                  : [ 14788, ],
    'CountWeighted_pt0to60'                                                          : [ 1.45402266e+04, 1.45525840e+04, 1.45239189e+04, ],
    'CountWeightedLHEWeightScale_pt0to60'                                            : [ 2.12015781e+04, 2.06424180e+04, 2.02291719e+04, 1.48630576e+04, 1.45402266e+04, 1.43144346e+04, 1.21899688e+04, 1.18974961e+04, 1.16824883e+04, ],
    'CountWeightedLHEEnvelope_pt0to60'                                               : [ 2.15096406e+04, 1.15082275e+04, ],
    'CountWeightedPSWeight_pt0to60'                                                  : [ 1.45441143e+04, 1.46699648e+04, 2.15450840e+04, 1.45355049e+04, 1.42444541e+04, 8.13005713e+03, ],
    'CountWeightedPSWeightOriginalXWGTUP_pt0to60'                                    : [ 7.94980811e+03, 8.02249609e+03, 1.17859443e+04, 7.94663916e+03, 7.78876611e+03, 4.44028223e+03, ],
    'CountWeightedFull_pt0to60'                                                      : [ 7.82243457e+03, 7.82913965e+03, 7.81367188e+03, ],
    'CountWeightedFullLHEWeightScale_pt0to60'                                        : [ 1.14059873e+04, 1.11052080e+04, 1.08828174e+04, 7.99598682e+03, 7.82243457e+03, 7.70086719e+03, 6.55797119e+03, 6.40058447e+03, 6.28491162e+03, ],
    'CountWeightedFullLHEEnvelope_pt0to60'                                           : [ 1.15717422e+04, 6.19115820e+03, ],
    'CountWeightedFullPSWeight_pt0to60'                                              : [ 7.82442969e+03, 7.89207178e+03, 1.15906270e+04, 7.81978174e+03, 7.66314258e+03, 4.37383398e+03, ],
    'CountWeightedFullPSWeightOriginalXWGTUP_pt0to60'                                : [ 4.27684180e+03, 4.31590869e+03, 6.34051221e+03, 4.27512109e+03, 4.19018506e+03, 2.38879419e+03, ],
    'Count_pt60to120'                                                                : [ 22285, ],
    'CountWeighted_pt60to120'                                                        : [ 2.19557031e+04, 2.19454199e+04, 2.19675918e+04, ],
    'CountWeightedLHEWeightScale_pt60to120'                                          : [ 3.19362598e+04, 3.10535410e+04, 3.03901895e+04, 2.24674316e+04, 2.19557031e+04, 2.15898145e+04, 1.84113301e+04, 1.79531172e+04, 1.76104336e+04, ],
    'CountWeightedLHEEnvelope_pt60to120'                                             : [ 3.23696582e+04, 1.73541406e+04, ],
    'CountWeightedPSWeight_pt60to120'                                                : [ 2.19648438e+04, 2.20273496e+04, 3.27540352e+04, 2.19429512e+04, 2.17975273e+04, 1.21901074e+04, ],
    'CountWeightedPSWeightOriginalXWGTUP_pt60to120'                                  : [ 1.20137344e+04, 1.20450537e+04, 1.79234336e+04, 1.20041963e+04, 1.19331201e+04, 6.66792285e+03, ],
    'CountWeightedFull_pt60to120'                                                    : [ 1.18123438e+04, 1.18066943e+04, 1.18176875e+04, ],
    'CountWeightedFullLHEWeightScale_pt60to120'                                      : [ 1.71811055e+04, 1.67062402e+04, 1.63494131e+04, 1.20870430e+04, 1.18123438e+04, 1.16148984e+04, 9.90489746e+03, 9.65837598e+03, 9.47406250e+03, ],
    'CountWeightedFullLHEEnvelope_pt60to120'                                         : [ 1.74142266e+04, 9.33615137e+03, ],
    'CountWeightedFullPSWeight_pt60to120'                                            : [ 1.18167510e+04, 1.18502617e+04, 1.76210508e+04, 1.18048848e+04, 1.17267363e+04, 6.55808350e+03, ],
    'CountWeightedFullPSWeightOriginalXWGTUP_pt60to120'                              : [ 6.46316504e+03, 6.47992871e+03, 9.64241113e+03, 6.45804688e+03, 6.41982617e+03, 3.58722339e+03, ],
    'Count_pt120to200'                                                               : [ 15608, ],
    'CountWeighted_pt120to200'                                                       : [ 1.52626592e+04, 1.52773086e+04, 1.52446377e+04, ],
    'CountWeightedLHEWeightScale_pt120to200'                                         : [ 2.22015605e+04, 2.15375273e+04, 2.10209883e+04, 1.56427012e+04, 1.52626592e+04, 1.49790088e+04, 1.28136602e+04, 1.24734053e+04, 1.22119912e+04, ],
    'CountWeightedLHEEnvelope_pt120to200'                                            : [ 2.24620781e+04, 1.20355498e+04, ],
    'CountWeightedPSWeight_pt120to200'                                               : [ 1.52730996e+04, 1.52254072e+04, 2.29157109e+04, 1.52537783e+04, 1.53563857e+04, 8.46324121e+03, ],
    'CountWeightedPSWeightOriginalXWGTUP_pt120to200'                                 : [ 8.38736914e+03, 8.38187402e+03, 1.26175537e+04, 8.37920703e+03, 8.43850488e+03, 4.64231836e+03, ],
    'CountWeightedFull_pt120to200'                                                   : [ 8.21109863e+03, 8.21900781e+03, 8.20144434e+03, ],
    'CountWeightedFullLHEWeightScale_pt120to200'                                     : [ 1.19438896e+04, 1.15866982e+04, 1.13088682e+04, 8.41533887e+03, 8.21109863e+03, 8.05836719e+03, 6.89340576e+03, 6.71039111e+03, 6.56977637e+03, ],
    'CountWeightedFullLHEEnvelope_pt120to200'                                        : [ 1.20840732e+04, 6.47484570e+03, ],
    'CountWeightedFullPSWeight_pt120to200'                                           : [ 8.21659668e+03, 8.19097412e+03, 1.23281289e+04, 8.20616797e+03, 8.26138184e+03, 4.55300000e+03, ],
    'CountWeightedFullPSWeightOriginalXWGTUP_pt120to200'                             : [ 4.51220898e+03, 4.50926514e+03, 6.78798828e+03, 4.50780371e+03, 4.53971387e+03, 2.49744995e+03, ],
    'Count_pt200to300'                                                               : [ 6386, ],
    'CountWeighted_pt200to300'                                                       : [ 6.23467529e+03, 6.21736230e+03, 6.25532080e+03, ],
    'CountWeightedLHEWeightScale_pt200to300'                                         : [ 9.06121094e+03, 8.73484277e+03, 8.47546094e+03, 6.42978418e+03, 6.23467529e+03, 6.08130615e+03, 5.25664746e+03, 5.08577979e+03, 4.94983008e+03, ],
    'CountWeightedLHEEnvelope_pt200to300'                                            : [ 9.15140039e+03, 4.88152588e+03, ],
    'CountWeightedPSWeight_pt200to300'                                               : [ 6.23457568e+03, 6.24865527e+03, 9.30816797e+03, 6.22964014e+03, 6.18098535e+03, 3.44104688e+03, ],
    'CountWeightedPSWeightOriginalXWGTUP_pt200to300'                                 : [ 3.45336255e+03, 3.45761621e+03, 5.14386426e+03, 3.44724194e+03, 3.41410449e+03, 1.90661475e+03, ],
    'CountWeightedFull_pt200to300'                                                   : [ 3.35384058e+03, 3.34463403e+03, 3.36502466e+03, ],
    'CountWeightedFullLHEWeightScale_pt200to300'                                     : [ 4.87444141e+03, 4.69886230e+03, 4.55932910e+03, 3.45886646e+03, 3.35384058e+03, 3.27140234e+03, 2.82778516e+03, 2.73586426e+03, 2.66273608e+03, ],
    'CountWeightedFullLHEEnvelope_pt200to300'                                        : [ 4.92294873e+03, 2.62599365e+03, ],
    'CountWeightedFullPSWeight_pt200to300'                                           : [ 3.35387622e+03, 3.36144507e+03, 5.00730615e+03, 3.35116064e+03, 3.32504639e+03, 1.85109448e+03, ],
    'CountWeightedFullPSWeightOriginalXWGTUP_pt200to300'                             : [ 1.85773596e+03, 1.86001660e+03, 2.76714453e+03, 1.85440759e+03, 1.83661853e+03, 1.02565515e+03, ],
    'Count_ptGt300'                                                                  : [ 2981, ],
    'CountWeighted_ptGt300'                                                          : [ 2.72514648e+03, 2.73719702e+03, 2.71166162e+03, ],
    'CountWeightedLHEWeightScale_ptGt300'                                            : [ 3.96811157e+03, 3.82225854e+03, 3.69984399e+03, 2.80651221e+03, 2.72514648e+03, 2.65579761e+03, 2.28276440e+03, 2.21124438e+03, 2.15037134e+03, ],
    'CountWeightedLHEEnvelope_ptGt300'                                               : [ 4.00391699e+03, 2.11685815e+03, ],
    'CountWeightedPSWeight_ptGt300'                                                  : [ 2.73040796e+03, 2.71317090e+03, 4.04115479e+03, 2.71757349e+03, 2.68279883e+03, 1.50355603e+03, ],
    'CountWeightedPSWeightOriginalXWGTUP_ptGt300'                                    : [ 1.58810596e+03, 1.58777173e+03, 2.35895361e+03, 1.58194739e+03, 1.55471472e+03, 8.69757080e+02, ],
    'CountWeightedFull_ptGt300'                                                      : [ 1.46594788e+03, 1.47243896e+03, 1.45871497e+03, ],
    'CountWeightedFullLHEWeightScale_ptGt300'                                        : [ 2.13462183e+03, 2.05616919e+03, 1.99031372e+03, 1.50972778e+03, 1.46594788e+03, 1.42866382e+03, 1.22798145e+03, 1.18951428e+03, 1.15676978e+03, ],
    'CountWeightedFullLHEEnvelope_ptGt300'                                           : [ 2.15388501e+03, 1.13874084e+03, ],
    'CountWeightedFullPSWeight_ptGt300'                                              : [ 1.46879822e+03, 1.45948291e+03, 2.17384424e+03, 1.46188513e+03, 1.44318652e+03, 8.08830933e+02, ],
    'CountWeightedFullPSWeightOriginalXWGTUP_ptGt300'                                : [ 8.54310669e+02, 8.54110474e+02, 1.26892346e+03, 8.50990112e+02, 8.36310852e+02, 4.67884796e+02, ],
    'Count_pt300to450'                                                               : [ 2266, ],
    'CountWeighted_pt300to450'                                                       : [ 2.10446313e+03, 2.11381934e+03, 2.09400073e+03, ],
    'CountWeightedLHEWeightScale_pt300to450'                                         : [ 3.05594482e+03, 2.94513135e+03, 2.85313599e+03, 2.16768042e+03, 2.10446313e+03, 2.05152100e+03, 1.76639514e+03, 1.71070117e+03, 1.66398608e+03, ],
    'CountWeightedLHEEnvelope_pt300to450'                                            : [ 3.08305005e+03, 1.64108630e+03, ],
    'CountWeightedPSWeight_pt300to450'                                               : [ 2.10947656e+03, 2.08278418e+03, 3.10656494e+03, 2.09755566e+03, 2.08067261e+03, 1.17334131e+03, ],
    'CountWeightedPSWeightOriginalXWGTUP_pt300to450'                                 : [ 1.20621777e+03, 1.19817554e+03, 1.78293079e+03, 1.20049207e+03, 1.18577637e+03, 6.67697693e+02, ],
    'CountWeightedFull_pt300to450'                                                   : [ 1.13202087e+03, 1.13705774e+03, 1.12640564e+03, ],
    'CountWeightedFullLHEWeightScale_pt300to450'                                     : [ 1.64385791e+03, 1.58426257e+03, 1.53477893e+03, 1.16602637e+03, 1.13202087e+03, 1.10355835e+03, 9.50169800e+02, 9.20218567e+02, 8.95095398e+02, ],
    'CountWeightedFullLHEEnvelope_pt300to450'                                        : [ 1.65844177e+03, 8.82775024e+02, ],
    'CountWeightedFullPSWeight_pt300to450'                                           : [ 1.13473193e+03, 1.12033655e+03, 1.67101477e+03, 1.12830542e+03, 1.11921509e+03, 6.31169983e+02, ],
    'CountWeightedFullPSWeightOriginalXWGTUP_pt300to450'                             : [ 6.48858276e+02, 6.44509766e+02, 9.59027466e+02, 6.45773865e+02, 6.37826721e+02, 3.59178650e+02, ],
    'Count_ptGt450'                                                                  : [ 715, ],
    'CountWeighted_ptGt450'                                                          : [ 6.20672607e+02, 6.23381104e+02, 6.17651978e+02, ],
    'CountWeightedLHEWeightScale_ptGt450'                                            : [ 9.12168579e+02, 8.77124268e+02, 8.46705383e+02, 6.38831970e+02, 6.20672607e+02, 6.04276917e+02, 5.16369080e+02, 5.00544464e+02, 4.86383820e+02, ],
    'CountWeightedLHEEnvelope_ptGt450'                                               : [ 9.20869141e+02, 4.75769623e+02, ],
    'CountWeightedPSWeight_ptGt450'                                                  : [ 6.20933105e+02, 6.30383057e+02, 9.34589844e+02, 6.20017212e+02, 6.02125610e+02, 3.30214539e+02, ],
    'CountWeightedPSWeightOriginalXWGTUP_ptGt450'                                    : [ 3.81888397e+02, 3.89596649e+02, 5.76024475e+02, 3.81454926e+02, 3.68938873e+02, 2.02059570e+02, ],
    'CountWeightedFull_ptGt450'                                                      : [ 3.33928345e+02, 3.35378754e+02, 3.32306458e+02, ],
    'CountWeightedFullLHEWeightScale_ptGt450'                                        : [ 4.90764435e+02, 4.71905182e+02, 4.55536072e+02, 3.43701263e+02, 3.33928345e+02, 3.25104919e+02, 2.77811005e+02, 2.69295197e+02, 2.61675049e+02, ],
    'CountWeightedFullLHEEnvelope_ptGt450'                                           : [ 4.95443115e+02, 2.55966125e+02, ],
    'CountWeightedFullPSWeight_ptGt450'                                              : [ 3.34066437e+02, 3.39144196e+02, 5.02828979e+02, 3.33579041e+02, 3.23971283e+02, 1.77660843e+02, ],
    'CountWeightedFullPSWeightOriginalXWGTUP_ptGt450'                                : [ 2.05452591e+02, 2.09600540e+02, 3.09895782e+02, 2.05216034e+02, 1.98484177e+02, 1.08706284e+02, ],
  }),
  ("nof_tree_events",                 13723),
  ("nof_db_events",                   7525991),
  ("fsize_local",                     90643363), # 90.64MB, avg file size 90.64MB
  ("fsize_db",                        469867184231), # 469.87GB, avg file size 2.10GB
  ("use_it",                          True),
  ("xsection",                        0.2118),
  ("genWeight",                       True),
  ("triggers",                        ['1e', '1mu', '2e', '2mu', '1e1mu', '3e', '3mu', '2e1mu', '1e2mu', '1e1tau', '1mu1tau', '2tau']),
  ("has_LHE",                         True),
  ("nof_PSweights",                   4),
  ("LHE_set",                         "LHA IDs 91400 - 91432 -> PDF4LHC15_nnlo_30_pdfas PDF set, expecting 33 weights (counted 33 weights)"),
  ("nof_reweighting",                 0),
  ("local_paths",
    [
      OD([
        ("path",      "/local/karl/ttHNtupleProduction/2018/2020Oct26_wPresel_nonNom_tth_sync/ntuples/ttHToNonbb_M125_powheg"),
        ("selection", "*"),
        ("blacklist", []),
      ]),
    ]
  ),
  ("missing_completely",           [
    # not computed
  ]),
  ("missing_from_superset",        [
    # not computed
  ]),
  ("missing_hlt_paths",            [

  ]),
  ("hlt_paths",                    [
    # not computed
  ]),
])

samples_2018["sum_events"] = [
]

