#include "tthAnalysis/HiggsToTauTau/interface/HadTopTagger.h" // HadTopTagger

#include "tthAnalysis/HiggsToTauTau/interface/HadTopKinFit.h" // HadTopKinFit
#include "tthAnalysis/HiggsToTauTau/interface/TMVAInterface.h" // TMVAInterface
#include "tthAnalysis/HiggsToTauTau/interface/XGBInterface.h" // XGBInterface
#include "tthAnalysis/HiggsToTauTau/interface/Particle.h" // Particle
#include "tthAnalysis/HiggsToTauTau/interface/hadTopTaggerAuxFunctions_internal.h" // isGenMatchedJetTriplet
#include "tthAnalysis/HiggsToTauTau/interface/hadTopTaggerAuxFunctions_geral.h" // kGenTop...
#include "DataFormats/Math/interface/deltaR.h" // deltaR
#include <TLorentzVector.h> // TLorentzVector

HadTopTagger::HadTopTagger(void)
  : kinFit_(new HadTopKinFit())
  , mva_hadTopTagger_xgb_withKinFit_(nullptr)
  , mva_hadTopTagger_xgb_HTT_multilep_(nullptr)
  , mva_xgb_HTT_CSVsort3rd_(nullptr)
  , mva_xgb_HTT_CSVsort3rd_withKinFit_(nullptr)
  , mva_xgb_HTT_highestCSV_(nullptr)
  , mva_xgb_HTT_highestCSV_withKinFit_(nullptr)
{

  std::string mvaFileNameWithKinFit = "tthAnalysis/HiggsToTauTau/data/HadTopTagger_resolved_XGB_CSV_sort_withKinFit.xml";
  std::string mvaFileNameHTT_multilep = "tthAnalysis/HiggsToTauTau/data/multilep_BDTs_2018/resTop_xgb_csv_order_qgl.xml";
  //std::string mvaFileNameHTT_CSVsort3rd = "tthAnalysis/HiggsToTauTau/data/hadTopTagger_2018/HTT_HadTopTagger_2017_boosted_final_XGB_ntrees_2500_deph_5_lr_0o03_btagSort3rd_nvar17_cat_3.pkl";
  //std::string mvaFileNameHTT_CSVsort3rd_withKinFit = "tthAnalysis/HiggsToTauTau/data/hadTopTagger_2018/HTT_HadTopTagger_2017_boosted_final_XGB_ntrees_2500_deph_5_lr_0o03_btagSort3rd_nvar20_cat_3_withKinFit_sel_5.pkl";
  //std::string mvaFileNameHTT_highestCSV = "tthAnalysis/HiggsToTauTau/data/hadTopTagger_2018/HTT_HadTopTagger_2017_boosted_final_XGB_ntrees_2500_deph_5_lr_0o010_btagHighest_nvar17_cat_3.pkl";
  //std::string mvaFileNameHTT_highestCSV_withKinFit = "tthAnalysis/HiggsToTauTau/data/hadTopTagger_2018/HTT_HadTopTagger_2017_boosted_final_XGB_ntrees_2500_deph_4_lr_0o010_btagHighest_nvar20_cat_3_withKinFit_sel_6.pkl";

  std::string mvaFileNameHTT_CSVsort3rd = "tthAnalysis/HiggsToTauTau/data/hadTopTagger_2018/no_mass_cut_on_new/HTT_HadTopTagger_2017_final_nomasscut_XGB_ntrees_2500_deph_5_lr_0o03_btagSort3rd_nvar17_cat_3.pkl";
  std::string mvaFileNameHTT_CSVsort3rd_withKinFit = "tthAnalysis/HiggsToTauTau/data/hadTopTagger_2018/no_mass_cut_on_new/HTT_HadTopTagger_2017_final_nomasscut_XGB_ntrees_2500_deph_5_lr_0o03_btagSort3rd_nvar20_cat_3_withKinFit_sel_9.pkl";
  std::string mvaFileNameHTT_highestCSV = "tthAnalysis/HiggsToTauTau/data/hadTopTagger_2018/no_mass_cut_on_new/HTT_HadTopTagger_2017_final_nomasscut_XGB_ntrees_2500_deph_5_lr_0o010_btagHighest_nvar17_cat_3.pkl";
  std::string mvaFileNameHTT_highestCSV_withKinFit = "tthAnalysis/HiggsToTauTau/data/hadTopTagger_2018/no_mass_cut_on_new/HTT_HadTopTagger_2017_final_nomasscut_XGB_ntrees_2500_deph_4_lr_0o010_btagHighest_nvar20_cat_3_withKinFit_sel_8.pkl";

  mvaInputsWithKinFitSort = {
    "CSV_b",
    "qg_Wj2",
    "pT_bWj1Wj2",
    "m_Wj1Wj2",
    "nllKinFit",
    "pT_b_o_kinFit_pT_b",
    "pT_Wj2"
  };
  mva_hadTopTagger_xgb_withKinFit_ = new TMVAInterface(
    mvaFileNameWithKinFit, mvaInputsWithKinFitSort
  );
  mva_hadTopTagger_xgb_withKinFit_->enableBDTTransform();

  mvaInputsHTT_multilepSort =  {
    "var_b_pt", "var_b_mass", "var_b_csv",
    "var_wj1_pt", "var_wj1_mass", "var_wj1_csv", "var_wj1_qgl",
    "var_wj2_pt", "var_wj2_mass", "var_wj2_csv", "var_wj2_qgl",
    "var_b_wj1_deltaR", "var_b_wj1_mass", "var_b_wj2_deltaR",
    "var_b_wj2_mass", "var_wcand_deltaR", "var_wcand_mass", "var_b_wcand_deltaR", "var_topcand_mass"
  };
  mva_hadTopTagger_xgb_HTT_multilep_ = new TMVAInterface(
    mvaFileNameHTT_multilep, mvaInputsHTT_multilepSort
  );

  mvaInputsHTTSort =  {
    "btagDisc_b", "btagDisc_Wj1", "btagDisc_Wj2", "qg_Wj1", "qg_Wj2",
    "m_Wj1Wj2_div_m_bWj1Wj2", "pT_Wj1Wj2", "dR_Wj1Wj2", "m_bWj1Wj2", "dR_bW", "m_bWj1", "m_bWj2",
    "mass_Wj1", "pT_Wj2", "mass_Wj2", "pT_b", "mass_b"
  };
  mva_xgb_HTT_CSVsort3rd_ = new XGBInterface(
    mvaFileNameHTT_CSVsort3rd, mvaInputsHTTSort
  );
  mva_xgb_HTT_highestCSV_ = new XGBInterface(
    mvaFileNameHTT_highestCSV, mvaInputsHTTSort
  );

  mvaInputsHTTSort =  {
    "btagDisc_b", "btagDisc_Wj1", "btagDisc_Wj2", "qg_Wj1", "qg_Wj2",
    "m_Wj1Wj2_div_m_bWj1Wj2", "pT_Wj1Wj2", "dR_Wj1Wj2", "m_bWj1Wj2", "dR_bW", "m_bWj1", "m_bWj2",
    "mass_Wj1", "pT_Wj2", "mass_Wj2", "pT_b", "mass_b",
    "kinFit_pT_b_o_pT_b", "kinFit_pT_Wj2_o_pT_Wj2", "nllKinFit"
  };
  mva_xgb_HTT_CSVsort3rd_withKinFit_ = new XGBInterface(
    mvaFileNameHTT_CSVsort3rd_withKinFit, mvaInputsHTTSort
  );
  mva_xgb_HTT_highestCSV_withKinFit_ = new XGBInterface(
    mvaFileNameHTT_highestCSV_withKinFit, mvaInputsHTTSort
  );

}

HadTopTagger::~HadTopTagger()
{
  delete kinFit_;
  delete mva_hadTopTagger_xgb_withKinFit_;
  delete mva_hadTopTagger_xgb_HTT_multilep_;
  delete mva_xgb_HTT_CSVsort3rd_;
  delete mva_xgb_HTT_CSVsort3rd_withKinFit_;
  delete mva_xgb_HTT_highestCSV_;
  delete mva_xgb_HTT_highestCSV_withKinFit_;
}

std::map<int, double>
HadTopTagger::operator()(const RecoJet & recBJet,
                         const RecoJet & recWJet1,
                         const RecoJet & recWJet2,
                         bool & calculate_matching, bool & isGenMatched,
                         double & genTopPt,
                         std::map<int, Particle::LorentzVector> genVar, std::map<int, Particle::LorentzVector> genVarAnti,
                         bool massCut = true
                       )
{
  std::map<int, double> result = {
    { kXGB_with_kinFit,  0. },
    { kXGB_multilep,  -1. },
    { kXGB_CSVsort3rd,  0. },
    { kXGB_CSVsort3rd_withKinFit,  0. },
    { kXGB_highestCSV,  0. },
    { kXGB_highestCSV_withKinFit,  0. }
  };

  const Particle::LorentzVector p4_bWj1Wj2 = recBJet.p4() + recWJet1.p4() + recWJet2.p4();
  const Particle::LorentzVector p4_Wj1Wj2  = recWJet1.p4() + recWJet2.p4();
  // apply the mass selections
  Particle::LorentzVector p4_Wj1  = recWJet1.p4();
  Particle::LorentzVector p4_Wj2  = recWJet2.p4();
  Particle::LorentzVector p4_b    = recBJet.p4();
  kinFit_->fit(p4_b, p4_Wj1, p4_Wj2);

  if ( calculate_matching ) {
    // calculate matching
    std::map<int, bool> genMatchingTop = isGenMatchedJetTriplet(
      recBJet.p4(), recWJet1.p4(), recWJet2.p4(),
      genVar[kGenTop], genVar[kGenTopB], genVar[kGenTopW], genVar[kGenTopWj1], genVar[kGenTopWj2],
      kGenTop
    );
    std::map<int, bool> genMatchingAntiTop = isGenMatchedJetTriplet(
      recBJet.p4(), recWJet1.p4(), recWJet2.p4(),
      genVarAnti[kGenTop], genVarAnti[kGenTopB], genVarAnti[kGenTopW], genVarAnti[kGenTopWj1], genVarAnti[kGenTopWj2],
      kGenAntiTop
    );
    if(genMatchingTop[kGenMatchedTriplet]) { genTopPt = genVar[kGenTop].pt(); }
    if(genMatchingAntiTop[kGenMatchedTriplet]) { genTopPt = genVarAnti[kGenTop].pt(); }
    isGenMatched = (genMatchingTop[kGenMatchedTriplet] || genMatchingAntiTop[kGenMatchedTriplet]);
  }

  mvaInputsWithKinFit["CSV_b"]              = recBJet.BtagCSV();
  mvaInputsWithKinFit["qg_Wj2"]             = recWJet2.QGDiscr();
  mvaInputsWithKinFit["pT_bWj1Wj2"]         = p4_bWj1Wj2.pt();
  mvaInputsWithKinFit["m_Wj1Wj2"]           = p4_Wj1Wj2.mass();
  mvaInputsWithKinFit["nllKinFit"]          = kinFit_->nll();
  mvaInputsWithKinFit["pT_b_o_kinFit_pT_b"] = recBJet.pt() / kinFit_->fittedBJet().pt();
  mvaInputsWithKinFit["pT_Wj2"]             = recWJet2.pt();
  const double HTT_WithKin_xgb = (*mva_hadTopTagger_xgb_withKinFit_)(mvaInputsWithKinFit);
  result[kXGB_with_kinFit] = HTT_WithKin_xgb;
  //std::cout << " HTT_2016 "<< HTT_WithKin_xgb << " " << kinFit_->nll() << std::endl;

  /*
  if ( massCut && !(p4_bWj1Wj2.mass() > 75. && p4_bWj1Wj2.mass() < 275. && p4_Wj1Wj2.mass() < 150.)) {
    result[kXGB_multilep] = -1.;
    result[kXGB_CSVsort3rd] = 0.;
    result[kXGB_highestCSV] = 0.;
    result[kXGB_CSVsort3rd_withKinFit] = 0.;
    result[kXGB_highestCSV_withKinFit] = 0.;
    return result;
  }
  */
  // "btagDisc_b", "btagDisc_Wj1", "btagDisc_Wj2", "qg_Wj1", "qg_Wj2",
  //"m_Wj1Wj2_div_m_bWj1Wj2", "pT_Wj1Wj2", "dR_Wj1Wj2", "m_bWj1Wj2", "dR_bW","m_bWj1", "m_bWj2",
  //"mass_Wj1", "pT_Wj2", "mass_Wj2", "pT_b", "mass_b"
  mvaInputsHTT["btagDisc_b"]         = recBJet.BtagCSV();
  mvaInputsHTT["btagDisc_Wj1"]       = recWJet1.BtagCSV();
  mvaInputsHTT["btagDisc_Wj2"]       = recWJet2.BtagCSV();
  mvaInputsHTT["qg_Wj1"]             = recWJet2.QGDiscr();
  mvaInputsHTT["qg_Wj2"]             = recWJet2.QGDiscr();
  mvaInputsHTT["m_Wj1Wj2_div_m_bWj1Wj2"]  = p4_Wj1Wj2.mass()/p4_bWj1Wj2.mass();
  mvaInputsHTT["pT_Wj1Wj2"]          = p4_Wj1Wj2.pt();
  mvaInputsHTT["dR_Wj1Wj2"]          = deltaR(recWJet1.p4(),recWJet2.p4());
  mvaInputsHTT["m_bWj1Wj2"]          = p4_bWj1Wj2.mass();
  mvaInputsHTT["dR_bW"]              = deltaR(recBJet.p4(), recWJet1.p4()+recWJet2.p4());
  mvaInputsHTT["m_bWj1"]            = deltaR(recBJet.p4(), recWJet1.p4());
  mvaInputsHTT["m_bWj2"]            = deltaR(recBJet.p4(), recWJet2.p4());
  mvaInputsHTT["mass_Wj1"]           = recWJet1.mass();
  mvaInputsHTT["pT_Wj2"]             = recWJet2.pt();
  mvaInputsHTT["mass_Wj2"]           = recWJet2.mass();
  mvaInputsHTT["pT_b"]               = recBJet.pt();
  mvaInputsHTT["mass_b"]             = recBJet.mass();
  const double HTT_CSVsort3rd = (*mva_xgb_HTT_CSVsort3rd_)(mvaInputsHTT);
  result[kXGB_CSVsort3rd] = HTT_CSVsort3rd;
  const double HTT_highestCSV = (*mva_xgb_HTT_highestCSV_)(mvaInputsHTT);
  result[kXGB_highestCSV] = HTT_highestCSV;
  //std::cout << " HTT_CSVsort3rd " << HTT_CSVsort3rd << std::endl;
  //std::cout << " HTT_highestCSV " << HTT_highestCSV << std::endl;
  //*/

  // "kinFit_pT_b_o_pT_b", "kinFit_pT_Wj2_o_pT_Wj2", "nllKinFit"
  mvaInputsHTT_withKinFit["btagDisc_b"]         = recBJet.BtagCSV();
  mvaInputsHTT_withKinFit["btagDisc_Wj1"]       = recWJet1.BtagCSV();
  mvaInputsHTT_withKinFit["btagDisc_Wj2"]       = recWJet2.BtagCSV();
  mvaInputsHTT_withKinFit["qg_Wj1"]             = recWJet2.QGDiscr();
  mvaInputsHTT_withKinFit["qg_Wj2"]             = recWJet2.QGDiscr();
  mvaInputsHTT_withKinFit["m_Wj1Wj2_div_m_bWj1Wj2"]  = p4_Wj1Wj2.mass()/p4_bWj1Wj2.mass();
  mvaInputsHTT_withKinFit["pT_Wj1Wj2"]          = p4_Wj1Wj2.pt();
  mvaInputsHTT_withKinFit["dR_Wj1Wj2"]          = deltaR(recWJet1.p4(),recWJet2.p4());
  mvaInputsHTT_withKinFit["m_bWj1Wj2"]          = p4_bWj1Wj2.mass();
  mvaInputsHTT_withKinFit["dR_bW"]              = deltaR(recBJet.p4(), recWJet1.p4()+recWJet2.p4());
  mvaInputsHTT_withKinFit["m_bWj1"]            = deltaR(recBJet.p4(), recWJet1.p4());
  mvaInputsHTT_withKinFit["m_bWj2"]            = deltaR(recBJet.p4(), recWJet2.p4());
  mvaInputsHTT_withKinFit["mass_Wj1"]           = recWJet1.mass();
  mvaInputsHTT_withKinFit["pT_Wj2"]             = recWJet2.pt();
  mvaInputsHTT_withKinFit["mass_Wj2"]           = recWJet2.mass();
  mvaInputsHTT_withKinFit["pT_b"]               = recBJet.pt();
  mvaInputsHTT_withKinFit["mass_b"]             = recBJet.mass();
  mvaInputsHTT_withKinFit["nllKinFit"]          = kinFit_->nll();
  mvaInputsHTT_withKinFit["kinFit_pT_b_o_pT_b"] =  kinFit_->fittedBJet().pt() / recBJet.pt();
  mvaInputsHTT_withKinFit["kinFit_pT_Wj2_o_pT_Wj2"] =  kinFit_->fittedWJet2().pt() / recWJet2.pt();
  const double HTT_CSVsort3rd_withKinFit = (*mva_xgb_HTT_CSVsort3rd_withKinFit_)(mvaInputsHTT_withKinFit);
  result[kXGB_CSVsort3rd_withKinFit] = HTT_CSVsort3rd_withKinFit;
  const double HTT_highestCSV_withKinFit = (*mva_xgb_HTT_highestCSV_withKinFit_)(mvaInputsHTT_withKinFit);
  result[kXGB_highestCSV_withKinFit] = HTT_highestCSV_withKinFit;
  //std::cout << " HTT_CSVsort3rd_withKinFit " << HTT_CSVsort3rd_withKinFit << std::endl;
  //std::cout << " HTT_highestCSV_withKinFit " << HTT_highestCSV_withKinFit << std::endl;
  //*/

  if ( massCut && !(p4_bWj1Wj2.mass() > 75. && p4_bWj1Wj2.mass() < 275.)) {
    result[kXGB_multilep] = -1.;
  } else {
    ///*
    mvaInputsHTT_multilep["var_b_pt"]             = recBJet.pt();
    mvaInputsHTT_multilep["var_b_mass"]           = recBJet.p4().mass();
    mvaInputsHTT_multilep["var_b_csv"]            = recBJet.BtagCSV();
    mvaInputsHTT_multilep["var_wj1_pt"]           = recWJet1.pt();
    mvaInputsHTT_multilep["var_wj1_mass"]         = recWJet1.p4().mass();
    mvaInputsHTT_multilep["var_wj1_csv"]          = recWJet1.BtagCSV();
    mvaInputsHTT_multilep["var_wj1_qgl"]          = recWJet1.QGDiscr();
    mvaInputsHTT_multilep["var_wj2_pt"]           = recWJet2.pt();
    mvaInputsHTT_multilep["var_wj2_mass"]         = recWJet2.p4().mass();
    mvaInputsHTT_multilep["var_wj2_csv"]          = recWJet2.BtagCSV();
    mvaInputsHTT_multilep["var_wj2_qgl"]          = recWJet2.QGDiscr();
    mvaInputsHTT_multilep["var_b_wj1_deltaR"]     = deltaR(recBJet.p4(), recWJet1.p4());
    mvaInputsHTT_multilep["var_b_wj1_mass"]       = (recBJet.p4()+recWJet1.p4()).mass();
    mvaInputsHTT_multilep["var_b_wj2_deltaR"]     = deltaR(recBJet.p4(), recWJet2.p4());
    mvaInputsHTT_multilep["var_b_wj2_mass"]       = (recBJet.p4()+recWJet2.p4()).mass();
    mvaInputsHTT_multilep["var_wcand_deltaR"]     = deltaR(recWJet2.p4(), recWJet1.p4());
    mvaInputsHTT_multilep["var_wcand_mass"]       = (recWJet2.p4()+recWJet1.p4()).mass();
    mvaInputsHTT_multilep["var_b_wcand_deltaR"]   = deltaR(recBJet.p4(), recWJet1.p4()+recWJet2.p4());
    mvaInputsHTT_multilep["var_topcand_mass"]     = (recBJet.p4()+ recWJet1.p4()+recWJet2.p4()).mass();
    const double HTT_multilep = (*mva_hadTopTagger_xgb_HTT_multilep_)(mvaInputsHTT_multilep);
    result[kXGB_multilep] = HTT_multilep;
    //std::cout << " HTT_multilep " << HTT_multilep << std::endl;
  }

  return result;
}

const std::map<std::string, double> &
HadTopTagger::mvaInputs() const
{
  return mvaInputsWithKinFit;
}

const HadTopKinFit *
HadTopTagger::kinFit() const
{
  return kinFit_;
}
