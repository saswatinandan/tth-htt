#include "tthAnalysis/HiggsToTauTau/interface/HadTopTagger_semi_boosted_AK8.h" // HadTopTagger_boosted

//#include "tthAnalysis/HiggsToTauTau/interface/HadTopKinFit.h" // HadTopKinFit
#include "tthAnalysis/HiggsToTauTau/interface/TMVAInterface.h" // TMVAInterface
#include "tthAnalysis/HiggsToTauTau/interface/XGBInterface.h" // XGBInterface
#include "tthAnalysis/HiggsToTauTau/interface/Particle.h" // Particle
#include "tthAnalysis/HiggsToTauTau/interface/hadTopTaggerAuxFunctions_internal.h" // isGenMatchedJetTriplet
#include "tthAnalysis/HiggsToTauTau/interface/hadTopTaggerAuxFunctions_geral.h" // kGenTop...
#include "DataFormats/Math/interface/deltaR.h" // deltaR
#include <TLorentzVector.h> // TLorentzVector


HadTopTagger_semi_boosted_AK8::HadTopTagger_semi_boosted_AK8(void)
  : mva_xgb_HTT_CSVsort3rd_AK8_(nullptr)
  , mva_xgb_HTT_CSVsort3rd_AK8_withKinFit_(nullptr)
{

  //std::string mvaFileNameHTT_CSVsort3rd = "tthAnalysis/HiggsToTauTau/data/hadTopTagger_2018/HTT_HadTopTagger_2017_boosted_final_XGB_ntrees_2500_deph_4_lr_0o03_btagSort3rd_nvar13_cat_2.pkl";
  //std::string mvaFileNameHTT_CSVsort3rd_withKinFit = "tthAnalysis/HiggsToTauTau/data/hadTopTagger_2018/HTT_HadTopTagger_2017_boosted_final_XGB_ntrees_2500_deph_4_lr_0o03_btagSort3rd_nvar16_cat_2_withKinFit_sel_5.pkl";

  std::string mvaFileNameHTT_CSVsort3rd = "tthAnalysis/HiggsToTauTau/data/hadTopTagger_2018/no_mass_cut_on_new/HTT_HadTopTagger_2017_final_nomasscut_AK8_XGB_ntrees_2500_deph_4_lr_0o03_btagSort3rd_nvar12_cat_2_AK8.pkl";
  std::string mvaFileNameHTT_CSVsort3rd_withKinFit = "tthAnalysis/HiggsToTauTau/data/hadTopTagger_2018/no_mass_cut_on_new/HTT_HadTopTagger_2017_final_nomasscut_AK8_XGB_ntrees_2500_deph_4_lr_0o03_btagSort3rd_nvar15_cat_2_withKinFit_sel_9.pkl";

  mvaInputsHTT_AK8Sort =  {
    "massW_SD", "tau21W", "btagDisc_b",
    "m_Wj1Wj2_div_m_bWj1Wj2", "dR_Wj1Wj2", "m_bWj1Wj2", "pT_bWj1Wj2", //"m_bWj2",
    "mass_Wj1", "pT_Wj2", "mass_Wj2", "pT_b", "mass_b"
  };
  mva_xgb_HTT_CSVsort3rd_AK8_ = new XGBInterface(
    mvaFileNameHTT_CSVsort3rd, mvaInputsHTT_AK8Sort
  );

  mvaInputsHTT_AK8Sort_withKinFit =  {
    "massW_SD", "tau21W", "btagDisc_b",
    "m_Wj1Wj2_div_m_bWj1Wj2", "dR_Wj1Wj2", "m_bWj1Wj2", "pT_bWj1Wj2", //"m_bWj2",
    "mass_Wj1", "pT_Wj2", "mass_Wj2", "pT_b", "mass_b",
    "kinFit_pT_b_o_pT_b", "kinFit_pT_Wj2_o_pT_Wj2", "nllKinFit"
  };
  mva_xgb_HTT_CSVsort3rd_AK8_withKinFit_ = new XGBInterface(
    mvaFileNameHTT_CSVsort3rd_withKinFit, mvaInputsHTT_AK8Sort_withKinFit
  );

}

HadTopTagger_semi_boosted_AK8::~HadTopTagger_semi_boosted_AK8()
{
  //delete kinFit_;
  delete mva_xgb_HTT_CSVsort3rd_AK8_;
  delete mva_xgb_HTT_CSVsort3rd_AK8_withKinFit_;
}

std::map<int, double>
HadTopTagger_semi_boosted_AK8::operator()(
  const RecoJetAK8 & jet_ptrsAK12, const RecoJet & b_jet_candidate,
  bool & calculate_matching, bool & isGenMatched,
  double & genTopPt,
  std::map<int, Particle::LorentzVector> genVar, std::map<int, Particle::LorentzVector> genVarAnti
)
{
  std::map<int, double> result = {
    { kXGB_semi_boosted_AK8_with_kinFit,  0. },
    { kXGB_semi_boosted_AK8_no_kinFit,  0. }
  };

  if ( !(jet_ptrsAK12.subJet1() && jet_ptrsAK12.subJet2()) ) {
    result[kXGB_semi_boosted_AK8_with_kinFit] = 0.;
    result[kXGB_semi_boosted_AK8_with_kinFit] = 0.;
    return result;
  }

  Particle::LorentzVector recBJet, recWJet1, recWJet2 ;
  recBJet = b_jet_candidate.p4();
  recWJet1 = jet_ptrsAK12.subJet1()->p4();
  recWJet2 = jet_ptrsAK12.subJet2()->p4();

  double m_bWj1Wj2 = (recWJet1 + recWJet2 + recBJet).mass();
  double  m_Wj1Wj2 = (recWJet1 + recWJet2).mass();
  //if ( !(m_bWj1Wj2 > 75. && m_bWj1Wj2 < 275. && m_Wj1Wj2 > 150.) ) {
  //  result[kXGB_semi_boosted_with_kinFit] = 0.;
  //  result[kXGB_semi_boosted_with_kinFit] = 0.;
  //  return result;
  //}
  //kinFit_->fit(recBJet, recWJet1, recWJet2);

  if ( calculate_matching ) {
    int typeTop = 2;
    std::map<int, bool> genMatchingTop = isGenMatchedJetTriplet(
      recBJet, recWJet1, recWJet2,
      genVar[kGenTop], genVar[kGenTopB], genVar[kGenTopW], genVar[kGenTopWj1], genVar[kGenTopWj2],
      kGenTop, typeTop, jet_ptrsAK12.p4()
    );
    std::map<int, bool> genMatchingAntiTop = isGenMatchedJetTriplet(
      recBJet, recWJet1, recWJet2,
      genVarAnti[kGenTop], genVarAnti[kGenTopB], genVarAnti[kGenTopW], genVarAnti[kGenTopWj1], genVarAnti[kGenTopWj2],
      kGenAntiTop, typeTop, jet_ptrsAK12.p4()
    );
    if(genMatchingTop[kGenMatchedTriplet]) { genTopPt = genVar[kGenTop].pt();}
    if(genMatchingAntiTop[kGenMatchedTriplet]) { genTopPt = genVarAnti[kGenTop].pt();}

    isGenMatched = (genMatchingTop[kGenMatchedTriplet] || genMatchingAntiTop[kGenMatchedTriplet]);
    //fatjet_isGenMatched = (genMatchingTop[kGenMatchedFatJet] || genMatchingAntiTop[kGenMatchedFatJet]);
  } // close gen matching


  // (*jetIter)->msoftdrop()
  //"massW_SD", "tau21W", "btagDisc_b",
  // "m_Wj1Wj2_div_m_bWj1Wj2", "dR_Wj1Wj2", "m_bWj1Wj2", "pT_bWj1Wj2", "m_bWj2",
  // "mass_Wj1", "pT_Wj2", "mass_Wj2", "pT_b", "mass_b",
  // "kinFit_pT_b_o_pT_b", "kinFit_pT_Wj2_o_pT_Wj2", "nllKinFit"
  mvaInputsHTT_AK8["massW_SD"]           = jet_ptrsAK12.msoftdrop();
  mvaInputsHTT_AK8["tau21W"]           = jet_ptrsAK12.tau2()/jet_ptrsAK12.tau1();
  mvaInputsHTT_AK8["btagDisc_b"]         = b_jet_candidate.BtagCSV();
  mvaInputsHTT_AK8["m_Wj1Wj2_div_m_bWj1Wj2"]  = m_Wj1Wj2/m_bWj1Wj2;
  mvaInputsHTT_AK8["dR_Wj1Wj2"]          = deltaR(recWJet1, recWJet2);
  mvaInputsHTT_AK8["m_bWj1Wj2"]          = (recWJet1 + recWJet2 + recBJet).mass();
  mvaInputsHTT_AK8["pT_bWj1Wj2"]          = (recWJet1 + recWJet2 + recBJet).pt();
  //mvaInputsHTT_AK8["m_bWj2"]             = (recBJet + recWJet2).mass();
  mvaInputsHTT_AK8["mass_Wj1"]           = recWJet1.mass();
  mvaInputsHTT_AK8["pT_Wj2"]             = recWJet2.pt();
  mvaInputsHTT_AK8["mass_Wj2"]           = recWJet2.mass();
  mvaInputsHTT_AK8["pT_b"]               = recBJet.pt();
  mvaInputsHTT_AK8["mass_b"]             = recBJet.mass();
  const double HTT_CSVsort3rd = (*mva_xgb_HTT_CSVsort3rd_AK8_)(mvaInputsHTT_AK8);
  result[kXGB_semi_boosted_AK8_no_kinFit] = HTT_CSVsort3rd;
  //std::cout << "semi-boosted HTT_CSVsort3rd " << HTT_CSVsort3rd << jet_ptrsAK12.tau2()/jet_ptrsAK12.tau1() << std::endl;

  /*
  // "kinFit_pT_b_o_pT_b", "kinFit_pT_Wj2_o_pT_Wj2", "nllKinFit"
  mvaInputsHTT_AK8_withKinFit["massW_SD"]           = jet_ptrsAK12.msoftdrop();
  mvaInputsHTT_AK8_withKinFit["tau21W"]           = jet_ptrsAK12.tau2()/jet_ptrsAK12.tau1();
  mvaInputsHTT_AK8_withKinFit["btagDisc_b"]         = b_jet_candidate.BtagCSV();
  mvaInputsHTT_AK8_withKinFit["m_Wj1Wj2_div_m_bWj1Wj2"]  = m_Wj1Wj2/m_bWj1Wj2;
  mvaInputsHTT_AK8_withKinFit["dR_Wj1Wj2"]          = deltaR(recWJet1, recWJet2);
  mvaInputsHTT_AK8_withKinFit["m_bWj1Wj2"]          = (recWJet1 + recWJet2 + recBJet).mass();
  mvaInputsHTT_AK8_withKinFit["pT_bWj1Wj2"]          = (recWJet1 + recWJet2 + recBJet).pt();
  mvaInputsHTT_AK8_withKinFit["m_bWj2"]             = (recBJet + recWJet2).mass();
  mvaInputsHTT_AK8_withKinFit["mass_Wj1"]           = recWJet1.mass();
  mvaInputsHTT_AK8_withKinFit["pT_Wj2"]             = recWJet2.pt();
  mvaInputsHTT_AK8_withKinFit["mass_Wj2"]           = recWJet2.mass();
  mvaInputsHTT_AK8_withKinFit["pT_b"]               = recBJet.pt();
  mvaInputsHTT_AK8_withKinFit["mass_b"]             = recBJet.mass();
  mvaInputsHTT_AK8_withKinFit["nllKinFit"]          = kinFit_->nll();
  mvaInputsHTT_AK8_withKinFit["kinFit_pT_b_o_pT_b"] =  kinFit_->fittedBJet().pt() / recBJet.pt();
  mvaInputsHTT_AK8_withKinFit["kinFit_pT_Wj2_o_pT_Wj2"] =  kinFit_->fittedWJet2().pt() / recWJet2.pt();
  const double HTT_CSVsort3rd_withKinFit = (*mva_xgb_HTT_CSVsort3rd_AK8_withKinFit_)(mvaInputsHTT_AK8_withKinFit);
  result[kXGB_semi_boosted_with_kinFit] = HTT_CSVsort3rd_withKinFit;
  //std::cout << "semi-boosted HTT_CSVsort3rd_withKinFit " << HTT_CSVsort3rd_withKinFit << kinFit_->nll() << std::endl;
  */
  result[kXGB_semi_boosted_AK8_with_kinFit] = 0.;

  return result;
}

const std::map<std::string, double> &
HadTopTagger_semi_boosted_AK8::mvaInputs() const
{
  return mvaInputsHTT_AK8;
}

//const HadTopKinFit *
//HadTopTagger_semi_boosted::kinFit() const
//{
//  return kinFit_;
//}
