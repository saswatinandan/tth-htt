#ifndef tthAnalysis_HiggsToTauTau_hadTopTaggerAuxFunctions_h
#define tthAnalysis_HiggsToTauTau_hadTopTaggerAuxFunctions_h

#include <TString.h> // TString, Form
#include <boost/math/special_functions/sign.hpp> // boost::math::sign()
#include <map>


#include "tthAnalysis/HiggsToTauTau/interface/RecoJet.h" // RecoJet

std::map<int, Particle::LorentzVector>
isGenMatchedJetTripletVar(const std::vector<GenParticle> & genTopQuarks,
                          const std::vector<GenParticle> & genBJets,
                          const std::vector<GenParticle> & genWBosons,
                          const std::vector<GenParticle> & genWJets,
                          int mode);

bool
passWbosonMassVeto(const GenParticle * genWJetFromTop_lead,
                   const GenParticle * genWJetFromTop_sublead,
                   const GenParticle * genWBosonFromTop);

int
getType(size_t sizeHTTv2, size_t sizeFatW, size_t sizeResolved);

std::vector<double>
getBdiscr(std::vector<const RecoJet*> selJetsIt);

std::vector<size_t>
calRank(std::vector<double> & btag_discEnter);

void get_two_maximal(
  double &  MVA_1_mvaOutput,
  bool &    MVA_1_truth,
  double &  MVA_1_genTopPt,
  double &  MVA_1_recTopPt,
  double & MVA_2_mvaOutput,
  bool &   MVA_2_truth,
  double & MVA_2_genTopPt,
  double & MVA_2_recTopPt,
  double & MVA_med_mvaOutput,
  bool   & MVA_med_truth,
  double & MVA_med_genTopPt,
  double & MVA_med_recTopPt,
  double & Wj1_pt_1, double & Wj2_pt_1, double & b_pt_1,
  double & Wj1_pt_2, double & Wj2_pt_2, double & b_pt_2,
  double & Wj1_pt_med, double & Wj2_pt_med, double & b_pt_med,
  double MVA, double recTopPt, double genTopPt_teste, bool isGenMatched,
  double Wj1_pt, double Wj2_pt, double b_pt
);

#endif
