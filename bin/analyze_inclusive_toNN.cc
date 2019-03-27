#include "FWCore/ParameterSet/interface/ParameterSet.h" // edm::ParameterSet
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h" // edm::readPSetsFrom()
#include "FWCore/Utilities/interface/Exception.h" // cms::Exception
#include "PhysicsTools/FWLite/interface/TFileService.h" // fwlite::TFileService
#include "DataFormats/FWLite/interface/InputSource.h" // fwlite::InputSource
#include "DataFormats/FWLite/interface/OutputFiles.h" // fwlite::OutputFiles
#include "DataFormats/Math/interface/LorentzVector.h" // math::PtEtaPhiMLorentzVector
#include "DataFormats/Math/interface/deltaR.h" // deltaR

#include <TLorentzVector.h> // TLorentzVector
#include <TVector3.h>
#include <TBenchmark.h> // TBenchmark
#include <TString.h> // TString, Form
#include <TError.h> // gErrorAbortLevel, kError

#include "tthAnalysis/HiggsToTauTau/interface/RecoLepton.h" // RecoLepton
#include "tthAnalysis/HiggsToTauTau/interface/RecoHadTau.h" // RecoHadTau
#include "tthAnalysis/HiggsToTauTau/interface/RecoJet.h" // RecoJet
#include "tthAnalysis/HiggsToTauTau/interface/RecoMEt.h" // RecoMEt
#include "tthAnalysis/HiggsToTauTau/interface/GenLepton.h" // GenLepton
#include "tthAnalysis/HiggsToTauTau/interface/GenJet.h" // GenJet
#include "tthAnalysis/HiggsToTauTau/interface/GenHadTau.h" // GenHadTau
#include "tthAnalysis/HiggsToTauTau/interface/EventInfo.h" // EventInfo
#include "tthAnalysis/HiggsToTauTau/interface/TMVAInterface.h" // TMVAInterface
#include "tthAnalysis/HiggsToTauTau/interface/XGBInterface.h" // XGBInterface
#include "tthAnalysis/HiggsToTauTau/interface/mvaAuxFunctions.h" // check_mvaInputs, get_mvaInputVariables
#include "tthAnalysis/HiggsToTauTau/interface/mvaInputVariables.h" // auxiliary functions for computing input variables of the MVA used for signal extraction in the inclusive_toNN category
#include "tthAnalysis/HiggsToTauTau/interface/LeptonFakeRateInterface.h" // LeptonFakeRateInterface
#include "tthAnalysis/HiggsToTauTau/interface/JetToTauFakeRateInterface.h" // JetToTauFakeRateInterface
#include "tthAnalysis/HiggsToTauTau/interface/RecoElectronReader.h" // RecoElectronReader
#include "tthAnalysis/HiggsToTauTau/interface/RecoMuonReader.h" // RecoMuonReader
#include "tthAnalysis/HiggsToTauTau/interface/RecoHadTauReader.h" // RecoHadTauReader
#include "tthAnalysis/HiggsToTauTau/interface/RecoJetReader.h" // RecoJetReader
#include "tthAnalysis/HiggsToTauTau/interface/RecoMEtReader.h" // RecoMEtReader
#include "tthAnalysis/HiggsToTauTau/interface/MEtFilterReader.h" // MEtFilterReader
#include "tthAnalysis/HiggsToTauTau/interface/GenLeptonReader.h" // GenLeptonReader
#include "tthAnalysis/HiggsToTauTau/interface/GenParticleReader.h" // GenParticleReader
#include "tthAnalysis/HiggsToTauTau/interface/GenHadTauReader.h" // GenHadTauReader
#include "tthAnalysis/HiggsToTauTau/interface/GenPhotonReader.h" // GenPhotonReader
#include "tthAnalysis/HiggsToTauTau/interface/GenJetReader.h" // GenJetReader
#include "tthAnalysis/HiggsToTauTau/interface/LHEInfoReader.h" // LHEInfoReader
#include "tthAnalysis/HiggsToTauTau/interface/EventInfoReader.h" // EventInfoReader
#include "tthAnalysis/HiggsToTauTau/interface/convert_to_ptrs.h" // convert_to_ptrs
#include "tthAnalysis/HiggsToTauTau/interface/ParticleCollectionCleaner.h" // RecoElectronCollectionCleaner, RecoMuonCollectionCleaner, RecoHadTauCollectionCleaner, RecoJetCollectionCleaner
#include "tthAnalysis/HiggsToTauTau/interface/ParticleCollectionGenMatcher.h" // RecoElectronCollectionGenMatcher, RecoMuonCollectionGenMatcher, RecoHadTauCollectionGenMatcher, RecoJetCollectionGenMatcher
#include "tthAnalysis/HiggsToTauTau/interface/RecoElectronCollectionSelectorLoose.h" // RecoElectronCollectionSelectorLoose
#include "tthAnalysis/HiggsToTauTau/interface/RecoElectronCollectionSelectorFakeable.h" // RecoElectronCollectionSelectorFakeable
#include "tthAnalysis/HiggsToTauTau/interface/RecoElectronCollectionSelectorTight.h" // RecoElectronCollectionSelectorTight
#include "tthAnalysis/HiggsToTauTau/interface/RecoMuonCollectionSelectorLoose.h" // RecoMuonCollectionSelectorLoose
#include "tthAnalysis/HiggsToTauTau/interface/RecoMuonCollectionSelectorFakeable.h" // RecoMuonCollectionSelectorFakeable
#include "tthAnalysis/HiggsToTauTau/interface/RecoMuonCollectionSelectorTight.h" // RecoMuonCollectionSelectorTight
#include "tthAnalysis/HiggsToTauTau/interface/RecoHadTauCollectionSelectorLoose.h" // RecoHadTauCollectionSelectorLoose
#include "tthAnalysis/HiggsToTauTau/interface/RecoHadTauCollectionSelectorFakeable.h" // RecoHadTauCollectionSelectorFakeable
#include "tthAnalysis/HiggsToTauTau/interface/RecoHadTauCollectionSelectorTight.h" // RecoHadTauCollectionSelectorTight
#include "tthAnalysis/HiggsToTauTau/interface/RecoJetCollectionSelector.h" // RecoJetCollectionSelector
#include "tthAnalysis/HiggsToTauTau/interface/RecoJetCollectionSelectorBtag.h" // RecoJetCollectionSelectorBtagLoose, RecoJetCollectionSelectorBtagMedium
#include "tthAnalysis/HiggsToTauTau/interface/RecoJetCollectionSelectorForward.h" // RecoJetSelectorForward
#include "tthAnalysis/HiggsToTauTau/interface/RunLumiEventSelector.h" // RunLumiEventSelector
#include "tthAnalysis/HiggsToTauTau/interface/MEtFilterSelector.h" // MEtFilterSelector
#include "tthAnalysis/HiggsToTauTau/interface/leptonTypes.h" // getLeptonType, kElectron, kMuon
#include "tthAnalysis/HiggsToTauTau/interface/analysisAuxFunctions.h" // getBTagWeight_option, getHadTau_genPdgId, isHigherPt, isMatched
#include "tthAnalysis/HiggsToTauTau/interface/leptonGenMatchingAuxFunctions.h" // getLeptonGenMatch_definitions_3lepton, getLeptonGenMatch_string, getLeptonGenMatch_int
#include "tthAnalysis/HiggsToTauTau/interface/hadTauGenMatchingAuxFunctions.h" // getHadTauGenMatch_definitions_3tau, getHadTauGenMatch_string, getHadTauGenMatch_int
#include "tthAnalysis/HiggsToTauTau/interface/fakeBackgroundAuxFunctions.h"
#include "tthAnalysis/HiggsToTauTau/interface/hltPath.h" // hltPath, create_hltPaths, hltPaths_isTriggered, hltPaths_delete
#include "tthAnalysis/HiggsToTauTau/interface/hltPathReader.h" // hltPathReader
#include "tthAnalysis/HiggsToTauTau/interface/Data_to_MC_CorrectionInterface_2016.h"
#include "tthAnalysis/HiggsToTauTau/interface/Data_to_MC_CorrectionInterface_2017.h"
#include "tthAnalysis/HiggsToTauTau/interface/Data_to_MC_CorrectionInterface_2018.h"
#include "tthAnalysis/HiggsToTauTau/interface/lutAuxFunctions.h" // loadTH2, get_sf_from_TH2
#include "tthAnalysis/HiggsToTauTau/interface/NtupleFillerBDT.h" // NtupleFillerBDT
#include "tthAnalysis/HiggsToTauTau/interface/TTreeWrapper.h" // TTreeWrapper
#include "tthAnalysis/HiggsToTauTau/interface/SyncNtupleManager.h" // SyncNtupleManager
#include "tthAnalysis/HiggsToTauTau/interface/hltFilter.h" // hltFilter()
#include "tthAnalysis/HiggsToTauTau/interface/EvtWeightManager.h" // EvtWeightManager
#include "tthAnalysis/HiggsToTauTau/interface/EvtYieldHistManager.h" // EvtYieldHistManager

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h" // ClassicSVfit
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h" // classic_svFit::MeasuredTauLepton
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h"

#include "tthAnalysis/HiggsToTauTau/interface/HadTopTagger.h" // HadTopTagger
#include "tthAnalysis/HiggsToTauTau/interface/HadTopTagger_boosted.h" // HadTopTagger_boosted
#include "tthAnalysis/HiggsToTauTau/interface/HadTopTagger_semi_boosted_AK8.h" // HadTopTagger_semi_boosted
#include "tthAnalysis/HiggsToTauTau/interface/hadTopTaggerAuxFunctions.h" // isGenMatchedJetTriplet
#include "tthAnalysis/HiggsToTauTau/interface/hadTopTaggerAuxFunctions_geral.h" // isGenMatchedJetTriplet tags
#include "tthAnalysis/HiggsToTauTau/interface/mvaAuxFunctions_Hj_and_Hjj_taggers.h" // comp_mvaOutput_Hj_tagger, comp_mvaOutput_Hjj_tagger
#include "tthAnalysis/HiggsToTauTau/interface/RecoJetCollectionSelectorHTTv2.h" // RecoJetSelectorHTTv2
#include "tthAnalysis/HiggsToTauTau/interface/RecoJetHTTv2.h"
#include "tthAnalysis/HiggsToTauTau/interface/RecoJetReaderHTTv2.h" // RecoJetReaderHTTv2
#include "tthAnalysis/HiggsToTauTau/interface/RecoJetReaderAK8.h" // RecoJetReaderAK8
#include "tthAnalysis/HiggsToTauTau/interface/JetHistManagerHTTv2.h" // JetHistManagerHTTv2
#include "tthAnalysis/HiggsToTauTau/interface/RecoJetCollectionSelectorAK8.h" // RecoJetSelectorAK8
#include "tthAnalysis/HiggsToTauTau/interface/ParticleCollectionCleanerSubJets.h" // RecoJetCollectionCleanerAK8SubJets
#include "tthAnalysis/HiggsToTauTau/interface/DYMCReweighting.h" // DYMCReweighting
#include "tthAnalysis/HiggsToTauTau/interface/DYMCNormScaleFactors.h" // DYMCNormScaleFactors

#include "tthAnalysis/HiggsToTauTau/interface/Data_to_MC_CorrectionInterface_1l_1tau_trigger.h" // Data_to_MC_CorrectionInterface_1l_1tau_trigger
#include "tthAnalysis/HiggsToTauTau/interface/Data_to_MC_CorrectionInterface_1l_2tau_trigger.h" // Data_to_MC_CorrectionInterface_1l_2tau_trigger
#include "tthAnalysis/HiggsToTauTau/interface/Data_to_MC_CorrectionInterface_0l_2tau_trigger.h" // Data_to_MC_CorrectionInterface_0l_2tau_trigger

#include <boost/math/special_functions/sign.hpp> // boost::math::sign()

#include <iostream> // std::cerr, std::fixed
#include <iomanip> // std::setprecision(), std::setw()
#include <string> // std::string
#include <vector> // std::vector<>
#include <cstdlib> // EXIT_SUCCESS, EXIT_FAILURE
#include <fstream> // std::ofstream
#include <assert.h> // assert
#include <array> // std::array<>
#include <tuple> // std::tuple<>, std::get<>(), std::make_tuple()

typedef math::PtEtaPhiMLorentzVector LV;
typedef std::vector<std::string> vstring;
typedef std::vector<double> vdouble;

enum { kFR_disabled, kFR_2lepton, kFR_4L, kFR_2tau };

//const int hadTauSelection_antiElectron = 1; // vLoose
//const int hadTauSelection_antiMuon = 1; // Loose
const int hadTauSelection_antiElectron = -1; // not applied
const int hadTauSelection_antiMuon = -1; // not applied


double comp_cosThetaS(const Particle::LorentzVector& hadTauP4_lead, const Particle::LorentzVector& hadTauP4_sublead)
{
  TLorentzVector hadTauP4tlv_lead;
  hadTauP4tlv_lead.SetPtEtaPhiM(hadTauP4_lead.pt(), hadTauP4_lead.eta(), hadTauP4_lead.phi(), hadTauP4_lead.mass());
  TLorentzVector hadTauP4tlv_sublead;
  hadTauP4tlv_sublead.SetPtEtaPhiM(hadTauP4_sublead.pt(), hadTauP4_sublead.eta(), hadTauP4_sublead.phi(), hadTauP4_sublead.mass());
  TLorentzVector hadTauBoost = hadTauP4tlv_lead;
  return std::fabs(hadTauBoost.CosTheta());
}


/**
 * @brief Produce datacard and control plots for inclusive_toNN category.
 */
int main(int argc, char* argv[])
{
//--- throw an exception in case ROOT encounters an error
  gErrorAbortLevel = kError;

//--- parse command-line arguments
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " [parameters.py]" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "<analyze_inclusive_toNN>:" << std::endl;

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("analyze_inclusive_toNN");

//--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") )
    throw cms::Exception("analyze_inclusive_toNN")
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::ParameterSet cfg_analyze = cfg.getParameter<edm::ParameterSet>("analyze_inclusive_toNN");

  std::string treeName = cfg_analyze.getParameter<std::string>("treeName");

  std::string process_string = cfg_analyze.getParameter<std::string>("process");
  bool isSignal = ( process_string == "signal" ) ? true : false;

  std::string era_string = cfg_analyze.getParameter<std::string>("era");
  const int era = get_era(era_string);

  vstring triggerNames_1e = cfg_analyze.getParameter<vstring>("triggers_1e");
  std::vector<hltPath*> triggers_1e = create_hltPaths(triggerNames_1e);
  bool use_triggers_1e = cfg_analyze.getParameter<bool>("use_triggers_1e");
  vstring triggerNames_2e = cfg_analyze.getParameter<vstring>("triggers_2e");
  std::vector<hltPath*> triggers_2e = create_hltPaths(triggerNames_2e);
  bool use_triggers_2e = cfg_analyze.getParameter<bool>("use_triggers_2e");
  vstring triggerNames_1mu = cfg_analyze.getParameter<vstring>("triggers_1mu");
  std::vector<hltPath*> triggers_1mu = create_hltPaths(triggerNames_1mu);
  bool use_triggers_1mu = cfg_analyze.getParameter<bool>("use_triggers_1mu");
  vstring triggerNames_2mu = cfg_analyze.getParameter<vstring>("triggers_2mu");
  std::vector<hltPath*> triggers_2mu = create_hltPaths(triggerNames_2mu);
  bool use_triggers_2mu = cfg_analyze.getParameter<bool>("use_triggers_2mu");
  vstring triggerNames_1e1mu = cfg_analyze.getParameter<vstring>("triggers_1e1mu");
  std::vector<hltPath*> triggers_1e1mu = create_hltPaths(triggerNames_1e1mu);
  bool use_triggers_1e1mu = cfg_analyze.getParameter<bool>("use_triggers_1e1mu");
  vstring triggerNames_2tau = cfg_analyze.getParameter<vstring>("triggers_2tau");
  std::vector<hltPath*> triggers_2tau = create_hltPaths(triggerNames_2tau);
  bool use_triggers_2tau = cfg_analyze.getParameter<bool>("use_triggers_2tau");
  // triple lepton triggers
  vstring triggerNames_3e = cfg_analyze.getParameter<vstring>("triggers_3e");
  std::vector<hltPath*> triggers_3e = create_hltPaths(triggerNames_3e, "triggers_3e");
  bool use_triggers_3e = cfg_analyze.getParameter<bool>("use_triggers_3e");
  vstring triggerNames_2e1mu = cfg_analyze.getParameter<vstring>("triggers_2e1mu");
  std::vector<hltPath*> triggers_2e1mu = create_hltPaths(triggerNames_2e1mu, "triggers_2e1mu");
  bool use_triggers_2e1mu = cfg_analyze.getParameter<bool>("use_triggers_2e1mu");
  vstring triggerNames_1e2mu = cfg_analyze.getParameter<vstring>("triggers_1e2mu");
  std::vector<hltPath*> triggers_1e2mu = create_hltPaths(triggerNames_1e2mu, "triggers_1e2mu");
  bool use_triggers_1e2mu = cfg_analyze.getParameter<bool>("use_triggers_1e2mu");
  vstring triggerNames_3mu = cfg_analyze.getParameter<vstring>("triggers_3mu");
  std::vector<hltPath*> triggers_3mu = create_hltPaths(triggerNames_3mu, "triggers_3mu");
  bool use_triggers_3mu = cfg_analyze.getParameter<bool>("use_triggers_3mu");
  vstring triggerNames_1e1tau = cfg_analyze.getParameter<vstring>("triggers_1e1tau");
  std::vector<hltPath*> triggers_1e1tau = create_hltPaths(triggerNames_1e1tau, "triggers_1e1tau");
  bool use_triggers_1e1tau = cfg_analyze.getParameter<bool>("use_triggers_1e1tau");
  vstring triggerNames_1mu1tau = cfg_analyze.getParameter<vstring>("triggers_1mu1tau");
  std::vector<hltPath*> triggers_1mu1tau = create_hltPaths(triggerNames_1mu1tau, "triggers_1mu1tau");
  bool use_triggers_1mu1tau = cfg_analyze.getParameter<bool>("use_triggers_1mu1tau");

  const std::string electronSelection_string = cfg_analyze.getParameter<std::string>("electronSelection");
  const std::string muonSelection_string     = cfg_analyze.getParameter<std::string>("muonSelection");
  std::cout << "electronSelection_string = " << electronSelection_string << "\n"
               "muonSelection_string     = " << muonSelection_string     << '\n'
  ;
  const int electronSelection = get_selection(electronSelection_string);
  const int muonSelection     = get_selection(muonSelection_string);
  double lep_mva_cut = cfg_analyze.getParameter<double>("lep_mva_cut"); // CV: used for tight lepton selection only

  bool apply_leptonGenMatching = cfg_analyze.getParameter<bool>("apply_leptonGenMatching");
  std::vector<leptonGenMatchEntry> leptonGenMatch_definitions = getLeptonGenMatch_definitions_2lepton(apply_leptonGenMatching);
  std::cout << "leptonGenMatch_definitions:" << std::endl;
  std::cout << leptonGenMatch_definitions;

  TString hadTauSelection_string = cfg_analyze.getParameter<std::string>("hadTauSelection").data();
  TObjArray* hadTauSelection_parts = hadTauSelection_string.Tokenize("|");
  assert(hadTauSelection_parts->GetEntries() >= 1);
  const std::string hadTauSelection_part1 = (dynamic_cast<TObjString*>(hadTauSelection_parts->At(0)))->GetString().Data();
  const int hadTauSelection = get_selection(hadTauSelection_part1);
  std::string hadTauSelection_part2 = ( hadTauSelection_parts->GetEntries() == 2 ) ? (dynamic_cast<TObjString*>(hadTauSelection_parts->At(1)))->GetString().Data() : "";
  delete hadTauSelection_parts;

  bool apply_hadTauGenMatching = cfg_analyze.getParameter<bool>("apply_hadTauGenMatching");
  std::vector<hadTauGenMatchEntry> hadTauGenMatch_definitions = getHadTauGenMatch_definitions_2tau(apply_hadTauGenMatching);
  std::cout << "hadTauGenMatch_definitions:" << std::endl;
  std::cout << hadTauGenMatch_definitions;

  bool isMC = cfg_analyze.getParameter<bool>("isMC");
  bool isMC_tH = ( process_string == "tHq" || process_string == "tHW" ) ? true : false;
  bool hasLHE = cfg_analyze.getParameter<bool>("hasLHE");
  std::string central_or_shift = cfg_analyze.getParameter<std::string>("central_or_shift");
  double lumiScale = ( process_string != "data_obs" ) ? cfg_analyze.getParameter<double>("lumiScale") : 1.;
  bool apply_genWeight = cfg_analyze.getParameter<bool>("apply_genWeight");
  edm::ParameterSet cfgMEtFilter = cfg_analyze.getParameter<edm::ParameterSet>("cfgMEtFilter");
  MEtFilterSelector metFilterSelector(cfgMEtFilter, isMC);
  bool apply_DYMCReweighting = cfg_analyze.getParameter<bool>("apply_DYMCReweighting");
  bool apply_DYMCNormScaleFactors = cfg_analyze.getParameter<bool>("apply_DYMCNormScaleFactors");
  //bool apply_hadTauFakeRateSF = cfg_analyze.getParameter<bool>("apply_hadTauFakeRateSF");
  const bool useNonNominal = cfg_analyze.getParameter<bool>("useNonNominal");
  const bool useNonNominal_jetmet = useNonNominal || ! isMC;

  const edm::ParameterSet syncNtuple_cfg = cfg_analyze.getParameter<edm::ParameterSet>("syncNtuple");
  const std::string syncNtuple_tree = syncNtuple_cfg.getParameter<std::string>("tree");
  const std::string syncNtuple_output = syncNtuple_cfg.getParameter<std::string>("output");

  const edm::ParameterSet additionalEvtWeight = cfg_analyze.getParameter<edm::ParameterSet>("evtWeight");
  const bool applyAdditionalEvtWeight = additionalEvtWeight.getParameter<bool>("apply");
  EvtWeightManager * eventWeightManager = nullptr;
  if(applyAdditionalEvtWeight)
  {
    eventWeightManager = new EvtWeightManager(additionalEvtWeight);
  }

  bool isDEBUG = cfg_analyze.getParameter<bool>("isDEBUG");
  if ( isDEBUG ) std::cout << "Warning: DEBUG mode enabled -> trigger selection will not be applied for data !!" << std::endl;

  checkOptionValidity(central_or_shift, isMC);
  const int jetToLeptonFakeRate_option = getJetToLeptonFR_option(central_or_shift);
  const int hadTauPt_option            = getHadTauPt_option     (central_or_shift);
  const int jetToTauFakeRate_option    = getJetToTauFR_option   (central_or_shift);
  const int lheScale_option            = getLHEscale_option     (central_or_shift);
  const int jetBtagSF_option           = getBTagWeight_option   (central_or_shift);
  const PUsys puSys_option             = getPUsys_option        (central_or_shift);
  const int dyMCReweighting_option  = getDYMCReweighting_option(central_or_shift);
  const int dyMCNormScaleFactors_option  = getDYMCNormScaleFactors_option(central_or_shift);

  const int met_option   = useNonNominal_jetmet ? kMEt_central_nonNominal : getMET_option(central_or_shift, isMC);
  const int jetPt_option = useNonNominal_jetmet ? kMEt_central_nonNominal : getJet_option(central_or_shift, isMC);

  std::cout
    << "central_or_shift = "               << central_or_shift           << "\n"
       " -> jetToLeptonFakeRate_option = " << jetToLeptonFakeRate_option << "\n"
       " -> hadTauPt_option            = " << hadTauPt_option            << "\n"
       " -> jetToTauFakeRate_option    = " << jetToTauFakeRate_option    << "\n"
       " -> lheScale_option            = " << lheScale_option            << "\n"
       " -> jetBtagSF_option           = " << jetBtagSF_option           << "\n"
       " -> met_option                 = " << met_option                 << "\n"
       " -> jetPt_option               = " << jetPt_option               << '\n'
  ;

  DYMCReweighting * dyReweighting = nullptr;
  if(apply_DYMCReweighting)
  {
    dyReweighting = new DYMCReweighting(era, dyMCReweighting_option);
  }
  DYMCNormScaleFactors* dyNormScaleFactors = nullptr;
  if ( apply_DYMCNormScaleFactors ) {
    dyNormScaleFactors = new DYMCNormScaleFactors(era, dyMCNormScaleFactors_option);
  }

  edm::ParameterSet cfg_dataToMCcorrectionInterface;
  cfg_dataToMCcorrectionInterface.addParameter<std::string>("era", era_string);
  cfg_dataToMCcorrectionInterface.addParameter<std::string>("hadTauSelection", hadTauSelection_part2);
  cfg_dataToMCcorrectionInterface.addParameter<int>("hadTauSelection_antiElectron", hadTauSelection_antiElectron);
  cfg_dataToMCcorrectionInterface.addParameter<int>("hadTauSelection_antiMuon", hadTauSelection_antiMuon);
  cfg_dataToMCcorrectionInterface.addParameter<std::string>("central_or_shift", central_or_shift);
  cfg_dataToMCcorrectionInterface.addParameter<bool>("isDEBUG", isDEBUG);
  Data_to_MC_CorrectionInterface_Base * dataToMCcorrectionInterface = nullptr;
  switch(era)
  {
    case kEra_2016: dataToMCcorrectionInterface = new Data_to_MC_CorrectionInterface_2016(cfg_dataToMCcorrectionInterface); break;
    case kEra_2017: dataToMCcorrectionInterface = new Data_to_MC_CorrectionInterface_2017(cfg_dataToMCcorrectionInterface); break;
    case kEra_2018: dataToMCcorrectionInterface = new Data_to_MC_CorrectionInterface_2018(cfg_dataToMCcorrectionInterface); break;
    default: throw cmsException("analyze_inclusive_toNN", __LINE__) << "Invalid era = " << era;
  }
  Data_to_MC_CorrectionInterface_1l_1tau_trigger* dataToMCcorrectionInterface_1l_1tau_trigger = new Data_to_MC_CorrectionInterface_1l_1tau_trigger(cfg_dataToMCcorrectionInterface);
  Data_to_MC_CorrectionInterface_0l_2tau_trigger* dataToMCcorrectionInterface_0l_2tau_trigger = new Data_to_MC_CorrectionInterface_0l_2tau_trigger(cfg_dataToMCcorrectionInterface);
  Data_to_MC_CorrectionInterface_1l_2tau_trigger* dataToMCcorrectionInterface_1l_2tau_trigger = new Data_to_MC_CorrectionInterface_1l_2tau_trigger(cfg_dataToMCcorrectionInterface);

  //--- initialize hadronic top tagger BDT
  HadTopTagger* hadTopTagger = new HadTopTagger();
  HadTopTagger_boosted* hadTopTagger_boosted = new HadTopTagger_boosted();
  HadTopTagger_semi_boosted_AK8* hadTopTagger_semi_boosted_fromAK8 = new HadTopTagger_semi_boosted_AK8();

  // Hj-tagger
  std::string mvaFileName_Hj_tagger = "tthAnalysis/HiggsToTauTau/data/Hj_deepcsv_BDTG_2017.weights.xml";
  std::vector<std::string> mvaInputVariables_Hj_tagger = {
    "Jet25_lepdrmin", "max(Jet25_bDiscriminator,0.)",
    "max(Jet25_qg,0.)", "Jet25_lepdrmax", "Jet25_pt" };
  TMVAInterface mva_Hj_tagger(mvaFileName_Hj_tagger, mvaInputVariables_Hj_tagger);

  std::cout << "start fake rates = " << std::endl;
  /*
  std::string applyFakeRateWeights_string = cfg_analyze.getParameter<std::string>("applyFakeRateWeights");
  int applyFakeRateWeights = -1;
  if      ( applyFakeRateWeights_string == "disabled" ) applyFakeRateWeights = kFR_disabled;
  else if ( applyFakeRateWeights_string == "2lepton"  ) applyFakeRateWeights = kFR_2lepton;
  else if ( applyFakeRateWeights_string == "4L"       ) applyFakeRateWeights = kFR_4L;
  else if ( applyFakeRateWeights_string == "2tau"     ) applyFakeRateWeights = kFR_2tau;
  else throw cms::Exception("analyze_inclusive_toNN")
    << "Invalid Configuration parameter 'applyFakeRateWeights' = " << applyFakeRateWeights_string << " !!\n";
  */

  LeptonFakeRateInterface* leptonFakeRateInterface = 0;
  edm::ParameterSet cfg_leptonFakeRateWeight = cfg_analyze.getParameter<edm::ParameterSet>("leptonFakeRateWeight");
  leptonFakeRateInterface = new LeptonFakeRateInterface(cfg_leptonFakeRateWeight, jetToLeptonFakeRate_option);

  JetToTauFakeRateInterface* jetToTauFakeRateInterface = 0;
  edm::ParameterSet cfg_hadTauFakeRateWeight = cfg_analyze.getParameter<edm::ParameterSet>("hadTauFakeRateWeight");
  cfg_hadTauFakeRateWeight.addParameter<std::string>("hadTauSelection", hadTauSelection_part2);
  jetToTauFakeRateInterface = new JetToTauFakeRateInterface(cfg_hadTauFakeRateWeight, jetToTauFakeRate_option);

  bool fillGenEvtHistograms = cfg_analyze.getParameter<bool>("fillGenEvtHistograms");

  std::string branchName_electrons = cfg_analyze.getParameter<std::string>("branchName_electrons");
  std::string branchName_muons = cfg_analyze.getParameter<std::string>("branchName_muons");
  std::string branchName_hadTaus = cfg_analyze.getParameter<std::string>("branchName_hadTaus");
  std::string branchName_jets = cfg_analyze.getParameter<std::string>("branchName_jets");
  std::string branchName_met = cfg_analyze.getParameter<std::string>("branchName_met");
  std::string branchName_jetsHTTv2 = cfg_analyze.getParameter<std::string>("branchName_jetsHTTv2");
  std::string branchName_subjetsHTTv2 = cfg_analyze.getParameter<std::string>("branchName_subjetsHTTv2");
  std::string branchName_jetsAK8 = cfg_analyze.getParameter<std::string>("branchName_jetsAK8");
  std::string branchName_subjetsAK8 = cfg_analyze.getParameter<std::string>("branchName_subjetsAK8");

  std::string branchName_genLeptons = cfg_analyze.getParameter<std::string>("branchName_genLeptons");
  std::string branchName_genHadTaus = cfg_analyze.getParameter<std::string>("branchName_genHadTaus");
  std::string branchName_genPhotons = cfg_analyze.getParameter<std::string>("branchName_genPhotons");
  std::string branchName_genJets = cfg_analyze.getParameter<std::string>("branchName_genJets");
  std::string branchName_genTauLeptons = cfg_analyze.getParameter<std::string>("branchName_genTauLeptons");
  std::string branchName_genTopQuarks = cfg_analyze.getParameter<std::string>("branchName_genTopQuarks");
  std::string branchName_genBJets = cfg_analyze.getParameter<std::string>("branchName_genBJets");
  std::string branchName_genWBosons = cfg_analyze.getParameter<std::string>("branchName_genWBosons");
  std::string branchName_genWJets = cfg_analyze.getParameter<std::string>("branchName_genWJets");
  std::string branchName_genQuarkFromTop = cfg_analyze.getParameter<std::string>("branchName_genQuarkFromTop");

  bool redoGenMatching = cfg_analyze.getParameter<bool>("redoGenMatching");

  std::string selEventsFileName_input = cfg_analyze.getParameter<std::string>("selEventsFileName_input");
  std::cout << "selEventsFileName_input = " << selEventsFileName_input << std::endl;
  RunLumiEventSelector* run_lumi_eventSelector = 0;
  if ( selEventsFileName_input != "" ) {
    edm::ParameterSet cfg_runLumiEventSelector;
    cfg_runLumiEventSelector.addParameter<std::string>("inputFileName", selEventsFileName_input);
    cfg_runLumiEventSelector.addParameter<std::string>("separator", ":");
    run_lumi_eventSelector = new RunLumiEventSelector(cfg_runLumiEventSelector);
  }

  std::string selEventsFileName_output = cfg_analyze.getParameter<std::string>("selEventsFileName_output");
  std::cout << "selEventsFileName_output = " << selEventsFileName_output << std::endl;

  fwlite::InputSource inputFiles(cfg);
  int maxEvents = inputFiles.maxEvents();
  std::cout << " maxEvents = " << maxEvents << std::endl;
  unsigned reportEvery = inputFiles.reportAfter();

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());

  TTreeWrapper* inputTree = new TTreeWrapper(treeName.data(), inputFiles.files(), maxEvents);
  std::cout << "Loaded " << inputTree->getFileCount() << " file(s)." << std::endl;

//--- declare event-level variables
  EventInfo eventInfo(isSignal, isMC, isMC_tH);
  EventInfoReader eventInfoReader(&eventInfo, puSys_option);
  inputTree->registerReader(&eventInfoReader);

  hltPathReader hltPathReader_instance({
    triggers_1e, triggers_2e, triggers_1mu, triggers_2mu, triggers_1e1mu,
    triggers_2tau,
    triggers_3e, triggers_2e1mu, triggers_1e2mu, triggers_3mu,
    triggers_1e1tau, triggers_1mu1tau
  });
  inputTree -> registerReader(&hltPathReader_instance);

  if(eventWeightManager)
  {
    inputTree->registerReader(eventWeightManager);
  }

//--- declare particle collections
  const bool readGenObjects = isMC && !redoGenMatching;
  RecoMuonReader* muonReader = new RecoMuonReader(era, branchName_muons, readGenObjects);
  inputTree->registerReader(muonReader);
  RecoMuonCollectionGenMatcher muonGenMatcher;
  RecoMuonCollectionSelectorLoose preselMuonSelector(era, -1, isDEBUG);
  RecoMuonCollectionSelectorFakeable fakeableMuonSelector(era, -1, isDEBUG);
  RecoMuonCollectionSelectorTight tightMuonSelector(era, -1, isDEBUG);
  tightMuonSelector.getSelector().set_min_mvaTTH(lep_mva_cut);

  RecoElectronReader* electronReader = new RecoElectronReader(era, branchName_electrons, readGenObjects);
  electronReader->readUncorrected(useNonNominal);
  inputTree->registerReader(electronReader);
  RecoElectronCollectionGenMatcher electronGenMatcher;
  RecoElectronCollectionCleaner electronCleaner(0.3, isDEBUG);
  RecoElectronCollectionSelectorLoose preselElectronSelector(era, -1, isDEBUG);
  RecoElectronCollectionSelectorFakeable fakeableElectronSelector(era, -1, isDEBUG);
  RecoElectronCollectionSelectorTight tightElectronSelector(era, -1, isDEBUG);
  tightElectronSelector.getSelector().set_min_mvaTTH(lep_mva_cut);

  RecoHadTauReader* hadTauReader = new RecoHadTauReader(era, branchName_hadTaus, readGenObjects);
  hadTauReader->setHadTauPt_central_or_shift(hadTauPt_option);
  inputTree->registerReader(hadTauReader);
  RecoHadTauCollectionGenMatcher hadTauGenMatcher;
  RecoHadTauCollectionCleaner hadTauCleaner(0.3, isDEBUG);
  RecoHadTauCollectionSelectorLoose preselHadTauSelector(era, -1, isDEBUG);
  preselHadTauSelector.set_if_looser(hadTauSelection_part2);
  preselHadTauSelector.set_min_antiElectron(hadTauSelection_antiElectron);
  preselHadTauSelector.set_min_antiMuon(hadTauSelection_antiMuon);
  RecoHadTauCollectionSelectorFakeable fakeableHadTauSelector(era, -1, isDEBUG);
  fakeableHadTauSelector.set_if_looser(hadTauSelection_part2);
  fakeableHadTauSelector.set_min_antiElectron(hadTauSelection_antiElectron);
  fakeableHadTauSelector.set_min_antiMuon(hadTauSelection_antiMuon);
  RecoHadTauCollectionSelectorTight tightHadTauSelector(era, -1, isDEBUG);
  tightHadTauSelector.set(hadTauSelection_part2);
  tightHadTauSelector.set_min_antiElectron(hadTauSelection_antiElectron);
  tightHadTauSelector.set_min_antiMuon(hadTauSelection_antiMuon);

  RecoJetReader* jetReader = new RecoJetReader(era, isMC, branchName_jets, readGenObjects);
  jetReader->setPtMass_central_or_shift(jetPt_option);
  jetReader->setBranchName_BtagWeight(jetBtagSF_option);
  inputTree->registerReader(jetReader);
  RecoJetCollectionGenMatcher jetGenMatcher;
  RecoJetCollectionCleaner jetCleaner(0.4, isDEBUG);
  RecoJetCollectionCleaner jetCleaner_large8(0.8, isDEBUG);
  RecoJetCollectionSelector jetSelector(era, -1, isDEBUG);
  RecoJetCollectionSelectorForward jetSelectorForward(era, -1, isDEBUG);
  RecoJetCollectionSelectorBtagLoose jetSelectorBtagLoose(era, -1, isDEBUG);
  RecoJetCollectionSelectorBtagMedium jetSelectorBtagMedium(era, -1, isDEBUG);

  RecoJetReaderHTTv2* jetReaderHTTv2 = new RecoJetReaderHTTv2(era, branchName_jetsHTTv2, branchName_subjetsHTTv2);
  inputTree -> registerReader(jetReaderHTTv2);
  RecoJetCollectionSelectorHTTv2 jetSelectorHTTv2(era);
  RecoJetCollectionCleanerHTTv2 jetCleanerHTTv2(1.5, isDEBUG); //to clean against leptons and hadronic taus
  RecoJetCollectionCleanerHTTv2SubJets jetCleanerHTTv2SubJets(0.4, isDEBUG); //to clean against leptons and hadronic taus

  RecoJetReaderAK8* jetReaderAK8 = new RecoJetReaderAK8(era, branchName_jetsAK8, branchName_subjetsAK8);
  inputTree -> registerReader(jetReaderAK8);
  RecoJetCollectionSelectorAK8 jetSelectorAK8(era);
  RecoJetCollectionCleanerAK8 jetCleanerAK8(0.8, isDEBUG); //to clean against leptons and hadronic taus
  RecoJetCollectionCleanerAK8SubJets jetCleanerAK8SubJets(0.4, isDEBUG); //to clean against leptons and hadronic taus

//--- declare missing transverse energy
  RecoMEtReader* metReader = new RecoMEtReader(era, isMC, branchName_met);
  metReader->setMEt_central_or_shift(met_option);
  inputTree->registerReader(metReader);

  MEtFilter metFilters;
  MEtFilterReader* metFilterReader = new MEtFilterReader(&metFilters, era);
  inputTree -> registerReader(metFilterReader);

  GenParticleReader* genTauLeptonReader = nullptr;
  if ( isMC && (apply_DYMCReweighting || apply_DYMCNormScaleFactors)) {
    genTauLeptonReader = new GenParticleReader(branchName_genTauLeptons);
    inputTree->registerReader(genTauLeptonReader);
  }
  GenParticleReader* genTopQuarkReader = new GenParticleReader(branchName_genTopQuarks);
  GenParticleReader* genBJetReader = new GenParticleReader(branchName_genBJets);
  GenParticleReader* genWBosonReader = new GenParticleReader(branchName_genWBosons);
  GenParticleReader* genWJetReader = new GenParticleReader(branchName_genWJets);
  GenParticleReader* genQuarkFromTopReader = new GenParticleReader(branchName_genQuarkFromTop);

  if ( isMC ) {
    inputTree -> registerReader(genTopQuarkReader);
    inputTree -> registerReader(genBJetReader);
    inputTree -> registerReader(genWBosonReader);
    inputTree -> registerReader(genWJetReader);
    inputTree -> registerReader(genQuarkFromTopReader);
  }

//--- declare generator level information
  GenLeptonReader* genLeptonReader = 0;
  GenHadTauReader* genHadTauReader = 0;
  GenPhotonReader* genPhotonReader = 0;
  GenJetReader* genJetReader = 0;
  LHEInfoReader* lheInfoReader = 0;
  if ( isMC ) {
    if(! readGenObjects)
    {
      if ( branchName_genLeptons != "" ) {
        genLeptonReader = new GenLeptonReader(branchName_genLeptons);
        inputTree->registerReader(genLeptonReader);
      }
      if ( branchName_genHadTaus != "" ) {
        genHadTauReader = new GenHadTauReader(branchName_genHadTaus);
        inputTree->registerReader(genHadTauReader);
      }
      if ( branchName_genPhotons != "" ) {
        genPhotonReader = new GenPhotonReader(branchName_genPhotons);
        inputTree -> registerReader(genPhotonReader);
      }
      if ( branchName_genJets != "" ) {
        genJetReader = new GenJetReader(branchName_genJets);
        inputTree->registerReader(genJetReader);
      }
    }
    lheInfoReader = new LHEInfoReader(hasLHE);
    inputTree->registerReader(lheInfoReader);
  }

//--- open output file containing run:lumi:event numbers of events passing final event selection criteria
  std::ostream* selEventsFile = ( selEventsFileName_output != "" ) ? new std::ofstream(selEventsFileName_output.data(), std::ios::out) : 0;
  std::cout << "selEventsFileName_output = " << selEventsFileName_output << std::endl;

  bool selectBDT = ( cfg_analyze.exists("selectBDT") ) ? cfg_analyze.getParameter<bool>("selectBDT") : false;

  NtupleFillerBDT<float, int>* bdt_filler = nullptr;
  typedef std::remove_pointer<decltype(bdt_filler)>::type::float_type float_type;
  typedef std::remove_pointer<decltype(bdt_filler)>::type::int_type   int_type;

  //--- initialize BDTs used to discriminate ttH vs. ttV and ttbar

  if ( selectBDT ) {
    bdt_filler = new std::remove_pointer<decltype(bdt_filler)>::type(
      makeHistManager_cfg(process_string, "inclusive/sel/evtntuple", central_or_shift)
    );
    bdt_filler -> register_variable<float_type>(
      "avg_dr_jet", "ptmiss",  "htmiss", "evtWeight",
      "mTauTauVis", "mTauTau_SVFit", "cosThetaS_hadTau",
      "mTauTauVis2", "mTauTau2_SVFit", "cosThetaS_hadTau2",
      "mbb_loose", "mbb_medium", "max_dr_jet", "met_LD",
      "jet1_pt", "jet1_eta", "jet1_phi", "jet1_mass",
      "jet2_pt", "jet2_eta", "jet2_phi", "jet2_mass",
      "jet3_pt", "jet3_eta", "jet3_phi", "jet3_mass",
      "jet4_pt", "jet4_eta", "jet4_phi", "jet4_mass",
      //
      "lep1_conept", "lep1_mT", "lep1_min_dr_jet", "lep1_pt", "lep1_eta", "lep1_phi", "lep1_mass",
      "lep2_conept", "lep2_mT", "lep2_min_dr_jet", "lep2_pt", "lep2_eta", "lep2_phi", "lep2_mass",
      "lep3_conept", "lep3_mT", "lep3_min_dr_jet", "lep3_pt", "lep3_eta", "lep3_phi", "lep3_mass",
      "lep4_conept", "lep4_mT", "lep4_min_dr_jet", "lep4_pt", "lep4_eta", "lep4_phi", "lep4_mass",
      //
      "tau1_min_dr_jet", "tau1_pt",  "tau1_eta", "tau1_phi", "tau1_mass",
      "tau2_min_dr_jet", "tau2_pt",  "tau2_eta", "tau2_phi", "tau2_mass",
      "tau3_min_dr_jet", "tau3_pt",  "tau3_eta", "tau3_phi", "tau3_mass",
      "tau4_min_dr_jet", "tau4_pt",  "tau4_eta", "tau4_phi", "tau4_mass",
      //
      "DijetForward_mass",
      "jetForward1_pt", "jetForward1_eta", "jetForward1_phi", "jetForward1_mass",
      "jetForward2_pt", "jetForward2_eta", "jetForward2_phi", "jetForward2_mass",
      //
      "res-HTT_CSVsort4rd", "res-HTT_CSVsort4rd_2",
      "HadTop_pt_CSVsort4rd", "HadTop_pt_CSVsort4rd_2", "HadTop_eta_CSVsort4rd",
      "genTopPt_CSVsort4rd", "genTopPt_CSVsort4rd_2",
      "HTTv2_lead_pt", "minDR_HTTv2_lep", "minDR_HTTv2subjets_lep",
      "HTT_boosted", "genTopPt_boosted", "HadTop_pt_boosted",
      "HTT_semi_boosted_fromAK8",
      "genTopPt_semi_boosted_fromAK8", "HadTop_pt_semi_boosted_fromAK8",
      "minDR_AK8_lep", "minDR_AK8subjets_lep",
      "weight_fakeRate", "weight_data_to_MC_correction", "weight_data_to_MC_correction_hadTau",
      "mvaOutput_Hj_tagger"
    );
    bdt_filler -> register_variable<int_type>(
      "nJet", "nJetForward", "nBJetLoose", "nBJetMedium", "nElectron", "nLepton", "nTau",
      "lep_match_jet", "lep_match_lep", "tau_match_jet", "tau_match_tau",
      // '1e', '1mu', '2e', '2mu', '1e1mu', '2tau', '3e', '3mu', '1e2mu', '2e1mu', '1e1tau', '1mu1tau'
      "selTrigger_1e",  "selTrigger_2e", "selTrigger_3e",
      "selTrigger_1mu", "selTrigger_2mu", "selTrigger_3mu",
      "selTrigger_1e1mu", "selTrigger_1e2mu","selTrigger_2e1mu",
      "selTrigger_1e1tau", "selTrigger_1mu1tau", "selTrigger_2tau",
      //
      "isCharge_hadTau_OS", "isCharge_Lepton_OS",
      "sum_lep_charge", "sum_tau_charge",
      //
      "pass_ttH",
      "pass_2lss_0tau", "pass_3l_0tau", "pass_4l_0tau",
      "pass_ttZctrl", "pass_ttWctrl", "pass_WZctrl", "pass_ZZctrl",
      "pass_1l_2tau",  "pass_2lss_1tau", "pass_2los_1tau", "pass_3l_1tau", "pass_2l_2tau",
      "pass_1l_1tau", "pass_1l_3tau", "pass_0l_4tau", "pass_2los_0tau", "pass_0l_3tau", "pass_0l_2tau",
      //
      "hadtruth", "hadtruth_2" , "hadtruth_boosted", "hadtruth_semi_boosted_fromAK8",
      "nHTTv2", "N_jetAK8",
      "bWj1Wj2_isGenMatched_CSVsort4rd", "bWj1Wj2_isGenMatched_CSVsort4rd_2",
      "bWj1Wj2_isGenMatched_boosted", "bWj1Wj2_isGenMatched_semi_boosted_fromAK8",
      "resolved_and_semi_AK8", "boosted_and_semi_AK8", "resolved_and_boosted",
      "failsZbosonMassVeto", "failsLowMassVeto", "failsHtoZZVeto"
    );
    bdt_filler -> bookTree(fs);
  }

  int analyzedEntries = 0;
  int selectedEntries = 0;
  double selectedEntries_weighted = 0.;
  TH1* histogram_analyzedEntries = fs.make<TH1D>("analyzedEntries", "analyzedEntries", 1, -0.5, +0.5);
  TH1* histogram_selectedEntries = fs.make<TH1D>("selectedEntries", "selectedEntries", 1, -0.5, +0.5);
  while ( inputTree->hasNextEvent() && (! run_lumi_eventSelector || (run_lumi_eventSelector && ! run_lumi_eventSelector -> areWeDone())) ) {
    if(inputTree -> canReport(reportEvery))
    {
      std::cout << "processing Entry " << inputTree -> getCurrentMaxEventIdx()
                << " or " << inputTree -> getCurrentEventIdx() << " entry in #"
                << (inputTree -> getProcessedFileCount() - 1)
                << " (" << eventInfo
                << ") file (" << selectedEntries << " Entries selected)\n";
    }
    ++analyzedEntries;
    //if ( analyzedEntries > 100000 ) break;
    histogram_analyzedEntries->Fill(0.);

    if ( isDEBUG ) {
      std::cout << "event #" << inputTree -> getCurrentMaxEventIdx() << ' ' << eventInfo << '\n';
    }

    if ( run_lumi_eventSelector && !(*run_lumi_eventSelector)(eventInfo) ) continue;

    if ( run_lumi_eventSelector ) {
      std::cout << "processing Entry #" << inputTree->getCumulativeMaxEventCount() << ": " << eventInfo << std::endl;
      if ( inputTree -> isOpen() ) {
	std::cout << "input File = " << inputTree->getCurrentFileName() << std::endl;
      }
    }

    bool pass_2lss_0tau = true;
    bool pass_2los_0tau = true;
    bool pass_3l_0tau = true;
    bool pass_4l_0tau = true;
    bool pass_2lss_1tau = true;
    bool pass_2los_1tau = true;
    bool pass_3l_1tau = true;
    bool pass_2l_2tau = true;
    bool pass_1l_1tau = true;
    bool pass_1l_2tau = true;
    bool pass_1l_3tau = true;
    bool pass_0l_4tau = true;
    bool pass_0l_3tau = true;
    bool pass_0l_2tau = true;
    bool pass_ttZctrl = true;
    bool pass_ttWctrl = true;
    bool pass_WZctrl = true;
    bool pass_ZZctrl = true;
    bool pass_ttH = true;

//--- build collections of generator level particles (before any cuts are applied, to check distributions in unbiased event samples)
    std::vector<GenLepton> genLeptons;
    std::vector<GenLepton> genElectrons;
    std::vector<GenLepton> genMuons;
    std::vector<GenHadTau> genHadTaus;
    std::vector<GenPhoton> genPhotons;
    std::vector<GenJet> genJets;
    if ( isMC && fillGenEvtHistograms ) {
      if ( genLeptonReader ) {
	genLeptons = genLeptonReader->read();
	for ( std::vector<GenLepton>::const_iterator genLepton = genLeptons.begin();
	      genLepton != genLeptons.end(); ++genLepton ) {
	  int abs_pdgId = std::abs(genLepton->pdgId());
	  if      ( abs_pdgId == 11 ) genElectrons.push_back(*genLepton);
	  else if ( abs_pdgId == 13 ) genMuons.push_back(*genLepton);
	}
      }
      if ( genHadTauReader ) {
	genHadTaus = genHadTauReader->read();
      }
      if ( genPhotonReader ) {
        genPhotons = genPhotonReader->read();
      }
      if ( genJetReader ) {
	genJets = genJetReader->read();
      }
      if ( isDEBUG ) {
        printCollection("genLeptons", genLeptons);
        printCollection("genHadTaus", genHadTaus);
        printCollection("genPhotons", genPhotons);
        printCollection("genJets", genJets);
      }
    }

    std::vector<GenParticle> genTauLeptons;
    if(isMC && (apply_DYMCReweighting || apply_DYMCNormScaleFactors))
    {
      genTauLeptons = genTauLeptonReader->read();
    }

    double evtWeight_inclusive = 1.;
    //std::cout<<" read LHE weights \n" <<
    //"apply_genWeight " << apply_genWeight << "\n" <<
    //"isMC_tH " << apply_genWeight << "\n" <<
    //"apply_DYMCReweighting " << apply_DYMCReweighting << "\n" <<
    //"eventWeightManager " << eventWeightManager << "\n" ;
    if(isMC)
    {
      if(apply_genWeight)    evtWeight_inclusive *= boost::math::sign(eventInfo.genWeight);
      if(isMC_tH)            evtWeight_inclusive *= eventInfo.genWeight_tH;
      if(apply_DYMCReweighting) evtWeight_inclusive *= dyReweighting->getWeight(genTauLeptons);
      if(eventWeightManager) evtWeight_inclusive *= eventWeightManager->getWeight();
      lheInfoReader->read();
      evtWeight_inclusive *= lheInfoReader->getWeight_scale(lheScale_option);
      evtWeight_inclusive *= eventInfo.pileupWeight;
      evtWeight_inclusive *= lumiScale;
    }

    std::vector<GenParticle> genTopQuarks;
    std::vector<GenParticle> genBJets;
    std::vector<GenParticle> genWBosons;
    std::vector<GenParticle> genWJets;
    std::vector<GenParticle> genQuarkFromTop;
    if ( isMC ) {
      genTopQuarks = genTopQuarkReader->read();
      genBJets = genBJetReader->read();
      genWBosons = genWBosonReader->read();
      genWJets = genWJetReader->read();
      genQuarkFromTop = genQuarkFromTopReader->read();
    }

    bool isTriggered_1e = hltPaths_isTriggered(triggers_1e);
    bool isTriggered_2e = hltPaths_isTriggered(triggers_2e);
    bool isTriggered_1mu = hltPaths_isTriggered(triggers_1mu);
    bool isTriggered_2mu = hltPaths_isTriggered(triggers_2mu);
    bool isTriggered_1e1mu = hltPaths_isTriggered(triggers_1e1mu);

    bool selTrigger_1e = use_triggers_1e && isTriggered_1e;
    bool selTrigger_2e = use_triggers_2e && isTriggered_2e;
    bool selTrigger_1mu = use_triggers_1mu && isTriggered_1mu;
    bool selTrigger_2mu = use_triggers_2mu && isTriggered_2mu;
    bool selTrigger_1e1mu = use_triggers_1e1mu && isTriggered_1e1mu;
    if ( !(selTrigger_1e || selTrigger_2e || selTrigger_1mu || selTrigger_2mu || selTrigger_1e1mu) ) {
      pass_2lss_0tau = false;
      pass_2los_0tau = false;
      pass_ttWctrl = false;
      pass_2lss_1tau = false;
      pass_2los_1tau = false;
      pass_2l_2tau = false;
    }

    bool isTriggered_3e = hltPaths_isTriggered(triggers_3e, isDEBUG);
    bool isTriggered_2e1mu = hltPaths_isTriggered(triggers_2e1mu, isDEBUG);
    bool isTriggered_1e2mu = hltPaths_isTriggered(triggers_1e2mu, isDEBUG);
    bool isTriggered_3mu = hltPaths_isTriggered(triggers_3mu, isDEBUG);

    bool selTrigger_3e = use_triggers_3e && isTriggered_3e;
    bool selTrigger_2e1mu = use_triggers_2e1mu && isTriggered_2e1mu;
    bool selTrigger_1e2mu = use_triggers_1e2mu && isTriggered_1e2mu;
    bool selTrigger_3mu = use_triggers_3mu && isTriggered_3mu;
    if ( !(selTrigger_1e || selTrigger_1mu   ||
     selTrigger_2e || selTrigger_1e1mu || selTrigger_2mu   ||
     selTrigger_3e || selTrigger_2e1mu || selTrigger_1e2mu || selTrigger_3mu) ) {
      pass_4l_0tau = false;
      pass_3l_0tau = false;
      pass_ttZctrl = false;
      pass_WZctrl = false;
      pass_ZZctrl = false;
    }

  bool isTriggered_2tau = hltPaths_isTriggered(triggers_2tau, isDEBUG);
  bool selTrigger_2tau = use_triggers_2tau && isTriggered_2tau;
  if ( !selTrigger_2tau ) {
    pass_0l_4tau = false;
    pass_0l_3tau = false;
    pass_0l_2tau = false;
  }

  bool isTriggered_1e1tau = hltPaths_isTriggered(triggers_1e1tau, isDEBUG);
  bool isTriggered_1mu1tau = hltPaths_isTriggered(triggers_1mu1tau, isDEBUG);

  bool selTrigger_1e1tau = use_triggers_1e1tau && isTriggered_1e1tau;
  bool selTrigger_1mu1tau = use_triggers_1mu1tau && isTriggered_1mu1tau;

  if ( !(selTrigger_1e || selTrigger_1e1tau || selTrigger_1mu || selTrigger_1mu1tau) ) {
    pass_1l_3tau = false;
    pass_1l_2tau = false;
  }
  if ( !(selTrigger_1e || selTrigger_1mu   ||
   selTrigger_2e || selTrigger_1e1mu || selTrigger_2mu   ||
   selTrigger_3e || selTrigger_2e1mu || selTrigger_1e2mu || selTrigger_3mu) ) pass_3l_1tau = false;

//--- build collections of electrons, muons and hadronic taus;
//    resolve overlaps in order of priority: muon, electron,
    std::vector<RecoMuon> muons = muonReader->read();
    std::vector<const RecoMuon*> muon_ptrs = convert_to_ptrs(muons);
    std::vector<const RecoMuon*> cleanedMuons = muon_ptrs; // CV: no cleaning needed for muons, as they have the highest priority in the overlap removal
    std::vector<const RecoMuon*> preselMuons = preselMuonSelector(cleanedMuons, isHigherConePt);
    std::vector<const RecoMuon*> fakeableMuons = fakeableMuonSelector(preselMuons, isHigherConePt);
    std::vector<const RecoMuon*> tightMuons = tightMuonSelector(fakeableMuons, isHigherConePt);
    if(isDEBUG || run_lumi_eventSelector)
    {
      printCollection("preselMuons",   preselMuons);
      printCollection("fakeableMuons", fakeableMuons);
      printCollection("tightMuons",    tightMuons);
    }

    std::vector<RecoElectron> electrons = electronReader->read();
    std::vector<const RecoElectron*> electron_ptrs = convert_to_ptrs(electrons);
    std::vector<const RecoElectron*> cleanedElectrons = electronCleaner(electron_ptrs, preselMuons);
    std::vector<const RecoElectron*> preselElectrons = preselElectronSelector(cleanedElectrons, isHigherConePt);
    std::vector<const RecoElectron*> fakeableElectrons = fakeableElectronSelector(preselElectrons, isHigherConePt);
    std::vector<const RecoElectron*> tightElectrons = tightElectronSelector(fakeableElectrons, isHigherConePt);
    if(isDEBUG || run_lumi_eventSelector)
    {
      printCollection("preselElectrons",   preselElectrons);
      printCollection("fakeableElectrons", fakeableElectrons);
      printCollection("tightElectrons",    tightElectrons);
    }

    std::vector<const RecoLepton*> preselLeptonsFull = mergeLeptonCollections(preselElectrons, preselMuons, isHigherConePt);
    std::vector<const RecoLepton*> fakeableLeptonsFull = mergeLeptonCollections(fakeableElectrons, fakeableMuons, isHigherConePt);
    std::vector<const RecoLepton*> tightLeptonsFull = mergeLeptonCollections(tightElectrons, tightMuons, isHigherConePt);

    std::vector<const RecoLepton*> preselLeptons = pickFirstNobjects(preselLeptonsFull, 4);
    std::vector<const RecoLepton*> fakeableLeptons = pickFirstNobjects(fakeableLeptonsFull, 4);
    std::vector<const RecoLepton*> tightLeptons = getIntersection(fakeableLeptons, tightLeptonsFull, isHigherConePt);

    std::vector<const RecoLepton*> selLeptons;
    std::vector<const RecoMuon*> selMuons;
    std::vector<const RecoElectron*> selElectrons;

    if(electronSelection == muonSelection)
    {
      // for SR, flip region and fake CR
      // doesn't matter if we supply electronSelection or muonSelection here
      selLeptons = selectObjects(muonSelection, preselLeptons, fakeableLeptons, tightLeptons);
      selMuons = getIntersection(preselMuons, selLeptons, isHigherConePt);
      selElectrons = getIntersection(preselElectrons, selLeptons, isHigherConePt);
    }
    else
    {
      // for MC closure
      // make sure that neither electron nor muon selections are loose
      assert(electronSelection != kLoose && muonSelection != kLoose);
      selMuons = selectObjects(muonSelection, preselMuons, fakeableMuons, tightMuons);
      selElectrons = selectObjects(electronSelection, preselElectrons, fakeableElectrons, tightElectrons);
      std::vector<const RecoLepton*> selLeptons_full = mergeLeptonCollections(selElectrons, selMuons, isHigherConePt);
      selLeptons = getIntersection(fakeableLeptons, selLeptons_full, isHigherConePt);
    }

    std::vector<RecoHadTau> hadTaus = hadTauReader->read();
    std::vector<const RecoHadTau*> hadTau_ptrs = convert_to_ptrs(hadTaus);
    std::vector<const RecoHadTau*> cleanedHadTaus = hadTauCleaner(hadTau_ptrs, preselMuons, preselElectrons);
    std::vector<const RecoHadTau*> preselHadTausFull = preselHadTauSelector(cleanedHadTaus, isHigherPt);
    std::vector<const RecoHadTau*> fakeableHadTausFull = fakeableHadTauSelector(preselHadTausFull, isHigherPt);
    std::vector<const RecoHadTau*> tightHadTausFull = tightHadTauSelector(fakeableHadTausFull, isHigherPt);

    std::vector<const RecoHadTau*> preselHadTaus = pickFirstNobjects(preselHadTausFull, 5);
    std::vector<const RecoHadTau*> fakeableHadTaus = pickFirstNobjects(fakeableHadTausFull, 5);
    std::vector<const RecoHadTau*> tightHadTaus = getIntersection(fakeableHadTaus, tightHadTausFull, isHigherPt);
    std::vector<const RecoHadTau*> selHadTaus = selectObjects(hadTauSelection, preselHadTaus, fakeableHadTaus, tightHadTaus);
    if(isDEBUG || run_lumi_eventSelector)
    {
      printCollection("preselHadTaus",   preselHadTaus);
      printCollection("fakeableHadTaus", fakeableHadTaus);
      printCollection("tightHadTaus",    tightHadTaus);
    }

    if(isDEBUG || run_lumi_eventSelector)
    {
      printCollection("selMuons", selMuons);
      printCollection("selElectrons", selElectrons);
      printCollection("selLeptons", selLeptons);
      printCollection("selHadTaus", selHadTaus);
    }

//--- build collections of jets and select subset of jets passing b-tagging criteria
    std::vector<RecoJet> jets = jetReader->read();
    std::vector<const RecoJet*> jet_ptrs = convert_to_ptrs(jets);
    std::vector<const RecoJet*> cleanedJets = jetCleaner(jet_ptrs, selMuons, selElectrons, selHadTaus);
    std::vector<const RecoJet*> selJets = jetSelector(cleanedJets, isHigherPt);
    std::vector<const RecoJet*> selJetsForward = jetSelectorForward(jet_ptrs, isHigherPt);
    std::vector<const RecoJet*> selBJets_loose = jetSelectorBtagLoose(cleanedJets, isHigherPt);
    std::vector<const RecoJet*> selBJets_medium = jetSelectorBtagMedium(cleanedJets, isHigherPt);
    if(isDEBUG || run_lumi_eventSelector)
    {
      printCollection("uncleanedJets", jet_ptrs);
      printCollection("selJets",       selJets);
    }

    //--- build collections of jets reconstructed by hep-top-tagger (HTTv2) algorithm
    std::vector<RecoJetHTTv2> jetsHTTv2 = jetReaderHTTv2->read();
    std::vector<const RecoJetHTTv2*> jet_ptrsHTTv2raw = convert_to_ptrs(jetsHTTv2);
    std::vector<const RecoJetHTTv2*> jet_ptrsHTTv2rawSel = jetCleanerHTTv2SubJets(jet_ptrsHTTv2raw, selMuons, selElectrons); // , selHadTaus //Xanda loose taus would overlap most of times
    std::vector<const RecoJetHTTv2*> sel_HTTv2 = jetSelectorHTTv2(jet_ptrsHTTv2raw, isHigherPt);
    //if (jetsHTTv2.size() > 0) std::cout << " sel_HTTv2.size() " << jetsHTTv2.size() << "\n"; //Xanda

//--- build collections of jets reconstructed by anti-kT algorithm with dR=0.8 (AK8)
    std::vector<RecoJetAK8> jetsAK8 = jetReaderAK8->read();
    std::vector<const RecoJetAK8*> jet_ptrsAK8raw1 = convert_to_ptrs(jetsAK8);
    std::vector<const RecoJetAK8*> jet_ptrsAK8raw = jetCleanerAK8SubJets(jet_ptrsAK8raw1, selMuons, selElectrons, selHadTaus);
    std::vector<const RecoJetAK8*> jet_ptrsAK8 = jetSelectorAK8(jet_ptrsAK8raw, isHigherPt);

    std::vector<const RecoJet*> cleanedJets_fromAK8;
    cleanedJets_fromAK8 = jetCleaner_large8(selJets, jet_ptrsAK8);

//--- build collections of generator level particles (after some cuts are applied, to save computing time)
    if ( isMC && redoGenMatching && !fillGenEvtHistograms ) {
      if ( genLeptonReader ) {
	genLeptons = genLeptonReader->read();
	for ( std::vector<GenLepton>::const_iterator genLepton = genLeptons.begin();
	      genLepton != genLeptons.end(); ++genLepton ) {
	  int abs_pdgId = std::abs(genLepton->pdgId());
	  if      ( abs_pdgId == 11 ) genElectrons.push_back(*genLepton);
	  else if ( abs_pdgId == 13 ) genMuons.push_back(*genLepton);
	}
      }
      if ( genHadTauReader ) {
	genHadTaus = genHadTauReader->read();
      }
      if ( genPhotonReader ) {
        genPhotons = genPhotonReader->read();
      }
      if ( genJetReader ) {
	genJets = genJetReader->read();
      }
    }

//--- match reconstructed to generator level particles
    if ( isMC && redoGenMatching ) {
      muonGenMatcher.addGenLeptonMatch(preselMuons, genLeptons, 0.2);
      muonGenMatcher.addGenHadTauMatch(preselMuons, genHadTaus, 0.2);
      muonGenMatcher.addGenJetMatch(preselMuons, genJets, 0.2);

      electronGenMatcher.addGenLeptonMatch(preselElectrons, genLeptons, 0.2);
      electronGenMatcher.addGenHadTauMatch(preselElectrons, genHadTaus, 0.2);
      electronGenMatcher.addGenPhotonMatch(preselElectrons, genPhotons, 0.2);
      electronGenMatcher.addGenJetMatch(preselElectrons, genJets, 0.2);

      hadTauGenMatcher.addGenLeptonMatch(preselHadTausFull, genLeptons, 0.2);
      hadTauGenMatcher.addGenHadTauMatch(preselHadTausFull, genHadTaus, 0.2);
      hadTauGenMatcher.addGenJetMatch(preselHadTausFull, genJets, 0.2);

      jetGenMatcher.addGenLeptonMatch(selJets, genLeptons, 0.2);
      jetGenMatcher.addGenHadTauMatch(selJets, genHadTaus, 0.2);
      jetGenMatcher.addGenJetMatch(selJets, genJets, 0.2);
    }

    const RecoJet* selJet_lead = nullptr;
    const RecoJet* selJet_sublead = nullptr;
    const RecoJet* selJet_third = nullptr;
    const RecoJet* selJet_fourth = nullptr;
    if (selJets.size() > 0) {
      selJet_lead = selJets[0];
      if (selJets.size() > 1) {
        selJet_sublead = selJets[1];
        if (selJets.size() > 2) {
          selJet_third = selJets[2];
          if (selJets.size() > 3) selJet_fourth = selJets[3];
        }
      }
    }

    // apply requirement on jets (incl. b-tagged jets) on preselection level
    if ( !((int)selJets.size() >= 2) ) {
      pass_2l_2tau = false;
      pass_3l_1tau = false;
      pass_4l_0tau = false;
      pass_0l_4tau = false;
      pass_1l_3tau = false;
      pass_3l_0tau = false;
      pass_ttZctrl = false;
      pass_WZctrl = false;
      pass_ZZctrl = false;
    }
    if ( !((int)selJets.size() == 3  || (int)selJets.size() == 2) ) pass_ttWctrl = false;
    if ( !((int)selJets.size() >= 3 - 2 ) ) {
      pass_1l_2tau = false;
      pass_2lss_1tau = false;
      pass_2los_1tau = false;
      pass_0l_3tau = false;
      //if ( !((int)selJets.size() >= 3 ) ) pass_ttH = false;
    }
    if ( !((int)selJets.size() >= 4 - 2) ) {
      pass_2los_0tau = false;
      pass_2lss_0tau = false;
      pass_0l_2tau = false;
      //if ( !((int)selJets.size() >= 4 ) ) pass_ttH = false;
    }

    if (
      !(
        (selBJets_loose.size() >= 2 || selBJets_medium.size() >= 1) || \
        (selJets.size() >= 2 && selBJets_medium.size() >= 1 && ((selJets.size() - selBJets_loose.size()) + selJetsForward.size() >= 1))
    )
    ) {
        pass_2los_0tau = false;
        pass_2lss_0tau = false;
        pass_3l_0tau = false;
        pass_ttZctrl = false;
        pass_1l_1tau = false;
        pass_0l_3tau = false;
        pass_0l_2tau = false;
    }

    if (
      !(
        (selBJets_loose.size() >= 2 || selBJets_medium.size() >= 1) || \
        (selJets.size() >= 1 && selBJets_medium.size() >= 1 && ((selJets.size() - selBJets_loose.size()) + selJetsForward.size() >= 1))
    )
    ) {
        pass_4l_0tau = false;
        pass_1l_2tau = false;
        pass_2lss_1tau = false;
        pass_2los_1tau = false;
        pass_3l_1tau = false;
        pass_2l_2tau = false;
        pass_1l_3tau = false;
        pass_0l_4tau = false;
    }


    if (
      !((int)selJets.size() >= 4) &&
      !(selJets.size() >= 2 && selBJets_medium.size() >= 1 && ((selJets.size() - selBJets_loose.size()) + selJetsForward.size()) >= 1)
      ) pass_2lss_0tau = false; // remove overlap with ttW

      if (
        !((int)selJets.size() >= 3) &&
        !(selJets.size() >= 2 && selBJets_medium.size() >= 1 && ((selJets.size() - selBJets_loose.size()) + selJetsForward.size()) >= 1)
        ) {
          pass_1l_2tau = false;
          pass_2lss_1tau = false;
          pass_2los_1tau = false;
          pass_0l_3tau = false;
        }

    if ( !(selBJets_loose.size() >= 2 || selBJets_medium.size() >= 1) ) {
    pass_ttWctrl = false; // Xanda: there would be the condition of 0j with |eta| > 2. to pass tH cat to avoid a large fraction of thq will enter the CR...
    }
    if ( !(selBJets_medium.size() == 0 && selBJets_loose.size() < 2) ) {
      pass_WZctrl = false;
      pass_ZZctrl = false;
    }

//--- compute MHT and linear MET discriminant (met_LD)
    RecoMEt met = metReader->read();
    Particle::LorentzVector mht_p4 = compMHT(fakeableLeptons, fakeableHadTaus, selJets);
    double met_LD = compMEt_LD(met.p4(), mht_p4);

//--- apply final event selection

    const RecoLepton* selLepton_lead = nullptr;
    const RecoLepton* selLepton_sublead = nullptr;
    const RecoLepton* selLepton_third = nullptr;
    const RecoLepton* selLepton_fourth = nullptr;
    int selLepton_lead_type = -1000;
    int selLepton_sublead_type = -100;
    int selLepton_third_type = -100;
    int selLepton_fourth_type = -100;
    bool isCharge_Lepton_OS = false;
    if (selLeptons.size() > 0) {
      selLepton_lead = selLeptons[0];
      selLepton_lead_type = getLeptonType(selLepton_lead->pdgId());
      const double minPt_lead = selLepton_lead->is_electron() ? 25. : 30.;
      if ( !(selLepton_lead->cone_pt() > minPt_lead && selLepton_lead->absEta() < 2.1) ) {
        pass_1l_2tau = false;
        pass_1l_1tau = false;
      }
      if (selLeptons.size() > 1) {
        selLepton_sublead = selLeptons[1];
        selLepton_sublead_type = getLeptonType(selLepton_sublead->pdgId());
        const double minPt_sublead = selLepton_sublead->is_electron() ? 15. : 10.;
        if ( !(selLepton_lead->cone_pt() > 25. && selLepton_sublead->cone_pt() > minPt_sublead) ) {
          pass_2lss_0tau = false;
          pass_2los_0tau = false;
          pass_ttWctrl = false;
          pass_2lss_1tau = false;
          pass_2los_1tau = false;
          pass_2l_2tau = false;
        } else {
          if ( selLepton_lead->charge()*selLepton_sublead->charge() < 0 ) {
            pass_2lss_0tau = false;
            pass_2lss_1tau = false;
            pass_ttWctrl = false;
          } else {
            pass_2los_0tau = false;
            pass_2los_1tau = false;
            isCharge_Lepton_OS = true;
          }
        }
        if (selLeptons.size() > 2) {
          selLepton_third = selLeptons[2];
          selLepton_third_type = getLeptonType(selLepton_third->pdgId());
          if ( !(selLepton_lead->cone_pt() > 25 && selLepton_sublead->cone_pt() > 15 && selLepton_third->cone_pt() > 10) )
          {
            pass_3l_1tau = false;
            pass_3l_0tau = false;
            pass_ttZctrl = false;
            pass_WZctrl = false;
          }
          if (selLeptons.size() > 3) {
            selLepton_fourth = selLeptons[3];
            selLepton_fourth_type = getLeptonType(selLepton_fourth->pdgId());
            if ( !(selLepton_lead->cone_pt() > 25 && selLepton_sublead->cone_pt() > 15 && selLepton_third->cone_pt() > 15 && selLepton_fourth->cone_pt() > 10) ) {
              pass_4l_0tau = false;
              pass_ZZctrl = false;
            }
          }
        }
      }
    }

    int lep_match_jet = 0; //selLepton_genMatch.numGenMatchedJets_)
    int lep_match_lep = 0;
    // require exactly two leptons passing tight selection criteria of final event selection
    //const leptonGenMatchEntry& selLepton_genMatch;
    int sum_lep_charge = -1;
    if ( !(selLeptons.size() == 1) ) {
      pass_1l_2tau = false;
      pass_1l_1tau = false;
      pass_1l_3tau = false;
    } else {
      sum_lep_charge = selLepton_lead->charge();
      const leptonGenMatchEntry& selLepton_genMatch = getLeptonGenMatch(
        leptonGenMatch_definitions,
        selLepton_lead
      );
      lep_match_jet = selLepton_genMatch.numGenMatchedJets_;
      lep_match_lep = selLepton_genMatch.numGenMatchedLeptons_;
    }
    if ( !(selLeptons.size() == 2) ) {
      pass_2lss_0tau = false;
      pass_2los_0tau = false;
      pass_2los_1tau = false;
      pass_2lss_1tau = false;
      pass_2l_2tau = false;
      pass_ttWctrl = false;
    } else {
      sum_lep_charge = selLepton_lead->charge() + selLepton_sublead->charge();
      if(
         (selLepton_lead->genLepton() && selLepton_lead->charge() != selLepton_lead->genLepton()->charge()) ||
         (selLepton_sublead->genLepton() && selLepton_sublead->charge() != selLepton_sublead->genLepton()->charge())
       ) {
           pass_2lss_0tau = false;
           pass_2lss_1tau = false;
           pass_2los_0tau = false;
           pass_2los_1tau = false;
         }
         const leptonGenMatchEntry& selLepton_genMatch = getLeptonGenMatch(
           leptonGenMatch_definitions,
           selLepton_lead,
           selLepton_sublead
         );
         lep_match_jet = selLepton_genMatch.numGenMatchedJets_;
         lep_match_lep = selLepton_genMatch.numGenMatchedLeptons_;
    }
    if ( !(selLeptons.size() == 3) ) {
      pass_3l_0tau = false;
      pass_3l_1tau = false;
      pass_WZctrl = false;
      pass_ttZctrl = false;
    } else {
      sum_lep_charge = selLepton_lead->charge() + selLepton_sublead->charge() + selLepton_third->charge();
      const leptonGenMatchEntry& selLepton_genMatch = getLeptonGenMatch(
        leptonGenMatch_definitions,
        selLepton_lead,
        selLepton_sublead,
        selLepton_third
      );
      lep_match_jet = selLepton_genMatch.numGenMatchedJets_;
      lep_match_lep = selLepton_genMatch.numGenMatchedLeptons_;
    }
    if ( !(selLeptons.size() == 4) ) {
      pass_4l_0tau = false;
      pass_ZZctrl = false;
    } else {
      sum_lep_charge = selLepton_lead->charge() + selLepton_sublead->charge() + selLepton_third->charge() + selLepton_fourth->charge();
      const leptonGenMatchEntry& selLepton_genMatch = getLeptonGenMatch(
        leptonGenMatch_definitions,
        selLepton_lead,
        selLepton_sublead,
        selLepton_third,
        selLepton_fourth
      );
      lep_match_jet = selLepton_genMatch.numGenMatchedJets_;
      lep_match_lep = selLepton_genMatch.numGenMatchedLeptons_;
    }
    if ( !(selLeptons.size() == 0) ) {
      pass_0l_4tau = false;
      pass_0l_3tau = false;
      pass_0l_2tau = false;
    }
    //int idxSelLepton_genMatch = selLepton_genMatch.idx_;
    //assert(idxSelLepton_genMatch != kGen_LeptonUndefined2);

//--- compute event-level weight for data/MC correction of b-tagging efficiency and mistag rate
//   (using the method "Event reweighting using scale factors calculated with a tag and probe method",
//    described on the BTV POG twiki https://twiki.cern.ch/twiki/bin/view/CMS/BTagShapeCalibration )
    double evtWeight = 1.;
    double btagWeight = 1.;
    if ( isMC ) {
      if ( apply_DYMCNormScaleFactors ) evtWeight_inclusive *= dyNormScaleFactors->getWeight(genTauLeptons, selJets.size(), selBJets_loose.size(), selBJets_medium.size());
      evtWeight *= evtWeight_inclusive;
      btagWeight = get_BtagWeight(selJets);
      evtWeight *= btagWeight;
      if ( isDEBUG ) {
	std::cout << "lumiScale = " << lumiScale << std::endl;
	if ( apply_genWeight ) std::cout << "genWeight = " << boost::math::sign(eventInfo.genWeight) << std::endl;
	std::cout << "pileupWeight = " << eventInfo.pileupWeight << std::endl;
	std::cout << "btagWeight = " << btagWeight << std::endl;
      }
    }

    // require that trigger paths match event category (with event category based on fakeableLeptons)
    if ( !((fakeableElectrons.size() >= 2 &&                              (selTrigger_2e    || selTrigger_1e                  )) ||
	   (fakeableElectrons.size() >= 1 && fakeableMuons.size() >= 1 && (selTrigger_1e1mu || selTrigger_1mu || selTrigger_1e)) ||
	   (                                 fakeableMuons.size() >= 2 && (selTrigger_2mu   || selTrigger_1mu                 ))) ) pass_2l_2tau = false;

    const RecoHadTau* selHadTau_lead = nullptr;
    const RecoHadTau* selHadTau_sublead = nullptr;
    const RecoHadTau* selHadTau_third = nullptr;
    const RecoHadTau* selHadTau_fourth = nullptr;
    bool isCharge_hadTau_OS = false;
    int sum_tau_charge = -1;
    double mTauTauVis       = -1000;
    double cosThetaS_hadTau = -1000;
    double mTauTauVis2       = -1000;
    double cosThetaS_hadTau2 = -1000;
    int selHadTau_lead_genPdgId = -1000;
    int selHadTau_sublead_genPdgId = -1000;
    int selHadTau_third_genPdgId = -1000;
    int selHadTau_fourth_genPdgId = -1000;
    if (selHadTaus.size() > 0) {
      selHadTau_lead = selHadTaus[0];
      selHadTau_lead_genPdgId  = getHadTau_genPdgId(selHadTau_lead);
      if ( !(selHadTau_lead->pt() > 30) ) {
        pass_1l_2tau = false;
        pass_1l_1tau = false;
        pass_0l_2tau = false;
        pass_0l_3tau = false;
      }
      if (selLepton_lead)
        if ( selLepton_lead->charge()*selHadTau_lead->charge() < 0 ) pass_1l_1tau = false; // isCharge_SS
      //if (!(selHadTau_lead->absEta() < 2.1)) {
      //  pass_0l_2tau = false;
      //  pass_0l_3tau = false;
      //}
      if (selHadTaus.size() > 1) {
        selHadTau_sublead = selHadTaus[1];
        selHadTau_sublead_genPdgId  = getHadTau_genPdgId(selHadTau_sublead);
        if ( !(selHadTau_sublead->pt() > 20) ) {
          pass_0l_2tau = false;
          pass_0l_3tau = false;
          pass_1l_2tau = false;
        }
        //if (!(selHadTau_sublead->absEta() < 2.1)) {
        //  pass_0l_2tau = false;
        //  pass_0l_3tau = false;
        //}
        if ( selHadTau_lead->charge()*selHadTau_sublead->charge() > 0 ) {
          pass_1l_2tau = false; // isCharge_OS
          pass_0l_2tau = false; // isCharge_OS
        } else {
          if (selLepton_lead && selLepton_sublead)
            if ( std::abs(selLepton_lead->charge() + selHadTau_lead->charge() + selLepton_sublead->charge() + selHadTau_sublead->charge()) > 0.1 ) pass_2l_2tau = false; // isCharge_SS
          if ( selLepton_lead ) if ( std::abs(selLepton_lead->charge() + selHadTau_lead->charge()+ selHadTau_sublead->charge()) != 1 ) pass_1l_2tau = false;
          isCharge_hadTau_OS = false;
        }
        if (selHadTaus.size() > 2) {
          selHadTau_third = selHadTaus[2];
          selHadTau_third_genPdgId  = getHadTau_genPdgId(selHadTau_third);
          if ( std::abs(selHadTau_lead->charge()+ selHadTau_sublead->charge() + selHadTau_third->charge()) != 1 ) pass_0l_3tau = false;
          if (selLeptons.size() > 0) if ( std::abs(selLepton_lead->charge() + selHadTau_lead->charge() + selHadTau_sublead->charge() + selHadTau_third->charge()) > 0.1 ) pass_1l_3tau = false;
          if (selHadTaus.size() > 3) {
            selHadTau_fourth = selHadTaus[3];
            selHadTau_fourth_genPdgId  = getHadTau_genPdgId(selHadTau_fourth);
            if ( std::abs(selHadTau_lead->charge() + selHadTau_sublead->charge() + selHadTau_third->charge() + selHadTau_fourth->charge()) > 0.1 ) pass_0l_4tau = false;
          }
        }
      }
    }

    int tau_match_jet = 0;
    int tau_match_tau = 0;
    if ( !(selHadTaus.size() == 1) ) {
      pass_3l_1tau = false;
      pass_2lss_1tau = false;
      pass_2los_1tau = false;
      pass_1l_1tau = false;
    } else {
      sum_tau_charge = selHadTau_lead->charge();
      const hadTauGenMatchEntry& selHadTau_genMatch = getHadTauGenMatch(
        hadTauGenMatch_definitions,
        selHadTau_lead
      );
      tau_match_jet = selHadTau_genMatch.numGenMatchedJets_;
      tau_match_tau = selHadTau_genMatch.numGenMatchedHadTaus_;
    }
    if ( !(selHadTaus.size() == 2) ) {
      pass_2l_2tau = false;
      pass_1l_2tau = false;
      pass_0l_2tau = false;
    } else {
      sum_tau_charge = selHadTau_lead->charge() + selHadTau_sublead->charge();
      const hadTauGenMatchEntry& selHadTau_genMatch = getHadTauGenMatch(
        hadTauGenMatch_definitions,
        selHadTau_lead,
        selHadTau_sublead
      );
      tau_match_jet = selHadTau_genMatch.numGenMatchedJets_;
      tau_match_tau = selHadTau_genMatch.numGenMatchedHadTaus_;
    }
    if ( !(selHadTaus.size() == 3) ) {
      pass_1l_3tau = false;
      pass_0l_3tau = false;
    } else {
      sum_tau_charge = selHadTau_lead->charge() + selHadTau_sublead->charge() + selHadTau_third->charge();
      const hadTauGenMatchEntry& selHadTau_genMatch = getHadTauGenMatch(
        hadTauGenMatch_definitions,
        selHadTau_lead,
        selHadTau_sublead,
        selHadTau_third
      );
      tau_match_jet = selHadTau_genMatch.numGenMatchedJets_;
      tau_match_tau = selHadTau_genMatch.numGenMatchedHadTaus_;
    }
    if ( !(selHadTaus.size() == 4) ) {
      pass_0l_4tau = false;
    } else {
      sum_tau_charge = selHadTau_lead->charge() + selHadTau_sublead->charge() + selHadTau_third->charge() + selHadTau_fourth->charge();
      const hadTauGenMatchEntry& selHadTau_genMatch = getHadTauGenMatch(
        hadTauGenMatch_definitions,
        selHadTau_lead,
        selHadTau_sublead,
        selHadTau_third,
        selHadTau_fourth
      );
      tau_match_jet = selHadTau_genMatch.numGenMatchedJets_;
      tau_match_tau = selHadTau_genMatch.numGenMatchedHadTaus_;
    }
    if ( !(selHadTaus.size() == 0) ) {
      pass_ttZctrl = false;
      pass_WZctrl = false;
      //pass_ZZctrl = false;
      pass_ttWctrl = false;
      pass_2lss_0tau = false;
      pass_2los_0tau = false;
      pass_3l_0tau = false;
      //pass_4l_0tau = false;
    }

    bool failsLowMassVeto = false;
    for ( std::vector<const RecoLepton*>::const_iterator lepton1 = preselLeptonsFull.begin();
    lepton1 != preselLeptonsFull.end(); ++lepton1 ) {
      for ( std::vector<const RecoLepton*>::const_iterator lepton2 = lepton1 + 1;
      lepton2 != preselLeptonsFull.end(); ++lepton2 ) {
	double mass = ((*lepton1)->p4() + (*lepton2)->p4()).mass();
	if ( mass < 12. ) {
	  failsLowMassVeto = true;
	}
      }
    }
    if ( failsLowMassVeto ) {
      pass_2lss_0tau = false;
      pass_2los_0tau = false;
      pass_3l_0tau = false;
      pass_4l_0tau = false;
      pass_ttZctrl = false;
      pass_ttWctrl = false;
      pass_2lss_1tau = false;
      pass_2los_1tau = false;
      pass_3l_1tau = false;
      pass_2l_2tau = false;
      pass_1l_2tau = false;
      pass_1l_1tau = false;
      pass_WZctrl = false;
      pass_ZZctrl = false;
    }

    bool isSameFlavor_OS = false;
    double massSameFlavor_OS = -1.;
    for ( std::vector<const RecoLepton*>::const_iterator lepton1 = preselLeptonsFull.begin();
    lepton1 != preselLeptonsFull.end(); ++lepton1 ) {
      for ( std::vector<const RecoLepton*>::const_iterator lepton2 = lepton1 + 1;
      lepton2 != preselLeptonsFull.end(); ++lepton2 ) {
	if ( (*lepton1)->pdgId() == -(*lepton2)->pdgId() ) { // pair of same flavor leptons of opposite charge
	  isSameFlavor_OS = true;
	  double mass = ((*lepton1)->p4() + (*lepton2)->p4()).mass();
	  if ( std::fabs(mass - z_mass) < std::fabs(massSameFlavor_OS - z_mass) ) massSameFlavor_OS = mass;
	}
      }
    }

    bool failsZbosonMassVeto = isSameFlavor_OS && std::fabs(massSameFlavor_OS - z_mass) < z_window;
    if ( failsZbosonMassVeto ) {
      pass_2lss_0tau = false;
      pass_2los_0tau = false;
      pass_3l_0tau = false;
      pass_4l_0tau = false;
      pass_ttWctrl = false;
      pass_2lss_1tau = false;
      pass_2los_1tau = false;
      pass_3l_1tau = false;
      pass_2l_2tau = false;
      pass_1l_2tau = false;
      pass_1l_1tau = false;
    } else {
      pass_ttZctrl = false;
      pass_WZctrl = false;
      pass_ZZctrl = false;
    }

    double met_LD_cut = 0.;
    if      ( selJets.size() >= 4 ) met_LD_cut = -1.; // MET LD cut not applied
    else if ( isSameFlavor_OS && selLeptons.size() > 1 ) met_LD_cut = 45.;
    else                            met_LD_cut = 30.;
    if ( met_LD_cut > 0 && met_LD < met_LD_cut ) {
      pass_3l_0tau = false;
      pass_4l_0tau = false;
      pass_ttWctrl = false;
      pass_2los_1tau = false;
      pass_3l_1tau = false;
      pass_2l_2tau = false;
      //pass_1l_2tau = false;
      //pass_1l_1tau = false;
      pass_WZctrl = false;
      pass_ZZctrl = false;
    }
    if (pass_2lss_0tau || pass_2lss_1tau || pass_ttZctrl) {
      if ( (selLepton_lead->is_electron() && selLepton_sublead->is_electron()) && met_LD < 30. ) {
      pass_2lss_0tau = false;
      pass_2lss_1tau = false;
      pass_ttZctrl = false;
    }
    }

    bool failsHtoZZVeto = false;
    for ( std::vector<const RecoLepton*>::const_iterator lepton1 = preselLeptonsFull.begin();
    lepton1 != preselLeptonsFull.end(); ++lepton1 ) {
      for ( std::vector<const RecoLepton*>::const_iterator lepton2 = lepton1 + 1;
      lepton2 != preselLeptonsFull.end(); ++lepton2 ) {
	if ( (*lepton1)->pdgId() == -(*lepton2)->pdgId() ) { // first pair of same flavor leptons of opposite charge
    for ( std::vector<const RecoLepton*>::const_iterator lepton3 = preselLeptonsFull.begin();
    lepton3 != preselLeptonsFull.end(); ++lepton3 ) {
	    if ( (*lepton3) == (*lepton1) || (*lepton3) == (*lepton2) ) continue;
	    for ( std::vector<const RecoLepton*>::const_iterator lepton4 = lepton3 + 1;
      lepton4 != preselLeptonsFull.end(); ++lepton4 ) {
	      if ( (*lepton4) == (*lepton1) || (*lepton4) == (*lepton2) ) continue;
	      if ( (*lepton3)->pdgId() == -(*lepton4)->pdgId() ) { // second pair of same flavor leptons of opposite charge
		double mass = ((*lepton1)->p4() + (*lepton2)->p4() + (*lepton3)->p4() + (*lepton4)->p4()).mass();
		if ( mass < 140. ) failsHtoZZVeto = true;
	      }
	    }
	  }
	}
      }
    }
    if ( failsHtoZZVeto ) {
      pass_3l_0tau = false;
      pass_4l_0tau = false;
      pass_ttZctrl = false;
      pass_WZctrl = false;
      pass_ZZctrl = false;
    }

    /// Calculate pair with largest mTTvis from tau pair with oposite charge
    // calculate as well SVFit mass
    double mTauTau = -1000;
    double mTauTau2 = -1000;
    ClassicSVfit svFitAlgo;
    svFitAlgo.addLogM_dynamic(false);
    svFitAlgo.addLogM_fixed(true, 4.);
    if ( pass_0l_2tau || pass_0l_3tau || pass_0l_4tau || pass_1l_3tau || pass_1l_2tau ) {
      const RecoHadTau* tau1 = nullptr;
      const RecoHadTau* tau2 = nullptr;
      for ( std::vector<const RecoHadTau*>::const_iterator lepton1 = selHadTaus.begin();
      lepton1 != selHadTaus.end(); ++lepton1 ) {
        for ( std::vector<const RecoHadTau*>::const_iterator lepton2 = lepton1 + 1;
        lepton2 != selHadTaus.end(); ++lepton2 ) {
        double mass = ((*lepton1)->p4() + (*lepton2)->p4()).mass();
        if ( (*lepton1)->charge()*(*lepton2)->charge() < 0 && mass > mTauTauVis) {
          mTauTauVis       = mass;
          cosThetaS_hadTau = comp_cosThetaS((*lepton1)->p4(), (*lepton2)->p4());
          tau1 = (*lepton1);
          tau2 = (*lepton2);
        }
      }
      }
      //--- reconstruct mass of tau-pair using SVfit algorithm
      //
      //    NOTE: SVfit needs to be run after all event selection cuts are applied,
      //          because the algorithm takes O(1 second per event) to run
      //
      std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
      classic_svFit::MeasuredTauLepton::kDecayType leg1Type = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
      double leg1Mass = tau1->mass();
      if ( leg1Mass < classic_svFit::chargedPionMass ) leg1Mass = classic_svFit::chargedPionMass;
      if ( leg1Mass > 1.5                            ) leg1Mass = 1.5;
      classic_svFit::MeasuredTauLepton::kDecayType leg2Type = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
      double leg2Mass = tau2->mass();
      if ( leg2Mass < classic_svFit::chargedPionMass ) leg2Mass = classic_svFit::chargedPionMass;
      if ( leg2Mass > 1.5                            ) leg2Mass = 1.5;
      measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(leg1Type, tau1->pt(), tau1->eta(), tau1->phi(), leg1Mass));
      measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(leg2Type, tau2->pt(), tau2->eta(), tau2->phi(), leg2Mass));
      svFitAlgo.integrate(measuredTauLeptons, met.p4().px(), met.p4().py(), met.cov());
      mTauTau   = ( svFitAlgo.isValidSolution() ) ? static_cast<classic_svFit::HistogramAdapterDiTau*>(svFitAlgo.getHistogramAdapter())->getMass() : -1.;

      if ( pass_0l_4tau )
      {
        const RecoHadTau* tau3 = nullptr;
        const RecoHadTau* tau4 = nullptr;
        for ( std::vector<const RecoHadTau*>::const_iterator lepton1 = selHadTaus.begin();
        lepton1 != selHadTaus.end(); ++lepton1 ) {
          if ( (*lepton1) == tau1 || (*lepton1) == tau2 ) continue;
          for ( std::vector<const RecoHadTau*>::const_iterator lepton2 = lepton1 + 1;
          lepton2 != selHadTaus.end(); ++lepton2 ) {
            if ( (*lepton1) == tau1 || (*lepton1) == tau2 ) continue;
            double mass = ((*lepton1)->p4() + (*lepton2)->p4()).mass();
            if ( (*lepton1)->charge()*(*lepton2)->charge() < 0 && mass > mTauTauVis2) {
            mTauTauVis2       = mass;
            cosThetaS_hadTau2 = comp_cosThetaS((*lepton1)->p4(), (*lepton2)->p4());
            tau3 = (*lepton1);
            tau4 = (*lepton2);
            }
        }
        }
        std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons2;
        double leg1Mass2 = tau1->mass();
        if ( leg1Mass2 < classic_svFit::chargedPionMass ) leg1Mass2 = classic_svFit::chargedPionMass;
        if ( leg1Mass2 > 1.5                            ) leg1Mass2 = 1.5;
        double leg2Mass2 = tau2->mass();
        if ( leg2Mass2 < classic_svFit::chargedPionMass ) leg2Mass2 = classic_svFit::chargedPionMass;
        if ( leg2Mass2 > 1.5                            ) leg2Mass2 = 1.5;
        measuredTauLeptons2.push_back(classic_svFit::MeasuredTauLepton(leg1Type, tau3->pt(), tau3->eta(), tau3->phi(), leg1Mass2));
        measuredTauLeptons2.push_back(classic_svFit::MeasuredTauLepton(leg2Type, tau4->pt(), tau4->eta(), tau4->phi(), leg2Mass2));
        svFitAlgo.integrate(measuredTauLeptons2, met.p4().px(), met.p4().py(), met.cov());
        //double mTauTau = -1.; // CV: temporarily comment-out the following line, to make code compile with "old" and "new" version of ClassicSVfit
        mTauTau   = ( svFitAlgo.isValidSolution() ) ? static_cast<classic_svFit::HistogramAdapterDiTau*>(svFitAlgo.getHistogramAdapter())->getMass() : -1.;
      }
    }

    if ( pass_1l_1tau ) {
      std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
      classic_svFit::MeasuredTauLepton::kDecayType leg1Type = classic_svFit::MeasuredTauLepton::kUndefinedDecayType;
      double leg1Mass = -1.;
      if ( selLepton_lead_type == kElectron ) {
        leg1Type = classic_svFit::MeasuredTauLepton::kTauToElecDecay;
        leg1Mass = classic_svFit::electronMass;
      } else if ( selLepton_lead_type == kMuon ) {
        leg1Type = classic_svFit::MeasuredTauLepton::kTauToMuDecay;
        leg1Mass = classic_svFit::muonMass;
      } else assert(0);
      classic_svFit::MeasuredTauLepton::kDecayType leg2Type = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
      double leg2Mass = selHadTau_lead->mass();
      if ( leg2Mass < classic_svFit::chargedPionMass ) leg2Mass = classic_svFit::chargedPionMass;
      if ( leg2Mass > 1.5                            ) leg2Mass = 1.5;
      measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(leg1Type, selLepton_lead->pt(), selLepton_lead->eta(), selLepton_lead->phi(), leg1Mass));
      measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(leg2Type, selHadTau_lead->pt(), selHadTau_lead->eta(), selHadTau_lead->phi(), leg2Mass));
      svFitAlgo.integrate(measuredTauLeptons, met.p4().px(), met.p4().py(), met.cov());
      mTauTau = ( svFitAlgo.isValidSolution() ) ? static_cast<classic_svFit::HistogramAdapterDiTau*>(svFitAlgo.getHistogramAdapter())->getMass() : -1.;
    }

//--- Declare the variables used as an input to the MVA/BDT in one place
//    so that there won't be any mismatches b/w the variables in the BDT Ntuple and
//    the variables used to evaluate the MVA/BDT scores.
//    Besides, we may want to use the said variables to fill sync Ntuple as well.

    if ( (pass_2los_0tau || pass_2lss_0tau || pass_0l_2tau) && !((int)selJets.size() >= 4 ) ) pass_ttH = false;
    if ( (pass_1l_2tau || pass_2lss_1tau || pass_2los_1tau || pass_0l_3tau) && !((int)selJets.size() >= 3 ) ) pass_ttH = false;
    if ( !((int)selJets.size() == 3) && pass_ttWctrl ) pass_ttH = false;
    if ( !(selBJets_loose.size() >= 2 || selBJets_medium.size() >= 1) && !(pass_WZctrl || pass_ZZctrl)) pass_ttH = false;

    if (! (
      pass_2lss_0tau || \
      pass_2los_0tau || \
      pass_3l_0tau || \
      pass_4l_0tau || \
      pass_ttZctrl || \
      pass_ttWctrl || \
      pass_1l_2tau || \
      pass_2lss_1tau || \
      pass_2los_1tau || \
      pass_3l_1tau || \
      pass_2l_2tau || \
      pass_1l_1tau  || \
      pass_1l_3tau || \
      pass_0l_4tau || \
      pass_0l_3tau || \
      pass_0l_2tau || \
      pass_WZctrl || \
      pass_ZZctrl)) continue;

    /* pass_2lss_0tau and pass_ttWctrl may overlap
    if (((pass_2lss_0tau) + \
    pass_2los_0tau + \
    pass_3l_0tau + \
    pass_4l_0tau + \
    pass_ttZctrl + \
    pass_ttWctrl + \
    pass_1l_2tau + \
    pass_2lss_1tau + \
    pass_2los_1tau + \
    pass_3l_1tau + \
    pass_2l_2tau + \
    pass_1l_1tau  + \
    pass_1l_3tau + \
    pass_0l_4tau + \
    pass_0l_3tau + \
    pass_0l_2tau + \
    pass_WZctrl + \
    pass_ZZctrl) > 1 ) {
      std::cout<<"event with categories overlap "<<eventInfo.run<<" "<< eventInfo.lumi<<" "<< eventInfo.event << "\n" << " selHadTaus.size() = "<< selHadTaus.size() << "\n" ;
      std::cout<<
      "(pass_2lss_0tau && pass_ttH) "<< (pass_2lss_0tau && pass_ttH) << "\n" <<
      "pass_2lss_0tau  "<< pass_2lss_0tau << "\n" <<
      "pass_2los_0tau  "<< pass_2los_0tau << "\n" <<
      "pass_3l_0tau   "<< pass_3l_0tau << "\n" <<
      "pass_4l_0tau   "<< pass_4l_0tau << "\n" <<
      "pass_ttZctrl   "<< pass_ttZctrl << "\n" <<
      "pass_ttWctrl   "<< pass_ttWctrl << "\n" <<
      "pass_1l_2tau   "<< pass_1l_2tau << "\n" <<
      "pass_2lss_1tau "<< pass_2lss_1tau << "\n" <<
      "pass_2los_1tau "<< pass_2los_0tau << "\n" <<
      "pass_3l_1tau   "<< pass_3l_1tau << "\n" <<
      "pass_2l_2tau   "<< pass_2l_2tau << "\n" <<
      "pass_1l_1tau   "<< pass_1l_1tau << "\n" <<
      "pass_1l_3tau   "<< pass_1l_3tau << "\n" <<
      "pass_0l_4tau   "<< pass_0l_4tau << "\n" <<
      "pass_0l_3tau   "<< pass_0l_3tau << "\n" <<
      "pass_0l_2tau   "<< pass_0l_2tau << "\n" <<
      "pass_WZctrl    "<< pass_WZctrl << "\n" <<
      "pass_ZZctrl    "<< pass_ZZctrl << "\n";
    }
    */

    // I only use tight guys, no fake rate is needed
    // we do apply dataToMC corrections
    double sf_triggerEff = 1.0;
    double weight_data_to_MC_correction = 1.;
    double triggerWeight = 1.;
    double weight_leptonEff = 1.;
    double weight_hadTauEff = 1.;
    double tauSF_weight = 1.;
    double sf_leptonEff = 1.;
    double sf_hadTauEff = 1.;
    double weight_data_to_MC_correction_hadTau = 1.;
    if ( selLeptons.size() > 0 ) {
      if (selLeptons.size() == 1) dataToMCcorrectionInterface->setLeptons(selLepton_lead_type, selLepton_lead->pt(), selLepton_lead->eta());
      if (selLeptons.size() == 2) dataToMCcorrectionInterface->setLeptons(selLepton_lead_type, selLepton_lead->pt(), selLepton_lead->eta(),
      selLepton_sublead_type, selLepton_sublead->pt(), selLepton_sublead->eta());
      if (selLeptons.size() == 3) dataToMCcorrectionInterface->setLeptons(
        selLepton_lead_type, selLepton_lead->pt(), selLepton_lead->eta(),
        selLepton_sublead_type, selLepton_sublead->pt(), selLepton_sublead->eta(),
        selLepton_third_type, selLepton_third->pt(), selLepton_third->eta());
      if (selLeptons.size() == 4) dataToMCcorrectionInterface->setLeptons(
        selLepton_lead_type, selLepton_lead->pt(), selLepton_lead->eta(),
        selLepton_sublead_type, selLepton_sublead->pt(), selLepton_sublead->eta(),
        selLepton_third_type, selLepton_third->pt(), selLepton_third->eta(),
        selLepton_fourth_type, selLepton_fourth->pt(), selLepton_fourth->eta());
      sf_triggerEff = dataToMCcorrectionInterface->getSF_leptonTriggerEff();
      //--- apply data/MC corrections for efficiencies for lepton to pass loose identification and isolation criteria
      sf_leptonEff *= dataToMCcorrectionInterface->getSF_leptonID_and_Iso_loose();
      //--- apply data/MC corrections for efficiencies of leptons passing the loose identification and isolation criteria
      //    to also pass the tight identification and isolation criteria
       if ( electronSelection == kFakeable && muonSelection == kFakeable ) {
        sf_leptonEff *= dataToMCcorrectionInterface->getSF_leptonID_and_Iso_fakeable_to_loose();
       } else if ( electronSelection >= kFakeable && muonSelection >= kFakeable ) {
        // apply loose-to-tight lepton ID SFs if either of the following is true:
        // 1) both electron and muon selections are tight -> corresponds to SR
        // 2) electron selection is relaxed to fakeable and muon selection is kept as tight -> corresponds to MC closure w/ relaxed e selection
        // 3) muon selection is relaxed to fakeable and electron selection is kept as tight -> corresponds to MC closure w/ relaxed mu selection
        // we allow (2) and (3) so that the MC closure regions would more compatible w/ the SR (1) in comparison
        sf_leptonEff *= dataToMCcorrectionInterface->getSF_leptonID_and_Iso_tight_to_loose_woTightCharge();
      }
    }

    if (selHadTaus.size() > 0 ) {
      //--- apply data/MC corrections for hadronic tau identification efficiency
      //    and for e->tau and mu->tau misidentification rates
      if (selHadTaus.size() == 1) dataToMCcorrectionInterface->setHadTaus(selHadTau_lead_genPdgId, selHadTau_lead->pt(), selHadTau_lead->eta());
      if (selHadTaus.size() == 2) dataToMCcorrectionInterface->setHadTaus(
        selHadTau_lead_genPdgId, selHadTau_lead->pt(), selHadTau_lead->eta(),
        selHadTau_sublead_genPdgId, selHadTau_sublead->pt(), selHadTau_sublead->eta());
      if (selHadTaus.size() == 3) dataToMCcorrectionInterface->setHadTaus(
        selHadTau_lead_genPdgId, selHadTau_lead->pt(), selHadTau_lead->eta(),
        selHadTau_sublead_genPdgId, selHadTau_sublead->pt(), selHadTau_sublead->eta(),
        selHadTau_third_genPdgId, selHadTau_third->pt(), selHadTau_third->eta()
      );
      if (selHadTaus.size() == 4) dataToMCcorrectionInterface->setHadTaus(
        selHadTau_lead_genPdgId, selHadTau_lead->pt(), selHadTau_lead->eta(),
        selHadTau_sublead_genPdgId, selHadTau_sublead->pt(), selHadTau_sublead->eta(),
        selHadTau_third_genPdgId, selHadTau_third->pt(), selHadTau_third->eta(),
        selHadTau_fourth_genPdgId, selHadTau_fourth->pt(), selHadTau_fourth->eta()
      );
      sf_hadTauEff *= dataToMCcorrectionInterface->getSF_hadTauID_and_Iso();
      sf_hadTauEff *= dataToMCcorrectionInterface->getSF_eToTauFakeRate();
      sf_hadTauEff *= dataToMCcorrectionInterface->getSF_muToTauFakeRate();
    }
    //--- apply data/MC corrections for trigger efficiency
    if (pass_1l_1tau) {
      dataToMCcorrectionInterface_1l_1tau_trigger->setLeptons(selLepton_lead_type, selLepton_lead->pt(), selLepton_lead->eta());
      dataToMCcorrectionInterface_1l_1tau_trigger->setHadTaus(selHadTau_lead->pt(), selHadTau_lead->eta(), selHadTau_lead->phi());
      dataToMCcorrectionInterface_1l_1tau_trigger->setTriggerBits(isTriggered_1e, isTriggered_1e1tau, isTriggered_1mu, isTriggered_1mu1tau);
      sf_triggerEff = dataToMCcorrectionInterface_1l_1tau_trigger->getSF_triggerEff();
    }
    if (pass_0l_2tau || pass_0l_3tau || pass_0l_4tau) {
      dataToMCcorrectionInterface_0l_2tau_trigger->setHadTaus(
        selHadTau_lead->pt(),    selHadTau_lead->eta(),    selHadTau_lead->phi(),
        selHadTau_sublead->pt(), selHadTau_sublead->eta(), selHadTau_sublead->phi()
      );
      dataToMCcorrectionInterface_0l_2tau_trigger->setTriggerBits(isTriggered_2tau);
      sf_triggerEff = dataToMCcorrectionInterface_0l_2tau_trigger->getSF_triggerEff();
    }
    if (pass_1l_2tau) {
      dataToMCcorrectionInterface_1l_2tau_trigger->setLeptons(selLepton_lead_type, selLepton_lead->pt(), selLepton_lead->eta());
      dataToMCcorrectionInterface_1l_2tau_trigger->setHadTaus(
        selHadTau_lead->pt(),    selHadTau_lead->eta(),    selHadTau_lead->phi(),
        selHadTau_sublead->pt(), selHadTau_sublead->eta(), selHadTau_sublead->phi()
      );
      dataToMCcorrectionInterface_1l_2tau_trigger->setTriggerBits(isTriggered_1e, isTriggered_1e1tau, isTriggered_1mu, isTriggered_1mu1tau);
      sf_triggerEff = dataToMCcorrectionInterface_1l_2tau_trigger->getSF_triggerEff();
    }
    if (pass_2l_2tau || pass_2los_1tau || pass_2lss_1tau || pass_2lss_0tau || pass_ttWctrl || pass_WZctrl || pass_3l_0tau || pass_3l_1tau || pass_ZZctrl || pass_4l_0tau) sf_triggerEff = dataToMCcorrectionInterface->getSF_leptonTriggerEff();

    triggerWeight *= sf_triggerEff;
    weight_data_to_MC_correction *= sf_triggerEff;
    weight_leptonEff *= sf_leptonEff;
    weight_data_to_MC_correction *= sf_leptonEff;
    weight_hadTauEff *= sf_hadTauEff;
    tauSF_weight *= weight_hadTauEff;
    weight_data_to_MC_correction *= sf_hadTauEff;
    tauSF_weight *= weight_data_to_MC_correction_hadTau;

    evtWeight *= weight_data_to_MC_correction;
    evtWeight *= weight_data_to_MC_correction_hadTau;

    // Fake rates
    double weight_fakeRate = 1.;
    double prob_fake_lepton_lead = 1.;
    double prob_fake_lepton_sublead = 1.;
    double prob_fake_lepton_third = 1.;
    double prob_fake_lepton_fourth = 1.;
    bool passesTight_lepton_lead = false;
    bool passesTight_lepton_sublead = false;
    bool passesTight_lepton_third = false;
    bool passesTight_lepton_fourth = false;
    if (selLeptons.size() > 0) {
      if      ( std::abs(selLepton_lead->pdgId()) == 11 ) prob_fake_lepton_lead = leptonFakeRateInterface->getWeight_e(selLepton_lead->cone_pt(), selLepton_lead->absEta());
      else if ( std::abs(selLepton_lead->pdgId()) == 13 ) prob_fake_lepton_lead = leptonFakeRateInterface->getWeight_mu(selLepton_lead->cone_pt(), selLepton_lead->absEta());
      else assert(0);
      passesTight_lepton_lead = isMatched(*selLepton_lead, tightElectrons) || isMatched(*selLepton_lead, tightMuons);
    }
    if (selLeptons.size() > 1) {
      if      ( std::abs(selLepton_sublead->pdgId()) == 11 ) prob_fake_lepton_sublead = leptonFakeRateInterface->getWeight_e(selLepton_sublead->cone_pt(), selLepton_sublead->absEta());
    	else if ( std::abs(selLepton_sublead->pdgId()) == 13 ) prob_fake_lepton_sublead = leptonFakeRateInterface->getWeight_mu(selLepton_sublead->cone_pt(), selLepton_sublead->absEta());
    	else assert(0);
      passesTight_lepton_sublead = isMatched(*selLepton_sublead, tightElectrons) || isMatched(*selLepton_sublead, tightMuons);
    }
    if (selLeptons.size() > 2) {
      if      ( std::abs(selLepton_third->pdgId()) == 11 ) prob_fake_lepton_third = leptonFakeRateInterface->getWeight_e(selLepton_third->cone_pt(), selLepton_third->absEta());
      else if ( std::abs(selLepton_third->pdgId()) == 13 ) prob_fake_lepton_third = leptonFakeRateInterface->getWeight_mu(selLepton_third->cone_pt(), selLepton_third->absEta());
      else assert(0);
      passesTight_lepton_third = isMatched(*selLepton_third, tightElectrons) || isMatched(*selLepton_third, tightMuons);
    }
    if (selLeptons.size() > 3) {
      if      ( std::abs(selLepton_fourth->pdgId()) == 11 ) prob_fake_lepton_fourth = leptonFakeRateInterface->getWeight_e(selLepton_fourth->cone_pt(), selLepton_fourth->absEta());
      else if ( std::abs(selLepton_fourth->pdgId()) == 13 ) prob_fake_lepton_fourth = leptonFakeRateInterface->getWeight_mu(selLepton_fourth->cone_pt(), selLepton_fourth->absEta());
      else assert(0);
      passesTight_lepton_fourth = isMatched(*selLepton_fourth, tightElectrons) || isMatched(*selLepton_fourth, tightMuons);
    }
    /////
    double prob_fake_hadTau_lead = 1.;
    double prob_fake_hadTau_sublead = 1.;
    double prob_fake_hadTau_third = 1.;
    double prob_fake_hadTau_fourth = 1.;
    bool passesTight_hadTau_lead = false;
    bool passesTight_hadTau_sublead = false;
    bool passesTight_hadTau_third = false;
    bool passesTight_hadTau_fourth = false;
    if (selHadTaus.size() > 0) {
      prob_fake_hadTau_lead = jetToTauFakeRateInterface->getWeight_lead(selHadTau_lead->pt(), selHadTau_lead->absEta());
      passesTight_hadTau_lead = isMatched(*selHadTau_lead, tightHadTausFull);
    }
    if (selHadTaus.size() > 1) {
      prob_fake_hadTau_sublead = jetToTauFakeRateInterface->getWeight_lead(selHadTau_sublead->pt(), selHadTau_sublead->absEta());
      passesTight_hadTau_sublead = isMatched(*selHadTau_sublead, tightHadTausFull);
    }
    if (selHadTaus.size() > 2) {
      prob_fake_hadTau_third = jetToTauFakeRateInterface->getWeight_lead(selHadTau_third->pt(), selHadTau_third->absEta());
    passesTight_hadTau_third = isMatched(*selHadTau_third, tightHadTausFull);
    }
    if (selHadTaus.size() > 3) {
      prob_fake_hadTau_fourth = jetToTauFakeRateInterface->getWeight_lead(selHadTau_fourth->pt(), selHadTau_fourth->absEta());
      passesTight_hadTau_fourth = isMatched(*selHadTau_fourth, tightHadTausFull);
    }

    if (pass_2lss_0tau || pass_2lss_1tau || pass_2los_0tau || pass_ttWctrl) weight_fakeRate = getWeight_2L(
            prob_fake_lepton_lead, passesTight_lepton_lead,
            prob_fake_lepton_sublead, passesTight_lepton_sublead);
    if (pass_2los_1tau) weight_fakeRate = getWeight_3L(
            prob_fake_lepton_lead, passesTight_lepton_lead,
            prob_fake_lepton_sublead, passesTight_lepton_sublead,
            prob_fake_hadTau_lead, passesTight_hadTau_lead
          );
    if (pass_3l_0tau || pass_ttZctrl || pass_WZctrl) weight_fakeRate = getWeight_3L(
            prob_fake_lepton_lead, passesTight_lepton_lead,
            prob_fake_lepton_sublead, passesTight_lepton_sublead,
            prob_fake_lepton_third, passesTight_lepton_third
          );
    if (pass_4l_0tau || pass_ZZctrl) weight_fakeRate = getWeight_4L(
            prob_fake_lepton_lead, passesTight_lepton_lead,
            prob_fake_lepton_sublead, passesTight_lepton_sublead,
            prob_fake_lepton_third, passesTight_lepton_third,
            prob_fake_lepton_fourth, passesTight_lepton_fourth
          );
    if (pass_3l_1tau) weight_fakeRate = getWeight_4L(
            prob_fake_lepton_lead, passesTight_lepton_lead,
            prob_fake_lepton_sublead, passesTight_lepton_sublead,
            prob_fake_lepton_third, passesTight_lepton_third,
            prob_fake_hadTau_lead, passesTight_hadTau_lead
          );
    if (pass_2l_2tau) weight_fakeRate = getWeight_4L(
      prob_fake_lepton_lead, passesTight_lepton_lead,
      prob_fake_lepton_sublead, passesTight_lepton_sublead,
      prob_fake_hadTau_lead, passesTight_hadTau_lead,
      prob_fake_hadTau_sublead, passesTight_hadTau_sublead
    );
    if (pass_1l_1tau) weight_fakeRate = getWeight_2L(
      prob_fake_lepton_lead, passesTight_lepton_lead,
      prob_fake_hadTau_lead, passesTight_hadTau_lead
    );
    if (pass_1l_2tau) weight_fakeRate = getWeight_3L(
      prob_fake_lepton_lead, passesTight_lepton_lead,
      prob_fake_hadTau_lead, passesTight_hadTau_lead,
      prob_fake_hadTau_sublead, passesTight_hadTau_sublead
    );
    if (pass_1l_3tau) weight_fakeRate = getWeight_4L(
      prob_fake_lepton_lead, passesTight_lepton_lead,
      prob_fake_hadTau_lead, passesTight_hadTau_lead,
      prob_fake_hadTau_sublead, passesTight_hadTau_sublead,
      prob_fake_hadTau_third, passesTight_hadTau_third
    );
    if (pass_0l_2tau) weight_fakeRate = getWeight_2L(
      prob_fake_hadTau_lead, passesTight_hadTau_lead,
      prob_fake_hadTau_sublead, passesTight_hadTau_sublead
    );
    if (pass_0l_3tau) weight_fakeRate = getWeight_3L(
      prob_fake_hadTau_lead, passesTight_hadTau_lead,
      prob_fake_hadTau_sublead, passesTight_hadTau_sublead,
      prob_fake_hadTau_third, passesTight_hadTau_third
    );
    if (pass_0l_4tau) weight_fakeRate = getWeight_4L(
      prob_fake_hadTau_lead, passesTight_hadTau_lead,
      prob_fake_hadTau_sublead, passesTight_hadTau_sublead,
      prob_fake_hadTau_third, passesTight_hadTau_third,
      prob_fake_hadTau_fourth, passesTight_hadTau_fourth
    );
    evtWeight *= weight_fakeRate;

    bool resolved_and_semi_AK8 = false;
    bool boosted_and_semi_AK8 = false;
    bool resolved_and_boosted = false;

    double max_mvaOutput_HTT_CSVsort4rd = 0.;
    bool max_truth_HTT_CSVsort4rd = false;
    double HadTop_pt_CSVsort4rd = 0.;
    double HadTop_eta_CSVsort4rd = 0.;
    double genTopPt_CSVsort4rd = 0.;
    double Wj1_pt_CSVsort4rd = -1000;
    double Wj2_pt_CSVsort4rd = -1000;
    double b_pt_CSVsort4rd   = -1000;

    double max_mvaOutput_HTT_CSVsort4rd_2 = 0.;
    bool max_truth_HTT_CSVsort4rd_2 = false;
    double HadTop_pt_CSVsort4rd_2 = 0.;
    double genTopPt_CSVsort4rd_2 = 0.;

    bool hadtruth = false;
    bool hadtruth_2 = false;

    bool calculate_matching = isMC && !applyAdditionalEvtWeight; // DY has not matching info
    std::map<int, Particle::LorentzVector> genVar;
    std::map<int, Particle::LorentzVector> genVarAnti;
    if (calculate_matching) {
      genVar = isGenMatchedJetTripletVar(genTopQuarks, genBJets, genWBosons, genQuarkFromTop, kGenTop);
      genVarAnti = isGenMatchedJetTripletVar(genTopQuarks, genBJets, genWBosons, genQuarkFromTop, kGenAntiTop);
    }

    if ((int)selJets.size() >= 3) {
      for ( std::vector<const RecoJet*>::const_iterator selBJet = selJets.begin(); selBJet != selJets.end(); ++selBJet ) {
        //btag_iterator++;
        for ( std::vector<const RecoJet*>::const_iterator selWJet1 = selBJet + 1; selWJet1 != selJets.end(); ++selWJet1 ) {
         if ( &(*selWJet1) == &(*selBJet) ) continue;
         for ( std::vector<const RecoJet*>::const_iterator selWJet2 = selWJet1 + 1; selWJet2 != selJets.end(); ++selWJet2 ) {
      if ( &(*selWJet2) == &(*selBJet) ) continue;
      if ( &(*selWJet2) == &(*selWJet1) ) continue;
      bool isGenMatched = false;
      double genTopPt_teste = 0.;
      const std::map<int, double> bdtResult = (*hadTopTagger)(**selBJet, **selWJet1, **selWJet2, calculate_matching, isGenMatched, genTopPt_teste, genVar, genVarAnti, true );
      // genTopPt_teste is filled with the result of gen-matching
      if ( isGenMatched ) hadtruth = true;
      // save genpt of all options
      double HadTop_pt = ((*selBJet)->p4() + (*selWJet1)->p4() + (*selWJet2)->p4()).pt();

      if ( bdtResult.at(kXGB_CSVsort4rd) > max_mvaOutput_HTT_CSVsort4rd ) {
        max_truth_HTT_CSVsort4rd = isGenMatched;
        max_mvaOutput_HTT_CSVsort4rd = bdtResult.at(kXGB_CSVsort4rd);
        HadTop_pt_CSVsort4rd = HadTop_pt;
        genTopPt_CSVsort4rd = genTopPt_teste;
        HadTop_eta_CSVsort4rd = std::fabs(((*selBJet)->p4() + (*selWJet1)->p4() + (*selWJet2)->p4()).eta());
        Wj1_pt_CSVsort4rd = (*selWJet1)->pt();
        Wj2_pt_CSVsort4rd = (*selWJet2)->pt();
        b_pt_CSVsort4rd   = (*selBJet)->pt();
      }
      }
        }
      }

      if ( selJets.size() > 5 ){
      for ( std::vector<const RecoJet*>::const_iterator selBJet = selJets.begin(); selBJet != selJets.end(); ++selBJet ) {
        if (b_pt_CSVsort4rd == (*selBJet)->pt() || Wj1_pt_CSVsort4rd == (*selBJet)->pt() || Wj2_pt_CSVsort4rd == (*selBJet)->pt()) continue;
        for ( std::vector<const RecoJet*>::const_iterator selWJet1 = selBJet + 1; selWJet1 != selJets.end(); ++selWJet1 ) {
         if ( &(*selWJet1) == &(*selBJet) ) continue;
         if (b_pt_CSVsort4rd == (*selWJet1)->pt() || Wj1_pt_CSVsort4rd == (*selWJet1)->pt() || Wj2_pt_CSVsort4rd  == (*selWJet1)->pt()) continue;
         for ( std::vector<const RecoJet*>::const_iterator selWJet2 = selWJet1 + 1; selWJet2 != selJets.end(); ++selWJet2 ) {
      if (b_pt_CSVsort4rd == (*selWJet2)->pt() || Wj1_pt_CSVsort4rd == (*selWJet2)->pt() || Wj2_pt_CSVsort4rd  == (*selWJet2)->pt()) continue;
      if ( &(*selWJet2) == &(*selBJet) ) continue;
      if ( &(*selWJet2) == &(*selWJet1) ) continue;
      bool isGenMatched = false;
      double genTopPt_teste = 0.;
      const std::map<int, double> bdtResult = (*hadTopTagger)(**selBJet, **selWJet1, **selWJet2, calculate_matching, isGenMatched, genTopPt_teste, genVar, genVarAnti, true );
      // genTopPt_teste is filled with the result of gen-matching
      if ( isGenMatched ) hadtruth_2 = true;
      // save genpt of all options
      double HadTop_pt = ((*selBJet)->p4() + (*selWJet1)->p4() + (*selWJet2)->p4()).pt();

      if ( bdtResult.at(kXGB_CSVsort4rd) > max_mvaOutput_HTT_CSVsort4rd_2 ) {
        max_truth_HTT_CSVsort4rd_2 = isGenMatched;
        max_mvaOutput_HTT_CSVsort4rd_2 = bdtResult.at(kXGB_CSVsort4rd);
        HadTop_pt_CSVsort4rd_2 = HadTop_pt;
        genTopPt_CSVsort4rd_2 = genTopPt_teste;
      }

      }
          }
        } // close loop on jetS
      } // close if 6 jets
    }

    std::map<std::string, double> mvaInputs_Hj_tagger;
    double mvaOutput_Hj_tagger = 0.;
    for ( std::vector<const RecoJet*>::const_iterator selJet = selJets.begin();
    selJet != selJets.end(); ++selJet ) {
      if ((*selJet)->pt()==Wj1_pt_CSVsort4rd || (*selJet)->pt()==Wj2_pt_CSVsort4rd || (*selJet)->pt()==b_pt_CSVsort4rd) continue;
      double mvaOutput = comp_mvaOutput_Hj_tagger(
        *selJet, fakeableLeptons, mvaInputs_Hj_tagger, mva_Hj_tagger,
        eventInfo);
      if ( mvaOutput > mvaOutput_Hj_tagger ) mvaOutput_Hj_tagger = mvaOutput;
    }

    //--- boosted hTT
    double HTT_boosted = 0.;
    bool bWj1Wj2_isGenMatched_boosted = false;
    double genTopPt_boosted = 0.;
    double HadTop_pt_boosted = 0.;
    bool hadtruth_boosted = false;
    double minDR_HTTv2_lep = -1.;
    double minDR_HTTv2subjets_lep = -1.;

    for ( std::vector<const RecoJetHTTv2*>::const_iterator jetIter = sel_HTTv2.begin();
    jetIter != sel_HTTv2.end(); ++jetIter ) {
    bool isGenMatched = false;
    double genTopPt_teste = 0.;
    const std::map<int, double> bdtResult = (*hadTopTagger_boosted)(**jetIter, calculate_matching, isGenMatched, genTopPt_teste, genVar, genVarAnti);
    if (isGenMatched) {hadtruth_boosted = true;}

    if ( bdtResult.at(kXGB_boosted_no_kinFit) > HTT_boosted ) {
      bWj1Wj2_isGenMatched_boosted = isGenMatched;
      HTT_boosted = bdtResult.at(kXGB_boosted_no_kinFit);
      HadTop_pt_boosted = (*jetIter)->pt();
      genTopPt_boosted = genTopPt_teste;

      if (selHadTau_lead && selHadTau_sublead) minDR_HTTv2_lep = std::min(
        deltaR(selHadTau_lead->p4(), (*jetIter)->p4()),
        deltaR(selHadTau_sublead->p4(), (*jetIter)->p4())
      );

      if (selHadTau_lead && selHadTau_sublead) minDR_HTTv2subjets_lep =
      std::min(
      std::min(
      std::min(
        deltaR(selHadTau_lead->p4(), (*jetIter)->subJet1()->p4()),
        deltaR(selHadTau_sublead->p4(), (*jetIter)->subJet1()->p4())
      ),
      std::min(
        deltaR(selHadTau_lead->p4(), (*jetIter)->subJet2()->p4()),
        deltaR(selHadTau_sublead->p4(), (*jetIter)->subJet2()->p4())
      )
     ),
     std::min(
       deltaR(selHadTau_lead->p4(), (*jetIter)->subJet3()->p4()),
       deltaR(selHadTau_sublead->p4(), (*jetIter)->subJet3()->p4())
     )
     );

    }
    }
    if (genTopPt_CSVsort4rd == genTopPt_boosted) resolved_and_boosted = true;

    // -- semi-boosted hTT -- AK8
    double HTT_semi_boosted_fromAK8 = 0.;
    bool bWj1Wj2_isGenMatched_semi_boosted_fromAK8 = false;
    double genTopPt_semi_boosted_fromAK8 = 0.;
    double HadTop_pt_semi_boosted_fromAK8 = 0.;
    bool hadtruth_semi_boosted_fromAK8 = false;
    double minDR_AK8_lep = -1.;
    double minDR_AK8subjets_lep = -1.;
    for ( std::vector<const RecoJet*>::const_iterator selBJet = cleanedJets_fromAK8.begin(); selBJet != cleanedJets_fromAK8.end(); ++selBJet )  { // cleanedJets.size()
    for ( std::vector<const RecoJetAK8*>::const_iterator jetIter = jet_ptrsAK8.begin();
          jetIter != jet_ptrsAK8.end(); ++jetIter ) {
        bool isGenMatched = false;
        double genTopPt_teste = 0.;
        double HadTop_pt = ((*jetIter)->p4() + (*selBJet)->p4()).pt();
        const std::map<int, double> bdtResult = (*hadTopTagger_semi_boosted_fromAK8)(**jetIter, **selBJet, calculate_matching, isGenMatched, genTopPt_teste, genVar, genVarAnti);
        if (isGenMatched) {hadtruth_semi_boosted_fromAK8 = true;}
        if ( bdtResult.at(kXGB_semi_boosted_AK8_no_kinFit) > HTT_semi_boosted_fromAK8 ) {
          bWj1Wj2_isGenMatched_semi_boosted_fromAK8 = isGenMatched;
          HTT_semi_boosted_fromAK8 = bdtResult.at(kXGB_semi_boosted_AK8_no_kinFit);
          HadTop_pt_semi_boosted_fromAK8 = HadTop_pt;
          genTopPt_semi_boosted_fromAK8 = genTopPt_teste;
          if (selHadTau_lead && selHadTau_sublead) minDR_AK8_lep = std::min(
            deltaR(selHadTau_lead->p4(), (*jetIter)->p4()),
            deltaR(selHadTau_sublead->p4(), (*jetIter)->p4())
          );
          if (selHadTau_lead && selHadTau_sublead) minDR_AK8subjets_lep =
          std::min(
          std::min(
            deltaR(selHadTau_lead->p4(), (*jetIter)->subJet1()->p4()),
            deltaR(selHadTau_sublead->p4(), (*jetIter)->subJet1()->p4())
          ),
          std::min(
            deltaR(selHadTau_lead->p4(), (*jetIter)->subJet2()->p4()),
            deltaR(selHadTau_sublead->p4(), (*jetIter)->subJet2()->p4())
          )
         );
        }
      }
    }
    if (genTopPt_CSVsort4rd == genTopPt_semi_boosted_fromAK8)  resolved_and_semi_AK8 = true;
    if (genTopPt_semi_boosted_fromAK8 == genTopPt_boosted)  boosted_and_semi_AK8 = true;

    if ( selEventsFile ) {
      (*selEventsFile) << eventInfo.run << ':' << eventInfo.lumi << ':' << eventInfo.event << '\n';
    }
    //double test1 = selLepton_lead ? comp_mindr_lep1_jet(*selLepton_lead, selJets) : -1000;
    //double test2 = selLepton_sublead ? comp_mindr_lep1_jet(*selLepton_sublead, selJets) : -1000;
    //if (
    //  std::abs(test1) < 0.3 || std::abs(test2) < 0.3
    //) std::cout<<"Isolation failed Xanda "<<test1 <<" "<<test2<<"\n";

    if ( bdt_filler ) {

      bdt_filler -> operator()({ eventInfo.run, eventInfo.lumi, eventInfo.event })
          ("avg_dr_jet",           comp_avg_dr_jet(selJets))
          ("ptmiss",               met.pt())
          ("htmiss",               mht_p4.pt())
          ("mTauTauVis",           mTauTauVis)
          ("mTauTau_SVFit",        mTauTau)
          ("cosThetaS_hadTau",     cosThetaS_hadTau)
          ("mTauTauVis2",           mTauTauVis2)
          ("mTauTau2_SVFit",        mTauTau2)
          ("cosThetaS_hadTau2",     cosThetaS_hadTau2)
          ("mbb_loose",            selBJets_loose.size() > 1  ? (selBJets_loose[0]->p4()  + selBJets_loose[1]->p4() ).mass() : -1000)
          ("mbb_medium",           selBJets_medium.size() > 1 ? (selBJets_medium[0]->p4() + selBJets_medium[1]->p4()).mass() : -1000)
          ("evtWeight",            evtWeight)
          // the next weights are already embebed on evtWeight, they are here just for debub purposes
          ("weight_fakeRate", weight_fakeRate)
          ("weight_data_to_MC_correction", weight_data_to_MC_correction)
          ("weight_data_to_MC_correction_hadTau", weight_data_to_MC_correction_hadTau)

          ("nJet",                 selJets.size())
          ("nBJetLoose",           selBJets_loose.size())
          ("nBJetMedium",          selBJets_medium.size())
          ("nElectron",                      selElectrons.size())
          ("nLepton",                      selLeptons.size())
          ("nTau", selHadTaus.size())
          //
          ("selTrigger_1e", selTrigger_1e)
          ("selTrigger_2e", selTrigger_2e)
          ("selTrigger_3e", selTrigger_3e)
          ("selTrigger_1mu", selTrigger_1mu)
          ("selTrigger_2mu", selTrigger_2mu)
          ("selTrigger_3mu", selTrigger_3mu)
          ("selTrigger_1e1mu", selTrigger_1e1mu)
          ("selTrigger_2e1mu", selTrigger_2e1mu)
          ("selTrigger_1e2mu", selTrigger_1e2mu)
          ("selTrigger_1e1tau", selTrigger_1e1tau)
          ("selTrigger_1mu1tau", selTrigger_1mu1tau)
          ("selTrigger_2tau", selTrigger_2tau)
          //
          ("pass_2lss_0tau", pass_2lss_0tau)
          ("pass_3l_0tau", pass_3l_0tau)
          ("pass_4l_0tau", pass_4l_0tau)
          ("pass_ttZctrl", pass_ttZctrl)
          ("pass_ttWctrl", pass_ttWctrl)
          ("pass_WZctrl", pass_WZctrl)
          ("pass_ZZctrl", pass_ZZctrl)
          ("pass_1l_2tau", pass_1l_2tau)
          ("pass_2lss_1tau", pass_2lss_1tau)
          ("pass_2los_1tau", pass_2los_1tau)
          ("pass_3l_1tau", pass_3l_1tau)
          ("pass_2l_2tau", pass_2l_2tau)
          ("pass_1l_1tau", pass_1l_1tau)
          ("pass_1l_3tau", pass_1l_3tau)
          ("pass_0l_4tau", pass_0l_4tau)
          ("pass_2los_0tau", pass_2los_0tau)
          ("pass_0l_3tau", pass_0l_3tau)
          ("pass_0l_2tau", pass_0l_2tau)
          ("pass_ttH", pass_ttH)
          ("isCharge_hadTau_OS", isCharge_hadTau_OS)
          ("isCharge_Lepton_OS", isCharge_Lepton_OS)
          ("jet1_pt",   selJet_lead ? selJet_lead->pt() : -1000)
          ("jet1_eta",  selJet_lead ? selJet_lead->eta() : -1000)
          ("jet1_phi",  selJet_lead ? selJet_lead->phi() : -1000)
          ("jet1_mass", selJet_lead ? selJet_lead->mass() : -1000)
          ("jet2_pt",   selJet_sublead ? selJet_sublead->pt() : -1000)
          ("jet2_eta",  selJet_sublead ? selJet_sublead->eta() : -1000)
          ("jet2_phi",  selJet_sublead ? selJet_sublead->phi() : -1000)
          ("jet2_mass", selJet_sublead ? selJet_sublead->mass() : -1000)
          ("jet3_pt",   selJet_third ? selJet_third->pt() : -1000)
          ("jet3_eta",  selJet_third ? selJet_third->eta() : -1000)
          ("jet3_phi",  selJet_third ? selJet_third->phi() : -1000)
          ("jet3_mass", selJet_third ? selJet_third->mass() : -1000)
          ("jet4_pt",   selJet_fourth ? selJet_fourth->pt() : -1000)
          ("jet4_eta",  selJet_fourth ? selJet_fourth->eta() : -1000)
          ("jet4_phi",  selJet_fourth ? selJet_fourth->phi() : -1000)
          ("jet4_mass", selJet_fourth ? selJet_fourth->mass() : -1000)
          //
          ("lep1_conept",     selLepton_lead ? comp_lep1_conePt(*selLepton_lead) : -1000)
          ("lep1_mT",         selLepton_lead ? comp_MT_met_lep1(*selLepton_lead,    met.pt(), met.phi()) : -1000)
          ("lep1_min_dr_jet", selLepton_lead ? comp_mindr_lep1_jet(*selLepton_lead, selJets) : -1000)
          ("lep1_pt",         selLepton_lead ? selLepton_lead->pt() : -1000)
          ("lep1_eta",        selLepton_lead ? selLepton_lead->eta() : -1000)
          ("lep1_phi",        selLepton_lead ? selLepton_lead->phi() : -1000)
          ("lep1_mass",       selLepton_lead ? selLepton_lead->mass() : -1000)
          ("lep2_conept",     selLepton_sublead ? comp_lep1_conePt(*selLepton_sublead) : -1000)
          ("lep2_mT",         selLepton_sublead ? comp_MT_met_lep1(*selLepton_sublead,    met.pt(), met.phi()) : -1000)
          ("lep2_min_dr_jet", selLepton_sublead ? comp_mindr_lep1_jet(*selLepton_sublead, selJets) : -1000)
          ("lep2_pt",         selLepton_sublead ? selLepton_sublead->pt() : -1000)
          ("lep2_eta",        selLepton_sublead ? selLepton_sublead->eta() : -1000)
          ("lep2_phi",        selLepton_sublead ? selLepton_sublead->phi() : -1000)
          ("lep2_mass",       selLepton_sublead ? selLepton_sublead->mass() : -1000)
          ("lep3_conept",     selLepton_third ? comp_lep1_conePt(*selLepton_third) : -1000)
          ("lep3_mT",         selLepton_third ? comp_MT_met_lep1(*selLepton_third,    met.pt(), met.phi()) : -1000)
          ("lep3_min_dr_jet", selLepton_third ? comp_mindr_lep1_jet(*selLepton_third, selJets) : -1000)
          ("lep3_pt",         selLepton_third ? selLepton_third->pt() : -1000)
          ("lep3_eta",        selLepton_third ? selLepton_third->eta() : -1000)
          ("lep3_phi",        selLepton_third ? selLepton_third->phi() : -1000)
          ("lep3_mass",       selLepton_third ? selLepton_third->mass() : -1000)
          ("lep4_conept",     selLepton_fourth ? comp_lep1_conePt(*selLepton_fourth) : -1000)
          ("lep4_mT",         selLepton_fourth ? comp_MT_met_lep1(*selLepton_fourth,    met.pt(), met.phi()) : -1000)
          ("lep4_min_dr_jet", selLepton_fourth ? comp_mindr_lep1_jet(*selLepton_fourth, selJets) : -1000)
          ("lep4_pt",         selLepton_fourth ? selLepton_fourth->pt() : -1000)
          ("lep4_eta",        selLepton_fourth ? selLepton_fourth->eta() : -1000)
          ("lep4_phi",        selLepton_fourth ? selLepton_fourth->phi() : -1000)
          ("lep4_mass",       selLepton_fourth ? selLepton_fourth->mass() : -1000)
          //
          ("tau1_min_dr_jet", selHadTau_lead ? comp_mindr_hadTau1_jet(*selHadTau_lead,    selJets) : -1000)
          //("tau1_min_dr_lep", selHadTau_lead ? comp_mindr_hadTau1_jet(*selHadTau_lead,    selLeptons) : -1000)
          ("tau1_pt",         selHadTau_lead ? selHadTau_lead->pt() : -1000)
          ("tau1_eta",        selHadTau_lead ? selHadTau_lead->eta() : -1000)
          ("tau1_phi",        selHadTau_lead ? selHadTau_lead->phi() : -1000)
          ("tau1_mass",       selHadTau_lead ? selHadTau_lead->mass() : -1000)
          ("tau2_min_dr_jet", selHadTau_sublead ? comp_mindr_hadTau1_jet(*selHadTau_sublead,    selJets) : -1000)
          //("tau2_min_dr_lep", selHadTau_sublead ? comp_mindr_hadTau1_jet(*selHadTau_sublead,    selLeptons) : -1000)
          ("tau2_pt",         selHadTau_sublead ? selHadTau_sublead->pt() : -1000)
          ("tau2_eta",        selHadTau_sublead ? selHadTau_sublead->eta() : -1000)
          ("tau2_phi",        selHadTau_sublead ? selHadTau_sublead->phi() : -1000)
          ("tau2_mass",       selHadTau_sublead ? selHadTau_sublead->mass() : -1000)
          ("tau3_min_dr_jet", selHadTau_third ? comp_mindr_hadTau1_jet(*selHadTau_third,    selJets) : -1000)
          //("tau3_min_dr_lep", selHadTau_third ? comp_mindr_hadTau1_jet(*selHadTau_third,    selLeptons) : -1000)
          ("tau3_pt",         selHadTau_third ? selHadTau_third->pt() : -1000)
          ("tau3_eta",        selHadTau_third ? selHadTau_third->eta() : -1000)
          ("tau3_phi",        selHadTau_third ? selHadTau_third->phi() : -1000)
          ("tau3_mass",       selHadTau_third ? selHadTau_third->mass() : -1000)
          ("tau4_min_dr_jet", selHadTau_fourth ? comp_mindr_hadTau1_jet(*selHadTau_fourth,    selJets) : -1000)
          //("tau4_min_dr_lep", selHadTau_fourth ? comp_mindr_hadTau1_jet(*selHadTau_fourth,    selLeptons) : -1000)
          ("tau4_pt",         selHadTau_fourth ? selHadTau_fourth->pt() : -1000)
          ("tau4_eta",        selHadTau_fourth ? selHadTau_fourth->eta() : -1000)
          ("tau4_phi",        selHadTau_fourth ? selHadTau_fourth->phi() : -1000)
          ("tau4_mass",       selHadTau_fourth ? selHadTau_fourth->mass() : -1000)
          //
          ("sum_lep_charge", sum_lep_charge)
          ("sum_tau_charge", sum_tau_charge)
          ("max_dr_jet",     comp_max_dr_jet(selJets))
          ("met_LD",         met_LD)
          //
          ("lep_match_jet", lep_match_jet)
          ("lep_match_lep", lep_match_lep)
          ("tau_match_jet", tau_match_jet)
          ("tau_match_tau", tau_match_tau)

          ("bWj1Wj2_isGenMatched_CSVsort4rd",              max_truth_HTT_CSVsort4rd)
          ("bWj1Wj2_isGenMatched_CSVsort4rd_2",              max_truth_HTT_CSVsort4rd_2)
          ("res-HTT_CSVsort4rd",                 max_mvaOutput_HTT_CSVsort4rd)
          ("res-HTT_CSVsort4rd_2",                 max_mvaOutput_HTT_CSVsort4rd_2)
          ("HadTop_pt_CSVsort4rd",            HadTop_pt_CSVsort4rd)
          ("HadTop_pt_CSVsort4rd_2",            HadTop_pt_CSVsort4rd_2)
          ("HadTop_eta_CSVsort4rd",            HadTop_eta_CSVsort4rd)
          ("genTopPt_CSVsort4rd",             genTopPt_CSVsort4rd)
          ("genTopPt_CSVsort4rd_2",             genTopPt_CSVsort4rd_2)
          ("hadtruth",                  hadtruth)
          ("hadtruth_2",                hadtruth_2)
          ("hadtruth_boosted",          hadtruth_boosted)
          ("hadtruth_semi_boosted_fromAK8",          hadtruth_semi_boosted_fromAK8)
          ("nHTTv2",                         sel_HTTv2.size())
          ("HTTv2_lead_pt",                  sel_HTTv2.size() > 0 ? sel_HTTv2[0]->pt() : -1 )
          ("minDR_HTTv2_lep",                minDR_HTTv2_lep)
          ("minDR_HTTv2subjets_lep",         minDR_HTTv2subjets_lep)
          ("HTT_boosted",                     HTT_boosted)
          ("bWj1Wj2_isGenMatched_boosted",    bWj1Wj2_isGenMatched_boosted)
          ("genTopPt_boosted",            genTopPt_boosted)
          ("HadTop_pt_boosted",           HadTop_pt_boosted)
          ("HTT_semi_boosted_fromAK8",                     HTT_semi_boosted_fromAK8)
          ("bWj1Wj2_isGenMatched_semi_boosted_fromAK8",    bWj1Wj2_isGenMatched_semi_boosted_fromAK8)
          ("genTopPt_semi_boosted_fromAK8",            genTopPt_semi_boosted_fromAK8)
          ("HadTop_pt_semi_boosted_fromAK8",           HadTop_pt_semi_boosted_fromAK8)
          ("N_jetAK8",     jet_ptrsAK8.size())
          ("minDR_AK8_lep",                minDR_AK8_lep)
          ("minDR_AK8subjets_lep",         minDR_AK8subjets_lep)
          ("resolved_and_semi_AK8",     resolved_and_semi_AK8)
          ("boosted_and_semi_AK8",      boosted_and_semi_AK8)
          ("resolved_and_boosted",      resolved_and_boosted)
          ("failsZbosonMassVeto",       failsZbosonMassVeto)
          ("failsLowMassVeto", failsLowMassVeto)
          ("failsHtoZZVeto", failsHtoZZVeto)

          ("nJetForward", selJetsForward.size())
          ("jetForward1_pt",   selJetsForward.size() > 0 ? selJetsForward[0]->pt() : -1000)
          ("jetForward1_eta",  selJetsForward.size() > 0 ? selJetsForward[0]->eta() : -1000)
          ("jetForward1_phi",  selJetsForward.size() > 0 ? selJetsForward[0]->phi() : -1000)
          ("jetForward1_mass", selJetsForward.size() > 0 ? selJetsForward[0]->mass() : -1000)
          ("jetForward2_pt",   selJetsForward.size() > 1 ? selJetsForward[1]->pt() : -1000)
          ("jetForward2_eta",  selJetsForward.size() > 1 ? selJetsForward[1]->eta() : -1000)
          ("jetForward2_phi",  selJetsForward.size() > 1 ? selJetsForward[1]->phi() : -1000)
          ("jetForward2_mass", selJetsForward.size() > 1 ? selJetsForward[1]->mass() : -1000)
          ("DijetForward_mass", selJetsForward.size() > 1 ? (selJetsForward[0]->p4() + selJetsForward[1]->p4()).mass() : -1000)

          ("mvaOutput_Hj_tagger", mvaOutput_Hj_tagger)

        .fill()
      ;
    }

    ++selectedEntries;
    selectedEntries_weighted += evtWeight;
    histogram_selectedEntries->Fill(0.);
  }

  std::cout << "max num. Entries = " << inputTree -> getCumulativeMaxEventCount()
            << " (limited by " << maxEvents << ") processed in "
            << inputTree -> getProcessedFileCount() << " file(s) (out of "
            << inputTree -> getFileCount() << ")\n"
            << " analyzed = " << analyzedEntries << '\n'
            << " selected = " << selectedEntries << " (weighted = " << selectedEntries_weighted << ")\n\n"
            << "cut-flow table" << std::endl;
  std::cout << std::endl;

  std::cout << "sel. Entries by gen. matching:" << std::endl;
  for ( std::vector<leptonGenMatchEntry>::const_iterator leptonGenMatch_definition = leptonGenMatch_definitions.begin();
	leptonGenMatch_definition != leptonGenMatch_definitions.end(); ++leptonGenMatch_definition ) {
    for ( std::vector<hadTauGenMatchEntry>::const_iterator hadTauGenMatch_definition = hadTauGenMatch_definitions.begin();
	  hadTauGenMatch_definition != hadTauGenMatch_definitions.end(); ++hadTauGenMatch_definition ) {

      std::string process_and_genMatch = process_string;
      if ( apply_leptonGenMatching ) process_and_genMatch += leptonGenMatch_definition->name_;
      if ( apply_leptonGenMatching && apply_hadTauGenMatching ) process_and_genMatch += "&";
      if ( apply_hadTauGenMatching ) process_and_genMatch += hadTauGenMatch_definition->name_;

    }
  }
  std::cout << std::endl;

  delete dataToMCcorrectionInterface;

  delete leptonFakeRateInterface;
  delete jetToTauFakeRateInterface;

  delete run_lumi_eventSelector;

  delete selEventsFile;

  delete muonReader;
  delete electronReader;
  delete hadTauReader;
  delete jetReader;
  delete metReader;
  delete metFilterReader;
  delete genLeptonReader;
  delete genHadTauReader;
  delete genPhotonReader;
  delete genJetReader;
  delete lheInfoReader;

  delete eventWeightManager;

  hltPaths_delete(triggers_1e);
  hltPaths_delete(triggers_2e);
  hltPaths_delete(triggers_1mu);
  hltPaths_delete(triggers_2mu);
  hltPaths_delete(triggers_1e1mu);

  delete inputTree;

  clock.Show("analyze_inclusive_toNN");

  return EXIT_SUCCESS;
}
