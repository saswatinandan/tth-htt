#include "FWCore/ParameterSet/interface/ParameterSet.h" // edm::ParameterSet
#include "FWCore/Utilities/interface/Exception.h" // cms::Exception
#include "PhysicsTools/FWLite/interface/TFileService.h" // fwlite::TFileService
#include "DataFormats/FWLite/interface/InputSource.h" // fwlite::InputSource
#include "DataFormats/FWLite/interface/OutputFiles.h" // fwlite::OutputFiles
#include "DataFormats/Math/interface/LorentzVector.h" // math::PtEtaPhiMLorentzVector
#include "DataFormats/Math/interface/deltaR.h" // deltaR

#if __has_include (<FWCore/ParameterSetReader/interface/ParameterSetReader.h>)
#  include <FWCore/ParameterSetReader/interface/ParameterSetReader.h> // edm::readPSetsFrom()
#else
#  include <FWCore/PythonParameterSet/interface/MakeParameterSets.h> // edm::readPSetsFrom()
#endif

#include <TBenchmark.h> // TBenchmark
#include <TString.h> // TString, Form
#include <TMatrixD.h> // TMatrixD
#include <TError.h> // gErrorAbortLevel, kError
#include <TMath.h> // TMath::
#include <TH2.h> // TH2
#include <TROOT.h> // TROOT

#include "tthAnalysis/HiggsToTauTau/interface/RecoLepton.h" // RecoLepton
#include "tthAnalysis/HiggsToTauTau/interface/RecoJet.h" // RecoJet
#include "tthAnalysis/HiggsToTauTau/interface/RecoHadTau.h" // RecoHadTau
#include "tthAnalysis/HiggsToTauTau/interface/GenLepton.h" // GenLepton
#include "tthAnalysis/HiggsToTauTau/interface/GenJet.h" // GenJet
#include "tthAnalysis/HiggsToTauTau/interface/GenHadTau.h" // GenHadTau
#include "tthAnalysis/HiggsToTauTau/interface/RecoMEt.h" // RecoMEt
#include "tthAnalysis/HiggsToTauTau/interface/EventInfo.h" // EventInfo
#include "tthAnalysis/HiggsToTauTau/interface/ObjectMultiplicity.h" // ObjectMultiplicity
#include "tthAnalysis/HiggsToTauTau/interface/TMVAInterface.h" // TMVAInterface
#include "tthAnalysis/HiggsToTauTau/interface/mvaAuxFunctions.h" // check_mvaInputs, get_mvaInputVariables
#include "tthAnalysis/HiggsToTauTau/interface/mvaAuxFunctions_Hj_and_Hjj_taggers.h" // comp_mvaOutput_Hj_tagger, comp_mvaOutput_Hjj_tagger
#include "tthAnalysis/HiggsToTauTau/interface/mvaInputVariables.h" // auxiliary functions for computing input variables of the MVA used for signal extraction in the 2lss category
#include "tthAnalysis/HiggsToTauTau/interface/LeptonFakeRateInterface.h" // LeptonFakeRateInterface
#include "tthAnalysis/HiggsToTauTau/interface/JetToTauFakeRateInterface.h" // JetToTauFakeRateInterface
#include "tthAnalysis/HiggsToTauTau/interface/RecoElectronReader.h" // RecoElectronReader
#include "tthAnalysis/HiggsToTauTau/interface/RecoMuonReader.h" // RecoMuonReader
#include "tthAnalysis/HiggsToTauTau/interface/RecoHadTauReader.h" // RecoHadTauReader
#include "tthAnalysis/HiggsToTauTau/interface/RecoJetReader.h" // RecoJetReader
#include "tthAnalysis/HiggsToTauTau/interface/RecoMEtReader.h" // RecoMEtReader
#include "tthAnalysis/HiggsToTauTau/interface/MEtFilterReader.h" // MEtFilterReader
#include "tthAnalysis/HiggsToTauTau/interface/GenLeptonReader.h" // GenLeptonReader
#include "tthAnalysis/HiggsToTauTau/interface/GenHadTauReader.h" // GenHadTauReader
#include "tthAnalysis/HiggsToTauTau/interface/GenPhotonReader.h" // GenPhotonReader
#include "tthAnalysis/HiggsToTauTau/interface/GenJetReader.h" // GenJetReader
#include "tthAnalysis/HiggsToTauTau/interface/GenParticleReader.h" // GenParticleReader
#include "tthAnalysis/HiggsToTauTau/interface/LHEInfoReader.h" // LHEInfoReader
#include "tthAnalysis/HiggsToTauTau/interface/PSWeightReader.h" // PSWeightReader
#include "tthAnalysis/HiggsToTauTau/interface/EventInfoReader.h" // EventInfoReader
#include "tthAnalysis/HiggsToTauTau/interface/ObjectMultiplicityReader.h" // ObjectMultiplicityReader
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
#include "tthAnalysis/HiggsToTauTau/interface/RecoJetCollectionSelector.h" // RecoJetCollectionSelector
#include "tthAnalysis/HiggsToTauTau/interface/RecoJetCollectionSelectorBtag.h" // RecoJetCollectionSelectorBtagLoose, RecoJetCollectionSelectorBtagMedium
#include "tthAnalysis/HiggsToTauTau/interface/RecoJetCollectionSelectorForward.h" // RecoJetSelectorForward
#include "tthAnalysis/HiggsToTauTau/interface/RunLumiEventSelector.h" // RunLumiEventSelector
#include "tthAnalysis/HiggsToTauTau/interface/MEtFilterSelector.h" // MEtFilterSelector
#include "tthAnalysis/HiggsToTauTau/interface/ElectronHistManager.h" // ElectronHistManager
#include "tthAnalysis/HiggsToTauTau/interface/MuonHistManager.h" // MuonHistManager
#include "tthAnalysis/HiggsToTauTau/interface/HadTauHistManager.h" // HadTauHistManager
#include "tthAnalysis/HiggsToTauTau/interface/JetHistManager.h" // JetHistManager
#include "tthAnalysis/HiggsToTauTau/interface/MEtHistManager.h" // MEtHistManager
#include "tthAnalysis/HiggsToTauTau/interface/MEtFilterHistManager.h" // MEtFilterHistManager
#include "tthAnalysis/HiggsToTauTau/interface/MVAInputVarHistManager.h" // MVAInputVarHistManager
#include "tthAnalysis/HiggsToTauTau/interface/EvtHistManager_2lss.h" // EvtHistManager_2lss
#include "tthAnalysis/HiggsToTauTau/interface/EvtYieldHistManager.h" // EvtYieldHistManager
#include "tthAnalysis/HiggsToTauTau/interface/CutFlowTableHistManager.h" // CutFlowTableHistManager
#include "tthAnalysis/HiggsToTauTau/interface/WeightHistManager.h" // WeightHistManager
#include "tthAnalysis/HiggsToTauTau/interface/GenEvtHistManager.h" // GenEvtHistManager
#include "tthAnalysis/HiggsToTauTau/interface/LHEInfoHistManager.h" // LHEInfoHistManager
#include "tthAnalysis/HiggsToTauTau/interface/leptonTypes.h" // getLeptonType, kElectron, kMuon
#include "tthAnalysis/HiggsToTauTau/interface/analysisAuxFunctions.h" // getBTagWeight_option, getHadTau_genPdgId, isHigherPt, isMatched
#include "tthAnalysis/HiggsToTauTau/interface/HistManagerBase.h" // HistManagerBase
#include "tthAnalysis/HiggsToTauTau/interface/histogramAuxFunctions.h" // fillWithOverFlow2d
#include "tthAnalysis/HiggsToTauTau/interface/hadTopTaggerAuxFunctions.h" // isGenMatchedJetTriplet
#include "tthAnalysis/HiggsToTauTau/interface/hadTopTaggerAuxFunctions_geral.h" // isGenMatchedJetTriplet tags
#include "tthAnalysis/HiggsToTauTau/interface/leptonGenMatchingAuxFunctions.h" // getLeptonGenMatch_definitions_1lepton, getLeptonGenMatch_string, getLeptonGenMatch_int
#include "tthAnalysis/HiggsToTauTau/interface/GenMatchInterface.h" // GenMatchInterface
#include "tthAnalysis/HiggsToTauTau/interface/fakeBackgroundAuxFunctions.h" // getWeight_2L
#include "tthAnalysis/HiggsToTauTau/interface/ChargeMisIdRate.h" // ChargeMisIdRate
#include "tthAnalysis/HiggsToTauTau/interface/hltPath.h" // hltPath, create_hltPaths, hltPaths_isTriggered, hltPaths_delete
#include "tthAnalysis/HiggsToTauTau/interface/hltPathReader.h" // hltPathReader
#include "tthAnalysis/HiggsToTauTau/interface/Data_to_MC_CorrectionInterface_2016.h"
#include "tthAnalysis/HiggsToTauTau/interface/Data_to_MC_CorrectionInterface_2017.h"
#include "tthAnalysis/HiggsToTauTau/interface/Data_to_MC_CorrectionInterface_2018.h"
#include "tthAnalysis/HiggsToTauTau/interface/lutAuxFunctions.h" // loadTH2, getSF_from_TH2
#include "tthAnalysis/HiggsToTauTau/interface/L1PreFiringWeightReader.h" // L1PreFiringWeightReader
#include "tthAnalysis/HiggsToTauTau/interface/cutFlowTable.h" // cutFlowTableType
#include "tthAnalysis/HiggsToTauTau/interface/NtupleFillerBDT.h" // NtupleFillerBDT
#include "tthAnalysis/HiggsToTauTau/interface/HadTopTagger.h" // HadTopTagger
#include "tthAnalysis/HiggsToTauTau/interface/XGBInterface.h" // XGBInterface
#include "tthAnalysis/HiggsToTauTau/interface/HadTopKinFit.h" // HadTopKinFit
#include "tthAnalysis/HiggsToTauTau/interface/TTreeWrapper.h" // TTreeWrapper
#include "tthAnalysis/HiggsToTauTau/interface/SyncNtupleManager.h" // SyncNtupleManager
#include "tthAnalysis/HiggsToTauTau/interface/hltFilter.h" // hltFilter()
#include "tthAnalysis/HiggsToTauTau/interface/EvtWeightManager.h" // EvtWeightManager
#include "tthAnalysis/HiggsToTauTau/interface/TensorFlowInterface.h"
#include "tthAnalysis/HiggsToTauTau/interface/mT2_2particle.h" // mT2_2particle
#include "tthAnalysis/HiggsToTauTau/interface/mT2_3particle.h" // mT2_3particle
#include "tthAnalysis/HiggsToTauTau/interface/weightAuxFunctions.h" // get_tH_weight_str()
#include "tthAnalysis/HiggsToTauTau/interface/EvtWeightRecorder.h" // EvtWeightRecorder
#include "tthAnalysis/HiggsToTauTau/interface/BtagSFRatioFacility.h" // BtagSFRatioFacility
#include "tthAnalysis/HiggsToTauTau/interface/HHWeightInterfaceLO.h" // HHWeightInterfaceLO
#include "tthAnalysis/HiggsToTauTau/interface/HHWeightInterfaceNLO.h" // HHWeightInterfaceNLO
#include "tthAnalysis/HiggsToTauTau/interface/HHWeightInterfaceCouplings.h" // HHWeightInterfaceCouplings
#include "tthAnalysis/HiggsToTauTau/interface/AnalysisConfig.h" // AnalysisConfig
#include "tthAnalysis/HiggsToTauTau/interface/RecoVertex.h" // RecoVertex
#include "tthAnalysis/HiggsToTauTau/interface/RecoVertexReader.h" // RecoVertexReader

#include <boost/algorithm/string/replace.hpp> // boost::replace_all_copy()
#include <boost/algorithm/string/predicate.hpp> // boost::starts_with()

#include <iostream> // std::cerr, std::fixed
#include <iomanip> // std::setprecision(), std::setw()
#include <string> // std::string
#include <vector> // std::vector<>
#include <cstdlib> // EXIT_SUCCESS, EXIT_FAILURE
#include <fstream> // std::ofstream
#include <assert.h> // assert

typedef math::PtEtaPhiMLorentzVector LV;
typedef std::vector<std::string> vstring;

enum { kFR_disabled, kFR_2lepton };

/**
 * @brief Produce datacard and control plots for 2lss categories.
 */
int main(int argc, char* argv[])
{
//--- throw an exception in case ROOT encounters an error
  gErrorAbortLevel = kError;

//--- stop ROOT from keeping track of all histograms
  TH1::AddDirectory(false);

//--- parse command-line arguments
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " [parameters.py]" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "<analyze_2lss>:" << std::endl;

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("analyze_2lss");

//--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") )
    throw cms::Exception("analyze_2lss")
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::ParameterSet cfg_analyze = cfg.getParameter<edm::ParameterSet>("analyze_2lss");
  AnalysisConfig analysisConfig("ttH->multilepton+tau", cfg_analyze);

  std::string treeName = cfg_analyze.getParameter<std::string>("treeName");

  std::string process_string = cfg_analyze.getParameter<std::string>("process");
  const bool isMC_tH = analysisConfig.isMC_tH();
  const bool isMC_VH = analysisConfig.isMC_VH();
  const bool isMC_H  = analysisConfig.isMC_H();
  const bool isMC_HH = analysisConfig.isMC_HH();
  const bool isMC_EWK = analysisConfig.isMC_EWK();
  const bool isMC_signal = analysisConfig.isMC_ttH();
  const bool isSignal = isMC_signal || isMC_tH || isMC_VH || isMC_HH || isMC_H;
  std::string histogramDir = cfg_analyze.getParameter<std::string>("histogramDir");
  bool isMCClosure_e = histogramDir.find("mcClosure_e") != std::string::npos;
  bool isMCClosure_m = histogramDir.find("mcClosure_m") != std::string::npos;

  std::string era_string = cfg_analyze.getParameter<std::string>("era");
  const Era era = get_era(era_string);
  const unsigned int skipEvery = cfg_analyze.getParameter<unsigned int>("skipEvery");

  vstring triggerNames_1e = cfg_analyze.getParameter<vstring>("triggers_1e");
  std::vector<hltPath*> triggers_1e = create_hltPaths(triggerNames_1e, "triggers_1e");
  bool use_triggers_1e = cfg_analyze.getParameter<bool>("use_triggers_1e");
  vstring triggerNames_2e = cfg_analyze.getParameter<vstring>("triggers_2e");
  std::vector<hltPath*> triggers_2e = create_hltPaths(triggerNames_2e, "triggers_2e");
  bool use_triggers_2e = cfg_analyze.getParameter<bool>("use_triggers_2e");
  vstring triggerNames_1mu = cfg_analyze.getParameter<vstring>("triggers_1mu");
  std::vector<hltPath*> triggers_1mu = create_hltPaths(triggerNames_1mu, "triggers_1mu");
  bool use_triggers_1mu = cfg_analyze.getParameter<bool>("use_triggers_1mu");
  vstring triggerNames_2mu = cfg_analyze.getParameter<vstring>("triggers_2mu");
  std::vector<hltPath*> triggers_2mu = create_hltPaths(triggerNames_2mu, "triggers_2mu");
  bool use_triggers_2mu = cfg_analyze.getParameter<bool>("use_triggers_2mu");
  vstring triggerNames_1e1mu = cfg_analyze.getParameter<vstring>("triggers_1e1mu");
  std::vector<hltPath*> triggers_1e1mu = create_hltPaths(triggerNames_1e1mu, "triggers_1e1mu");
  bool use_triggers_1e1mu = cfg_analyze.getParameter<bool>("use_triggers_1e1mu");

  bool apply_offline_e_trigger_cuts_1e = cfg_analyze.getParameter<bool>("apply_offline_e_trigger_cuts_1e");
  bool apply_offline_e_trigger_cuts_2e = cfg_analyze.getParameter<bool>("apply_offline_e_trigger_cuts_2e");
  bool apply_offline_e_trigger_cuts_1mu = cfg_analyze.getParameter<bool>("apply_offline_e_trigger_cuts_1mu");
  bool apply_offline_e_trigger_cuts_2mu = cfg_analyze.getParameter<bool>("apply_offline_e_trigger_cuts_2mu");
  bool apply_offline_e_trigger_cuts_1e1mu = cfg_analyze.getParameter<bool>("apply_offline_e_trigger_cuts_1e1mu");

  enum { kOS, kSS };
  std::string leptonChargeSelection_string = cfg_analyze.getParameter<std::string>("leptonChargeSelection");
  int leptonChargeSelection = -1;
  if      ( leptonChargeSelection_string == "OS" ) leptonChargeSelection = kOS;
  else if ( leptonChargeSelection_string == "SS" ) leptonChargeSelection = kSS;
  else throw cms::Exception("analyze_2lss")
    << "Invalid Configuration parameter 'leptonChargeSelection' = " << leptonChargeSelection_string << " !!\n";

  const std::string electronSelection_string = cfg_analyze.getParameter<std::string>("electronSelection");
  const std::string muonSelection_string     = cfg_analyze.getParameter<std::string>("muonSelection");
  std::cout << "electronSelection_string = " << electronSelection_string << "\n"
               "muonSelection_string     = " << muonSelection_string     << '\n'
  ;
  const int electronSelection = get_selection(electronSelection_string);
  const int muonSelection     = get_selection(muonSelection_string);

  bool apply_leptonGenMatching = cfg_analyze.getParameter<bool>("apply_leptonGenMatching");
  std::vector<leptonChargeFlipGenMatchEntry> leptonGenMatch_definitions = getLeptonChargeFlipGenMatch_definitions_2lepton(true);
  std::cout << "leptonGenMatch_definitions:" << std::endl;
  std::cout << leptonGenMatch_definitions;

  GenMatchInterface genMatchInterface(2, apply_leptonGenMatching, true);

  TString hadTauSelection_string = cfg_analyze.getParameter<std::string>("hadTauSelection").data();
  TObjArray* hadTauSelection_parts = hadTauSelection_string.Tokenize("|");
  assert(hadTauSelection_parts->GetEntries() >= 1);
  std::string hadTauSelection_part2 = ( hadTauSelection_parts->GetEntries() == 2 ) ? (dynamic_cast<TObjString*>(hadTauSelection_parts->At(1)))->GetString().Data() : "";
  delete hadTauSelection_parts;

  bool isMC = cfg_analyze.getParameter<bool>("isMC");
  bool hasLHE = cfg_analyze.getParameter<bool>("hasLHE");
  bool hasPS = cfg_analyze.getParameter<bool>("hasPS");
  bool apply_LHE_nom = cfg_analyze.getParameter<bool>("apply_LHE_nom");
  bool useObjectMultiplicity = cfg_analyze.getParameter<bool>("useObjectMultiplicity");
  std::string central_or_shift_main = cfg_analyze.getParameter<std::string>("central_or_shift");
  std::vector<std::string> central_or_shifts_local = cfg_analyze.getParameter<std::vector<std::string>>("central_or_shifts_local");
  edm::VParameterSet lumiScale = cfg_analyze.getParameter<edm::VParameterSet>("lumiScale");
  bool apply_genWeight = cfg_analyze.getParameter<bool>("apply_genWeight");
  std::string apply_topPtReweighting_str = cfg_analyze.getParameter<std::string>("apply_topPtReweighting");
  bool apply_topPtReweighting = ! apply_topPtReweighting_str.empty();
  bool apply_l1PreFireWeight = cfg_analyze.getParameter<bool>("apply_l1PreFireWeight");
  bool apply_btagSFRatio = cfg_analyze.getParameter<bool>("applyBtagSFRatio");
  bool apply_hlt_filter = cfg_analyze.getParameter<bool>("apply_hlt_filter");
  bool apply_met_filters = cfg_analyze.getParameter<bool>("apply_met_filters");
  edm::ParameterSet cfgMEtFilter = cfg_analyze.getParameter<edm::ParameterSet>("cfgMEtFilter");
  MEtFilterSelector metFilterSelector(cfgMEtFilter, isMC);
  const bool useNonNominal = cfg_analyze.getParameter<bool>("useNonNominal");
  const bool useNonNominal_jetmet = useNonNominal || ! isMC;

  if(! central_or_shifts_local.empty())
  {
    assert(central_or_shift_main == "central");
    assert(std::find(central_or_shifts_local.cbegin(), central_or_shifts_local.cend(), "central") != central_or_shifts_local.cend());
  }
  else
  {
    central_or_shifts_local = { central_or_shift_main };
  }

  edm::ParameterSet triggerWhiteList;
  if(! isMC)
  {
    triggerWhiteList = cfg_analyze.getParameter<edm::ParameterSet>("triggerWhiteList");
  }

  const edm::ParameterSet syncNtuple_cfg = cfg_analyze.getParameter<edm::ParameterSet>("syncNtuple");
  const std::string syncNtuple_tree = syncNtuple_cfg.getParameter<std::string>("tree");
  const std::string syncNtuple_output = syncNtuple_cfg.getParameter<std::string>("output");
  const bool jetCleaningByIndex = cfg_analyze.getParameter<bool>("jetCleaningByIndex");
  const bool do_sync = ! syncNtuple_tree.empty() && ! syncNtuple_output.empty();

  const edm::ParameterSet additionalEvtWeight = cfg_analyze.getParameter<edm::ParameterSet>("evtWeight");
  const bool applyAdditionalEvtWeight = additionalEvtWeight.getParameter<bool>("apply");
  EvtWeightManager * eventWeightManager = nullptr;
  if(applyAdditionalEvtWeight)
  {
    eventWeightManager = new EvtWeightManager(additionalEvtWeight);
    eventWeightManager->set_central_or_shift(central_or_shift_main);
  }

  bool isDEBUG = cfg_analyze.getParameter<bool>("isDEBUG");
  if ( isDEBUG ) std::cout << "Warning: DEBUG mode enabled -> trigger selection will not be applied for data !!" << std::endl;

  checkOptionValidity(central_or_shift_main, isMC);
  const int met_option      = useNonNominal_jetmet ? kJetMET_central_nonNominal : getMET_option(central_or_shift_main, isMC);
  const int jetPt_option    = useNonNominal_jetmet ? kJetMET_central_nonNominal : getJet_option(central_or_shift_main, isMC);
  const int hadTauPt_option = useNonNominal_jetmet ? kHadTauPt_uncorrected      : getHadTauPt_option(central_or_shift_main);

  std::cout
    << "central_or_shift = "    << central_or_shift_main << "\n"
       " -> hadTauPt_option = " << hadTauPt_option       << "\n"
       " -> met_option      = " << met_option            << "\n"
       " -> jetPt_option    = " << jetPt_option          << '\n'
  ;

  edm::ParameterSet cfg_dataToMCcorrectionInterface;
  cfg_dataToMCcorrectionInterface.addParameter<std::string>("era", era_string);
  cfg_dataToMCcorrectionInterface.addParameter<std::string>("hadTauSelection", hadTauSelection_part2);
  cfg_dataToMCcorrectionInterface.addParameter<bool>("isDEBUG", isDEBUG);
  Data_to_MC_CorrectionInterface_Base * dataToMCcorrectionInterface = nullptr;
  switch(era)
  {
    case Era::k2016: dataToMCcorrectionInterface = new Data_to_MC_CorrectionInterface_2016(cfg_dataToMCcorrectionInterface); break;
    case Era::k2017: dataToMCcorrectionInterface = new Data_to_MC_CorrectionInterface_2017(cfg_dataToMCcorrectionInterface); break;
    case Era::k2018: dataToMCcorrectionInterface = new Data_to_MC_CorrectionInterface_2018(cfg_dataToMCcorrectionInterface); break;
    default: throw cmsException("analyze_2lss", __LINE__) << "Invalid era = " << static_cast<int>(era);
  }
  const ChargeMisIdRate chargeMisIdRate(era);

  std::string applyFakeRateWeights_string = cfg_analyze.getParameter<std::string>("applyFakeRateWeights");
  int applyFakeRateWeights = -1;
  if      ( applyFakeRateWeights_string == "disabled" ) applyFakeRateWeights = kFR_disabled;
  else if ( applyFakeRateWeights_string == "2lepton"  ) applyFakeRateWeights = kFR_2lepton;
  else throw cms::Exception("analyze_2lss")
    << "Invalid Configuration parameter 'applyFakeRateWeights' = " << applyFakeRateWeights_string << " !!\n";
  std::cout << "Applying fake rate weights = " << applyFakeRateWeights_string << " (" << applyFakeRateWeights << ")\n";

  LeptonFakeRateInterface* leptonFakeRateInterface = 0;
  if ( applyFakeRateWeights == kFR_2lepton ) {
    edm::ParameterSet cfg_leptonFakeRateWeight = cfg_analyze.getParameter<edm::ParameterSet>("leptonFakeRateWeight");
    cfg_leptonFakeRateWeight.addParameter<std::string>("era", era_string);
    leptonFakeRateInterface = new LeptonFakeRateInterface(cfg_leptonFakeRateWeight);
  }

  bool fillGenEvtHistograms = cfg_analyze.getParameter<bool>("fillGenEvtHistograms");
  edm::ParameterSet cfg_EvtYieldHistManager = cfg_analyze.getParameter<edm::ParameterSet>("cfgEvtYieldHistManager");

  std::string branchName_electrons = cfg_analyze.getParameter<std::string>("branchName_electrons");
  std::string branchName_muons = cfg_analyze.getParameter<std::string>("branchName_muons");
  std::string branchName_hadTaus = cfg_analyze.getParameter<std::string>("branchName_hadTaus");
  std::string branchName_jets = cfg_analyze.getParameter<std::string>("branchName_jets");
  std::string branchName_met = cfg_analyze.getParameter<std::string>("branchName_met");
  std::string branchName_vertex = cfg_analyze.getParameter<std::string>("branchName_vertex");

  std::string branchName_genLeptons = cfg_analyze.getParameter<std::string>("branchName_genLeptons");
  std::string branchName_genHadTaus = cfg_analyze.getParameter<std::string>("branchName_genHadTaus");
  std::string branchName_genPhotons = cfg_analyze.getParameter<std::string>("branchName_genPhotons");
  std::string branchName_genJets = cfg_analyze.getParameter<std::string>("branchName_genJets");

  std::string branchName_muonGenMatch     = cfg_analyze.getParameter<std::string>("branchName_muonGenMatch");
  std::string branchName_electronGenMatch = cfg_analyze.getParameter<std::string>("branchName_electronGenMatch");
  std::string branchName_hadTauGenMatch   = cfg_analyze.getParameter<std::string>("branchName_hadTauGenMatch");
  std::string branchName_jetGenMatch      = cfg_analyze.getParameter<std::string>("branchName_jetGenMatch");

  std::string branchName_genTopQuarks = cfg_analyze.getParameter<std::string>("branchName_genTopQuarks");
  std::string branchName_genBJets = cfg_analyze.getParameter<std::string>("branchName_genBJets");
  std::string branchName_genWBosons = cfg_analyze.getParameter<std::string>("branchName_genWBosons");
  std::string branchName_genWJets = cfg_analyze.getParameter<std::string>("branchName_genWJets");
  std::string branchName_genQuarkFromTop = cfg_analyze.getParameter<std::string>("branchName_genQuarkFromTop");

  bool redoGenMatching = cfg_analyze.getParameter<bool>("redoGenMatching");
  bool genMatchingByIndex = cfg_analyze.getParameter<bool>("genMatchingByIndex");

  bool selectBDT = ( cfg_analyze.exists("selectBDT") ) ? cfg_analyze.getParameter<bool>("selectBDT") : false;

  std::string selEventsFileName_input = cfg_analyze.getParameter<std::string>("selEventsFileName_input");
  std::cout << "selEventsFileName_input = " << selEventsFileName_input << std::endl;
  RunLumiEventSelector* run_lumi_eventSelector = 0;
  if ( selEventsFileName_input != "" ) {
    edm::ParameterSet cfgRunLumiEventSelector;
    cfgRunLumiEventSelector.addParameter<std::string>("inputFileName", selEventsFileName_input);
    cfgRunLumiEventSelector.addParameter<std::string>("separator", ":");
    run_lumi_eventSelector = new RunLumiEventSelector(cfgRunLumiEventSelector);
  }

  std::string selEventsFileName_output = cfg_analyze.getParameter<std::string>("selEventsFileName_output");

  fwlite::InputSource inputFiles(cfg);
  int maxEvents = inputFiles.maxEvents();
  std::cout << " maxEvents = " << maxEvents << std::endl;
  unsigned reportEvery = inputFiles.reportAfter();

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());

  TTreeWrapper * inputTree = new TTreeWrapper(treeName.data(), inputFiles.files(), maxEvents);

  std::cout << "Loaded " << inputTree -> getFileCount() << " file(s).\n";

//--- prepare sync Ntuple
  SyncNtupleManager * snm = nullptr;
  if(do_sync)
  {
    snm = new SyncNtupleManager(syncNtuple_output, syncNtuple_tree, SyncGenMatchCharge::kLepton);
    snm->initializeBranches();
    snm->initializeHLTBranches({ triggers_1e, triggers_2e, triggers_1mu, triggers_2mu, triggers_1e1mu });
  }

//--- declare event-level variables
  EventInfo eventInfo(analysisConfig);
  if(isMC)
  {
    const double ref_genWeight = cfg_analyze.getParameter<double>("ref_genWeight");
    eventInfo.set_refGetWeight(ref_genWeight);
  }
  const std::string default_cat_str = "default";
  std::vector<std::string> evt_cat_strs = { default_cat_str };
  const std::vector<std::pair<std::string, int>> evt_htxs_binning = get_htxs_binning(isMC_signal);
  eventInfo.read_htxs(! evt_htxs_binning.empty());

  //--- HH scan
  const edm::ParameterSet hhWeight_cfg = cfg_analyze.getParameterSet("hhWeight_cfg");
  const bool apply_HH_rwgt_lo = analysisConfig.isHH_rwgt_allowed() && hhWeight_cfg.getParameter<bool>("apply_rwgt_lo");
  const bool apply_HH_rwgt_nlo = analysisConfig.isHH_rwgt_allowed() && hhWeight_cfg.getParameter<bool>("apply_rwgt_nlo");
  const HHWeightInterfaceCouplings * hhWeight_couplings = nullptr;
  const HHWeightInterfaceLO * HHWeightLO_calc = nullptr;
  const HHWeightInterfaceNLO * HHWeightNLO_calc = nullptr;
  if(apply_HH_rwgt_lo || apply_HH_rwgt_nlo)
  {
    hhWeight_couplings = new HHWeightInterfaceCouplings(hhWeight_cfg);

    if(apply_HH_rwgt_lo)
    {
      HHWeightLO_calc = new HHWeightInterfaceLO(hhWeight_couplings, hhWeight_cfg);
      evt_cat_strs = hhWeight_couplings->get_bm_names();
    }

    if(apply_HH_rwgt_nlo)
    {
      HHWeightNLO_calc = new HHWeightInterfaceNLO(hhWeight_couplings, era);
    }
  }

  const std::vector<edm::ParameterSet> tHweights = cfg_analyze.getParameterSetVector("tHweights");
  if((isMC_tH || isMC_signal) && ! tHweights.empty())
  {
    eventInfo.set_central_or_shift(central_or_shift_main);
    eventInfo.loadWeight_tH(tHweights);
    const std::vector<std::string> evt_cat_tH_strs = eventInfo.getWeight_tH_str();
    evt_cat_strs.insert(evt_cat_strs.end(), evt_cat_tH_strs.begin(), evt_cat_tH_strs.end());
  }
  EventInfoReader eventInfoReader(&eventInfo);
  if(apply_topPtReweighting)
  {
    eventInfoReader.setTopPtRwgtBranchName(apply_topPtReweighting_str);
  }
  inputTree -> registerReader(&eventInfoReader);

  RecoVertex vertex;
  RecoVertexReader vertexReader(&vertex, branchName_vertex);
  inputTree -> registerReader(&vertexReader);

  ObjectMultiplicity objectMultiplicity;
  ObjectMultiplicityReader objectMultiplicityReader(&objectMultiplicity);
  if(useObjectMultiplicity)
  {
    inputTree -> registerReader(&objectMultiplicityReader);
  }
  const int minLeptonSelection = std::min(electronSelection, muonSelection);

  hltPathReader hltPathReader_instance({ triggers_1e, triggers_2e, triggers_1mu, triggers_2mu, triggers_1e1mu });
  inputTree -> registerReader(&hltPathReader_instance);

  if(eventWeightManager)
  {
    inputTree->registerReader(eventWeightManager);
  }

  L1PreFiringWeightReader * l1PreFiringWeightReader = nullptr;
  if(apply_l1PreFireWeight)
  {
    l1PreFiringWeightReader = new L1PreFiringWeightReader(era);
    inputTree->registerReader(l1PreFiringWeightReader);
  }

  BtagSFRatioFacility * btagSFRatioFacility = nullptr;
  if(apply_btagSFRatio)
  {
    const edm::ParameterSet btagSFRatio = cfg_analyze.getParameterSet("btagSFRatio");
    btagSFRatioFacility = new BtagSFRatioFacility(btagSFRatio);
  }

//--- declare particle collections
  const bool readGenObjects = isMC && !redoGenMatching;
  RecoMuonReader* muonReader = new RecoMuonReader(era, branchName_muons, isMC, readGenObjects);
  inputTree -> registerReader(muonReader);
  RecoMuonCollectionGenMatcher muonGenMatcher;
  RecoMuonCollectionSelectorLoose preselMuonSelector(era, -1, isDEBUG);
  RecoMuonCollectionSelectorFakeable fakeableMuonSelector(era, -1, isDEBUG);
  RecoMuonCollectionSelectorTight tightMuonSelector(era, -1, isDEBUG);

  RecoElectronReader* electronReader = new RecoElectronReader(era, branchName_electrons, isMC, readGenObjects);
  inputTree -> registerReader(electronReader);
  RecoElectronCollectionGenMatcher electronGenMatcher;
  RecoElectronCollectionCleaner electronCleaner(0.3, isDEBUG);
  RecoElectronCollectionSelectorLoose preselElectronSelector(era, -1, isDEBUG);
  RecoElectronCollectionSelectorFakeable fakeableElectronSelector(era, -1, isDEBUG);
  RecoElectronCollectionSelectorTight tightElectronSelector(era, -1, isDEBUG);

  RecoHadTauReader* hadTauReader = new RecoHadTauReader(era, branchName_hadTaus, isMC, readGenObjects);
  hadTauReader->setHadTauPt_central_or_shift(hadTauPt_option);
  inputTree -> registerReader(hadTauReader);
  RecoHadTauCollectionGenMatcher hadTauGenMatcher;
  RecoHadTauCollectionCleaner hadTauCleaner(0.3, isDEBUG);
  RecoHadTauCollectionSelectorLoose looseHadTauSelector(era, -1, isDEBUG);
  looseHadTauSelector.set_if_looser(hadTauSelection_part2);
  RecoHadTauCollectionSelectorLoose tightHadTauSelector(era, -1, isDEBUG);
  tightHadTauSelector.set(hadTauSelection_part2);

  RecoJetReader* jetReader = new RecoJetReader(era, isMC, branchName_jets, readGenObjects);
  jetReader->setPtMass_central_or_shift(jetPt_option);
  jetReader->read_btag_systematics((central_or_shifts_local.size() > 1 || central_or_shift_main != "central") && isMC);
  inputTree -> registerReader(jetReader);
  RecoJetCollectionGenMatcher jetGenMatcher;
  RecoJetCollectionCleaner jetCleaner(0.4, isDEBUG);
  RecoJetCollectionCleanerByIndex jetCleanerByIndex(isDEBUG);
  RecoJetCollectionSelector jetSelector(era, -1, isDEBUG);
  RecoJetCollectionSelectorBtagLoose jetSelectorBtagLoose(era, -1, isDEBUG);
  RecoJetCollectionSelectorBtagMedium jetSelectorBtagMedium(era, -1, isDEBUG);
  const RecoJetCollectionSelectorForward jetSelectorForward(era, -1, isDEBUG);

//--- declare missing transverse energy
  RecoMEtReader* metReader = new RecoMEtReader(era, isMC, branchName_met);
  metReader->setMEt_central_or_shift(met_option);
  metReader->set_phiModulationCorrDetails(&eventInfo, &vertex);
  inputTree -> registerReader(metReader);

  MEtFilter metFilters;
  MEtFilterReader* metFilterReader = new MEtFilterReader(&metFilters, era);
  inputTree -> registerReader(metFilterReader);

//--- declare generator level information
  GenLeptonReader * genLeptonReader = nullptr;
  GenHadTauReader * genHadTauReader = nullptr;
  GenPhotonReader * genPhotonReader = nullptr;
  GenJetReader * genJetReader = nullptr;
  LHEInfoReader * lheInfoReader = nullptr;
  PSWeightReader * psWeightReader = nullptr;

  GenParticleReader * genMatchToMuonReader     = nullptr;
  GenParticleReader * genMatchToElectronReader = nullptr;
  GenParticleReader * genMatchToHadTauReader   = nullptr;
  GenParticleReader * genMatchToJetReader      = nullptr;
  if(isMC)
  {
    if(! readGenObjects)
    {
      genLeptonReader = new GenLeptonReader(branchName_genLeptons);
      inputTree -> registerReader(genLeptonReader);
      genHadTauReader = new GenHadTauReader(branchName_genHadTaus);
      inputTree -> registerReader(genHadTauReader);
      genJetReader = new GenJetReader(branchName_genJets);
      inputTree -> registerReader(genJetReader);

      if(genMatchingByIndex)
      {
        genMatchToMuonReader = new GenParticleReader(branchName_muonGenMatch);
        genMatchToMuonReader -> readGenPartFlav(true);
        inputTree -> registerReader(genMatchToMuonReader);

        genMatchToElectronReader = new GenParticleReader(branchName_electronGenMatch);
        genMatchToElectronReader -> readGenPartFlav(true);
        inputTree -> registerReader(genMatchToElectronReader);

        genMatchToHadTauReader = new GenParticleReader(branchName_hadTauGenMatch);
        genMatchToHadTauReader -> readGenPartFlav(true);
        inputTree -> registerReader(genMatchToHadTauReader);

        genMatchToJetReader = new GenParticleReader(branchName_jetGenMatch);
        genMatchToJetReader -> readGenPartFlav(true);
        inputTree -> registerReader(genMatchToJetReader);
      }
      else
      {
        genPhotonReader = new GenPhotonReader(branchName_genPhotons);
        inputTree -> registerReader(genPhotonReader);
      }
    }
    lheInfoReader = new LHEInfoReader(hasLHE);
    inputTree -> registerReader(lheInfoReader);
    psWeightReader = new PSWeightReader(hasPS, apply_LHE_nom);
    inputTree -> registerReader(psWeightReader);
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

  double evtWeightSum=0; // to devbug
  //--- initialize hadronic top tagger BDT
  HadTopTagger* hadTopTagger = new HadTopTagger();

//--- initialize BDTs used to discriminate ttH vs. ttV and ttH vs. ttbar
//    in 2lss category of ttH multilepton analysis
  std::string mvaFileName_2lss_ttV = "tthAnalysis/HiggsToTauTau/data/2lss_ttV_BDTG.weights.xml";
  std::vector<std::string> mvaInputVariables_2lss_ttV;
  mvaInputVariables_2lss_ttV.push_back("max(abs(LepGood_eta[iF_Recl[0]]),abs(LepGood_eta[iF_Recl[1]]))");
  mvaInputVariables_2lss_ttV.push_back("MT_met_lep1");
  mvaInputVariables_2lss_ttV.push_back("nJet25_Recl");
  mvaInputVariables_2lss_ttV.push_back("mindr_lep1_jet");
  mvaInputVariables_2lss_ttV.push_back("mindr_lep2_jet");
  mvaInputVariables_2lss_ttV.push_back("LepGood_conePt[iF_Recl[0]]");
  mvaInputVariables_2lss_ttV.push_back("LepGood_conePt[iF_Recl[1]]");
  TMVAInterface mva_2lss_ttV(mvaFileName_2lss_ttV, mvaInputVariables_2lss_ttV,
    { "iF_Recl[0]", "iF_Recl[1]", "iF_Recl[2]" });

  std::string mvaFileName_2lss_ttbar = "tthAnalysis/HiggsToTauTau/data/2lss_ttbar_BDTG.weights.xml";
  std::vector<std::string> mvaInputVariables_2lss_ttbar;
  mvaInputVariables_2lss_ttbar.push_back("max(abs(LepGood_eta[iF_Recl[0]]),abs(LepGood_eta[iF_Recl[1]]))");
  mvaInputVariables_2lss_ttbar.push_back("nJet25_Recl");
  mvaInputVariables_2lss_ttbar.push_back("mindr_lep1_jet");
  mvaInputVariables_2lss_ttbar.push_back("mindr_lep2_jet");
  mvaInputVariables_2lss_ttbar.push_back("min(met_pt,400)");
  mvaInputVariables_2lss_ttbar.push_back("avg_dr_jet");
  mvaInputVariables_2lss_ttbar.push_back("MT_met_lep1");
  TMVAInterface mva_2lss_ttbar(mvaFileName_2lss_ttbar, mvaInputVariables_2lss_ttbar,
    { "iF_Recl[0]", "iF_Recl[1]", "iF_Recl[2]" });
  // the same vector will be used to XGB training

  std::vector<std::string> mvaInputVariables_2lss = get_mvaInputVariables(mvaInputVariables_2lss_ttV, mvaInputVariables_2lss_ttbar);
  std::map<std::string, double> mvaInputs_2lss;

  std::string mvaFileName_Hj_tagger = "tthAnalysis/HiggsToTauTau/data/NN_for_legacy_opt/Hjtagger_legacy_xgboost_v1.weights.xml";
  std::vector<std::string> mvaInputVariables_Hj_tagger = {
    //"Jet25_lepdrmin", "max(Jet25_bDiscriminator,0.)",
    //"max(Jet25_qg,0.)", "Jet25_lepdrmax", "Jet25_pt"
    "Jet25_bDiscriminator",
    "Jet25_pt",
    "Jet25_lepdrmin",
    "Jet25_lepdrmax",
    "Jet25_qg"
  };
  TMVAInterface mva_Hj_tagger(mvaFileName_Hj_tagger, mvaInputVariables_Hj_tagger);
  mva_Hj_tagger.enableBDTTransform();

  std::string mvaFileName_Hjj_tagger = "tthAnalysis/HiggsToTauTau/data/Hjj_csv_BDTG.weights.xml";
  std::vector<std::string> mvaInputVariables_Hjj_tagger;
  mvaInputVariables_Hjj_tagger.push_back("bdtJetPair_minlepmass");
  mvaInputVariables_Hjj_tagger.push_back("bdtJetPair_sumbdt");
  mvaInputVariables_Hjj_tagger.push_back("bdtJetPair_dr");
  mvaInputVariables_Hjj_tagger.push_back("bdtJetPair_minjdr");
  mvaInputVariables_Hjj_tagger.push_back("bdtJetPair_mass");
  mvaInputVariables_Hjj_tagger.push_back("bdtJetPair_minjOvermaxjdr");
  TMVAInterface mva_Hjj_tagger(mvaFileName_Hjj_tagger, mvaInputVariables_Hjj_tagger);

  // enter the pb
  /*
  if there are >3 jets: run HTT first, then run Hj for the rest;
  if n_jet ==3, run HTT, Hj = -9; if n_jet <=2, HTT=-9, run Hj.
  */
  std::string mvaFileName_TensorFlow_2lss_ttH_tH_4cat_onlyTHQ_v4 = "tthAnalysis/HiggsToTauTau/data/NN_for_legacy_opt/2017tautag2p1samples_xsecrwonly_oldvars_tH_selection.pb";
  std::vector<std::string> mvaInputVariables_NN = {
    "jet3_pt","jet3_eta","lep1_eta",
    "jet2_pt","jet1_pt","jetFwd1_eta",
    "mT_lep1","mT_lep2","jet4_phi",
    "lep2_conePt","hadTop_BDT",
    "jet1_phi","jet2_eta",
    "n_presel_jetFwd","n_presel_jet",
    "lep1_charge","avg_dr_jet","lep1_phi",
    "Hj_tagger_hadTop","nBJetLoose",
    "jet4_pt","mindr_lep1_jet","lep1_conePt",
    "jetFwd1_pt","lep2_phi","jet2_phi","lep2_eta","mbb",
    "mindr_lep2_jet","jet4_eta","nBJetMedium",
    "Dilep_pdgId","metLD","jet3_phi","maxeta","jet1_eta"
  };
  // the order also matters
  std::vector<std::string> classes_TensorFlow_2lss_ttH_tH_4cat = {"predictions_ttH","predictions_Rest","predictions_ttW","predictions_tHQ"};
  TensorFlowInterface mva_NN(
    mvaFileName_TensorFlow_2lss_ttH_tH_4cat_onlyTHQ_v4,
    mvaInputVariables_NN,
    classes_TensorFlow_2lss_ttH_tH_4cat
  );

//--- open output file containing run:lumi:event numbers of events passing final event selection criteria
  std::ostream* selEventsFile = ( selEventsFileName_output != "" ) ? new std::ofstream(selEventsFileName_output.data(), std::ios::out) : 0;
  std::cout << "selEventsFileName_output = " << selEventsFileName_output << std::endl;

//--- declare histograms
  struct selHistManagerType
  {
    ElectronHistManager* electrons_;
    MuonHistManager* muons_;
    HadTauHistManager* hadTaus_;
    JetHistManager* jets_;
    JetHistManager* leadJet_;
    JetHistManager* subleadJet_;
    JetHistManager* BJets_loose_;
    JetHistManager* leadBJet_loose_;
    JetHistManager* subleadBJet_loose_;
    JetHistManager* BJets_medium_;
    MEtHistManager* met_;
    MEtFilterHistManager* metFilters_;
    MVAInputVarHistManager* mvaInputVariables_2lss_;
    std::map<std::string, EvtHistManager_2lss*> evt_;
    std::map<std::string, std::map<std::string, EvtHistManager_2lss*>> evt_htxs_;
    std::map<std::string, std::map<std::string, EvtHistManager_2lss*>> evt_in_decayModes_;
    std::map<std::string, std::map<std::string, std::map<std::string, EvtHistManager_2lss*>>> evt_htxs_in_decayModes_;
    EvtYieldHistManager* evtYield_;
    WeightHistManager* weights_;
  };

  std::map<std::string, GenEvtHistManager*> genEvtHistManager_beforeCuts;
  std::map<std::string, GenEvtHistManager*> genEvtHistManager_afterCuts;
  std::map<std::string, LHEInfoHistManager*> lheInfoHistManager;
  std::map<std::string, std::map<int, selHistManagerType*>> selHistManagers;

  const std::map<std::string, std::vector<double>> categories_NN = {
    {"output_NN_ttH_ee",      {0.0, 0.372851725029, 0.417126148611, 0.472538609775, 0.552644373227, 1.0}},
    {"output_NN_ttH_em",      {0.0, 0.336146706937, 0.363275650194, 0.38587797283, 0.405287896913, 0.42485721155, 0.445131216975, 0.468393114181, 0.495677675653, 0.523271652099, 0.556861718693, 0.591809523852, 0.650997335872, 1.0}},
    {"output_NN_ttH_mm",      {0.0, 0.33979852243, 0.364827437727, 0.3867203125, 0.40887020513, 0.426840842152, 0.44465645747, 0.465767926997, 0.488303020496, 0.51055658437, 0.549262538136, 0.588036217014, 0.650446923062, 1.0}},
    {"output_NN_ttW_ee",      {0.0, 0.371067397487, 0.427817254056, 0.479343688012, 0.543130393748, 0.628246753114, 1.0}},
    {"output_NN_ttW_em",      {0.0, 0.327841819271, 0.351888068965, 0.373128247902, 0.392545542021, 0.408310026608, 0.424846061502, 0.443887374249, 0.460469462678, 0.477238622854, 0.494281337826, 0.513386657952, 0.533711848197, 0.553009906357, 0.577164811518, 0.604601394495, 0.634911463675, 0.678598395425, 0.741082054423, 1.0}},
    {"output_NN_ttW_mm",      {0.0, 0.339439167603, 0.366807793797, 0.388418572766, 0.414270450605, 0.435699430178, 0.459671303678, 0.481839427108, 0.505013485008, 0.529061151818, 0.556865531897, 0.586840543306, 0.624634673301, 0.668919938188, 0.737786937571, 1.0}},
    {"output_NN_rest_ee",     {0.0, 0.390692225186, 0.44178120355, 0.48292513672, 0.522541789168, 0.564936947403, 0.610694555772, 0.654897513907, 1.0}},
    {"output_NN_rest_em",     {0.0, 0.391022040358, 0.435295534005, 0.470411315061, 0.51049663999, 0.555209638571, 0.601641720069, 0.6452841146, 1.0}},
    {"output_NN_rest_mm",     {0.0, 0.365608908907, 0.403828483817, 0.433243927402, 0.458660511153, 0.483372514182, 0.506262119418, 0.541374669081, 0.570458505259, 0.61026722902, 0.645195516132, 1.0}},
    {"output_NN_tH_ee",       {0.0, 0.448709812819, 0.533888778589, 0.651762437729, 1.0}},
    {"output_NN_tH_em",       {0.0, 0.372455420213, 0.414555824854, 0.448523285159, 0.479447825462, 0.510106000464, 0.540665579769, 0.582643520639, 0.627764586846, 0.68690681219, 0.774823730915, 1.0}},
    {"output_NN_tH_mm",       {0.0, 0.401323590943, 0.466079891197, 0.504747341883, 0.549550724173, 0.613479399002, 0.707290599558, 1.0}}
   };
   const std::map<std::string, std::vector<double>> categories_list_SVA =
   {
     {"mass_2L_ee_lj_pos", {}},
     {"mass_2L_ee_lj_neg", {}},
     {"mass_2L_ee_hj_pos", {}},
     {"mass_2L_ee_hj_neg", {}},
     {"mass_2L_em_lj_pos", {}},
     {"mass_2L_em_lj_neg", {}},
     {"mass_2L_em_hj_pos", {}},
     {"mass_2L_em_hj_neg", {}},
     {"mass_2L_mm_lj_pos", {}},
     {"mass_2L_mm_lj_neg", {}},
     {"mass_2L_mm_hj_pos", {}},
     {"mass_2L_mm_hj_neg", {}},
     {"mass_2L_cr",        {}},
     {"mass_2L_rest",        {0,1}}
   };

  for(const std::string & central_or_shift: central_or_shifts_local)
  {
    const bool skipBooking = central_or_shift != central_or_shift_main;
    std::vector<const GenMatchEntry*> genMatchDefinitions = genMatchInterface.getGenMatchDefinitions();
    for (const GenMatchEntry * genMatchDefinition : genMatchDefinitions)
    {
      std::string process_and_genMatch = process_string;
      process_and_genMatch += genMatchDefinition->getName();

      int idxLepton = genMatchDefinition->getIdx();

      selHistManagerType* selHistManager = new selHistManagerType();
      if(! skipBooking)
      {
        selHistManager->electrons_ = new ElectronHistManager(makeHistManager_cfg(process_and_genMatch,
          Form("%s/sel/electrons", histogramDir.data()), era_string, central_or_shift, "allHistograms"));
        selHistManager->electrons_->bookHistograms(fs);
        selHistManager->muons_ = new MuonHistManager(makeHistManager_cfg(process_and_genMatch,
          Form("%s/sel/muons", histogramDir.data()), era_string, central_or_shift, "allHistograms"));
        selHistManager->muons_->bookHistograms(fs);
        selHistManager->hadTaus_ = new HadTauHistManager(makeHistManager_cfg(process_and_genMatch,
          Form("%s/sel/hadTaus", histogramDir.data()), era_string, central_or_shift, "allHistograms"));
        selHistManager->hadTaus_->bookHistograms(fs);
        selHistManager->jets_ = new JetHistManager(makeHistManager_cfg(process_and_genMatch,
          Form("%s/sel/jets", histogramDir.data()), era_string, central_or_shift, "allHistograms"));
        selHistManager->jets_->bookHistograms(fs);
        selHistManager->leadJet_ = new JetHistManager(makeHistManager_cfg(process_and_genMatch,
          Form("%s/sel/leadJet", histogramDir.data()), era_string, central_or_shift, "minimalHistograms", 0));
        selHistManager->leadJet_->bookHistograms(fs);
        selHistManager->subleadJet_ = new JetHistManager(makeHistManager_cfg(process_and_genMatch,
          Form("%s/sel/subleadJet", histogramDir.data()), era_string, central_or_shift, "minimalHistograms", 1));
        selHistManager->subleadJet_->bookHistograms(fs);
        selHistManager->BJets_loose_ = new JetHistManager(makeHistManager_cfg(process_and_genMatch,
          Form("%s/sel/BJets_loose", histogramDir.data()), era_string, central_or_shift, "allHistograms"));
        selHistManager->BJets_loose_->bookHistograms(fs);
        selHistManager->leadBJet_loose_ = new JetHistManager(makeHistManager_cfg(process_and_genMatch,
          Form("%s/sel/leadBJet_loose", histogramDir.data()), era_string, central_or_shift, "minimalHistograms", 0));
        selHistManager->leadBJet_loose_->bookHistograms(fs);
        selHistManager->subleadBJet_loose_ = new JetHistManager(makeHistManager_cfg(process_and_genMatch,
          Form("%s/sel/subleadBJet_loose", histogramDir.data()), era_string, central_or_shift, "minimalHistograms", 1));
        selHistManager->subleadBJet_loose_->bookHistograms(fs);
        selHistManager->BJets_medium_ = new JetHistManager(makeHistManager_cfg(process_and_genMatch,
          Form("%s/sel/BJets_medium", histogramDir.data()), era_string, central_or_shift, "allHistograms"));
        selHistManager->BJets_medium_->bookHistograms(fs);
        selHistManager->met_ = new MEtHistManager(makeHistManager_cfg(process_and_genMatch,
          Form("%s/sel/met", histogramDir.data()), era_string, central_or_shift));
        selHistManager->met_->bookHistograms(fs);
        selHistManager->metFilters_ = new MEtFilterHistManager(makeHistManager_cfg(process_and_genMatch,
          Form("%s/sel/metFilters", histogramDir.data()), era_string, central_or_shift));
        selHistManager->metFilters_->bookHistograms(fs);
        selHistManager->mvaInputVariables_2lss_ = new MVAInputVarHistManager(makeHistManager_cfg(process_and_genMatch,
          Form("%s/sel/mvaInputs_2lss", histogramDir.data()), era_string, central_or_shift));
        selHistManager->mvaInputVariables_2lss_->bookHistograms(fs, mvaInputVariables_NN);
      }

      for(const std::string & evt_cat_str: evt_cat_strs)
      {
        if(skipBooking && evt_cat_str != default_cat_str)
        {
          continue;
        }
        const std::string process_string_new = evt_cat_str == default_cat_str ?
          process_string :
          process_string + "_" + evt_cat_str
        ;
        const std::string process_and_genMatchName = boost::replace_all_copy(
          process_and_genMatch, process_string, process_string_new
        );

        selHistManager->evt_[evt_cat_str] = new EvtHistManager_2lss(makeHistManager_cfg(
          process_and_genMatchName, Form("%s/sel/evt", histogramDir.data()), era_string, central_or_shift,
          skipBooking ? "minimalHistograms" : "allHistograms"
        ));
        selHistManager->evt_[evt_cat_str]->bookCategories(fs, categories_NN, categories_list_SVA);
        selHistManager->evt_[evt_cat_str]->bookHistograms(fs);

        for(const auto & kv: evt_htxs_binning)
        {
          const std::string htxs_process_and_genMatchName = Form("htxs_%s_%s", kv.first.data(), process_and_genMatchName.data());
          selHistManager->evt_htxs_[evt_cat_str][kv.first] = new EvtHistManager_2lss(makeHistManager_cfg(
            htxs_process_and_genMatchName, Form("%s/sel/evt", histogramDir.data()), era_string, central_or_shift,
            "minimalHistograms"
          ));
          selHistManager->evt_htxs_[evt_cat_str][kv.first]->bookCategories(fs, categories_NN, categories_list_SVA);
          selHistManager->evt_htxs_[evt_cat_str][kv.first]->bookHistograms(fs);
        }

      }

      if(isSignal)
      {
        const vstring decayModes_evt = get_key_list_hist(eventInfo, isMC_HH, isMC_VH);
        for(const std::string & decayMode_evt: decayModes_evt)
        {
          if((isMC_tH || isMC_H ) && (decayMode_evt == "hzg" || decayMode_evt == "hmm"))
          {
            continue;
          }

          std::string decayMode_and_genMatch = get_prefix(process_string, isMC_tH,  isMC_HH, isMC_H, isMC_VH);
          decayMode_and_genMatch += decayMode_evt;
	  decayMode_and_genMatch += genMatchDefinition->getName();

          for(const std::string & evt_cat_str: evt_cat_strs)
          {
            if(skipBooking && evt_cat_str != default_cat_str)
            {
              continue;
            }
            const std::string process_string_new = evt_cat_str == default_cat_str ?
              process_string:
              process_string + "_" + evt_cat_str
            ;
            const std::string decayMode_and_genMatchName = boost::replace_all_copy(
              decayMode_and_genMatch, process_string, process_string_new
            );

            selHistManager -> evt_in_decayModes_[evt_cat_str][decayMode_evt] = new EvtHistManager_2lss(makeHistManager_cfg(
              decayMode_and_genMatchName, Form("%s/sel/evt", histogramDir.data()), era_string, central_or_shift,
              skipBooking ? "minimalHistograms" : "allHistograms"
            ));
            selHistManager -> evt_in_decayModes_[evt_cat_str][decayMode_evt] -> bookCategories(fs, categories_NN, categories_list_SVA);
            selHistManager -> evt_in_decayModes_[evt_cat_str][decayMode_evt] -> bookHistograms(fs);

            for(const auto & kv: evt_htxs_binning)
            {
              const std::string htxs_decayMode_and_genMatchName = Form("htxs_%s_%s", kv.first.data(), decayMode_and_genMatchName.data());
              selHistManager->evt_htxs_in_decayModes_[evt_cat_str][decayMode_evt][kv.first] = new EvtHistManager_2lss(makeHistManager_cfg(
                htxs_decayMode_and_genMatchName, Form("%s/sel/evt", histogramDir.data()), era_string, central_or_shift,
                "minimalHistograms"
              ));
              selHistManager->evt_htxs_in_decayModes_[evt_cat_str][decayMode_evt][kv.first]->bookCategories(fs, categories_NN, categories_list_SVA);
              selHistManager->evt_htxs_in_decayModes_[evt_cat_str][decayMode_evt][kv.first]->bookHistograms(fs);
            }
          }
        }
      }

      if(! skipBooking)
      {
        edm::ParameterSet cfg_EvtYieldHistManager_sel = makeHistManager_cfg(process_and_genMatch,
          Form("%s/sel/evtYield", histogramDir.data()), era_string, central_or_shift);
        cfg_EvtYieldHistManager_sel.addParameter<edm::ParameterSet>("runPeriods", cfg_EvtYieldHistManager);
        cfg_EvtYieldHistManager_sel.addParameter<bool>("isMC", isMC);
        selHistManager->evtYield_ = new EvtYieldHistManager(cfg_EvtYieldHistManager_sel);
        selHistManager->evtYield_->bookHistograms(fs);
        selHistManager->weights_ = new WeightHistManager(makeHistManager_cfg(process_and_genMatch,
          Form("%s/sel/weights", histogramDir.data()), era_string, central_or_shift));
        selHistManager->weights_->bookHistograms(fs, { "genWeight", "pileupWeight", "triggerWeight", "data_to_MC_correction", "fakeRate" });
      }
      selHistManagers[central_or_shift][idxLepton] = selHistManager;
    }

    if(isMC && ! skipBooking)
    {
      genEvtHistManager_beforeCuts[central_or_shift] = new GenEvtHistManager(makeHistManager_cfg(process_string,
        Form("%s/unbiased/genEvt", histogramDir.data()), era_string, central_or_shift));
      genEvtHistManager_beforeCuts[central_or_shift]->bookHistograms(fs);
      genEvtHistManager_afterCuts[central_or_shift] = new GenEvtHistManager(makeHistManager_cfg(process_string,
        Form("%s/sel/genEvt", histogramDir.data()), era_string, central_or_shift));
      genEvtHistManager_afterCuts[central_or_shift]->bookHistograms(fs);
      lheInfoHistManager[central_or_shift] = new LHEInfoHistManager(makeHistManager_cfg(process_string,
        Form("%s/sel/lheInfo", histogramDir.data()), era_string, central_or_shift));
      lheInfoHistManager[central_or_shift]->bookHistograms(fs);

      if(eventWeightManager)
      {
        genEvtHistManager_beforeCuts[central_or_shift]->bookHistograms(fs, eventWeightManager);
        genEvtHistManager_afterCuts[central_or_shift]->bookHistograms(fs, eventWeightManager);
      }
    }
  }

  NtupleFillerBDT<float, int>* bdt_filler = nullptr;
  typedef std::remove_pointer<decltype(bdt_filler)>::type::float_type float_type;
  typedef std::remove_pointer<decltype(bdt_filler)>::type::int_type int_type;

  if ( selectBDT ) {
    bdt_filler = new std::remove_pointer<decltype(bdt_filler)>::type(
      makeHistManager_cfg(process_string, Form("%s/sel/evtntuple", histogramDir.data()), era_string, central_or_shift_main)
    );
    bdt_filler->register_variable<float_type>(
      "lep1_pt", "lep1_conePt", "lep1_eta", "lep1_phi", "lep1_tth_mva", "mindr_lep1_jet",
      "mindr_lep2_jet", "mT_lep1",  "MT_met_lep1",
      "lep2_pt", "lep2_conePt", "lep2_eta", "lep2_phi", "max_Lep_eta", "avg_dr_lep",
      "lep2_tth_mva", "mT_lep2",
      "avg_dr_jet",  "nJet25_Recl", "ptmiss", "htmiss",
       "dr_leps",
      "lep1_genLepPt", "lep2_genLepPt",
      "lep1_frWeight", "lep2_frWeight",
      "mvaOutput_2lss_ttV",  "mvaOutput_2lss_ttbar", "mvaDiscr_2lss",
      "Hj_tagger",
      "lumiScale", "genWeight", "evtWeight", "min(met_pt,400)",
      "mbb_loose", "mbb_medium",
      "res_HTT", "HadTop_pt", "genTopPt_CSVsort4rd",
      "massLT", "massL", "min_Deta_mostfwdJet_jet", "min_Deta_leadfwdJet_jet",
      "met_LD", "jet1_pt", "jet1_eta",
      "jet1_pt", "jet1_eta", "jet1_phi", "jet1_E",
      "jet2_pt", "jet2_eta", "jet2_phi", "jet2_E",
      "jet3_pt", "jet3_eta", "jet3_phi", "jet3_E",
      "jet4_pt", "jet4_eta", "jet4_phi", "jet4_E",
      "mostFwdJet_eta", "mostFwdJet_pt", "mostFwdJet_phi", "mostFwdJet_E",
      "leadFwdJet_eta", "leadFwdJet_pt", "leadFwdJet_phi", "leadFwdJet_E"
    );
    for(const std::string & evt_cat_str: evt_cat_strs)
    {
      bdt_filler->register_variable<float_type>(evt_cat_str);
    }
    bdt_filler->register_variable<int_type>(
      "nJet", "nBJetLoose", "nBJetMedium", "nLep", "nJetForward",
      "lep1_isTight", "lep2_isTight", "failsTightChargeCut",
      "nElectron", "sum_Lep_charge",
      "hadtruth", "bWj1Wj2_isGenMatched_CSVsort4rd"
    );
    bdt_filler->bookTree(fs);
  }

  int analyzedEntries = 0;
  int selectedEntries = 0;
  double selectedEntries_weighted = 0.;
  std::map<std::string, int> selectedEntries_byGenMatchType;                                    // key = process_and_genMatch
  std::map<std::string, std::map<std::string, double>> selectedEntries_weighted_byGenMatchType; // key = central_or_shift, process_and_genMatch
  TH1* histogram_analyzedEntries = fs.make<TH1D>("analyzedEntries", "analyzedEntries", 1, -0.5, +0.5);
  TH1* histogram_selectedEntries = fs.make<TH1D>("selectedEntries", "selectedEntries", 1, -0.5, +0.5);
  //---------------------------------------------------------------------------
  // CV: book a few extra histogram for the HIG-19-008 paper
  TH1* histogram_genBJetPt_fake = fs.make<TH1D>("genBJetPt_fake", "genBJetPt_fake", 250, 0., 250.);
  TH1* histogram_genBJetEta_fake = fs.make<TH1D>("genBJetEta_fake", "genBJetEta_fake", 100, -5., +5.);
  TH1* histogram_numGenBJets_fake = fs.make<TH1D>("numGenBJets_fake", "numGenBJets_fake", 10, -0.5, 9.5);
  TH1* histogram_numRecBJets_loose_fake = fs.make<TH1D>("numRecBJetsLoose_fake", "numRecBJetsLoose_fake", 10, -0.5, 9.5);
  TH1* histogram_numRecBJets_medium_fake = fs.make<TH1D>("numRecBJetsMedium_fake", "numRecBJetsMedium_fake", 10, -0.5, 9.5);
  TH1* histogram_genBJetPt_nonfake = fs.make<TH1D>("genBJetPt_nonfake", "genBJetPt_nonfake", 250, 0., 250.);
  TH1* histogram_genBJetEta_nonfake = fs.make<TH1D>("genBJetEta_nonfake", "genBJetEta_nonfake", 100, -5., +5.);
  TH1* histogram_numGenBJets_nonfake = fs.make<TH1D>("numGenBJets_nonfake", "numGenBJets_nonfake", 10, -0.5, 9.5);
  TH1* histogram_numRecBJets_loose_nonfake = fs.make<TH1D>("numRecBJetsLoose_nonfake", "numRecBJetsLoose_nonfake", 10, -0.5, 9.5);
  TH1* histogram_numRecBJets_medium_nonfake = fs.make<TH1D>("numRecBJetsMedium_nonfake", "numRecBJetsMedium_nonfake", 10, -0.5, 9.5);
  //---------------------------------------------------------------------------
  cutFlowTableType cutFlowTable;
  const edm::ParameterSet cutFlowTableCfg = makeHistManager_cfg(
    process_string, Form("%s/sel/cutFlow", histogramDir.data()), era_string, central_or_shift_main
  );
  const std::vector<std::string> cuts = {
    "run:ls:event selection",
    "object multiplicity",
    "trigger",
    ">= 2 presel leptons",
    ">= 2 jets",
    ">= 2 sel leptons",
    "<= 2 tight leptons",
    "HLT filter matching",
    "Hadronic selection",
    "sel tau veto",
    "m(ll) > 12 GeV",
    "lead lepton pT > 25 GeV && sublead lepton pT > 15 GeV",
    "tight lepton charge",
    "sel lepton-pair OS/SS charge",
    "Z-boson mass veto",
    "met LD > 30 GeV",
    "MEt filters",
    "sel lepton-pair gen=rec charge match",
    "signal region veto",
  };
  CutFlowTableHistManager * cutFlowHistManager = new CutFlowTableHistManager(cutFlowTableCfg, cuts);
  cutFlowHistManager->bookHistograms(fs);

  bool isDEBUG_TF = false;
  while(inputTree -> hasNextEvent() && (! run_lumi_eventSelector || (run_lumi_eventSelector && ! run_lumi_eventSelector -> areWeDone())))
  {
    if(inputTree -> canReport(reportEvery))
    {
      std::cout << "processing Entry " << inputTree -> getCurrentMaxEventIdx()
                << " or " << inputTree -> getCurrentEventIdx() << " entry in #"
                << (inputTree -> getProcessedFileCount() - 1)
                << " (" << eventInfo
                << ") file (" << selectedEntries << " Entries selected)\n";
    }
    ++analyzedEntries;
    histogram_analyzedEntries->Fill(0.);

    if (run_lumi_eventSelector && !(*run_lumi_eventSelector)(eventInfo))
    {
      continue;
    }
    EvtWeightRecorder evtWeightRecorder(central_or_shifts_local, central_or_shift_main, isMC);
    if(skipEvery > 1)
    {
      if(selectBDT)
      {
        if(eventInfo.event % skipEvery == 0)
        {
          // skip every N-th event when running in BDT mode
          continue;
        }
        else
        {
          // rescale event weight by N / (N - 1)
          evtWeightRecorder.record_rescaling(skipEvery / (skipEvery - 1.));
        }
      }
      else
      {
        if(eventInfo.event % skipEvery != 0)
        {
          // we enter here (N-1) times -> select every N-th event when running the analysis regularly
          continue;
        }
        else
        {
          // rescale event weight by N
          evtWeightRecorder.record_rescaling(skipEvery);
        }
      }
    }
    cutFlowTable.update("run:ls:event selection", evtWeightRecorder.get(central_or_shift_main));
    cutFlowHistManager->fillHistograms("run:ls:event selection", evtWeightRecorder.get(central_or_shift_main));

    if ( isDEBUG || run_lumi_eventSelector) {
      std::cout << "event #" << inputTree -> getCurrentMaxEventIdx() << ' ' << eventInfo << '\n';
    }

    if(run_lumi_eventSelector)
    {
      std::cout << "processing Entry " << inputTree -> getCurrentMaxEventIdx() << ": " << eventInfo << '\n';
      if(inputTree -> isOpen())
      {
        std::cout << "input File = " << inputTree -> getCurrentFileName() << '\n';
      }
    }

    if(useObjectMultiplicity)
    {
      if(objectMultiplicity.getNRecoLepton(minLeptonSelection) < 2 ||
         objectMultiplicity.getNRecoLepton(kTight)             > 2  )
      {
        if(isDEBUG || run_lumi_eventSelector)
        {
          std::cout << "event " << eventInfo.str() << " FAILS preliminary object multiplicity cuts\n";
        }
        continue;
      }
    }
    cutFlowTable.update("object multiplicity", evtWeightRecorder.get(central_or_shift_main));
    cutFlowHistManager->fillHistograms("object multiplicity", evtWeightRecorder.get(central_or_shift_main));

//--- build collections of generator level particles (before any cuts are applied, to check distributions in unbiased event samples)
    std::vector<GenLepton> genLeptons;
    std::vector<GenLepton> genElectrons;
    std::vector<GenLepton> genMuons;
    std::vector<GenHadTau> genHadTaus;
    std::vector<GenPhoton> genPhotons;
    std::vector<GenJet> genJets;

    std::vector<GenParticle> muonGenMatch;
    std::vector<GenParticle> electronGenMatch;
    std::vector<GenParticle> hadTauGenMatch;
    std::vector<GenParticle> jetGenMatch;
    if(isMC && fillGenEvtHistograms)
    {
      if(genLeptonReader)
      {
        genLeptons = genLeptonReader->read();
        for(const GenLepton & genLepton: genLeptons)
        {
          const int abs_pdgId = std::abs(genLepton.pdgId());
          switch(abs_pdgId)
          {
            case 11: genElectrons.push_back(genLepton); break;
            case 13: genMuons.push_back(genLepton);     break;
            default: assert(0);
          }
        }
      }
      if(genHadTauReader) genHadTaus = genHadTauReader->read();
      if(genPhotonReader) genPhotons = genPhotonReader->read();
      if(genJetReader)    genJets = genJetReader->read();

      if(genMatchToMuonReader)     muonGenMatch = genMatchToMuonReader->read();
      if(genMatchToElectronReader) electronGenMatch = genMatchToElectronReader->read();
      if(genMatchToHadTauReader)   hadTauGenMatch = genMatchToHadTauReader->read();
      if(genMatchToJetReader)      jetGenMatch = genMatchToJetReader->read();
    }

    if(isMC)
    {
      if(apply_genWeight)         evtWeightRecorder.record_genWeight(eventInfo);
      if(eventWeightManager)      evtWeightRecorder.record_auxWeight(eventWeightManager);
      if(l1PreFiringWeightReader) evtWeightRecorder.record_l1PrefireWeight(l1PreFiringWeightReader);
      if(apply_topPtReweighting)  evtWeightRecorder.record_toppt_rwgt(eventInfo.topPtRwgtSF);
      lheInfoReader->read();
      psWeightReader->read();
      evtWeightRecorder.record_lheScaleWeight(lheInfoReader);
      evtWeightRecorder.record_psWeight(psWeightReader);
      evtWeightRecorder.record_puWeight(&eventInfo);
      evtWeightRecorder.record_nom_tH_weight(&eventInfo);
      evtWeightRecorder.record_lumiScale(lumiScale);
      for(const std::string & central_or_shift: central_or_shifts_local)
      {
        if(central_or_shift != central_or_shift_main)
        {
          continue;
        }
        genEvtHistManager_beforeCuts[central_or_shift]->fillHistograms(
          genElectrons, genMuons, genHadTaus, genPhotons, genJets, evtWeightRecorder.get_inclusive(central_or_shift)
        );
        if(eventWeightManager)
        {
          genEvtHistManager_beforeCuts[central_or_shift]->fillHistograms(
            eventWeightManager, evtWeightRecorder.get_inclusive(central_or_shift)
          );
        }
      }
    }

    bool isTriggered_1e = hltPaths_isTriggered(triggers_1e, triggerWhiteList, eventInfo, isMC, isDEBUG);
    bool isTriggered_2e = hltPaths_isTriggered(triggers_2e, triggerWhiteList, eventInfo, isMC, isDEBUG);
    bool isTriggered_1mu = hltPaths_isTriggered(triggers_1mu, triggerWhiteList, eventInfo, isMC, isDEBUG);
    bool isTriggered_2mu = hltPaths_isTriggered(triggers_2mu, triggerWhiteList, eventInfo, isMC, isDEBUG);
    bool isTriggered_1e1mu = hltPaths_isTriggered(triggers_1e1mu, triggerWhiteList, eventInfo, isMC, isDEBUG);

    bool selTrigger_1e = use_triggers_1e && isTriggered_1e;
    bool selTrigger_2e = use_triggers_2e && isTriggered_2e;
    bool selTrigger_1mu = use_triggers_1mu && isTriggered_1mu;
    bool selTrigger_2mu = use_triggers_2mu && isTriggered_2mu;
    bool selTrigger_1e1mu = use_triggers_1e1mu && isTriggered_1e1mu;

    //--- build collections of generator level particles
    //std::cout << "built gen variable" << std::endl;
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

    if ( !(selTrigger_1e || selTrigger_2e || selTrigger_1mu || selTrigger_2mu || selTrigger_1e1mu) ) {
      if ( run_lumi_eventSelector ) {
    std::cout << "event " << eventInfo.str() << " FAILS trigger selection." << std::endl;
	std::cout << " (selTrigger_1e = " << selTrigger_1e
		  << ", selTrigger_2e = " << selTrigger_2e
		  << ", selTrigger_1mu = " << selTrigger_1mu
		  << ", selTrigger_2mu = " << selTrigger_2mu
		  << ", selTrigger_1e1mu = " << selTrigger_1e1mu << ")" << std::endl;
      }
      continue;
    }

//--- rank triggers by priority and ignore triggers of lower priority if a trigger of higher priority has fired for given event;
//    the ranking of the triggers is as follows: 2mu, 1e1mu, 2e, 1mu, 1e
// CV: this logic is necessary to avoid that the same event is selected multiple times when processing different primary datasets
    if ( !isMC && !isDEBUG ) {
      if ( selTrigger_1e && (isTriggered_1mu || isTriggered_2mu || isTriggered_1e1mu) ) {
	if ( run_lumi_eventSelector ) {
          std::cout << "event " << eventInfo.str() << " FAILS trigger selection." << std::endl;
	  std::cout << " (selTrigger_1e = " << selTrigger_1e
		    << ", isTriggered_1mu = " << isTriggered_1mu
		    << ", isTriggered_2mu = " << isTriggered_2mu
		    << ", isTriggered_1e1mu = " << isTriggered_1e1mu << ")" << std::endl;
	}
	continue;
      }
      if ( selTrigger_1e && isTriggered_2e && era != Era::k2018 ) {
        if ( run_lumi_eventSelector ) {
          std::cout << "event " << eventInfo.str() << " FAILS trigger selection." << std::endl;
          std::cout << " (selTrigger_1e = " << selTrigger_1e
                    << ", isTriggered_2e = " << isTriggered_2e << ")" << std::endl;
        }
        continue;
      }
      if ( selTrigger_2e && (isTriggered_2mu || isTriggered_1e1mu) ) {
	if ( run_lumi_eventSelector ) {
      std::cout << "event " << eventInfo.str() << " FAILS trigger selection." << std::endl;
	  std::cout << " (selTrigger_2e = " << selTrigger_2e
		    << ", isTriggered_2mu = " << isTriggered_2mu
		    << ", isTriggered_1e1mu = " << isTriggered_1e1mu << ")" << std::endl;
	}
	continue;
      }
      if ( selTrigger_1mu && (isTriggered_2e || isTriggered_2mu || isTriggered_1e1mu) ) {
	if ( run_lumi_eventSelector ) {
      std::cout << "event " << eventInfo.str() << " FAILS trigger selection." << std::endl;
	  std::cout << " (selTrigger_1mu = " << selTrigger_1mu
		    << ", isTriggered_2e = " << isTriggered_2e
		    << ", isTriggered_2mu = " << isTriggered_2mu
		    << ", isTriggered_1e1mu = " << isTriggered_1e1mu << ")" << std::endl;
	}
	continue;
      }
      if ( selTrigger_1e1mu && isTriggered_2mu ) {
	if ( run_lumi_eventSelector ) {
      std::cout << "event " << eventInfo.str() << " FAILS trigger selection." << std::endl;
	  std::cout << " (selTrigger_1e1mu = " << selTrigger_1e1mu
		    << ", isTriggered_2mu = " << isTriggered_2mu << ")" << std::endl;
	}
	continue;
      }
    }
    cutFlowTable.update("trigger", evtWeightRecorder.get(central_or_shift_main));
    cutFlowHistManager->fillHistograms("trigger", evtWeightRecorder.get(central_or_shift_main));

    if ( (selTrigger_2mu   && !apply_offline_e_trigger_cuts_2mu)   ||
	 (selTrigger_1mu   && !apply_offline_e_trigger_cuts_1mu)   ||
	 (selTrigger_2e    && !apply_offline_e_trigger_cuts_2e)    ||
	 (selTrigger_1e1mu && !apply_offline_e_trigger_cuts_1e1mu) ||
	 (selTrigger_1e    && !apply_offline_e_trigger_cuts_1e)    ) {
      fakeableElectronSelector.disable_offline_e_trigger_cuts();
      tightElectronSelector.disable_offline_e_trigger_cuts();
    } else {
      fakeableElectronSelector.enable_offline_e_trigger_cuts();
      tightElectronSelector.enable_offline_e_trigger_cuts();
    }

//--- build collections of electrons, muons and hadronic taus;
//    resolve overlaps in order of priority: muon, electron,
    const std::vector<RecoMuon> muons = muonReader->read();
    const std::vector<const RecoMuon*> muon_ptrs = convert_to_ptrs(muons);
    const std::vector<const RecoMuon*> cleanedMuons = muon_ptrs; // CV: no cleaning needed for muons, as they have the highest priority in the overlap removal
    const std::vector<const RecoMuon*> preselMuons = preselMuonSelector(cleanedMuons, isHigherConePt);
    const std::vector<const RecoMuon*> fakeableMuons = fakeableMuonSelector(preselMuons, isHigherConePt);
    const std::vector<const RecoMuon*> tightMuons = tightMuonSelector(fakeableMuons, isHigherConePt);
    if(isDEBUG || run_lumi_eventSelector)
    {
      printCollection("preselMuons",   preselMuons);
      printCollection("fakeableMuons", fakeableMuons);
      printCollection("tightMuons",    tightMuons);
    }

    const std::vector<RecoElectron> electrons = electronReader->read();
    const std::vector<const RecoElectron*> electron_ptrs = convert_to_ptrs(electrons);
    const std::vector<const RecoElectron*> cleanedElectrons = electronCleaner(electron_ptrs, preselMuons);
    const std::vector<const RecoElectron*> preselElectrons = preselElectronSelector(cleanedElectrons, isHigherConePt);
    const std::vector<const RecoElectron*> preselElectronsUncleaned = preselElectronSelector(electron_ptrs, isHigherConePt);
    const std::vector<const RecoElectron*> fakeableElectrons = fakeableElectronSelector(preselElectrons, isHigherConePt);
    const std::vector<const RecoElectron*> tightElectrons = tightElectronSelector(fakeableElectrons, isHigherConePt);
    if(isDEBUG || run_lumi_eventSelector)
    {
      printCollection("preselElectrons",   preselElectrons);
      printCollection("preselElectronsUncleaned", preselElectronsUncleaned);
      printCollection("preselElectronsUncleaned", preselElectronsUncleaned);
      printCollection("fakeableElectrons", fakeableElectrons);
      printCollection("tightElectrons",    tightElectrons);
    }

    const std::vector<const RecoLepton*> preselLeptonsFull = mergeLeptonCollections(preselElectrons, preselMuons, isHigherConePt);
    const std::vector<const RecoLepton*> preselLeptonsFullUncleaned = mergeLeptonCollections(preselElectronsUncleaned, preselMuons, isHigherConePt);
    const std::vector<const RecoLepton*> fakeableLeptonsFull = mergeLeptonCollections(fakeableElectrons, fakeableMuons, isHigherConePt);
    const std::vector<const RecoLepton*> tightLeptonsFull = mergeLeptonCollections(tightElectrons, tightMuons, isHigherConePt);

    const std::vector<const RecoLepton*> preselLeptons = pickFirstNobjects(preselLeptonsFull, 2);
    const std::vector<const RecoLepton*> fakeableLeptons = pickFirstNobjects(fakeableLeptonsFull, 2);
    const std::vector<const RecoLepton*> tightLeptons = getIntersection(fakeableLeptons, tightLeptonsFull, isHigherConePt);

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
    }
    const std::vector<const RecoLepton*> selLeptons_full = mergeLeptonCollections(selElectrons, selMuons, isHigherConePt);
    if(!(electronSelection == muonSelection)) selLeptons = getIntersection(fakeableLeptons, selLeptons_full, isHigherConePt);

    const std::vector<RecoHadTau> hadTaus = hadTauReader->read();
    const std::vector<const RecoHadTau*> hadTau_ptrs = convert_to_ptrs(hadTaus);
    const std::vector<const RecoHadTau*> cleanedHadTaus = hadTauCleaner(hadTau_ptrs, preselMuons, preselElectrons);
    const std::vector<const RecoHadTau*> looseHadTaus = looseHadTauSelector(cleanedHadTaus, isHigherPt);
    const std::vector<const RecoHadTau*> selHadTaus = tightHadTauSelector(looseHadTaus, isHigherPt);

    if(isDEBUG || run_lumi_eventSelector)
    {
      printCollection("selMuons", selMuons);
      printCollection("selElectrons", selElectrons);
      printCollection("selLeptons", selLeptons);
      printCollection("looseHadTaus", looseHadTaus);
      printCollection("selHadTaus", selHadTaus);
    }

//--- build collections of jets and select subset of jets passing b-tagging criteria
    const std::vector<RecoJet> jets = jetReader->read();
    const std::vector<const RecoJet*> jet_ptrs = convert_to_ptrs(jets);
    const std::vector<const RecoJet*> cleanedJets = jetCleaningByIndex ?
      jetCleanerByIndex(jet_ptrs, selectBDT ? selLeptons_full : fakeableLeptonsFull, looseHadTaus) :
      jetCleaner       (jet_ptrs, selectBDT ? selLeptons_full : fakeableLeptonsFull, looseHadTaus)
    ;
    const std::vector<const RecoJet*> selJets = jetSelector(cleanedJets, isHigherPt);
    const std::vector<const RecoJet*> selBJets_loose = jetSelectorBtagLoose(cleanedJets, isHigherPt);
    const std::vector<const RecoJet*> selBJets_medium = jetSelectorBtagMedium(cleanedJets, isHigherPt);
    const std::vector<const RecoJet *> selJetsForward = jetSelectorForward(cleanedJets, isHigherPt);
    if(isDEBUG || run_lumi_eventSelector)
    {
      printCollection("uncleanedJets", jet_ptrs);
      printCollection("cleanedJets", cleanedJets);
      printCollection("selJets",       selJets);
      printCollection("selJetsForward", selJetsForward);
    }

//--- build collections of generator level particles (after some cuts are applied, to safe computing time)
    if(isMC && redoGenMatching && ! fillGenEvtHistograms)
    {
      if(genLeptonReader)
      {
        genLeptons = genLeptonReader->read();
        for(const GenLepton & genLepton: genLeptons)
        {
          const int abs_pdgId = std::abs(genLepton.pdgId());
          switch(abs_pdgId)
          {
            case 11: genElectrons.push_back(genLepton); break;
            case 13: genMuons.push_back(genLepton);     break;
            default: assert(0);
          }
        }
      }
      if(genHadTauReader) genHadTaus = genHadTauReader->read();
      if(genPhotonReader) genPhotons = genPhotonReader->read();
      if(genJetReader)    genJets = genJetReader->read();

      if(genMatchToMuonReader)     muonGenMatch = genMatchToMuonReader->read();
      if(genMatchToElectronReader) electronGenMatch = genMatchToElectronReader->read();
      if(genMatchToHadTauReader)   hadTauGenMatch = genMatchToHadTauReader->read();
      if(genMatchToJetReader)      jetGenMatch = genMatchToJetReader->read();
    }

//--- match reconstructed to generator level particles
    if(isMC && redoGenMatching)
    {
      if(genMatchingByIndex)
      {
        muonGenMatcher.addGenLeptonMatchByIndex(preselMuons, muonGenMatch, GenParticleType::kGenMuon);
        muonGenMatcher.addGenHadTauMatch       (preselMuons, genHadTaus);
        muonGenMatcher.addGenJetMatch          (preselMuons, genJets);

        electronGenMatcher.addGenLeptonMatchByIndex(preselElectrons, electronGenMatch, GenParticleType::kGenElectron);
        electronGenMatcher.addGenPhotonMatchByIndex(preselElectrons, electronGenMatch);
        electronGenMatcher.addGenHadTauMatch       (preselElectrons, genHadTaus);
        electronGenMatcher.addGenJetMatch          (preselElectrons, genJets);

        hadTauGenMatcher.addGenLeptonMatchByIndex(selHadTaus, hadTauGenMatch, GenParticleType::kGenAnyLepton);
        hadTauGenMatcher.addGenHadTauMatch       (selHadTaus, genHadTaus);
        hadTauGenMatcher.addGenJetMatch          (selHadTaus, genJets);

        jetGenMatcher.addGenLeptonMatch    (selJets, genLeptons);
        jetGenMatcher.addGenHadTauMatch    (selJets, genHadTaus);
        jetGenMatcher.addGenJetMatchByIndex(selJets, jetGenMatch);
      }
      else
      {
        muonGenMatcher.addGenLeptonMatch(preselMuons, genMuons);
        muonGenMatcher.addGenHadTauMatch(preselMuons, genHadTaus);
        muonGenMatcher.addGenJetMatch   (preselMuons, genJets);

        electronGenMatcher.addGenLeptonMatch(preselElectrons, genElectrons);
        electronGenMatcher.addGenPhotonMatch(preselElectrons, genPhotons);
        electronGenMatcher.addGenHadTauMatch(preselElectrons, genHadTaus);
        electronGenMatcher.addGenJetMatch   (preselElectrons, genJets);

        hadTauGenMatcher.addGenLeptonMatch(selHadTaus, genLeptons);
        hadTauGenMatcher.addGenHadTauMatch(selHadTaus, genHadTaus);
        hadTauGenMatcher.addGenJetMatch   (selHadTaus, genJets);

        jetGenMatcher.addGenLeptonMatch(selJets, genLeptons);
        jetGenMatcher.addGenHadTauMatch(selJets, genHadTaus);
        jetGenMatcher.addGenJetMatch   (selJets, genJets);
      }
    }

//--- apply preselection
    // require at least two leptons passing loose preselection criteria
    if ( !(preselLeptonsFull.size() >= 2) ) {
      if ( run_lumi_eventSelector ) {
    std::cout << "event " << eventInfo.str() << " FAILS preselLeptons selection." << std::endl;
  printCollection("preselLeptons", preselLeptonsFull);
      }
      continue;
    }
    cutFlowTable.update(">= 2 presel leptons", evtWeightRecorder.get(central_or_shift_main));
    cutFlowHistManager->fillHistograms(">= 2 presel leptons", evtWeightRecorder.get(central_or_shift_main));

//--- compute MHT and linear MET discriminant (met_LD)
    const RecoMEt met = metReader->read();
    const Particle::LorentzVector mht_p4 = compMHT(fakeableLeptonsFull, looseHadTaus, selJets);
    const double met_LD = compMEt_LD(met.p4(), mht_p4);

//--- apply final event selection
    if ( !(selLeptons.size() >= 2) ) {
      if ( run_lumi_eventSelector ) {
    std::cout << "event " << eventInfo.str() << " FAILS selLeptons selection." << std::endl;
	printCollection("selLeptons", selLeptons);
      }
      continue;
    }
    cutFlowTable.update(">= 2 sel leptons", evtWeightRecorder.get(central_or_shift_main));
    cutFlowHistManager->fillHistograms(">= 2 sel leptons", evtWeightRecorder.get(central_or_shift_main));
    const RecoLepton* selLepton_lead = selLeptons[0];
    int selLepton_lead_type = getLeptonType(selLepton_lead->pdgId());
    const RecoLepton* selLepton_sublead = selLeptons[1];
    int selLepton_sublead_type = getLeptonType(selLepton_sublead->pdgId());
    const leptonChargeFlipGenMatchEntry& selLepton_genMatch = getLeptonChargeFlipGenMatch(leptonGenMatch_definitions, selLepton_lead, selLepton_sublead);

    if(isMC)
    {
//--- compute event-level weight for data/MC correction of b-tagging efficiency and mistag rate
//   (using the method "Event reweighting using scale factors calculated with a tag and probe method",
//    described on the BTV POG twiki https://twiki.cern.ch/twiki/bin/view/CMS/BTagShapeCalibration )
      evtWeightRecorder.record_btagWeight(selJets);
      if(btagSFRatioFacility)
      {
        evtWeightRecorder.record_btagSFRatio(btagSFRatioFacility, selJets.size());
      }

      if(isMC_EWK)
      {
        evtWeightRecorder.record_ewk_jet(selJets);
        evtWeightRecorder.record_ewk_bjet(selBJets_medium);
      }

      dataToMCcorrectionInterface->setLeptons({ selLepton_lead, selLepton_sublead }, true);

//--- apply data/MC corrections for trigger efficiency
      evtWeightRecorder.record_leptonTriggerEff(dataToMCcorrectionInterface);

//--- apply data/MC corrections for efficiencies for lepton to pass loose identification and isolation criteria
      evtWeightRecorder.record_leptonIDSF_recoToLoose(dataToMCcorrectionInterface);

//--- apply data/MC corrections for efficiencies of leptons passing the loose identification and isolation criteria
//    to also pass the tight identification and isolation criteria
      if(electronSelection >= kFakeable && muonSelection >= kFakeable)
      {
        // apply looseToTight SF to leptons matched to generator-level prompt leptons and passing Tight selection conditions
        evtWeightRecorder.record_leptonIDSF_looseToTight(dataToMCcorrectionInterface, false);
      }
    }

    const bool passesTight_lepton_lead = isMatched(*selLepton_lead, tightElectrons) || isMatched(*selLepton_lead, tightMuons);
    const bool passesTight_lepton_sublead = isMatched(*selLepton_sublead, tightElectrons) || isMatched(*selLepton_sublead, tightMuons);

    if(leptonFakeRateInterface)
    {
      evtWeightRecorder.record_jetToLepton_FR_lead(leptonFakeRateInterface, selLepton_lead);
      evtWeightRecorder.record_jetToLepton_FR_sublead(leptonFakeRateInterface, selLepton_sublead);
    }

    if(! selectBDT)
    {
      if(applyFakeRateWeights == kFR_2lepton)
      {
        evtWeightRecorder.compute_FR_2l(passesTight_lepton_lead, passesTight_lepton_sublead);
      }
    }

    // require exactly two leptons passing tight selection criteria, to avoid overlap with other channels
    if ( !(tightLeptonsFull.size() <= 2) ) {
      if ( run_lumi_eventSelector ) {
        std::cout << "event " << eventInfo.str() << " FAILS tightLeptons selection.\n";
        printCollection("tightLeptonsFull", tightLeptonsFull);
      }
      continue;
    }
    cutFlowTable.update("<= 2 tight leptons", evtWeightRecorder.get(central_or_shift_main));
    cutFlowHistManager->fillHistograms("<= 2 tight leptons", evtWeightRecorder.get(central_or_shift_main));

//--- apply HLT filter
    if(apply_hlt_filter)
    {
      const std::map<hltPathsE, bool> trigger_bits = {
        { hltPathsE::trigger_1e,    selTrigger_1e    },
        { hltPathsE::trigger_1mu,   selTrigger_1mu   },
        { hltPathsE::trigger_2e,    selTrigger_2e    },
        { hltPathsE::trigger_2mu,   selTrigger_2mu   },
        { hltPathsE::trigger_1e1mu, selTrigger_1e1mu },
      };
      if(! hltFilter(trigger_bits, fakeableLeptons, {}))
      {
        if(run_lumi_eventSelector || isDEBUG)
        {
          std::cout << "event " << eventInfo.str() << " FAILS HLT filter matching\n";
        }
        continue;
      }
    }
    cutFlowTable.update("HLT filter matching", evtWeightRecorder.get(central_or_shift_main));
    cutFlowHistManager->fillHistograms("HLT filter matching", evtWeightRecorder.get(central_or_shift_main));

    bool ttH_like = false;
    bool ttW_like = false;
    bool tH_like = false;
    if ((selBJets_loose.size() >= 2 || selBJets_medium.size() >= 1) && (selJets.size() >= 4)) ttH_like = true;
    if ((selBJets_loose.size() >= 2 || selBJets_medium.size() >= 1) && (selJets.size() == 3)) ttW_like = true;
    if (
      (selBJets_medium.size() >= 1 && ((selJets.size() - selBJets_loose.size()) + selJetsForward.size()) >= 1)
    ) tH_like = true;

    bool is_tH_like_and_not_ttH_like = false;
    if ( tH_like && !ttH_like ) is_tH_like_and_not_ttH_like = true;

    // apply requirement on jets (incl. b-tagged jets) and hadronic taus on level of final event selection
    bool passEvents = ttH_like || ttW_like || tH_like;
    if(do_sync) passEvents = ttH_like || tH_like;
    if ( !(passEvents) ) {
      if ( run_lumi_eventSelector ) {
        std::cout << "event " << eventInfo.str() << " FAILS Hadronic selection." << std::endl;
        //printCollection("selJets", selJets);
        std::cout << "event " << eventInfo.str()
          << "\n is_tH_like_and_not_ttH_like = " << is_tH_like_and_not_ttH_like
          << "\n selJets.size() = " << selJets.size()
          << "\n selBJets_medium.size() = " << selBJets_medium.size()
          << "\n selBJets_loose.size() = "<< selBJets_loose.size()
          << "\n selJetsForward.size()" << selJetsForward.size()
          << std::endl;
        printCollection("selJetsForward", selJetsForward);
      }
      continue;
    }
    cutFlowTable.update("Hadronic selection", evtWeightRecorder.get(central_or_shift_main));
    cutFlowHistManager->fillHistograms("Hadronic selection", evtWeightRecorder.get(central_or_shift_main));

    if ( selHadTaus.size() > 0 ) {
      if ( run_lumi_eventSelector ) {
        std::cout << "event " << eventInfo.str() << " FAILS selHadTaus veto." << std::endl;
        printCollection("selHadTaus", selHadTaus);
      }
      continue;
    }
    cutFlowTable.update("sel tau veto", evtWeightRecorder.get(central_or_shift_main));
    cutFlowHistManager->fillHistograms("sel tau veto", evtWeightRecorder.get(central_or_shift_main));

    const bool failsLowMassVeto = isfailsLowMassVeto(preselLeptonsFullUncleaned);
    if ( failsLowMassVeto  ) {
      if ( run_lumi_eventSelector || isDEBUG ) {
    std::cout << "event " << eventInfo.str() << " FAILS low mass lepton pair veto. " << std::endl;
      }
      continue;
    }
    cutFlowTable.update("m(ll) > 12 GeV", evtWeightRecorder.get(central_or_shift_main));
    cutFlowHistManager->fillHistograms("m(ll) > 12 GeV", evtWeightRecorder.get(central_or_shift_main));

    const double minPt_lead = 25.;
    const double minPt_sublead = 15.;
    // CV: according to Giovanni, the pT cuts should be applied on cone_pt
    //    (combined efficiency of single lepton and double lepton triggers assumed to be high,
    //     even if one or two leptons and fakes and hence cone_pt may be significantly smaller than lepton_pt,
    //     on which pT cuts are applied on trigger level)
    if ( !(selLepton_lead->cone_pt() > minPt_lead && selLepton_sublead->cone_pt() > minPt_sublead) ) {
      if ( run_lumi_eventSelector ) {
    std::cout << "event " << eventInfo.str() << " FAILS lepton pT selection." << std::endl;
	std::cout << " (leading selLepton pT = " << selLepton_lead->pt() << ", minPt_lead = " << minPt_lead
		  << ", subleading selLepton pT = " << selLepton_sublead->pt() << ", minPt_sublead = " << minPt_sublead << ")" << std::endl;
      }
      continue;
    }
    cutFlowTable.update("lead lepton pT > 25 GeV && sublead lepton pT > 15 GeV", evtWeightRecorder.get(central_or_shift_main));
    cutFlowHistManager->fillHistograms("lead lepton pT > 25 GeV && sublead lepton pT > 15 GeV", evtWeightRecorder.get(central_or_shift_main));

    bool failsTightChargeCut = false;
    for ( std::vector<const RecoLepton*>::const_iterator lepton = fakeableLeptons.begin();
	  lepton != fakeableLeptons.end(); ++lepton ) {
      if ( (*lepton)->is_electron() ) {
	const RecoElectron* electron = dynamic_cast<const RecoElectron*>(*lepton);
	assert(electron);
	if ( electron->tightCharge() < 2 ) failsTightChargeCut = true;
      }
      if ( (*lepton)->is_muon() ) {
	const RecoMuon* muon = dynamic_cast<const RecoMuon*>(*lepton);
	assert(muon);
	if ( muon->tightCharge() < 2 ) failsTightChargeCut = true;
      }
    }
    if ( failsTightChargeCut ) {
      if ( run_lumi_eventSelector ) {
    std::cout << "event " << eventInfo.str() << " FAILS tight lepton charge requirement." << std::endl;
      }
      //if (!selectBDT)
      continue;
    }
    cutFlowTable.update("tight lepton charge", evtWeightRecorder.get(central_or_shift_main));
    cutFlowHistManager->fillHistograms("tight lepton charge", evtWeightRecorder.get(central_or_shift_main));

    bool isLeptonCharge_SS = selLepton_lead->charge()*selLepton_sublead->charge() > 0;
    bool isLeptonCharge_OS = selLepton_lead->charge()*selLepton_sublead->charge() < 0;
    if ( leptonChargeSelection == kOS && isLeptonCharge_SS ) {
      if ( run_lumi_eventSelector ) {
    std::cout << "event " << eventInfo.str() << " FAILS lepton charge selection." << std::endl;
	std::cout << " (leading selLepton charge = " << selLepton_lead->charge()
		  << ", subleading selLepton charge = " << selLepton_sublead->charge() << ", leptonChargeSelection = OS)" << std::endl;
      }
      continue;
    }
    if ( leptonChargeSelection == kSS && isLeptonCharge_OS ) {
      if ( run_lumi_eventSelector ) {
    std::cout << "event " << eventInfo.str() << " FAILS lepton charge selection." << std::endl;
	std::cout << " (leading selLepton charge = " << selLepton_lead->charge()
		  << ", subleading selLepton charge = " << selLepton_sublead->charge() << ", leptonChargeSelection = SS)" << std::endl;
      }
      continue;
    }
    if(leptonChargeSelection == kOS)
    {
      const double prob_chargeMisID_sum = chargeMisIdRate.get(selLepton_lead, selLepton_sublead);
      // Karl: reject the event, if the applied probability of charge misidentification is 0;
      //       note that this can happen only if both selected leptons are muons (their misId prob is 0).
      if(prob_chargeMisID_sum == 0.)
      {
        if(run_lumi_eventSelector)
        {
          std::cout << "event " << eventInfo.str() << " FAILS charge flip selection\n"
                       "(leading lepton charge (pdgId) = " << selLepton_lead->charge() << " (" << selLepton_lead->pdgId()
                    << "); "
                       "subleading lepton charge (pdgId) = " << selLepton_sublead->charge() << " (" << selLepton_sublead->pdgId()
                    << "))\n"
          ;
        }
        continue;
      }
      evtWeightRecorder.record_chargeMisIdProb(prob_chargeMisID_sum);
    }
    cutFlowTable.update(Form("sel lepton-pair %s charge", leptonChargeSelection_string.data()), evtWeightRecorder.get(central_or_shift_main));
    cutFlowHistManager->fillHistograms("sel lepton-pair OS/SS charge", evtWeightRecorder.get(central_or_shift_main));

    const bool failsZbosonMassVeto = isfailsZbosonMassVeto(preselLeptonsFull) || (
        selLepton_lead->is_electron() &&
        selLepton_sublead->is_electron() &&
        isfailsZbosonMassVeto({ selLepton_lead, selLepton_sublead }, true)
      )
    ;
    if ( failsZbosonMassVeto ) {
      if ( run_lumi_eventSelector ) {
        std::cout << "event " << eventInfo.str() << " FAILS Z-boson veto. " << std::endl;
      }
      continue;
    }
    cutFlowTable.update("Z-boson mass veto", evtWeightRecorder.get(central_or_shift_main));
    cutFlowHistManager->fillHistograms("Z-boson mass veto", evtWeightRecorder.get(central_or_shift_main));
    if ( !(selLepton_lead->is_muon() || selLepton_sublead->is_muon() || met_LD >= 30.) ) {
      if ( run_lumi_eventSelector ) {
        std::cout << "event " << eventInfo.str() << " FAILS MET LD selection.\n"
                     " (LD = " << met_LD << ")\n"
        ;
      }
      continue;
    }
    cutFlowTable.update("met LD > 0.2", evtWeightRecorder.get(central_or_shift_main));
    cutFlowHistManager->fillHistograms("met LD > 0.2", evtWeightRecorder.get(central_or_shift_main));

    if ( apply_met_filters ) {
      if ( !metFilterSelector(metFilters) ) {
	if ( run_lumi_eventSelector ) {
      std::cout << "event " << eventInfo.str() << " FAILS MEt filters." << std::endl;
	}
	continue;
      }
    }
    cutFlowTable.update("MEt filters", evtWeightRecorder.get(central_or_shift_main));
    cutFlowHistManager->fillHistograms("MEt filters", evtWeightRecorder.get(central_or_shift_main));

    if (isMC) {
      if((selLepton_lead->genLepton() && selLepton_lead->charge() != selLepton_lead->genLepton()->charge()) ||
         (selLepton_sublead->genLepton() && selLepton_sublead->charge() != selLepton_sublead->genLepton()->charge())){
        if(run_lumi_eventSelector )
          {
	    std::cout << "event " << eventInfo.str() << " FAILS lepton-par gen=rec charge matching\n"
              "(leading lepton charge = " << selLepton_lead->charge() << " genlepton charge = " << selLepton_lead->genLepton()->charge()<< "; "
              "subleading lepton charge = " << selLepton_sublead->charge() << " genlepton charge = " << selLepton_sublead->genLepton()->charge()<< "\n"
              ;
          }
        continue;
      }
    }
    cutFlowTable.update("sel lepton-pair gen=rec charge match", evtWeightRecorder.get(central_or_shift_main));
    cutFlowHistManager->fillHistograms("sel lepton-pair gen=rec charge match", evtWeightRecorder.get(central_or_shift_main));



    bool failsSignalRegionVeto = false;
    if ( isMCClosure_e || isMCClosure_m ) {
      bool applySignalRegionVeto = (isMCClosure_e && countFakeElectrons(selLeptons) >= 1) || (isMCClosure_m && countFakeMuons(selLeptons) >= 1);
      if ( applySignalRegionVeto && tightLeptons.size() >= 2 ) failsSignalRegionVeto = true;
    } else if ( electronSelection == kFakeable || muonSelection == kFakeable ) {
      if ( tightLeptons.size() >= 2 ) failsSignalRegionVeto = true;
    }
    if ( failsSignalRegionVeto ) {
      if ( run_lumi_eventSelector ) {
	std::cout << "event " << eventInfo.str() << " FAILS overlap w/ the SR: "
	             "# tight leptons = " << tightLeptons.size() << " >= 2\n"
        ;
	printCollection("tightLeptons", tightLeptons);
      }
      continue; // CV: avoid overlap with signal region
    }

    cutFlowTable.update("signal region veto", evtWeightRecorder.get(central_or_shift_main));
    cutFlowHistManager->fillHistograms("signal region veto", evtWeightRecorder.get(central_or_shift_main));

    std::map<std::string, double> Weight_ktScan; // weights to do histograms for BMs
    if(apply_HH_rwgt_lo)
    {
      assert(HHWeightLO_calc);
      evtWeightRecorder.record_hhWeight_lo(HHWeightLO_calc, eventInfo, isDEBUG);
      for(const std::string & HHWeightName: evt_cat_strs)
      {
        Weight_ktScan[HHWeightName] = HHWeightLO_calc->getRelativeWeight(HHWeightName, eventInfo.gen_mHH, eventInfo.gen_cosThetaStar, isDEBUG);
      }
    }
    if(apply_HH_rwgt_nlo)
    {
      assert(HHWeightNLO_calc);
      evtWeightRecorder.record_hhWeight_nlo(HHWeightNLO_calc, eventInfo, isDEBUG);
      for(const std::string & HHWeightName: evt_cat_strs)
      {
        Weight_ktScan[HHWeightName] *= HHWeightNLO_calc->getRelativeWeight_LOtoNLO(HHWeightName, eventInfo.gen_mHH, eventInfo.gen_cosThetaStar, isDEBUG);
      }
    }

//--- compute variables BDTs used to discriminate ttH vs. ttV and ttH vs. ttbar -- they will be used more than once -- Xanda
    double mindr_lep1_jet=comp_mindr_jet(*selLepton_lead, selJets);
    double mindr_lep2_jet=comp_mindr_jet(*selLepton_sublead, selJets);
    const double max_lep_eta=TMath::Max(std::abs(selLepton_lead -> eta()), std::abs(selLepton_sublead -> eta()));
    double avg_dr_jet= comp_avg_dr_jet(selJets);
    double nJet25_Recl=comp_n_jet25_recl(selJets);
    double lep1_conePt=comp_lep_conePt(*selLepton_lead);
    double lep2_conePt=comp_lep_conePt(*selLepton_sublead);
    double minMET400=std::min(met.pt(), (Double_t)400.);
    double min_Deta_mostfwdJet_jet = 0;
    double min_Deta_leadfwdJet_jet = 0;
    // take the highest eta selJetsForward
    Particle::LorentzVector mostFwdJet = HighestEtaFwdJet(selJetsForward);
    if (selJetsForward.size() > 0 && selJets.size() > 0)
    {
      min_Deta_mostfwdJet_jet = min_Deta_fwdJet_jet(mostFwdJet, selJets);
      Particle::LorentzVector leadFwdJet = selJetsForward[0]-> p4();
      min_Deta_leadfwdJet_jet = min_Deta_fwdJet_jet(leadFwdJet, selJets);
    }

    double mT2_top_3particle = -1.;
    double mT2_top_2particle = -1.;
    double mT2_W = -1.;
    const Particle::LorentzVector & selLeptonP4_lead = selLepton_lead->cone_p4();
    const Particle::LorentzVector & selLeptonP4_sublead = selLepton_sublead->cone_p4();

    if(selJets.size() >= 2)
    {
      const Particle::LorentzVector & selJetP4_Hbb_lead = selJets[0]->p4();
      const Particle::LorentzVector & selJetP4_Hbb_sublead = selJets[1]->p4();
      mT2_2particle mT2Algo_2particle;
      const Particle::LorentzVector& metP4 = met.p4();
      mT2Algo_2particle(
        selLeptonP4_lead.px(), selLeptonP4_lead.py(), selLeptonP4_lead.mass(),
        selLeptonP4_sublead.px(), selLeptonP4_sublead.py(), selLeptonP4_sublead.mass(),
        metP4.px(), metP4.py(), 0.
      );
      mT2_W = mT2Algo_2particle.get_min_mT2();

      const double cSumPx = selLeptonP4_lead.px() + selLeptonP4_sublead.px() + metP4.px();
      const double cSumPy = selLeptonP4_lead.py() + selLeptonP4_sublead.py() + metP4.py();
      mT2Algo_2particle(
        selJetP4_Hbb_lead.px(), selJetP4_Hbb_lead.py(), selJetP4_Hbb_lead.mass(),
        selJetP4_Hbb_sublead.px(), selJetP4_Hbb_sublead.py(), selJetP4_Hbb_sublead.mass(),
        cSumPx, cSumPy, wBosonMass
      );
      mT2_top_2particle = mT2Algo_2particle.get_min_mT2();
      mT2_3particle mT2Algo_3particle;
      mT2Algo_3particle(
        selJetP4_Hbb_lead.px(), selJetP4_Hbb_lead.py(), selJetP4_Hbb_lead.mass(),
        selJetP4_Hbb_sublead.px(), selJetP4_Hbb_sublead.py(), selJetP4_Hbb_sublead.mass(),
        selLeptonP4_lead.px(), selLeptonP4_lead.py(), selLeptonP4_lead.mass(),
        selLeptonP4_sublead.px(), selLeptonP4_sublead.py(), selLeptonP4_sublead.mass(),
        metP4.px(), metP4.py(), 0.
      );
      const double mT2_top_3particle_permutation1 = mT2Algo_3particle.get_min_mT2();
      mT2Algo_3particle(
       selJetP4_Hbb_lead.px(), selJetP4_Hbb_lead.py(), selJetP4_Hbb_lead.mass(),
       selJetP4_Hbb_sublead.px(), selJetP4_Hbb_sublead.py(), selJetP4_Hbb_sublead.mass(),
       selLeptonP4_sublead.px(), selLeptonP4_sublead.py(), selLeptonP4_sublead.mass(),
       selLeptonP4_lead.px(), selLeptonP4_lead.py(), selLeptonP4_lead.mass(),
       metP4.px(), metP4.py(), 0.
      );
      const double mT2_top_3particle_permutation2 = mT2Algo_3particle.get_min_mT2();

      if(mT2_top_3particle_permutation1 <= mT2_top_3particle_permutation2)
      {
        mT2_top_3particle = mT2_top_3particle_permutation1;
      }
      else
      {
        mT2_top_3particle = mT2_top_3particle_permutation2;
      }
    }

//--- compute output of BDTs used to discriminate ttH vs. ttV and ttH vs. ttbar
//    in 2lss category of ttH multilepton analysis
    mvaInputs_2lss["max(abs(LepGood_eta[iF_Recl[0]]),abs(LepGood_eta[iF_Recl[1]]))"] = std::max(std::fabs(selLepton_lead->eta()), std::fabs(selLepton_sublead->eta()));
    mvaInputs_2lss["MT_met_lep1"]                = comp_MT_met(selLepton_lead, met.pt(), met.phi());
    mvaInputs_2lss["nJet25_Recl"]                = nJet25_Recl;
    mvaInputs_2lss["mindr_lep1_jet"]             = TMath::Min(10., mindr_lep1_jet);
    mvaInputs_2lss["mindr_lep2_jet"]             = TMath::Min(10., mindr_lep2_jet);
    mvaInputs_2lss["LepGood_conePt[iF_Recl[0]]"] = lep1_conePt;
    mvaInputs_2lss["LepGood_conePt[iF_Recl[1]]"] = lep2_conePt;
    mvaInputs_2lss["min(met_pt,400)"]            = minMET400;
    mvaInputs_2lss["avg_dr_jet"]                 = avg_dr_jet;
    check_mvaInputs(mvaInputs_2lss, eventInfo);
    double mvaOutput_2lss_ttV = mva_2lss_ttV(mvaInputs_2lss);
    double mvaOutput_2lss_ttbar = mva_2lss_ttbar(mvaInputs_2lss);
    //std::cout << "mvaOutput_2lss_ttbar = " << mvaOutput_2lss_ttbar << std::endl;

//--- compute integer discriminant based on both BDT outputs,
//    as defined in Table 16 () of AN-2015/321 (AN-2016/211) for analysis of 2015 (2016) data
    Double_t mvaDiscr_2lss = -1;
    if      ( mvaOutput_2lss_ttbar > +0.4 && mvaOutput_2lss_ttV >  +0.4 ) mvaDiscr_2lss = 7.;
    else if ( mvaOutput_2lss_ttbar > +0.4 && mvaOutput_2lss_ttV >  +0.1 ) mvaDiscr_2lss = 6.;
    else if ( mvaOutput_2lss_ttbar > +0.4 && mvaOutput_2lss_ttV <= +0.1 ) mvaDiscr_2lss = 5.;
    else if ( mvaOutput_2lss_ttbar > +0.1 && mvaOutput_2lss_ttV >  +0.3 ) mvaDiscr_2lss = 4.;
    else if ( mvaOutput_2lss_ttbar > +0.1 && mvaOutput_2lss_ttV <= +0.3 ) mvaDiscr_2lss = 3.;
    else if ( mvaOutput_2lss_ttbar > -0.2                               ) mvaDiscr_2lss = 2.;
    else                                                                  mvaDiscr_2lss = 1.;

    //--- compute output of hadronic top tagger BDT
    // it returns the gen-triplets organized in top/anti-top
    bool calculate_matching = isMC && selectBDT && !applyAdditionalEvtWeight; // DY has not matching info
    std::map<int, Particle::LorentzVector> genVar;
    std::map<int, Particle::LorentzVector> genVarAnti;
    if (calculate_matching) {
      genVar = isGenMatchedJetTripletVar(genTopQuarks, genBJets, genWBosons, genQuarkFromTop, kGenTop);
      genVarAnti = isGenMatchedJetTripletVar(genTopQuarks, genBJets, genWBosons, genQuarkFromTop, kGenAntiTop);
    }

    // resolved HTT

    double max_mvaOutput_HTT_CSVsort4rd = -9.;
    bool max_truth_HTT_CSVsort4rd = false;
    double HadTop_pt_CSVsort4rd = 0.;
    //double HadTop_eta_CSVsort4rd = 0.;
    double genTopPt_CSVsort4rd = 0.;
    double b_pt_CSVsort4rd_1 = 0.1;
    double Wj1_pt_CSVsort4rd_1 = 0.1;
    double Wj2_pt_CSVsort4rd_1 = 0.1;
    bool hadtruth = false;
    for ( std::vector<const RecoJet*>::const_iterator selBJet = selJets.begin(); selBJet != selJets.end(); ++selBJet ) {
      for ( std::vector<const RecoJet*>::const_iterator selWJet1 = selJets.begin(); selWJet1 != selJets.end(); ++selWJet1 ) {
        if ( &(*selWJet1) == &(*selBJet) ) continue;
        for ( std::vector<const RecoJet*>::const_iterator selWJet2 = selJets.begin(); selWJet2 != selJets.end(); ++selWJet2 ) {
          if ( &(*selWJet2) == &(*selBJet) ) continue;
          if ( &(*selWJet2) == &(*selWJet1) ) continue;
          bool isGenMatched = false;
          double genTopPt_teste = 0.;
          const std::map<int, double> bdtResult = (*hadTopTagger)(**selBJet, **selWJet1, **selWJet2, calculate_matching, isGenMatched, genTopPt_teste, genVar, genVarAnti );
          // genTopPt_teste is filled with the result of gen-matching
          if ( isGenMatched ) hadtruth = true;
          // save genpt of all options
          double HadTop_pt = ((*selBJet)->p4() + (*selWJet1)->p4() + (*selWJet2)->p4()).pt();

          if ( bdtResult.at(kXGB_CSVsort4rd) > max_mvaOutput_HTT_CSVsort4rd ) {
            max_truth_HTT_CSVsort4rd = isGenMatched;
            max_mvaOutput_HTT_CSVsort4rd = bdtResult.at(kXGB_CSVsort4rd);
            HadTop_pt_CSVsort4rd = HadTop_pt;
            genTopPt_CSVsort4rd = genTopPt_teste;
            //HadTop_eta_CSVsort4rd = std::fabs(((*selBJet)->p4() + (*selWJet1)->p4() + (*selWJet2)->p4()).eta());
            Wj1_pt_CSVsort4rd_1 = (*selWJet1)->pt();
            Wj2_pt_CSVsort4rd_1 = (*selWJet2)->pt();
            b_pt_CSVsort4rd_1   = (*selBJet)->pt();
          }
        }
      }
    }

    //std::map<std::string, double> mvaOutput_Hj_tagger;
    std::map<std::string, double> mvaInputs_Hj_tagger;
    double mvaOutput_Hj_tagger = -9.;
    for ( std::vector<const RecoJet*>::const_iterator selJet = selJets.begin();
	  selJet != selJets.end(); ++selJet ) {
      if ((*selJet)->pt()==Wj1_pt_CSVsort4rd_1 || (*selJet)->pt()==Wj2_pt_CSVsort4rd_1 || (*selJet)->pt()==b_pt_CSVsort4rd_1) continue;
      double mvaOutput = comp_mvaOutput_Hj_tagger(
        *selJet, fakeableLeptons, mvaInputs_Hj_tagger, mva_Hj_tagger,
        eventInfo);
      if ( mvaOutput > mvaOutput_Hj_tagger ) mvaOutput_Hj_tagger = mvaOutput;
    }

    evtWeightSum += evtWeightRecorder.get(central_or_shift_main);

    /*
    "lep1_conePt", "lep1_eta", "lep1_phi", "mT_lep1", "lep1_charge",
    "lep2_conePt", "lep2_eta", "lep2_phi", "mT_lep2",
    "mindr_lep1_jet", "mindr_lep2_jet",
    "Dilep_pdgId", "maxeta",
    "jet1_pt", "jet1_eta", "jet1_phi",
    "jet2_pt", "jet2_eta", "jet2_phi",
    "jet3_pt", "jet3_eta", "jet3_phi",
    "jet4_pt", "jet4_eta", "jet4_phi",
    "jetFwd1_pt", "jetFwd1_eta",
    "hadTop_BDT",
    "n_presel_jet",
    "avg_dr_jet",
    "Hj_tagger_hadTop","nBJetLoose",
    "mbb",
    "nBJetMedium", "n_presel_jetFwd",
    "metLD",
    */
    std::map<std::string, double> mvaInputVariables_NN_list = {
      {"lep1_conePt",     lep1_conePt},
      {"lep1_eta",        selLepton_lead -> eta()},
      {"lep1_phi",        selLepton_lead -> phi()},
      {"mT_lep1",         comp_MT_met(selLepton_lead, met.pt(), met.phi())}, //
      {"mindr_lep1_jet",  mindr_lep1_jet},
      {"lep1_charge",     selLepton_lead -> charge()},
      {"lep2_conePt",     lep2_conePt},
      {"lep2_eta",        selLepton_sublead -> eta()},
      {"lep2_phi",        selLepton_sublead -> phi()},
      {"mT_lep2",         comp_MT_met(selLepton_sublead, met.pt(), met.phi())}, //
      {"mindr_lep2_jet",  mindr_lep2_jet},
      {"Dilep_pdgId",     selElectrons.size() + 1},
      {"maxeta",          TMath::Max(selLepton_lead -> absEta(), selLepton_sublead -> absEta())},
      {"jetFwd1_eta",     selJetsForward.size() > 0 ? selJetsForward[0] -> absEta() : 9.},
      {"jetFwd1_pt",      selJetsForward.size() > 0 ? selJetsForward[0] -> pt()     : -9.},
      {"mbb",             selBJets_loose.size()>1 ?  (selBJets_loose[0]->p4()+selBJets_loose[1]->p4()).mass() : -9},
      {"avg_dr_jet",      selJets.size() > 1 ?  avg_dr_jet : -9},
      {"metLD",           met_LD > 0 ? met_LD : -9},
      {"hadTop_BDT",      max_mvaOutput_HTT_CSVsort4rd},
      {"n_presel_jet",    selJets.size()},
      {"n_presel_jetFwd", selJetsForward.size()},
      {"jet1_pt",         selJets[0]->pt()},
      {"jet1_eta",        selJets[0]->absEta()},
      {"jet1_phi",        selJets[0]->phi()},
      {"jet2_pt",         selJets.size() > 1 ?  selJets[1]->pt()     : -9.},
      {"jet2_eta",        selJets.size() > 1 ?  selJets[1]->absEta() :  9.},
      {"jet2_phi",        selJets.size() > 1 ?  selJets[1]->phi()    : -9.},
      {"jet3_pt",         selJets.size() > 2 ?  selJets[2]->pt()     : -9.},
      {"jet3_eta",        selJets.size() > 2 ?  selJets[2]->absEta() :  9.},
      {"jet3_phi",        selJets.size() > 2 ?  selJets[2]->phi()    : -9.},
      {"jet4_pt",         selJets.size() > 3 ?  selJets[3]->pt()     : -9.},
      {"jet4_eta",        selJets.size() > 3 ?  selJets[3]->absEta() :  9.},
      {"jet4_phi",        selJets.size() > 3 ?  selJets[3]->phi()    : -9.},
      {"nBJetLoose",      selBJets_loose.size()},
      {"nBJetMedium",     selBJets_medium.size()},
      {"Hj_tagger_hadTop", mvaOutput_Hj_tagger}
      };
    std::map<std::string, double> mvaOutput_NN = mva_NN(mvaInputVariables_NN_list);
    if ( isDEBUG_TF ) {
      std::cout << "event " << eventInfo.str() << "\n";
      std::cout << "Variables :\n";
      for(auto elem : mvaInputVariables_NN_list) std::cout << elem.first << " " << elem.second << "\n";
      std::cout << "\n";
      std::cout << "result :";
      for(auto elem : mvaOutput_NN) std::cout << elem.first << " " << elem.second << " ";
      std::cout << "\n";
    }

//--- do NN categories
    std::string category_NN = "output_NN_";
    double output_NN = -10;
    // "predictions_ttH","predictions_Rest","predictions_ttW","predictions_tHQ"
    if (passEvents) {
      if (
        mvaOutput_NN["predictions_ttH"] >= mvaOutput_NN["predictions_tHQ"] &&
        mvaOutput_NN["predictions_ttH"] >= mvaOutput_NN["predictions_ttW"] &&
        mvaOutput_NN["predictions_ttH"] >= mvaOutput_NN["predictions_Rest"]
      ) {
        category_NN += "ttH";
        output_NN = mvaOutput_NN["predictions_ttH"];
      }
      if (
        mvaOutput_NN["predictions_tHQ"] > mvaOutput_NN["predictions_ttH"] &&
        mvaOutput_NN["predictions_tHQ"] > mvaOutput_NN["predictions_ttW"] &&
        mvaOutput_NN["predictions_tHQ"] >= mvaOutput_NN["predictions_Rest"]
      ) {
        category_NN += "tH";
        output_NN = mvaOutput_NN["predictions_tHQ"];
      }
      if (
        mvaOutput_NN["predictions_ttW"] > mvaOutput_NN["predictions_ttH"] &&
        mvaOutput_NN["predictions_ttW"] >= mvaOutput_NN["predictions_tHQ"] &&
        mvaOutput_NN["predictions_ttW"] >= mvaOutput_NN["predictions_Rest"]
      ) {
        category_NN += "ttW";
        output_NN = mvaOutput_NN["predictions_ttW"];
      }
      if (
        mvaOutput_NN["predictions_Rest"] > mvaOutput_NN["predictions_ttH"] &&
        mvaOutput_NN["predictions_Rest"] > mvaOutput_NN["predictions_tHQ"] &&
        mvaOutput_NN["predictions_Rest"] > mvaOutput_NN["predictions_ttW"]
      ) {
        category_NN += "rest";
        output_NN = mvaOutput_NN["predictions_Rest"];
      }

      if  ( ( selLepton_lead_type == kElectron && selLepton_sublead_type == kElectron ) ) {
        category_NN += "_ee";
      } else if (  selLepton_lead_type == kMuon     && selLepton_sublead_type == kMuon      ) {
        category_NN += "_mm";
      } else if ( (selLepton_lead_type == kElectron && selLepton_sublead_type == kMuon    ) ||
      (selLepton_lead_type == kMuon     && selLepton_sublead_type == kElectron) ) {
        category_NN += "_em";
      }

    }

    ///////////////////////////////
    // SVA variables
    const double mass_2L           = (selLepton_lead->cone_p4() + selLepton_sublead->cone_p4()).mass();
    const int    sum_Lep_charge    = selLepton_lead -> charge() + selLepton_sublead -> charge();
    std::string category_SVA = "mass_2L_";
    if (! is_tH_like_and_not_ttH_like) {
      if ( selJets.size() > 3 )
        {
          if  ( ( selLepton_lead_type == kElectron && selLepton_sublead_type == kElectron ) ) {
            category_SVA += "ee";
          } else if (  selLepton_lead_type == kMuon     && selLepton_sublead_type == kMuon      ) {
            category_SVA += "mm";
          } else if ( (selLepton_lead_type == kElectron && selLepton_sublead_type == kMuon    ) ||
          (selLepton_lead_type == kMuon     && selLepton_sublead_type == kElectron) ) {
            category_SVA += "em";
          }
          if (selJets.size() < 6) category_SVA += "_lj";
          else category_SVA += "_hj";
          if (sum_Lep_charge > 0 ) category_SVA += "_pos";
          else category_SVA += "_neg";
        } else category_SVA += "cr";
    } else category_SVA += "rest";

//--- retrieve gen-matching flags
    std::vector<const GenMatchEntry*> genMatches = genMatchInterface.getGenMatch(selLeptons);

//--- fill histograms with events passing final selection
    std::map<std::string, double> tH_weight_map_main;
    for(const std::string & central_or_shift: central_or_shifts_local)
    {
      const double evtWeight = evtWeightRecorder.get(central_or_shift);
      const bool skipFilling = central_or_shift != central_or_shift_main;
      for (const GenMatchEntry* genMatch : genMatches)
      {
        selHistManagerType* selHistManager = selHistManagers[central_or_shift][genMatch->getIdx()];
        assert(selHistManager);
        if(! skipFilling)
        {
          selHistManager->electrons_->fillHistograms(selElectrons, evtWeight);
          selHistManager->muons_->fillHistograms(selMuons, evtWeight);
          selHistManager->hadTaus_->fillHistograms(selHadTaus, evtWeight);
          selHistManager->jets_->fillHistograms(selJets, evtWeight);
          selHistManager->leadJet_->fillHistograms(selJets, evtWeight);
          selHistManager->subleadJet_->fillHistograms(selJets, evtWeight);
          selHistManager->BJets_loose_->fillHistograms(selBJets_loose, evtWeight);
          selHistManager->leadBJet_loose_->fillHistograms(selBJets_loose, evtWeight);
          selHistManager->subleadBJet_loose_->fillHistograms(selBJets_loose, evtWeight);
          selHistManager->BJets_medium_->fillHistograms(selBJets_medium, evtWeight);
          selHistManager->met_->fillHistograms(met, mht_p4, met_LD, evtWeight);
          selHistManager->metFilters_->fillHistograms(metFilters, evtWeight);
          selHistManager->mvaInputVariables_2lss_->fillHistograms(mvaInputVariables_NN_list, evtWeight);
        }

        const std::string central_or_shift_tH = eventInfo.has_central_or_shift(central_or_shift) ? central_or_shift : central_or_shift_main;
        const double evtWeight_tH_nom = evtWeightRecorder.get_nom_tH_weight(central_or_shift_tH);
        std::map<std::string, double> tH_weight_map;
        for(const std::string & evt_cat_str: evt_cat_strs)
        {
          if(skipFilling && evt_cat_str != default_cat_str)
          {
            continue;
          }
          if(isMC_tH)
          {
            const std::string evt_cat_str_query = evt_cat_str == default_cat_str ? get_tH_SM_str() : evt_cat_str;
            tH_weight_map[evt_cat_str] = evtWeight / evtWeight_tH_nom * eventInfo.genWeight_tH(central_or_shift_tH, evt_cat_str_query);
          }
          else if(apply_HH_rwgt_lo)
          {
            tH_weight_map[evt_cat_str] = evtWeight * Weight_ktScan[evt_cat_str];
          }
          else
          {
            tH_weight_map[evt_cat_str] = evtWeight;
          }
        }
        if(central_or_shift == central_or_shift_main)
        {
          tH_weight_map_main = tH_weight_map;
        }

        const int htxs_category = eventInfo.get_htxs_category();
        for(const auto & kv: tH_weight_map)
        {
          const EvtHistManager_2lss_Input fillVariables {
            selElectrons.size(),
            selMuons.size(),
            selHadTaus.size(),
            selJets.size(),
            selBJets_loose.size(),
            selBJets_medium.size(),
            kv.second,
            mvaOutput_2lss_ttV,
            mvaOutput_2lss_ttbar,
            mvaDiscr_2lss,
            mvaOutput_Hj_tagger,
            output_NN,
            category_NN,
            mass_2L,
            category_SVA,
          };
          EvtHistManager_2lss* selHistManager_evt = selHistManager->evt_[kv.first];
          if ( selHistManager_evt )
          {
            selHistManager_evt->fillHistograms(fillVariables);
          }
          for(const auto & kw: evt_htxs_binning)
          {
            if(htxs_category & kw.second)
            {
              EvtHistManager_2lss_Input fillVariables_htxs = fillVariables;
              fillVariables_htxs.evtWeight = evtWeightRecorder.get(central_or_shift, kw.first);
              selHistManager->evt_htxs_[kv.first][kw.first]->fillHistograms(fillVariables_htxs);
            }
          }

          if ( isSignal ) {
            std::string decayModeStr = get_key_hist(eventInfo, genWBosons, isMC_HH, isMC_VH);
            if ( ( isMC_tH || isMC_H ) && ( decayModeStr == "hzg" || decayModeStr == "hmm" ) ) continue;
            if(! decayModeStr.empty())
            {
              EvtHistManager_2lss* selHistManager_evt_decay = selHistManager->evt_in_decayModes_[kv.first][decayModeStr];
              if ( selHistManager_evt_decay )
              {
                selHistManager_evt_decay->fillHistograms(fillVariables);
              }
              for(const auto & kw: evt_htxs_binning)
              {
                if(htxs_category & kw.second)
                {
                  EvtHistManager_2lss_Input fillVariables_htxs = fillVariables;
                  fillVariables_htxs.evtWeight = evtWeightRecorder.get(central_or_shift, kw.first);
                  selHistManager->evt_htxs_in_decayModes_[kv.first][decayModeStr][kw.first]->fillHistograms(fillVariables_htxs);
                }
              }
            }
          }
        }
        if(! skipFilling)
        {
          selHistManager->evtYield_->fillHistograms(eventInfo, evtWeight);
          selHistManager->weights_->fillHistograms("genWeight", eventInfo.genWeight);
          selHistManager->weights_->fillHistograms("pileupWeight", evtWeightRecorder.get_puWeight(central_or_shift));
          selHistManager->weights_->fillHistograms("triggerWeight", evtWeightRecorder.get_sf_triggerEff(central_or_shift));
          selHistManager->weights_->fillHistograms("data_to_MC_correction", evtWeightRecorder.get_data_to_MC_correction(central_or_shift));
          selHistManager->weights_->fillHistograms("fakeRate", evtWeightRecorder.get_FR(central_or_shift));
        }
      }

      if(isMC && ! skipFilling)
      {
        genEvtHistManager_afterCuts[central_or_shift]->fillHistograms(
          genElectrons, genMuons, genHadTaus, genPhotons, genJets, evtWeightRecorder.get_inclusive(central_or_shift)
        );
        lheInfoHistManager[central_or_shift]->fillHistograms(*lheInfoReader, evtWeight);
        if(eventWeightManager)
        {
          genEvtHistManager_afterCuts[central_or_shift]->fillHistograms(
            eventWeightManager, evtWeightRecorder.get_inclusive(central_or_shift)
          );
        }
      }
    }

    if ( selEventsFile ) {
      (*selEventsFile) << eventInfo.str() << '\n';
    }

    bool isGenMatched = false;
    if ( isMC )
    {
      for (const GenMatchEntry* genMatch : genMatches)
      {
	if ( genMatch->getName() == "" ) isGenMatched = true; // non-fake
      }
    }

    if ( bdt_filler ) {

      const double lep1_frWeight = selLepton_lead->isGenMatched(true) ? 1. : evtWeightRecorder.get_jetToLepton_FR_lead(central_or_shift_main);
      const double lep2_frWeight = selLepton_sublead->isGenMatched(true) ? 1. : evtWeightRecorder.get_jetToLepton_FR_sublead(central_or_shift_main);
      const double evt_frWeight = lep1_frWeight * lep2_frWeight;

      double lep1_genLepPt = selLepton_lead->genLepton() ? selLepton_lead->genLepton()->pt() : 0.;
      double lep2_genLepPt = selLepton_sublead->genLepton() ? selLepton_sublead->genLepton()->pt() : 0.;

      bdt_filler -> operator()({ eventInfo.run, eventInfo.lumi, eventInfo.event })
          ("lep1_pt",                selLepton_lead -> pt())
          ("lep1_conePt",            comp_lep_conePt(*selLepton_lead))
          ("lep1_eta",               selLepton_lead -> eta())
          ("lep1_phi",               selLepton_lead -> phi())
          ("lep1_tth_mva",           selLepton_lead -> mvaRawTTH())
          ("mindr_lep1_jet",         TMath::Min(10., mindr_lep1_jet) )
          ("mindr_lep2_jet",         TMath::Min(10., mindr_lep2_jet) )
          ("mT_lep1",                comp_MT_met(selLepton_lead, met.pt(), met.phi()))
          ("MT_met_lep1",            comp_MT_met(selLepton_lead, met.pt(), met.phi()))
          ("lep2_pt",                selLepton_sublead -> pt())
          ("lep2_conePt",            comp_lep_conePt(*selLepton_sublead))
          ("lep2_eta",               selLepton_sublead -> eta())
          ("lep2_phi",               selLepton_sublead -> phi())
          ("max_Lep_eta",            TMath::Max(std::abs(selLepton_lead -> eta()), std::abs(selLepton_sublead -> eta())))
          ("avg_dr_lep",             1.0) // comp_avg_dr_jet(selLeptons))
          ("lep2_tth_mva",           selLepton_sublead -> mvaRawTTH())
          ("mT_lep2",                comp_MT_met(selLepton_sublead, met.pt(), met.phi()))
          ("avg_dr_jet",             comp_avg_dr_jet(selJets))
          ("nJet25_Recl",            comp_n_jet25_recl(selJets))
          ("ptmiss",                 met.pt())
          ("htmiss",                 mht_p4.pt())
          ("dr_leps",                deltaR(selLepton_lead -> p4(), selLepton_sublead -> p4()))
          ("lep1_genLepPt",          lep1_genLepPt)
          ("lep2_genLepPt",          lep2_genLepPt)
          ("lep1_frWeight",          lep1_frWeight)
          ("lep2_frWeight",          lep2_frWeight)
          ("mvaOutput_2lss_ttV",     mvaOutput_2lss_ttV)
          ("mvaOutput_2lss_ttbar",   mvaOutput_2lss_ttbar)
          ("mvaDiscr_2lss",          mvaDiscr_2lss)
          ("Hj_tagger",              mvaOutput_Hj_tagger)
          ("lumiScale",              evtWeightRecorder.get_lumiScale(central_or_shift_main))
          ("genWeight",              eventInfo.genWeight)
          ("evtWeight",              evtWeightRecorder.get(central_or_shift_main) * evt_frWeight)
          ("nJet",                   selJets.size())
          ("nLep",                   selLeptons.size())
          ("nBJetLoose",             selBJets_loose.size())
          ("nBJetMedium",            selBJets_medium.size())
          ("lep1_isTight",           int(selLepton_lead -> isTight()))
          ("lep2_isTight",           int(selLepton_sublead -> isTight()))

          ("hadtruth",               hadtruth)
          ("bWj1Wj2_isGenMatched_CSVsort4rd",              max_truth_HTT_CSVsort4rd)
          ("res_HTT",                 max_mvaOutput_HTT_CSVsort4rd)
          ("HadTop_pt",               HadTop_pt_CSVsort4rd)
          ("genTopPt_CSVsort4rd",             genTopPt_CSVsort4rd)

          ("min(met_pt,400)",            std::min(met.pt(), (Double_t)400.))
          ("mbb_loose",       selBJets_loose.size()>1 ?  (selBJets_loose[0]->p4()+selBJets_loose[1]->p4()).mass() : 0  )
          ("failsTightChargeCut",          failsTightChargeCut)

          ("mbb_medium",          selBJets_medium.size()>1 ?  (selBJets_medium[0]->p4()+selBJets_medium[1]->p4()).mass() : 0 )
          ("nElectron",                      selElectrons.size())
          ("sum_Lep_charge", selLepton_lead -> charge() + selLepton_sublead -> charge())
          ("massLT",          selLeptons.size() > 1 ? comp_massL2(selLeptons[0], selLeptons[1], met.pt(), met.phi())  : 0.)
          ("massL",           massL(fakeableLeptons))
          ("min_Deta_mostfwdJet_jet", min_Deta_mostfwdJet_jet)
          ("min_Deta_leadfwdJet_jet", min_Deta_leadfwdJet_jet)
          ("met_LD",              met_LD)
          ("jet1_pt",   selJets.size() > 0 ? selJets[0]->pt() : -1000)
          ("jet1_eta",  selJets.size() > 0 ? selJets[0]->eta() : -1000)
          ("jet1_phi",  selJets.size() > 0 ? selJets[0]->phi() : -1000)
          ("jet1_E",    selJets.size() > 0 ? selJets[0]->p4().energy() : -1000)
          ("jet2_pt",   selJets.size() > 1 ? selJets[1]->pt() : -1000)
          ("jet2_eta",  selJets.size() > 1 ? selJets[1]->eta() : -1000)
          ("jet2_phi",  selJets.size() > 1 ? selJets[1]->phi() : -1000)
          ("jet2_E",    selJets.size() > 1 ? selJets[1]->p4().energy() : -1000)
          ("jet3_pt",   selJets.size() > 2 ? selJets[2]->pt() : -1000)
          ("jet3_eta",  selJets.size() > 2 ? selJets[2]->eta() : -1000)
          ("jet3_phi",  selJets.size() > 2 ? selJets[2]->phi() : -1000)
          ("jet3_E",    selJets.size() > 2 ? selJets[2]->p4().energy() : -1000)
          ("jet4_pt",   selJets.size() > 3 ? selJets[3]->pt() : -1000)
          ("jet4_eta",  selJets.size() > 3 ? selJets[3]->eta() : -1000)
          ("jet4_phi",  selJets.size() > 3 ? selJets[3]->phi() : -1000)
          ("jet4_E",    selJets.size() > 3 ? selJets[3]->p4().energy() : -1000)
          ("mostFwdJet_eta",      selJetsForward.size() > 0 ? std::abs(mostFwdJet.Eta()) : -1000)
          ("mostFwdJet_pt",       selJetsForward.size() > 0 ? mostFwdJet.pt() : -1000)
          ("mostFwdJet_phi",      selJetsForward.size() > 0 ? mostFwdJet.phi() : -1000)
          ("mostFwdJet_E",        selJetsForward.size() > 0 ? mostFwdJet.energy() : -1000)
          ("leadFwdJet_eta",      selJetsForward.size() > 0 ? selJetsForward[0] -> absEta() : -1000)
          ("leadFwdJet_pt",       selJetsForward.size() > 0 ? selJetsForward[0] -> pt() : -1000)
          ("leadFwdJet_phi",      selJetsForward.size() > 0 ? selJetsForward[0] -> phi() : -1000)
          ("leadFwdJet_E",        selJetsForward.size() > 0 ? selJetsForward[0] -> p4().energy() : -1000)
          ("nJetForward",         selJetsForward.size())
          (tH_weight_map_main)
        .fill()
      ;
    }

    if(snm)
    {
      const double mT_lep1        = comp_MT_met(selLepton_lead, met.pt(), met.phi());
      const double mT_lep2        = comp_MT_met(selLepton_sublead, met.pt(), met.phi());
      const double max_dr_jet     = comp_max_dr_jet(selJets);
      const double mbb            = selBJets_medium.size() > 1 ? (selBJets_medium[0]->p4() + selBJets_medium[1]->p4()).mass() : -1.;
      const double mbb_loose      = selBJets_loose.size() > 1 ? (selBJets_loose[0]->p4() + selBJets_loose[1]->p4()).mass() : -1.;
      const double min_dr_lep_jet = std::min(mindr_lep1_jet, mindr_lep2_jet);
      const double dr_leps        = deltaR(selLepton_lead->p4(), selLepton_sublead->p4());
      const int nLightJet         = selJets.size() - selBJets_loose.size() + selJetsForward.size();

      snm->read(eventInfo);
      snm->read(selLeptons);
      snm->read(preselMuons,     fakeableMuons,     tightMuons);
      snm->read(preselElectrons, fakeableElectrons, tightElectrons);
      snm->read(looseHadTaus);
      snm->read(selJets, false);
      snm->read(selJetsForward, true);

      snm->read({ triggers_1e, triggers_2e, triggers_1mu, triggers_2mu, triggers_1e1mu });
      snm->read(
        isGenMatched, selBJets_medium.size(), selBJets_loose.size(), nLightJet,
        SyncNtupleManager::placeholder_value, is_tH_like_and_not_ttH_like
      );

      snm->read(met.pt(),                               FloatVariableType::PFMET);
      snm->read(met.phi(),                              FloatVariableType::PFMETphi);
      snm->read(mht_p4.pt(),                            FloatVariableType::MHT);
      snm->read(met_LD,                                 FloatVariableType::metLD);

      snm->read(mindr_lep1_jet,                         FloatVariableType::mindr_lep1_jet);
      snm->read(mindr_lep2_jet,                         FloatVariableType::mindr_lep2_jet);
      // mindr_lep3_jet not filled
      // mindr_lep4_jet not filled

      // mindr_tau1_jet not filled
      // mindr_tau2_jet not filled

      snm->read(avg_dr_jet,                             FloatVariableType::avg_dr_jet);
      // avr_dr_lep_tau not filled
      snm->read(max_dr_jet,                             FloatVariableType::max_dr_jet);
      // max_dr_lep_tau not filled
      // min_dr_tau_jet not filled
      // min_dr_lep_tau not filled
      snm->read(min_dr_lep_jet,                         FloatVariableType::min_dr_lep_jet);

      snm->read(dr_leps,                                FloatVariableType::dr_leps);
      // dr_taus not filled

      // dr_lep_tau_ss not filled
      // dr_lep1_tau1 not filled
      // dr_lep1_tau2 not filled
      // dr_lep2_tau1 not filled
      // dr_lep3_tau1 not filled
      // dr_lep2_tau2 not filled

      snm->read(max_lep_eta,                            FloatVariableType::max_lep_eta);

      snm->read(mT_lep1,                                FloatVariableType::mT_met_lep1);
      snm->read(mT_lep2,                                FloatVariableType::mT_met_lep2);
      // mT_met_lep3 not filled
      // mT_met_lep4 not filled
      snm->read(massL(fakeableLeptons),                 FloatVariableType::massL);

      snm->read(mT2_W,                                  FloatVariableType::mT2_W);
      snm->read(mT2_top_2particle,                      FloatVariableType::mT2_top_2particle);
      snm->read(mT2_top_3particle,                      FloatVariableType::mT2_top_3particle);

      // mTauTauVis not filled
      // mvis_l1tau not filled
      // mvis_l2tau not filled

      snm->read(mbb,                                    FloatVariableType::mbb);
      snm->read(mbb_loose,                              FloatVariableType::mbb_loose);

      // cosThetaS_hadTau not filled
      // HTT not filled
      // HadTop_pt not filled
      // Hj_tagger not filled

      // mvaOutput_plainKin_ttV not filled
      // mvaOutput_plainKin_tt not filled
      // mvaOutput_plainKin_1B_VT not filled
      // mvaOutput_HTT_SUM_VT not filled

      // mvaOutput_plainKin_SUM_VT not filled

      // mvaOutput_2lss_1tau_plainKin_tt not filled
      // mvaOutput_2lss_1tau_plainKin_ttV not filled
      // mvaOutput_2lss_1tau_plainKin_1B_M not filled
      // mvaOutput_2lss_1tau_plainKin_SUM_M not filled
      // mvaOutput_2lss_1tau_HTT_SUM_M not filled
      // mvaOutput_2lss_1tau_HTTMEM_SUM_M not filled

      // mvaOutput_3l_ttV not filled
      // mvaOutput_3l_ttbar not filled
      // mvaOutput_plainKin_SUM_M not filled
      // mvaOutput_plainKin_1B_M not filled

      snm->read(evtWeightRecorder.get_FR(central_or_shift_main),             FloatVariableType::FR_weight);
      snm->read(evtWeightRecorder.get_sf_triggerEff(central_or_shift_main),  FloatVariableType::triggerSF_weight);
      snm->read(evtWeightRecorder.get_leptonSF(),                            FloatVariableType::leptonSF_weight);
      // tauSF_weight not filled
      snm->read(evtWeightRecorder.get_btag(central_or_shift_main),           FloatVariableType::bTagSF_weight);
      snm->read(evtWeightRecorder.get_puWeight(central_or_shift_main),       FloatVariableType::PU_weight);
      snm->read(evtWeightRecorder.get_genWeight(),                           FloatVariableType::MC_weight);

      // mvaOutput_2lss_ttV not filled
      // mvaOutput_2lss_tt not filled
      // mvaOutput_2lss_1tau_plainKin_tt not filled
      // mvaOutput_2lss_1tau_plainKin_ttV not filled
      // mvaOutput_2lss_1tau_plainKin_1B_M not filled
      // mvaOutput_2lss_1tau_plainKin_SUM_M not filled
      // mvaOutput_2lss_1tau_HTT_SUM_M not filled
      // mvaOutput_2lss_1tau_HTTMEM_SUM_M not filled

      snm->read(mvaOutput_NN["predictions_ttH"],  FloatVariableType::mvaOutput_2lss_0tau_ttH);
      snm->read(mvaOutput_NN["predictions_tHQ"],   FloatVariableType::mvaOutput_2lss_0tau_tH);
      snm->read(mvaOutput_NN["predictions_ttW"],  FloatVariableType::mvaOutput_2lss_0tau_ttW);
      snm->read(mvaOutput_NN["predictions_Rest"], FloatVariableType::mvaOutput_2lss_0tau_rest);

      // mvaOutput_3l_ttH_tH_3cat_v8_ttH not filled
      // mvaOutput_3l_ttH_tH_3cat_v8_tH not filled
      // mvaOutput_3l_ttH_tH_3cat_v8_rest not filled

      snm->read(eventInfo.genWeight,                    FloatVariableType::genWeight);

      if(isGenMatched)
      {
        snm->fill();
      }
      else
      {
        snm->resetBranches();
      }
    }

    ++selectedEntries;
    selectedEntries_weighted += evtWeightRecorder.get(central_or_shift_main);
    std::string process_and_genMatch = process_string;
    process_and_genMatch += selLepton_genMatch.name_;
    ++selectedEntries_byGenMatchType[process_and_genMatch];
    for(const std::string & central_or_shift: central_or_shifts_local)
    {
      selectedEntries_weighted_byGenMatchType[central_or_shift][process_and_genMatch] += evtWeightRecorder.get(central_or_shift);
    }
    histogram_selectedEntries->Fill(0.);
    if(isDEBUG)
    {
      std::cout << evtWeightRecorder << '\n';
    }

    //---------------------------------------------------------------------------
    // CV: fill a few extra histogram for the HIG-19-008 paper
    if ( genMatches.size() > 0 )
    {
      bool isFake = false;
      for ( auto genMatch : genMatches )
      {
        if ( genMatch->getName() == "_fake" ) isFake = true;
      }
      TH1* histogram_genBJetPt = nullptr;
      TH1* histogram_genBJetEta = nullptr;
      TH1* histogram_numGenBJets = nullptr;
      TH1* histogram_numRecBJets_loose = nullptr;
      TH1* histogram_numRecBJets_medium = nullptr;
      if ( isFake )
      {
        histogram_genBJetPt = histogram_genBJetPt_fake;
        histogram_genBJetEta = histogram_genBJetEta_fake;
        histogram_numGenBJets = histogram_numGenBJets_fake;
        histogram_numRecBJets_loose = histogram_numRecBJets_loose_fake;
        histogram_numRecBJets_medium = histogram_numRecBJets_medium_fake;
      }
      else
      {
        histogram_genBJetPt = histogram_genBJetPt_nonfake;
        histogram_genBJetEta = histogram_genBJetEta_nonfake;
        histogram_numGenBJets = histogram_numGenBJets_nonfake;
        histogram_numRecBJets_loose = histogram_numRecBJets_loose_nonfake;
        histogram_numRecBJets_medium = histogram_numRecBJets_medium_nonfake;
      }
      double evtWeight = evtWeightRecorder.get(central_or_shift_main);
      int numGenBJets = 0;
      for ( auto genBJet : genBJets )
      {
        if ( std::fabs(genBJet.eta()) <  2.4 ) histogram_genBJetPt->Fill(genBJet.pt(), evtWeight);
        if (           genBJet.pt()   > 25.  ) histogram_genBJetEta->Fill(genBJet.eta(), evtWeight);
        if ( std::fabs(genBJet.eta()) <  2.4 &&
                       genBJet.pt()   > 25.  ) ++numGenBJets;
      }
      histogram_numGenBJets->Fill(numGenBJets, evtWeight);
      histogram_numRecBJets_loose->Fill(selBJets_loose.size(), evtWeight);
      histogram_numRecBJets_medium->Fill(selBJets_medium.size(), evtWeight);
    }
    //---------------------------------------------------------------------------
  }

  if(snm)
  {
    snm->write();
  }

  std::cout << "max num. Entries = " << inputTree -> getCumulativeMaxEventCount()
            << " (limited by " << maxEvents << ") processed in "
            << inputTree -> getProcessedFileCount() << " file(s) (out of "
            << inputTree -> getFileCount() << ")\n"
            << " analyzed = " << analyzedEntries << '\n'
            << " selected = " << selectedEntries << " (weighted = " << selectedEntries_weighted << ")\n\n"
            << "cut-flow table" << std::endl;
  cutFlowTable.print(std::cout);
  std::cout << std::endl;

  std::cout << "sel. Entries by gen. matching:" << std::endl;
  for(const std::string & central_or_shift: central_or_shifts_local)
  {
    std::cout << "central_or_shift = " << central_or_shift << '\n';
    for(const leptonChargeFlipGenMatchEntry & leptonGenMatch_definition: leptonGenMatch_definitions)
    {
      std::string process_and_genMatch = process_string;
      process_and_genMatch += leptonGenMatch_definition.name_;
      std::cout << " " << process_and_genMatch << " = " << selectedEntries_byGenMatchType[process_and_genMatch]
                << " (weighted = " << selectedEntries_weighted_byGenMatchType[central_or_shift][process_and_genMatch] << ")\n";
    }
  }
  std::cout << std::endl;
  std::cout << "Sum of weights "<< evtWeightSum << std::endl;

//--- manually write histograms to output file
  fs.file().cd();
  //histogram_analyzedEntries->Write();
  //histogram_selectedEntries->Write();
  HistManagerBase::writeHistograms();

//--- memory clean-up
  delete dataToMCcorrectionInterface;

  delete leptonFakeRateInterface;

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
  delete psWeightReader;

  delete hadTopTagger;

  for(auto & kv: genEvtHistManager_beforeCuts)
  {
    delete kv.second;
  }
  for(auto & kv: genEvtHistManager_afterCuts)
  {
    delete kv.second;
  }
  for(auto & kv: lheInfoHistManager)
  {
    delete kv.second;
  }
  delete l1PreFiringWeightReader;
  delete cutFlowHistManager;
  delete eventWeightManager;

  delete inputTree;
  delete snm;

  clock.Show("analyze_2lss");

  return EXIT_SUCCESS;
}
