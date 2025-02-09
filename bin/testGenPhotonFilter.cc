#include "tthAnalysis/HiggsToTauTau/interface/TTreeWrapper.h" // TTreeWrapper
#include "tthAnalysis/HiggsToTauTau/interface/GenPhotonReader.h" // GenPhotonReader, GenPhoton, GenParticle
#include "tthAnalysis/HiggsToTauTau/interface/GenParticleReader.h" // GenParticleReader
#include "tthAnalysis/HiggsToTauTau/interface/GenLeptonReader.h" // GenLeptonReader, GenLepton
#include "tthAnalysis/HiggsToTauTau/interface/EventInfo.h" // EventInfo
#include "tthAnalysis/HiggsToTauTau/interface/EventInfoReader.h" // EventInfoReader
#include "tthAnalysis/HiggsToTauTau/interface/LHEParticleReader.h" // LHEParticleReader, LHEParticle
#include "tthAnalysis/HiggsToTauTau/interface/RunLumiEventSelector.h" // RunLumiEventSelector
#include "tthAnalysis/HiggsToTauTau/interface/GenPhotonFilter.h" // GenPhotonFilter
#include "tthAnalysis/HiggsToTauTau/interface/EvtWeightRecorder.h" // EvtWeightRecorder
#include "tthAnalysis/HiggsToTauTau/interface/EvtWeightManager.h" // EvtWeightManager
#include "tthAnalysis/HiggsToTauTau/interface/HistManagerBase.h" // HistManagerBase
#include "tthAnalysis/HiggsToTauTau/interface/histogramAuxFunctions.h" // fillWithOverFlow()

#include "PhysicsTools/FWLite/interface/TFileService.h" // fwlite::TFileService
#include "DataFormats/FWLite/interface/InputSource.h" // fwlite::InputSource
#include "DataFormats/FWLite/interface/OutputFiles.h" // fwlite::OutputFiles

#if __has_include (<FWCore/ParameterSetReader/interface/ParameterSetReader.h>)
#  include <FWCore/ParameterSetReader/interface/ParameterSetReader.h> // edm::readPSetsFrom()
#else
#  include <FWCore/PythonParameterSet/interface/MakeParameterSets.h> // edm::readPSetsFrom()
#endif

#include <TError.h> // gErrorAbortLevel, kError
#include <TBenchmark.h> // TBenchmark

class GenPhotonFilterHistManager
  : public HistManagerBase
{
public:
  GenPhotonFilterHistManager(const edm::ParameterSet & cfg)
    : HistManagerBase(cfg)
    , histogram_analyzed_(nullptr)
    , histogram_selected_subleading_(nullptr)
    , histogram_subleading_lepton_pt_(nullptr)
    , histogram_subleading_electron_pt_(nullptr)
    , histogram_subleading_muon_pt_(nullptr)
    , histogram_selected_third_(nullptr)
    , histogram_third_lepton_pt_(nullptr)
    , histogram_third_electron_pt_(nullptr)
    , histogram_third_muon_pt_(nullptr)
  {
    const std::vector<std::string> sysOpts = {
      "analyzed",
      "selected_subleading",
      "subleading_lepton_pt",
      "subleading_electron_pt",
      "subleading_muon_pt",
      "selected_third",
      "third_lepton_pt",
      "third_electron_pt",
      "third_muon_pt",
    };
    for(const std::string & sysOpt: sysOpts)
    {
      central_or_shiftOptions_[sysOpt] = { "central" };
    }
  }
  ~GenPhotonFilterHistManager() {}

  void
  bookHistograms(TFileDirectory & dir) override
  {
    histogram_analyzed_ = book1D(dir, "analyzed", "analyzed", 1, -0.5, +0.5);

    histogram_selected_subleading_    = book1D(dir, "selected_subleading",    "selected_subleading",      1,  -0.5,   +0.5);
    histogram_subleading_lepton_pt_   = book1D(dir, "subleading_lepton_pt",   "subleading_lepton_pt",   100,   0.,   100.);
    histogram_subleading_electron_pt_ = book1D(dir, "subleading_electron_pt", "subleading_electron_pt", 100,   0.,   100.);
    histogram_subleading_muon_pt_     = book1D(dir, "subleading_muon_pt",     "subleading_muon_pt",     100,   0.,   100.);

    histogram_selected_third_    = book1D(dir, "selected_third",    "selected_third",      1,  -0.5,   +0.5);
    histogram_third_lepton_pt_   = book1D(dir, "third_lepton_pt",   "third_lepton_pt",   100,   0.,   100.);
    histogram_third_electron_pt_ = book1D(dir, "third_electron_pt", "third_electron_pt", 100,   0.,   100.);
    histogram_third_muon_pt_     = book1D(dir, "third_muon_pt",     "third_muon_pt",     100,   0.,   100.);
  }

  void
  fillHistograms(const std::vector<GenLepton> & genLeptons,
                 double evtWeight)
  {
    const double evtWeightErr = 0.;

    fillWithOverFlow(histogram_analyzed_, 0., evtWeight, evtWeightErr);

    if(genLeptons.size() < 2)
    {
      return;
    }
    assert(std::all_of(
      genLeptons.cbegin(), genLeptons.cend(),
      [](const GenLepton & genLepton) -> bool
      {
        return
          (genLepton.checkStatusFlag(StatusFlag::isPrompt) || genLepton.checkStatusFlag(StatusFlag::isPromptTauDecayProduct)) &&
          genLepton.status() == 1
        ;
      }
    ));

    fillWithOverFlow(histogram_selected_subleading_, 0., evtWeight, evtWeightErr);

    const GenLepton subleading_lepton = genLeptons.at(1);
    const int subleading_lepton_absPdgId = std::abs(subleading_lepton.pdgId());
    const double subleading_lepton_pt = subleading_lepton.pt();

    fillWithOverFlow(histogram_subleading_lepton_pt_, subleading_lepton_pt, evtWeight, evtWeightErr);

    if(subleading_lepton_absPdgId == 11)
    {
      fillWithOverFlow(histogram_subleading_electron_pt_, subleading_lepton_pt, evtWeight, evtWeightErr);
    }
    else if(subleading_lepton_absPdgId == 13)
    {
      fillWithOverFlow(histogram_subleading_muon_pt_, subleading_lepton_pt, evtWeight, evtWeightErr);
    }
    else
    {
      throw cmsException(this, __func__, __LINE__)
        << "Invalid PDG ID = " << subleading_lepton_absPdgId << " of generator-level lepton: " << subleading_lepton
      ;
    }

    if(genLeptons.size() < 3)
    {
      return;
    }

    fillWithOverFlow(histogram_selected_third_, 0., evtWeight, evtWeightErr);

    const GenLepton third_lepton = genLeptons.at(2);
    const int third_lepton_absPdgId = std::abs(third_lepton.pdgId());
    const double third_lepton_pt = third_lepton.pt();

    fillWithOverFlow(histogram_third_lepton_pt_, third_lepton_pt, evtWeight, evtWeightErr);

    if(third_lepton_absPdgId == 11)
    {
      fillWithOverFlow(histogram_third_electron_pt_, third_lepton_pt, evtWeight, evtWeightErr);
    }
    else if(third_lepton_absPdgId == 13)
    {
      fillWithOverFlow(histogram_third_muon_pt_, third_lepton_pt, evtWeight, evtWeightErr);
    }
    else
    {
      throw cmsException(this, __func__, __LINE__)
        << "Invalid PDG ID = " << third_lepton_absPdgId << " of generator-level lepton: " << third_lepton
      ;
    }
  }

private:
  TH1 * histogram_analyzed_;

  TH1 * histogram_selected_subleading_;
  TH1 * histogram_subleading_lepton_pt_;
  TH1 * histogram_subleading_electron_pt_;
  TH1 * histogram_subleading_muon_pt_;

  TH1 * histogram_selected_third_;
  TH1 * histogram_third_lepton_pt_;
  TH1 * histogram_third_electron_pt_;
  TH1 * histogram_third_muon_pt_;
};

std::string
getCategory(const LHEParticleReader * const lheParticleReader)
{
  assert(lheParticleReader);
  const std::vector<LHEParticle> lheParticles = lheParticleReader->read();
  const std::size_t nofLeptons = std::count_if(
    lheParticles.cbegin(), lheParticles.cend(),
    [](const LHEParticle & lheParticle) -> bool
    {
      const int lheParticleAbsPdgId = std::abs(lheParticle.pdgId());
      return lheParticleAbsPdgId == 11 || lheParticleAbsPdgId == 13 || lheParticleAbsPdgId == 15;
    }
  );
  if(nofLeptons == 0)
  {
    return "Hadronic";
  }
  else if(nofLeptons == 1)
  {
    return "SemiLept";
  }
  else if(nofLeptons == 2)
  {
    return "DiLept";
  }
  throw cmsException(__func__, __LINE__) << "Unexpected number of leptons found in the event: " << nofLeptons;
}

int
main(int argc,
     char ** argv)
{
//--- throw an exception in case ROOT encounters an error
  gErrorAbortLevel = kError;

//--- parse command-line arguments
  if(argc < 2)
  {
    std::cout << "Usage: " << argv[0] << " [parameters.py]\n";
    return EXIT_FAILURE;
  }

  std::cout << '<' << argv[0] << ">:\b";

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start(argv[0]);

//--- read python configuration parameters
  if(! edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process"))
  {
    throw cmsException(argv[0]) << "No ParameterSet 'process' found in configuration file = " << argv[1];
  }

  const edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");
  const edm::ParameterSet cfg_analyze = cfg.getParameter<edm::ParameterSet>("analyze_testGenPhotonFilter");
  const std::string process_string = cfg_analyze.getParameter<std::string>("process");
  const std::string era_string = cfg_analyze.getParameter<std::string>("era");
  const std::string treeName = cfg_analyze.getParameter<std::string>("treeName");
  const AnalysisConfig analysisConfig("testGenPhotonFilter", cfg_analyze);

  const bool isMC = cfg_analyze.getParameter<bool>("isMC");
  const bool isDEBUG = cfg_analyze.getParameter<bool>("isDEBUG");
  const bool apply_genWeight = cfg_analyze.getParameter<bool>("apply_genWeight");
  const double ref_genWeight = cfg_analyze.getParameter<double>("ref_genWeight");

  const std::string branchName_genLeptons = cfg_analyze.getParameter<std::string>("branchName_genLeptons");
  const std::string branchName_genPhotons = cfg_analyze.getParameter<std::string>("branchName_genPhotons");
  const std::string branchName_genProxyPhotons = cfg_analyze.getParameter<std::string>("branchName_genProxyPhotons");
  const std::string branchName_genFromHardProcess = cfg_analyze.getParameter<std::string>("branchName_genFromHardProcess");

  const std::string central_or_shift_main = cfg_analyze.getParameter<std::string>("central_or_shift");
  const std::vector<std::string> central_or_shifts_local = { central_or_shift_main };
  const edm::VParameterSet lumiScale = cfg_analyze.getParameter<edm::VParameterSet>("lumiScale");

  const GenPhotonFilter genPhotonFilter("enabled");

  const std::string selEventsFileName_input = cfg_analyze.getParameter<std::string>("selEventsFileName_input");
  std::cout << "selEventsFileName_input = " << selEventsFileName_input << '\n';
  const RunLumiEventSelector * run_lumi_eventSelector = nullptr;
  if(! selEventsFileName_input.empty())
  {
    edm::ParameterSet cfgRunLumiEventSelector;
    cfgRunLumiEventSelector.addParameter<std::string>("inputFileName", selEventsFileName_input);
    cfgRunLumiEventSelector.addParameter<std::string>("separator", ":");
    run_lumi_eventSelector = new RunLumiEventSelector(cfgRunLumiEventSelector);
  }

  const fwlite::InputSource inputFiles(cfg);
  const int maxEvents = inputFiles.maxEvents();
  std::cout << " maxEvents = " << maxEvents << '\n';
  const unsigned reportEvery = inputFiles.reportAfter();

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());

  TTreeWrapper * const inputTree = new TTreeWrapper(treeName.data(), inputFiles.files(), maxEvents);
  std::cout << "Loaded " << inputTree -> getFileCount() << " file(s).\n";

  assert(isMC);

  EventInfo eventInfo(analysisConfig);
  eventInfo.set_refGetWeight(ref_genWeight);
  EventInfoReader eventInfoReader(&eventInfo);
  inputTree -> registerReader(&eventInfoReader);

  const edm::ParameterSet additionalEvtWeight = cfg_analyze.getParameter<edm::ParameterSet>("evtWeight");
  const bool applyAdditionalEvtWeight = additionalEvtWeight.getParameter<bool>("apply");
  EvtWeightManager * eventWeightManager = nullptr;
  if(applyAdditionalEvtWeight)
  {
    eventWeightManager = new EvtWeightManager(additionalEvtWeight);
    eventWeightManager->set_central_or_shift(central_or_shift_main);
    inputTree->registerReader(eventWeightManager);
  }

  GenLeptonReader * const genLeptonReader = new GenLeptonReader(branchName_genLeptons);
  inputTree -> registerReader(genLeptonReader);

  GenPhotonReader * const genPhotonReader = new GenPhotonReader(branchName_genPhotons);
  inputTree -> registerReader(genPhotonReader);

  GenPhotonReader * const genProxyPhotonReader = new GenPhotonReader(branchName_genProxyPhotons);
  inputTree -> registerReader(genProxyPhotonReader);

  GenParticleReader * const genFromHardProcessReader = new GenParticleReader(branchName_genFromHardProcess);
  inputTree -> registerReader(genFromHardProcessReader);

  LHEParticleReader * lheParticleReader = nullptr;
  std::vector<std::string> process_strings = { process_string };
  if(process_string == "TTGamma")
  {
    // split TTGamma into 3 subcategories
    for(const std::string & suffix: { "DiLept", "SemiLept", "Hadronic" })
    {
      process_strings.push_back(Form("%s_%s", process_string.data(), suffix.data()));
    }
    // initialize LHEParticleReader so that we can find out if the top-antitop pair decays DL, SL or hadronically
    lheParticleReader = new LHEParticleReader();
    inputTree -> registerReader(lheParticleReader);
  }

  std::map<std::string, GenPhotonFilterHistManager *> genPhotonFilterHistManager_selected;
  std::map<std::string, GenPhotonFilterHistManager *> genPhotonFilterHistManager_rejected;
  for(const std::string & cat: process_strings)
  {
    genPhotonFilterHistManager_selected[cat] = new GenPhotonFilterHistManager(makeHistManager_cfg(
      cat, "genPhotonFilter/selected", era_string, central_or_shift_main
    ));
    genPhotonFilterHistManager_selected[cat]->bookHistograms(fs);

    genPhotonFilterHistManager_rejected[cat] = new GenPhotonFilterHistManager(makeHistManager_cfg(
      cat, "genPhotonFilter/rejected", era_string, central_or_shift_main
    ));
    genPhotonFilterHistManager_rejected[cat]->bookHistograms(fs);
  }
  int analyzedEntries = 0;
  int selectedEntries = 0;
  double selectedEntries_weighted = 0.;

  while(inputTree->hasNextEvent() && (! run_lumi_eventSelector || (run_lumi_eventSelector && ! run_lumi_eventSelector -> areWeDone())))
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

    if(isDEBUG)
    {
      std::cout << "event #" << inputTree -> getCurrentMaxEventIdx() << ' ' << eventInfo << '\n';
    }

    if( run_lumi_eventSelector && !(*run_lumi_eventSelector)(eventInfo))
    {
      continue;
    }

    if(run_lumi_eventSelector)
    {
      std::cout << "processing Entry #" << inputTree->getCumulativeMaxEventCount() << ": " << eventInfo << '\n';
      if(inputTree -> isOpen())
      {
        std::cout << "input File = " << inputTree->getCurrentFileName() << '\n';
      }
    }

    EvtWeightRecorder evtWeightRecorder(central_or_shifts_local, central_or_shift_main, isMC);
    if(apply_genWeight)
    {
      evtWeightRecorder.record_genWeight(eventInfo);
    }
    if(eventWeightManager)
    {
      evtWeightRecorder.record_auxWeight(eventWeightManager);
    }
    evtWeightRecorder.record_lumiScale(lumiScale);

    const std::vector<GenLepton> genLeptons = genLeptonReader->read();
    const std::vector<GenPhoton> genPhotons = genPhotonReader->read(true);
    const std::vector<GenPhoton> genProxyPhotons = genProxyPhotonReader->read(true);
    const std::vector<GenParticle> genFromHardProcess = genFromHardProcessReader->read();

    std::vector<std::string> catToFill = { process_string };
    if(lheParticleReader)
    {
      const std::string cat = Form("%s_%s", process_string.data(), getCategory(lheParticleReader).data());
      catToFill.push_back(cat);
    }
    for(const std::string & cat: catToFill)
    {
      if(genPhotonFilter(genPhotons, genProxyPhotons, genFromHardProcess))
      {
        genPhotonFilterHistManager_selected[cat]->fillHistograms(genLeptons, evtWeightRecorder.get(central_or_shift_main));
      }
      else
      {
        genPhotonFilterHistManager_rejected[cat]->fillHistograms(genLeptons, evtWeightRecorder.get(central_or_shift_main));
      }
    }

    ++selectedEntries;
    selectedEntries_weighted += evtWeightRecorder.get(central_or_shift_main);
    if(isDEBUG)
    {
      std::cout << evtWeightRecorder << '\n';
    }
  }

  std::cout
    << "max num. Entries = " << inputTree -> getCumulativeMaxEventCount()
    << " (limited by " << maxEvents << ") processed in "
    << inputTree -> getProcessedFileCount() << " file(s) (out of "
    << inputTree -> getFileCount() << ")\n"
    << " analyzed = " << analyzedEntries << '\n'
    << " selected = " << selectedEntries << " (weighted = " << selectedEntries_weighted << ")\n"
  ;

  delete inputTree;
  delete run_lumi_eventSelector;
  delete genPhotonReader;
  delete genProxyPhotonReader;
  delete genFromHardProcessReader;
  delete lheParticleReader;
  delete eventWeightManager;

  clock.Show(argv[0]);

  return EXIT_SUCCESS;
}
