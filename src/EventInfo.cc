#include "tthAnalysis/HiggsToTauTau/interface/EventInfo.h" // EventInfo

#include "tthAnalysis/HiggsToTauTau/interface/weightAuxFunctions.h" // get_tH_weight_str()
#include "tthAnalysis/HiggsToTauTau/interface/cmsException.h" // cmsException()

#include <boost/range/adaptor/map.hpp> // boost::adaptors::map_keys
#include <boost/range/algorithm/copy.hpp> // boost::copy()

#include <cmath> // std::fabs()
#include <cstring> // std::memcpy()
#include <cassert> // assert()

#define EPS 1E-2

const std::map<std::string, Int_t> EventInfo::decayMode_idString_singleHiggs =
{
  { "hww",     24 },
  { "hzz",     23 },
  { "htt",     15 },
  { "hzg", 220023 },
  { "hmm",     13 },
  { "hbb",      5 },
};

const std::map<std::string, Int_t> EventInfo::decayMode_idString_diHiggs =
{
  { "tttt",       15 },
  { "zzzz",       23 },
  { "wwww",       24 },
  { "bbtt",  5000015 },
  { "bbzz",  5000023 },
  { "bbww",  5000024 },
  { "ttzz", 15000023 },
  { "ttww", 15000024 },
  { "zzww", 23000024 },
};

const std::map<Int_t, std::string> EventInfo::productionMode_idString_singleHiggs = {
  { 23000025, "ZH" },
  { 24000025, "WH" },
};

EventInfo::EventInfo()
  : EventInfo(false, false, false, false, false)
{}

EventInfo::EventInfo(const AnalysisConfig & analysisConfig)
  : EventInfo(analysisConfig.isMC(),
              analysisConfig.isMC_H(),
              analysisConfig.isMC_HH(),
              analysisConfig.isHH_rwgt_allowed(),
              analysisConfig.apply_topPtReweighting())
{
  process_string_ = analysisConfig.process();
  std::cout << "EventInfo::EventInfo():: process_string_: " << process_string_ << "\n";
}

EventInfo::EventInfo(bool isMC,
                     bool isMC_H,
                     bool isMC_HH,
                     bool isHH_rwgt_allowed,
                     bool apply_topPtRwgt)
  : run(0)
  , lumi(0)
  , event(0)
  , genHiggsDecayMode(-1)
  , genWeight(1.)
  , pileupWeight(1.)
  , genDiHiggsDecayMode(-1)
  , gen_mHH(0.)
  , gen_cosThetaStar(-2.)
  , topPtRwgtSF(-1.)
  , isMC_(isMC)
  , isMC_H_(isMC_H)
  , isMC_HH_(isMC_HH)
  , isHH_rwgt_allowed_(isHH_rwgt_allowed)
  , apply_topPtRwgt_(apply_topPtRwgt)
  , central_or_shift_("central")
  , nLHEReweightingWeight(0)
  , LHEReweightingWeight(nullptr)
  , LHEReweightingWeight_max(69)
  , is_owner(false)
  , read_htxs_(false)
  , refGenWeight_(0.)
  , productionMode_(-1)
{
  int checksum = 0;
  if ( isMC_H_  ) ++checksum;
  if ( isMC_HH_ ) ++checksum;
  if ( !(checksum == 0 || (checksum == 1 && isMC_)) )
    throw cmsException(this, __func__, __LINE__)
      << "Invalid combination of Configuration parameters:\n" 
      << " isMC = " << isMC_ << "\n" 
      << " isMC_H = " << isMC_H_ << "\n"  
      << " isMC_HH = " << isMC_HH_ << "\n";
}

EventInfo::EventInfo(const EventInfo & eventInfo)
  : EventInfo()
{
  copy(eventInfo);
}

EventInfo &
EventInfo::operator=(const EventInfo & eventInfo)
{
  copy(eventInfo);
  return *this;
}

void
EventInfo::copy(const EventInfo & eventInfo)
{
  run = eventInfo.run;
  lumi = eventInfo.lumi;
  event = eventInfo.event;
  genHiggsDecayMode = eventInfo.genHiggsDecayMode;
  genWeight = eventInfo.genWeight;
  isMC_H_ = eventInfo.isMC_H_;
  isMC_ = eventInfo.isMC_;
  genDiHiggsDecayMode = eventInfo.genDiHiggsDecayMode;
  isMC_HH_ = eventInfo.isMC_HH_;
  isHH_rwgt_allowed_ = eventInfo.isHH_rwgt_allowed_;
  gen_mHH = eventInfo.gen_mHH;
  gen_cosThetaStar = eventInfo.gen_cosThetaStar;
  apply_topPtRwgt_ = eventInfo.apply_topPtRwgt_;
  topPtRwgtSF = eventInfo.topPtRwgtSF;

  nLHEReweightingWeight = eventInfo.nLHEReweightingWeight;
  if(eventInfo.LHEReweightingWeight)
  {
    if(LHEReweightingWeight)
    {
      if(is_owner)
      {
        delete[] LHEReweightingWeight;
        LHEReweightingWeight = nullptr;
      }
      else
      {
        throw cmsException(this, __func__, __LINE__)
          << "Cannot overwrite with instance" << eventInfo << " because this instance doesn't own LHEReweightingWeight"
        ;
      }
    }
    LHEReweightingWeight = new Float_t[LHEReweightingWeight_max];
    is_owner = true;
    std::memcpy(LHEReweightingWeight, eventInfo.LHEReweightingWeight, nLHEReweightingWeight);
  }

  htxs_ = eventInfo.htxs_;
}

EventInfo::~EventInfo()
{}

bool
EventInfo::isMC() const
{
  return isMC_;
}

bool
EventInfo::isMC_H() const
{
  return isMC_H_;
}

bool
EventInfo::isMC_HH() const
{
  return isMC_HH_;
}

double
EventInfo::genWeight_tH() const
{
  return genWeight_tH(central_or_shift_);
}

double
EventInfo::genWeight_tH(double kv,
                        double kt,
                        double cosa) const
{
  return genWeight_tH(central_or_shift_, get_tH_weight_str(kv, kt, cosa));
}

double
EventInfo::genWeight_tH(const std::string & central_or_shift) const
{
  return genWeight_tH(central_or_shift, get_tH_SM_str());
}

double
EventInfo::genWeight_tH(const std::string & central_or_shfit,
                        const std::string & name) const
{
  if(tH_sf.empty())
  {
    return 1.;
  }
  if(! tH_sf.count(central_or_shfit))
  {
    throw cmsException(this, __func__, __LINE__) << "No such central or shift option: " << central_or_shfit;
  }
  if(! tH_sf.at(central_or_shfit).count(name))
  {
    throw cmsException(this, __func__, __LINE__)
      << "Invalid coupling requested for central or shift option " << central_or_shfit << ": " << name
    ;
  }
  const std::pair<int, double> & settings = tH_sf.at(central_or_shfit).at(name);
  assert(settings.first < static_cast<int>(nLHEReweightingWeight));

  return (settings.first >= 0 ? LHEReweightingWeight[settings.first] : 1.) * settings.second;
}

void
EventInfo::loadWeight_tH(const std::vector<edm::ParameterSet> & cfg)
{
  std::map<std::string, int> main_map;
  for(const edm::ParameterSet & cfg_entry: cfg)
  {
    const double weight = cfg_entry.getParameter<double>("weight");
    const double kt = cfg_entry.getParameter<double>("kt");
    const double kv = cfg_entry.getParameter<double>("kv");
    const double cosa = cfg_entry.exists("cosa") ?
      cfg_entry.getParameter<double>("cosa")     :
      std::numeric_limits<double>::quiet_NaN()
    ;
    const std::string central_or_shift = cfg_entry.getParameter<std::string>("central_or_shift");
    const int idx = cfg_entry.getParameter<int>("idx");
    const std::string name = get_tH_weight_str(kt, kv, cosa);

    // make sure that we can recover the paramters encoded by the string
    const std::tuple<double, double, double> kt_kv_cosa = get_tH_params(name);
    assert(std::get<0>(kt_kv_cosa) == kt);
    assert(std::get<1>(kt_kv_cosa) == kv);
    if(! std::isnan(cosa))
    {
      assert(std::get<2>(kt_kv_cosa) == cosa);
    }
    if(! has_central_or_shift(central_or_shift))
    {
      tH_sf[central_or_shift] = {};
    }

    assert(! tH_sf.at(central_or_shift).count(name));
    tH_sf[central_or_shift][name] = std::pair<int, double>(idx, weight);
    std::cout
      << "Loaded weight '" << name << "': kt = " << kt << " & kv = " << kv << " & cosa = " << cosa
      << " for systematic option " << central_or_shift
      << " -> weight = " << weight <<" @ idx = " << idx << '\n'
    ;
  }
}

std::vector<std::string>
EventInfo::getWeight_tH_str(bool include_sm) const
{
  return getWeight_tH_str(central_or_shift_, include_sm);
}

std::vector<std::string>
EventInfo::getWeight_tH_str(const std::string & central_or_shift,
                            bool include_sm) const
{
  std::vector<std::string> names;
  const std::string sm_str = get_tH_SM_str();
  if(! has_central_or_shift(central_or_shift))
  {
    throw cmsException(this, __func__, __LINE__) << "No such central or shift option: " << central_or_shift;
  }
  for(const auto & kv: tH_sf.at(central_or_shift))
  {
    if(! include_sm && kv.first == sm_str)
    {
      continue;
    }
    names.push_back(kv.first);
  }
  return names;
}

void
EventInfo::set_refGetWeight(double refGenWeight)
{
  refGenWeight_ = std::fabs(refGenWeight);
}

double
EventInfo::get_refGenWeight() const
{
  return refGenWeight_;
}

void
EventInfo::set_central_or_shift(const std::string & central_or_shift)
{
  central_or_shift_ = central_or_shift;
}

bool
EventInfo::has_central_or_shift(const std::string & central_or_shift) const
{
  return tH_sf.count(central_or_shift);
}

bool
EventInfo::read_htxs() const
{
  return read_htxs_;
}

void
EventInfo::read_htxs(bool flag)
{
  read_htxs_ = flag;
}

int
EventInfo::get_htxs_category() const
{
  return read_htxs_ ? htxs_.get_category() : HTXSCategory::kUnknown;
}

bool
EventInfo::is_initialized() const
{
  return run != 0 && lumi != 0 && event != 0;
}

std::string
EventInfo::getDecayModeString() const
{
  if(! isMC_H_)
  {
    throw cmsException(this, __func__, __LINE__)
      << "The event " << *this << " is not a H signal event => request "
         "for a decay mode as a string is not applicable\n"
    ;
  }
  return EventInfo::getDecayModeString(decayMode_idString_singleHiggs);
}

std::string
EventInfo::getDiHiggsDecayModeString() const
{
  if(! (isMC_HH_))
  {
    throw cmsException(this, __func__, __LINE__)
      << "The event " << *this << " is not a HH signal event => request "
         "for a decay mode as a string is not applicable\n"
    ;
  }
  return EventInfo::getDecayModeString(decayMode_idString_diHiggs);
}

std::string
EventInfo::getProductionModeString() const
{
  if(! isMC_H_)
  {
    throw cmsException(this, __func__, __LINE__)
      << "The event " << *this << " is not a H signal event => request "
         "for a production mode as a string is not applicable\n"
    ;
  }
  return
    productionMode_idString_singleHiggs.count(productionMode_) ?
    productionMode_idString_singleHiggs.at(productionMode_)    :
    ""
  ;
}

std::string
EventInfo::getDecayModeString(const std::map<std::string, Int_t> & decayMode_idString) const
{
  for(const auto & kv: decayMode_idString)
  {
    if(genHiggsDecayMode == kv.second)
    {
      return kv.first;
    }
  }
  return ""; // unrecognizable decay mode
}

std::vector<std::string>
EventInfo::getDecayModes()
{
  return getDecayModes(decayMode_idString_singleHiggs);
}

std::vector<std::string>
EventInfo::getDiHiggsDecayModes()
{
  return getDecayModes(decayMode_idString_diHiggs);
}

std::vector<std::string>
EventInfo::getDecayModes(const std::map<std::string, Int_t> & decayMode_idString)
{
  std::vector<std::string> decayModes;
  decayModes.reserve(decayMode_idString.size());
  boost::copy(decayMode_idString | boost::adaptors::map_keys, std::back_inserter(decayModes));
  return decayModes;
}

void
EventInfo::reset_productionMode()
{
  set_productionMode(-1);
}

void
EventInfo::set_productionMode(int productionMode)
{
  productionMode_ = productionMode;
}

std::string
EventInfo::str() const
{
  std::stringstream ss;
  ss << run << ':' << lumi << ':' << event;
  return ss.str();
}

std::ostream &
operator<<(std::ostream & os,
           const EventInfo & info)
{
  os << "run = " << info.run << ", ls = " << info.lumi << ", event = " << info.event;
  return os;
}
