#include "tthAnalysis/HiggsToTauTau/interface/RecoMEtReader.h"

#include "tthAnalysis/HiggsToTauTau/interface/analysisAuxFunctions.h"     // getBranchName_MEt()
#include "tthAnalysis/HiggsToTauTau/interface/cmsException.h"             // cmsException()
#include "tthAnalysis/HiggsToTauTau/interface/BranchAddressInitializer.h" // BranchAddressInitializer, TTree, Form()
#include "tthAnalysis/HiggsToTauTau/interface/sysUncertOptions.h"         // kJetMET_*
#include "tthAnalysis/HiggsToTauTau/interface/metPhiModulation.h"         // METXYCorr_Met_MetPhi()

std::map<std::string, int> RecoMEtReader::numInstances_;
std::map<std::string, RecoMEtReader *> RecoMEtReader::instances_;

RecoMEtReader::RecoMEtReader(Era era,
                             bool isMC)
  : RecoMEtReader(era, isMC, "MET")
{}

RecoMEtReader::RecoMEtReader(Era era,
                             bool isMC,
                             const std::string & branchName_obj,
                             const std::string & branchName_cov)
  : era_(era)
  , isMC_(isMC)
  , branchName_obj_(branchName_obj)
  , branchName_cov_(branchName_cov.empty() ? branchName_obj_ : branchName_cov)
  , eventInfo_(nullptr)
  , recoVertex_(nullptr)
  , enable_phiModulationCorr_(false)
  , ptPhiOption_(isMC_ ? kJetMET_central : kJetMET_central_nonNominal)
  , read_ptPhi_systematics_(false)
{
  setBranchNames();
}

RecoMEtReader::~RecoMEtReader()
{
  --numInstances_[branchName_obj_];
  assert(numInstances_[branchName_obj_] >= 0);
  if(numInstances_[branchName_obj_] == 0)
  {
    RecoMEtReader * const gInstance = instances_[branchName_obj_];
    assert(gInstance);
    instances_[branchName_obj_] = nullptr;
  }
}

void
RecoMEtReader::setMEt_central_or_shift(int central_or_shift)
{
  if(! isMC_ && central_or_shift != kJetMET_central_nonNominal)
  {
    throw cmsException(this, __func__, __LINE__) << "Nominal MET available only in MC";
  }
  if(! isValidJESsource(era_, central_or_shift))
  {
    throw cmsException(this, __func__, __LINE__) << "Invalid option for the era = " << static_cast<int>(era_) << ": " << central_or_shift;
  }
  if(central_or_shift <= kJetMET_UnclusteredEnDown)
  {
    ptPhiOption_ = central_or_shift;
  }
  else
  {
    std::cout
        << get_human_line(this, __func__, __LINE__)
        << "Not setting the systematics option to " << central_or_shift
        << " but keeping it at " << ptPhiOption_ << std::endl;
    ;
  }
}

void
RecoMEtReader::read_ptPhi_systematics(bool flag)
{
  if(! isMC_ && flag)
  {
    throw cmsException(this, __func__, __LINE__)
      << "Cannot read MET systematics from data"
    ;
  }
  read_ptPhi_systematics_ = flag;
}

void
RecoMEtReader::setBranchNames()
{
  if(numInstances_[branchName_obj_] == 0)
  {
    for(int met_option = kJetMET_central_nonNominal; met_option <= kJetMET_UnclusteredEnDown; ++met_option)
    {
      if(! isValidJESsource(era_, met_option))
      {
        continue;
      }
      branchName_pt_[met_option]  = getBranchName_jetMET(branchName_obj_, era_, met_option, true);
      branchName_phi_[met_option] = getBranchName_jetMET(branchName_obj_, era_, met_option, false);
    }
    branchName_sumEt_ = Form("%s_%s", branchName_obj_.data(), "sumEt");
    branchName_covXX_ = Form("%s_%s", branchName_cov_.data(), "covXX");
    branchName_covXY_ = Form("%s_%s", branchName_cov_.data(), "covXY");
    branchName_covYY_ = Form("%s_%s", branchName_cov_.data(), "covYY");
    instances_[branchName_obj_] = this;
  }
  ++numInstances_[branchName_obj_];
}

std::vector<std::string>
RecoMEtReader::setBranchAddresses(TTree * tree)
{
  if(instances_[branchName_obj_] == this)
  {
    BranchAddressInitializer bai(tree);
    bai.setBranchAddress(met_.systematics_[ptPhiOption_].pt_,  branchName_pt_[ptPhiOption_]);
    bai.setBranchAddress(met_.systematics_[ptPhiOption_].phi_, branchName_phi_[ptPhiOption_]);
    if(read_ptPhi_systematics_)
    {
      for(int met_option = kJetMET_central_nonNominal; met_option <= kJetMET_UnclusteredEnDown; ++met_option)
      {
        if(! isValidJESsource(era_, met_option))
        {
          continue;
        }
        if(met_option == ptPhiOption_)
        {
          continue; // do not bind the same branch twice
        }
        met_.systematics_[met_option] = {0., 0.};
        bai.setBranchAddress(met_.systematics_[met_option].pt_,  branchName_pt_[met_option]);
        bai.setBranchAddress(met_.systematics_[met_option].phi_, branchName_phi_[met_option]);
      }
    }
    bai.setBranchAddress(met_.sumEt_, branchName_sumEt_);
    bai.setBranchAddress(met_.covXX_, branchName_covXX_);
    bai.setBranchAddress(met_.covXY_, branchName_covXY_);
    bai.setBranchAddress(met_.covYY_, branchName_covYY_);

    return bai.getBoundBranchNames();
  }
  return {};
}

void
RecoMEtReader::set_phiModulationCorrDetails(const EventInfo * const eventInfo,
                                            const RecoVertex * const recoVertex,
                                            bool enable_phiModulationCorr)
{
  eventInfo_ = eventInfo;
  recoVertex_ = recoVertex;
  enable_phiModulationCorr_ = enable_phiModulationCorr;
}

RecoMEt
RecoMEtReader::read() const
{
  const RecoMEtReader * const gInstance = instances_[branchName_obj_];
  assert(gInstance);
  RecoMEt met = met_;
  if(enable_phiModulationCorr_)
  {
    const std::pair TheXYCorr_Met_MetPhi = METXYCorr_Met_MetPhi(eventInfo_, recoVertex_, era_);
    met.shift_PxPy(TheXYCorr_Met_MetPhi);
  }
  met.set_default(ptPhiOption_);
  return met;
}
