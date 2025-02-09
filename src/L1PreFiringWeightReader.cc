#include "tthAnalysis/HiggsToTauTau/interface/L1PreFiringWeightReader.h" // L1PreFiringWeightReader

#include "tthAnalysis/HiggsToTauTau/interface/BranchAddressInitializer.h" // BranchAddressInitializer, TTree, Form()
#include "tthAnalysis/HiggsToTauTau/interface/analysisAuxFunctions.h" // Era::k2018
#include "tthAnalysis/HiggsToTauTau/interface/sysUncertOptions.h" // L1PreFiringWeightSys

#include <cassert> // assert()

std::map<std::string, int> L1PreFiringWeightReader::numInstances_;
std::map<std::string, L1PreFiringWeightReader*> L1PreFiringWeightReader::instances_;

L1PreFiringWeightReader::L1PreFiringWeightReader(Era era)
  : era_(era)
  , branchName_l1PreFiringWeight_("L1PreFiringWeight")
  , l1PreFiringWeight_nominal_(1.)
{
  setBranchNames();
}

L1PreFiringWeightReader::~L1PreFiringWeightReader()
{
  --numInstances_[branchName_l1PreFiringWeight_];
  assert(numInstances_[branchName_l1PreFiringWeight_] >= 0);

  if(numInstances_[branchName_l1PreFiringWeight_] == 0)
  {
    L1PreFiringWeightReader * const gInstance = instances_[branchName_l1PreFiringWeight_];
    assert(gInstance);
    instances_[branchName_l1PreFiringWeight_] = nullptr;
  }
}

void
L1PreFiringWeightReader::setBranchNames()
{
  if(numInstances_[branchName_l1PreFiringWeight_] == 0)
  {
    instances_[branchName_l1PreFiringWeight_] = this;
  }
  ++numInstances_[branchName_l1PreFiringWeight_];
}

std::vector<std::string>
L1PreFiringWeightReader::setBranchAddresses(TTree * tree)
{
  if(era_ != Era::k2018)
  {
    BranchAddressInitializer bai(tree);
    bai.setBranchAddress(l1PreFiringWeight_nominal_, Form("%s_Nom", branchName_l1PreFiringWeight_.data()), 1.);
    bai.setBranchAddress(l1PreFiringWeight_up_,      Form("%s_Up",  branchName_l1PreFiringWeight_.data()), 1.);
    bai.setBranchAddress(l1PreFiringWeight_down_,    Form("%s_Dn",  branchName_l1PreFiringWeight_.data()), 1.);
    return bai.getBoundBranchNames();
  }
  return {};
}

double
L1PreFiringWeightReader::getWeight(L1PreFiringWeightSys option) const
{
  switch(option)
  {
    case L1PreFiringWeightSys::nominal: return l1PreFiringWeight_nominal_; break;
    case L1PreFiringWeightSys::up:      return l1PreFiringWeight_up_;      break;
    case L1PreFiringWeightSys::down:    return l1PreFiringWeight_down_;    break;
    default: throw cmsException(this, __func__, __LINE__) << "Invalid option";
  }
}
