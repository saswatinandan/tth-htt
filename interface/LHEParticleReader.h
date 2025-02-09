#ifndef tthAnalysis_HiggsToTauTau_LHEParticleReader_h
#define tthAnalysis_HiggsToTauTau_LHEParticleReader_h

#include "tthAnalysis/HiggsToTauTau/interface/ReaderBase.h" // ReaderBase
#include "tthAnalysis/HiggsToTauTau/interface/LHEParticle.h" // LHEParticle

#include <map> // std::map<,>

// forward declarations
class TTree;

class LHEParticleReader
  : public ReaderBase
{
public:
  LHEParticleReader();
  LHEParticleReader(const std::string & branchName_particles);
  ~LHEParticleReader();

  /**
   * @brief Call tree->SetBranchAddress for all LHEParticle branches
   */
  std::vector<std::string>
  setBranchAddresses(TTree * tree) override;

  /**
   * @brief Read branches from tree and use information to fill collection of LHEParticle objects
   * @return Collection of LHEParticle objects
   */
  std::vector<LHEParticle>
  read() const;

protected:
 /**
   * @brief Initialize names of branches to be read from tree
   */
  void
  setBranchNames();

  const unsigned int max_nParticles_;
  std::string branchName_nParticles_;
  std::string branchName_particles_;

  std::string branchName_particle_pt_;
  std::string branchName_particle_eta_;
  std::string branchName_particle_phi_;
  std::string branchName_particle_mass_;
  std::string branchName_particle_pdgId_;
  std::string branchName_particle_status_;

  UInt_t nParticles_;
  Float_t * particle_pt_;
  Float_t * particle_eta_;
  Float_t * particle_phi_;
  Float_t * particle_mass_;
  Int_t * particle_pdgId_;
  Int_t * particle_status_;

  // CV: make sure that only one LHEParticleReader instance exists for a given branchName,
  //     as ROOT cannot handle multiple TTree::SetBranchAddress calls for the same branch.
  static std::map<std::string, int> numInstances_;
  static std::map<std::string, LHEParticleReader*> instances_;
};

#endif // tthAnalysis_HiggsToTauTau_LHEParticleReader_h

