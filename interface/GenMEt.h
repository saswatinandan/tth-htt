#ifndef tthAnalysis_HiggsToTauTau_GenMEt_h
#define tthAnalysis_HiggsToTauTau_GenMEt_h

#include "tthAnalysis/HiggsToTauTau/interface/GenParticle.h" // Particle::LorentzVector

#include <TMatrixD.h> // TMatrixD

#include <map> // std::map<,>

class GenMEt
{
public:
  GenMEt();

  GenMEt(Float_t pt,
	 Float_t phi);

  GenMEt(const math::PtEtaPhiMLorentzVector & p4);

  GenMEt &
  operator=(const GenMEt & other);

  /**
   * @brief Funtions to access data-members
   * @return Values of data-members
   */
  Double_t pt() const;
  Double_t phi() const;
  
  const Particle::LorentzVector & p4() const;

  friend class GenMEtReader;
  // friend class GenMEtWriter; // since not needed 

protected:

  Float_t pt_;  ///< pT of missing transverse momentum vector
  Float_t phi_; ///< phi of missing transverse momentum vector

  ///< (default) 4-momentum constructed from pT and phi, assuming eta and mass to be equal to zero
  math::PtEtaPhiMLorentzVector p4_;

  ///< Update cov and p4 (needed by GenMEtReader)
  void update();
  void update_p4();
};

std::ostream& operator<<(std::ostream& stream, const GenMEt& met);

#endif // tthAnalysis_HiggsToTauTau_GenMEt_h  
