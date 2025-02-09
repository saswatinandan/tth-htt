#include "tthAnalysis/HiggsToTauTau/interface/RecoElectron.h"

#include "tthAnalysis/HiggsToTauTau/interface/GenLepton.h" // GenLepton
#include "tthAnalysis/HiggsToTauTau/interface/analysisAuxFunctions.h" // EGammaID, EGammaWP, as_integer(), cmsException()

RecoElectron::RecoElectron(const RecoLepton & lepton,
                           Double_t eCorr,
                           Double_t sigmaEtaEta,
                           Double_t HoE,
                           Double_t deltaEta,
                           Double_t deltaPhi,
                           Double_t OoEminusOoP,
                           Int_t    nLostHits,
                           Bool_t   passesConversionVeto,
                           Int_t    cutbasedID_HLT)
  : RecoLepton(lepton)
  , eCorr_(eCorr)
  , sigmaEtaEta_(sigmaEtaEta)
  , HoE_(HoE)
  , deltaEta_(deltaEta) 
  , deltaPhi_(deltaPhi) 
  , OoEminusOoP_(OoEminusOoP)
  , nLostHits_(nLostHits)
  , passesConversionVeto_(passesConversionVeto)
  , cutbasedID_HLT_(cutbasedID_HLT)
{
  set_mvaRawTTH_cut(0.80);
}

Double_t
RecoElectron::eCorr() const
{
  return eCorr_;
}

Double_t
RecoElectron::mvaRaw_POG() const
{
  return mvaRaw_POG(EGammaID::Fall17V2noIso);
}

Double_t
RecoElectron::mvaRaw_POG(EGammaID id) const
{
  if(! egammaID_raws_.count(id))
  {
    throw cmsException(this, __func__, __LINE__) << "No such EGamma ID available: " << as_integer(id);
  }
  return egammaID_raws_.at(id);
}

Bool_t
RecoElectron::mvaID_POG(EGammaWP wp) const
{
  return mvaID_POG(EGammaID::Fall17V2noIso, wp);
}

Bool_t
RecoElectron::mvaID_POG(EGammaID id,
                        EGammaWP wp) const
{
  if(! egammaID_ids_.count(id))
  {
    throw cmsException(this, __func__, __LINE__) << "No such EGamma ID available: " << as_integer(id);
  }
  const auto & egammaID_id = egammaID_ids_.at(id);
  if(! egammaID_id.count(wp))
  {
    throw cmsException(this, __func__, __LINE__)
      << "No such EGamma WP available for EGamma ID " << as_integer(id) << ": " << as_integer(wp)
    ;
  }
  return egammaID_id.at(wp);
}

Double_t
RecoElectron::sigmaEtaEta() const
{
  return sigmaEtaEta_;
}

Double_t
RecoElectron::HoE() const
{
  return HoE_;
}

Double_t
RecoElectron::deltaEta() const
{
  return deltaEta_; 
}

Double_t
RecoElectron::deltaPhi() const
{
  return deltaPhi_; 
}

Double_t
RecoElectron::etaSC() const
{
  return deltaEta_ + eta_;
}

Double_t
RecoElectron::phiSC() const
{
  return deltaPhi_ + phi_;
}

Double_t
RecoElectron::absEtaSC() const
{
  return std::fabs(etaSC());
}

Double_t
RecoElectron::OoEminusOoP() const
{
  return OoEminusOoP_;
}

Int_t
RecoElectron::nLostHits() const
{
  return nLostHits_;
}

Bool_t
RecoElectron::passesConversionVeto() const
{
  return passesConversionVeto_;
}

Int_t
RecoElectron::cutbasedID_HLT() const
{
  return cutbasedID_HLT_;
}

bool
RecoElectron::is_electron() const
{
  return true;
}

bool
RecoElectron::is_muon() const
{
  return false;
}

bool
RecoElectron::isGenMatched(bool requireChargeMatch) const
{
  return requireChargeMatch                                    ?
    (!! genLepton_ ? charge() == genLepton_->charge() : false) :
    (!! genLepton_ || !! genPhoton_)
  ;
}

std::ostream &
operator<<(std::ostream & stream,
           const RecoElectron & electron)
{
  stream << static_cast<const RecoLepton & >(electron)                      << ",\n "
            "eCorr = "                << electron.eCorr()                   << ", "
            "mvaRaw_POG = "           << electron.mvaRaw_POG()              << " ("
            "mvaID_POG (loose) = "    << electron.mvaID_POG(EGammaWP::WPL)  << ", "
            "mvaID_POG (80%) = "      << electron.mvaID_POG(EGammaWP::WP80) << "), "
            "nLostHits = "            << electron.nLostHits()               << ",\n "
            "passesConversionVeto = " << electron.passesConversionVeto()    << ", "
            "sigmaEtaEta = "          << electron.sigmaEtaEta()             << ", "
            "deltaEta = "             << electron.deltaEta()                << ",\n "
            "deltaPhi = "             << electron.deltaPhi()                << ", "
            "HoE = "                  << electron.HoE()                     << ", "
            "OoEminusOoP = "          << electron.OoEminusOoP()             << ",\n"
            "cutbasedID_HLT = "       << electron.cutbasedID_HLT()          << '\n'
  ;
  return stream;
}
