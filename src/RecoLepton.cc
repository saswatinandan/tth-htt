#include "tthAnalysis/HiggsToTauTau/interface/RecoLepton.h" // RecoLepton, ChargedParticle
#include "tthAnalysis/HiggsToTauTau/interface/GenLepton.h" // GenLepton
#include "tthAnalysis/HiggsToTauTau/interface/GenHadTau.h" // GenHadTau
#include "tthAnalysis/HiggsToTauTau/interface/GenPhoton.h" // GenPhoton
#include "tthAnalysis/HiggsToTauTau/interface/GenJet.h" // GenJet
#include "tthAnalysis/HiggsToTauTau/interface/analysisAuxFunctions.h" // cmsException(), as_integer()

Double_t
RecoLepton::get_assocJet_pt(Double_t reco_pt,
                            Double_t jetPtRatio)
{
  return jetPtRatio > 1.e-3 ? 0.90 * reco_pt / jetPtRatio : reco_pt;
}

RecoLepton::RecoLepton(const ChargedParticle & lepton,
                       Double_t dxy,
                       Double_t dz,
                       Double_t relIso,
                       Double_t pfRelIso04All,
                       Double_t miniRelIsoCharged,
                       Double_t miniRelIsoNeutral,
                       Double_t sip3d,
                       Double_t mvaRawTTH,
                       Double_t jetPtRatio,
                       Double_t jetPtRel,
                       Int_t    jetNDauChargedMVASel,
                       Int_t    tightCharge,
                       UInt_t   filterBits,
                       Int_t    jetIdx,
                       UChar_t  genPartFlav,
                       Int_t    genMatchIdx)
  : ChargedParticle(lepton)
  , dxy_(dxy)
  , dz_(dz)
  , relIso_(relIso)
  , pfRelIso04All_(pfRelIso04All)
  , miniRelIsoCharged_(miniRelIsoCharged)
  , miniRelIsoNeutral_(miniRelIsoNeutral)
  , sip3d_(sip3d)
  , mvaRawTTH_(mvaRawTTH)
  , jetPtRatio_(jetPtRatio)
  , jetPtRel_(jetPtRel)
  , jetNDauChargedMVASel_(jetNDauChargedMVASel)
  , tightCharge_(tightCharge)
  , filterBits_(filterBits)
  , jetIdx_(jetIdx)
  , genPartFlav_(genPartFlav)
  , genMatchIdx_(genMatchIdx)
  , mvaRawTTH_cut_(-1.)
  , assocJet_pt_(get_assocJet_pt(pt_, jetPtRatio_))
  , assocJet_p4_(assocJet_pt_, eta_, phi_, mass_)
  , genLepton_(nullptr)
  , genHadTau_(nullptr)
  , genJet_(nullptr)
  , isCMSPOG_(false)
  , isLoose_(false)
  , isFakeable_(false)
  , isTight_(false)
{}

RecoLepton::~RecoLepton()
{}

void
RecoLepton::set_isCMSPOG() const
{
  isCMSPOG_ = true;
}

void
RecoLepton::set_isLoose() const
{
  isLoose_ = true;
}

void
RecoLepton::set_isFakeable() const
{
  isFakeable_ = true;
}

void
RecoLepton::set_isTight() const
{
  isTight_ = true;
}

void
RecoLepton::set_genLepton(const GenLepton * genLepton)
{
  genLepton_.reset(genLepton);
}

void
RecoLepton::set_genHadTau(const GenHadTau * genHadTau)
{
  genHadTau_.reset(genHadTau);
}

void
RecoLepton::set_genPhoton(const GenPhoton * genPhoton)
{
  genPhoton_.reset(genPhoton);
}

void
RecoLepton::set_genJet(const GenJet * genJet)
{
  genJet_.reset(genJet);
}

bool
RecoLepton::is_electron() const
{
  return false;
}

bool
RecoLepton::is_muon() const
{
  return false;
}

Double_t
RecoLepton::lepton_pt() const
{
  return pt_;
}

const Particle::LorentzVector &
RecoLepton::lepton_p4() const
{
  return p4_;
}

Double_t
RecoLepton::cone_pt() const
{
  assert(mvaRawTTH_cut_ > 0.);
  return mvaRawTTH_ >= mvaRawTTH_cut_ ? pt_ : assocJet_pt_;
}

const Particle::LorentzVector &
RecoLepton::cone_p4() const
{
  assert(mvaRawTTH_cut_ > 0.);
  return mvaRawTTH_ >= mvaRawTTH_cut_ ? p4_ : assocJet_p4_;
}

Double_t
RecoLepton::assocJet_pt() const
{
  return assocJet_pt_;
}

const Particle::LorentzVector &
RecoLepton::assocJet_p4() const
{
  return assocJet_p4_;
}

Double_t
RecoLepton::dxy() const
{
  return dxy_;
}

Double_t
RecoLepton::dz() const
{
  return dz_;
}

Double_t
RecoLepton::relIso() const
{
  return relIso_;
}

Double_t
RecoLepton::pfRelIso04All() const
{
  return pfRelIso04All_;
}

Double_t
RecoLepton::miniRelIsoCharged() const
{
  return miniRelIsoCharged_;
}

Double_t
RecoLepton::miniRelIsoNeutral() const
{
  return miniRelIsoNeutral_;
}

Double_t
RecoLepton::sip3d() const
{
  return sip3d_;
}

Double_t
RecoLepton::mvaRawTTH() const
{
  return mvaRawTTH_;
}

Double_t
RecoLepton::mvaRawTTH_cut() const
{
  assert(mvaRawTTH_cut_ > 0.);
  return mvaRawTTH_cut_;
}

Double_t
RecoLepton::jetPtRatio() const
{
  return jetPtRatio_;
}

Double_t
RecoLepton::jetRelIso() const
{
  return 1. / jetPtRatio_ - 1.;
}

Double_t
RecoLepton::jetPtRel() const
{
  return jetPtRel_;
}

Double_t
RecoLepton::jetBtagCSV(bool doAssoc) const
{
  return jetBtagCSV(Btag::kDeepJet, doAssoc);
}

Double_t
RecoLepton::jetBtagCSV(Btag btag,
                       bool doAssoc) const
{
  if(! hasJetBtagCSV(btag, doAssoc))
  {
    throw cmsException(this, __func__, __LINE__)
      << "b-tagging discriminator not available: " << as_integer(btag)
    ;
  }
  return doAssoc ? assocJetBtagCSVs_.at(btag) : jetBtagCSVs_.at(btag);
}

Int_t
RecoLepton::jetNDauChargedMVASel() const
{
  return jetNDauChargedMVASel_;
}

Int_t
RecoLepton::tightCharge() const
{
  return tightCharge_;
}

UInt_t
RecoLepton::filterBits() const
{
  return filterBits_;
}

Int_t
RecoLepton::jetIdx() const
{
  return jetIdx_;
}

UChar_t
RecoLepton::genPartFlav() const
{
  return genPartFlav_;
}

Int_t
RecoLepton::genMatchIdx() const
{
  return genMatchIdx_;
}

const GenLepton *
RecoLepton::genLepton() const
{
  return genLepton_.get();
}

const GenHadTau *
RecoLepton::genHadTau() const
{
  return genHadTau_.get();
}

const GenPhoton *
RecoLepton::genPhoton() const
{
  return genPhoton_.get();
}

const GenJet *
RecoLepton::genJet() const
{
  return genJet_.get();
}

bool
RecoLepton::hasJetBtagCSV(Btag btag,
                          bool doAssoc) const
{
  return doAssoc ? assocJetBtagCSVs_.count(btag) : jetBtagCSVs_.count(btag);
}

bool
RecoLepton::isGenMatched(bool requireChargeMatch) const
{
  return requireChargeMatch ? (!! genLepton_ ? charge() == genLepton_->charge() : false) : !! genLepton_;
}

bool
RecoLepton::hasAnyGenMatch() const
{
  return !! genLepton_ || !! genHadTau_ || !! genPhoton_ || !! genJet_;
}

bool
RecoLepton::isCMSPOG() const
{
  return isCMSPOG_;
}

bool
RecoLepton::isLoose() const
{
  return isLoose_;
}

bool
RecoLepton::isFakeable() const
{
  return isFakeable_;
}

bool
RecoLepton::isTight() const
{
  return isTight_;
}

void
RecoLepton::set_mvaRawTTH_cut(Double_t mvaRawTTH_cut)
{
  mvaRawTTH_cut_ = mvaRawTTH_cut;
}

void
RecoLepton::set_p4(const Particle::LorentzVector & p4)
{
  Particle::set_p4(p4);
  set_assocJet_p4();
}

void
RecoLepton::set_ptEtaPhiMass(Double_t pt,
                             Double_t eta,
                             Double_t phi,
                             Double_t mass)
{
  Particle::set_ptEtaPhiMass(pt, eta, phi, mass);
  set_assocJet_p4();
}

void
RecoLepton::set_assocJet_p4()
{
  assocJet_pt_ = get_assocJet_pt(pt_, jetPtRatio_);
  assocJet_p4_ = { assocJet_pt_, eta_, phi_, mass_ };
}

std::ostream &
operator<<(std::ostream & stream,
           const RecoLepton & lepton)
{
  stream << static_cast<const ChargedParticle &>(lepton)  << ",\n"
            " cone_pT = "             << lepton.cone_pt()                       << ","
            " assocJet pT = "         << lepton.assocJet_pt()                   << ","
            " dxy = "                 << lepton.dxy()                           << ","
            " dz = "                  << lepton.dz()                            << ",\n"
            " sip3d = "               << lepton.sip3d()                         << ","
            " relIso = "              << lepton.relIso()                        << ","
            " pfRelIso04All = "       << lepton.pfRelIso04All()                 << ",\n"
            " genPartFlav = "         << static_cast<int>(lepton.genPartFlav()) << ", "
            " tightCharge = "         << lepton.tightCharge()                   << ","
            " jetPtRatio = "          << lepton.jetPtRatio()                    << ","
            " jetPtRel = "            << lepton.jetPtRel()                      << ",\n"
            " jetBtagCSV(default) = " << lepton.jetBtagCSV()
  ;
  std::cout << ", jetBtagCSV(DeepJet) = ";
  if(lepton.hasJetBtagCSV(Btag::kDeepJet)) stream << lepton.jetBtagCSV(Btag::kDeepJet); else stream << "N/A";
  std::cout << ", jetBtagCSV(DeepCSV) = ";
  if(lepton.hasJetBtagCSV(Btag::kDeepCSV)) stream << lepton.jetBtagCSV(Btag::kDeepCSV); else stream << "N/A";
  std::cout << ", jetBtagCSV(CSVv2) = ";
  if(lepton.hasJetBtagCSV(Btag::kCSVv2))   stream << lepton.jetBtagCSV(Btag::kCSVv2);   else stream << "N/A";
  stream << ", mvaRawTTH = "               << lepton.mvaRawTTH()
         << ", mvaRawTTH_cut = "           << lepton.mvaRawTTH_cut()
         << ", passes CMS POG ID = "       << lepton.isCMSPOG()
         << ", is loose/fakeable/tight = " << lepton.isLoose()       << '/'
                                           << lepton.isFakeable()    << '/'
                                           << lepton.isTight()       << ",\n"
            " gen. matching:";
  stream << ",\n  lepton = " << lepton.genLepton();
  if(lepton.genLepton())
  {
    stream << ": " << *(lepton.genLepton());
  }
  stream << ",\n  hadTau = " << lepton.genHadTau();
  if(lepton.genHadTau())
  {
    stream << ": " << *(lepton.genHadTau());
  }
  stream << ",\n  photon = " << lepton.genPhoton();
  if(lepton.genPhoton())
  {
    stream << ": " << *(lepton.genPhoton());
  }
  stream << ",\n  jet    = " << lepton.genJet();
  if(lepton.genJet())
  {
    stream << ": " << *(lepton.genJet());
  }
  stream << '\n';
  return stream;
}
