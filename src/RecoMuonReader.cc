#include "tthAnalysis/HiggsToTauTau/interface/RecoMuonReader.h"

#include "tthAnalysis/HiggsToTauTau/interface/RecoLeptonReader.h"         // RecoLeptonReader
#include "tthAnalysis/HiggsToTauTau/interface/cmsException.h"             // cmsException()
#include "tthAnalysis/HiggsToTauTau/interface/BranchAddressInitializer.h" // BranchAddressInitializer, TTree, Form()

std::map<std::string, int> RecoMuonReader::numInstances_;
std::map<std::string, RecoMuonReader *> RecoMuonReader::instances_;

RecoMuonReader::RecoMuonReader(Era era,
                               bool isMC,
                               bool readGenMatching)
  : RecoMuonReader(era, "Muon", isMC, readGenMatching)
{}

RecoMuonReader::RecoMuonReader(Era era,
                               const std::string & branchName_obj,
                               bool isMC,
                               bool readGenMatching)
  : era_(era)
  , branchName_num_(Form("n%s", branchName_obj.data()))
  , branchName_obj_(branchName_obj)
  , leptonReader_(new RecoLeptonReader(branchName_obj_, era, isMC, readGenMatching))
  , mediumIdPOG_(nullptr)
  , segmentCompatibility_(nullptr)
  , ptErr_(nullptr)
  , mvaTTH_wp_(-1.)
{
  setBranchNames();
}

RecoMuonReader::~RecoMuonReader()
{
  --numInstances_[branchName_obj_];
  assert(numInstances_[branchName_obj_] >= 0);

  if(numInstances_[branchName_obj_] == 0)
  {
    RecoMuonReader * gInstance = instances_[branchName_obj_];
    assert(gInstance);
    delete gInstance->leptonReader_;
    delete[] gInstance->mediumIdPOG_;
    delete[] gInstance->segmentCompatibility_;
    delete[] gInstance->ptErr_;
    instances_[branchName_obj_] = nullptr;
  }
}

void
RecoMuonReader::setBranchNames()
{
  if(numInstances_[branchName_obj_] == 0)
  {
    // Karl: already includes HIP mitigation for 2016 B-F
    branchName_mediumIdPOG_ = Form("%s_%s", branchName_obj_.data(), "mediumId");
    branchName_segmentCompatibility_ = Form("%s_%s", branchName_obj_.data(), "segmentComp");
    branchName_ptErr_ = Form("%s_%s", branchName_obj_.data(), "ptErr");
    instances_[branchName_obj_] = this;
  }
  else
  {
    if(branchName_num_ != instances_[branchName_obj_]->branchName_num_)
    {
      throw cmsException(this)
        << "Association between configuration parameters 'branchName_num' and 'branchName_obj' must be unique:"
        << " present association 'branchName_num' = " << branchName_num_ << " with 'branchName_obj' = " << branchName_obj_
        << " does not match previous association 'branchName_num' = " << instances_[branchName_obj_]->branchName_num_
        << " with 'branchName_obj' = " << instances_[branchName_obj_]->branchName_obj_ << " !!\n";
    }
  }
  ++numInstances_[branchName_obj_];
}

std::vector<std::string>
RecoMuonReader::setBranchAddresses(TTree * tree)
{
  std::vector<std::string> bound_branches;
  if(instances_[branchName_obj_] == this)
  {
    const std::vector<std::string> leptonBranches = leptonReader_->setBranchAddresses(tree);
    bound_branches.insert(bound_branches.end(), leptonBranches.begin(), leptonBranches.end());

    const unsigned int max_nLeptons = leptonReader_->max_nLeptons_;
    BranchAddressInitializer bai(tree, max_nLeptons);
    bai.setBranchAddress(mediumIdPOG_, branchName_mediumIdPOG_);
    bai.setBranchAddress(segmentCompatibility_, branchName_segmentCompatibility_);
    bai.setBranchAddress(ptErr_, branchName_ptErr_, -1.);

    const std::vector<std::string> recoMuonBranches = bai.getBoundBranchNames();
    bound_branches.insert(bound_branches.end(), recoMuonBranches.begin(), recoMuonBranches.end());
  }
  return bound_branches;
}

std::vector<RecoMuon>
RecoMuonReader::read() const
{
  const RecoLeptonReader * const gLeptonReader = leptonReader_->instances_[branchName_obj_];
  assert(gLeptonReader);
  const RecoMuonReader * const gMuonReader = instances_[branchName_obj_];
  assert(gMuonReader);
  std::vector<RecoMuon> muons;
  const UInt_t nLeptons = gLeptonReader->nLeptons_;
  
  if(nLeptons > leptonReader_->max_nLeptons_)
  {
    throw cmsException(this)
      << "Number of leptons stored in Ntuple = " << nLeptons << ", exceeds max_nLeptons = "
      << leptonReader_->max_nLeptons_;
  }
  if(nLeptons > 0)
  {
    muons.reserve(nLeptons);
    for(UInt_t idxLepton = 0; idxLepton < nLeptons; ++idxLepton)
    {
      if(std::abs(gLeptonReader->pdgId_[idxLepton]) == 13)
      {
        // Karl: For *some* leptons that don't have an associated jet,
        //       the deepCSV score is nan (and not -1) for an unknown reason.
        //       Adding a safeguard for these instances.
        muons.push_back(RecoMuon({
          {
            {
              gLeptonReader->pt_[idxLepton],
              gLeptonReader->eta_[idxLepton],
              gLeptonReader->phi_[idxLepton],
              gLeptonReader->mass_[idxLepton],
              gLeptonReader->pdgId_[idxLepton],
              gLeptonReader->charge_[idxLepton],
            },
            gLeptonReader->dxy_[idxLepton],
            gLeptonReader->dz_[idxLepton],
            gLeptonReader->relIso_all_[idxLepton],
            gLeptonReader->pfRelIso04_all_[idxLepton],
            gLeptonReader->relIso_chg_[idxLepton],
            gLeptonReader->relIso_neu_[idxLepton],
            gLeptonReader->sip3d_[idxLepton],
            gLeptonReader->mvaRawTTH_[idxLepton],
            gLeptonReader->jetPtRatio_[idxLepton],
            gLeptonReader->jetPtRel_[idxLepton],
            gLeptonReader->jetNDauChargedMVASel_[idxLepton],
            gLeptonReader->tightCharge_[idxLepton],
            gLeptonReader->filterBits_[idxLepton],
            gLeptonReader->jetIdx_[idxLepton],
            gLeptonReader->genPartFlav_[idxLepton],
            gLeptonReader->genMatchIdx_[idxLepton],
          },
          true, // Karl: all muon objects pass Muon POG's loose definition at the nanoAOD prodction level
          gMuonReader->mediumIdPOG_[idxLepton],
          gMuonReader->segmentCompatibility_[idxLepton],
          gMuonReader->ptErr_[idxLepton]
        }));

        RecoMuon & muon = muons.back();
        for(const auto & kv: gLeptonReader->jetBtagCSVs_)
        {
          const double val = kv.second[idxLepton];
          muon.jetBtagCSVs_[kv.first] = std::isnan(val) ? -2. : val;
        }
        for(const auto & kv: gLeptonReader->assocJetBtagCSVs_)
        {
          const double val = kv.second[idxLepton];
          muon.assocJetBtagCSVs_[kv.first] = std::isnan(val) ? -2. : val;
        }
        if(mvaTTH_wp_ > 0.)
        {
          muon.set_mvaRawTTH_cut(mvaTTH_wp_);
        }
      }
    }
    gLeptonReader->readGenMatching(muons);
  }
  return muons;
}

void
RecoMuonReader::set_mvaTTH_wp(double mvaTTH_wp)
{
  mvaTTH_wp_ = mvaTTH_wp;
}
