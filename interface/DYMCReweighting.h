#ifndef tthAnalysis_HiggsToTauTau_DYMCReweighting_h
#define tthAnalysis_HiggsToTauTau_DYMCReweighting_h

#include "tthAnalysis/HiggsToTauTau/interface/GenParticle.h"          // GenParticle
#include "tthAnalysis/HiggsToTauTau/interface/analysisAuxFunctions.h" // Era

#include <map> // std::map<,>

// forward declarations
class TFile;
class lutWrapperBase;

class DYMCReweighting
{
public:
  DYMCReweighting(Era era,
                  bool debug = false);
  ~DYMCReweighting();

  // reweight (LO) Drell-Yan MC produced by Madgraph in order to improve modelling of dilepton mass and pT distribution,
  // cf. slides 7-11 of presentation by Alexei Raspereza in HTT meeting on October 10th 2018
  // (https://indico.cern.ch/event/762837/contributions/3172618/attachments/1731302/2798220/Recoils_20181010.pdf )
  double
  getWeight(const std::vector<GenParticle> & genTauLeptons,
            int central_or_shift) const;

protected:
  Era era_;
  bool debug_;
  std::map<std::string, TFile *> inputFiles_;
  lutWrapperBase * weights_;
};

#endif // tthAnalysis_HiggsToTauTau_DYMCReweighting_h
