#include "tthAnalysis/HiggsToTauTau/interface/hadTopTaggerAuxFunctions.h"

#include "tthAnalysis/HiggsToTauTau/interface/analysisAuxFunctions.h" // isHigherPt()
#include "tthAnalysis/HiggsToTauTau/interface/cmsException.h" // cmsException()
#include "tthAnalysis/HiggsToTauTau/interface/hadTopTaggerAuxFunctions_geral.h" // kGenTop...

#include <algorithm> // std::sort()
#include <numeric> // iota
#include <map>

std::map<int, Particle::LorentzVector>
isGenMatchedJetTripletVar(const std::vector<GenParticle> & genTopQuarks,
                          const std::vector<GenParticle> & genBJets,
                          const std::vector<GenParticle> & genWBosons,
                          const std::vector<GenParticle> & genQuarkFromTop,
                          int mode)
{
  std::map<int, Particle::LorentzVector> genMatchValues = {
    { kGenTop,        Particle::LorentzVector(0., 0., 0., 0.) },
    { kGenTopB,       Particle::LorentzVector(0., 0., 0., 0.) },
    { kGenTopW,       Particle::LorentzVector(0., 0., 0., 0.) },
    { kGenTopWj1,     Particle::LorentzVector(0., 0., 0., 0.) },
    { kGenTopWj2,     Particle::LorentzVector(0., 0., 0., 0.) }
  };

  int pddIdTop, pddIdWBoson, pdgIdBJet;
  if(mode == kGenTop)
  {
    pddIdTop    =  +6;
    pddIdWBoson = +24;
    pdgIdBJet   =  +5;
  }
  else if(mode == kGenAntiTop)
  {
    pddIdTop    =  -6;
    pddIdWBoson = -24;
    pdgIdBJet   =  -5;
  }
  else
  {
    throw cmsException(__func__, __LINE__) << "Invalid parameter mode = " << mode;
  }

  const GenParticle * genTop = 0;
  for(const GenParticle & genTopQuark: genTopQuarks)
  {
    if(genTopQuark.pdgId() == pddIdTop && ! genTop)
    {
      genTop = &genTopQuark;
      genMatchValues[kGenTop].SetPxPyPzE(genTopQuark.p4().x(), genTopQuark.p4().y(), genTopQuark.p4().z(), genTopQuark.p4().energy());
    }
  }

  const GenParticle * genBJetFromTop = nullptr;
  for(const GenParticle & genBJet: genBJets)
  {
    if(genBJet.pdgId() == pdgIdBJet && ! genBJetFromTop )
    {
      genBJetFromTop = &genBJet;
      genMatchValues[kGenTopB].SetPxPyPzE(genBJet.p4().x(), genBJet.p4().y(), genBJet.p4().z(), genBJet.p4().energy());
    }
  }

  std::vector<const GenParticle *> genWBosonsFromTop; // = nullptr;
  for(const GenParticle & genWBoson: genWBosons)
  {
    if(genWBoson.pdgId() == pddIdWBoson ) // && ! genWBosonFromTop
    {
      genWBosonsFromTop.push_back(&genWBoson);
    }
  }

  double mass_diff = 1000000.;
  const GenParticle * genWBosonFromTopFinal = nullptr;
  //for(auto itW = genWBosons.cbegin(); itW != genWBosons.cend(); ++itW)
  for(auto itW = genWBosonsFromTop.cbegin(); itW != genWBosonsFromTop.cend(); ++itW) // Edit Siddhesh
  {
    const GenParticle * genWBosonFromTop = (*itW);
    if(! genBJetFromTop || ! genTop)
    {
      continue;
    }
    if ( std::fabs((genWBosonFromTop->p4() + genBJetFromTop->p4()).mass() - genTop->p4().mass()) < mass_diff)
    {
      genWBosonFromTopFinal = genWBosonFromTop;
      genMatchValues[kGenTopW].SetPxPyPzE(genWBosonFromTop->p4().x(), genWBosonFromTop->p4().y(), genWBosonFromTop->p4().z(), genWBosonFromTop->p4().energy());
      mass_diff = std::fabs((genWBosonFromTop->p4() + genBJetFromTop->p4()).mass() - genTop->p4().mass());
    }
  }

  //std::cout<<" genWBosonsFromTop.size() "<< genWBosonsFromTop.size() << std::endl;
  std::vector<const GenParticle *> genWJetsFromTop;
  if (genWBosonFromTopFinal) {
  double genWJetsFromTop_mass = -1.;
  double genWJetsFromTop_massW = -1.;
  for(auto it1 = genQuarkFromTop.cbegin(); it1 != genQuarkFromTop.cend(); ++it1)
  {
    const GenParticle * genWJet1 = &(*it1);
    for(auto it2 = it1 + 1; it2 != genQuarkFromTop.cend(); ++it2)
    {
      const GenParticle * genWJet2 = &(*it2);
      if ( boost::math::sign(genWJet1->charge()) == boost::math::sign(genWBosonFromTopFinal->charge()) &&
	     boost::math::sign(genWJet2->charge()) == boost::math::sign(genWBosonFromTopFinal->charge()) ) // Edit Siddhesh
      {
        const double genWJetsFromTop_massCurrent = (genWJet1->p4() + genWJet2->p4() + genBJetFromTop->p4()).mass();
        const double diff_massCurrent = std::fabs((genWJet1->p4() + genWJet2->p4() + genBJetFromTop->p4()).mass() - genTop->p4().mass());
        const double diff_mass        = std::fabs(genWJetsFromTop_mass        - genTop->p4().mass());

        const double genWJetsFromTop_massWCurrent = (genWJet1->p4() + genWJet2->p4()).mass();
        const double diff_massWCurrent = std::fabs(genWJetsFromTop_massWCurrent - genWBosonFromTopFinal->p4().mass());
        const double diff_massW        = std::fabs(genWJetsFromTop_massW        - genWBosonFromTopFinal->p4().mass());

        if(genWJetsFromTop_mass == -1. || ((diff_massCurrent < diff_mass) || (diff_massWCurrent < diff_massW)))
        {
          const bool passWbosonMassVeto_top = passWbosonMassVeto(
            genWJet1, genWJet2, genWBosonFromTopFinal
          );
          if(passWbosonMassVeto_top)
          {
            genWJetsFromTop = { genWJet1, genWJet2 };
            genWJetsFromTop_mass = genWJetsFromTop_massCurrent;
            genWJetsFromTop_massW = genWJetsFromTop_massWCurrent;
          }
        }
      } // close if charge
    }
  }
  }

  if(genWJetsFromTop.size() != 2)
  {
    return genMatchValues;
  }

  std::sort(genWJetsFromTop.begin(), genWJetsFromTop.end(), isHigherPt);
  const GenParticle * genWJetFromTop_lead = genWJetsFromTop[0];
  const GenParticle * genWJetFromTop_sublead = genWJetsFromTop[1];

  genMatchValues[kGenTopWj1].SetPxPyPzE(genWJetFromTop_lead->p4().x(), genWJetFromTop_lead->p4().y(), genWJetFromTop_lead->p4().z(), genWJetFromTop_lead->p4().energy());
  genMatchValues[kGenTopWj2].SetPxPyPzE(genWJetFromTop_sublead->p4().x(), genWJetFromTop_sublead->p4().y(), genWJetFromTop_sublead->p4().z(), genWJetFromTop_sublead->p4().energy());


  return genMatchValues;

}

bool
passWbosonMassVeto(const GenParticle * genWJetFromTop_lead,
                   const GenParticle * genWJetFromTop_sublead,
                   const GenParticle * genWBosonFromTop)
{
  const double genWJetFromTop_mass = (genWJetFromTop_lead->p4() + genWJetFromTop_sublead->p4()).mass();
  return std::fabs(genWJetFromTop_mass - genWBosonFromTop->mass()) < 15. && genWJetFromTop_mass > 60.;
}

int
getType(size_t sizeHTTv2, size_t sizeFatW, size_t sizeResolved){
  int typeTop = -1;
  if (sizeHTTv2 > 0) typeTop = 1;
  else if (sizeFatW >0) typeTop = 2;
  else if (sizeResolved >0) typeTop = 3;
  return typeTop;
}

std::vector<double>
getBdiscr(std::vector<const RecoJet*> selJetsIt){
  std::vector<double> btag_disc;
  for ( std::vector<const RecoJet*>::const_iterator jetIterB = selJetsIt.begin();
  jetIterB != selJetsIt.end(); ++jetIterB ) {
    btag_disc.push_back((*jetIterB)->BtagCSV());
  }
  return btag_disc;
}

std::vector<size_t>
calRank( std::vector<double> & btag_disc ) {
    std::vector<size_t> result(btag_disc.size(),0);
    //sorted index
    std::vector<size_t> indx(btag_disc.size());
    iota(indx.begin(),indx.end(),0);
    sort(indx.begin(),indx.end(),[&btag_disc](int i1, int i2){return btag_disc[i1]>btag_disc[i2];});
    //return ranking
    for(size_t iter=0;iter<btag_disc.size();++iter) result[indx[iter]]=iter+1;
    //std::cout<<"btag discriminant  ";
    //for (auto i: btag_disc) std::cout << i << " ";
    //std::cout<<std::endl;
    return result;
}

void get_two_maximal(
  double &  MVA_1_mvaOutput,
  bool &    MVA_1_truth,
  double &  MVA_1_genTopPt,
  double &  MVA_1_recTopPt,
  double & MVA_2_mvaOutput,
  bool &   MVA_2_truth,
  double & MVA_2_genTopPt,
  double & MVA_2_recTopPt,
  double & MVA_med_mvaOutput,
  bool   & MVA_med_truth,
  double & MVA_med_genTopPt,
  double & MVA_med_recTopPt,
  double & Wj1_pt_1, double & Wj2_pt_1, double & b_pt_1,
  double & Wj1_pt_2, double & Wj2_pt_2, double & b_pt_2,
  double & Wj1_pt_med, double & Wj2_pt_med, double & b_pt_med,
  double MVA, double recTopPt, double genTopPt_teste, bool isGenMatched,
  double Wj1_pt, double Wj2_pt, double b_pt
  // remove combinations in which jets overlap among first and second

) {

  if ( MVA == 0. ) return;

  if ( MVA_1_mvaOutput == 0. || MVA_1_mvaOutput == -1. ) {
    MVA_1_mvaOutput = MVA;
    MVA_1_recTopPt = recTopPt;
    MVA_1_genTopPt = genTopPt_teste;
    Wj1_pt_1 = Wj1_pt;
    Wj2_pt_1 = Wj2_pt;
    b_pt_1 = b_pt;
    //std::cout << "start 1" << std::endl;
    return;
  }

  if ( MVA_2_mvaOutput == 0. || MVA_2_mvaOutput == -1.) {
    if (  MVA < MVA_1_mvaOutput &&
        (
        Wj1_pt == -1.0 || (
        Wj1_pt != Wj1_pt_1 && Wj1_pt != Wj2_pt_1 && Wj1_pt != b_pt_1 &&
        Wj2_pt != Wj1_pt_1 && Wj2_pt != Wj2_pt_1 && Wj2_pt != b_pt_1 &&
        b_pt   != Wj1_pt_1 && b_pt   != Wj2_pt_1 && b_pt   != b_pt_1
        )
        )
    ) {
      MVA_2_mvaOutput = MVA;
      MVA_2_recTopPt = recTopPt;
      MVA_2_genTopPt = genTopPt_teste;
      Wj1_pt_2 = Wj1_pt;
      Wj2_pt_2 = Wj2_pt;
      b_pt_2 = b_pt;
      //std::cout << "start 2" << std::endl;
      return;
    } else if ( MVA > MVA_1_mvaOutput ) {
        //std::cout << "start 1A" << std::endl;
        if (
        Wj1_pt == -1.0 || (
        Wj1_pt != Wj1_pt_1 && Wj1_pt != Wj2_pt_1 && Wj1_pt != b_pt_1 &&
        Wj2_pt != Wj1_pt_1 && Wj2_pt != Wj2_pt_1 && Wj2_pt != b_pt_1 &&
        b_pt   != Wj1_pt_1 && b_pt   != Wj2_pt_1 && b_pt   != b_pt_1
        )
        ) {
          MVA_2_mvaOutput = MVA_1_mvaOutput;
          MVA_2_truth     = MVA_1_truth;
          MVA_2_genTopPt  = MVA_1_genTopPt;
          MVA_2_recTopPt  = MVA_1_recTopPt;
          Wj1_pt_2 = Wj1_pt_1;
          Wj2_pt_2 = Wj2_pt_1;
          b_pt_2 = b_pt_1;

          MVA_1_mvaOutput = MVA;
          MVA_1_recTopPt = recTopPt;
          MVA_1_genTopPt = genTopPt_teste;
          Wj1_pt_1 = Wj1_pt;
          Wj2_pt_1 = Wj2_pt;
          b_pt_1 = b_pt;
          //std::cout << "start 2A" << std::endl;
          return;
        } else {
          MVA_1_mvaOutput = MVA;
          MVA_1_recTopPt = recTopPt;
          MVA_1_genTopPt = genTopPt_teste;
          Wj1_pt_1 = Wj1_pt;
          Wj2_pt_1 = Wj2_pt;
          b_pt_1 = b_pt;
        }

      }
      //std::cout << "could not fill 2 yet -- > try again" << std::endl;
      return;
  } //  close if MV2 not filled

  /// -- first stage
  if ( ( MVA >  MVA_2_mvaOutput && MVA < MVA_1_mvaOutput ) &&
      (
      Wj1_pt == -1.0 || (
      Wj1_pt != Wj1_pt_1 && Wj1_pt != Wj2_pt_1 && Wj1_pt != b_pt_1 &&
      Wj2_pt != Wj1_pt_1 && Wj2_pt != Wj2_pt_1 && Wj2_pt != b_pt_1 &&
      b_pt   != Wj1_pt_1 && b_pt   != Wj2_pt_1 && b_pt   != b_pt_1
      )
      )
    ) {
      MVA_med_mvaOutput = MVA_2_mvaOutput;
      MVA_med_truth     = MVA_2_truth;
      MVA_med_genTopPt  = MVA_2_genTopPt;
      MVA_med_recTopPt  = MVA_2_recTopPt;
      Wj1_pt_med = Wj1_pt_2;
      Wj2_pt_med = Wj2_pt_2;
      b_pt_med = b_pt_2;
      //
      MVA_2_truth = isGenMatched;
      MVA_2_mvaOutput = MVA;
      MVA_2_recTopPt = recTopPt;
      MVA_2_genTopPt = genTopPt_teste;
      Wj1_pt_2 = Wj1_pt;
      Wj2_pt_2 = Wj2_pt;
      b_pt_2 = b_pt;
      //std::cout << "change A (MVA_2_mvaOutput) "<< MVA_2_mvaOutput << std::endl;
  } else if ( MVA >  MVA_2_mvaOutput &&  MVA > MVA_1_mvaOutput && MVA > MVA_med_mvaOutput ) {
      if ( ( MVA_2_mvaOutput < MVA_1_mvaOutput && MVA_med_mvaOutput <  MVA_1_mvaOutput )
      &&
        (Wj1_pt == -1.0 || (
        Wj1_pt != Wj1_pt_1 && Wj1_pt != Wj2_pt_1 && Wj1_pt != b_pt_1 &&
        Wj2_pt != Wj1_pt_1 && Wj2_pt != Wj2_pt_1 && Wj2_pt != b_pt_1 &&
        b_pt   != Wj1_pt_1 && b_pt   != Wj2_pt_1 && b_pt   != b_pt_1
        )
        )
      ) {
        if ( MVA_med_mvaOutput < MVA_2_mvaOutput  ){
          MVA_med_mvaOutput = MVA_2_mvaOutput;
          MVA_med_truth     = MVA_2_truth;
          MVA_med_genTopPt  = MVA_2_genTopPt;
          MVA_med_recTopPt  = MVA_2_recTopPt;
          Wj1_pt_med = Wj1_pt_2;
          Wj2_pt_med = Wj2_pt_2;
          b_pt_med = b_pt_2;
          //std::cout << "change 8" << std::endl;
        }
        ///
        MVA_2_mvaOutput = MVA_1_mvaOutput;
        MVA_2_truth     = MVA_1_truth;
        MVA_2_genTopPt  = MVA_1_genTopPt;
        MVA_2_recTopPt  = MVA_1_recTopPt;
        Wj1_pt_2 = Wj1_pt_1;
        Wj2_pt_2 = Wj2_pt_1;
        b_pt_2 = b_pt_1;
        //std::cout << "change B1" << std::endl;
    } else if ( MVA_2_mvaOutput < MVA_1_mvaOutput ) {
        // chck med
        if ( Wj1_pt == -1.0 ||
          ( //
        Wj1_pt_med != Wj1_pt && Wj1_pt_med != Wj2_pt && Wj1_pt_med != b_pt_1 &&
        Wj2_pt_med != Wj1_pt && Wj2_pt_med != Wj2_pt && Wj2_pt_med != b_pt_1 &&
        b_pt_med   != Wj1_pt && b_pt_med   != Wj2_pt && b_pt_med   != b_pt_1
        )
        ) {
          MVA_2_mvaOutput = MVA_med_mvaOutput;
          MVA_2_truth     = MVA_med_truth;
          MVA_2_genTopPt  = MVA_med_genTopPt;
          MVA_2_recTopPt  = MVA_med_recTopPt;
          Wj1_pt_2 = Wj1_pt_med;
          Wj2_pt_2 = Wj2_pt_med;
          b_pt_2 = b_pt_med;
          ///
          MVA_med_mvaOutput = MVA_1_mvaOutput;
          MVA_med_truth     = MVA_1_truth;
          MVA_med_genTopPt  = MVA_1_genTopPt;
          MVA_med_recTopPt  = MVA_1_recTopPt;
          Wj1_pt_med = Wj1_pt_1;
          Wj2_pt_med = Wj2_pt_1;
          b_pt_med = b_pt_1;
          //std::cout << "change 7" << std::endl;
        } else if ( !( Wj1_pt == -1.0 ||
          ( // not 2 +! new
        Wj1_pt_2 != Wj1_pt && Wj1_pt_2 != Wj2_pt && Wj1_pt != b_pt_1 &&
        Wj2_pt_2 != Wj1_pt && Wj2_pt_2 != Wj2_pt && Wj2_pt != b_pt_1 &&
        b_pt_2   != Wj1_pt && b_pt_2   != Wj2_pt && b_pt   != b_pt_1
        ))
      ) {
        MVA_2_mvaOutput = 0.;
        MVA_2_truth     = false;
        MVA_2_genTopPt  = 0.;
        MVA_2_recTopPt  = 0.;
        Wj1_pt_2 = 0.;
        Wj2_pt_2 = 0.;
        b_pt_2 = 0.;
        // the second HTT does not pass sorting + mass cuts if keeping the first HTT as the highest
        //std::cout << "change 15" << std::endl;
      } //else std::cout << "change B fail : " << MVA << " " << MVA_1_mvaOutput << " " << MVA_2_mvaOutput << " " << MVA_med_mvaOutput << std::endl;
    } //
    MVA_1_truth = isGenMatched;
    MVA_1_mvaOutput = MVA;
    MVA_1_recTopPt = recTopPt;
    MVA_1_genTopPt = genTopPt_teste;
    Wj1_pt_1 = Wj1_pt;
    Wj2_pt_1 = Wj2_pt;
    b_pt_1 = b_pt; //--> new also has to be different of 2
    //std::cout << "change B" << std::endl;
  } else if ( MVA > MVA_1_mvaOutput && MVA_1_mvaOutput < MVA_med_mvaOutput && MVA < MVA_med_mvaOutput )
  {
    //std::cout << "change 6 : " << MVA << " " << MVA_1_mvaOutput << " " << MVA_2_mvaOutput << " " << MVA_med_mvaOutput << std::endl;
      double MVA_ex_mvaOutput = MVA_1_mvaOutput;
      double MVA_ex_truth     = MVA_1_truth;
      double MVA_ex_genTopPt  = MVA_1_genTopPt;
      double MVA_ex_recTopPt  = MVA_1_recTopPt;
      double Wj1_pt_ex = Wj1_pt_1;
      double Wj2_pt_ex = Wj2_pt_1;
      double b_pt_ex = b_pt_1;
      // check med

      MVA_1_mvaOutput = MVA_med_mvaOutput;
      MVA_1_truth     = MVA_med_truth;
      MVA_1_genTopPt  = MVA_med_genTopPt;
      MVA_1_recTopPt  = MVA_med_recTopPt;
      Wj1_pt_1 = Wj1_pt_med;
      Wj2_pt_1 = Wj2_pt_med;
      b_pt_1 = b_pt_med;

    if ( ( MVA_2_mvaOutput < MVA ) &&
      (
        Wj1_pt == -1.0 || (
        Wj1_pt != Wj1_pt_1 && Wj1_pt != Wj2_pt_1 && Wj1_pt != b_pt_1 &&
        Wj2_pt != Wj1_pt_1 && Wj2_pt != Wj2_pt_1 && Wj2_pt != b_pt_1 &&
        b_pt   != Wj1_pt_1 && b_pt   != Wj2_pt_1 && b_pt   != b_pt_1
        )
        )
      ) {

        if ( MVA_2_mvaOutput > MVA_ex_mvaOutput ) {
          MVA_med_mvaOutput = MVA_2_mvaOutput;
          MVA_med_truth     = MVA_2_truth;
          MVA_med_genTopPt  = MVA_2_genTopPt;
          MVA_med_recTopPt  = MVA_2_recTopPt;
          Wj1_pt_med = Wj1_pt_2;
          Wj2_pt_med = Wj2_pt_2;
          b_pt_med = b_pt_2;
          //std::cout << "change 4A" << std::endl;
        } else {
          MVA_med_mvaOutput = MVA_ex_mvaOutput;
          MVA_med_truth     = MVA_ex_truth;
          MVA_med_genTopPt  = MVA_ex_genTopPt;
          MVA_med_recTopPt  = MVA_ex_recTopPt;
          Wj1_pt_med = Wj1_pt_ex;
          Wj2_pt_med = Wj2_pt_ex;
          b_pt_med = b_pt_ex;
          //std::cout << "change 4B" << std::endl;
        }
        MVA_2_truth = isGenMatched;
        MVA_2_mvaOutput = MVA;
        MVA_2_recTopPt = recTopPt;
        MVA_2_genTopPt = genTopPt_teste;
        Wj1_pt_2 = Wj1_pt;
        Wj2_pt_2 = Wj2_pt;
        b_pt_2 = b_pt;
        //std::cout << "change 9" << std::endl;
    } else  { //(mva2 > mva)
      if ( MVA <  MVA_2_mvaOutput ){
        MVA_med_truth = isGenMatched;
        MVA_med_mvaOutput = MVA;
        MVA_med_recTopPt = recTopPt;
        MVA_med_genTopPt = genTopPt_teste;
        Wj1_pt_med = Wj1_pt;
        Wj2_pt_med = Wj2_pt;
        b_pt_med = b_pt;
        //std::cout << "change 5" << std::endl;
      } else if ( MVA >  MVA_ex_mvaOutput && MVA >  MVA_2_mvaOutput &&
          (Wj1_pt == -1.0 || (
            Wj1_pt != Wj1_pt_1 && Wj1_pt != Wj2_pt_1 && Wj1_pt != b_pt_1 &&
            Wj2_pt != Wj1_pt_1 && Wj2_pt != Wj2_pt_1 && Wj2_pt != b_pt_1 &&
            b_pt   != Wj1_pt_1 && b_pt   != Wj2_pt_1 && b_pt   != b_pt_1
            )
            )
          ) {
              MVA_2_truth = isGenMatched;
              MVA_2_mvaOutput = MVA;
              MVA_2_recTopPt = recTopPt;
              MVA_2_genTopPt = genTopPt_teste;
              Wj1_pt_2 = Wj1_pt;
              Wj2_pt_2 = Wj2_pt;
              b_pt_2 = b_pt;
              //std::cout << "change 10" << std::endl;
              if ( MVA_ex_mvaOutput > MVA_2_mvaOutput ) {
                MVA_med_mvaOutput = MVA_ex_mvaOutput;
                MVA_med_truth     = MVA_ex_truth;
                MVA_med_genTopPt  = MVA_ex_genTopPt;
                MVA_med_recTopPt  = MVA_ex_recTopPt;
                Wj1_pt_med = Wj1_pt_ex;
                Wj2_pt_med = Wj2_pt_ex;
                b_pt_med = b_pt_ex;
                //std::cout << "change 11" << std::endl;
              } else {
                MVA_med_mvaOutput = MVA_2_mvaOutput;
                MVA_med_truth     = MVA_2_truth;
                MVA_med_genTopPt  = MVA_2_genTopPt;
                MVA_med_recTopPt  = MVA_2_recTopPt;
                Wj1_pt_med = Wj1_pt_2;
                Wj2_pt_med = Wj2_pt_2;
                b_pt_med = b_pt_2;
                //std::cout << "change 4B" << std::endl;
              }

      }
    } //else std::cout << MVA_1_mvaOutput << " " << MVA_2_mvaOutput << " " << MVA_med_mvaOutput << std::endl;

  }

    return;

}
