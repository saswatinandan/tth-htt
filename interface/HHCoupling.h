#ifndef tthAnalysis_HiggsToTauTau_HHCoupling_h
#define tthAnalysis_HiggsToTauTau_HHCoupling_h

#include <string> // std::string
#include <vector> // std::vector<>
#include <map> // std::map<,>
#include <iosfwd> // std::ostream

class HHCoupling
{
public:
  HHCoupling(double kl = 1.,
             double kt = 1.,
             double c2 = 0.,
             double cg = 0.,
             double c2g = 0.,
             const std::string & name = "SM",
             const std::string & training = "SM");

  HHCoupling(const std::vector<double> & couplings,
             const std::string & name = "SM",
             const std::string & training = "SM");

  HHCoupling(const std::map<std::string, double> & couplings,
             const std::string & name = "SM",
             const std::string & training = "SM");

  double
  kl() const;

  double
  kt() const;

  double
  c2() const;

  double
  cg() const;

  double
  c2g() const;

  std::string
  name() const;

  std::string
  weightName() const;

  std::string
  training() const;

  friend std::ostream &
  operator<<(std::ostream & os,
             const HHCoupling & coupling);

protected:
  double kl_;
  double kt_;
  double c2_;
  double cg_;
  double c2g_;

  std::string name_;
  std::string training_;
};

#endif // tthAnalysis_HiggsToTauTau_HHCoupling_h
