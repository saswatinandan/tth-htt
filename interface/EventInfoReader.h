#ifndef EventInfoReader_H
#define EventInfoReader_H

#include "tthAnalysis/HiggsToTauTau/interface/ReaderBase.h" // ReaderBase
#include "tthAnalysis/HiggsToTauTau/interface/sysUncertOptions.h" // PUsys

// forward declarations
class TTree;
class EventInfo;

class EventInfoReader
  : public ReaderBase
{
public:
  explicit EventInfoReader(PUsys puSys_option = PUsys::central,
                           bool read_genHiggsDecayMode = true,
                           bool read_puWeight = true);
  explicit EventInfoReader(EventInfo * info,
                           PUsys puSys_option = PUsys::central,
                           bool read_genHiggsDecayMode = true,
                           bool read_puWeight = true);
  ~EventInfoReader() {}

  void
  setBranchAddresses(TTree * tree) override;

  void
  setEventInfo(EventInfo * info);

protected:
  bool read_genHiggsDecayMode_;
  bool read_puWeight_;

  EventInfo * info_;

public:
  const std::string branchName_run;
  const std::string branchName_lumi;
  const std::string branchName_event;
  const std::string branchName_genHiggsDecayMode;
  const std::string branchName_genWeight;
  const std::string branchName_puWeight;
  const std::string branchName_lheWeightSM;

  const std::string branchName_lheWeight_kt_m3p0_kv_1p0;
  const std::string branchName_lheWeight_kt_m2p0_kv_1p0;
  const std::string branchName_lheWeight_kt_m1p5_kv_1p0;
  const std::string branchName_lheWeight_kt_m1p25_kv_1p0;
  const std::string branchName_lheWeight_kt_m0p75_kv_1p0;
  const std::string branchName_lheWeight_kt_m0p5_kv_1p0;
  const std::string branchName_lheWeight_kt_m0p25_kv_1p0;
  const std::string branchName_lheWeight_kt_0p0_kv_1p0;
  const std::string branchName_lheWeight_kt_0p25_kv_1p0;
  const std::string branchName_lheWeight_kt_0p5_kv_1p0;
  const std::string branchName_lheWeight_kt_0p75_kv_1p0;
  const std::string branchName_lheWeight_kt_1p0_kv_1p0;
  const std::string branchName_lheWeight_kt_1p25_kv_1p0;
  const std::string branchName_lheWeight_kt_1p5_kv_1p0;
  const std::string branchName_lheWeight_kt_2p0_kv_1p0;
  const std::string branchName_lheWeight_kt_3p0_kv_1p0;
  //
  const std::string branchName_lheWeight_kt_m3p0_kv_1p5;
  const std::string branchName_lheWeight_kt_m2p0_kv_1p5;
  const std::string branchName_lheWeight_kt_m1p5_kv_1p5;
  const std::string branchName_lheWeight_kt_m1p25_kv_1p5;
  const std::string branchName_lheWeight_kt_m0p75_kv_1p5;
  const std::string branchName_lheWeight_kt_m0p5_kv_1p5;
  const std::string branchName_lheWeight_kt_m0p25_kv_1p5;
  const std::string branchName_lheWeight_kt_0p0_kv_1p5;
  const std::string branchName_lheWeight_kt_0p25_kv_1p5;
  const std::string branchName_lheWeight_kt_0p5_kv_1p5;
  const std::string branchName_lheWeight_kt_0p75_kv_1p5;
  const std::string branchName_lheWeight_kt_1p0_kv_1p5;
  const std::string branchName_lheWeight_kt_1p25_kv_1p5;
  const std::string branchName_lheWeight_kt_1p5_kv_1p5;
  const std::string branchName_lheWeight_kt_2p0_kv_1p5;
  const std::string branchName_lheWeight_kt_3p0_kv_1p5;
  //
  const std::string branchName_lheWeight_kt_m3p0_kv_0p5;
  const std::string branchName_lheWeight_kt_m2p0_kv_0p5;
  const std::string branchName_lheWeight_kt_m1p5_kv_0p5;
  const std::string branchName_lheWeight_kt_m1p25_kv_0p5;
  const std::string branchName_lheWeight_kt_m0p75_kv_0p5;
  const std::string branchName_lheWeight_kt_m0p5_kv_0p5;
  const std::string branchName_lheWeight_kt_m0p25_kv_0p5;
  const std::string branchName_lheWeight_kt_0p0_kv_0p5;
  const std::string branchName_lheWeight_kt_0p25_kv_0p5;
  const std::string branchName_lheWeight_kt_0p5_kv_0p5;
  const std::string branchName_lheWeight_kt_0p75_kv_0p5;
  const std::string branchName_lheWeight_kt_1p0_kv_0p5;
  const std::string branchName_lheWeight_kt_1p25_kv_0p5;
  const std::string branchName_lheWeight_kt_1p5_kv_0p5;
  const std::string branchName_lheWeight_kt_2p0_kv_0p5;
  const std::string branchName_lheWeight_kt_3p0_kv_0p5;

};

#endif // EventInfoReader_H
