#include "tthAnalysis/HiggsToTauTau/interface/EventInfoReader.h"

#include "tthAnalysis/HiggsToTauTau/interface/EventInfo.h" // EventInfo
#include "tthAnalysis/HiggsToTauTau/interface/BranchAddressInitializer.h" // BranchAddressInitializer, TTree, Form()

EventInfoReader::EventInfoReader(PUsys puSys_option,
                                 bool read_genHiggsDecayMode,
                                 bool read_puWeight)
  : EventInfoReader(nullptr, puSys_option, read_genHiggsDecayMode, read_puWeight)
{}

EventInfoReader::EventInfoReader(EventInfo * info,
                                 PUsys puSys_option,
                                 bool read_genHiggsDecayMode,
                                 bool read_puWeight)
  : read_genHiggsDecayMode_(read_genHiggsDecayMode)
  , read_puWeight_(read_puWeight)
  , info_(info)
  , branchName_run("run")
  , branchName_lumi("luminosityBlock")
  , branchName_event("event")
  , branchName_genHiggsDecayMode("genHiggsDecayMode")
  , branchName_genWeight("genWeight")
  , branchName_puWeight(getBranchName_pileup(puSys_option))
  , branchName_lheWeightSM("LHEWeight_rwgt_12")
  , branchName_lheWeight_kt_m3p0_kv_1p0("LHEWeight_rwgt_1")
  , branchName_lheWeight_kt_m2p0_kv_1p0("LHEWeight_rwgt_2")
  , branchName_lheWeight_kt_m1p5_kv_1p0("LHEWeight_rwgt_3")
  , branchName_lheWeight_kt_m1p25_kv_1p0("LHEWeight_rwgt_4")
  , branchName_lheWeight_kt_m0p75_kv_1p0("LHEWeight_rwgt_5")
  , branchName_lheWeight_kt_m0p5_kv_1p0("LHEWeight_rwgt_6")
  , branchName_lheWeight_kt_m0p25_kv_1p0("LHEWeight_rwgt_7")
  , branchName_lheWeight_kt_0p0_kv_1p0("LHEWeight_rwgt_8")
  , branchName_lheWeight_kt_0p25_kv_1p0("LHEWeight_rwgt_9")
  , branchName_lheWeight_kt_0p5_kv_1p0("LHEWeight_rwgt_10")
  , branchName_lheWeight_kt_0p75_kv_1p0("LHEWeight_rwgt_11")
  , branchName_lheWeight_kt_1p0_kv_1p0("LHEWeight_rwgt_12")
  , branchName_lheWeight_kt_1p25_kv_1p0("LHEWeight_rwgt_13")
  , branchName_lheWeight_kt_1p5_kv_1p0("LHEWeight_rwgt_14")
  , branchName_lheWeight_kt_2p0_kv_1p0("LHEWeight_rwgt_15")
  , branchName_lheWeight_kt_3p0_kv_1p0("LHEWeight_rwgt_16")
  //
  , branchName_lheWeight_kt_m3p0_kv_1p5("LHEWeight_rwgt_17")
  , branchName_lheWeight_kt_m2p0_kv_1p5("LHEWeight_rwgt_18")
  , branchName_lheWeight_kt_m1p5_kv_1p5("LHEWeight_rwgt_19")
  , branchName_lheWeight_kt_m1p25_kv_1p5("LHEWeight_rwgt_20")
  , branchName_lheWeight_kt_m0p75_kv_1p5("LHEWeight_rwgt_21")
  , branchName_lheWeight_kt_m0p5_kv_1p5("LHEWeight_rwgt_22")
  , branchName_lheWeight_kt_m0p25_kv_1p5("LHEWeight_rwgt_23")
  , branchName_lheWeight_kt_0p0_kv_1p5("LHEWeight_rwgt_24")
  , branchName_lheWeight_kt_0p25_kv_1p5("LHEWeight_rwgt_25")
  , branchName_lheWeight_kt_0p5_kv_1p5("LHEWeight_rwgt_26")
  , branchName_lheWeight_kt_0p75_kv_1p5("LHEWeight_rwgt_27")
  , branchName_lheWeight_kt_1p0_kv_1p5("LHEWeight_rwgt_28")
  , branchName_lheWeight_kt_1p25_kv_1p5("LHEWeight_rwgt_29")
  , branchName_lheWeight_kt_1p5_kv_1p5("LHEWeight_rwgt_30")
  , branchName_lheWeight_kt_2p0_kv_1p5("LHEWeight_rwgt_31")
  , branchName_lheWeight_kt_3p0_kv_1p5("LHEWeight_rwgt_32")
  //
  , branchName_lheWeight_kt_m3p0_kv_0p5("LHEWeight_rwgt_33")
  , branchName_lheWeight_kt_m2p0_kv_0p5("LHEWeight_rwgt_34")
  , branchName_lheWeight_kt_m1p5_kv_0p5("LHEWeight_rwgt_35")
  , branchName_lheWeight_kt_m1p25_kv_0p5("LHEWeight_rwgt_36")
  , branchName_lheWeight_kt_m0p75_kv_0p5("LHEWeight_rwgt_37")
  , branchName_lheWeight_kt_m0p5_kv_0p5("LHEWeight_rwgt_38")
  , branchName_lheWeight_kt_m0p25_kv_0p5("LHEWeight_rwgt_39")
  , branchName_lheWeight_kt_0p0_kv_0p5("LHEWeight_rwgt_40")
  , branchName_lheWeight_kt_0p25_kv_0p5("LHEWeight_rwgt_41")
  , branchName_lheWeight_kt_0p5_kv_0p5("LHEWeight_rwgt_42")
  , branchName_lheWeight_kt_0p75_kv_0p5("LHEWeight_rwgt_43")
  , branchName_lheWeight_kt_1p0_kv_0p5("LHEWeight_rwgt_44")
  , branchName_lheWeight_kt_1p25_kv_0p5("LHEWeight_rwgt_45")
  , branchName_lheWeight_kt_1p5_kv_0p5("LHEWeight_rwgt_46")
  , branchName_lheWeight_kt_2p0_kv_0p5("LHEWeight_rwgt_47")
  , branchName_lheWeight_kt_3p0_kv_0p5("LHEWeight_rwgt_48")
{}

void
EventInfoReader::setBranchAddresses(TTree * tree)
{
  BranchAddressInitializer bai(tree);
  bai.setBranchAddress(info_ -> run, branchName_run);
  bai.setBranchAddress(info_ -> lumi, branchName_lumi);
  bai.setBranchAddress(info_ -> event, branchName_event);
  if(info_ -> is_signal() && read_genHiggsDecayMode_)
  {
    bai.setBranchAddress(info_ -> genHiggsDecayMode, branchName_genHiggsDecayMode);
  }
  if(info_ -> is_mc())
  {
    bai.setBranchAddress(info_ -> genWeight, branchName_genWeight);
    if (read_puWeight_)
    {
      bai.setBranchAddress(info_ -> pileupWeight, branchName_puWeight);
    }
  }
  if(info_ -> is_mc_th())
  {
    bai.setBranchAddress(info_ -> genWeight_tH, branchName_lheWeightSM);

    bai.setBranchAddress(info_ -> genWeight_tH_kt_m3p0_kv_1p0, branchName_lheWeight_kt_m3p0_kv_1p0);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_m2p0_kv_1p0, branchName_lheWeight_kt_m2p0_kv_1p0);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_m1p5_kv_1p0, branchName_lheWeight_kt_m1p5_kv_1p0);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_m1p25_kv_1p0, branchName_lheWeight_kt_m1p25_kv_1p0);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_m0p75_kv_1p0, branchName_lheWeight_kt_m0p75_kv_1p0);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_m0p5_kv_1p0, branchName_lheWeight_kt_m0p5_kv_1p0);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_m0p25_kv_1p0, branchName_lheWeight_kt_m0p25_kv_1p0);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_0p0_kv_1p0, branchName_lheWeight_kt_0p0_kv_1p0);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_0p25_kv_1p0, branchName_lheWeight_kt_0p25_kv_1p0);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_0p5_kv_1p0, branchName_lheWeight_kt_0p5_kv_1p0);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_0p75_kv_1p0, branchName_lheWeight_kt_0p75_kv_1p0);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_1p0_kv_1p0, branchName_lheWeight_kt_1p0_kv_1p0);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_1p25_kv_1p0, branchName_lheWeight_kt_1p25_kv_1p0);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_1p5_kv_1p0, branchName_lheWeight_kt_1p5_kv_1p0);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_2p0_kv_1p0, branchName_lheWeight_kt_2p0_kv_1p0);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_3p0_kv_1p0, branchName_lheWeight_kt_3p0_kv_1p0);
    //
    bai.setBranchAddress(info_ -> genWeight_tH_kt_m3p0_kv_1p5, branchName_lheWeight_kt_m3p0_kv_1p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_m2p0_kv_1p5, branchName_lheWeight_kt_m2p0_kv_1p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_m1p5_kv_1p5, branchName_lheWeight_kt_m1p5_kv_1p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_m1p25_kv_1p5, branchName_lheWeight_kt_m1p25_kv_1p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_m0p75_kv_1p5, branchName_lheWeight_kt_m0p75_kv_1p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_m0p5_kv_1p5, branchName_lheWeight_kt_m0p5_kv_1p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_m0p25_kv_1p5, branchName_lheWeight_kt_m0p25_kv_1p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_0p0_kv_1p5, branchName_lheWeight_kt_0p0_kv_1p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_0p25_kv_1p5, branchName_lheWeight_kt_0p25_kv_1p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_0p5_kv_1p5, branchName_lheWeight_kt_0p5_kv_1p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_0p75_kv_1p5, branchName_lheWeight_kt_0p75_kv_1p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_1p0_kv_1p5, branchName_lheWeight_kt_1p0_kv_1p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_1p25_kv_1p5, branchName_lheWeight_kt_1p25_kv_1p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_1p5_kv_1p5, branchName_lheWeight_kt_1p5_kv_1p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_2p0_kv_1p5, branchName_lheWeight_kt_2p0_kv_1p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_3p0_kv_1p5, branchName_lheWeight_kt_3p0_kv_1p5);
    //
    bai.setBranchAddress(info_ -> genWeight_tH_kt_m3p0_kv_0p5, branchName_lheWeight_kt_m3p0_kv_0p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_m2p0_kv_0p5, branchName_lheWeight_kt_m2p0_kv_0p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_m1p5_kv_0p5, branchName_lheWeight_kt_m1p5_kv_0p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_m1p25_kv_0p5, branchName_lheWeight_kt_m1p25_kv_0p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_m0p75_kv_0p5, branchName_lheWeight_kt_m0p75_kv_0p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_m0p5_kv_0p5, branchName_lheWeight_kt_m0p5_kv_0p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_m0p25_kv_0p5, branchName_lheWeight_kt_m0p25_kv_0p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_0p0_kv_0p5, branchName_lheWeight_kt_0p0_kv_0p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_0p25_kv_0p5, branchName_lheWeight_kt_0p25_kv_0p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_0p5_kv_0p5, branchName_lheWeight_kt_0p5_kv_0p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_0p75_kv_0p5, branchName_lheWeight_kt_0p75_kv_0p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_1p0_kv_0p5, branchName_lheWeight_kt_1p0_kv_0p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_1p25_kv_0p5, branchName_lheWeight_kt_1p25_kv_0p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_1p5_kv_0p5, branchName_lheWeight_kt_1p5_kv_0p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_2p0_kv_0p5, branchName_lheWeight_kt_2p0_kv_0p5);
    bai.setBranchAddress(info_ -> genWeight_tH_kt_3p0_kv_0p5, branchName_lheWeight_kt_3p0_kv_0p5);
  }
}

void
EventInfoReader::setEventInfo(EventInfo * info)
{
  info_ = info;
}
