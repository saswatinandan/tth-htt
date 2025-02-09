#ifndef tthAnalysis_HiggsToTauTau_RunLumiEventSelector_h
#define tthAnalysis_HiggsToTauTau_RunLumiEventSelector_h

/** \class RunLumiEventSelector
 *
 * Select events based on run + luminosity section + event number pairs
 * written (a three columns separated by white-space character) into an ASCII file
 * 
 * \author Christian Veelken, Tallinn
 *
 */

#include "tthAnalysis/HiggsToTauTau/interface/EventInfo.h" // EventInfo

#include <FWCore/ParameterSet/interface/ParameterSet.h> // edm::ParameterSet

#include <set> // std::set<>

class RunLumiEventSelector 
{
public:

  // constructor 
  explicit RunLumiEventSelector(const edm::ParameterSet & cfg);
  explicit RunLumiEventSelector(const std::string & inputFileName,
                                const std::string & separator = ":");

  // destructor
  virtual ~RunLumiEventSelector();

  bool
  operator()(ULong_t run,
             ULong_t ls,
             ULong_t event) const;

  bool
  operator()(const EventInfo & info) const;

  bool
  areWeDone() const;

private:

//--- read ASCII file containing run and event numbers
  void readInputFile();
  
  std::string inputFileName_;
  std::string separator_;

  typedef ULong_t   RunType;
  typedef ULong_t   LumiSectionType;
  typedef ULong64_t EventType;

  typedef std::set<EventType> eventNumberSet;
  typedef std::map<LumiSectionType, eventNumberSet> lumiSectionEventNumberMap;
  std::map<RunType, lumiSectionEventNumberMap> runLumiSectionEventNumbers_;
  
  typedef std::map<EventType, int> matchedEventNumbersMap;
  typedef std::map<LumiSectionType, matchedEventNumbersMap> matchedLumiSectionEventNumberMap;
  mutable std::map<RunType, matchedLumiSectionEventNumberMap> matchedRunLumiSectionEventNumbers_;

  mutable long numEventsProcessed_;
  mutable long numEventsToBeSelected_;
  mutable long numEventsSelected_;
};

#endif // tthAnalysis_HiggsToTauTau_RunLumiEventSelector_h
