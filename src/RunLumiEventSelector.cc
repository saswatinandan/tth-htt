#include "tthAnalysis/HiggsToTauTau/interface/RunLumiEventSelector.h"

#include "tthAnalysis/HiggsToTauTau/interface/cmsException.h" // cmsException()

#include <boost/lexical_cast.hpp> // boost::lexical_cast<>()
#include <boost/algorithm/string/join.hpp> // boost::algorithm::join()

#include <fstream> // std::ifstream
#include <regex> // std::regex, std::regex_search(), std::smatch

RunLumiEventSelector::RunLumiEventSelector(const std::string & inputFileName,
                                           const std::string & separator)
  : inputFileName_(inputFileName)
  , separator_(separator)
  , numEventsProcessed_(0)
  , numEventsToBeSelected_(0)
  , numEventsSelected_(0)
{
  if(inputFileName_.empty())
  {
    throw cmsException(this) << "inputFileName empty";
  }

  readInputFile();
}

RunLumiEventSelector::RunLumiEventSelector(const edm::ParameterSet & cfg)
  : RunLumiEventSelector(cfg.getParameter<std::string>("inputFileName"),
                         cfg.exists("separator") ? cfg.getParameter<std::string>("separator") : "[[:space:]]+")
{}

RunLumiEventSelector::~RunLumiEventSelector()
{
  const std::string matchRemark = numEventsSelected_ == numEventsToBeSelected_ ? "matches" : "does NOT match";
  std::cout << "<RunLumiEventSelector::~RunLumiEventSelector>:" 
               " Number of Events processed = " << numEventsProcessed_ << "\n"
               " Number of Events selected = " << numEventsSelected_ << ", " << matchRemark
            << " Number of Events to be selected = " << numEventsToBeSelected_ << '\n';
 
//--- check for events specified by run + event number in ASCII file
//    and not found in EDM input .root file
  int numRunLumiSectionEventNumbersUnmatched = 0;
  for(const auto & run: matchedRunLumiSectionEventNumbers_)
  {
    for(const auto & lumiSection: run.second)
    {
      for(const auto & event: lumiSection.second)
      {
        if(event.second < 1)
        {
          if(numRunLumiSectionEventNumbersUnmatched == 0)
          {
            std::cout << "Events not found:\n";
          }
          std::cout << " run# = " << run.first << ","
                       " ls# "    << lumiSection.first << ","
                       " event# " << event.first << '\n';
          ++numRunLumiSectionEventNumbersUnmatched;
        } // if(event.second < 1)
      } // event
    } // lumiSection
  } // run

  if(numRunLumiSectionEventNumbersUnmatched > 0)
  {
    std::cout << "--> Number of unmatched Events = "
              << numRunLumiSectionEventNumbersUnmatched << '\n';
  }

//--- check for events specified by run + event number in ASCII file
//    and found more than once in EDM input .root file
  int numRunLumiSectionEventNumbersAmbiguousMatch = 0;
  for(const auto & run: matchedRunLumiSectionEventNumbers_)
  {
    for(const auto & lumiSection: run.second)
    {
      for(const auto & event: lumiSection.second)
      {
        if(event.second > 1)
        {
          if(numRunLumiSectionEventNumbersAmbiguousMatch == 0)
          {
            std::cout << "Events found more than once:\n";
          }
          std::cout << " run# = " << run.first << ","
                       " ls# "    << lumiSection.first << ","
                       " event# " << event.first << '\n';
          ++numRunLumiSectionEventNumbersAmbiguousMatch;
        } // if(event.second > 1)
      } // event
    } // lumiSection
  } // run

  if(numRunLumiSectionEventNumbersAmbiguousMatch > 0)
  {
    std::cout << "--> Number of ambiguously matched Events = "
              << numRunLumiSectionEventNumbersAmbiguousMatch << '\n';
  }
}

void
RunLumiEventSelector::readInputFile()
{
//--- read run + luminosity section + event number pairs from ASCII file
  std::ifstream inputFile(inputFileName_);
  if(inputFile)
  {
    const std::regex pattern(boost::algorithm::join(std::vector<std::string>(3, "([[:digit:]]+)"), separator_));
    std::smatch match;
    int iLine = 0;

    for(std::string line; std::getline(inputFile, line); )
    {
      ++iLine;
      if(line.empty()) continue;

      bool parseError = false;
      try
      {
        if(std::regex_search(line, match, pattern))
        {
          const RunType         runNumber         = boost::lexical_cast<RunType>        (match[1]);
          const LumiSectionType lumiSectionNumber = boost::lexical_cast<LumiSectionType>(match[2]);
          const EventType       eventNumber       = boost::lexical_cast<EventType>      (match[3]);
          std::cout << "--> adding "
                       "run# = " << runNumber         << ", "
                       "ls# "    << lumiSectionNumber << ", "
                       "event# " << eventNumber       << '\n';

          runLumiSectionEventNumbers_[runNumber][lumiSectionNumber].insert(eventNumber);
          matchedRunLumiSectionEventNumbers_[runNumber][lumiSectionNumber][eventNumber] = 0;
          ++numEventsToBeSelected_;
        }
        else
        {
          parseError = true;
        }
      }
      catch(const boost::bad_lexical_cast & e)
      {
        parseError = true;
      }

      if(parseError)
      {
        std::cerr << "<RunLumiEventSelector::readInputFile>: Error in parsing line "
                  << iLine << " = '" << line << "'" << " of input file = " << inputFileName_ << '\n';
      }
    } // line
  } // if(inputFile)
  else
  {
    throw cmsException(this, __func__) << "Could not open file " << inputFileName_ << " for reading";
  }

  if(numEventsToBeSelected_ == 0)
  {
    throw cmsException(this, __func__)
      << "Failed to read any run, lumi and event numbers from input file " << inputFileName_;
  }
}

bool
RunLumiEventSelector::operator()(ULong_t run,
                                 ULong_t ls,
                                 ULong_t event) const
{
//--- check if run number matches any of the runs containing events to be selected
  bool isSelected = false;
  if(runLumiSectionEventNumbers_.count(run))
  {
    const lumiSectionEventNumberMap & lumiSectionEventNumbers = runLumiSectionEventNumbers_.at(run);
    if(lumiSectionEventNumbers.count(ls))
    {
      const eventNumberSet & eventNumbers = lumiSectionEventNumbers.at(ls);
      if(eventNumbers.count(event))
      {
        isSelected = true;
      }
    }
  }

  ++numEventsProcessed_;
  if(isSelected)
  {
    std::cout << "<RunLumiEventSelector::operator>: selecting "
                 "run# = " << run << ", ls# " << ls << ", event# " << event << '\n';

    ++matchedRunLumiSectionEventNumbers_[run][ls][event];
    if ( matchedRunLumiSectionEventNumbers_[run][ls][event] == 1 ) { // CV: select only first match
      ++numEventsSelected_;
      return true;
    } else {
      return false;
    }
  }

  return false;
}

bool
RunLumiEventSelector::areWeDone() const
{
  return numEventsToBeSelected_ == numEventsSelected_;
}

bool
RunLumiEventSelector::operator()(const EventInfo & info) const
{
  return RunLumiEventSelector::operator()(info.run, info.lumi, info.event);
}
