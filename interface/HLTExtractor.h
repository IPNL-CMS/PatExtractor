#ifndef HLTEXTRACTOR_H
#define HLTEXTRACTOR_H

/**
 * HLTExtractor
 * \brief: Base class for extracting HLT info
 */

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include <FWCore/ParameterSet/interface/ParameterSet.h>

#include "Extractors/PatExtractor/interface/BaseExtractor.h"
#include "Extractors/PatExtractor/interface/MCExtractor.h"

//Include std C++
#include <iostream>
#include <vector>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

#include <boost/regex.hpp>

namespace tinyxml2 {
  class XMLElement;
}

template<typename T>
class Range {
  public:
    Range(T from, T to):
      mFrom(from), mTo(to) {}

    T from() const {
      return mFrom;
    }

    T to() const {
      return mTo;
    }

    bool in(T value) const {
      /*if (mFrom < 0)
        return value <= mTo;
      else if (mTo < 0)
        return value >= mFrom;
      else*/
        return value >= mFrom && value <= mTo;
    }

    bool operator<(const Range<T>& other) const {
      return mFrom < other.mFrom;
    }

    template<typename U>
    friend std::ostream& operator<<(std::ostream& stream, const Range<U>& range);

  private:
    T mFrom;
    T mTo;
};

template<typename T>
std::ostream& operator<<(std::ostream& stream, const Range<T>& range)
{
  stream << "[" << range.from() << ", " << range.to() << "]";

  return stream;
}

typedef boost::regex PathData;
typedef std::vector<PathData> PathVector;

class Triggers {
  public:
    Triggers(const std::string& xmlContent):
      mCachedRange(NULL), mCachedVector(NULL) {
        parse(xmlContent);
      }

    void print();

    //const boost::regex& getHLTPath(unsigned int run, float pt);
    const PathVector& getTriggers(unsigned int run);

  private:
    std::string mXmlFile;
    std::map<Range<unsigned int>, PathVector> mTriggers;

    const Range<unsigned int>* mCachedRange;
    const PathVector* mCachedVector;

    bool parseRunsElement(const tinyxml2::XMLElement* runs);

    bool parse(const std::string& xmlContent);
};

class HLTExtractor: public SuperBaseExtractor
{

 public:

  HLTExtractor(const std::string& name, bool doTree, const edm::ParameterSet& config);
  HLTExtractor(const std::string& name, TFile *a_file);
  ~HLTExtractor();


  virtual void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor* mcExtractor);
  void getInfo(int ievt); 

  void reset();
  void fillTree(); 

  int getSize() const;

  // Setters/Getters

  std::string paths(int i) {return m_HLT_vector->at(i);}

  std::vector<std::string>* getPaths() {
    return m_HLT_vector;
  }

  bool isTriggerFired() const {
    return m_passed;
  }

 private:
  
  TTree* m_tree_HLT;

  int m_n_HLTs;
  std::vector< std::string > *m_HLT_vector;

  bool m_filterHLT;
  std::string m_triggersXML;
  std::string* m_mustPass;
  bool m_passed;

  std::shared_ptr<Triggers> m_triggersService; 
};

#endif 
