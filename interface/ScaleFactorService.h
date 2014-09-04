#pragma once

#include <map>
#include <vector>
#include <memory>
#include <tuple>

#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>

#include <CondFormats/PhysicsToolsObjects/interface/BinningPointByMap.h>

#include <FWCore/Framework/interface/EventSetup.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>

#include <RecoBTag/PerformanceDB/interface/BtagPerformance.h>
#include <RecoBTag/Records/interface/BTagPerformanceRecord.h>

#include <Extractors/PatExtractor/interface/ScaleFactor.h>

#include <FWCore/Utilities/interface/Exception.h>


#include <TString.h>
#include <TF1.h>

class BTagScaleFactorService {
  public:
    BTagScaleFactorService(const std::string& jet_type) {
      m_jet_type = jet_type;
    }

    virtual void init(const edm::EventSetup& iSetup) {
      // Do nothing
    }

    virtual ScaleFactor getScaleFactor(float et, float eta) = 0;

  protected:
    std::string m_jet_type;
};

class BTagScaleFactorFromGlobalTagService: public BTagScaleFactorService {
  public:
    BTagScaleFactorFromGlobalTagService(const std::string& jet_type, const std::string& payload):
      BTagScaleFactorService(jet_type) {
      m_payload = payload;
      m_init = false;

      std::cout << "Reading scale factors from Global Tag for " << jet_type << " jets" << std::endl;
    }

    virtual void init(const edm::EventSetup& iSetup) {
      try {
        iSetup.get<BTagPerformanceRecord>().get(m_payload, mBTagPerf);
        m_init = mBTagPerf.isValid();
      } catch (const cms::Exception& ex) {
        m_init = false;
        std::cout << "WARNING: an error occured when loading scale factors from Global Tag. Please check that '" << m_payload << "' is a correct payload" << std::endl;
        std::cout << ex << std::endl;
      }

      if (m_init) {
        std::cout << "B-tagging scale factors correctly loaded from Global Tag" << std::endl;
      } else {
      }
    }

    virtual ScaleFactor getScaleFactor(float et, float eta) {

      if (!m_init)
        return ScaleFactor();

      eta = fabs(eta);

      // Double errors if pt < 20 || pt > 800
      bool doubleErrors = et < 20 || et > 800;
      if (et < 20)
        et = 20;

      if (et > 800)
        et = 800;

      BinningPointByMap measurePoint;
      measurePoint.insert(BinningVariables::JetEt, et);
      measurePoint.insert(BinningVariables::JetEta, eta);

      if (mBTagPerf->isResultOk(PerformanceResult::BTAGBEFFCORR, measurePoint) && mBTagPerf->isResultOk(PerformanceResult::BTAGBERRCORR, measurePoint)) {
        ScaleFactor sf(mBTagPerf->getResult(PerformanceResult::BTAGBEFFCORR, measurePoint), mBTagPerf->getResult(PerformanceResult::BTAGBERRCORR, measurePoint));

        if (doubleErrors)
          return {sf.getValue(), 2 * sf.getErrorHigh(), 2 * sf.getErrorLow()};

        return sf;
      }

      return ScaleFactor();
    }

  private:
    std::string m_payload;
    bool m_init;
    edm::ESHandle<BtagPerformance> mBTagPerf;
};

class BTagScaleFactorFromVSetService: public BTagScaleFactorService {
  public:
    BTagScaleFactorFromVSetService(const std::string& jet_type, const edm::VParameterSet& sfs):
      BTagScaleFactorService(jet_type) {

      int nEta = 0;
      for (const edm::ParameterSet& etaBin: sfs) {
        nEta++;
        const std::vector<double> etaBinVector = etaBin.getParameter<std::vector<double>>("eta");
        auto eta = std::make_pair(etaBinVector[0], etaBinVector[1]);
        float ptMax = ((eta.first - 1.5 < 1e-6) || (eta.first - 1.6 < 1e-6)) ? 850. : 1000;
        auto functions = std::make_tuple(
            new TF1(TString::Format("%s_nominal_%f_%f", jet_type.c_str(), eta.first, eta.second), etaBin.getParameter<std::string>("value").c_str(), 20, ptMax),
            new TF1(TString::Format("%s_max_%f_%f", jet_type.c_str(), eta.first, eta.second), etaBin.getParameter<std::string>("error_high").c_str(), 20, ptMax),
            new TF1(TString::Format("%s_min_%f_%f", jet_type.c_str(), eta.first, eta.second), etaBin.getParameter<std::string>("error_low").c_str(), 20, ptMax)
            );

        m_scaleFactors[eta] = functions;
      }

      std::cout << nEta << " eta bins loaded for " << jet_type << " jets." << std::endl;
    }

    virtual ~BTagScaleFactorFromVSetService() {
      for (auto& scaleFactor: m_scaleFactors) {
        delete std::get<0>(scaleFactor.second);
        delete std::get<1>(scaleFactor.second);
        delete std::get<2>(scaleFactor.second);
      }
    }

    virtual ScaleFactor getScaleFactor(float et, float eta) {
      eta = fabs(eta);
      float ptMax = getPtMax(eta);

      // Double errors if pt < 20 or pt > ptMax
      float errorFactor = (et < 20 || et > ptMax) ? 2. : 1.;
      if (et > ptMax)
        et = ptMax;
      else if (et < 20)
        et = 20;

      for (const auto& etaBin: m_scaleFactors) {
        if (eta >= etaBin.first.first && eta < etaBin.first.second) {
          TF1* nominal = std::get<0>(etaBin.second);
          TF1* error_high = std::get<1>(etaBin.second);
          TF1* error_low = std::get<2>(etaBin.second);

          float nominal_value = nominal->Eval(et);
          float error_high_value = error_high->Eval(et);
          float error_low_value = error_low->Eval(et);

          return ScaleFactor(nominal_value, errorFactor * (error_high_value - nominal_value), errorFactor * (nominal_value - error_low_value));
        }
      }

      return ScaleFactor();
    }

  private:
    std::map<
      std::pair<float, float>, // Name
      std::tuple<TF1*, TF1*, TF1*> // Functions (nominal, up, down)
    > m_scaleFactors;

    float getPtMax(float eta) const {
      if (fabs(eta) >= 1.6)
        return 850;
      else
        return 1000;
    }
};

class ScaleFactorService {
  public:
    enum WorkingPoint {
      LOOSE,
      TIGHT
    };

    enum Flavor {
      B,
      C,
      LIGHT
    };

  private:
    friend class PatExtractor;
    static std::shared_ptr<ScaleFactorService> m_instance;
    static void createInstance(const edm::ParameterSet& settings) {
      m_instance.reset( new ScaleFactorService(settings) );
    }

    typedef std::map<
      std::string,
      std::map<
        // Eta binning
        std::pair<double, double>,
        std::map<
          // Pt binning
          std::pair<double, double>,
          ScaleFactor
        >
      >
    > ScaleFactorMap;

    typedef std::map<
      std::string, // Name
      std::shared_ptr<BTagScaleFactorService>
    > BTagScaleFactorMap;

    ScaleFactorService(const edm::ParameterSet& settings) {
      parseMuonScaleFactors(settings);
      parseElectronScaleFactors(settings);
      parseBTagScaleFactors(settings);
    }

    ScaleFactorService(ScaleFactorService const&); // Don't Implement
    void operator=(ScaleFactorService const&); // Don't implement

  public:

    static ScaleFactorService& getInstance() {
      return *m_instance;
    }

    void prepareBTaggingScaleFactors(const edm::EventSetup& iSetup) {
      for (auto& btagSF: mBTagScaleFactors) {
        btagSF.second->init(iSetup);
      }
    }

    ScaleFactor getElectronScaleFactor(WorkingPoint effWp, WorkingPoint isoWp, double pt, double eta) {
      return getScaleFactorFromMap(effWp, isoWp, pt, eta, mElectronScaleFactors);
    }

    ScaleFactor getMuonScaleFactor(WorkingPoint effWp, WorkingPoint isoWp, double pt, double eta) const {
      return getScaleFactorFromMap(effWp,  isoWp, pt, eta, mMuonScaleFactors);
    }

    ScaleFactor getBTaggingScaleFactor(Flavor flavor, double et, double eta) const {

      const std::string flavorName = flavorToString(flavor);
      if (mBTagScaleFactors.count(flavorName) == 0)
        return ScaleFactor();

      ScaleFactor sf = mBTagScaleFactors.at(flavorName)->getScaleFactor(et, eta);

      if (flavor == C) {
        // For C jets, SFc = SFb, but errors are twice larger
        return {sf.getValue(), 2 * sf.getErrorHigh(), 2 * sf.getErrorLow() };
      } else {
        return sf;
      }
    }

    ScaleFactor getHLTScaleFactor(double pt, double eta) {
      return {1., 0, 0};
    }

    // First selection working point, Second isolation working pointc
    std::vector<std::pair<WorkingPoint, WorkingPoint>> getElectronScaleFactorWorkingPoints() const {
      return getWorkingPointsFromMap(mElectronScaleFactors);
    }

    // First selection working point, Second isolation working pointc
    std::vector<std::pair<WorkingPoint, WorkingPoint>> getMuonScaleFactorWorkingPoints() const {
      return getWorkingPointsFromMap(mMuonScaleFactors);
    }

    static std::string nameFromWorkingPoints(WorkingPoint eff, WorkingPoint iso) {
      return workingPointToString(eff) + "eff_" + workingPointToString(iso) + "iso";
    }

    static std::string workingPointToString(WorkingPoint wp) {
      switch (wp) {
        case LOOSE:
        return "loose";

        case TIGHT:
        return "tight";
      }

      return "";
    }

    static std::pair<WorkingPoint, WorkingPoint> getWorkingPointFromName(const std::string& name) {
      static boost::regex regex("([a-z]+)eff_([a-z]+)iso");
      boost::smatch result;

      if (boost::regex_search(name, result, regex)) {
        return std::make_pair(stringToWorkingPoint(std::string(result[1].first, result[1].second)), stringToWorkingPoint(std::string(result[2].first, result[2].second)));
      }

      assert(false);
      return std::make_pair(LOOSE, LOOSE);
    }

  private:
    ScaleFactor getScaleFactorFromMap(WorkingPoint effWp, WorkingPoint isoWp, double pt, double eta, const ScaleFactorMap& map) const {
      const std::string wpName = nameFromWorkingPoints(effWp, isoWp);
      if (map.count(wpName) == 0)
        return ScaleFactor();

      eta = fabs(eta);
      for (const auto& etaBin: map.at(wpName)) {
        if (eta >= etaBin.first.first && eta < etaBin.first.second) {
          // Look for pt bin
          for (const auto& ptBin: etaBin.second) {
            if (pt >= ptBin.first.first && pt < ptBin.first.second) {
              return ptBin.second;
            }
          }
        }
      }

      return ScaleFactor();
    }

    void parseScaleFactors(const edm::ParameterSet& settings, const std::string& setName, ScaleFactorMap& map, const std::string& wpName) {
      uint32_t nEta = 0, nPt = 0;
      const edm::VParameterSet& sfs = settings.getParameterSetVector(setName);

      for (const edm::ParameterSet& etaBin: sfs) {
        nEta++;
        nPt = 0;
        const std::vector<double> etaBinVector = etaBin.getParameter<std::vector<double>>("eta");
        auto eta = std::make_pair(etaBinVector[0], etaBinVector[1]);
        for (const edm::ParameterSet& ptBin: etaBin.getParameterSetVector("SF")) {
          nPt++;
          const std::vector<double> ptBinVector = ptBin.getParameter<std::vector<double>>("pt");
          auto pt = std::make_pair(ptBinVector[0], ptBinVector[1]);

          ScaleFactor sf(ptBin.getParameter<double>("value"), ptBin.getParameter<double>("error_low"), ptBin.getParameter<double>("error_high"));
          map[wpName][eta][pt] = sf;
        }
      }

      std::cout << nEta << " eta bins and " << nPt << " pt bins loaded for " << wpName << " working points." << std::endl;
    }

    void parseMuonScaleFactors(const edm::ParameterSet& settings) {
      std::cout << "Loading muon scale factors..." << std::endl;
      const std::vector<std::string>& scaleFactors = settings.getParameter<std::vector<std::string>>("muon_scale_factors");
      for (const std::string& name: scaleFactors) {
        std::string wpName = name;
        boost::replace_all(wpName, "muon_scale_factors_", "");
        parseScaleFactors(settings, name, mMuonScaleFactors, wpName);
      }
      std::cout << std::endl;
    }

    void parseElectronScaleFactors(const edm::ParameterSet& settings) {
      std::cout << "Loading electron scale factors..." << std::endl;
      const std::vector<std::string>& scaleFactors = settings.getParameter<std::vector<std::string>>("electron_scale_factors");
      for (const std::string& name: scaleFactors) {
        std::string wpName = name;
        boost::replace_all(wpName, "electron_scale_factors_", "");
        parseScaleFactors(settings, name, mElectronScaleFactors, wpName);
      }
      std::cout << std::endl;
    }

    void parseBTagScaleFactors(const edm::ParameterSet& settings) {
      std::cout << "Loading btag scale factors..." << std::endl;
      const std::vector<std::string>& scaleFactors = settings.getParameter<std::vector<std::string>>("b_tagging_scale_factors");
      for (const std::string& name: scaleFactors) {
        const edm::ParameterSet& sf = settings.getParameterSet(name);
        std::string jet_type = sf.getParameter<std::string>("jet_type");
        bool from_globaltag = sf.getParameter<bool>("from_globaltag");
        if (from_globaltag) {
          std::string payload = sf.getParameter<std::string>("payload");
          mBTagScaleFactors[jet_type] = std::make_shared<BTagScaleFactorFromGlobalTagService>(jet_type, payload);
        } else {
          const edm::VParameterSet& sfs = sf.getParameterSetVector("scale_factors");
          mBTagScaleFactors[jet_type] = std::make_shared<BTagScaleFactorFromVSetService>(jet_type, sfs);
        }
      }

      std::cout << std::endl;
    }

    static WorkingPoint stringToWorkingPoint(const std::string& name) {

      if (name == "loose")
        return LOOSE;

      if (name == "tight")
        return TIGHT;

      return LOOSE;
    }

    static std::string flavorToString(Flavor flavor) {
      switch (flavor) {
        case B:
          return "b";

        case C:
          return "c";

        case LIGHT:
          return "light";
      }

      return "";
    }

    static std::vector<std::pair<WorkingPoint, WorkingPoint>> getWorkingPointsFromMap(const ScaleFactorMap& map) {
      std::vector<std::pair<WorkingPoint, WorkingPoint>> result;

      for (const auto& it: map) {
        result.push_back(getWorkingPointFromName(it.first));
      }

      return result;
    }

    ScaleFactorMap mMuonScaleFactors;
    ScaleFactorMap mElectronScaleFactors;
    BTagScaleFactorMap mBTagScaleFactors;

    // BTag
    float mBTagEfficiency;
};
