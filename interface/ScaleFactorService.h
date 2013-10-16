#pragma once

#include <map>
#include <vector>

#include <boost/algorithm/string.hpp>

#include <CondFormats/PhysicsToolsObjects/interface/BinningPointByMap.h>

#include <FWCore/Framework/interface/EventSetup.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>

#include <RecoBTag/PerformanceDB/interface/BtagPerformance.h>
#include <RecoBTag/Records/interface/BTagPerformanceRecord.h>

#include <Extractors/PatExtractor/interface/ScaleFactor.h>

#include <FWCore/Utilities/interface/Exception.h>

class ScaleFactorService {
  public:
    enum WorkingPoint {
      LOOSE,
      TIGHT
    };

  private:
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

  public:
    ScaleFactorService(const edm::ParameterSet& settings) {
      parseMuonScaleFactors(settings);
      parseElectronScaleFactors(settings);
      mBTagPerfInit = false;
    }

    void prepareBTaggingScaleFactors(const edm::EventSetup& iSetup) {
      if ( mBTagPerfInit)
        return;

      try {
        iSetup.get<BTagPerformanceRecord>().get("MUJETSWPBTAGCSVM", mBTagPerf);
        mBTagPerfInit = mBTagPerf.isValid();
      } catch (const cms::Exception& ex) {
        mBTagPerfInit = false;
      }
    }

    ScaleFactor getElectronScaleFactor(WorkingPoint effWp, WorkingPoint isoWp, double pt, double eta) {
      return getScaleFactorFromMap(effWp, isoWp, pt, eta, mElectronScaleFactors);
    }

    ScaleFactor getMuonScaleFactor(WorkingPoint effWp, WorkingPoint isoWp, double pt, double eta) const {
      return getScaleFactorFromMap(effWp,  isoWp, pt, eta, mMuonScaleFactors);
    }

    ScaleFactor getBTaggingScaleFactor(double et, double eta) const {

      if (!mBTagPerfInit)
        return ScaleFactor();

      eta = fabs(eta);

      BinningPointByMap measurePoint;
      measurePoint.insert(BinningVariables::JetEt, et);
      measurePoint.insert(BinningVariables::JetEta, eta);

      if (mBTagPerf->isResultOk(PerformanceResult::BTAGBEFFCORR, measurePoint) && mBTagPerf->isResultOk(PerformanceResult::BTAGBERRCORR, measurePoint)) {
        ScaleFactor sf(mBTagPerf->getResult(PerformanceResult::BTAGBEFFCORR, measurePoint), mBTagPerf->getResult(PerformanceResult::BTAGBERRCORR, measurePoint));

        return sf;
      }

      return ScaleFactor();
    }

    ScaleFactor getHLTScaleFactor(double pt, double eta) {
      return {1., 0, 0};
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

    std::string workingPointToString(WorkingPoint wp) const {
      switch (wp) {
        case LOOSE:
        return "loose";

        case TIGHT:
        return "tight";
      }

      return "";
    }

    std::string nameFromWorkingPoints(WorkingPoint eff, WorkingPoint iso) const {
      return workingPointToString(eff) + "eff_" + workingPointToString(iso) + "iso";
    }

    ScaleFactorMap mMuonScaleFactors;
    ScaleFactorMap mElectronScaleFactors;

    // BTag
    bool mBTagPerfInit;
    edm::ESHandle<BtagPerformance> mBTagPerf;
    float mBTagEfficiency;
};
