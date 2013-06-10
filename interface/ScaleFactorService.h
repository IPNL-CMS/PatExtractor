#pragma once

#include <map>
#include <vector>

#include <CondFormats/PhysicsToolsObjects/interface/BinningPointByMap.h>

#include <FWCore/Framework/interface/EventSetup.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>

#include <RecoBTag/PerformanceDB/interface/BtagPerformance.h>
#include <RecoBTag/Records/interface/BTagPerformanceRecord.h>

#include <Extractors/PatExtractor/interface/ScaleFactor.h>

#include <FWCore/Utilities/interface/Exception.h>

class ScaleFactorService {
  private:
    typedef std::map<
      // Eta binning
      std::pair<double, double>,
      std::map<
        // Pt binning
        std::pair<double, double>,
        ScaleFactor
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
	mBTagPerfInit =  mBTagPerf.isValid();
      } catch (const cms::Exception& ex) {
        mBTagPerfInit = false;
      }
    }

    ScaleFactor getElectronScaleFactor(double pt, double eta) {
      return getScaleFactorFromMap(pt, eta, mElectronScaleFactors);
    }

    ScaleFactor getMuonScaleFactor(double pt, double eta) const {
      return getScaleFactorFromMap(pt, eta, mMuonScaleFactors);
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
    ScaleFactor getScaleFactorFromMap(double pt, double eta, const ScaleFactorMap& map) const {
      eta = fabs(eta);
      for (const auto& etaBin: map) {
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

    void parseScaleFactors(const edm::ParameterSet& settings, const std::string& name, ScaleFactorMap& map) {
      uint32_t nEta = 0, nPt = 0;
      const edm::VParameterSet& sfs = settings.getParameterSetVector(name);

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
          map[eta][pt] = sf;
        }
      }

      std::cout << "Done. " << nEta << " eta bins and " << nPt << " pt bins loaded." << std::endl;
      std::cout << std::endl;
    }

    void parseMuonScaleFactors(const edm::ParameterSet& settings) {
      std::cout << "Loading muon scale factors..." << std::endl;
      parseScaleFactors(settings, "muon_scale_factors", mMuonScaleFactors);
    }

    void parseElectronScaleFactors(const edm::ParameterSet& settings) {
      std::cout << "Loading electron scale factors..." << std::endl;
      parseScaleFactors(settings, "electron_scale_factors", mElectronScaleFactors);
    }

    ScaleFactorMap mMuonScaleFactors;
    ScaleFactorMap mElectronScaleFactors;

    // BTag
    bool mBTagPerfInit;
    edm::ESHandle<BtagPerformance> mBTagPerf;
    float mBTagEfficiency;
};
