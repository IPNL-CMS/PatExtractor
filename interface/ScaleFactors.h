#pragma once

#include <map>
#include <vector>

#include <FWCore/ParameterSet/interface/ParameterSet.h>

struct ScaleFactor {
  double value;
  double error_low;
  double error_high;
};

class ScaleFactors {
  public:
    ScaleFactors(const edm::ParameterSet& settings) {
      parseMuonScaleFactors(settings);
    }

    ScaleFactor getElectronScaleFactor(double pt, double eta) {
      return {1., 0, 0};
    }

    ScaleFactor getMuonScaleFactor(double pt, double eta) {

      eta = fabs(eta);
      for (const auto& etaBin: mMuonScaleFactors) {
        if (eta >= etaBin.first.first && eta < etaBin.first.second) {
          // Look for pt bin
          for (const auto& ptBin: etaBin.second) {
            if (pt >= ptBin.first.first && pt < ptBin.first.second) {
              return ptBin.second;
            }
          }
        }
      }

      return {1., 0, 0};
    }

    ScaleFactor getBTaggingScaleFactor(double pt, double eta) {
      return {1., 0, 0};
    }

    ScaleFactor getHLTScaleFactor(double pt, double eta) {
      return {1., 0, 0};
    }


  private:
    void parseMuonScaleFactors(const edm::ParameterSet& settings) {
      std::cout << "Loading muon scale factors..." << std::endl;
      uint32_t nEta = 0, nPt = 0;
      const edm::VParameterSet& sfs = settings.getParameterSetVector("muon_scale_factors");

      for (const edm::ParameterSet& etaBin: sfs) {
        nEta++;
        nPt = 0;
        const std::vector<double> etaBinVector = etaBin.getParameter<std::vector<double>>("eta");
        auto eta = std::make_pair(etaBinVector[0], etaBinVector[1]);

        for (const edm::ParameterSet& ptBin: etaBin.getParameterSetVector("SF")) {
          nPt++;
          const std::vector<double> ptBinVector = ptBin.getParameter<std::vector<double>>("pt");
          auto pt = std::make_pair(ptBinVector[0], ptBinVector[1]);

          ScaleFactor sf = {
            ptBin.getParameter<double>("value"),
            ptBin.getParameter<double>("error_low"),
            ptBin.getParameter<double>("error_high")
          };

          mMuonScaleFactors[eta][pt] = sf;
        }
      }

      std::cout << "Done. " << nEta << " eta bins and " << nPt << " pt bins loaded." << std::endl;
      std::cout << std::endl;
    }


    std::map<
      // Eta binning
      std::pair<double, double>,
      std::map<
        // Pt binning
        std::pair<double, double>,
        ScaleFactor
      >
    > mMuonScaleFactors;
};
