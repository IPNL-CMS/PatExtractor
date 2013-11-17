#pragma once

#include <memory>

#include "Extractors/PatExtractor/interface/ScaleFactorService.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class BTaggingEfficiencyProvider {
  public:
    BTaggingEfficiencyProvider(const edm::ParameterSet& settings) {
      TH1::AddDirectory(false);

      const edm::ParameterSet& pset = settings.getParameterSet("b_tagging_efficiency");
      std::string filename = edm::FileInPath(pset.getParameter<std::string>("filename")).fullPath();

      TFile* f = TFile::Open(filename.c_str());

      // Get histo from file
      m_b_tagging_efficiency.reset(static_cast<TH2*>(f->Get(pset.getParameter<std::string>("b_eff_histo_name").c_str())->Clone()));
      m_cjets_fakerate.reset(static_cast<TH2*>(f->Get(pset.getParameter<std::string>("cjets_fakerate_histo_name").c_str())->Clone()));
      m_lightjets_fakerate.reset(static_cast<TH2*>(f->Get(pset.getParameter<std::string>("lightjets_fakerate_histo_name").c_str())->Clone()));
      
      f->Close();
      delete f;
    }

    float getEfficiency(ScaleFactorService::Flavor flavor, float pt, float eta) {
      if (pt > 800)
        pt = 800;

      eta = fabs(eta);
      TH2* histo = nullptr;

      switch (flavor) {
        case ScaleFactorService::B:
          histo = m_b_tagging_efficiency.get();
          break;

        case ScaleFactorService::C:
          histo = m_cjets_fakerate.get();
          break;

        case ScaleFactorService::LIGHT:
          histo = m_lightjets_fakerate.get();
          break;
      }

      if (! histo)
        return 0.;

      return histo->GetBinContent(histo->FindBin(pt, eta));
    }

    float getEfficiencyError(ScaleFactorService::Flavor flavor, float pt, float eta) {
      if (pt > 800)
        pt = 800;

      eta = fabs(eta);
      TH2* histo = nullptr;

      switch (flavor) {
        case ScaleFactorService::B:
          histo = m_b_tagging_efficiency.get();
          break;

        case ScaleFactorService::C:
          histo = m_cjets_fakerate.get();
          break;

        case ScaleFactorService::LIGHT:
          histo = m_lightjets_fakerate.get();
          break;
      }

      if (! histo)
        return 0.;

      return histo->GetBinError(histo->FindBin(pt, eta));
    }

  private:
    std::shared_ptr<TH2> m_b_tagging_efficiency;
    std::shared_ptr<TH2> m_cjets_fakerate;
    std::shared_ptr<TH2> m_lightjets_fakerate;
};
