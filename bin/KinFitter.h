#pragma once

#include <vector>
#include <boost/shared_ptr.hpp>

#include "TopQuarkAnalysis/TopKinFitter/interface/TtSemiLepKinFitter.h"

namespace edm {
  class ParameterSet;
}

class TLorentzVector;

class KinFitter {
  public:
    KinFitter(const edm::ParameterSet& cfg);

    TtSemiLepKinFitter::Param param(unsigned val) const;
    TtSemiLepKinFitter::Constraint constraint(unsigned val) const;
    std::vector<TtSemiLepKinFitter::Constraint> constraints(std::vector<unsigned>& val) const;

    TtSemiLepKinFitter& getFitter();

    bool PzNeutrino(const TLorentzVector& lepton, TLorentzVector& neutrino, const TLorentzVector& bJet) const;

  private:
    /// maximal number of iterations to be performed for the fit
    unsigned int maxNrIter_;
    /// maximal chi2 equivalent
    double maxDeltaS_;
    /// maximal deviation for contstraints
    double maxF_;
    unsigned int jetParam_;
    unsigned int lepParam_;
    unsigned int metParam_;
    /// constrains
    std::vector<unsigned> constraints_;
    double mW_;
    double mTop_;
    /// scale factors for jet energy resolution
    std::vector<double> jetEnergyResolutionScaleFactors_;
    std::vector<double> jetEnergyResolutionEtaBinning_;
    /// config-file-based object resolutions
    std::vector<edm::ParameterSet> udscResolutions_;
    std::vector<edm::ParameterSet> bResolutions_;
    std::vector<edm::ParameterSet> lepResolutions_;
    std::vector<edm::ParameterSet> metResolutions_;

    boost::shared_ptr<TtSemiLepKinFitter> m_fitter;
};
