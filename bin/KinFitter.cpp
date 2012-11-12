#include "KinFitter.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

KinFitter::KinFitter(const edm::ParameterSet& cfg):
  maxNrIter_               (cfg.getParameter<unsigned>     ("maxNrIter"           )),
  maxDeltaS_               (cfg.getParameter<double>       ("maxDeltaS"           )),
  maxF_                    (cfg.getParameter<double>       ("maxF"                )),
  jetParam_                (cfg.getParameter<unsigned>     ("jetParametrisation"  )),
  lepParam_                (cfg.getParameter<unsigned>     ("lepParametrisation"  )),
  metParam_                (cfg.getParameter<unsigned>     ("metParametrisation"  )),
  constraints_             (cfg.getParameter<std::vector<unsigned> >("constraints")),
  mW_                      (cfg.getParameter<double>       ("mW"                  )),
  mTop_                    (cfg.getParameter<double>       ("mTop"                )),
  jetEnergyResolutionScaleFactors_(cfg.getParameter<std::vector<double> >("jetEnergyResolutionScaleFactors")),
  jetEnergyResolutionEtaBinning_  (cfg.getParameter<std::vector<double> >("jetEnergyResolutionEtaBinning")),
  udscResolutions_(0), bResolutions_(0), lepResolutions_(0), metResolutions_(0)
{

  if (cfg.exists("udscResolutions") && cfg.exists("bResolutions") && cfg.exists("lepResolutions") && cfg.exists("metResolutions")) {
    udscResolutions_ = cfg.getParameter<std::vector<edm::ParameterSet> >("udscResolutions");
    bResolutions_    = cfg.getParameter<std::vector<edm::ParameterSet> >("bResolutions"   );
    lepResolutions_  = cfg.getParameter<std::vector<edm::ParameterSet> >("lepResolutions" );
    metResolutions_  = cfg.getParameter<std::vector<edm::ParameterSet> >("metResolutions" );
  }
  else if (cfg.exists("udscResolutions") || cfg.exists("bResolutions") || cfg.exists("lepResolutions") || cfg.exists("metResolutions")) {
    throw cms::Exception("Configuration") << "Parameters 'udscResolutions', 'bResolutions', 'lepResolutions', 'metResolutions' should be used together.\n";
  }

  m_fitter.reset(new TtSemiLepKinFitter(param(jetParam_), param(lepParam_), param(metParam_), maxNrIter_, maxDeltaS_, maxF_,
        constraints(constraints_), mW_, mTop_, &udscResolutions_, &bResolutions_, &lepResolutions_, &metResolutions_,
        &jetEnergyResolutionScaleFactors_, &jetEnergyResolutionEtaBinning_));
}

TtSemiLepKinFitter::Param KinFitter::param(unsigned val) const {
  TtSemiLepKinFitter::Param result;
  switch (val) {
    case TtSemiLepKinFitter::kEMom:
      result = TtSemiLepKinFitter::kEMom;
      break;
    case TtSemiLepKinFitter::kEtEtaPhi:
      result = TtSemiLepKinFitter::kEtEtaPhi;
      break;
    case TtSemiLepKinFitter::kEtThetaPhi:
      result = TtSemiLepKinFitter::kEtThetaPhi;
      break;
    default:
      throw cms::Exception("Configuration")  << "Chosen jet parametrization is not supported: " << val << std::endl;
      break;
  }

  return result;
}

TtSemiLepKinFitter::Constraint KinFitter::constraint(unsigned val) const {
  TtSemiLepKinFitter::Constraint result;
  switch (val) {
    case TtSemiLepKinFitter::kWHadMass:
      result = TtSemiLepKinFitter::kWHadMass;
      break;
    case TtSemiLepKinFitter::kWLepMass:
      result = TtSemiLepKinFitter::kWLepMass;
      break;
    case TtSemiLepKinFitter::kTopHadMass:
      result = TtSemiLepKinFitter::kTopHadMass;
      break;
    case TtSemiLepKinFitter::kTopLepMass:
      result = TtSemiLepKinFitter::kTopLepMass;
      break;
    case TtSemiLepKinFitter::kNeutrinoMass:
      result = TtSemiLepKinFitter::kNeutrinoMass;
      break;
    case TtSemiLepKinFitter::kEqualTopMasses:
      result = TtSemiLepKinFitter::kEqualTopMasses;
      break;
    case TtSemiLepKinFitter::kSumPt:
      result = TtSemiLepKinFitter::kSumPt;
      break;

    default:
      throw cms::Exception("Configuration")  << "Chosen fit constraint is not supported: " << val << std::endl;
      break;
  }

  return result;
}

std::vector<TtSemiLepKinFitter::Constraint> KinFitter::constraints(std::vector<unsigned>& val) const
{
  std::vector<TtSemiLepKinFitter::Constraint> result;
  for (unsigned i = 0; i < val.size(); ++i) {
    result.push_back(constraint(val[i]));
  }

  return result; 
}

TtSemiLepKinFitter& KinFitter::getFitter() {
  return *m_fitter;
}

/**
 * Try to compute neutrino Pz using W mass constraint
 * @return True if the neutrino is corrected, false otherwise
 */
bool KinFitter::PzNeutrino(const TLorentzVector& lepton, TLorentzVector& neutrino, const TLorentzVector& bJet) const
{
  if (!lepton.E())
    return 0;

  double x = (mW_ * mW_ - lepton.M() * lepton.M() + 2. * (neutrino.Px() * lepton.Px() + neutrino.Py() * lepton.Py())) / (2 * lepton.E());
  double a = 1 - (lepton.Pz() * lepton.Pz()) / (lepton.E() * lepton.E());
  double b = -2. * (lepton.Pz() / lepton.E()) * x;
  double c = neutrino.Pt() * neutrino.Pt() - x * x;

  if (!a && !b) return 0;

  if (!a)
  {     
    neutrino.SetPz(-1 * c / b);
    neutrino.SetE(sqrt(neutrino.Px() * neutrino.Px() + neutrino.Py() * neutrino.Py() + neutrino.Pz() * neutrino.Pz()));
    return 1;
  }


  double delta = b * b - 4 * a *c;


  if (delta < 0) // No solution, try to correct MET
  {
    double rat = neutrino.Py() / neutrino.Px();

    double u = 4. / (lepton.E() * lepton.E()) * ((lepton.Px() + rat * lepton.Py()) * (lepton.Px() + rat * lepton.Py()) / (1 + rat * rat)
        - (lepton.E() * lepton.E()) + (lepton.Pz() * lepton.Pz()));

    double v = 4. / (lepton.E() * lepton.E()) * (mW_ * mW_ - lepton.M() *lepton.M())
      * (lepton.Px() + rat * lepton.Py()) / sqrt(1 + rat * rat);

    double w = (mW_ * mW_ - lepton.M() * lepton.M()) * (mW_ * mW_ - lepton.M() *lepton.M()) / (lepton.E() * lepton.E());

    double deltan = v * v - 4 * u * w;

    if (deltan < 0)
      return false; // Hopeless, MET can't be corrected

    double pt      = 0.;
    double corfact = 0.;

    if (u == 0)
    {
      pt = -w/v;
      if (pt <= 0)
        return false; // There is no way out...

      corfact = pt / neutrino.Pt();
    }
    else // Deltan>=0 and u!=0
    {
      double pt2 = (v - (sqrt(deltan))) / (2 * u);
      double pt1 = (v + (sqrt(deltan))) / (2 * u);

      // Pas de correction car negative
      if (pt1 <= 0 && pt2 <= 0)
        return false;

      if (pt1 > 0 && pt2 < 0)
        pt = pt1;

      if (pt2 > 0 && pt1 < 0)
        pt = pt2;

      if (pt1 > 0 && pt2 > 0)
      {
        (fabs(pt1 - neutrino.Pt()) <= fabs(pt2 - neutrino.Pt()))
          ? pt = pt1
          : pt = pt2;     
      }

      corfact = pt / neutrino.Pt();
    }

    // Now we have the correction factor

    neutrino.SetPx(corfact*neutrino.Px());
    neutrino.SetPy(corfact*neutrino.Py());

    // Recompute the new parameters

    x = (mW_ * mW_ - lepton.M() * lepton.M() + 2. * (neutrino.Px() * lepton.Px() + neutrino.Py() * lepton.Py())) / (2 * lepton.E());
    a = 1 - (lepton.Pz() * lepton.Pz()) / (lepton.E() * lepton.E());
    b = -2. * (lepton.Pz() / lepton.E()) * x;
    c = neutrino.Px() * neutrino.Px() + neutrino.Py() * neutrino.Py() - x * x;

//         std::cout << "We have rescaled the MET " << lepton->E() << " / " << corfact <<  " , now delta should be null:" << std::endl;
//         std::cout << "Previous delta: " << delta<< std::endl;


    delta = b * b - 4 * a * c;

    if (fabs(delta) < 0.000001)
      delta = 0.;

//    std::cout << "New delta     : " << delta << std::endl;

    if (delta != 0)
      return false; // This should not happen, but who knows...
  }


  // We can go back to the normal path: 

  TLorentzVector TopCand1 = lepton + bJet;
  TLorentzVector TopCand2 = lepton + bJet;

  neutrino.SetPz((-b - (sqrt(delta))) / (2 * a));
  neutrino.SetE(sqrt(neutrino.Px() * neutrino.Px() + neutrino.Py() * neutrino.Py() + neutrino.Pz() * neutrino.Pz()));
  TopCand1 += neutrino;

  neutrino.SetPz((-b + (sqrt(delta))) / (2 * a));
  neutrino.SetE(sqrt(neutrino.Px() * neutrino.Px() + neutrino.Py() * neutrino.Py() + neutrino.Pz() * neutrino.Pz()));
  TopCand2 += neutrino;

  double mtt_1 = sqrt(std::max(0., TopCand1.M2()));
  double mtt_2 = sqrt(std::max(0., TopCand2.M2()));

  if(fabs(mtt_1 - mTop_) <= fabs(mtt_2 - mTop_)) // Otherwise it's OK
  {
    neutrino.SetPz((-b - (sqrt(delta))) / (2 * a));
    neutrino.SetE(sqrt(neutrino.Px() * neutrino.Px() + neutrino.Py() * neutrino.Py() + neutrino.Pz() * neutrino.Pz()));
  }

  return true;
} 

