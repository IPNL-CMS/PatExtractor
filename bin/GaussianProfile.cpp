#include "GaussianProfile.h"

#include <sstream>
#include <TH1D.h>
#include <TF1.h>


void GaussianProfile::createProfiles() {
  for (int i = 0; i < m_nXBins; i++) {

    std::stringstream ss;
    ss << m_name << "_" << m_prefix << "_" << (int) m_XBins[i] << "_" << (int) m_XBins[i + 1];

    int nBins = m_nYBins;
    double min = m_YMin, max = m_YMax;

    if (m_autoBinning) {
      nBins = 100;
      min = m_XBins[i] - m_autoBinningLowPercent * m_XBins[i];
      max = m_XBins[i + 1] + m_autoBinningHighPercent * m_XBins[i + 1];
    }

    TH1* object = new TH1D(ss.str().c_str(), ss.str().c_str(), nBins, min, max);
    m_profiles.push_back(object);
  }
}

void GaussianProfile::createGraph() {

  if (m_profiles.size() == 0 || (m_graph.get() && !m_dirty))
    return;

  std::stringstream ss;
  ss << m_name << "_graph";

  m_graph.reset(new TGraphErrors(m_nXBins));
  m_graph->SetName(ss.str().c_str());
  
  // Create gaussian for fitting
  TF1* gauss = new TF1("g", "gaus");

  for (int i = 0; i < m_nXBins; i++) {
    TH1* hist = m_profiles[i];
    double mean = (m_XBins[i] + m_XBins[i + 1]) / 2.;
    //double min = hist->GetXaxis()->GetBinLowEdge(1);
    //double max = hist->GetXaxis()->GetBinUpEdge(hist->GetXaxis()->GetLast());
    
    double min = hist->GetXaxis()->GetXmin();
    double max = hist->GetXaxis()->GetXmax();
    if (m_autoBinning) {
      min = m_XBins[i] * 0.90;
      max = m_XBins[i + 1] * 1.10;
    }

    gauss->SetRange(min, max);
    gauss->SetParameter(1, (min + max) / 2.);
    gauss->SetParLimits(1, min, max);
    gauss->SetParameter(2, 20);

    hist->Fit(gauss, "QR");

    m_graph->SetPoint(i, mean, gauss->GetParameter(1));
    m_graph->SetPointError(i, 0, gauss->GetParameter(2));
  }

  m_dirty = false;
}
