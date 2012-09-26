#include "../interface/dimuon_analysis.h"

dimuon_analysis::dimuon_analysis(AnalysisSettings *settings)
{
  std::cout << "Entering dimuon analysis" << std::endl;

  /// Tree definition containing your analysis results

  m_tree_dimuon = new TTree("dimuon","dimuon Analysis info");  
  
  /// Branches definition

  m_tree_dimuon->Branch("evt",         &m_evt    ,"evt/I");      // Simple evt number or event ID
  m_tree_dimuon->Branch("dimuon_mass", &m_mdimuon,"dimuon_mass/F");
  m_tree_dimuon->Branch("dimuon_mom",  &m_totp   ,"dimuon_mom/F");
  m_tree_dimuon->Branch("mu1_pt",      &m_mu1pt  ,"mu1_pt/F");
  m_tree_dimuon->Branch("mu2_pt",      &m_mu2pt  ,"mu2_pt/F");


  /// Analysis settings (you define them in your python script)

  // If the setting is not defined the method returns -1. In this 
  // case we set a default value for the cut, in order to avoid 
  // unwanted crash

  settings->getSetting("pTmu_min", m_mu_Ptcut);
//    ? m_mu_Ptcut = settings->getSetting("pTmu_min") // Value from the joboption
//    : m_mu_Ptcut = 0.5;                               // Default val


  settings->getSetting("nmu_min", m_mu_Mult);
//    ? m_mu_Mult = settings->getSetting("nmu_min") // Value from the joboption
//    : m_mu_Mult = 2;                                // Default val

}

dimuon_analysis::~dimuon_analysis(){;}

void dimuon_analysis::reset()
{
  m_evt = 0;

  m_mdimuon = 0.;
  m_totp    = 0.;
  m_mu1pt   = 0.;
  m_mu2pt   = 0.;
}



int dimuon_analysis::dimuon_Sel(MuonExtractor *muon, int evtnum) 
{
  int n_mus = muon->getSize();


  // First check the number of muons

  if (n_mus < m_mu_Mult) return 0; 


  // Then loop on them and compute their mass

  for (int i=0;i<n_mus-1;++i)
  {
    for (int j=i;j<n_mus;++j)
    {



      TLorentzVector *mui = muon->getMuLorentzVector(i);
      TLorentzVector *muj = muon->getMuLorentzVector(j);

      if (!dimuon_analysis::isMuSel(mui)) continue; // apply the pt cut
      if (!dimuon_analysis::isMuSel(muj)) continue;


      double E  = mui->E()+muj->E();
      double px = mui->Px()+muj->Px();    
      double py = mui->Py()+muj->Py();    
      double pz = mui->Pz()+muj->Pz();    
      
      m_evt     = evtnum;
      m_mdimuon = sqrt(E*E-px*px-py*py-pz*pz);
      m_totp    = sqrt(px*px+py*py+pz*pz);
      m_mu1pt   = mui->Pt();
      m_mu2pt   = muj->Pt();

      // Finally fill the analysis tree 
      dimuon_analysis::fillTree();
    }
  }

  return 0;
}


// Here we do the muon selection, using the cuts defined in the settings
// In this example this is a simple pt cut, but you can add whatever you want


bool dimuon_analysis::isMuSel(TLorentzVector *muon)
{  
  double mupt = muon->Pt();

  if (mupt<m_mu_Ptcut)            return false;

  return true;
}


// Fill the root tree containing analysis results

void dimuon_analysis::fillTree()
{
  m_tree_dimuon->Fill(); 
}


