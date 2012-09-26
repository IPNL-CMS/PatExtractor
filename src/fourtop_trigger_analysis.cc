#include "../interface/fourtop_trigger_analysis.h"

fourtop_trigger_analysis::fourtop_trigger_analysis(AnalysisSettings *settings)
{
  std::cout << "Entering fourtop_trigger analysis" << std::endl;

  /// Tree definition containing your analysis results

  //  m_tree_fourtop_trigger = new TTree("fourtop_trigger","fourtop_trigger Analysis info");  
  
  /// Branches definition

  //  m_tree_fourtop_trigger->Branch("evt",         &m_evt    ,"evt/I");      // Simple evt number or event ID
  //  m_tree_fourtop_trigger->Branch("fourtop_trigger_mass", &m_mfourtop_trigger,"fourtop_trigger_mass/F");

  /// Analysis settings (you define them in your python script)

  // If the setting is not defined the method returns -1. In this 
  // case we set a default value for the cut, in order to avoid 
  // unwanted crash

  settings->getSetting("n_SSlept", m_SSlept);
//    ? m_SSlept = settings->getSetting("n_SSlept") // Value from the joboption
//    : m_SSlept = 0;                               // Default val

  fourtop_trigger_analysis::reset();
}

fourtop_trigger_analysis::~fourtop_trigger_analysis(){;}

void fourtop_trigger_analysis::reset()
{
  m_ndiff_triggers=0;
  m_HLT_vector_tot.clear();
  n_tot_evt = 0;

  for (int i=0;i<500;++i) m_tot[i]=0;   
}



int fourtop_trigger_analysis::fourtop_trigger_Sel(HLTExtractor *hlt, MCExtractor *mc, int evtnum) 
{
  // First check the number of SS lepton in the original events

  if (fourtop_trigger_analysis::nSSlept(mc) < m_SSlept) return 0;

  ++n_tot_evt;

  int npaths = hlt->getSize();

  for (int i=0;i<npaths;++i)
  {
    // Is this path already known of not?
    m_iter = m_HLT_vector_tot.find(hlt->paths(i));
    int indext = -1;

    if (m_iter == m_HLT_vector_tot.end()) // New path, add it...
    {
      m_HLT_vector_tot.insert(std::make_pair(hlt->paths(i),m_ndiff_triggers));
      indext = m_ndiff_triggers;
      ++m_ndiff_triggers;
    }
    else
    {
      indext = m_iter->second;
    }

    ++m_tot[indext];         
  }

  return 0;
}


// Here we do the muon selection, using the cuts defined in the settings
// In this example this is a simple pt cut, but you can add whatever you want


int fourtop_trigger_analysis::nSSlept(MCExtractor *mc)
{  
  int n_part =  mc->getSize();

  int n_m_lept=0;
  int n_p_lept=0;

  for (int i=0;i<n_part;++i)
  {
    if (mc->getType(i)==11 || mc->getType(i)==13) ++n_p_lept;
    if (mc->getType(i)==-11 || mc->getType(i)==-13) ++n_m_lept;
  }

  //  std::cout << n_p_lept << " / " << n_m_lept << " / " << std::max(n_p_lept,n_m_lept) << std::endl;

  return std::max(n_p_lept,n_m_lept);
}


// Fill the root tree containing analysis results

void fourtop_trigger_analysis::fillTree()
{
  //  m_tree_fourtop_trigger->Fill(); 
}


void fourtop_trigger_analysis::fourtop_trigger_finalize(int evtnum)
{
  for (std::multimap< std::string, int>::const_iterator ii=m_HLT_vector_tot.begin(); ii!=m_HLT_vector_tot.end(); ++ii)
    print_results((*ii).second,(*ii).first,evtnum);
}

void fourtop_trigger_analysis::print_results(int index,std::string pathname,int evtnum)
{
  std::cout << pathname << " / " << m_tot[index] << " / " << n_tot_evt << " / = / " 
	    << 100.*double(m_tot[index])/double(n_tot_evt)  << std::endl; 
}
