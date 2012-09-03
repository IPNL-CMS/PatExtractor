/* 
############################################################
#
# AnalysisSettings.cc
#
############################################################
#
# Author: Seb Viret <viret@in2p3.fr>, inspired from ATLAS BH generator
#         written by W.Bell
#
# May 27th, 2010
#
# Goal: 
# Get a set of cuts used in any analysis from python joboption
#
# Each parameter is defined as a string and a cut value
#
#############################################################
*/

#include "../interface/AnalysisSettings.h"
 
AnalysisSettings::AnalysisSettings(std::vector<std::string> *settings): m_analysisSettings(settings),	       		 
									m_settingsParsed(false) 
{}
 
//---------------------------------------------------------------------
//
// The settings appear as a text line in the job options, one has to 
// parse and decode them. This is the role of this method
//
//---------------------------------------------------------------------

int AnalysisSettings::parseSettings() 
{ 
  std::vector<std::string> strVector;
  std::vector<std::string>::iterator row_itr;
  std::vector<std::string>::iterator row_itr_end;
  std::vector<std::string>::iterator col_itr;
  std::vector<std::string>::iterator col_itr_end;
  
  // Loop over all settings.
  row_itr     = m_analysisSettings->begin();
  row_itr_end = m_analysisSettings->end();
 
  for(;row_itr!=row_itr_end;++row_itr) 
  {
    strVector.clear();
    strVector = AsciiInput::strToStrVec(*row_itr);
    
    if(strVector.size() == 0) continue;
    if(strVector.size() == 1) 
    {
      std::cerr << "A generator setting must be followed by a value" << std::endl;
      continue;
    }
 
       
    // The basic method: lower and upper limit for the considered value 
    AnalysisSettings::parseLimitSetting(&strVector); 
      
  }
 
  AnalysisSettings::printSettings();  
  m_settingsParsed = true;
 
  return 0;
}
 
//---------------------------------------------------------------------
//
// Method setting upper and lower limit for a given parameter contained in
// a pre-defined list
//
//---------------------------------------------------------------------

int AnalysisSettings::parseLimitSetting(std::vector<std::string> *commandVector) 
{
  double limit = -1.0;
    
  std::vector<std::string>::iterator itr     = commandVector->begin();
  std::vector<std::string>::iterator itr_end = commandVector->end();
 
  std::string limit_name = (*itr);
  itr++;

  if(itr!=itr_end) limit = AsciiInput::strToDouble(*itr);

  m_settings.insert(std::make_pair(limit_name,limit));
  
  return 0;
}
 
  

 
//---------------------------------------------------------------------

bool AnalysisSettings::checkSetting(std::string key, double value) 
{
  std::map<std::string,float>::iterator itr = m_settings.find(key);
  if(itr == m_settings.end()) 
  {
    std::cerr << "Error: the limit " << key << " is not defined." << std::endl;
    return false;
  }
  
  if((*itr).second >= 0. && value < (*itr).second) return false;
  
  return true;
}

float AnalysisSettings::getSetting(std::string key)
{
  std::map<std::string,float>::iterator itr = m_settings.find(key);
  if(itr == m_settings.end()) 
  {
    std::cerr << "Error: the limit " << key << " is not defined." << std::endl;
    return -1;
  }

  return (*itr).second;
}

//---------------------------------------------------------------------

void AnalysisSettings::printSettings() 
{
  std::cout << "##################################################" << std::endl;

  std::cout << "Description of analysis cuts:" << std::endl;


 
  std::map<std::string, float >::iterator limit_itr = m_settings.begin();
  std::map<std::string, float >::iterator limit_itr_end = m_settings.end();

  for(;limit_itr!=limit_itr_end;++limit_itr) 
  {
    std::cout << " " << (*limit_itr).first << ": ";
    
    if((*limit_itr).second <= 0.) std::cout << "disabled, ";
    else std::cout << (*limit_itr).second;

    std::cout << std::endl;
  }
}
