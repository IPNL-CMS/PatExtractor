#ifndef ANALYSISSETTINGS_H
#define ANALYSISSETTINGS_H

#include "../interface/AsciiInput.h"

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <cmath>

class AnalysisSettings 
{
 public:
 
  AnalysisSettings(std::vector<std::string> *settings);
 
  // Returns 0 if successful and a status code otherwise.
  int parseSettings(void);

  int parseLimitSetting(std::vector<std::string> *commandVector); 

  void printSettings(void);

  template<class T> bool getSetting(std::string key, T& value) {

    std::map<std::string, float>::iterator itr = m_numericSettings.find(key);
    if(itr == m_numericSettings.end())
    {
      std::cerr << "Error: the limit " << key << " is not defined." << std::endl;
      return false;
    }

    value = itr->second;
    return true;
  }


 private:


  
  
  // A vector of strings to configure the analysis settings.

  std::vector<std::string>* m_analysisSettings;  

  std::map<std::string, float> m_numericSettings;
  std::map<std::string, std::string> m_stringSettings;
 

  // Flag
  bool m_settingsParsed;
};

template<> inline bool AnalysisSettings::getSetting<std::string>(std::string key, std::string& value) {

  std::map<std::string, std::string>::iterator itr = m_stringSettings.find(key);
  if(itr == m_stringSettings.end()) 
  {
    std::cerr << "Error: the limit " << key << " is not defined." << std::endl;
    return false;
  }

  value = itr->second;
  return true;
}

#endif
