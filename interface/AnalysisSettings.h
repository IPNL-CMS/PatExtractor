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

  bool checkSetting(std::string key, double value);
  float getSetting(std::string key);

 private:


  
  
  // A vector of strings to configure the analysis settings.

  std::vector<std::string>* m_analysisSettings;  
  std::map<std::string, float > m_settings;
 

  // Flag
  bool m_settingsParsed;
};
 
#endif
