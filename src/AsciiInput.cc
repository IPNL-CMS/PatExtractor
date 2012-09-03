/* 
############################################################
#
# AsciiInput.cc
#
############################################################
#
# Author: Seb Viret <viret@in2p3.fr>, inspired from ATLAS BH generator
#         written by W.Bell
#
# May 27th, 2010
#
# Goal: 
# Handling the text files used in specific joboptions
#
#############################################################
*/

#include "../interface/AsciiInput.h"


AsciiInput::AsciiInput(std::string fileName): m_fileName(fileName),
                                               m_file() 
{}
 
//-------------------------------------------------------------------------
 
AsciiInput::~AsciiInput() 
{}
 
//-------------------------------------------------------------------------

int AsciiInput::open() 
{
  m_file.open(m_fileName.c_str());
  if(!m_file) 
  {
    std::cerr << "Error: could not open " << m_fileName << std::endl;
    return 1;
  }  
  return 0;
}
 
int AsciiInput::close() 
{
  if(m_file.is_open()) 
  {
    m_file.close();
  }
  return 0;
}


//-------------------------------------------------------------------------

std::vector<std::string> AsciiInput::readRow() 
{
  std::vector<std::string> row;
  
  // Check if the end of file has been reached.
  if(m_file.eof()) 
  {
    std::cout << "End of file reached."  << std::endl;
    return row;
  }
 
  // Read one line from the input file.  (Need one additional character
  // for string termination of a word.)
  
  m_file.getline(m_lineBuffer,MAX_LINE_LENGTH-1);
  if(m_file.eof()) { // True when normally reaching the end of the file.
    return row;
  }
   
  // Convert the cstring to a string.
  std::string tmpString(m_lineBuffer);
 
  // Split the string into substrings.
  return strToStrVec(tmpString);
}
 
//-------------------------------------------------------------------------

std::vector<std::string> AsciiInput::strToStrVec(std::string inputString) 
{
  std::vector<std::string> strVec;
  std::string tmpString;
  size_t stringLength, i;
  
  stringLength = inputString.length();
  
  // Loop over all characters in the string.
  i=0;
  while(i<stringLength) 
  {
    
    // Skip over any white spaces at the start of the string.
    while(inputString[i] == ' ' || 
          inputString[i] == '\t') {
      i++;
      if(i >= stringLength) break;
    }
    
    if(i>=stringLength) continue;
    
    // Copy all non-white space characters until a white space, end of string or comment character.
    tmpString.clear();
    while(inputString[i] != ' ' &&
	  inputString[i] != '\t' && 
          inputString[i] != '#') 
      {
	tmpString.push_back(inputString[i]);
	i++;
	if(i >= stringLength) break;
      }
    
    // Push the string back if its length is greater than zero.
    if(tmpString.length() > 0) 
    {
      strVec.push_back(tmpString);
    }

    // Check if a comment was found.
    if(i < stringLength) 
    {
      if(inputString[i] == '#') break;
    }
  }
  
  return strVec;
}


//-------------------------------------------------------------------------

long AsciiInput::strToLong(std::string inputString) 
{
  long longValue;
  std::istringstream inStr(inputString);
  inStr >> longValue;
  return longValue;
}
 
//-------------------------------------------------------------------------
 
double AsciiInput::strToDouble(std::string inputString) 
{
  double doubleValue;
  std::istringstream inStr(inputString);
  inStr >> doubleValue;
  return doubleValue;
}
