#pragma once

#include <string>
#include <vector>

// Custom corrections
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

#include <FWCore/ParameterSet/interface/FileInPath.h>

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>

using namespace xercesc;

class XMLSimpleStr {
  public:
    XMLSimpleStr(XMLCh *str) :
      string(XMLString::transcode(str)), ch(str) {}

    XMLSimpleStr(const XMLCh *str) :
      string(XMLString::transcode(str)), ch(NULL) {}

    ~XMLSimpleStr() {
      XMLString::release(&string);
      if (ch) {
        XMLString::release(&ch);
      }
    }

    inline operator const char *() const {
      return string;
    }

    inline std::string get() const {
      return std::string(string);
    }

  private:
    char    *string;
    XMLCh   *ch;
};


inline FactorizedJetCorrector* makeFactorizedJetCorrectorFromXML(const std::string& xmlfile, const std::string& jetAlgo, const bool isMC)
{

  try {
    XMLPlatformUtils::Initialize();
  }
  catch (const XMLException& toCatch) {
    char* message = XMLString::transcode(toCatch.getMessage());
    std::cout << "Error during initialization! :\n" << message << "\n";
    XMLString::release(&message);

    return NULL;
  }

  const std::string mcDataText = (isMC) ? "MC" : "DATA";

  XercesDOMParser parser;
  parser.setValidationScheme(XercesDOMParser::Val_Auto);
  parser.setDoNamespaces(false);
  parser.setDoSchema(false);
  parser.setValidationSchemaFullChecking(false);

  try {
    parser.parse(xmlfile.c_str());

    DOMDocument* xmlDoc = parser.getDocument();
    if (! xmlDoc) {
      throw std::runtime_error("Failed to open '" + xmlfile + "'");
    }
    DOMElement* xmlRoot = xmlDoc->getDocumentElement();

    XMLCh* correctionsStr = XMLString::transcode("corrections");
    DOMElement *corrections = dynamic_cast<DOMElement*>(xmlRoot->getElementsByTagName(correctionsStr)->item(0));
    XMLString::release(&correctionsStr);

    if (! corrections) {
      throw std::runtime_error("Missing corrections tag");
    }

    XMLCh* prefixStr = XMLString::transcode("prefix");
    XMLCh* pathStr = XMLString::transcode("path");

    std::vector<JetCorrectorParameters> correctors;
    XMLSimpleStr prefix(corrections->getAttribute(prefixStr));
    XMLSimpleStr path(corrections->getAttribute(pathStr));

    XMLString::release(&prefixStr);
    XMLString::release(&pathStr);

    DOMNodeList* children = corrections->getChildNodes();
    const XMLSize_t nodeCount = children->getLength();

    XMLCh* nameStr = XMLString::transcode("name");
    XMLCh* onlyOnDataStr = XMLString::transcode("onlyondata");
    XMLCh* trueStr = XMLString::transcode("true");

    for (XMLSize_t i = 0; i < nodeCount; i++) {
      DOMNode * node = children->item(i);

      if (node->getNodeType() && node->getNodeType() == DOMNode::ELEMENT_NODE) {
        DOMElement* element = dynamic_cast<DOMElement*>(node);

        XMLSimpleStr name(element->getAttribute(nameStr));
        bool onlyOnData = false;
        if (element->hasAttribute(onlyOnDataStr)) {
          const XMLCh* foo = element->getAttribute(onlyOnDataStr);
          onlyOnData = XMLString::equals(foo, trueStr);
        }

        std::string filename = path.get() + prefix.get() + "_" + mcDataText + "_" + name.get() + "_" + jetAlgo + ".txt";
        //std::string filename = path.get() + prefix.get() + "_" + name.get() + "_" + jetAlgo + ".txt";
        if (!onlyOnData || (onlyOnData && !isMC)) {
          filename = edm::FileInPath(filename).fullPath();
          std::cout << "Using payload '" << filename << "'" << std::endl;
          correctors.push_back(JetCorrectorParameters(filename));
        }
      }
    }

    XMLString::release(&nameStr);
    XMLString::release(&onlyOnDataStr);
    XMLString::release(&trueStr);

    return new FactorizedJetCorrector(correctors);

  } catch (const XMLException& toCatch) {
    char* message = XMLString::transcode(toCatch.getMessage());
    std::cout << "Exception message is: \n" << message << "\n";
    XMLString::release(&message);
    return NULL;
  } catch (const DOMException& toCatch) {
    char* message = XMLString::transcode(toCatch.msg);
    std::cout << "Exception message is: \n" << message << "\n";
    XMLString::release(&message);
    return NULL;
  }

  XMLPlatformUtils::Terminate();
  return NULL;
}

