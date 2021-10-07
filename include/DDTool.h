#ifndef DDTool_h
#define DDTool_h

/*
 * @class TauIDSFTool
 *
 * Class to retrieve tau ID SFs.
 *  - pT-dependent SFs for MVAoldDM2017v2
 *  - DM-dependent SFs for MVAoldDM2017v2
 *  - eta-dependent SFs for anti-lepton discriminators
 * Source: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV
 * Inspiration from TauTriggerSFs/src/TauTriggerSFs2017.cc
 *
 * @author Izaak Neutelings
 * @date July 2019
 *
 */

#include <TFile.h>   // TFile
#include <TH1.h>     // TH1
#include <TH2D.h>     // TH1
#include <TF1.h>     // TF1
#include <TString.h> // Form
#include <string>    // std::string
#include <vector>    // std::vector
#include <map>       // std::map
#include <stdlib.h>  // getenv
#include <functional>

class DDTool {
    
  protected:
    
    std::map<const std::string,const TF1*> func;
    TH2F hist;
    
  public:
    
    DDTool(const std::string& year);
    ~DDTool() { }
    
    float getFF( float taupt, float taupt_o_jetpt, const std::string& year, const std::string& unc="") ;
  
};

#endif // DDTool_h

