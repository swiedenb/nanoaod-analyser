#ifndef MuonSFTool_h
#define MuonSFTool_h

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

class MuonSFTool {
    
  protected:
    
    std::map<const std::string,const TF1*> func;
    TH2D hist;
    TH2D hist_stat;
    TH2D hist_syst;
    TH2D hist_iso;
    TH2D hist_iso_stat;
    TH2D hist_iso_syst;
    
  public:
    
    std::string ID;
    std::string ISO;
    MuonSFTool(const std::string& year, const std::string& wp="Tight", const std::string& iso="LooseRelTkIso");
    ~MuonSFTool() { }
    
    float getSFID( double pt, double eta, const std::string& unc="") ;
    float getSFISO( double pt, double eta, const std::string& unc="") ;
  
};

#endif // MuonSFTool_h

