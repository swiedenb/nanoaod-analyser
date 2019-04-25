#include "config.hh"
#include <fstream>
#include <sstream>
#include <TFile.h>
#include <TROOT.h>

namespace config {
	// set default values for all used quantities
	
	// minimum tau pt
    float tau_pt = 80;
    
    // maximum tau eta
    float tau_eta = 2.3;
    
    // working points of tau identification
    char tau_iso_WP[10] = "2";
    char tau_antiMu_WP[10] = "1";
    char tau_antiE_WP[10] = "1";
    
    // minimum missing transverse momentum
    float met_pt = 150;
    
    // gen particle cut procedure
    int gen_pdgID = 0;
    std::string cut_type = "";
    float cut_value_min = 0;
    float cut_value_max = 999999;
    
    // pileup hist
    TH1D* pileup_hist = NULL;
}

bool config::load_config_file(json cfg)
{
	// load all numbers from the config file
	json globalcfg;
    std::ifstream cfg_json;
    cfg_json.open("cfg/analysis_config.cfg");
    cfg_json >> globalcfg;
    
    
    tau_pt = globalcfg["tau_pt"];
    tau_eta = globalcfg["tau_eta"];
    met_pt = globalcfg["met_pt"];
    
    std::string tmp1 = globalcfg["tau_isolation_WP"];
    std::string tmp2 = globalcfg["tau_antiE_WP"];
    std::string tmp3 = globalcfg["tau_antiMu_WP"];
    strcpy(tau_iso_WP, tmp1.c_str());
    strcpy(tau_antiE_WP, tmp2.c_str());
    strcpy(tau_antiMu_WP, tmp3.c_str());
    
    gen_pdgID = cfg["gen_part_pdg_id"];
    cut_type = cfg["gen_cut_type"];
    cut_value_min = cfg["gen_cut_min"];
    cut_value_max = cfg["gen_cut_max"];
    
    TFile* pileup_file = new TFile(((std::string) cfg["pileup_file"]).c_str(), "READ");
    pileup_hist = (TH1D*) pileup_file->Get("pileup");
    
}
