#include "config.hh"
#include <fstream>
#include <sstream>
#include <TFile.h>
#include <TROOT.h>
#include "TH2D.h"

namespace config {
	// set default values for all used quantities
	
	// minimum tau pt
    float tau_pt = 80;
    
    // maximum tau eta
    float tau_eta = 2.3;
    
	// minimum muon pt
    float muon_pt = 20;
    
    // maximum tau eta
    float muon_eta = 2.4;
    
    // working points of tau identification
    uint tau_iso_WP = 1;
    uint tau_antiMu_WP = 1;
    uint tau_antiE_WP = 1;
    
    // minimum missing transverse momentum
    float met_pt = 150;
    
    // gen particle cut procedure
    int gen_pdgID = 0;
    std::string cut_type = "";
    float cut_value_min = 0;
    float cut_value_max = 999999;
    
    // pileup hist
    TH1D* pileup_hist = NULL;
    
    // runOnData?
    bool runOnData = false;
    
    // which era?
    int era = 2018;
    
    // this checks systematics later
    std::string run_type = "";
    
    // lets start with systematics
    
    // tau scale factors
    double tau_scale = 1.0;
    double tau_scale_up = 1.0;
    double tau_scale_down = 1.0;
}

bool config::load_config_file(json cfg)
{
	// load all numbers from the config file
	json globalcfg;
    std::ifstream cfg_json;
    cfg_json.open("cfg/analysis_config.cfg");
    cfg_json >> globalcfg;
    
    // always check, if the key is actually in the config file
    if (globalcfg.find("tau_pt") != globalcfg.end()) 	tau_pt = globalcfg["tau_pt"];
    if (globalcfg.find("tau_eta") != globalcfg.end()) 	tau_eta = globalcfg["tau_eta"];
    if (globalcfg.find("muon_pt") != globalcfg.end()) 	muon_pt = globalcfg["muon_pt"];
    if (globalcfg.find("muon_eta") != globalcfg.end()) 	muon_eta = globalcfg["muon_eta"];
    if (globalcfg.find("met_pt") != globalcfg.end()) 	met_pt = globalcfg["met_pt"];
    
    if (globalcfg.find("tau_isolation_WP") != globalcfg.end())
		tau_iso_WP = globalcfg["tau_isolation_WP"];
    if (globalcfg.find("tau_antiE_WP") != globalcfg.end())
		tau_antiE_WP = globalcfg["tau_antiE_WP"];
    if (globalcfg.find("tau_antiMu_WP") != globalcfg.end())
		tau_antiMu_WP = globalcfg["tau_antiMu_WP"];
		
    if (cfg.find("gen_part_pdg_id") != cfg.end()) 	gen_pdgID = cfg["gen_part_pdg_id"];
    if (cfg.find("gen_cut_type") != cfg.end())		cut_type = cfg["gen_cut_type"];
    if (cfg.find("gen_cut_min") != cfg.end())		cut_value_min = cfg["gen_cut_min"];
    if (cfg.find("gen_cut_max") != cfg.end())		cut_value_max = cfg["gen_cut_max"];
    
    
    if (cfg.find("runOnData") != cfg.end())			runOnData = cfg["runOnData"];
    if (cfg.find("era") != cfg.end())				era = cfg["era"];
    
    if (!runOnData) { 
		// read in pileup file
		TFile* pileup_file = new TFile(((std::string) cfg["pileup_file"]).c_str(), "READ");
		pileup_hist = (TH1D*) pileup_file->Get("pileup");
		
		// read in tau scale factor file
		TFile* tau_scale_file = new TFile(((std::string) cfg["tau_scale_factor_file"]).c_str(), "READ");
		TH2D* tau_scale_hist = (TH2D*) tau_scale_file->Get("tau_scale_factor");
		tau_scale = tau_scale_hist->GetBinContent( tau_scale_hist->GetXaxis()->FindBin( tau_iso_WP ),
												  tau_scale_hist->GetYaxis()->FindBin( era ) );
		tau_scale_up = tau_scale_hist->GetBinContent( tau_scale_hist->GetXaxis()->FindBin( tau_iso_WP ), 
													 tau_scale_hist->GetYaxis()->FindBin( era ) ) 
							+ tau_scale_hist->GetBinErrorUp( tau_scale_hist->GetXaxis()->FindBin( tau_iso_WP ),
															tau_scale_hist->GetYaxis()->FindBin( era ) );
		
		tau_scale_down = tau_scale_hist->GetBinContent( tau_scale_hist->GetXaxis()->FindBin( tau_iso_WP ), 
													 tau_scale_hist->GetYaxis()->FindBin( era ) ) 
							- tau_scale_hist->GetBinErrorLow( tau_scale_hist->GetXaxis()->FindBin( tau_iso_WP ),
															tau_scale_hist->GetYaxis()->FindBin( era ) );
		delete tau_scale_hist;
		delete tau_scale_file;    
	}
}
