#include "config.hh"
#include <fstream>
#include <sstream>
#include <TFile.h>
#include <TROOT.h>
#include "TH2D.h"

namespace config {
	// set default values for all used quantities
	
	// trigger string
	std::string trigger = "HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr || HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1";
	
	// minimum tau pt
    float tau_pt = 80;
    
    // maximum tau eta
    float tau_eta = 2.3;
    
	// minimum muon pt
    float muon_pt = 20;
    
    // maximum tau eta
    float muon_eta = 2.4;
    
    // working points of tau identification
    std::string tau_dm = "Tau_idDecayModeNewDMs";
    std::string tau_iso = "Tau_idMVAnewDM2017v2";
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
    
    // W kfactor hist
    TH1D* W_kfactor_hist = NULL;
    
    // tau electron fake scale factor hist
    TH1D* tau_ele_fake_hist = NULL;
    
    // tau muon fake scale factor hist
    TH1D* tau_muo_fake_hist = NULL;
    
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
    
    // tau energy scale
    std::map< std::string, std::vector< double > > tau_energy_scale = {
							{"h+-", 			{1.0, 1.0, 1.0} },
							{"h+- pi0s", 		{1.0, 1.0, 1.0} },
							{"h+-h+-h+-", 		{1.0, 1.0, 1.0} },
							{"h+-h+-h+- pi0s", 	{1.0, 1.0, 1.0} }
						};
	
}

bool config::load_config_file(json cfg)
{
	// load all numbers from the config file
	json globalcfg;
    std::ifstream cfg_json;
    cfg_json.open("cfg/analysis_config.cfg");
    cfg_json >> globalcfg;
    
    // always check, if the key is actually in the config file
    if (globalcfg.find("tau_pt") != globalcfg.end()) 			tau_pt = globalcfg["tau_pt"];
    if (globalcfg.find("tau_eta") != globalcfg.end()) 			tau_eta = globalcfg["tau_eta"];
    if (globalcfg.find("muon_pt") != globalcfg.end()) 			muon_pt = globalcfg["muon_pt"];
    if (globalcfg.find("muon_eta") != globalcfg.end()) 			muon_eta = globalcfg["muon_eta"];
    if (globalcfg.find("met_pt") != globalcfg.end()) 			met_pt = globalcfg["met_pt"];
    if (globalcfg.find("tau_decayMode") != globalcfg.end()) 	tau_dm = globalcfg["tau_decayMode"];
    if (tau_dm == "Tau_idDecayModeNewDMs") {
		tau_iso = "Tau_idMVAnewDM2017v2";
	} else {
		tau_iso = "Tau_idMVAoldDM2017v2";
	}
    
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
    
    if (cfg.find("trigger") != cfg.end())			trigger = cfg["trigger"];
    if (cfg.find("runOnData") != cfg.end())			runOnData = cfg["runOnData"];
    if (cfg.find("era") != cfg.end())				era = cfg["era"];
    
    if (!runOnData) { 
		auto getrightbincontent = [](TH2D* hist, const int binX, const int binY, std::string type) {
			if (type == "") {
				return hist->GetBinContent(	hist->GetXaxis()->FindBin( binX ),
											hist->GetYaxis()->FindBin( binY ) );
			} else if (type == "Up") {
				return hist->GetBinContent(	hist->GetXaxis()->FindBin( binX ),
											hist->GetYaxis()->FindBin( binY ) )
							+ hist->GetBinErrorUp(	hist->GetXaxis()->FindBin( binX ),
													hist->GetYaxis()->FindBin( binY ) );
			} else if (type == "Down") {
				return hist->GetBinContent(	hist->GetXaxis()->FindBin( binX ),
											hist->GetYaxis()->FindBin( binY ) )
							- hist->GetBinErrorLow(	hist->GetXaxis()->FindBin( binX ),
													hist->GetYaxis()->FindBin( binY ) );
			}
			return 1.0;
		};
		// read in pileup file
		TFile* pileup_file = new TFile(((std::string) cfg["pileup_file"]).c_str(), "READ");
		pileup_hist = (TH1D*) pileup_file->Get("pileup");
		
		// read in tau scale factor file
		TFile* tau_scale_file = new TFile(((std::string) cfg["tau_scale_factor_file"]).c_str(), "READ");
		TH2D* tau_scale_hist = (TH2D*) tau_scale_file->Get("tau_scale_factor");		
		
		tau_scale = getrightbincontent(tau_scale_hist, tau_iso_WP, era, "");											
		tau_scale_up = getrightbincontent(tau_scale_hist, tau_iso_WP, era, "Up");												
		tau_scale_down = getrightbincontent(tau_scale_hist, tau_iso_WP, era, "Down");													
		delete tau_scale_hist;
		delete tau_scale_file;  
		
		TFile* tau_energy_scale_file = new TFile(((std::string) cfg["tau_energy_scale_file"]).c_str(), "READ");
		TH2D* tau_energy_scale_hist = (TH2D*) tau_energy_scale_file->Get("tau_energy_scale");
		
		tau_energy_scale["h+-"][0] = getrightbincontent(tau_energy_scale_hist, 1, era, "Up");
		tau_energy_scale["h+-"][1] = getrightbincontent(tau_energy_scale_hist, 1, era, "");
		tau_energy_scale["h+-"][2] = getrightbincontent(tau_energy_scale_hist, 1, era, "Down");
		
		tau_energy_scale["h+- pi0s"][0] = getrightbincontent(tau_energy_scale_hist, 2, era, "Up");
		tau_energy_scale["h+- pi0s"][1] = getrightbincontent(tau_energy_scale_hist, 2, era, "");
		tau_energy_scale["h+- pi0s"][2] = getrightbincontent(tau_energy_scale_hist, 2, era, "Down");
		
		tau_energy_scale["h+-h+-h+-"][0] = getrightbincontent(tau_energy_scale_hist, 3, era, "Up");
		tau_energy_scale["h+-h+-h+-"][1] = getrightbincontent(tau_energy_scale_hist, 3, era, "");
		tau_energy_scale["h+-h+-h+-"][2] = getrightbincontent(tau_energy_scale_hist, 3, era, "Down");
		
		tau_energy_scale["h+-h+-h+- pi0s"][0] = getrightbincontent(tau_energy_scale_hist, 4, era, "Up");
		tau_energy_scale["h+-h+-h+- pi0s"][1] = getrightbincontent(tau_energy_scale_hist, 4, era, "");
		tau_energy_scale["h+-h+-h+- pi0s"][2] = getrightbincontent(tau_energy_scale_hist, 4, era, "Down");
		
		delete tau_energy_scale_hist; 
		delete tau_energy_scale_file;
		
		if (cfg.find("W_kfactor_file") != cfg.end() ) {
			TFile* W_kfactor_file = new TFile(((std::string) cfg["W_kfactor_file"]).c_str(), "READ");
			W_kfactor_hist = (TH1D*) W_kfactor_file->Get("h_t_kfac_add");
		}
		
		if (cfg.find("tau_ele_fake_scale_file") != cfg.end() ) {
			TFile* tau_ele_fake_scale_file = new TFile(((std::string) cfg["tau_ele_fake_scale_file"]).c_str(), "READ");
			tau_ele_fake_hist = (TH1D*) tau_ele_fake_scale_file->Get("tau_ele_fake_rate");
		}
		
		if (cfg.find("tau_muon_fake_scale_file") != cfg.end() ) {
			TFile* tau_muon_fake_scale_file = new TFile(((std::string) cfg["tau_muon_fake_scale_file"]).c_str(), "READ");
			tau_muo_fake_hist = (TH1D*) tau_ele_fake_scale_file->Get("tau_muon_fake_rate");
		}
	}
}
