#include "config.hh"
#include <fstream>
#include <sstream>
#include <TFile.h>
#include <TROOT.h>
#include "TH2D.h"
#include "TF1.h"

namespace config {
	// set default values for all used quantities
	
	// trigger string
	std::string trigger = "HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr || HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1";

	// primary vertex cuts
    float pv_z = 24;
    float pv_d = 2;
    float pv_ndof = 4;
	
	// minimum tau pt
    float tau_pt = 80;
    
    // maximum tau eta
    float tau_eta = 2.3;
    
	// minimum muon pt
    float muon_pt = 20;
    
    // maximum muon eta
    float muon_eta = 2.4;
    
	// minimum ele pt
    float ele_pt = 20;
    
    // maximum ele eta
    float ele_eta = 2.5;
    
    // working points of tau identification
    std::string tau_dm = "Tau_idDecayModeNewDMs";
    std::string tau_iso = "Tau_idMVAnewDM2017v2";
	std::string tau_antiEle = "Tau_idAntiEle";
	std::string tau_antiMuon = "Tau_idAntiMu";
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
    
    // Prefire histograms
    TH2D* prefire_photon_hist = NULL;
    TH2D* prefire_jet_hist = NULL;
    
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
    
    // lets start with systematics and scale factors
    
    // tau scale factors
    TF1* tau_scale = NULL;
    TF1* tau_scale_up = NULL;
    TF1* tau_scale_down = NULL;
    
    // tau energy scale
    TH1D* tau_energy_scale_hist = NULL;
}

bool config::load_config_file(json cfg)
{
	// load all numbers from the config file
	json globalcfg;
    std::ifstream cfg_json;
    cfg_json.open("cfg/analysis_config.cfg");
    cfg_json >> globalcfg;

    auto int_to_WP = [] (const int WP_int) {
        switch (WP_int) {
            case 1:
                return (std::string) "VVVLoose";
            case 2:
                return (std::string) "VVLoose";
            case 4:
                return (std::string) "VLoose";
            case 8:
                return (std::string) "Loose";
            case 16:
                return (std::string) "Medium";
            case 32:
                return (std::string) "Tight";
            case 64:
                return (std::string) "VTight";
            case 128:
                return (std::string) "VVTight";
        };
    };
    
    // always check, if the key is actually in the config file
    if (globalcfg.find("PV_Z") != globalcfg.end())
  pv_z = globalcfg["PV_Z"];
    if (globalcfg.find("PV_D") != globalcfg.end())
  pv_d = globalcfg["PV_D"];
    if (globalcfg.find("PV_NDOF") != globalcfg.end())
  pv_ndof = globalcfg["PV_NDOF"];
    if (globalcfg.find("tau_pt") != globalcfg.end()) 			tau_pt = globalcfg["tau_pt"];
    if (globalcfg.find("tau_eta") != globalcfg.end()) 			tau_eta = globalcfg["tau_eta"];
    if (globalcfg.find("muon_pt") != globalcfg.end()) 			muon_pt = globalcfg["muon_pt"];
    if (globalcfg.find("muon_eta") != globalcfg.end()) 			muon_eta = globalcfg["muon_eta"];
    if (globalcfg.find("electron_pt") != globalcfg.end()) 		ele_pt = globalcfg["electron_pt"];
    if (globalcfg.find("electron_eta") != globalcfg.end()) 		ele_eta = globalcfg["electron_eta"];
    if (globalcfg.find("met_pt") != globalcfg.end()) 			met_pt = globalcfg["met_pt"];
    
    if (globalcfg.find("tau_decayMode") != globalcfg.end()) 	tau_dm = globalcfg["tau_decayMode"];
    if (globalcfg.find("tau_isolation") != globalcfg.end()) 	tau_iso = globalcfg["tau_isolation"];
    if (globalcfg.find("tau_antiEle") != globalcfg.end()) 		tau_antiEle = globalcfg["tau_antiEle"];
    if (globalcfg.find("tau_antiMuon") != globalcfg.end()) 		tau_antiMuon = globalcfg["tau_antiMuon"];
    
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
		// read in pileup file
		TFile* pileup_file = new TFile(((std::string) cfg["pileup_file"]).c_str(), "READ");
		pileup_hist = (TH1D*) pileup_file->Get("pileup");
		
		// read in prefiring jet file
		if (cfg.find("prefire_jet_hist") != cfg.end()) {
			TFile* prefire_jet_file = new TFile(((std::string) cfg["prefire_jet_hist"]).c_str(), "READ");
			std::string pj_hist_name = (std::string) cfg["prefire_jet_hist"];
			pj_hist_name.resize( pj_hist_name.size() - 5);
			prefire_jet_hist = (TH2D*) prefire_jet_file->Get(pj_hist_name.c_str());
		}
		
		// read in prefiring photon file
		if (cfg.find("prefire_photon_hist") != cfg.end()) {
			TFile* prefire_photon_file = new TFile(((std::string) cfg["prefire_photon_hist"]).c_str(), "READ");
			std::string pp_hist_name = (std::string) cfg["prefire_jet_hist"];
			pp_hist_name.resize( pp_hist_name.size() - 5);
			prefire_photon_hist = (TH2D*) prefire_photon_file->Get(pp_hist_name.c_str());
		}
			
		// read in tau scale factor file
		if (globalcfg.find("tau_scale_factor_file") != globalcfg.end()) {
            TFile* tau_scale_file = new TFile(((std::string) globalcfg["tau_scale_factor_file"]).c_str(), "READ");
            tau_scale = (TF1*) tau_scale_file->Get( (int_to_WP(tau_iso_WP) + "_cent").c_str() );		
            tau_scale_down = (TF1*) tau_scale_file->Get( (int_to_WP(tau_iso_WP) + "_down").c_str() );		
            tau_scale_up = (TF1*) tau_scale_file->Get( (int_to_WP(tau_iso_WP) + "_up").c_str() );
        }
		
        // read in tau energy scale factor (decay mode dependent)
        if (globalcfg.find("tau_energy_scale_factor_file") != globalcfg.end()) {
            TFile* tau_energy_scale_file = new TFile(((std::string) globalcfg["tau_energy_scale_factor_file"]).c_str(), "READ");
            tau_energy_scale_hist = (TH1D*) tau_energy_scale_file->Get("tes");
        }
		
		if (globalcfg.find("tau_ele_fake_scale_file") != globalcfg.end() ) {
			TFile* tau_ele_fake_scale_file = new TFile(((std::string) globalcfg["tau_ele_fake_scale_file"]).c_str(), "READ");
			tau_ele_fake_hist = (TH1D*) tau_ele_fake_scale_file->Get("VLoose");
		}
		
		if (globalcfg.find("tau_muon_fake_scale_file") != globalcfg.end() ) {
			TFile* tau_muon_fake_scale_file = new TFile(((std::string) globalcfg["tau_muon_fake_scale_file"]).c_str(), "READ");
			tau_muo_fake_hist = (TH1D*) tau_muon_fake_scale_file->Get("Loose");
		}
    
		if (cfg.find("W_kfactor_file") != cfg.end() ) {
			TFile* W_kfactor_file = new TFile(((std::string) cfg["W_kfactor_file"]).c_str(), "READ");
			W_kfactor_hist = (TH1D*) W_kfactor_file->Get("h_t_kfac_add");
		}

	} //!runOnData
}
