#ifndef CONFIG_H
#define CONFIG_H

#include "include/json.hpp"
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <map>
using json = nlohmann::json;

namespace config {
	extern std::string trigger;
	
	extern float pv_z;
	extern float pv_d;
	extern float pv_ndof;

	extern float tau_pt;
	extern float tau_eta;
	extern std::string tau_dm;
	extern std::string tau_iso;
	extern std::string tau_antiEle;
	extern std::string tau_antiMuon;
    extern uint tau_iso_WP;
    extern uint tau_antiE_WP;
    extern uint tau_antiMu_WP;
    
	extern float muon_pt;
	extern float muon_eta;
	extern float ele_pt;
	extern float ele_eta;
	extern float met_pt;
    
    extern int gen_pdgID;
    extern std::string cut_type;
    extern float cut_value_min;
    extern float cut_value_max;
    
    extern bool runOnData;
    extern int era;
    
    extern TH1D* pileup_hist;
    extern TH2D* prefire_photon_hist;
    extern TH2D* prefire_jet_hist;
    
    extern std::string run_type;
    
    extern TF1* tau_scale;
    extern TF1* tau_scale_up;
    extern TF1* tau_scale_down;
    
    extern TH1D* tau_energy_scale_hist;
    
    extern TH1D* W_kfactor_hist;
    
    extern TH1D* tau_ele_fake_hist;
    extern TH1D* tau_muo_fake_hist;

    
    bool load_config_file(json cfg);
    void clean_memory(); 
}

#endif
