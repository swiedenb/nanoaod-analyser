#ifndef CONFIG_H
#define CONFIG_H

#include "json.hpp"
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <map>
#include "TauIDSFTool.h"
#include "MuonSFTool.h"
#include "EleSFTool.h"
#include "LHAPDF/LHAPDF.h"
using json = nlohmann::json;

namespace config {
	extern std::string trigger;
	extern std::string metfilters;
	
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
    extern uint tau_iso_WP_datadriven;
    extern uint tau_antiE_WP_datadriven;
    extern uint tau_antiMu_WP_datadriven;

    extern std::string muon_id;
    extern std::string muon_iso;
    extern uint muon_id_WP;
    extern uint muon_iso_WP;
    
	extern float muon_pt;
	extern float muon_eta;
	extern float ele_pt;
	extern float ele_eta;
	extern float met_pt;
    
    extern std::vector<int> gen_pdgID;
    extern int gen_motherpdgID;
    extern std::string cut_type;
    extern float cut_value_min;
    extern float cut_value_max;
    extern bool cut_double;
    extern bool cut_single;
    
    extern bool TT;
    extern bool doXY;
    extern bool WW;
    extern bool DY;
    extern bool runOnData;
    extern bool runOnSignal;
    extern bool wwuncertainty;
    extern bool runDataDriven;
    extern bool use2017XY;
    extern bool doSnapshot;
    extern bool calcDataDriven;
    extern bool closure;
    extern int era;

    extern bool use_EEMET;
    
    extern TH1D* pileup_hist;
    extern TH1D* pileup_hist_up;
    extern TH1D* pileup_hist_down;
    extern TH2D* trigger_hist;
    extern TH2D* trigger_hist_up;
    extern TH2D* trigger_hist_down;
    extern TH2D* ff_hist;
    extern TH2D* ff_hist_low;
    extern TH2D* ff_hist_high;
    extern TH2D* ff_hist_barrel_low;
    extern TH2D* ff_hist_barrel_high;
    extern TH2D* ff_hist_endcap_low;
    extern TH2D* ff_hist_endcap_high;
    extern TH2D* ff_closure_hist;
    extern TH2D* ff_closure_hist_low;
    extern TH2D* ff_closure_hist_high;
    extern TH2D* ff_closure_hist_barrel_low;
    extern TH2D* ff_closure_hist_barrel_high;
    extern TH2D* ff_closure_hist_endcap_low;
    extern TH2D* ff_closure_hist_endcap_high;
    extern TH2D* prefire_photon_hist;
    extern TH2D* prefire_jet_hist;
    
    extern std::string run_type;
    
    
    extern TH1D* tau_energy_scale_hist;
    
    extern TH1D* W_kfactor_hist;
    
    extern TauESTool* tau_dm_scale;
    extern TauIDSFTool* tau_vsjet_SF;
    extern TauIDSFTool* tau_vsmu_SF;
    extern TauIDSFTool* tau_vse_SF;

    extern MuonSFTool* muon_SF;
    extern EleSFTool* ele_SF;

    extern std::string pdf_set_name;
    extern std::string pdf_prod_set_name;
    extern int pdf_nweights;
    extern int pdf_setid;
    extern bool pdf_is_initialized;
    extern std::vector< LHAPDF::PDF* > init_pdf_sets;

    
    bool load_config_file(json cfg);
    void clean_memory(); 
}

#endif
