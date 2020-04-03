#include "config.hh"
#include <fstream>
#include <sstream>
#include <iostream>
#include <TFile.h>
#include <TROOT.h>
#include "TH2D.h"
#include "TF1.h"

namespace config {
	// set default values for all used quantities
	
	// trigger string
	std::string trigger = "HLT_Mu50 || HLT_TkMu100 || HLT_OldMu100";

	// primary vertex cuts
    float pv_z = 24;
    float pv_d = 2;
    float pv_ndof = 4;
	
	// minimum tau pt
    float tau_pt = 30;
    
    // maximum tau eta
    float tau_eta = 2.3;
    
	// minimum muon pt
    float muon_pt = 53;
    
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
    uint tau_iso_WP = 16;
    uint tau_antiMu_WP = 1;
    uint tau_antiE_WP = 2;
    
    // muon id and iso 
    std::string muon_id = "HighPtID";
    std::string muon_iso = "LooseRelTkIso";
    uint muon_id_WP = 2;
    uint muon_iso_WP = 1;

    // ele id
    std::string ele_id = "heep";

    // minimum missing transverse momentum
    float met_pt = 150;
    
    // gen particle cut procedure
    std::vector<int> gen_pdgID;
    std::string cut_type = "";
    float cut_value_min = 0;
    float cut_value_max = 999999;
    
    // pileup hist
    TH1D* pileup_hist = NULL;
    TH1D* pileup_hist_up = NULL;
    TH1D* pileup_hist_down = NULL;
    
    // Prefire histograms
    TH2D* prefire_photon_hist = NULL;
    TH2D* prefire_jet_hist = NULL;
    
    // W kfactor hist
    TH1D* W_kfactor_hist = NULL;

    // runOnData?
    bool runOnData = false;
    
    // which era?
    int era = 2018;
    
    // this checks systematics later
    std::string run_type = "";
    
    // lets start with systematics and scale factors
    
    // tau energy scale
    TauESTool* tau_dm_scale;
    TauIDSFTool* tau_vsjet_SF;
    TauIDSFTool* tau_vsmu_SF;
    TauIDSFTool* tau_vse_SF;

    // ele id scale factor
    EleSFTool* ele_SF; 
    // muon scale factor
    MuonSFTool* muon_SF; 

    // pdf setup
    std::string pdf_set_name = "";
    int pdf_nweights = 0;
    int pdf_setid = 0;
    bool pdf_is_initialized = false;
    std::vector< LHAPDF::PDF* > init_pdf_sets;

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
    if (globalcfg.find("trigger") != globalcfg.end())			trigger = globalcfg["trigger"];

    if (globalcfg.find("PV_Z") != globalcfg.end())              pv_z = globalcfg["PV_Z"];
    if (globalcfg.find("PV_D") != globalcfg.end())              pv_d = globalcfg["PV_D"];
    if (globalcfg.find("PV_NDOF") != globalcfg.end())           pv_ndof = globalcfg["PV_NDOF"];

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
    
    if (globalcfg.find("tau_isolation_WP") != globalcfg.end())  tau_iso_WP = globalcfg["tau_isolation_WP"];
    if (globalcfg.find("tau_antiE_WP") != globalcfg.end())      tau_antiE_WP = globalcfg["tau_antiE_WP"];
    if (globalcfg.find("tau_antiMu_WP") != globalcfg.end())     tau_antiMu_WP = globalcfg["tau_antiMu_WP"];
		
    if (globalcfg.find("muon_id") != globalcfg.end())           muon_id = globalcfg["muon_id"];
    if (globalcfg.find("muon_iso") != globalcfg.end())          muon_iso = globalcfg["muon_iso"];
    if (globalcfg.find("muon_id_WP") != globalcfg.end())        muon_id_WP = globalcfg["muon_id_WP"];
    if (globalcfg.find("muon_iso_WP") != globalcfg.end())       muon_iso_WP = globalcfg["muon_iso_WP"];

    if (cfg.find("gen_part_pdg_id") != cfg.end()){
        std::string gen_pdgID_str = cfg["gen_part_pdg_id"];
        std::stringstream ss(gen_pdgID_str);
        for (int i; ss >> i;) {
              gen_pdgID.push_back(i);    
              if (ss.peek() == ',')
                  ss.ignore();
         }
    }
    if (cfg.find("gen_cut_type") != cfg.end())		cut_type = cfg["gen_cut_type"];
    if (cfg.find("gen_cut_min") != cfg.end())		cut_value_min = cfg["gen_cut_min"];
    if (cfg.find("gen_cut_max") != cfg.end())		cut_value_max = cfg["gen_cut_max"];
    
    if (cfg.find("runOnData") != cfg.end())			runOnData = cfg["runOnData"];
    if (cfg.find("era") != cfg.end())				era = cfg["era"];
    
    if (cfg.find("PDF_set_name") != cfg.end())      pdf_set_name = cfg["PDF_set_name"];
    if (cfg.find("PDF_nWeights") != cfg.end())      pdf_nweights = cfg["PDF_nWeights"];
    if (cfg.find("PDF_SetID") != cfg.end())         pdf_setid = cfg["PDF_SetID"];
    
    if (!runOnData) { 
		// read in pileup file
		TFile* pileup_file = new TFile(((std::string) cfg["pileup_file"]).c_str(), "READ");
		pileup_hist = (TH1D*) pileup_file->Get("pileup");
		pileup_hist_up = (TH1D*) pileup_file->Get("pileup_up");
		pileup_hist_down = (TH1D*) pileup_file->Get("pileup_down");
		
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
			
		if (cfg.find("W_kfactor_file") != cfg.end() ) {
			TFile* W_kfactor_file = new TFile(((std::string) cfg["W_kfactor_file"]).c_str(), "READ");
			W_kfactor_hist = (TH1D*) W_kfactor_file->Get("h_t_kfac_add");
		}
        auto era_to_SFToolYear = [](int era) {
            if (era == 2016)
                return (std::string) "2016Legacy";
            if (era == 2017)
                return (std::string) "2017ReReco";
            if (era == 2018)
                return (std::string) "2018ReReco";
        };
        
        
        
        
        std::string tau_iso_short = tau_iso;
        tau_iso_short.erase(0,6);
        
        std::string tau_vsmu_short = tau_antiMuon;
        tau_vsmu_short.erase(0,6);
        
        std::string tau_vse_short = tau_antiEle;
        tau_vse_short.erase(0,6);
        
        tau_dm_scale = new TauESTool(
                                era_to_SFToolYear(era),     // year
                                tau_iso_short);             // embedded sample?
        
        tau_vsjet_SF = new TauIDSFTool(
                                era_to_SFToolYear(era),     // year
                                tau_iso_short,              // ID
                                int_to_WP(tau_iso_WP),      // WP
                                false,                      // decay mode?
                                false);                     // embedded sample?
                                
        tau_vsmu_SF = new TauIDSFTool(
                                era_to_SFToolYear(era),     // year
                                tau_vsmu_short,             // ID
                                "Loose",                    // WP
                                false,                      // decay mode?
                                false);                     // embedded sample?


        tau_vse_SF = new TauIDSFTool(
                                era_to_SFToolYear(era),     // year
                                tau_vse_short,              // ID
                                int_to_WP(tau_antiE_WP),    // WP
                                false,                      // decay mode?
                                false);                     // embedded sample?
        muon_SF = new MuonSFTool(
                                era_to_SFToolYear(era),     // year
                                muon_id,
                                muon_iso);

        ele_SF = new EleSFTool(
                                era_to_SFToolYear(era),     // year
                                ele_id);

	} //!runOnData
}

// reset config parameters, delete pointers etc.
void config::clean_memory() {
	// trigger string
	trigger = "";

	// primary vertex cuts
    pv_z = 24;
    pv_d = 2;
    pv_ndof = 4;
	
	// minimum tau pt
    tau_pt = 80;
    
    // maximum tau eta
    tau_eta = 2.3;
    
	// minimum muon pt
    muon_pt = 20;
    
    // maximum muon eta
    muon_eta = 2.4;
    
	// minimum ele pt
    ele_pt = 20;
    
    // maximum ele eta
    ele_eta = 2.5;
    
    // working points of tau identification
    tau_dm = "Tau_idDecayModeNewDMs";
    tau_iso = "Tau_idMVAnewDM2017v2";
    tau_antiEle = "Tau_idAntiEle";
	tau_antiMuon = "Tau_idAntiMu";
    tau_iso_WP = 1;
    tau_antiMu_WP = 1;
    tau_antiE_WP = 1;
    
    // muon id 
    muon_id = "HighPtID";
    muon_id_WP = 2;

    // minimum missing transverse momentum
    met_pt = 150;
    
    // gen particle cut procedure
    gen_pdgID.clear();
    cut_type = "";
    cut_value_min = 0;
    cut_value_max = 999999;
    
    // pileup hist
    delete pileup_hist;
    delete pileup_hist_up;
    delete pileup_hist_down;
    pileup_hist = NULL;
    pileup_hist_up = NULL;
    pileup_hist_down = NULL;
    
    // Prefire histograms
    delete prefire_photon_hist;
    delete prefire_jet_hist;
    prefire_photon_hist = NULL;
    prefire_jet_hist = NULL;
    
    // W kfactor hist
    delete W_kfactor_hist;
    W_kfactor_hist = NULL;
    
    // runOnData?
    runOnData = false;
    
    // which era?
    era = 2018;
    
    // this checks systematics later
    run_type = "";
    
    // lets start with systematics and scale factors
    
    // tau energy scale    
    delete tau_dm_scale;
    delete tau_vsjet_SF;
    delete tau_vsmu_SF;
    delete tau_vse_SF;
    tau_dm_scale = NULL;
    tau_vsjet_SF = NULL;
    tau_vsmu_SF = NULL;
    tau_vse_SF = NULL;
    // ele id scale factor
    delete ele_SF;
    ele_SF = NULL;

    // muon id scale factor
    delete muon_SF;
    muon_SF = NULL;
    
    // pdf setup
    pdf_set_name = "";
    pdf_nweights = 0;
    pdf_setid = 0;
    pdf_is_initialized = false;
    init_pdf_sets.clear();
}
    

