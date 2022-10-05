#include "config.hh"
#include <fstream>
#include <sstream>
#include <iostream>
#include <TFile.h>
#include <TROOT.h>
#include "TH2D.h"
#include "TF1.h"
#include "TPRegexp.h"

namespace config {
	// set default values for all used quantities
	
	// trigger string
	std::string trigger = "HLT_Mu50 || HLT_TkMu100 || HLT_OldMu100";
	// metfilters string
	std::string metfilters = "Flag_goodVertices & Flag_globalSuperTightHalo2016Filter & Flag_HBHENoiseFilter & Flag_HBHENoiseIsoFilter & Flag_EcalDeadCellTriggerPrimitiveFilter & Flag_BadPFMuonFilter & Flag_ecalBadCalibFilterV2";

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
    uint tau_iso_WP_datadriven = 16;
    uint tau_antiMu_WP_datadriven = 1;
    uint tau_antiE_WP_datadriven = 2;
    
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
    int gen_motherpdgID;
    std::string cut_type = "";
    float cut_value_min = 0;
    float cut_value_max = 999999;
    bool cut_double = false;
    bool cut_single = false;
    
    // pileup hist
    TH1D* pileup_hist = NULL;
    TH1D* pileup_hist_up = NULL;
    TH1D* pileup_hist_down = NULL;
    
    // trigger histograms
    TH2D* trigger_hist = NULL;
    TH2D* trigger_hist_up = NULL;
    TH2D* trigger_hist_down = NULL;
    // fakerate histograms
    TH2D* ff_hist = NULL;
    TH2D* ff_hist_low = NULL;
    TH2D* ff_hist_high = NULL;
    TH2D* ff_hist_barrel_low = NULL;
    TH2D* ff_hist_barrel_high = NULL;
    TH2D* ff_hist_endcap_low = NULL;
    TH2D* ff_hist_endcap_high = NULL;
    TH2D* ff_closure_hist = NULL;
    TH2D* ff_closure_hist_low = NULL;
    TH2D* ff_closure_hist_high = NULL;
    TH2D* ff_closure_hist_barrel_low = NULL;
    TH2D* ff_closure_hist_barrel_high = NULL;
    TH2D* ff_closure_hist_endcap_low = NULL;
    TH2D* ff_closure_hist_endcap_high = NULL;
    // Prefire histograms
    TH2D* prefire_photon_hist = NULL;
    TH2D* prefire_jet_hist = NULL;
    
    //WW shape uncertainty switch

    bool wwuncertainty = false;
    
    // W kfactor hist
    TH1D* W_kfactor_hist = NULL;

    //do xy
    bool doXY = false;

    // TT?
    bool TT = false;
    bool DY = false;
    bool WW = false;
    // runOnData?
    bool runOnData = false;
    // signal?
    bool runOnSignal = false;
    
    //do Snapshot
    bool doSnapshot = false;
    // runDataDriven
    bool runDataDriven = false;

    bool use2017XY = false;

    // calcDataDriven
    bool calcDataDriven = false;

    // closure
    bool closure = false;

    // which era?
    int era = 2018;
    bool use_EEMET = false;
    
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

    // electron uncertainty tool
    EnergyScaleCorrection *eleCorr;
    // pdf setup
    std::string pdf_set_name = "";
    std::string pdf_prod_set_name = "NNPDF31_nnlo_as_0118";
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
    if (cfg.find("globalconfig") == cfg.end()) {
        std::cout << "No global Configuration file given!" << std::endl;
        cfg_json.open("cfg/analysis_config.cfg");
        cfg_json >> globalcfg;
    }
	// load all numbers from the config file
    else{
        cfg_json.open(cfg["globalconfig"]);
        cfg_json >> globalcfg;
    }

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
    if (globalcfg.find("tau_isolation_WP_datadriven") != globalcfg.end())  tau_iso_WP_datadriven = globalcfg["tau_isolation_WP_datadriven"];
    if (globalcfg.find("tau_antiE_WP_datadriven") != globalcfg.end())      tau_antiE_WP_datadriven = globalcfg["tau_antiE_WP_datadriven"];
    if (globalcfg.find("tau_antiMu_WP_datadriven") != globalcfg.end())     tau_antiMu_WP_datadriven = globalcfg["tau_antiMu_WP_datadriven"];
		
    if (globalcfg.find("muon_id") != globalcfg.end())           muon_id = globalcfg["muon_id"];
    if (globalcfg.find("muon_iso") != globalcfg.end())          muon_iso = globalcfg["muon_iso"];
    if (globalcfg.find("muon_id_WP") != globalcfg.end())        muon_id_WP = globalcfg["muon_id_WP"];
    if (globalcfg.find("muon_iso_WP") != globalcfg.end())       muon_iso_WP = globalcfg["muon_iso_WP"];
    if (globalcfg.find("doXY") != globalcfg.end())		    	doXY = globalcfg["doXY"];
    if (globalcfg.find("doSnapshot") != globalcfg.end())		doSnapshot = globalcfg["doSnapshot"];
    if (globalcfg.find("closure") != globalcfg.end())	        closure = globalcfg["closure"];

    if (cfg.find("gen_part_pdg_id") != cfg.end()){
        std::string gen_pdgID_str = cfg["gen_part_pdg_id"];
        std::stringstream ss(gen_pdgID_str);
        for (int i; ss >> i;) {
              gen_pdgID.push_back(i);    
              if (ss.peek() == ',')
                  ss.ignore();
         }
    }
    if (cfg.find("gen_part_mother_pdg_id") != cfg.end())		gen_motherpdgID = cfg["gen_part_mother_pdg_id"];
    if (cfg.find("gen_cut_type") != cfg.end())		cut_type = cfg["gen_cut_type"];
    if (cfg.find("gen_cut_double") != cfg.end())	cut_double = cfg["gen_cut_double"];
    if (cfg.find("gen_cut_single") != cfg.end())	cut_single = cfg["gen_cut_single"];
    if (cfg.find("gen_cut_min") != cfg.end())		cut_value_min = cfg["gen_cut_min"];
    if (cfg.find("gen_cut_max") != cfg.end())		cut_value_max = cfg["gen_cut_max"];
    
    if (cfg.find("wwuncertainty") != cfg.end())	    wwuncertainty = cfg["wwuncertainty"];
    if (cfg.find("runOnData") != cfg.end())			runOnData = cfg["runOnData"];
    if (cfg.find("runOnSignal") != cfg.end())	    runOnSignal = cfg["runOnSignal"];
    if (cfg.find("TT") != cfg.end())		    	TT = cfg["TT"];
    if (cfg.find("DY") != cfg.end())		    	DY = cfg["DY"];
    if (cfg.find("WW") != cfg.end())		    	WW = cfg["WW"];
    if (cfg.find("runDataDriven") != cfg.end())		runDataDriven = cfg["runDataDriven"];
    if (cfg.find("use2017XY") != cfg.end())		    use2017XY = cfg["use2017XY"];
    if (cfg.find("calcDataDriven") != cfg.end())	calcDataDriven = cfg["calcDataDriven"];
    if (cfg.find("era") != cfg.end())				era = cfg["era"];
    if (cfg.find("use_EEMET") != cfg.end())			use_EEMET = cfg["use_EEMET"];
    
    if (cfg.find("PDF_set_name") != cfg.end())      pdf_set_name = cfg["PDF_set_name"];
    if (cfg.find("PDF_prod_set_name") != cfg.end())      pdf_prod_set_name = cfg["PDF_prod_set_name"];
    if (cfg.find("PDF_nWeights") != cfg.end())      pdf_nweights = cfg["PDF_nWeights"];
    if (cfg.find("PDF_SetID") != cfg.end())         pdf_setid = cfg["PDF_SetID"];
    if (cfg.find("metfilters") != cfg.end())        metfilters = cfg["metfilters"];
    if (cfg.find("trigger") != cfg.end())			trigger = cfg["trigger"];
    if (cfg.find("closure") != cfg.end())	        closure = cfg["closure"];
    TFile* ff_file = new TFile(("cfg/fakerate/fakerate_" + std::to_string(era)+ "_etau" + ".root").c_str(), "READ");
    TFile* ff_closure_file = new TFile(("cfg/fakerate/fakerate_" + std::to_string(era) + "_closure" + "_etau" + ".root").c_str(), "READ");
    TString ff_hist_name = "iso_fake_rate";
    TString ff_hist_name_low = "iso_fake_rate_low";
    TString ff_hist_name_high = "iso_fake_rate_high";
    TString ff_hist_name_barrel_low = "iso_fake_rate_barrel_low";
    TString ff_hist_name_barrel_high = "iso_fake_rate_barrel_high";
    TString ff_hist_name_endcap_low = "iso_fake_rate_endcap_low";
    TString ff_hist_name_endcap_high = "iso_fake_rate_endcap_high";
    std::cout << "Reading in hist " << ff_hist_name << std::endl;
    ff_hist = (TH2D*) ff_file->Get(ff_hist_name);
    ff_hist_low = (TH2D*) ff_file->Get(ff_hist_name_low);
    ff_hist_high = (TH2D*) ff_file->Get(ff_hist_name_high);
    ff_hist_barrel_low = (TH2D*) ff_file->Get(ff_hist_name_barrel_low);
    ff_hist_barrel_high = (TH2D*) ff_file->Get(ff_hist_name_barrel_high);
    ff_hist_endcap_low = (TH2D*) ff_file->Get(ff_hist_name_endcap_low);
    ff_hist_endcap_high = (TH2D*) ff_file->Get(ff_hist_name_endcap_high);
    ff_closure_hist = (TH2D*) ff_closure_file->Get(ff_hist_name);
    ff_closure_hist_low = (TH2D*) ff_closure_file->Get(ff_hist_name_low);
    ff_closure_hist_high = (TH2D*) ff_closure_file->Get(ff_hist_name_high);
    ff_closure_hist_barrel_low = (TH2D*) ff_closure_file->Get(ff_hist_name_barrel_low);
    ff_closure_hist_barrel_high = (TH2D*) ff_closure_file->Get(ff_hist_name_barrel_high);
    ff_closure_hist_endcap_low = (TH2D*) ff_closure_file->Get(ff_hist_name_endcap_low);
    ff_closure_hist_endcap_high = (TH2D*) ff_closure_file->Get(ff_hist_name_endcap_high);
    
    if (!runOnData) { 
		// read in pileup file
		TFile* pileup_file = new TFile(((std::string) cfg["pileup_file"]).c_str(), "READ");
		pileup_hist = (TH1D*) pileup_file->Get("pileup");
		pileup_hist_up = (TH1D*) pileup_file->Get("pileup_up");
		pileup_hist_down = (TH1D*) pileup_file->Get("pileup_down");
		

        // read in prefiring jet file
        if (cfg.find("prefire_jet_hist") != cfg.end()) {
            TPRegexp r1("[^/]+(?=/$|$)");
            TFile* prefire_jet_file = new TFile(((std::string) cfg["prefire_jet_hist"]).c_str(), "READ");
            TString pj_hist_name = (std::string) cfg["prefire_jet_hist"];
            pj_hist_name.Resize( pj_hist_name.Length() - 5);
            pj_hist_name = pj_hist_name(r1);

            std::cout << "Reading in hist " << pj_hist_name << std::endl;
            prefire_jet_hist = (TH2D*) prefire_jet_file->Get(pj_hist_name);
        }

        // read in prefiring photon file
        if (cfg.find("prefire_photon_hist") != cfg.end()) {
            TPRegexp r1("[^/]+(?=/$|$)");
            TFile* prefire_photon_file = new TFile(((std::string) cfg["prefire_photon_hist"]).c_str(), "READ");
            TString pp_hist_name = (std::string) cfg["prefire_photon_hist"];
            pp_hist_name.Resize( pp_hist_name.Length() - 5);
            pp_hist_name = pp_hist_name(r1);
            std::cout << "Reading in hist " << pp_hist_name << std::endl;
            prefire_photon_hist = (TH2D*) prefire_photon_file->Get(pp_hist_name);
        }
        if (cfg.find("trigger_hist") != cfg.end() && cfg.find("trigger_file") != cfg.end()) {
            TPRegexp r1("[^/]+(?=/$|$)");
            TFile* trigger_file = new TFile(((std::string) cfg["trigger_file"]).c_str(), "READ");
            TString trigger_hist_name = (std::string) cfg["trigger_hist"];
            std::cout << "Reading in hist " << trigger_hist_name << std::endl;
            trigger_hist = (TH2D*) trigger_file->Get(trigger_hist_name);
            if (cfg.find("trigger_hist_up") != cfg.end() ) {
                TString trigger_hist_name = (std::string) cfg["trigger_hist_up"];
                std::cout << "Reading in hist " << trigger_hist_name << std::endl;
                trigger_hist_up = (TH2D*) trigger_file->Get(trigger_hist_name);
            }
            if (cfg.find("trigger_hist_down") != cfg.end() ) {
                TString trigger_hist_name = (std::string) cfg["trigger_hist_down"];
                trigger_hist_down = (TH2D*) trigger_file->Get(trigger_hist_name);
            }
        }

//		// read in prefiring jet file
//		if (cfg.find("prefire_jet_hist") != cfg.end() and cfg.find("prefire_jet_hist_name") != cfg.end()) {
//			TFile* prefire_jet_file = new TFile(((std::string) cfg["prefire_jet_hist"]).c_str(), "READ");
//			std::string pj_hist_name = (std::string) cfg["prefire_jet_hist_name"];
//			prefire_jet_hist = (TH2D*) prefire_jet_file->Get(pj_hist_name.c_str());
//		}
//		
//		// read in prefiring photon file
//		if (cfg.find("prefire_photon_hist") != cfg.end() and cfg.find("prefire_photon_hist_name") != cfg.end()) {
//			TFile* prefire_photon_file = new TFile(((std::string) cfg["prefire_photon_hist"]).c_str(), "READ");
//			std::string pp_hist_name = (std::string) cfg["prefire_jet_hist_name"];
//			prefire_photon_hist = (TH2D*) prefire_photon_file->Get(pp_hist_name.c_str());
//		}
			
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
        if(era == 2016){
            eleCorr= new EnergyScaleCorrection((std::string) getenv("MY_ANALYSIS_PATH") + "/cfg/SF/EleSFs/Legacy2016_07Aug2017_FineEtaR9_v3_ele",EnergyScaleCorrection::ECALELF);
        }
        else if(era == 2017){
            eleCorr= new EnergyScaleCorrection((std::string) getenv("MY_ANALYSIS_PATH") + "/cfg/SF/EleSFs/Run2017_17Nov2017_v1_ele_unc",EnergyScaleCorrection::ECALELF);
        }
        else if(era == 2018){
            eleCorr= new EnergyScaleCorrection((std::string) getenv("MY_ANALYSIS_PATH") + "/cfg/SF/EleSFs/Run2018_Step2Closure_CoarseEtaR9Gain_v2",EnergyScaleCorrection::ECALELF);
        }

	} //!runOnData
    return true;
}

// reset config parameters, delete pointers etc.
void config::clean_memory() {
	//// trigger string
	//trigger = "";
	//// METFilter string
	//metfilters = "";

	//// primary vertex cuts
    //pv_z = 24;
    //pv_d = 2;
    //pv_ndof = 4;
	//
	//// minimum tau pt
    //tau_pt = 80;
    //
    //// maximum tau eta
    //tau_eta = 2.3;
    //
	//// minimum muon pt
    //muon_pt = 20;
    //
    //// maximum muon eta
    //muon_eta = 2.4;
    //
	//// minimum ele pt
    //ele_pt = 20;
    //
    //// maximum ele eta
    //ele_eta = 2.5;
    //
    //// working points of tau identification
    //tau_dm = "Tau_idDecayModeNewDMs";
    //tau_iso = "Tau_idMVAnewDM2017v2";
    //tau_antiEle = "Tau_idAntiEle";
	//tau_antiMuon = "Tau_idAntiMu";
    //tau_iso_WP = 1;
    //tau_antiMu_WP = 1;
    //tau_antiE_WP = 1;
    //
    //// muon id 
    //muon_id = "HighPtID";
    //muon_id_WP = 2;

    //// minimum missing transverse momentum
    //met_pt = 150;
    //
    //// gen particle cut procedure
    //gen_pdgID.clear();
    //cut_type = "";
    //cut_value_min = 0;
    //cut_value_max = 999999;
    //cut_double = false;
    //cut_single = false;
    //
    //// pileup hist
    //delete pileup_hist;
    //delete pileup_hist_up;
    //delete pileup_hist_down;
    //pileup_hist = NULL;
    //pileup_hist_up = NULL;
    //pileup_hist_down = NULL;
    //
    //// Prefire histograms
    //delete prefire_photon_hist;
    //delete prefire_jet_hist;
    //prefire_photon_hist = NULL;
    //prefire_jet_hist = NULL;
    //
    //// W kfactor hist
    //delete W_kfactor_hist;
    //W_kfactor_hist = NULL;
    //
    //// runOnData?
    //runOnData = false;

    //// runDataDriven
    //runDataDriven = false;
    //
    //// calcDataDriven
    //calcDataDriven = false;

    //// which era?
    //era = 2018;
    //
    //use_EEMET = false;
    //
    //// this checks systematics later
    //run_type = "";
    //
    //// lets start with systematics and scale factors
    //
    //// tau energy scale    
    //delete tau_dm_scale;
    //delete tau_vsjet_SF;
    //delete tau_vsmu_SF;
    //delete tau_vse_SF;
    //tau_dm_scale = NULL;
    //tau_vsjet_SF = NULL;
    //tau_vsmu_SF = NULL;
    //tau_vse_SF = NULL;
    //// ele id scale factor
    //delete ele_SF;
    //ele_SF = NULL;

    //// muon id scale factor
    //delete muon_SF;
    //muon_SF = NULL;
    //
    //// pdf setup
    //pdf_set_name = "";
    //pdf_nweights = 0;
    //pdf_setid = 0;
    //pdf_is_initialized = false;
    //init_pdf_sets.clear();
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

    // minimum missing transverse momentum
    met_pt = 150;

    // gen particle cut procedure
    gen_pdgID.clear();
    gen_motherpdgID = 0;
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

    // trigger histograms
    delete trigger_hist;
    trigger_hist = NULL;
    delete trigger_hist_up;
    trigger_hist_up = NULL;
    delete trigger_hist_down;
    trigger_hist_down = NULL;
    // fakerate histogram
    delete ff_hist;
    ff_hist = NULL;
    delete ff_hist_low;
    ff_hist_low = NULL;
    delete ff_hist_high;
    ff_hist_high = NULL;
    delete ff_hist_barrel_high;
    ff_hist_barrel_high = NULL;
    delete ff_hist_endcap_high;
    ff_hist_endcap_high = NULL;
    delete ff_closure_hist;
    ff_closure_hist = NULL;
    delete ff_closure_hist_low;
    ff_closure_hist_low = NULL;
    delete ff_closure_hist_high;
    ff_closure_hist_high = NULL;
    delete ff_closure_hist_barrel_high;
    ff_closure_hist_barrel_high = NULL;
    delete ff_closure_hist_endcap_high;
    ff_closure_hist_endcap_high = NULL;
    // Prefire histograms
    delete prefire_photon_hist;
    delete prefire_jet_hist;
    prefire_photon_hist = NULL;
    prefire_jet_hist = NULL;

    doXY = false;
    // TT?
    TT = false;
    DY = false;
    WW = false;
    // runOnData?
    runOnData = false;
    runOnSignal = false;
    calcDataDriven = false;
    runDataDriven = false;
    use2017XY = false;
    doSnapshot = false;
    closure = false;
    
    //ww uncertainty switch
    wwuncertainty = false;

    // which era?
    era = 2018;
    use_EEMET = false;

    // this checks systematics later
    metfilters = "";
    

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

    delete ele_SF;
    delete muon_SF;
    delete eleCorr;
    ele_SF = NULL;
    muon_SF = NULL;
    eleCorr = NULL;

    // pdf setup
    pdf_set_name = "";
    pdf_prod_set_name = "";
    pdf_nweights = 0;
    pdf_setid = 0;
    pdf_is_initialized = false;
    init_pdf_sets.clear();

}
    

