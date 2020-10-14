#include "include/PhysicalQuantities.hh"
#include "include/Particles.hh"
#include "include/weights.hh"
#include "include/event_cleaner.hh"
#include "include/PDFTool.hh"
#include "include/TauIDSFTool.h"
#include "include/MuonSFTool.h"
#include "include/DDTool.h"
#include <iostream>
#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RVec.hxx>
#include <math.h>
#include <TSystemDirectory.h>
#include "TPRegexp.h"

// reduce amount to write for each vector by a lot
template < typename T >
using rvec = ROOT::VecOps::RVec<T>;
using RNode = ROOT::RDF::RNode;

json goldenjson;
json cfg;

std::multimap< std::string, ROOT::RDF::RResultPtr<TH1D> > hist_dict;
std::multimap< std::string, ROOT::RDF::RResultPtr<TH2D> > hist_dict_2d;

// Di muon pair number cut: how many "trues" in mask
bool deltar_tau_muon_cut(
                   const rvec<float>& tau_eta,
                   const rvec<float>& tau_phi,
                   const rvec<float>& muon_eta,
                   const rvec<float>& muon_phi,
                   const rvec<bool>& tau_mask,
                   const rvec<bool>& muon_mask,
                   const rvec<int>& col_idx) {
               double delR = sqrt(pow((tau_eta[tau_mask][col_idx[0]] - muon_eta[muon_mask][col_idx[1]]),2) + pow(tau_phi[tau_mask][col_idx[0]]-muon_phi[muon_mask][col_idx[1]],2));
               if (delR < 0.5){
                    return false;
               }
	return true;
};
// Di muon pair number cut: how many "trues" in mask
bool dimuonpair_cut(const rvec<float>& pt,
                   const rvec<float>& eta,
                   const rvec<float>& phi,
                   const rvec<bool>& mask) {
    if ( pt[mask].size()> 2){
        for (uint i = 0; i<pt[mask].size(); i++){
            for( uint j = i + 1; j < pt[mask].size(); j++){
               double delR = sqrt(pow((eta[mask][i] - eta[mask][j]),2) + pow(phi[mask][i]-phi[mask][j],2));
               if (delR > 0.2){
                    return false;
               }
            }
        }
    }
	return true;
};

// Particle number cut: how many "trues" in mask
bool nparticle_cut(const rvec<bool>& mask) {
	return std::count(mask.begin(), mask.end(), true) > 0.;
};
// Particle number veto: how many "trues" in mask
bool nparticle_veto(const rvec<bool>& mask) {
	return std::count(mask.begin(), mask.end(), true) == 0;
};

// Returns first mask_true particle for quantity
float selected_jet_quant(const rvec<int>& quantity, 
								const rvec<bool>& mask) {
	return quantity[mask][0];
};
// Returns first mask_true particle for quantity
float selected_part_col_idx_tau(const rvec<float>& quantity, 
								const rvec<bool>& mask,
                            const rvec<int>& col_idx) {
	return quantity[mask][col_idx[0]];
};
float selected_part_col_idx_muon(const rvec<float>& quantity, 
								const rvec<bool>& mask,
                            const rvec<int>& col_idx) {
	return quantity[mask][col_idx[1]];
};
// Returns first mask_true particle for quantity
float selected_part_quant(const rvec<float>& quantity, 
								const rvec<bool>& mask) {
	return quantity[mask][0];
};

// Returns first mask_true particle for quantity
rvec <float> select_part_quants(const rvec<float>& quantity, 
								const rvec<bool>& mask) {
	return quantity[mask];
};

// fill trigger plots
RNode trigger(	RNode df) {
	auto triggered = df.Filter(config::trigger, "Trigger requirement");
	if (!config::runOnData) {
		auto trigger_eff_tau_all = df.Histo1D( 	{"Tau_pt_all", "GenVisTau_pt", 600, 0, 6000}, "GenVisTau_pt");
		auto trigger_eff_met_all = df.Histo1D(	{"MET_pt_all", "MET_pt", 600, 0, 6000}, "MET_pt");
		auto trigger_eff_tau_passed = triggered.Histo1D( 	{"Tau_pt_passed", "Tau_pt", 600, 0, 6000}, "GenVisTau_pt");
		auto trigger_eff_met_passed = triggered.Histo1D(	{"MET_pt_passed", "MET_pt", 600, 0, 6000}, "MET_pt");
		hist_dict.emplace("Trigger", trigger_eff_tau_all);
		hist_dict.emplace("Trigger", trigger_eff_met_all);
		hist_dict.emplace("Trigger", trigger_eff_tau_passed);
		hist_dict.emplace("Trigger", trigger_eff_met_passed);
	};
	return triggered;
}

// fill preselection histograms
void fill_preselection(	RNode df) {
    if( !config::runOnData){
	    auto gentau_pt = df.Histo1D(	{"GenVisTau_pt", "", 			6000u, 0, 6000}, 					"GenVisTau_pt");
	    auto gen_pdgid = df.Histo1D(	{"GenPart_pdgId", "", 			200u, -100, 100}, 					"GenPart_pdgId");
	    auto gen_mass = df.Histo1D(	{"GenPart_mass", "", 			1000u, 0, 1000}, 					"GenPart_mass");
	    hist_dict.emplace("Preselection", gentau_pt);
	    hist_dict.emplace("Preselection", gen_pdgid);
	    hist_dict.emplace("Preselection", gen_mass);
    }
	auto tau_pt = df.Histo1D(	{"Tau_pt", "", 			6000u, 0, 6000}, 					"Tau_pt_ES");
	auto tau_eta = df.Histo1D(	{"Tau_eta", "", 		100u, -5, 5}, 						"Tau_eta_ES");
	auto tau_phi = df.Histo1D(	{"Tau_phi", "", 		100u, -3.2, 3.2}, 					"Tau_phi_ES");	
	auto muon_tP_pt = df.Histo1D(	{"Muon_tP_pt", "", 		6000u, 0, 6000}, 					"Muon_tP_pt");
	auto muon_pt = df.Histo1D(	{"Muon_pt", "", 		6000u, 0, 6000}, 					"Muon_pt");
	auto muon_eta = df.Histo1D(	{"Muon_eta", "", 		100u, -5, 5}, 						"Muon_eta");
	auto muon_phi = df.Histo1D(	{"Muon_phi", "", 		100u, -3.2, 3.2}, 					"Muon_phi");	
	auto met_pt = df.Histo1D(	{"MET_pt", "MET_pt", 	6000u, 0, 6000}, 					"MET_pt");
	auto met_phi = df.Histo1D(	{"MET_phi", "MET_phi", 	100u, -3.2, 3.2}, 					"MET_phi");
	hist_dict.emplace("Preselection", muon_pt);
	hist_dict.emplace("Preselection", muon_tP_pt);
	hist_dict.emplace("Preselection", muon_eta);
	hist_dict.emplace("Preselection", muon_phi);
	hist_dict.emplace("Preselection", tau_pt);
	hist_dict.emplace("Preselection", tau_eta);
	hist_dict.emplace("Preselection", tau_phi);
	hist_dict.emplace("Preselection", met_pt);
	hist_dict.emplace("Preselection", met_phi);
};


// fill datadriven histograms
void fill_datadriven(	RNode df,
                        std::string name = "") {
	auto tau_pt = df.Histo1D(	{"Tau_pt", "", 			6000u, 0, 6000}, 					"Tau_pt_new");
	auto tau_eta = df.Histo1D(	{"Tau_eta", "", 		100u, -5, 5}, 						"Tau_eta_new");
	auto tau_phi = df.Histo1D(	{"Tau_phi", "", 		100u, -3.2, 3.2}, 					"Tau_phi_new");	
	auto met_pt = df.Histo1D(	{"MET_pt", "MET_pt", 	6000u, 0, 6000}, 					"MET_pt");
	auto met_phi = df.Histo1D(	{"MET_phi", "MET_phi", 	100u, -3.2, 3.2}, 					"MET_phi");
	hist_dict.emplace("Datadriven_" + name, tau_pt);
	hist_dict.emplace("Datadriven_" + name, tau_eta);
	hist_dict.emplace("Datadriven_" + name, tau_phi);
	hist_dict.emplace("Datadriven_" + name, met_pt);
	hist_dict.emplace("Datadriven_" + name, met_phi);
};


// fill hist function - fills hists for each stage
void create_datadriven_hists(	RNode df,
                                std::string name,
                                std::string weight_column) {
    std::string stringcopy = weight_column;
    stringcopy.erase(0,12);
    //if( !config::runOnData){
	//    auto gentau_pt = df.Histo1D(	{((TString) "GenVisTau_pt" + stringcopy), "", 			6000u, 0, 6000}, 					"GenVisTau_pt",     weight_column);
	//    hist_dict.emplace(name, gentau_pt);
    //    
    //}
	auto tau_pt_ov_jet_pt_2d = df.Histo2D(	{"Tau_pt_vs_Tau_pt_over_Jet_pt", "", 	6000u, 0, 6000, 600u, 0,6}, 	"sel_Tau_pt_alt",	"Tau_pt_over_Jet_pt", 				weight_column);
	auto tau_pt_ov_jet_pt = df.Histo1D(	{((TString) "Tau_pt_over_Jet_pt" + stringcopy), "", 					6000u, 0, 1}, 		"Tau_pt_over_Jet_pt", 				weight_column);
	auto tau_pt = df.Histo1D(	{((TString) "Tau_pt" + stringcopy), "", 					6000u, 0, 6000}, 		"sel_Tau_pt", 				weight_column);
	auto tau_eta = df.Histo1D(	{((TString) "Tau_eta" + stringcopy), "", 					100u, -5, 5}, 			"sel_Tau_eta", 				weight_column);
	auto tau_phi = df.Histo1D(	{((TString) "Tau_phi" + stringcopy), "", 					100u, -3.2, 3.2}, 		"sel_Tau_phi", 				weight_column);
	auto muon_pt = df.Histo1D(	{((TString) "Muon_pt" + stringcopy), "", 					6000u, 0, 6000}, 		"sel_Muon_pt", 				weight_column);
	auto muon_eta = df.Histo1D(	{((TString) "Muon_eta" + stringcopy), "", 				100u, -5, 5}, 			"sel_Muon_eta", 			weight_column);
	auto muon_phi = df.Histo1D(	{((TString) "Muon_phi" + stringcopy), "", 				100u, -3.2, 3.2}, 		"sel_Muon_phi", 			weight_column);
	auto met_pt = df.Histo1D(	{((TString) "MET_pt" + stringcopy), "", 					6000u, 0, 6000}, 		"MET_pt", 					weight_column);
	auto met_phi = df.Histo1D(	{((TString) "MET_phi" + stringcopy), "", 					100u, -3.2, 3.2}, 		"MET_phi", 					weight_column);
	auto MT = df.Histo1D(		{((TString) "MT" + stringcopy), "", 						6000u, 0, 6000},		"MT", 						weight_column);
	auto coll_mass = df.Histo1D({((TString) "CollMass" + stringcopy), "",					6000u, 0, 6000},		"CollMass",					weight_column);
	auto coll_mass_alt = df.Histo1D({((TString) "CollMass_alt" + stringcopy), "",					6000u, 0, 6000},		"CollMass_alt",					weight_column);
	auto nvtx = df.Histo1D(		{((TString) "nvtx" + stringcopy), "", 					80u, 0, 80}, 			"PV_npvs", 					weight_column);
	auto njets = df.Histo1D(		{((TString) "nJets" + stringcopy), "", 					40u, 0, 40}, 			"nJet", 					weight_column);
	auto nvtxgood = df.Histo1D(	{((TString) "nvtx_good" + stringcopy), "", 				80u, 0, 80}, 			"PV_npvsGood", 				weight_column);
	//hist_dict.emplace(name, tau_genpartflav);
	hist_dict_2d.emplace(name, tau_pt_ov_jet_pt_2d);
	hist_dict.emplace(name, tau_pt_ov_jet_pt);
	hist_dict.emplace(name, njets);
	hist_dict.emplace(name, tau_pt);
	hist_dict.emplace(name, tau_eta);
	hist_dict.emplace(name, tau_phi);
	hist_dict.emplace(name, muon_pt);
	hist_dict.emplace(name, muon_eta);
	hist_dict.emplace(name, muon_phi);
	hist_dict.emplace(name, met_pt);
	hist_dict.emplace(name, met_phi);
	hist_dict.emplace(name, MT);
	hist_dict.emplace(name, coll_mass);
	hist_dict.emplace(name, coll_mass_alt);
	hist_dict.emplace(name, nvtx);
	hist_dict.emplace(name, nvtxgood);
};
void create_fr_hists(	RNode df,
                                std::string name,
                                std::string weight_column) {
    std::string stringcopy = weight_column;
    stringcopy.erase(0,12);
    //if( !config::runOnData){
	//    auto gentau_pt = df.Histo1D(	{((TString) "GenVisTau_pt" + stringcopy), "", 			6000u, 0, 6000}, 					"GenVisTau_pt",     weight_column);
	//    hist_dict.emplace(name, gentau_pt);
    //    
    //}
	auto tau_pt_ov_jet_pt = df.Histo1D(	{((TString) "Tau_pt_over_Jet_pt" + stringcopy), "", 					6000u, 0, 1}, 		"Tau_pt_over_Jet_pt", 				weight_column);
	auto tau_pt = df.Histo1D(	{((TString) "Tau_pt" + stringcopy), "", 					6000u, 0, 6000}, 		"sel_Tau_pt", 				weight_column);
	auto tau_eta = df.Histo1D(	{((TString) "Tau_eta" + stringcopy), "", 					100u, -5, 5}, 			"sel_Tau_eta", 				weight_column);
	auto tau_phi = df.Histo1D(	{((TString) "Tau_phi" + stringcopy), "", 					100u, -3.2, 3.2}, 		"sel_Tau_phi", 				weight_column);
	auto muon_pt = df.Histo1D(	{((TString) "Muon_pt" + stringcopy), "", 					6000u, 0, 6000}, 		"sel_Muon_pt", 				weight_column);
	auto muon_eta = df.Histo1D(	{((TString) "Muon_eta" + stringcopy), "", 				100u, -5, 5}, 			"sel_Muon_eta", 			weight_column);
	auto muon_phi = df.Histo1D(	{((TString) "Muon_phi" + stringcopy), "", 				100u, -3.2, 3.2}, 		"sel_Muon_phi", 			weight_column);
	auto met_pt = df.Histo1D(	{((TString) "MET_pt" + stringcopy), "", 					6000u, 0, 6000}, 		"MET_pt", 					weight_column);
	auto met_phi = df.Histo1D(	{((TString) "MET_phi" + stringcopy), "", 					100u, -3.2, 3.2}, 		"MET_phi", 					weight_column);
	auto MT = df.Histo1D(		{((TString) "MT" + stringcopy), "", 						6000u, 0, 6000},		"MT", 						weight_column);
	auto coll_mass = df.Histo1D({((TString) "CollMass" + stringcopy), "",					6000u, 0, 6000},		"CollMass",					weight_column);
	auto coll_mass_alt = df.Histo1D({((TString) "CollMass_alt" + stringcopy), "",					6000u, 0, 6000},		"CollMass_alt",					weight_column);
	auto nvtx = df.Histo1D(		{((TString) "nvtx" + stringcopy), "", 					80u, 0, 80}, 			"PV_npvs", 					weight_column);
	auto njets = df.Histo1D(		{((TString) "nJets" + stringcopy), "", 					40u, 0, 40}, 			"nJet", 					weight_column);
	auto nvtxgood = df.Histo1D(	{((TString) "nvtx_good" + stringcopy), "", 				80u, 0, 80}, 			"PV_npvsGood", 				weight_column);
	//hist_dict.emplace(name, tau_genpartflav);
	hist_dict.emplace(name, tau_pt_ov_jet_pt);
	hist_dict.emplace(name, njets);
	hist_dict.emplace(name, tau_pt);
	hist_dict.emplace(name, tau_eta);
	hist_dict.emplace(name, tau_phi);
	hist_dict.emplace(name, muon_pt);
	hist_dict.emplace(name, muon_eta);
	hist_dict.emplace(name, muon_phi);
	hist_dict.emplace(name, met_pt);
	hist_dict.emplace(name, met_phi);
	hist_dict.emplace(name, MT);
	hist_dict.emplace(name, coll_mass);
	hist_dict.emplace(name, coll_mass_alt);
	hist_dict.emplace(name, nvtx);
	hist_dict.emplace(name, nvtxgood);
};
void create_reso_hist( RNode df,
                       std::string name,
                       std::string weight_column) {
    std::string stringcopy = weight_column;
    if( !config::runOnData){
	    auto reco_reso = df.Histo1D(	{"Reco_coll_mass_reso", "", 			200u, -100, 100}, 					"Reco_coll_mass_reso",     weight_column);
	    hist_dict.emplace("Reso", reco_reso);
        
    }

                       }
// fill hist function - fills hists for each stage
void create_hists(	RNode df,
					std::string name,
                    std::string weight_column) {
    std::string stringcopy = weight_column;
    stringcopy.erase(0,12);
    //if( !config::runOnData){
	//    auto gentau_pt = df.Histo1D(	{((TString) "GenVisTau_pt" + stringcopy), "", 			6000u, 0, 6000}, 					"GenVisTau_pt",     weight_column);
	//    hist_dict.emplace(name, gentau_pt);
    //    
    //}
	auto tau_pt = df.Histo1D(	{((TString) "Tau_pt" + stringcopy), "", 					6000u, 0, 6000}, 		"sel_Tau_pt", 				weight_column);
	auto tau_eta = df.Histo1D(	{((TString) "Tau_eta" + stringcopy), "", 					100u, -5, 5}, 			"sel_Tau_eta", 				weight_column);
	auto tau_phi = df.Histo1D(	{((TString) "Tau_phi" + stringcopy), "", 					100u, -3.2, 3.2}, 		"sel_Tau_phi", 				weight_column);
	auto muon_pt = df.Histo1D(	{((TString) "Muon_pt" + stringcopy), "", 					6000u, 0, 6000}, 		"sel_Muon_pt", 				weight_column);
	auto muon_eta = df.Histo1D(	{((TString) "Muon_eta" + stringcopy), "", 				100u, -5, 5}, 			"sel_Muon_eta", 			weight_column);
	auto muon_phi = df.Histo1D(	{((TString) "Muon_phi" + stringcopy), "", 				100u, -3.2, 3.2}, 		"sel_Muon_phi", 			weight_column);
	auto met_pt = df.Histo1D(	{((TString) "MET_pt" + stringcopy), "", 					6000u, 0, 6000}, 		"MET_pt", 					weight_column);
	auto met_phi = df.Histo1D(	{((TString) "MET_phi" + stringcopy), "", 					100u, -3.2, 3.2}, 		"MET_phi", 					weight_column);
	auto MT = df.Histo1D(		{((TString) "MT" + stringcopy), "", 						6000u, 0, 6000},		"MT", 						weight_column);
	auto coll_mass = df.Histo1D({((TString) "CollMass" + stringcopy), "",					6000u, 0, 6000},		"CollMass",					weight_column);
	auto coll_mass_alt = df.Histo1D({((TString) "CollMass_alt" + stringcopy), "",					6000u, 0, 6000},		"CollMass_alt",					weight_column);
	auto nvtx = df.Histo1D(		{((TString) "nvtx" + stringcopy), "", 					80u, 0, 80}, 			"PV_npvs", 					weight_column);
	auto njets = df.Histo1D(		{((TString) "nJets" + stringcopy), "", 					40u, 0, 40}, 			"nJet", 					weight_column);
	auto nvtxgood = df.Histo1D(	{((TString) "nvtx_good" + stringcopy), "", 				80u, 0, 80}, 			"PV_npvsGood", 				weight_column);
	//hist_dict.emplace(name, tau_genpartflav);
	hist_dict.emplace(name, njets);
	hist_dict.emplace(name, tau_pt);
	hist_dict.emplace(name, tau_eta);
	hist_dict.emplace(name, tau_phi);
	hist_dict.emplace(name, muon_pt);
	hist_dict.emplace(name, muon_eta);
	hist_dict.emplace(name, muon_phi);
	hist_dict.emplace(name, met_pt);
	hist_dict.emplace(name, met_phi);
	hist_dict.emplace(name, MT);
	hist_dict.emplace(name, coll_mass);
	hist_dict.emplace(name, coll_mass_alt);
	hist_dict.emplace(name, nvtx);
	hist_dict.emplace(name, nvtxgood);
};

bool gen_match(	const rvec<int>& gen_pdgId,
				const rvec<float>& gen_part_pt,
				const rvec<float>& gen_part_eta,
				const rvec<float>& gen_part_phi,
				const rvec<float>& gen_part_mass,
				const int pdgId,
				const float& part_pt,
				const float& part_eta,
				const float& part_phi,
				const float& part_mass
				) {
	TLorentzVector p1, p2; 
	p1.SetPtEtaPhiM(part_pt, part_eta, part_phi, part_mass);
	for (uint i = 0; i < gen_pdgId.size(); i++) {
		if (gen_pdgId[i] == pdgId) {
			p2.SetPtEtaPhiM(gen_part_pt[i], gen_part_eta[i], gen_part_phi[i], gen_part_mass[i]);
			if (p1.DeltaR(p2) < 0.3) 
				return true;
		}
	}
	
	return false;
};

// Calculate datadriven distributions
RNode calc_datadriven(RNode df) {
	
	// begin evaluating the branch with non isolated taus for datadriven background estimation in QCD and ZToNuNu						 
	auto noniso_tau = df.Define("NonIso_Tau", tau_acceptance_and_id_and_dm_noniso,	{"Tau_pt_ES", "Tau_eta_ES", "Tau_dz", config::tau_dm, "Tau_decayMode", config::tau_iso, config::tau_antiEle, config::tau_antiMuon});
											
	auto iso_tau = df.Define("Iso_Tau", tau_acceptance_and_id_and_dm,	{"Tau_pt_ES", "Tau_eta_ES", "Tau_dz", config::tau_dm, "Tau_decayMode", config::tau_iso, config::tau_antiEle, config::tau_antiMuon});
											
	auto evntselect_noniso_tau = noniso_tau.Filter(nparticle_cut, {"NonIso_Tau"}, "select_taus")
                                .Filter(nparticle_cut, {"Muon_mask"}, "select_muons")
                                .Filter(nparticle_veto, {"Electron_mask"}, "select_eles");
	auto evntselect_iso_tau = iso_tau.Filter(nparticle_cut, {"Iso_Tau"}, "select_taus")
                                .Filter(nparticle_cut, {"Muon_mask"}, "select_muons")
                                .Filter(nparticle_veto, {"Electron_mask"}, "select_eles");
	
	// Now that we have our real tau, write it to special columns to handle it more easily
	auto defquants = evntselect_noniso_tau.Define("sel_Tau_pt", selected_part_quant, {"Tau_pt_ES", "NonIso_Tau"})
							   .Define("sel_Tau_eta", selected_part_quant, {"Tau_eta_ES", "NonIso_Tau"})
							   .Define("sel_Tau_phi", selected_part_quant, {"Tau_phi_ES", "NonIso_Tau"})
							   .Define("sel_Tau_jetIdx", selected_jet_quant, {"Tau_jetIdx", "NonIso_Tau"})
                               .Define("sel_Muon_pt", selected_part_quant, {"Muon_tP_pt", "Muon_mask"})
                               .Define("sel_Muon_eta", selected_part_quant, {"Muon_eta", "Muon_mask"})
                               .Define("sel_Muon_phi", selected_part_quant, {"Muon_phi", "Muon_mask"});

    auto dimuon_cut = defquants.Filter(dimuonpair_cut,{"Muon_tP_pt","Muon_eta","Muon_phi", "DiMuon_mask"},"Di Muon Pair Cut");
	// Calculate MT distribution
	auto mtcalc = dimuon_cut.Define("MT", mass_transv, {"sel_Muon_pt", "sel_Muon_phi", "MET_pt", "MET_phi"}); //.Define("sel_Tau_pt_ov_Jet_pt_NonIso", [](const rvec<float> tau_pt, 
                                                                                                             //                                  const rvec<Int_t> jetidx,
                                                                                                               //                                const rvec<float> jet_pt){
//        return tau_pt/jet_pt[jetidx];                                             
  //  },
    //                {"sel_Tau_pt_NonIso","sel_Tau_jetIdx_NonIso","Jet_pt"});
	RNode define_weights = mtcalc;
    // Define weights
    if (!config::runOnData) {
        auto defineScaleFactors = mtcalc.Define(   "TauEleFakeScaleFactor", 
                                                    []( const rvec<float>& tau_pt,
                                                        const rvec<float>& tau_eta,
                                                        const rvec<bool>& tau_mask,
                                                        const rvec<UChar_t>& tau_genPartFlav
                                                        )
                                                    {
                                                        return tau_fake_scale_factor(tau_pt, tau_eta, tau_mask, tau_genPartFlav, "Ele", "");
                                                    }, {"Tau_pt", "Tau_eta", "NonIso_Tau", "Tau_genPartFlav"})
                                        
                                        .Define(    "TauJetFakeScaleFactor", 
                                                    []( const rvec<float>& tau_pt,
                                                        const rvec<float>& tau_eta,
                                                        const rvec<bool>& tau_mask,
                                                        const rvec<UChar_t>& tau_genPartFlav
                                                        )
                                                    {
                                                        return tau_fake_scale_factor(tau_pt, tau_eta, tau_mask, tau_genPartFlav, "Jet", "");
                                                    }, {"Tau_pt", "Tau_eta", "NonIso_Tau", "Tau_genPartFlav"})
                                        
                                        .Define(    "TauMuonFakeScaleFactor", 
                                                    []( const rvec<float>& tau_pt,
                                                        const rvec<float>& tau_eta,
                                                        const rvec<bool>& tau_mask,
                                                        const rvec<UChar_t>& tau_genPartFlav
                                                        )
                                                    {
                                                        return tau_fake_scale_factor(tau_pt, tau_eta, tau_mask, tau_genPartFlav, "Muon", "");
                                                    }, {"Tau_pt", "Tau_eta", "NonIso_Tau", "Tau_genPartFlav"})
                                        
                                        .Define("PrefiringWeight", prefire_factor, {"Jet_pt", "Jet_eta", "Photon_pt", "Photon_eta"})
                                        .Define("MuonISOScaleFactor", [](const rvec<float>& muon_pt,
                                                                         const rvec<float>& muon_eta,
                                                                         const rvec<bool>& muon_mask
                                                                         )
                                                                         {
                                                                            return muon_iso_scale_factor(muon_pt, muon_eta, muon_mask, "");
                                                                         }, {"Muon_pt", "Muon_eta", "Muon_mask"})
                                        .Define("MuonIDScaleFactor", []( const rvec<float>& muon_pt,
                                                                         const rvec<float>& muon_eta,
                                                                         const rvec<bool>& muon_mask
                                                                         )
                                                                         {
                                                                            return muon_id_scale_factor(muon_pt, muon_eta, muon_mask, "");
                                                                         }, {"Muon_pt", "Muon_eta", "Muon_mask"});
        
        define_weights = defineScaleFactors;
    }
	RNode total_weights = define_weights; 
    if (config::runOnData) {
        total_weights = define_weights.Define("total_weight_fr", "1.0");
    } else {
        auto def_weights = [](const rvec<float>& weights){
            float mult = 1.0;
            for (const auto& it : weights)
                mult *= it;
            return mult;
        };
        total_weights = define_weights.Define("total_weight_fr", ROOT::RDF::PassAsVec<8, float>(def_weights), {
                                                                                                "pileup_weight",
                                                                                                "genWeight",
                                                                                                "TauJetFakeScaleFactor",
                                                                                                "TauEleFakeScaleFactor",
                                                                                                "TauMuonFakeScaleFactor",
                                                                                                "MuonIDScaleFactor",
                                                                                                "MuonISOScaleFactor",
                                                                                                "PrefiringWeight"});
    }
    auto df_col_idx = total_weights.Define("col_idx", col_idx, {"Tau_pt_ES", "Muon_pt", "Tau_eta_ES", "Muon_eta", "Tau_phi_ES", "Muon_phi", "Tau_mass_ES", "Muon_mass", "MET_pt", "MET_phi", "NonIso_Tau", "Muon_mask"});   
    auto delta_r_cut = df_col_idx.Filter(deltar_tau_muon_cut,{"Tau_eta_ES", "Tau_phi_ES", "Muon_eta","Muon_phi", "NonIso_Tau","Muon_mask", "col_idx"},"tau deltaR");
    auto df_tau_pt_ov_jet_pt = delta_r_cut.Define("Tau_pt_over_Jet_pt", tau_pt_over_jet_pt,{"col_idx","Tau_pt_ES","Tau_jetIdx","Jet_pt","NonIso_Tau"})
                                          .Define("sel_Tau_pt_alt",sel_tau_pt,{"Tau_pt_ES","NonIso_Tau","col_idx"});
    auto df_coll = df_tau_pt_ov_jet_pt.Define("CollMass", collinear_mass, {"col_idx","Tau_pt_ES", "Muon_pt", "Tau_eta_ES", "Muon_eta", "Tau_phi_ES", "Muon_phi", "Tau_mass_ES", "Muon_mass", "MET_pt", "MET_phi","NonIso_Tau", "Muon_mask"});  
    auto df_coll_alt = df_coll.Define("CollMass_alt", collinear_mass_alt, {"col_idx","Tau_pt_ES", "Muon_pt", "Tau_eta_ES", "Muon_eta", "Tau_phi_ES", "Muon_phi", "Tau_mass_ES", "Muon_mass", "MET_pt", "MET_phi","NonIso_Tau", "Muon_mask"});  
	// Now that we have our real tau, write it to special columns to handle it more easily
	auto defquants_iso = evntselect_iso_tau.Define("sel_Tau_pt", selected_part_quant, {"Tau_pt_ES", "Iso_Tau"})
							   .Define("sel_Tau_eta", selected_part_quant, {"Tau_eta_ES", "Iso_Tau"})
							   .Define("sel_Tau_phi", selected_part_quant, {"Tau_phi_ES", "Iso_Tau"})
							   .Define("sel_Tau_jetIdx", selected_jet_quant, {"Tau_jetIdx", "Iso_Tau"})
                               .Define("sel_Muon_pt", selected_part_quant, {"Muon_pt", "Muon_mask"})
                               .Define("sel_Muon_eta", selected_part_quant, {"Muon_eta", "Muon_mask"})
                               .Define("sel_Muon_phi", selected_part_quant, {"Muon_phi", "Muon_mask"});

    auto dimuon_cut_iso = defquants_iso.Filter(dimuonpair_cut,{"Muon_pt","Muon_eta","Muon_phi", "DiMuon_mask"},"Di Muon Pair Cut");
	// Calculate MT distribution
	auto mtcalc_iso = dimuon_cut_iso.Define("MT", mass_transv, {"sel_Muon_pt", "sel_Muon_phi", "MET_pt", "MET_phi"}); //.Define("sel_Tau_pt_ov_Jet_pt_NonIso", [](const rvec<float> tau_pt, 
                                                                                                             //                                  const rvec<Int_t> jetidx,
                                                                                                               //                                const rvec<float> jet_pt){
//        return tau_pt/jet_pt[jetidx];                                             
  //  },
    //                {"sel_Tau_pt_NonIso","sel_Tau_jetIdx_NonIso","Jet_pt"});
	RNode define_weights_iso = mtcalc_iso;
    if (!config::runOnData) {
        auto defineScaleFactors_iso = mtcalc_iso.Define(   "TauEleFakeScaleFactor", 
                                                    []( const rvec<float>& tau_pt,
                                                        const rvec<float>& tau_eta,
                                                        const rvec<bool>& tau_mask,
                                                        const rvec<UChar_t>& tau_genPartFlav
                                                        )
                                                    {
                                                        return tau_fake_scale_factor(tau_pt, tau_eta, tau_mask, tau_genPartFlav, "Ele", "");
                                                    }, {"Tau_pt", "Tau_eta", "Iso_Tau", "Tau_genPartFlav"})
                                        
                                        .Define(    "TauJetFakeScaleFactor", 
                                                    []( const rvec<float>& tau_pt,
                                                        const rvec<float>& tau_eta,
                                                        const rvec<bool>& tau_mask,
                                                        const rvec<UChar_t>& tau_genPartFlav
                                                        )
                                                    {
                                                        return tau_fake_scale_factor(tau_pt, tau_eta, tau_mask, tau_genPartFlav, "Jet", "");
                                                    }, {"Tau_pt", "Tau_eta", "Iso_Tau", "Tau_genPartFlav"})
                                        
                                        .Define(    "TauMuonFakeScaleFactor", 
                                                    []( const rvec<float>& tau_pt,
                                                        const rvec<float>& tau_eta,
                                                        const rvec<bool>& tau_mask,
                                                        const rvec<UChar_t>& tau_genPartFlav
                                                        )
                                                    {
                                                        return tau_fake_scale_factor(tau_pt, tau_eta, tau_mask, tau_genPartFlav, "Muon", "");
                                                    }, {"Tau_pt", "Tau_eta", "Iso_Tau", "Tau_genPartFlav"})
                                        
                                        .Define("PrefiringWeight", prefire_factor, {"Jet_pt", "Jet_eta", "Photon_pt", "Photon_eta"})
                                        .Define("MuonISOScaleFactor", [](const rvec<float>& muon_pt,
                                                                         const rvec<float>& muon_eta,
                                                                         const rvec<bool>& muon_mask
                                                                         )
                                                                         {
                                                                            return muon_iso_scale_factor(muon_pt, muon_eta, muon_mask, "");
                                                                         }, {"Muon_pt", "Muon_eta", "Muon_mask"})
                                        .Define("MuonIDScaleFactor", []( const rvec<float>& muon_pt,
                                                                         const rvec<float>& muon_eta,
                                                                         const rvec<bool>& muon_mask
                                                                         )
                                                                         {
                                                                            return muon_id_scale_factor(muon_pt, muon_eta, muon_mask, "");
                                                                         }, {"Muon_pt", "Muon_eta", "Muon_mask"});
        
        define_weights_iso = defineScaleFactors_iso;
    }
    // Define weights
	RNode total_weights_iso = define_weights_iso; 
    if (config::runOnData) {
        total_weights_iso = define_weights_iso.Define("total_weight_fr", "1.0");
    } else {
        auto def_weights = [](const rvec<float>& weights){
            float mult = 1.0;
            for (const auto& it : weights)
                mult *= it;
            return mult;
        };
        total_weights_iso = define_weights_iso.Define("total_weight_fr", ROOT::RDF::PassAsVec<8, float>(def_weights), {
                                                                                                "pileup_weight",
                                                                                                "genWeight",
                                                                                                "TauJetFakeScaleFactor",
                                                                                                "TauEleFakeScaleFactor",
                                                                                                "TauMuonFakeScaleFactor",
                                                                                                "MuonIDScaleFactor",
                                                                                                "MuonISOScaleFactor",
                                                                                                "PrefiringWeight"});
    }

    auto df_col_idx_iso = total_weights_iso.Define("col_idx", col_idx, {"Tau_pt_ES", "Muon_pt", "Tau_eta_ES", "Muon_eta", "Tau_phi_ES", "Muon_phi", "Tau_mass_ES", "Muon_mass", "MET_pt", "MET_phi", "Iso_Tau", "Muon_mask"});   
    auto delta_r_cut_iso = df_col_idx_iso.Filter(deltar_tau_muon_cut,{"Tau_eta_ES", "Tau_phi_ES", "Muon_eta","Muon_phi","Iso_Tau","Muon_mask", "col_idx"},"tau deltaR");
    auto df_tau_pt_ov_jet_pt_iso = df_col_idx_iso.Define("Tau_pt_over_Jet_pt", tau_pt_over_jet_pt,{"col_idx","Tau_pt_ES","Tau_jetIdx","Jet_pt","Iso_Tau"})
                                                  .Define("sel_Tau_pt_alt",sel_tau_pt,{"Tau_pt_ES","Iso_Tau","col_idx"});
    auto df_coll_iso = df_tau_pt_ov_jet_pt_iso.Define("CollMass", collinear_mass, {"col_idx","Tau_pt_ES", "Muon_pt", "Tau_eta_ES", "Muon_eta", "Tau_phi_ES", "Muon_phi", "Tau_mass_ES", "Muon_mass", "MET_pt", "MET_phi","Iso_Tau", "Muon_mask"});  
    auto df_coll_alt_iso = df_coll_iso.Define("CollMass_alt", collinear_mass_alt, {"col_idx","Tau_pt_ES", "Muon_pt", "Tau_eta_ES", "Muon_eta", "Tau_phi_ES", "Muon_phi", "Tau_mass_ES", "Muon_mass", "MET_pt", "MET_phi","Iso_Tau", "Muon_mask"});  
    // Define Stage 0: any event with one tau fulfilling acceptance and id	
	create_datadriven_hists(df_coll_alt, "NonIso", "total_weight_fr");
	create_datadriven_hists(df_coll_alt_iso, "Iso", "total_weight_fr");
	
	// actual analysis cut  -- MT > 120. GeV 
	auto df_mt = df_coll_alt.Filter("MT < 120.", "mt_cut");
	auto df_mt_iso = df_coll_alt_iso.Filter("MT < 120.", "mt_cut");
	auto df_mt_signal = df_coll_alt.Filter("MT > 120.", "mt_cut");
	auto df_mt_signal_iso = df_coll_alt_iso.Filter("MT > 120.", "mt_cut");
	create_datadriven_hists(df_mt, "NonIso_FakeRegion", "total_weight_fr");
	create_datadriven_hists(df_mt_iso, "Iso_FakeRegion", "total_weight_fr");
	create_datadriven_hists(df_mt_signal, "NonIso_SignalRegion", "total_weight_fr");
	create_datadriven_hists(df_mt_signal_iso, "Iso_SignalRegion", "total_weight_fr");
    if (not config::runOnData){
        auto df_mc_tau_veto_noniso = df_mt.Filter([](const rvec<int> col_idx, const rvec<UChar_t> genpartflav, const rvec<bool> tau_mask){
            if (genpartflav[tau_mask][col_idx[0]] == 5) return false;
            else return true;}, {"col_idx", "Tau_genPartFlav", "NonIso_Tau"},"Data Driven MC real Tau veto");
        auto df_mc_tau_veto_iso = df_mt_iso.Filter([](const rvec<int> col_idx, const rvec<UChar_t> genpartflav, const rvec<bool> tau_mask){
            if (genpartflav[tau_mask][col_idx[0]] == 5) return false;
            else return true;}, {"col_idx", "Tau_genPartFlav", "Iso_Tau"},"Data Driven MC real Tau veto");
        auto df_mc_tau_veto_signal_noniso = df_mt_signal.Filter([](const rvec<int> col_idx, const rvec<UChar_t> genpartflav, const rvec<bool> tau_mask){
            if (genpartflav[tau_mask][col_idx[0]] == 5) return false;
            else return true;}, {"col_idx", "Tau_genPartFlav", "NonIso_Tau"},"Data Driven MC real Tau veto");
        auto df_mc_tau_veto_signal_iso = df_mt_signal_iso.Filter([](const rvec<int> col_idx, const rvec<UChar_t> genpartflav, const rvec<bool> tau_mask){
            if (genpartflav[tau_mask][col_idx[0]] == 5) return false;
            else return true;}, {"col_idx", "Tau_genPartFlav", "Iso_Tau"},"Data Driven MC real Tau veto");
	    create_datadriven_hists(df_mc_tau_veto_noniso, "NonIso_FakeRegion_MC", "total_weight_fr");
	    create_datadriven_hists(df_mc_tau_veto_iso, "Iso_FakeRegion_MC", "total_weight_fr");
	    create_datadriven_hists(df_mc_tau_veto_signal_noniso, "NonIso_SignalRegion_MC", "total_weight_fr");
	    create_datadriven_hists(df_mc_tau_veto_signal_iso, "Iso_SignalRegion_MC", "total_weight_fr");

    }
    return df_mt;
//
//    for (const auto& it : list_of_weights) {
//        std::string folder_name = "Stage1";
//        if ( it.find("total_weight_") != std::string::npos )
//            folder_name = "Stage1/Systematics";
//        if ( it.find("pdf_weight") != std::string::npos )
//            folder_name = "Stage1/PDF";
//        create_datadriven_hists(df_mt, folder_name, it);
//    }
//    return df_mt;
};

// Apply datadriven distributions
RNode apply_datadriven(RNode df) {
	
	// begin evaluating the branch with non isolated taus for datadriven background estimation in QCD and ZToNuNu						 
	auto noniso_tau = df.Define("NonIso_Tau", tau_acceptance_and_id_and_dm_noniso,	{"Tau_pt_ES", "Tau_eta_ES", "Tau_dz", config::tau_dm, "Tau_decayMode", config::tau_iso, config::tau_antiEle, config::tau_antiMuon});
											
	auto iso_tau = df.Define("Iso_Tau", tau_acceptance_and_id_and_dm,	{"Tau_pt_ES", "Tau_eta_ES", "Tau_dz", config::tau_dm, "Tau_decayMode", config::tau_iso, config::tau_antiEle, config::tau_antiMuon});
											
	auto evntselect_noniso_tau = noniso_tau.Filter(nparticle_cut, {"NonIso_Tau"}, "select_taus")
                                .Filter(nparticle_cut, {"Muon_mask"}, "select_muons")
                                .Filter(nparticle_veto, {"Electron_mask"}, "select_eles");
	auto evntselect_iso_tau = iso_tau.Filter(nparticle_cut, {"Iso_Tau"}, "select_taus")
                                .Filter(nparticle_cut, {"Muon_mask"}, "select_muons")
                                .Filter(nparticle_veto, {"Electron_mask"}, "select_eles");
	
	// Now that we have our real tau, write it to special columns to handle it more easily
	auto defquants_iso = evntselect_iso_tau.Define("sel_Tau_pt", selected_part_quant, {"Tau_pt_ES", "Iso_Tau"})
							   .Define("sel_Tau_eta", selected_part_quant, {"Tau_eta_ES", "Iso_Tau"})
							   .Define("sel_Tau_phi", selected_part_quant, {"Tau_phi_ES", "Iso_Tau"})
							   .Define("sel_Tau_jetIdx", selected_jet_quant, {"Tau_jetIdx", "Iso_Tau"})
                               .Define("sel_Muon_pt", selected_part_quant, {"Muon_pt", "Muon_mask"})
                               .Define("sel_Muon_eta", selected_part_quant, {"Muon_eta", "Muon_mask"})
                               .Define("sel_Muon_phi", selected_part_quant, {"Muon_phi", "Muon_mask"});

    auto dimuon_cut_iso = defquants_iso.Filter(dimuonpair_cut,{"Muon_pt","Muon_eta","Muon_phi", "DiMuon_mask"},"Di Muon Pair Cut");
	// Calculate MT distribution
	auto mtcalc_iso = dimuon_cut_iso.Define("MT", mass_transv, {"sel_Muon_pt", "sel_Muon_phi", "MET_pt", "MET_phi"}); //.Define("sel_Tau_pt_ov_Jet_pt_NonIso", [](const rvec<float> tau_pt, 
                                                                                                             //                                  const rvec<Int_t> jetidx,
                                                                                                               //                                const rvec<float> jet_pt){
//        return tau_pt/jet_pt[jetidx];                                             
  //  },
    //                {"sel_Tau_pt_NonIso","sel_Tau_jetIdx_NonIso","Jet_pt"});

    auto df_col_idx_iso = mtcalc_iso.Define("col_idx", col_idx, {"Tau_pt_ES", "Muon_pt", "Tau_eta_ES", "Muon_eta", "Tau_phi_ES", "Muon_phi", "Tau_mass_ES", "Muon_mass", "MET_pt", "MET_phi", "Iso_Tau", "Muon_mask"});   
    auto delta_r_cut_iso = df_col_idx_iso.Filter(deltar_tau_muon_cut,{"Tau_eta_ES", "Tau_phi_ES", "Muon_eta","Muon_phi","Iso_Tau","Muon_mask", "col_idx"},"tau deltaR");
    auto df_tau_pt_ov_jet_pt_iso = df_col_idx_iso.Define("Tau_pt_over_Jet_pt", tau_pt_over_jet_pt,{"col_idx","Tau_pt_ES","Tau_jetIdx","Jet_pt","Iso_Tau"})
                                                  .Define("sel_Tau_pt_alt",sel_tau_pt,{"Tau_pt_ES","Iso_Tau","col_idx"});
    auto df_coll_iso = df_tau_pt_ov_jet_pt_iso.Define("CollMass", collinear_mass, {"col_idx","Tau_pt_ES", "Muon_pt", "Tau_eta_ES", "Muon_eta", "Tau_phi_ES", "Muon_phi", "Tau_mass_ES", "Muon_mass", "MET_pt", "MET_phi","Iso_Tau", "Muon_mask"});  
    auto df_coll_alt_iso = df_coll_iso.Define("CollMass_alt", collinear_mass_alt, {"col_idx","Tau_pt_ES", "Muon_pt", "Tau_eta_ES", "Muon_eta", "Tau_phi_ES", "Muon_phi", "Tau_mass_ES", "Muon_mass", "MET_pt", "MET_phi","Iso_Tau", "Muon_mask"});  
    auto define_weights_iso = df_coll_alt_iso;
    if (config::runOnData){
        define_weights_iso  = df_coll_alt_iso.Define("FakeRate", 
                                                    []( const rvec<float>& tau_pt,
                                                           const float & taupt_o_jetpt,
                                                           const rvec<int>& col_idx, 
                                                           const rvec<bool>& tau_mask
                                                        )
                                                    {
                                                        return dd_fakerate(tau_pt, taupt_o_jetpt, col_idx, tau_mask, "");
                                                    }, {"Tau_pt_ES","Tau_pt_over_Jet_pt","col_idx","Iso_Tau"});
    }
    else
    {
        define_weights_iso = df_coll_alt_iso.Define("FakeRate", 
                                                    []( const rvec<float>& tau_pt,
                                                           const float & taupt_o_jetpt,
                                                           const rvec<int>& col_idx, 
                                                           const rvec<bool>& tau_mask
                                                        )
                                                    {
                                                        return dd_fakerate(tau_pt, taupt_o_jetpt, col_idx, tau_mask, "");
                                                    }, {"Tau_pt_ES","Tau_pt_over_Jet_pt","col_idx","Iso_Tau"})
                                                    .Define(   "TauEleFakeScaleFactor", 
                                                    []( const rvec<float>& tau_pt,
                                                        const rvec<float>& tau_eta,
                                                        const rvec<bool>& tau_mask,
                                                        const rvec<UChar_t>& tau_genPartFlav
                                                        )
                                                    {
                                                        return tau_fake_scale_factor(tau_pt, tau_eta, tau_mask, tau_genPartFlav, "Ele", "");
                                                    }, {"Tau_pt", "Tau_eta", "Iso_Tau", "Tau_genPartFlav"})
                                        
                                        .Define(    "TauJetFakeScaleFactor", 
                                                    []( const rvec<float>& tau_pt,
                                                        const rvec<float>& tau_eta,
                                                        const rvec<bool>& tau_mask,
                                                        const rvec<UChar_t>& tau_genPartFlav
                                                        )
                                                    {
                                                        return tau_fake_scale_factor(tau_pt, tau_eta, tau_mask, tau_genPartFlav, "Jet", "");
                                                    }, {"Tau_pt", "Tau_eta", "Iso_Tau", "Tau_genPartFlav"})
                                        
                                        .Define(    "TauMuonFakeScaleFactor", 
                                                    []( const rvec<float>& tau_pt,
                                                        const rvec<float>& tau_eta,
                                                        const rvec<bool>& tau_mask,
                                                        const rvec<UChar_t>& tau_genPartFlav
                                                        )
                                                    {
                                                        return tau_fake_scale_factor(tau_pt, tau_eta, tau_mask, tau_genPartFlav, "Muon", "");
                                                    }, {"Tau_pt", "Tau_eta", "Iso_Tau", "Tau_genPartFlav"})
                                        
                                        .Define("PrefiringWeight", prefire_factor, {"Jet_pt", "Jet_eta", "Photon_pt", "Photon_eta"})
                                        .Define("MuonISOScaleFactor", [](const rvec<float>& muon_pt,
                                                                         const rvec<float>& muon_eta,
                                                                         const rvec<bool>& muon_mask
                                                                         )
                                                                         {
                                                                            return muon_iso_scale_factor(muon_pt, muon_eta, muon_mask, "");
                                                                         }, {"Muon_pt", "Muon_eta", "Muon_mask"})
                                        .Define("MuonIDScaleFactor", []( const rvec<float>& muon_pt,
                                                                         const rvec<float>& muon_eta,
                                                                         const rvec<bool>& muon_mask
                                                                         )
                                                                         {
                                                                            return muon_id_scale_factor(muon_pt, muon_eta, muon_mask, "");
                                                                         }, {"Muon_pt", "Muon_eta", "Muon_mask"});
        
    }
    //                                
    //
    // Define weights
	RNode total_weights_iso = define_weights_iso; 
    if (config::runOnData) {
        total_weights_iso = define_weights_iso.Define("total_weight", "FakeRate");
    } else {
        auto def_weights = [](const rvec<float>& weights){
            float mult = 1.0;
            for (const auto& it : weights)
                mult *= it;
            return mult;
        };
        total_weights_iso = define_weights_iso.Define("total_weight", ROOT::RDF::PassAsVec<9, float>(def_weights), {
                                                                                                "FakeRate",
                                                                                                "pileup_weight",
                                                                                                "genWeight",
                                                                                                "TauJetFakeScaleFactor",
                                                                                                "TauEleFakeScaleFactor",
                                                                                                "TauMuonFakeScaleFactor",
                                                                                                "MuonIDScaleFactor",
                                                                                                "MuonISOScaleFactor",
                                                                                                "PrefiringWeight"});
    }
    // Define Stage 0: any event with one tau fulfilling acceptance and id	
	create_fr_hists(total_weights_iso, "Iso", "total_weight");
	
	// actual analysis cut  -- MT > 120. GeV 
	auto df_mt = total_weights_iso.Filter("MT < 120.", "mt_cut");
	//auto df_mt_iso = total_weights_iso.Filter("MT < 120.", "mt_cut");
	//auto df_mt_signal = total_weights_iso.Filter("MT > 120.", "mt_cut");
	//auto df_mt_signal_iso = total_weights_iso.Filter("MT > 120.", "mt_cut");
	//create_fr_hists(df_mt_iso, "Iso_FakeRegion", "FakeRate");
	//create_fr_hists(df_mt_signal_iso, "Iso_SignalRegion", "FakeRate");
    return df_mt;
//
//    for (const auto& it : list_of_weights) {
//        std::string folder_name = "Stage1";
//        if ( it.find("total_weight_") != std::string::npos )
//            folder_name = "Stage1/Systematics";
//        if ( it.find("pdf_weight") != std::string::npos )
//            folder_name = "Stage1/PDF";
//        create_datadriven_hists(df_mt, folder_name, it);
//    }
//    return df_mt;
};

// analyse function - gets called for each systematic
void analyse(	RNode df, 
				TFile* outFile) {	
	
	// init counter
	auto genWeightSum = df.Sum("genWeight");
    auto eventcounter = df.Count();
    
    if (!config::runOnData)
        df = init_PDFs(df);
	
	// Select trigger requirements
	auto triggered = trigger(df);
								
	//auto met_filter = triggered.Filter("Flag_METFilters", "MET filters");
    auto met_filter = triggered.Filter(config::metfilters , "MET filters");

		
	// Trigger turn on cut
	//auto trigger_obj1 = met_filter.Filter("Muon_pt > " + std::to_string(config::muon_pt), "Avoid trigger turn on with muon object, pT>" + std::to_string(config::muon_pt));
	

	//auto good_primary_vertex = triggered.Filter([](	const float& pvx, 
	auto good_primary_vertex = met_filter.Filter([](	const float& pvx, 
								const float& pvy, 
								const float& pvz,
								const float& pvndof ){
									if (std::abs(pvz) > config::pv_z)
										return false;
									if (pow(pvx, 2) + pow(pvy, 2) > config::pv_d)
										return false;
									if (pvndof < config::pv_ndof)
										return false;
									return true;
						}, {"PV_x", "PV_y", "PV_z", "PV_ndof"}, "Ensure good primary vertex selection");
	
	RNode fixed = good_primary_vertex;
	// handle 2018 HEM problem
	/* if (config::era == 2018) {
		fixed = trigger_obj3.Filter([]         (const rvec<float>& jet_pt,
												const rvec<float>& jet_eta,
												const rvec<float>& jet_phi){
													for (uint i = 0; i < jet_pt.size(); i++) {
														if (jet_pt[i] < 30)
															continue;
														if ( jet_eta[i] < -3.0 || jet_eta[i] > -1.3 )
															continue;
														if ( jet_phi[i] < -1.57 || jet_phi[i] > -0.87 )
															continue;
														return false;
												    }
													return true;
												}, {"Jet_pt", "Jet_eta", "Jet_phi"}, "2018 HEM problem");
	} */
	
	// fill preselection
	if (config::run_type == "") {
		fill_preselection(good_primary_vertex);
	}
	
	// creates mask, which fills bool tags for taus which fulfil id and are in acceptance
	auto masked = good_primary_vertex.Define("Tau_mask", tau_acceptance_and_id_and_dm,	{"Tau_pt_ES", "Tau_eta_ES", "Tau_dz", config::tau_dm, "Tau_decayMode", config::tau_iso, config::tau_antiEle, config::tau_antiMuon})
							  .Define("Muon_mask", muon_acceptance_and_id, {"Muon_tP_pt", "Muon_eta", "Muon_highPtId", "Muon_tkIsoId"})
							  .Define("DiMuon_mask", di_muon_id, {"Muon_tP_pt", "Muon_eta", "Muon_highPtId", "Muon_tkRelIso"})
							  .Define("Electron_mask", ele_acceptance_and_simpleid, {"Electron_pt", "Electron_eta", "Electron_cutBased_HEEP"});
	
	
	
	// this is for datadriven part of analysis
    auto df_final = [&]() -> RNode{
        if (config::calcDataDriven){
            return calc_datadriven(masked);
        }
        if (config::runDataDriven ){
            return apply_datadriven(masked);
        }
	
        
        // Select events with certain number of taus, which fulfil all acceptance & id
        // also, cut any event with electrons or muons
        auto evntselect = masked.Filter(nparticle_cut, {"Tau_mask"}, "select_taus")
                                .Filter(nparticle_cut, {"Muon_mask"}, "select_muons")
                                .Filter(nparticle_veto, {"Electron_mask"}, "select_eles");
        
        
        auto df_col_idx = evntselect.Define("col_idx", col_idx, {"Tau_pt_ES", "Muon_tP_pt", "Tau_eta_ES", "Muon_eta", "Tau_phi_ES", "Muon_phi", "Tau_mass_ES", "Muon_mass", "MET_pt", "MET_phi", "Tau_mask", "Muon_mask"});   
        // Now that we have our real tau, write it to special columns to handle it more easily
        auto defquants = df_col_idx.Define("sel_Tau_pt", selected_part_col_idx_tau, {"Tau_pt_ES", "Tau_mask", "col_idx"})
                                   .Define("sel_Tau_eta", selected_part_col_idx_tau, {"Tau_eta_ES", "Tau_mask", "col_idx"})
                                   .Define("sel_Tau_phi", selected_part_col_idx_tau, {"Tau_phi_ES", "Tau_mask", "col_idx"})
                                   .Define("sel_Muon_pt", selected_part_col_idx_muon, {"Muon_tP_pt", "Muon_mask", "col_idx"})
                                   .Define("sel_Muon_eta", selected_part_col_idx_muon, {"Muon_eta", "Muon_mask", "col_idx"})
                                   .Define("sel_Muon_phi", selected_part_col_idx_muon, {"Muon_phi", "Muon_mask", "col_idx"});
        
        //auto dimuon_cut = defquants.Filter(dimuonpair_cut,{"Muon_pt","Muon_eta","Muon_phi", "DiMuon_mask"},"Di Muon Pair Cut").Define("sel_tau_genPartFlav_int",[](const rvec<UChar_t>& genPartFlav){ 
        //    rvec<float> copy;
        //    for (const auto it: genPartFlav){
        //        copy.push_back((float) it);
        //    }
        //return copy;},{"Tau_genPartFlav"});
        
        auto dimuon_cut = defquants.Filter(dimuonpair_cut,{"Muon_tP_pt","Muon_eta","Muon_phi", "DiMuon_mask"},"Di Muon Pair Cut");
        
        // Calculate MT distribution
        //auto mtcalc = dimuon_cut.Define("MT", mass_transv, {"sel_Muon_pt", "sel_Muon_phi", "MET_pt", "MET_phi"}).Define("sel_Tau_genPartFlav", selected_part_quant, {"sel_tau_genPartFlav_int", "Tau_mask"});
        auto mtcalc = dimuon_cut.Define("MT", mass_transv, {"sel_Muon_pt", "sel_Muon_phi", "MET_pt", "MET_phi"});
        
        RNode define_weights = mtcalc;
        // Define and calculate weights for monte carlo
        if (!config::runOnData) {
            auto defineScaleFactors = mtcalc.Define(   "TauEleFakeScaleFactor", 
                                                        []( const rvec<float>& tau_pt,
                                                            const rvec<float>& tau_eta,
                                                            const rvec<bool>& tau_mask,
                                                            const rvec<UChar_t>& tau_genPartFlav
                                                            )
                                                        {
                                                            return tau_fake_scale_factor(tau_pt, tau_eta, tau_mask, tau_genPartFlav, "Ele", "");
                                                        }, {"Tau_pt", "Tau_eta", "Tau_mask", "Tau_genPartFlav"})
                                            
                                            .Define(    "TauEleFakeScaleFactorUp", 
                                                        []( const rvec<float>& tau_pt,
                                                            const rvec<float>& tau_eta,
                                                            const rvec<bool>& tau_mask,
                                                            const rvec<UChar_t>& tau_genPartFlav
                                                            )
                                                        {
                                                            return tau_fake_scale_factor(tau_pt, tau_eta, tau_mask, tau_genPartFlav, "Ele", "Up");
                                                        }, {"Tau_pt", "Tau_eta", "Tau_mask", "Tau_genPartFlav"})
                                            
                                            .Define(    "TauEleFakeScaleFactorDown", 
                                                        []( const rvec<float>& tau_pt,
                                                            const rvec<float>& tau_eta,
                                                            const rvec<bool>& tau_mask,
                                                            const rvec<UChar_t>& tau_genPartFlav
                                                            )
                                                        {
                                                            return tau_fake_scale_factor(tau_pt, tau_eta, tau_mask, tau_genPartFlav, "Ele", "Down");
                                                        }, {"Tau_pt", "Tau_eta", "Tau_mask", "Tau_genPartFlav"})
                                            
                                            .Define(    "TauJetFakeScaleFactor", 
                                                        []( const rvec<float>& tau_pt,
                                                            const rvec<float>& tau_eta,
                                                            const rvec<bool>& tau_mask,
                                                            const rvec<UChar_t>& tau_genPartFlav
                                                            )
                                                        {
                                                            return tau_fake_scale_factor(tau_pt, tau_eta, tau_mask, tau_genPartFlav, "Jet", "");
                                                        }, {"Tau_pt", "Tau_eta", "Tau_mask", "Tau_genPartFlav"})
                                            
                                            .Define(    "TauJetFakeScaleFactorUp", 
                                                        []( const rvec<float>& tau_pt,
                                                            const rvec<float>& tau_eta,
                                                            const rvec<bool>& tau_mask,
                                                            const rvec<UChar_t>& tau_genPartFlav
                                                            )
                                                        {
                                                            return tau_fake_scale_factor(tau_pt, tau_eta, tau_mask, tau_genPartFlav, "Jet", "Up");
                                                        }, {"Tau_pt", "Tau_eta", "Tau_mask", "Tau_genPartFlav"})
                                            
                                            .Define(    "TauJetFakeScaleFactorDown", 
                                                        []( const rvec<float>& tau_pt,
                                                            const rvec<float>& tau_eta,
                                                            const rvec<bool>& tau_mask,
                                                            const rvec<UChar_t>& tau_genPartFlav
                                                            )
                                                        {
                                                            return tau_fake_scale_factor(tau_pt, tau_eta, tau_mask, tau_genPartFlav, "Jet", "Down");
                                                        }, {"Tau_pt", "Tau_eta", "Tau_mask", "Tau_genPartFlav"})
                                            
                                            .Define(    "TauMuonFakeScaleFactor", 
                                                        []( const rvec<float>& tau_pt,
                                                            const rvec<float>& tau_eta,
                                                            const rvec<bool>& tau_mask,
                                                            const rvec<UChar_t>& tau_genPartFlav
                                                            )
                                                        {
                                                            return tau_fake_scale_factor(tau_pt, tau_eta, tau_mask, tau_genPartFlav, "Muon", "");
                                                        }, {"Tau_pt", "Tau_eta", "Tau_mask", "Tau_genPartFlav"})
                                            
                                            .Define(    "TauMuonFakeScaleFactorUp", 
                                                        []( const rvec<float>& tau_pt,
                                                            const rvec<float>& tau_eta,
                                                            const rvec<bool>& tau_mask,
                                                            const rvec<UChar_t>& tau_genPartFlav
                                                            )
                                                        {
                                                            return tau_fake_scale_factor(tau_pt, tau_eta, tau_mask, tau_genPartFlav, "Muon", "Up");
                                                        }, {"Tau_pt", "Tau_eta", "Tau_mask", "Tau_genPartFlav"})
                                            
                                            .Define(    "TauMuonFakeScaleFactorDown", 
                                                        []( const rvec<float>& tau_pt,
                                                            const rvec<float>& tau_eta,
                                                            const rvec<bool>& tau_mask,
                                                            const rvec<UChar_t>& tau_genPartFlav
                                                            )
                                                        {
                                                            return tau_fake_scale_factor(tau_pt, tau_eta, tau_mask, tau_genPartFlav, "Muon", "Down");
                                                        }, {"Tau_pt", "Tau_eta", "Tau_mask", "Tau_genPartFlav"})
                                            

                                            .Define("PrefiringWeight", prefire_factor, {"Jet_pt", "Jet_eta", "Photon_pt", "Photon_eta"})
                                            .Define("MuonISOScaleFactor", [](const rvec<float>& muon_pt,
                                                                             const rvec<float>& muon_eta,
                                                                             const rvec<bool>& muon_mask
                                                                             )
                                                                             {
                                                                                return muon_iso_scale_factor(muon_pt, muon_eta, muon_mask, "");
                                                                             }, {"Muon_tP_pt", "Muon_eta", "Muon_mask"})
                                            .Define("MuonISOScaleFactorUp", [](const rvec<float>& muon_pt,
                                                                             const rvec<float>& muon_eta,
                                                                             const rvec<bool>& muon_mask
                                                                             )
                                                                             {
                                                                                return muon_iso_scale_factor(muon_pt, muon_eta, muon_mask, "Up");
                                                                             }, {"Muon_tP_pt", "Muon_eta", "Muon_mask"})
                                            .Define("MuonISOScaleFactorDown", [](const rvec<float>& muon_pt,
                                                                             const rvec<float>& muon_eta,
                                                                             const rvec<bool>& muon_mask
                                                                             )
                                                                             {
                                                                                return muon_iso_scale_factor(muon_pt, muon_eta, muon_mask, "Down");
                                                                             }, {"Muon_tP_pt", "Muon_eta", "Muon_mask"})
                                            .Define("MuonIDScaleFactor", []( const rvec<float>& muon_pt,
                                                                             const rvec<float>& muon_eta,
                                                                             const rvec<bool>& muon_mask
                                                                             )
                                                                             {
                                                                                return muon_id_scale_factor(muon_pt, muon_eta, muon_mask, "");
                                                                             }, {"Muon_tP_pt", "Muon_eta", "Muon_mask"})
            
                                            .Define("MuonIDScaleFactorUp", []( const rvec<float>& muon_pt,
                                                                             const rvec<float>& muon_eta,
                                                                             const rvec<bool>& muon_mask
                                                                             )
                                                                             {
                                                                                return muon_id_scale_factor(muon_pt, muon_eta, muon_mask, "Up");
                                                                             }, {"Muon_tP_pt", "Muon_eta", "Muon_mask"})
                                            .Define("MuonIDScaleFactorDown", []( const rvec<float>& muon_pt,
                                                                             const rvec<float>& muon_eta,
                                                                             const rvec<bool>& muon_mask
                                                                             )
                                                                             {
                                                                                return muon_id_scale_factor(muon_pt, muon_eta, muon_mask, "Down");
                                                                             }, {"Muon_tP_pt", "Muon_eta", "Muon_mask"});
            define_weights = defineScaleFactors;
        }
        

        // Define weights
        RNode total_weights = define_weights; 

        std::vector < std::string > list_of_weights;
        if (config::runOnData) {
            total_weights = define_weights.Define("total_weight", "1.0");
            list_of_weights.push_back("total_weight");
        } else {
            auto def_weights = [](const rvec<float>& weights){
                float mult = 1.0;
                for (const auto& it : weights)
                    mult *= it;
                return mult;
            };
            total_weights = define_weights.Define("total_weight", ROOT::RDF::PassAsVec<8, float>(def_weights), {
                                                                                                    "pileup_weight",
                                                                                                    "genWeight",
                                                                                                    "TauJetFakeScaleFactor",
                                                                                                    "TauEleFakeScaleFactor",
                                                                                                    "TauMuonFakeScaleFactor",
                                                                                                    "MuonIDScaleFactor",
                                                                                                    "MuonISOScaleFactor",
                                                                                                    "PrefiringWeight"});
            total_weights = total_weights.Define("total_weight_MuonISOScaleFactorUp", ROOT::RDF::PassAsVec<8, float>(def_weights), {
                                                                                                    "pileup_weight",
                                                                                                    "genWeight",
                                                                                                    "TauJetFakeScaleFactor",
                                                                                                    "TauEleFakeScaleFactor",
                                                                                                    "TauMuonFakeScaleFactor",
                                                                                                    "MuonIDScaleFactor",
                                                                                                    "MuonISOScaleFactorUp",
                                                                                                    "PrefiringWeight"});
            total_weights = total_weights.Define("total_weight_MuonISOScaleFactorDown", ROOT::RDF::PassAsVec<8, float>(def_weights), {
                                                                                                    "pileup_weight",
                                                                                                    "genWeight",
                                                                                                    "TauJetFakeScaleFactor",
                                                                                                    "TauEleFakeScaleFactor",
                                                                                                    "TauMuonFakeScaleFactor",
                                                                                                    "MuonIDScaleFactor",
                                                                                                    "MuonISOScaleFactorDown",
                                                                                                    "PrefiringWeight"});
            total_weights = total_weights.Define("total_weight_MuonIDScaleFactorUp", ROOT::RDF::PassAsVec<8, float>(def_weights), {
                                                                                                    "pileup_weight",
                                                                                                    "genWeight",
                                                                                                    "TauJetFakeScaleFactor",
                                                                                                    "TauEleFakeScaleFactor",
                                                                                                    "TauMuonFakeScaleFactor",
                                                                                                    "MuonIDScaleFactorUp",
                                                                                                    "MuonISOScaleFactor",
                                                                                                    "PrefiringWeight"});
            total_weights = total_weights.Define("total_weight_MuonIDScaleFactorDown", ROOT::RDF::PassAsVec<8, float>(def_weights), {
                                                                                                    "pileup_weight",
                                                                                                    "genWeight",
                                                                                                    "TauJetFakeScaleFactor",
                                                                                                    "TauEleFakeScaleFactor",
                                                                                                    "TauMuonFakeScaleFactor",
                                                                                                    "MuonIDScaleFactorDown",
                                                                                                    "MuonISOScaleFactor",
                                                                                                    "PrefiringWeight"});
            total_weights = total_weights.Define("total_weight_TauJetFakeScaleFactorUp", ROOT::RDF::PassAsVec<8, float>(def_weights), {
                                                                                                    "pileup_weight",
                                                                                                    "genWeight",
                                                                                                    "TauJetFakeScaleFactorUp",
                                                                                                    "TauEleFakeScaleFactor",
                                                                                                    "TauMuonFakeScaleFactor",
                                                                                                    "MuonIDScaleFactor",
                                                                                                    "MuonISOScaleFactor",
                                                                                                    "PrefiringWeight"});
            total_weights = total_weights.Define("total_weight_TauJetFakeScaleFactorDown", ROOT::RDF::PassAsVec<8, float>(def_weights), {
                                                                                                    "pileup_weight",
                                                                                                    "genWeight",
                                                                                                    "TauJetFakeScaleFactorDown",
                                                                                                    "TauEleFakeScaleFactor",
                                                                                                    "TauMuonFakeScaleFactor",
                                                                                                    "MuonIDScaleFactor",
                                                                                                    "MuonISOScaleFactor",
                                                                                                    "PrefiringWeight"});
            total_weights = total_weights.Define("total_weight_TauEleFakeScaleFactorUp", ROOT::RDF::PassAsVec<8, float>(def_weights), {
                                                                                                    "pileup_weight",
                                                                                                    "genWeight",
                                                                                                    "TauJetFakeScaleFactor",
                                                                                                    "TauEleFakeScaleFactorUp",
                                                                                                    "TauMuonFakeScaleFactor",
                                                                                                    "MuonIDScaleFactor",
                                                                                                    "MuonISOScaleFactor",
                                                                                                    "PrefiringWeight"});
            total_weights = total_weights.Define("total_weight_TauEleFakeScaleFactorDown", ROOT::RDF::PassAsVec<8, float>(def_weights), {
                                                                                                    "pileup_weight",
                                                                                                    "genWeight",
                                                                                                    "TauJetFakeScaleFactor",
                                                                                                    "TauEleFakeScaleFactorDown",
                                                                                                    "TauMuonFakeScaleFactor",
                                                                                                    "MuonIDScaleFactor",
                                                                                                    "MuonISOScaleFactor",
                                                                                                    "PrefiringWeight"});
            total_weights = total_weights.Define("total_weight_TauMuonFakeScaleFactorUp", ROOT::RDF::PassAsVec<8, float>(def_weights), {
                                                                                                    "pileup_weight",
                                                                                                    "TauJetFakeScaleFactor",
                                                                                                    "genWeight",
                                                                                                    "TauEleFakeScaleFactor",
                                                                                                    "TauMuonFakeScaleFactorUp",
                                                                                                    "MuonIDScaleFactor",
                                                                                                    "MuonISOScaleFactor",
                                                                                                    "PrefiringWeight"});
            total_weights = total_weights.Define("total_weight_TauMuonFakeScaleFactorDown", ROOT::RDF::PassAsVec<8, float>(def_weights), {
                                                                                                    "pileup_weight",
                                                                                                    "TauJetFakeScaleFactor",
                                                                                                    "genWeight",
                                                                                                    "TauEleFakeScaleFactor",
                                                                                                    "TauMuonFakeScaleFactorDown",
                                                                                                    "MuonIDScaleFactor",
                                                                                                    "MuonISOScaleFactor",
                                                                                                    "PrefiringWeight"});
            //~ total_weights = total_weights.Define("total_weight_pileupUp", ROOT::RDF::PassAsVec<8, float>(def_weights), {
                                                                                                    //~ "pileup_weight_Up",
                                                                                                    //~ "TauJetFakeScaleFactor",
                                                                                                    //~ "genWeight",
                                                                                                    //~ "W_kfactor",
                                                                                                    //~ "TauEleFakeScaleFactor",
                                                                                                    //~ "TauMuonFakeScaleFactor",
                                                                                                    //~ "PrefiringWeight"});
            //~ total_weights = total_weights.Define("total_weight_pileupDown", ROOT::RDF::PassAsVec<8, float>(def_weights), {
                                                                                                    //~ "pileup_weight_Down",
                                                                                                    //~ "TauJetFakeScaleFactor",
                                                                                                    //~ "genWeight",
                                                                                                    //~ "W_kfactor",
                                                                                                    //~ "TauEleFakeScaleFactor",
                                                                                                    //~ "TauMuonFakeScaleFactor",
                                                                                                    //~ "PrefiringWeight"});
        
            list_of_weights.push_back("total_weight");
            list_of_weights.push_back("total_weight_TauJetFakeScaleFactorUp");
            list_of_weights.push_back("total_weight_TauJetFakeScaleFactorDown");
            list_of_weights.push_back("total_weight_TauEleFakeScaleFactorUp");
            list_of_weights.push_back("total_weight_TauEleFakeScaleFactorDown");
            list_of_weights.push_back("total_weight_TauMuonFakeScaleFactorUp");
            list_of_weights.push_back("total_weight_TauMuonFakeScaleFactorDown");
            list_of_weights.push_back("total_weight_MuonIDScaleFactorUp");
            list_of_weights.push_back("total_weight_MuonIDScaleFactorDown");
            list_of_weights.push_back("total_weight_MuonISOScaleFactorUp");
            list_of_weights.push_back("total_weight_MuonISOScaleFactorDown");
            
            
            // pdf weights
            std::vector<std::string> colNames = total_weights.GetColumnNames();
            if (std::find(colNames.begin(), colNames.end(), "LHEPdfWeight") != colNames.end()) {
                for (uint i = 1; i < config::pdf_nweights; i++) {
                    total_weights = total_weights.Define("total_weight_pdf_weight_" + std::to_string(i), [i](const float& weight, const rvec<float>& pdf_weights)
                        {
                            if ( pdf_weights.size() == 0 ) return weight;
                            return weight*pdf_weights[i];
                        }, {"total_weight", "LHEPdfWeight"});
                    list_of_weights.push_back("total_weight_pdf_weight_" + std::to_string(i));
                }
            }
        }

            
        auto df_coll = total_weights.Define("CollMass", collinear_mass, {"col_idx","Tau_pt_ES", "Muon_tP_pt", "Tau_eta_ES", "Muon_eta", "Tau_phi_ES", "Muon_phi", "Tau_mass_ES", "Muon_mass", "MET_pt", "MET_phi","Tau_mask", "Muon_mask"});  
        //auto df_coll_alt = df_coll.Define("CollMass_alt", collinear_mass_alt, {"col_idx","Tau_pt_ES", "Muon_tP_pt", "Tau_eta_ES", "Muon_eta", "Tau_phi_ES", "Muon_phi", "Tau_mass_ES", "Muon_mass", "MET_pt", "MET_phi","Tau_mask", "Muon_mask"});  
        // Define Stage 0: any event with one tau fulfilling acceptance and id	
        //create_hists(df_coll_alt, "Stage0", "total_weight");
        
        // actual analysis cut  -- MT > 120. GeV 
        //auto df_mt = df_coll_alt.Filter("MT > 120.", "mt_cut");
        auto df_mt = df_coll.Filter("MT > 120.", "mt_cut");

        auto df_coll_alt = df_mt.Define("CollMass_alt", collinear_mass_alt, {"col_idx","Tau_pt_ES", "Muon_tP_pt", "Tau_eta_ES", "Muon_eta", "Tau_phi_ES", "Muon_phi", "Tau_mass_ES", "Muon_mass", "MET_pt", "MET_phi","Tau_mask", "Muon_mask"});  
        auto df_delphi_taumet = df_coll_alt.Define("DeltaPhi_tau_met", delta_phi, {"sel_Tau_phi", "MET_phi"});  
        auto df_delphi_cut = df_delphi_taumet.Filter("abs(DeltaPhi_tau_met) < 1.5707963267948966","delphi cut");
        create_hists(df_coll_alt, "Stage0", "total_weight");
        for (const auto& it : list_of_weights) {
            std::string folder_name = "Stage1";
            if ( it.find("total_weight_") != std::string::npos )
                folder_name = "Stage1/Systematics";
            if ( it.find("pdf_weight") != std::string::npos )
                folder_name = "Stage1/PDF";
            //create_hists(df_mt, folder_name, it);
            create_hists(df_delphi_cut, folder_name, it);
        }
        return df_delphi_cut;
   }();
   auto df_report = df_final.Report();
   // write all hists
   outFile->cd();
   for (const auto& it : hist_dict) {
       const auto& folder_name = it.first;
       auto hist = it.second;
       if (outFile->GetDirectory(folder_name.c_str()) == 0)
           outFile->mkdir(folder_name.c_str());
       outFile->cd(folder_name.c_str());
       hist->Write();
       outFile->cd();
   }
   for (const auto& it : hist_dict_2d) {
       const auto& folder_name = it.first;
       auto hist = it.second;
       if (outFile->GetDirectory(folder_name.c_str()) == 0)
           outFile->mkdir(folder_name.c_str());
       outFile->cd(folder_name.c_str());
       hist->Write();
       outFile->cd();
   }
   
   auto counter = TH1D("counter", "counter", 10u, 0, 10);
   counter.SetBinContent(1, *genWeightSum );
   counter.SetEntries(*eventcounter);
   counter.Write();
   df_report->Print();
    
}

std::vector < std::string > readfiles(const char* file_directory) {
	TSystemDirectory dir(file_directory, file_directory);
	TList *files = dir.GetListOfFiles();
	std::vector < std::string > rootfiles;
	if (files) {
		TSystemFile *file;
		TString fname;
		TIter next(files);
		while ((file=(TSystemFile*)next())) {
			fname = file->GetName();
			if (!file->IsDirectory() and fname.EndsWith(".root")) {
				std::string full_path = file_directory;
				full_path += file->GetName();
				rootfiles.push_back(full_path);
			}
		}
	}
	return rootfiles;
}



int main (int argc, char* argv[]) {
	
	// read in files to run over
	json files;
	std::ifstream files_json;
	files_json.open(argv[1]);
	files_json >> files;
	
	// loop over all types of backgrounds
	for (auto& it : files.items()) {
		
		
		// Create outfile name from path
		TPRegexp r1("[^/]+(?=/$|$)");
		TString tmp = (std::string) it.key();
		tmp = "output/" + tmp(r1)+ ".root";
		
		std::cout << "\t \t Running on: " << tmp(r1) << std::endl;
		
		// create outputfile to write hists
		TFile* outFile = new TFile(tmp, "RECREATE");
		
		// Files to run over
		std::string path = it.key();
		std::vector< std::string > names = readfiles(path.c_str());
		
		// read in config file for the corresponding file
		std::ifstream cfg_json;
		cfg_json.open(it.value());
		cfg_json >> cfg;
		
		// load analysis setup config file (see config.cc)
		config::load_config_file(cfg);
			
		// open root tree
		RNode df = ROOT::RDataFrame("Events", names);
	
        //// enable multithreading
        //ROOT::EnableImplicitMT();
		

        std::vector<std::string> colNames = df.GetColumnNames();
		if (std::find(colNames.begin(), colNames.end(), "genWeight") == colNames.end()) {
			df = df.Define("genWeight", "1.0");
		}
        auto tunep_pt = df.Define("Muon_tP_pt" , []( const rvec<float>& muon_pt, rvec<float>& muon_relTP_pt){
            return muon_relTP_pt * muon_pt;}, {"Muon_pt", "Muon_tunepRelPt"});
		auto loopcounter = tunep_pt.Filter([](ULong64_t e){if (0ULL == e) std::cout << "Running evtloop" << std::endl; return true;},{"rdfentry_"});
		
	
	
	
		
		if (config::runOnData) {
			std::ifstream goldenjson_file;
			goldenjson_file.open(cfg["json_file"]);
			goldenjson_file >> goldenjson;

			auto jsoncleaned = loopcounter.Filter(json_check, {"run", "luminosityBlock"}, "json cleaning");
			
			auto catch_up = jsoncleaned.Define("Tau_pt_ES", [](const rvec<float>& tau_pt){return tau_pt;}, {"Tau_pt"})
										.Define("Tau_eta_ES", [](const rvec<float>& tau_eta){return tau_eta;}, {"Tau_eta"})
										.Define("Tau_phi_ES", [](const rvec<float>& tau_phi){return tau_phi;}, {"Tau_phi"})
										.Define("Tau_mass_ES", [](const rvec<float>& tau_mass){return tau_mass;}, {"Tau_mass"});
		
			// this function does all analysis steps
			analyse(catch_up, outFile);
		
		} else {
			// counter - PSWeight is filled with ones
			auto definecounts = loopcounter.Define("abs_gen_weight", [](const float& x){return std::abs(x);}, {"genWeight"})
										   .Define("top_pt_weight", calc_top_pt_reweighting, {"GenPart_pdgId", "GenPart_pt"});
					
			// clean gen files
			auto gencleaned = definecounts.Filter(clean_gen_file, {"GenPart_pdgId", "GenPart_mass", "GenPart_pt", "GenPart_eta", "GenPart_phi"}, "gen cleaning");
			
			// calculate pileup weight
			auto df_pu_weight = gencleaned.Define("pileup_weight", [](const float& nvtx_true){ return pu_weight(nvtx_true, ""); }, {"Pileup_nTrueInt"});
			//auto df_pu_weight = definecounts.Define("pileup_weight", [](const float& nvtx_true){ return pu_weight(nvtx_true, ""); }, {"Pileup_nTrueInt"});
            
            auto calc_tau_es = df_pu_weight.Define("Tau_ES_vectors", calc_tau_energy_scale, {"Tau_decayMode", "Tau_genPartFlav", "Tau_pt", "Tau_eta", "Tau_phi", "Tau_mass"});
            			
			auto tau_es_applies = calc_tau_es.Define("Tau_pt_ES", [](const rvec< rvec < float > >& es_vector){return apply_tau_energy_scale(es_vector, 0);}, {"Tau_ES_vectors"})
											 .Define("Tau_eta_ES", [](const rvec< rvec < float > >& es_vector){return apply_tau_energy_scale(es_vector, 1);}, {"Tau_ES_vectors"})
											 .Define("Tau_phi_ES", [](const rvec< rvec < float > >& es_vector){return apply_tau_energy_scale(es_vector, 2);}, {"Tau_ES_vectors"})
											 .Define("Tau_mass_ES", [](const rvec< rvec < float > >& es_vector){return apply_tau_energy_scale(es_vector, 3);}, {"Tau_ES_vectors"});
						
			
			// this function does all analysis steps
			analyse(tau_es_applies, outFile);
		};
		
        config::clean_memory();
		hist_dict.clear();
		hist_dict_2d.clear();
		delete outFile;
		
	//	// disable multithreading
	//	ROOT::DisableImplicitMT();
	};
	return 0;
}
