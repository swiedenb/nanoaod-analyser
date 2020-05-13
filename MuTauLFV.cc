#include "include/PhysicalQuantities.hh"
#include "include/Particles.hh"
#include "include/weights.hh"
#include "include/event_cleaner.hh"
#include "include/PDFTool.hh"
#include "include/TauIDSFTool.h"
#include "include/MuonSFTool.h"
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
	    hist_dict.emplace("Preselection", gentau_pt);
    }
	auto tau_pt = df.Histo1D(	{"Tau_pt", "", 			6000u, 0, 6000}, 					"Tau_pt_ES");
	auto tau_eta = df.Histo1D(	{"Tau_eta", "", 		100u, -5, 5}, 						"Tau_eta_ES");
	auto tau_phi = df.Histo1D(	{"Tau_phi", "", 		100u, -3.2, 3.2}, 					"Tau_phi_ES");	
	auto muon_pt = df.Histo1D(	{"Muon_pt", "", 		6000u, 0, 6000}, 					"Muon_pt");
	auto muon_eta = df.Histo1D(	{"Muon_eta", "", 		100u, -5, 5}, 						"Muon_eta");
	auto muon_phi = df.Histo1D(	{"Muon_phi", "", 		100u, -3.2, 3.2}, 					"Muon_phi");	
	auto met_pt = df.Histo1D(	{"MET_pt", "MET_pt", 	6000u, 0, 6000}, 					"MET_pt");
	auto met_phi = df.Histo1D(	{"MET_phi", "MET_phi", 	100u, -3.2, 3.2}, 					"MET_phi");
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
void calc_datadriven(RNode df) {
	// calc pt_ratio vector
	auto pt_ratio = df.Define("pt_ratio", ratio_vector, {"Tau_pt_ES", "MET_pt"});
	
	auto no_light_lepton = pt_ratio.Filter(nparticle_veto, {"Muon_mask"}, "muon_veto")
									.Filter(nparticle_veto, {"Electron_mask"}, "electron_veto");
	
	// begin evaluating the branch with non isolated taus for datadriven background estimation in QCD and ZToNuNu						 
	auto noniso_tau = no_light_lepton.Define("NonIso_Tau", [](	const rvec<float>& pt, 
															    const rvec<float>& eta,
														    	const rvec<bool>& dm,
														    	const rvec<UChar_t>& iso, 
														    	const rvec<UChar_t>& antiEle_disc, 
														    	const rvec<UChar_t>& antiMu_disc) {
		// check if tau is in acceptance and fulfils id requirements BUT NOT ISO
		auto mask_pt = pt > config::tau_pt;
	
		// tau eta cut
		auto mask_eta = abs(eta) < config::tau_eta;
												
		// tau decay mode
		auto mask_dm = dm;
													
		// tau iso requirement (1: VVLoose, 2: VLoose, 4: Loose, 8: Medium, 16: Tight, 32: VTight, 64: VVTight)
		auto mask_iso = !((iso & config::tau_iso_WP) == config::tau_iso_WP);
												
	    // Anti Electron discriminator (1: VLoose, 2: Loose, 4: Medium, 8: Tight, 16: VTight)
		auto mask_antiele = (config::tau_antiE_WP & antiEle_disc) == config::tau_antiE_WP;
												
		// Anti Muon discriminator (1: Loose, 2: Tight)
		auto mask_antimu = (config::tau_antiMu_WP & antiMu_disc) == config::tau_antiMu_WP;
												
		// return vector with true, if particle fulfils all requirements - else false
		rvec <bool> mask = mask_pt & mask_eta & mask_dm & mask_iso & mask_antiele & mask_antimu;
		return mask;
	}, 
											{"Tau_pt_ES", "Tau_eta_ES", config::tau_dm, config::tau_iso, "Tau_idAntiEle", "Tau_idAntiMu"});
											
		
	
	auto iso_tau = no_light_lepton.Define("Iso_Tau", tau_acceptance_and_id, {"Tau_pt_ES", "Tau_eta_ES", config::tau_dm, config::tau_iso, config::tau_antiEle, config::tau_antiMuon});
	
	
	auto noniso_select = noniso_tau.Filter([](const rvec<bool>& mask){return std::count(mask.begin(), mask.end(), true) >= 1;}, {"NonIso_Tau"}, "select at least one nonIso tau");
	auto noniso_redefine = noniso_select.Define("Tau_pt_new", select_part_quants, {"Tau_pt_ES", "NonIso_Tau"})
										.Define("Tau_eta_new", select_part_quants, {"Tau_eta_ES", "NonIso_Tau"})
										.Define("Tau_phi_new", select_part_quants, {"Tau_phi_ES", "NonIso_Tau"});
	
	
	
	auto iso_select = iso_tau.Filter(nparticle_cut, {"Iso_Tau"}, "select exactly one tau");
	auto iso_redefine = iso_select.Define("Tau_pt_new", select_part_quants, {"Tau_pt_ES", "Iso_Tau"})
									.Define("Tau_eta_new", select_part_quants, {"Tau_eta_ES", "Iso_Tau"})
									.Define("Tau_phi_new", select_part_quants, {"Tau_phi_ES", "Iso_Tau"});
	
	// select regions with pt/ptmiss > 1.5 and evaluate tight to loose ratio there
	auto noniso_nonsignal_region = noniso_redefine.Filter([](const rvec<float>& pt_ratio){
																	for (auto& it : pt_ratio) {
																		if (it >= 1.5) return true;
																	}
																	return false;
																}, {"pt_ratio"}, "Non-signal region - nonIsolated");
																
	auto noniso_signal_region = noniso_redefine.Filter([](const rvec<float>& pt_ratio){
																	for (auto& it : pt_ratio) {
																		if (it >= 1.5) return false;
																	}
																	return true;
																}, {"pt_ratio"}, "Signal region - nonIsolated");
																
																
	auto iso_nonsignal_region = iso_redefine.Filter([](const rvec<float>& pt_ratio){
																	for (auto& it : pt_ratio) {
																		if (it >= 1.5) return true;
																	}
																	return false;
																}, {"pt_ratio"}, "Non-signal region - Isolated");
																
	auto iso_signal_region = iso_redefine.Filter([](const rvec<float>& pt_ratio){
																	for (auto& it : pt_ratio) {
																		if (it >= 1.5) return false;
																	}
																	return true;
																}, {"pt_ratio"}, "Signal region - Isolated");
	
	fill_datadriven(noniso_nonsignal_region, "NonIsoNonSignal");
	fill_datadriven(noniso_signal_region, "NonIsoSignal");
	fill_datadriven(iso_nonsignal_region, "IsoNonSignal");
	fill_datadriven(iso_signal_region, "IsoSignal");
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
								
	auto met_filter = triggered.Filter("Flag_METFilters", "MET filters");
		
	// Trigger turn on cut
	//auto trigger_obj1 = met_filter.Filter("Muon_pt > " + std::to_string(config::muon_pt), "Avoid trigger turn on with muon object, pT>" + std::to_string(config::muon_pt));
	

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
							  .Define("Muon_mask", muon_acceptance_and_id, {"Muon_pt", "Muon_eta", "Muon_highPtId", "Muon_tkIsoId"})
							  .Define("DiMuon_mask", di_muon_id, {"Muon_pt", "Muon_eta", "Muon_highPtId", "Muon_tkRelIso"})
							  .Define("Electron_mask", ele_acceptance_and_simpleid, {"Electron_pt", "Electron_eta", "Electron_cutBased_HEEP"});
	
	
	
	// this is for datadriven part of analysis
	//calc_datadriven(masked);
	
	
	// Select events with certain number of taus, which fulfil all acceptance & id
	// also, cut any event with electrons or muons
	auto evntselect = masked.Filter(nparticle_cut, {"Tau_mask"}, "select_taus")
							.Filter(nparticle_cut, {"Muon_mask"}, "select_muons")
							.Filter(nparticle_veto, {"Electron_mask"}, "select_eles");
	
	
	// Now that we have our real tau, write it to special columns to handle it more easily
	auto defquants = evntselect.Define("sel_Tau_pt", selected_part_quant, {"Tau_pt_ES", "Tau_mask"})
							   .Define("sel_Tau_eta", selected_part_quant, {"Tau_eta_ES", "Tau_mask"})
							   .Define("sel_Tau_phi", selected_part_quant, {"Tau_phi_ES", "Tau_mask"})
                               .Define("sel_Muon_pt", selected_part_quant, {"Muon_pt", "Muon_mask"})
                               .Define("sel_Muon_eta", selected_part_quant, {"Muon_eta", "Muon_mask"})
                               .Define("sel_Muon_phi", selected_part_quant, {"Muon_phi", "Muon_mask"});
	
    //auto dimuon_cut = defquants.Filter(dimuonpair_cut,{"Muon_pt","Muon_eta","Muon_phi", "DiMuon_mask"},"Di Muon Pair Cut").Define("sel_tau_genPartFlav_int",[](const rvec<UChar_t>& genPartFlav){ 
    //    rvec<float> copy;
    //    for (const auto it: genPartFlav){
    //        copy.push_back((float) it);
    //    }
    //return copy;},{"Tau_genPartFlav"});
    auto dimuon_cut = defquants.Filter(dimuonpair_cut,{"Muon_pt","Muon_eta","Muon_phi", "DiMuon_mask"},"Di Muon Pair Cut");
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
                                                                         }, {"Muon_pt", "Muon_eta", "Muon_mask"})
										.Define("MuonISOScaleFactorUp", [](const rvec<float>& muon_pt,
                                                                         const rvec<float>& muon_eta,
                                                                         const rvec<bool>& muon_mask
                                                                         )
                                                                         {
                                                                            return muon_iso_scale_factor(muon_pt, muon_eta, muon_mask, "Up");
                                                                         }, {"Muon_pt", "Muon_eta", "Muon_mask"})
										.Define("MuonISOScaleFactorDown", [](const rvec<float>& muon_pt,
                                                                         const rvec<float>& muon_eta,
                                                                         const rvec<bool>& muon_mask
                                                                         )
                                                                         {
                                                                            return muon_iso_scale_factor(muon_pt, muon_eta, muon_mask, "Down");
                                                                         }, {"Muon_pt", "Muon_eta", "Muon_mask"})
										.Define("MuonIDScaleFactor", []( const rvec<float>& muon_pt,
                                                                         const rvec<float>& muon_eta,
                                                                         const rvec<bool>& muon_mask
                                                                         )
                                                                         {
                                                                            return muon_id_scale_factor(muon_pt, muon_eta, muon_mask, "");
                                                                         }, {"Muon_pt", "Muon_eta", "Muon_mask"})
		
										.Define("MuonIDScaleFactorUp", []( const rvec<float>& muon_pt,
                                                                         const rvec<float>& muon_eta,
                                                                         const rvec<bool>& muon_mask
                                                                         )
                                                                         {
                                                                            return muon_id_scale_factor(muon_pt, muon_eta, muon_mask, "Up");
                                                                         }, {"Muon_pt", "Muon_eta", "Muon_mask"})
										.Define("MuonIDScaleFactorDown", []( const rvec<float>& muon_pt,
                                                                         const rvec<float>& muon_eta,
                                                                         const rvec<bool>& muon_mask
                                                                         )
                                                                         {
                                                                            return muon_id_scale_factor(muon_pt, muon_eta, muon_mask, "Down");
                                                                         }, {"Muon_pt", "Muon_eta", "Muon_mask"});
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

	    
    auto df_col_idx = total_weights.Define("col_idx", col_idx, {"Tau_pt_ES", "Muon_pt", "Tau_eta_ES", "Muon_eta", "Tau_phi_ES", "Muon_phi", "Tau_mass_ES", "Muon_mass", "MET_pt", "MET_phi", "Tau_mask", "Muon_mask"});   
    auto df_coll = df_col_idx.Define("CollMass", collinear_mass, {"col_idx","Tau_pt_ES", "Muon_pt", "Tau_eta_ES", "Muon_eta", "Tau_phi_ES", "Muon_phi", "Tau_mass_ES", "Muon_mass", "MET_pt", "MET_phi","Tau_mask", "Muon_mask"});  
    auto df_coll_alt = df_coll.Define("CollMass_alt", collinear_mass_alt, {"col_idx","Tau_pt_ES", "Muon_pt", "Tau_eta_ES", "Muon_eta", "Tau_phi_ES", "Muon_phi", "Tau_mass_ES", "Muon_mass", "MET_pt", "MET_phi","Tau_mask", "Muon_mask"});  
    // Define Stage 0: any event with one tau fulfilling acceptance and id	
	create_hists(df_coll_alt, "Stage0", "total_weight");
	
	// actual analysis cut  -- MT > 120. GeV 
	auto df_mt = df_coll_alt.Filter("MT > 120.", "mt_cut");

    for (const auto& it : list_of_weights) {
        std::string folder_name = "Stage1";
        if ( it.find("total_weight_") != std::string::npos )
            folder_name = "Stage1/Systematics";
        if ( it.find("pdf_weight") != std::string::npos )
            folder_name = "Stage1/PDF";
        create_hists(df_mt, folder_name, it);
    }
	auto df_report = df_mt.Report();
	
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
		auto loopcounter = df.Filter([](ULong64_t e){if (0ULL == e) std::cout << "Running evtloop" << std::endl; return true;},{"rdfentry_"});
		
	
	
	
		
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
			auto gencleaned = definecounts.Filter(clean_gen_file, {"GenPart_pdgId", "GenPart_mass"}, "gen cleaning");
			
			// calculate pileup weight
			auto df_pu_weight = gencleaned.Define("pileup_weight", [](const float& nvtx_true){ return pu_weight(nvtx_true, ""); }, {"Pileup_nTrueInt"});
            
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
		delete outFile;
		
	//	// disable multithreading
	//	ROOT::DisableImplicitMT();
	};
	return 0;
}
