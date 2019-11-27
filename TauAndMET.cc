#include "PhysicalQuantities.hh"
#include "Particles.hh"
#include "weights.hh"
#include "include/event_cleaner.hh"
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

// Particle number cut: how many "trues" in mask
bool nparticle_cut(const rvec<bool>& mask) {
	return std::count(mask.begin(), mask.end(), true) == 1;
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
RNode trigger(	RNode df,
				TFile* outFile) {
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
void fill_preselection(	RNode df,
						TFile* outFile) {
	auto tau_pt = df.Histo1D(	{"Tau_pt", "", 			6000u, 0, 6000}, 					"Tau_pt_ES");
	auto tau_eta = df.Histo1D(	{"Tau_eta", "", 		100u, -5, 5}, 						"Tau_eta_ES");
	auto tau_phi = df.Histo1D(	{"Tau_phi", "", 		100u, -3.2, 3.2}, 					"Tau_phi_ES");	
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
						TFile* outFile) {
	auto tau_pt = df.Histo1D(	{"Tau_pt", "", 			6000u, 0, 6000}, 					"Tau_pt_new");
	auto tau_eta = df.Histo1D(	{"Tau_eta", "", 		100u, -5, 5}, 						"Tau_eta_new");
	auto tau_phi = df.Histo1D(	{"Tau_phi", "", 		100u, -3.2, 3.2}, 					"Tau_phi_new");	
	auto met_pt = df.Histo1D(	{"MET_pt", "MET_pt", 	6000u, 0, 6000}, 					"MET_pt");
	auto met_phi = df.Histo1D(	{"MET_phi", "MET_phi", 	100u, -3.2, 3.2}, 					"MET_phi");
	hist_dict.emplace("Datadriven", tau_pt);
	hist_dict.emplace("Datadriven", tau_eta);
	hist_dict.emplace("Datadriven", tau_phi);
	hist_dict.emplace("Datadriven", met_pt);
	hist_dict.emplace("Datadriven", met_phi);
};


// fill hist function - fills hists for each stage
void create_hists(	RNode df,
					std::string name = "") {
					
	RNode dummy = df;
	if (config::runOnData) {
		dummy = df.Define("total_weight", "1.0");
	} else {
		if (config::W_kfactor_hist != NULL) {
			dummy = df.Define("total_weight", "pileup_weight*TauScaleFactor*top_pt_weight*genWeight*W_kfactor*TauEleFakeScaleFactor*TauMuonFakeScaleFactor*PrefiringWeight*1.0");
		} else {
			dummy = df.Define("total_weight", "pileup_weight*TauScaleFactor*top_pt_weight*genWeight*TauEleFakeScaleFactor*TauMuonFakeScaleFactor*PrefiringWeight*1.0");
		}
	}
	auto tau_pt = dummy.Histo1D(	{((TString) "Tau_pt" + config::run_type), "", 					6000u, 0, 6000}, 		"sel_Tau_pt", 				"total_weight");
	auto tau_eta = dummy.Histo1D(	{((TString) "Tau_eta" + config::run_type), "", 					100u, -5, 5}, 			"sel_Tau_eta", 				"total_weight");
	auto tau_phi = dummy.Histo1D(	{((TString) "Tau_phi" + config::run_type), "", 					100u, -3.2, 3.2}, 		"sel_Tau_phi", 				"total_weight");
	auto met_pt = dummy.Histo1D(	{((TString) "MET_pt" + config::run_type), "", 					6000u, 0, 6000}, 		"MET_pt", 					"total_weight");
	auto met_phi = dummy.Histo1D(	{((TString) "MET_phi" + config::run_type), "", 					100u, -3.2, 3.2}, 		"MET_phi", 					"total_weight");
	auto pt_ratio = dummy.Histo1D(	{((TString) "pT_ratio" + config::run_type), "", 				100u, 0, 10}, 			"pt_o_ptmiss", 				"total_weight");
	auto dphi = dummy.Histo1D(		{((TString) "DeltaPhi" + config::run_type), "", 				100u, 0, 3.2}, 			"dphi", 					"total_weight");
	auto MT = dummy.Histo1D(		{((TString) "MT" + config::run_type), "", 						6000u, 0, 6000},		"MT", 						"total_weight");
	auto nvtx = dummy.Histo1D(		{((TString) "nvtx" + config::run_type), "", 					80u, 0, 80}, 			"PV_npvs", 					"total_weight");
	auto nvtxgood = dummy.Histo1D(	{((TString) "nvtx_good" + config::run_type), "", 				80u, 0, 80}, 			"PV_npvsGood", 				"total_weight");
	hist_dict.emplace(name, tau_pt);
	hist_dict.emplace(name, tau_eta);
	hist_dict.emplace(name, tau_phi);
	hist_dict.emplace(name, met_pt);
	hist_dict.emplace(name, met_phi);
	hist_dict.emplace(name, pt_ratio);
	hist_dict.emplace(name, dphi);
	hist_dict.emplace(name, MT);
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
void calc_datadriven(TFile* outFile, RNode df) {
	// calc pt_ratio vector
	auto pt_ratio = df.Define("pt_ratio", ratio_vector, {"Tau_pt_ES", "MET_pt"});
	
	auto no_light_lepton = pt_ratio.Filter(nparticle_veto, {"Muon_mask"}, "muon_veto")
									.Filter(nparticle_veto, {"Electron_mask"}, "electron_veto");
	
	// begin evaluating the branch with non isolated taus for datadriven background estimation in QCD and ZToNuNu						 
	auto noniso_tau = no_light_lepton.Define("NonIso_Tau", [](		const rvec<float>& pt, 
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
	
	outFile->cd();
	outFile->mkdir("Datadriven");
	outFile->mkdir("Datadriven/NonSignalNonIso");
	outFile->mkdir("Datadriven/SignalNonIso");
	outFile->mkdir("Datadriven/NonSignalIso");
	outFile->mkdir("Datadriven/SignalIso");
	
	outFile->cd("Datadriven/NonSignalNonIso");
	fill_datadriven(noniso_nonsignal_region, outFile);
	outFile->cd("Datadriven/SignalNonIso");
	fill_datadriven(noniso_signal_region, outFile);
	outFile->cd("Datadriven/NonSignalIso");
	fill_datadriven(iso_nonsignal_region, outFile);
	outFile->cd("Datadriven/SignalIso");
	fill_datadriven(iso_signal_region, outFile);
	outFile->cd();
};

// analyse function - gets called for each systematic
void analyse(	RNode df, 
				TFile* outFile) {	
	
	// init counter
	auto eventcounter = df.Sum("genWeight");
	
	// Select trigger requirements
	auto triggered = trigger(df, outFile);
								
	auto met_filter = triggered.Filter("Flag_METFilters", "MET filters");
		
	// Trigger turn on cut
	auto trigger_obj1 = met_filter.Filter("MET_pt > " + std::to_string(config::met_pt), "Avoid trigger turn on with MET object, pT>" + std::to_string(config::met_pt));
	
	auto trigger_obj2 = trigger_obj1.Filter([]  (const rvec<float>& tau_pt) {	if (tau_pt.size() == 0) return false;
																				if (tau_pt[0] < config::tau_pt) return false;
																				return true;
																			}, {"Tau_pt"}, "Avoid trigger turn on with tau object");
																		
	auto trigger_obj3 = trigger_obj2.Filter([] (const rvec<float>& jet_pt,
											    const rvec<float>& jet_eta){	
												    if (jet_pt.size() == 0) return false;
													for (uint i = 0; i < jet_pt.size(); i++) {
														if (abs(jet_eta[i]) > 2.3)
															continue;
														if (jet_pt[i] > 100)
															return true;
														}
	                    								return true;
	                    							}, {"Jet_pt", "Jet_eta"}, "Avoid trigger turn on with jet object, pT>100");
	

	auto good_primary_vertex = trigger_obj3.Filter([](	const float& pvx, 
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
	if (config::era == 2018) {
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
	}
	
	// fill preselection
	if (config::run_type == "") {
		fill_preselection(trigger_obj3, outFile);
	}
	
	// creates mask, which fills bool tags for taus which fulfil id and are in acceptance
	auto masked = trigger_obj3.Define("Tau_mask", tau_acceptance_and_id,	{"Tau_pt_ES", "Tau_eta_ES", config::tau_dm, config::tau_iso, config::tau_antiEle, config::tau_antiMuon})
							  .Define("Muon_mask", muon_acceptance_and_id, {"Muon_pt", "Muon_eta", "Muon_softId", "Muon_pfRelIso03_all"})
							  .Define("Electron_mask", ele_acceptance_and_id, {"Electron_pt", "Electron_eta", "Electron_cutBased"});
	
	
	
	// this is for datadriven part of analysis
	//~ calc_datadriven(outFile, masked);
	
	
	// Select events with certain number of taus, which fulfil all acceptance & id
	// also, cut any event with electrons or muons
	auto evntselect = masked.Filter(nparticle_cut, {"Tau_mask"}, "select_one_tau")
							.Filter(nparticle_veto, {"Muon_mask"}, "muon_veto")
							.Filter(nparticle_veto, {"Electron_mask"}, "electron_veto");
	
	
	// Now that we have our real tau, write it to special columns to handle it more easily
	auto defquants = evntselect.Define("sel_Tau_pt", selected_part_quant, {"Tau_pt_ES", "Tau_mask"})
							   .Define("sel_Tau_eta", selected_part_quant, {"Tau_eta_ES", "Tau_mask"})
							   .Define("sel_Tau_phi", selected_part_quant, {"Tau_phi_ES", "Tau_mask"});
	
	// Calculate MT distribution
	auto mtcalc = defquants.Define("MT", mass_transv, {"sel_Tau_pt", "sel_Tau_phi", "MET_pt", "MET_phi"})
						   .Define("pt_o_ptmiss", ratio, {"sel_Tau_pt", "MET_pt"})
						   .Define("dphi", delta_phi, {"sel_Tau_phi", "MET_phi"});
	
	
	RNode saviour = mtcalc;
	// Define and calculate weights for monte carlo
	if (!config::runOnData) {
		auto defineScaleFactors = mtcalc.Define("TauScaleFactor", apply_scale_factor, {})
										.Define("TauEleFakeScaleFactor", [](const float& tau_eta,
										                                    const float& tau_phi,
										                                    const rvec<int>& gen_pdgID,
										                                    const rvec<float>& gen_eta,
										                                    const rvec<float>& gen_phi) {
																				return tau_fake_scale_factor(tau_eta, tau_phi, gen_pdgID, gen_eta, gen_phi, 11);
																			}, {"sel_Tau_eta", "sel_Tau_phi", "GenPart_pdgId", "GenPart_eta", "GenPart_phi"})
										.Define("TauMuonFakeScaleFactor", [](const float& tau_eta,
										                                    const float& tau_phi,
										                                    const rvec<int>& gen_pdgID,
										                                    const rvec<float>& gen_eta,
										                                    const rvec<float>& gen_phi) {
																				return tau_fake_scale_factor(tau_eta, tau_phi, gen_pdgID, gen_eta, gen_phi, 13);
																			}, {"sel_Tau_eta", "sel_Tau_phi", "GenPart_pdgId", "GenPart_eta", "GenPart_phi"})
										.Define("PrefiringWeight", prefire_factor, {"Jet_pt", "Jet_eta", "Photon_pt", "Photon_eta"});
		
		if (config::W_kfactor_hist != NULL) {
			saviour = defineScaleFactors.Define("W_kfactor", get_kfactor, {"GenPart_pdgId", "GenPart_mass"});
		}
		else {
			saviour = defineScaleFactors;
		}
	}
	
	// Define Stage 0: any event with one tau fulfilling acceptance and id	
	create_hists(saviour, "Stage0");
	
	// actual analysis cut  -- 0.7 < pt/ptmiss < 1.3 
	auto df_ptmiss = saviour.Filter("(pt_o_ptmiss > 0.7) && (pt_o_ptmiss < 1.3)", "pt_miss_cut");
	create_hists(df_ptmiss, "Stage1");
	
	// dphi cut
	auto df_dphi = df_ptmiss.Filter("(dphi > 2.4) || (dphi < -2.4)", "deltaPhi_cut");
	create_hists(df_dphi, "Stage2");
	
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
	counter.SetBinContent(1, *eventcounter);
	counter.SetEntries(*eventcounter);
	counter.Write();


	//~ if (!config::runOnData) {
		//~ for (uint pdf_weight_i = 0; pdf_weight_i < (*df_dphi.Sum( "nLHEPdfWeight" ) / *df_dphi.Count()); pdf_weight_i++) {
			//~ std::string pdf_column_title = "pdf_weight" + std::to_string(pdf_weight_i);
			
			//~ auto pdf_weights = df_dphi.Define(pdf_column_title, [pdf_weight_i](rvec<float> pdf_weights){ return get_pdf_weight(pdf_weight_i, pdf_weights); }, {"LHEPdfWeight"});
			
			//~ std::string final_weights_title = "final_weights" + std::to_string(pdf_weight_i);
			
			//~ auto final_weights = pdf_weights.Define(final_weights_title, ( (std::string) "pileup_weight*TauScaleFactor*top_pt_weight*genWeight*" + pdf_column_title ) );

			//~ auto pdf_hist = final_weights.Histo1D<double>( { ( (std::string) ("MT_pdf_" + std::to_string(pdf_weight_i)) ).c_str() , "MT", 5000u, 0 , 5000}, "MT", final_weights_title);

			//~ pdf_hist->Write();
		//~ }
	//~ }
		
	
	// Print information about cut efficiencies !! WARNING; THIS TRIGGERS A NEW EVENT LOOP !!!
	//~ df_dphi.Report()->Print();
	//~ std::cout << "This was run_type: " << config::run_type << std::endl;
	
	//~ // Prints the graph to the rd1.dot file in the current directory
	//~ ROOT::RDF::SaveGraph(df, "./mydot.dot");
	//~ // Prints the graph to standard output
	//~ ROOT::RDF::SaveGraph(df);
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
		std::vector<std::string> colNames = df.GetColumnNames();
		if (std::find(colNames.begin(), colNames.end(), "genWeight") == colNames.end()) {
			df = df.Define("genWeight", "1.0");
		}
		auto loopcounter = df.Filter([](ULong64_t e){if (0ULL == e) std::cout << "Running evtloop" << std::endl; return true;},{"rdfentry_"});
		
		
		// enable multithreading
		ROOT::EnableImplicitMT();
		
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
			auto df_pu_weight = gencleaned.Define("pileup_weight", pu_weight, {"Pileup_nTrueInt"});
			
			auto tau_es_applies = df_pu_weight.Define("Tau_pt_ES", apply_tau_energy_scale_on_pt, {"Tau_decayMode", "Tau_pt", "Tau_eta", "Tau_phi", "Tau_mass"})
												.Define("Tau_eta_ES", apply_tau_energy_scale_on_eta, {"Tau_decayMode", "Tau_pt", "Tau_eta", "Tau_phi", "Tau_mass"})
												.Define("Tau_phi_ES", apply_tau_energy_scale_on_phi, {"Tau_decayMode", "Tau_pt", "Tau_eta", "Tau_phi", "Tau_mass"})
												.Define("Tau_mass_ES", apply_tau_energy_scale_on_mass, {"Tau_decayMode", "Tau_pt", "Tau_eta", "Tau_phi", "Tau_mass"});
			
			
			// this function does all analysis steps
			analyse(tau_es_applies, outFile);
		};
		
		hist_dict.clear();
		delete outFile;
		
		// disable multithreading
		ROOT::DisableImplicitMT();
	};
	return 0;
}
