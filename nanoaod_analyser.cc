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

// Particle number cut: how many 'trues' in mask
bool nparticle_cut(const rvec<bool>& mask) {
	return std::count(mask.begin(), mask.end(), true) == 1;
};

// Particle number veto: how many 'trues' in mask
bool nparticle_veto(const rvec<bool>& mask) {
	return std::count(mask.begin(), mask.end(), true) == 0;
};

// Returns first mask_true particle for quantity
float selected_part_quant(const rvec<float>& quantity, 
								const rvec<bool>& mask) {
	return quantity[mask][0];
};

// fill preselection histograms
void fill_preselection(	RNode df,
						TFile* outFile) {
	auto tau_pt = df.Histo1D(	{"Tau_pt", "", 			6000u, 0, 6000}, 					"Tau_pt");
	auto tau_eta = df.Histo1D(	{"Tau_eta", "", 		100u, -5, 5}, 						"Tau_eta");
	auto tau_phi = df.Histo1D(	{"Tau_phi", "", 		100u, -3.2, 3.2}, 					"Tau_phi");	
	auto met_pt = df.Histo1D(	{"MET_pt", "MET_pt", 	6000u, 0, 6000}, 					"MET_pt");
	auto met_phi = df.Histo1D(	{"MET_phi", "MET_phi", 	100u, -3.2, 3.2}, 					"MET_phi");
	tau_pt->Write();
	tau_eta->Write();
	tau_phi->Write();
	met_pt->Write();
	met_phi->Write();
};


// fill hist function - fills hists for each stage
void fill_hists(RNode df, 
				TFile* outFile) {
					
	RNode dummy = df;
	if (config::runOnData) {
		dummy = df.Define("total_weight", "1.0");
	} else {
		dummy = df.Define("total_weight", "pileup_weight*TauScaleFactor*top_pt_weight*genWeight*1.0");
	}
	auto tau_pt = dummy.Histo1D(	{((TString) "Tau_pt" + config::run_type), "", 					6000u, 0, 6000}, 		"sel_Tau_pt", 				"total_weight");
	auto tau_eta = dummy.Histo1D(	{((TString) "Tau_eta" + config::run_type), "", 					100u, -5, 5}, 			"sel_Tau_eta", 				"total_weight");
	auto tau_phi = dummy.Histo1D(	{((TString) "Tau_phi" + config::run_type), "", 					100u, -3.2, 3.2}, 		"sel_Tau_phi", 				"total_weight");
	auto met_pt = dummy.Histo1D(	{((TString) "MET_pt" + config::run_type), "", 					6000u, 0, 6000}, 		"MET_pt", 					"total_weight");
	auto met_phi = dummy.Histo1D(	{((TString) "MET_phi" + config::run_type), "", 					100u, -3.2, 3.2}, 		"MET_phi", 					"total_weight");
	auto pt_ratio = dummy.Histo1D(	{((TString) "pT_ratio" + config::run_type), "", 				100u, 0, 10}, 			"pt_o_ptmiss", 				"total_weight");
	auto dphi = dummy.Histo1D(		{((TString) "DeltaPhi" + config::run_type), "", 				100u, 0, 3.2}, 			"dphi", 					"total_weight");
	auto MT = dummy.Histo1D(		{((TString) "MT" + config::run_type), "", 						6000u, 0, 6000},		"MT", 						"total_weight");
	tau_pt->Write();
	tau_eta->Write();
	tau_phi->Write();
	met_pt->Write();
	met_phi->Write();
	pt_ratio->Write();
	dphi->Write();
	MT->Write();
};

// analyse function - gets called for each systematic
void analyse(	RNode df, 
				TFile* outFile) {
	TString folder = "";
	if (config::run_type == "") {
		folder = "Selection";
	} else {
		folder = "Systs";
	};
	
	
	// Select trigger requirements
	auto triggered = df.Filter(cfg["trigger"],	"Trigger has to fire");
								
	auto met_filter = triggered.Filter("Flag_METFilters", "Remove events, which do not fulfil MET filters");
		
	// Trigger turn on cut, minimum tau requirement
	auto trigger_obj = met_filter.Filter("MET_pt > " + std::to_string(config::met_pt), "Avoid trigger turn on with trigger object")
								 .Filter("nTau >= 1", "At least one tau candidate");
	
	
	// fill preselection
	if (config::run_type == "") {
		outFile->cd();
		outFile->cd("Preselection");
		fill_preselection(triggered, outFile);
	}
	
	// creates mask, which fills bool tags for taus which fulfil id and are in acceptance
	auto masked = trigger_obj.Define("Tau_mask", tau_acceptance_and_id,	{"Tau_pt_ES", "Tau_eta_ES", "Tau_idMVAoldDM2017v2", "Tau_idAntiEle", "Tau_idAntiMu"})
							 .Define("Muon_mask", muon_acceptance_and_id, {"Muon_pt", "Muon_eta", "Muon_softId", "Muon_pfRelIso03_all"})
							 .Define("Electron_mask", ele_acceptance_and_id, {"Electron_pt", "Electron_eta", "Electron_cutBased"});
	
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
		auto defineScaleFactors = mtcalc.Define("TauScaleFactor", apply_scale_factor, {});
		
		
		saviour = defineScaleFactors;
	}
	
	// Define Stage 0: any event with one tau fulfilling acceptance and id
	outFile->cd();
	if (outFile->GetDirectory(folder + "/Stage0") == 0)
		outFile->mkdir(folder + "/Stage0");
	outFile->cd(folder + "/Stage0");
	
	fill_hists(saviour, outFile);
	
	
	
	
	// actual analysis cut  -- 0.7 < pt/ptmiss < 1.4 
	auto df_ptmiss = saviour.Filter("(pt_o_ptmiss > 0.7) && (pt_o_ptmiss < 1.3)", "pt_miss_cut");
	
	// Define Stage 1: fulfils ptmiss cut
	outFile->cd();
	if (outFile->GetDirectory(folder + "/Stage1") == 0)
		outFile->mkdir(folder + "/Stage1");
	outFile->cd(folder + "/Stage1");
	
	fill_hists(df_ptmiss, outFile);
	
	
	
	// dphi cut
	auto df_dphi = df_ptmiss.Filter("(dphi > 2.4) || (dphi < -2.4)", "deltaPhi_cut");
	
	// Define Stage 1: fulfils deltaphi cut
	outFile->cd();
	if (outFile->GetDirectory(folder + "/Stage2") == 0)
		outFile->mkdir(folder + "/Stage2");
	outFile->cd(folder + "/Stage2");
	
	fill_hists(df_dphi, outFile);
	
	if (!config::runOnData) {
		for (uint pdf_weight_i = 0; pdf_weight_i < (*df_dphi.Sum( "nLHEPdfWeight" ) / *df_dphi.Count()); pdf_weight_i++) {
			std::string pdf_column_title = "pdf_weight" + std::to_string(pdf_weight_i);
			
			auto pdf_weights = df_dphi.Define(pdf_column_title, [pdf_weight_i](rvec<float> pdf_weights){ return get_pdf_weight(pdf_weight_i, pdf_weights); }, {"LHEPdfWeight"});
			
			std::string final_weights_title = "final_weights" + std::to_string(pdf_weight_i);
			
			auto final_weights = pdf_weights.Define(final_weights_title, ( (std::string) "total_weight*" + pdf_column_title ) );

			auto pdf_hist = final_weights.Histo1D<double>( { ( (std::string) ("MT_pdf_" + std::to_string(pdf_weight_i)) ).c_str() , "MT", 5000u, 0 , 5000}, "MT", final_weights_title);

			pdf_hist->Write();
		}
	}
		
	
	// Print information about cut efficiencies
	df_dphi.Report()->Print();
	std::cout << "This was run_type: " << config::run_type << std::endl;
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
		
		// create outputfile to write hists
		TFile* outFile = new TFile(tmp, "RECREATE");
		
		// create wanted folder structure
		outFile->mkdir("Custom");
		outFile->mkdir("Preselection");
		outFile->mkdir("Selection");
		outFile->mkdir("Systs");
		
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
		auto df = ROOT::RDataFrame("Events", names);
		
		
		// enable multithreading
		ROOT::EnableImplicitMT();
		
		if (config::runOnData) {
			std::ifstream goldenjson_file;
			goldenjson_file.open(cfg["json_file"]);
			goldenjson_file >> goldenjson;

			auto jsoncleaned = df.Filter(json_check, {"run", "luminosityBlock"}, "json cleaning");
			// create counter
			auto counter = TH1I("counter", "counter", 10u, 0, 10);
			counter.SetBinContent(1, *jsoncleaned.Count());
			counter.SetEntries(*jsoncleaned.Count());
		
			// save counter
			outFile->cd();
			counter.Write();
			
			auto catch_up = jsoncleaned.Define("Tau_pt_ES", [](rvec<float>& tau_pt){return tau_pt;}, {"Tau_pt"})
										.Define("Tau_eta_ES", [](rvec<float>& tau_eta){return tau_eta;}, {"Tau_eta"})
										.Define("Tau_phi_ES", [](rvec<float>& tau_phi){return tau_phi;}, {"Tau_phi"})
										.Define("Tau_mass_ES", [](rvec<float>& tau_mass){return tau_mass;}, {"Tau_mass"});
		
			// this function does all analysis steps
			analyse(catch_up, outFile);
		
		} else {
			// counter - PSWeight is filled with ones
			auto definecounts = df.Define("abs_gen_weight", [](const float& x){return std::abs(x);}, {"genWeight"})
								  .Define("top_pt_weight", calc_top_pt_reweighting, {"GenPart_pdgId", "GenPart_pt"});
			auto counter = TH1D("counter", "counter", 10u, 0, 10);
			counter.SetBinContent(1, *df.Count());
			counter.SetEntries(*df.Count());
			counter.SetBinContent(2, *definecounts.Sum("abs_gen_weight"));
			// save counter
			outFile->cd();
			counter.Write();
		
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
			
			// start handling systematics
			//~ config::run_type = "_TauScaleUp";
			//~ analyse(df_pu_weight, outFile);
			
			//~ config::run_type = "_TauScaleDown";
			//~ analyse(df_pu_weight, outFile);
			
			
		};
		delete outFile;
		
		
		// enable multithreading
		ROOT::DisableImplicitMT();
	};
	return 0;
}
