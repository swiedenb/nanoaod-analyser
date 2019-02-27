#include "PhysicalQuantities.hh"
#include "include/event_cleaner.hh"
#include <iostream>
#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RVec.hxx>
#include <math.h>

// reduce amount to write for each vector by a lot
template < typename T >
using rvec = ROOT::VecOps::RVec<T>;
using RNode = ROOT::RDF::RNode;

json goldenjson;
std::ifstream goldenjson_file;

// check if tau is in acceptance and fulfils id requirements
rvec<bool> tau_acceptance_and_id(	const rvec<float>& pt, 
									const rvec<float>& eta, 
									const rvec<Bool_t>& id, 
									const rvec<UChar_t>& iso, 
									const rvec<UChar_t>& antiEle_disc, 
									const rvec<UChar_t>& antiMu_disc) {
	// tau pt cut
	auto mask_pt = pt > 80;
	
	// tau eta cut
	auto mask_eta = abs(eta) < 2.3;
	
	// tau id requirement 
	auto mask_id = id;
	
	// tau iso requirement (1: VVLoose, 2: VLoose, 4: Loose, 8: Medium, 16: Tight, 32: VTight, 64: VVTight)
	auto mask_iso = iso > 1;
	
	// Anti Electron discriminator (1: VLoose, 2: Loose, 4: Medium, 8: Tight, 16: VTight)
	auto mask_antiele = antiEle_disc > 0;
	
	// Anti Muon discriminatot (1: Loose, 2: Tight)
	auto mask_antimu = antiMu_disc > 0;
	
	// return vector with true, if particle fulfils all requirements - else false
	return mask_pt & mask_eta & mask_id & mask_iso & mask_antiele & mask_antimu;
};


// Particle number cut: how many 'trues' in mask
bool nparticle_cut(const rvec<bool>& mask) {
	return std::count(mask.begin(), mask.end(), true) == 1;
};

// returns vector of true particles from all candidates
float selected_part_quant(const rvec<float>& quantity, 
								const rvec<bool>& mask) {
	return quantity[mask][0];
};


// fill preselection histograms
void fill_preselection(	RNode df,
						TFile* outFile,
						TString& endung) {
	auto tau_pt = df.Histo1D(	{"Tau_pt", "", 			6000u, 0, 6000}, 					"Tau_pt");
	auto tau_eta = df.Histo1D(	{"Tau_eta", "", 		100u, -5, 5}, 						"Tau_eta");
	auto tau_phi = df.Histo1D(	{"Tau_phi", "", 		100u, -3.2, 3.2}, 					"Tau_phi");	
	auto met_pt = df.Histo1D(	{"MET_pt", "MET_pt", 	6000u, 0, 6000}, 					"MET_pt");
	auto met_phi = df.Histo1D(	{"MET_phi", "MET_phi", 	100u, -3.2, 3.2}, 					"MET_phi");
};


// fill hist function - fills hists for each stage
void fill_hists(RNode df, 
				TFile* outFile,
				TString& endung) {
	auto tau_pt = df.Histo1D(	{"Tau_pt", "", 						6000u, 0, 6000}, 		"sel_Tau_pt");
	auto tau_eta = df.Histo1D(	{"Tau_eta", "", 					100u, -5, 5}, 			"sel_Tau_eta");
	auto tau_phi = df.Histo1D(	{"Tau_phi", "", 					100u, -3.2, 3.2}, 		"sel_Tau_phi");
	auto pt_ratio = df.Histo1D(	{"p_{T}^{tau} / p_{T}^{miss}", "", 	100u, 0, 10}, 			"pt_o_ptmiss");
	auto dphi = df.Histo1D(		{"DeltaPhi", "", 					100u, 0, 3.2}, 			"dphi");
	auto MT = df.Histo1D(		{"MT", "", 							6000u, 0, 6000},		"MT");
	tau_pt->Write();
	tau_eta->Write();
	tau_phi->Write();
	pt_ratio->Write();
	dphi->Write();
	MT->Write();
};

// analyse function - gets called for each systematic
void analyse_dataframe(RNode df, 
						TFile* outFile, 
						TString endung="") {
	TString folder = "";
	if (endung == "") {
		folder = "Selection";
	} else {
		folder = "Systs";
	};
	
	
	// Select trigger requirements
	auto triggered = df.Filter("HLT_CaloMET90_HBHECleaned || \
								HLT_CaloMET100_HBHECleaned || \
								HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100 || \
								HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110	|| \
								HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120	|| \
								HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130	|| \
								HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140	|| \
								HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90", 
								"Trigger has to fire");
		
	// Trigger turn on cut, minimum tau requirement
	auto trigger_obj = triggered.Filter("MET_pt > 150", "Avoid trigger turn on with trigger object")
								.Filter("nTau >= 1", "At least one tau candidate");
	
	
	// this need to be preselection
	// fill_preselection
	fill_preselection(triggered, outFile, endung);
	
	
	// creates mask, which tags taus which fulfil id and are in acceptance
	auto masked = trigger_obj.Define("Tau_mask", tau_acceptance_and_id, {"Tau_pt", "Tau_eta", "Tau_idDecayMode", "Tau_idMVAoldDM2017v2", "Tau_idAntiEle", "Tau_idAntiMu"});
	
	// Select events with certain number of taus, which fulfil all acceptance & id
	auto evntselect = masked.Filter(nparticle_cut, {"Tau_mask"});
	
	
	// Now that we have our real tau, write it to special columns to handle it more easily
	auto defquants = evntselect.Define("sel_Tau_pt", selected_part_quant, {"Tau_pt", "Tau_mask"})
							   .Define("sel_Tau_eta", selected_part_quant, {"Tau_eta", "Tau_mask"})
							   .Define("sel_Tau_phi", selected_part_quant, {"Tau_phi", "Tau_mask"});
	
	// Calculate MT distribution
	auto mtcalc = defquants.Define("MT", transverse_mass, {"sel_Tau_pt", "sel_Tau_phi", "MET_pt", "MET_phi"})
						   .Define("pt_o_ptmiss", ratio, {"sel_Tau_pt", "MET_pt"})
						   .Define("dphi", delta_phi, {"sel_Tau_phi", "MET_phi"});
	
	// Define Stage 0: any event with one tau fulfilling acceptance and id
	outFile->cd();
	outFile->mkdir(folder + "/Stage0");
	outFile->cd(folder + "/Stage0");
	
	fill_hists(mtcalc, outFile, endung);
	
	// actual analysis cut  -- 0.7 < pt/ptmiss < 1.4 
	auto df_ptmiss = mtcalc.Filter("(pt_o_ptmiss > 0.7) && (pt_o_ptmiss < 1.4)", "pt_miss_cut");
	
	// Define Stage 1: fulfils ptmiss cut
	outFile->cd();
	outFile->mkdir(folder + "/Stage1");
	outFile->cd(folder + "/Stage1");
	
	fill_hists(df_ptmiss, outFile, endung);
	
	// dphi cut
	auto df_dphi = df_ptmiss.Filter("(dphi > 2.4)", "deltaPhi_cut");
	
	// Define Stage 1: fulfils ptmiss cut
	outFile->cd();
	outFile->mkdir(folder + "/Stage2");
	outFile->cd(folder + "/Stage2");
	
	fill_hists(df_dphi, outFile, endung);
}


int main (int argc, char* argv[]) {
	// argc = argument count; argv = argument 'vector'
	// Execute analyser with a list of files given in command line
	std::vector< std::string > names;
	if(argc < 2) {
        throw std::invalid_argument("Please provide exactly one input argument: The path to the input file.");
    } else {
        names.assign(argv + 1, argv + argc);
    }
	
	// enable multithreading
	ROOT::EnableImplicitMT();
	
	
	// create outputfile to write hists
	TFile* outFile = new TFile("histograms.root", "RECREATE");

	// create wanted folder structure
	outFile->mkdir("Custom");
	outFile->mkdir("Preselection");
	outFile->mkdir("Selection");
	outFile->mkdir("Systs");
	
	bool runOnData = true;
	
	
	// Let's Go - first, read in tree
	auto df = ROOT::RDataFrame("Events", names);
	
	// initialize counter
	//~ auto df1 = df.Define("absGenWeight", calc_abs, {"genWeight"});
	//~ auto counter = df1.Histo1D({"counter", "", 10u, 0, 10}, "absGenWeight");
	

	goldenjson_file.open("cfg/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt");
	goldenjson_file >> goldenjson;
	
	// if we are working on datafiles, throw out all non-golden json events
	if (runOnData) {
		
		auto jsoncleaned = df.Filter(json_check, {"run", "luminosityBlock"});
	} else {
		// consider "genWeight"
	};
	
	// this function does all analysis steps
	analyse_dataframe(df, outFile);
	
	// save counter
	outFile->cd();
	//~ counter->Write();
	
	delete outFile;
	
	// creates logic graph of analysis
	ROOT::RDF::SaveGraph(df, "./mydot.dot");
	
	return 0;
}
