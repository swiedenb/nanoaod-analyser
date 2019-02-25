#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include "TLorentzVector.h"

// reduce amount to write for each vector by a lot
template < typename T >
using rvec = ROOT::VecOps::RVec<T>;
using RNode = ROOT::RDF::RNode;

static std::unordered_map<std::string, TH1D * > histo;

// check if tau is in acceptance and fulfils id requirements
rvec<bool> tau_acceptance_and_id(	const rvec<float>& pt, 
									const rvec<float>& eta, 
									const rvec<Bool_t>& id, 
									const rvec<UChar_t>& iso, 
									const rvec<UChar_t>& antiEle_disc, 
									const rvec<UChar_t>& antiMu_disc) {
	// tau pt cut
	auto mask_pt = pt > 30;
	
	// tau eta cut
	auto mask_eta = abs(eta) < 2.3;
	
    // tau id decay mode
    auto mask_id = id;
    //
	// tau id iso (1: VVLoose, 2: VLoose, 4: Loose, 8: Medium, 16: Tight, 32: VTight, 64: VVTight)
	auto mask_iso = iso > 15 ;

	// Anti Electron discriminator (1: VLoose, 2: Loose, 4: Medium, 8: Tight, 16: VTight)
	auto mask_antiele = antiEle_disc > 0;
	
	// Anti Muon discriminatot (1: Loose, 2: Tight)
	auto mask_antimu = antiMu_disc > 1;
	
	// return vector with true, if particle fulfils all requirements - else false
	return mask_pt & mask_eta & mask_id & mask_iso & mask_antiele & mask_antimu;
};

// check if muon is in acceptance and fulfils id requirements
rvec<bool> muon_acceptance_and_id(	const rvec<float>& pt, 
									const rvec<float>& eta, 
									const rvec<UChar_t>& id, 
									const rvec<UChar_t>& iso
									) {
	// muon pt cut
	auto mask_pt = pt > 53;
	
	// muon eta cut
	auto mask_eta = abs(eta) < 2.4;
	
    // muon high pt id
    auto mask_id = id > 1;
    //
	// muon tk iso (1=TkIsoLoose, 2=TkIsoTight)
	auto mask_iso = iso == 1 ;

	
	// return vector with true, if particle fulfils all requirements - else false
	return mask_pt & mask_eta & mask_id & mask_iso;
};

// Calculate MT from one particle (which fulfils mask) and MET
rvec<float> transverse_mass(const rvec<float>& part_pt, 
							const rvec<float>& part_phi, 
							const float& met_pt, 
							const float& met_phi, 
							const rvec<bool>& mask) {
								
	return sqrt(2*part_pt[mask]*met_pt*cos(part_phi[mask]*met_phi));
};
	
// Calculate invariant mass from two particles
float invariant_mass(const rvec<float>& part1_pt,
                           const rvec<float>& part2_pt,
                           const rvec<float>& part1_eta,
                           const rvec<float>& part2_eta,
                           const rvec<float>& part1_phi,
                           const rvec<float>& part2_phi,
                           const rvec<float>& part1_mass,
                           const rvec<float>& part2_mass){
    TLorentzVector p1, p2;
    p1.SetPtEtaPhiM(part1_pt[0], part1_eta[0], part1_phi[0], part1_mass[0]);
    p2.SetPtEtaPhiM(part2_pt[0], part2_eta[0], part2_phi[0], part2_mass[0]);
    return (p1 + p2).M();
};

	
// Define pt over pt_miss Cut between particle and MET
bool pt_o_ptmiss(	const rvec<float>& part_pt, 
					const float& met_pt, 
					const rvec<bool>& mask) {
						
	if ((part_pt[mask][0] / met_pt) < 0.7)
		return false;
	else if ((part_pt[mask][0] / met_pt) > 1.4)
		return false;
	else
		return true;
};


// Particle number cut: how many 'trues' in mask
bool nparticle_cut(const rvec<bool>& mask) {
	return std::count(mask.begin(), mask.end(), true) > 0;
};

float calc_abs(const float& numb) {
	return std::abs(numb);
};

// Check if run number and lumi ranges are in json.txt - still need to write
bool json_check(const float& runnumber, const float& luminumber) {
	/* PSEUDOCODE
	 * if (runnumber in json) {
	 * 		for (luminumber1_json, luminumber2_json in json[runnumber]) {
	 * 			if (luminumber > luminumber1_json && luminumber < luminumber2_json)
	 * 				return true;
	 *			else
	 * 				return false;
	 * 		};
	 * };
	 * else
	 * 		return false;
	 */
	return true;
}

// returns vector of true particles from all candidates
rvec<float> selected_part_quant(const rvec<float>& quantity, const rvec<bool>& mask) {
	return quantity[mask];
};


// fill hist function - fills hists for each stage
void fill_hists(RNode df, TString& endung, TFile* outFile) {
	auto tau_pt = df.Histo1D("sel_Tau_pt");
	auto tau_eta = df.Histo1D("sel_Tau_eta");
	auto tau_phi = df.Histo1D("sel_Tau_phi");
	auto InvMass = df.Histo1D("InvMass");
	tau_pt->Write();
	tau_eta->Write();
	tau_phi->Write();
	InvMass->Write();
};

// analyse function - gets called for each systematic
void analyse_dataframe(RNode df, TFile* outFile, TString endung="") {
	TString folder = "";
	if (endung == "") {
		folder = "Selection";
	} else {
		folder = "Systs";
	};
	
	
	// Select trigger requirements
	auto triggered = df.Filter("HLT_Mu50 || HLT_TkMu100", "Trigger has to fire");
		
	
	
	// this need to be preselection
	// fill_preselection
	
	
	// creates mask, which tags taus which fulfil id and are in acceptance
	auto masked = triggered.Define("Tau_mask", tau_acceptance_and_id, {"Tau_pt", "Tau_eta", "Tau_idDecayMode" , "Tau_idMVAoldDM2017v2",  "Tau_idAntiEle", "Tau_idAntiMu"});
	
	// creates mask, which tags muons which fulfil id and are in acceptance
	auto masked2 = masked.Define("Muon_mask", muon_acceptance_and_id, {"Muon_pt", "Muon_eta", "Muon_highPtId", "Muon_tkIsoId"});

	// Select events with certain number of taus, which fulfil all acceptance & id
	auto evntselect = masked2.Filter(nparticle_cut, {"Tau_mask"});
	
	// Select events with certain number of muons, which fulfil all acceptance & id
	auto evntselect2 = evntselect.Filter(nparticle_cut, {"Muon_mask"});
	
	// Now that we have our real tau, write it to special columns to handle it more easily
	auto defquants = evntselect2.Define("sel_Tau_pt", selected_part_quant, {"Tau_pt", "Tau_mask"})
							   .Define("sel_Tau_eta", selected_part_quant, {"Tau_eta", "Tau_mask"})
							   .Define("sel_Tau_phi", selected_part_quant, {"Tau_phi", "Tau_mask"})
							   .Define("sel_Tau_mass", selected_part_quant, {"Tau_mass", "Tau_mask"})
                               .Define("sel_Muon_pt", selected_part_quant, {"Muon_pt", "Muon_mask"})
                               .Define("sel_Muon_eta", selected_part_quant, {"Muon_eta", "Muon_mask"})
                               .Define("sel_Muon_phi", selected_part_quant, {"Muon_phi", "Muon_mask"})
                               .Define("sel_Muon_mass", selected_part_quant, {"Muon_mass", "Muon_mask"});
	
	// Calculate MT distribution
	//auto mtcalc = defquants.Define("MT", transverse_mass, {"Tau_pt", "Tau_phi", "MET_pt", "MET_phi", "Tau_mask"});
	
    // Calculate invariant mass distribution
    auto invmasscalc = defquants.Define("InvMass", invariant_mass, {"Tau_pt", "Muon_pt", "Tau_eta", "Muon_eta", "Tau_phi", "Muon_phi", "Tau_mass", "Muon_mass"});
	// Define Stage 0: any event with one tau fulfilling acceptance and id
	outFile->cd();
	outFile->mkdir(folder + "/Stage0");
	outFile->cd(folder + "/Stage0");
	
	fill_hists(invmasscalc, endung, outFile);
	
	// actual analysis cut  -- 0.7 < pt/ptmiss < 1.4 
	//auto ptmiss = mtcalc.Filter(pt_o_ptmiss, {"Tau_phi", "MET_phi", "Tau_mask"});
	
	
	
	// Define Stage 1: fulfils ptmiss cut
//	outFile->cd();
//	outFile->mkdir(folder + "/Stage1");
//	outFile->cd(folder + "/Stage1");
//	
//	fill_hists(ptmiss, endung, outFile);
}


int main () {
	
	// enable multithreading
	ROOT::EnableImplicitMT();
	
	// define input files
	std::vector < std::string > names = {
		//~ "/net/scratch_cms3a/wiedenbeck/nanoaodfiles/ttbar/034EB431-E3E1-234D-A5FD-5D09B4BEC4DF.root ",
		//~ "/net/scratch_cms3a/wiedenbeck/nanoaodfiles/ttbar/1C2995B4-BEE6-FF43-B6C6-D00A7B165E09.root ",
		//~ "/net/scratch_cms3a/wiedenbeck/nanoaodfiles/ttbar/262A894B-282D-D04E-8CFA-97C2CCCC27E9.root ",
		//~ "/net/scratch_cms3a/wiedenbeck/nanoaodfiles/ttbar/2ABDDA8A-AE9A-C74D-97E1-FB6FB65FD1E6.root ",
		//~ "/net/scratch_cms3a/wiedenbeck/nanoaodfiles/ttbar/4A57A884-28D4-E04C-94F6-CFD905FC0C14.root ",
		//~ "/net/scratch_cms3a/wiedenbeck/nanoaodfiles/ttbar/79097697-485B-9542-8E6B-43A747EA7F4B.root ",
		//~ "/net/scratch_cms3a/wiedenbeck/nanoaodfiles/ttbar/8FBAA9CE-C37A-5E4C-95BB-7ACD29EB5B94.root ",
		//~ "/net/scratch_cms3a/wiedenbeck/nanoaodfiles/ttbar/A0F8D607-6B79-C348-B583-0BF3786455FF.root ",
		//~ "/net/scratch_cms3a/wiedenbeck/nanoaodfiles/ttbar/B0D91BD7-FD5C-F14F-9A7F-098F0BEE6658.root ",
		//~ "/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon_RunD_2017.root ",
		
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/0726FD7F-D51D-4348-8B48-FA2B7D2B403D.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/3E9394FA-EC63-114B-8AB1-18BD46BC5D43.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/7E78E240-4195-3C4F-976B-FEACA862289B.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/B727A5D5-B2C6-8345-9172-CACBB164CAB7.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/D377CDD4-D7E9-F84A-9AB9-AF9735FB5B46.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/0EF2170E-CD00-244D-8B6A-204A343D5947.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/46665B58-B1DF-4641-9CFB-B72822DDD495.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/7E8E6ACD-CDD3-7D43-9C0A-26D837ABF96A.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/B955E5C5-0EFB-9A41-A9A8-1836AC6E06E3.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/F477C46B-4D91-984E-9920-C29C70024778.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/1592DE0E-1474-8E43-85B6-C3BB3989B7F5.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/49F968DC-170C-F643-BE68-852DC4B59C93.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/815D3755-9B3F-734B-B305-1374A378D473.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/B98E474E-C2A7-394B-8401-BE169A390542.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/F60F986C-6A03-6A4C-ADFF-8C22BFE7F8CE.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/21A3C033-D062-7745-9892-A1194A4547BF.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/557920E8-7A80-C545-8ACF-3A83270949A7.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/81B82AF9-7AC2-BF44-8FB4-AEC8E0D3B66D.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/BA9E59EC-CFA3-9C4D-9369-0713CEA57FA8.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/FAE0760D-5E8D-954D-9821-5A195B165DFF.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/27EB22D8-139C-B841-BCC8-355820D75667.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/63CEC7BB-A7D4-D649-A3FB-DC8ABADFD4A7.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/83B8092C-126D-6E43-9644-252ED585D4A3.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/BF4124A8-53D5-704A-87FD-6B715507C1CD.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/FDBC726F-4376-D04D-89B6-43280B2642DB.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/2B117BC0-213F-EC4D-B1E3-A65C5D343F65.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/64468164-1F34-5542-82DB-337556EF5ECA.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/9FBB1A04-E382-6046-BD2B-3DE3FEE890BE.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/C859575A-CDB9-D94D-9382-711852F1CD7A.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/2F8304C8-2251-AD4C-A4E8-95C27ABD41EC.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/7A848068-4C62-1543-83D7-1A9B7BC4D871.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/ACB53DBE-22F8-0647-B7AC-DFF455950D4E.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/342194DB-2816-3146-9B7E-D6BA4452F6D8.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/7CF45BAA-3D37-E940-B9CC-0802987903C4.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/B562A8FC-B2C3-5046-A29A-C10E74A8D476.root",
		"/net/scratch_cms3a/wiedenbeck/nanoaodfiles/SingleMuon/D0E5F84C-796D-494E-8BCD-5FE856986950.root",
	};
	
	
	// create outputfile to write hists
	TFile* outFile = new TFile("histograms.root", "RECREATE");

	// create wanted folder structure
	outFile->mkdir("Custom");
	outFile->mkdir("Systs");
	outFile->mkdir("Selection");
	
	bool runOnData = false;
	
	
	// Let's Go - first, read in tree
	auto df = ROOT::RDataFrame("Events", names);
	
	// initialize counter
	//~ auto df1 = df.Define("absGenWeight", calc_abs, {"genWeight"});
	//~ auto counter = df1.Histo1D({"counter", "", 10u, 0, 10}, "absGenWeight");
	
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
