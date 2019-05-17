#include "Particles.hh"


// check if tau is in acceptance and fulfils id requirements
rvec<bool> tau_acceptance_and_id(	const rvec<float>& pt, 
									const rvec<float>& eta, 
									const rvec<UChar_t>& iso, 
									const rvec<UChar_t>& antiEle_disc, 
									const rvec<UChar_t>& antiMu_disc) {
	// tau pt cut
	auto mask_pt = pt > config::tau_pt;
	
	// tau eta cut
	auto mask_eta = abs(eta) < config::tau_eta;
	
	// tau iso requirement (1: VVLoose, 2: VLoose, 4: Loose, 8: Medium, 16: Tight, 32: VTight, 64: VVTight)
	auto mask_iso = (iso & config::tau_iso_WP) == config::tau_iso_WP;
	
	// Anti Electron discriminator (1: VLoose, 2: Loose, 4: Medium, 8: Tight, 16: VTight)
	auto mask_antiele = (config::tau_antiE_WP & antiEle_disc) == config::tau_antiE_WP;
	
	// Anti Muon discriminator (1: Loose, 2: Tight)
	auto mask_antimu = (config::tau_antiMu_WP & antiMu_disc) == config::tau_antiMu_WP;
	
	// return vector with true, if particle fulfils all requirements - else false
	return mask_pt & mask_eta & mask_iso & mask_antiele & mask_antimu;
};

// check if muon is in acceptance and fulfils id requirements
rvec<bool> muon_acceptance_and_id(	const rvec<float>& pt, 
									const rvec<float>& eta, 
									const rvec<bool>& id, 
									const rvec<float>& iso) {
	// muon pt cut
	auto mask_pt = pt > config::muon_pt;
	
	// muon eta cut
	auto mask_eta = abs(eta) < config::muon_eta;
	
    // muon id
    auto mask_id = id;
    
	// muon relative iso 
	auto mask_iso = iso / pt < 0.1;

	
	// return vector with true, if particle fulfils all requirements - else false
	return mask_pt & mask_eta & mask_id & mask_iso;
};

// check if electron is in acceptance and fulfils id requirements
rvec<bool> ele_acceptance_and_id(	const rvec<float>& pt,
									const rvec<float>& eta,
									const rvec<Int_t>& id) {
	// electron pt cut
	auto mask_pt = pt > 20;
	
	// electron eta cut
	auto mask_eta = abs(eta) < 2.5;
	
	// electron id cut
	auto mask_id = id >= 2;
	
	// return vector with true, if particle fulfils all requirements - else false
	return mask_pt & mask_eta & mask_id;
};


// apply tau energy scale
RNode apply_tau_energy_scale(	RNode df,
								const rvec<float>& tau_decayMode,
								const rvec<float>& tau_pt,
								const rvec<float>& tau_eta,
								const rvec<float>& tau_phi,
								const rvec<float>& tau_mass) {	
	rvec<float> new_tau_pt;				
	rvec<float> new_tau_eta;				
	rvec<float> new_tau_phi;				
	rvec<float> new_tau_mass;				
	
	for (uint i = 0; i<tau_pt.size(); i++) {
		if (tau_pt[i] > 400)
			new_tau_pt.push_back(tau_pt[i]);
			new_tau_eta.push_back(tau_eta[i]);
			new_tau_phi.push_back(tau_phi[i]);
			new_tau_mass.push_back(tau_mass[i]);
			continue;
		
		TLorentzVector p1, p2; 
		p1.SetPtEtaPhiM(tau_pt[i], tau_eta[i], tau_phi[i], tau_mass[i]);
		
		if (tau_decayMode[i] == 0) {										// this is one prong decay
			if (config::run_type == "")
				p1 *= (config::tau_energy_scale["h+-"])[1];
			if (config::run_type == "_TauEnergyScaleUp")
				p1 *= (config::tau_energy_scale["h+-"])[0];
			if (config::run_type == "_TauEnergyScaleDown")
				p1 *= (config::tau_energy_scale["h+-"])[2];
		} else if ((tau_decayMode[i] == 1) or (tau_decayMode[i] == 2)) {	// this is one prong decay with additional pi0s
			if (config::run_type == "")
				p1 *= (config::tau_energy_scale["h+- pi0s"])[1];
			if (config::run_type == "_TauEnergyScaleUp")
				p1 *= (config::tau_energy_scale["h+- pi0s"])[0];
			if (config::run_type == "_TauEnergyScaleDown")
				p1 *= (config::tau_energy_scale["h+- pi0s"])[2];
		} else if (tau_decayMode[i] == 10) {								// this is three prong decay
			if (config::run_type == "")
				p1 *= (config::tau_energy_scale["h+-h+-h+-"])[1];
			if (config::run_type == "_TauEnergyScaleUp")
				p1 *= (config::tau_energy_scale["h+-h+-h+-"])[0];
			if (config::run_type == "_TauEnergyScaleDown")
				p1 *= (config::tau_energy_scale["h+-h+-h+-"])[2];
		} else if (tau_decayMode[i] == 11) {								// this is one prong decay with additional pi0s
			if (config::run_type == "")
				p1 *= (config::tau_energy_scale["h+-h+-h+- pi0s"])[1];
			if (config::run_type == "_TauEnergyScaleUp")
				p1 *= (config::tau_energy_scale["h+-h+-h+- pi0s"])[0];
			if (config::run_type == "_TauEnergyScaleDown")
				p1 *= (config::tau_energy_scale["h+-h+-h+- pi0s"])[2];
		}
		
		new_tau_pt.push_back(p1.Pt());
		new_tau_eta.push_back(p1.Eta());
		new_tau_phi.push_back(p1.Phi());
		new_tau_mass.push_back(p1.M());
	}
	
	df = df.Define("Tau_pt_ES", [new_tau_pt](){return new_tau_pt;})
				.Define("Tau_eta_ES", [new_tau_eta](){return new_tau_eta;})
				.Define("Tau_phi_ES", [new_tau_phi](){return new_tau_phi;})
				.Define("Tau_mass_ES", [new_tau_mass](){return new_tau_mass;});
									
	
	return df;
};
