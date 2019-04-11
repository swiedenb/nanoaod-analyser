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
	auto mask_iso = (iso & *config::tau_iso_WP) != 0;
	
	// Anti Electron discriminator (1: VLoose, 2: Loose, 4: Medium, 8: Tight, 16: VTight)
	auto mask_antiele = (*config::tau_antiE_WP & antiEle_disc) != 0;
	
	// Anti Muon discriminator (1: Loose, 2: Tight)
	auto mask_antimu = (*config::tau_antiMu_WP & antiMu_disc) != 0;
	
	// return vector with true, if particle fulfils all requirements - else false
	return mask_pt & mask_eta & mask_iso & mask_antiele & mask_antimu;
};

// check if muon is in acceptance and fulfils id requirements
rvec<bool> muon_acceptance_and_id(	const rvec<float>& pt, 
									const rvec<float>& eta, 
									const rvec<bool>& id, 
									const rvec<float>& iso) {
	// muon pt cut
	auto mask_pt = pt > 20;
	
	// muon eta cut
	auto mask_eta = abs(eta) < 2.4;
	
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
