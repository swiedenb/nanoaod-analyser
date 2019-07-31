#include "weights.hh"

// Get Pileup weight from file
float pu_weight( const float nvtx_true){
    auto weight = config::pileup_hist->GetBinContent(nvtx_true);
    return weight;
};


// get tau scale factor
float apply_scale_factor() {
    if (config::run_type == "_TauScaleUp") {
        return config::tau_scale_up;
    } else if (config::run_type == "_TauScaleDown") {
        return config::tau_scale_down;
    } else {
        return config::tau_scale;
	}
};

// get tau fake rate scale factor (for both e an mu)
float tau_fake_scale_factor(	const float& tau_eta,
								const float& tau_phi,
								const rvec<int>& gen_pdgID,
								const rvec<float>& gen_eta,
								const rvec<float>& gen_phi,
								const int& required_pdgID) {
	
	auto deltaR = [](const float& deltaPhi, const float& deltaEta){
		return sqrt( pow(deltaPhi, 2) + pow(deltaEta, 2) );
	};
	
	bool match = false;
	for (uint i = 0; i < gen_pdgID.size(); i++) {
		if (match) break;
		if ( std::abs(gen_pdgID[i]) != required_pdgID )
			continue;
		if ( deltaR( delta_phi(tau_phi, gen_phi[i]), delta_eta(tau_eta, gen_eta[i]) < 0.3 ) ) {
			match = true;
		}
	}
	
	float scale_factor = 1.0;
	
	if (match) {
		if (required_pdgID == 11) {
			scale_factor *= config::tau_ele_fake_hist->GetBinContent( config::tau_ele_fake_hist->GetXaxis()->FindBin( std::abs(tau_eta) ), config::tau_ele_fake_hist->GetYaxis()->FindBin( config::tau_antiE_WP ) );
		} 
		else if (required_pdgID == 13)  {
			scale_factor *= config::tau_muo_fake_hist->GetBinContent( config::tau_muo_fake_hist->GetXaxis()->FindBin( std::abs(tau_eta) ), config::tau_ele_fake_hist->GetYaxis()->FindBin( config::tau_antiMu_WP ) );
		}
	}
	
	return scale_factor;
};


// calculate top pt reweighting
float calc_top_pt_reweighting( const rvec<int>& gen_pdg,
                               const rvec<float>& gen_pt) {
   auto func = [] ( const double pt ) { return std::exp( 0.0615 - 0.0005 * pt ); };
   if (gen_pdg[2] == 6 && gen_pdg[3] == -6)
        return std::sqrt(func(gen_pt[2]) * func(gen_pt[3]));
   return 1.0;
};

// function for pdf weights
float get_pdf_weight(const unsigned int n_pdf_weight,
                     const rvec<float>& pdf_weight) {
    return pdf_weight[n_pdf_weight];
};

// apply kfactor for W background
float get_kfactor(	const rvec<int>& gen_pdg,
					const rvec<float>& gen_mass) {
	auto mass = 0.0;
	for (int i = 0; i < gen_pdg.size(); i++) {
		if (abs(gen_pdg[i]) == 24) {
			mass = gen_mass[i];
			break;
		}
	}
	
	return config::W_kfactor_hist->GetBinContent(config::W_kfactor_hist->GetXaxis()->FindBin(mass));
}
