#include "weights.hh"

// Get Pileup weight from file
float pu_weight( const float& nvtx_true, const std::string& unc = ""){
    if (unc == "")
        return config::pileup_hist->GetBinContent( config::pileup_hist->FindBin(nvtx_true) );
    else if (unc == "Up")
        return config::pileup_hist_up->GetBinContent( config::pileup_hist->FindBin(nvtx_true) );
    else if (unc == "Down")
        return config::pileup_hist_down->GetBinContent( config::pileup_hist->FindBin(nvtx_true) );
};


// calculate prefiring weight for the corresponding years
float prefire_factor(	const rvec<float>& jet_pt,
						const rvec<float>& jet_eta,
						const rvec<float>& photon_pt,
						const rvec<float>& photon_eta) {
	if (config::era == 2018)
		return 1.0;
	float weight = 1.0;
	for (uint i = 0; i < photon_pt.size(); i++) {
		weight *= 1 - config::prefire_photon_hist->GetBinContent(
								config::prefire_photon_hist->GetXaxis()->FindBin(photon_eta[i]), 
								config::prefire_photon_hist->GetYaxis()->FindBin(photon_pt[i]) );
	}
	for (uint i = 0; i < photon_pt.size(); i++) {
		weight *= 1 - config::prefire_jet_hist->GetBinContent(
								config::prefire_jet_hist->GetXaxis()->FindBin(jet_eta[i]), 
								config::prefire_jet_hist->GetYaxis()->FindBin(jet_pt[i]) );
	}
	return weight;		
};

// get ele id scale factor
float ele_id_scale_factor( const rvec<float>& ele_pt,
                            const rvec<float>& ele_eta,
                            const rvec<float>& ele_mask,
                            const std::string& run_type) {

    float scale_factor = 1.0;
    for (uint i = 0; i < ele_pt[ele_mask].size(); i++) {
        scale_factor *= config::ele_SF->getSFID(ele_pt[ele_mask][i],ele_eta[ele_mask][i], run_type);
    return scale_factor;
    }

}
// get muon id scale factor
float muon_id_scale_factor( const rvec<float>& muon_pt,
                            const rvec<float>& muon_eta,
                            const rvec<float>& muon_mask,
                            const std::string& run_type) {

    float scale_factor = 1.0;
    for (uint i = 0; i < muon_pt[muon_mask].size(); i++) {
        scale_factor *= config::muon_SF->getSFID(muon_pt[muon_mask][i],muon_eta[muon_mask][i], run_type);
    return scale_factor;
    }

}
float muon_iso_scale_factor( const rvec<float>& muon_pt,
                            const rvec<float>& muon_eta,
                            const rvec<float>& muon_mask,
                            const std::string& run_type) {

    float scale_factor = 1.0;
    for (uint i = 0; i < muon_pt[muon_mask].size(); i++) {
        scale_factor *= config::muon_SF->getSFISO(muon_pt[muon_mask][i],muon_eta[muon_mask][i], run_type);
    return scale_factor;
    }

}


// get tau fake rate scale factor (for both e an mu)
float tau_fake_scale_factor( const rvec<float>& tau_pt,
                                const rvec<float>& tau_eta,
                                const rvec<bool>& tau_mask,
                                const rvec<UChar_t>& tau_genPartFlav,
                                const std::string& part,
                                const std::string& run_type) {
    
    float scale_factor = 1.0;
    
    
    for (uint i = 0; i < tau_genPartFlav[tau_mask].size(); i++) {
        if (part == "Ele") {                                      // there are outgoing gen_electrons
            if ((tau_genPartFlav[tau_mask][i] & 1 == 1) || (tau_genPartFlav[tau_mask][i] & 3 == 3)) {
                scale_factor *= config::tau_vse_SF->getSFvsEta(tau_eta[tau_mask][i], tau_genPartFlav[tau_mask][i], run_type);
            }
        }
        else if (part == "Muon") {                                 // these are outgoing gen_muons
            if ((tau_genPartFlav[tau_mask][i] & 2 == 2) || (tau_genPartFlav[tau_mask][i] & 4 == 4)) {
                scale_factor *= config::tau_vsmu_SF->getSFvsEta(tau_eta[tau_mask][i], tau_genPartFlav[tau_mask][i], run_type);
            }
        }
        else if (part == "Jet") {                                  // these are outgoing hadronic taus
            if ((tau_genPartFlav[tau_mask][i] & 5 == 5)) {
                scale_factor *= config::tau_vsjet_SF->getSFvsPT(tau_pt[tau_mask][i], run_type);
            }
        }
        else {                                                                                          // all the rest
            continue;
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
                        
    if (config::W_kfactor_hist != NULL)  {
        auto mass = 0.0;
        for (int i = 0; i < gen_pdg.size(); i++) {
            if (abs(gen_pdg[i]) == 24) {
                mass = gen_mass[i];
                break;
            }
        }    
        return config::W_kfactor_hist->GetBinContent(config::W_kfactor_hist->GetXaxis()->FindBin(mass));
    }
    else
        return 1.0;
}
