#include "Particles.hh"



// check if tau is in acceptance and fulfils id requirements
rvec<bool> tau_acceptance_and_id(const rvec<float>& pt, 
                                 const rvec<float>& eta,
                                 const rvec<bool>& dm,
                                 const rvec<UChar_t>& iso, 
                                 const rvec<UChar_t>& antiEle_disc, 
                                 const rvec<UChar_t>& antiMu_disc) {
        // tau pt cut
        auto mask_pt = pt > config::tau_pt;
        
        // tau eta cut
        auto mask_eta = abs(eta) < config::tau_eta;
        
        // tau decay mode
        auto mask_dm = dm;
                
        // tau iso requirement (1: VVLoose, 2: VLoose, 4: Loose, 8: Medium, 16: Tight, 32: VTight, 64: VVTight)
        auto mask_iso = (iso & config::tau_iso_WP) == config::tau_iso_WP;
        
        // Anti Electron discriminator (1: VLoose, 2: Loose, 4: Medium, 8: Tight, 16: VTight)
        auto mask_antiele = (config::tau_antiE_WP & antiEle_disc) == config::tau_antiE_WP;
        
        // Anti Muon discriminator (1: Loose, 2: Tight)
        auto mask_antimu = (config::tau_antiMu_WP & antiMu_disc) == config::tau_antiMu_WP;
        
        // return vector with true, if particle fulfils all requirements - else false
        return mask_pt & mask_eta & mask_dm & mask_iso & mask_antiele & mask_antimu;
};

// check if tau is in acceptance and fulfils id requirements with certain dm vetos
rvec<bool> tau_acceptance_and_id_and_dm(const rvec<float>& pt, 
                                        const rvec<float>& eta,
                                        const rvec<float>& dz,
                                        const rvec<bool>& dm,
                                        const rvec<int>& dm_number,
                                        const rvec<UChar_t>& iso, 
                                        const rvec<UChar_t>& antiEle_disc, 
                                        const rvec<UChar_t>& antiMu_disc) {
        // tau pt cut
        auto mask_pt = pt > config::tau_pt;
        
        // tau eta cut
        auto mask_eta = abs(eta) < config::tau_eta;
        
        // tau dz cut
        auto mask_dz = abs(dz) < 0.2;
        // tau decay mode
        auto mask_dm = dm;

        // tau decay mode as integer
        auto mask = not (dm_number == 5 || dm_number == 6);
                
        // tau iso requirement (1: VVLoose, 2: VLoose, 4: Loose, 8: Medium, 16: Tight, 32: VTight, 64: VVTight)
        auto mask_iso = (iso & config::tau_iso_WP) == config::tau_iso_WP;
        
        // Anti Electron discriminator (1: VLoose, 2: Loose, 4: Medium, 8: Tight, 16: VTight)
        auto mask_antiele = (config::tau_antiE_WP & antiEle_disc) == config::tau_antiE_WP;
        
        // Anti Muon discriminator (1: Loose, 2: Tight)
        auto mask_antimu = (config::tau_antiMu_WP & antiMu_disc) == config::tau_antiMu_WP;
        
        // return vector with true, if particle fulfils all requirements - else false
        return mask_pt & mask_eta & mask_dz & mask & mask_dm & mask_iso & mask_antiele & mask_antimu;
};

// check if event has a well seperated di muon pair
rvec<bool> di_muon_id(const rvec<float>& pt, 
                                  const rvec<float>& eta, 
                                  const rvec<UChar_t>& id, 
                                  const rvec<float>& iso) {
        // muon pt cut
        auto mask_pt = pt > 10.;
        
        // muon eta cut
        auto mask_eta = abs(eta) < config::muon_eta;
        
        // muon id
        auto mask_id = (id == config::muon_id_WP);
    
        // muon tracker iso 
        auto mask_iso = iso / pt < 0.15;

        
        // return vector with true, if particle fulfils all requirements - else false
        return mask_pt & mask_eta & mask_id & mask_iso;
};

// check if muon is in acceptance and fulfils id requirements
rvec<bool> muon_acceptance_and_id(const rvec<float>& pt, 
                                  const rvec<float>& eta, 
                                  const rvec<UChar_t>& id, 
                                  const rvec<float>& iso) {
        // muon pt cut
        auto mask_pt = pt > config::muon_pt;
        
        // muon eta cut
        auto mask_eta = abs(eta) < config::muon_eta;
        
        // muon id
        auto mask_id = (id == config::muon_id_WP);
    
        // muon relative iso 
        auto mask_iso = (1 <= config::muon_iso_WP);

        
        // return vector with true, if particle fulfils all requirements - else false
        return mask_pt & mask_eta & mask_id & mask_iso;
};

// check if electron is in acceptance and fulfils id requirements
rvec<bool> ele_acceptance_and_id(const rvec<float>& pt,
                                 const rvec<float>& eta,
                                 const rvec<Int_t>& id) {
        // electron pt cut
        auto mask_pt = pt > config::ele_pt;
        
        // electron eta cut
        auto mask_eta = abs(eta) < config::ele_eta;
        
        // electron id cut
        auto mask_id = id >= 2;
        
        // return vector with true, if particle fulfils all requirements - else false
        return mask_pt & mask_eta & mask_id;
};


rvec<bool> ele_acceptance_and_simpleid(const rvec<float>& pt,
                                       const rvec<float>& eta,
                                       const rvec<bool>& id) {
        // electron pt cut
        auto mask_pt = pt > config::ele_pt;
        
        // electron eta cut
        auto mask_eta = abs(eta) < config::ele_eta;
        
        // electron id cut
        auto mask_id = id;
        
        // return vector with true, if particle fulfils all requirements - else false
        return mask_pt & mask_eta & mask_id;
};


// calc tau energy scale
rvec< rvec< float > > calc_tau_energy_scale(const rvec<int>& tau_decayMode,
                                            const rvec<UChar_t>& tau_genPartFlav,
                                            const rvec<float>& tau_pt,
                                            const rvec<float>& tau_eta,
                                            const rvec<float>& tau_phi,
                                            const rvec<float>& tau_mass) {  
        rvec<rvec<float>> binder(4);
        
        rvec<float> new_tau_pt(tau_pt.size());                          
        rvec<float> new_tau_eta(tau_eta.size());                                
        rvec<float> new_tau_phi(tau_phi.size());                                
        rvec<float> new_tau_mass(tau_mass.size());                              
        
        for (uint i = 0; i<tau_pt.size(); i++) {

            TLorentzVector p1; 
            p1.SetPtEtaPhiM(tau_pt[i], tau_eta[i], tau_phi[i], tau_mass[i]);

            if (config::run_type == "_TauEnergyScaleUp")
                p1 *= config::tau_dm_scale->getTES_highpt(tau_pt[i], tau_decayMode[i], (int) tau_genPartFlav[i], "Up");
            else if (config::run_type == "Tau_EnergyScaleDown")
                p1 *= config::tau_dm_scale->getTES_highpt(tau_pt[i], tau_decayMode[i], (int) tau_genPartFlav[i], "Down");
            else
                p1 *= config::tau_dm_scale->getTES_highpt(tau_pt[i], tau_decayMode[i], (int) tau_genPartFlav[i], "");
                 
            new_tau_pt[i] = p1.Pt();
            new_tau_eta[i] = p1.Eta();
            new_tau_phi[i] = p1.Phi();
            new_tau_mass[i] = p1.M();
        }
        
        binder[0] = new_tau_pt;                         
        binder[1] = new_tau_eta;                                
        binder[2] = new_tau_phi;                                
        binder[3] = new_tau_mass;
        
        return binder;
}


// apply tau energy scale
rvec <float> apply_tau_energy_scale(const rvec< rvec<float> > binder, const int& index) {
    return binder[index];
}
