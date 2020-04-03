#include "PhysicalQuantities.hh"
#include "TLorentzVector.h"

#include <math.h>

// Calculate MT from one particle and MET
float mass_transv(	const float& part_pt, 
					const float& part_phi, 
					const float& met_pt, 
					const float& met_phi) {
	return sqrt(2*part_pt*met_pt*(1 - cos( delta_phi(part_phi, met_phi) ) ) );
};

// Calculate MT from first (highest pt) particle, which fulfils mask, -  and MET
float mass_transv_masked(	const rvec<float>& part_pt, 
							const rvec<float>& part_phi, 
							const rvec<bool>& part_mask, 
							const float& met_pt, 
							const float& met_phi) {
	return sqrt(2*part_pt[part_mask][0]*met_pt*(1 - cos( delta_phi(part_phi[part_mask][0], met_phi) ) ) );
};

// Calculate invariant mass from two particles
float mass_inv(	const float& part1_pt,
				const float& part2_pt,
				const float& part1_eta,
				const float& part2_eta,
				const float& part1_phi,
				const float& part2_phi,
				const float& part1_mass,
				const float& part2_mass) {
    TLorentzVector p1, p2;
    p1.SetPtEtaPhiM(part1_pt, part1_eta, part1_phi, part1_mass);
    p2.SetPtEtaPhiM(part2_pt, part2_eta, part2_phi, part2_mass);
    return (p1 + p2).M();
};

// Calculate invariant mass from two highest pt particles, which fulfil mask
float mass_inv_masked(	const rvec<float>& p1_pt,
						const rvec<float>& p2_pt,
						const rvec<float>& p1_eta,
						const rvec<float>& p2_eta,
						const rvec<float>& p1_phi,
						const rvec<float>& p2_phi,
						const rvec<float>& p1_mass,
						const rvec<float>& p2_mass,
						const rvec<bool>& p1_mask,
						const rvec<bool>& p2_mask) {
    TLorentzVector p1, p2;
    p1.SetPtEtaPhiM(p1_pt[p1_mask][0], p1_eta[p1_mask][0], p1_phi[p1_mask][0], p1_mass[p1_mask][0]);
    p2.SetPtEtaPhiM(p2_pt[p2_mask][0], p2_eta[p2_mask][0], p2_phi[p2_mask][0], p2_mass[p2_mask][0]);
    return (p1 + p2).M();
};

// Calculate collinear mass from two highest pt particles, which fulfil mask and MET with the eta from the first particle
float collinear_mass(           const rvec<int>& col_idx,
                                const rvec<float>& p1_pt,
                                const rvec<float>& p2_pt,
                                const rvec<float>& p1_eta,
                                const rvec<float>& p2_eta,
                                const rvec<float>& p1_phi,
                                const rvec<float>& p2_phi,
                                const rvec<float>& p1_mass,
                                const rvec<float>& p2_mass,
                                const float& met_pt,
                                const float& met_phi,
                                const rvec<bool>& p1_mask,
                                const rvec<bool>& p2_mask){
    TLorentzVector p1, p2, p3, p_best;
    p1.SetPtEtaPhiM(p1_pt[p1_mask][col_idx[0]], p1_eta[p1_mask][col_idx[0]], p1_phi[p1_mask][col_idx[0]], p1_mass[p1_mask][col_idx[0]]);
    p2.SetPtEtaPhiM(p2_pt[p2_mask][col_idx[1]], p2_eta[p2_mask][col_idx[1]], p2_phi[p2_mask][col_idx[1]], p2_mass[p2_mask][col_idx[1]]);
    p3.SetPtEtaPhiM(met_pt, p1_eta[p1_mask][col_idx[0]], met_phi, 0.);
      
    return (p1 + p2 + p3).M(); 
};
// Calculate collinear mass from two highest pt particles, which fulfil mask and MET is projected onto the tau pT vector to calculate a rescale ratio for the tau 4 momentum
float collinear_mass_alt(       const rvec<int>& col_idx,
                                const rvec<float>& p1_pt,
                                const rvec<float>& p2_pt,
                                const rvec<float>& p1_eta,
                                const rvec<float>& p2_eta,
                                const rvec<float>& p1_phi,
                                const rvec<float>& p2_phi,
                                const rvec<float>& p1_mass,
                                const rvec<float>& p2_mass,
                                const float& met_pt,
                                const float& met_phi,
                                const rvec<bool>& p1_mask,
                                const rvec<bool>& p2_mask){
    TLorentzVector p_tau, p2, p_met;
    p_tau.SetPtEtaPhiM(p1_pt[p1_mask][col_idx[0]], p1_eta[p1_mask][col_idx[0]], p1_phi[p1_mask][col_idx[0]], p1_mass[p1_mask][col_idx[0]]);
    p2.SetPtEtaPhiM(p2_pt[p2_mask][col_idx[1]], p2_eta[p2_mask][col_idx[1]], p2_phi[p2_mask][col_idx[1]], p2_mass[p2_mask][col_idx[1]]);
    p_met.SetPtEtaPhiM(met_pt, p1_eta[p1_mask][col_idx[0]], met_phi, 0.);
    double x = p_tau.Pt() / (p_tau.Pt() + met_pt*cos(delta_phi(p_tau.Phi(),met_phi))); 

    p_tau.SetPxPyPzE(p_tau.Px()/x,p_tau.Py()/x,p_tau.Pz()/x,p_tau.E()/x);
      
    return (p_tau + p2).M(); 
};
rvec<int> col_idx(                    const rvec<float>& p1_pt,
                                      const rvec<float>& p2_pt,
                                      const rvec<float>& p1_eta,
                                      const rvec<float>& p2_eta,
                                      const rvec<float>& p1_phi,
                                      const rvec<float>& p2_phi,
                                      const rvec<float>& p1_mass,
                                      const rvec<float>& p2_mass,
                                      const float& met_pt,
                                      const float& met_phi,
                                      const rvec<bool>& p1_mask,
                                      const rvec<bool>& p2_mask) {
    TLorentzVector p1, p2, p3, p_best;
    rvec<int> idx(2);
    float highest_mass = 0;
    int best_i = 0;
    int best_j = 0;
    for (unsigned int i = 0; i < p1_pt[p1_mask].size(); i++){
        p3.SetPtEtaPhiM(met_pt, p1_eta[p1_mask][i], met_phi, 0.);
        p1.SetPtEtaPhiM(p1_pt[p1_mask][i], p1_eta[p1_mask][i], p1_phi[p1_mask][i], p1_mass[p1_mask][i]);
        for(unsigned int j = 0; j < p2_pt[p2_mask].size(); j++){
            p2.SetPtEtaPhiM(p2_pt[p2_mask][j], p2_eta[p2_mask][j], p2_phi[p2_mask][j], p2_mass[p2_mask][j]);
            float mass = (p1 + p2 + p3).M();
            if ( highest_mass < mass){
                highest_mass = mass;
                best_i = i;
                best_j = j;
            }
        }
    }
    idx[0] = best_i;
    idx[1] = best_j;
    return idx; 
};

// Returns ratio of two given quants
float ratio(	const float& part1_qt, 
				const float& part2_qt) {
	return part1_qt / part2_qt;
};

// Returns ratio of all valeus in vector divided by quant
rvec <float> ratio_vector(	const rvec<float> & part1_qt, 
							const float& part2_qt) {
	rvec<float> result_vec(part1_qt.size());
	for (uint i = 0; i < part1_qt.size(); i++) {
		result_vec[i] = part1_qt[i] / part2_qt;
	}
	return result_vec;
};

// Returns Delta Phi between two particles
float delta_phi(	const float& part1_phi,
					const float& part2_phi) {
	double dphi = part1_phi - part2_phi;
	if (dphi > M_PI) dphi -= 2*M_PI;
	else if (dphi < - M_PI) dphi += 2*M_PI;
	return dphi;
};

// Returns Delta Eta between two particles
float delta_eta(	const float& part1_eta,
					const float& part2_eta) {
	return abs(part1_eta - part2_eta);
};

// Returns abs of given number
float calc_abs(const float& numb) {
	return std::abs(numb);
};
