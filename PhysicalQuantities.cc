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

// Returns ratio of two given quants
float ratio(	const float& part1_qt, 
				const float& part2_qt) {
	return part1_qt / part2_qt;
};

// Returns Delta Phi between two particles
float delta_phi(	const float& part1_phi,
					const float& part2_phi) {
	double dphi = part1_phi - part2_phi;
	if (dphi > M_PI) dphi -= M_PI;
	else if (dphi < - M_PI) dphi += M_PI;
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
