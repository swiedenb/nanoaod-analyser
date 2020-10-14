#ifndef PhysicalQuantities_hh
#define PhysicalQuantities_hh

#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RVec.hxx>


template < typename T >
using rvec = ROOT::VecOps::RVec<T>;

// Calculate MT from one particle and MET
float mass_transv(	const float& part_pt, 
					const float& part_phi, 
					const float& met_pt, 
					const float& met_phi);

// Calculate MT from first (highest pt) particle, which fulfils mask, -  and MET
float mass_transv_masked(	const rvec<float>& part_pt, 
							const rvec<float>& part_phi, 
							const rvec<bool>& part_mask, 
							const float& met_pt, 
							const float& met_phi);

// Calculate invariant mass from two particles
float mass_inv(	const float& part1_pt,
				const float& part2_pt,
				const float& part1_eta,
				const float& part2_eta,
				const float& part1_phi,
				const float& part2_phi,
				const float& part1_mass,
				const float& part2_mass);

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
						const rvec<bool>& p2_mask);
						
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
                                const rvec<bool>& p2_mask);
// Calculate collinear mass from two highest pt particles, which fulfil mask and MET is projected onto the tau pT vector to calculate a rescale ratio for the tau 4 momentum
float collinear_mass_alt(           const rvec<int>& col_idx,
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
                                const rvec<bool>& p2_mask);
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
                                      const rvec<bool>& p2_mask) ;
// Returns ratio of two given quants
float ratio(const float& part1_pt, 
			const float& part2_pt);

// Returns ratio of two given quants
rvec<float> ratio_vector(	const rvec<float> & part1_qt, 
							const float& part2_qt);


// Returns Delta Phi between two particles
float delta_phi(const float& part1_phi,
				const float& part2_phi);

// Returns Delta Eta between two particles				
float delta_eta(	const float& part1_eta,
					const float& part2_eta);

// Returns abs of given number
float calc_abs(const float& numb);

#endif