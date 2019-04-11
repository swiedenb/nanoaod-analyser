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
						
// Returns ratio of two given quants
float ratio(const float& part1_pt, 
			const float& part2_pt);


// Returns Delta Phi between two particles
float delta_phi(const float& part1_phi,
				const float& part2_phi);
				
// Returns abs of given number
float calc_abs(const float& numb);

#endif
