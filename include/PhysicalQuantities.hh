#ifndef PhysicalQuantities_hh
#define PhysicalQuantities_hh

#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RVec.hxx>


#include "config.hh"
template < typename T >
using rvec = ROOT::VecOps::RVec<T>;

// Calculate sphericity from one particle class
float sphericity(	const rvec<float>& part_pt, 
					const rvec<float>& part_eta, 
					const rvec<float>& part_phi, 
					const rvec<float>& part_mass,
                    const rvec<bool>& mask); 
float sphericity_full(	const rvec<float>& part1_pt, 
					const rvec<float>& part1_eta, 
					const rvec<float>& part1_phi, 
					const rvec<float>& part1_mass,
                    const rvec<float>& part2_pt, 
					const rvec<float>& part2_eta, 
					const rvec<float>& part2_phi, 
					const rvec<float>& part2_mass,
                    const rvec<float>& part3_pt, 
					const rvec<float>& part3_eta, 
					const rvec<float>& part3_phi, 
					const rvec<float>& part3_mass,
                    const rvec<bool>&  mask,
                    const rvec<bool>&  mask2);
float sphericity_leptons(	const rvec<float>& part_pt, 
					const rvec<float>& part_eta, 
					const rvec<float>& part_phi, 
					const rvec<float>& part_mass,
					const rvec<float>& part2_pt, 
					const rvec<float>& part2_eta, 
					const rvec<float>& part2_phi, 
					const rvec<float>& part2_mass,
					const rvec<bool>& mask,
                    const rvec<bool>& mask2);
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

float lepton_inv_mass_dy(const rvec<float>& lhept,
                         const rvec<float>& lheeta,
                         const rvec<float>& lhephi,
                         const rvec<float>& lhemass,
                         const rvec<int>& lhepdgid);
float lepton_inv_mass(const rvec<int>& pdgID,
                      const rvec<float>& mass,
                      const rvec<float>& pt,
                      const rvec<float>& eta,
                      const rvec<float>& phi,
                      const rvec<int>& mother_idx);
float lepton_inv_mass_ww(const rvec<int>& pdgID,
                      const rvec<float>& mass,
                      const rvec<float>& pt,
                      const rvec<float>& eta,
                      const rvec<float>& phi,
                      const rvec<int>& status);
// Calculate invariant mass from two particles
float mass_inv(	const float& part1_pt,
				const float& part2_pt,
				const float& part1_eta,
				const float& part2_eta,
				const float& part1_phi,
				const float& part2_phi,
				const float& part1_mass,
				const float& part2_mass);

float mass_inv_idx(	const rvec<float>& pt,
						const rvec<float>& eta,
						const rvec<float>& phi,
						const rvec<float>& mass,
						const rvec<int>& idx); 
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
						
float sel_tau_pt(               const rvec<float>& tau_pt,
                                const rvec<bool>& tau_mask,
                                const rvec<int>& col_idx);
float tau_pt_over_jet_pt(       const rvec<int>& col_idx,
                                const rvec<float>& tau_pt,
                                const rvec<int>& jet_idx,
                                const rvec<float>& jet_pt,
                                const rvec<bool>& tau_mask);
float tau_pt_over_jet_pt_barrel(       const rvec<int>& col_idx,
                                const rvec<float>& tau_pt,
                                const rvec<int>& jet_idx,
                                const rvec<float>& jet_pt,
                                const rvec<bool>& tau_mask,
                                const rvec<float>& tau_eta);
float tau_pt_over_jet_pt_endcap(       const rvec<int>& col_idx,
                                const rvec<float>& tau_pt,
                                const rvec<int>& jet_idx,
                                const rvec<float>& jet_pt,
                                const rvec<bool>& tau_mask,
                                const rvec<float>& tau_eta);
// Calculate collinear mass from two highest pt particles, which fulfil mask and MET with the eta from the first particle
float collinear_mass(           const float& p1_pt,
                                const float& p2_pt,
                                const float& p1_eta,
                                const float& p2_eta,
                                const float& p1_phi,
                                const float& p2_phi,
                                const float& p1_mass,
                                const float& p2_mass,
                                const float& met_pt,
                                const float& met_phi);
//float collinear_mass(           const rvec<int>& col_idx,
//                                const rvec<float>& p1_pt,
//                                const rvec<float>& p2_pt,
//                                const rvec<float>& p1_eta,
//                                const rvec<float>& p2_eta,
//                                const rvec<float>& p1_phi,
//                                const rvec<float>& p2_phi,
//                                const rvec<float>& p1_mass,
//                                const rvec<float>& p2_mass,
//                                const float& met_pt,
//                                const float& met_phi,
//                                const rvec<bool>& p1_mask,
//                                const rvec<bool>& p2_mask);
// Calculate collinear mass from two highest pt particles, which fulfil mask and MET is projected onto the tau pT vector to calculate a rescale ratio for the tau 4 momentum
float collinear_mass_alt_2(       const float& p1_pt,
                                const float& p2_pt,
                                const float& p1_eta,
                                const float& p2_eta,
                                const float& p1_phi,
                                const float& p2_phi,
                                const float& p1_mass,
                                const float& p2_mass,
                                const float& met_pt,
                                const float& met_phi);
float collinear_mass_alt(       const float& p1_pt,
                                const float& p2_pt,
                                const float& p1_eta,
                                const float& p2_eta,
                                const float& p1_phi,
                                const float& p2_phi,
                                const float& p1_mass,
                                const float& p2_mass,
                                const float& met_pt,
                                const float& met_phi);
//float collinear_mass_alt(           const rvec<int>& col_idx,
//                                const rvec<float>& p1_pt,
//                                const rvec<float>& p2_pt,
//                                const rvec<float>& p1_eta,
//                                const rvec<float>& p2_eta,
//                                const rvec<float>& p1_phi,
//                                const rvec<float>& p2_phi,
//                                const rvec<float>& p1_mass,
//                                const rvec<float>& p2_mass,
//                                const float& met_pt,
//                                const float& met_phi,
//                                const rvec<bool>& p1_mask,
//                                const rvec<bool>& p2_mask);
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
