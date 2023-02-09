#ifndef WEIGHTS_H
#define WEIGHTS_H

#include "config.hh"
#include "PhysicalQuantities.hh"

#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RVec.hxx>

template < typename T >
using rvec = ROOT::VecOps::RVec<T>;

float trigger_sf(   const rvec<float>& pt,
                        const rvec<float>& eta,
                        const rvec<bool>& mask,
                        const std::string& unc); 

float GetTopQscale(
                    const float& Mll,
                    const std::string& var);
float GetTopPDF(
                    const float& Mll,
                    const std::string& var);
float pu_weight( const int& nvtx_true,
                 const std::string& unc);

float muon_reco_eff( const rvec<float>& muon_pt,
                               const rvec<float>& muon_eta,
                               const rvec<float>& muon_phi,
                               const rvec<float>& muon_mass,
                               const rvec<bool>& muon_mask,
                               const std::string runtype); 

float prefire_factor(	const rvec<float>& jet_pt,
						const rvec<float>& jet_eta,
						const rvec<float>& photon_pt,
						const rvec<float>& photon_eta,
                        const std::string& unc);

float tau_fake_scale_factor(	const float& tau_pt,
                                const float& tau_eta,
								const UChar_t& tau_genPartFlav,
								const std::string& part,
                                const std::string& run_type);

float ele_id_scale_factor( const rvec<float>& ele_pt,
                            const rvec<float>& ele_eta,
                            const rvec<float>& ele_mask,
                            const std::string& run_type);
float dd_fakerate( const rvec<float>& tau_pt,
                            const rvec<float>& tau_eta,
                            const float & tau_pt_over_jet_pt,
                            const rvec<int>& col_idx,
                            const rvec<bool>& tau_mask,
                            const std::string& run_type);

float muon_id_scale_factor( const rvec<float>& muon_pt,
                            const rvec<float>& muon_eta,
                            const rvec<float>& muon_mask,
                            const std::string& run_type);

float muon_iso_scale_factor( const rvec<float>& muon_pt,
                            const rvec<float>& muon_eta,
                            const rvec<float>& muon_mask,
                            const std::string& run_type);
float calc_top_pt_reweighting( const rvec<int>& gen_pdg,
                               const rvec<float>& gen_pt);

float get_pdf_weight(const unsigned int n_pdf_weight,
                     const rvec<float>& pdf_weight);

float get_kfactor(	const rvec<int>& gen_pdg,
					const rvec<float>& gen_mass);                     
                     
#endif
