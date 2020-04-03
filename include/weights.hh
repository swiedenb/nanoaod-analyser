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


float pu_weight( const float& nvtx_true,
                 const std::string& unc);


float prefire_factor(	const rvec<float>& jet_pt,
						const rvec<float>& jet_eta,
						const rvec<float>& photon_pt,
						const rvec<float>& photon_eta);

float tau_fake_scale_factor(    const rvec<float>& tau_pt,
                                const rvec<float>& tau_eta,
                                const rvec<bool>& tau_mask,
                                const rvec<UChar_t>& tau_genPartFlav,
                                const std::string& part,
                                const std::string& run_type);

float ele_id_scale_factor( const rvec<float>& ele_pt,
                            const rvec<float>& ele_eta,
                            const rvec<float>& ele_mask,
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
