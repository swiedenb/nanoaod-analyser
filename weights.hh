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


float pu_weight( const float nvtx_true);

float apply_scale_factor();

float tau_fake_scale_factor(	const float& tau_eta,
								const float& tau_phi,
								const rvec<int>& gen_pdgID,
								const rvec<float>& gen_eta,
								const rvec<float>& gen_phi,
								const int& required_pdgID);
								
float calc_top_pt_reweighting( const rvec<int>& gen_pdg,
                               const rvec<float>& gen_pt);

float get_pdf_weight(const unsigned int n_pdf_weight,
                     const rvec<float>& pdf_weight);

float get_kfactor(	const rvec<int>& gen_pdg,
					const rvec<float>& gen_mass);                     
                     
#endif
