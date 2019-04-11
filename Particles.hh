#ifndef Particles_hh
#define Particles_hh

#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RVec.hxx>

#include "config.hh"

template < typename T >
using rvec = ROOT::VecOps::RVec<T>;

rvec<bool> tau_acceptance_and_id(	const rvec<float>& pt, 
									const rvec<float>& eta, 
									const rvec<UChar_t>& iso, 
									const rvec<UChar_t>& antiEle_disc, 
									const rvec<UChar_t>& antiMu_disc);


rvec<bool> muon_acceptance_and_id(	const rvec<float>& pt, 
									const rvec<float>& eta, 
									const rvec<bool>& id, 
									const rvec<float>& iso
									);

rvec<bool> ele_acceptance_and_id(	const rvec<float>& pt,
									const rvec<float>& eta,
									const rvec<Int_t>& id);

#endif
