#ifndef Particles_hh
#define Particles_hh

#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TLorentzVector.h>

#include "config.hh"

template < typename T >
using rvec = ROOT::VecOps::RVec<T>;
using RNode = ROOT::RDF::RNode;


rvec<bool> tau_acceptance_and_id(	const rvec<float>& pt, 
                                    const rvec<float>& eta, 
                                    const rvec<float>& dz,
                                    const rvec<int>& charge, 
                                    const rvec<bool> & dm,
                                    const rvec<UChar_t>& iso, 
                                    const rvec<UChar_t>& antiEle_disc, 
                                    const rvec<UChar_t>& antiMu_disc);



rvec<bool> tau_acceptance_and_id_and_dm(	const rvec<float>& pt, 
                                            const rvec<float>& eta, 
                                            const rvec<float>& dz,
                                            const rvec<int>& charge,
                                            const rvec<bool>& dm,
                                            const rvec<int>& dm_number, 
                                            const rvec<UChar_t>& iso, 
                                            const rvec<UChar_t>& antiEle_disc, 
                                            const rvec<UChar_t>& antiMu_disc);



rvec<bool> muon_acceptance_and_id(	const rvec<float>& pt, 
									const rvec<float>& eta, 
									const rvec<UChar_t>& id, 
									const rvec<UChar_t>& iso
									);

rvec<bool> di_muon_id(const rvec<float>& pt, 
                                  const rvec<float>& eta, 
                                  const rvec<UChar_t>& id, 
                                  const rvec<float>& iso) ;
rvec<bool> ele_acceptance_and_id(	const rvec<float>& pt,
									const rvec<float>& eta,
									const rvec<Int_t>& id);
									
rvec<bool> ele_acceptance_and_simpleid(	const rvec<float>& pt,
										const rvec<float>& eta,
										const rvec<bool>& id);
									
rvec< rvec< float > > calc_tau_energy_scale(const rvec<int>& tau_decayMode,
                                            const rvec<UChar_t>& tau_genPartFlav,
                                            const rvec<float>& tau_pt,
                                            const rvec<float>& tau_eta,
                                            const rvec<float>& tau_phi,
                                            const rvec<float>& tau_mass,
                                            std::string type);

rvec <float> apply_tau_energy_scale(const rvec< rvec<float> > binder, 
                                    const int& index);

#endif
