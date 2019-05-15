#ifndef WEIGHTS_H
#define WEIGHTS_H

#include "config.hh"

#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RVec.hxx>

template < typename T >
using rvec = ROOT::VecOps::RVec<T>;
float pu_weight( const float nvtx_true);
float apply_scale_factor(); 
float calc_top_pt_reweighting( const rvec<int>& gen_pdg,
                               const rvec<float>& gen_pt); 
#endif
