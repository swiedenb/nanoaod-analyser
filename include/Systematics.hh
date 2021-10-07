#ifndef Systematics_hh
#define Systematics_hh

#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RVec.hxx>


template < typename T >
using rvec = ROOT::VecOps::RVec<T>;


rvec<float> muonResolutionSmearing(const rvec<float> pt, const rvec<float> eta , const std::string muon_shift);
float applyResolutionSmearing(const float pt, const float eta);

#endif
