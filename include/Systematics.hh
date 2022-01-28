#ifndef Systematics_hh
#define Systematics_hh

#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RVec.hxx>
#include "config.hh"


template < typename T >
using rvec = ROOT::VecOps::RVec<T>;


rvec<float> muonResolutionSmearing(const rvec<float> pt, const rvec<float> eta , const std::string muon_shift);
float applyResolutionSmearing(const float pt, const float eta);
rvec<float> electronScaleUncertainty(const rvec<float>& pt, 
                               const rvec<float>& eta, 
                               const rvec<float>& phi, 
                               const rvec<float>& mass, 
                               const rvec<UChar_t>& gain, 
                               const rvec<float>& eCorr, 
                               const rvec<float>& r9, 
                               const rvec<float>& deltaSCEta,
                               const UInt_t& runnb,
                               const std::string& type,
                               const std::string& direction);

#endif
