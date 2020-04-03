#ifndef PDFTOOL_h
#define PDFTOOL_h

#include "LHAPDF/LHAPDF.h"

#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RVec.hxx>

#include "config.hh"

template < typename T >
using rvec = ROOT::VecOps::RVec<T>;
using RNode = ROOT::RDF::RNode;


rvec<float> calc_pdf_weights(   const float& q2scale,
                                const float& x1,
                                const float& x2,
                                const int& id1,
                                const int& id2);
                                    
std::pair<float, float> calc_as_weights(   const float& q2scale,
                                           const float& x1,
                                           const float& x2,
                                           const int& id1,
                                           const int& id2);

       
RNode init_PDFs(RNode& df);  
    

#endif
