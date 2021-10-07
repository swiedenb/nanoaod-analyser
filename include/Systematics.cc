#include "Systematics.hh"
#include "TLorentzVector.h"
#include "TRandom3.h"

#include <math.h>


rvec<float> muonResolutionSmearing(const rvec<float> pt, const rvec<float> eta, const std::string muon_shift){
    rvec<float> pt_smeared;
    for(unsigned int i = 0; i < pt.size(); i++){
        float ratio = applyResolutionSmearing(pt[i],eta[i]);
        if (muon_shift == "_muonResolutionUp"){
            pt_smeared.push_back(pt[i]*(1+ratio));
        }
        else if (muon_shift == "_muonResolutionDown"){
            pt_smeared.push_back(pt[i]*(1-ratio));
        }
    }
    return pt_smeared;

}

float applyResolutionSmearing(const float pt, const float eta){
// Apply Gaussian smearing to transverse momentum
     // Reference and explanation: https://indico.cern.ch/event/601852/
     // Determine value of smearing based on muon pT
     TRandom3* gauss = new TRandom3();
     double smearing = 0.0;
     if (pt < 200.) {
       smearing = 0.003;
     } else if (pt < 500.) {
       smearing = 0.005;
     } else {
       smearing = 0.01;
     }   
     // Double smearing for muons in the endcaps
     if (std::abs(eta) > 1.2) smearing *= 2;
     // Return Gaussian smearing
     double ratio = gauss->Gaus(0, smearing);
     return ratio;
    
}
