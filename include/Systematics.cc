#include "Systematics.hh"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "EnergyScaleCorrection.h"

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
     delete gauss;
     return ratio;
    
}
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
                               const std::string& direction){
    TRandom3* gauss = new TRandom3();
    double nrandom = gauss->Gaus(0, 1);
    TLorentzVector p, praw;
    rvec<float> pt_unc;
    rvec<float> mass_unc;
    for(int i = 0; i < pt.size() ; i ++){
        p.SetPtEtaPhiM(pt[i], eta[i], phi[i], mass[i]);
        praw = 1./eCorr[i] * p;
        auto et = praw.Et();
        auto abseta = abs(praw.Eta() + deltaSCEta[i]);
        auto eleSmear    = config::eleCorr->smearingSigma(runnb, et, abseta, r9[i], 12, 0, 0.);
        auto escaleErr   = config::eleCorr->scaleCorrUncert(runnb, et, abseta, r9[i]);
        auto eleSmearUp  = config::eleCorr->smearingSigma(runnb, et, abseta, r9[i], 12,  1, 0.);
        auto eleSmearDo  = config::eleCorr->smearingSigma(runnb, et, abseta, r9[i], 12, -1, 0.);
        auto eleSmearUnc = nrandom * sqrt( (eleSmearUp-eleSmear)*(eleSmearUp-eleSmear) + (eleSmearDo-eleSmear)*(eleSmearDo-eleSmear) );
        auto vEleUp = (1+escaleErr+eleSmearUnc)* p;
        auto vEleDo = (1-escaleErr-eleSmearUnc)* p;
        if (direction == "Up"){
            if (type == "mass") mass_unc.push_back(vEleUp.M());
            if (type == "pt")   pt_unc.push_back(vEleUp.Pt());
        }
        else if (direction == "Down"){
            if (type == "mass") mass_unc.push_back(vEleDo.M());
            if (type == "pt")  pt_unc.push_back(vEleDo.Pt());
        }
    }
    delete gauss;
    if (type == "mass") return mass_unc;
    if (type == "pt") return pt_unc;
    return pt;
                               }
