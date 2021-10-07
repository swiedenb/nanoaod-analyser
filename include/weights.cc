#include "weights.hh"
#include "TLorentzVector.h"

// Get Pileup weight from file
float pu_weight( const int& nvtx_true, const std::string& unc = ""){
    if (unc == "")
        return config::pileup_hist->GetBinContent( config::pileup_hist->FindBin(nvtx_true) );
    else if (unc == "Up")
        return config::pileup_hist_up->GetBinContent( config::pileup_hist->FindBin(nvtx_true) );
    else if (unc == "Down")
        return config::pileup_hist_down->GetBinContent( config::pileup_hist->FindBin(nvtx_true) );
};


// calculate trigger weight for the corresponding years
float trigger_sf(   const rvec<float>& pt,
                    const rvec<float>& eta,
                    const rvec<bool>& mask,
                    const std::string& unc = "") {
    float weight = 1.0;
    
    auto get_weight_from_hist = []( TH2D* hist, 
                                    TH2D* hist_up,
                                    TH2D* hist_down,
                                    const rvec<float>& pt, 
                                    const rvec<float>& eta, 
                                    const rvec<bool>& mask,
                                    const std::string& unc ) {
        float weight = 1.0;
        for (uint i = 0; i < pt[mask].size(); i++) {
            int eta_bin = hist->GetYaxis()->FindBin(abs(eta[mask][i]));
            int pt_bin = hist->GetXaxis()->FindBin(pt[mask][i]);
            int lasteta_bin = hist->GetYaxis()->FindBin(500);
            int lastpt_bin = hist->GetXaxis()->FindBin(999.);
            if (std::abs(eta[i]) < 5) {
                if (unc == "") {
                        if (pt[mask][i] > 1000.){
                            weight *= hist->GetBinContent(lastpt_bin, eta_bin);
                        }
                        else{
                            weight *= hist->GetBinContent(pt_bin, eta_bin);
                        }
                } else if (unc == "Up") {
                        if (pt[mask][i] > 1000.){
                            weight *= ( hist->GetBinContent(lastpt_bin,eta_bin) + hist_up->GetBinContent(lastpt_bin,eta_bin));
                        }
                        else{
                            weight *= ( hist->GetBinContent(pt_bin, eta_bin)
                                        + hist_up->GetBinContent(pt_bin, eta_bin) );
                        }
                } else if (unc == "Down") {
                        if (pt[mask][i] > 1000.){
                            weight *=  ( hist->GetBinContent(lastpt_bin, eta_bin) - hist_down->GetBinContent(lastpt_bin, eta_bin));
                        }
                        else{
                            weight *= ( hist->GetBinContent(pt_bin, eta_bin)
                                            - hist_down->GetBinContent(pt_bin, eta_bin) );
                        }
                }
            }
        }
        return weight;
    };
    
    weight *= get_weight_from_hist(config::trigger_hist, config::trigger_hist_up, config::trigger_hist_down, pt, eta, mask, unc);
    return weight;      
};
// calculate prefiring weight for the corresponding years
float prefire_factor(   const rvec<float>& jet_pt,
                        const rvec<float>& jet_eta,
                        const rvec<float>& photon_pt,
                        const rvec<float>& photon_eta,
                        const std::string& unc = "") {
    if (config::era == 2018) {
        return 1.0;
    }
    float weight = 1.0;
    
    auto get_weight_from_hist = []( TH2D* hist, 
                                    const rvec<float>& pt, 
                                    const rvec<float>& eta, 
                                    const std::string& unc ) {
        float weight = 1.0;
        for (uint i = 0; i < pt.size(); i++) {
            int eta_bin = hist->GetXaxis()->FindBin(eta[i]);
            int pt_bin = hist->GetYaxis()->FindBin(pt[i]);
            int lastpt_bin = hist->GetYaxis()->FindBin(500);
            if (std::abs(eta[i]) < 5) {
                if (unc == "") {
                    if (pt[i] <= 500) {
                        weight *= 1 - hist->GetBinContent(eta_bin, pt_bin);
                    } else {
                        weight *= 1 - hist->GetBinContent(eta_bin, lastpt_bin);
                    }
                } else if (unc == "Up") {
                    if (pt[i] <= 500) {
                        weight *= 1 - ( hist->GetBinContent(eta_bin, pt_bin)
                                        + hist->GetBinErrorUp(eta_bin, pt_bin) );
                    } else {
                        weight *= 1 - ( hist->GetBinContent(eta_bin, lastpt_bin)
                                        + hist->GetBinErrorUp(eta_bin, lastpt_bin) );
                    }
                } else if (unc == "Down") {
                    if (pt[i] <= 500) {
                        weight *= 1 - ( hist->GetBinContent(eta_bin, pt_bin)
                                        - hist->GetBinErrorLow(eta_bin, pt_bin) );
                    } else {
                        weight *= 1 - ( hist->GetBinContent(eta_bin, lastpt_bin)
                                        - hist->GetBinErrorLow(eta_bin, lastpt_bin) );
                    }
                }
            }
        }
        return weight;
    };
    
    weight *= get_weight_from_hist(config::prefire_photon_hist, photon_pt, photon_eta, unc);
    weight *= get_weight_from_hist(config::prefire_jet_hist, jet_pt, jet_eta, unc);
    return weight;      
};

// get DD fake rate 
float dd_fakerate( const rvec<float>& tau_pt,
                            const float & tau_pt_over_jet_pt,
                            const rvec<int>& col_idx,
                            const rvec<bool>& tau_mask,
                            const std::string& run_type) {

  Int_t tauptbin = config::ff_hist->GetXaxis()->FindBin(tau_pt[tau_mask][col_idx[0]]);
  Int_t tauptojetptbin = config::ff_hist->GetYaxis()->FindBin(tau_pt_over_jet_pt);
  float FF  = config::ff_hist->GetBinContent(tauptbin,tauptojetptbin);
  return FF;
    

}

// get ele id scale factor
float ele_id_scale_factor( const rvec<float>& ele_pt,
                            const rvec<float>& ele_eta,
                            const rvec<float>& ele_mask,
                            const std::string& run_type) {

    float scale_factor = 1.0;
    for (uint i = 0; i < ele_pt[ele_mask].size(); i++) {
        scale_factor *= config::ele_SF->getSFID(ele_pt[ele_mask][i],ele_eta[ele_mask][i], run_type);
    }
    return scale_factor;
    

}
// get muon id scale factor
float muon_id_scale_factor( const rvec<float>& muon_pt,
                            const rvec<float>& muon_eta,
                            const rvec<float>& muon_mask,
                            const std::string& run_type) {

    float scale_factor = 1.0;
    for (uint i = 0; i < muon_pt[muon_mask].size(); i++) {
        scale_factor *= config::muon_SF->getSFID(muon_pt[muon_mask][i],muon_eta[muon_mask][i], run_type);
    }
    return scale_factor;
    

}
float muon_iso_scale_factor( const rvec<float>& muon_pt,
                            const rvec<float>& muon_eta,
                            const rvec<float>& muon_mask,
                            const std::string& run_type) {

    float scale_factor = 1.0;
    for (uint i = 0; i < muon_pt[muon_mask].size(); i++) {
        scale_factor *= config::muon_SF->getSFISO(muon_pt[muon_mask][i],muon_eta[muon_mask][i], run_type);
    }
    return scale_factor;
    

}


// get tau fake rate scale factor (for both e an mu)
float tau_fake_scale_factor(	const float& tau_pt,
                                const float& tau_eta,
								const UChar_t& tau_genPartFlav,
								const std::string& part,
                                const std::string& run_type) {
	
	float scale_factor = 1.0;
	
    if (part == "Ele") {                                      // there are outgoing gen_electrons
        if ((tau_genPartFlav & 1 == 1) || (tau_genPartFlav & 3 == 3)) {
            scale_factor *= config::tau_vse_SF->getSFvsEta(tau_eta, tau_genPartFlav, run_type);
        }
    }
    else if (part == "Muon") {                                 // these are outgoing gen_muons
        if ((tau_genPartFlav & 2 == 2) || (tau_genPartFlav & 4 == 4)) {
            scale_factor *= config::tau_vsmu_SF->getSFvsEta(tau_eta, tau_genPartFlav, run_type);
        }
    }
    else if (part == "Jet") {                                  // these are outgoing hadronic taus
        if ((tau_genPartFlav & 5 == 5)) {
            scale_factor *= config::tau_vsjet_SF->getSFvsPT(tau_pt, run_type);
        }
    }
	
	return scale_factor;
};

float GetTopQscale(
                    const float& Mll,
                    const std::string& var){
  if (Mll <0 || var=="nom" || !config::TT) {
    return 1;
  }
  else {
    //double unc = 0.135 - 5.981*pow(10,-5)*Mll + 1.807*pow(10,-7)*pow(Mll,2) - 1.815*pow(10,-10)*pow(Mll,3) + 7.875*pow(10,-14)*pow(Mll,4) - 1.229*pow(10,-17)*pow(Mll,5);
    double unc = .007 - 1.238*pow(10,-5)*Mll + 9.69*pow(10,-9)*pow(Mll,2);

    double weight = 0;
    if (var=="nom") weight = 1;
    else if (var=="up") weight = 1+fabs(unc);
    else if (var=="down") weight = 1-fabs(unc);

    return weight;
  }
}


float GetTopPDF(
                 const float& Mll,
                 const std::string& var){
  if (Mll <0 || var=="nom" || !config::TT) {
    return 1;
  }
  else {
    //double unc = 0.49 - 0.0007795*Mll + 1.59*pow(10,-6)*pow(Mll,2) - 1.166*pow(10,-9)*pow(Mll,3) + 3.93*pow(10,-13)*pow(Mll,4) - 4.72*pow(10,-17)*pow(Mll,5);
    double unc = .07 - 0.0001739*Mll + 1.383*pow(10,-7)*pow(Mll,2);

    double weight = 0;
    if (var=="nom") weight = 1;
    else if (var=="up") weight = 1+fabs(unc);
    else if (var=="down") weight = 1-fabs(unc);

    return weight;
  }
}

// calculate top pt reweighting
float calc_top_pt_reweighting( const rvec<int>& gen_pdg,
                               const rvec<float>& gen_pt) {
    auto func = [] ( const double pt ) { return std::exp( 0.0615 - 0.0005 * pt ); };
    float top1_pt = -10.;
    float top2_pt = -10.;
    bool find_top1 = false;
    bool find_top2 = false;
    for(uint i = 0; i < gen_pt.size(); i++){
        if(gen_pdg[i] == 6){
            top1_pt = gen_pt[i];
        }
        else if(gen_pdg[i] == -6){
            top2_pt = gen_pt[i];
        }
        if (find_top1 && find_top2) break;
    }
    if (find_top1 && find_top2){
        if (top1_pt > 500.) top1_pt = 500.;
        if (top2_pt > 500.) top2_pt = 500.;
        return std::sqrt(func(gen_pt[2]) * func(gen_pt[3]));
    }
    return 1.0;
};
// calculate muon reco eff
float muon_reco_eff( const rvec<float>& muon_pt,
                               const rvec<float>& muon_eta,
                               const rvec<float>& muon_phi,
                               const rvec<float>& muon_mass,
                               const rvec<bool>& muon_mask,
                               const std::string runtype) {
    double weight = 1.;
    TLorentzVector muon;
    for(uint i = 0; i < muon_pt[muon_mask].size(); i++){
        muon.SetPtEtaPhiM(muon_pt[muon_mask][i],muon_eta[muon_mask][i],muon_phi[muon_mask][i],muon_mass[muon_mask][i]);
        if (runtype == ""){
            if( fabs(muon.Eta()) < 1.6){
                if( muon.Pt() <= 100.){
                    weight *= 0.9943;
                }
                else if ( 100. < muon.Pt() <= 150.){
                    weight *= 0.9948;
                }
                else if ( 150. < muon.Pt() <= 200.){
                    weight *= 0.9950;
                }
                else if ( 200. < muon.Pt() <= 300.){
                    weight *= 0.994;
                }
                else if ( 300. < muon.Pt() <= 400.){
                    weight *= 0.9914;
                }
                else if ( 400. < muon.Pt() <= 600.){
                    weight *= 0.993;
                }
                else if ( 600. < muon.Pt() <= 1500.){
                    weight *= 0.991;
                }
                else if ( 1500. < muon.Pt() <= 3500.){
                    weight *= 1.;
                }
            }
            if( fabs(muon.Eta()) >= 1.6){
                if ( 100. < muon.Pt() <= 150.){
                    weight *= 0.993;
                }
                else if ( 150. < muon.Pt() <= 200.){
                    weight *= 0.99;
                }
                else if ( 200. < muon.Pt() <= 300.){
                    weight *= 0.988;
                }
                else if ( 300. < muon.Pt() <= 400.){
                    weight *= 0.981;
                }
                else if ( 400. < muon.Pt() <= 600.){
                    weight *= 0.983;
                }
                else if ( 600. < muon.Pt() <= 1500.){
                    weight *= 0.978;
                }
                else if ( 1500. < muon.Pt() <= 3500.){
                    weight *= 0.98;
                }
            }
        }
        else if (runtype == "Up"){
            if( fabs(muon.Eta()) < 1.6){
                if( muon.Pt() <= 100.){
                    weight *= (0.9943 + 0.0007);
                }
                else if ( 100. < muon.Pt() <= 150.){
                    weight *= (0.9948 + 0.0007);
                }
                else if ( 150. < muon.Pt() <= 200.){
                    weight *= (0.9950 + 0.0009); 
                }
                else if ( 200. < muon.Pt() <= 300.){
                    weight *= (0.994 + 0.001);
                }
                else if ( 300. < muon.Pt() <= 400.){
                    weight *= (0.9914 + 0.0009); 
                }
                else if ( 400. < muon.Pt() <= 600.){
                    weight *= (0.993 + 0.002);
                }
                else if ( 600. < muon.Pt() <= 1500.){
                    weight *= (0.991 + 0.004);
                }
                else if ( 1500. < muon.Pt() <= 3500.){
                    weight *= 1.;
                }
            }
            if( fabs(muon.Eta()) >= 1.6){
                if ( 100. < muon.Pt() <= 150.){
                    weight *=( 0.993 + 0.001);
                }
                else if ( 150. < muon.Pt() <= 200.){
                    weight *= (0.99 + 0.001); 
                }
                else if ( 200. < muon.Pt() <= 300.){
                    weight *= (0.988 + 0.001);
                }
                else if ( 300. < muon.Pt() <= 400.){
                    weight *= (0.981 + 0.002);
                }
                else if ( 400. < muon.Pt() <= 600.){
                    weight *= (0.983 + 0.003);
                }
                else if ( 600. < muon.Pt() <= 1500.){
                    weight *= (0.978 + 0.006);
                }
                else if ( 1500. < muon.Pt() <= 3500.){
                    weight *= (0.98 + 0.03);
                }
            }
        }
        else if (runtype == "Down"){
            if( fabs(muon.Eta()) < 1.6){
                if( muon.Pt() <= 100.){
                    weight *= (0.9943 - 0.0007);
                }
                else if ( 100. < muon.Pt() <= 150.){
                    weight *= (0.9948 - 0.0007);
                }
                else if ( 150. < muon.Pt() <= 200.){
                    weight *= (0.9950 - 0.0009); 
                }
                else if ( 200. < muon.Pt() <= 300.){
                    weight *= (0.994 - 0.001);
                }
                else if ( 300. < muon.Pt() <= 400.){
                    weight *= (0.9914 - 0.0009); 
                }
                else if ( 400. < muon.Pt() <= 600.){
                    weight *= (0.993 - 0.002);
                }
                else if ( 600. < muon.Pt() <= 1500.){
                    weight *= (0.991 - 0.004);
                }
                else if ( 1500. < muon.Pt() <= 3500.){
                    weight *= (1. - 0.1);
                }
            }
            if( fabs(muon.Eta()) >= 1.6){
                if ( 100. < muon.Pt() <= 150.){
                    weight *=( 0.993 - 0.001);
                }
                else if ( 150. < muon.Pt() <= 200.){
                    weight *= (0.99 - 0.001); 
                }
                else if ( 200. < muon.Pt() <= 300.){
                    weight *= (0.988 - 0.001);
                }
                else if ( 300. < muon.Pt() <= 400.){
                    weight *= (0.981 - 0.002);
                }
                else if ( 400. < muon.Pt() <= 600.){
                    weight *= (0.983 - 0.003);
                }
                else if ( 600. < muon.Pt() <= 1500.){
                    weight *= (0.978 - 0.006);
                }
                else if ( 1500. < muon.Pt() <= 3500.){
                    weight *= (0.98 - 0.03);
                }
            }
        }
    }
    return 1.0;
};

// function for pdf weights
float get_pdf_weight(const unsigned int n_pdf_weight,
                     const rvec<float>& pdf_weight) {
    return pdf_weight[n_pdf_weight];
};

// apply kfactor for W background
float get_kfactor(	const rvec<int>& gen_pdg,
					const rvec<float>& gen_mass) {
                        
    if (config::W_kfactor_hist != NULL)  {
        auto mass = 0.0;
        for (int i = 0; i < gen_pdg.size(); i++) {
            if (abs(gen_pdg[i]) == 24) {
                mass = gen_mass[i];
                break;
            }
        }    
        return config::W_kfactor_hist->GetBinContent(config::W_kfactor_hist->GetXaxis()->FindBin(mass));
    }
    else
        return 1.0;
}
