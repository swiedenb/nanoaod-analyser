#include "include/PhysicalQuantities.hh"
#include "include/Particles.hh"
#include "include/weights.hh"
#include "include/event_cleaner.hh"
#include "include/PDFTool.hh"
#include "include/TauIDSFTool.h"
#include "include/MuonSFTool.h"
#include "include/DDTool.h"
#include "include/Systematics.hh"
#include "include/MuonScale/GEScaleSyst.h"
#include "include/MuonScale/GEScaleSyst.cc"
#include <iostream>
#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RVec.hxx>
#include <math.h>
#include <TRandom3.h>
#include <TSystemDirectory.h>
#include "TPRegexp.h"
#include <variant>


// reduce amount to write for each vector by a lot
template < typename T >
using rvec = ROOT::VecOps::RVec<T>;
using RNode = ROOT::RDF::RNode;

json goldenjson;
json cfg;

std::multimap< std::string, ROOT::RDF::RResultPtr<TH2D> > hist_dict_2d;

template<typename ... Base>
struct Visitor : Base...
{
    using Base::operator()...;
};
template<typename ... T> Visitor(T...) -> Visitor<T...>;

std::multimap< std::string, std::variant< ROOT::RDF::RResultPtr<TH1D>, ROOT::RDF::RResultPtr<TH2D> > > hist_dict;


TString runName;
TRandom3* gauss = new TRandom3();

std::string met_branch_name = "MET";




void define_weight(RNode& df,  
                    const std::vector< std::string >& def_weights,
                    std::vector< std::string >& list_of_weights,
                    const std::string& weight_shift = "",
                    const std::string& weight_replace = "") {
    std::string weight_string = "";
    for (const auto& it: def_weights) {
       // df = df.Filter([it](const float& weight){
       //                     if (weight > 10.){
       //                         std::cout<<it<< "    ";
       //                         std::cout<< weight << std::endl;
       //                     }
       //                     return true;
       //                     },{it});
        if (it == weight_replace) {
            weight_string += weight_shift + "*";
        } else {
            weight_string += it + "*";
        }
        //std::cout<<"Weight: "<< it << " " << *df.Sum(it) <<std::endl;

    }
    //std::cout<<std::endl;
    weight_string.pop_back();
    
    auto defColNames = df.GetDefinedColumnNames();
    if (weight_shift == "") {
        if (std::find( defColNames.begin(), defColNames.end(), "total_weight" ) == defColNames.end()) {
            df = df.Define("total_weight", weight_string);
            list_of_weights.push_back("total_weight");
        }
    } else {
        if (std::find( defColNames.begin(), defColNames.end(), "total_weight_" + weight_shift ) == defColNames.end()) {
            df = df.Define("total_weight_" + weight_shift, weight_string);
            list_of_weights.push_back("total_weight_" + weight_shift);
        }
    }
}


void calc_tau_uncertainty ( RNode& df,
                            const std::string& shift_part,
                            const std::string& shift_dir,
                            const std::string& bin_type,
                            const std::string& bin) {    
    df = df.Define( "Tau" + shift_part + "FakeScaleFactor" + shift_dir + "_" + bin_type + '_' + bin, 
                    [shift_part, shift_dir, bin, bin_type]( 
                        const rvec<float>& tau_pt,
                        const rvec<float>& tau_eta,
                        const rvec<int>& tau_dm,
                        const rvec<bool>& tau_mask,
                        const rvec<UChar_t>& tau_genPartFlav
                        )
                        {
                            float weight = 1.0;
                            if (bin_type == "" && bin == "") {
                                for (uint i = 0; i < tau_pt[tau_mask].size(); i++) {
                                    weight *= tau_fake_scale_factor(tau_pt[tau_mask][i], tau_eta[tau_mask][i], tau_genPartFlav[tau_mask][i], shift_part, shift_dir);
                                }
                            }
                            if (bin_type == "dm") {
                                for (uint i = 0; i < tau_dm[tau_mask].size(); i++) {
                                    int cdm = tau_dm[tau_mask][i];
                                    if (cdm == std::stoi(bin)) {
                                        weight *= tau_fake_scale_factor(tau_pt[tau_mask][i], tau_eta[tau_mask][i], tau_genPartFlav[tau_mask][i], shift_part, shift_dir);
                                    }
                                }
                            } else if (bin_type == "eta") {
                                for (uint i = 0; i < tau_eta[tau_mask].size(); i++) {
                                    float ceta = std::abs( tau_eta[tau_mask][i] );
                                    if (bin == "eta0p4") {
                                        if (ceta < 0.4) {
                                            weight *= tau_fake_scale_factor(tau_pt[tau_mask][i], tau_eta[tau_mask][i], tau_genPartFlav[tau_mask][i], shift_part, shift_dir);
                                        }
                                    } else if (bin == "eta0p4to0p8") {
                                        if (ceta >= 0.4 && ceta < 0.8) {
                                            weight *= tau_fake_scale_factor(tau_pt[tau_mask][i], tau_eta[tau_mask][i], tau_genPartFlav[tau_mask][i], shift_part, shift_dir);
                                        }
                                    } else if (bin == "eta0p8to1p2") {
                                        if (ceta >= 0.8 && ceta < 1.2) {
                                            weight *= tau_fake_scale_factor(tau_pt[tau_mask][i], tau_eta[tau_mask][i], tau_genPartFlav[tau_mask][i], shift_part, shift_dir);
                                        }
                                        
                                    } else if (bin == "eta1p2to1p7") {
                                        if (ceta >= 1.2 && ceta < 1.7) {
                                            weight *= tau_fake_scale_factor(tau_pt[tau_mask][i], tau_eta[tau_mask][i], tau_genPartFlav[tau_mask][i], shift_part, shift_dir);
                                        }
                                        
                                    } else if (bin == "eta1p7") {
                                        if (ceta >= 1.7) {
                                            weight *= tau_fake_scale_factor(tau_pt[tau_mask][i], tau_eta[tau_mask][i], tau_genPartFlav[tau_mask][i], shift_part, shift_dir);
                                        }
                                    } else if (bin == "barrel") {
                                        if (ceta < 1.4442) {
                                            weight *= tau_fake_scale_factor(tau_pt[tau_mask][i], tau_eta[tau_mask][i], tau_genPartFlav[tau_mask][i], shift_part, shift_dir);
                                        }
                                    } else if (bin == "endcap") {
                                        if (ceta < 1.566) {
                                            weight *= tau_fake_scale_factor(tau_pt[tau_mask][i], tau_eta[tau_mask][i], tau_genPartFlav[tau_mask][i], shift_part, shift_dir);
                                        }
                                    }
                                }
                            }
                            return weight;
                        },
                        {"Tau_pt", "Tau_eta", "Tau_decayMode", "Tau_mask", "Tau_genPartFlav"});
}
rvec<float> muon_reso_smearing(
                        const rvec<float>& pt,
                        const rvec<float>& eta,
                        const rvec<float>& phi,
                        const rvec<float>& mass,
                        const std::string& direction){
          TLorentzVector muon;
          rvec<float> pt_smeared;
          if(!config::runOnData){
              for( unsigned int i = 0; i < pt.size(); i++){
                  muon.SetPtEtaPhiM(pt[i],eta[i],phi[i],mass[i]);
                  auto sigma = 0.;
                  auto p = muon.P();
                  if(0. <= abs(eta[i]) < 1.2){
                    if(config::era == 2016){
                        sigma = 0.0062 + 0.0001*p - 1.0*pow(10,-7)*p*p + 5.7*pow(10,-11)*p*p*p - 1.1*pow(10,-14)*p*p*p*p;
                    }
                    else if(config::era == 2017){
                        sigma = 0.0053 + 0.00011*p - 1.3*pow(10,-7)*p*p + 6.9*pow(10,-11)*p*p*p - 1.3*pow(10,-14)*p*p*p*p;
                    }
                    else if(config::era == 2018){
                        sigma = 0.0062 + 0.000096*p - 9.7*pow(10,-8)*p*p + 4.9*pow(10,-11)*p*p*p - 9.*pow(10,-15)*p*p*p*p;
                    }
                  }
                  else if(1.2 <= abs(eta[i]) < 2.1){
                    if(config::era == 2016){
                        sigma = 0.0134 + 0.000063*p - 4.7*pow(10,-8)*p*p + 2.6*pow(10,-11)*p*p*p - 5*pow(10,-15)*p*p*p*p;
                    }
                    else if(config::era == 2017){
                        sigma = 0.0136 + 0.000063*p - 2.6*pow(10,-8)*p*p + 1.3*pow(10,-12)*p*p*p + 3.*pow(10,-15)*p*p*p*p;
                    }
                    else if(config::era == 2018){
                        sigma = 0.0136 + 0.000052*p - 2.4*pow(10,-8)*p*p + 5.0*pow(10,-12)*p*p*p;
                    }
                  }
                  else if(2.1 <= abs(eta[i]) < 2.4){
                    if(config::era == 2016){
                        sigma = 0.0151 + 0.000114*p - 3.7*pow(10,-8)*p*p - 3.9*pow(10,-12)*p*p*p + 1.*pow(10,-15)*p*p*p*p;
                    }
                    else if(config::era == 2017){
                        sigma = 0.0170 + 0.000084*p - 2.6*pow(10,-9)*p*p - 2.3*pow(10,-11)*p*p*p + 8.*pow(10,-15)*p*p*p*p;
                    }
                    else if(config::era == 2018){
                        sigma = 0.0174 + 0.000087*p - 3.3*pow(10,-9)*p*p - 1.6*pow(10,-11)*p*p*p + 5.*pow(10,-15)*p*p*p*p;
                    }
                  }
                  auto p_new = p;
                  if(abs(eta[i]) >= 1.2){
                      p_new = p * ( 1 + gauss->Gaus(0,sigma * 0.57));
                  }
                  if(direction == "Up"){
                    p_new = p * ( 1 + gauss->Gaus(0,sigma * 0.46));
                  }
                  else if(direction == "Down"){
                    p_new = p * ( 1 - gauss->Gaus(0,sigma * 0.46));
                  }
                  auto pt_new = p_new * 1./cosh(eta[i]);
                  pt_smeared.push_back(std::max(0.,pt_new));
              
              }
          }
          else {
            pt_smeared = pt;
          }
          return pt_smeared;
                        }
bool jet_matching(
                        const float& tau_pt,
                        const float& tau_eta,
                        const float& tau_phi,
                        const float& tau_mass,
                        const rvec<float>& jet_pt,
                        const rvec<float>& jet_eta,
                        const rvec<float>& jet_phi,
                        const rvec<float>& jet_mass){
          TLorentzVector tau;
          TLorentzVector jet;
          bool matched_jet = false;
          tau.SetPtEtaPhiM(tau_pt,tau_eta,tau_phi, tau_mass);
          for( unsigned int i = 0; i < jet_pt.size(); i++){
            jet.SetPtEtaPhiM(jet_pt[i],jet_eta[i],jet_phi[i],jet_mass[i]);
            if(tau.DeltaR(jet) < 0.2){
                if(abs(jet_eta[i]) >= 2.5){
                    continue;
                }
                else{
                    matched_jet = true;
                    break;
                }
            }
            else{
                continue;
            }
          }
          return matched_jet;
 }


rvec<float> scale_correct_muon(
                        const rvec<float>& pt,
                        const rvec<float>& eta,
                        const rvec<float>& phi,
                        const rvec<int>& q,
                        const std::string& direction){
          GEScaleSyst *GE = new GEScaleSyst();
          GE->SetVerbose(0);
          rvec<float> pt_corr;
          int matrix_number;
          if(config::era == 2016){
            matrix_number = 160000;
          }
          else if(config::era == 2017){
            matrix_number = 170000;
          }
          else if(config::era == 2018){
            matrix_number = 180000;
          }
          
          for( unsigned int i = 0; i < pt.size(); i++){
            float pt_sum;
            if(abs(eta[i])>2.4 or abs(phi[i])>=3.1416){
                pt_corr.push_back(pt[i]);
                continue;
            }
            for( unsigned int j = 0; j < 50; j++){
                pt_sum+=GE->GEScaleCorrPt(matrix_number + j,pt[i],eta[i],phi[i],q[i],false);
            }
            pt_sum = pt_sum/50;
            float pt_ratio = 1 - pt_sum/pt[i];
            if(direction == "_muonScaleUp"){
                pt_corr.push_back(pt[i] * (1 + pt_ratio));
            }
            else if(direction == "_muonScaleDown"){
                pt_corr.push_back(pt[i]*(1- pt_ratio));
            }
            else{
                pt_corr.push_back(pt[i]);
            }
          }

          delete GE;
          return pt_corr;
}
float correct_met_single(
                   const float& muon_pt,
                   const float& muon_tuneP_pt,
                   const float& muon_eta,
                   const float& muon_phi,
                   const float& muon_mass,
                   const bool& PFCand,
                   const float& tau_pt,
                   const float& tau_eta,
                   const float& tau_phi,
                   const float& tau_mass,
                   const float& tau_pt_norm,
                   const float& tau_eta_norm,
                   const float& tau_phi_norm,
                   const float& tau_mass_norm,
                   const float& met,
                   const float& met_phi,
                   const UInt_t& runnb,
                   const int& npv) {

        float met_corr = met;
        float met_phi_corr = met_phi;
        if(config::doXY){
             met_corr = METXYCorr_Met_MetPhi(met,met_phi,runnb,config::era, config::runOnData, npv).first;
             met_phi_corr = METXYCorr_Met_MetPhi(met,met_phi,runnb,config::era, config::runOnData, npv).second;
        }
        TLorentzVector muon_pf_p4, muon_tune_p4,met_p4,tau_p4,tau_es_p4;
        met_p4.SetPtEtaPhiE(met_corr,0,met_phi_corr,met_corr);
        muon_tune_p4.SetPtEtaPhiM(muon_tuneP_pt,0,muon_phi,muon_mass);
        muon_pf_p4.SetPtEtaPhiM(muon_pt,0,muon_phi,muon_mass);
        tau_es_p4.SetPtEtaPhiM(tau_pt,tau_eta,tau_phi,tau_mass);
        tau_p4.SetPtEtaPhiM(tau_pt_norm,tau_eta_norm,tau_phi_norm,tau_mass_norm);
        if(PFCand){
            float met_px = met_p4.Px() + muon_pf_p4.Px() - muon_tune_p4.Px();
            float met_py = met_p4.Py() + muon_pf_p4.Py() - muon_tune_p4.Py();
            float met_pt = sqrt(pow(met_px,2) + pow(met_py,2)); 
            met_p4.SetPxPyPzE(met_px,met_py,0,met_pt);
        }

        //float met_px = met_p4.Px() + tau_p4.Px() - tau_es_p4.Px();
        //float met_py = met_p4.Py() + tau_p4.Py() - tau_es_p4.Py();
        //float met_pt = sqrt(pow(met_px,2) + pow(met_py,2)); 
        //met_p4.SetPxPyPzE(met_px,met_py,0,met_pt);

        return met_p4.Pt();
        

};
float correct_met_single_phi(
                   const float& muon_pt,
                   const float& muon_tuneP_pt,
                   const float& muon_eta,
                   const float& muon_phi,
                   const float& muon_mass,
                   const bool& PFCand,
                   const float& tau_pt,
                   const float& tau_eta,
                   const float& tau_phi,
                   const float& tau_mass,
                   const float& tau_pt_norm,
                   const float& tau_eta_norm,
                   const float& tau_phi_norm,
                   const float& tau_mass_norm,
                   const float& met,
                   const float& met_phi,
                   const UInt_t& runnb,
                   const int& npv) {

        float met_corr = met;
        float met_phi_corr = met_phi;
        if(config::doXY){
             met_corr = METXYCorr_Met_MetPhi(met,met_phi,runnb,config::era, config::runOnData, npv).first;
             met_phi_corr = METXYCorr_Met_MetPhi(met,met_phi,runnb,config::era, config::runOnData, npv).second;
        }
        TLorentzVector muon_pf_p4, muon_tune_p4,met_p4,tau_p4,tau_es_p4;
        met_p4.SetPtEtaPhiE(met_corr,0,met_phi_corr,met_corr);
        muon_tune_p4.SetPtEtaPhiM(muon_tuneP_pt,0,muon_phi,muon_mass);
        muon_pf_p4.SetPtEtaPhiM(muon_pt,0,muon_phi,muon_mass);
        tau_es_p4.SetPtEtaPhiM(tau_pt,tau_eta,tau_phi,tau_mass);
        tau_p4.SetPtEtaPhiM(tau_pt_norm,tau_eta_norm,tau_phi_norm,tau_mass_norm);
        if(PFCand){
            float met_px = met_p4.Px() + muon_pf_p4.Px() - muon_tune_p4.Px();
            float met_py = met_p4.Py() + muon_pf_p4.Py() - muon_tune_p4.Py();
            float met_pt = sqrt(pow(met_px,2) + pow(met_py,2)); 
            met_p4.SetPxPyPzE(met_px,met_py,0,met_pt);
        }
        //float met_px = met_p4.Px() + tau_p4.Px() - tau_es_p4.Px();
        //float met_py = met_p4.Py() + tau_p4.Py() - tau_es_p4.Py();
        //float met_pt = sqrt(pow(met_px,2) + pow(met_py,2)); 
        //met_p4.SetPxPyPzE(met_px,met_py,0,met_pt);

        return met_p4.Phi();
        

};
float correct_met(
                   const rvec<float>& muon_pt,
                   const rvec<float>& muon_tuneP_pt,
                   const rvec<float>& muon_eta,
                   const rvec<float>& muon_phi,
                   const rvec<float>& muon_mass,
                   const rvec<UChar_t>& id,
                   const rvec<bool>& PFCand,
                   const float& met,
                   const float& met_phi,
                   const UInt_t& runnb,
                   const int& npv) {

        float met_corr = met;
        float met_phi_corr = met_phi;
        if(config::doXY){
             met_corr = METXYCorr_Met_MetPhi(met,met_phi,runnb,config::era, config::runOnData, npv).first;
             met_phi_corr = METXYCorr_Met_MetPhi(met,met_phi,runnb,config::era, config::runOnData, npv).second;
        }
        TLorentzVector muon_pf_p4, muon_tune_p4,met_p4;
        met_p4.SetPtEtaPhiE(met_corr,0,met_phi_corr,met_corr);
        for( unsigned int i = 0; i < muon_pt.size(); i++){
            muon_tune_p4.SetPtEtaPhiM(muon_tuneP_pt[i],0,muon_phi[i],muon_mass[i]);
            muon_pf_p4.SetPtEtaPhiM(muon_pt[i],0,muon_phi[i],muon_mass[i]);

            if (id[i] == 2){
                if(PFCand[i]){
                    float met_px = met_p4.Px() + muon_pf_p4.Px() - muon_tune_p4.Px();
                    float met_py = met_p4.Py() + muon_pf_p4.Py() - muon_tune_p4.Py();
                    float met_pt = sqrt(pow(met_px,2) + pow(met_py,2)); 
                    met_p4.SetPxPyPzE(met_px,met_py,0,met_pt);

                }
                else{
                   continue;
                   // float met_px = met_p4.Px() - muon_tune_p4.Px();
                   // float met_py = met_p4.Py() - muon_tune_p4.Py();
                   // float met_pt = sqrt(pow(met_px,2) + pow(met_py,2)); 
                   // met_p4.SetPxPyPzE(met_px,met_py,0,met_pt);
                }
            }
        }
        return met_p4.Pt();
        

};
float correct_met_phi(
                   const rvec<float>& muon_pt,
                   const rvec<float>& muon_tuneP_pt,
                   const rvec<float>& muon_eta,
                   const rvec<float>& muon_phi,
                   const rvec<float>& muon_mass,
                   const rvec<UChar_t>& id,
                   const rvec<bool>& PFCand,
                   const float& met,
                   const float& met_phi,
                   const UInt_t& runnb,
                   const int& npv) {

        float met_corr = met;
        float met_phi_corr = met_phi;
        if(config::doXY){
             met_corr = METXYCorr_Met_MetPhi(met,met_phi,runnb,config::era, config::runOnData, npv).first;
             met_phi_corr = METXYCorr_Met_MetPhi(met,met_phi,runnb,config::era, config::runOnData, npv).second;
        }
        TLorentzVector muon_pf_p4, muon_tune_p4,met_p4;
        met_p4.SetPtEtaPhiE(met_corr,0,met_phi_corr,met_corr);
        for( unsigned int i = 0; i < muon_pt.size(); i++){
            muon_tune_p4.SetPtEtaPhiM(muon_tuneP_pt[i],0,muon_phi[i],muon_mass[i]);
            muon_pf_p4.SetPtEtaPhiM(muon_pt[i],0,muon_phi[i],muon_mass[i]);

            if (id[i] == 2){
                if(PFCand[i]){
                    float met_px = met_p4.Px() + muon_pf_p4.Px() - muon_tune_p4.Px();
                    float met_py = met_p4.Py() + muon_pf_p4.Py() - muon_tune_p4.Py();
                    float met_pt = sqrt(pow(met_px,2) + pow(met_py,2)); 
                    met_p4.SetPxPyPzE(met_px,met_py,0,met_pt);

                }
                else{
                    //float met_px = met_p4.Px() - muon_tune_p4.Px();
                    //float met_py = met_p4.Py() - muon_tune_p4.Py();
                    //float met_pt = sqrt(pow(met_px,2) + pow(met_py,2)); 
                    //met_p4.SetPxPyPzE(met_px,met_py,0,met_pt);
                    continue;
                }
            }
        }
        return met_p4.Phi();
        

};

float closest_pf_muo_pt(
                   const rvec<float>& muon_pt,
                   const rvec<float>& muon_eta,
                   const rvec<float>& muon_phi,
                   const rvec<bool>& muon_PFcand,
                   const float& sel_muon_pt,
                   const float& sel_muon_eta,
                   const float& sel_muon_phi
                   ) {
            int best_i = 0;
            double min_delR = 0.2;
            for (unsigned int i = 0; i < muon_pt.size(); i++){
                if (muon_PFcand[i]){
                    double delR = sqrt(pow((sel_muon_eta - muon_eta[i]),2) + pow(sel_muon_phi - muon_phi[i],2));
                    if (delR < min_delR){
                        min_delR = delR;
                        best_i = i;
                    }
                        
                }
            }
	return muon_pt[best_i];
};
// Di muon pair number cut: how many "trues" in mask
bool deltar_tau_muon_cut(
                   const rvec<float>& tau_eta,
                   const rvec<float>& tau_phi,
                   const rvec<float>& muon_eta,
                   const rvec<float>& muon_phi,
                   const rvec<bool>& tau_mask,
                   const rvec<bool>& muon_mask,
                   const rvec<int>& col_idx) {
               double delR = sqrt(pow((tau_eta[tau_mask][col_idx[0]] - muon_eta[muon_mask][col_idx[1]]),2) + pow(tau_phi[tau_mask][col_idx[0]]-muon_phi[muon_mask][col_idx[1]],2));
               if (delR < 0.5){
                    return false;
               }
	return true;
};
// Di muon pair number cut: how many "trues" in mask
bool dimuonpair_cut(const rvec<float>& pt,
                   const rvec<float>& eta,
                   const rvec<float>& phi,
                   const rvec<bool>& mask) {
    if ( pt[mask].size()>= 2){
        for (uint i = 0; i<pt[mask].size(); i++){
            for( uint j = i + 1; j < pt[mask].size(); j++){
               double delphi = delta_phi(phi[mask][i],phi[mask][j]);
               double deleta = delta_eta(eta[mask][i],eta[mask][j]);
               double delR = sqrt(pow(deleta,2) + pow(delphi,2));
               if (delR > 0.2){
                    return false;
               }
            }
        }
    }
	return true;
};
// mu tau candidate delR cut:
bool mutaudelR_cut(const float& eta_mu,
                   const float& phi_mu,
                   const float& eta_tau,
                   const float& phi_tau) {
    double delEta = delta_eta(eta_mu,eta_tau);
    double delPhi = delta_phi(phi_mu,phi_tau);
    double delR = sqrt(pow(delEta,2) + pow(delPhi,2));
    if (delR < 0.5){
        return false;
    }
    return true;
};
float wwuncertainty(const float& coll_mass){
    float weight = 1.;
    if(config::wwuncertainty){
        weight = (1 - (0.993-0.00002001 * coll_mass + 0.00000000283 * coll_mass * coll_mass));
    }
    return weight;
}
// Particle number cut: how many "trues" in mask
bool nparticle_cut(const rvec<bool>& mask) {
	return std::count(mask.begin(), mask.end(), true) > 0.;
};
// Particle number veto: how many "trues" in mask
bool nparticle_veto(const rvec<bool>& mask) {
	return std::count(mask.begin(), mask.end(), true) == 0;
};

// Returns first mask_true particle for quantity
float selected_jet_quant(const rvec<int>& quantity, 
								const rvec<bool>& mask,
                                const rvec<int>& col_idx) {
	return quantity[mask][col_idx[0]];
};
// Returns first mask_true particle for quantity
float selected_part_col_idx_tau_flav(const rvec<UChar_t>& quantity, 
								const rvec<bool>& mask,
                            const rvec<int>& col_idx) {
	return quantity[mask][col_idx[0]];
};
float selected_part_col_idx_tau(const rvec<float>& quantity, 
								const rvec<bool>& mask,
                            const rvec<int>& col_idx) {
	return quantity[mask][col_idx[0]];
};
int selected_part_col_idx_tau_idx(const rvec<int>& quantity, 
								const rvec<bool>& mask,
                            const rvec<int>& col_idx) {
	return quantity[mask][col_idx[0]];
};
float selected_part_col_idx_tau_barrel(const rvec<float>& quantity, 
								const rvec<bool>& mask,
                            const rvec<int>& col_idx,
                            const rvec<float>& eta) {
    if(abs(eta[mask][col_idx[0]]) < 1.46){
	    return quantity[mask][col_idx[0]];
    }
    return 0.;
};
float selected_part_col_idx_tau_endcap(const rvec<float>& quantity, 
								const rvec<bool>& mask,
                            const rvec<int>& col_idx,
                            const rvec<float>& eta) {
    if(abs(eta[mask][col_idx[0]]) > 1.56){
	    return quantity[mask][col_idx[0]];
    }
    return 0.;
};
bool selected_part_col_idx_muon_bool(const rvec<bool>& quantity, 
								const rvec<bool>& mask,
                            const rvec<int>& col_idx) {
	return quantity[mask][col_idx[1]];
};
float selected_part_col_idx_muon(const rvec<float>& quantity, 
								const rvec<bool>& mask,
                            const rvec<int>& col_idx) {
	return quantity[mask][col_idx[1]];
};
int selected_part_col_idx_muon_idx(const rvec<int>& quantity, 
								const rvec<bool>& mask,
                            const rvec<int>& col_idx) {
	return quantity[mask][col_idx[1]];
};
// Returns first mask_true particle for quantity
float selected_part_quant(const rvec<float>& quantity, 
								const rvec<bool>& mask) {
	return quantity[mask][0];
};

// Returns first mask_true particle for quantity
rvec <float> select_part_quants(const rvec<float>& quantity, 
								const rvec<bool>& mask) {
	return quantity[mask];
};

//Find two tops

rvec <int> lepton_pair_idx(const rvec<int>& pdgID, rvec<int>& mother_idx){
        rvec<int> idx(2);

        int idx_1;
        int idx_2;

        bool first{false};
        bool second{false};

        for (const auto gen_pdg : config::gen_pdgID){
           auto iter = pdgID.begin();
           while (true)
           {
              iter = std::find_if(iter, pdgID.end(), [rhs = gen_pdg](const auto lhs){return abs(lhs) == rhs;});
              const bool found = iter != pdgID.end();
              if (not found)
              {
                 break;
              }
              else
              {
                 if (not first)
                 {
                    idx_1 = std::distance(pdgID.begin(), iter++);
                    if (abs(pdgID[mother_idx[idx_1]]) == 24){
                        first = true;
                        idx[0] = idx_1;
                    }
                 }
                 else if (not second)
                 {
                    idx_2 = std::distance(pdgID.begin(), iter++);
                    if (abs(pdgID[mother_idx[idx_2]]) == 24){
                        second = true;
                        idx[1] = idx_2;
                        break;
                    }
                 }
              }
           }
           if (second)
           {
              break;
           }
        }
        return idx;
}


// fill trigger plots
RNode trigger(	RNode df , bool default_run) {
	auto triggered = df.Filter(config::trigger, "Trigger requirement");
	if (!config::runOnData && default_run) {
		auto trigger_eff_tau_all = df.Histo1D( 	{"Tau_pt_all", "GenVisTau_pt", 600, 0, 6000}, "GenVisTau_pt");
		auto trigger_eff_tau_passed = triggered.Histo1D( 	{"Tau_pt_passed", "Tau_pt", 600, 0, 6000}, "GenVisTau_pt");
		hist_dict.emplace("Trigger", trigger_eff_tau_all);
		hist_dict.emplace("Trigger", trigger_eff_tau_passed);
	};
	return triggered;
}

void create_preselection_hists(	RNode df) {
    if( !config::runOnData){
	    auto gentau_pt = df.Histo1D(	{"GenVisTau_pt", "", 			6000u, 0, 6000}, 					"GenVisTau_pt");
	    auto gen_pdgid = df.Histo1D(	{"GenPart_pdgId", "", 			200u, -100, 100}, 					"GenPart_pdgId");
	    auto gen_mass = df.Histo1D(	{"GenPart_mass", "", 			1000u, 0, 1000}, 					"GenPart_mass");
	    auto mll = df.Histo1D(	{"mll", "", 			1000u, 0, 1000}, 					"Mll");
	    auto ht = df.Histo1D(	{"LHE HT", "", 			1000u, 0, 1000}, 					"LHE_HT");
	    //auto top_pair_invmass = df.Histo1D(	{"top_pair_invmass", "", 			6000u, 0, 6000}, 					"top_pair_invmass");
	    hist_dict.emplace("Preselection", gentau_pt);
	    hist_dict.emplace("Preselection", gen_pdgid);
	    hist_dict.emplace("Preselection", gen_mass);
        hist_dict.emplace("Preselection", mll);
        hist_dict.emplace("Preselection", ht);
	    //hist_dict.emplace("Preselection", top_pair_invmass);
    }
	auto tau_pt = df.Histo1D(	{"Tau_pt", "", 			6000u, 0, 6000}, 								"Tau_pt_ES",                "presel_weight");
	auto tau_eta = df.Histo1D(	{"Tau_eta", "", 		100u, -5, 5}, 									"Tau_eta_ES",               "presel_weight");
	auto tau_phi = df.Histo1D(	{"Tau_phi", "", 		100u, -3.2, 3.2}, 								"Tau_phi_ES",               "presel_weight");
	auto met_pt = df.Histo1D(	{"MET_pt", "MET_pt", 	6000u, 0, 6000}, 							    met_branch_name + "_pt",    "presel_weight");
	auto met_phi = df.Histo1D(	{"MET_phi", "MET_phi", 	100u, -3.2, 3.2}, 					            met_branch_name + "_phi",   "presel_weight");
	auto multiplicity = df.Histo1D(	{"multiplicity", "multiplicity", 	50u, 0., 50.}, 					           "multiplicity"   ,   "presel_weight");
	hist_dict.emplace("Preselection", tau_pt);
	hist_dict.emplace("Preselection", tau_eta);
	hist_dict.emplace("Preselection", tau_phi);
	hist_dict.emplace("Preselection", met_pt);
	hist_dict.emplace("Preselection", met_phi);
	hist_dict.emplace("Preselection", multiplicity);
};
// fill preselection histograms
void fill_preselection(	RNode df) {
    if( !config::runOnData){
	    auto gentau_pt = df.Histo1D(	{"GenVisTau_pt", "", 			6000u, 0, 6000}, 					"GenVisTau_pt");
	    auto gen_pdgid = df.Histo1D(	{"GenPart_pdgId", "", 			200u, -100, 100}, 					"GenPart_pdgId");
	    auto gen_mass = df.Histo1D(	{"GenPart_mass", "", 			1000u, 0, 1000}, 					"GenPart_mass");
	    //auto top_pair_invmass = df.Histo1D(	{"top_pair_invmass", "", 			6000u, 0, 6000}, 					"top_pair_invmass");
	    hist_dict.emplace("Preselection", gentau_pt);
	    hist_dict.emplace("Preselection", gen_pdgid);
	    hist_dict.emplace("Preselection", gen_mass);
	    //hist_dict.emplace("Preselection", top_pair_invmass);
    }
	auto tau_pt = df.Histo1D(	{"Tau_pt", "", 			6000u, 0, 6000}, 					"Tau_pt_ES");
	auto tau_eta = df.Histo1D(	{"Tau_eta", "", 		100u, -5, 5}, 						"Tau_eta_ES");
	auto tau_phi = df.Histo1D(	{"Tau_phi", "", 		100u, -3.2, 3.2}, 					"Tau_phi_ES");	
	auto muon_tP_pt = df.Histo1D(	{"Muon_tP_pt", "", 		6000u, 0, 6000}, 					"Muon_tP_pt");
	auto muon_pt = df.Histo1D(	{"Muon_pt", "", 		6000u, 0, 6000}, 					"Muon_pt");
	auto muon_eta = df.Histo1D(	{"Muon_eta", "", 		100u, -5, 5}, 						"Muon_eta");
	auto muon_phi = df.Histo1D(	{"Muon_phi", "", 		100u, -3.2, 3.2}, 					"Muon_phi");	
	auto met_pt = df.Histo1D(	{"MET_pt", "MET_pt", 	6000u, 0, 6000}, 	    met_branch_name + "_pt" );
	auto met_phi = df.Histo1D(	{"MET_phi", "MET_phi", 	100u, -3.2, 3.2}, 		met_branch_name	+ "_phi" );
	hist_dict.emplace("Preselection", muon_pt);
	hist_dict.emplace("Preselection", muon_tP_pt);
	hist_dict.emplace("Preselection", muon_eta);
	hist_dict.emplace("Preselection", muon_phi);
	hist_dict.emplace("Preselection", tau_pt);
	hist_dict.emplace("Preselection", tau_eta);
	hist_dict.emplace("Preselection", tau_phi);
	hist_dict.emplace("Preselection", met_pt);
	hist_dict.emplace("Preselection", met_phi);
};


// fill datadriven histograms
void fill_datadriven(	RNode df,
                        std::string name = "") {
	auto tau_pt = df.Histo1D(	{"Tau_pt", "", 			6000u, 0, 6000}, 					"Tau_pt_new");
	auto tau_eta = df.Histo1D(	{"Tau_eta", "", 		100u, -5, 5}, 						"Tau_eta_new");
	auto tau_phi = df.Histo1D(	{"Tau_phi", "", 		100u, -3.2, 3.2}, 					"Tau_phi_new");	
	auto met_pt = df.Histo1D(	{"MET_pt", "MET_pt", 	6000u, 0, 6000}, 			met_branch_name + "_pt");
	auto met_phi = df.Histo1D(	{"MET_phi", "MET_phi", 	100u, -3.2, 3.2}, 			met_branch_name + "_phi");
	hist_dict.emplace("Datadriven_" + name, tau_pt);
	hist_dict.emplace("Datadriven_" + name, tau_eta);
	hist_dict.emplace("Datadriven_" + name, tau_phi);
	hist_dict.emplace("Datadriven_" + name, met_pt);
	hist_dict.emplace("Datadriven_" + name, met_phi);
};

// fill hist function - fills hists for each stage
void create_datadriven_hists(	RNode df,
                                std::string name,
                                std::string muon_shift,
                                std::string tau_shift, 
                                std::string met_shift,
                                std::string weight_column) {
    std::string stringcopy = weight_column + muon_shift + tau_shift + met_shift;
    if (tau_shift == "_ES" && (met_shift == "_jer" || met_shift == "_nom") && muon_shift == "") {
        stringcopy = weight_column;
    }
    else if(tau_shift == "_ES" && (met_shift == "_jer" || met_shift == "_nom") && muon_shift != ""){
        stringcopy = weight_column + muon_shift;
    }
    stringcopy.erase(0,12);
    //if( !config::runOnData){
	//    auto gentau_pt = df.Histo1D(	{((TString) "GenVisTau_pt" + stringcopy), "", 			6000u, 0, 6000}, 					"GenVisTau_pt",     weight_column);
	//    hist_dict.emplace(name, gentau_pt);
    //    
    //}
	auto tau_pt_ov_jet_pt_2d = df.Histo2D(	{((TString) "Tau_pt_vs_Tau_pt_over_Jet_pt" + stringcopy), "", 	600u, 0, 3000, 60u, 0,3.}, 	"sel_Tau_pt",	"Tau_pt_over_Jet_pt", 				weight_column);
	auto tau_pt_ov_jet_pt_2d_barrel = df.Histo2D(	{((TString) "Tau_pt_vs_Tau_pt_over_Jet_pt_barrel" + stringcopy), "", 	600u, 0, 3000, 60u, 0,3.}, 	"sel_Tau_pt_barrel",	"Tau_pt_over_Jet_pt_barrel", 				weight_column);
	auto tau_pt_ov_jet_pt_2d_endcap = df.Histo2D(	{((TString) "Tau_pt_vs_Tau_pt_over_Jet_pt_endcap" + stringcopy), "", 	600u, 0, 3000, 60u, 0,3.}, 	"sel_Tau_pt_endcap",	"Tau_pt_over_Jet_pt_endcap", 				weight_column);
	auto tau_pt_ov_jet_pt = df.Histo1D(	{((TString) "Tau_pt_over_Jet_pt" + stringcopy), "", 					6000u, 0, 1}, 		"Tau_pt_over_Jet_pt", 				weight_column);
	auto tau_pt = df.Histo1D(	{((TString) "Tau_pt" + stringcopy), "", 					6000u, 0, 6000}, 		"sel_Tau_pt", 				weight_column);
	auto tau_eta = df.Histo1D(	{((TString) "Tau_eta" + stringcopy), "", 					100u, -5, 5}, 			"sel_Tau_eta", 				weight_column);
	auto tau_phi = df.Histo1D(	{((TString) "Tau_phi" + stringcopy), "", 					100u, -3.2, 3.2}, 		"sel_Tau_phi", 				weight_column);
	auto muon_pt = df.Histo1D(	{((TString) "Muon_pt" + stringcopy), "", 					6000u, 0, 6000}, 		"sel_Muon_pt", 				weight_column);
	auto muon_eta = df.Histo1D(	{((TString) "Muon_eta" + stringcopy), "", 				100u, -5, 5}, 			"sel_Muon_eta", 			weight_column);
	auto muon_phi = df.Histo1D(	{((TString) "Muon_phi" + stringcopy), "", 				100u, -3.2, 3.2}, 		"sel_Muon_phi", 			weight_column);
	auto met_pt = df.Histo1D(	{((TString) "MET_pt" + stringcopy), "", 					6000u, 0, 6000}, "sel_MET_pt" + met_shift, 					weight_column);
	auto met_phi = df.Histo1D(	{((TString) "MET_phi" + stringcopy), "", 					100u, -3.2, 3.2},  "sel_MET_phi" + met_shift, 					weight_column);
	auto MT = df.Histo1D(		{((TString) "MT" + stringcopy), "", 						6000u, 0, 6000},		"MT", 						weight_column);
	auto coll_mass = df.Histo1D({((TString) "CollMass" + stringcopy), "",					10000u, 0, 10000},		"CollMass",					weight_column);
	auto nvtx = df.Histo1D(		{((TString) "nvtx" + stringcopy), "", 					80u, 0, 80}, 			"PV_npvs", 					weight_column);
	auto njets = df.Histo1D(		{((TString) "nJets" + stringcopy), "", 					40u, 0, 40}, 			"nJet", 					weight_column);
	auto nvtxgood = df.Histo1D(	{((TString) "nvtx_good" + stringcopy), "", 				80u, 0, 80}, 			"PV_npvsGood", 				weight_column);
    if( !config::runOnData){
        df = df.Define("no_pileupweight",[](const float genweight,
                                            const float prefire,
                                            const float trigger,
                                            const float muoniso,
                                            const float muonid,
                                            const float ww,
                                            const float topq,
                                            const float toppdf,
                                            const float toppt,
                                            const float fr,
                                            const float muonreco){
                                            return genweight *  prefire * trigger * muoniso * muonid * ww * topq * toppdf * toppt * muonreco * fr;
                                                                        },{"genWeight","PrefiringWeight","TriggerScaleFactor","MuonISOScaleFactor","MuonIDScaleFactor","WWshape","TopQ","TopPDF","top_pt_weight","FakeRate","MuonReco"});
        auto nvtxgood_nopu = df.Histo1D(	{((TString) "nvtx_good_nopu" + stringcopy), "", 				80u, 0, 80}, 			"PV_npvsGood", 				"no_pileupweight");
        hist_dict.emplace(name, nvtxgood_nopu);
    }
    else{
        auto nvtxgood_nopu = df.Histo1D(	{((TString) "nvtx_good_nopu" + stringcopy), "", 				80u, 0, 80}, 			"PV_npvsGood", 			weight_column);
        hist_dict.emplace(name, nvtxgood_nopu);
    }
	//hist_dict.emplace(name, tau_genpartflav);
	hist_dict.emplace(name, tau_pt_ov_jet_pt_2d);
	hist_dict.emplace(name, tau_pt_ov_jet_pt_2d_barrel);
	hist_dict.emplace(name, tau_pt_ov_jet_pt_2d_endcap);
	hist_dict.emplace(name, tau_pt_ov_jet_pt);
	hist_dict.emplace(name, njets);
	hist_dict.emplace(name, tau_pt);
	hist_dict.emplace(name, tau_eta);
	hist_dict.emplace(name, tau_phi);
	hist_dict.emplace(name, muon_pt);
	hist_dict.emplace(name, muon_eta);
	hist_dict.emplace(name, muon_phi);
	hist_dict.emplace(name, met_pt);
	hist_dict.emplace(name, met_phi);
	hist_dict.emplace(name, MT);
	hist_dict.emplace(name, coll_mass);
	hist_dict.emplace(name, nvtx);
	hist_dict.emplace(name, nvtxgood);
};
void create_fr_hists(	RNode df,
                                std::string name,
                                std::string weight_column) {
    std::string stringcopy = weight_column;
    stringcopy.erase(0,12);
    //if( !config::runOnData){
	//    auto gentau_pt = df.Histo1D(	{((TString) "GenVisTau_pt" + stringcopy), "", 			6000u, 0, 6000}, 					"GenVisTau_pt",     weight_column);
	//    hist_dict.emplace(name, gentau_pt);
    //    
    //}
	auto tau_pt_ov_jet_pt = df.Histo1D(	{((TString) "Tau_pt_over_Jet_pt" + stringcopy), "", 					6000u, 0, 1}, 		"Tau_pt_over_Jet_pt", 				weight_column);
	auto tau_pt = df.Histo1D(	{((TString) "Tau_pt" + stringcopy), "", 					6000u, 0, 6000}, 		"sel_Tau_pt", 				weight_column);
	auto tau_eta = df.Histo1D(	{((TString) "Tau_eta" + stringcopy), "", 					100u, -5, 5}, 			"sel_Tau_eta", 				weight_column);
	auto tau_phi = df.Histo1D(	{((TString) "Tau_phi" + stringcopy), "", 					100u, -3.2, 3.2}, 		"sel_Tau_phi", 				weight_column);
	auto muon_pt = df.Histo1D(	{((TString) "Muon_pt" + stringcopy), "", 					6000u, 0, 6000}, 		"sel_Muon_pt", 				weight_column);
	auto muon_eta = df.Histo1D(	{((TString) "Muon_eta" + stringcopy), "", 				100u, -5, 5}, 			"sel_Muon_eta", 			weight_column);
	auto muon_phi = df.Histo1D(	{((TString) "Muon_phi" + stringcopy), "", 				100u, -3.2, 3.2}, 		"sel_Muon_phi", 			weight_column);
	auto met_pt = df.Histo1D(	{((TString) "MET_pt" + stringcopy), "", 					6000u, 0, 6000},  	"sel_MET_pt" , 					weight_column);
	auto met_phi = df.Histo1D(	{((TString) "MET_phi" + stringcopy), "", 					100u, -3.2, 3.2}, "sel_MET_phi", 					weight_column);
	auto MT = df.Histo1D(		{((TString) "MT" + stringcopy), "", 						6000u, 0, 6000},		"MT", 						weight_column);
	auto coll_mass = df.Histo1D({((TString) "CollMass" + stringcopy), "",					10000u, 0, 10000},		"CollMass",					weight_column);
	auto nvtx = df.Histo1D(		{((TString) "nvtx" + stringcopy), "", 					80u, 0, 80}, 			"PV_npvs", 					weight_column);
	auto njets = df.Histo1D(		{((TString) "nJets" + stringcopy), "", 					40u, 0, 40}, 			"nJet", 					weight_column);
	auto nvtxgood = df.Histo1D(	{((TString) "nvtx_good" + stringcopy), "", 				80u, 0, 80}, 			"PV_npvsGood", 				weight_column);
	auto sphericity = df.Histo1D(	{((TString) "Sphericity" + stringcopy), "", 					100u, 0, 1.}, 		"Sphericity", 				weight_column);
	auto sphericity_leptons = df.Histo1D(	{((TString) "Sphericity_leptons" + stringcopy), "", 					100u, 0, 1.}, 		"Sphericity_leptons", 				weight_column);
	//hist_dict.emplace(name, tau_genpartflav);
	hist_dict.emplace(name, tau_pt_ov_jet_pt);
	hist_dict.emplace(name, njets);
	hist_dict.emplace(name, tau_pt);
	hist_dict.emplace(name, tau_eta);
	hist_dict.emplace(name, tau_phi);
	hist_dict.emplace(name, muon_pt);
	hist_dict.emplace(name, muon_eta);
	hist_dict.emplace(name, muon_phi);
	hist_dict.emplace(name, met_pt);
	hist_dict.emplace(name, met_phi);
	hist_dict.emplace(name, MT);
	hist_dict.emplace(name, coll_mass);
	hist_dict.emplace(name, nvtx);
	hist_dict.emplace(name, nvtxgood);
	hist_dict.emplace(name, sphericity);
	hist_dict.emplace(name, sphericity_leptons);
};
void create_reso_hist( RNode df,
                       std::string name,
                       std::string weight_column) {
    std::string stringcopy = weight_column;
    if( !config::runOnData){
	    auto reco_reso = df.Histo1D(	{"Reco_coll_mass_reso", "", 			200u, -100, 100}, 					"Reco_coll_mass_reso",     weight_column);
	    hist_dict.emplace("Reso", reco_reso);
        
    }

                       }
// fill hist function - fills hists for each stage
void create_gen_hists( RNode df,
                       std::string name
					   ) {
    auto coll_mass_gen = df.Histo1D({((TString) "CollMass_gen" ), "",					6000u, 0, 6000},		"CollMass_gen",					"genWeight");
    auto coll_mass_gen_sel = df.Histo1D({((TString) "CollMass_gen_sel" ), "",					6000u, 0, 6000},		"CollMass_gen_sel",					"genWeight");
   // auto tau_pdgid = df.Histo1D({((TString) "sel_Tau_pdgId" ), "",					100, -50, 50},		"sel_Tau_pdgId",					"genWeight");
   // auto tau_genpt = df.Histo1D({((TString) "sel_Tau_genPt" ), "",					3000u, 0, 3000},		"sel_Tau_genPt",					"genWeight");
    hist_dict.emplace(name, coll_mass_gen);
    hist_dict.emplace(name, coll_mass_gen_sel);
   // hist_dict.emplace("Generation", tau_pdgid);
   // hist_dict.emplace("Generation", tau_genpt);
}
void fill_gen_hists( RNode df){
    auto split_true = df.Filter([](const float & genPartFlav){
        return (genPartFlav == 5);
    }, {"sel_Tau_genPartFlav"});
    
    auto split_false = df.Filter([](const float & genPartFlav){
        return !(genPartFlav == 5);
    }, {"sel_Tau_genPartFlav"});
    create_gen_hists(df,"Generation_nosplit");
    create_gen_hists(split_true, "Generation");
    create_gen_hists(split_false, "Generation_split");

}
                        
void create_hists(	RNode df,
					std::string name,
                    std::string muon_shift,
					std::string tau_shift, 
					std::string met_shift,
                    std::string weight_column) {
    std::string stringcopy = weight_column + muon_shift + tau_shift + met_shift;
    //if (tau_shift == "_ES" && (met_shift == "_jer" || met_shift == "_nom") && muon_shift == "") {
    if (tau_shift == "_ES" && (met_shift == "_jer" || met_shift == "_nom" || met_shift == "") && muon_shift == "") {
        stringcopy = weight_column;
    }
    //else if(tau_shift == "_ES" && (met_shift == "_jer" || met_shift == "_nom") && muon_shift != ""){
    else if(tau_shift == "_ES" && (met_shift == "_jer" || met_shift == "_nom" || met_shift == "") && muon_shift != ""){
        stringcopy = weight_column + muon_shift;
    }
    stringcopy.erase(0,12);
    std::string nopuweight = weight_column + "/ pileup_weight";
    //if( !config::runOnData){
	//    auto gentau_pt = df.Histo1D(	{((TString) "GenVisTau_pt" + stringcopy), "", 			6000u, 0, 6000}, 					"GenVisTau_pt",     weight_column);
	//    hist_dict.emplace(name, gentau_pt);
    //    
    //}
	auto tau_pt = df.Histo1D(	{((TString) "Tau_pt" + stringcopy), "", 					6000u, 0, 6000}, 		"sel_Tau_pt", 				weight_column);
	auto tau_eta = df.Histo1D(	{((TString) "Tau_eta" + stringcopy), "", 					100u, -5, 5}, 			"sel_Tau_eta", 				weight_column);
	auto tau_phi = df.Histo1D(	{((TString) "Tau_phi" + stringcopy), "", 					100u, -3.2, 3.2}, 		"sel_Tau_phi", 				weight_column);
	auto sphericity = df.Histo1D(	{((TString) "Sphericity" + stringcopy), "", 					100u, 0, 1.}, 		"Sphericity", 				weight_column);
	auto sphericity_leptons = df.Histo1D(	{((TString) "Sphericity_leptons" + stringcopy), "", 					100u, 0, 1.}, 		"Sphericity_leptons", 				weight_column);
	auto muon_pt = df.Histo1D(	{((TString) "Muon_pt" + stringcopy), "", 					6000u, 0, 6000}, 		"sel_Muon_pt", 				weight_column);
	auto muon_pt_PF = df.Histo1D(	{((TString) "Muon_pt_PF" + stringcopy), "", 					6000u, 0, 6000}, 		"sel_Muon_PF_pt", 				weight_column);
	auto muon_PFcand = df.Histo1D(	{((TString) "Muon_PFcand" + stringcopy), "", 					5u, 0, 5}, 		"sel_Muon_isPF", 				weight_column);
	auto muon_eta = df.Histo1D(	{((TString) "Muon_eta" + stringcopy), "", 				100u, -5, 5}, 			"sel_Muon_eta", 			weight_column);
	auto muon_phi = df.Histo1D(	{((TString) "Muon_phi" + stringcopy), "", 				100u, -3.2, 3.2}, 		"sel_Muon_phi", 			weight_column);
	auto met_pt = df.Histo1D(	{((TString) "MET_pt" + stringcopy), "", 					6000u, 0, 6000}, met_branch_name  + "_pt" + met_shift, 					weight_column);
	auto met_phi = df.Histo1D(	{((TString) "MET_phi" + stringcopy), "", 					100u, -3.2, 3.2}, 		met_branch_name + "_phi" + met_shift, 					weight_column);
	auto met_pt_uncorrected = df.Histo1D(	{((TString) "MET_pt_uncorrected" + stringcopy), "", 					6000u, 0, 6000}, "MET_pt", 					weight_column);
	auto met_phi_uncorrected = df.Histo1D(	{((TString) "MET_phi_uncorrected" + stringcopy), "", 					100u, -3.2, 3.2}, 		"MET_phi", 					weight_column);
	auto MT = df.Histo1D(		{((TString) "MT" + stringcopy), "", 						6000u, 0, 6000},		"MT", 						weight_column);
	auto coll_mass = df.Histo1D({((TString) "CollMass" + stringcopy), "",					10000u, 0, 10000},		"CollMass",					weight_column);
	auto coll_mass_full = df.Histo1D({((TString) "CollMass_full" + stringcopy), "",					10000u, 0, 10000},		"CollMass_full",					weight_column);
	auto nvtx = df.Histo1D(		{((TString) "nvtx" + stringcopy), "", 					80u, 0, 80}, 			"PV_npvs", 					weight_column);
	auto njets = df.Histo1D(		{((TString) "nJets" + stringcopy), "", 					40u, 0, 40}, 			"nJet", 					weight_column);
    if( !config::runOnData){
        df = df.Define("no_pileupweight",[](const float genweight,
                                            const float taujet,
                                            const float tauele,
                                            const float taumuon,
                                            const float prefire,
                                            const float trigger,
                                            const float muoniso,
                                            const float muonid,
                                            const float ww,
                                            const float topq,
                                            const float toppdf,
                                            const float toppt,
                                            const float muonreco){
                                            return genweight * taujet * tauele * taumuon * prefire * trigger * muoniso * muonid * ww * topq * toppdf * toppt * muonreco;
                                                                        },{"genWeight","TauJetFakeScaleFactor__","TauEleFakeScaleFactor__","TauMuonFakeScaleFactor__","PrefiringWeight","TriggerScaleFactor","MuonISOScaleFactor","MuonIDScaleFactor","WWshape","TopQ","TopPDF","top_pt_weight","MuonReco"});
        auto nvtxgood_nopu = df.Histo1D(	{((TString) "nvtx_good_nopu" + stringcopy), "", 				80u, 0, 80}, 			"PV_npvsGood", 				"no_pileupweight");
        hist_dict.emplace(name, nvtxgood_nopu);
    }
    else{
        auto nvtxgood_nopu = df.Histo1D(	{((TString) "nvtx_good_nopu" + stringcopy), "", 				80u, 0, 80}, 			"PV_npvsGood", 			weight_column);
        hist_dict.emplace(name, nvtxgood_nopu);
    }
	auto nvtxgood = df.Histo1D(	{((TString) "nvtx_good" + stringcopy), "", 				80u, 0, 80}, 			"PV_npvsGood", 				weight_column);
	auto multiplicity = df.Histo1D(	{((TString) "multiplicity" + stringcopy), "", 				30u, 0, 30}, 			"multiplicity", 				weight_column);
	//hist_dict.emplace(name, tau_genpartflav);
	hist_dict.emplace(name, njets);
	hist_dict.emplace(name, tau_pt);
	hist_dict.emplace(name, tau_eta);
	hist_dict.emplace(name, tau_phi);
	hist_dict.emplace(name, muon_pt);
	hist_dict.emplace(name, muon_pt_PF);
	hist_dict.emplace(name, muon_PFcand);
	hist_dict.emplace(name, muon_eta);
	hist_dict.emplace(name, muon_phi);
	hist_dict.emplace(name, met_pt);
	hist_dict.emplace(name, met_phi);
	hist_dict.emplace(name, met_pt_uncorrected);
	hist_dict.emplace(name, met_phi_uncorrected);
	hist_dict.emplace(name, MT);
	hist_dict.emplace(name, coll_mass);
	hist_dict.emplace(name, coll_mass_full);
	hist_dict.emplace(name, nvtx);
	hist_dict.emplace(name, nvtxgood);
	hist_dict.emplace(name, multiplicity);
	hist_dict.emplace(name, sphericity);
	hist_dict.emplace(name, sphericity_leptons);
};

void fill_stage_with_syst(  RNode df, 
                            std::string stage_name,
                            std::string muon_shift,
                            std::string tau_shift,
                            std::string met_shift,
                            std::vector< std::string > list_of_weights) {
    
    bool default_run = (tau_shift == "_ES" && (met_shift == "_jer" || met_shift == "_nom") && muon_shift == "");
    
    
    if (default_run) {                           
        for (const auto& it : list_of_weights) {
            std::string folder_name = stage_name + "_nosplit";
            if ( it.find("total_weight_") != std::string::npos ) {
                folder_name += "/Systematics";
            }
            if ( it.find("pdf_weight") != std::string::npos ) {
                folder_name += "/PDF";
            }
            create_hists(df, folder_name, muon_shift, tau_shift, met_shift, it);
        }
    } else {
        std::string folder_name = stage_name + "_nosplit/Systematics";
        create_hists(df, folder_name, muon_shift, tau_shift, met_shift, "total_weight");
    }
    
    if (!config::runOnData) {
        
        auto split_true = df.Filter([](const float & genPartFlav){
            return !(genPartFlav == 0);
        }, {"sel_Tau_genPartFlav"});
        
        auto split_false = df.Filter([](const float & genPartFlav){
            return (genPartFlav == 0);
        }, {"sel_Tau_genPartFlav"});

        auto split_genIdx = df.Filter([](const int & genIdx,
                                         const float & genPartFlav){
            return (genIdx == -1) and !(genPartFlav == 0);
        }, {"sel_Muon_genIdx", "sel_Tau_genPartFlav"});
        
        if (default_run) {                           
            for (const auto& it : list_of_weights) {
                std::string folder_name_tau = stage_name;
                std::string folder_name_notau = stage_name + "_inv";
                std::string folder_name_nofakes = stage_name + "_nofakes";
                if ( it.find("total_weight_") != std::string::npos ) {
                    folder_name_tau += "/Systematics";
                    folder_name_notau += "/Systematics";
                }
                if ( it.find("pdf_weight") != std::string::npos ) {
                    folder_name_tau += "/PDF";
                    folder_name_notau += "/PDF";
                }
                create_hists(split_true, folder_name_tau, muon_shift, tau_shift, met_shift, it);
                create_hists(split_false, folder_name_notau, muon_shift, tau_shift, met_shift, it);
                create_hists(split_genIdx, folder_name_nofakes, muon_shift, tau_shift, met_shift, it);
            }
        } else {
            std::string folder_name_tau = stage_name + "/Systematics";
            std::string folder_name_notau = stage_name + "_inv/Systematics";
            create_hists(split_true, folder_name_tau, muon_shift, tau_shift, met_shift, "total_weight");
            create_hists(split_false, folder_name_notau, muon_shift, tau_shift, met_shift, "total_weight");
        }
        
    } else {
        create_hists(df, stage_name, muon_shift, tau_shift, met_shift, "total_weight");
    }
};
bool gen_match(	const rvec<int>& gen_pdgId,
				const rvec<float>& gen_part_pt,
				const rvec<float>& gen_part_eta,
				const rvec<float>& gen_part_phi,
				const rvec<float>& gen_part_mass,
				const int pdgId,
				const float& part_pt,
				const float& part_eta,
				const float& part_phi,
				const float& part_mass
				) {
	TLorentzVector p1, p2; 
	p1.SetPtEtaPhiM(part_pt, part_eta, part_phi, part_mass);
	for (uint i = 0; i < gen_pdgId.size(); i++) {
		if (gen_pdgId[i] == pdgId) {
			p2.SetPtEtaPhiM(gen_part_pt[i], gen_part_eta[i], gen_part_phi[i], gen_part_mass[i]);
			if (p1.DeltaR(p2) < 0.3) 
				return true;
		}
	}
	
	return false;
};
int gen_match_idx(	
                    const rvec<float>& gen_part_pt,
                    const rvec<float>& gen_part_eta,
                    const rvec<float>& gen_part_phi,
                    const rvec<float>& gen_part_mass,
                    const float& part_pt,
                    const float& part_eta,
                    const float& part_phi,
                    const float& part_mass
				) {
	TLorentzVector p1, p2; 
	p1.SetPtEtaPhiM(part_pt, part_eta, part_phi, part_mass);
	for (uint i = 0; i < gen_part_pt.size(); i++) {
			p2.SetPtEtaPhiM(gen_part_pt[i], gen_part_eta[i], gen_part_phi[i], gen_part_mass[i]);
			if (p1.DeltaR(p2) < 0.3) 
				return i;
		}
	
	
	return 0;
};

void fill_datadriven_with_syst(  RNode df, 
                            std::string stage_name,
                            std::string muon_shift,
                            std::string tau_shift,
                            std::string met_shift,
                            std::vector< std::string > list_of_weights) {
    
    bool default_run = (tau_shift == "_ES" && (met_shift == "_jer" || met_shift == "_nom") && muon_shift == "");
    
    
    if (default_run) {                           
        for (const auto& it : list_of_weights) {
            std::string folder_name = stage_name + "_nosplit";
            if ( it.find("total_weight_") != std::string::npos ) {
                folder_name += "/Systematics";
            }
            if ( it.find("pdf_weight") != std::string::npos ) {
                folder_name += "/PDF";
            }
            create_datadriven_hists(df, folder_name, muon_shift, tau_shift, met_shift, it);
        }
    } else {
        std::string folder_name = stage_name + "_nosplit/Systematics";
        create_datadriven_hists(df, folder_name, muon_shift, tau_shift, met_shift, "total_weight");
    }
    
    if (!config::runOnData) {
        
        auto split_true = df.Filter([](const float & genPartFlav){
            return (genPartFlav == 5);
        }, {"sel_Tau_genPartFlav"});
        
        auto split_false = df.Filter([](const float & genPartFlav){
            return !(genPartFlav == 5);
        }, {"sel_Tau_genPartFlav"});
        
        if (default_run) {                           
            for (const auto& it : list_of_weights) {
                std::string folder_name_tau = stage_name;
                std::string folder_name_notau = stage_name + "_inv";
                if ( it.find("total_weight_") != std::string::npos ) {
                    folder_name_tau += "/Systematics";
                    folder_name_notau += "/Systematics";
                }
                if ( it.find("pdf_weight") != std::string::npos ) {
                    folder_name_tau += "/PDF";
                    folder_name_notau += "/PDF";
                }
                create_datadriven_hists(split_true, folder_name_tau, muon_shift, tau_shift, met_shift, it);
                create_datadriven_hists(split_false, folder_name_notau, muon_shift, tau_shift, met_shift, it);
            }
        } else {
            std::string folder_name_tau = stage_name + "/Systematics";
            std::string folder_name_notau = stage_name + "_inv/Systematics";
            create_datadriven_hists(split_true, folder_name_tau, muon_shift, tau_shift, met_shift, "total_weight");
            create_datadriven_hists(split_false, folder_name_notau, muon_shift, tau_shift, met_shift, "total_weight");
        }
        
    } else {
        create_datadriven_hists(df, stage_name, muon_shift, tau_shift, met_shift, "total_weight");
    }
};
// Calculate datadriven distributions
RNode calc_datadriven(RNode df,
                      std::string muon_shift,
                      std::string tau_shift,
                      std::string met_shift) {
	//// begin evaluating the branch with non isolated taus for datadriven background estimation in QCD and ZToNuNu						 
	auto noniso_tau = df.Define("NonIso_Tau", tau_acceptance_and_id_and_dm_noniso,	{"Tau_pt" + tau_shift, "Tau_eta" + tau_shift, "Tau_dz","Tau_charge", config::tau_dm, "Tau_decayMode", "Tau_jetIdx", config::tau_iso, config::tau_antiEle, config::tau_antiMuon});
	//										
	auto iso_tau = df.Define("Iso_Tau", tau_acceptance_and_id_and_dm,	{"Tau_pt" + tau_shift, "Tau_eta" + tau_shift, "Tau_dz", "Tau_charge", config::tau_dm, "Tau_decayMode", "Tau_jetIdx", config::tau_iso, config::tau_antiEle, config::tau_antiMuon});
	//										
	auto evntselect_noniso_tau = noniso_tau.Filter(nparticle_cut, {"NonIso_Tau"}, "select_taus")
                                .Filter(nparticle_cut, {"Muon_mask"}, "select_muons")
                                .Filter(nparticle_veto, {"Electron_mask"}, "select_eles");
	auto evntselect_iso_tau = iso_tau.Filter(nparticle_cut, {"Iso_Tau"}, "select_taus")
                                .Filter(nparticle_cut, {"Muon_mask"}, "select_muons")
                                .Filter(nparticle_veto, {"Electron_mask"}, "select_eles");
    auto df_col_idx = evntselect_noniso_tau.Define("col_idx", col_idx, {"Tau_pt" + tau_shift, "Muon_pt_corr" + muon_shift, "Tau_eta" + tau_shift, "Muon_eta", "Tau_phi" + tau_shift, "Muon_phi", "Tau_mass" + tau_shift, "Muon_mass", met_branch_name + "_pt" + met_shift, met_branch_name + "_phi" + met_shift , "NonIso_Tau", "Muon_mask"});   
	//
	//// Now that we have our real tau, write it to special columns to handle it more easily
        auto defquants = df_col_idx.Define("sel_Tau_pt", selected_part_col_idx_tau, {"Tau_pt" + tau_shift, "NonIso_Tau", "col_idx"})
                                   .Define("sel_Tau_pt_barrel", selected_part_col_idx_tau_barrel, {"Tau_pt" + tau_shift, "NonIso_Tau", "col_idx","Tau_eta" + tau_shift})
                                   .Define("sel_Tau_pt_endcap", selected_part_col_idx_tau_endcap, {"Tau_pt" + tau_shift, "NonIso_Tau", "col_idx","Tau_eta" + tau_shift})
                                   .Define("sel_Tau_eta", selected_part_col_idx_tau, {"Tau_eta" + tau_shift, "NonIso_Tau", "col_idx"})
                                   .Define("sel_Tau_phi", selected_part_col_idx_tau, {"Tau_phi" + tau_shift, "NonIso_Tau", "col_idx"})
                                   .Define("sel_Tau_mass", selected_part_col_idx_tau, {"Tau_mass" + tau_shift, "NonIso_Tau", "col_idx"})
                                   .Define("sel_Tau_jetIdx", selected_jet_quant, {"Tau_jetIdx", "NonIso_Tau", "col_idx"})
                                   .Define("sel_Tau_pt_norm", selected_part_col_idx_tau, {"Tau_pt" , "Tau_mask", "col_idx"})
                                   .Define("sel_Tau_eta_norm", selected_part_col_idx_tau, {"Tau_eta" , "Tau_mask", "col_idx"})
                                   .Define("sel_Tau_phi_norm", selected_part_col_idx_tau, {"Tau_phi" , "Tau_mask", "col_idx"})
                                   .Define("sel_Tau_mass_norm", selected_part_col_idx_tau, {"Tau_mass", "Tau_mask", "col_idx"})
                                   .Define("sel_Muon_PF_pt", selected_part_col_idx_muon, {"Muon_pt", "Muon_mask", "col_idx"})
                                   .Define("sel_Muon_isPF", selected_part_col_idx_muon_bool, {"Muon_isPFcand", "Muon_mask", "col_idx"})
                                   .Define("sel_Muon_pt", selected_part_col_idx_muon, {"Muon_pt_corr" + muon_shift, "Muon_mask", "col_idx"})
                                   .Define("sel_Muon_eta", selected_part_col_idx_muon, {"Muon_eta", "Muon_mask", "col_idx"})
                                   .Define("sel_Muon_phi", selected_part_col_idx_muon, {"Muon_phi", "Muon_mask", "col_idx"})
                                   .Define("sel_Muon_mass", selected_part_col_idx_muon, {"Muon_mass", "Muon_mask", "col_idx"})
                                   .Define("sel_MET_pt" + met_shift, correct_met_single,{"sel_Muon_PF_pt","sel_Muon_pt","sel_Muon_eta","sel_Muon_phi","sel_Muon_mass","sel_Muon_isPF","sel_Tau_pt","sel_Tau_eta","sel_Tau_phi","sel_Tau_mass","sel_Tau_pt_norm","sel_Tau_eta_norm","sel_Tau_phi_norm","sel_Tau_mass_norm",met_branch_name + "_pt" + met_shift,met_branch_name + "_phi" + met_shift,"run","PV_npvsGood"})
                                   .Define("sel_MET_phi" + met_shift, correct_met_single_phi,{"sel_Muon_PF_pt","sel_Muon_pt","sel_Muon_eta","sel_Muon_phi","sel_Muon_mass","sel_Muon_isPF","sel_Tau_pt","sel_Tau_eta","sel_Tau_phi","sel_Tau_mass","sel_Tau_pt_norm","sel_Tau_eta_norm","sel_Tau_phi_norm","sel_Tau_mass_norm",met_branch_name + "_pt" + met_shift,met_branch_name + "_phi" + met_shift,"run","PV_npvsGood"});
    if (!config::runOnData) {
        defquants = defquants.Define("sel_Tau_genPartFlav", selected_part_col_idx_tau_flav, {"Tau_genPartFlav", "NonIso_Tau", "col_idx"});
    }
	//auto defquants = evntselect_noniso_tau.Define("sel_Tau_pt", selected_part_quant, {"Tau_pt_ES", "NonIso_Tau"})
	//						   .Define("sel_Tau_eta", selected_part_quant, {"Tau_eta_ES", "NonIso_Tau"})
	//						   .Define("sel_Tau_phi", selected_part_quant, {"Tau_phi_ES", "NonIso_Tau"})
	//						   .Define("sel_Tau_jetIdx", selected_jet_quant, {"Tau_jetIdx", "NonIso_Tau"})
    //                           .Define("sel_Muon_pt", selected_part_quant, {"Muon_tP_pt", "Muon_mask"})
    //                           .Define("sel_Muon_eta", selected_part_quant, {"Muon_eta", "Muon_mask"})
    //                           .Define("sel_Muon_phi", selected_part_quant, {"Muon_phi", "Muon_mask"});

    //auto dimuon_cut = defquants.Filter(dimuonpair_cut,{"Muon_tP_pt","Muon_eta","Muon_phi", "DiMuon_mask"},"Di Muon Pair Cut");
	//// Calculate MT distribution
    auto dimuon_cut = defquants.Filter(dimuonpair_cut,{"Muon_pt_corr" + muon_shift,"Muon_eta","Muon_phi", "DiMuon_mask"},"Di Muon Pair Cut");
    auto delR_cut = dimuon_cut.Filter(mutaudelR_cut,{"sel_Muon_eta","sel_Muon_phi", "sel_Tau_eta", "sel_Tau_phi"},"Mu Tau Del R Cut");
    auto df_jet_matching = delR_cut.Define("matched_jet",jet_matching, {"sel_Tau_pt","sel_Tau_eta","sel_Tau_phi","sel_Tau_mass","Jet_pt","Jet_eta","Jet_phi","Jet_mass"});
    auto jet_matching_cut = df_jet_matching.Filter("matched_jet == true","jet matching cut");


    auto corr_met_df  = jet_matching_cut;
	auto mtcalc = corr_met_df.Define("MT", mass_transv, {"sel_Muon_pt", "sel_Muon_phi", "sel_MET_pt" + met_shift, "sel_MET_phi" + met_shift});// .Define("sel_Tau_pt_ov_Jet_pt_NonIso", [](const rvec<float> tau_pt, 
                                                                                                           //                                    const rvec<Int_t> jetidx,
                                                                                                           //                                    const rvec<float> jet_pt){
        //return tau_pt/jet_pt[jetidx];                                             
//    },
  //                  {"sel_Tau_pt","sel_Tau_jetIdx","Jet_pt"});
    auto df_coll = mtcalc.Define("CollMass", collinear_mass, {"sel_Tau_pt", "sel_Muon_pt", "sel_Tau_eta", "sel_Muon_eta", "sel_Tau_phi", "sel_Muon_phi", "sel_Tau_mass", "sel_Muon_mass", "sel_MET_pt" + met_shift, "sel_MET_phi" + met_shift}).Define("CollMass_alt", collinear_mass_alt, {"sel_Tau_pt", "sel_Muon_pt", "sel_Tau_eta", "sel_Muon_eta", "sel_Tau_phi", "sel_Muon_phi", "sel_Tau_mass", "sel_Muon_mass", "sel_MET_pt" + met_shift, "sel_MET_phi" + met_shift}).Define("Tau_pt_over_Jet_pt", tau_pt_over_jet_pt,{"col_idx","Tau_pt" + tau_shift,"Tau_jetIdx","Jet_pt_nom","NonIso_Tau"}).Define("Tau_pt_over_Jet_pt_barrel", tau_pt_over_jet_pt_barrel,{"col_idx","Tau_pt" + tau_shift, "Tau_jetIdx", "Jet_pt_nom", "NonIso_Tau", "Tau_eta" + tau_shift}).Define("Tau_pt_over_Jet_pt_endcap", tau_pt_over_jet_pt_endcap,{"col_idx", "Tau_pt" + tau_shift, "Tau_jetIdx", "Jet_pt_nom", "NonIso_Tau", "Tau_eta" + tau_shift});
    RNode total_weights = df_coll;
    std::vector < std::string > list_of_weights;

    if (config::runOnData) {
        total_weights = total_weights.Define("total_weight", "1.0");
        list_of_weights.push_back("total_weight");
    } else {    
        total_weights = init_PDFs(total_weights);
        
        std::vector < std::string > default_weights = {
            "pileup_weight",
            "genWeight",
            "PrefiringWeight",
            "TriggerScaleFactor",
            "MuonISOScaleFactor",
            "MuonIDScaleFactor",
            "WWshape",
            "TopQ",
            "TopPDF",
            "top_pt_weight",
            "MuonReco"
        };


    total_weights = total_weights.Define("MuonISOScaleFactor", [](const rvec<float>& muon_pt,
                                                                  const rvec<float>& muon_eta,
                                                                  const rvec<bool>& muon_mask
                                                                  )
                                                                  {
                                                                     return muon_iso_scale_factor(muon_pt, muon_eta, muon_mask, "");
                                                                  }, {"Muon_pt_corr", "Muon_eta", "Muon_mask"})
                                 .Define("MuonISOScaleFactorUp", [](const rvec<float>& muon_pt,
                                                                  const rvec<float>& muon_eta,
                                                                  const rvec<bool>& muon_mask
                                                                  )
                                                                  {
                                                                     return muon_iso_scale_factor(muon_pt, muon_eta, muon_mask, "Up");
                                                                  }, {"Muon_pt_corr", "Muon_eta", "Muon_mask"})
                                 .Define("MuonISOScaleFactorDown", [](const rvec<float>& muon_pt,
                                                                      const rvec<float>& muon_eta,
                                                                      const rvec<bool>& muon_mask
                                                                      )
                                                                      {
                                                                         return muon_iso_scale_factor(muon_pt, muon_eta, muon_mask, "Down");
                                                                      }, {"Muon_pt_corr", "Muon_eta", "Muon_mask"})
                                 .Define("MuonIDScaleFactor", []( const rvec<float>& muon_pt,
                                                                  const rvec<float>& muon_eta,
                                                                  const rvec<bool>& muon_mask
                                                                  )
                                                                  {
                                                                     return muon_id_scale_factor(muon_pt, muon_eta, muon_mask, "");
                                                                  }, {"Muon_pt_corr", "Muon_eta", "Muon_mask"})
                                 .Define("MuonIDScaleFactorUp", []( const rvec<float>& muon_pt,
                                                                  const rvec<float>& muon_eta,
                                                                  const rvec<bool>& muon_mask
                                                                  )
                                                                  {
                                                                     return muon_id_scale_factor(muon_pt, muon_eta, muon_mask, "Up");
                                                                  }, {"Muon_pt_corr", "Muon_eta", "Muon_mask"})
                                 .Define("MuonIDScaleFactorDown", []( const rvec<float>& muon_pt,
                                                                  const rvec<float>& muon_eta,
                                                                  const rvec<bool>& muon_mask
                                                                  )
                                                                  {
                                                                     return muon_id_scale_factor(muon_pt, muon_eta, muon_mask, "Down");
                                                                  }, {"Muon_pt_corr", "Muon_eta", "Muon_mask"})
                                 .Define("WWshape", [](){ return (float) 1.;})
                                 .Define("WWshapeUp", []( const float& mll
                                                                  )
                                                                  {
                                                                     if (config::wwuncertainty){
                                                                         float weight = 1 + 1 - (0.993-0.00002001 * mll + 0.00000000283 * mll * mll);
                                                                         return weight ;
                                                                     }
                                                                     return (float) 1.;
                                                                  }, {"Mll"})
                                 .Define("WWshapeDown", []( const float& mll
                                                                  )
                                                                  {
                                                                     if (config::wwuncertainty){
                                                                         float weight = (0.993-0.00002001 * mll + 0.00000000283 * mll * mll);
                                                                         return weight;
                                                                     }
                                                                     return (float) 1.;
                                                                  }, {"Mll"})
                                 .Define("TopQ", []( const float& Mll
                                                                  )
                                                                  {
                                                                     if (config::TT){
                                                                         return GetTopQscale(Mll,"nom");
                                                                     }
                                                                     return (float) 1.;
                                                                  }, {"Mll"})
                                 .Define("TopQUp", []( const float& Mll
                                                                  )
                                                                  {
                                                                     if (config::TT){
                                                                         return GetTopQscale(Mll,"up");
                                                                     }
                                                                     return (float) 1.;
                                                                  }, {"Mll"})
                                 .Define("TopQDown", []( const float& Mll
                                                                  )
                                                                  {
                                                                     if (config::TT){
                                                                         return GetTopQscale(Mll,"down");
                                                                     }
                                                                     return (float) 1.;
                                                                  }, {"Mll"})
                                 .Define("TopPDF", []( const float& Mll
                                                                  )
                                                                  {
                                                                     if (config::TT){
                                                                         return GetTopPDF(Mll,"nom");
                                                                     }
                                                                     return (float) 1.;
                                                                  }, {"Mll"})
                                 .Define("TopPDFUp", []( const float& Mll
                                                                  )
                                                                  {
                                                                     if (config::TT){
                                                                         return GetTopPDF(Mll,"up");
                                                                     }
                                                                     return (float) 1.;
                                                                  }, {"Mll"})
                                 .Define("TopPDFDown", []( const float& Mll
                                                                  )
                                                                  {
                                                                     if (config::TT){
                                                                         return GetTopPDF(Mll,"down");
                                                                     }
                                                                     return (float) 1.;
                                                                  }, {"Mll"})
                            .Define("TriggerScaleFactor", [](const rvec<float>& pt,
                                                             const rvec<float>& eta,
                                                             const rvec<bool>& mask
                                                             ) { 
                                                                 return trigger_sf(pt, eta,mask, "");
                                                            }, {"Muon_pt_corr", "Muon_eta", "Muon_mask"})
                            .Define("TriggerScaleFactorUp", [](const rvec<float>& pt,
                                                             const rvec<float>& eta,
                                                             const rvec<bool>& mask
                                                             ) { 
                                                                 return trigger_sf(pt, eta, mask, "Up");
                                                            }, {"Muon_pt_corr", "Muon_eta", "Muon_mask"})
                            .Define("TriggerScaleFactorDown", [](const rvec<float>& pt,
                                                             const rvec<float>& eta,
                                                             const rvec<bool>& mask
                                                             ) { 
                                                                 return trigger_sf(pt, eta, mask, "Down");
                                                            }, {"Muon_pt_corr", "Muon_eta", "Muon_mask"})
                            .Define("MuonReco", [](const rvec<float>& pt,
                                                             const rvec<float>& eta,
                                                             const rvec<float>& phi,
                                                             const rvec<float>& mass,
                                                             const rvec<bool>& mask
                                                             ) { 
                                                                 return muon_reco_eff(pt, eta, phi, mass,mask,"");
                                                            }, {"Muon_pt_corr", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_mask"});
    //define_weight(total_weights, default_weights, list_of_weights);
            //DEFINE NON-ENERGY RELATED UNCERTAINTIES
            total_weights = total_weights.Define("pileup_weightUp", [](const int& nvtx_true){ return pu_weight(nvtx_true, "Up"); }, {"Pileup_nPU"})
                                         .Define("pileup_weightDown", [](const int& nvtx_true){ return pu_weight(nvtx_true, "Down"); }, {"Pileup_nPU"})
                                         .Define("prefiring_weight_Up", [](const rvec<float>& jet_pt,
                                                                 const rvec<float>& jet_eta,
                                                                 const rvec<float>& photon_pt,
                                                                 const rvec<float>& photon_eta) { 
                                                                     return prefire_factor(jet_pt, jet_eta, photon_pt, photon_eta, "Up");
                                                                }, {"Jet_pt", "Jet_eta", "Photon_pt", "Photon_eta"})
                                         .Define("prefiring_weight_Down", [](const rvec<float>& jet_pt,
                                                                 const rvec<float>& jet_eta,
                                                                 const rvec<float>& photon_pt,
                                                                 const rvec<float>& photon_eta) { 
                                                                     return prefire_factor(jet_pt, jet_eta, photon_pt, photon_eta, "Down");
                                                            }, {"Jet_pt", "Jet_eta", "Photon_pt", "Photon_eta"});

            
        
            define_weight(total_weights, default_weights, list_of_weights);
            define_weight(total_weights, default_weights, list_of_weights, "pileup_weightUp", "pileup_weight");
            define_weight(total_weights, default_weights, list_of_weights, "pileup_weightDown", "pileup_weight");
            define_weight(total_weights, default_weights, list_of_weights, "TriggerScaleFactorUp", "TriggerScaleFactor");
            define_weight(total_weights, default_weights, list_of_weights, "TriggerScaleFactorDown", "TriggerScaleFactor");
            define_weight(total_weights, default_weights, list_of_weights, "MuonISOScaleFactorUp", "MuonISOScaleFactor");
            define_weight(total_weights, default_weights, list_of_weights, "MuonISOScaleFactorDown", "MuonISOScaleFactor");
            define_weight(total_weights, default_weights, list_of_weights, "MuonIDScaleFactorUp", "MuonIDScaleFactor");
            define_weight(total_weights, default_weights, list_of_weights, "MuonIDScaleFactorDown", "MuonIDScaleFactor");
            define_weight(total_weights, default_weights, list_of_weights, "WWshapeUp", "WWshape");
            define_weight(total_weights, default_weights, list_of_weights, "WWshapeDown", "WWshape");
            define_weight(total_weights, default_weights, list_of_weights, "TopQUp", "TopQ");
            define_weight(total_weights, default_weights, list_of_weights, "TopQDown", "TopQ");
            define_weight(total_weights, default_weights, list_of_weights, "TopPDFUp", "TopPDF");
            define_weight(total_weights, default_weights, list_of_weights, "TopPDFDown", "TopPDF");

           
            // factorization & renormalization scale uncertainty
            std::vector<std::string> columnNames = total_weights.GetColumnNames();
            if (std::find(columnNames.begin(), columnNames.end(), "LHEScaleWeight") != columnNames.end()) {
                total_weights = total_weights.Define("total_weight_ScaleWeightUp", [](const rvec<float>& weights, const float& tot_weight) {
                                                                                                    float max_val = 1.0;
                                                                                                    for (uint i = 0; i < weights.size(); i++) {
                                                                                                        if ((i==5) || (i==7)) {
                                                                                                            continue;
                                                                                                        }
                                                                                                        if (max_val < weights[i]) {
                                                                                                            max_val = weights[i];
                                                                                                        }
                                                                                                    }
                                                                                                    return max_val*tot_weight;
                                                                                    }, {"LHEScaleWeight", "total_weight"})
                                                                            .Define("total_weight_ScaleWeightDown", [](const rvec<float>& weights, const float& tot_weight) {
                                                                                                    float min_val = 1.0;
                                                                                                    for (uint i = 0; i < weights.size(); i++) {
                                                                                                        if ((i==5) || (i==7)) {
                                                                                                            continue;
                                                                                                        }
                                                                                                        if (min_val > weights[i]) {
                                                                                                            min_val = weights[i];
                                                                                                        }
                                                                                                    }
                                                                                                    return min_val*tot_weight;
                                                                                    }, {"LHEScaleWeight", "total_weight"});
                
            } else {
                total_weights = total_weights.Define("total_weight_ScaleWeightUp", "total_weight")
                                             .Define("total_weight_ScaleWeightDown", "total_weight");
            }
            list_of_weights.push_back("total_weight_ScaleWeightUp");
            list_of_weights.push_back("total_weight_ScaleWeightDown");
            // pdf weights
            std::vector<std::string> colNames = total_weights.GetColumnNames();
            if (std::find(colNames.begin(), colNames.end(), "LHEPdfWeight_def") != colNames.end()) {
                if (!config::runOnSignal) {
                    if (std::find(colNames.begin(), colNames.end(), "LHEWeight_originalXWGTUP") == colNames.end()) {
                        total_weights = total_weights.Define("LHEWeight_originalXWGTUP", [](const rvec<float>& pdf_weights){return pdf_weights[0];}, {"LHEPdfWeight_def"});
                    }
                    for (uint i = 0; i < config::pdf_nweights; i++) {
                        total_weights = total_weights.Define("total_weight_pdf_weight_" + std::to_string(i), [i](
                                const float& weight, 
                                const rvec<float>& pdf_weights, 
                                const float& orig_weight)
                            {
                                if (pdf_weights.size() != config::pdf_nweights) {
                                    return weight;
                                }
                                //~ return weight*pdf_weights[i]/orig_weight;
                                return weight*pdf_weights[i]/pdf_weights[0];
                            }, {"total_weight", "LHEPdfWeight_def", "LHEWeight_originalXWGTUP"});
                        list_of_weights.push_back("total_weight_pdf_weight_" + std::to_string(i));
                    }
                }
            }
        }
    auto df_col_idx_iso = evntselect_iso_tau.Define("col_idx", col_idx, {"Tau_pt" + tau_shift, "Muon_pt_corr" + muon_shift, "Tau_eta" + tau_shift, "Muon_eta", "Tau_phi" + tau_shift, "Muon_phi", "Tau_mass" + tau_shift, "Muon_mass", met_branch_name + "_pt" + met_shift, met_branch_name + "_phi" + met_shift , "Iso_Tau", "Muon_mask"});   
        auto defquants_iso = df_col_idx_iso.Define("sel_Tau_pt", selected_part_col_idx_tau, {"Tau_pt" + tau_shift, "Iso_Tau", "col_idx"})
                                   .Define("sel_Tau_pt_barrel", selected_part_col_idx_tau_barrel, {"Tau_pt" + tau_shift, "Iso_Tau", "col_idx","Tau_eta" + tau_shift})
                                   .Define("sel_Tau_pt_endcap", selected_part_col_idx_tau_endcap, {"Tau_pt" + tau_shift, "Iso_Tau", "col_idx","Tau_eta" + tau_shift})
                                   .Define("sel_Tau_eta", selected_part_col_idx_tau, {"Tau_eta" + tau_shift, "Iso_Tau", "col_idx"})
                                   .Define("sel_Tau_phi", selected_part_col_idx_tau, {"Tau_phi" + tau_shift, "Iso_Tau", "col_idx"})
                                   .Define("sel_Tau_mass", selected_part_col_idx_tau, {"Tau_mass" + tau_shift, "Iso_Tau", "col_idx"})
                                   .Define("sel_Tau_jetIdx", selected_jet_quant, {"Tau_jetIdx", "Iso_Tau", "col_idx"})
                                   .Define("sel_Tau_pt_norm", selected_part_col_idx_tau, {"Tau_pt" , "Tau_mask", "col_idx"})
                                   .Define("sel_Tau_eta_norm", selected_part_col_idx_tau, {"Tau_eta" , "Tau_mask", "col_idx"})
                                   .Define("sel_Tau_phi_norm", selected_part_col_idx_tau, {"Tau_phi" , "Tau_mask", "col_idx"})
                                   .Define("sel_Tau_mass_norm", selected_part_col_idx_tau, {"Tau_mass", "Tau_mask", "col_idx"})
                                   .Define("sel_Muon_PF_pt", selected_part_col_idx_muon, {"Muon_pt", "Muon_mask", "col_idx"})
                                   .Define("sel_Muon_isPF", selected_part_col_idx_muon_bool, {"Muon_isPFcand", "Muon_mask", "col_idx"})
                                   .Define("sel_Muon_pt", selected_part_col_idx_muon, {"Muon_pt_corr" + muon_shift, "Muon_mask", "col_idx"})
                                   .Define("sel_Muon_eta", selected_part_col_idx_muon, {"Muon_eta", "Muon_mask", "col_idx"})
                                   .Define("sel_Muon_phi", selected_part_col_idx_muon, {"Muon_phi", "Muon_mask", "col_idx"})
                                   .Define("sel_Muon_mass", selected_part_col_idx_muon, {"Muon_mass", "Muon_mask", "col_idx"})
                                   .Define("sel_MET_pt" + met_shift, correct_met_single,{"sel_Muon_PF_pt","sel_Muon_pt","sel_Muon_eta","sel_Muon_phi","sel_Muon_mass","sel_Muon_isPF","sel_Tau_pt","sel_Tau_eta","sel_Tau_phi","sel_Tau_mass","sel_Tau_pt_norm","sel_Tau_eta_norm","sel_Tau_phi_norm","sel_Tau_mass_norm",met_branch_name + "_pt" + met_shift,met_branch_name + "_phi" + met_shift,"run","PV_npvsGood"})
                                   .Define("sel_MET_phi" + met_shift, correct_met_single_phi,{"sel_Muon_PF_pt","sel_Muon_pt","sel_Muon_eta","sel_Muon_phi","sel_Muon_mass","sel_Muon_isPF","sel_Tau_pt","sel_Tau_eta","sel_Tau_phi","sel_Tau_mass","sel_Tau_pt_norm","sel_Tau_eta_norm","sel_Tau_phi_norm","sel_Tau_mass_norm",met_branch_name + "_pt" + met_shift,met_branch_name + "_phi" + met_shift,"run","PV_npvsGood"});
    if (!config::runOnData) {
        defquants_iso = defquants_iso.Define("sel_Tau_genPartFlav", selected_part_col_idx_tau_flav, {"Tau_genPartFlav", "Iso_Tau", "col_idx"});
    }
    auto dimuon_cut_iso = defquants_iso.Filter(dimuonpair_cut,{"Muon_pt_corr" + muon_shift,"Muon_eta","Muon_phi", "DiMuon_mask"},"Di Muon Pair Cut");
    auto delR_cut_iso = dimuon_cut_iso.Filter(mutaudelR_cut,{"sel_Muon_eta","sel_Muon_phi", "sel_Tau_eta", "sel_Tau_phi"},"Mu Tau Del R Cut");
    auto df_jet_matching_iso = delR_cut_iso.Define("matched_jet",jet_matching, {"sel_Tau_pt","sel_Tau_eta","sel_Tau_phi","sel_Tau_mass","Jet_pt","Jet_eta","Jet_phi","Jet_mass"});
    auto jet_matching_cut_iso = df_jet_matching_iso.Filter("matched_jet == true","jet matching cut");


    auto corr_met_df_iso  = jet_matching_cut_iso;
	auto mtcalc_iso = corr_met_df_iso.Define("MT", mass_transv, {"sel_Muon_pt", "sel_Muon_phi", "sel_MET_pt" + met_shift, "sel_MET_phi" + met_shift});//.Define("sel_Tau_pt_ov_Jet_pt_Iso", [](const rvec<float> tau_pt, 
                                                                                                                     //                          const rvec<Int_t> jetidx,
                                                                                                                       //                        const rvec<float> jet_pt){
//        return tau_pt/jet_pt[jetidx];                                             
  //  },
    //                {"sel_Tau_pt","sel_Tau_jetIdx","Jet_pt"});
    auto df_coll_iso = mtcalc_iso.Define("CollMass", collinear_mass, {"sel_Tau_pt", "sel_Muon_pt", "sel_Tau_eta", "sel_Muon_eta", "sel_Tau_phi", "sel_Muon_phi", "sel_Tau_mass", "sel_Muon_mass", "sel_MET_pt" + met_shift, "sel_MET_phi" + met_shift}).Define("CollMass_alt", collinear_mass_alt, {"sel_Tau_pt", "sel_Muon_pt", "sel_Tau_eta", "sel_Muon_eta", "sel_Tau_phi", "sel_Muon_phi", "sel_Tau_mass", "sel_Muon_mass", "sel_MET_pt" + met_shift, "sel_MET_phi" + met_shift}).Define("Tau_pt_over_Jet_pt", tau_pt_over_jet_pt,{"col_idx","Tau_pt" + tau_shift,"Tau_jetIdx","Jet_pt_nom","Iso_Tau"}).Define("Tau_pt_over_Jet_pt_endcap", tau_pt_over_jet_pt_endcap,{"col_idx","Tau_pt" + tau_shift, "Tau_jetIdx", "Jet_pt_nom", "Iso_Tau", "Tau_eta" + tau_shift}).Define("Tau_pt_over_Jet_pt_barrel", tau_pt_over_jet_pt_barrel,{"col_idx","Tau_pt" + tau_shift, "Tau_jetIdx", "Jet_pt_nom", "Iso_Tau", "Tau_eta" + tau_shift});  
	RNode total_weights_iso = df_coll_iso;
    std::vector < std::string > list_of_weights_iso;

    if (config::runOnData) {
        total_weights_iso = total_weights_iso.Define("total_weight", "1.0");
        list_of_weights_iso.push_back("total_weight");
    } else {    
    
    total_weights_iso = init_PDFs(total_weights_iso);
    
    
    std::vector < std::string > default_weights_iso = {
        "pileup_weight",
        "genWeight",
        "PrefiringWeight",
        "TriggerScaleFactor",
        "MuonISOScaleFactor",
        "MuonIDScaleFactor",
        "WWshape",
        "TopQ",
        "TopPDF",
        "top_pt_weight",
        "MuonReco"
    };


total_weights_iso = total_weights_iso.Define("MuonISOScaleFactor", [](const rvec<float>& muon_pt,
                                                              const rvec<float>& muon_eta,
                                                              const rvec<bool>& muon_mask
                                                              )
                                                              {
                                                                 return muon_iso_scale_factor(muon_pt, muon_eta, muon_mask, "");
                                                              }, {"Muon_pt_corr", "Muon_eta", "Muon_mask"})
                             .Define("MuonISOScaleFactorUp", [](const rvec<float>& muon_pt,
                                                              const rvec<float>& muon_eta,
                                                              const rvec<bool>& muon_mask
                                                              )
                                                              {
                                                                 return muon_iso_scale_factor(muon_pt, muon_eta, muon_mask, "Up");
                                                              }, {"Muon_pt_corr", "Muon_eta", "Muon_mask"})
                             .Define("MuonISOScaleFactorDown", [](const rvec<float>& muon_pt,
                                                                  const rvec<float>& muon_eta,
                                                                  const rvec<bool>& muon_mask
                                                                  )
                                                                  {
                                                                     return muon_iso_scale_factor(muon_pt, muon_eta, muon_mask, "Down");
                                                                  }, {"Muon_pt_corr", "Muon_eta", "Muon_mask"})
                             .Define("MuonIDScaleFactor", []( const rvec<float>& muon_pt,
                                                              const rvec<float>& muon_eta,
                                                              const rvec<bool>& muon_mask
                                                              )
                                                              {
                                                                 return muon_id_scale_factor(muon_pt, muon_eta, muon_mask, "");
                                                              }, {"Muon_pt_corr", "Muon_eta", "Muon_mask"})
                             .Define("MuonIDScaleFactorUp", []( const rvec<float>& muon_pt,
                                                              const rvec<float>& muon_eta,
                                                              const rvec<bool>& muon_mask
                                                              )
                                                              {
                                                                 return muon_id_scale_factor(muon_pt, muon_eta, muon_mask, "Up");
                                                              }, {"Muon_pt_corr", "Muon_eta", "Muon_mask"})
                             .Define("MuonIDScaleFactorDown", []( const rvec<float>& muon_pt,
                                                              const rvec<float>& muon_eta,
                                                              const rvec<bool>& muon_mask
                                                              )
                                                              {
                                                                 return muon_id_scale_factor(muon_pt, muon_eta, muon_mask, "Down");
                                                              }, {"Muon_pt_corr", "Muon_eta", "Muon_mask"})
                             .Define("WWshape", [](){ return (float) 1.;})
                             .Define("WWshapeUp", []( const float& coll_mass
                                                              )
                                                              {
                                                                 if (config::wwuncertainty){
                                                                     float weight = 1 + 1 - (0.993-0.00002001 * coll_mass + 0.00000000283 * coll_mass * coll_mass);
                                                                     return weight ;
                                                                 }
                                                                 return (float) 1.;
                                                              }, {"CollMass_alt"})
                             .Define("WWshapeDown", []( const float& coll_mass
                                                              )
                                                              {
                                                                 if (config::wwuncertainty){
                                                                     float weight = (0.993-0.00002001 * coll_mass + 0.00000000283 * coll_mass * coll_mass);
                                                                     return weight;
                                                                 }
                                                                 return (float) 1.;
                                                              }, {"CollMass_alt"})
                             .Define("TopQ", []( const float& Mll
                                                              )
                                                              {
                                                                 if (config::TT){
                                                                     return GetTopQscale(Mll,"nom");
                                                                 }
                                                                 return (float) 1.;
                                                              }, {"Mll"})
                             .Define("TopQUp", []( const float& Mll
                                                              )
                                                              {
                                                                 if (config::TT){
                                                                     return GetTopQscale(Mll,"up");
                                                                 }
                                                                 return (float) 1.;
                                                              }, {"Mll"})
                             .Define("TopQDown", []( const float& Mll
                                                              )
                                                              {
                                                                 if (config::TT){
                                                                     return GetTopQscale(Mll,"down");
                                                                 }
                                                                 return (float) 1.;
                                                              }, {"Mll"})
                             .Define("TopPDF", []( const float& Mll
                                                              )
                                                              {
                                                                 if (config::TT){
                                                                     return GetTopPDF(Mll,"nom");
                                                                 }
                                                                 return (float) 1.;
                                                              }, {"Mll"})
                             .Define("TopPDFUp", []( const float& Mll
                                                              )
                                                              {
                                                                 if (config::TT){
                                                                     return GetTopPDF(Mll,"up");
                                                                 }
                                                                 return (float) 1.;
                                                              }, {"Mll"})
                             .Define("TopPDFDown", []( const float& Mll
                                                              )
                                                              {
                                                                 if (config::TT){
                                                                     return GetTopPDF(Mll,"down");
                                                                 }
                                                                 return (float) 1.;
                                                              }, {"Mll"})
                        .Define("TriggerScaleFactor", [](const rvec<float>& pt,
                                                         const rvec<float>& eta,
                                                         const rvec<bool>& mask
                                                         ) { 
                                                             return trigger_sf(pt, eta,mask, "");
                                                        }, {"Muon_pt_corr", "Muon_eta", "Muon_mask"})
                        .Define("TriggerScaleFactorUp", [](const rvec<float>& pt,
                                                         const rvec<float>& eta,
                                                         const rvec<bool>& mask
                                                         ) { 
                                                             return trigger_sf(pt, eta, mask, "Up");
                                                        }, {"Muon_pt_corr", "Muon_eta", "Muon_mask"})
                        .Define("TriggerScaleFactorDown", [](const rvec<float>& pt,
                                                         const rvec<float>& eta,
                                                         const rvec<bool>& mask
                                                         ) { 
                                                             return trigger_sf(pt, eta, mask, "Down");
                                                        }, {"Muon_pt_corr", "Muon_eta", "Muon_mask"})
                        .Define("MuonReco", [](const rvec<float>& pt,
                                                         const rvec<float>& eta,
                                                         const rvec<float>& phi,
                                                         const rvec<float>& mass,
                                                         const rvec<bool>& mask
                                                         ) { 
                                                             return muon_reco_eff(pt, eta, phi, mass,mask,"");
                                                        }, {"Muon_pt_corr", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_mask"});
        //define_weight(total_weights_iso, default_weights_iso, list_of_weights_iso);

        //DEFINE NON-ENERGY RELATED UNCERTAINTIES
        total_weights_iso = total_weights_iso.Define("pileup_weightUp", [](const int& nvtx_true){ return pu_weight(nvtx_true, "Up"); }, {"Pileup_nPU"})
                                     .Define("pileup_weightDown", [](const int& nvtx_true){ return pu_weight(nvtx_true, "Down"); }, {"Pileup_nPU"})
                                     .Define("prefiring_weight_Up", [](const rvec<float>& jet_pt,
                                                             const rvec<float>& jet_eta,
                                                             const rvec<float>& photon_pt,
                                                             const rvec<float>& photon_eta) { 
                                                                 return prefire_factor(jet_pt, jet_eta, photon_pt, photon_eta, "Up");
                                                            }, {"Jet_pt", "Jet_eta", "Photon_pt", "Photon_eta"})
                                     .Define("prefiring_weight_Down", [](const rvec<float>& jet_pt,
                                                             const rvec<float>& jet_eta,
                                                             const rvec<float>& photon_pt,
                                                             const rvec<float>& photon_eta) { 
                                                                 return prefire_factor(jet_pt, jet_eta, photon_pt, photon_eta, "Down");
                                                        }, {"Jet_pt", "Jet_eta", "Photon_pt", "Photon_eta"});

    
        define_weight(total_weights_iso, default_weights_iso, list_of_weights_iso);
        define_weight(total_weights_iso, default_weights_iso, list_of_weights_iso, "pileup_weightUp", "pileup_weight");
        define_weight(total_weights_iso, default_weights_iso, list_of_weights_iso, "pileup_weightDown", "pileup_weight");
        define_weight(total_weights_iso, default_weights_iso, list_of_weights_iso, "TriggerScaleFactorUp", "TriggerScaleFactor");
        define_weight(total_weights_iso, default_weights_iso, list_of_weights_iso, "TriggerScaleFactorDown", "TriggerScaleFactor");
        define_weight(total_weights_iso, default_weights_iso, list_of_weights_iso, "MuonISOScaleFactorUp", "MuonISOScaleFactor");
        define_weight(total_weights_iso, default_weights_iso, list_of_weights_iso, "MuonISOScaleFactorDown", "MuonISOScaleFactor");
        define_weight(total_weights_iso, default_weights_iso, list_of_weights_iso, "MuonIDScaleFactorUp", "MuonIDScaleFactor");
        define_weight(total_weights_iso, default_weights_iso, list_of_weights_iso, "MuonIDScaleFactorDown", "MuonIDScaleFactor");
        define_weight(total_weights_iso, default_weights_iso, list_of_weights_iso, "WWshapeUp", "WWshape");
        define_weight(total_weights_iso, default_weights_iso, list_of_weights_iso, "WWshapeDown", "WWshape");
        define_weight(total_weights_iso, default_weights_iso, list_of_weights_iso, "TopQUp", "TopQ");
        define_weight(total_weights_iso, default_weights_iso, list_of_weights_iso, "TopQDown", "TopQ");
        define_weight(total_weights_iso, default_weights_iso, list_of_weights_iso, "TopPDFUp", "TopPDF");
        define_weight(total_weights_iso, default_weights_iso, list_of_weights_iso, "TopPDFDown", "TopPDF");

        
        // factorization & renormalization scale uncertainty
        std::vector<std::string> columnNames_iso = total_weights_iso.GetColumnNames();
        if (std::find(columnNames_iso.begin(), columnNames_iso.end(), "LHEScaleWeight") != columnNames_iso.end()) {
            total_weights_iso = total_weights_iso.Define("total_weight_ScaleWeightUp", [](const rvec<float>& weights, const float& tot_weight) {
                                                                                                if (weights.size() == 0)
                                                                                                    return (float) 1.0;
                                                                                                float max_val = *std::max_element(weights.begin(), weights.end());
                                                                                                return max_val*tot_weight;
                                                                                }, {"LHEScaleWeight", "total_weight"})
                                                                        .Define("total_weight_ScaleWeightDown", [](const rvec<float>& weights, const float& tot_weight) {
                                                                                                if (weights.size() == 0)
                                                                                                    return (float) 1.0;
                                                                                                float min_val = *std::min_element(weights.begin(), weights.end());
                                                                                                return min_val*tot_weight;
                                                                                }, {"LHEScaleWeight", "total_weight"});
            list_of_weights_iso.push_back("total_weight_ScaleWeightUp");
            list_of_weights_iso.push_back("total_weight_ScaleWeightDown");
        }
        
        // pdf weights
        std::vector<std::string> colNames_iso = total_weights_iso.GetColumnNames();
        if (std::find(colNames_iso.begin(), colNames_iso.end(), "LHEPdfWeight_def") != colNames_iso.end()) {
            for (uint i = 1; i < config::pdf_nweights; i++) {
                total_weights_iso = total_weights_iso.Define("total_weight_pdf_weight_" + std::to_string(i), [i](const float& weight, const rvec<float>& pdf_weights)
                    {
                        if (pdf_weights.size() != config::pdf_nweights) {
                            return weight;
                        }
                        return weight*pdf_weights[i]/pdf_weights[0];
                    }, {"total_weight", "LHEPdfWeight_def"});
                list_of_weights_iso.push_back("total_weight_pdf_weight_" + std::to_string(i));
            }
        }
    }

    //auto df_coll_iso = df_tau_pt_ov_jet_pt_iso.Define("CollMass", collinear_mass, {"col_idx","Tau_pt_ES", "Muon_pt", "Tau_eta_ES", "Muon_eta", "Tau_phi_ES", "Muon_phi", "Tau_mass_ES", "Muon_mass", "MET_pt", "MET_phi","Iso_Tau", "Muon_mask"});  
    //auto df_coll_alt_iso = df_coll_iso.Define("CollMass_alt", collinear_mass_alt, {"col_idx","Tau_pt_ES", "Muon_pt", "Tau_eta_ES", "Muon_eta", "Tau_phi_ES", "Muon_phi", "Tau_mass_ES", "Muon_mass", "MET_pt", "MET_phi","Iso_Tau", "Muon_mask"});  
    //// Define Stage 0: any event with one tau fulfilling acceptance and id	
	auto df_mt = total_weights.Filter("MT < 120.", "mt_cut");
	auto df_mt_iso = total_weights_iso.Filter("MT < 120.", "mt_cut");
	auto df_mt_signal = total_weights.Filter("MT > 120.", "mt_cut");
	auto df_mt_signal_iso = total_weights_iso.Filter("MT > 120.", "mt_cut");
    if (config::runOnData){
        create_datadriven_hists(df_mt, "NonIso_FakeRegion", muon_shift, tau_shift, met_shift, "total_weight");
        create_datadriven_hists(df_mt_iso, "Iso_FakeRegion", muon_shift, tau_shift, met_shift, "total_weight");
        create_datadriven_hists(df_mt_signal, "NonIso_SignalRegion", muon_shift, tau_shift, met_shift, "total_weight");
        create_datadriven_hists(df_mt_signal_iso, "Iso_SignalRegion", muon_shift, tau_shift, met_shift, "total_weight");
    }
    //fill_datadriven_with_syst(total_weights,"NonIso",muon_shift,tau_shift,met_shift,list_of_weights);
//	create_datadriven_hists(df_coll_iso, "Iso", "total_weight");
//	//
//	//// actual analysis cut  -- MT > 120. GeV 
//	auto df_mt = df_coll.Filter("MT < 120.", "mt_cut");
//	auto df_mt_iso = df_coll_iso.Filter("MT < 120.", "mt_cut");
//	auto df_mt_signal = df_coll.Filter("MT > 120.", "mt_cut");
//	auto df_mt_signal_iso = df_coll_iso.Filter("MT > 120.", "mt_cut");
//	create_datadriven_hists(df_mt, "NonIso_FakeRegion", "total_weight");
//	create_datadriven_hists(df_mt_iso, "Iso_FakeRegion", "total_weight");
//	create_datadriven_hists(df_mt_signal, "NonIso_SignalRegion", "total_weight");
//	create_datadriven_hists(df_mt_signal_iso, "Iso_SignalRegion", "total_weight");
    if (not config::runOnData){
        auto df_mc_tau_veto_noniso_closure = df_mt.Filter([](const rvec<int> col_idx, const rvec<UChar_t> genpartflav, const rvec<bool> tau_mask){
            if (genpartflav[tau_mask][col_idx[0]] != 0) return false;
            else return true;}, {"col_idx", "Tau_genPartFlav", "NonIso_Tau"},"Data Driven MC real Tau veto");
        auto df_mc_tau_veto_iso_closure = df_mt_iso.Filter([](const rvec<int> col_idx, const rvec<UChar_t> genpartflav, const rvec<bool> tau_mask){
            if (genpartflav[tau_mask][col_idx[0]] != 0) return false;
            else return true;}, {"col_idx", "Tau_genPartFlav", "Iso_Tau"},"Data Driven MC real Tau veto");
        auto df_mc_tau_veto_signal_noniso_closure = df_mt_signal.Filter([](const rvec<int> col_idx, const rvec<UChar_t> genpartflav, const rvec<bool> tau_mask){
            if (genpartflav[tau_mask][col_idx[0]] != 0) return false;
            else return true;}, {"col_idx", "Tau_genPartFlav", "NonIso_Tau"},"Data Driven MC real Tau veto");
        auto df_mc_tau_veto_signal_iso_closure = df_mt_signal_iso.Filter([](const rvec<int> col_idx, const rvec<UChar_t> genpartflav, const rvec<bool> tau_mask){
            if (genpartflav[tau_mask][col_idx[0]] != 0) return false;
            else return true;}, {"col_idx", "Tau_genPartFlav", "Iso_Tau"},"Data Driven MC real Tau veto");
	    create_datadriven_hists(df_mc_tau_veto_noniso_closure, "NonIso_FakeRegion_MC_Closure", muon_shift, tau_shift, met_shift,"total_weight");
	    create_datadriven_hists(df_mc_tau_veto_iso_closure, "Iso_FakeRegion_MC_Closure", muon_shift, tau_shift, met_shift,"total_weight");
	    create_datadriven_hists(df_mc_tau_veto_signal_noniso_closure, "NonIso_SignalRegion_MC_Closure", muon_shift, tau_shift, met_shift,"total_weight");
	    create_datadriven_hists(df_mc_tau_veto_signal_iso_closure, "Iso_SignalRegion_MC_Closure", muon_shift, tau_shift, met_shift,"total_weight");
        auto df_mc_tau_veto_noniso = df_mt.Filter([](const rvec<int> col_idx, const rvec<UChar_t> genpartflav, const rvec<bool> tau_mask){
            if (genpartflav[tau_mask][col_idx[0]] == 0) return false;
            else return true;}, {"col_idx", "Tau_genPartFlav", "NonIso_Tau"},"Data Driven MC real Tau veto");
        auto df_mc_tau_veto_iso = df_mt_iso.Filter([](const rvec<int> col_idx, const rvec<UChar_t> genpartflav, const rvec<bool> tau_mask){
            if (genpartflav[tau_mask][col_idx[0]] == 0) return false;
            else return true;}, {"col_idx", "Tau_genPartFlav", "Iso_Tau"},"Data Driven MC real Tau veto");
        auto df_mc_tau_veto_signal_noniso = df_mt_signal.Filter([](const rvec<int> col_idx, const rvec<UChar_t> genpartflav, const rvec<bool> tau_mask){
            if (genpartflav[tau_mask][col_idx[0]] == 0) return false;
            else return true;}, {"col_idx", "Tau_genPartFlav", "NonIso_Tau"},"Data Driven MC real Tau veto");
        auto df_mc_tau_veto_signal_iso = df_mt_signal_iso.Filter([](const rvec<int> col_idx, const rvec<UChar_t> genpartflav, const rvec<bool> tau_mask){
            if (genpartflav[tau_mask][col_idx[0]] == 0) return false;
            else return true;}, {"col_idx", "Tau_genPartFlav", "Iso_Tau"},"Data Driven MC real Tau veto");
	    create_datadriven_hists(df_mc_tau_veto_noniso, "NonIso_FakeRegion_MC", muon_shift, tau_shift, met_shift,"total_weight");
	    create_datadriven_hists(df_mc_tau_veto_iso, "Iso_FakeRegion_MC", muon_shift, tau_shift, met_shift,"total_weight");
	    create_datadriven_hists(df_mc_tau_veto_signal_noniso, "NonIso_SignalRegion_MC", muon_shift, tau_shift, met_shift,"total_weight");
	    create_datadriven_hists(df_mc_tau_veto_signal_iso, "Iso_SignalRegion_MC", muon_shift, tau_shift, met_shift,"total_weight");
//
    }
    return df_coll;
    }
//
//    for (const auto& it : list_of_weights) {
//        std::string folder_name = "Stage1";
//        if ( it.find("total_weight_") != std::string::npos )
//            folder_name = "Stage1/Systematics";
//        if ( it.find("pdf_weight") != std::string::npos )
//            folder_name = "Stage1/PDF";
//        create_datadriven_hists(df_mt, folder_name, it);
//    }
//    return df_mt;
//};

// Apply datadriven distributions
RNode apply_datadriven(RNode df, std::string muon_shift, std::string tau_shift, std::string met_shift) {
//	
	auto noniso_tau = df.Define("NonIso_Tau", tau_acceptance_and_id_and_dm_noniso,	{"Tau_pt" + tau_shift, "Tau_eta" + tau_shift, "Tau_dz","Tau_charge", config::tau_dm, "Tau_decayMode", "Tau_jetIdx", config::tau_iso, config::tau_antiEle, config::tau_antiMuon});
	//										
	auto iso_tau = df.Define("Iso_Tau", tau_acceptance_and_id_and_dm,	{"Tau_pt" + tau_shift, "Tau_eta" + tau_shift, "Tau_dz", "Tau_charge", config::tau_dm, "Tau_decayMode", "Tau_jetIdx", config::tau_iso, config::tau_antiEle, config::tau_antiMuon});
	//										
	auto evntselect_noniso_tau = noniso_tau.Filter(nparticle_cut, {"NonIso_Tau"}, "select_taus")
                                .Filter(nparticle_cut, {"Muon_mask"}, "select_muons")
                                .Filter(nparticle_veto, {"Electron_mask"}, "select_eles");
	auto evntselect_iso_tau = iso_tau.Filter(nparticle_cut, {"Iso_Tau"}, "select_taus")
                                .Filter(nparticle_cut, {"Muon_mask"}, "select_muons")
                                .Filter(nparticle_veto, {"Electron_mask"}, "select_eles");
//	
    auto df_col_idx_noniso = evntselect_noniso_tau.Define("col_idx", col_idx, {"Tau_pt" + tau_shift, "Muon_pt_corr" + muon_shift, "Tau_eta" + tau_shift, "Muon_eta", "Tau_phi" + tau_shift, "Muon_phi", "Tau_mass" + tau_shift, "Muon_mass", met_branch_name + "_pt" + met_shift, met_branch_name + "_phi" + met_shift , "NonIso_Tau", "Muon_mask"});   
        auto defquants_noniso = df_col_idx_noniso.Define("sel_Tau_pt", selected_part_col_idx_tau, {"Tau_pt" + tau_shift, "NonIso_Tau", "col_idx"})
                                   .Define("sel_Tau_pt_barrel", selected_part_col_idx_tau_barrel, {"Tau_pt" + tau_shift, "NonIso_Tau", "col_idx", "Tau_eta" + tau_shift})
                                   .Define("sel_Tau_pt_endcap", selected_part_col_idx_tau_endcap, {"Tau_pt" + tau_shift, "NonIso_Tau", "col_idx", "Tau_eta" + tau_shift})
                                   .Define("sel_Tau_eta", selected_part_col_idx_tau, {"Tau_eta" + tau_shift, "NonIso_Tau", "col_idx"})
                                   .Define("sel_Tau_phi", selected_part_col_idx_tau, {"Tau_phi" + tau_shift, "NonIso_Tau", "col_idx"})
                                   .Define("sel_Tau_mass", selected_part_col_idx_tau, {"Tau_mass" + tau_shift, "NonIso_Tau", "col_idx"})
                                   .Define("sel_Tau_jetIdx", selected_jet_quant, {"Tau_jetIdx", "NonIso_Tau", "col_idx"})
                                   .Define("sel_Tau_pt_norm", selected_part_col_idx_tau, {"Tau_pt" , "Tau_mask", "col_idx"})
                                   .Define("sel_Tau_eta_norm", selected_part_col_idx_tau, {"Tau_eta" , "Tau_mask", "col_idx"})
                                   .Define("sel_Tau_phi_norm", selected_part_col_idx_tau, {"Tau_phi" , "Tau_mask", "col_idx"})
                                   .Define("sel_Tau_mass_norm", selected_part_col_idx_tau, {"Tau_mass", "Tau_mask", "col_idx"})
                                   .Define("sel_Muon_PF_pt", selected_part_col_idx_muon, {"Muon_pt", "Muon_mask", "col_idx"})
                                   .Define("sel_Muon_isPF", selected_part_col_idx_muon_bool, {"Muon_isPFcand", "Muon_mask", "col_idx"})
                                   .Define("sel_Muon_pt", selected_part_col_idx_muon, {"Muon_pt_corr" + muon_shift, "Muon_mask", "col_idx"})
                                   .Define("sel_Muon_eta", selected_part_col_idx_muon, {"Muon_eta", "Muon_mask", "col_idx"})
                                   .Define("sel_Muon_phi", selected_part_col_idx_muon, {"Muon_phi", "Muon_mask", "col_idx"})
                                   .Define("sel_Muon_mass", selected_part_col_idx_muon, {"Muon_mass", "Muon_mask", "col_idx"})
                                   .Define("sel_MET_pt" + met_shift, correct_met_single,{"sel_Muon_PF_pt","sel_Muon_pt","sel_Muon_eta","sel_Muon_phi","sel_Muon_mass","sel_Muon_isPF","sel_Tau_pt","sel_Tau_eta","sel_Tau_phi","sel_Tau_mass","sel_Tau_pt_norm","sel_Tau_eta_norm","sel_Tau_phi_norm","sel_Tau_mass_norm",met_branch_name + "_pt" + met_shift,met_branch_name + "_phi" + met_shift,"run","PV_npvsGood"})
                                   .Define("sel_MET_phi" + met_shift, correct_met_single_phi,{"sel_Muon_PF_pt","sel_Muon_pt","sel_Muon_eta","sel_Muon_phi","sel_Muon_mass","sel_Muon_isPF","sel_Tau_pt","sel_Tau_eta","sel_Tau_phi","sel_Tau_mass","sel_Tau_pt_norm","sel_Tau_eta_norm","sel_Tau_phi_norm","sel_Tau_mass_norm",met_branch_name + "_pt" + met_shift,met_branch_name + "_phi" + met_shift,"run","PV_npvsGood"});
    if (!config::runOnData) {
        defquants_noniso = defquants_noniso.Define("sel_Tau_genPartFlav", selected_part_col_idx_tau_flav, {"Tau_genPartFlav", "NonIso_Tau", "col_idx"});
    }
    auto dimuon_cut_noniso = defquants_noniso.Filter(dimuonpair_cut,{"Muon_pt_corr" + muon_shift,"Muon_eta","Muon_phi", "DiMuon_mask"},"Di Muon Pair Cut");
    auto delR_cut_noniso = dimuon_cut_noniso.Filter(mutaudelR_cut,{"sel_Muon_eta","sel_Muon_phi", "sel_Tau_eta", "sel_Tau_phi"},"Mu Tau Del R Cut");
    auto df_jet_matching_noniso = delR_cut_noniso.Define("matched_jet",jet_matching, {"sel_Tau_pt","sel_Tau_eta","sel_Tau_phi","sel_Tau_mass","Jet_pt","Jet_eta","Jet_phi","Jet_mass"});
    auto jet_matching_cut_noniso = df_jet_matching_noniso.Filter("matched_jet == true","jet matching cut");


    auto corr_met_df_noniso  = jet_matching_cut_noniso;
	auto mtcalc_noniso = corr_met_df_noniso.Define("MT", mass_transv, {"sel_Muon_pt", "sel_Muon_phi", "sel_MET_pt" + met_shift, "sel_MET_phi" + met_shift});//.Define("sel_Tau_pt_ov_Jet_pt_Iso", [](const rvec<float> tau_pt, 
    auto df_coll_noniso = mtcalc_noniso.Define("CollMass", collinear_mass, {"sel_Tau_pt", "sel_Muon_pt", "sel_Tau_eta", "sel_Muon_eta", "sel_Tau_phi", "sel_Muon_phi", "sel_Tau_mass", "sel_Muon_mass", "sel_MET_pt" + met_shift, "sel_MET_phi" + met_shift}).Define("CollMass_alt", collinear_mass_alt, {"sel_Tau_pt", "sel_Muon_pt", "sel_Tau_eta", "sel_Muon_eta", "sel_Tau_phi", "sel_Muon_phi", "sel_Tau_mass", "sel_Muon_mass", "sel_MET_pt" + met_shift, "sel_MET_phi" + met_shift}).Define("Tau_pt_over_Jet_pt", tau_pt_over_jet_pt,{"col_idx","Tau_pt" + tau_shift,"Tau_jetIdx","Jet_pt_nom","NonIso_Tau"}).Define("Tau_pt_over_Jet_pt_barrel", tau_pt_over_jet_pt_barrel,{"col_idx","Tau_pt" + tau_shift, "Tau_jetIdx","Jet_pt_nom","NonIso_Tau","Tau_eta" + tau_shift}).Define("Tau_pt_over_Jet_pt_endcap", tau_pt_over_jet_pt_endcap,{"col_idx","Tau_pt" + tau_shift, "Tau_jetIdx", "Jet_pt_nom", "NonIso_Tau", "Tau_eta" + tau_shift});  
	RNode total_weights_noniso = df_coll_noniso;
    std::vector < std::string > list_of_weights_noniso;
    if (config::runOnData) {
        total_weights_noniso  = df_coll_noniso.Define("total_weight", 
                                                    []( const rvec<float>& tau_pt,
                                                        const rvec<float>& tau_eta,
                                                           const float & taupt_o_jetpt,
                                                           const rvec<int>& col_idx, 
                                                           const rvec<bool>& tau_mask
                                                        )
                                                    {
                                                        return dd_fakerate(tau_pt, tau_eta, taupt_o_jetpt, col_idx, tau_mask, "");
                                                    }, {"Tau_pt" + tau_shift,"Tau_eta" + tau_shift,"Tau_pt_over_Jet_pt","col_idx","NonIso_Tau"});
        list_of_weights_noniso.push_back("total_weight");
    } else {    
        std::vector < std::string > default_weights_noniso = {
            "pileup_weight",
            "genWeight",
            "PrefiringWeight",
            "TriggerScaleFactor",
            "MuonISOScaleFactor",
            "MuonIDScaleFactor",
            "WWshape",
            "TopQ",
            "TopPDF",
            "top_pt_weight",
            "MuonReco",
            "FakeRate"
        };
        total_weights_noniso = df_coll_noniso.Define("FakeRate", 
                                                    []( const rvec<float>& tau_pt,
                                                        const rvec<float>& tau_eta,
                                                           const float & taupt_o_jetpt,
                                                           const rvec<int>& col_idx, 
                                                           const rvec<bool>& tau_mask
                                                        )
                                                    {
                                                        return dd_fakerate(tau_pt, tau_eta, taupt_o_jetpt, col_idx, tau_mask, "");
                                                    }, {"Tau_pt" + tau_shift, "Tau_eta" + tau_shift, "Tau_pt_over_Jet_pt","col_idx","NonIso_Tau"}).Define("MuonISOScaleFactor", [](const rvec<float>& muon_pt,
                                                              const rvec<float>& muon_eta,
                                                              const rvec<bool>& muon_mask
                                                              )
                                                              {
                                                                 return muon_iso_scale_factor(muon_pt, muon_eta, muon_mask, "");
                                                              }, {"Muon_pt_corr", "Muon_eta", "Muon_mask"})
                             .Define("MuonIDScaleFactor", []( const rvec<float>& muon_pt,
                                                              const rvec<float>& muon_eta,
                                                              const rvec<bool>& muon_mask
                                                              )
                                                              {
                                                                 return muon_id_scale_factor(muon_pt, muon_eta, muon_mask, "");
                                                              }, {"Muon_pt_corr", "Muon_eta", "Muon_mask"})
                             .Define("WWshape", [](){ return (float) 1.;})
                             .Define("TopQ", []( const float& Mll
                                                              )
                                                              {
                                                                 if (config::TT){
                                                                     return GetTopQscale(Mll,"nom");
                                                                 }
                                                                 return (float) 1.;
                                                              }, {"Mll"})
                             .Define("TopPDF", []( const float& Mll
                                                              )
                                                              {
                                                                 if (config::TT){
                                                                     return GetTopPDF(Mll,"nom");
                                                                 }
                                                                 return (float) 1.;
                                                              }, {"Mll"})
                        .Define("TriggerScaleFactor", [](const rvec<float>& pt,
                                                         const rvec<float>& eta,
                                                         const rvec<bool>& mask
                                                         ) { 
                                                             return trigger_sf(pt, eta,mask, "");
                                                        }, {"Muon_pt_corr", "Muon_eta", "Muon_mask"})
                        .Define("MuonReco", [](const rvec<float>& pt,
                                                         const rvec<float>& eta,
                                                         const rvec<float>& phi,
                                                         const rvec<float>& mass,
                                                         const rvec<bool>& mask
                                                         ) { 
                                                             return muon_reco_eff(pt, eta, phi, mass,mask,"");
                                                        }, {"Muon_pt_corr", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_mask"});
        define_weight(total_weights_noniso, default_weights_noniso, list_of_weights_noniso);
        }
	auto df_mt_noniso = total_weights_noniso.Filter("MT < 120.", "mt_cut").Define("Sphericity",sphericity_full,{"Tau_pt","Tau_eta","Tau_phi","Tau_mass","Muon_pt_corr" + muon_shift,"Muon_eta","Muon_phi","Muon_mass","Jet_pt_nom","Jet_eta", "Jet_phi", "Jet_mass","NonIso_Tau", "Muon_mask"}).Define("Sphericity_leptons",sphericity_leptons,{"Tau_pt","Tau_eta","Tau_phi","Tau_mass","Muon_pt_corr" + muon_shift,"Muon_eta","Muon_phi","Muon_mass","NonIso_Tau", "Muon_mask"});
	auto df_mt_signal_noniso = total_weights_noniso.Filter("MT > 120.", "mt_cut").Define("Sphericity",sphericity_full,{"Tau_pt","Tau_eta","Tau_phi","Tau_mass","Muon_pt_corr" + muon_shift,"Muon_eta","Muon_phi","Muon_mass","Jet_pt_nom","Jet_eta", "Jet_phi", "Jet_mass","NonIso_Tau", "Muon_mask"}).Define("Sphericity_leptons",sphericity_leptons,{"Tau_pt","Tau_eta","Tau_phi","Tau_mass","Muon_pt_corr" + muon_shift,"Muon_eta","Muon_phi","Muon_mass","NonIso_Tau", "Muon_mask"});
    if (config::runOnData){
        create_datadriven_hists(df_mt_signal_noniso, "Fakes_SignalRegion", muon_shift, tau_shift, met_shift, "total_weight");
    }
    if (not config::runOnData){
        auto df_mc_tau_veto_signal_noniso = df_mt_signal_noniso.Filter([](const rvec<int> col_idx, const rvec<UChar_t> genpartflav, const rvec<bool> tau_mask){
            if (genpartflav[tau_mask][col_idx[0]] != 0) return false;
            else return true;}, {"col_idx", "Tau_genPartFlav", "NonIso_Tau"},"Data Driven MC real Tau veto");
        auto df_mc_no_fake_signal_noniso = df_mt_signal_noniso.Filter([](const rvec<int> col_idx, const rvec<UChar_t> genpartflav, const rvec<bool> tau_mask){
            if (genpartflav[tau_mask][col_idx[0]] == 0) return false;
            else return true;}, {"col_idx", "Tau_genPartFlav", "NonIso_Tau"},"Data Driven MC real Tau veto");
	    //create_datadriven_hists(df_mc_tau_veto_signal_noniso, "Fakes_SignalRegion_MC", muon_shift, tau_shift, met_shift,"total_weight");
	    //create_datadriven_hists(df_mt_signal_noniso, "Fakes_SignalRegion_MC", muon_shift, tau_shift, met_shift,"total_weight");
	    create_datadriven_hists(df_mc_tau_veto_signal_noniso, "Fakes_SignalRegion_MC", muon_shift, tau_shift, met_shift,"total_weight");
	    create_datadriven_hists(df_mc_no_fake_signal_noniso, "No_Fakes_SignalRegion_MC", muon_shift, tau_shift, met_shift,"total_weight");
//
    }
    return df_mt_noniso;
    }
//    auto define_weights_iso = df_coll_alt_iso;
//    if (config::runOnData){
//        define_weights_iso  = df_coll_alt_iso.Define("FakeRate", 
//                                                    []( const rvec<float>& tau_pt,
//                                                           const float & taupt_o_jetpt,
//                                                           const rvec<int>& col_idx, 
//                                                           const rvec<bool>& tau_mask
//                                                        )
//                                                    {
//                                                        return dd_fakerate(tau_pt, taupt_o_jetpt, col_idx, tau_mask, "");
//                                                    }, {"Tau_pt_ES","Tau_pt_over_Jet_pt","col_idx","Iso_Tau"});
//    }
//    else
//    {
//        define_weights_iso = df_coll_alt_iso.Define("FakeRate", 
//                                                    []( const rvec<float>& tau_pt,
//                                                           const float & taupt_o_jetpt,
//                                                           const rvec<int>& col_idx, 
//                                                           const rvec<bool>& tau_mask
//                                                        )
//                                                    {
//                                                        return dd_fakerate(tau_pt, taupt_o_jetpt, col_idx, tau_mask, "");
//                                                    }, {"Tau_pt_ES","Tau_pt_over_Jet_pt","col_idx","Iso_Tau"})
//                                                    .Define(   "TauEleFakeScaleFactor", 
//                                                    []( const rvec<float>& tau_pt,
//                                                        const rvec<float>& tau_eta,
//                                                        const rvec<bool>& tau_mask,
//                                                        const rvec<UChar_t>& tau_genPartFlav
//                                                        )
//                                                    {
//                                                        return tau_fake_scale_factor(tau_pt, tau_eta, tau_mask, tau_genPartFlav, "Ele", "");
//                                                    }, {"Tau_pt", "Tau_eta", "Iso_Tau", "Tau_genPartFlav"})
//                                        
//                                        .Define(    "TauJetFakeScaleFactor", 
//                                                    []( const rvec<float>& tau_pt,
//                                                        const rvec<float>& tau_eta,
//                                                        const rvec<bool>& tau_mask,
//                                                        const rvec<UChar_t>& tau_genPartFlav
//                                                        )
//                                                    {
//                                                        return tau_fake_scale_factor(tau_pt, tau_eta, tau_mask, tau_genPartFlav, "Jet", "");
//                                                    }, {"Tau_pt", "Tau_eta", "Iso_Tau", "Tau_genPartFlav"})
//                                        
//                                        .Define(    "TauMuonFakeScaleFactor", 
//                                                    []( const rvec<float>& tau_pt,
//                                                        const rvec<float>& tau_eta,
//                                                        const rvec<bool>& tau_mask,
//                                                        const rvec<UChar_t>& tau_genPartFlav
//                                                        )
//                                                    {
//                                                        return tau_fake_scale_factor(tau_pt, tau_eta, tau_mask, tau_genPartFlav, "Muon", "");
//                                                    }, {"Tau_pt", "Tau_eta", "Iso_Tau", "Tau_genPartFlav"})
//                                        
//                                        .Define("PrefiringWeight", prefire_factor, {"Jet_pt", "Jet_eta", "Photon_pt", "Photon_eta"})
//                                        .Define("MuonISOScaleFactor", [](const rvec<float>& muon_pt,
//                                                                         const rvec<float>& muon_eta,
//                                                                         const rvec<bool>& muon_mask
//                                                                         )
//                                                                         {
//                                                                            return muon_iso_scale_factor(muon_pt, muon_eta, muon_mask, "");
//                                                                         }, {"Muon_pt", "Muon_eta", "Muon_mask"})
//                                        .Define("MuonIDScaleFactor", []( const rvec<float>& muon_pt,
//                                                                         const rvec<float>& muon_eta,
//                                                                         const rvec<bool>& muon_mask
//                                                                         )
//                                                                         {
//                                                                            return muon_id_scale_factor(muon_pt, muon_eta, muon_mask, "");
//                                                                         }, {"Muon_pt", "Muon_eta", "Muon_mask"});
//        
//    }
//    //                                
//    //
//    // Define weights
//	RNode total_weights_iso = define_weights_iso; 
//    if (config::runOnData) {
//        total_weights_iso = define_weights_iso.Define("total_weight", "FakeRate");
//    } else {
//        auto def_weights = [](const rvec<float>& weights){
//            float mult = 1.0;
//            for (const auto& it : weights)
//                mult *= it;
//            return mult;
//        };
//        total_weights_iso = define_weights_iso.Define("total_weight", ROOT::RDF::PassAsVec<9, float>(def_weights), {
//                                                                                                "FakeRate",
//                                                                                                "pileup_weight",
//                                                                                                "genWeight",
//                                                                                                "TauJetFakeScaleFactor",
//                                                                                                "TauEleFakeScaleFactor",
//                                                                                                "TauMuonFakeScaleFactor",
//                                                                                                "MuonIDScaleFactor",
//                                                                                                "MuonISOScaleFactor",
//                                                                                                "PrefiringWeight"});
//    }
//    // Define Stage 0: any event with one tau fulfilling acceptance and id	
//	create_fr_hists(total_weights_iso, "Iso", "total_weight");
//	
//	// actual analysis cut  -- MT > 120. GeV 
//	auto df_mt = total_weights_iso.Filter("MT < 120.", "mt_cut");
//	//auto df_mt_iso = total_weights_iso.Filter("MT < 120.", "mt_cut");
//	//auto df_mt_signal = total_weights_iso.Filter("MT > 120.", "mt_cut");
//	//auto df_mt_signal_iso = total_weights_iso.Filter("MT > 120.", "mt_cut");
//	//create_fr_hists(df_mt_iso, "Iso_FakeRegion", "FakeRate");
//	//create_fr_hists(df_mt_signal_iso, "Iso_SignalRegion", "FakeRate");
//    return df_mt;
////
////    for (const auto& it : list_of_weights) {
////        std::string folder_name = "Stage1";
////        if ( it.find("total_weight_") != std::string::npos )
////            folder_name = "Stage1/Systematics";
////        if ( it.find("pdf_weight") != std::string::npos )
////            folder_name = "Stage1/PDF";
////        create_datadriven_hists(df_mt, folder_name, it);
////    }
////    return df_mt;

// analyse function - gets called for each systematic
RNode run_analyser(  RNode df,
                    std::string muon_shift,
                    std::string tau_shift,
                    std::string met_shift,
                    std::vector< ROOT::RDF::RResultPtr< double >  >& genWeightSums,
                    TFile* outFile) {
    
    // should only happen on default run (non scale systs and others)
    bool default_run = (tau_shift == "_ES" && (met_shift == "_jer" || met_shift == "_nom") && muon_shift == "");

   // if (met_shift == "_nom"){
   //     met_shift = "";
   // }

	// init counter
//	auto genWeightSum = df.Sum("genWeight");
//  auto eventcounter = df.Count();
    
	// Select trigger requirements
								
	// Trigger turn on cut
	//auto trigger_obj1 = met_filter.Filter("Muon_pt > " + std::to_string(config::muon_pt), "Avoid trigger turn on with muon object, pT>" + std::to_string(config::muon_pt));
	

    auto preselection = df;
    if (config::runOnData) {
        preselection = preselection.Define("presel_weight", "1.0");
    } else {
        // calculate pileup weight
        preselection = preselection.Define("pileup_weight", [](const int& nvtx_true){ return pu_weight(nvtx_true, ""); }, {"Pileup_nPU"})
                                   .Define("PrefiringWeight", [](const rvec<float>& jet_pt,
                                                                     const rvec<float>& jet_eta,
                                                                     const rvec<float>& photon_pt,
                                                                     const rvec<float>& photon_eta) { 
                                                                         return prefire_factor(jet_pt, jet_eta, photon_pt, photon_eta, "");
                                                                    }, {"Jet_pt", "Jet_eta", "Photon_pt", "Photon_eta"});
        preselection = preselection.Define("W_kfactor", get_kfactor, {"GenPart_pdgId", "GenPart_mass"});
        preselection = preselection.Define("presel_weight", "pileup_weight*PrefiringWeight*W_kfactor*genWeight");
    }
	auto triggered = trigger(preselection,default_run);
	
	
	// creates mask, which fills bool tags for taus which fulfil id and are in acceptance
	auto masked = triggered.Define("Tau_mask", tau_acceptance_and_id_and_dm, {"Tau_pt" + tau_shift, "Tau_eta" + tau_shift, "Tau_dz", "Tau_charge", config::tau_dm, "Tau_decayMode", "Tau_jetIdx", config::tau_iso, config::tau_antiEle, config::tau_antiMuon})
							  .Define("Muon_mask", muon_acceptance_and_id, {"Muon_pt_corr" + muon_shift, "Muon_eta", "Muon_isPFcand", "Muon_highPtId", "Muon_tkIsoId"})
							  .Define("DiMuon_mask", di_muon_id, {"Muon_pt_corr" + muon_shift, "Muon_eta", "Muon_highPtId", "Muon_tkRelIso"})
							  .Define("Electron_mask", ele_acceptance_and_simpleid, {"Electron_pt", "Electron_eta", "Electron_cutBased_HEEP"});
	
	
    //auto df_report = masked.Report();
    //df_report->Print();
	
	// this is for datadriven part of analysis
    auto df_final = [&]() -> RNode{
        if (config::calcDataDriven){
            if(default_run){
            return calc_datadriven(masked,
                     muon_shift,
                     tau_shift,
                     met_shift);
            }
        }
        if (config::runDataDriven ){
            return apply_datadriven(masked,
                                    muon_shift,
                                    tau_shift,
                                    met_shift);
        }
        else{	
        
        // Select events with certain number of taus, which fulfil all acceptance & id
        // also, cut any event with electrons or muons
        auto evntselect = masked.Filter(nparticle_cut, {"Tau_mask"}, "select_taus")
                                .Filter(nparticle_cut, {"Muon_mask"}, "select_muons")
                                .Filter(nparticle_veto, {"Electron_mask"}, "select_eles");
        
        
        auto df_col_idx = evntselect.Define("col_idx", col_idx, {"Tau_pt" + tau_shift, "Muon_pt_corr" + muon_shift, "Tau_eta" + tau_shift, "Muon_eta", "Tau_phi" + tau_shift, "Muon_phi", "Tau_mass" + tau_shift, "Muon_mass", met_branch_name + "_pt" + met_shift, met_branch_name + "_phi" + met_shift , "Tau_mask", "Muon_mask"});   
        // Now that we have our real tau, write it to special columns to handle it more easily
        auto defquants = df_col_idx.Define("sel_Tau_pt", selected_part_col_idx_tau, {"Tau_pt" + tau_shift, "Tau_mask", "col_idx"})
                                   .Define("sel_Tau_pt_barrel", selected_part_col_idx_tau_barrel, {"Tau_pt" + tau_shift, "Tau_mask", "col_idx", "Tau_eta" + tau_shift})
                                   .Define("sel_Tau_pt_endcap", selected_part_col_idx_tau_endcap, {"Tau_pt" + tau_shift, "Tau_mask", "col_idx", "Tau_eta" + tau_shift})
                                   .Define("sel_Tau_eta", selected_part_col_idx_tau, {"Tau_eta" + tau_shift, "Tau_mask", "col_idx"})
                                   .Define("sel_Tau_phi", selected_part_col_idx_tau, {"Tau_phi" + tau_shift, "Tau_mask", "col_idx"})
                                   .Define("sel_Tau_mass", selected_part_col_idx_tau, {"Tau_mass" + tau_shift, "Tau_mask", "col_idx"})
                                   .Define("sel_Tau_pt_norm", selected_part_col_idx_tau, {"Tau_pt" , "Tau_mask", "col_idx"})
                                   .Define("sel_Tau_eta_norm", selected_part_col_idx_tau, {"Tau_eta" , "Tau_mask", "col_idx"})
                                   .Define("sel_Tau_phi_norm", selected_part_col_idx_tau, {"Tau_phi" , "Tau_mask", "col_idx"})
                                   .Define("sel_Tau_mass_norm", selected_part_col_idx_tau, {"Tau_mass", "Tau_mask", "col_idx"})
                                   .Define("sel_Muon_PF_pt", selected_part_col_idx_muon, {"Muon_pt", "Muon_mask", "col_idx"})
                                   .Define("sel_Muon_isPF", selected_part_col_idx_muon_bool, {"Muon_isPFcand", "Muon_mask", "col_idx"})
                                   .Define("sel_Muon_pt", selected_part_col_idx_muon, {"Muon_pt_corr" + muon_shift, "Muon_mask", "col_idx"})
                                   .Define("sel_Muon_eta", selected_part_col_idx_muon, {"Muon_eta", "Muon_mask", "col_idx"})
                                   .Define("sel_Muon_phi", selected_part_col_idx_muon, {"Muon_phi", "Muon_mask", "col_idx"})
                                   .Define("sel_Muon_mass", selected_part_col_idx_muon, {"Muon_mass", "Muon_mask", "col_idx"})
                                   .Define("sel_MET_pt" + met_shift, correct_met_single,{"sel_Muon_PF_pt","sel_Muon_pt","sel_Muon_eta","sel_Muon_phi","sel_Muon_mass","sel_Muon_isPF","sel_Tau_pt","sel_Tau_eta","sel_Tau_phi","sel_Tau_mass","sel_Tau_pt_norm","sel_Tau_eta_norm","sel_Tau_phi_norm","sel_Tau_mass_norm",met_branch_name + "_pt" + met_shift,met_branch_name + "_phi" + met_shift,"run","PV_npvsGood"})
                                   .Define("sel_MET_phi" + met_shift, correct_met_single_phi,{"sel_Muon_PF_pt","sel_Muon_pt","sel_Muon_eta","sel_Muon_phi","sel_Muon_mass","sel_Muon_isPF","sel_Tau_pt","sel_Tau_eta","sel_Tau_phi","sel_Tau_mass","sel_Tau_pt_norm","sel_Tau_eta_norm","sel_Tau_phi_norm","sel_Tau_mass_norm",met_branch_name + "_pt" + met_shift,met_branch_name + "_phi" + met_shift,"run","PV_npvsGood"})
                                   .Define("sel_MET_pt_full" + met_shift, correct_met,{"Muon_pt","Muon_pt_corr","Muon_eta","Muon_phi","Muon_mass","Muon_highPtId","Muon_isPFcand",met_branch_name + "_pt",met_branch_name + "_phi","run","PV_npvs"})
                                   .Define("sel_MET_phi_full" + met_shift, correct_met_phi,{"Muon_pt","Muon_pt_corr","Muon_eta","Muon_phi","Muon_mass","Muon_highPtId","Muon_isPFcand",met_branch_name + "_pt_nom",met_branch_name + "_phi","run","PV_npvs"});
                if (!config::runOnData) {
            defquants = defquants.Define("sel_Tau_genPartFlav", selected_part_col_idx_tau_flav, {"Tau_genPartFlav", "Tau_mask", "col_idx"})
                                 .Define("sel_Tau_genIdx", selected_part_col_idx_tau_idx,{"Tau_genPartIdx","Tau_mask","col_idx"})
                                 .Define("sel_Tau_genVisIdx", gen_match_idx,{"GenVisTau_pt","GenVisTau_eta","GenVisTau_phi","GenVisTau_mass", "sel_Tau_pt","sel_Tau_eta","sel_Tau_phi", "sel_Tau_mass", })
                                 .Define("sel_Muon_genIdx", selected_part_col_idx_muon_idx,{"Muon_genPartIdx","Muon_mask","col_idx"});
            defquants = defquants.Define("CollMass_gen", [](const rvec<float>& genpt,
                                                            const rvec<float>& geneta, 
                                                            const rvec<float>& genphi, 
                                                            const rvec<float>& genmass, 
                                                            const float& genPartFlav,
                                                            const rvec<float>& genmuonpt,
                                                            const rvec<float>& genmuoneta,
                                                            const rvec<float>& genmuonphi,
                                                            const rvec<float>& genmuonmass,
                                                            const int& Tau_genIdx,
                                                            const int& Tau_genVisIdx,
                                                            const int& Muon_genIdx,
                                                            const float& met_pt,
                                                            const float& met_phi){
                                            TLorentzVector tau, muon, met;
                                            if(genPartFlav == 5){
                                                tau.SetPtEtaPhiM(genpt[Tau_genVisIdx],geneta[Tau_genVisIdx],genphi[Tau_genVisIdx],genmass[Tau_genVisIdx]);
                                                muon.SetPtEtaPhiM(genmuonpt[Muon_genIdx],genmuoneta[Muon_genIdx],genmuonphi[Muon_genIdx],genmuonmass[Muon_genIdx]);
                                                met.SetPtEtaPhiM(met_pt,0.,met_phi,0.);
                                                double METproj=(met.Px()*tau.Px()+met.Py()*tau.Py())/tau.Pt();
                                                double xth = 1;
                                                if(METproj>0) xth=tau.Pt()/(tau.Pt() + METproj);
                                                else xth = 1;
                                                double mass_vis = (tau + muon).M();
                                                double mcol = 0;
                                                if (mass_vis != mass_vis) mass_vis=0;
                                                if (mass_vis <= 0.) mass_vis = 0.;
                                                mcol = mass_vis/sqrt(xth);
                                                return mcol;
                                            }
                                            else{
                                                tau.SetPtEtaPhiM(genmuonpt[Tau_genIdx],genmuoneta[Tau_genIdx],genmuonphi[Tau_genIdx],genmuonmass[Tau_genIdx]);
                                                muon.SetPtEtaPhiM(genmuonpt[Muon_genIdx],genmuoneta[Muon_genIdx],genmuonphi[Muon_genIdx],genmuonmass[Muon_genIdx]);
                                                met.SetPtEtaPhiM(met_pt,0.,met_phi,0.);
                                                double METproj=(met.Px()*tau.Px()+met.Py()*tau.Py())/tau.Pt();
                                                double xth = 1;
                                                if(METproj>0) xth=tau.Pt()/(tau.Pt() + METproj);
                                                else xth = 1;
                                                double mass_vis = (tau + muon).M();
                                                double mcol = 0;
                                                if (mass_vis != mass_vis) mass_vis=0;
                                                if (mass_vis <= 0.) mass_vis = 0.;
                                                mcol = mass_vis/sqrt(xth);
                                                return mcol;
                                                
                                            }
                                            return -1.;

                                 
                                 },{"GenVisTau_pt","GenVisTau_eta","GenVisTau_phi","GenVisTau_mass","sel_Tau_genPartFlav", "GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass","sel_Tau_genIdx", "sel_Tau_genVisIdx","sel_Muon_genIdx","GenMET_pt","GenMET_phi"})
                                           .Define("CollMass_gen_sel", [](const rvec<float>& genpt,
                                                                        const rvec<float>& geneta, 
                                                                        const rvec<float>& genphi, 
                                                                        const rvec<float>& genmass, 
                                                                        const rvec<int>& genmotheridx,
                                                                        const rvec<float>& genmuonpt,
                                                                        const rvec<float>& genmuoneta,
                                                                        const rvec<float>& genmuonphi,
                                                                        const rvec<float>& genmuonmass,
                                                                        const rvec<int>& genmuonpdgid,
                                                                        const rvec<int>& genmuonmotheridx,
                                                                        const rvec<int>& genmuonstatus,
                                                                        const float& met_pt,
                                                                        const float& met_phi){
                                                        TLorentzVector tau, muon, met, highest_tau, highest_muon;
                                                        float highest_mass = 0;
                                                        int best_i = 0;
                                                        int best_j = 0;
                                                        for (unsigned int i = 0; i < genmuonpt.size(); i++){
                                                                if( abs(genmuonpdgid[i]) == 15 and abs(genmuonpdgid[genmuonmotheridx[i]]) == 32){
                                                                    tau.SetPtEtaPhiM(genmuonpt[i], genmuoneta[i], genmuonphi[i], genmuonmass[i]);
                                                                    double mass_vis = 0;
                                                                    for(unsigned int j = 0; j < genmuonpt.size(); j++){
                                                                        if( abs(genmuonpdgid[j]) == 13 and abs(genmuonpdgid[genmuonmotheridx[j]]) == 32){
                                                                            muon.SetPtEtaPhiM(genmuonpt[j], genmuoneta[j], genmuonphi[j], genmuonmass[j]);
                                                                            mass_vis = (tau + muon).M();
                                                                            if ( highest_mass < mass_vis){
                                                                                highest_mass = mass_vis;
                                                                                highest_tau = tau;
                                                                                highest_muon = muon;
                                                                                best_i = i;
                                                                                best_j = j;
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                        }
                                                        return highest_mass;
                                             },{"GenVisTau_pt","GenVisTau_eta","GenVisTau_phi","GenVisTau_mass","GenVisTau_genPartIdxMother","GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass","GenPart_pdgId","GenPart_genPartIdxMother","GenPart_statusFlags","GenMET_pt","GenMET_phi"})
           
                                 .Define("sel_Tau_genPt", [](const rvec<float>& genpt,
                                                             const int& Tau_genIdx){
                                                                return genpt[Tau_genIdx];
                                                             },{"GenVisTau_pt","sel_Tau_genIdx"});
        }
        
        auto dimuon_cut = defquants.Filter(dimuonpair_cut,{"Muon_pt_corr" + muon_shift,"Muon_eta","Muon_phi", "DiMuon_mask"},"Di Muon Pair Cut");
        auto delR_cut = dimuon_cut.Filter(mutaudelR_cut,{"sel_Muon_eta","sel_Muon_phi", "sel_Tau_eta", "sel_Tau_phi"},"Mu Tau Del R Cut");
        auto df_jet_matching = delR_cut.Define("matched_jet",jet_matching, {"sel_Tau_pt","sel_Tau_eta","sel_Tau_phi","sel_Tau_mass","Jet_pt","Jet_eta","Jet_phi","Jet_mass"});
        auto jet_matching_cut = df_jet_matching.Filter("matched_jet == true","jet matching cut");


        auto corr_met_df  = jet_matching_cut;
        
        // Calculate MT distribution
        //auto mtcalc = corr_met_df.Define("MT", mass_transv, {"sel_Muon_pt", "sel_Muon_phi", met_branch_name + "_pt" + met_shift, met_branch_name + "_phi" + met_shift});
        auto mtcalc = jet_matching_cut.Define("MT", mass_transv, {"sel_Muon_pt", "sel_Muon_phi", "sel_MET_pt" + met_shift ,"sel_MET_phi" + met_shift});
        auto df_coll = mtcalc.Define("CollMass", collinear_mass, {"sel_Tau_pt", "sel_Muon_pt", "sel_Tau_eta", "sel_Muon_eta", "sel_Tau_phi", "sel_Muon_phi", "sel_Tau_mass", "sel_Muon_mass",  "sel_MET_pt" + met_shift,  "sel_MET_phi" + met_shift})  
                             .Define("CollMass_full", collinear_mass, {"sel_Tau_pt", "sel_Muon_pt", "sel_Tau_eta", "sel_Muon_eta", "sel_Tau_phi", "sel_Muon_phi", "sel_Tau_mass", "sel_Muon_mass", "sel_MET_pt_full" + met_shift, "sel_MET_phi_full" + met_shift});  
        auto df_mll = df_coll;
        // Define and calculate weights for monte carlo
        RNode total_weights = df_mll;
        std::vector < std::string > list_of_weights;
        if (config::runOnData) {
            total_weights = total_weights.Define("total_weight", "1.0");
            list_of_weights.push_back("total_weight");
        } else {    
                total_weights = init_PDFs(total_weights);
                
                
                calc_tau_uncertainty(total_weights, "Jet", "", "", "");
                calc_tau_uncertainty(total_weights, "Ele", "", "", "");
                calc_tau_uncertainty(total_weights, "Muon", "", "", "");
                
                
                std::vector < std::string > default_weights = {
                    "pileup_weight",
                    "genWeight",
                    "TauJetFakeScaleFactor__",
                    "TauEleFakeScaleFactor__",
                    "TauMuonFakeScaleFactor__",
                    "PrefiringWeight",
                    "TriggerScaleFactor",
                    "MuonISOScaleFactor",
                    "MuonIDScaleFactor",
                    "WWshape",
                    "TopQ",
                    "TopPDF",
                    "top_pt_weight",
                    "MuonReco"
                };


            total_weights = total_weights.Define("MuonISOScaleFactor", [](const rvec<float>& muon_pt,
                                                                          const rvec<float>& muon_eta,
                                                                          const rvec<bool>& muon_mask
                                                                          )
                                                                          {
                                                                             return muon_iso_scale_factor(muon_pt, muon_eta, muon_mask, "");
                                                                          }, {"Muon_pt_corr" + muon_shift, "Muon_eta", "Muon_mask"})
                                         .Define("MuonISOScaleFactorUp", [](const rvec<float>& muon_pt,
                                                                          const rvec<float>& muon_eta,
                                                                          const rvec<bool>& muon_mask
                                                                          )
                                                                          {
                                                                             return muon_iso_scale_factor(muon_pt, muon_eta, muon_mask, "Up");
                                                                          }, {"Muon_pt_corr" + muon_shift, "Muon_eta", "Muon_mask"})
                                         .Define("MuonISOScaleFactorDown", [](const rvec<float>& muon_pt,
                                                                              const rvec<float>& muon_eta,
                                                                              const rvec<bool>& muon_mask
                                                                              )
                                                                              {
                                                                                 return muon_iso_scale_factor(muon_pt, muon_eta, muon_mask, "Down");
                                                                              }, {"Muon_pt_corr" + muon_shift, "Muon_eta", "Muon_mask"})
                                         .Define("MuonIDScaleFactor", []( const rvec<float>& muon_pt,
                                                                          const rvec<float>& muon_eta,
                                                                          const rvec<bool>& muon_mask
                                                                          )
                                                                          {
                                                                             return muon_id_scale_factor(muon_pt, muon_eta, muon_mask, "");
                                                                          }, {"Muon_pt_corr" + muon_shift, "Muon_eta", "Muon_mask"})
                                         .Define("MuonIDScaleFactorUp", []( const rvec<float>& muon_pt,
                                                                          const rvec<float>& muon_eta,
                                                                          const rvec<bool>& muon_mask
                                                                          )
                                                                          {
                                                                             return muon_id_scale_factor(muon_pt, muon_eta, muon_mask, "Up");
                                                                          }, {"Muon_pt_corr" + muon_shift, "Muon_eta", "Muon_mask"})
                                         .Define("MuonIDScaleFactorDown", []( const rvec<float>& muon_pt,
                                                                          const rvec<float>& muon_eta,
                                                                          const rvec<bool>& muon_mask
                                                                          )
                                                                          {
                                                                             return muon_id_scale_factor(muon_pt, muon_eta, muon_mask, "Down");
                                                                          }, {"Muon_pt_corr" + muon_shift, "Muon_eta", "Muon_mask"})
                                         .Define("WWshape", [](){ return (float) 1.;})
                                         .Define("WWshapeUp", []( const float& coll_mass
                                                                          )
                                                                          {
                                                                             if (config::wwuncertainty){
                                                                                 float weight = 1 + 1 - (0.993-0.00002001 * coll_mass + 0.00000000283 * coll_mass * coll_mass);
                                                                                 return weight ;
                                                                             }
                                                                             return (float) 1.;
                                                                          }, {"CollMass"})
                                         .Define("WWshapeDown", []( const float& coll_mass
                                                                          )
                                                                          {
                                                                             if (config::wwuncertainty){
                                                                                 float weight = (0.993-0.00002001 * coll_mass + 0.00000000283 * coll_mass * coll_mass);
                                                                                 return weight;
                                                                             }
                                                                             return (float) 1.;
                                                                          }, {"CollMass"})
                                         .Define("TopQ", []( const float& Mll
                                                                          )
                                                                          {
                                                                             if (config::TT){
                                                                                 return GetTopQscale(Mll,"nom");
                                                                             }
                                                                             return (float) 1.;
                                                                          }, {"Mll"})
                                         .Define("TopQUp", []( const float& Mll
                                                                          )
                                                                          {
                                                                             if (config::TT){
                                                                                 return GetTopQscale(Mll,"up");
                                                                             }
                                                                             return (float) 1.;
                                                                          }, {"Mll"})
                                         .Define("TopQDown", []( const float& Mll
                                                                          )
                                                                          {
                                                                             if (config::TT){
                                                                                 return GetTopQscale(Mll,"down");
                                                                             }
                                                                             return (float) 1.;
                                                                          }, {"Mll"})
                                         .Define("TopPDF", []( const float& Mll
                                                                          )
                                                                          {
                                                                             if (config::TT){
                                                                                 return GetTopPDF(Mll,"nom");
                                                                             }
                                                                             return (float) 1.;
                                                                         }, {"Mll"})
                                         .Define("TopPDFUp", []( const float& Mll
                                                                          )
                                                                          {
                                                                             if (config::TT){
                                                                                 return GetTopPDF(Mll,"up");
                                                                             }
                                                                             return (float) 1.;
                                                                          }, {"Mll"})
                                         .Define("TopPDFDown", []( const float& Mll
                                                                          )
                                                                          {
                                                                             if (config::TT){
                                                                                 return GetTopPDF(Mll,"down");
                                                                             }
                                                                             return (float) 1.;
                                                                          }, {"Mll"})
                                    .Define("TriggerScaleFactor", [](const rvec<float>& pt,
                                                                     const rvec<float>& eta,
                                                                     const rvec<bool>& mask
                                                                     ) { 
                                                                         return trigger_sf(pt, eta,mask, "");
                                                                    }, {"Muon_pt_corr"+ muon_shift, "Muon_eta", "Muon_mask"})
                                    .Define("TriggerScaleFactorUp", [](const rvec<float>& pt,
                                                                     const rvec<float>& eta,
                                                                     const rvec<bool>& mask
                                                                     ) { 
                                                                         return trigger_sf(pt, eta, mask, "Up");
                                                                    }, {"Muon_pt_corr", "Muon_eta", "Muon_mask"})
                                    .Define("TriggerScaleFactorDown", [](const rvec<float>& pt,
                                                                     const rvec<float>& eta,
                                                                     const rvec<bool>& mask
                                                                     ) { 
                                                                         return trigger_sf(pt, eta, mask, "Down");
                                                                    }, {"Muon_pt_corr" + muon_shift, "Muon_eta", "Muon_mask"})
                                    .Define("MuonReco", [](const rvec<float>& pt,
                                                                     const rvec<float>& eta,
                                                                     const rvec<float>& phi,
                                                                     const rvec<float>& mass,
                                                                     const rvec<bool>& mask
                                                                     ) { 
                                                                         return muon_reco_eff(pt, eta, phi, mass,mask,"");
                                                                    }, {"Muon_pt_corr" + muon_shift, "Muon_eta", "Muon_phi", "Muon_mass", "Muon_mask"})
                                    .Define("MuonRecoUp", [](const rvec<float>& pt,
                                                                     const rvec<float>& eta,
                                                                     const rvec<float>& phi,
                                                                     const rvec<float>& mass,
                                                                     const rvec<bool>& mask
                                                                     ) { 
                                                                         return muon_reco_eff(pt, eta, phi, mass,mask,"Up");
                                                                    }, {"Muon_pt_corr" + muon_shift, "Muon_eta", "Muon_phi", "Muon_mass", "Muon_mask"})
                                    .Define("MuonRecoDown", [](const rvec<float>& pt,
                                                                     const rvec<float>& eta,
                                                                     const rvec<float>& phi,
                                                                     const rvec<float>& mass,
                                                                     const rvec<bool>& mask
                                                                     ) { 
                                                                         return muon_reco_eff(pt, eta, phi, mass,mask,"Down");
                                                                    }, {"Muon_pt_corr"+ muon_shift , "Muon_eta", "Muon_phi", "Muon_mass", "Muon_mask"});
            define_weight(total_weights, default_weights, list_of_weights);

            if (default_run) {
                    //DEFINE NON-ENERGY RELATED UNCERTAINTIES
                    total_weights = total_weights.Define("pileup_weightUp", [](const int& nvtx_true){ return pu_weight(nvtx_true, "Up"); }, {"Pileup_nPU"})
                                                 .Define("pileup_weightDown", [](const int& nvtx_true){ return pu_weight(nvtx_true, "Down"); }, {"Pileup_nPU"})
                                                 .Define("prefiring_weight_Up", [](const rvec<float>& jet_pt,
                                                                         const rvec<float>& jet_eta,
                                                                         const rvec<float>& photon_pt,
                                                                         const rvec<float>& photon_eta) { 
                                                                             return prefire_factor(jet_pt, jet_eta, photon_pt, photon_eta, "Up");
                                                                        }, {"Jet_pt", "Jet_eta", "Photon_pt", "Photon_eta"})
                                                 .Define("prefiring_weight_Down", [](const rvec<float>& jet_pt,
                                                                         const rvec<float>& jet_eta,
                                                                         const rvec<float>& photon_pt,
                                                                         const rvec<float>& photon_eta) { 
                                                                             return prefire_factor(jet_pt, jet_eta, photon_pt, photon_eta, "Down");
                                                                    }, {"Jet_pt", "Jet_eta", "Photon_pt", "Photon_eta"});

                    calc_tau_uncertainty(total_weights, "Jet", "Up", "", "");
                    calc_tau_uncertainty(total_weights, "Jet", "Down", "", "");
                    
                    calc_tau_uncertainty(total_weights, "Ele", "Up", "eta", "barrel");
                    calc_tau_uncertainty(total_weights, "Ele", "Down", "eta", "barrel");
                    calc_tau_uncertainty(total_weights, "Ele", "Up", "eta", "endcap");
                    calc_tau_uncertainty(total_weights, "Ele", "Down", "eta", "endcap");
                    
                    calc_tau_uncertainty(total_weights, "Muon", "Up", "eta", "eta0p4");
                    calc_tau_uncertainty(total_weights, "Muon", "Up", "eta", "eta0p4to0p8");
                    calc_tau_uncertainty(total_weights, "Muon", "Up", "eta", "eta0p8to1p2");
                    calc_tau_uncertainty(total_weights, "Muon", "Up", "eta", "eta1p2to1p7");
                    calc_tau_uncertainty(total_weights, "Muon", "Up", "eta", "eta1p7");
                    calc_tau_uncertainty(total_weights, "Muon", "Down", "eta", "eta0p4");
                    calc_tau_uncertainty(total_weights, "Muon", "Down", "eta", "eta0p4to0p8");
                    calc_tau_uncertainty(total_weights, "Muon", "Down", "eta", "eta0p8to1p2");
                    calc_tau_uncertainty(total_weights, "Muon", "Down", "eta", "eta1p2to1p7");
                    calc_tau_uncertainty(total_weights, "Muon", "Down", "eta", "eta1p7");
                    
                
                    define_weight(total_weights, default_weights, list_of_weights);
                    define_weight(total_weights, default_weights, list_of_weights, "pileup_weightUp", "pileup_weight");
                    define_weight(total_weights, default_weights, list_of_weights, "pileup_weightDown", "pileup_weight");
                    define_weight(total_weights, default_weights, list_of_weights, "TriggerScaleFactorUp", "TriggerScaleFactor");
                    define_weight(total_weights, default_weights, list_of_weights, "TriggerScaleFactorDown", "TriggerScaleFactor");
                    define_weight(total_weights, default_weights, list_of_weights, "MuonISOScaleFactorUp", "MuonISOScaleFactor");
                    define_weight(total_weights, default_weights, list_of_weights, "MuonISOScaleFactorDown", "MuonISOScaleFactor");
                    define_weight(total_weights, default_weights, list_of_weights, "MuonIDScaleFactorUp", "MuonIDScaleFactor");
                    define_weight(total_weights, default_weights, list_of_weights, "MuonIDScaleFactorDown", "MuonIDScaleFactor");
                    define_weight(total_weights, default_weights, list_of_weights, "WWshapeUp", "WWshape");
                    define_weight(total_weights, default_weights, list_of_weights, "WWshapeDown", "WWshape");
                    define_weight(total_weights, default_weights, list_of_weights, "TopQUp", "TopQ");
                    define_weight(total_weights, default_weights, list_of_weights, "TopQDown", "TopQ");
                    define_weight(total_weights, default_weights, list_of_weights, "TopPDFUp", "TopPDF");
                    define_weight(total_weights, default_weights, list_of_weights, "TopPDFDown", "TopPDF");
                    define_weight(total_weights, default_weights, list_of_weights, "MuonRecoUp", "MuonReco");
                    define_weight(total_weights, default_weights, list_of_weights, "MuonRecoDown", "MuonReco");

                    define_weight(total_weights, default_weights, list_of_weights, "TauJetFakeScaleFactorUp__", "TauJetFakeScaleFactor__");
                    define_weight(total_weights, default_weights, list_of_weights, "TauJetFakeScaleFactorDown__", "TauJetFakeScaleFactor__");

                    define_weight(total_weights, default_weights, list_of_weights, "TauEleFakeScaleFactorUp_eta_barrel", "TauEleFakeScaleFactor__");
                    define_weight(total_weights, default_weights, list_of_weights, "TauEleFakeScaleFactorDown_eta_barrel", "TauEleFakeScaleFactor__");
                    define_weight(total_weights, default_weights, list_of_weights, "TauEleFakeScaleFactorUp_eta_endcap", "TauEleFakeScaleFactor__");
                    define_weight(total_weights, default_weights, list_of_weights, "TauEleFakeScaleFactorDown_eta_endcap", "TauEleFakeScaleFactor__");

                    define_weight(total_weights, default_weights, list_of_weights, "TauMuonFakeScaleFactorUp_eta_eta0p4", "TauMuonFakeScaleFactor__");
                    define_weight(total_weights, default_weights, list_of_weights, "TauMuonFakeScaleFactorDown_eta_eta0p4", "TauMuonFakeScaleFactor__");
                    define_weight(total_weights, default_weights, list_of_weights, "TauMuonFakeScaleFactorUp_eta_eta0p4to0p8", "TauMuonFakeScaleFactor__");
                    define_weight(total_weights, default_weights, list_of_weights, "TauMuonFakeScaleFactorDown_eta_eta0p4to0p8", "TauMuonFakeScaleFactor__");
                    define_weight(total_weights, default_weights, list_of_weights, "TauMuonFakeScaleFactorUp_eta_eta0p8to1p2", "TauMuonFakeScaleFactor__");
                    define_weight(total_weights, default_weights, list_of_weights, "TauMuonFakeScaleFactorDown_eta_eta0p8to1p2", "TauMuonFakeScaleFactor__");
                    define_weight(total_weights, default_weights, list_of_weights, "TauMuonFakeScaleFactorUp_eta_eta1p2to1p7", "TauMuonFakeScaleFactor__");
                    define_weight(total_weights, default_weights, list_of_weights, "TauMuonFakeScaleFactorDown_eta_eta1p2to1p7", "TauMuonFakeScaleFactor__");
                    define_weight(total_weights, default_weights, list_of_weights, "TauMuonFakeScaleFactorUp_eta_eta1p7", "TauMuonFakeScaleFactor__");
                    define_weight(total_weights, default_weights, list_of_weights, "TauMuonFakeScaleFactorDown_eta_eta1p7", "TauMuonFakeScaleFactor__");
                    define_weight(total_weights, default_weights, list_of_weights, "prefiring_weight_Up", "PrefiringWeight");
                    define_weight(total_weights, default_weights, list_of_weights, "prefiring_weight_Down", "PrefiringWeight");

                
                                // factorization & renormalization scale uncertainty
                    std::vector<std::string> columnNames = total_weights.GetColumnNames();
                    if (std::find(columnNames.begin(), columnNames.end(), "LHEScaleWeight") != columnNames.end()) {
                        total_weights = total_weights.Define("total_weight_ScaleWeightUp", [](const rvec<float>& weights, const float& tot_weight) {
                                                                                                            float max_val = 1.0;
                                                                                                            for (uint i = 0; i < weights.size(); i++) {
                                                                                                                if ((i==5) || (i==7)) {
                                                                                                                    continue;
                                                                                                                }
                                                                                                                if (max_val < weights[i]) {
                                                                                                                    max_val = weights[i];
                                                                                                                }
                                                                                                            }
                                                                                                            return max_val*tot_weight;
                                                                                            }, {"LHEScaleWeight", "total_weight"})
                                                                                    .Define("total_weight_ScaleWeightDown", [](const rvec<float>& weights, const float& tot_weight) {
                                                                                                            float min_val = 1.0;
                                                                                                            for (uint i = 0; i < weights.size(); i++) {
                                                                                                                if ((i==5) || (i==7)) {
                                                                                                                    continue;
                                                                                                                }
                                                                                                                if (min_val > weights[i]) {
                                                                                                                    min_val = weights[i];
                                                                                                                }
                                                                                                            }
                                                                                                            return min_val*tot_weight;
                                                                                            }, {"LHEScaleWeight", "total_weight"});
                        
                    } else {
                        total_weights = total_weights.Define("total_weight_ScaleWeightUp", "total_weight")
                                                     .Define("total_weight_ScaleWeightDown", "total_weight");
                    }
                    list_of_weights.push_back("total_weight_ScaleWeightUp");
                    list_of_weights.push_back("total_weight_ScaleWeightDown");
                    
                    // pdf weights
                    std::vector<std::string> colNames = total_weights.GetColumnNames();
                    if (std::find(colNames.begin(), colNames.end(), "LHEPdfWeight_def") != colNames.end()) {
                        if (!config::runOnSignal) {
                            if (std::find(colNames.begin(), colNames.end(), "LHEWeight_originalXWGTUP") == colNames.end()) {
                                total_weights = total_weights.Define("LHEWeight_originalXWGTUP", [](const rvec<float>& pdf_weights){return pdf_weights[0];}, {"LHEPdfWeight_def"});
                            }
                            for (uint i = 0; i < config::pdf_nweights; i++) {
                                total_weights = total_weights.Define("total_weight_pdf_weight_" + std::to_string(i), [i](
                                        const float& weight, 
                                        const rvec<float>& pdf_weights, 
                                        const float& orig_weight)
                                    {
                                        if (pdf_weights.size() != config::pdf_nweights) {
                                            return weight;
                                        }
                                        //~ return weight*pdf_weights[i]/orig_weight;
                                        return weight*pdf_weights[i]/pdf_weights[0];
                                    }, {"total_weight", "LHEPdfWeight_def", "LHEWeight_originalXWGTUP"});
                                list_of_weights.push_back("total_weight_pdf_weight_" + std::to_string(i));
                            }
                        }
                    }
                } // default_run
            } // !runOnData

        

            
        //auto df_coll_alt = df_coll.Define("CollMass_alt", collinear_mass_alt, {"col_idx","Tau_pt_ES", "Muon_tP_pt", "Tau_eta_ES", "Muon_eta", "Tau_phi_ES", "Muon_phi", "Tau_mass_ES", "Muon_mass", "MET_pt", "MET_phi","Tau_mask", "Muon_mask"});  
        // Define Stage 0: any event with one tau fulfilling acceptance and id	
        //create_hists(df_coll_alt, "Stage0", "total_weight");
        
        // actual analysis cut  -- MT > 120. GeV 
        //auto df_mt = df_coll_alt.Filter("MT > 120.", "mt_cut");
        auto df_mt = total_weights.Filter("MT > 120.", "mt_cut").Define("Sphericity",sphericity_full,{"Tau_pt","Tau_eta","Tau_phi","Tau_mass","Muon_pt_corr" + muon_shift,"Muon_eta","Muon_phi","Muon_mass","Jet_pt_nom","Jet_eta", "Jet_phi", "Jet_mass","Tau_mask", "Muon_mask"}).Define("Sphericity_leptons",sphericity_leptons,{"Tau_pt","Tau_eta","Tau_phi","Tau_mass","Muon_pt_corr" + muon_shift,"Muon_eta","Muon_phi","Muon_mass","Tau_mask", "Muon_mask"});
        auto df_dip = df_mt.Filter("CollMass > 800. && CollMass < 1000.", "dip");

        auto df_delphi_taumet = df_mt.Define("DeltaPhi_tau_met", delta_phi, {"sel_Tau_phi", met_branch_name + "_phi"});  
        //auto df_njets = df_mt.Filter("nJet > 9");  
        auto df_delphi_cut_test = df_delphi_taumet.Filter("DeltaPhi_tau_met < 1.5707963267948966","delphi cut");
        auto df_delphi_cut = df_delphi_cut_test.Filter("abs(DeltaPhi_tau_met) < 1.5707963267948966","delphi cut");
        //create_hists(df_coll_alt, "Stage0", "total_weight");
        if ( default_run ) {
            // fill preselection
            create_preselection_hists(preselection);
            
            //// tau efficiency study
            //create_tau_efficiency_hists(masked);
            //create_tau_efficiency_hists(df_dphi, "_PostSelection");
            
                
            //create_hists(total_weights, "Stage0", tau_shift, met_shift, "total_weight");
            
            //create_hists(df_mt, "Stage0", muon_shift, tau_shift, met_shift, "total_weight");
            if(!config::runOnData && config::doSnapshot){
                ROOT::RDF::RSnapshotOptions opts;
                opts.fLazy = false;
                std::string snapshot;
                total_weights.Snapshot
                //total_weights.Snapshot<
                //              UInt_t,
                //              UInt_t,
                //               rvec<float>,
                //               rvec<float>,
                //               rvec<float>,
                //               rvec<float>,
                //               rvec<float>,
                //               rvec<float>,
                //               rvec<float>,
                //               rvec<float>,
                //               rvec<float>,
                //               rvec<float>,
                //               rvec<float>,
                //               rvec<float>,
                //               rvec<float>,
                //               rvec<float>,
                //               rvec<float>,
                //               rvec<float>,
                //               rvec<float>,
                //               rvec<float>,
                //               rvec<float>,
                //               rvec<float>,
                //               float,
                //               float,
                //               float,
                //               float,
                //               float,
                //               float,
                //               float,
                //               float,
                //               float,
                //               float,
                //               float,
                //               float,
                //               float,
                //               float,
                //              // float,
                //              // float,
                //              // float,
                //              // float,
                //              // float,
                //              // float,
                //              // float,
                //              // float,
                //               float,
                //               double,
                //               float,
                //               float,
                //               rvec<bool>,
                //               rvec<bool>,
                //               rvec<int>,
                //               rvec<UChar_t>,
                //               float,
                //               rvec<float>,
                //               rvec<int>,
                //               rvec<bool>,
                //               rvec<int>,
                //               rvec<int>,
                //               rvec<UChar_t>,
                //               rvec<UChar_t>,
                //               rvec<UChar_t>,
                //               rvec<bool>,
                //               rvec<UChar_t>,
                //               rvec<UChar_t>
                //               >
                               ("Events","output/" + std::to_string(config::era) + "/" + runName + "_snapshot.root",
                               {
                               "nTau",
                               "nMuon",
                               "nGenPart",
                               "nGenVisTau",
                               "nJet",
                               "Muon_pt_corr",
                                "Muon_pt_corr_muonResolutionUp",
                                "Muon_pt_corr_muonResolutionDown",
                                "Muon_pt_corr_muonScaleUp",
                                "Muon_pt_corr_muonScaleDown",
                                "Muon_eta",
                                "Muon_phi",
                                "Muon_mass", 
                                "Tau_pt_ES",
                                "Tau_eta_ES",
                                "Tau_phi_ES", 
                                "Tau_mass_ES", 
                                "Tau_pt_ESUp",
                                "Tau_eta_ESUp",
                                "Tau_phi_ESUp", 
                                "Tau_mass_ESUp",
                                "Tau_pt_ESDown",
                                "Tau_eta_ESDown",
                                "Tau_phi_ESDown", 
                                "Tau_mass_ESDown",
                                met_branch_name + "_pt_jer" ,
                                met_branch_name + "_phi_jer",
                                met_branch_name + "_pt_jerUp",
                                met_branch_name + "_pt_jerDown",
                                met_branch_name + "_phi_jerUp",
                                met_branch_name + "_phi_jerDown",
                                met_branch_name + "_pt_jesTotalUp" ,
                                met_branch_name + "_pt_jesTotalDown",
                                met_branch_name + "_phi_jesTotalUp" ,
                                met_branch_name + "_phi_jesTotalDown" ,
                                met_branch_name + "_pt_unclustEnUp",
                                met_branch_name + "_pt_unclustEnDown",
                                met_branch_name + "_phi_unclustEnUp",
                                met_branch_name + "_phi_unclustEnDown",
                                "CollMass",
                               // "CollMass_ESUp",
                               // "CollMass_ESDown",
                               // "CollMass_muonResolutionUp",
                               // "CollMass_muonResolutionDown",
                               // "CollMass_muonScaleUp",
                               // "CollMass_muonScaleDown",
                               // "CollMass_jerUp",
                               // "CollMass_jerDown",
                               // "CollMass_jesTotalUp",
                               // "CollMass_jesTotalDown",
                               // "CollMass_unclustEnUp",
                               // "CollMass_unclustEnDown",
                                "CollMass_gen",
                                "genWeight",
                                "total_weight",
                                "Tau_mask",
                                "Muon_mask",
                                "col_idx",
                                "Tau_genPartIdx",
                                "Tau_genPartFlav",
                                "sel_Tau_genPartFlav",
                                "Tau_dz",
                                "Tau_charge",
                                config::tau_dm,
                                "Tau_decayMode",
                                "Tau_jetIdx",
                                "Tau_idDeepTau2017v2p1VSjet",
                                "Tau_idDeepTau2017v2p1VSe",
                                "Tau_idDeepTau2017v2p1VSmu",
                                "Muon_isPFcand",
                                "Muon_highPtId",
                                "Muon_tkIsoId",
                                "Muon_genPartIdx",
                                "MT",
                                "GenVisTau_pt",
                                "GenVisTau_eta",
                                "GenVisTau_phi",
                                "GenVisTau_mass",
                                "GenVisTau_status",
                                "GenPart_pt",
                                "GenPart_phi",
                                "GenPart_eta",
                                "GenPart_mass",
                                "GenPart_statusFlags",
                                "GenPart_genPartIdxMother",
                                "sel_Tau_genIdx",
                                "sel_Tau_genVisIdx",
                                "sel_Muon_genIdx",
                                "sel_Muon_pt",
                                "sel_Muon_eta",
                                "sel_Muon_phi",
                                "sel_Muon_mass",
                                "GenMET_pt",
                                "GenMET_phi",
                                "GenPart_pdgId",
                                "GenPart_status",
                                "Tau_jetIdx",
                                "Jet_pt",
                                "Jet_pt_nom",
                                "Jet_eta",
                                "Jet_phi",
                                "Jet_mass"

                                },opts);
                fill_gen_hists(df_mt);
            }
         } 
        fill_stage_with_syst(df_mt, "Stage1", muon_shift, tau_shift, met_shift, list_of_weights);
        fill_stage_with_syst( df_delphi_cut, "Stage2",muon_shift, tau_shift, met_shift, list_of_weights);
        return df_delphi_cut;
        }
   }();
   //auto df_report = df_final.Report();
   //df_report->Print();
   return df_final;

}

std::vector < std::string > readfiles(const char* file_directory) {
	TSystemDirectory dir(file_directory, file_directory);
	TList *files = dir.GetListOfFiles();
	std::vector < std::string > rootfiles;
	if (files) {
		TSystemFile *file;
		TString fname;
		TIter next(files);
		while ((file=(TSystemFile*)next())) {
			fname = file->GetName();
			if (!file->IsDirectory() and fname.EndsWith(".root")) {
				std::string full_path = file_directory;
				full_path += file->GetName();
                rootfiles.push_back(full_path);
			}
		}
	}
	return rootfiles;
}
// analyse function - gets called for each systematic
void analyse(   RNode df,
                TFile* outFile) {
    // init counter
    df = df.Define("GenWeight",[](const float genWeight){
                                            return (double) genWeight;
                                            },{"genWeight"});
    std::vector< ROOT::RDF::RResultPtr< double >  > genWeightsSum;
    genWeightsSum.push_back(df.Sum< double >("GenWeight"));
    auto eventcounter = df.Count();
    auto genweighthist = df.Histo1D({"Genweight","Genweight",1000u,-10,10},"genWeight");
    
    hist_dict.emplace("",genweighthist);
    // generic analysis part -------------------------------------------
    auto met_filter = df.Filter(config::metfilters, "MET filters");
    
    auto good_primary_vertex = met_filter.Filter("PV_npvsGood > 1", "has good PV");
    //auto good_primary_vertex = met_filter;
    
    // handle 2018 HEM problem
    RNode fixed = good_primary_vertex;
    //if (config::era == 2018) {
    //    fixed = good_primary_vertex.Filter([](   const rvec<float>& jet_pt,
    //                                    const rvec<float>& jet_eta,
    //                                    const rvec<float>& jet_phi){
    //                                        for (uint i = 0; i < jet_pt.size(); i++) {
    //                                            if (jet_pt[i] < 30)
    //                                                continue;
    //                                            if ( jet_eta[i] < -3.0 || jet_eta[i] > -1.3 )
    //                                                continue;
    //                                            if ( jet_phi[i] < -1.57 || jet_phi[i] > -0.87 )
    //                                                continue;
    //                                            return false;
    //                                        }
    //                                        return true;
    //                                    }, {"Jet_pt_nom", "Jet_eta", "Jet_phi"}, "2018 HEM problem");
    //}
    if (!config::runOnData) {
            // Signal Study
            fixed = init_PDFs(fixed);
        }
        
        //std::vector<double> pdf_all( config::pdf_nweights, 0);
        //std::vector<double> pdf_passed( config::pdf_nweights, 0);
        if (config::runOnSignal) {
            fixed = fixed.Define("hasTauInAcceptance", [](
                        const rvec<float>& genVisTau_pt,
                        const rvec<float>& genVisTau_eta){
                            bool hasGenTauInAcceptance = false;
                            for (uint i = 0; i < genVisTau_pt.size(); i++) {
                                if (genVisTau_pt[i] < config::tau_pt) {
                                    continue;
                                }
                                if (std::abs(genVisTau_eta[i]) > 2.3) {
                                    continue;
                                }
                                hasGenTauInAcceptance = true;
                                break;
                            }
                            return hasGenTauInAcceptance;
                            }, {"GenVisTau_pt", "GenVisTau_eta"})
                          .Define("hasMuonInAcceptance", [](
                                  const rvec<float>& GenPart_pt,
                                  const rvec<float>& GenPart_eta,
                                  const rvec<int>& GenPart_pdgId){
                                  bool hasGenMuonInAcceptance = false;
                                    for (uint i = 0; i < GenPart_pt.size(); i++) {
                                        if(abs(GenPart_pdgId[i]) == 13 && GenPart_pt[i] > config::muon_pt && abs(GenPart_eta[i]) < config::muon_eta){
                                            hasGenMuonInAcceptance = true;
                                            break;
                                        }
                                    }
                                    return hasGenMuonInAcceptance;
                                  
                                  }, {"GenPart_pt","GenPart_eta","GenPart_pdgId"});
            fixed = fixed.Define("pdf_passed", [](
                                        const rvec<float>& pdf_weights,
                                        const bool hasTauAcceptance,
                                        const bool hasMuonAcceptance)
                                    {
                                        rvec<float> tmp( pdf_weights.size(), 0);
                                        for (uint i = 0; i < pdf_weights.size(); i++) {
                                            tmp[i] = 0;
                                            if (hasTauAcceptance && hasMuonAcceptance) {
                                                //if(pdf_weights[i] >= 0.){
                                                    tmp[i] = pdf_weights[i];
                                               // }
                                               // else{
                                               //     tmp[i] = 1.;
                                               // }
                                            }
                                        }
                                        return tmp;
                                    }, {"LHEPdfWeight_def", "hasTauInAcceptance", "hasMuonInAcceptance"});
            }
    
    // END generic analysis part ---------------------------------------
    
        
    std::vector < std::string > tau_shifts = 
    { 
        "_ES",
        "_ESUp",
        "_ESDown"
    };
    
    std::vector < std::string > met_shifts = 
    { 
        "_jer",
        "_jerUp",
        "_jerDown",
        "_jesTotalUp",
        "_jesTotalDown",
        "_unclustEnUp",
        "_unclustEnDown"
    };
    std::vector < std::string > muon_shifts =
    {
        "",
        "_muonResolutionUp",
        "_muonResolutionDown",
        "_muonScaleUp",
        "_muonScaleDown"

    };
    //auto df_analyzed = fixed;
    
    auto final_note = fixed;
    auto base = fixed;
    for (const auto& tau_shift: tau_shifts) {
        if (!config::runOnData) {
            base = base.Define("Tau" + tau_shift + "_vectors", [tau_shift](const rvec<int>& tau_decayMode,
                                            const rvec<UChar_t>& tau_genPartFlav,
                                            const rvec<float>& tau_pt,
                                            const rvec<float>& tau_eta,
                                            const rvec<float>& tau_phi,
                                            const rvec<float>& tau_mass) {
                                                return calc_tau_energy_scale(tau_decayMode, tau_genPartFlav, tau_pt, tau_eta, tau_phi, tau_mass, tau_shift);
                                            }, {"Tau_decayMode", "Tau_genPartFlav", "Tau_pt", "Tau_eta", "Tau_phi", "Tau_mass"});
            
            base = base.Define("Tau_pt" + tau_shift, [](const rvec< rvec < float > >& es_vector){return apply_tau_energy_scale(es_vector, 0);}, {"Tau" + tau_shift + "_vectors"})
                        .Define("Tau_eta" + tau_shift, [](const rvec< rvec < float > >& es_vector){return apply_tau_energy_scale(es_vector, 1);}, {"Tau" + tau_shift + "_vectors"})
                        .Define("Tau_phi" + tau_shift, [](const rvec< rvec < float > >& es_vector){return apply_tau_energy_scale(es_vector, 2);}, {"Tau" + tau_shift + "_vectors"})
                        .Define("Tau_mass" + tau_shift, [](const rvec< rvec < float > >& es_vector){return apply_tau_energy_scale(es_vector, 3);}, {"Tau" + tau_shift + "_vectors"});
        }
        
        if (tau_shift == "_ES") {
            if (config::runOnData) {
                base = base.Define("Tau_pt_ES", [](const rvec<float>& tau_pt){return tau_pt;}, {"Tau_pt"})
                            .Define("Tau_eta_ES", [](const rvec<float>& tau_eta){return tau_eta;}, {"Tau_eta"})
                            .Define("Tau_phi_ES", [](const rvec<float>& tau_phi){return tau_phi;}, {"Tau_phi"})
                            .Define("Tau_mass_ES", [](const rvec<float>& tau_mass){return tau_mass;}, {"Tau_mass"});
            }
            for (const auto& muon_shift: muon_shifts) {
                if (config::runOnData) {
                    break;
                }
                else{
                    if( muon_shift != ""){
                        if(muon_shift == "_muonResolutionUp"){
                            base = base.Define("Muon_pt_corr" + muon_shift, [](const rvec<float>& muon_pt, const rvec<float>& muon_eta, const rvec<float>& muon_phi, const rvec<float>& muon_mass){return muon_reso_smearing(muon_pt,muon_eta,muon_phi,muon_mass,"Up");},{"Muon_pt_corr","Muon_eta","Muon_phi","Muon_mass"});
                        }
                        else if( muon_shift == "_muonResolutionDown"){
                            base = base.Define("Muon_pt_corr" + muon_shift, [](const rvec<float>& muon_pt, const rvec<float>& muon_eta, const rvec<float>& muon_phi, const rvec<float>& muon_mass){return muon_reso_smearing(muon_pt,muon_eta,muon_phi, muon_mass,"Down");},{"Muon_pt_corr","Muon_eta","Muon_phi","Muon_mass"});
                        }
                        else if (muon_shift == "_muonScaleUp" or muon_shift == "_muonScaleDown"){
                            base = base.Define("Muon_pt_corr" + muon_shift, [muon_shift](const rvec<float>& muon_pt, const rvec<float>& muon_eta, const rvec<float>& muon_phi, const rvec<int> q){return scale_correct_muon(muon_pt,muon_eta,muon_phi,q,muon_shift);},{"Muon_pt_corr","Muon_eta", "Muon_phi","Muon_charge"});
                        }
                    }
                }
            }
        } 
    }
    for (const auto& tau_shift: tau_shifts) {
        if(config::runOnData){
            if (tau_shift == "_ES") final_note = run_analyser(base, "", tau_shift, "_nom", genWeightsSum, outFile);
        }
        else{
            final_note = run_analyser(base, "", tau_shift, "_jer", genWeightsSum, outFile);
        }
    
    }
    for (const auto& muon_shift: muon_shifts) {
        if (config::runOnData) {
            break;
        } 
        if( muon_shift != ""){
            final_note = run_analyser(base, muon_shift, "_ES", "_jer", genWeightsSum, outFile);
        }
    }
    for (const auto& met_shift: met_shifts) {
        if (config::runOnData) {
            break;
        }
        else {
           if(met_shift == "_jer"){
               continue;
           }
           final_note = run_analyser(base, "", "_ES", met_shift, genWeightsSum, outFile);
        }
    }
    
    auto df_report = final_note.Report();
    df_report->Print();
    // write all hists
    outFile->cd();
    for (const auto& [folder_name, hist] : hist_dict) {
        if (outFile->GetDirectory(folder_name.c_str()) == 0) {
            outFile->mkdir(folder_name.c_str());
        }
        outFile->cd(folder_name.c_str());
        std::visit( Visitor{[](auto hist){hist->Write();}}, hist);
        outFile->cd();
    }
    if (config::runOnSignal) {
        rvec<float> base(config::pdf_nweights, 0);
        auto passed_sum = fixed.Reduce< rvec<float> (const rvec<float> &, rvec<float>)>([](const rvec<float>& weights, rvec<float> base) {
            rvec<float> result(weights.size(), 0);
            for (uint i = 0; i < weights.size(); i++) {
                result[i] = weights[i] + base[i];
            }
            return result;
            }, {"pdf_passed"}, base);
        auto all_sum = fixed.Reduce< rvec<float> (const rvec<float> &, rvec<float> )>([](const rvec<float>& weights, rvec<float> base) {
            rvec<float> result(weights.size(), 0);
            for (uint i = 0; i < weights.size(); i++) {
                //if (weights[i] > 0.){
                    result[i] = weights[i] + base[i];
               // }
                //else result[i] = 1. + base[i];

            }
            return result;
            }, {"LHEPdfWeight_def"}, base);
        
        double sum = 0; 
        double sum2 = 0; 
        for (uint i = 0; i < config::pdf_nweights; i++) {
            double ratio = passed_sum->at(i) / all_sum->at(i);
            sum += ratio;
            sum2 += ratio * ratio;
        }
        
        double pdf_signal_scale = 1./(config::pdf_nweights) * sqrt( (config::pdf_nweights) * sum2 - sum*sum  );
        //double pdf_signal_scale = 1./(npdfweights) * sqrt( (npdfweights) * sum2 - sum*sum  );

        auto pdf_copyhist = (TH1F*) outFile->Get("Stage1_nosplit/CollMass")->Clone();
        pdf_copyhist->Scale( 1+pdf_signal_scale );
        outFile->cd("Stage1_nosplit/Systematics");
        pdf_copyhist->SetName("CollMass_pdf_total");
        pdf_copyhist->Write();
        outFile->cd();
    }

    //if (config::runOnSignal) {
    //    TH1F* pdf_copyhist = (TH1F*) outFile->Get("Stage1/CollMass")->Clone();

    //    double sum = 0;
    //    double sum2 = 0;
    //    for (uint i = 0; i < pdf_all.size(); i++) {
    //        double ratio = pdf_passed[i] / pdf_all[i];
    //        sum += ratio;
    //        sum2 += ratio * ratio;
    //    }
    //    double pdf_signal_scale = 1./pdf_all.size() * sqrt( pdf_all.size()* sum2 - sum*sum  );
    //    pdf_copyhist->Scale( 1.+pdf_signal_scale );
    //    outFile->cd("Stage1/Systematics");
    //    pdf_copyhist->SetName("CollMass_pdf_total");
    //    pdf_copyhist->Write();
    //    outFile->cd();
    //}
    
    auto counter = TH1D("counter", "counter", genWeightsSum.size(), 0, genWeightsSum.size());
    for (uint i = 0; i < genWeightsSum.size(); i++) {
        counter.SetBinContent(i+1, genWeightsSum[i].GetValue() );
    }
    counter.SetEntries(*eventcounter );
    counter.Write();
}






int main (int argc, char* argv[]) {
//	// read in files to run over
//	json files;
//	std::ifstream files_json;
//	files_json.open(argv[1]);
//	files_json >> files;
//	
//	// loop over all types of backgrounds
//	for (auto& it : files.items()) {
//		
//		
//		// Create outfile name from path
//		TPRegexp r1("[^/]+(?=/$|$)");
//		TString tmp = (std::string) it.key();
//		tmp = "output/" + tmp(r1)+ ".root";
//		
//		std::cout << "\t \t Running on: " << tmp(r1) << std::endl;
//		
//		// create outputfile to write hists
//		TFile* outFile = new TFile(tmp, "RECREATE");
//		
//		// Files to run over
//		std::string path = it.key();
//		std::vector< std::string > names = readfiles(path.c_str());
//		
//		// read in config file for the corresponding file
//		std::ifstream cfg_json;
//		cfg_json.open(it.value());
//		cfg_json >> cfg;
//		
//		// load analysis setup config file (see config.cc)
//		config::load_config_file(cfg);
//			
//		// open root tree
//		RNode df = ROOT::RDataFrame("Events", names);
//	
//        //// enable multithreading
//        //ROOT::EnableImplicitMT();
//		
//
//        std::vector<std::string> colNames = df.GetColumnNames();
//		if (std::find(colNames.begin(), colNames.end(), "genWeight") == colNames.end()) {
//            df = df.Define("genWeight", [](){return (float)1.0;});
//		}
//        auto tunep_pt = df.Define("Muon_tP_pt" , []( const rvec<float>& muon_pt, rvec<float>& muon_relTP_pt){
//           return muon_relTP_pt * muon_pt;}, {"Muon_pt", "Muon_tunepRelPt"});
//            //return (float) 1. * muon_pt;}, {"Muon_pt", "Muon_tunepRelPt"});
//		auto loopcounter = tunep_pt.Filter([](ULong64_t e){if (0ULL == e) std::cout << "Running evtloop" << std::endl; return true;},{"rdfentry_"});
		
	
// read in files to run over
    json files;
    std::ifstream files_json;
    files_json.open(argv[1]);
    files_json >> files;

    // loop over all types of backgrounds
    std::string sample = files["sample"];
    std::string config = files["config"];
    std::string folder = files["folder"];
    std::vector <std::string > filenames = files["files"];
    std::vector <std::string > filepaths;

    for (auto& it : filenames){
        filepaths.push_back(folder + it);
    }

    std::ifstream cfg_json;
    cfg_json.open(config);
    cfg_json >> cfg;

    config::load_config_file(cfg);

    TString tmp = "output/" + std::to_string(config::era) + "/" + sample + ".root";

    std::cout << "\t \t Running on: " << sample << std::endl;

    runName = sample;

    // create outputfile to write hists
    TFile* outFile = new TFile(tmp, "RECREATE", "", 101);

    // open root tree
    RNode df = ROOT::RDataFrame("Events", filepaths);
    df = df.Define("Muon_tP_pt" , []( const rvec<float>& muon_pt, rvec<float>& muon_relTP_pt){
       return muon_relTP_pt * muon_pt;}, {"Muon_pt", "Muon_tunepRelPt"})
       //return muon_pt;}, {"Muon_pt", "Muon_tunepRelPt"})
           .Define("multiplicity", []( const unsigned int& nele, const unsigned int& nmuon, const unsigned int& njets){ return nele + nmuon + njets;}, {"nElectron","nMuon","nJet"});
    if(!config::runOnData){
        df = df.Define("Muon_pt_corr", [](const rvec<float>& muon_pt, const rvec<float>& muon_eta, const rvec<float>& muon_phi, const rvec<float>& muon_mass){return muon_reso_smearing(muon_pt,muon_eta,muon_phi,muon_mass,"");},{"Muon_tP_pt", "Muon_eta", "Muon_phi", "Muon_mass"});
    }
    else{
        df = df.Define("Muon_pt_corr", "Muon_tP_pt");
    }
    
    if (argc == 3 && std::string(argv[2]) == "test") {
        std::cout << "Running in test mode." << std::endl;
        df = df.Range(500000);
    }
    else if (argc == 4 && std::string(argv[2]) == "test") {
        std::cout << "Running in test mode with "<< argv[3]<< " events." << std::endl;
        const char* value = argv[3];
        stringstream strValue;
        strValue << value;

        unsigned int intValue;
        strValue >> intValue;
        df = df.Range(intValue);
        
    }

    std::vector<std::string> colNames = df.GetColumnNames();
    if (std::find(colNames.begin(), colNames.end(), "genWeight") == colNames.end()) {
        df = df.Define("genWeight", [](){return (float)1.0;});
    }
    if (!config::runOnData){
        if (std::find(colNames.begin(), colNames.end(), "LHE_HT") == colNames.end()){
            df = df.Define("LHE_HT", [](){return (float)0.0;}).Define("LHEPart_pt", [](){return (float)0.0;}).Define("LHEPart_eta", [](){return (float)0.0;}).Define("LHEPart_phi", [](){return (float)0.0;}).Define("LHEPart_mass", [](){return (float)0.0;}).Define("LHEPart_status", [](){return (float)0.0;}).Define("LHEPart_pdgId", [](){return (float)0.0;});
        }
    }
    auto loopcounter = df.Filter([](ULong64_t e){if (0ULL == e) std::cout << "Running evtloop" << std::endl; return true;},{"rdfentry_"});

    if (std::find(colNames.begin(), colNames.end(), "METFixEE2017_pt") != colNames.end()) {
        if (config::era == 2017) {
            met_branch_name = "METFixEE2017";
        }
        else if (config::use_EEMET){
            met_branch_name = "METFixEE2017";
        }
    }
    if(config::era == 2022){
        loopcounter = loopcounter.Define("MET_pt_nom",[](const float& met_pt){return met_pt;},{"MET_pt"})
                                 .Define("MET_phi_nom",[](const float& met_pt){return met_pt;},{"MET_phi"})
                                 .Define("Jet_eta_nom",[](const rvec<float>& jet_pt){return jet_pt;},{"Jet_eta"})
                                 .Define("Jet_phi_nom",[](const rvec<float>& jet_pt){return jet_pt;},{"Jet_phi"})
                                 .Define("Jet_pt_nom",[](const rvec<float>& jet_pt){return jet_pt;},{"Jet_pt"});
                                 //.Define("Electron_cutBased_HEEP",[](const rvec<float>& ele_pt){
                                 //       rvec<bool> heep;
                                 //       for(unsigned int i; i < ele_pt.size(); i++){
                                 //               heep.push_back(false);
                                 //           }
                                 //       return heep;
                                 //       },{"Electron_pt"});
    }
//std::vector < std::string > met_shifts = 
//{ 
//    "_jer",
//    "_jerUp",
//    "_jerDown",
//    "_jesTotalUp",
//    "_jesTotalDown",
//    "_unclustEnUp",
//    "_unclustEnDown"
//};
    auto corr_met_df = loopcounter;


    

    if (config::runOnData) {
        std::ifstream goldenjson_file;
        goldenjson_file.open(cfg["json_file"]);
        goldenjson_file >> goldenjson;

        auto jsoncleaned = corr_met_df.Filter(json_check, {"run", "luminosityBlock"}, "json cleaning");
    
        // this function does all analysis steps
        analyse(jsoncleaned, outFile);
    
    } else {
        // counter - PSWeight is filled with ones
        auto definecounts = corr_met_df.Define("abs_gen_weight", [](const float& x, const UInt_t& run, const UInt_t& lumi, const unsigned long long& event){
        
                                                return std::abs(x);
                                                }, {"genWeight","run","luminosityBlock","event"});
          if(config::runOnSignal){
                    definecounts = definecounts.Define("CollMass_gen_presel", [](const rvec<float>& genpt,
                                                                    const rvec<float>& geneta, 
                                                                    const rvec<float>& genphi, 
                                                                    const rvec<float>& genmass, 
                                                                    const rvec<int>& genmotheridx,
                                                                    const rvec<float>& genmuonpt,
                                                                    const rvec<float>& genmuoneta,
                                                                    const rvec<float>& genmuonphi,
                                                                    const rvec<float>& genmuonmass,
                                                                    const rvec<int>& genmuonpdgid,
                                                                    const rvec<int>& genmuonmotheridx,
                                                                    const rvec<int>& genmuonstatus,
                                                                    const float& met_pt,
                                                                    const float& met_phi){
                                                    TLorentzVector tau, muon, met, highest_tau, highest_muon;
                                                    float highest_mass = 0;
                                                    int best_i = 0;
                                                    int best_j = 0;
                                                    for (unsigned int i = 0; i < genmuonpt.size(); i++){
                                                            if( (abs(genmuonpdgid[i]) == 15 and abs(genmuonpdgid[genmuonmotheridx[i]]) == 32) and config::runOnZprime){
                                                                tau.SetPtEtaPhiM(genmuonpt[i], genmuoneta[i], genmuonphi[i], genmuonmass[i]);
                                                                double mass_vis = 0;
                                                                for(unsigned int j = 0; j < genmuonpt.size(); j++){
                                                                    if( (abs(genmuonpdgid[j]) == 13 and abs(genmuonpdgid[genmuonmotheridx[j]]) == 32) ){
                                                                        muon.SetPtEtaPhiM(genmuonpt[j], genmuoneta[j], genmuonphi[j], genmuonmass[j]);
                                                                        mass_vis = (tau + muon).M();
                                                                        if ( highest_mass < mass_vis){
                                                                            highest_mass = mass_vis;
                                                                            highest_tau = tau;
                                                                            highest_muon = muon;
                                                                            best_i = i;
                                                                            best_j = j;
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                            else if( abs(genmuonpdgid[i] == 15) and not config::runOnZprime){
                                                                tau.SetPtEtaPhiM(genmuonpt[i], genmuoneta[i], genmuonphi[i], genmuonmass[i]);
                                                                double mass_vis = 0;
                                                                for(unsigned int j = 0; j < genmuonpt.size(); j++){
                                                                    if( (abs(genmuonpdgid[j]) == 13 and genmuonmotheridx[j] == genmuonmotheridx[i]) and genmuonstatus[i] == genmuonstatus[j]){
                                                                        muon.SetPtEtaPhiM(genmuonpt[j], genmuoneta[j], genmuonphi[j], genmuonmass[j]);
                                                                        mass_vis = (tau + muon).M();
                                                                        if ( highest_mass < mass_vis){
                                                                            highest_mass = mass_vis;
                                                                            highest_tau = tau;
                                                                            highest_muon = muon;
                                                                            best_i = i;
                                                                            best_j = j;
                                                                        }
                                                                    }
                                                                }
                                                                
                                                            }
                                                    }
                                                    return highest_mass;
                                       
                                         },{"GenVisTau_pt","GenVisTau_eta","GenVisTau_phi","GenVisTau_mass","GenVisTau_genPartIdxMother","GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass","GenPart_pdgId","GenPart_genPartIdxMother","GenPart_statusFlags","GenMET_pt","GenMET_phi"});
                }
               // else{
               //     definecounts = definecounts.Define("CollMass_gen_presel", [](const rvec<int>& pdgid,
               //                                                          const rvec<float>& mass,
               //                                                          const rvec<float>& pt,
               //                                                          const rvec<float>& eta,
               //                                                          const rvec<float>& phi,
               //                                                          const rvec<int>& status
               //                                                     ){
               //                                     TLorentzVector tau, muon, met, highest_tau, highest_muon;
               //                                     float highest_mass = 0;
               //                                     int best_i = 0;
               //                                     int best_j = 0;
               //                                     for (unsigned int i = 0; i < pt.size(); i++){
               //                                             if( (abs(pdgid[i]) == 15) and status[i] == 1){
               //                                                 tau.SetPtEtaPhiM(pt[i], eta[i], phi[i], mass[i]);
               //                                                 double mass_vis = 0;
               //                                                 for(unsigned int j = 0; j < pt.size(); j++){
               //                                                     if( abs(pdgid[j]) == 13 and status[j] == 1){
               //                                                         muon.SetPtEtaPhiM(pt[j], eta[j], phi[j], mass[j]);
               //                                                         mass_vis = (tau + muon).M();
               //                                                         if ( highest_mass < mass_vis){
               //                                                             highest_mass = mass_vis;
               //                                                             highest_tau = tau;
               //                                                             highest_muon = muon;
               //                                                             best_i = i;
               //                                                             best_j = j;
               //                                                         }
               //                                                     }
               //                                                 }
               //                                             }
               //                                     }
               //                                     //for (unsigned int i = 0; i < genmuonpt.size(); i++){
               //                                     //        if( (abs(genmuonpdgid[i]) == 15 and abs(genmuonpdgid[genmuonmotheridx[i]]) == 32)){
               //                                     //            tau.SetPtEtaPhiM(genmuonpt[i], genmuoneta[i], genmuonphi[i], genmuonmass[i]);
               //                                     //            double mass_vis = 0;
               //                                     //            for(unsigned int j = 0; j < genmuonpt.size(); j++){
               //                                     //                if( (abs(genmuonpdgid[j]) == 13 and abs(genmuonpdgid[genmuonmotheridx[j]]) == 32) ){
               //                                     //                    muon.SetPtEtaPhiM(genmuonpt[j], genmuoneta[j], genmuonphi[j], genmuonmass[j]);
               //                                     //                    mass_vis = (tau + muon).M();
               //                                     //                    if ( highest_mass < mass_vis){
               //                                     //                        highest_mass = mass_vis;
               //                                     //                        highest_tau = tau;
               //                                     //                        highest_muon = muon;
               //                                     //                        best_i = i;
               //                                     //                        best_j = j;
               //                                     //                    }
               //                                     //                }
               //                                     //            }
               //                                     //        }
               //                                     //}
               //                                     return highest_mass;
               //                          },{"LHEPart_pdgId","LHEPart_mass","LHEPart_pt","LHEPart_eta","LHEPart_phi","LHEPart_status"});
               //                          
               //                          //},{"GenVisTau_pt","GenVisTau_eta","GenVisTau_phi","GenVisTau_mass","GenVisTau_genPartIdxMother","GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass","GenPart_pdgId","GenPart_genPartIdxMother","GenPart_statusFlags","GenMET_pt","GenMET_phi"});
               //     }
                
        else{
            definecounts = definecounts.Define("CollMass_gen_presel",[](){ return (float) 0.;},{});
        }
        auto collmass_gen_presel_hist = definecounts.Histo1D(	{"CollMass_gen_preselection", "", 			6000u, 0, 6000}, 					"CollMass_gen_presel", "genWeight");
        //auto top_pair_invmass = df.Histo1D(	{"top_pair_invmass", "", 			6000u, 0, 6000}, 					"top_pair_invmass");
        hist_dict.emplace("Preselection", collmass_gen_presel_hist);
        if(config::TT){
            definecounts = definecounts.Define("top_pt_weight", calc_top_pt_reweighting, {"GenPart_pdgId", "GenPart_pt"});
        }
        else{

            definecounts = definecounts.Define("top_pt_weight", [](){return (float) 1.;},{});
        }
        //if(config::TT ){
        //    //definecounts = definecounts.Define("Mll", lepton_inv_mass, {"GenPart_pdgId", "GenPart_mass", "GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_genPartIdxMother"});
        //}
        if(config::WW or config::TT ){
            definecounts = definecounts.Define("Mll", lepton_inv_mass_ww, {"LHEPart_pdgId", "LHEPart_mass","LHEPart_pt","LHEPart_eta","LHEPart_phi","LHEPart_status"});
        }
        else if(config::DY){
            definecounts = definecounts.Define("Mll", lepton_inv_mass_dy, {"LHEPart_pt","LHEPart_eta","LHEPart_phi","LHEPart_mass","LHEPart_pdgId"});
            
        }
        else{
            
            definecounts = definecounts.Define("Mll",[](){ return (float) 0.;},{});
        }
                
        // clean gen files
        RNode gencleaned = definecounts;
        if(config::cut_type != ""){
            gencleaned = definecounts.Filter(clean_gen_file, {"LHE_HT", "LHEPart_pt", "LHEPart_eta", "LHEPart_phi", "LHEPart_mass", "LHEPart_pdgId", "LHEPart_status"}, "gen cleaning");
        }
        //auto lepton_pair_idx = gencleaned.Define("lepton_pair_idx",lepton_pair_idx, {"GenPart_pdgId","GenPart_genPartIdxMother"});
        //auto lepton_pair_mass = lepton_pair_idx.Define("lepton_pair_invmass", mass_inv_idx, {"GenPart_pt", "GenPart_eta","GenPart_phi","GenPart_mass","top_pair_idx"});

        
        
        // this function does all analysis steps
        analyse(gencleaned, outFile);
    };

    
    config::clean_memory();
    hist_dict.clear();
    hist_dict_2d.clear();
    delete outFile;
    
    //	// disable multithreading
    //	ROOT::DisableImplicitMT();
	return 0;
}
