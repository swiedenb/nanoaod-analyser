#include "PhysicalQuantities.hh"
#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"
#include <math.h>
#include "Math/LorentzVector.h"

// Calculate sphericity from one particle class
float sphericity(	const rvec<float>& part_pt, 
					const rvec<float>& part_eta, 
					const rvec<float>& part_phi, 
					const rvec<float>& part_mass,
                    const rvec<bool>& mask) {
    TMatrixD h(2,2);
    float sum_px = 0.;
    float sum_py = 0.;
    float sum_pxpy = 0.;
    float sum_pypx = 0.;
    float sum_pt = 0;
    for ( int i = 0; i < part_pt[mask].size(); i++){
        TLorentzVector part;
        part.SetPtEtaPhiM(part_pt[mask][i],part_eta[mask][i],part_phi[mask][i],part_mass[mask][i]);
        sum_px += part.Px()*part.Px()/part.Pt();
        sum_py += part.Py()*part.Py()/part.Pt();
        sum_pxpy += part.Px()*part.Py()/part.Pt();
        sum_pypx += part.Py()*part.Px()/part.Pt();
        sum_pt += part.Pt();
    }
    h[0][0] = sum_px/sum_pt;
    h[0][1] = sum_pxpy/sum_pt;
    h[1][0] = sum_pypx/sum_pt;
    h[1][1] = sum_py/sum_pt;
    const TMatrixDEigen eigen(h);
    const TMatrixD eigenVal = eigen.GetEigenValues();
    float St = (2*eigenVal[1][1])/(eigenVal[1][1]+eigenVal[0][0]);
	return St;
};
float sphericity_leptons(	const rvec<float>& part_pt, 
					const rvec<float>& part_eta, 
					const rvec<float>& part_phi, 
					const rvec<float>& part_mass,
					const rvec<float>& part2_pt, 
					const rvec<float>& part2_eta, 
					const rvec<float>& part2_phi, 
					const rvec<float>& part2_mass,
					const rvec<bool>& mask,
                    const rvec<bool>& mask2) {
    TMatrixD h(2,2);
    float sum_px = 0.;
    float sum_py = 0.;
    float sum_pxpy = 0.;
    float sum_pypx = 0.;
    float sum_pt = 0;
    for ( int i = 0; i < part_pt[mask].size(); i++){
        TLorentzVector part;
        part.SetPtEtaPhiM(part_pt[mask][i],part_eta[mask][i],part_phi[mask][i],part_mass[mask][i]);
        sum_px += part.Px()*part.Px()/part.Pt();
        sum_py += part.Py()*part.Py()/part.Pt();
        sum_pxpy += part.Px()*part.Py()/part.Pt();
        sum_pypx += part.Py()*part.Px()/part.Pt();
        sum_pt += part.Pt();
    }
    for ( int i = 0; i < part2_pt[mask2].size(); i++){
        TLorentzVector part;
        part.SetPtEtaPhiM(part2_pt[mask2][i],part2_eta[mask2][i],part2_phi[mask2][i],part2_mass[mask2][i]);
        sum_px += part.Px()*part.Px()/part.Pt();
        sum_py += part.Py()*part.Py()/part.Pt();
        sum_pxpy += part.Px()*part.Py()/part.Pt();
        sum_pypx += part.Py()*part.Px()/part.Pt();
        sum_pt += part.Pt();
    }
    h[0][0] = sum_px/sum_pt;
    h[0][1] = sum_pxpy/sum_pt;
    h[1][0] = sum_pypx/sum_pt;
    h[1][1] = sum_py/sum_pt;
    const TMatrixDEigen eigen(h);
    const TMatrixD eigenVal = eigen.GetEigenValues();
    float St = (2*eigenVal[1][1])/(eigenVal[1][1]+eigenVal[0][0]);
	return St;
};
float sphericity_full(	const rvec<float>& part1_pt, 
					const rvec<float>& part1_eta, 
					const rvec<float>& part1_phi, 
					const rvec<float>& part1_mass,
                    const rvec<float>& part2_pt, 
					const rvec<float>& part2_eta, 
					const rvec<float>& part2_phi, 
					const rvec<float>& part2_mass,
                    const rvec<float>& part3_pt, 
					const rvec<float>& part3_eta, 
					const rvec<float>& part3_phi, 
					const rvec<float>& part3_mass,
                    const rvec<bool>&  mask,
                    const rvec<bool>&  mask2) {
    TMatrixD h(2,2);
    float sum_px = 0.;
    float sum_py = 0.;
    float sum_pxpy = 0.;
    float sum_pypx = 0.;
    float sum_pt = 0;
    for ( int i = 0; i < part1_pt[mask].size(); i++){
        TLorentzVector part;
        part.SetPtEtaPhiM(part1_pt[mask][i],part1_eta[mask][i],part1_phi[mask][i],part1_mass[mask][i]);
        sum_px += part.Px()*part.Px()/part.Pt();
        sum_py += part.Py()*part.Py()/part.Pt();
        sum_pxpy += part.Px()*part.Py()/part.Pt();
        sum_pypx += part.Py()*part.Px()/part.Pt();
        sum_pt += part.Pt();
    }
    for ( int i = 0; i < part2_pt[mask2].size(); i++){
        TLorentzVector part;
        part.SetPtEtaPhiM(part2_pt[mask2][i],part2_eta[mask2][i],part2_phi[mask2][i],part2_mass[mask2][i]);
        sum_px += part.Px()*part.Px()/part.Pt();
        sum_py += part.Py()*part.Py()/part.Pt();
        sum_pxpy += part.Px()*part.Py()/part.Pt();
        sum_pypx += part.Py()*part.Px()/part.Pt();
        sum_pt += part.Pt();
    }
    for ( int i = 0; i < part3_pt.size(); i++){
        TLorentzVector part;
        part.SetPtEtaPhiM(part3_pt[i],part3_eta[i],part3_phi[i],part3_mass[i]);
        sum_px += part.Px()*part.Px()/part.Pt();
        sum_py += part.Py()*part.Py()/part.Pt();
        sum_pxpy += part.Px()*part.Py()/part.Pt();
        sum_pypx += part.Py()*part.Px()/part.Pt();
        sum_pt += part.Pt();
    }
    h[0][0] = sum_px/sum_pt;
    h[0][1] = sum_pxpy/sum_pt;
    h[1][0] = sum_pypx/sum_pt;
    h[1][1] = sum_py/sum_pt;
    const TMatrixDEigen eigen(h);
    const TMatrixD eigenVal = eigen.GetEigenValues();
    float St = (2*eigenVal[1][1])/(eigenVal[1][1]+eigenVal[0][0]);
	return St;
};
// Calculate MT from one particle and MET
float mass_transv(	const float& part_pt, 
					const float& part_phi, 
					const float& met_pt, 
					const float& met_phi) {
	return sqrt(2*part_pt*met_pt*(1 - cos( delta_phi(part_phi, met_phi) ) ) );
};

// Calculate MT from first (highest pt) particle, which fulfils mask, -  and MET
float mass_transv_masked(	const rvec<float>& part_pt, 
							const rvec<float>& part_phi, 
							const rvec<bool>& part_mask, 
							const float& met_pt, 
							const float& met_phi) {
	return sqrt(2*part_pt[part_mask][0]*met_pt*(1 - cos( delta_phi(part_phi[part_mask][0], met_phi) ) ) );
};
float lepton_inv_mass_dy(const rvec<float>& lhept,
                         const rvec<float>& lheeta,
                         const rvec<float>& lhephi,
                         const rvec<float>& lhemass,
                         const rvec<int>& lhepdgid){
           
                TLorentzVector l1,l2,ll;
                bool found_1 = false;
                bool found_2 = false;
                int l1_pdg;
                int l2_pdg;
                for (uint i = 0; i < lhepdgid.size(); i++){
                    if(lhepdgid[i] == 11 || lhepdgid[i] == 13 || lhepdgid[i] == 15){
                        l1_pdg = lhepdgid[i];
                        l1.SetPtEtaPhiM(lhept[i],lheeta[i],lhephi[i],lhemass[i]);
                        found_1 = true;
                    }
                    else if(lhepdgid[i] == -11 || lhepdgid[i] == -13 || lhepdgid[i] == -15){
                        l2_pdg = lhepdgid[i];
                        l2.SetPtEtaPhiM(lhept[i],lheeta[i],lhephi[i],lhemass[i]);
                        found_2 = true;
                        
                    }
                    if( l1_pdg == -l2_pdg && found_1 && found_2){
                        auto mass = (l1 + l2).M();
                        return mass;
                    }
                }
                return 0.;
          }

float lepton_inv_mass_ww(const rvec<int>& pdgID,
                      const rvec<float>& mass,
                      const rvec<float>& pt,
                      const rvec<float>& eta,
                      const rvec<float>& phi,
                      const rvec<int>& status){
    
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
                       if(status[idx_1] == 1){
                           first = true;
                       }
                    }
                    else if (not second)
                    {
                       idx_2 = std::distance(pdgID.begin(), iter++);
                       if(status[idx_2] == 1){
                           second = true;
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
           if (first and second){
               TLorentzVector part1;
               TLorentzVector part2;
               part1.SetPtEtaPhiM(pt[idx_1],eta[idx_1],phi[idx_1],mass[idx_1]);
               part2.SetPtEtaPhiM(pt[idx_2],eta[idx_2],phi[idx_2],mass[idx_2]);
               auto inv_mass = (part1 + part2).M();
               return inv_mass;
           }
           return 0.;
};
float lepton_inv_mass(const rvec<int>& pdgID,
                      const rvec<float>& mass,
                      const rvec<float>& pt,
                      const rvec<float>& eta,
                      const rvec<float>& phi,
                      const rvec<int>& mother_idx){
    
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
                       }
                    }
                    else if (not second)
                    {
                       idx_2 = std::distance(pdgID.begin(), iter++);
                       if (abs(pdgID[mother_idx[idx_2]]) == 24){
                           second = true;
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
           if (first and second){
               TLorentzVector part1;
               TLorentzVector part2;
               part1.SetPtEtaPhiM(pt[idx_1],eta[idx_1],phi[idx_1],mass[idx_1]);
               part2.SetPtEtaPhiM(pt[idx_2],eta[idx_2],phi[idx_2],mass[idx_2]);
               auto inv_mass = (part1 + part2).M();
               return inv_mass;
           }
           return 0.;
};
// Calculate invariant mass from two particles with idx
float mass_inv(	const float& part1_pt,
				const float& part2_pt,
				const float& part1_eta,
				const float& part2_eta,
				const float& part1_phi,
				const float& part2_phi,
				const float& part1_mass,
				const float& part2_mass) {
    TLorentzVector p1, p2;
    p1.SetPtEtaPhiM(part1_pt, part1_eta, part1_phi, part1_mass);
    p2.SetPtEtaPhiM(part2_pt, part2_eta, part2_phi, part2_mass);
    return (p1 + p2).M();
};

// Calculate invariant mass for to gen particles
float mass_inv_idx(	const rvec<float>& pt,
						const rvec<float>& eta,
						const rvec<float>& phi,
						const rvec<float>& mass,
						const rvec<int>& idx) {
    TLorentzVector p1, p2;
    p1.SetPtEtaPhiM(pt[idx[0]], eta[idx[0]], phi[idx[0]], mass[idx[0]]);
    p2.SetPtEtaPhiM(pt[idx[1]], eta[idx[1]], phi[idx[1]], mass[idx[1]]);
    return (p1 + p2).M();
};
// Calculate invariant mass from two highest pt particles, which fulfil mask
float mass_inv_masked(	const rvec<float>& p1_pt,
						const rvec<float>& p2_pt,
						const rvec<float>& p1_eta,
						const rvec<float>& p2_eta,
						const rvec<float>& p1_phi,
						const rvec<float>& p2_phi,
						const rvec<float>& p1_mass,
						const rvec<float>& p2_mass,
						const rvec<bool>& p1_mask,
						const rvec<bool>& p2_mask) {
    TLorentzVector p1, p2;
    p1.SetPtEtaPhiM(p1_pt[p1_mask][0], p1_eta[p1_mask][0], p1_phi[p1_mask][0], p1_mass[p1_mask][0]);
    p2.SetPtEtaPhiM(p2_pt[p2_mask][0], p2_eta[p2_mask][0], p2_phi[p2_mask][0], p2_mass[p2_mask][0]);
    return (p1 + p2).M();
};

float sel_tau_pt(               const rvec<float>& tau_pt,
                                const rvec<bool>& tau_mask,    
                                const rvec<int>& col_idx){
    return tau_pt[tau_mask][col_idx[0]];
};
// Calculate collinear mass from two highest pt particles, which fulfil mask and MET with the eta from the first particle
float collinear_mass_mask(      const rvec<float>& p1_pt_vec,
                                const rvec<float>& p2_pt_vec,
                                const rvec<float>& p1_eta_vec,
                                const rvec<float>& p2_eta_vec,
                                const rvec<float>& p1_phi_vec,
                                const rvec<float>& p2_phi_vec,
                                const rvec<float>& p1_mass_vec,
                                const rvec<float>& p2_mass_vec,
                                const float& met_pt,
                                const float& met_phi,
                                const rvec<bool>& mask_1,
                                const rvec<bool>& mask_2,
                                const rvec<int>& idx){
    TLorentzVector p1, p2, p3, p_best;
    auto p1_pt = p1_pt_vec[mask_1][idx[0]];
    auto p1_eta = p1_eta_vec[mask_1][idx[0]];
    auto p1_phi = p1_phi_vec[mask_1][idx[0]];
    auto p1_mass = p1_mass_vec[mask_1][idx[0]];
    auto p2_pt = p2_pt_vec[mask_2][idx[1]];
    auto p2_eta = p2_eta_vec[mask_2][idx[1]];
    auto p2_phi = p2_phi_vec[mask_2][idx[1]];
    auto p2_mass = p2_mass_vec[mask_2][idx[1]];
    
    p1.SetPtEtaPhiM(p1_pt, p1_eta, p1_phi, p1_mass);
    p2.SetPtEtaPhiM(p2_pt, p2_eta, p2_phi, p2_mass);
    p3.SetPtEtaPhiM(met_pt, 0., met_phi, 0.);
	double METproj=(p3.Px()*p1.Px()+p3.Py()*p1.Py())/p1.Pt();
	//if (METproj < 0) METproj = 0;
	double xth=1;
	if(METproj>0) xth=p1.Pt()/(p1.Pt()+METproj);
	else xth = 1;
	double mass_vis=(p1+p2).M();
	double mcol = 0;
	if (mass_vis != mass_vis) mass_vis=0;
	if (mass_vis <= 0) mass_vis = 0;

	mcol=mass_vis/sqrt(xth);
	return mcol;
};
float collinear_mass(           const float& p1_pt,
                                const float& p2_pt,
                                const float& p1_eta,
                                const float& p2_eta,
                                const float& p1_phi,
                                const float& p2_phi,
                                const float& p1_mass,
                                const float& p2_mass,
                                const float& met_pt,
                                const float& met_phi){
    TLorentzVector p1, p2, p3, p_best;
    p1.SetPtEtaPhiM(p1_pt, p1_eta, p1_phi, p1_mass);
    p2.SetPtEtaPhiM(p2_pt, p2_eta, p2_phi, p2_mass);
    p3.SetPtEtaPhiM(met_pt, 0., met_phi, 0.);
	double METproj=(p3.Px()*p1.Px()+p3.Py()*p1.Py())/p1.Pt();
	//if (METproj < 0) METproj = 0;
	double xth=1;
	if(METproj>0) xth=p1.Pt()/(p1.Pt()+METproj);
	else xth = 1;
	double mass_vis=(p1+p2).M();
	double mcol = 0;
	if (mass_vis != mass_vis) mass_vis=0;
	if (mass_vis <= 0) mass_vis = 0;

	mcol=mass_vis/sqrt(xth);

	return mcol;
};
//float collinear_mass(           const rvec<int>& col_idx,
//                                const rvec<float>& p1_pt,
//                                const rvec<float>& p2_pt,
//                                const rvec<float>& p1_eta,
//                                const rvec<float>& p2_eta,
//                                const rvec<float>& p1_phi,
//                                const rvec<float>& p2_phi,
//                                const rvec<float>& p1_mass,
//                                const rvec<float>& p2_mass,
//                                const float& met_pt,
//                                const float& met_phi,
//                                const rvec<bool>& p1_mask,
//                                const rvec<bool>& p2_mask){
//    TLorentzVector p1, p2, p3, p_best;
//    p1.SetPtEtaPhiM(p1_pt[p1_mask][col_idx[0]], p1_eta[p1_mask][col_idx[0]], p1_phi[p1_mask][col_idx[0]], p1_mass[p1_mask][col_idx[0]]);
//    p2.SetPtEtaPhiM(p2_pt[p2_mask][col_idx[1]], p2_eta[p2_mask][col_idx[1]], p2_phi[p2_mask][col_idx[1]], p2_mass[p2_mask][col_idx[1]]);
//    p3.SetPtEtaPhiM(met_pt, p1_eta[p1_mask][col_idx[0]], met_phi, 0.);
//      
//    return (p1 + p2 + p3).M(); 
//};
float tau_pt_over_jet_pt(       const rvec<int>& col_idx,
                                const rvec<float>& tau_pt,
                                const rvec<int>& jet_idx,
                                const rvec<float>& jet_pt,
                                const rvec<bool>& tau_mask){
    return tau_pt[tau_mask][col_idx[0]]/jet_pt[jet_idx[col_idx[0]]];
};
float tau_pt_over_jet_pt_barrel(       const rvec<int>& col_idx,
                                       const rvec<float>& tau_pt,
                                       const rvec<int>& jet_idx,
                                       const rvec<float>& jet_pt,
                                       const rvec<bool>& tau_mask,
                                       const rvec<float>& tau_eta){
    if (abs(tau_eta[tau_mask][col_idx[0]]) < 1.46)  return tau_pt[tau_mask][col_idx[0]]/jet_pt[jet_idx[col_idx[0]]];
};
float tau_pt_over_jet_pt_endcap(       const rvec<int>& col_idx,
                                       const rvec<float>& tau_pt,
                                       const rvec<int>& jet_idx,
                                       const rvec<float>& jet_pt,
                                       const rvec<bool>& tau_mask,
                                       const rvec<float>& tau_eta){
    if (abs(tau_eta[tau_mask][col_idx[0]]) > 1.56)  return tau_pt[tau_mask][col_idx[0]]/jet_pt[jet_idx[col_idx[0]]];
};
// Calculate collinear mass from two highest pt particles, which fulfil mask and MET is projected onto the tau pT vector to calculate a rescale ratio for the tau 4 momentum
float collinear_mass_alt_2(       const float& p1_pt,
                                const float& p2_pt,
                                const float& p1_eta,
                                const float& p2_eta,
                                const float& p1_phi,
                                const float& p2_phi,
                                const float& p1_mass,
                                const float& p2_mass,
                                const float& met_pt,
                                const float& met_phi){
    TLorentzVector p_tau, p2, p_met, p_tau2,p_tau3;
    p_tau.SetPtEtaPhiM(p1_pt, p1_eta, p1_phi, p1_mass);
    p2.SetPtEtaPhiM(p2_pt, p2_eta, p2_phi, p2_mass);
    p_met.SetPtEtaPhiM(met_pt, 0. , met_phi, 0.);
    double x = 0;
    double p_met_proj;
    double mass_vis = 0;
    double mass_coll = 4000.;
    //p_met_proj = std::max((p_met.Px()*p_tau.Px() + p_met.Py()*p_tau.Py())/p_tau.Pt(),0.);
    p_met_proj = (p_met.Px()*p_tau.Px() + p_met.Py()*p_tau.Py())/p_tau.Pt();
    x = p_tau.Pt() / (p_tau.Pt() + p_met_proj); 
    
    mass_vis = (p_tau + p2).M();
    double mass_coll_alt = mass_vis/sqrt(x);
    p_tau3.SetPxPyPzE(p_tau.Px()/x,p_tau.Py()/x,p_tau.Pz(),p_tau.E()/x);
    p_tau2.SetPxPyPzE(p_tau.Px()/x,p_tau.Py()/x,p_tau.Pz()/x,p_tau.E()/x);
    mass_coll = (p_tau2 + p2).M();
    double mass_coll_alt2 = (p_tau3 + p2).M();
    return mass_coll_alt2;
};
float collinear_mass_alt(       const float& p1_pt,
                                const float& p2_pt,
                                const float& p1_eta,
                                const float& p2_eta,
                                const float& p1_phi,
                                const float& p2_phi,
                                const float& p1_mass,
                                const float& p2_mass,
                                const float& met_pt,
                                const float& met_phi){
    TLorentzVector p_tau, p2, p_met, p_tau2,p_tau3;
    p_tau.SetPtEtaPhiM(p1_pt, p1_eta, p1_phi, p1_mass);
    p2.SetPtEtaPhiM(p2_pt, p2_eta, p2_phi, p2_mass);
    //p_met.SetPtEtaPhiM(met_pt, 0., met_phi, 0.);
    p_met.SetPtEtaPhiM(met_pt, p1_eta, met_phi, 0.);
    double x = 0;
    double p_met_proj;
    double mass_vis = 0;
    double mass_coll = 4000.;
    p_met_proj = std::max((p_met.Px()*p_tau.Px() + p_met.Py()*p_tau.Py())/p_tau.Pt(),0.);
    //p_met_proj = (p_met.Px()*p_tau.Px() + p_met.Py()*p_tau.Py())/p_tau.Pt();
    x = p_tau.Pt() / (p_tau.Pt() + p_met_proj); 
    
    mass_vis = (p_tau + p2).M();
    if (mass_vis < 20.){
        std::cout<<"mass vis: "<<mass_vis<<std::endl;
        std::cout<<"x > 1 "<<x<<std::endl;
        std::cout<<"proj met: "<<p_met_proj<<std::endl;
        std::cout<<"tau pt: "<<p_tau.Pt()<<std::endl;
        std::cout<<"mass vis: "<<mass_vis<<std::endl;
        std::cout<<"MET: " << p_met.Pt()<<std::endl;
        std::cout<<"MET phi: " << met_phi<<std::endl;
        std::cout<<"tau px: "<< p_tau.Px()<<std::endl;
        std::cout<<"tau py: "<< p_tau.Py()<<std::endl;
        std::cout<<"tau pz: "<< p_tau.Pz()<<std::endl;
        std::cout<<"tau E: "<< p_tau.E()<<std::endl;
        std::cout<<"muon px: "<< p2.Px()<<std::endl;
        std::cout<<"muon py: "<< p2.Py()<<std::endl;
        std::cout<<"muon pz: "<< p2.Pz()<<std::endl;
        std::cout<<"muon E: "<< p2.E()<<std::endl;
    }
    double mass_coll_alt = mass_vis/sqrt(x);
    p_tau3.SetPxPyPzE(p_tau.Px()/x,p_tau.Py()/x,p_tau.Pz(),p_tau.E()/x);
    p_tau2.SetPxPyPzE(p_tau.Px()/x,p_tau.Py()/x,p_tau.Pz()/x,p_tau.E()/x);
    mass_coll = (p_tau2 + p2).M();
    double mass_coll_alt2 = (p_tau3 + p2).M();
    //std::cout<<mass_coll-mass_coll_alt<<std::endl;
    //std::cout<<"sqrt x: "<<sqrt(x)<<std::endl;
    //std::cout<<"mass coll tau weighted: "<<mass_coll<<std::endl;
    //std::cout<<"mass coll mass weighted: "<<mass_coll_alt<<std::endl;
    //std::cout<<"mass coll tau pt weighted: "<<mass_coll_alt2<<std::endl;
    if (mass_coll_alt < 20.){
        std::cout<<"                                            x: "<<x<<std::endl;
        std::cout<<"tauPt: "<<p_tau.Pt()<<std::endl;
        std::cout<<"tauPx: "<<p_tau.Px()<<std::endl;
        std::cout<<"tauPy: "<<p_tau.Py()<<std::endl;
        std::cout<<"tauPz: "<<p_tau.Pz()<<std::endl;
        std::cout<<"tauE: "<<p_tau.E()<<std::endl;
        std::cout<<"muonPt: "<<p2.Pt()<<std::endl;
        std::cout<<"muonPx: "<<p2.Px()<<std::endl;
        std::cout<<"muonPy: "<<p2.Py()<<std::endl;
        std::cout<<"muonPz: "<<p2.Pz()<<std::endl;
        std::cout<<"muonE: "<<p2.E()<<std::endl;
        std::cout<<"muon phi: "<<p2.Phi()<<std::endl;
        std::cout<<"met: "<<met_pt<<"    "<< p_met.Pt()<<std::endl;
        std::cout<<"mass: "<<mass_vis<<std::endl;
        std::cout<<"mass coll: "<<mass_coll<<std::endl;
        std::cout<<"mass coll_alt: "<<mass_coll_alt<<std::endl;
        std::cout<<"met phi: "<<met_phi<<std::endl;
        std::cout<<"tau phi: "<<p_tau.Phi()<<std::endl;
        std::cout<<"projected met: "<<p_met_proj<<std::endl;
    }


    //if (abs(delta_phi(p_tau.Phi(),met_phi)) > M_PI/2.){
    //    std::cout<<"sum: "<<(p_tau.Pt() + met_pt*cos(M_PI - abs(delta_phi(p_tau.Phi(),met_phi))))<<std::endl;
    //    std::cout<<"tau pt: "<<p_tau.Pt()<<" met: "<<met_pt << " tau_phi: "<<p_tau.Phi() << " met phi: "<<met_phi<<std::endl;
    //    x = p_tau.Pt() / (p_tau.Pt() - met_pt*cos(M_PI - abs(delta_phi(p_tau.Phi(),met_phi)))); 
    //    std::cout<< "x: "<< x<< std::endl;

    //}
    //else{
    //    x = p_tau.Pt() / (p_tau.Pt() + met_pt*cos(delta_phi(p_tau.Phi(),met_phi))); 
    //}

    //p_tau.SetPxPyPzE(p_tau.Px()/x,p_tau.Py()/x,p_tau.Pz()/x,p_tau.E()/x);
    //if( (p_tau + p2).M()< 5.){
    //    std::cout<< "px: "<< p_tau.Px() << " py: " << p_tau.Py() << " pz: "<<p_tau.Pz() << " E: " << p_tau.E() << std::endl;
    //    std::cout<< "M: "<<(p_tau + p2).M()<<std::endl;
    //    std::cout<< "x: "<< x<< std::endl;
    //}
    //return (p_tau + p2).M(); 
    return mass_coll_alt; 
};
//float collinear_mass_alt(       const rvec<int>& col_idx,
//                                const rvec<float>& p1_pt,
//                                const rvec<float>& p2_pt,
//                                const rvec<float>& p1_eta,
//                                const rvec<float>& p2_eta,
//                                const rvec<float>& p1_phi,
//                                const rvec<float>& p2_phi,
//                                const rvec<float>& p1_mass,
//                                const rvec<float>& p2_mass,
//                                const float& met_pt,
//                                const float& met_phi,
//                                const rvec<bool>& p1_mask,
//                                const rvec<bool>& p2_mask){
//    TLorentzVector p_tau, p2, p_met, p_tau2,p_tau3;
//    p_tau.SetPtEtaPhiM(p1_pt[p1_mask][col_idx[0]], p1_eta[p1_mask][col_idx[0]], p1_phi[p1_mask][col_idx[0]], p1_mass[p1_mask][col_idx[0]]);
//    p2.SetPtEtaPhiM(p2_pt[p2_mask][col_idx[1]], p2_eta[p2_mask][col_idx[1]], p2_phi[p2_mask][col_idx[1]], p2_mass[p2_mask][col_idx[1]]);
//    p_met.SetPtEtaPhiM(met_pt, 0, met_phi, 0.);
//    double x = 0;
//    double p_met_proj;
//    double mass_vis = 0;
//    double mass_coll = 4000.;
//    p_met_proj = (p_met.Px()*p_tau.Px() + p_met.Py()*p_tau.Py())/p_tau.Pt();
//    x = p_tau.Pt() / (p_tau.Pt() + p_met_proj); 
//    mass_vis = (p_tau + p2).M();
//    double mass_coll_alt = mass_vis/sqrt(x);
//    p_tau3.SetPxPyPzE(p_tau.Px()/x,p_tau.Py()/x,p_tau.Pz(),p_tau.E()/x);
//    p_tau2.SetPxPyPzE(p_tau.Px()/x,p_tau.Py()/x,p_tau.Pz()/x,p_tau.E()/x);
//    mass_coll = (p_tau2 + p2).M();
//    double mass_coll_alt2 = (p_tau3 + p2).M();
//    //std::cout<<mass_coll-mass_coll_alt<<std::endl;
//    //std::cout<<"sqrt x: "<<sqrt(x)<<std::endl;
//    //std::cout<<"mass coll tau weighted: "<<mass_coll<<std::endl;
//    //std::cout<<"mass coll mass weighted: "<<mass_coll_alt<<std::endl;
//    //std::cout<<"mass coll tau pt weighted: "<<mass_coll_alt2<<std::endl;
//    //if (mass_coll < 5.){
//    //    std::cout<<"                                            x: "<<x<<std::endl;
//    //    std::cout<<"sumPt: "<<p_sum.Pt()<<std::endl;
//    //    std::cout<<"sumPx: "<<p_sum.Px()<<std::endl;
//    //    std::cout<<"sumPy: "<<p_sum.Py()<<std::endl;
//    //    std::cout<<"sumPz: "<<p_sum.Pz()<<std::endl;
//    //    std::cout<<"sumE: "<<p_sum.E()<<std::endl;
//    //    std::cout<<"tauPt: "<<p_tau.Pt()<<std::endl;
//    //    std::cout<<"tauPx: "<<p_tau.Px()<<std::endl;
//    //    std::cout<<"tauPy: "<<p_tau.Py()<<std::endl;
//    //    std::cout<<"tauPz: "<<p_tau.Pz()<<std::endl;
//    //    std::cout<<"tauE: "<<p_tau.E()<<std::endl;
//    //    std::cout<<"muonPt: "<<p2.Pt()<<std::endl;
//    //    std::cout<<"muonPx: "<<p2.Px()<<std::endl;
//    //    std::cout<<"muonPy: "<<p2.Py()<<std::endl;
//    //    std::cout<<"muonPz: "<<p2.Pz()<<std::endl;
//    //    std::cout<<"muonE: "<<p2.E()<<std::endl;
//    //    std::cout<<"muon phi: "<<p2.Phi()<<std::endl;
//    //    std::cout<<"met: "<<met_pt<<"    "<< p_met.Pt()<<std::endl;
//    //    std::cout<<"mass: "<<mass_vis<<std::endl;
//    //    std::cout<<"mass coll: "<<mass_coll<<std::endl;
//    //    std::cout<<"met phi: "<<met_phi<<std::endl;
//    //    std::cout<<"tau phi: "<<p_tau.Phi()<<std::endl;
//    //    std::cout<<"projected met: "<<p_met_proj<<std::endl;
//    //}
//
//
//    //if (abs(delta_phi(p_tau.Phi(),met_phi)) > M_PI/2.){
//    //    std::cout<<"sum: "<<(p_tau.Pt() + met_pt*cos(M_PI - abs(delta_phi(p_tau.Phi(),met_phi))))<<std::endl;
//    //    std::cout<<"tau pt: "<<p_tau.Pt()<<" met: "<<met_pt << " tau_phi: "<<p_tau.Phi() << " met phi: "<<met_phi<<std::endl;
//    //    x = p_tau.Pt() / (p_tau.Pt() - met_pt*cos(M_PI - abs(delta_phi(p_tau.Phi(),met_phi)))); 
//    //    std::cout<< "x: "<< x<< std::endl;
//
//    //}
//    //else{
//    //    x = p_tau.Pt() / (p_tau.Pt() + met_pt*cos(delta_phi(p_tau.Phi(),met_phi))); 
//    //}
//
//    //p_tau.SetPxPyPzE(p_tau.Px()/x,p_tau.Py()/x,p_tau.Pz()/x,p_tau.E()/x);
//    //if( (p_tau + p2).M()< 5.){
//    //    std::cout<< "px: "<< p_tau.Px() << " py: " << p_tau.Py() << " pz: "<<p_tau.Pz() << " E: " << p_tau.E() << std::endl;
//    //    std::cout<< "M: "<<(p_tau + p2).M()<<std::endl;
//    //    std::cout<< "x: "<< x<< std::endl;
//    //}
//    //return (p_tau + p2).M(); 
//    return mass_coll_alt2; 
//};
rvec<int> col_idx(                    const rvec<float>& p1_pt,
                                      const rvec<float>& p2_pt,
                                      const rvec<float>& p1_eta,
                                      const rvec<float>& p2_eta,
                                      const rvec<float>& p1_phi,
                                      const rvec<float>& p2_phi,
                                      const rvec<float>& p1_mass,
                                      const rvec<float>& p2_mass,
                                      const float& met_pt,
                                      const float& met_phi,
                                      const rvec<bool>& p1_mask,
                                      const rvec<bool>& p2_mask) {
    TLorentzVector p1, p2, p_met;
    rvec<int> idx(2);
    float highest_mass = 0;
    int best_i = 0;
    int best_j = 0;
    for (unsigned int i = 0; i < p1_pt[p1_mask].size(); i++){
        p1.SetPtEtaPhiM(p1_pt[p1_mask][i], p1_eta[p1_mask][i], p1_phi[p1_mask][i], p1_mass[p1_mask][i]);
        //double x = 0;
        //double p_met_proj;
        double mass_vis = 0;
        //double mass_coll = 0;
        //p_met_proj = (p_met.Px()*p1.Px() + p_met.Py()*p1.Py())/p1.Pt();
        //x = p1.Pt() / (p1.Pt() + p_met_proj); 
        for(unsigned int j = 0; j < p2_pt[p2_mask].size(); j++){
            p2.SetPtEtaPhiM(p2_pt[p2_mask][j], p2_eta[p2_mask][j], p2_phi[p2_mask][j], p2_mass[p2_mask][j]);
            mass_vis = (p1 + p2).M();
            //mass_coll = mass_vis/sqrt(x);
            if ( highest_mass < mass_vis){
                highest_mass = mass_vis;
                best_i = i;
                best_j = j;
            }
        }
    }
    idx[0] = best_i;
    idx[1] = best_j;
    return idx; 
};

// Returns ratio of two given quants
float ratio(	const float& part1_qt, 
				const float& part2_qt) {
	return part1_qt / part2_qt;
};

// Returns ratio of all valeus in vector divided by quant
rvec <float> ratio_vector(	const rvec<float> & part1_qt, 
							const float& part2_qt) {
	rvec<float> result_vec(part1_qt.size());
	for (uint i = 0; i < part1_qt.size(); i++) {
		result_vec[i] = part1_qt[i] / part2_qt;
	}
	return result_vec;
};

// Returns Delta Phi between two particles
float delta_phi(	const float& part1_phi,
					const float& part2_phi) {
	double dphi = part1_phi - part2_phi;
	if (dphi > M_PI) dphi -= 2*M_PI;
	else if (dphi < - M_PI) dphi += 2*M_PI;
	return dphi;
};

// Returns Delta Eta between two particles
float delta_eta(	const float& part1_eta,
					const float& part2_eta) {
	return abs(part1_eta - part2_eta);
};

// Returns abs of given number
float calc_abs(const float& numb) {
	return std::abs(numb);
};
