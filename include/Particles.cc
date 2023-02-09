#include "Particles.hh"

// MET xy correction


std::pair<double,double> METXYCorr_Met_MetPhi(double originalMet, double originalMet_phi, int runnb, int year, bool isData, int npv){

  enum TheRunEra{y2016B,y2016C,y2016D,y2016E,y2016F,y2016G,y2016H,y2017B,y2017C,y2017D,y2017E,y2017F,y2018A,y2018B,y2018C,y2018D,y2016MC,y2017MC,y2018MC};
  std::pair<double,double>  TheXYCorr_Met_MetPhi(originalMet,originalMet_phi);
  
  if(npv>100) npv=100;
  int runera =-1;
  bool usemetv2 =false;
  //bool usemetv2 =true;
  if(!isData && year == 2016) runera = y2016MC;
  else if(!isData && (year == 2017 or config::use2017XY == true)) {runera = y2017MC; usemetv2 =true;}
  else if(!isData && year == 2018) runera = y2018MC;
  
  else if(isData && runnb >=272007 &&runnb<=275376  ) runera = y2016B;
  else if(isData && runnb >=275657 &&runnb<=276283  ) runera = y2016C;
  else if(isData && runnb >=276315 &&runnb<=276811  ) runera = y2016D;
  else if(isData && runnb >=276831 &&runnb<=277420  ) runera = y2016E;
  else if(isData && runnb >=277772 &&runnb<=278808  ) runera = y2016F;
  else if(isData && runnb >=278820 &&runnb<=280385  ) runera = y2016G;
  else if(isData && runnb >=280919 &&runnb<=284044  ) runera = y2016H;
  
  else if(isData && runnb >=297020 &&runnb<=299329 ){ runera = y2017B; usemetv2 =true;}
  else if(isData && runnb >=299337 &&runnb<=302029 ){ runera = y2017C; usemetv2 =true;}
  else if(isData && runnb >=302030 &&runnb<=303434 ){ runera = y2017D; usemetv2 =true;}
  else if(isData && runnb >=303435 &&runnb<=304826 ){ runera = y2017E; usemetv2 =true;}
  else if(isData && runnb >=304911 &&runnb<=306462 ){ runera = y2017F; usemetv2 =true;}
  
  else if(isData && runnb >=315252 &&runnb<=316995 ) runera = y2018A;
  else if(isData && runnb >=316998 &&runnb<=319312 ) runera = y2018B;
  else if(isData && runnb >=319313 &&runnb<=320393 ) runera = y2018C;
  else if(isData && runnb >=320394 &&runnb<=325273 ) runera = y2018D;

  else {
    //Couldn't find data/MC era => no correction applied
    return TheXYCorr_Met_MetPhi;
  }
  
  double METxcorr(0.),METycorr(0.);

  if(!usemetv2){//Current recommendation for 2016 and 2018
    if(runera==y2016B) METxcorr = -(-0.0478335*npv -0.108032);
    if(runera==y2016B) METycorr = -(0.125148*npv +0.355672);
    if(runera==y2016C) METxcorr = -(-0.0916985*npv +0.393247);
    if(runera==y2016C) METycorr = -(0.151445*npv +0.114491);
    if(runera==y2016D) METxcorr = -(-0.0581169*npv +0.567316);
    if(runera==y2016D) METycorr = -(0.147549*npv +0.403088);
    if(runera==y2016E) METxcorr = -(-0.065622*npv +0.536856);
    if(runera==y2016E) METycorr = -(0.188532*npv +0.495346);
    if(runera==y2016F) METxcorr = -(-0.0313322*npv +0.39866);
    if(runera==y2016F) METycorr = -(0.16081*npv +0.960177);
    if(runera==y2016G) METxcorr = -(0.040803*npv -0.290384);
    if(runera==y2016G) METycorr = -(0.0961935*npv +0.666096);
    if(runera==y2016H) METxcorr = -(0.0330868*npv -0.209534);
    if(runera==y2016H) METycorr = -(0.141513*npv +0.816732);
    if(runera==y2017B) METxcorr = -(-0.259456*npv +1.95372);
    if(runera==y2017B) METycorr = -(0.353928*npv -2.46685);
    if(runera==y2017C) METxcorr = -(-0.232763*npv +1.08318);
    if(runera==y2017C) METycorr = -(0.257719*npv -1.1745);
    if(runera==y2017D) METxcorr = -(-0.238067*npv +1.80541);
    if(runera==y2017D) METycorr = -(0.235989*npv -1.44354);
    if(runera==y2017E) METxcorr = -(-0.212352*npv +1.851);
    if(runera==y2017E) METycorr = -(0.157759*npv -0.478139);
    if(runera==y2017F) METxcorr = -(-0.232733*npv +2.24134);
    if(runera==y2017F) METycorr = -(0.213341*npv +0.684588);
    if(runera==y2018A) METxcorr = -(0.362865*npv -1.94505);
    if(runera==y2018A) METycorr = -(0.0709085*npv -0.307365);
    if(runera==y2018B) METxcorr = -(0.492083*npv -2.93552);
    if(runera==y2018B) METycorr = -(0.17874*npv -0.786844);
    if(runera==y2018C) METxcorr = -(0.521349*npv -1.44544);
    if(runera==y2018C) METycorr = -(0.118956*npv -1.96434);
    if(runera==y2018D) METxcorr = -(0.531151*npv -1.37568);
    if(runera==y2018D) METycorr = -(0.0884639*npv -1.57089);
    if(runera==y2016MC) METxcorr = -(-0.195191*npv -0.170948);
    if(runera==y2016MC) METycorr = -(-0.0311891*npv +0.787627);
    if(runera==y2017MC) METxcorr = -(-0.217714*npv +0.493361);
    if(runera==y2017MC) METycorr = -(0.177058*npv -0.336648);
    if(runera==y2018MC) METxcorr = -(0.296713*npv -0.141506);
    if(runera==y2018MC) METycorr = -(0.115685*npv +0.0128193);
  }
  else {//these are the corrections for v2 MET recipe (currently recommended for 2017)
    if(runera==y2016B) METxcorr = -(-0.0374977*npv +0.00488262);
    if(runera==y2016B) METycorr = -(0.107373*npv +-0.00732239);
    if(runera==y2016C) METxcorr = -(-0.0832562*npv +0.550742);
    if(runera==y2016C) METycorr = -(0.142469*npv +-0.153718);
    if(runera==y2016D) METxcorr = -(-0.0400931*npv +0.753734);
    if(runera==y2016D) METycorr = -(0.127154*npv +0.0175228);
    if(runera==y2016E) METxcorr = -(-0.0409231*npv +0.755128);
    if(runera==y2016E) METycorr = -(0.168407*npv +0.126755);
    if(runera==y2016F) METxcorr = -(-0.0161259*npv +0.516919);
    if(runera==y2016F) METycorr = -(0.141176*npv +0.544062);
    if(runera==y2016G) METxcorr = -(0.0583851*npv +-0.0987447);
    if(runera==y2016G) METycorr = -(0.0641427*npv +0.319112);
    if(runera==y2016H) METxcorr = -(0.0706267*npv +-0.13118);
    if(runera==y2016H) METycorr = -(0.127481*npv +0.370786);
    if(runera==y2017B) METxcorr = -(-0.19563*npv +1.51859);
    if(runera==y2017B) METycorr = -(0.306987*npv +-1.84713);
    if(runera==y2017C) METxcorr = -(-0.161661*npv +0.589933);
    if(runera==y2017C) METycorr = -(0.233569*npv +-0.995546);
    if(runera==y2017D) METxcorr = -(-0.180911*npv +1.23553);
    if(runera==y2017D) METycorr = -(0.240155*npv +-1.27449);
    if(runera==y2017E) METxcorr = -(-0.149494*npv +0.901305);
    if(runera==y2017E) METycorr = -(0.178212*npv +-0.535537);
    if(runera==y2017F) METxcorr = -(-0.165154*npv +1.02018);
    if(runera==y2017F) METycorr = -(0.253794*npv +0.75776);
    if(runera==y2018A) METxcorr = -(0.362642*npv +-1.55094);
    if(runera==y2018A) METycorr = -(0.0737842*npv +-0.677209);
    if(runera==y2018B) METxcorr = -(0.485614*npv +-2.45706);
    if(runera==y2018B) METycorr = -(0.181619*npv +-1.00636);
    if(runera==y2018C) METxcorr = -(0.503638*npv +-1.01281);
    if(runera==y2018C) METycorr = -(0.147811*npv +-1.48941);
    if(runera==y2018D) METxcorr = -(0.520265*npv +-1.20322);
    if(runera==y2018D) METycorr = -(0.143919*npv +-0.979328);
    if(runera==y2016MC) METxcorr = -(-0.159469*npv +-0.407022);
    if(runera==y2016MC) METycorr = -(-0.0405812*npv +0.570415);
    if(runera==y2017MC) METxcorr = -(-0.182569*npv +0.276542);
    if(runera==y2017MC) METycorr = -(0.155652*npv +-0.417633);
    if(runera==y2018MC) METxcorr = -(0.299448*npv +-0.13866);
    if(runera==y2018MC) METycorr = -(0.118785*npv +0.0889588);
  }

  double CorrectedMET_x = originalMet *cos( originalMet_phi)+METxcorr;
  double CorrectedMET_y = originalMet *sin( originalMet_phi)+METycorr;

  double CorrectedMET = sqrt(CorrectedMET_x*CorrectedMET_x+CorrectedMET_y*CorrectedMET_y);
  double CorrectedMETPhi;
  if(CorrectedMET_x==0 && CorrectedMET_y>0) CorrectedMETPhi = TMath::Pi();
  else if(CorrectedMET_x==0 && CorrectedMET_y<0 )CorrectedMETPhi = -TMath::Pi();
  else if(CorrectedMET_x >0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x);
  else if(CorrectedMET_x <0&& CorrectedMET_y>0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x) + TMath::Pi();
  else if(CorrectedMET_x <0&& CorrectedMET_y<0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x) - TMath::Pi();
  else CorrectedMETPhi =0;

  TheXYCorr_Met_MetPhi.first= CorrectedMET;
  TheXYCorr_Met_MetPhi.second= CorrectedMETPhi;
  return TheXYCorr_Met_MetPhi;

};

// check if tau is in acceptance and fulfils id requirements
rvec<bool> tau_acceptance_and_id(const rvec<float>& pt, 
                                 const rvec<float>& eta,
                                 const rvec<float>& dz,
                                 const rvec<int>& charge,
                                 const rvec<bool>& dm,
                                 const rvec<UChar_t>& iso, 
                                 const rvec<UChar_t>& antiEle_disc, 
                                 const rvec<UChar_t>& antiMu_disc) {
        // tau pt cut
        auto mask_pt = pt > config::tau_pt;
        
        // tau eta cut
        auto mask_eta = abs(eta) < config::tau_eta;
        
        // tau dz cut
        auto mask_dz = abs(dz) < 0.2;
        
        // tau charge
        auto mask_charge = charge != 0;
        
        // tau decay mode
        auto mask_dm = dm;
                
        // tau iso requirement (1: VVLoose, 2: VLoose, 4: Loose, 8: Medium, 16: Tight, 32: VTight, 64: VVTight)
        auto mask_iso = (iso & config::tau_iso_WP) == config::tau_iso_WP;
        
        // Anti Electron discriminator (1: VLoose, 2: Loose, 4: Medium, 8: Tight, 16: VTight)
        auto mask_antiele = (config::tau_antiE_WP & antiEle_disc) == config::tau_antiE_WP;
        
        // Anti Muon discriminator (1: Loose, 2: Tight)
        auto mask_antimu = (config::tau_antiMu_WP & antiMu_disc) == config::tau_antiMu_WP;
        
        // return vector with true, if particle fulfils all requirements - else false
        return mask_pt & mask_eta & mask_charge & mask_dz & mask_dm & mask_iso & mask_antiele & mask_antimu;
};

rvec<bool> tau_acceptance_and_id_and_dm_noniso(const rvec<float>& pt, 
                                        const rvec<float>& eta,
                                        const rvec<float>& dz,
                                        const rvec<int>& charge,
                                        const rvec<bool>& dm,
                                        const rvec<int>& dm_number,
                                        const rvec<int>& jetIdx, 
                                        const rvec<UChar_t>& iso, 
                                        const rvec<UChar_t>& antiEle_disc, 
                                        const rvec<UChar_t>& antiMu_disc) {
        // tau pt cut
        auto mask_pt = pt > 50.;
        
        // tau eta cut
        auto mask_eta = abs(eta) < 2.3;
        
        // tau dz cut
        auto mask_dz = abs(dz) < 0.2;
        
        // tau charge
        auto mask_charge = charge != 0;
        
        // tau decay mode
        auto mask_dm = dm;

        // tau decay mode as integer
        auto mask = not (dm_number == 5 || dm_number == 6);
                
        // tau iso requirement
        auto mask_iso = ((1 & iso) == 1) && ((32 & iso) != 32);
        
        // Anti Electron discriminator
        auto mask_antiele = (config::tau_antiE_WP & antiEle_disc) == config::tau_antiE_WP;
        
        // Anti Muon discriminator
        auto mask_antimu = (config::tau_antiMu_WP & antiMu_disc) == config::tau_antiMu_WP;
        
        // Jet matching
        auto mask_jetmatching = jetIdx != -1;
        // return vector with true, if particle fulfils all requirements - else false
        return mask_pt & mask_eta & mask_dz & mask_charge & mask & mask_dm & mask_iso & mask_antiele & mask_antimu & mask_jetmatching;
        //return mask_pt & mask_eta & mask_dz & mask & mask_dm & mask_iso & mask_antiele & mask_antimu;
};
rvec<bool> tau_acceptance_and_id_and_dm(const rvec<float>& pt, 
                                        const rvec<float>& eta,
                                        const rvec<float>& dz,
                                        const rvec<int>& charge,
                                        const rvec<bool>& dm,
                                        const rvec<int>& dm_number,
                                        const rvec<int>& jetIdx, 
                                        const rvec<UChar_t>& iso, 
                                        const rvec<UChar_t>& antiEle_disc, 
                                        const rvec<UChar_t>& antiMu_disc) {
        // tau pt cut
        auto mask_pt = pt > config::tau_pt;
        
        // tau eta cut
        auto mask_eta = abs(eta) < config::tau_eta;
        
        // tau dz cut
        auto mask_dz = abs(dz) < 0.2;
        
        // tau charge
        auto mask_charge = charge != 0;
        
        // tau decay mode
        auto mask_dm = dm;

        // tau decay mode as integer
        auto mask = not (dm_number == 5 || dm_number == 6);
                
        // tau iso requirement
        auto mask_iso = (config::tau_iso_WP & iso) == config::tau_iso_WP;
        
        // Anti Electron discriminator
        auto mask_antiele = (config::tau_antiE_WP & antiEle_disc) == config::tau_antiE_WP;
        
        // Anti Muon discriminator
        auto mask_antimu = (config::tau_antiMu_WP & antiMu_disc) == config::tau_antiMu_WP;

        // Jet Matching
        auto mask_jetmatching = jetIdx != -1;
        
        // return vector with true, if particle fulfils all requirements - else false
        return mask_pt & mask_eta & mask_dz & mask_charge & mask & mask_dm & mask_iso & mask_antiele & mask_antimu & mask_jetmatching;
        //return mask_pt & mask_eta & mask_dz & mask & mask_dm & mask_iso & mask_antiele & mask_antimu;
};

// check if event has a well seperated di muon pair
rvec<bool> di_muon_id(const rvec<float>& pt, 
                                  const rvec<float>& eta, 
                                  const rvec<UChar_t>& id, 
                                  const rvec<float>& iso){ 
        // muon pt cut
        auto mask_pt = pt > 10.;
        
        // muon eta cut
        auto mask_eta = abs(eta) < config::muon_eta;
        
        // muon id
        auto mask_id = (id == config::muon_id_WP);
    
        // muon tracker iso 
        auto mask_iso = iso < 0.15;

        
        // return vector with true, if particle fulfils all requirements - else false
        return mask_pt & mask_eta & mask_id & mask_iso;
};
rvec<bool> muon_eta_mask(const rvec<float>&eta){
        auto mask_eta = abs(eta) < config::muon_eta;
        return mask_eta;
}
// check if muon is in acceptance and fulfils id requirements
rvec<bool> muon_acceptance_and_id(const rvec<float>& pt, 
                                  const rvec<float>& eta, 
                                  const rvec<bool>& pfcand,
                                  const rvec<UChar_t>& id, 
                                  const rvec<UChar_t>& iso) {
        // muon pt cut
        auto mask_pt = pt > config::muon_pt;
        
        // muon eta cut
        auto mask_eta = abs(eta) < config::muon_eta;
        
        // muon id
        auto mask_id = (id == config::muon_id_WP);
    
        // muon relative iso id
        auto mask_iso = (iso >= 1);

        // muon PF cand
        auto mask_pf = pfcand == true;
        
        // return vector with true, if particle fulfils all requirements - else false
        return mask_pt & mask_eta & mask_id & mask_iso;
};

// check if electron is in acceptance and fulfils id requirements
rvec<bool> ele_acceptance_and_id(const rvec<float>& pt,
                                 const rvec<float>& eta,
                                 const rvec<Int_t>& id) {
        // electron pt cut
        auto mask_pt = pt > config::ele_pt;
        
        // electron eta cut
        auto mask_eta = abs(eta) < config::ele_eta;
        
        // electron id cut
        auto mask_id = id >= 2;
        
        // return vector with true, if particle fulfils all requirements - else false
        return mask_pt & mask_eta & mask_id;
};


rvec<bool> ele_acceptance_and_simpleid(const rvec<float>& pt,
                                       const rvec<float>& eta,
                                       const rvec<bool>& id) {
        // electron pt cut
        auto mask_pt = pt > config::ele_pt;
        
        // electron eta cut
        auto mask_eta = abs(eta) < config::ele_eta;
        
        // electron id cut
        auto mask_id = id;
        
        // return vector with true, if particle fulfils all requirements - else false
        return mask_pt & mask_eta & mask_id;
};


// calc tau energy scale
rvec< rvec< float > > calc_tau_energy_scale(const rvec<int>& tau_decayMode,
                                            const rvec<UChar_t>& tau_genPartFlav,
                                            const rvec<float>& tau_pt,
                                            const rvec<float>& tau_eta,
                                            const rvec<float>& tau_phi,
                                            const rvec<float>& tau_mass,
                                            std::string type = "") {  
        rvec<rvec<float>> binder(4);
        
        rvec<float> new_tau_pt(tau_pt.size());                          
        rvec<float> new_tau_eta(tau_eta.size());                                
        rvec<float> new_tau_phi(tau_phi.size());                                
        rvec<float> new_tau_mass(tau_mass.size());                              
        
        for (uint i = 0; i<tau_pt.size(); i++) {

            TLorentzVector p1; 
            p1.SetPtEtaPhiM(tau_pt[i], tau_eta[i], tau_phi[i], tau_mass[i]);

            if (type == "_ESUp") {
                p1 *= config::tau_dm_scale->getTES(tau_pt[i], tau_decayMode[i], (int) tau_genPartFlav[i], "Up");
            } else if (type == "_ESDown") {
                p1 *= config::tau_dm_scale->getTES(tau_pt[i], tau_decayMode[i], (int) tau_genPartFlav[i], "Down");
            } else {
                p1 *= config::tau_dm_scale->getTES(tau_pt[i], tau_decayMode[i], (int) tau_genPartFlav[i], "");
            }

            new_tau_pt[i] = p1.Pt();
            new_tau_eta[i] = p1.Eta();
            new_tau_phi[i] = p1.Phi();
            new_tau_mass[i] = p1.M();
        }
        
        binder[0] = new_tau_pt;                         
        binder[1] = new_tau_eta;                                
        binder[2] = new_tau_phi;                                
        binder[3] = new_tau_mass;
        
        return binder;
}


// apply tau energy scale
rvec <float> apply_tau_energy_scale(const rvec< rvec<float> > binder, const int& index) {
    return binder[index];
}
