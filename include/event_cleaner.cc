#include "event_cleaner.hh"
#include "TLorentzVector.h"

// Json Cleaning
bool json_check(const UInt_t& runnumber, const UInt_t& luminumber) {
    std::string runnumber_string = std::to_string(runnumber);
    char const *runnumber_char = runnumber_string.c_str();

    if (goldenjson.find(runnumber_char) != goldenjson.end()){
        for ( auto it = goldenjson[runnumber_char].begin(); it != goldenjson[runnumber_char].end(); ++it){
            auto lumirange = *it;
            if (luminumber >= lumirange[0] && luminumber <= lumirange[1]) return true; 
        }
    }
    return false;
};

rvec <int> pair_idx(const rvec<int>& pdgID, const rvec<int>& mother_idx) {
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
            //const auto idx = std::distance(pdgID.begin(), iter++);
            //pts.push_back(pt[idx]);
            //etas.push_back(eta[idx]);
            //phis.push_back(phi[idx]);
            //masses.push_back(mass[idx]);
            //indices.push_back(std::distance(pdgID.begin(), iter++));
            //indices.push_back(idx);
         }
      }
      if (second)
      {
         break;
      }
   }
   //std::exit(1);
   if (first and second){
       idx[0] = idx_1;
       idx[1] = idx_2;
       return idx;
   }
   else{
        std::cout<<"Could not find the defined particles! return true"<<std::endl;
        return idx;
   }
   return idx;
};

// Gen Cleaning
bool clean_gen_file( const float& lheht, const rvec<float>& lhept , const rvec<float>& lheeta, const rvec<float>& lhephi, const rvec<float>& lhemass, const rvec<int>& lhepdgid,const rvec<int>& status) {
    
    if (config::cut_type == "HT"){
        if ((lheht >= config::cut_value_min) && (lheht < config::cut_value_max)){
            return true;

        }
        else
            return false;
    }
    else if (config::gen_pdgID.empty()){
        return true;
    }
    else {
        if ( config::gen_pdgID.size() == 1 and config::gen_pdgID.at(0) == 0) return true;
        // mass gen cleaning
        if (config::cut_type == "mass") {
            if (config::DY == true){
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
                        if ((mass >= config::cut_value_min) && (mass < config::cut_value_max))
                            return true;
                        else
                            return false;
                    }
                }
                
            }
           //     TLorentzVector l1,l2;
           //     int l1_pdg;
           //     int l2_pdg;
           //     bool reject = false;
           //     for (uint i = 0; i < lhepdgid.size();i++){
           //         if (abs(lhepdgid[i]) == 11 || abs(lhepdgid[i]) == 13 || abs(lhepdgid[i]) == 15){
           //             l1.SetPtEtaPhiM(lhept[i],lheeta[i],lhephi[i],lhemass[i]);
           //             l1_pdg = lhepdgid[i];
           //         }
           //         for (uint j = 0; j < i; j++){
           //             if (abs(lhepdgid[j]) == 11 || abs(lhepdgid[j]) == 13 || abs(lhepdgid[j]) == 15){
           //                 l2.SetPtEtaPhiM(lhept[j],lheeta[j],lhephi[j],lhemass[j]);
           //                 l2_pdg = lhepdgid[j];
           //             }
           //             if((l1+l2).M() > config::cut_value_max or (l1+l2).M() < config::cut_value_min) reject = true;
           //         if (reject) break;
           //         }
           //         if (reject) break;
           //     }
           //     if (reject) return false;
           // }
           // if(config::cut_single == 1){
           //     for (uint i = 0; i < pdgID.size(); i++) {
           //         for (uint j = 0; j < config::gen_pdgID.size(); j++) {
           //             if ( abs(pdgID[i]) == config::gen_pdgID.at(j)) {
           //                 if ((mass[i] > config::cut_value_min) && (mass[i] < config::cut_value_max)) 
           //                     return true;
           //                 else
           //                     return false;
           //             }
           //         }
           //     }
           // }
            if(config::TT or config::WW){
               int idx_1;
               int idx_2;
               bool first{false};
               bool second{false};

               for (const auto gen_pdg : config::gen_pdgID){
                  auto iter = lhepdgid.begin();
                  while (true)
                  {
                     iter = std::find_if(iter, lhepdgid.end(), [rhs = gen_pdg](const auto lhs){return abs(lhs) == rhs;});
                     const bool found = iter != lhepdgid.end();
                     if (not found)
                     {
                        break;
                     }
                     else
                     {
                        if (not first)
                        {
                           idx_1 = std::distance(lhepdgid.begin(), iter++);
                           if (status[idx_1] == 1){
                               first = true;
                           }
                        }
                        else if (not second)
                        {
                           idx_2 = std::distance(lhepdgid.begin(), iter++);
                           if (status[idx_2] == 1){
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
                   part1.SetPtEtaPhiM(lhept[idx_1],lheeta[idx_1],lhephi[idx_1],lhemass[idx_1]);
                   part2.SetPtEtaPhiM(lhept[idx_2],lheeta[idx_2],lhephi[idx_2],lhemass[idx_2]);
                   auto inv_mass = (part1 + part2).M();
                   if( inv_mass >= config::cut_value_min and inv_mass < config::cut_value_max){
                        return true;
                   }
                   else{
                        return false;
                   }
               }
               else{
                    std::cout<<"Could not find the defined particles! return true"<<std::endl;
                    return true;
               }

            }
        }
    }
    return true;
};

