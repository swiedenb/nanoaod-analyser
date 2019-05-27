#include "weights.hh"

// Get Pileup weight from file
float pu_weight( const float nvtx_true){
    auto weight = config::pileup_hist->GetBinContent(nvtx_true);
    return weight;
};


// get tau scale factor
float apply_scale_factor() {
    if (config::run_type == "") {
        return config::tau_scale;
    } else if (config::run_type == "_TauScaleUp") {
        return config::tau_scale_up;
    } else if (config::run_type == "_TauScaleDown") {
        return config::tau_scale_down;
    }
    return 1.0;
};


// calculate top pt reweighting
float calc_top_pt_reweighting( const rvec<int>& gen_pdg,
                               const rvec<float>& gen_pt) {
   auto func = [] ( const double pt ) { return std::exp( 0.0615 - 0.0005 * pt ); };
   if (gen_pdg[2] == 6 && gen_pdg[3] == -6)
        return std::sqrt(func(gen_pt[2]) * func(gen_pt[3]));
   return 1.0;
};
