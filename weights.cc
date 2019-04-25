#include "weights.hh"
// Get Pileup weight from file
float pu_weight(const float nvtx_true){
    auto weight = config::pileup_hist->GetBinContent(nvtx_true);
    return weight;
};

float tau_scale_factor(const float tau_pt) {
	
}
