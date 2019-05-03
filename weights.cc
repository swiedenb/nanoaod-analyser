#include "weights.hh"
// Get Pileup weight from file
float pu_weight(const float& nvtx_true){
    auto weight = config::pileup_hist->GetBinContent(nvtx_true);
    return weight;
};

float apply_scale_factor() {
	if (config::run_type == "") {
		return config::tau_scale;
	} else if (config::run_type == "_TauScaleUp") {
		return config::tau_scale_up;
	} else if (config::run_type == "_TauScaleDown") {
		return config::tau_scale_down;
	}
	return 1.0;
}
