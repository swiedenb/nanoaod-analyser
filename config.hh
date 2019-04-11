#ifndef CONFIG_H
#define CONFIG_H

#include "include/json.hpp"
using json = nlohmann::json;

namespace config {
	extern float tau_pt;
	extern float tau_eta;
	extern float met_pt;
    extern char tau_iso_WP[10];
    extern char tau_antiE_WP[10];
    extern char tau_antiMu_WP[10];
    
    extern int gen_pdgID;
    extern std::string cut_type;
    extern float cut_value_min;
    extern float cut_value_max;
    
    bool load_config_file(json cfg);
}

#endif
