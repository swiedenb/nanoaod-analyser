#include "event_cleaner.hh"

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


// Gen Cleaning
bool clean_gen_file(const rvec<int>& pdgID, const rvec<float>& mass) {
	if (config::gen_pdgID.size() == 0)
		return true;
		
    else {
        if ( config::gen_pdgID.size() == 1 and config::gen_pdgID.at(0) == 0) return true;
	    // mass gen cleaning
	    if (config::cut_type == "mass") {
	    	for (uint i = 0; i < pdgID.size(); i++) {
                for (uint j = 0; j < config::gen_pdgID.size(); j++) {
	    		    if ( abs(pdgID[i]) == config::gen_pdgID.at(j)) {
	    		    	if ((mass[i] > config::cut_value_min) && (mass[i] < config::cut_value_max)) 
	    		    		return true;
	    		    	else
	    		    		return false;
	    		    }
                }
	    	}
	    }
    }
	return true;
};

