#include "event_cleaner.hh"

// Check if run number and lumi ranges are in json.txt - still need to write
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
