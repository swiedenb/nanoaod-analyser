#ifndef event_cleaner_hh
#define event_cleaner_hh

#include "json.hpp"
#include <fstream>
#include <sstream>
#include <TROOT.h>
using json = nlohmann::json;


extern json goldenjson;
extern std::ifstream goldenjson_file;

bool json_check(const UInt_t& runnumber, const UInt_t& luminumber); 

#endif
