#ifndef event_cleaner_hh
#define event_cleaner_hh

#include "json.hpp"
#include <fstream>
#include <sstream>
#include <TROOT.h>
#include <ROOT/RVec.hxx>

#include "../config.hh"
using json = nlohmann::json;

template < typename T >
using rvec = ROOT::VecOps::RVec<T>;

extern json goldenjson;

bool json_check(const UInt_t& runnumber, const UInt_t& luminumber); 


bool clean_gen_file(const rvec<int>& pdgID, const rvec<float>& mass);

#endif
