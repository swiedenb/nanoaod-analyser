#ifndef event_cleaner_hh
#define event_cleaner_hh

#include "json.hpp"
#include <fstream>
#include <sstream>
#include <TROOT.h>
#include <ROOT/RVec.hxx>

#include "config.hh"
using json = nlohmann::json;

template < typename T >
using rvec = ROOT::VecOps::RVec<T>;

extern json goldenjson;

bool json_check(const UInt_t& runnumber, const UInt_t& luminumber); 

rvec <int> pair_idx(const rvec<int>& pdgID, const rvec<int>& mother_idx);

bool clean_gen_file(const float& lheht, const rvec<float>& lhept , const rvec<float>& lheeta, const rvec<float>& lhephi, const rvec<float>& lhemass, const rvec<int>& lhepdgid, const rvec<int>& status);

#endif
