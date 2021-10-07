#include "EleSFTool.h"
#include <iostream> // std::cerr, std::endl
#include <iomanip> 
#include <assert.h> // assert
#include "TauIDSFTool.h"

TFile* ensureTFileEle(const TString filename, bool verbose=false){
  if(verbose)
    std::cout << "Opening " << filename << std::endl;
  TFile* file = new TFile(filename);
  if(!file or file->IsZombie()) {
    std::cerr << std::endl << "ERROR! Failed to open input file = '" << filename << "'!" << std::endl;
    assert(0);
  }
  return file;
}

TH2F extractTH1Ele(const TFile* file, const std::string& histname){
  auto *tmp = dynamic_cast<TH2F*>((const_cast<TFile*>(file))->Get(histname.data()));
  if(!tmp){
    std::cerr << std::endl << "ERROR! Failed to load histogram = '" << histname << "' from input file!" << std::endl;
    assert(0);
  }
  auto hist = TH2F(*tmp);
  delete tmp;
  return hist;
}



EleSFTool::EleSFTool(const std::string& year, const std::string& id): ID(id){
  
  bool verbose = false;
  std::string datapath                = (std::string) getenv("MY_ANALYSIS_PATH") + "/cfg/SF/EleSFs/";
  std::vector<std::string> years      = {"2016Legacy","2017ReReco","2018ReReco"};
  std::vector<std::string> EleIDs  = {"heep"};
  
  if(std::find(years.begin(),years.end(),year)==years.end()){
    std::cerr << std::endl << "ERROR! '"<<year<<"' is not a valid year! Please choose from ";
    std::vector<std::string>::iterator it = years.begin();
    for(it=years.begin(); it!=years.end(); it++){
      if(it!=years.begin()) std::cerr << ", ";
      std::cerr << *it;
    }
    std::cerr << std::endl;
    assert(0);
  }
  TString filename_reco;
  std::string histname_reco;
  filename_reco = Form("%s/sf_ele_reco_%s.root",datapath.data(),year.data());
  TFile* file_reco = ensureTFileEle(filename_reco,verbose);
  histname_reco = "EGamma_SF2D";
  hist_reco = extractTH1Ele(file_reco,histname_reco);
  hist_reco.SetDirectory(0);
  file_reco->Close();
  delete file_reco;
  
//  if(std::find(EleIDs.begin(),EleIDs.end(),ID)!=EleIDs.end()){
//      TString filename;
//      std::string histname;
//      filename = Form("%s/sf_ele_%s_%s.root",datapath.data(),ID.data(),year.data());
//      TFile* file = ensureTFileEle(filename,verbose);
//      histname = "SF_MC";
//      hist = extractTH1Ele(file,histname);
//      hist.SetDirectory(0);
//      file->Close();
//      delete file;
//    }else{
//        std::cerr << "Did not recognize ele ID '" << ID << "'!" << std::endl;
//        assert(0);
//    }
}

float EleSFTool::getSFID(double pt, double eta, const std::string& unc) {
  Int_t ptbin = hist.GetXaxis()->FindBin(pt);
  Int_t etabin = hist.GetYaxis()->FindBin(fabs(eta));
  float SF  = hist.GetBinContent(etabin,ptbin);
  if(unc=="Up")
    if( eta < 1.4442){
        if (pt < 90) SF += 0.01;
        else SF += std::min(1+(pt-90.)*0.0022,3.)/100.;
    }
    else if (1.566< eta < 2.5){
        if (pt < 90) SF += 0.02;
        else SF += std::min(1+(pt-90.)*0.0143,5.)/100.;
   }
  else if(unc=="Down")
    if (eta < 1.4442){
        if (pt < 90) SF -= 0.01;
        else SF -= std::min(1+(pt-90.)*0.0022,3.)/100.;
    }
    else if (1.566< eta < 2.5){
        if (pt < 90) SF -= 0.02;
        else SF -= std::min(1+(pt-90.)*0.0143,5.)/100.;
   }
  return SF;
}

float EleSFTool::getSFReco(double pt, double eta, const std::string& unc) {
  Int_t ptbin = hist.GetXaxis()->FindBin(pt);
  Int_t etabin = hist.GetYaxis()->FindBin(fabs(eta));
  float SF  = hist.GetBinContent(etabin,ptbin);
  return SF;
}



