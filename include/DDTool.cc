#include "DDTool.h"
#include <iostream> // std::cerr, std::endl
#include <iomanip> 
#include <assert.h> // assert
#include "TauIDSFTool.h"

TFile* ensureTFileDD(const TString filename, bool verbose=false){
  if(verbose)
    std::cout << "Opening " << filename << std::endl;
  TFile* file = new TFile(filename);
  if(!file or file->IsZombie()) {
    std::cerr << std::endl << "ERROR! Failed to open input file = '" << filename << "'!" << std::endl;
    assert(0);
  }
  return file;
}

TH2F extractTH1DD(const TFile* file, const std::string& histname){
  auto *tmp = dynamic_cast<TH2F*>((const_cast<TFile*>(file))->Get(histname.data()));
  if(!tmp){
    std::cerr << std::endl << "ERROR! Failed to load histogram = '" << histname << "' from input file!" << std::endl;
    assert(0);
  }
  auto hist = TH2F(*tmp);
  delete tmp;
  return hist;
}



DDTool::DDTool(const std::string& year){
  
}


float DDTool::getFF(float taupt, float taupt_o_jetpt, const std::string& year, const std::string& unc) {
  //auto filein = new TFile("/home/home1/institut_3a/wiedenbeck/PhD/analysis/nanoaod-analyser/cfg/fakerate/fakerate_2018.root","READ");

  TString datapath                = (TString) getenv("MY_ANALYSIS_PATH") + "/cfg/fakerate/fakerate_" + year + ".root";
  std::cout<<"Reading in FF from: "<<datapath<<std::endl;
  auto filein = new TFile(datapath,"READ");
  TH2D* histin =(TH2D*) filein->Get("iso_fake_rate");
  Int_t tauptbin = histin->GetXaxis()->FindBin(taupt);
  Int_t tauptojetptbin = histin->GetYaxis()->FindBin(taupt_o_jetpt);
  float FF  = histin->GetBinContent(tauptbin,tauptojetptbin);
  filein->Close();
  delete filein;
  return FF;
}



