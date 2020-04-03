#include "MuonSFTool.h"
#include <iostream> // std::cerr, std::endl
#include <iomanip> 
#include <assert.h> // assert
#include "TauIDSFTool.h"

TFile* ensureTFileMuon(const TString filename, bool verbose=false){
  if(verbose)
    std::cout << "Opening " << filename << std::endl;
  TFile* file = new TFile(filename);
  if(!file or file->IsZombie()) {
    std::cerr << std::endl << "ERROR! Failed to open input file = '" << filename << "'!" << std::endl;
    assert(0);
  }
  return file;
}

TH2D extractTH1Muon(const TFile* file, const std::string& histname){
  auto *tmp = dynamic_cast<TH2D*>((const_cast<TFile*>(file))->Get(histname.data()));
  if(!tmp){
    std::cerr << std::endl << "ERROR! Failed to load histogram = '" << histname << "' from input file!" << std::endl;
    assert(0);
  }
  auto hist = TH2D(*tmp);
  delete tmp;
  return hist;
}



MuonSFTool::MuonSFTool(const std::string& year, const std::string& id, const std::string& iso): ID(id), ISO(iso){
  
  bool verbose = false;
  std::string datapath                = (std::string) getenv("MY_ANALYSIS_PATH") + "/cfg/SF/MuonSFs/";
  std::vector<std::string> years      = {"2016Legacy","2017ReReco","2018ReReco"};
  std::vector<std::string> MuIDs  = {"LooseID", "MediumID", "MediumPromptID", "TightID", "SoftID", "HighPtID", "TrkHighPtID"};
  std::vector<std::string> MuISOs  = {"LooseRelIso", "LooseRelTkIso", "TightRelIso", "TightRelTkIso"};
  
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
  
  if(std::find(MuIDs.begin(),MuIDs.end(),ID)!=MuIDs.end()){
      TString filename;
      std::string histname;
      std::string histname_syst;
      std::string histname_stat;
      filename = Form("%s/RunABCD_SF_ID_%s.root",datapath.data(),year.data());
      TFile* file = ensureTFileMuon(filename,verbose);
      if( id == "HighPtID" or id == "TrkHighPtID"){
        if ( id == "HighPtID"){
            histname = "NUM_HighPtID_DEN_TrackerMuons_pair_newTuneP_probe_pt_abseta";
            histname_syst = "NUM_HighPtID_DEN_TrackerMuons_pair_newTuneP_probe_pt_abseta_syst";
            histname_stat = "NUM_HighPtID_DEN_TrackerMuons_pair_newTuneP_probe_pt_abseta_stat";
        }
        if ( id == "TrkHighPtID"){
            histname = "NUM_TrikHighPtID_DEN_TrackerMuons_pair_newTuneP_probe_pt_abseta";
            histname_syst = "NUM_TrikHighPtID_DEN_TrackerMuons_pair_newTuneP_probe_pt_abseta_syst";
            histname_stat = "NUM_TrikHighPtID_DEN_TrackerMuons_pair_newTuneP_probe_pt_abseta_stat";
        }
      }
      else{
        histname = Form("NUM_%s_DEN_TrackerMuons_pt_abseta",ID.data());
        histname_syst = Form("NUM_%s_DEN_TrackerMuons_pt_abseta_syst",ID.data());
        histname_stat = Form("NUM_%s_DEN_TrackerMuons_pt_abseta_stat",ID.data());
      }
      hist = extractTH1Muon(file,histname);
      hist_syst = extractTH1Muon(file,histname_syst);
      hist_stat = extractTH1Muon(file,histname_stat);
      hist.SetDirectory(0);
      hist_syst.SetDirectory(0);
      hist_stat.SetDirectory(0);
      file->Close();
      delete file;
      if(std::find(MuISOs.begin(),MuISOs.end(),iso)!=MuISOs.end()){
          TString filename_iso;
          std::string histname_iso;
          std::string histname_iso_syst;
          std::string histname_iso_stat;
          filename_iso = Form("%s/RunABCD_SF_ISO_%s.root",datapath.data(),year.data());
          TFile* file_iso = ensureTFileMuon(filename_iso,verbose);
          if( id == "HighPtID" or id == "TrkHighPtID" or "TightID"){
            if ( id == "HighPtID"){
                histname_iso = Form("NUM_%s_DEN_%sandIPCut_pair_newTuneP_probe_pt_abseta",ISO.data(),ID.data());
                histname_iso_syst = Form("NUM_%s_DEN_%sandIPCut_pair_newTuneP_probe_pt_abseta_syst",ISO.data(),ID.data());
                histname_iso_stat = Form("NUM_%s_DEN_%sandIPCut_pair_newTuneP_probe_pt_abseta_stat",ISO.data(),ID.data());
            }
            if ( id == "TrkHighPtID"){
                histname_iso = Form("NUM_%s_DEN_%s_pair_newTuneP_probe_pt_abseta",ISO.data(),ID.data());
                histname_iso_syst = Form("NUM_%s_DEN_%s_pair_newTuneP_probe_pt_abseta_syst",ISO.data(),ID.data());
                histname_iso_stat = Form("NUM_%s_DEN_%s_pair_newTuneP_probe_pt_abseta_stat",ISO.data(),ID.data());
            }
            if ( id == "TightID"){
                histname_iso = Form("NUM_%s_DEN_%sandIPCut_pt_abseta",ISO.data(),ID.data());
                histname_iso_syst = Form("NUM_%s_DEN_%sandIPCut_pt_abseta_syst",ISO.data(),ID.data());
                histname_iso_stat = Form("NUM_%s_DEN_%sandIPCut_pt_abseta_stat",ISO.data(),ID.data());
            }
          }
          else{
            histname_iso = Form("NUM_%s_DEN_%s_pt_abseta",ISO.data(),ID.data());
            histname_iso_syst = Form("NUM_%s_DEN_%s_pt_abseta_syst",ISO.data(),ID.data());
            histname_iso_stat = Form("NUM_%s_DEN_%s_pt_abseta_stat",ISO.data(),ID.data());
          }
          hist_iso = extractTH1Muon(file_iso,histname_iso);
          hist_iso_syst = extractTH1Muon(file_iso,histname_iso_syst);
          hist_iso_stat = extractTH1Muon(file_iso,histname_iso_stat);
          hist_iso.SetDirectory(0);
          hist_iso_syst.SetDirectory(0);
          hist_iso_stat.SetDirectory(0);
          file_iso->Close();
          delete file_iso;
      }else{
          std::cerr << "Did not recognize muon ISO '" << ID << "'!" << std::endl;
          assert(0);
      }
    }else{
        std::cerr << "Did not recognize muon ID '" << ID << "'!" << std::endl;
        assert(0);
    }
}

float MuonSFTool::getSFID(double pt, double eta, const std::string& unc) {
  if (pt > 100.) pt = 100.;
  Int_t ptbin = hist.GetXaxis()->FindBin(pt);
  Int_t etabin = hist.GetYaxis()->FindBin(fabs(eta));
  float SF  = hist.GetBinContent(ptbin,etabin);
  if(unc=="Up")
    SF += hist_syst.GetBinError(hist_syst.FindBin(pt,fabs(eta)));
  else if(unc=="Down")
    SF -= hist_syst.GetBinError(hist_syst.FindBin(pt,fabs(eta)));
  return SF;
}
float MuonSFTool::getSFISO(double pt, double eta, const std::string& unc) {
  if (pt > 100.) pt = 100.;
  Int_t ptbin = hist.GetXaxis()->FindBin(pt);
  Int_t etabin = hist.GetYaxis()->FindBin(fabs(eta));
  float SF  = hist.GetBinContent(ptbin,etabin);
  if(unc=="Up")
    SF += hist_syst.GetBinError(hist_syst.FindBin(pt,fabs(eta)));
  else if(unc=="Down")
    SF -= hist_syst.GetBinError(hist_syst.FindBin(pt,fabs(eta)));
  return SF;
}




