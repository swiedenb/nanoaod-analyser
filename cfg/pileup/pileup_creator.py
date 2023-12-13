import ROOT

mc_file = ROOT.TFile("mc_pileup_Summer16.root","READ")
#mc_file = ROOT.TFile("mc_pileup_2018.root","READ")
#mc_file = ROOT.TFile("mc_pileup_2017.root","READ")
#data_file = ROOT.TFile("data_pileup_2016_69p2.root","READ")
#data_file = ROOT.TFile("Data_PileUp_may_2018_25ns_69.2.root","READ")
data_file = ROOT.TFile("PileupData_GoldenJSON_Full2016.root","READ")
#data_file = ROOT.TFile("data_pileup_2017_69p2.root","READ")
#data_file_up = ROOT.TFile("Data_PileUp_may_2018_25ns_69.2_up5.root", "READ")
#data_file_down = ROOT.TFile("Data_PileUp_may_2018_25ns_69.2_down5.root", "READ")
#data_file_up = ROOT.TFile("Data_PileUp_2018_25ns_69.2_up5.root", "READ")
#data_file_down = ROOT.TFile("Data_PileUp_2018_25ns_69.2_down5.root", "READ")

mc_hist = mc_file.Get("pileup")
data_hist = data_file.Get("pileup")
#data_hist_up = data_file_up.Get("pileup")
#data_hist_down = data_file_down.Get("pileup")
data_hist_up = data_file.Get("pileup_plus")
data_hist_down = data_file.Get("pileup_minus")


#mc_hist.Scale(1./mc_hist.Integral())
data_hist.Scale(1./data_hist.Integral())
data_hist_up.Scale(1./data_hist.Integral())
data_hist_down.Scale(1./data_hist.Integral())
new_mc_hist = ROOT.TH1D("mc_pileup","mc_pileup",100,0,100)
new_data_up_hist = ROOT.TH1D("pileup_up","pileup_up",100,0,100)
new_data_down_hist = ROOT.TH1D("pileup_down","pileup_down",100,0,100)
new_data_hist = ROOT.TH1D("pileup","pileup",100,0,100)

for i in range(0,data_hist.GetNbinsX() + 1):
    new_mc_hist.SetBinContent(i,mc_hist.GetBinContent(i))
    new_data_hist.SetBinContent(i,data_hist.GetBinContent(i))
    new_data_up_hist.SetBinContent(i,data_hist_up.GetBinContent(i))
    new_data_down_hist.SetBinContent(i,data_hist_down.GetBinContent(i))


#mc_hist.Sumw2()
#data_hist.Sumw2()
new_data_hist.Scale(1./new_data_hist.Integral())
new_data_up_hist.Scale(1./new_data_up_hist.Integral())
new_data_down_hist.Scale(1./new_data_down_hist.Integral())
#new_data_hist.Sumw2()
#new_mc_hist.Sumw2()
new_data_hist.Divide(new_mc_hist)

new_data_up_hist.Divide(new_mc_hist)
new_data_down_hist.Divide(new_mc_hist)

for i in range(new_data_hist.GetNbinsX() + 1):
    new_data_hist.SetBinError(i,new_data_hist.GetBinContent(i)*0.049)
    new_data_up_hist.SetBinError(i,new_data_up_hist.GetBinContent(i)*0.049)
    new_data_down_hist.SetBinError(i,new_data_down_hist.GetBinContent(i)*0.049)
#
#data_hist_up.SetName("pileup_up")
#data_hist_down.SetName("pileup_down")

outfile = ROOT.TFile("pileup_weight_2016.root","RECREATE")
#outfile = ROOT.TFile("pileup_weight_2018.root","RECREATE")
#outfile = ROOT.TFile("pileup_weight_2017.root","RECREATE")
#data_hist.Sumw2()
#data_hist.Write()
#data_hist_up.Write()
#data_hist_down.Write()
new_data_hist.Write()
new_data_up_hist.Write()
new_data_down_hist.Write()
outfile.Close()
