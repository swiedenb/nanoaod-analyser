import ROOT

mc_file = ROOT.TFile("mc_pileup_2018.root","READ")
data_file = ROOT.TFile("Data_PileUp_2018_25ns_69.2.root","READ")
data_file_up = ROOT.TFile("Data_PileUp_2018_25ns_69.2_up5.root", "READ")
data_file_down = ROOT.TFile("Data_PileUp_2018_25ns_69.2_down5.root", "READ")

mc_hist = mc_file.Get("pileup")
data_hist = data_file.Get("pileup")
data_hist_up = data_file_up.Get("pileup")
data_hist_down = data_file_down.Get("pileup")

mc_hist.Sumw2()
data_hist.Sumw2()

mc_hist.Scale(1./mc_hist.Integral())
data_hist.Scale(1./data_hist.Integral())
data_hist_up.Scale(1./data_hist.Integral())
data_hist_down.Scale(1./data_hist.Integral())
mc_hist.Sumw2()
data_hist.Sumw2()

data_hist.Divide(mc_hist)
data_hist_up.Divide(mc_hist)
data_hist_down.Divide(mc_hist)

data_hist_up.SetName("pileup_up")
data_hist_down.SetName("pileup_down")

outfile = ROOT.TFile("pileup_weight_2018.root","RECREATE")
data_hist.Write()
data_hist_up.Write()
data_hist_down.Write()
outfile.Close()
