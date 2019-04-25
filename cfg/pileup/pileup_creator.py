import ROOT

mc_file = ROOT.TFile("mc_pileup_2017.root","READ")
data_file = ROOT.TFile("data_pileup_2018.root","READ")

mc_hist = mc_file.Get("pileup")
data_hist = data_file.Get("pileup")

mc_hist.Scale(1./mc_hist.Integral())
data_hist.Scale(1./data_hist.Integral())

data_hist.Divide(mc_hist)

outfile = ROOT.TFile("pileup_weight.root","RECREATE")
data_hist.Write()
outfile.Close()
