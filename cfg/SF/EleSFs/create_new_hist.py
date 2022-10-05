from ROOT import *

infile = TFile("sf_ele_heep_2017ReReco.root","READ")

hist = infile.Get("SF_MC")


hist_new = TH2F("SF_MC_2017","SF_MC_2017",50,-2.5,2.5,1,0,10000)

for i in range(hist.GetNbinsX()+1):
    hist_new.SetBinContent(i+1,1,hist.GetBinContent(i+1,1))

outfile = TFile("sf_ele_heep_2017_test.root","RECREATE")
hist_new.Write()
outfile.Close()

