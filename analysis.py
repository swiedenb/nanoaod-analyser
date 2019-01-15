from particles.Muon import *
from particles.Electron import *
from particles.Tau import *
from cfg.cfg import *

import uproot
import ROOT as r
import numpy as np
import sys


def progress(count, total, status=''):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
    sys.stdout.flush() # As suggested by Rom Ruben (see: http://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console/27871113#comment50529068_27871113)

def read_in_arrays(datapath):
    tree = uproot.open(datapath)["Events"]
    muonVar = {}
    eleVar = {}
    tauVar = {}
    trigger_array = []
    nMuons = 0.
    nElectrons = 0.
    nTaus = 0.
    for trigger in trigger_list:
        trigger_array.append(tree.array(trigger))
    if useMuon:
        nMuons = tree.array("nMuon")
        muonVar["Muon_charge"] = tree.array("Muon_charge")
        muonVar["Muon_pt"] = tree.array("Muon_pt")
        muonVar["Muon_eta"] = tree.array("Muon_eta")
        muonVar["Muon_phi"] = tree.array("Muon_phi")
        muonVar["Muon_highPtId"] = tree.array("Muon_highPtId")
        muonVar["Muon_mediumId"] = tree.array("Muon_mediumId")
        #muonVar["Muon_mediumPromptId"] = tree.array("Muon_mediumPromptId")
        muonVar["Muon_tightId"] = tree.array("Muon_tightId")
        muonVar["Muon_tkIsoId"] = tree.array("Muon_tkIsoId")
    if useEle:
        eleVar["Electron_charge"] = tree.array("Electron_charge")
        eleVar["Electron_pt"] = tree.array("Electron_pt")
        eleVar["Electron_eta"] = tree.array("Electron_eta")
        eleVar["Electron_phi"] = tree.array("Electron_phi")
        nElectrons = tree.array("nElectron")
    if useTau:
        tauVar["Tau_charge"] = tree.array("Tau_charge")
        tauVar["Tau_pt"] = tree.array("Tau_pt")
        tauVar["Tau_eta"] = tree.array("Tau_eta")
        tauVar["Tau_phi"] = tree.array("Tau_phi")
        tauVar["Tau_idDecayMode"] = tree.array("Tau_idDecayMode")
        #tauVa(tree.array("Tau_idMVAoldDM"))
        tauVar["Tau_idAntiEle"] = tree.array("Tau_idAntiEle")
        tauVar["Tau_idAntiMu"] = tree.array("Tau_idAntiMu")
        nTaus = tree.array("nTau")

    return muonVar, nMuons, eleVar, nElectrons, tauVar, nTaus, trigger_array

def invMass(par1,par2):
    return np.sqrt(2*par1.pt*par2.pt*(np.cosh(par1.eta-par2.eta)-np.cos(par1.phi-par2.phi)))

muons,nMuons, electrons, nElectrons,taus, nTaus, trigger_array = read_in_arrays("/net/scratch_cms3a/wiedenbeck/nanoaodfiles/TTbar.root")
muon_pt_hist = r.TH1D("muon_pt", "muon_pt", 2000,0,2000)
tau_pt_hist = r.TH1D("tau_pt", "tau_pt", 2000,0,2000)
nMuons_hist = r.TH1D("nMuons", "nMuons", 10,0,10)
for event in range(len(nMuons)):
    progress(event, len(nMuons), status='Doing very long job')
    muonCollection = np.array([])
    eleCollection = np.array([])
    tauCollection = np.array([])
    trigger_bool = False
    for trigger in trigger_array:
        if trigger[event]:trigger_bool = True
    if not trigger_bool: continue
    if useMuon:
        for muon in range(nMuons[event]):
            muonDict = {}
            for key in muons:
                muonDict[key] = muons[key][event][muon]
            muonCollection =  np.append(muonCollection,Muon(muonDict))
    if useEle:
         for electron in range(nElectrons[event]):
            eleDict = {}
            for key in electrons:
                eleDict[key] = electrons[key][event][tau]
            eleCollection = np.append(eleCollection,Electron(eleDict))
    if useTau:
        for tau in range(nTaus[event]):
            tauDict = {}
            for key in taus:
                tauDict[key] = taus[key][event][tau]
            tauCollection = np.append(tauCollection,Tau(tauDict))
    for muon in muonCollection: 
        if muon.pt >53. and muon.highPtId and abs(muon.eta)<2.4 and muon.tkIsoId == 1:
                muon_pt_hist.Fill(muon.pt)
    nMuons_hist.Fill(nMuons[event])
    #for tau in tauCollection:
    #    if tau.decay_mode and tau.anti_muon and tau.anti_ele: tau_pt_hist.Fill(tau.pt)
outfile = r.TFile("analysisOut.root","RECREATE")
muon_pt_hist.Write()
tau_pt_hist.Write()
nMuons_hist.Write()
outfile.Close()

