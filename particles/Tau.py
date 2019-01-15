class Tau():
    def __init__(self,physVar):
        self.charge = physVar["Tau_charge"]
        self.pt = physVar["Tau_pt"]
        self.eta = physVar["Tau_eta"]
        self.phi = physVar["Tau_phi"]
        self.decay_mode = physVar["Tau_idDecayMode"]
        #self.iso = physVar[5]
        self.anti_ele = physVar["Tau_idAntiEle"]
        self.anti_muon = physVar["Tau_idAntiMu"]

       # self.dxy = detVal[0]
       # self.dxyErr = detVal[1]
       # self.dz = detVal[2]
       # self.dzErr = detVal[3]
