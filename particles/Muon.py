class Muon():
    def __init__(self,physVar):
        self.charge = physVar["Muon_charge"]
        self.pt = physVar["Muon_pt"]
        self.eta = physVar["Muon_eta"]
        self.phi = physVar["Muon_phi"]
        
        self.highPtId = physVar["Muon_highPtId"]
        self.mediumId = physVar["Muon_mediumId"]
        #self.mediumPromtId = physVar["Muon_mediumPromptId"]
        self.tightId = physVar["Muon_tightId"]
        self.tkIsoId = physVar["Muon_tkIsoId"]

    def setPt(self,pt):
        self.pt = pt

