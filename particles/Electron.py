class Electron():
    #def __init__(self,physVar,detVal,idDict):
    def __init__(self,physVar):
        self.charge = physVar["Electron_charge"]
        self.pt = physVar["Electron_pt"]
        self.eta = physVar["Electron_eta"]
        self.phi = physVar["Electron_phi"]

       # self.dxy = detVal[0]
       # self.dxyErr = detVal[1]
       # self.dz = detVal[2]
       # self.dzErr = detVal[3]
       # 
       # self.idDict = idDict 
    
    def getPt(self):
        return self.pt
