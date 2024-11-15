#provides a class to calculate coaxial jets with aerosandbox as its basis
#therefore it must be used within an optimizer, especially to solve the implicit equations

class CoaxJet:
    def __init__(self,Ui,Ua,Ti,Ta,d):
        #Ui is velocity of the Jet in [m/s]
        #Ua is velocity of the freestream [m/s]
        #Ti is temperature of the Jet [K]
        #Ta is temperature of the freestream [K]
        #d is diameter of the Jet
        self.Ui = Ui
        self.Ua = Ua
        self.Ti = Ti
        self.Ta = Ta
    
    def calcualte_lam(self,U):
        #calcualtes the lam value that corresponds to the given velocity
        return self.Ua / (U - self.Ua)

    def calculate_theta(self): #calculates the theta for this jet