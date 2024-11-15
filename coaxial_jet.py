#provides a class to calculate coaxial jets with aerosandbox as its basis
#therefore it must be used within an optimizer, especially to solve the implicit equations
import aerosandbox.numpy as np
import aerosandbox as asb

class CoaxJet:
    def __init__(self,Ui,Ua,Ti,Ta,d,opti: asb.Opti):
        #Ui is velocity of the Jet in [m/s]
        #Ua is velocity of the freestream [m/s]
        #Ti is temperature of the Jet [K]
        #Ta is temperature of the freestream [K]
        #d is diameter of the Jet
        #opti is the optimizer class these calculations will be "appended to"
        self.Ui = Ui
        self.Ua = Ua
        self.Ti = Ti
        self.Ta = Ta
        self.d = d
        self.lam_0 = self.calculate_lam(Ui)
        self.theta = self.calculate_theta()
        self.opti = opti 
        self.x_core = self.calculate_x_core()
    
    def calculate_lam(self,Uc):
        #calcualtes the lam value that corresponds to the given velocity
        return self.Ua / (Uc - self.Ua)

    def calculate_Uc_from_lam(self,lam):
        return self.Ua *((1/lam) + 1)
    
    def calculate_x_core(self):
        #calculates the xv at which the core region (aka. constant Uc) ends
        return self.calculate_xv(self.Ui)

    def calculate_theta(self): #calculates the theta for this jet
         # Calculates theta for this jet
        lam_0 = self.lam_0  # Calculate lambda_0 for the jet velocity Ui
        theta_d_ratio = np.sqrt((np.pi / 4) * (self.Ti / self.Ta) * ((lam_0 + 1) / (lam_0 ** 2)))
        theta = self.d * theta_d_ratio
        return theta
    
    def calculate_xv(self, Uc):
        # Calculate lambda for the given U
        lam = self.calculate_lam(Uc)
        # Calculate Lambda
        Lambda = lam + 0.468
        # Calculate xv/theta term
        xv_theta_ratio = (
            6.39 * Lambda ** (3 / 2) +
            5.38 * Lambda ** (1 / 2) +
            1.76 * Lambda ** (-1 / 2) -
            8.3
        )
        # Calculate theta
        theta = self.calculate_theta()
        # Calculate xv
        xv = theta * xv_theta_ratio
        return xv
    
    def calculate_Uc_from_xv(self,x_v):
        #implicitly calculates the responding centerline U_c value to a x_v value
        l_U = self.opti.variable(self.Ua + 0.5 * (self.Ui - self.Ua))

        #if we are within the core, setting l_U to Ui must return x_v
        #if wer are smaller than x_core blend should be very negative and we should get the "first as output",
        #where the equation is fullfiled if l_U equals Ui
        self.opti.subject_to(x_v == np.blend(9999**4 * (x_v - self.x_core) ,self.calculate_xv(l_U)  ,x_v * (self.Ui / l_U)))
        return l_U