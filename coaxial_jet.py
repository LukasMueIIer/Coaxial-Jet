#provides a class to calculate coaxial jets with aerosandbox as its basis
#therefore it must be used within an optimizer, especially to solve the implicit equations
import aerosandbox.numpy as np
import aerosandbox as asb

class CoaxJet:
    def __init__(self,Ui,Ua,Ti,Ta,rho_c_i,rho_c_a,d,opti: asb.Opti):
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
        self.rho_c_i = rho_c_i
        self.rho_c_a = rho_c_a
        self.d = d
        self.lam_0 = self.calculate_lam(Ui)
        self.theta = self.calculate_theta()
        self.opti = opti 
        self.x_core = self.calculate_x_core()
        self.delta_U_core = self.calculate_delta_U_full_Virtual(self.x_core)
    
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
        l_U = self.opti.variable(self.Ua + 0.5 * (self.Ui - self.Ua),log_transform=True,upper_bound=self.Ui,lower_bound=self.Ua)
        #if we are within the core, setting l_U to Ui must return x_v
        #if wer are smaller than x_core blend should be very negative and we should get the "first as output",
        #where the equation is fullfiled if l_U equals Ui
        self.opti.subject_to(x_v == np.blend(9999**4 * (x_v - self.x_core) ,self.calculate_xv(l_U)  ,x_v + l_U - self.Ui))
        return l_U
    
    def calculate_delta_U_full_Virtual(self, x_v): #no blending for the core
        # Calculates delta for a given x_v
        Phi_41 = 0.0950
        Phi_51 = 0.0445

        # Calculate Uc from x_v
        Uc = self.calculate_Uc_from_xv(x_v)
        # Calculate lambda
        lam = self.calculate_lam(Uc)
        # Calculate delta/theta
        delta_theta_ratio = lam / np.sqrt(2 * np.pi * (lam * Phi_41 + Phi_51))
        # Calculate theta
        theta = self.calculate_theta()
        # Calculate delta
        delta = delta_theta_ratio * theta
        return delta
    
    def calculate_delta_U(self,x_v):
        return np.blend(9999**4 * (self.x_core - x_v), self.d/2 + x_v / self.x_core * (self.delta_U_core - 0.5 * self.d), self.calculate_delta_U_full_Virtual(x_v))
    
    def calculate_dTc_ratio_from_Uc(self,U_c):

        E = 2  # Constant value of E
        
        # Calculate terms
        ratio_Uc = (U_c - self.Ua) / (self.Ui - self.Ua)
        ln_term_numerator = ratio_Uc * (1 + 2 * (self.Ua / (self.Ui - self.Ua)))
        ln_term_denominator = ratio_Uc + 2 * (self.Ua / (self.Ui - self.Ua))
        
        # Calculate the main equation
        dTc_ratio = (1 / E) * (1 / (
            (1 / ratio_Uc - 1) +
            (1 / 6) * ((self.Ui - self.Ua) / self.Ua) * np.log(ln_term_numerator / ln_term_denominator) +
            (1 / E)
        ))
        return dTc_ratio
    
    def calculate_dTc_ratio(self,x_v):
        U_c = self.calculate_Uc_from_xv(x_v)
        return self.calculate_dTc_ratio_from_Uc(U_c)
    
    def calculate_Tc(self,x_v):
        #calculates the core Temperature for a given x_v
        deltaTi = self.Ti - self.Ta
        l_deltaTi = self.calculate_dTc_ratio(x_v) * deltaTi
        return self.Ta + l_deltaTi

    def calculate_dRho_ratio(self,x_v):
        return self.calculate_dTc_ratio(x_v)
    
    def calculate_rho_c_c(self,x_v):
        deltaRho_c = self.rho_c_i - self.rho_c_a
        l_deltaRho_c = self.calculate_dRho_ratio(x_v)
        return self.rho_c_a + l_deltaRho_c
    
    def calculate_delat_T(self,x_v):
        l_U = self.calculate_Uc_from_xv(x_v)
        dTc_ratio = self.calculate_dTc_ratio_from_Uc(l_U)
        ratio_Uc = (l_U - self.Ua) / (self.Ui - self.Ua)
        delta_ratio = np.sqrt(ratio_Uc) / np.sqrt(dTc_ratio)
        return delta_ratio * self.calculate_delta_U(x_v)
