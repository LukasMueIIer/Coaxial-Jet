#usage of the coaxial_jet.py class
import aerosandbox as asb
from coaxial_jet import *

Ui = 10 #[m/s]
Ua = 5 #[m/s]
d = 0.1 #[m]

opti = asb.Opti()

jet = CoaxJet(Ui,Ua,1,1,d,opti) #creating the Jet class

x_c = jet.calculate_x_core()
print(x_c)
U = jet.calculate_Uc_from_xv(1)

sol = opti.solve()
print(sol(U))