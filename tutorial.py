#usage of the coaxial_jet.py class
import aerosandbox as asb
from coaxial_jet import *

#parameters of the Jet
Ui = 10     #Jet velocity [m/s]
Ua = 5      #Freestream Velocity [m/s]
Ti = 450    #Temperature of the Jet [K]
Ta = 225    #Surrounding Temperature [K]
d = 0.1     #Jet diameter [m]

opti = asb.Opti() #Creating our optimizer

jet = CoaxJet(Ui,Ua,Ti,Ta,d,opti) #creating the Jet class, passing it the jet data and optimizer

## Center Line calculations
x_c = jet.calculate_x_core()
print("Core Length of the Jet: "+ str(x_c))

U_c_1 = jet.calculate_Uc_from_xv(1)     #Get the center line velocity, this is an implicit calculation
sol = opti.solve()                      #therefore result is only avaliable after running the solver
print("Centerline Velicity at 1m: " + str(sol(U_c_1)))
print(sol(jet.delta_U_core))



