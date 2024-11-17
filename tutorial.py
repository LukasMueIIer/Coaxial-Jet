#usage of the coaxial_jet.py class
import aerosandbox as asb
import aerosandbox.numpy as np
from coaxial_jet import *
import matplotlib.pyplot as plt
import copy

#parameters of the Jet
Ui = 10     #Jet velocity [m/s]
Ua = 5      #Freestream Velocity [m/s]
Ti = 450    #Temperature of the Jet [K]
Ta = 225    #Surrounding Temperature [K]
d = 0.1     #Jet diameter [m]
rho_c_i = 0.01  #mass concentration of the contaminant in jet [kg/m^3]
rho_c_a = 0     #mass concentration of the contaminant on freestream [kg/m^3]

opti = asb.Opti() #Creating our optimizer

jet = CoaxJet(Ui,Ua,Ti,Ta,rho_c_i,rho_c_a,d,opti) #creating the Jet class, passing it the jet data and optimizer

## Center Line calculations
x_c = jet.calculate_x_core()
print("Core Length of the Jet: "+ str(x_c))

X = np.linspace(0,20,50)
U_c = []
T_c = []
for i in X:
    #U_c.append(jet.calculate_Uc_from_xv(i))
    #T_c.append(jet.calculate_dTc_ratio(i))
    U_c.append(jet.calculate_delta_U(i))
    T_c.append(jet.calculate_delat_T(i))



sol = opti.solve()

U_c = sol(U_c)
T_c = sol(T_c)
print(U_c)
print(T_c)

# Create the plot
plt.figure(figsize=(10, 6))

# Plot U_c
plt.plot(X, U_c, label=r"$U_c$", color="blue", linestyle="-", marker="o")

# Plot T_c
plt.plot(X, T_c, label=r"$T_c$", color="red", linestyle="--", marker="x")

# Add labels and legend
plt.xlabel("X (m)")
plt.ylabel("Values")
plt.title("U_c and T_c Profiles Along X")
plt.legend()
plt.grid(True)

# Show the plot
plt.show()