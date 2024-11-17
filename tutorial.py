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

jet = CoaxJet(Ui,Ua,Ti,Ta,d,opti) #creating the Jet class, passing it the jet data and optimizer

## Center Line calculations
x_c = jet.calculate_x_core()
print("Core Length of the Jet: "+ str(x_c))

X = np.linspace(0,20,50)
U_c = []
T_c = []
for i in X:
    U_c.append(jet.calculate_Uc_from_xv(i))
    T_c.append(jet.calculate_dTc_ratio(i))



sol = opti.solve()

U_c = sol(U_c)
T_c = sol(T_c)
print(U_c)
print(T_c)
# Create the plot
fig, ax1 = plt.subplots(figsize=(10, 6))  # Create figure and first axis

# Plot U_c on the primary y-axis (left)
ax1.plot(X, U_c, label="Centerline Velocity $U_c$", color="blue", marker="o")
ax1.set_xlabel("x (m)")
ax1.set_ylabel("$U_c$ (m/s)", color="blue")
ax1.tick_params(axis="y", labelcolor="blue")
ax1.set_ylim(jet.Ua,jet.Ui)
ax1.grid(True)

# Add secondary y-axis (right) for T_c
#ax2 = ax1.twinx()  # Create a second y-axis sharing the same x-axis
#ax2.plot(x_stats, T_c, label="Centerline Temperature $T_c$", color="red", linestyle="--", marker="x")
#ax2.set_ylabel("$T_c$ (K)", color="red")
#ax2.tick_params(axis="y", labelcolor="red")
#ax2.set_ylim(290, 360)  # Adjust limits as needed

# Add legends
fig.legend(loc="upper center", bbox_to_anchor=(0.5, 1.15), ncol=2)

# Add a title
plt.title("Centerline Velocity and Temperature Profile")

# Show the plot
plt.show()
