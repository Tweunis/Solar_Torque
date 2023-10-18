import numpy as np
import math as m
from sympy import sin, cos, tan, pi
import matplotlib.pyplot as plt

#%%
CR = 2
GscoverR2 = 590 # /m^2
c = 3*10**6 # m/s

dtheta = np.pi/10
theta = np.arange(0, 2*np.pi + dtheta, dtheta) # Angle in orbit

A_bus = 2.045 # m^2
r_b = 0.480 # m
r_a = 1.04938 # m
h = 1.050 # m

def cross_prod(a, b):
    result = [a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]]

    return result

def rounding(A):
    for i in range(len(A)):
        for j in range(len(A[i])):
            if abs(A[i][j]) < 10**-15 :
                A[i][j] = 0

    return A

#%% Rotating Force unit vector 

n_F2 = []
for i in theta:
    if  0 < i < pi/2:
        for i in theta:
            n_F2.append([sin(i), 0, cos(i)])
    elif pi/2 < i < pi:
        for i in theta:
            n_F2.append([sin(i), 0, -cos(i)])
    elif pi < i < 3/2*pi:
        for i in theta:
            n_F2.append([sin(i), 0, -cos(i)])
    elif 3/2*pi < i < 2*pi:
        for i in theta:
            n_F2.append([sin(i), 0, cos(i)])

n_F = np.array(n_F2)


# %% Bus
phiQ2 = []
phiQ3 = []
phiQ4 = []
for i in theta:
    phiQ2.append(np.pi-i)
    phiQ3.append(i-np.pi)
    phiQ4.append(2*np.pi-i)
            
r_bus = np.array([0.028367, -0.015365, -0.232788])

A_bus2 = []
r_bus2 = []
for i in theta:
    if 0 < i < pi/2 or pi < i < 3/2*pi:
        for i in theta:
            y = r_a - r_b/cos(i) - h*sin(i) + r_b*tan(i)*sin(i)
            x = y / sin(i)
            if A_bus*sin(i) - 2*r_b* y >= 0:
                A_bus2.append(A_bus*sin(i) - 2*r_b* y)
                x_COA = -x/2 # delta z in catia
                deltaRbus = np.array([0, 0, x_COA])
                r_bus2.append(np.add(r_bus, deltaRbus))
            else:
                A_bus2.append(0)
                r_bus2.append(r_bus)
        
    elif pi/2 < i < pi:
        for i in phiQ2:
            y = r_a - r_b/cos(i) - h*sin(i) + r_b*tan(i)*sin(i)
            x = y / sin(i)
            if A_bus*sin(i) - 2*r_b* y >= 0:
                A_bus2.append(A_bus*sin(i) - 2*r_b* y)
                x_COA = -x/2 # delta z in catia
                deltaRbus = np.array([0, 0, x_COA])
                r_bus2.append(np.add(r_bus, deltaRbus))
            else:
                A_bus2.append(0)
                r_bus2.append(r_bus)
        
    elif pi < i < 3/2*pi:
        for i in phiQ3:
            y = r_a - r_b/cos(i) - h*sin(i) + r_b*tan(i)*sin(i)
            x = y / sin(i)
            if A_bus*sin(i) - 2*r_b* y >= 0:
                A_bus2.append(A_bus*sin(i) - 2*r_b* y)
                x_COA = -x/2 # delta z in catia
                deltaRbus = np.array([0, 0, x_COA])
                r_bus2.append(np.add(r_bus, deltaRbus))
            else:
                A_bus2.append(0)
                r_bus2.append(r_bus)
        
    elif 3/2*pi < i < 2*pi:
        for i in phiQ4:
            y = r_a - r_b/cos(i) - h*sin(i) + r_b*tan(i)*sin(i)
            x = y / sin(i)
            if A_bus*sin(i) - 2*r_b* y >= 0:
                A_bus2.append(A_bus*sin(i) - 2*r_b* y)
                x_COA = -x/2 # delta z in catia
                deltaRbus = np.array([0, 0, x_COA])
                r_bus2.append(np.add(r_bus, deltaRbus))
            else:
                A_bus2.append(0)
                r_bus2.append(r_bus)

F_bus = []
for i in range(len(theta)):
    F_bus.append(n_F[i] * CR * GscoverR2 / c * A_bus2[i])

tau_bus = []
for i in range(len(theta)):
    tau_bus.append(cross_prod(r_bus2[i], F_bus[i]))

#%% Antenna

A_ant = 3.47 # m^2
r_ant = np.array([0.028367, -0.015365, 1.827212])
F_ant0 = np.array(CR * GscoverR2 / c * A_ant)

F_ant = []
for i in range(len(theta)):
    F_ant.append(n_F[i] * F_ant0)

tau_ant = []
for i in range(len(theta)):
    tau_ant.append(cross_prod(r_ant, F_ant[i]))

#%% Arrays

A_array = 2.015 # m^2
r_array1 = np.array([0.028367, 1.765072, -0.344661])
r_array2 = np.array([0.028367, -1.795802, -0.344661])

F_array0 = np.array(CR * GscoverR2 / c * A_array)

F_array = []
for i in range(len(theta)):
    F_array.append(n_F[i] * F_array0)

tau_array1 = []
tau_array2 = []
for i in range(len(theta)):
    tau_array1.append(cross_prod(r_array1, F_array[i]))
    tau_array2.append(cross_prod(r_array2, F_array[i]))

#%% Total

tau_total_arrays = []
for i in tau_array1:
    for j in tau_array2:
        sum = [x + y for x, y in zip(i, j)]
        tau_total_arrays.append(sum)

tau_total_bus_ant = []
for i in tau_bus:
    for j in tau_ant:
        sum = [x + y for x, y in zip(i, j)]
        tau_total_bus_ant.append(sum)

tau_total = []
for i in tau_total_arrays:
    for j in tau_total_bus_ant:
        sum = [x + y for x, y in zip(i, j)]
        tau_total.append(sum)

tau_bus = rounding(tau_bus)
tau_ant = rounding(tau_ant)
tau_array1 = rounding(tau_array1)
tau_array2 = rounding(tau_array2)
tau_total = rounding(tau_total)

print(tau_total[0], tau_total[-1])

# Finding max values
lstx = []
lsty = []
lstz = []
for i in tau_total:
    lstx.append(i[0])
    lsty.append(i[1])
    lstz.append(i[2])
tau_x_max = max(lstx)
tau_y_max = max(lsty)
tau_z_max = max(lstz)

indexx = int(0)
for i in theta:
    if i - dtheta/2 < pi/2 < i + dtheta/2:
        indexy = int(i/dtheta)
        indexz = int(i/dtheta)

print("Maximum $𝜏_x$ =", round(lstx[indexx]*1000, 5), "mNm")    
print("Maximum $𝜏_y$ =", round(lsty[indexy]*1000, 5), "mNm")
print("Maximum $𝜏_z$ =", round(lstz[indexz]*1000, 5), "mNm")

#%% Plotting

tau_total_x = []
tau_total_y = []
tau_total_z = []
for i in range(len(theta)):
    tau_total_x.append(tau_total[i][0]*1000)
    tau_total_y.append(tau_total[i][1]*1000)
    tau_total_z.append(tau_total[i][2]*1000)

fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.plot(theta, tau_total_x, c="r", linestyle = "dashed", label="$𝜏_x$")
ax1.plot(theta, tau_total_z, c="b", linestyle = "dashed", label="$𝜏_z$")
ax2.plot(theta, tau_total_y, c="g", linestyle = "dashed", label="$𝜏_y$")

ax1.set_xlabel('True anomaly angle $θ$ [rad]')
ax1.set_ylabel('Torque $𝜏$ [mN$\cdot$m]')
ax2.set_xlabel('True anomaly angle $θ$ [rad]')
ax2.set_ylabel('Torque $𝜏$ [mN$\cdot$m]')
ax1.legend(loc='lower left')
ax2.legend(loc='lower left')

x_ticks = [0, m.pi/2, m.pi, 3/2*m.pi, 2*m.pi]
x_label = [0, '$\pi$/2', '$\pi$', '3/2$\cdot \pi$', '2$\cdot \pi$']

ax1.set_xticks(x_ticks) 
ax1.set_xticklabels(x_label, fontsize=10)
ax2.set_xticks(x_ticks) 
ax2.set_xticklabels(x_label, fontsize=10)

fig.tight_layout(pad=1.0)
ax1.grid()
ax2.grid()
plt.show()
plt.figure(dpi = 1200)

#fig.savefig('Torque_Plots.pdf')













