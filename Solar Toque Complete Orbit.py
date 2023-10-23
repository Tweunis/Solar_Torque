import numpy as np
import math as m
from sympy import sin, cos, pi
import matplotlib.pyplot as plt

#%%
CR = 2
GscoverR2 = 590 # /m^2
c = 3*10**6 # m/s

dtheta = np.pi/1000
theta = np.arange(0, 2*np.pi, dtheta) # Angle in orbit

A_bus = 2.045 # m^2
r_b = 0.480 # m
r_a = 1.04938 # m
h = 1.050 # m
L_bus = 2.13 # m

def cross_prod(a, b):
    result = [a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]]

    return result

def add(a, b):
    result = [a[0]+b[0],
            a[1]+b[1],
            a[2]+b[2]]

    return result

#%% Rotating Force unit vector 

n_F2 = []
            
for i in theta:
    if  0 <= i < pi/2:
        n_F2.append([sin(i), 0, -cos(i)])
    elif pi/2 <= i < pi:
        n_F2.append([sin(i), 0, cos(i)])
    elif pi <= i < 3/2*pi:
        n_F2.append([sin(i), 0, cos(i)])
    elif 3/2*pi <= i <= 2*pi:
        n_F2.append([sin(i), 0, -cos(i)])
            
n_F = np.array(n_F2)


# %% Bus

r_bus = np.array([0.028367, -0.015365, -0.232788])

A_bus2 = []
r_bus2 = []
for i in theta:
    if 0 <= i < pi/2:
            i = -i
            y = sin(i)**2*(r_a-r_b/cos(i))+sin(i)**2*r_b/cos(i)-h*sin(i)
            x = sin(i)*(r_a-r_b/cos(i))+sin(i)*r_b/cos(i)-h
            if A_bus*sin(i) > 2*r_b*y:
                A_bus2.append(A_bus*sin(i) - 2*r_b* y)
                x_COA = -x/2 # delta z in catia
                if pi/2 - pi/100 < i < pi/2 + pi/100:
                    x_COA = 0
                if x_COA <= L_bus:
                    deltaRbus = np.array([0, 0, x_COA])
                    r_bus2.append(add(r_bus, deltaRbus))
                else:
                    r_bus2.append(r_bus)
            else:
                A_bus2.append(0)
                r_bus2.append(r_bus)
        
    elif pi/2 < i < pi:
            i = i-np.pi
            y = sin(i)**2*(r_a-r_b/cos(i))+sin(i)**2*r_b/cos(i)-h*sin(i)
            x = sin(i)*(r_a-r_b/cos(i))+sin(i)*r_b/cos(i)-h
            if A_bus*sin(i) > 2*r_b*y:
                A_bus2.append(A_bus*sin(i) - 2*r_b* y)
                x_COA = -x/2 # delta z in catia
                if pi/2 - pi/100 < i < pi/2 + pi/100:
                    x_COA = 0
                if x_COA <= L_bus:
                    deltaRbus = np.array([0, 0, x_COA])
                    r_bus2.append(add(r_bus, deltaRbus))
                else:
                    r_bus2.append(r_bus)
            else:
                A_bus2.append(0)
                r_bus2.append(r_bus)
        
    elif pi < i < 3/2*pi:
            i = np.pi-i
            y = sin(i)**2*(r_a-r_b/cos(i))+sin(i)**2*r_b/cos(i)-h*sin(i)
            x = sin(i)*(r_a-r_b/cos(i))+sin(i)*r_b/cos(i)-h
            if A_bus*sin(i) > 2*r_b*y:
                A_bus2.append(A_bus*sin(i) - 2*r_b* y)
                x_COA = -x/2 # delta z in catia
                if pi/2 - pi/100 < i < pi/2 + pi/100:
                    x_COA = 0
                if x_COA <= L_bus:
                    deltaRbus = np.array([0, 0, x_COA])
                    r_bus2.append(add(r_bus, deltaRbus))
                else:
                    r_bus2.append(r_bus)
            else:
                A_bus2.append(0)
                r_bus2.append(r_bus)
        
    elif 3/2*pi < i <= 2*pi:
            i = i-2*np.pi
            y = sin(i)**2*(r_a-r_b/cos(i))+sin(i)**2*r_b/cos(i)-h*sin(i)
            x = sin(i)*(r_a-r_b/cos(i))+sin(i)*r_b/cos(i)-h
            if A_bus*sin(i) > 2*r_b*y:
                A_bus2.append(A_bus*sin(i) - 2*r_b* y)
                x_COA = -x/2 # delta z in catia
                if pi/2 - pi/100 < i < pi/2 + pi/100:
                    x_COA = 0
                if x_COA <= L_bus:
                    deltaRbus = np.array([0, 0, x_COA])
                    r_bus2.append(add(r_bus, deltaRbus))
                else:
                    r_bus2.append(r_bus)
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

#%% Total torque per axis

tau_total_x = []
tau_total_y = []
tau_total_z = []
for i in range(len(theta)):
    tau_total_x.append((tau_array1[i][0] + tau_array2[i][0] + tau_bus[i][0] + tau_ant[i][0])*1000)
    tau_total_y.append((tau_array1[i][1] + tau_array2[i][1] + tau_bus[i][1] + tau_ant[i][1])*1000)
    tau_total_z.append((tau_array1[i][2] + tau_array2[i][2] + tau_bus[i][2] + tau_ant[i][2])*1000)

#%% Plotting

fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.plot(theta, tau_total_x, c="r",  label="$T_x$")
ax1.plot(theta, tau_total_z, c="b",  label="$T_z$")
ax2.plot(theta, tau_total_y, c="g",  label="$T_y$")

ax1.set_xlabel('True anomaly angle $θ$ [rad]')
ax1.set_ylabel('Solar Radiation Torque $T_s$ [mN$\cdot$m]')
ax2.set_xlabel('True anomaly angle $θ$ [rad]')
ax2.set_ylabel('Solar Radiation Torque $T_s$ [mN$\cdot$m]')
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

fig.savefig('Torque_Plots.pdf')













