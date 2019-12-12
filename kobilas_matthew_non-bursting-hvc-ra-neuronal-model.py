# sources from research author:
# http://www.personal.psu.edu/dzj2/PDFs/JinPRE2009.pdf
# http://www.personal.psu.edu/dzj2/PDFs/LongJinFeeNature2010.pdf
# http://www.personal.psu.edu/dzj2/PDFs/LongJinFeeNature2010SI.pdf

import math
import matplotlib.pyplot as plt
import numpy as np

C_m = 1 # microF/cm^2 membrane capacitance
A = 5000 # microm^2
# leak current for soma
def I_L(V):
    G_L = 0.1 # mS/cm^2 conductance
    E_L = -80 # mV reversal potential
    return -G_L * (V - E_L) # leak current
# sodium current
def I_Na(V, h):
    G_Na = 60 # mS/cm^2 conductance
    E_Na = 55 # mV reversal potential
    return -G_Na * m_infinity(V)**3 * h * (V - E_Na) # sodium current
# potassium current
def I_Kdr(V, n):
    G_Kdr = 8 # mS/cm^2 conductance
    E_K = -90 # mV reversal potential
    return -G_Kdr * n**4 * (V - E_K) # potassium current
def I_exc(V, g_exc):
    # g_exc is the total excitatory synaptic conductance
    return -g_exc * V
def I_inh(V, g_inh):
    E_I = -80
    return -g_inh * (V - E_I)

# the general equation for the gating variables m, h, and n is
# dx/dt = alpha_x(V)(1 - x) - beta_x(V)x
# where x=m, h, n
# the voltage dependent coefficients of the gating variables are
def gatingVarHN(x_infinity, tau_x, x, V):
    return (x_infinity(V) - x) / tau_x(V)
def m_infinity(V):
        return 1 / (1 + math.exp(-(V + 30) / 9.5))
def h_infinity(V):
    return 1 / (1 + math.exp((V + 45) / 7))
def tau_h(V):
    return 0.1 + 0.75 / (1 + math.exp((V + 40.5) / 6))
def n_infinity(V):
    return 1 / (1 + math.exp(-(V + 35) / 10))
def tau_n(V):
    return 0.1 + 0.5 / (1 + math.exp((V + 27) / 15))

def getTotalCurrent(V, h, n, g_exc, g_inh, I_ext):
    return I_L(V) + I_Na(V, h) + I_Kdr(V, n) + I_exc(V, g_exc) + I_inh(V, g_inh) + I_ext / A
def getConductance(g):
    return -g / tau

# 500 ms is also too long but more clear
# 40 ms is good, but incorrect behavior
T_max = 50
dt = 0.01
t = np.arange(0, T_max, dt).tolist()
V = [0 for i in range(len(t))]
h = [0 for i in range(len(t))]
n = [0 for i in range(len(t))]
g_exc = [0 for i in range(len(t))]
g_inh = [0 for i in range(len(t))]
I_ext = [0 for i in range(len(t))]
start = 1000
end = 3000
I_ext[start:end] = [50000] * (end - start)
#I_ext[start:end] = [26540] * (end - start)
# generates 1 spike at 18778 nA
# generates 2 spikes at 26540 nA
G = 0.1
tau_exc = 5 # ms
tau_inh = 5 # ms
tau = tau_exc
Vrst = -80
Vth = -55
Vspk = 60
V[0] = Vrst

for i in range(len(t) - 1):
    kV0 = getTotalCurrent(V[i], h[i], n[i], g_exc[i], g_inh[i], I_ext[i]) / C_m
    aV = V[i] + kV0 * dt
    kV1 = getTotalCurrent(aV, h[i+1], n[i+1], g_exc[i+1], g_inh[i+1], I_ext[i+1]) / C_m
    V[i+1] = V[i] + (kV0 + kV1) * dt / 2
    
    kgexc0 = getConductance(g_exc[i])
    agexc = g_exc[i] + kgexc0 * dt
    kgexc1 = getConductance(agexc)
    g_exc[i+1] = g_exc[i] + (kgexc0 + kgexc1) * dt / 2
    
    kginh0 = getConductance(g_inh[i])
    aginh = g_inh[i] + kginh0 * dt
    kginh1 = getConductance(aginh)
    g_inh[i+1] = g_inh[i] + (kginh0 + kginh1) * dt / 2
    
    kh0 = gatingVarHN(h_infinity, tau_h, h[i], V[i])
    ah = h[i] + kh0 * dt
    kh1 = gatingVarHN(h_infinity, tau_h, ah, V[i+1])
    h[i+1] = h[i] + (kh0 + kh1) * dt / 2
    
    kn0 = gatingVarHN(n_infinity, tau_n, n[i], V[i])
    an = n[i] + kn0 * dt
    kn1 = gatingVarHN(n_infinity, tau_n, an, V[i])
    n[i+1] = n[i] + (kn0 + kn1) * dt / 2

    if V[i+1] >= Vth:
        V[i] = Vspk
        V[i+1] = Vrst
        g_exc[i+1] = g_exc[i] + G

plt.figure(1)
plt.plot(t, V, 'r--')
plt.axis([0, 50, -80, 60])
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mV)')
plt.title('Non-bursting model of HVC(RA) neurons')
plt.show()

plt.figure(2)
plt.plot(t, h, 'b--', t, n, 'g--')
plt.xlabel('Time (ms)')
plt.ylabel('Gating Variables')
plt.title('Time vs. Gating Variable Value')
plt.show()