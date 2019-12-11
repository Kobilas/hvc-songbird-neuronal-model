import math
import matplotlib.pyplot as plt
import numpy as np

C_m = 1 # microF/cm^2 membrane capacitance
A_s = 100 # microm^2 surface area of soma
A_d = 50000 # microm^2 surface area of dendrite
# leak current for soma
def I_Ls(V_s):
    g_Ls = 0.05 # mS/cm^2 conductance
    E_r = -85 # mV reversal potential
    return g_Ls * (E_r - V_s) # leak current
# sodium current
def I_Na(V_s, m, h):
    g_Na = 100 # mS/cm^2 conductance
    E_Na = 55 # mV reversal potential
    return g_Na * m**3 * h * (E_Na - V_s) # sodium current
# potassium current
def I_K(V_s, n):
    g_K = 2 # mS/cm^2 conductance
    E_K = -90 # mV reversal potential
    return g_K * n**4 * (E_K - V_s) # potassium current
# high threshold potassium current
def I_KHT(V_s, w):
    g_KHT = 300 # mS/cm^2 conductance
    E_K = -90 # mV reversal potential
    return g_KHT * w * (E_K - V_s) # high threshold potassium current
# low threshold potassium current
def I_KLT(V_s, l):
    g_KLT = 25 # mS/cm^2 conductance
    E_K = -90 # mV reversal potential
    return g_KLT * l * (E_K - V_s) # low threshold potassium current
# feedforward excitatory input to soma
def I_FFs(V_s, g_FFs):
    return -g_FFs * V_s
R_c = 250 # MOhms resistance of connection between soma and dendrite
# leak current for dendrite
def I_Ld(V_d):
    g_Ld = 0.1 # mS/cm^2 conductance
    E_r = -85 # mV reversal potential
    return g_Ld * (E_r - V_d) # leak current
# high threshold calcium current
def I_Ca(V_d, m_infinity):
    g_Ca = 200 # mS/cm^2 conductance
    E_Ca = 120 # mV reversal potential
    return g_Ca * m_infinity**2 * (E_Ca - V_d) # high threshold calcium current
# calcium dependent potassium current
def I_CaK(V_d, q):
    g_CaK = 100 # mS/cm^2 conductance
    E_K = -90 # mV reversal potential
    return g_CaK * q * (E_K - V_d) # calcium dependent potassium current
# excitatory synaptic current
def I_syn(V_d, g_syn):
    return -g_syn * V_d # excitatory synaptic current
def getG_syn(g_syn):
    return -g_syn / tau_syn
# calcium concentration follows first-order kinetics:
# dConc_Ca/dt = 0.1 * I_Ca - Conc_Ca/tau_Ca;
def getCaConc(V_d, m_infinity, conc_Ca):
    tau_Ca = 100 # ms decay time constant
    return 0.1 * I_Ca(V_d, m_infinity) - conc_Ca / tau_Ca
tau_Ca = 100 # ms decay time constant
# feedforward excitatory input to dendrite
def I_FFd(V_d, g_FFd):
    return -g_FFd * V_d # feedforward excitatory input to dendrite

# the general equation for the gating variables m, h, and n is
# dx/dt = alpha_x(V)(1 - x) - beta_x(V)x
# where x=m, h, n
# the voltage dependent coefficients of the gating variables are
def gatingVarMHN(alpha_x, beta_x, x, V):
    return alpha_x(V) * (1 - x) - beta_x(V) * x
def alpha_m(V):
    try:
        res = -0.5 * (V + 22) / (math.exp(-(V + 22) / 10) - 1)
        return res
    except OverflowError:
        print('*******************' + str(V))
def beta_m(V):
    return 20 * math.exp(-(V + 47) / 18)
def m_infinity(V_d):
    try:
        res = 1 / (1 + math.exp(-(V_d - 20) / 15))
        return res
    except OverflowError:
        print('*****************' + str(V_d))
def alpha_h(V):
    return 0.35 * math.exp(-(V + 34) / 20)
def beta_h(V):
    return 5 / (math.exp(-(V + 4) / 10) + 1)
def alpha_n(V):
    return -0.075 * (V + 30) / (math.exp(-(V + 30) / 10) - 1)
def beta_n(V):
    return 0.1 * math.exp(-(V + 40) / 80)
# the general equation for the gating variables w, l, and q is
# dx/dt = (x_infinity(V) - x) / tau_x
def gatingVarWLQ(x_infinity, tau_x, x, V):
    return (x_infinity(V) - x) / tau_x
def w_infinity(V):
    return 1 / (math.exp(-V / 5) + 1)
tau_w = 1 # ms
def l_infinity(V):
    return 1 / (math.exp(-(V + 40) / 5) + 1)
tau_l = 10 # ms
def q_infinity(conc_Ca):
    return (0.0005 * conc_Ca)**2
def tau_q(conc_Ca):
    return 0.0338 / (min(0.0001 * conc_Ca,  0.01) +  0.001)

def getTotalSomaCurrent(V_s, m, h, n, w, l, g_FFs):
    ils = I_Ls(V_s)
    ina = I_Na(V_s, m, h)
    ik = I_K(V_s, n)
    ikht = I_KHT(V_s, w)
    iklt = I_KLT(V_s, l)
    iffs = I_FFs(V_s, g_FFs)
    #print('ils: ' + str(ils) + '; ina: ' + str(ina) + '; ik: ' + str(ik) + '; ikht: ' + str(ikht) + '; iklt: ' + str(iklt) + '; iffs: ' + str(iffs))
    return I_Ls(V_s) + I_Na(V_s, m, h) + I_K(V_s, n) + I_KHT(V_s, w) + I_KLT(V_s, l) + I_FFs(V_s, g_FFs)
def getTotalDendriteCurrent(V_d, m_infinity, q, g_syn, g_FFd):
    return I_Ld(V_d) + I_Ca(V_d, m_infinity) + I_CaK(V_d, q) + I_syn(V_d, g_syn) + I_FFd(V_d, g_FFd)

T_max = 1000
dt = 0.01
t = np.arange(0, T_max, dt).tolist()
Vs = [0 for i in range(len(t))]
Vd = [0 for i in range(len(t))]
m = [0 for i in range(len(t))]
h = [0 for i in range(len(t))]
n = [0 for i in range(len(t))]
w = [0 for i in range(len(t))]
l = [0 for i in range(len(t))]
q = [0 for i in range(len(t))]
conc_ca = [0 for i in range(len(t))]
g_syn = [0 for i in range(len(t))]
g_ffd = 1 # mS/cm^2
g_ffs = 1 # mS/cm^2
I_ext = 0 # nA
G = 0.1
#g_syn = 1
tau_syn = 5 # ms
Vrst = -85
Vd[0] = Vrst
Vs[0] = Vrst
Vth = -60

for i in range(len(t) - 1):
    kg0 = getG_syn(g_syn[i])
    kVs0 = ((A_s * getTotalSomaCurrent(Vs[i], m[i], h[i], n[i], w[i], l[i], g_ffs)) + I_ext + (Vd[i] - Vs[i]) / R_c) / (C_m * A_s)
    ag = g_syn[i] + kg0 * dt
    aVs = Vs[i] + kVs0 * dt
    kg1 = getG_syn(ag)
    kVs1 = ((A_s * getTotalSomaCurrent(aVs, m[i+1], h[i+1], n[i+1], w[i+1], l[i+1], g_ffs)) + I_ext + (Vd[i+1] - aVs) / R_c) / (C_m * A_s)
    g_syn[i+1] = g_syn[i] + (kg0 + kg1) * dt / 2
    Vs[i+1] = Vs[i] + (kVs0 + kVs1) * dt / 2
    
    #print('loop i: ' + str(i) + '; Vs[i]: ' + str(Vs[i]) + '; Vs[i+1]: ' + str(Vs[i+1]))
    if Vs[i] >= Vth:
        #print('hit')
        g_syn[i+1] = g_syn[i] + G
        Vs[i+1] = Vrst
    
    kVd0 = ((A_d * getTotalDendriteCurrent(Vd[i], m_infinity(Vd[i]), q[i], g_syn[i], g_ffd)) + (Vs[i] - Vd[i]) / R_c) / (C_m * A_d)
    aVd = Vd[i] + kVd0 * dt
    kVd1 = ((A_d * getTotalDendriteCurrent(aVd, m_infinity(aVd), q[i+1], g_syn[i+1], g_ffd)) + (Vs[i+1] - aVd) / R_c) / (C_m * A_d)
    Vd[i+1] = Vd[i] + (kVd0 + kVd1) * dt / 2
    
    kCa0 = getCaConc(Vd[i], m_infinity(Vd[i]), conc_ca[i])
    aCa = conc_ca[i] + kCa0 * dt
    kCa1 = getCaConc(Vd[i+1], m_infinity(Vd[i+1]), aCa)
    conc_ca[i+1] = conc_ca[i] + (kCa0 + kCa1) * dtertte / 2
    
    km0 = gatingVarMHN(alpha_m, beta_m, m[i], Vs[i])
    am = m[i] + km0 * dt
    km1 = gatingVarMHN(alpha_m, beta_m, am, Vs[i+1])
    m[i+1] = m[i] + (km0 + km1) * dt / 2
    
    kh0 = gatingVarMHN(alpha_h, beta_h, h[i], Vs[i])
    ah = h[i] + kh0 * dt
    kh1 = gatingVarMHN(alpha_h, beta_h, ah, Vs[i+1])
    h[i+1] = h[i] + (kh0 + kh1) * dt / 2
    
    kn0 = gatingVarMHN(alpha_n, beta_n, n[i], Vs[i])
    an = n[i] + kn0 * dt
    kn1 = gatingVarMHN(alpha_n, beta_n, an, Vs[i+1])
    n[i+1] = n[i] + (kn0 + kn1) * dt / 2
    
    kw0 = gatingVarWLQ(w_infinity, tau_w, w[i], Vs[i])
    aw = w[i] + kw0 * dt
    kw1 = gatingVarWLQ(w_infinity, tau_w, aw, Vs[i])
    w[i+1] = w[i] + (kw0 + kw1) * dt / 2
    
    kl0 = gatingVarWLQ(l_infinity, tau_l, l[i], Vs[i])
    al = l[i] + kl0 * dt
    kl1 = gatingVarWLQ(l_infinity, tau_l, al, Vs[i+1])
    l[i+1] = l[i] + (kl0 + kl1) * dt / 2
    
    kq0 = gatingVarWLQ(q_infinity, tau_q(conc_ca[i]), q[i], Vd[i])
    aq = q[i] + kq0 * dt
    kq1 = gatingVarWLQ(q_infinity, tau_q(conc_ca[i+1]), aq, Vd[i+1])
    q[i+1] = q[i] + (kq0 + kq1) * dt / 2
        
plt.plot(t, Vs, 'r--', t, Vd, 'b--')
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mV)')