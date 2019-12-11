% membrane potentials of V_s and V_d obey following equations
% C_m * A_s * (dV_s/dt) = A_s(I_Ls + I_Na + I_K + I_KHT + I_KLT + I_FFs)
%   + I_ext + (V_d - V_s) / R_c
% C_m * A_d * (dV_d/dt) = A_d(I_Ld + I_Ca + I_CaK + I_syn + I_FFd)
%   + (V_s - V_d) / R_c

% soma constants/equations
C_m = 1; % ?F/cm^2, membrance capacitance
A_s = 100; % ?m^2, soma surface area
A_d = 50000; % ?m^2, dendrite surface area
I_Ls = @(V_s) 0.05 .* (-85 - V_s); % leak current for soma, constants below
% g_Ls = 0.05 mS/cm^2, conductance
% E_r = -85 mV, reversal potential
I_Na = @(V_s, m, h) 100 .* m^3 .* h .* (55 - V_s); % sodium current, constants below
% g_Na = 100 mS/cm^2, conductance
% E_Na = 55 mV, reversal potential
% gating variables m and h
I_K = @(V_s, n) 2 .* n^4 .* (-90 - V_s); % potassium current, constants below
% g_K = 2 mS/cm^2, conductance
% E_K = -90 mV, reversal potential
% gating variable n
I_KHT = @(V_s, w) 300 .* w .* (-90 - V_s); % high threshold potassium current, constants below
% g_KHT = 300 mS/cm^2, conductance
% E_K = -90 mV, reversal potential
% gating variable w
I_KLT = @(V_s, l) 25 .* l .* (-90 - V_s); % low threshold potassium current, constants below
% g_KLT = 25 mS/cm^2, conductance
% E_K = -90 mV, reversal potential
% gating variable l
I_FFs = @(V_s, g_FFs) -g_FFs .* V_s; % feedforward excitatory input to soma, constants below
% g_FFs = constant conductance
% I_ext, external current injection to soma
R_c = 250; % M?, resistance of connection between soma and dendrite

% dendrite constants/equations
I_Ld = @(V_d) 0.1 .* (-85 - V_d); % leak current for dendrite, constants below
% g_Ld = 0.1 mS/cm^2, conductance
% E_r = -85 mV, reversal potential
I_Ca = @(V_d, m_infin) 200 .* m_infin^2 .* (120 - V_d); % high threshold calcium current, constants below
% g_Ca = 200 mS/cm^2, conductance
% E_Ca = 120 mV, reversal potential
% voltage dependant factor m_infin
I_CaK = @(V_d, q) 100 .* q .* (-90 - V_d); % calcium dependent potassium current, constants below
% g_CaK = 100 mS/cm^2, conductance
% E_K = -90 mV, reversal potential
% calcium dependent variable q
I_syn = @(V_d, g_syn) -g_syn .* V_d; % excitatory synaptic current, constants below
% Calcium concentration follows a first order kinetics:
% dConc_Ca/dt = 0.1 * I_Ca - Conc_Ca/tau_Ca;
tau_Ca = 100; % ms, decay time constant
I_FFd = @(V_d, g_FFd) -g_FFd .* V_d; % feedforward excitatory to the dendrite, constants below
% g_FFd = constant conductance

% Synaptic conductance follows a "kick-and-decay" kinetics: g_syn -> g_syn + G
% when a spike arrives from another HVC(RA) neuron at a synapse with conductance G
% and dg_syn/dt = -g_syn/tau_syn in between spikes with synaptic time constant
% tau_syn = 5ms

% Gating variable equations/constants
% General equation for gating variables x = m, h, n:
% dx/dt = alpha_x(V)(1 - x) - beta_x(V)x
alpha_m = @(V) -0.5 .* (V + 22) ./ (exp(-(V + 22) ./ 10) - 1);
beta_m = @(V) 20 .* exp(-(V + 47) ./ 18);
m_infinity = @(V) 1 ./ (1 + exp(-(V - 20) ./ 15));
alpha_h = @(V) 0.35 .* exp(-(V + 34) ./ 20);
beta_h = @(V) 5 ./ (exp(-(V + 4) ./ 10) + 1);
alpha_n = @(V) (0.032.*(V+52))./(1-(exp(-(V+52)./5)));
beta_n = @(V) (0.5.*(exp(-(57+V)./(40))));
w_infinity = @(V) 1 ./ (exp(-V ./ 5) + 1);
tau_w = 1; % ms
l_infinity = @(V) 1 ./ (exp(-(V + 40) ./ 5) + 1);
tau_l = 10; % ms
q_infinity = @(conc_Ca) (0.0005 .* conc_Ca)^2;
tau_q = @(conc_Ca) (0.0338) ./ (min([(0.0001 .* conc_Ca) 0.01]) + 0.0001);
% General equation for gating variables x = w, l, q
% dx/dt = (x_infinity(V) - x)/tau_x

Tmax = 1000;
dt = 0.1;
t = 0:dt:Tmax;
Vs = zeros(0, length(t));
Vd = zeros(0, length(t));
m = zeros(0, length(t));
h = zeros(0, length(t));
n = zeros(0, length(t));
w = zeros(0, length(t));
l = zeros(0, length(t));
q = zeros(0, length(t));
ca_conc = zeros(0, length(t));
g_ffd = 1; % mS/cm^2
g_ffs = 1; % mS/cm^2
I_ext = 1; % nA

% TODO: add I_syn, look at homework 3 g_sra (spike rate adaptation)
for i=1:length(t) - 1
    kVs0 = (A_s .* (I_Ls(Vs(i)) + I_Na(Vs(i), m(i), h(i)) ...
        + I_K(Vs(i), n(i)) + I_KHT(Vs(i), w(i)) ...
        + I_KLT(Vs(i), l(i)) + I_FFs(Vs(i), g_ffs)) ...
        + I_ext + (Vd(i) - Vs(i)) ./ R_c) ./ (C_m .* A_s);
    aVs = Vs(i) + kVs0 * dt;
    kVs1 = (A_s .* (I_Ls(aVs) + I_Na(aVs, m(i + 1), h(i + 1)) ...
        + I_K(aVs, n(i + 1)) + I_KHT(aVs, w(i + 1)) ...
        + I_KLT(aVs, l(i + 1)) + I_FFs(aVs, g_ffs)) ...
        + I_ext + (Vd(i + 1) - aVs) ./ R_c)/(C_m .* A_s);
    Vs(i + 1) = Vs(i) + (kVs0 + kVs1)*dt/2;
    kVd0 = (A_d .* (I_Ld(Vd(i)) + I_Ca(Vd(i), m_infinity(Vd(i))) ...
        + I_CaK(Vd(i), q(i)) + I_syn(Vd(i), ) ...
        + I_FFd(Vd(i), g_ffd)) + (Vs(i) - Vd(i)) ./ R_c) ./ (C_m .* A_d);
    aVd = Vd(i) + kVd0 * dt;
    kVd1 = (A_d .* (I_Ld(aVd) + I_Ca(aVd, m_infinity(aVd)) ...
        + I_CaK(aVd, q(i + 1)) + I_syn(aVd, ) ...
        + I_FFd(aVd, g_ffd)) + (Vs(i + 1) - aVd) ./ R_c) ./ (C_m .* A_d);
    Vd(i + 1) = Vd(i) + (kVd0 + kVd1)*dt/2;
    kCa0 = 0.1 .* I_Ca(Vd(i), m_infinity(Vd(i))) - ca_conc(i) ./ tau_Ca;
    aCa = ca_conc(i) + kCa0 * dt;
    kCa1 = 0.1 .* I_Ca(Vd(i + 1), m_infinity(Vd(i + 1))) - aCa ./ tau_Ca;
    ca_conc(i + 1) = ca_conc(i) + (kCa0 + kCa1) * dt/2;
    km0 = alpha_m(Vs(i)) .* (1 - m(i)) - beta_m(Vs(i)) .* m(i);
    am = m(i) + km0 * dt;
    km1 = alpha_m(Vs(i + 1)) .* (1 - am) - beta_m(Vs(i + 1)) .* am;
    m(i + 1) = m(i) + (km0 + km1) * dt/2;
    kh0 = alpha_h(Vs(i)) .* (1 - h(i)) - beta_h(Vs(i)) .* h(i);
    ah = h(i) + kh0 * dt;
    kh1 = alpha_h(Vs(i + 1)) .* (1 - ah) - beta_h(Vs(i + 1)) .* ah;
    h(i + 1) = h(i) + (kh0 + kh1) * dt/2;
    kn0 = alpha_n(Vs(i)) .* (1 - n(i)) - beta_n(Vs(i)) .* n(i);
    an = n(i) + kn0 * dt;
    kn1 = alpha_n(Vs(i + 1)) .* (1 - an) - beta_n(Vs(i + 1)) .* an;
    n(i + 1) = n(i) + (kn0 + kn1) * dt/2;
    kw0 = (w_infinity(Vs(i)) - w(i)) ./ tau_w;
    aw = w(i) + kw0 * dt;
    kw1 = (w_infinity(Vs(i + 1)) - aw) ./ tau_w;
    w(i + 1) = w(i) + (kw0 + kw1) * dt/2;
    kl0 = (l_infinity(Vs(i)) - l(i)) ./ tau_l;
    al = l(i) + kl0 * dt;
    kl1 = (l_infinity(Vs(i + 1)) - al) ./ tau_l;
    l(i + 1) = l(i) + (kl0 + kl1) * dt/2;
    kq0 = (q_infinity(Vd(i)) - q(i)) ./ tau_q(ca_conc(i));
    aq = q(i) + kq0 * dt;
    kq1 = (q_infinity(Vd(i + 1)) - aq) ./ tau_q(ca_conc(i + 1));
    q(i + 1) = q(i) + (kq0 + kq1) * dt/2;
end
