%COURTEMANCHE CELL MODEL
%INCLUDING THREE FACTORS TO MODEL BRUGADA SYNDROME IN EPICARDIAL CELLS 
%(THEY SHOULD BE ALL 1 FOR NORMAL AP): decreaseFacgNa, increaseFacgto, rateOfInact

%Simulation-dependent stuff
dt = .005; %ms
tFinal = 1200;%ms
t = 0:dt:tFinal;%ms
I_begin = [50 1000 1000 1500 2100];
I_stim = 0; %pA/pF
I_mag = -15; %pA/pF
I_dur = 2; %ms
Icounter = 1;
inAP = 0;

%normal parameters
R = 8.3143; %J/(K*mol) gas constant
T = 310; %K normal body temperature
F = 96.4867; %C/mmol Faraday's constant
Cm = 1; %pF membrane capacitance
Vol_cell = 20100; %um^3  cell volume
Vol_i = 13668; %um^3 intracellular volume
Vol_up = 1109.52; %um^3 Sarcoplasmic Reticulum (SR) uptake compartment volume
Vol_rel = 96.48; %um^3 SR release compartment volume
K_o = 5.4; %mM extracellular potassium concentration
Na_o = 140; %mM extracellular sodium concentration
Ca_o = 1.8; %mM extracellular calcium concentration
g_Na = 7.8; %nS/pF max I_Na conductance
g_K1 = .09; %nS/pF max I_K1 conductance
g_to = .1652; %nS/pF max I_to conductance
g_Kr = .0294; %nS/pF max I_Kr conductance
g_Ks = .129; %nS/pF max I_Ks conductance
g_CaL = .1238; %nS/pF max I_CaL conductance
g_bCa = .00113; %nS/pF max I_bCa conductance
g_bNa = .000674; %nS/pF max I_bNa conductance
I_NaK_max = .60; %pA/pF max I_NaK
I_NaCa_max = 1600; %pA/pF max I_NaCa
I_pCa_max = .275; %pA/pF max I_pCa calcium pump current
I_up_max = .005; %mM/ms max I_up
K_Q10 = 3; %Temperature scaling factor for I_Kur and I_to kinetics
gamma = .35; %Voltage-dependence parameter for I_NaCa
K_m_Na_i = 10; %mM  Na_intra concentration half-saturation constant for I_NaK
K_m_K_o = 1.5; %mM  K_extra concentration half-saturation constant for I_NaK
K_m_Na = 87.5; %mM   Na_extra concentration half-saturation constant for I_NaCa
K_m_Ca = 1.38; %mM   Ca_extra concentration half-saturation constant for I_NaCa
k_sat = .1; %saturation factor for I_NaCa
k_rel = 30; %ms^-1    maximal release rate for I_rel
K_up = .00092; %mM   Ca_intra concentration half-saturation constant for I_up
Ca_up_max = 15; %mM    max Ca concentration in uptake compartment
Cmdn_max = .05; %mM   total Calmodulin concentration in myoplasm
Trpn_max = .07; %mM   total Troponin concentration in myoplasm
Csqn_max = 10; %mM   total Calsequestrin concentration in SR release compartment
Km_Cmdn = .00238; %mM  Ca_intra concentration half-saturation constant for Calmodulin
Km_Trpn = .0005; %mM  Ca_intra concentration half-saturation constant for Troponin
Km_Csqn = .8; %mM  Ca_rel concentration half-saturation constant for I_up
cell_length = 100; %um
cell_diameter = 16; %um

%Brugada Syndrome stuff
decreaseFacgNa = 1.; %factor to decrease conductance of Na, 1 if normal conductance
increaseFacgto = 1.; %factor to increase conductance of transient outwards, 1 if normal conductance
rateOfInact = 1.;%factor to control speed of fast inactivation for INa, 1 if normal speed
percBacNav = 0.;%factor to control how much of BacNav current will be added, 0 if no BacNav at all
g_Na_BacNav = 21.6399; %nS/pF max I_Na_BacNav conductance   (21.6399 is original)  instead of 10 to affect peak of notch in AP (the higher it is, the higher the peak is)
tauhParam = 84.6609;%parameter within calculation of tau_h_BN to affect AP shape, the higher it is the shorter the APD is... original is 84.6609 instead of 90.
hBN = zeros(1,length(t));
mBN = zeros(1,length(t));
hBN(1) = 0.8231; %IC
mBN(1) = 0.0000094; %IC

%Calculated parameters
cell_volume = pi * (cell_diameter/20000) * (cell_diameter/20000) * (cell_length/10000);% cm^3
cell_surface_area = 2.0 * pi * (cell_diameter/20000) * (cell_diameter/20000) + 2.0 * pi * (cell_diameter/20000) * (cell_length/10000);%cm^2
M_surface_to_volume_ratio = cell_surface_area /  cell_volume; %cm^-1



%Updating gating variables y(t+1) = y_inf + ( y(t) - y_inf )*exp(-dt/tau_y)

%initialization
Vm = zeros(1,length(t)); %mV
h = zeros(1,length(t));
d = zeros(1,length(t));
x_r = zeros(1,length(t));
Na_i = zeros(1,length(t));
K_i = zeros(1,length(t));
Ca_rel = zeros(1,length(t));
o_i = zeros(1,length(t));
u_i = zeros(1,length(t));
Cmdn_Ca_i = zeros(1,length(t));
Csqn_Ca_i = zeros(1,length(t));
v = zeros(1,length(t));
m = zeros(1,length(t));
j = zeros(1,length(t));
f = zeros(1,length(t));
x_s = zeros(1,length(t));
Ca_i = zeros(1,length(t));
Ca_up = zeros(1,length(t));
o_a = zeros(1,length(t));
u_a = zeros(1,length(t));
f_Ca = zeros(1,length(t));
Trpn_Ca_i = zeros(1,length(t));
u = zeros(1,length(t));
w = zeros(1,length(t));

I_BacNav = zeros(1,length(t));
I_Na = zeros(1,length(t));
I_K1 = zeros(1,length(t));
I_to = zeros(1,length(t));
I_Kur = zeros(1,length(t));
I_Kr = zeros(1,length(t));
I_Ks = zeros(1,length(t));
I_CaL = zeros(1,length(t));
I_pCa = zeros(1,length(t));
I_NaK = zeros(1,length(t));
I_NaCa = zeros(1,length(t));
I_bNa = zeros(1,length(t));
I_bCa = zeros(1,length(t));
I_tr = zeros(1,length(t));
I_rel = zeros(1,length(t));
I_up = zeros(1,length(t));
I_up_leak = zeros(1,length(t));

%IC
Vm(1) = -81.2; %mV
h(1) = .965;
d(1) = .000137;
x_r(1) = .0000329;
Na_i(1) = 11.2;
K_i(1) = 139;
Ca_rel(1) = 1.49;
o_i(1) = .999;
u_i(1) = .999;
Cmdn_Ca_i(1) = .00205;
Csqn_Ca_i(1) = 6.51;
v(1) = 1;
m(1) = .00291;
j(1) = .978;
f(1) = .999;
x_s(1) = .0187;
Ca_i(1) = .000102;
Ca_up(1) = 1.49;
o_a(1) = .0304;
u_a(1) = .00496;
f_Ca(1) = .775;
Trpn_Ca_i(1) = .0118;
u(1) = 0.0;
w(1) = .999;

I_ion = [];
I_ion(1) = 0;

tic
for i=1:length(t)-1
    
    
%     %Equilibrium potential CHECK
    E_Ca = (R*T/(2*F))*log(Ca_o/Ca_i(i));
    E_Na = (R*T/(1*F))*log(Na_o/Na_i(i));
    E_K = (R*T/(1*F))*log(K_o/K_i(i));
    %Equilibrium potential CHECK
%     E_Ca = 134;%(R*T/(2*F))*log10(Ca_o/Ca_i(i));   mV
%     E_Na = 52;%(R*T/(1*F))*log10(Na_o/Na_i(i));    mV
%     E_K = -96;%(R*T/(1*F))*log10(K_o/K_i(i));      mV

    %BacNav current
    I_BacNav(i+1) = percBacNav.*g_Na_BacNav.*(mBN(i).^3).*hBN(i).*( Vm(i) - E_Na );

    tau_m_BN = (4.2451./(exp((Vm(i) - -38.3561)/11.43877) + exp(-(Vm(i) - -34.4288)/1.)) + 0.14824);
    tau_h_BN = tauhParam+(0.01 - 84.6609)*(1.0/(1.0+exp((-18.9945-Vm(i))/2.4304))) + 0.01 + (12.47004 - 0.01)*(1.0/(1.0+exp((-40. - Vm(i) )/1.))); %((107.8)./(exp((Vm(i) + 27.15)./.1281) + exp((Vm(i) + 25.63)./(-25.19)))) + 9.593;
    %84.6609 instead of 90. above to go back to original (only change is in
    %tauh in bacnav)
    m_inf_BN = 1.0/(1.0 + exp((-22.1573 - Vm(i))/8.1769));%1./(1 + exp((-22.5-Vm(i))./2.704));
    h_inf_BN = 1.-1.0/(1.0+exp((-76.7507-Vm(i) )/10.4215)); %1./(1 + exp((77.05+Vm(i))./10.64));
    
    %Fast Na current
    I_Na(i+1) = decreaseFacgNa.*g_Na*(m(i)^3)*h(i)*j(i)*( Vm(i) - E_Na );
    
    if Vm(i) == -47.13
        alpha_m = 3.2;
    else
        alpha_m = .32*(Vm(i) + 47.13)/(1 - exp(-.1*(Vm(i) + 47.13)));
    end
    beta_m = .08*exp(-Vm(i)/11);
    
    if Vm(i) >= -40
        alpha_h = 0;
    else
        alpha_h = .135*exp((Vm(i) + 80)/(-6.8));
    end
    if Vm(i) >= -40
        beta_h = 1/(.13*(1 + exp((Vm(i) + 10.66)/-11.1)));
    else
        beta_h = 3.56*exp(.079*Vm(i)) + 310000*exp(.35*Vm(i));
    end
    
    if Vm(i) >= -40
        alpha_j = 0;
    else
        alpha_j = (-127140*exp(.2444*Vm(i)) - .00003474*exp(-.04391*Vm(i)))*((Vm(i) + 37.78)/(1 + exp(.311*(Vm(i) + 79.23))));
    end
    if Vm(i) >= -40
        beta_j = (.3*exp(-.0000002535*Vm(i)))/(1 + exp(-.1*(Vm(i) + 32)));
    else
        beta_j = (.1212*exp(-.01052*Vm(i)))/(1 + exp(-.1378*(Vm(i) + 40.14)));
    end
    
    tau_m = 1/(alpha_m + beta_m);
    tau_h = 1/(alpha_h + beta_h);
    tau_j = 1/(alpha_j + beta_j);
    m_inf = alpha_m/(alpha_m + beta_m);
    h_inf = alpha_h/(alpha_h + beta_h);
    j_inf = alpha_j/(alpha_j + beta_j);
    
    
    %Time-independent K current
    I_K1(i+1) = (g_K1*(Vm(i) - E_K))/(1 + exp(.07*(Vm(i) + 80)));
    
    
    %Transient outward K current
    I_to(i+1) = increaseFacgto.*g_to*(o_a(i)^3)*(o_i(i))*(Vm(i) - E_K);
    
    alpha_oa = .65/(exp((Vm(i) + 10)/-8.5) + exp((Vm(i) - 30)/-59));
    beta_oa = .65/(2.5 + exp((Vm(i) + 82)/17));
    tau_oa = 1/(K_Q10*(alpha_oa + beta_oa));
    oa_inf = 1/(1 + exp((Vm(i) + 20.47)/-17.54));
    alpha_oi = 1/(18.53 + exp((Vm(i) + 113.7)/10.95));
    beta_oi = 1/(35.56 + exp((Vm(i) + 1.26)/-7.44));
    tau_oi = 1/(K_Q10*(alpha_oi + beta_oi));
    oi_inf = 1/(1 + exp((Vm(i) + 43.1)/5.3));
    
    
    %Ultrarapid rapid delayed rectifier K current
    g_Kur = .005 + (.05)/(1 + exp((Vm(i) - 15)/-13));
    I_Kur(i+1) = g_Kur*(u_a(i)^3)*(u_i(i))*(Vm(i) - E_K);
    
    alpha_ua = .65/(exp((Vm(i) + 10)/-8.5) + exp((Vm(i) - 30)/-59));
    beta_ua = .65/(2.5 + exp((Vm(i) + 82)/17));
    tau_ua = 1/(K_Q10*(alpha_ua + beta_ua));
    ua_inf = 1/(1 + exp((Vm(i) + 30.3)/-9.6));
    alpha_ui = 1/(21 + exp((Vm(i) - 185)/-28));
    beta_ui = exp((Vm(i) - 158)/16);
    tau_ui = 1/(K_Q10*(alpha_ui + beta_ui));
    ui_inf = 1/(1 + exp((Vm(i) - 99.45)/27.48));
    
    
    %Rapid delayed outward rectifier K current
    I_Kr(i+1) = (g_Kr*x_r(i)*(Vm(i) - E_K))/(1 + exp((Vm(i) + 15)/22.4));
    alpha_xr = (.0003*(Vm(i) + 14.1))/(1 - exp((Vm(i) + 14.1)/-5));
    beta_xr = (.000073898*(Vm(i) - 3.3328))/(-1 + exp((Vm(i) - 3.3328)/5.1237));
    tau_xr = 1/(alpha_xr + beta_xr);
    xr_inf = 1/(1 + exp((Vm(i) + 14.1)/-6.5));
    
    
    %Slow delayed outward rectifier K current
    I_Ks(i+1) = g_Ks*(x_s(i)^2)*(Vm(i) - E_K);
    alpha_xs = (.00004*(Vm(i) - 19.9))/(1 - exp((Vm(i) - 19.9)/-17));
    beta_xs = (.000035*(Vm(i) - 19.9))/(-1 + exp((Vm(i) - 19.9)/9));
    tau_xs = 1/(2*(alpha_xs + beta_xs));
    xs_inf = 1/(sqrt(1 + exp((Vm(i) - 19.9)/-12.7)));
    
    
    %L-type Ca current
    I_CaL(i+1) = g_CaL*d(i)*f(i)*f_Ca(i)*(Vm(i) - 65);
    tau_d = (1 - exp((Vm(i) + 10)/-6.24))/(.035*(Vm(i) + 10)*(1 + exp((Vm(i) + 10)/-6.24)));
    d_inf = 1/(1 + exp((Vm(i) + 10)/-8));
    tau_f = 9/(.02 + .0197*exp(-(.0337^2)*(Vm(i) + 10)^2));
    f_inf = 1/(1 + exp((Vm(i) + 28)/6.9));
    tau_fCa = 2;
    fCa_inf = 1/(1 + Ca_i(i)/.00035);
    
    
    %NaK pump current
    sigma = (1/7)*(exp(Na_o/67.3) - 1);
    f_NaK = 1/(1 + (.1245*exp(-.1*F*Vm(i)/(R*T))) + (.0365*sigma*exp(-F*Vm(i)/(R*T))));
    I_NaK(i+1) = I_NaK_max*f_NaK*(1/(1 + ((K_m_Na_i/Na_i(i))^(1.5))))*(K_o/(K_o + K_m_K_o));
    
    
    %NaCa exchanger current
    I_NaCa(i+1) = ( I_NaCa_max*(exp(gamma*F*Vm(i)/(R*T))*(Na_i(i)^3)*(Ca_o) - exp((gamma-1)*F*Vm(i)/(R*T))*(Na_o^3)*(Ca_i(i))) )/( (K_m_Na^3 + Na_o^3)*(K_m_Ca + Ca_o)*(1 + k_sat*exp((gamma - 1)*F*Vm(i)/(R*T))) );
    
    
    %Background currents
    I_bCa(i+1) = g_bCa*(Vm(i) - E_Ca);
    I_bNa(i+1) = g_bNa*(Vm(i) - E_Na);
    I_bK = 0;
    
    
    %Ca pump current
    I_pCa(i+1) = (I_pCa_max*Ca_i(i))/(.0005 + Ca_i(i));
    
    
    %Ca release current from JSR
    I_rel(i+1) = k_rel*(u(i)^2)*v(i)*w(i)*(Ca_rel(i) - Ca_i(i));
    
    F_n = (10^(-12))*Vol_rel*I_rel(i+1) - ((5E-13)/F)*(.5*I_CaL(i+1) - .2*I_NaCa(i+1));
    tau_u = 8;
    u_inf = 1/(1 + exp((F_n - (3.4175E-13))/(-(13.67E-16))));
    tau_v = 1.91 + 2.09/(1 + exp((F_n - (3.4175E-13))/(-(13.67E-16))));
    v_inf = 1 - 1/(1 + exp((F_n - (6.835E-14))/(-(13.67E-16))));
    tau_w = 6*(1 - exp((Vm(i) - 7.9)/-5))/((1 + .3*exp((Vm(i) - 7.9)/-5))*(Vm(i) - 7.9));
    w_inf = 1 - 1/(1 + exp((Vm(i) - 40)/(-(17))));
    
    
    %Transfer current from NSR to JSR
    tau_tr = 180;
    I_tr(i+1) = (Ca_up(i) - Ca_rel(i))/(tau_tr);
    
    
    %Ca uptake current by NSR
    I_up(i+1) = I_up_max/(1 + (K_up/Ca_i(i)));
    
    
    %Ca leak current by NSR
    I_up_leak(i+1) = I_up_max*Ca_up(i)/(Ca_up_max);
    
    
    %Ca Buffer
    Cmdn_Ca_i(i+1) = Cmdn_max*(Ca_i(i))/(Ca_i(i) + Km_Cmdn); 
    Trpn_Ca_i(i+1) = Trpn_max*(Ca_i(i))/(Ca_i(i) + Km_Trpn); 
    Csqn_Ca_i(i+1) = Csqn_max*(Ca_rel(i))/(Ca_rel(i) + Km_Csqn); 
    
    
    %Differential equations (RHS)
    if t(i) >= I_begin(Icounter) && t(i) <= (I_begin(Icounter)+I_dur)
        I_stim = I_mag;
        inAP = 1;
%         display("here")
    else
        I_stim = 0;
        if inAP == 1
            Icounter = Icounter+1;
            inAP = 0;
        end
        
    end
    
    I_ion(i+1) = I_BacNav(i+1) + I_Na(i+1) + I_K1(i+1) + I_to(i+1) + I_Kur(i+1) + I_Kr(i+1) + I_Ks(i+1) + I_CaL(i+1) + I_pCa(i+1) + I_NaK(i+1) + I_NaCa(i+1) + I_bNa(i+1) + I_bCa(i+1);
    dVm_dt = -(I_ion(i+1) + I_stim)/Cm;
    dNai_dt = (-3*I_NaK(i+1) - 3*I_NaCa(i+1) - I_bNa(i+1) - I_Na(i+1) - I_BacNav(i+1))/(F*Vol_i);
    dKi_dt = (2*I_NaK(i+1) - I_K1(i+1) - I_to(i+1) - I_Kur(i+1) - I_Kr(i+1) - I_Ks(i+1) - I_bK)/(F*Vol_i);
    B1 = ((2*I_NaCa(i+1) - I_pCa(i+1) - I_CaL(i+1) - I_bCa(i+1))/(2*F*Vol_i)) + ((Vol_up*(I_up_leak(i+1) - I_up(i+1)) + I_rel(i+1)*Vol_rel)/Vol_i);
    B2 = 1 + ((Trpn_max*Km_Trpn)/((Ca_i(i) + Km_Trpn)^2)) + ((Cmdn_max*Km_Cmdn)/((Ca_i(i) + Km_Cmdn)^2));
    dCai_dt = B1/B2;
    dCaup_dt = I_up(i+1) - I_up_leak(i+1) -I_tr(i+1)*Vol_rel/Vol_up;
    dCarel_dt = (I_tr(i+1) - I_rel(i+1))/(1 + (((Csqn_max*Km_Csqn)/((Ca_rel(i) + Km_Csqn)^2))));
    
    
    %Integrating numerically
    Vm(i+1) = Vm(i) + dt*dVm_dt;
    Na_i(i+1) = Na_i(i) + dt*dNai_dt;
    K_i(i+1) = K_i(i) + dt*dKi_dt;
    Ca_i(i+1) = Ca_i(i) + dt*dCai_dt;
    Ca_up(i+1) = Ca_up(i) + dt*dCaup_dt;
    Ca_rel(i+1) = Ca_rel(i) + dt*dCarel_dt;
    
    h(i+1) =  h(i) + dt.*(rateOfInact).*((alpha_h).*(1 - h(i)) - (beta_h).*h(i));% h_inf + ( h(i) - h_inf )*exp(-dt/tau_h);
    d(i+1) = d_inf + ( d(i) - d_inf )*exp(-dt/tau_d);%d(i) + dt.*((d_inf - d(i))/tau_d);%
    x_r(i+1) = x_r(i) + dt.*((alpha_xr).*(1 - x_r(i)) - (beta_xr).*x_r(i));%xr_inf + ( x_r(i) - xr_inf )*exp(-dt/tau_xr);
    o_i(i+1) = o_i(i) + dt.*((alpha_oi).*(1 - o_i(i)) - (beta_oi).*o_i(i));%oi_inf + ( o_i(i) - oi_inf )*exp(-dt/tau_oi);
    u_i(i+1) = u_i(i) + dt.*((alpha_ui).*(1 - u_i(i)) - (beta_ui).*u_i(i));%ui_inf + ( u_i(i) - ui_inf )*exp(-dt/tau_ui);
    v(i+1) = v_inf + ( v(i) - v_inf )*exp(-dt/tau_v);
    m(i+1) = m(i) + dt.*((alpha_m).*(1 - m(i)) - (beta_m).*m(i));%m_inf + ( m(i) - m_inf )*exp(-dt/tau_m);
    j(i+1) = j(i) + dt.*((alpha_j).*(1 - j(i)) - (beta_j).*j(i));%j_inf + ( j(i) - j_inf )*exp(-dt/tau_j);
    f(i+1) = f_inf + ( f(i) - f_inf )*exp(-dt/tau_f);%f(i) + dt.*((f_inf - f(i))/tau_f);%
    x_s(i+1) = xs_inf + ( x_s(i) - xs_inf )*exp(-dt/tau_xs);
    o_a(i+1) = oa_inf + ( o_a(i) - oa_inf )*exp(-dt/tau_oa);
    u_a(i+1) = ua_inf + ( u_a(i) - ua_inf )*exp(-dt/tau_ua);
    f_Ca(i+1) = fCa_inf + ( f_Ca(i) - fCa_inf )*exp(-dt/tau_fCa);%f_Ca(i) + dt.*((fCa_inf - f_Ca(i))/tau_fCa);%
    u(i+1) = u_inf + ( u(i) - u_inf )*exp(-dt/tau_u);
    w(i+1) = w_inf + ( w(i) - w_inf )*exp(-dt/tau_w);

    mBN(i+1) = mBN(i) + dt.*((m_inf_BN - mBN(i))/tau_m_BN);% m_inf_BN + ( mBN(i) - m_inf_BN )*exp(-dt/tau_m_BN);
    hBN(i+1) = hBN(i) + dt.*((h_inf_BN - hBN(i))/tau_h_BN);%h_inf_BN + ( hBN(i) - h_inf_BN )*exp(-dt/tau_h_BN);
    
%     if (t(i+1) == I_begin(2))
%         %Vm(i+1) = -95;
%     end
    
end
toc

%%
figure
plot(t, Vm, 'r','LineWidth',4)
xlabel('Time (ms)')
ylabel('Potential (mV)')
% title('Courtemanche AP')
set(gca,'FontSize',35)


% plot(t1, Vm1, 'r--')
% hold on
% plot(t2, Vm2, 'b--')
% xlabel('Time (ms)')
% ylabel('Potential (mV)')
% title('Courtemanche AP')
% plot(t3, Vm3, 'k--')
% legend('dt of .001', 'dt of .01', 'dt of .05')

% %Test stuff
% testVm = -100:dt*10:50;
% [alpha_m, alpha_h, alpha_j, beta_m, beta_h, beta_j] = I_Na_stuff(testVm);
% tau_m = 1./(alpha_m + beta_m);
% tau_h = 1./(alpha_h + beta_h);
% tau_j = 1./(alpha_j + beta_j);
% m_inf = alpha_m./(alpha_m + beta_m);
% h_inf = alpha_h./(alpha_h + beta_h);
% j_inf = alpha_j./(alpha_j + beta_j);
% 
% % I_Na_inf = g_Na.*(m_inf.^3).*h_inf.*j_inf.*( testVm - E_Na );
% % 
% % E_Ca = (R.*T./(2.*F)).*log10(Ca_o./Ca_i);
% % E_Na = (R.*T./(1.*F)).*log10(Na_o./Na_i);
% % E_K = (R.*T./(1.*F)).*log10(K_o./K_i);
% 
% 
% 
% figure
% subplot(2,1,1)
% plot(testVm, m_inf.^3, 'r')
% hold on
% plot(testVm, h_inf, 'k')
% legend('m_inf','h_inf')
% hold off
% subplot(2,1,2)
% plot(testVm, tau_m.*500, 'r')
% hold on
% plot(testVm, tau_h, 'k')
% plot(testVm, tau_j, 'g')
% hold off
% legend('tau_m','tau_h','tau_j')
% 
% 
% E_Ca = (R*T/(2*F))*log10(Ca_o/Ca_i(1));
% E_Na = (R*T/(1*F))*log10(Na_o/Na_i(1));
% E_K = (R*T/(1*F))*log10(K_o/K_i(1));
% 
% I_K1_inf = (g_K1.*(testVm - E_K))./(1 + exp(.07.*(testVm + 80)));
% f_NaK_inf = 1./(1 + (.1245.*exp(-.1.*F.*testVm./(R.*T))) + (.0365.*sigma.*exp(-F.*testVm./(R.*T))));
% I_NaK_inf = I_NaK_max.*f_NaK_inf.*(1./(1 + ((K_m_Na_i./Na_i(1)).^(1.5)))).*(K_o./(K_o + K_m_K_o));
% I_NaCa_inf = ( I_NaCa_max.*(exp(gamma.*F.*testVm./(R.*T)).*(Na_i(1).^3).*(Ca_o) - exp((gamma-1).*F.*testVm./(R.*T)).*(Na_o.^3).*(Ca_i(1))) )./( (K_m_Na.^3 + Na_o.^3).*(K_m_Ca + Ca_o).*(1 + k_sat.*exp((gamma - 1).*F.*testVm./(R.*T))) );
% I_bCa_inf = g_bCa.*(testVm - E_Ca);
% I_bNa_inf = g_bNa.*(testVm - E_Na);
% I_pCa_inf = (I_pCa_max.*Ca_i(1))./(.0005 + Ca_i(1));
% 
% 
% figure   
% plot(testVm, I_K1_inf, 'r')
% hold on
% plot(testVm, I_NaK_inf, 'k')
% plot(testVm, I_NaK_inf+I_K1_inf+I_NaCa_inf+I_bCa_inf+I_bNa_inf+I_pCa_inf, 'g')
% legend('I_(K1)','I_(NaK)','I_v')
% hold off

% phiE = [];
% phiE(1) = 0;
% for i=2:length(t)
%     I_m = (Vm(i)-Vm(i-1))/dt - I_ion(i);
%     phiE(i) = (1./(4.*pi.*20))*I_m/.1;
% end
% 
% plot(t,phiE)
