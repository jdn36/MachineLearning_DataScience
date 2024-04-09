function [dVm_dt, statesCMpo, I_BacNav] = CourtemancheFunction(statesCMpr, I_stim, dt, boolBrugada, percBacNav, g_Na_BacNav, tauhParam, decreaseFacgNa, increaseFacgto, rateOfInact )
%COURTEMANCHE CELL MODEL
%INCLUDING THREE FACTORS TO MODEL BRUGADA SYNDROME IN EPICARDIAL CELLS 
%(THEY SHOULD BE ALL 1 FOR NORMAL AP): decreaseFacgNa, increaseFacgto, rateOfInact

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
% decreaseFacgNa = 1.;
% increaseFacgto = 1.;
% rateOfInact = 1.;
% if boolBrugada
%     decreaseFacgNa = .5; %factor to decrease conductance of Na, 1 if normal conductance
%     increaseFacgto = 7.; %factor to increase conductance of transient outwards, 1 if normal conductance
%     rateOfInact = 3.5;%factor to control speed of fast inactivation for INa, 1 if normal speed
% end
% g_Na_BacNav = 10.; %nS/pF max I_Na_BacNav conductance   (21.6399 is original)  instead of 10 to affect peak of notch in AP (the higher it is, the higher the peak is)
% tauhParam = 90.;%parameter within calculation of tau_h_BN to affect AP shape, the higher it is the shorter the APD is... original is 84.6609 instead of 90.

%Calculated parameters
cell_volume = pi * (cell_diameter/20000) * (cell_diameter/20000) * (cell_length/10000);% cm^3
cell_surface_area = 2.0 * pi * (cell_diameter/20000) * (cell_diameter/20000) + 2.0 * pi * (cell_diameter/20000) * (cell_length/10000);%cm^2
M_surface_to_volume_ratio = cell_surface_area /  cell_volume; %cm^-1

%Updating gating variables y(t+1) = y_inf + ( y(t) - y_inf )*exp(-dt/tau_y)
%BacNav variables
hBN = statesCMpr(1); mBN = statesCMpr(2);

%Rest of variables
Vm = statesCMpr(3); h = statesCMpr(4); d = statesCMpr(5);
x_r = statesCMpr(6); Na_i = statesCMpr(7); K_i = statesCMpr(8); Ca_rel = statesCMpr(9); o_i = statesCMpr(10);
u_i = statesCMpr(11); Cmdn_Ca_i = statesCMpr(12); Csqn_Ca_i = statesCMpr(13); v = statesCMpr(14); m = statesCMpr(15);
j = statesCMpr(16); f = statesCMpr(17); x_s = statesCMpr(18); Ca_i = statesCMpr(19); Ca_up = statesCMpr(20);
o_a=statesCMpr(21);u_a=statesCMpr(22);f_Ca=statesCMpr(23);Trpn_Ca_i=statesCMpr(24);u=statesCMpr(25);w=statesCMpr(26);


% for i=1:length(t)-1    
%     %Equilibrium potential CHECK
    E_Ca = (R*T/(2*F))*log(Ca_o./Ca_i);
    E_Na = (R*T/(1*F))*log(Na_o./Na_i);
    E_K = (R*T/(1*F))*log(K_o./K_i);

    %BacNav current
    I_BacNav = percBacNav.*g_Na_BacNav.*(mBN.^3).*hBN.*( Vm - E_Na );

    tau_m_BN = (4.2451./(exp((Vm - -38.3561)/11.43877) + exp(-(Vm - -34.4288)/1.)) + 0.14824);
    tau_h_BN = tauhParam+(0.01 - tauhParam).*(1.0./(1.0+exp((-18.9945-Vm)./2.4304))) + 0.01 + (12.47004 - 0.01)*(1.0./(1.0+exp((-40. - Vm )./1.)));
    %84.6609 instead of 90. above to go back to original (only change is in
    %tauh in bacnav)
    m_inf_BN = 1.0./(1.0 + exp((-22.1573 - Vm)./8.1769));
    h_inf_BN = 1.-1.0./(1.0+exp((-76.7507-Vm )./10.4215));
    
    %Fast Na current
    I_Na = decreaseFacgNa.*g_Na.*(m.^3).*h.*j.*( Vm - E_Na );
    
    if Vm == -47.13
        alpha_m = 3.2;
    else
        alpha_m = .32.*(Vm + 47.13)./(1 - exp(-.1.*(Vm + 47.13)));
    end
    beta_m = .08.*exp(-Vm./11);
    
    if Vm >= -40
        alpha_h = 0;
    else
        alpha_h = .135.*exp((Vm + 80)./(-6.8));
    end
    if Vm >= -40
        beta_h = 1./(.13.*(1 + exp((Vm + 10.66)./-11.1)));
    else
        beta_h = 3.56.*exp(.079.*Vm) + 310000.*exp(.35.*Vm);
    end
    
    if Vm >= -40
        alpha_j = 0;
    else
        alpha_j = (-127140.*exp(.2444.*Vm) - .00003474.*exp(-.04391.*Vm)).*((Vm + 37.78)./(1 + exp(.311.*(Vm + 79.23))));
    end
    if Vm >= -40
        beta_j = (.3.*exp(-.0000002535.*Vm))./(1 + exp(-.1.*(Vm + 32)));
    else
        beta_j = (.1212.*exp(-.01052.*Vm))./(1 + exp(-.1378.*(Vm + 40.14)));
    end
    
    tau_m = 1./(alpha_m + beta_m);
    tau_h = 1./(alpha_h + beta_h);
    tau_j = 1./(alpha_j + beta_j);
    m_inf = alpha_m./(alpha_m + beta_m);
    h_inf = alpha_h./(alpha_h + beta_h);
    j_inf = alpha_j./(alpha_j + beta_j);
    
    
    %Time-independent K current
    I_K1 = (g_K1.*(Vm - E_K))./(1 + exp(.07.*(Vm + 80)));
    
    
    %Transient outward K current
    I_to = increaseFacgto.*g_to.*(o_a.^3).*(o_i).*(Vm - E_K);
    
    alpha_oa = .65./(exp((Vm + 10)./-8.5) + exp((Vm - 30)./-59));
    beta_oa = .65./(2.5 + exp((Vm + 82)./17));
    tau_oa = 1./(K_Q10.*(alpha_oa + beta_oa));
    oa_inf = 1./(1 + exp((Vm + 20.47)./-17.54));
    alpha_oi = 1./(18.53 + exp((Vm + 113.7)./10.95));
    beta_oi = 1./(35.56 + exp((Vm + 1.26)./-7.44));
    tau_oi = 1./(K_Q10.*(alpha_oi + beta_oi));
    oi_inf = 1./(1 + exp((Vm + 43.1)./5.3));
    
    
    %Ultrarapid rapid delayed rectifier K current
    g_Kur = .005 + (.05)./(1 + exp((Vm - 15)./-13));
    I_Kur = g_Kur.*(u_a.^3).*(u_i).*(Vm - E_K);
    
    alpha_ua = .65./(exp((Vm + 10)./-8.5) + exp((Vm - 30)./-59));
    beta_ua = .65./(2.5 + exp((Vm + 82)./17));
    tau_ua = 1./(K_Q10.*(alpha_ua + beta_ua));
    ua_inf = 1./(1 + exp((Vm + 30.3)./-9.6));
    alpha_ui = 1./(21 + exp((Vm - 185)./-28));
    beta_ui = exp((Vm - 158)./16);
    tau_ui = 1./(K_Q10.*(alpha_ui + beta_ui));
    ui_inf = 1./(1 + exp((Vm - 99.45)./27.48));
    
    
    %Rapid delayed outward rectifier K current
    I_Kr = (g_Kr.*x_r.*(Vm - E_K))./(1 + exp((Vm + 15)./22.4));
    alpha_xr = (.0003.*(Vm + 14.1))./(1 - exp((Vm + 14.1)./-5));
    beta_xr = (.000073898.*(Vm - 3.3328))./(-1 + exp((Vm - 3.3328)./5.1237));
    tau_xr = 1./(alpha_xr + beta_xr);
    xr_inf = 1./(1 + exp((Vm + 14.1)./-6.5));
    
    
    %Slow delayed outward rectifier K current
    I_Ks = g_Ks.*(x_s.^2).*(Vm - E_K);
    alpha_xs = (.00004.*(Vm - 19.9))./(1 - exp((Vm - 19.9)./-17));
    beta_xs = (.000035.*(Vm - 19.9))./(-1 + exp((Vm - 19.9)./9));
    tau_xs = 1./(2.*(alpha_xs + beta_xs));
    xs_inf = 1./(sqrt(1 + exp((Vm - 19.9)./-12.7)));
    
    
    %L-type Ca current
    I_CaL = g_CaL.*d.*f.*f_Ca.*(Vm - 65);
    tau_d = (1 - exp((Vm + 10)./-6.24))./(.035.*(Vm + 10).*(1 + exp((Vm + 10)./-6.24)));
    d_inf = 1./(1 + exp((Vm + 10)./-8));
    tau_f = 9./(.02 + .0197.*exp(-(.0337.^2).*(Vm + 10).^2));
    f_inf = 1./(1 + exp((Vm + 28)./6.9));
    tau_fCa = 2;
    fCa_inf = 1./(1 + Ca_i./.00035);
    
    
    %NaK pump current
    sigma = (1./7).*(exp(Na_o./67.3) - 1);
    f_NaK = 1./(1 + (.1245.*exp(-.1.*F.*Vm./(R.*T))) + (.0365.*sigma.*exp(-F.*Vm./(R.*T))));
    I_NaK = I_NaK_max.*f_NaK.*(1./(1 + ((K_m_Na_i./Na_i).^(1.5)))).*(K_o./(K_o + K_m_K_o));
    
    
    %NaCa exchanger current
    I_NaCa = ( I_NaCa_max.*(exp(gamma.*F.*Vm./(R.*T)).*(Na_i.^3).*(Ca_o) - exp((gamma-1).*F.*Vm./(R.*T)).*(Na_o.^3).*(Ca_i)) )./( (K_m_Na.^3 + Na_o.^3).*(K_m_Ca + Ca_o).*(1 + k_sat.*exp((gamma - 1).*F.*Vm./(R.*T))) );
    
    
    %Background currents
    I_bCa = g_bCa.*(Vm - E_Ca);
    I_bNa = g_bNa.*(Vm - E_Na);
    I_bK = 0;
    
    
    %Ca pump current
    I_pCa = (I_pCa_max.*Ca_i)./(.0005 + Ca_i);
    
    
    %Ca release current from JSR
    I_rel = k_rel.*(u.^2).*v.*w.*(Ca_rel - Ca_i);
    
    F_n = (10.^(-12)).*Vol_rel.*I_rel - ((5E-13)./F).*(.5.*I_CaL - .2.*I_NaCa);
    tau_u = 8;
    u_inf = 1./(1 + exp((F_n - (3.4175E-13))./(-(13.67E-16))));
    tau_v = 1.91 + 2.09./(1 + exp((F_n - (3.4175E-13))./(-(13.67E-16))));
    v_inf = 1 - 1./(1 + exp((F_n - (6.835E-14))./(-(13.67E-16))));
    tau_w = 6.*(1 - exp((Vm - 7.9)./-5))./((1 + .3.*exp((Vm - 7.9)./-5)).*(Vm - 7.9));
    w_inf = 1 - 1./(1 + exp((Vm - 40)./(-(17))));
    
    
    %Transfer current from NSR to JSR
    tau_tr = 180;
    I_tr = (Ca_up - Ca_rel)./(tau_tr);
    
    
    %Ca uptake current by NSR
    I_up = I_up_max./(1 + (K_up./Ca_i));
    
    
    %Ca leak current by NSR
    I_up_leak = I_up_max.*Ca_up./(Ca_up_max);
    
    
    %Ca Buffer
    Cmdn_Ca_i = Cmdn_max.*(Ca_i)./(Ca_i + Km_Cmdn); 
    Trpn_Ca_i = Trpn_max.*(Ca_i)./(Ca_i + Km_Trpn); 
    Csqn_Ca_i = Csqn_max.*(Ca_rel)./(Ca_rel + Km_Csqn); 
    
%     stimu = I_mag;
    %Differential equations (RHS)
    
    I_ion = I_BacNav + I_Na + I_K1 + I_to + I_Kur + I_Kr + I_Ks + I_CaL + I_pCa + I_NaK + I_NaCa + I_bNa + I_bCa;
    dVm_dt = -(I_ion + I_stim)./Cm;
%     dVm_dt(1) = -(I_ion(1) + stimu)./Cm;
    dNai_dt = (-3.*I_NaK - 3.*I_NaCa - I_bNa - I_Na - I_BacNav)./(F.*Vol_i);
    dKi_dt = (2.*I_NaK - I_K1 - I_to - I_Kur - I_Kr - I_Ks - I_bK)./(F.*Vol_i);
    B1 = ((2.*I_NaCa - I_pCa - I_CaL - I_bCa)./(2.*F.*Vol_i)) + ((Vol_up.*(I_up_leak - I_up) + I_rel.*Vol_rel)./Vol_i);
    B2 = 1 + ((Trpn_max.*Km_Trpn)./((Ca_i + Km_Trpn).^2)) + ((Cmdn_max.*Km_Cmdn)./((Ca_i + Km_Cmdn).^2));
    dCai_dt = B1./B2;
    dCaup_dt = I_up - I_up_leak - I_tr.*Vol_rel./Vol_up;
    dCarel_dt = (I_tr - I_rel)./(1 + (((Csqn_max.*Km_Csqn)./((Ca_rel + Km_Csqn).^2))));
    
    
    %Integrating numerically
    Vm = Vm + dt.*dVm_dt;
    Na_i = Na_i + dt.*dNai_dt;
    K_i = K_i + dt.*dKi_dt;
    Ca_i = Ca_i + dt.*dCai_dt;
    Ca_up = Ca_up + dt.*dCaup_dt;
    Ca_rel = Ca_rel + dt.*dCarel_dt;
    
    h =  h + dt.*(rateOfInact).*((alpha_h).*(1 - h) - (beta_h).*h);% h_inf + ( h(i) - h_inf )*exp(-dt/tau_h);
    d = d_inf + ( d - d_inf ).*exp(-dt./tau_d);%d(i) + dt.*((d_inf - d(i))/tau_d);%
    x_r = x_r + dt.*((alpha_xr).*(1 - x_r) - (beta_xr).*x_r);%xr_inf + ( x_r(i) - xr_inf )*exp(-dt/tau_xr);
    o_i = o_i + dt.*((alpha_oi).*(1 - o_i) - (beta_oi).*o_i);%oi_inf + ( o_i(i) - oi_inf )*exp(-dt/tau_oi);
    u_i = u_i + dt.*((alpha_ui).*(1 - u_i) - (beta_ui).*u_i);%ui_inf + ( u_i(i) - ui_inf )*exp(-dt/tau_ui);
    v = v_inf + ( v - v_inf ).*exp(-dt./tau_v);
    m = m + dt.*((alpha_m).*(1 - m) - (beta_m).*m);%m_inf + ( m(i) - m_inf )*exp(-dt/tau_m);
    j = j + dt.*((alpha_j).*(1 - j) - (beta_j).*j);%j_inf + ( j(i) - j_inf )*exp(-dt/tau_j);
    f = f_inf + ( f - f_inf ).*exp(-dt./tau_f);%f(i) + dt.*((f_inf - f(i))/tau_f);%
    x_s = xs_inf + ( x_s - xs_inf ).*exp(-dt./tau_xs);
    o_a = oa_inf + ( o_a - oa_inf ).*exp(-dt./tau_oa);
    u_a = ua_inf + ( u_a - ua_inf ).*exp(-dt./tau_ua);
    f_Ca = fCa_inf + ( f_Ca - fCa_inf ).*exp(-dt./tau_fCa);%f_Ca(i) + dt.*((fCa_inf - f_Ca(i))/tau_fCa);%
    u = u_inf + ( u - u_inf ).*exp(-dt./tau_u);
    w = w_inf + ( w - w_inf ).*exp(-dt./tau_w);

    %BacNav update variables
    mBN = mBN + dt.*((m_inf_BN - mBN)./tau_m_BN);% m_inf_BN + ( mBN(i) - m_inf_BN )*exp(-dt/tau_m_BN);
    hBN = hBN + dt.*((h_inf_BN - hBN)./tau_h_BN);%h_inf_BN + ( hBN(i) - h_inf_BN )*exp(-dt/tau_h_BN);
    

statesCMpo = [hBN, mBN, Vm, h, d, x_r, Na_i, K_i, Ca_rel, o_i, u_i, Cmdn_Ca_i, Csqn_Ca_i, v, m, j, f,...
    x_s, Ca_i, Ca_up, o_a, u_a, f_Ca, Trpn_Ca_i, u, w];
    
% end
end
