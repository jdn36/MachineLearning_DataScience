function [finalyimportant, finaltimportant, finalcurrentsimportant, timeImport, yImport, currentsImport] = ShannonBers_rabbit(percBacNavG, HFbool, numBeats, CycleLength, param, diffRateBool)
% global PoUIC3 PoUIC2 PoUIF PoUIM1 PoUIM2 PoUC3 PoUC2 PoUC1 PoUO PoLC3 PoLC2 PoLC1 PoLO
% ParbGal
% ParWT_insDold
% Rabbit E-C coupling model.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%NB: yfinal after 150s 1Hz%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reference:
% Thomas R. Shannon, Fei Wang, Jose Puglisi, Christopher Weber, Donald M. Bers
% "A Mathematical Treatment of Integrated Ca Dynamics Within the Ventricular Myocyte",
% Biophys J., 87:3351-3371 (2004).
%

%BacNav IC's
hBNIC = 0.8231; %IC for BacNav
mBNIC = 0.0000094; %IC for BacNav
% percBacNavG = 0.6;
% HFbool = 1;
g_Na_BacNavG = 14.1722;%6.229541985593778;%10.229541985593778; %nS/pF max I_Na_BacNav conductance   (21.6399 is original)  instead of 10 to affect peak of notch in AP (the higher it is, the higher the peak is)
% numBeats = 100;
% CycleLength = 1000;

p = 0;  % Parameter array for passing nondefault conditions
%% Initial conditions
mo=1.405627e-3;
ho= 9.867005e-1;
jo=9.915620e-1; 
do=7.175662e-6; 
fo=1.000681; 
fcaBjo=2.421991e-2;
fcaBslo=1.452605e-2;
xtoso=4.051574e-3;
ytoso=9.945511e-1; 
xtofo=4.051574e-3; 
ytofo= 9.945511e-1; 
xkro=8.641386e-3; 
xkso= 5.412034e-3;
RyRro=8.884332e-1;
RyRoo=8.156628e-7; 
RyRio=1.024274e-7; 
NaBjo=3.539892;
NaBslo=7.720854e-1; 
TnCLo=8.773191e-3; 
TnCHco=1.078283e-1; 
TnCHmo=1.524002e-2; 
CaMo=2.911916e-4; 
Myoco=1.298754e-3; 
Myomo=1.381982e-1;
SRBo=2.143165e-3; 
SLLjo=9.566355e-3; 
SLLslo=1.110363e-1; 
SLHjo=7.347888e-3; 
SLHslo=7.297378e-2; 
Csqnbo= 1.242988;
Ca_sro=5.545201e-1; 
Najo=8.80329; 
Naslo=8.80733; 
Naio=8.80853;
%Naio = 8.80853 in normal models
%Naio = 35 in DAD models
Kio=135; 
Cajo=1.737475e-4; 
Caslo= 1.031812e-4; 
Caio=8.597401e-5; 
Vmo=-8.556885e+1; 
rtoso=0.9946; 
ICajuncinto=0; 
ICaslinto=0;


% Gating variables      
%   1       2       3       4       5       6       7       8       9       10      11      12      13
%%   m       h       j       d       f       fcaBj   fcaBsl   xtos    ytos    xtof    ytof    xkr     xks   
%y10=[1.2e-3;0.99;   0.99;   0.0;    1.0;    0.0141; 0.0141;     0;      1;      0.0;    1.0;    0.0;    6e-3;];
y10=[mo; ho; jo; do; fo; fcaBjo; fcaBslo; xtoso; ytoso; xtofo; ytofo; xkro; xkso;];   
% RyR and Buffering variables
%   14      15      16      17      18      19      20      21      22      23      24
%%   RyRr    RyRo    RyRi    NaBj    NaBsl   TnCL    TnCHc   TnCHm   CaM     Myoc    Myom  
y20=[RyRro; RyRoo; RyRio; NaBjo; NaBslo; TnCLo; TnCHco; TnCHmo; CaMo; Myoco; Myomo;];           
%y20=[1;     0;      0;      1.8;   0.8;    0.012;   0.112;  0.01;   0.4e-3; 1.9e-3; 0.135;];
% More buffering variables
%   25      26      27      28      29      30
%%   SRB     SLLj   SLLsl    SLHj    SLHsl  Csqnb
y30=[SRBo; SLLjo; SLLslo; SLHjo; SLHslo; Csqnbo];
%y30=[3.3e-3; 0.012; 0.012; 0.13;  0.13;  1.5;];
%   Intracellular concentrations/ Membrane voltage
%    31      32      33      34      35      36      37     38       39    40       41      42
%%    Ca_sr   Naj     Nasl    Nai     Ki      Caj     Casl    Cai    Vm    rtos    hBN      mBN
y40=[Ca_sro; Najo; Naslo; Naio;      Kio;    Cajo;    Caslo;  Caio; Vmo; rtoso; hBNIC; mBNIC];    
%y40=[0.9;    8.8;    8.8;    8.8;    135;    0.1e-3; 0.1e-3; 0.1e-3; -88;  0.89; 0;          0;];
  
% Put everything together
%                           43            44         45
y0  = [y10;y20;y30;y40; percBacNavG; g_Na_BacNavG; HFbool; 0; 0; 0; 0 ];    %Caio; Caio*.78; Caio*.21; Caio*.01
% load('yfinal'); % load output of previous simulation saved as yfinal.mat
% y0 = yfinal;

%% Single Run Simulation
% tspan = [0;50e3];
options = odeset('RelTol',1e-5,'MaxStep',1,'Stats','off'); 
% [t,y] = ode15s(@f,tspan,y0,options,p);


h4 = waitbar(0,' Matlab is working hard, Please wait ...');

if param.GlycosidesBool
    display('Glycosides simulation start')
end

x0 = y0;
y=[];t=[];

y2 =[]; t2=[];

param.caffeineBoolCurr = false;
param.DADboolLocal =  false;



for p=1:numBeats
%     [Ti,X]=ode15s(@f,[0  CycleLength],x0,options,p);
    %[t,X]=ode15s('cell_alg',[0  data.bcl],x0,opts,data);


    if diffRateBool
        if p <= 15
            param.ratePace = .00025;
            CycleLength = 1/param.ratePace;
        else
            param.ratePace = .00381;%.0073  .007
            CycleLength = 1/param.ratePace;
        end
    end

    if param.DADbool && p > numBeats - 10 %the last ten beats are the ones that can have DADs
        % x0(32) = 35;
        % x0(33) = 35;
        % x0(34) = 35;
        % CycleLength = 2000.;
        % param.ratePace = (1.e-3/2);
        
        param.DADboolLocal =  true;
    end

    [Ti,X] = ode15s(@(t, y) f(t, y, p, 'ydot', param), [0  CycleLength], x0, options);

    waitbar(p/numBeats,h4)
    t=[t; Ti+(CycleLength*(p-1))];
    y=[y; X];
    x0=[y(end,1:end)]';

    if diffRateBool
        if p >= 16 && p <= 30
            t2=[t2; Ti+(CycleLength*(p-1))];
            y2=[y2; X];
        elseif p >= numBeats - 10
            t2=[t2; Ti+(CycleLength*(p-1))];
            y2=[y2; X];
        end
    end

    if param.DADbool && p >= numBeats - 20
        t2=[t2; Ti+(CycleLength*(p-1))];
        y2=[y2; X];
    end


end

close(h4)

finalyimportant = X;
finaltimportant = Ti;
finalcurrentsimportant = calcCurrents(Ti,X,p,param);

currentsImport = [];
if diffRateBool || param.DADbool
    currentsImport = calcCurrents(t2,y2,p,param);
    yImport = y2;
    timeImport = t2;
else
    currentsImport = finalcurrentsimportant;
    yImport = finalyimportant;
    timeImport = finaltimportant;
end

size(finalcurrentsimportant)


if param.caffeineBool
    x0=[X(end,1:end)]';
    x0(14) = 1.;
    x0(15) = 1.;
    x0(16) = 1.;

    param.caffeineBoolCurr = true;
    param.ratePace = .00002;

    [Ti2,X2] = ode15s(@(t, y) f(t, y, p, 'ydot', param), [CycleLength (CycleLength + 12*CycleLength)], x0, options);

    finalyimportant = [finalyimportant; X2];
    finaltimportant = [finaltimportant; Ti2];
    finalcurrentsimportant = [finalcurrentsimportant; calcCurrents(Ti2,X2,p,param)];
end





yfinal = y(end,:);
output = yfinal;

if param.GlycosidesBool
    display('Glycosides simulation end')
end


% save 'yfinal'
% voltage=y(:,39);
% save 'voltage'
% save 't'

% currents = calcCurrents(t,y,p, param);
% figure(1);
% subplot(4,1,1); plot(t,y(:,39)); title('Voltage');
% subplot(4,1,2); plot(t,y(:,38)); title('Cyto calcium');
% subplot(4,1,3); plot(t,y(:,32),t,y(:,33),t,y(:,34)); title('Na concs');
% subplot(4,1,4); plot(t,y(:,36),t,y(:,37),t,y(:,31)); title('Ca concs');
% 




function output = f(t,y,p,runType, param)

ydot = zeros(size(y));

%% Model Parameters

%BacNav properties
g_Na_BacNav = y(44); %nS/pF max I_Na_BacNav conductance   (21.6399 is original)  instead of 10 to affect peak of notch in AP (the higher it is, the higher the peak is)
tauhParam = 84.6609;%parameter within calculation of tau_h_BN to affect AP shape, the higher it is the shorter the APD is... original is 84.6609 instead of 90.
percBacNav = y(43);

%HF
HFbool = y(45);
Itos_Multiplier = 1.;
Itof_Multiplier = 1.;
IK1_Multiplier = 1.;
INaCa_Multiplier = 1.;
INaK_Multiplier = 1.;
Jup_Multiplier = 1.;
Jupleak_Multiplier = 1.;
ICab_Multiplier = 1.;
INab_Multiplier = 1.;

IKs_Multiplier = 1.;
IKoCa_Multiplier = 1.;

if HFbool == 1
    Itos_Multiplier = .65;
    Itof_Multiplier = 0.6667;
    IK1_Multiplier = .5;
    INaCa_Multiplier = 2.;
%     INaK_Multiplier = .58;
    Jup_Multiplier = 0.6014;
    Jupleak_Multiplier = 3.0991;
%     ICab_Multiplier = 1.5294;
%     INab_Multiplier = 0.;

    IKs_Multiplier = .6;
    IKoCa_Multiplier = 3.;
end

if param.GlycosidesBool
    INaK_Multiplier = .5;
end

if param.caffeineBoolCurr
    %Open the gates for RyR activity
    y(14) = 1.;
    y(15) = 1.;
    y(16) = 1.;
end


if param.GKir ~= 1.
    IK1_Multiplier = param.GKir;
end


% %AP CLAMP BOOLEAN
% APClamp = y(62);

% Constants
R = 8314;       % [J/kmol*K]  
Frdy = 96485;   % [C/mol]  
Temp = 310;     % [K]
FoRT = Frdy/R/Temp;
Cmem = 1.3810e-10;   % [F] membrane capacitance
Qpow = (Temp-310)/10;

% Cell geometry
cellLength = 100;     % cell length [um]
cellRadius = 10.25;   % cell radius [um]
junctionLength = 160e-3;  % junc length [um]
junctionRadius = 15e-3;   % junc radius [um]
distSLcyto = 0.45;    % dist. SL to cytosol [um]
distJuncSL = 0.5;  % dist. junc to SL [um]
DcaJuncSL = 1.64e-6;  % Dca junc to SL [cm^2/sec]
DcaSLcyto = 1.22e-6; % Dca SL to cyto [cm^2/sec]
DnaJuncSL = 1.09e-5;  % Dna junc to SL [cm^2/sec]
DnaSLcyto = 1.79e-5;  % Dna SL to cyto [cm^2/sec] 
Vcell = pi*cellRadius^2*cellLength*1e-15;    % [L]
Vmyo = 0.65*Vcell; Vsr = 0.035*Vcell; Vsl = 0.02*Vcell; Vjunc = 0.0539*.01*Vcell; 
SAjunc = 20150*pi*2*junctionLength*junctionRadius;  % [um^2]
SAsl = pi*2*cellRadius*cellLength;          % [um^2]

J_ca_juncsl = 1/1.2134e12; % [L/msec] = 8.2413e-13
J_ca_slmyo = 1/2.68510e11; % [L/msec] = 3.2743e-12
J_na_juncsl = 1/(1.6382e12/3*100); % [L/msec] = 6.1043e-13
J_na_slmyo = 1/(1.8308e10/3*100);  % [L/msec] = 5.4621e-11

% Fractional currents in compartments
Fjunc = 0.11;   Fsl = 1-Fjunc;
Fjunc_CaL = 0.9; Fsl_CaL = 1-Fjunc_CaL;

% Fixed ion concentrations     
Cli = 15;   % Intracellular Cl  [mM]
Clo = 150;  % Extracellular Cl  [mM]

if param.DADboolLocal
    Ko = 3.;
else
    Ko = 5.4;   % Extracellular K   [mM]
end

Nao = 140;  % Extracellular Na  [mM]
Cao = 1.8;  % Extracellular Ca  [mM]
Mgi = 1;    % Intracellular Mg  [mM]

% Nernst Potentials
ena_junc = (1/FoRT)*log(Nao/y(32));     % [mV]
ena_sl = (1/FoRT)*log(Nao/y(33));       % [mV]
ek = (1/FoRT)*log(Ko/y(35));	        % [mV]
eca_junc = (1/FoRT/2)*log(Cao/y(36));   % [mV]
eca_sl = (1/FoRT/2)*log(Cao/y(37));     % [mV]
ecl = (1/FoRT)*log(Cli/Clo);            % [mV]

% Na transport parameters

GNa=16;
GNaB = 0.297e-3;    % [mS/uF] 
IbarNaK = 1.90719;     % [uA/uF]
KmNaip = 11;         % [mM]
KmKo = 1.5;         % [mM]
Q10NaK = 1.63;  
Q10KmNai = 1.39;

%% K current parameters
pNaK = 0.01833;      
GtoSlow = 0.06*1;     % [mS/uF] %0.09 CaMKII
GtoFast = 0.02*1;     % [mS/uF] 
gkp = 0.001;

% Cl current parameters
GClCa = 0.109625;   % [mS/uF]
GClB = 9e-3;        % [mS/uF]
KdClCa = 100e-3;    % [mM]

% I_Ca parameters
pNa = 1.5e-8;       % [cm/sec]
pCa = 5.4e-4;       % [cm/sec]
pK = 2.7e-7;        % [cm/sec]
Q10CaL = 1.8;       

% Ca transport parameters
IbarNCX = 9.0;      % [uA/uF]
KmCai = 3.59e-3;    % [mM]
KmCao = 1.3;        % [mM]
KmNai = 12.29;      % [mM]
KmNao = 87.5;       % [mM]
ksat = 0.27;        % [none]  
nu = 0.35;          % [none]
Kdact = 0.256e-3;   % [mM] 
Q10NCX = 1.57;      % [none]
IbarSLCaP = 0.0673; % [uA/uF](2.2 umol/L cytosol/sec) 
KmPCa = 0.5e-3;     % [mM] 
GCaB = 2.513e-4;    % [uA/uF] 
Q10SLCaP = 2.35;    % [none]

% SR flux parameters
Q10SRCaP = 2.6;          % [none]
Vmax_SRCaP = 5.3114e-3;  % [mM/msec] (286 umol/L cytosol/sec)
Kmf = 0.246e-3;          % [mM] default
%Kmf = 0.175e-3;          % [mM]
Kmr = 1.7;               % [mM]L cytosol
hillSRCaP = 1.787;       % [mM]
ks = 25;                 % [1/ms]      
koCa = 10;               % [mM^-2 1/ms]   %default 10   modified 20
kom = 0.06;              % [1/ms]     
kiCa = 0.5;              % [1/mM/ms]
kim = 0.005;             % [1/ms]
ec50SR = 0.45;           % [mM]

% Buffering parameters
% Note: we are using [1/ms] and [1/mM/ms], which differs from that in the paper 
% koff: [1/s] = 1e-3*[1/ms];  kon: [1/uM/s] = [1/mM/ms]
Bmax_Naj = 7.561;       % [mM]  % Na buffering
Bmax_Nasl = 1.65;       % [mM]
koff_na = 1e-3;         % [1/ms]
kon_na = 0.1e-3;        % [1/mM/ms]
Bmax_TnClow = 70e-3;    % [mM]                      % TnC low affinity
koff_tncl = 19.6e-3;    % [1/ms] 
kon_tncl = 32.7;        % [1/mM/ms]
Bmax_TnChigh = 140e-3;  % [mM]                      % TnC high affinity 
koff_tnchca = 0.032e-3; % [1/ms] 
kon_tnchca = 2.37;      % [1/mM/ms]
koff_tnchmg = 3.33e-3;  % [1/ms] 
kon_tnchmg = 3e-3;      % [1/mM/ms]
Bmax_CaM = 24e-3;       % [mM]  % CaM buffering
koff_cam = 238e-3;      % [1/ms] 
kon_cam = 34;           % [1/mM/ms]
Bmax_myosin = 140e-3;   % [mM]                      % Myosin buffering
koff_myoca = 0.46e-3;   % [1/ms]
kon_myoca = 13.8;       % [1/mM/ms]
koff_myomg = 0.057e-3;  % [1/ms]
kon_myomg = 0.0157;     % [1/mM/ms]
Bmax_SR = 19*.9e-3;     % [mM] (Bers text says 47e-3) 19e-3
koff_sr = 60e-3;        % [1/ms]
kon_sr = 100;           % [1/mM/ms]
Bmax_SLlowsl = 37.4e-3*Vmyo/Vsl;        % [mM]    % SL buffering
Bmax_SLlowj = 4.6e-3*Vmyo/Vjunc*0.1;    % [mM]    
koff_sll = 1300e-3;     % [1/ms]
kon_sll = 100;          % [1/mM/ms]
Bmax_SLhighsl = 13.4e-3*Vmyo/Vsl;       % [mM] 
Bmax_SLhighj = 1.65e-3*Vmyo/Vjunc*0.1;  % [mM] 
koff_slh = 30e-3;       % [1/ms]
kon_slh = 100;          % [1/mM/ms]
Bmax_Csqn = 140e-3*Vmyo/Vsr;            % [mM] % Bmax_Csqn = 2.6;      % Csqn buffering
koff_csqn = 65;         % [1/ms] 
kon_csqn = 100;         % [1/mM/ms] 


%% Membrane Currents

if param.DADboolLocal
    pscale_Ca = 1.6;
    pscale_K = .5;
else
    pscale_Ca = 1.; 
    pscale_K = 1.;
end

taum_scale = param.taum_scale;
tauh_scale = param.tauh_scale;

%BacNav current
%use fractions of space for E_Na to have two different bacnav currents and
%add them into one
I_BacNav_junc = Fjunc*percBacNav.*g_Na_BacNav.*(y(42).^3).*y(41).*( y(39) - ena_junc );
I_BacNav_sl = Fsl*percBacNav.*g_Na_BacNav.*(y(42).^3).*y(41).*( y(39) - ena_sl );
I_BacNav = I_BacNav_sl + I_BacNav_junc;

tau_m_BN = taum_scale*((4.2451./(exp((y(39) - -38.3561)/11.43877) + exp(-(y(39) - -34.4288)/1.)) + 0.14824));
tau_h_BN = tauh_scale*(tauhParam + (0.01 - tauhParam).*(1.0./(1.0+exp((-18.9945-y(39))./2.4304))) + 0.01 + (12.47004 - 0.01)*(1.0./(1.0+exp((-40. - y(39) )./1.))));
%84.6609 instead of 90. above to go back to original (only change is in
%tauh in bacnav)
m_inf_BN = 1.0./(1.0 + exp((-22.1573 - y(39))./8.1769));
h_inf_BN = 1.-1.0./(1.0+exp((-76.7507-y(39) )./10.4215));


ydot(41) = (h_inf_BN - y(41)) / tau_h_BN;
ydot(42) = (m_inf_BN - y(42)) / tau_m_BN;


% I_Na: Fast Na Current
am = 0.32*(y(39)+47.13)/(1-exp(-0.1*(y(39)+47.13)));
bm = 0.08*exp(-y(39)/11);
if y(39) >= -40
    ah = 0; aj = 0;
    bh = 1/(0.13*(1+exp(-(y(39)+10.66)/11.1)));
    bj = 0.3*exp(-2.535e-7*y(39))/(1+exp(-0.1*(y(39)+32)));
else
    ah = 0.135*exp((80+y(39))/-6.8);
    bh = 3.56*exp(0.079*y(39))+3.1e5*exp(0.35*y(39));
    aj = (-1.2714e5*exp(0.2444*y(39))-3.474e-5*exp(-0.04391*y(39)))*(y(39)+37.78)/(1+exp(0.311*(y(39)+79.23)));
    bj = 0.1212*exp(-0.01052*y(39))/(1+exp(-0.1378*(y(39)+40.14)));
end

ydot(1) = am*(1-y(1))-bm*y(1);
ydot(2) = ah*(1-y(2))-bh*y(2);
ydot(3) = aj*(1-y(3))-bj*y(3);

% if param.DADboolLocal && y(3) <= .02
%     % ydot(3) = 0.;
% end
% if param.DADboolLocal && y(2) <= .02 %.03 - .02 worked
%     % ydot(2) = 0.;
%     % y(2) = 1.;
% end


I_Na_junc = Fjunc*GNa*y(1)^3*y(2)*y(3)*(y(39)-ena_junc);
I_Na_sl = Fsl*GNa*y(1)^3*y(2)*y(3)*(y(39)-ena_sl);

% I_Na_junc= I_Na_junc1*(1-flag)+I_Na_junc2*flag;
% I_Na_sl= I_Na_sl1*(1-flag)+I_Na_sl2*flag;
I_Na = I_Na_junc+I_Na_sl;



% I_nabk: Na Background Current
I_nabk_junc = INab_Multiplier.*(Fjunc*GNaB*(y(39)-ena_junc));
I_nabk_sl = INab_Multiplier.*(Fsl*GNaB*(y(39)-ena_sl));
I_nabk = I_nabk_junc+I_nabk_sl;

% I_nak: Na/K Pump Current
sigma = (exp(Nao/67.3)-1)/7;
fnak = 1/(1+0.1245*exp(-0.1*y(39)*FoRT)+0.0365*sigma*exp(-y(39)*FoRT));
I_nak_junc = INaK_Multiplier.*(Fjunc*IbarNaK*fnak*Ko /(1+(KmNaip/y(32))^4) /(Ko+KmKo));
I_nak_sl = INaK_Multiplier.*(Fsl*IbarNaK*fnak*Ko /(1+(KmNaip/y(33))^4) /(Ko+KmKo));
I_nak = I_nak_junc+I_nak_sl;

% I_kr: Rapidly Activating K Current
gkr = 0.03*sqrt(Ko/5.4);
xrss = 1/(1+exp(-(y(39)+50)/7.5));
tauxr = 1/(1.38e-3*(y(39)+7)/(1-exp(-0.123*(y(39)+7)))+6.1e-4*(y(39)+10)/(exp(0.145*(y(39)+10))-1));
ydot(12) = (xrss-y(12))/tauxr;
rkr = 1/(1+exp((y(39)+33)/22.4));
I_kr = gkr*y(12)*rkr*(y(39)-ek);

% I_ks: Slowly Activating K Current
pcaks_junc = -log10(y(36))+3.0; 
pcaks_sl = -log10(y(37))+3.0;  
gks_junc = 0.07*(0.057 +0.19/(1+ exp((-7.2+pcaks_junc)/0.6)));
gks_sl = 0.07*(0.057 +0.19/(1+ exp((-7.2+pcaks_sl)/0.6))); 
eks = (1/FoRT)*log((Ko+pNaK*Nao)/(y(35)+pNaK*y(34)));	
xsss = 1/(1+exp(-(y(39)-1.5)/16.7));
tauxs = 1/(7.19e-5*(y(39)+30)/(1-exp(-0.148*(y(39)+30)))+1.31e-4*(y(39)+30)/(exp(0.0687*(y(39)+30))-1)); 
ydot(13) = (xsss-y(13))/tauxs;
I_ks_junc = IKs_Multiplier.*Fjunc*gks_junc*y(13)^2*(y(39)-eks);
I_ks_sl = IKs_Multiplier.*Fsl*gks_sl*y(13)^2*(y(39)-eks);
I_ks = I_ks_junc+I_ks_sl;

%I_kp: Plateau K current
kp_kp = 1/(1+exp(7.488-y(39)/5.98));
I_kp_junc = pscale_K*(Fjunc*gkp*kp_kp*(y(39)-ek));
I_kp_sl = pscale_K*(Fsl*gkp*kp_kp*(y(39)-ek));
I_kp = I_kp_junc+I_kp_sl;

%% I_to: Transient Outward K Current (slow and fast components)
xtoss = 1/(1+exp(-(y(39)+3.0)/15));
ytoss = 1/(1+exp((y(39)+33.5)/10));
rtoss = 1/(1+exp((y(39)+33.5)/10));
tauxtos = 9/(1+exp((y(39)+3.0)/15))+0.5;
tauytos = 3e3/(1+exp((y(39)+60.0)/10))+30;
%tauytos = 182/(1+exp((y(39)+33.5)/10))+1;
taurtos = 2.8e3/(1+exp((y(39)+60.0)/10))+220; %Fei changed here!! time-dependent gating variable
%taurtos =8085/(1+exp((y(39)+33.5)/10))+313;
ydot(8) = (xtoss-y(8))/tauxtos;
ydot(9) = (ytoss-y(9))/tauytos;
ydot(40)= (rtoss-y(40))/taurtos; %Fei changed here!! time-dependent gating variable
I_tos = pscale_K*(Itos_Multiplier.*GtoSlow*y(8)*(y(9)+0.5*y(40))*(y(39)-ek));    % [uA/uF]

tauxtof = 3.5*exp(-y(39)*y(39)/30/30)+1.5;
%tauxtof = 3.5*exp(-((y(39)+3)/30)^2)+1.5;
tauytof = 20.0/(1+exp((y(39)+33.5)/10))+20.0;
%tauytof = 20.0/(1+exp((y(39)+33.5)/10))+20.0;
ydot(10) = (xtoss-y(10))/tauxtof;
ydot(11) = (ytoss-y(11))/tauytof;
I_tof = pscale_K*(Itof_Multiplier.*GtoFast*y(10)*y(11)*(y(39)-ek));
I_to = I_tos + I_tof;

% I_ki: Time-Independent K Current
aki = 1.02/(1+exp(0.2385*(y(39)-ek-59.215)));
bki =(0.49124*exp(0.08032*(y(39)+5.476-ek)) + exp(0.06175*(y(39)-ek-594.31))) /(1 + exp(-0.5143*(y(39)-ek+4.753)));
kiss = aki/(aki+bki);
I_ki = IK1_Multiplier.*0.9*sqrt(Ko/5.4)*kiss*(y(39)-ek);

% I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current
I_ClCa_junc = Fjunc*GClCa/(1+KdClCa/y(36))*(y(39)-ecl);
I_ClCa_sl = Fsl*GClCa/(1+KdClCa/y(37))*(y(39)-ecl);
I_ClCa = I_ClCa_junc+I_ClCa_sl;
I_Clbk = GClB*(y(39)-ecl);

%% I_Ca: L-type Calcium Current
dss = 1/(1+exp(-(y(39)+14.5)/6.0));
taud = dss*(1-exp(-(y(39)+14.5)/6.0))/(0.035*(y(39)+14.5));
fss = 1/(1+exp((y(39)+35.06)/3.6))+0.6/(1+exp((50-y(39))/20));
tauf = 1/(0.0197*exp( -(0.0337*(y(39)+14.5))^2 )+0.02);

ydot(4) = (dss-y(4))/taud;

if param.DADboolLocal
    ydot(4) = 1.5*((dss-y(4))/taud);
end

ydot(5) = (fss-y(5))/tauf;
ydot(6) = 1.7*y(36)*(1-y(6))-11.9e-3*y(6); % fCa_junc   koff!!!!!!!!
ydot(7) = 1.7*y(37)*(1-y(7))-11.9e-3*y(7); % fCa_sl
fcaCaMSL= 0.1/(1+(0.01/y(37)));
fcaCaj= 0.1/(1+(0.01/y(36)));
fcaCaMSL=0;
fcaCaj= 0;
%y(6)=0;
%y(7)=0;
ibarca_j = pCa*4*(y(39)*Frdy*FoRT) * (0.341*y(36)*exp(2*y(39)*FoRT)-0.341*Cao) /(exp(2*y(39)*FoRT)-1);
ibarca_sl = pCa*4*(y(39)*Frdy*FoRT) * (0.341*y(37)*exp(2*y(39)*FoRT)-0.341*Cao) /(exp(2*y(39)*FoRT)-1);
ibark = pK*(y(39)*Frdy*FoRT)*(0.75*y(35)*exp(y(39)*FoRT)-0.75*Ko) /(exp(y(39)*FoRT)-1);
ibarna_j = pNa*(y(39)*Frdy*FoRT) *(0.75*y(32)*exp(y(39)*FoRT)-0.75*Nao)  /(exp(y(39)*FoRT)-1);
ibarna_sl = pNa*(y(39)*Frdy*FoRT) *(0.75*y(33)*exp(y(39)*FoRT)-0.75*Nao)  /(exp(y(39)*FoRT)-1);
I_Ca_junc = pscale_Ca*((Fjunc_CaL*ibarca_j*y(4)*y(5)*((1-y(6))+fcaCaj)*Q10CaL^Qpow)*0.45*1);
I_Ca_sl = pscale_Ca*((Fsl_CaL*ibarca_sl*y(4)*y(5)*((1-y(7))+fcaCaMSL)*Q10CaL^Qpow)*0.45*1);
I_Ca = I_Ca_junc+I_Ca_sl;
I_CaK = (ibark*y(4)*y(5)*(Fjunc_CaL*(fcaCaj+(1-y(6)))+Fsl_CaL*(fcaCaMSL+(1-y(7))))*Q10CaL^Qpow)*0.45*1;
I_CaNa_junc = 1.*(Fjunc_CaL*ibarna_j*y(4)*y(5)*((1-y(6))+fcaCaj)*Q10CaL^Qpow)*0.45*1;
I_CaNa_sl = 1.*(Fsl_CaL*ibarna_sl*y(4)*y(5)*((1-y(7))+fcaCaMSL)*Q10CaL^Qpow)*.45*1;
I_CaNa = I_CaNa_junc+I_CaNa_sl;
I_Catot = I_Ca+I_CaK+I_CaNa;

%CaNa is the current due to the permeability of Sodium in L-type channel
%CaK is the current due to the permeability of Potassium in L-type channel

% I_ncx: Na/Ca Exchanger flux
Ka_junc = 1/(1+(Kdact/y(36))^3);
Ka_sl = 1/(1+(Kdact/y(37))^3);
s1_junc = exp(nu*y(39)*FoRT)*y(32)^3*Cao;
s1_sl = exp(nu*y(39)*FoRT)*y(33)^3*Cao;
s2_junc = exp((nu-1)*y(39)*FoRT)*Nao^3*y(36);
s3_junc = KmCai*Nao^3*(1+(y(32)/KmNai)^3) + KmNao^3*y(36)*(1+y(36)/KmCai)+KmCao*y(32)^3+y(32)^3*Cao+Nao^3*y(36);
s2_sl = exp((nu-1)*y(39)*FoRT)*Nao^3*y(37);
s3_sl = KmCai*Nao^3*(1+(y(33)/KmNai)^3) + KmNao^3*y(37)*(1+y(37)/KmCai)+KmCao*y(33)^3+y(33)^3*Cao+Nao^3*y(37);
I_ncx_junc = INaCa_Multiplier.*(Fjunc*(param.NCXMult*IbarNCX)*Q10NCX^Qpow*Ka_junc*(s1_junc-s2_junc)/s3_junc/(1+ksat*exp((nu-1)*y(39)*FoRT)));
I_ncx_sl = INaCa_Multiplier.*(Fsl*(param.NCXMult*IbarNCX)*Q10NCX^Qpow*Ka_sl*(s1_sl-s2_sl)/s3_sl/(1+ksat*exp((nu-1)*y(39)*FoRT)));
I_ncx = I_ncx_junc+I_ncx_sl;

% I_pca: Sarcolemmal Ca Pump Current
I_pca_junc = Fjunc*Q10SLCaP^Qpow*IbarSLCaP*y(36)^1.6/(KmPCa^1.6+y(36)^1.6);
I_pca_sl = Fsl*Q10SLCaP^Qpow*IbarSLCaP*y(37)^1.6/(KmPCa^1.6+y(37)^1.6);
I_pca = I_pca_junc+I_pca_sl;

% I_cabk: Ca Background Current
I_cabk_junc = ICab_Multiplier.*(Fjunc*GCaB*(y(39)-eca_junc));
I_cabk_sl = ICab_Multiplier.*(Fsl*GCaB*(y(39)-eca_sl));
I_cabk = I_cabk_junc+I_cabk_sl;

%% SR fluxes: Calcium Release, SR Ca pump, SR Ca leak
MaxSR = 15; MinSR = 1;
kCaSR = MaxSR - (MaxSR-MinSR)/(1+(ec50SR/y(31))^2.5);
koSRCa = IKoCa_Multiplier.*koCa/kCaSR;
kiSRCa = kiCa*kCaSR;
RI = 1-y(14)-y(15)-y(16);
ydot(14) = (kim*RI-kiSRCa*y(36)*y(14))-(koSRCa*y(36)^2*y(14)-kom*y(15));   % R
ydot(15) = (koSRCa*y(36)^2*y(14)-kom*y(15))-(kiSRCa*y(36)*y(15)-kim*y(16));% O
ydot(16) = (kiSRCa*y(36)*y(15)-kim*y(16))-(kom*y(16)-koSRCa*y(36)^2*RI);   % I
J_SRCarel = ks*y(15)*(y(31)-y(36));          % [mM/ms]

%SR Ca pump
J_serca = Jup_Multiplier.*(Q10SRCaP^Qpow*Vmax_SRCaP*((y(38)/Kmf)^hillSRCaP-(y(31)/Kmr)^hillSRCaP)...
    /(1+(y(38)/Kmf)^hillSRCaP+(y(31)/Kmr)^hillSRCaP));
J_SRleak = Jupleak_Multiplier.*(5.348e-6*(y(31)-y(36)));           %   [mM/ms]

%Fluxes for Ca total calculation
% J_SRcaPump = (Q*Vmax)./(1 + (Km./ Cai).^H);
J_SLCaPump_Catot = Fsl*(((2.35^1)*.002)./(1 + (.5/(y(37)*1000))^1.6)) + Fjunc*((2.35*.002)./(1 + (.5/(y(36)*1000))^1.6));%(7.3)./(1 + (.42./ (y(38)*1000)).^4.5);%Q10SLCaP^Qpow*IbarSLCaP*y(38)^1.6/(KmPCa^1.6+y(38)^1.6);%( Q10SLCaP .* 2.2 )./( 1 + ( KmPCa./ y(38) ).^1.6 );
J_SRCaPump_Catot = (Fsl*(2.6^1)*.286.*( ((y(38)*1000)./.246).^1.787 - ((y(31))./1.7).^1.787 ))./( 1 + ((y(38)*1000)./.246).^1.787 + ((y(31))./1.7).^1.787 );%Fsl_CaL*Jup_Multiplier.*(Q10SRCaP^Qpow*Vmax_SRCaP*((y(38)/Kmf)^hillSRCaP-(y(31)/Kmr)^hillSRCaP)...
    %/(1+(y(38)/Kmf)^hillSRCaP+(y(31)/Kmr)^hillSRCaP));%Fsl_CaL*(Q10SRCaP).*( ( (Vmax_SRCaP*1000*( y(38)./Kmf ).^hillSRCaP) - (286*( y(31)./Kmr ).^hillSRCaP) )./( 1 + ( y(38)./Kmf ).^hillSRCaP + ( y(31)./Kmr ).^hillSRCaP ) );

Ka_Catot = 1/(1+(Kdact/y(38))^3);
s1_Catot = exp(nu*y(39)*FoRT)*y(34)^3*Cao;
s2_Catot = exp((nu-1)*y(39)*FoRT)*Nao^3*y(38);
s3_Catot = KmCai*Nao^3*(1+(y(34)/KmNai)^3) + KmNao^3*y(38)*(1+y(38)/KmCai)+KmCao*y(34)^3+y(34)^3*Cao+Nao^3*y(38);
% J_NCX_Catot = INaCa_Multiplier.*(IbarNCX*Q10NCX^Qpow*Ka_Catot*(s1_Catot-s2_Catot)/s3_Catot/(1+ksat*exp((nu-1)*y(39)*FoRT)));%(47)./(1 + (.42./ (y(38)*1000) ).^4.1);%



Ka_junc_NCX = 1/(1+(.256/(y(36)*1000))^3);
Ka_sl_NCX =  1/(1+(.256/(y(37)*1000))^3);
s1_junc_NCX = exp(nu*y(39)*FoRT)*(y(32)*1000)^3*Cao*1000;
s1_sl_NCX = exp(nu*y(39)*FoRT)*(y(33)*1000)^3*Cao*1000;
s2_junc_NCX = exp((nu-1)*y(39)*FoRT)*(Nao*1000)^3*y(36)*1000;
s3_junc_NCX = 3.59*(Nao*1000)^3*(1+(y(32)/12.29)^3) + (87.5*1000)^3*(y(36)*1000)*(1+(y(36)*1000)/3.59)+(1.3*1000)*(y(32)*1000)^3+(y(32)*1000)^3*(1000*Cao)+(1000*Nao)^3*1000*y(36);
s2_sl_NCX = exp((nu-1)*y(39)*FoRT)*(Nao*1000)^3*y(37)*1000;
s3_sl_NCX = 3.59*(Nao*1000)^3*(1+(y(33)/12.29)^3) + (87.5*1000)^3*(y(37)*1000)*(1+(y(37)*1000)/3.59)+(1.3*1000)*(y(33)*1000)^3+(y(33)*1000)^3*(1000*Cao)+(1000*Nao)^3*1000*y(37);
IbarNCX_junc = IbarNCX*Cmem/(Frdy*Vjunc);
IbarNCX_sl = IbarNCX*Cmem/(Frdy*Vsl);
J_NCX_Catot_junc = INaCa_Multiplier.*(Fjunc*IbarNCX_junc*(Q10NCX^1)*Ka_junc_NCX*(s1_junc_NCX-s2_junc_NCX)/s3_junc_NCX/(1+ksat*exp((nu-1)*y(39)*FoRT)));
J_NCX_Catot_sl = INaCa_Multiplier.*(Fsl*IbarNCX_sl*(Q10NCX^1)*Ka_sl_NCX*(s1_sl_NCX-s2_sl_NCX)/s3_sl_NCX/(1+ksat*exp((nu-1)*y(39)*FoRT)));
J_NCX_Catot = J_NCX_Catot_junc+J_NCX_Catot_sl;

dCa_SLCaPump = J_SLCaPump_Catot;%I_pca;%J_SLCaPump_Catot;
dCa_SRCaPump = J_SRCaPump_Catot;%J_serca;%J_SRCaPump_Catot;
dCa_NCX = -2*J_NCX_Catot;% -2*I_ncx;%-2*J_NCX_Catot; (47)./(1 + (.42./ (y(38)*1000) ).^4.1);%

dCa_Tot = dCa_SLCaPump + dCa_SRCaPump + dCa_NCX;

ydot(46) = dCa_Tot;
ydot(47) = dCa_SRCaPump;
ydot(48) = dCa_NCX;
ydot(49) = dCa_SLCaPump;


%% Sodium and Calcium Buffering
ydot(17) = kon_na*y(32)*(Bmax_Naj-y(17))-koff_na*y(17);        % NaBj      [mM/ms]
ydot(18) = kon_na*y(33)*(Bmax_Nasl-y(18))-koff_na*y(18);       % NaBsl     [mM/ms]

% Cytosolic Ca Buffers
ydot(19) = kon_tncl*y(38)*(Bmax_TnClow-y(19))-koff_tncl*y(19);            % TnCL      [mM/ms]
ydot(20) = kon_tnchca*y(38)*(Bmax_TnChigh-y(20)-y(21))-koff_tnchca*y(20); % TnCHc     [mM/ms]
ydot(21) = kon_tnchmg*Mgi*(Bmax_TnChigh-y(20)-y(21))-koff_tnchmg*y(21);   % TnCHm     [mM/ms]
ydot(22) = kon_cam*y(38)*(Bmax_CaM-y(22))-koff_cam*y(22);                 % CaM       [mM/ms]
ydot(23) = kon_myoca*y(38)*(Bmax_myosin-y(23)-y(24))-koff_myoca*y(23);    % Myosin_ca [mM/ms]
ydot(24) = kon_myomg*Mgi*(Bmax_myosin-y(23)-y(24))-koff_myomg*y(24);      % Myosin_mg [mM/ms]
ydot(25) = kon_sr*y(38)*(Bmax_SR-y(25))-koff_sr*y(25);                    % SRB       [mM/ms]
J_CaB_cytosol = sum(ydot(19:25));

% Junctional and SL Ca Buffers
ydot(26) = kon_sll*y(36)*(Bmax_SLlowj-y(26))-koff_sll*y(26);       % SLLj      [mM/ms]
ydot(27) = kon_sll*y(37)*(Bmax_SLlowsl-y(27))-koff_sll*y(27);      % SLLsl     [mM/ms]
ydot(28) = kon_slh*y(36)*(Bmax_SLhighj-y(28))-koff_slh*y(28);      % SLHj      [mM/ms]
ydot(29) = kon_slh*y(37)*(Bmax_SLhighsl-y(29))-koff_slh*y(29);     % SLHsl     [mM/ms]
J_CaB_junction = ydot(26)+ydot(28);
J_CaB_sl = ydot(27)+ydot(29);

%% Ion concentrations
% SR Ca Concentrations
ydot(30) = kon_csqn*y(31)*(Bmax_Csqn-y(30))-koff_csqn*y(30);       % Csqn      [mM/ms]
ydot(31) = J_serca-(J_SRleak*Vmyo/Vsr+J_SRCarel)-ydot(30);         % Ca_sr     [mM/ms] %Ratio 3 leak current

% Sodium Concentrations
I_Na_tot_junc = I_BacNav_junc + I_Na_junc+I_nabk_junc+3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc;   % [uA/uF]
I_Na_tot_sl = I_BacNav_sl + I_Na_sl+I_nabk_sl+3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl;   % [uA/uF]


ydot(32) = -I_Na_tot_junc*Cmem/(Vjunc*Frdy)+J_na_juncsl/Vjunc*(y(33)-y(32))-ydot(17);
ydot(33) = -I_Na_tot_sl*Cmem/(Vsl*Frdy)+J_na_juncsl/Vsl*(y(32)-y(33))...
   +J_na_slmyo/Vsl*(y(34)-y(33))-ydot(18);
%ydot(32) = 0;
%ydot(33) = 0;
ydot(34) = J_na_slmyo/Vmyo*(y(33)-y(34));             % [mM/msec] 
%ydot(34)=0;

% Potassium Concentration
I_K_tot = I_to+I_kr+I_ks+I_ki-2*I_nak+I_CaK+I_kp;     % [uA/uF]
% ydot(35) = 0; %-I_K_tot*Cmem/(Vmyo*Frdy);           % [mM/msec]
ydot(35) = 0;%-I_K_tot*Cmem/(Vmyo*Frdy);

% Calcium Concentrations
I_Ca_tot_junc = I_Ca_junc+I_cabk_junc+I_pca_junc-2*I_ncx_junc;                   % [uA/uF]
I_Ca_tot_sl = I_Ca_sl+I_cabk_sl+I_pca_sl-2*I_ncx_sl;            % [uA/uF]
ydot(36) = -I_Ca_tot_junc*Cmem/(Vjunc*2*Frdy)+J_ca_juncsl/Vjunc*(y(37)-y(36))...
    -J_CaB_junction+(J_SRCarel)*Vsr/Vjunc+J_SRleak*Vmyo/Vjunc;  % Ca_j
ydot(37) = -I_Ca_tot_sl*Cmem/(Vsl*2*Frdy)+J_ca_juncsl/Vsl*(y(36)-y(37))...
    + J_ca_slmyo/Vsl*(y(38)-y(37))-J_CaB_sl;   % Ca_sl
%ydot(36)=0;
%ydot(37)=0;
% ydot(38) = -J_serca*Vsr/Vmyo-J_CaB_cytosol;%+J_ca_slmyo/Vmyo*(y(37)-y(38));    % [mM/msec]
ydot(38) = -J_serca*Vsr/Vmyo-J_CaB_cytosol +J_ca_slmyo/Vmyo*(y(37)-y(38));
%ydot(38)=0;
%if (t<15000)
%    ydot(41) = 0;
%    ydot(42) = 0;
%else
%ydot(41) = -I_Na*Cmem/(Vmyo*Frdy);
%ydot(42) = -I_Ca_sl*Cmem/(Vjunc*2*Frdy)*Vjunc/Vmyo;
%end


%% Simulation type
protocol = param.protocol;
rate = param.ratePace;

switch lower(protocol)
    case {'none',''},
        I_app = 0;
    case 'pace1',        % pace w/ current injection at rate 'rate' 9.5 nA/nF USUAL
% 		rate = 1e-3;
		if mod(t,1/rate) <= 5
            I_app = 15.;
        else
            I_app = 0.0;
        end
    
    case 'pace2',        % pace w/ current injection at rate 'rate'
		factor = 2;
        rate = factor*1e-3;
		if (mod(t+900,1/rate) <= 5) & ((t > 5000) & (t < 10000))
            I_app = 12.0;
        elseif (mod(t+900,1/rate*factor*2) <= 5)  & ((t <= 5000) | (t >= 10000))
            I_app = 12.0;
        else
            I_app = 0;
        end        
    case 'pace3',
        rate = 40e-3;
        if (t > 1000) & (t < 2000) & (mod(t+900,1/rate) <= 5) 
                    I_app = 10;    
        elseif (t > 2000) & (mod(t+900,1/rate*10) <= 5)
                    I_app = 10;
                    else
            I_app = 0;
        end
    case 'vclamp',      
		V_hold = param.vcParameters.V_hold;           
        V_test = param.vcParameters.V_test;
		if (t > param.vcParameters.tstartclamp & t < param.vcParameters.tendclamp)
		    V_clamp = V_test;
		else
		    V_clamp = V_hold;
		end
		R_clamp = param.vcParameters.R_clamp;
		I_app = (V_clamp-y(39))/R_clamp;


    case 'vclampspecial'
        if t <= param.vcParameters.tstartclamp
            V_clamp = param.vcParameters.V_hold;
        elseif (t > param.vcParameters.tstartclamp & t < param.vcParameters.tendclamp)
            V_clamp = param.vcParameters.topVoltage + (t - param.vcParameters.tstartclamp)*param.vcParameters.voltageRate;
        elseif t >= param.vcParameters.tendclamp
            V_clamp = param.vcParameters.V_hold;
        end
		R_clamp = param.vcParameters.R_clamp;
		I_app = (V_clamp-y(39))/R_clamp;

    case 'fig5',
        rate = 0.5e-3;
        %if t<=60000
         %   V_clamp=-90;
        %else
        if mod(t,1/rate)<=200
            V_clamp = -20;
		else
		    V_clamp = -90;
        end
        %end
  
    R_clamp = 0.02;
		I_app = (V_clamp-y(39))/R_clamp;
            
end  




%




%% Membrane Potential
I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;          % [uA/uF]
I_Cl_tot = I_ClCa+I_Clbk;                        % [uA/uF]
I_Ca_tot = I_Ca_tot_junc+I_Ca_tot_sl;
I_tot = I_Na_tot+I_Cl_tot+I_Ca_tot+I_K_tot;
%ydot(39) = -(I_Ca_tot+I_K_tot+I_Na_tot-I_app);
ydot(39) = -(I_tot-I_app);
vmax = ydot(39);
% ----- END EC COUPLING MODEL ---------------
% adjust output depending on the function call

if strcmp(runType,'ydot')
    output = ydot;
elseif strcmp(runType,'currents')
    currents= [ I_Na I_to I_kr I_ks I_Catot I_ncx I_BacNav I_nabk 3*I_ncx 3*I_nak I_CaNa I_Ca_junc I_Ca_sl ];
    %currents = [];
    output = currents;
end

if (nargin == 3)
    output = ydot;
elseif (nargin == 4) & strcmp(runType,'ydot')
    output = ydot;
elseif (nargin == 4) & strcmp(runType,'rates')
    output = r;
elseif (nargin == 4) & strcmp(runType,'currents')
    %currents = [I_Na I_nabk I_nak I_kr I_ks I_kp I_tos I_tof I_ki I_ClCa I_Clbk I_Catot I_ncx I_pca I_cabk J_serca*Vmyo/Vsr];
%     currents = [I_Na I_tof I_tos I_kr I_ks I_ClCa I_Catot J_SRCarel*Vsr/Vmyo J_SRleak RI I_ncx]; 
%                 1    2   3     4    5       6     7        8      9       10      11     12        13      
    currents= [ I_Na I_to I_kr I_ks I_Catot I_ncx I_BacNav I_nabk 3*I_ncx 3*I_nak I_CaNa I_Ca_junc I_Ca_sl ];
    output = currents;
end

%% Calculate timecourse for currents and other intermediates
function currents = calcCurrents(t,y,p, param)
% After running a simulation, feed the time vector and state variables into
% this function to compute ionic currents, etc.
% currents: [I_Na,I_Catot];
currents=[];
for i=1:size(t)
%     if ceil(i/1000)==i/1000
%         disp(['t = ',num2str(ceil(t(i)))]);
%     end
    currents=[currents;f(t(i),y(i,:),p,'currents', param)];
end

% end calcCurrents