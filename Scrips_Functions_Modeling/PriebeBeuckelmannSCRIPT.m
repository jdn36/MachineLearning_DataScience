%% Setting parameters
% clear 

tic

param.DADbool = false;
param.GKir = 1.;
param.verbose = true; % printing numbers of beats simulated.
options = []; % parameters for ode15s - usually empty

param.bcl = 1000; % basic cycle length in ms
beats = 1000; % number of beats
beatsToShow = 1; %in case simulation runs for long and interested only in last X number of beats
ignoreFirst = beats - beatsToShow; % this many beats at the start of the simulations are ignored when extracting the structure of simulation outputs (i.e., beats - 1 keeps the last beat).

%BACNAV PARAMETERS
param.g_Na_BacNav = 6.229541985593778;%8.8659;%10.229541985593778;%21.6399 original; modified based on IV curve peak matching
param.tauhParam = 84.6609; %parameter relevant to tau_h calculation
param.percBacNav = 0.; %multiplier of BacNav current

param.INa_Multiplier = 1.;%.3;%This can be changed to affect excitability of the cell
param.IKr_Multiplier = 1.;%1.25;%1.5;
param.rateOfInact = 1.;%3.5;
param.Jup_Multiplier = 1.;%.7; it's actually .65 in the HCM model, already low
param.ICaL_Multiplier = 1.;%0.575;%1.;%.625;

param.NCXMult = 1.;


%HF Alterations
param.HFboolean = 1;
% param.Ito_Multiplier = .64;
% param.IK1_Multiplier = .75;
% param.INaCa_Multiplier = 1.65;
% param.INaK_Multiplier = .58;
% param.Jup_Multiplier = .33333;
% param.Jupleak_Multiplier = 0.6538;
% param.ICab_Multiplier = 1.5294;
% param.INab_Multiplier = 0.;

% 
%Normal case
param.HFboolean = 0;
param.Ito_Multiplier = 1.;%.64;
param.IK1_Multiplier = 1.;%.75;
param.INaCa_Multiplier = 1.;%1.65;
param.INaK_Multiplier = 1.;%.58;
param.Jup_Multiplier = 1.;%.33333;
param.Jupleak_Multiplier = 1.;%0.6538;
param.ICab_Multiplier = 1.;%1.5294;
param.INab_Multiplier = 1.;%0.;

%CURRENT PULSE
param.stimAmp = -15.;%negative value of strength of stimulus; zero if in voltage or AP clamp; -15 for 2ms
param.stimDur = 2.; %stimulus duration






tic

%[VOI, STATES, ALGEBRAIC, CONSTANTS] = PriebeBeuckelmannFun(param, t0, tend);
param.BCL = 425.;%Every how many ms does it get paced?
param.protocol = 'pace1';%'vclamp';
param.ratePace = 4.e-3;%.001 is 1Hz since 1/.001=1000 so every 1000ms
param.g_Na_BacNav = 12.6927; %PB BacNav conductance

param.HFboolean = 0; param.percBacNav = 0.;
[VOI, STATES, ALGEBRAIC, CONSTANTS] = PriebeBeuckelmannFun(param, 0, 180000);
param.HFboolean = 0; param.percBacNav = 0.;
[VOIHF, STATESHF, ALGEBRAICHF, CONSTANTSHF] = PriebeBeuckelmannFun(param, 0, 20000);
param.HFboolean = 0; param.percBacNav = 0.05;
[VOIHF_BNd05X, STATESHF_BNd05X, ALGEBRAICHF_BNd05X, CONSTANTSHF_BNd05X] = PriebeBeuckelmannFun(param, 0, 20000);
param.HFboolean = 0; param.percBacNav = 0.2;
[VOIHF_BNd2X, STATESHF_BNd2X, ALGEBRAICHF_BNd2X, CONSTANTSHF_BNd2X] = PriebeBeuckelmannFun(param, 0, 20000);
param.HFboolean = 0; param.percBacNav = 0.4;
[VOIHF_BNd4X, STATESHF_BNd4X, ALGEBRAICHF_BNd4X, CONSTANTSHF_BNd4X] = PriebeBeuckelmannFun(param, 0, 20000);
param.HFboolean = 0; param.percBacNav = 0.6;
[VOIHF_BNd6X, STATESHF_BNd6X, ALGEBRAICHF_BNd6X, CONSTANTSHF_BNd6X] = PriebeBeuckelmannFun(param, 0, 20000);
param.HFboolean = 0; param.percBacNav = 0.8;
[VOIHF_BNd8X, STATESHF_BNd8X, ALGEBRAICHF_BNd8X, CONSTANTSHF_BNd8X] = PriebeBeuckelmannFun(param, 0, 20000);
param.HFboolean = 0; param.percBacNav = 1.;
[VOIHF_BN1X, STATESHF_BN1X, ALGEBRAICHF_BN1X, CONSTANTSHF_BN1X] = PriebeBeuckelmannFun(param, 0, 20000);

toc
%% Simulation and extraction of outputs
%     CL = parameters.bcl;
%     time = cell(beats,1);
%     X = cell(beats, 1);
%     parameters.isFailed = 0; % We assume the simulation does not fail (which can change later).
%     for n=1:beats
%         if (verbose)
%             disp(['Beat = ' num2str(n)]);
%         end
%         
%         if (maxTimePerBeat < Inf) % before each beat, restart the time counter, if time is to be measured
%             timerVar = tic;
%         end
%         
%         [time{n}, X{n}]=ode15s(parameters.model,[0 CL],X0,options,1,  cellType, ICaL_Multiplier, ...
%             INa_Multiplier, Ito_Multiplier, INaL_Multiplier, IKr_Multiplier, IKs_Multiplier, IK1_Multiplier, IKb_Multiplier,INaCa_Multiplier,...
%             INaK_Multiplier, INab_Multiplier, ICab_Multiplier, IpCa_Multiplier, ICaCl_Multiplier, IClb_Multiplier, Jrel_Multiplier,Jup_Multiplier, Jupleak_Multiplier, nao,cao,...
%             ko,ICaL_fractionSS,INaCa_fractionSS, stimAmp, stimDur, vcParameters, apClamp, extraParams, g_Na_BacNav, tauhParam, percBacNav, rateOfInact);
%         
%         
%         if ~isequal(time{n}(end), parameters.bcl) % If simulation was killed prematurely in ode15sTimed, it is handled separately, setting dummy values.
%             toc(timerVar);
%             parameters.isFailed = 1; % It is noted the simulation was killed prematurely.
%             try
%                 time(1:ignoreFirst) = [];
%                 X(1:ignoreFirst) = [];
%                 
%             catch
%                 time = [];
%                 X = [];
%             end
%             
%             return;
%         end
%         
%         X0=X{n}(size(X{n},1),:);
%     end
%   
% toc

%% Plotting everything...

indi = find( VOI >= max(VOI) - param.bcl );
indiHF = find( VOIHF >= max(VOIHF) - param.bcl );
% indiHF_BNd1X = find( VOIHF_BNd1X >= max(VOIHF_BNd1X) - param.bcl );
indiHF_BNd2X = find( VOIHF_BNd2X >= max(VOIHF_BNd2X) - param.bcl );
indiHF_BNd4X = find( VOIHF_BNd4X >= max(VOIHF_BNd4X) - param.bcl );
% indiHF_BNd5X = find( VOIHF_BNd5X >= max(VOIHF_BNd5X) - param.bcl );
indiHF_BNd6X = find( VOIHF_BNd6X >= max(VOIHF_BNd6X) - param.bcl );
indiHF_BNd8X = find( VOIHF_BNd8X >= max(VOIHF_BNd8X) - param.bcl );
indiHF_BN1X = find( VOIHF_BN1X >= max(VOIHF_BN1X) - param.bcl );

indi2 = find( VOI >= 10000 & VOI <= 11000 );
indiHF2 = find( VOIHF >= 10000 & VOIHF <= 11000 );
% indiHF2_BNd1X = find( VOIHF_BNd1X >= 10000 & VOIHF_BNd1X <= 11000 );
indiHF2_BNd2X = find( VOIHF_BNd2X >= 10000 & VOIHF_BNd2X <= 11000 );
indiHF2_BNd05X = find( VOIHF_BNd05X >= 10000 & VOIHF_BNd05X <= 11000 );
indiHF2_BNd4X = find( VOIHF_BNd4X >= 10000 & VOIHF_BNd4X <= 11000 );
% indiHF2_BNd5X = find( VOIHF_BNd5X >= 10000 & VOIHF_BNd5X <= 11000 );
indiHF2_BNd6X = find( VOIHF_BNd6X >= 10000 & VOIHF_BNd6X <= 11000 );
indiHF2_BNd8X = find( VOIHF_BNd8X >= 10000 & VOIHF_BNd8X <= 11000 );
indiHF2_BN1X = find( VOIHF_BN1X >= 10000 & VOIHF_BN1X <= 11000 );

%Vm
figure;
% plot(VOIHF(indiHF2),STATESHF(indiHF2,1), 'DisplayName','HF');
hold on;
plot(VOI(indi2),STATES(indi2,1), 'DisplayName','Healthy');
% plot(VOIHF_BNd1X(indiHF2_BNd1X),STATESHF_BNd1X(indiHF2_BNd1X,1), 'DisplayName','HF BacNav .1X');
plot(VOIHF_BNd05X(indiHF2_BNd05X),STATESHF_BNd05X(indiHF2_BNd05X,1), 'DisplayName','HF BacNav .05X');
plot(VOIHF_BNd2X(indiHF2_BNd2X),STATESHF_BNd2X(indiHF2_BNd2X,1), 'DisplayName','HF BacNav .2X');
plot(VOIHF_BNd4X(indiHF2_BNd4X),STATESHF_BNd4X(indiHF2_BNd4X,1), 'DisplayName','HF BacNav .4X');
% plot(VOIHF_BNd5X(indiHF2_BNd5X),STATESHF_BNd5X(indiHF2_BNd5X,1), 'DisplayName','HF BacNav .5X');
plot(VOIHF_BNd6X(indiHF2_BNd6X),STATESHF_BNd6X(indiHF2_BNd6X,1), 'DisplayName','HF BacNav .6X');
plot(VOIHF_BNd8X(indiHF2_BNd8X),STATESHF_BNd8X(indiHF2_BNd8X,1), 'DisplayName','HF BacNav .8X');
plot(VOIHF_BN1X(indiHF2_BN1X),STATESHF_BN1X(indiHF2_BN1X,1),'DisplayName','HF BacNav 1X');
title('Vm in HF with and without BacNav')
xlabel('time (ms)')
ylabel('Vm (mV)')
legend()


%Cai
figure;
plot(VOIHF(indiHF2),STATESHF(indiHF2,6), 'DisplayName','HF');
hold on;
plot(VOI(indi2),STATES(indi2,6), 'DisplayName','Healthy');
% plot(VOIHF_BNd1X(indiHF2_BNd1X),STATESHF_BNd1X(indiHF2_BNd1X,6), 'DisplayName','HF BacNav .1X');
plot(VOIHF_BNd2X(indiHF2_BNd2X),STATESHF_BNd2X(indiHF2_BNd2X,6), 'DisplayName','HF BacNav .2X');
plot(VOIHF_BNd4X(indiHF2_BNd4X),STATESHF_BNd4X(indiHF2_BNd4X,6), 'DisplayName','HF BacNav .4X');
% plot(VOIHF_BNd5X(indiHF2_BNd5X),STATESHF_BNd5X(indiHF2_BNd5X,6), 'DisplayName','HF BacNav .5X');
plot(VOIHF_BNd6X(indiHF2_BNd6X),STATESHF_BNd6X(indiHF2_BNd6X,6), 'DisplayName','HF BacNav .6X');
plot(VOIHF_BNd8X(indiHF2_BNd8X),STATESHF_BNd8X(indiHF2_BNd8X,6), 'DisplayName','HF BacNav .8X');
plot(VOIHF_BN1X(indiHF2_BN1X),STATESHF_BN1X(indiHF2_BN1X,6), 'DisplayName','HF BacNav 1X');
title('Ca_i concentration in HF with and without BacNav')
xlabel('time (ms)')
ylabel('Ca_i (mmol/mL)')
legend()


%ICaL
figure;
plot(VOIHF(indiHF2),ALGEBRAICHF(indiHF2,25), 'DisplayName','HF');
hold on;
plot(VOI(indi2),ALGEBRAIC(indi2,25), 'DisplayName','Healthy');
% plot(VOIHF_BNd1X(indiHF2_BNd1X),STATESHF_BNd1X(indiHF2_BNd1X,6), 'DisplayName','HF BacNav .1X');
plot(VOIHF_BNd2X(indiHF2_BNd2X),ALGEBRAICHF_BNd2X(indiHF2_BNd2X,25), 'DisplayName','HF BacNav .2X');
plot(VOIHF_BNd4X(indiHF2_BNd4X),ALGEBRAICHF_BNd4X(indiHF2_BNd4X,25), 'DisplayName','HF BacNav .4X');
% plot(VOIHF_BNd5X(indiHF2_BNd5X),STATESHF_BNd5X(indiHF2_BNd5X,6), 'DisplayName','HF BacNav .5X');
plot(VOIHF_BNd6X(indiHF2_BNd6X),ALGEBRAICHF_BNd6X(indiHF2_BNd6X,25), 'DisplayName','HF BacNav .6X');
plot(VOIHF_BNd8X(indiHF2_BNd8X),ALGEBRAICHF_BNd8X(indiHF2_BNd8X,25), 'DisplayName','HF BacNav .8X');
plot(VOIHF_BN1X(indiHF2_BN1X),ALGEBRAICHF_BN1X(indiHF2_BN1X,25), 'DisplayName','HF BacNav 1X');
title('ICa current in HF with and without BacNav')
xlabel('time (ms)')
ylabel('ICa (uA/uF)')
legend()


%INaCa
figure;
plot(VOIHF(indiHF2),ALGEBRAICHF(indiHF2,42), 'DisplayName','HF');
hold on;
plot(VOI(indi2),ALGEBRAIC(indi2,42), 'DisplayName','Healthy');
% plot(VOIHF_BNd1X(indiHF2_BNd1X),STATESHF_BNd1X(indiHF2_BNd1X,6), 'DisplayName','HF BacNav .1X');
plot(VOIHF_BNd2X(indiHF2_BNd2X),ALGEBRAICHF_BNd2X(indiHF2_BNd2X,42), 'DisplayName','HF BacNav .2X');
plot(VOIHF_BNd4X(indiHF2_BNd4X),ALGEBRAICHF_BNd4X(indiHF2_BNd4X,42), 'DisplayName','HF BacNav .4X');
% plot(VOIHF_BNd5X(indiHF2_BNd5X),STATESHF_BNd5X(indiHF2_BNd5X,6), 'DisplayName','HF BacNav .5X');
plot(VOIHF_BNd6X(indiHF2_BNd6X),ALGEBRAICHF_BNd6X(indiHF2_BNd6X,42), 'DisplayName','HF BacNav .6X');
plot(VOIHF_BNd8X(indiHF2_BNd8X),ALGEBRAICHF_BNd8X(indiHF2_BNd8X,42), 'DisplayName','HF BacNav .8X');
plot(VOIHF_BN1X(indiHF2_BN1X),ALGEBRAICHF_BN1X(indiHF2_BN1X,42), 'DisplayName','HF BacNav 1X');
title('INaCa current in HF with and without BacNav')
xlabel('time (ms)')
ylabel('INaCa (uA/uF)')
legend()


%CaSR
figure;
subplot(3,1,1)
plot(VOIHF(indiHF2),STATESHF(indiHF2,14), 'DisplayName','HF');
hold on;
plot(VOI(indi2),STATES(indi2,14), 'DisplayName','Healthy');
% plot(VOIHF_BNd1X(indiHF2_BNd1X),STATESHF_BNd1X(indiHF2_BNd1X,6), 'DisplayName','HF BacNav .1X');
plot(VOIHF_BNd2X(indiHF2_BNd2X),STATESHF_BNd2X(indiHF2_BNd2X,14), 'DisplayName','HF BacNav .2X');
plot(VOIHF_BNd4X(indiHF2_BNd4X),STATESHF_BNd4X(indiHF2_BNd4X,14), 'DisplayName','HF BacNav .4X');
% plot(VOIHF_BNd5X(indiHF2_BNd5X),STATESHF_BNd5X(indiHF2_BNd5X,6), 'DisplayName','HF BacNav .5X');
plot(VOIHF_BNd6X(indiHF2_BNd6X),STATESHF_BNd6X(indiHF2_BNd6X,14), 'DisplayName','HF BacNav .6X');
plot(VOIHF_BNd8X(indiHF2_BNd8X),STATESHF_BNd8X(indiHF2_BNd8X,14), 'DisplayName','HF BacNav .8X');
plot(VOIHF_BN1X(indiHF2_BN1X),STATESHF_BN1X(indiHF2_BN1X,14), 'DisplayName','HF BacNav 1X');
title('CaJSR concentration in HF with and without BacNav')
xlabel('time (ms)')
ylabel('CaJSR (mmol/mL)')
legend()
subplot(3,1,2)
plot(VOIHF(indiHF2),STATESHF(indiHF2,22), 'DisplayName','HF');
hold on;
plot(VOI(indi2),STATES(indi2,22), 'DisplayName','Healthy');
% plot(VOIHF_BNd1X(indiHF2_BNd1X),STATESHF_BNd1X(indiHF2_BNd1X,6), 'DisplayName','HF BacNav .1X');
plot(VOIHF_BNd2X(indiHF2_BNd2X),STATESHF_BNd2X(indiHF2_BNd2X,22), 'DisplayName','HF BacNav .2X');
plot(VOIHF_BNd4X(indiHF2_BNd4X),STATESHF_BNd4X(indiHF2_BNd4X,22), 'DisplayName','HF BacNav .4X');
% plot(VOIHF_BNd5X(indiHF2_BNd5X),STATESHF_BNd5X(indiHF2_BNd5X,6), 'DisplayName','HF BacNav .5X');
plot(VOIHF_BNd6X(indiHF2_BNd6X),STATESHF_BNd6X(indiHF2_BNd6X,22), 'DisplayName','HF BacNav .6X');
plot(VOIHF_BNd8X(indiHF2_BNd8X),STATESHF_BNd8X(indiHF2_BNd8X,22), 'DisplayName','HF BacNav .8X');
plot(VOIHF_BN1X(indiHF2_BN1X),STATESHF_BN1X(indiHF2_BN1X,22), 'DisplayName','HF BacNav 1X');
title('CaNSR concentration in HF with and without BacNav')
xlabel('time (ms)')
ylabel('CaNSR (mmol/mL)')
legend()
subplot(3,1,3)
plot(VOIHF(indiHF2),STATESHF(indiHF2,14) + STATESHF(indiHF2,22), 'DisplayName','HF');
hold on;
plot(VOI(indi2),STATES(indi2,14) + STATES(indi2,22), 'DisplayName','Healthy');
% plot(VOIHF_BNd1X(indiHF2_BNd1X),STATESHF_BNd1X(indiHF2_BNd1X,6), 'DisplayName','HF BacNav .1X');
plot(VOIHF_BNd2X(indiHF2_BNd2X),STATESHF_BNd2X(indiHF2_BNd2X,14) + STATESHF_BNd2X(indiHF2_BNd2X,22), 'DisplayName','HF BacNav .2X');
plot(VOIHF_BNd4X(indiHF2_BNd4X),STATESHF_BNd4X(indiHF2_BNd4X,14) + STATESHF_BNd4X(indiHF2_BNd4X,22), 'DisplayName','HF BacNav .4X');
% plot(VOIHF_BNd5X(indiHF2_BNd5X),STATESHF_BNd5X(indiHF2_BNd5X,6), 'DisplayName','HF BacNav .5X');
plot(VOIHF_BNd6X(indiHF2_BNd6X),STATESHF_BNd6X(indiHF2_BNd6X,14) + STATESHF_BNd6X(indiHF2_BNd6X,22), 'DisplayName','HF BacNav .6X');
plot(VOIHF_BNd8X(indiHF2_BNd8X),STATESHF_BNd8X(indiHF2_BNd8X,14) + STATESHF_BNd8X(indiHF2_BNd8X,22), 'DisplayName','HF BacNav .8X');
plot(VOIHF_BN1X(indiHF2_BN1X),STATESHF_BN1X(indiHF2_BN1X,14) + STATESHF_BN1X(indiHF2_BN1X,22), 'DisplayName','HF BacNav 1X');
title('Total Ca in SR concentration in HF with and without BacNav')
xlabel('time (ms)')
ylabel('Ca in SR (mmol/mL)')
legend()
% %I_Na vs I_BacNav
% figure
% plot(VOIHF_BN1X(indiHF_BN1X),ALGEBRAICHF_BN1X(indiHF_BN1X,52), '--','DisplayName','HF BacNav 1X');
% hold on
% plot(VOIHF(indiHF),ALGEBRAICHF(indiHF,22), 'DisplayName','HF');
% legend()

%% GRANDI
tic

param.protocol = 'pace1';%'vclamp';
param.ratePace = 1.e-3;%.001 is 1Hz since 1/.001=1000 so every 1000ms

[finalyimportant, finaltimportant, finalcurrentsimportant] = Grandi_human_vclamp2(0., 0, 1000, 1000, param);
[finalyimportantHF, finaltimportantHF, finalcurrentsimportantHF] = Grandi_human_vclamp2(0., 1, 1000, 1000, param);
[finalyimportantHF_BNd025X, finaltimportantHF_BNd025X, finalcurrentsimportantHF_BNd025X] = Grandi_human_vclamp2(0.025, 1, 1000, 1000, param);
[finalyimportantHF_BNd075X, finaltimportantHF_BNd075X, finalcurrentsimportantHF_BNd075X] = Grandi_human_vclamp2(0.075, 1, 1000, 1000, param);
[finalyimportantHF_BNd15X, finaltimportantHF_BNd15X, finalcurrentsimportantHF_BNd15X] = Grandi_human_vclamp2(0.15, 1, 1000, 1000, param);
[finalyimportantHF_BNd05X, finaltimportantHF_BNd05X, finalcurrentsimportantHF_BNd05X] = Grandi_human_vclamp2(0.05, 1, 1000, 1000, param);
[finalyimportantHF_BNd1X, finaltimportantHF_BNd1X, finalcurrentsimportantHF_BNd1X] = Grandi_human_vclamp2(0.1, 1, 1000, 1000, param);
[finalyimportantHF_BNd2X, finaltimportantHF_BNd2X, finalcurrentsimportantHF_BNd2X] = Grandi_human_vclamp2(0.2, 1, 1000, 1000, param);
[finalyimportantHF_BNd4X, finaltimportantHF_BNd4X, finalcurrentsimportantHF_BNd4X] = Grandi_human_vclamp2(0.4, 1, 1000, 1000, param);
[finalyimportantHF_BNd6X, finaltimportantHF_BNd6X, finalcurrentsimportantHF_BNd6X] = Grandi_human_vclamp2(0.6, 1, 1000, 1000, param);
[finalyimportantHF_BNd8X, finaltimportantHF_BNd8X, finalcurrentsimportantHF_BNd8X] = Grandi_human_vclamp2(0.8, 1, 1000, 1000, param);
[finalyimportantHF_BN1X, finaltimportantHF_BN1X, finalcurrentsimportantHF_BN1X] = Grandi_human_vclamp2(1., 1, 1000, 1000, param);

toc

%Vm
figure
plot(finaltimportantHF, finalyimportantHF(:,39),'DisplayName','HF')
hold on
plot(finaltimportant, finalyimportant(:,39),'DisplayName','Healthy')
plot(finaltimportantHF_BNd025X, finalyimportantHF_BNd025X(:,39),'DisplayName','HF BacNav .025X')
plot(finaltimportantHF_BNd05X, finalyimportantHF_BNd05X(:,39),'DisplayName','HF BacNav .05X')
plot(finaltimportantHF_BNd075X, finalyimportantHF_BNd075X(:,39),'DisplayName','HF BacNav .075X')
plot(finaltimportantHF_BNd1X, finalyimportantHF_BNd1X(:,39),'DisplayName','HF BacNav .1X')
plot(finaltimportantHF_BNd15X, finalyimportantHF_BNd15X(:,39),'DisplayName','HF BacNav .15X')
plot(finaltimportantHF_BNd2X, finalyimportantHF_BNd2X(:,39),'DisplayName','HF BacNav .2X')
% plot(finaltimportantHF_BNd4X, finalyimportantHF_BNd4X(:,39),'DisplayName','HF BacNav .4X')
% plot(finaltimportantHF_BNd6X, finalyimportantHF_BNd6X(:,39),'DisplayName','HF BacNav .6X')
% plot(finaltimportantHF_BNd8X, finalyimportantHF_BNd8X(:,39),'DisplayName','HF BacNav .8X')
% plot(finaltimportantHF_BN1X, finalyimportantHF_BN1X(:,39),'DisplayName','HF BacNav 1X')
title('Grandi Vm in HF with and without BacNav')
xlabel('time (ms)')
ylabel('Vm (mV)')
legend()


%Cai
figure
plot(finaltimportantHF, finalyimportantHF(:,38),'DisplayName','HF')
hold on
plot(finaltimportant, finalyimportant(:,38),'DisplayName','Healthy')
plot(finaltimportantHF_BNd025X, finalyimportantHF_BNd025X(:,38),'DisplayName','HF BacNav .025X')
plot(finaltimportantHF_BNd05X, finalyimportantHF_BNd05X(:,38),'DisplayName','HF BacNav .05X')
plot(finaltimportantHF_BNd075X, finalyimportantHF_BNd075X(:,38),'DisplayName','HF BacNav .075X')
plot(finaltimportantHF_BNd1X, finalyimportantHF_BNd1X(:,38),'DisplayName','HF BacNav .1X')
plot(finaltimportantHF_BNd15X, finalyimportantHF_BNd15X(:,38),'DisplayName','HF BacNav .15X')
plot(finaltimportantHF_BNd2X, finalyimportantHF_BNd2X(:,38),'DisplayName','HF BacNav .2X')
% plot(finaltimportantHF_BNd4X, finalyimportantHF_BNd4X(:,38),'DisplayName','HF BacNav .4X')
% plot(finaltimportantHF_BNd6X, finalyimportantHF_BNd6X(:,38),'DisplayName','HF BacNav .6X')
% plot(finaltimportantHF_BNd8X, finalyimportantHF_BNd8X(:,38),'DisplayName','HF BacNav .8X')
% plot(finaltimportantHF_BN1X, finalyimportantHF_BN1X(:,38),'DisplayName','HF BacNav 1X')
title('Grandi Ca_i in HF with and without BacNav')
xlabel('time (ms)')
ylabel('Ca_i (mmol/mL)')
legend()


%ICa
figure
plot(finaltimportantHF, finalcurrentsimportantHF(:,2),'DisplayName','HF')
hold on
plot(finaltimportant, finalcurrentsimportant(:,2),'DisplayName','Healthy')
plot(finaltimportantHF_BNd025X, finalcurrentsimportantHF_BNd025X(:,2),'DisplayName','HF BacNav .025X')
plot(finaltimportantHF_BNd05X, finalcurrentsimportantHF_BNd05X(:,2),'DisplayName','HF BacNav .05X')
plot(finaltimportantHF_BNd075X, finalcurrentsimportantHF_BNd075X(:,2),'DisplayName','HF BacNav .075X')
plot(finaltimportantHF_BNd1X, finalcurrentsimportantHF_BNd1X(:,2),'DisplayName','HF BacNav .1X')
plot(finaltimportantHF_BNd15X, finalcurrentsimportantHF_BNd15X(:,2),'DisplayName','HF BacNav .15X')
plot(finaltimportantHF_BNd2X, finalcurrentsimportantHF_BNd2X(:,2),'DisplayName','HF BacNav .2X')
% plot(finaltimportantHF_BNd4X, finalcurrentsimportantHF_BNd4X(:,2),'DisplayName','HF BacNav .4X')
% plot(finaltimportantHF_BNd6X, finalcurrentsimportantHF_BNd6X(:,2),'DisplayName','HF BacNav .6X')
% plot(finaltimportantHF_BNd8X, finalcurrentsimportantHF_BNd8X(:,2),'DisplayName','HF BacNav .8X')
% plot(finaltimportantHF_BN1X, finalcurrentsimportantHF_BN1X(:,2),'DisplayName','HF BacNav 1X')
title('Grandi ICa in HF with and without BacNav')
xlabel('time (ms)')
ylabel('ICa (uA/uF)')
legend()


%INaCa
figure
plot(finaltimportantHF, finalcurrentsimportantHF(:,3),'DisplayName','HF')
hold on
plot(finaltimportant, finalcurrentsimportant(:,3),'DisplayName','Healthy')
plot(finaltimportantHF_BNd025X, finalcurrentsimportantHF_BNd025X(:,3),'DisplayName','HF BacNav .025X')
plot(finaltimportantHF_BNd05X, finalcurrentsimportantHF_BNd05X(:,3),'DisplayName','HF BacNav .05X')
plot(finaltimportantHF_BNd075X, finalcurrentsimportantHF_BNd075X(:,3),'DisplayName','HF BacNav .075X')
plot(finaltimportantHF_BNd1X, finalcurrentsimportantHF_BNd1X(:,3),'DisplayName','HF BacNav .1X')
plot(finaltimportantHF_BNd15X, finalcurrentsimportantHF_BNd15X(:,3),'DisplayName','HF BacNav .15X')
plot(finaltimportantHF_BNd2X, finalcurrentsimportantHF_BNd2X(:,3),'DisplayName','HF BacNav .2X')
% plot(finaltimportantHF_BNd4X, finalcurrentsimportantHF_BNd4X(:,3),'DisplayName','HF BacNav .4X')
% plot(finaltimportantHF_BNd6X, finalcurrentsimportantHF_BNd6X(:,3),'DisplayName','HF BacNav .6X')
% plot(finaltimportantHF_BNd8X, finalcurrentsimportantHF_BNd8X(:,3),'DisplayName','HF BacNav .8X')
% plot(finaltimportantHF_BN1X, finalcurrentsimportantHF_BN1X(:,3),'DisplayName','HF BacNav 1X')
title('Grandi INaCa in HF with and without BacNav')
xlabel('time (ms)')
ylabel('INaCa (uA/uF)')
legend()



%CaSR
figure
subplot(3,1,1)
plot(finaltimportantHF, finalyimportantHF(:,31),'DisplayName','HF')
hold on
plot(finaltimportant, finalyimportant(:,31),'DisplayName','Healthy')
plot(finaltimportantHF_BNd025X, finalyimportantHF_BNd025X(:,31),'DisplayName','HF BacNav .025X')
plot(finaltimportantHF_BNd05X, finalyimportantHF_BNd05X(:,31),'DisplayName','HF BacNav .05X')
plot(finaltimportantHF_BNd075X, finalyimportantHF_BNd075X(:,31),'DisplayName','HF BacNav .075X')
plot(finaltimportantHF_BNd1X, finalyimportantHF_BNd1X(:,31),'DisplayName','HF BacNav .1X')
plot(finaltimportantHF_BNd15X, finalyimportantHF_BNd15X(:,31),'DisplayName','HF BacNav .15X')
plot(finaltimportantHF_BNd2X, finalyimportantHF_BNd2X(:,31),'DisplayName','HF BacNav .2X')
% plot(finaltimportantHF_BNd4X, finalyimportantHF_BNd4X(:,31),'DisplayName','HF BacNav .4X')
% plot(finaltimportantHF_BNd6X, finalyimportantHF_BNd6X(:,31),'DisplayName','HF BacNav .6X')
% plot(finaltimportantHF_BNd8X, finalyimportantHF_BNd8X(:,31),'DisplayName','HF BacNav .8X')
% plot(finaltimportantHF_BN1X, finalyimportantHF_BN1X(:,31),'DisplayName','HF BacNav 1X')
title('Grandi Ca in SR in HF with and without BacNav')
xlabel('time (ms)')
ylabel('SR Ca (mmol/mL)')
legend()
subplot(3,1,2)
plot(finaltimportantHF, finalyimportantHF(:,36),'DisplayName','HF')
hold on
plot(finaltimportant, finalyimportant(:,36),'DisplayName','Healthy')
plot(finaltimportantHF_BNd025X, finalyimportantHF_BNd025X(:,36),'DisplayName','HF BacNav .025X')
plot(finaltimportantHF_BNd05X, finalyimportantHF_BNd05X(:,36),'DisplayName','HF BacNav .05X')
plot(finaltimportantHF_BNd075X, finalyimportantHF_BNd075X(:,36),'DisplayName','HF BacNav .075X')
plot(finaltimportantHF_BNd1X, finalyimportantHF_BNd1X(:,36),'DisplayName','HF BacNav .1X')
plot(finaltimportantHF_BNd15X, finalyimportantHF_BNd15X(:,36),'DisplayName','HF BacNav .15X')
plot(finaltimportantHF_BNd2X, finalyimportantHF_BNd2X(:,36),'DisplayName','HF BacNav .2X')
% plot(finaltimportantHF_BNd4X, finalyimportantHF_BNd4X(:,36),'DisplayName','HF BacNav .4X')
% plot(finaltimportantHF_BNd6X, finalyimportantHF_BNd6X(:,36),'DisplayName','HF BacNav .6X')
% plot(finaltimportantHF_BNd8X, finalyimportantHF_BNd8X(:,36),'DisplayName','HF BacNav .8X')
% plot(finaltimportantHF_BN1X, finalyimportantHF_BN1X(:,36),'DisplayName','HF BacNav 1X')
title('Grandi Junctional Ca in HF with and without BacNav')
xlabel('time (ms)')
ylabel('Junctional Ca (mmol/mL)')
legend()
subplot(3,1,3)
plot(finaltimportantHF, finalyimportantHF(:,37),'DisplayName','HF')
hold on
plot(finaltimportant, finalyimportant(:,37),'DisplayName','Healthy')
plot(finaltimportantHF_BNd025X, finalyimportantHF_BNd025X(:,37),'DisplayName','HF BacNav .025X')
plot(finaltimportantHF_BNd05X, finalyimportantHF_BNd05X(:,37),'DisplayName','HF BacNav .05X')
plot(finaltimportantHF_BNd075X, finalyimportantHF_BNd075X(:,37),'DisplayName','HF BacNav .075X')
plot(finaltimportantHF_BNd1X, finalyimportantHF_BNd1X(:,37),'DisplayName','HF BacNav .1X')
plot(finaltimportantHF_BNd15X, finalyimportantHF_BNd15X(:,37),'DisplayName','HF BacNav .15X')
plot(finaltimportantHF_BNd2X, finalyimportantHF_BNd2X(:,37),'DisplayName','HF BacNav .2X')
% plot(finaltimportantHF_BNd4X, finalyimportantHF_BNd4X(:,37),'DisplayName','HF BacNav .4X')
% plot(finaltimportantHF_BNd6X, finalyimportantHF_BNd6X(:,37),'DisplayName','HF BacNav .6X')
% plot(finaltimportantHF_BNd8X, finalyimportantHF_BNd8X(:,37),'DisplayName','HF BacNav .8X')
% plot(finaltimportantHF_BN1X, finalyimportantHF_BN1X(:,37),'DisplayName','HF BacNav 1X')
title('Grandi Subsarcolemmal Ca in HF with and without BacNav')
xlabel('time (ms)')
ylabel('Subsarcolemmal Ca (mmol/mL)')
legend()



%% SHANNON-BERS
tic

param.protocol = 'pace1';%'vclamp';
param.ratePace = 4.e-3;%.001 is 1Hz since 1/.001=1000 so every 1000ms

param.GlycosidesBool = false;
numBeats = 1000;
CycleLength = 1000;
diffRateBool = true;

[SBfinaly, SBfinalt, SBfinalcurrents, timeImport, yImport, currentsImport] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);
% [SBfinalyHF, SBfinaltHF, SBfinalcurrentsHF] = ShannonBers_rabbit(0., 1, numBeats, 1000, param);
% figure
% plot(timeImport, yImport(:,39))
% title('Healthy')

param.GlycosidesBool = true;
[SBfinalyGLY, SBfinaltGLY, SBfinalcurrentsGLY, timeImportGLY, yImportGLY, currentsImportGLY] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);


param.GlycosidesBool = false;
[SBfinalyHF_BNd1X, SBfinaltHF_BNd1X, SBfinalcurrentsHF_BNd1X, timeImportHF_BNd1X, yImportHF_BNd1X, currentsImportHF_BNd1X] = ShannonBers_rabbit(0.1, 1, numBeats, CycleLength, param, diffRateBool);
[SBfinalyHF_BNd2X, SBfinaltHF_BNd2X, SBfinalcurrentsHF_BNd2X, timeImportHF_BNd2X, yImportHF_BNd2X, currentsImportHF_BNd2X] = ShannonBers_rabbit(0.2, 1, numBeats, CycleLength, param, diffRateBool);
[SBfinalyHF_BNd4X, SBfinaltHF_BNd4X, SBfinalcurrentsHF_BNd4X, timeImportHF_BNd4X, yImportHF_BNd4X, currentsImportHF_BNd4X] = ShannonBers_rabbit(0.4, 1, numBeats, CycleLength, param, diffRateBool);
[SBfinalyHF_BNd6X, SBfinaltHF_BNd6X, SBfinalcurrentsHF_BNd6X, timeImportHF_BNd6X, yImportHF_BNd6X, currentsImportHF_BNd6X] = ShannonBers_rabbit(0.6, 1, numBeats, CycleLength, param, diffRateBool);
[SBfinalyHF_BNd8X, SBfinaltHF_BNd8X, SBfinalcurrentsHF_BNd8X, timeImportHF_BNd8X, yImportHF_BNd8X, currentsImportHF_BNd8X] = ShannonBers_rabbit(0.8, 1, numBeats, CycleLength, param, diffRateBool);
[SBfinalyHF_BN1X, SBfinaltHF_BN1X, SBfinalcurrentsHF_BN1X, timeImportHF_BN1X, yImportHF_BN1X, currentsImportHF_BN1X] = ShannonBers_rabbit(1., 1, numBeats, CycleLength, param, diffRateBool);





tic

param.protocol = 'pace1';%'vclamp';
param.ratePace = 1.e-3;%.001 is 1Hz since 1/.001=1000 so every 1000ms

param.GlycosidesBool = false;
numBeats = 1000;
CycleLength = 1000;
diffRateBool = false;
%Healthy
[SBfinaly2, SBfinalt2, SBfinalcurrents2, timeImport2, yImport2, currentsImport2] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);

%HF
[SBfinaly2HF, SBfinalt2HF, SBfinalcurrents2HF, timeImport2HF, yImport2HF, currentsImport2HF] = ShannonBers_rabbit(0., 1, numBeats, CycleLength, param, diffRateBool);
% 
% %Glycosides + HF
% param.GlycosidesBool = true;
% [SBfinaly_glycosides2, SBfinalt_glycosides2, SBfinalcurrents_glycosides2, timeImport_glycosides2, yImport_glycosides2, currentsImport_glycosides2] = ShannonBers_rabbit(0., 1, numBeats, CycleLength, param, diffRateBool);
% 
% %Glycosides
% [SBfinaly_glycosides3, SBfinalt_glycosides3, SBfinalcurrents_glycosides3, timeImport_glycosides3, yImport_glycosides3, currentsImport_glycosides3] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);








%FINAL THING WORKING FOR PAPER
%BacNav .25X and .5X and 1X + HF
param.GlycosidesBool = false;
[SBfinalyHF_BNd1X2, SBfinaltHF_BNd1X2, SBfinalcurrentsHF_BNd1X2, timeImportHF_BNd1X2, yImportHF_BNd1X2, currentsImportHF_BNd1X2] = ShannonBers_rabbit(0.1, 1, numBeats, CycleLength, param, diffRateBool);
[SBfinalyHF_BNd2X2, SBfinaltHF_BNd2X2, SBfinalcurrentsHF_BNd2X2, timeImportHF_BNd2X2, yImportHF_BNd2X2, currentsImportHF_BNd2X2] = ShannonBers_rabbit(0.2, 1, numBeats, CycleLength, param, diffRateBool);
[SBfinalyHF_BNd4X2, SBfinaltHF_BNd4X2, SBfinalcurrentsHF_BNd4X2, timeImportHF_BNd4X2, yImportHF_BNd4X2, currentsImportHF_BNd4X2] = ShannonBers_rabbit(0.4, 1, numBeats, CycleLength, param, diffRateBool);
[SBfinalyHF_BNd6X2, SBfinaltHF_BNd6X2, SBfinalcurrentsHF_BNd6X2, timeImportHF_BNd6X2, yImportHF_BNd6X2, currentsImportHF_BNd6X2] = ShannonBers_rabbit(0.6, 1, numBeats, CycleLength, param, diffRateBool);
[SBfinalyHF_BNd8X2, SBfinaltHF_BNd8X2, SBfinalcurrentsHF_BNd8X2, timeImportHF_BNd8X2, yImportHF_BNd8X2, currentsImportHF_BNd8X2] = ShannonBers_rabbit(0.8, 1, numBeats, CycleLength, param, diffRateBool);
[SBfinalyHF_BN1X2, SBfinaltHF_BN1X2, SBfinalcurrentsHF_BN1X2, timeImportHF_BN1X2, yImportHF_BN1X2, currentsImportHF_BN1X2] = ShannonBers_rabbit(1., 1, numBeats, CycleLength, param, diffRateBool);


%BacNav .25X and .5X and 1X
[SBfinaly_BNd1X2, SBfinalt_BNd1X2, SBfinalcurrents_BNd1X2, timeImport_BNd1X2, yImport_BNd1X2, currentsImport_BNd1X2] = ShannonBers_rabbit(0.1, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd2X2, SBfinalt_BNd2X2, SBfinalcurrents_BNd2X2, timeImport_BNd2X2, yImport_BNd2X2, currentsImport_BNd2X2] = ShannonBers_rabbit(0.2, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd4X2, SBfinalt_BNd4X2, SBfinalcurrents_BNd4X2, timeImport_BNd4X2, yImport_BNd4X2, currentsImport_BNd4X2] = ShannonBers_rabbit(0.4, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd6X2, SBfinalt_BNd6X2, SBfinalcurrents_BNd6X2, timeImport_BNd6X2, yImport_BNd6X2, currentsImport_BNd6X2] = ShannonBers_rabbit(0.6, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd8X2, SBfinalt_BNd8X2, SBfinalcurrents_BNd8X2, timeImport_BNd8X2, yImport_BNd8X2, currentsImport_BNd8X2] = ShannonBers_rabbit(0.8, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BN1X2, SBfinalt_BN1X2, SBfinalcurrents_BN1X2, timeImport_BN1X2, yImport_BN1X2, currentsImport_BN1X2] = ShannonBers_rabbit(1., 0, numBeats, CycleLength, param, diffRateBool);









%%
tic

%CAFFEINE EFFECT
param.caffeineBool = true;


%Healthy
[SBfinaly22, SBfinalt22, SBfinalcurrents22, timeImport22, yImport22, currentsImport22] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);



% normCaTot = ((SBfinaly22(:,46) - min(SBfinaly22(:,46)))./(max(SBfinaly22(:,46)) - min(SBfinaly22(:,46))))*100;
% normCaSR = ((SBfinaly22(:,47) - min(SBfinaly22(:,47)))./(max(SBfinaly22(:,46)) - min(SBfinaly22(:,46))))*100;
% normCaNCX = ((SBfinaly22(:,48) - min(SBfinaly22(:,48)))./(max(SBfinaly22(:,46)) - min(SBfinaly22(:,46))))*100;
% normCaSL = ((SBfinaly22(:,49) - min(SBfinaly22(:,49)))./(max(SBfinaly22(:,46)) - min(SBfinaly22(:,46))))*100;
% 
% figure;
% plot(SBfinalt22, normCaTot);
% hold on;
% plot(SBfinalt22, normCaSR);
% plot(SBfinalt22, normCaNCX);
% plot(SBfinalt22, normCaSL);
% % legend( ['Total ' num2str(round(max(normCaTot))) ' Perc. '],'SERCA','INCX','Slow')
% legend(['Total ' num2str(round(max(normCaTot))) '%'],['SERCA ' num2str(round(max(normCaSR))) '%'],['INCX ' num2str(round(max(normCaNCX))) '%'],['Slow ' num2str(round(max(normCaSL))) '%'])
% axis([1000 3000 0 100])
% ylabel("Integrated Ca flux normalized")
% xlabel("time (ms)")


%HF
% [SBfinaly22HF, SBfinalt22HF, SBfinalcurrents22HF, timeImport22HF, yImport22HF, currentsImport22HF] = ShannonBers_rabbit(0., 1, numBeats, CycleLength, param, diffRateBool);
% 

% %BacNav .25X and .5X and 1X + HF
% param.GlycosidesBool = false;
% [SBfinalyHF_BNd1X22, SBfinaltHF_BNd1X22, SBfinalcurrentsHF_BNd1X22, timeImportHF_BNd1X22, yImportHF_BNd1X22, currentsImportHF_BNd1X22] = ShannonBers_rabbit(0.1, 1, numBeats, CycleLength, param, diffRateBool);
% [SBfinalyHF_BNd2X22, SBfinaltHF_BNd2X22, SBfinalcurrentsHF_BNd2X22, timeImportHF_BNd2X22, yImportHF_BNd2X22, currentsImportHF_BNd2X22] = ShannonBers_rabbit(0.2, 1, numBeats, CycleLength, param, diffRateBool);
% [SBfinalyHF_BNd4X22, SBfinaltHF_BNd4X22, SBfinalcurrentsHF_BNd4X22, timeImportHF_BNd4X22, yImportHF_BNd4X22, currentsImportHF_BNd4X22] = ShannonBers_rabbit(0.4, 1, numBeats, CycleLength, param, diffRateBool);
% [SBfinalyHF_BNd6X22, SBfinaltHF_BNd6X22, SBfinalcurrentsHF_BNd6X22, timeImportHF_BNd6X22, yImportHF_BNd6X22, currentsImportHF_BNd6X22] = ShannonBers_rabbit(0.6, 1, numBeats, CycleLength, param, diffRateBool);
% [SBfinalyHF_BNd8X22, SBfinaltHF_BNd8X22, SBfinalcurrentsHF_BNd8X22, timeImportHF_BNd8X22, yImportHF_BNd8X22, currentsImportHF_BNd8X22] = ShannonBers_rabbit(0.8, 1, numBeats, CycleLength, param, diffRateBool);
% [SBfinalyHF_BN1X22, SBfinaltHF_BN1X22, SBfinalcurrentsHF_BN1X22, timeImportHF_BN1X22, yImportHF_BN1X22, currentsImportHF_BN1X22] = ShannonBers_rabbit(1., 1, numBeats, CycleLength, param, diffRateBool);


%BacNav .25X and .5X and 1X
[SBfinaly_BNd1X22, SBfinalt_BNd1X22, SBfinalcurrents_BNd1X22, timeImport_BNd1X22, yImport_BNd1X22, currentsImport_BNd1X22] = ShannonBers_rabbit(0.1, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd2X22, SBfinalt_BNd2X22, SBfinalcurrents_BNd2X22, timeImport_BNd2X22, yImport_BNd2X22, currentsImport_BNd2X22] = ShannonBers_rabbit(0.2, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd4X22, SBfinalt_BNd4X22, SBfinalcurrents_BNd4X22, timeImport_BNd4X22, yImport_BNd4X22, currentsImport_BNd4X22] = ShannonBers_rabbit(0.4, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd6X22, SBfinalt_BNd6X22, SBfinalcurrents_BNd6X22, timeImport_BNd6X22, yImport_BNd6X22, currentsImport_BNd6X22] = ShannonBers_rabbit(0.6, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd8X22, SBfinalt_BNd8X22, SBfinalcurrents_BNd8X22, timeImport_BNd8X22, yImport_BNd8X22, currentsImport_BNd8X22] = ShannonBers_rabbit(0.8, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BN1X22, SBfinalt_BN1X22, SBfinalcurrents_BN1X22, timeImport_BN1X22, yImport_BN1X22, currentsImport_BN1X22] = ShannonBers_rabbit(1., 0, numBeats, CycleLength, param, diffRateBool);








normCaTot = ((SBfinaly22(:,46) - min(SBfinaly22(:,46)))./(max(SBfinaly22(:,46)) - min(SBfinaly22(:,46))))*100;
normCaSR = ((SBfinaly22(:,47) - min(SBfinaly22(:,47)))./(max(SBfinaly22(:,46)) - min(SBfinaly22(:,46))))*100;
normCaNCX = ((SBfinaly22(:,48) - min(SBfinaly22(:,48)))./(max(SBfinaly22(:,46)) - min(SBfinaly22(:,46))))*100;
normCaSL = ((SBfinaly22(:,49) - min(SBfinaly22(:,49)))./(max(SBfinaly22(:,46)) - min(SBfinaly22(:,46))))*100;

normCaTot_BNd1X = ((SBfinaly_BNd1X22(:,46) - min(SBfinaly_BNd1X22(:,46)))./(max(SBfinaly_BNd1X22(:,46)) - min(SBfinaly_BNd1X22(:,46))))*100;
normCaSR_BNd1X = ((SBfinaly_BNd1X22(:,47) - min(SBfinaly_BNd1X22(:,47)))./(max(SBfinaly_BNd1X22(:,46)) - min(SBfinaly_BNd1X22(:,46))))*100;
normCaNCX_BNd1X = ((SBfinaly_BNd1X22(:,48) - min(SBfinaly_BNd1X22(:,48)))./(max(SBfinaly_BNd1X22(:,46)) - min(SBfinaly_BNd1X22(:,46))))*100;
normCaSL_BNd1X = ((SBfinaly_BNd1X22(:,49) - min(SBfinaly_BNd1X22(:,49)))./(max(SBfinaly_BNd1X22(:,46)) - min(SBfinaly_BNd1X22(:,46))))*100;

normCaTot_BNd2X = ((SBfinaly_BNd2X22(:,46) - min(SBfinaly_BNd2X22(:,46)))./(max(SBfinaly_BNd2X22(:,46)) - min(SBfinaly_BNd2X22(:,46))))*100;
normCaSR_BNd2X = ((SBfinaly_BNd2X22(:,47) - min(SBfinaly_BNd2X22(:,47)))./(max(SBfinaly_BNd2X22(:,46)) - min(SBfinaly_BNd2X22(:,46))))*100;
normCaNCX_BNd2X = ((SBfinaly_BNd2X22(:,48) - min(SBfinaly_BNd2X22(:,48)))./(max(SBfinaly_BNd2X22(:,46)) - min(SBfinaly_BNd2X22(:,46))))*100;
normCaSL_BNd2X = ((SBfinaly_BNd2X22(:,49) - min(SBfinaly_BNd2X22(:,49)))./(max(SBfinaly_BNd2X22(:,46)) - min(SBfinaly_BNd2X22(:,46))))*100;

normCaTot_BNd4X = ((SBfinaly_BNd4X22(:,46) - min(SBfinaly_BNd4X22(:,46)))./(max(SBfinaly_BNd4X22(:,46)) - min(SBfinaly_BNd4X22(:,46))))*100;
normCaSR_BNd4X = ((SBfinaly_BNd4X22(:,47) - min(SBfinaly_BNd4X22(:,47)))./(max(SBfinaly_BNd4X22(:,46)) - min(SBfinaly_BNd4X22(:,46))))*100;
normCaNCX_BNd4X = ((SBfinaly_BNd4X22(:,48) - min(SBfinaly_BNd4X22(:,48)))./(max(SBfinaly_BNd4X22(:,46)) - min(SBfinaly_BNd4X22(:,46))))*100;
normCaSL_BNd4X = ((SBfinaly_BNd4X22(:,49) - min(SBfinaly_BNd4X22(:,49)))./(max(SBfinaly_BNd4X22(:,46)) - min(SBfinaly_BNd4X22(:,46))))*100;

normCaTot_BNd6X = ((SBfinaly_BNd6X22(:,46) - min(SBfinaly_BNd6X22(:,46)))./(max(SBfinaly_BNd6X22(:,46)) - min(SBfinaly_BNd6X22(:,46))))*100;
normCaSR_BNd6X = ((SBfinaly_BNd6X22(:,47) - min(SBfinaly_BNd6X22(:,47)))./(max(SBfinaly_BNd6X22(:,46)) - min(SBfinaly_BNd6X22(:,46))))*100;
normCaNCX_BNd6X = ((SBfinaly_BNd6X22(:,48) - min(SBfinaly_BNd6X22(:,48)))./(max(SBfinaly_BNd6X22(:,46)) - min(SBfinaly_BNd6X22(:,46))))*100;
normCaSL_BNd6X = ((SBfinaly_BNd6X22(:,49) - min(SBfinaly_BNd6X22(:,49)))./(max(SBfinaly_BNd6X22(:,46)) - min(SBfinaly_BNd6X22(:,46))))*100;

normCaTot_BNd8X = ((SBfinaly_BNd8X22(:,46) - min(SBfinaly_BNd8X22(:,46)))./(max(SBfinaly_BNd8X22(:,46)) - min(SBfinaly_BNd8X22(:,46))))*100;
normCaSR_BNd8X = ((SBfinaly_BNd8X22(:,47) - min(SBfinaly_BNd8X22(:,47)))./(max(SBfinaly_BNd8X22(:,46)) - min(SBfinaly_BNd8X22(:,46))))*100;
normCaNCX_BNd8X = ((SBfinaly_BNd8X22(:,48) - min(SBfinaly_BNd8X22(:,48)))./(max(SBfinaly_BNd8X22(:,46)) - min(SBfinaly_BNd8X22(:,46))))*100;
normCaSL_BNd8X = ((SBfinaly_BNd8X22(:,49) - min(SBfinaly_BNd8X22(:,49)))./(max(SBfinaly_BNd8X22(:,46)) - min(SBfinaly_BNd8X22(:,46))))*100;

normCaTot_BN1X = ((SBfinaly_BN1X22(:,46) - min(SBfinaly_BN1X22(:,46)))./(max(SBfinaly_BN1X22(:,46)) - min(SBfinaly_BN1X22(:,46))))*100;
normCaSR_BN1X = ((SBfinaly_BN1X22(:,47) - min(SBfinaly_BN1X22(:,47)))./(max(SBfinaly_BN1X22(:,46)) - min(SBfinaly_BN1X22(:,46))))*100;
normCaNCX_BN1X = ((SBfinaly_BN1X22(:,48) - min(SBfinaly_BN1X22(:,48)))./(max(SBfinaly_BN1X22(:,46)) - min(SBfinaly_BN1X22(:,46))))*100;
normCaSL_BN1X = ((SBfinaly_BN1X22(:,49) - min(SBfinaly_BN1X22(:,49)))./(max(SBfinaly_BN1X22(:,46)) - min(SBfinaly_BN1X22(:,46))))*100;




figure;
plot(SBfinalt22, normCaTot);
hold on;
plot(SBfinalt22, normCaSR);
plot(SBfinalt22, normCaNCX);
plot(SBfinalt22, normCaSL);
% legend( ['Total ' num2str(round(max(normCaTot))) ' Perc. '],'SERCA','INCX','Slow')
legend(['Total ' num2str(round(max(normCaTot))) '%'],['SERCA ' num2str(round(max(normCaSR))) '%'],['INCX ' num2str(round(max(normCaNCX))) '%'],['Slow ' num2str(round(max(normCaSL))) '%'])
% axis([1000 13000 0 100])
ylabel("Integrated Ca flux normalized")
xlabel("time (ms)")
title("Healthy Cell")


figure;
plot(SBfinalt_BNd1X22, normCaTot_BNd1X);
hold on;
plot(SBfinalt_BNd1X22, normCaSR_BNd1X);
plot(SBfinalt_BNd1X22, normCaNCX_BNd1X);
plot(SBfinalt_BNd1X22, normCaSL_BNd1X);
% legend( ['Total ' num2str(round(max(normCaTot))) ' Perc. '],'SERCA','INCX','Slow')
legend(['Total ' num2str(round(max(normCaTot_BNd1X))) '%'],['SERCA ' num2str(round(max(normCaSR_BNd1X))) '%'],['INCX ' num2str(round(max(normCaNCX_BNd1X))) '%'],['Slow ' num2str(round(max(normCaSL_BNd1X))) '%'])
% axis([1000 13000 0 100])
ylabel("Integrated Ca flux normalized")
xlabel("time (ms)")
title("Healthy Cell + .1X BacNav")

figure;
plot(SBfinalt_BNd2X22, normCaTot_BNd2X);
hold on;
plot(SBfinalt_BNd2X22, normCaSR_BNd2X);
plot(SBfinalt_BNd2X22, normCaNCX_BNd2X);
plot(SBfinalt_BNd2X22, normCaSL_BNd2X);
% legend( ['Total ' num2str(round(max(normCaTot))) ' Perc. '],'SERCA','INCX','Slow')
legend(['Total ' num2str(round(max(normCaTot_BNd2X))) '%'],['SERCA ' num2str(round(max(normCaSR_BNd2X))) '%'],['INCX ' num2str(round(max(normCaNCX_BNd2X))) '%'],['Slow ' num2str(round(max(normCaSL_BNd2X))) '%'])
% axis([1000 13000 0 100])
ylabel("Integrated Ca flux normalized")
xlabel("time (ms)")
title("Healthy Cell + .2X BacNav")

figure;
plot(SBfinalt_BNd4X22, normCaTot_BNd4X);
hold on;
plot(SBfinalt_BNd4X22, normCaSR_BNd4X);
plot(SBfinalt_BNd4X22, normCaNCX_BNd4X);
plot(SBfinalt_BNd4X22, normCaSL_BNd4X);
% legend( ['Total ' num2str(round(max(normCaTot))) ' Perc. '],'SERCA','INCX','Slow')
legend(['Total ' num2str(round(max(normCaTot_BNd4X))) '%'],['SERCA ' num2str(round(max(normCaSR_BNd4X))) '%'],['INCX ' num2str(round(max(normCaNCX_BNd4X))) '%'],['Slow ' num2str(round(max(normCaSL_BNd4X))) '%'])
% axis([1000 13000 0 100])
ylabel("Integrated Ca flux normalized")
xlabel("time (ms)")
title("Healthy Cell + .4X BacNav")

figure;
plot(SBfinalt_BNd6X22, normCaTot_BNd6X);
hold on;
plot(SBfinalt_BNd6X22, normCaSR_BNd6X);
plot(SBfinalt_BNd6X22, normCaNCX_BNd6X);
plot(SBfinalt_BNd6X22, normCaSL_BNd6X);
% legend( ['Total ' num2str(round(max(normCaTot))) ' Perc. '],'SERCA','INCX','Slow')
legend(['Total ' num2str(round(max(normCaTot_BNd6X))) '%'],['SERCA ' num2str(round(max(normCaSR_BNd6X))) '%'],['INCX ' num2str(round(max(normCaNCX_BNd6X))) '%'],['Slow ' num2str(round(max(normCaSL_BNd6X))) '%'])
% axis([1000 13000 0 100])
ylabel("Integrated Ca flux normalized")
xlabel("time (ms)")
title("Healthy Cell + .6X BacNav")

figure;
plot(SBfinalt_BNd8X22, normCaTot_BNd8X);
hold on;
plot(SBfinalt_BNd8X22, normCaSR_BNd8X);
plot(SBfinalt_BNd8X22, normCaNCX_BNd8X);
plot(SBfinalt_BNd8X22, normCaSL_BNd8X);
% legend( ['Total ' num2str(round(max(normCaTot))) ' Perc. '],'SERCA','INCX','Slow')
legend(['Total ' num2str(round(max(normCaTot_BNd8X))) '%'],['SERCA ' num2str(round(max(normCaSR_BNd8X))) '%'],['INCX ' num2str(round(max(normCaNCX_BNd8X))) '%'],['Slow ' num2str(round(max(normCaSL_BNd8X))) '%'])
% axis([1000 13000 0 100])
ylabel("Integrated Ca flux normalized")
xlabel("time (ms)")
title("Healthy Cell + .8X BacNav")

figure;
plot(SBfinalt_BN1X22, normCaTot_BN1X);
hold on;
plot(SBfinalt_BN1X22, normCaSR_BN1X);
plot(SBfinalt_BN1X22, normCaNCX_BN1X);
plot(SBfinalt_BN1X22, normCaSL_BN1X);
% legend( ['Total ' num2str(round(max(normCaTot))) ' Perc. '],'SERCA','INCX','Slow')
legend(['Total ' num2str(round(max(normCaTot_BN1X))) '%'],['SERCA ' num2str(round(max(normCaSR_BN1X))) '%'],['INCX ' num2str(round(max(normCaNCX_BN1X))) '%'],['Slow ' num2str(round(max(normCaSL_BN1X))) '%'])
% axis([1000 13000 0 100])
ylabel("Integrated Ca flux normalized")
xlabel("time (ms)")
title("Healthy Cell + 1X BacNav")






A_Caffeine = { SBfinalt22, normCaTot, SBfinalt22, normCaSR, SBfinalt22, normCaNCX, SBfinalt22, normCaSL, zeros(1000,1),...
               SBfinalt_BNd1X22, normCaTot_BNd1X, SBfinalt_BNd1X22, normCaSR_BNd1X, SBfinalt_BNd1X22, normCaNCX_BNd1X, SBfinalt_BNd1X22, normCaSL_BNd1X, zeros(1000,1),...
               SBfinalt_BNd2X22, normCaTot_BNd2X, SBfinalt_BNd2X22, normCaSR_BNd2X, SBfinalt_BNd2X22, normCaNCX_BNd2X, SBfinalt_BNd2X22, normCaSL_BNd2X, zeros(1000,1),...
               SBfinalt_BNd4X22, normCaTot_BNd4X, SBfinalt_BNd4X22, normCaSR_BNd4X, SBfinalt_BNd4X22, normCaNCX_BNd4X, SBfinalt_BNd4X22, normCaSL_BNd4X, zeros(1000,1),...
               SBfinalt_BNd6X22, normCaTot_BNd6X, SBfinalt_BNd6X22, normCaSR_BNd6X, SBfinalt_BNd6X22, normCaNCX_BNd6X, SBfinalt_BNd6X22, normCaSL_BNd6X, zeros(1000,1),...
               SBfinalt_BNd8X22, normCaTot_BNd8X, SBfinalt_BNd8X22, normCaSR_BNd8X, SBfinalt_BNd8X22, normCaNCX_BNd8X, SBfinalt_BNd8X22, normCaSL_BNd8X, zeros(1000,1),...
               SBfinalt_BN1X22, normCaTot_BN1X, SBfinalt_BN1X22, normCaSR_BN1X, SBfinalt_BN1X22, normCaNCX_BN1X, SBfinalt_BN1X22, normCaSL_BN1X, zeros(1000,1)};

filename = 'CaffeineEffect_CaIntegratedFlux.xlsx';
alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

% Write the arrays to the Excel file, each array in a separate column
for i = 1:length(A_Caffeine)

    if i > 130
        range = [ 'E' alphabet(i-130) num2str(1) ':E' alphabet(i-130) num2str(length(A_Caffeine{i})) ];
    elseif i > 104
        range = [ 'D' alphabet(i-104) num2str(1) ':D' alphabet(i-104) num2str(length(A_Caffeine{i})) ];
    elseif i > 78
        range = [ 'C' alphabet(i-78) num2str(1) ':C' alphabet(i-78) num2str(length(A_Caffeine{i})) ];
    elseif i > 52
        range = [ 'B' alphabet(i-52) num2str(1) ':B' alphabet(i-52) num2str(length(A_Caffeine{i})) ];
    elseif i > 26
        range = [ 'A' alphabet(i-26) num2str(1) ':A' alphabet(i-26) num2str(length(A_Caffeine{i})) ];
    else
        range = [ alphabet(i) num2str(1) ':' alphabet(i) num2str(length(A_Caffeine{i})) ]; 
    end
        

    writematrix(A_Caffeine{i}, filename, 'Range', range );
    % 'Range' specifies the column (A, B, C, etc.) where data is written
    % 'WriteMode' is set to 'append' to add columns for each array
end



toc
%%









% [finalyimportant, finaltimportant, finalcurrentsimportant] = Grandi_human_vclamp_BacNav(0., 0, 300);
% [finalyimportantHF, finaltimportantHF, finalcurrentsimportantHF] = Grandi_human_vclamp_BacNav(0., 1, 300);
% [finalyimportantHF_BNd2X, finaltimportantHF_BNd2X, finalcurrentsimportantHF_BNd2X] = Grandi_human_vclamp_BacNav(0.2, 1, 300);
% [finalyimportantHF_BNd4X, finaltimportantHF_BNd4X, finalcurrentsimportantHF_BNd4X] = Grandi_human_vclamp_BacNav(0.4, 1, 300);
% [finalyimportantHF_BNd6X, finaltimportantHF_BNd6X, finalcurrentsimportantHF_BNd6X] = Grandi_human_vclamp_BacNav(0.6, 1, 300);
% [finalyimportantHF_BNd8X, finaltimportantHF_BNd8X, finalcurrentsimportantHF_BNd8X] = Grandi_human_vclamp_BacNav(0.8, 1, 300);
% [finalyimportantHF_BN1X, finaltimportantHF_BN1X, finalcurrentsimportantHF_BN1X] = Grandi_human_vclamp_BacNav(1., 1, 300);


 
% toc

%Vm
figure
% plot(SBfinaltHF, SBfinalyHF(:,39),'DisplayName','HF')
plot(SBfinalt, SBfinaly(:,39),'DisplayName','Healthy')
hold on
plot(SBfinaltHF_BNd05X, SBfinalyHF_BNd1X(:,39),'DisplayName','HF BacNav .1X')
plot(SBfinaltHF_BNd2X, SBfinalyHF_BNd2X(:,39),'DisplayName','HF BacNav .2X')
plot(SBfinaltHF_BNd4X, SBfinalyHF_BNd4X(:,39),'DisplayName','HF BacNav .4X')
plot(SBfinaltHF_BNd6X, SBfinalyHF_BNd6X(:,39),'DisplayName','HF BacNav .6X')
plot(SBfinaltHF_BNd8X, SBfinalyHF_BNd8X(:,39),'DisplayName','HF BacNav .8X')
plot(SBfinaltHF_BN1X, SBfinalyHF_BN1X(:,39),'DisplayName','HF BacNav 1X')
title('Shannon-Bers Vm in HF with and without BacNav')
xlabel('time (ms)')
ylabel('Vm (mV)')
legend()


%Cai
figure
% plot(SBfinaltHF, SBfinalyHF(:,38),'DisplayName','HF')
hold on
plot(SBfinalt, SBfinaly(:,38),'DisplayName','Healthy')
plot(SBfinaltHF_BNd2X, SBfinalyHF_BNd2X(:,38),'DisplayName','HF BacNav .2X')
plot(SBfinaltHF_BNd4X, SBfinalyHF_BNd4X(:,38),'DisplayName','HF BacNav .4X')
plot(SBfinaltHF_BNd6X, SBfinalyHF_BNd6X(:,38),'DisplayName','HF BacNav .6X')
plot(SBfinaltHF_BNd8X, SBfinalyHF_BNd8X(:,38),'DisplayName','HF BacNav .8X')
plot(SBfinaltHF_BN1X, SBfinalyHF_BN1X(:,38),'DisplayName','HF BacNav 1X')
title('Shannon-Bers Ca_i in HF with and without BacNav')
xlabel('time (ms)')
ylabel('Ca_i (mmol/mL)')
legend()


%ICa
figure
plot(SBfinaltHF, SBfinalcurrentsHF(:,5),'DisplayName','HF')
hold on
plot(SBfinalt, SBfinalcurrents(:,5),'DisplayName','Healthy')
plot(SBfinaltHF_BNd2X, SBfinalcurrentsHF_BNd2X(:,5),'DisplayName','HF BacNav .2X')
plot(SBfinaltHF_BNd4X, SBfinalcurrentsHF_BNd4X(:,5),'DisplayName','HF BacNav .4X')
plot(SBfinaltHF_BNd6X, SBfinalcurrentsHF_BNd6X(:,5),'DisplayName','HF BacNav .6X')
plot(SBfinaltHF_BNd8X, SBfinalcurrentsHF_BNd8X(:,5),'DisplayName','HF BacNav .8X')
plot(SBfinaltHF_BN1X, SBfinalcurrentsHF_BN1X(:,5),'DisplayName','HF BacNav 1X')
title('Shannon-Bers ICa in HF with and without BacNav')
xlabel('time (ms)')
ylabel('ICa (uA/uF)')
legend()


%INaCa
figure
% plot(SBfinaltHF, SBfinalcurrentsHF(:,6),'DisplayName','HF')
plot(SBfinaltHF_BNd05X, SBfinalcurrentsHF_BNd05X(:,39),'DisplayName','HF BacNav .05X')
hold on
plot(SBfinalt, SBfinalcurrents(:,6),'DisplayName','Healthy')
plot(SBfinaltHF_BNd2X, SBfinalcurrentsHF_BNd2X(:,6),'DisplayName','HF BacNav .2X')
plot(SBfinaltHF_BNd4X, SBfinalcurrentsHF_BNd4X(:,6),'DisplayName','HF BacNav .4X')
plot(SBfinaltHF_BNd6X, SBfinalcurrentsHF_BNd6X(:,6),'DisplayName','HF BacNav .6X')
plot(SBfinaltHF_BNd8X, SBfinalcurrentsHF_BNd8X(:,6),'DisplayName','HF BacNav .8X')
% plot(SBfinaltHF_BN1X, SBfinalcurrentsHF_BN1X(:,6),'DisplayName','HF BacNav 1X')
title('Shannon-Bers INaCa in HF with and without BacNav')
xlabel('time (ms)')
ylabel('INaCa (uA/uF)')
legend()



%CaSR
figure
subplot(3,1,1)
plot(SBfinaltHF, SBfinalyHF(:,31),'DisplayName','HF')
hold on
plot(SBfinalt, SBfinaly(:,31),'DisplayName','Healthy')
plot(SBfinaltHF_BNd2X, SBfinalyHF_BNd2X(:,31),'DisplayName','HF BacNav .2X')
plot(SBfinaltHF_BNd4X, SBfinalyHF_BNd4X(:,31),'DisplayName','HF BacNav .4X')
plot(SBfinaltHF_BNd6X, SBfinalyHF_BNd6X(:,31),'DisplayName','HF BacNav .6X')
plot(SBfinaltHF_BNd8X, SBfinalyHF_BNd8X(:,31),'DisplayName','HF BacNav .8X')
plot(SBfinaltHF_BN1X, SBfinalyHF_BN1X(:,31),'DisplayName','HF BacNav 1X')
title('Shannon-Bers Ca in SR in HF with and without BacNav')
xlabel('time (ms)')
ylabel('SR Ca (mmol/mL)')
legend()
subplot(3,1,2)
plot(SBfinaltHF, SBfinalyHF(:,36),'DisplayName','HF')
hold on
plot(SBfinalt, SBfinaly(:,36),'DisplayName','Healthy')
plot(SBfinaltHF_BNd2X, SBfinalyHF_BNd2X(:,36),'DisplayName','HF BacNav .2X')
plot(SBfinaltHF_BNd4X, SBfinalyHF_BNd4X(:,36),'DisplayName','HF BacNav .4X')
plot(SBfinaltHF_BNd6X, SBfinalyHF_BNd6X(:,36),'DisplayName','HF BacNav .6X')
plot(SBfinaltHF_BNd8X, SBfinalyHF_BNd8X(:,36),'DisplayName','HF BacNav .8X')
plot(SBfinaltHF_BN1X, SBfinalyHF_BN1X(:,36),'DisplayName','HF BacNav 1X')
title('Shannon-Bers Junctional Ca in HF with and without BacNav')
xlabel('time (ms)')
ylabel('Junctional Ca (mmol/mL)')
legend()
subplot(3,1,3)
plot(SBfinaltHF, SBfinalyHF(:,37),'DisplayName','HF')
hold on
plot(SBfinalt, SBfinaly(:,37),'DisplayName','Healthy')
plot(SBfinaltHF_BNd2X, SBfinalyHF_BNd2X(:,37),'DisplayName','HF BacNav .2X')
plot(SBfinaltHF_BNd4X, SBfinalyHF_BNd4X(:,37),'DisplayName','HF BacNav .4X')
plot(SBfinaltHF_BNd6X, SBfinalyHF_BNd6X(:,37),'DisplayName','HF BacNav .6X')
plot(SBfinaltHF_BNd8X, SBfinalyHF_BNd8X(:,37),'DisplayName','HF BacNav .8X')
plot(SBfinaltHF_BN1X, SBfinalyHF_BN1X(:,37),'DisplayName','HF BacNav 1X')
title('Shannon-Bers Subsarcolemmal Ca in HF with and without BacNav')
xlabel('time (ms)')
ylabel('Subsarcolemmal Ca (mmol/mL)')
legend()


%IBacNav
figure
% plot(SBfinaltHF, SBfinalcurrentsHF(:,7),'DisplayName','HF')
plot(SBfinaltHF_BNd05X, SBfinalcurrentsHF_BNd05X(:,7),'DisplayName','HF BacNav .05X')
hold on
plot(SBfinalt, SBfinalcurrents(:,7),'DisplayName','Healthy')
plot(SBfinaltHF_BNd2X, SBfinalcurrentsHF_BNd2X(:,7),'DisplayName','HF BacNav .2X')
plot(SBfinaltHF_BNd4X, SBfinalcurrentsHF_BNd4X(:,7),'DisplayName','HF BacNav .4X')
plot(SBfinaltHF_BNd6X, SBfinalcurrentsHF_BNd6X(:,7),'DisplayName','HF BacNav .6X')
plot(SBfinaltHF_BNd8X, SBfinalcurrentsHF_BNd8X(:,7),'DisplayName','HF BacNav .8X')
% plot(SBfinaltHF_BN1X, SBfinalcurrentsHF_BN1X(:,7),'DisplayName','HF BacNav 1X')
title('Shannon-Bers IBacNav in HF with and without BacNav')
xlabel('time (ms)')
ylabel('IBacNav (uA/uF)')
legend()
% 
% %IBacNav
% figure
% plot(SBfinaltHF, SBfinalcurrentsHF(:,7),'DisplayName','HF')
% hold on
% plot(SBfinalt, SBfinalcurrents(:,7),'DisplayName','Healthy')
% plot(SBfinaltHF_BNd2X, SBfinalcurrentsHF_BNd2X(:,7),'DisplayName','HF BacNav .2X')
% plot(SBfinaltHF_BNd4X, SBfinalcurrentsHF_BNd4X(:,7),'DisplayName','HF BacNav .4X')
% plot(SBfinaltHF_BNd6X, SBfinalcurrentsHF_BNd6X(:,7),'DisplayName','HF BacNav .6X')
% plot(SBfinaltHF_BNd8X, SBfinalcurrentsHF_BNd8X(:,7),'DisplayName','HF BacNav .8X')
% plot(SBfinaltHF_BN1X, SBfinalcurrentsHF_BN1X(:,7),'DisplayName','HF BacNav 1X')
% title('Shannon-Bers IBacNav in HF with and without BacNav')
% xlabel('time (ms)')
% ylabel('IBacNav (uA/uF)')
% legend()


%% VOLTAGE-CLAMP CONDUCTIVITIES


vecTestVolt = -120:2:80;

VmoGrandi=-8.09763e+1; 
VmoPB=-90.7796417483135;
VmoShannonBers=-8.556885e+1; 

vecPeakNaCurr_Grandi = [];
vecPeakNaCurr_GrandiBN = [];
vecPeakBacNavNaCurr_Grandi = [];

vecPeakNaCurr_ShannBers = [];
vecPeakNaCurr_ShannBersBN = [];
vecPeakBacNavNaCurr_ShannBers = [];

vecPeakNaCurr_PB = [];
vecPeakNaCurr_PBBN = [];
vecPeakBacNavNaCurr_PB = [];

% [finalyimportant, finaltimportant, finalcurrentsimportant] = Grandi_human_vclamp2(0., 0, 10, 1000, param); %to have 10 beats with cycle length 1000 and define rate to 1Hz
param.HFboolean = 0;

for i = 1:length(vecTestVolt)

    param.vcParameters.V_test = vecTestVolt(i);

    param.vcParameters.tstartclamp = .1;
    param.vcParameters.tendclamp = 2222.1;
    param.vcParameters.R_clamp = .02;
    % param.protocol = 'vclampspecial';%'pace1';%'vclamp';
    param.ratePace = 1.e-3;%.001 is 1Hz since 1/.001=1000 so every 1000ms
    


    % 
    % %Grandi VOLTAGE CLAMP
    % param.vcParameters.V_hold = VmoGrandi;
    % %[finalyimportant, finaltimportant, finalcurrentsimportant] = Grandi_human_vclamp2(percBacNavG, HFbool, numBeats, CycleLength, param2)
    % [Grandiyimportant, Granditimportant, Grandicurrentsimportant] = Grandi_human_vclamp2(0., 0, 1, 3000, param);
    % [GrandiyimportantBN, GranditimportantBN, GrandicurrentsimportantBN] = Grandi_human_vclamp2(1., 0, 1, 3000, param);
    % 
    % vecPeakNaCurr_Grandi(i) = min(Grandicurrentsimportant(:,1));% + max(finalcurrentsimportant(:,1)); %peak value of INa
    % vecPeakNaCurr_GrandiBN(i) = min(GrandicurrentsimportantBN(:,1));% + max(finalcurrentsimportant(:,1)); %peak value of INa
    % vecPeakBacNavNaCurr_Grandi(i) = min(GrandicurrentsimportantBN(:,10));% + max(finalcurrentsimportantBN(:,10)); %peak value of I_BacNav
    % 



    %Shannon-Bers VOLTAGE CLAMP
    param.vcParameters.tstartclamp = .1;
    param.vcParameters.tendclamp = 2000;%909.02;%1111.1;
    param.vcParameters.R_clamp = .02;
    param.protocol = 'vclamp';%'vclampspecial';%'pace1';%'vclamp';
    param.vcParameters.V_hold = VmoShannonBers;%-40.;%VmoShannonBers;
    param.vcParameters.voltageRate = -.22;%-.18;
    param.vcParameters.topVoltage = 80;
    

    [ShannBersyimportant, ShannBerstimportant, ShannBerscurrentsimportant, ShannBersyimportant2, ShannBerstimportant2, ShannBerscurrentsimportant2] = ShannonBers_rabbit(0., 0, 1, param.vcParameters.tendclamp+.1, param, diffRateBool);
    % [ShannBersyimportant_d1XBN, ShannBerstimportant_d1XBN, ShannBerscurrentsimportant_d1XBN, ShannBersyimportant2_d1XBN, ShannBerstimportant2_d1XBN, ShannBerscurrentsimportant2_d1XBN] = ShannonBers_rabbit(.1, 0, 1, param.vcParameters.tendclamp+.1, param, diffRateBool);
    % [ShannBersyimportant_d2XBN, ShannBerstimportant_d2XBN, ShannBerscurrentsimportant_d2XBN, ShannBersyimportant2_d2XBN, ShannBerstimportant2_d2XBN, ShannBerscurrentsimportant2_d2XBN] = ShannonBers_rabbit(.2, 0, 1, param.vcParameters.tendclamp+.1, param, diffRateBool);
    % [ShannBersyimportant_d4XBN, ShannBerstimportant_d4XBN, ShannBerscurrentsimportant_d4XBN, ShannBersyimportant2_d4XBN, ShannBerstimportant2_d4XBN, ShannBerscurrentsimportant2_d4XBN] = ShannonBers_rabbit(.4, 0, 1, param.vcParameters.tendclamp+.1, param, diffRateBool);
    % [ShannBersyimportant_d6XBN, ShannBerstimportant_d6XBN, ShannBerscurrentsimportant_d6XBN, ShannBersyimportant2_d6XBN, ShannBerstimportant2_d6XBN, ShannBerscurrentsimportant2_d6XBN] = ShannonBers_rabbit(.6, 0, 1, param.vcParameters.tendclamp+.1, param, diffRateBool);
    % [ShannBersyimportant_d8XBN, ShannBerstimportant_d8XBN, ShannBerscurrentsimportant_d8XBN, ShannBersyimportant2_d8XBN, ShannBerstimportant2_d8XBN, ShannBerscurrentsimportant2_d8XBN] = ShannonBers_rabbit(.8, 0, 1, param.vcParameters.tendclamp+.1, param, diffRateBool);
    % [ShannBersyimportant_1XBN, ShannBerstimportant_1XBN, ShannBerscurrentsimportant_1XBN, ShannBersyimportant2_1XBN, ShannBerstimportant2_1XBN, ShannBerscurrentsimportant2_1XBN] = ShannonBers_rabbit(1.0, 0, 1, param.vcParameters.tendclamp+.1, param, diffRateBool);
    % 
    vecPeakNaCurr_ShannBers(i) = min(ShannBerscurrentsimportant(:,6))+ max(ShannBerscurrentsimportant(:,6)); %peak value of INCX
    % vecPeakNaCurr_ShannBersBN(i) = min(ShannBerscurrentsimportantBN(:,1));% + max(finalcurrentsimportant(:,1)); %peak value of INa
    % vecPeakBacNavNaCurr_ShannBers(i) = min(ShannBerscurrentsimportantBN(:,7));% + max(finalcurrentsimportantBN(:,10)); %peak value of I_BacNav

    % 



    % %TESTING START
    % Vmtesting = [];
    % for i = 1:length(t)
    %     if t(i) < 10
    %         Vmtesting(i) = -40;
    %     elseif t(i) >= 10 & t(i) <= 2232
    %         Vmtesting(i) = 80 - (t(i) - 10)*.09;
    %     elseif t(i) > 2232
    %         Vmtesting(i) = -40;
    %     end
    % end
    % 
    % indiShB = find( ShannBerstimportant > param.vcParameters.tstartclamp & ShannBerstimportant < param.vcParameters.tendclamp );
    % indiShB_d1XBN = find( ShannBerstimportant_d1XBN > param.vcParameters.tstartclamp & ShannBerstimportant_d1XBN < param.vcParameters.tendclamp );
    % indiShB_d2XBN = find( ShannBerstimportant_d2XBN > param.vcParameters.tstartclamp & ShannBerstimportant_d2XBN < param.vcParameters.tendclamp );
    % indiShB_d4XBN = find( ShannBerstimportant_d4XBN > param.vcParameters.tstartclamp & ShannBerstimportant_d4XBN < param.vcParameters.tendclamp );
    % indiShB_d6XBN = find( ShannBerstimportant_d6XBN > param.vcParameters.tstartclamp & ShannBerstimportant_d6XBN < param.vcParameters.tendclamp );
    % indiShB_d8XBN = find( ShannBerstimportant_d8XBN > param.vcParameters.tstartclamp & ShannBerstimportant_d8XBN < param.vcParameters.tendclamp );
    % indiShB_1XBN = find( ShannBerstimportant_1XBN > param.vcParameters.tstartclamp & ShannBerstimportant_1XBN < param.vcParameters.tendclamp );

    % 
    % figure;
    % subplot(3,1,1);
    % plot(ShannBerstimportant(indiShB), ShannBersyimportant(indiShB,39), 'DisplayName','Healthy');
    % hold on
    % plot(ShannBerstimportant_d1XBN(indiShB_d1XBN), ShannBersyimportant_d1XBN(indiShB_d1XBN,39), 'DisplayName','BN .1X');
    % plot(ShannBerstimportant_d2XBN(indiShB_d2XBN), ShannBersyimportant_d2XBN(indiShB_d2XBN,39), 'DisplayName','BN .2X');
    % plot(ShannBerstimportant_d4XBN(indiShB_d4XBN), ShannBersyimportant_d4XBN(indiShB_d4XBN,39), 'DisplayName','BN .4X');
    % plot(ShannBerstimportant_d6XBN(indiShB_d6XBN), ShannBersyimportant_d6XBN(indiShB_d6XBN,39), 'DisplayName','BN .6X');
    % plot(ShannBerstimportant_d8XBN(indiShB_d8XBN), ShannBersyimportant_d8XBN(indiShB_d8XBN,39), 'DisplayName','BN .8X');
    % plot(ShannBerstimportant_1XBN(indiShB_1XBN), ShannBersyimportant_1XBN(indiShB_1XBN,39), 'DisplayName','BN 1X');
    % legend()
    % title('Vm');
    % xlabel('time (ms)')
    % ylabel('Vm')
    % subplot(3,1,2);
    % plot(ShannBerstimportant(indiShB), ShannBerscurrentsimportant(indiShB,6), 'DisplayName','Healthy');
    % hold on
    % plot(ShannBerstimportant_d1XBN(indiShB_d1XBN), ShannBerscurrentsimportant_d1XBN(indiShB_d1XBN,6), 'DisplayName','BN .1X');
    % plot(ShannBerstimportant_d2XBN(indiShB_d2XBN), ShannBerscurrentsimportant_d2XBN(indiShB_d2XBN,6), 'DisplayName','BN .2X');
    % plot(ShannBerstimportant_d4XBN(indiShB_d4XBN), ShannBerscurrentsimportant_d4XBN(indiShB_d4XBN,6), 'DisplayName','BN .4X');
    % plot(ShannBerstimportant_d6XBN(indiShB_d6XBN), ShannBerscurrentsimportant_d6XBN(indiShB_d6XBN,6), 'DisplayName','BN .6X');
    % plot(ShannBerstimportant_d8XBN(indiShB_d8XBN), ShannBerscurrentsimportant_d8XBN(indiShB_d8XBN,6), 'DisplayName','BN .8X');
    % plot(ShannBerstimportant_1XBN(indiShB_1XBN), ShannBerscurrentsimportant_1XBN(indiShB_1XBN,6), 'DisplayName','BN 1X');
    % legend()
    % title('INCX');
    % xlabel('time (ms)')
    % ylabel('I_NCX')
    % subplot(3,1,3);
    % plot(ShannBersyimportant(indiShB,39), ShannBerscurrentsimportant(indiShB,6), 'DisplayName','Healthy');
    % hold on
    % plot(ShannBersyimportant_d1XBN(indiShB_d1XBN,39), ShannBerscurrentsimportant_d1XBN(indiShB_d1XBN,6), 'DisplayName','BN .1X');
    % plot(ShannBersyimportant_d2XBN(indiShB_d2XBN,39), ShannBerscurrentsimportant_d2XBN(indiShB_d2XBN,6), 'DisplayName','BN .2X');
    % plot(ShannBersyimportant_d4XBN(indiShB_d4XBN,39), ShannBerscurrentsimportant_d4XBN(indiShB_d4XBN,6), 'DisplayName','BN .4X');
    % plot(ShannBersyimportant_d6XBN(indiShB_d6XBN,39), ShannBerscurrentsimportant_d6XBN(indiShB_d6XBN,6), 'DisplayName','BN .6X');
    % plot(ShannBersyimportant_d8XBN(indiShB_d8XBN,39), ShannBerscurrentsimportant_d8XBN(indiShB_d8XBN,6), 'DisplayName','BN .8X');
    % plot(ShannBersyimportant_1XBN(indiShB_1XBN,39), ShannBerscurrentsimportant_1XBN(indiShB_1XBN,6), 'DisplayName','BN 1X');
    % legend()
    % title('I-V');
    % xlabel('Vm')
    % ylabel('I_NCX')

    %TESTING END



    % 
    % 
    % %PB VOLTAGE CLAMP
    % param.vcParameters.V_hold = VmoPB;
    % param.BCL = 3000.;
    % param.percBacNav = 0.;
    % [VOI, STATES, ALGEBRAIC, CONSTANTS] = PriebeBeuckelmannFun(param, 0, 3000);
    % param.g_Na_BacNav = 12.6927; %PB BacNav conductance
    % param.percBacNav = 1.;
    % [VOIBN, STATESBN, ALGEBRAICBN, CONSTANTSBN] = PriebeBeuckelmannFun(param, 0, 3000);
    % 
    % vecPeakNaCurr_PB(i) = min(ALGEBRAIC(:,22));% + max(finalcurrentsimportant(:,1)); %peak value of INa
    % vecPeakNaCurr_PBBN(i) = min(ALGEBRAICBN(:,22));% + max(finalcurrentsimportant(:,1)); %peak value of INa
    % vecPeakBacNavNaCurr_PB(i) = min(ALGEBRAICBN(:,52));% + max(finalcurrentsimportantBN(:,10)); %peak value of I_BacNav
    % 
end


figure
plot(vecTestVolt, vecPeakNaCurr_Grandi)
hold on
plot(vecTestVolt, vecPeakBacNavNaCurr_Grandi)
plot(vecTestVolt, vecPeakNaCurr_GrandiBN)
legend('INa','IBacNav','INa with BN')
title('Grandi I-V curve with BacNav max conductance of 32.3936')
xlabel('Voltage (mV)')
ylabel('Current (uA/uF)')

figure
plot(vecTestVolt, vecPeakNaCurr_ShannBers)
hold on
plot(vecTestVolt, vecPeakBacNavNaCurr_ShannBers)
plot(vecTestVolt, vecPeakNaCurr_ShannBersBN)
legend('INa','IBacNav','INa with BN')
title('Shannon-Bers I-V curve with BacNav max conductance of 14.1722')
xlabel('Voltage (mV)')
ylabel('Current (uA/uF)')

figure
plot(vecTestVolt, vecPeakNaCurr_PB)
hold on
plot(vecTestVolt, vecPeakNaCurr_PBBN)
plot(vecTestVolt, vecPeakBacNavNaCurr_PB)
legend('INa','IBacNav','INa with BN')
title('PB I-V curve with BacNav max conductance of 12.6927')
xlabel('Voltage (mV)')
ylabel('Current (uA/uF)')















%Vm
figure;
% plot(VOIHF(indiHF2),STATESHF(indiHF2,1), 'DisplayName','HF');
hold on;
plot(VOI(indi2),STATES(indi2,1), 'DisplayName','Healthy');
% plot(VOIHF_BNd1X(indiHF2_BNd1X),STATESHF_BNd1X(indiHF2_BNd1X,1), 'DisplayName','HF BacNav .1X');
plot(VOIHF_BNd05X(indiHF2_BNd05X),STATESHF_BNd05X(indiHF2_BNd05X,1), 'DisplayName','BacNav .05X');
plot(VOIHF_BNd2X(indiHF2_BNd2X),STATESHF_BNd2X(indiHF2_BNd2X,1), 'DisplayName','BacNav .2X');
plot(VOIHF_BNd4X(indiHF2_BNd4X),STATESHF_BNd4X(indiHF2_BNd4X,1), 'DisplayName','BacNav .4X');
% plot(VOIHF_BNd5X(indiHF2_BNd5X),STATESHF_BNd5X(indiHF2_BNd5X,1), 'DisplayName','HF BacNav .5X');
plot(VOIHF_BNd6X(indiHF2_BNd6X),STATESHF_BNd6X(indiHF2_BNd6X,1), 'DisplayName','BacNav .6X');
plot(VOIHF_BNd8X(indiHF2_BNd8X),STATESHF_BNd8X(indiHF2_BNd8X,1), 'DisplayName','BacNav .8X');
plot(VOIHF_BN1X(indiHF2_BN1X),STATESHF_BN1X(indiHF2_BN1X,1),'DisplayName','BacNav 1X');
title('Vm with and without BacNav')
xlabel('time (ms)')
ylabel('Vm (mV)')
legend()


%Cai
figure;
% plot(VOIHF(indiHF2),STATESHF(indiHF2,6), 'DisplayName','HF');
plot(VOIHF_BNd05X(indiHF2_BNd05X),STATESHF_BNd05X(indiHF2_BNd05X,6), 'DisplayName','BacNav .05X');
hold on;
plot(VOI(indi2),STATES(indi2,6), 'DisplayName','Healthy');
% plot(VOIHF_BNd1X(indiHF2_BNd1X),STATESHF_BNd1X(indiHF2_BNd1X,6), 'DisplayName','HF BacNav .1X');
plot(VOIHF_BNd2X(indiHF2_BNd2X),STATESHF_BNd2X(indiHF2_BNd2X,6), 'DisplayName','BacNav .2X');
plot(VOIHF_BNd4X(indiHF2_BNd4X),STATESHF_BNd4X(indiHF2_BNd4X,6), 'DisplayName','BacNav .4X');
% plot(VOIHF_BNd5X(indiHF2_BNd5X),STATESHF_BNd5X(indiHF2_BNd5X,6), 'DisplayName','HF BacNav .5X');
plot(VOIHF_BNd6X(indiHF2_BNd6X),STATESHF_BNd6X(indiHF2_BNd6X,6), 'DisplayName','BacNav .6X');
plot(VOIHF_BNd8X(indiHF2_BNd8X),STATESHF_BNd8X(indiHF2_BNd8X,6), 'DisplayName','BacNav .8X');
plot(VOIHF_BN1X(indiHF2_BN1X),STATESHF_BN1X(indiHF2_BN1X,6), 'DisplayName','BacNav 1X');
title('Ca_i concentration with and without BacNav')
xlabel('time (ms)')
ylabel('Ca_i (mmol/mL)')
legend()









%Vm
figure
% plot(SBfinaltHF, SBfinalyHF(:,39),'DisplayName','HF')
plot(SBfinalt, SBfinaly(:,39),'DisplayName','Healthy')
hold on
plot(SBfinaltHF_BNd05X, SBfinalyHF_BNd05X(:,39),'DisplayName','BacNav .1X')
plot(SBfinaltHF_BNd2X, SBfinalyHF_BNd2X(:,39),'DisplayName','BacNav .2X')
plot(SBfinaltHF_BNd4X, SBfinalyHF_BNd4X(:,39),'DisplayName','BacNav .4X')
plot(SBfinaltHF_BNd6X, SBfinalyHF_BNd6X(:,39),'DisplayName','BacNav .6X')
plot(SBfinaltHF_BNd8X, SBfinalyHF_BNd8X(:,39),'DisplayName','BacNav .8X')
% plot(SBfinaltHF_BN1X, SBfinalyHF_BN1X(:,39),'DisplayName','HF BacNav 1X')
title('Shannon-Bers Vm with and without BacNav')
xlabel('time (ms)')
ylabel('Vm (mV)')
legend()


%Cai
figure
% plot(SBfinaltHF, SBfinalyHF(:,38),'DisplayName','HF')
plot(SBfinaltHF_BNd05X, SBfinalyHF_BNd05X(:,38),'DisplayName','BacNav .05X')
hold on
plot(SBfinalt, SBfinaly(:,38),'DisplayName','Healthy')
plot(SBfinaltHF_BNd2X, SBfinalyHF_BNd2X(:,38),'DisplayName','BacNav .2X')
plot(SBfinaltHF_BNd4X, SBfinalyHF_BNd4X(:,38),'DisplayName','BacNav .4X')
plot(SBfinaltHF_BNd6X, SBfinalyHF_BNd6X(:,38),'DisplayName','BacNav .6X')
plot(SBfinaltHF_BNd8X, SBfinalyHF_BNd8X(:,38),'DisplayName','BacNav .8X')
% plot(SBfinaltHF_BN1X, SBfinalyHF_BN1X(:,38),'DisplayName','BacNav 1X')
title('Shannon-Bers Ca_i with and without BacNav')
xlabel('time (ms)')
ylabel('Ca_i (mmol/mL)')
legend()











%Vm
figure
hold on
plot(finaltimportant, finalyimportant(:,39),'DisplayName','Healthy')
plot(finaltimportantHF_BNd025X, finalyimportantHF_BNd025X(:,39),'DisplayName','BacNav .025X')
plot(finaltimportantHF_BNd05X, finalyimportantHF_BNd05X(:,39),'DisplayName','BacNav .05X')
plot(finaltimportantHF_BNd075X, finalyimportantHF_BNd075X(:,39),'DisplayName','BacNav .075X')
plot(finaltimportantHF_BNd1X, finalyimportantHF_BNd1X(:,39),'DisplayName','BacNav .1X')
plot(finaltimportantHF_BNd15X, finalyimportantHF_BNd15X(:,39),'DisplayName','BacNav .15X')
plot(finaltimportantHF_BNd2X, finalyimportantHF_BNd2X(:,39),'DisplayName','BacNav .2X')
% plot(finaltimportantHF_BNd4X, finalyimportantHF_BNd4X(:,39),'DisplayName','BacNav .4X')
% plot(finaltimportantHF_BNd6X, finalyimportantHF_BNd6X(:,39),'DisplayName','HF BacNav .6X')
% plot(finaltimportantHF_BNd8X, finalyimportantHF_BNd8X(:,39),'DisplayName','HF BacNav .8X')
% plot(finaltimportantHF_BN1X, finalyimportantHF_BN1X(:,39),'DisplayName','HF BacNav 1X')
title('Grandi Vm with and without BacNav')
xlabel('time (ms)')
ylabel('Vm (mV)')
legend()


%Cai
figure
hold on
plot(finaltimportant, finalyimportant(:,38),'DisplayName','Healthy')
plot(finaltimportantHF_BNd025X, finalyimportantHF_BNd025X(:,38),'DisplayName','BacNav .025X')
plot(finaltimportantHF_BNd05X, finalyimportantHF_BNd05X(:,38),'DisplayName','BacNav .05X')
plot(finaltimportantHF_BNd075X, finalyimportantHF_BNd075X(:,38),'DisplayName','BacNav .075X')
plot(finaltimportantHF_BNd1X, finalyimportantHF_BNd1X(:,38),'DisplayName','BacNav .1X')
plot(finaltimportantHF_BNd15X, finalyimportantHF_BNd15X(:,38),'DisplayName','BacNav .15X')
plot(finaltimportantHF_BNd2X, finalyimportantHF_BNd2X(:,38),'DisplayName','BacNav .2X')
% plot(finaltimportantHF_BNd4X, finalyimportantHF_BNd4X(:,38),'DisplayName','BacNav .4X')
% plot(finaltimportantHF_BNd6X, finalyimportantHF_BNd6X(:,38),'DisplayName','HF BacNav .6X')
% plot(finaltimportantHF_BNd8X, finalyimportantHF_BNd8X(:,38),'DisplayName','HF BacNav .8X')
% plot(finaltimportantHF_BN1X, finalyimportantHF_BN1X(:,38),'DisplayName','HF BacNav 1X')
title('Grandi Ca_i with and without BacNav')
xlabel('time (ms)')
ylabel('Ca_i (mmol/mL)')
legend()








%% CHANGING K1 CONDUCTANCE

tic

%CAFFEINE EFFECT
param.caffeineBool = false;

param.protocol = 'pace1';%'vclamp';
param.ratePace = 1.e-3;%.001 is 1Hz since 1/.001=1000 so every 1000ms

param.GlycosidesBool = false;
numBeats = 1015;
CycleLength = 1000;
diffRateBool = false;

param.DADbool = true;

%Multiplier = 1.
param.GKir = 1.;
%Healthy
[SBfinalyKir1, SBfinaltKir1, SBfinalcurrentsKir1, timeImportKir1, yImportKir1, currentsImportKir1] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);
%BacNav .25X and .5X and 1X
[SBfinaly_BNd1XKir1, SBfinalt_BNd1XKir1, SBfinalcurrents_BNd1XKir1, timeImport_BNd1XKir1, yImport_BNd1XKir1, currentsImport_BNd1XKir1] = ShannonBers_rabbit(0.1, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd2XKir1, SBfinalt_BNd2XKir1, SBfinalcurrents_BNd2XKir1, timeImport_BNd2XKir1, yImport_BNd2XKir1, currentsImport_BNd2XKir1] = ShannonBers_rabbit(0.2, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd4XKir1, SBfinalt_BNd4XKir1, SBfinalcurrents_BNd4XKir1, timeImport_BNd4XKir1, yImport_BNd4XKir1, currentsImport_BNd4XKir1] = ShannonBers_rabbit(0.4, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd6XKir1, SBfinalt_BNd6XKir1, SBfinalcurrents_BNd6XKir1, timeImport_BNd6XKir1, yImport_BNd6XKir1, currentsImport_BNd6XKir1] = ShannonBers_rabbit(0.6, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd8XKir1, SBfinalt_BNd8XKir1, SBfinalcurrents_BNd8XKir1, timeImport_BNd8XKir1, yImport_BNd8XKir1, currentsImport_BNd8XKir1] = ShannonBers_rabbit(0.8, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BN1XKir1, SBfinalt_BN1XKir1, SBfinalcurrents_BN1XKir1, timeImport_BN1XKir1, yImport_BN1XKir1, currentsImport_BN1XKir1] = ShannonBers_rabbit(1., 0, numBeats, CycleLength, param, diffRateBool);



%Multiplier = .2
param.GKir = .2;
%Healthy
[SBfinalyKird2, SBfinaltKird2, SBfinalcurrentsKird2, timeImportKird2, yImportKird2, currentsImportKird2] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);
%BacNav .25X and .5X and 1X
[SBfinaly_BNd1XKird2, SBfinalt_BNd1XKird2, SBfinalcurrents_BNd1XKird2, timeImport_BNd1XKird2, yImport_BNd1XKird2, currentsImport_BNd1XKird2] = ShannonBers_rabbit(0.1, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd2XKird2, SBfinalt_BNd2XKird2, SBfinalcurrents_BNd2XKird2, timeImport_BNd2XKird2, yImport_BNd2XKird2, currentsImport_BNd2XKird2] = ShannonBers_rabbit(0.2, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd4XKird2, SBfinalt_BNd4XKird2, SBfinalcurrents_BNd4XKird2, timeImport_BNd4XKird2, yImport_BNd4XKird2, currentsImport_BNd4XKird2] = ShannonBers_rabbit(0.4, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd6XKird2, SBfinalt_BNd6XKird2, SBfinalcurrents_BNd6XKird2, timeImport_BNd6XKird2, yImport_BNd6XKird2, currentsImport_BNd6XKird2] = ShannonBers_rabbit(0.6, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd8XKird2, SBfinalt_BNd8XKird2, SBfinalcurrents_BNd8XKird2, timeImport_BNd8XKird2, yImport_BNd8XKird2, currentsImport_BNd8XKird2] = ShannonBers_rabbit(0.8, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BN1XKird2, SBfinalt_BN1XKird2, SBfinalcurrents_BN1XKird2, timeImport_BN1XKird2, yImport_BN1XKird2, currentsImport_BN1XKird2] = ShannonBers_rabbit(1., 0, numBeats, CycleLength, param, diffRateBool);




%Multiplier = .4
param.GKir = .4;
%Healthy
[SBfinalyKird4, SBfinaltKird4, SBfinalcurrentsKird4, timeImportKird4, yImportKird4, currentsImportKird4] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);
%BacNav .22X and .2X and 1X
[SBfinaly_BNd1XKird4, SBfinalt_BNd1XKird4, SBfinalcurrents_BNd1XKird4, timeImport_BNd1XKird4, yImport_BNd1XKird4, currentsImport_BNd1XKird4] = ShannonBers_rabbit(0.1, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd2XKird4, SBfinalt_BNd2XKird4, SBfinalcurrents_BNd2XKird4, timeImport_BNd2XKird4, yImport_BNd2XKird4, currentsImport_BNd2XKird4] = ShannonBers_rabbit(0.2, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd4XKird4, SBfinalt_BNd4XKird4, SBfinalcurrents_BNd4XKird4, timeImport_BNd4XKird4, yImport_BNd4XKird4, currentsImport_BNd4XKird4] = ShannonBers_rabbit(0.4, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd6XKird4, SBfinalt_BNd6XKird4, SBfinalcurrents_BNd6XKird4, timeImport_BNd6XKird4, yImport_BNd6XKird4, currentsImport_BNd6XKird4] = ShannonBers_rabbit(0.6, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd8XKird4, SBfinalt_BNd8XKird4, SBfinalcurrents_BNd8XKird4, timeImport_BNd8XKird4, yImport_BNd8XKird4, currentsImport_BNd8XKird4] = ShannonBers_rabbit(0.8, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BN1XKird4, SBfinalt_BN1XKird4, SBfinalcurrents_BN1XKird4, timeImport_BN1XKird4, yImport_BN1XKird4, currentsImport_BN1XKird4] = ShannonBers_rabbit(1., 0, numBeats, CycleLength, param, diffRateBool);




%Multiplier = .6
param.GKir = .6;
%Healthy
[SBfinalyKird6, SBfinaltKird6, SBfinalcurrentsKird6, timeImportKird6, yImportKird6, currentsImportKird6] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);
%BacNav .3X and .5X and 1X
[SBfinaly_BNd1XKird6, SBfinalt_BNd1XKird6, SBfinalcurrents_BNd1XKird6, timeImport_BNd1XKird6, yImport_BNd1XKird6, currentsImport_BNd1XKird6] = ShannonBers_rabbit(0.1, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd2XKird6, SBfinalt_BNd2XKird6, SBfinalcurrents_BNd2XKird6, timeImport_BNd2XKird6, yImport_BNd2XKird6, currentsImport_BNd2XKird6] = ShannonBers_rabbit(0.2, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd4XKird6, SBfinalt_BNd4XKird6, SBfinalcurrents_BNd4XKird6, timeImport_BNd4XKird6, yImport_BNd4XKird6, currentsImport_BNd4XKird6] = ShannonBers_rabbit(0.4, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd6XKird6, SBfinalt_BNd6XKird6, SBfinalcurrents_BNd6XKird6, timeImport_BNd6XKird6, yImport_BNd6XKird6, currentsImport_BNd6XKird6] = ShannonBers_rabbit(0.6, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd8XKird6, SBfinalt_BNd8XKird6, SBfinalcurrents_BNd8XKird6, timeImport_BNd8XKird6, yImport_BNd8XKird6, currentsImport_BNd8XKird6] = ShannonBers_rabbit(0.8, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BN1XKird6, SBfinalt_BN1XKird6, SBfinalcurrents_BN1XKird6, timeImport_BN1XKird6, yImport_BN1XKird6, currentsImport_BN1XKird6] = ShannonBers_rabbit(1., 0, numBeats, CycleLength, param, diffRateBool);




%Multiplier = .8
param.GKir = .8;
%Healthy
[SBfinalyKird8, SBfinaltKird8, SBfinalcurrentsKird8, timeImportKird8, yImportKird8, currentsImportKird8] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);
%BacNav .25X and .5X and 1X
[SBfinaly_BNd1XKird8, SBfinalt_BNd1XKird8, SBfinalcurrents_BNd1XKird8, timeImport_BNd1XKird8, yImport_BNd1XKird8, currentsImport_BNd1XKird8] = ShannonBers_rabbit(0.1, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd2XKird8, SBfinalt_BNd2XKird8, SBfinalcurrents_BNd2XKird8, timeImport_BNd2XKird8, yImport_BNd2XKird8, currentsImport_BNd2XKird8] = ShannonBers_rabbit(0.2, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd4XKird8, SBfinalt_BNd4XKird8, SBfinalcurrents_BNd4XKird8, timeImport_BNd4XKird8, yImport_BNd4XKird8, currentsImport_BNd4XKird8] = ShannonBers_rabbit(0.4, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd6XKird8, SBfinalt_BNd6XKird8, SBfinalcurrents_BNd6XKird8, timeImport_BNd6XKird8, yImport_BNd6XKird8, currentsImport_BNd6XKird8] = ShannonBers_rabbit(0.6, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd8XKird8, SBfinalt_BNd8XKird8, SBfinalcurrents_BNd8XKird8, timeImport_BNd8XKird8, yImport_BNd8XKird8, currentsImport_BNd8XKird8] = ShannonBers_rabbit(0.8, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BN1XKird8, SBfinalt_BN1XKird8, SBfinalcurrents_BN1XKird8, timeImport_BN1XKird8, yImport_BN1XKird8, currentsImport_BN1XKird8] = ShannonBers_rabbit(1., 0, numBeats, CycleLength, param, diffRateBool);




%Multiplier = 1.2
param.GKir = 1.2;
%Healthy
[SBfinalyKir1d2, SBfinaltKir1d2, SBfinalcurrentsKir1d2, timeImportKir1d2, yImportKir1d2, currentsImportKir1d2] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);
%BacNav .25X and .5X and 1X
[SBfinaly_BNd1XKir1d2, SBfinalt_BNd1XKir1d2, SBfinalcurrents_BNd1XKir1d2, timeImport_BNd1XKir1d2, yImport_BNd1XKir1d2, currentsImport_BNd1XKir1d2] = ShannonBers_rabbit(0.1, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd2XKir1d2, SBfinalt_BNd2XKir1d2, SBfinalcurrents_BNd2XKir1d2, timeImport_BNd2XKir1d2, yImport_BNd2XKir1d2, currentsImport_BNd2XKir1d2] = ShannonBers_rabbit(0.2, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd4XKir1d2, SBfinalt_BNd4XKir1d2, SBfinalcurrents_BNd4XKir1d2, timeImport_BNd4XKir1d2, yImport_BNd4XKir1d2, currentsImport_BNd4XKir1d2] = ShannonBers_rabbit(0.4, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd6XKir1d2, SBfinalt_BNd6XKir1d2, SBfinalcurrents_BNd6XKir1d2, timeImport_BNd6XKir1d2, yImport_BNd6XKir1d2, currentsImport_BNd6XKir1d2] = ShannonBers_rabbit(0.6, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd8XKir1d2, SBfinalt_BNd8XKir1d2, SBfinalcurrents_BNd8XKir1d2, timeImport_BNd8XKir1d2, yImport_BNd8XKir1d2, currentsImport_BNd8XKir1d2] = ShannonBers_rabbit(0.8, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BN1XKir1d2, SBfinalt_BN1XKir1d2, SBfinalcurrents_BN1XKir1d2, timeImport_BN1XKir1d2, yImport_BN1XKir1d2, currentsImport_BN1XKir1d2] = ShannonBers_rabbit(1., 0, numBeats, CycleLength, param, diffRateBool);




%Multiplier = 1.4
param.GKir = 1.4;
%Healthy
[SBfinalyKir1d4, SBfinaltKir1d4, SBfinalcurrentsKir1d4, timeImportKir1d4, yImportKir1d4, currentsImportKir1d4] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);
%BacNav .25X and .5X and 1X
[SBfinaly_BNd1XKir1d4, SBfinalt_BNd1XKir1d4, SBfinalcurrents_BNd1XKir1d4, timeImport_BNd1XKir1d4, yImport_BNd1XKir1d4, currentsImport_BNd1XKir1d4] = ShannonBers_rabbit(0.1, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd2XKir1d4, SBfinalt_BNd2XKir1d4, SBfinalcurrents_BNd2XKir1d4, timeImport_BNd2XKir1d4, yImport_BNd2XKir1d4, currentsImport_BNd2XKir1d4] = ShannonBers_rabbit(0.2, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd4XKir1d4, SBfinalt_BNd4XKir1d4, SBfinalcurrents_BNd4XKir1d4, timeImport_BNd4XKir1d4, yImport_BNd4XKir1d4, currentsImport_BNd4XKir1d4] = ShannonBers_rabbit(0.4, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd6XKir1d4, SBfinalt_BNd6XKir1d4, SBfinalcurrents_BNd6XKir1d4, timeImport_BNd6XKir1d4, yImport_BNd6XKir1d4, currentsImport_BNd6XKir1d4] = ShannonBers_rabbit(0.6, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd8XKir1d4, SBfinalt_BNd8XKir1d4, SBfinalcurrents_BNd8XKir1d4, timeImport_BNd8XKir1d4, yImport_BNd8XKir1d4, currentsImport_BNd8XKir1d4] = ShannonBers_rabbit(0.8, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BN1XKir1d4, SBfinalt_BN1XKir1d4, SBfinalcurrents_BN1XKir1d4, timeImport_BN1XKir1d4, yImport_BN1XKir1d4, currentsImport_BN1XKir1d4] = ShannonBers_rabbit(1., 0, numBeats, CycleLength, param, diffRateBool);




%Multiplier = 1.6
param.GKir = 1.6;
%Healthy
[SBfinalyKir1d6, SBfinaltKir1d6, SBfinalcurrentsKir1d6, timeImportKir1d6, yImportKir1d6, currentsImportKir1d6] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);
%BacNav .25X and .5X and 1X
[SBfinaly_BNd1XKir1d6, SBfinalt_BNd1XKir1d6, SBfinalcurrents_BNd1XKir1d6, timeImport_BNd1XKir1d6, yImport_BNd1XKir1d6, currentsImport_BNd1XKir1d6] = ShannonBers_rabbit(0.1, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd2XKir1d6, SBfinalt_BNd2XKir1d6, SBfinalcurrents_BNd2XKir1d6, timeImport_BNd2XKir1d6, yImport_BNd2XKir1d6, currentsImport_BNd2XKir1d6] = ShannonBers_rabbit(0.2, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd4XKir1d6, SBfinalt_BNd4XKir1d6, SBfinalcurrents_BNd4XKir1d6, timeImport_BNd4XKir1d6, yImport_BNd4XKir1d6, currentsImport_BNd4XKir1d6] = ShannonBers_rabbit(0.4, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd6XKir1d6, SBfinalt_BNd6XKir1d6, SBfinalcurrents_BNd6XKir1d6, timeImport_BNd6XKir1d6, yImport_BNd6XKir1d6, currentsImport_BNd6XKir1d6] = ShannonBers_rabbit(0.6, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd8XKir1d6, SBfinalt_BNd8XKir1d6, SBfinalcurrents_BNd8XKir1d6, timeImport_BNd8XKir1d6, yImport_BNd8XKir1d6, currentsImport_BNd8XKir1d6] = ShannonBers_rabbit(0.8, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BN1XKir1d6, SBfinalt_BN1XKir1d6, SBfinalcurrents_BN1XKir1d6, timeImport_BN1XKir1d6, yImport_BN1XKir1d6, currentsImport_BN1XKir1d6] = ShannonBers_rabbit(1., 0, numBeats, CycleLength, param, diffRateBool);




%Multiplier = 1.8
param.GKir = 1.8;
%Healthy
[SBfinalyKir1d8, SBfinaltKir1d8, SBfinalcurrentsKir1d8, timeImportKir1d8, yImportKir1d8, currentsImportKir1d8] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);
%BacNav .25X and .5X and 1X
[SBfinaly_BNd1XKir1d8, SBfinalt_BNd1XKir1d8, SBfinalcurrents_BNd1XKir1d8, timeImport_BNd1XKir1d8, yImport_BNd1XKir1d8, currentsImport_BNd1XKir1d8] = ShannonBers_rabbit(0.1, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd2XKir1d8, SBfinalt_BNd2XKir1d8, SBfinalcurrents_BNd2XKir1d8, timeImport_BNd2XKir1d8, yImport_BNd2XKir1d8, currentsImport_BNd2XKir1d8] = ShannonBers_rabbit(0.2, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd4XKir1d8, SBfinalt_BNd4XKir1d8, SBfinalcurrents_BNd4XKir1d8, timeImport_BNd4XKir1d8, yImport_BNd4XKir1d8, currentsImport_BNd4XKir1d8] = ShannonBers_rabbit(0.4, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd6XKir1d8, SBfinalt_BNd6XKir1d8, SBfinalcurrents_BNd6XKir1d8, timeImport_BNd6XKir1d8, yImport_BNd6XKir1d8, currentsImport_BNd6XKir1d8] = ShannonBers_rabbit(0.6, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd8XKir1d8, SBfinalt_BNd8XKir1d8, SBfinalcurrents_BNd8XKir1d8, timeImport_BNd8XKir1d8, yImport_BNd8XKir1d8, currentsImport_BNd8XKir1d8] = ShannonBers_rabbit(0.8, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BN1XKir1d8, SBfinalt_BN1XKir1d8, SBfinalcurrents_BN1XKir1d8, timeImport_BN1XKir1d8, yImport_BN1XKir1d8, currentsImport_BN1XKir1d8] = ShannonBers_rabbit(1., 0, numBeats, CycleLength, param, diffRateBool);




%Multiplier = 2.
param.GKir = 2.;
%Healthy
[SBfinalyKir2, SBfinaltKir2, SBfinalcurrentsKir2, timeImportKir2, yImportKir2, currentsImportKir2] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);
%BacNav .25X and .5X and 1X
[SBfinaly_BNd1XKir2, SBfinalt_BNd1XKir2, SBfinalcurrents_BNd1XKir2, timeImport_BNd1XKir2, yImport_BNd1XKir2, currentsImport_BNd1XKir2] = ShannonBers_rabbit(0.1, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd2XKir2, SBfinalt_BNd2XKir2, SBfinalcurrents_BNd2XKir2, timeImport_BNd2XKir2, yImport_BNd2XKir2, currentsImport_BNd2XKir2] = ShannonBers_rabbit(0.2, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd4XKir2, SBfinalt_BNd4XKir2, SBfinalcurrents_BNd4XKir2, timeImport_BNd4XKir2, yImport_BNd4XKir2, currentsImport_BNd4XKir2] = ShannonBers_rabbit(0.4, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd6XKir2, SBfinalt_BNd6XKir2, SBfinalcurrents_BNd6XKir2, timeImport_BNd6XKir2, yImport_BNd6XKir2, currentsImport_BNd6XKir2] = ShannonBers_rabbit(0.6, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd8XKir2, SBfinalt_BNd8XKir2, SBfinalcurrents_BNd8XKir2, timeImport_BNd8XKir2, yImport_BNd8XKir2, currentsImport_BNd8XKir2] = ShannonBers_rabbit(0.8, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BN1XKir2, SBfinalt_BN1XKir2, SBfinalcurrents_BN1XKir2, timeImport_BN1XKir2, yImport_BN1XKir2, currentsImport_BN1XKir2] = ShannonBers_rabbit(1., 0, numBeats, CycleLength, param, diffRateBool);




% 
% %Multiplier = 1.1
% param.GKir = 1.1;
% %Healthy
% [SBfinalyKir1d1, SBfinaltKir1d1, SBfinalcurrentsKir1d1, timeImportKir1d1, yImportKir1d1, currentsImportKir1d1] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);
% %BacNav .25X and .5X and 1X
% [SBfinaly_BNd1XKir1d1, SBfinalt_BNd1XKir1d1, SBfinalcurrents_BNd1XKir1d1, timeImport_BNd1XKir1d1, yImport_BNd1XKir1d1, currentsImport_BNd1XKir1d1] = ShannonBers_rabbit(0.1, 0, numBeats, CycleLength, param, diffRateBool);
% [SBfinaly_BNd2XKir1d1, SBfinalt_BNd2XKir1d1, SBfinalcurrents_BNd2XKir1d1, timeImport_BNd2XKir1d1, yImport_BNd2XKir1d1, currentsImport_BNd2XKir1d1] = ShannonBers_rabbit(0.2, 0, numBeats, CycleLength, param, diffRateBool);
% [SBfinaly_BNd4XKir1d1, SBfinalt_BNd4XKir1d1, SBfinalcurrents_BNd4XKir1d1, timeImport_BNd4XKir1d1, yImport_BNd4XKir1d1, currentsImport_BNd4XKir1d1] = ShannonBers_rabbit(0.4, 0, numBeats, CycleLength, param, diffRateBool);
% [SBfinaly_BNd6XKir1d1, SBfinalt_BNd6XKir1d1, SBfinalcurrents_BNd6XKir1d1, timeImport_BNd6XKir1d1, yImport_BNd6XKir1d1, currentsImport_BNd6XKir1d1] = ShannonBers_rabbit(0.6, 0, numBeats, CycleLength, param, diffRateBool);
% [SBfinaly_BNd8XKir1d1, SBfinalt_BNd8XKir1d1, SBfinalcurrents_BNd8XKir1d1, timeImport_BNd8XKir1d1, yImport_BNd8XKir1d1, currentsImport_BNd8XKir1d1] = ShannonBers_rabbit(0.8, 0, numBeats, CycleLength, param, diffRateBool);
% [SBfinaly_BN1XKir1d1, SBfinalt_BN1XKir1d1, SBfinalcurrents_BN1XKir1d1, timeImport_BN1XKir1d1, yImport_BN1XKir1d1, currentsImport_BN1XKir1d1] = ShannonBers_rabbit(1., 0, numBeats, CycleLength, param, diffRateBool);
% 



toc



% multiplierFactors = [.1, .2, .3, .4, .5, .6, .7, .8, .9, 1., 1.1];
multiplierFactors = [ .2, .4, .6, .8, 1., 1.2, 1.4, 1.6, 1.8, 2. 1.1 ];

peakCai_BN0X = [ max(SBfinalyKird2(:,38)), max(SBfinalyKird4(:,38)), max(SBfinalyKird6(:,38)), max(SBfinalyKird8(:,38)), max(SBfinalyKir1(:,38)),...
    max(SBfinalyKir1d2(:,38)), max(SBfinalyKir1d4(:,38)), max(SBfinalyKir1d6(:,38)), max(SBfinalyKir1d8(:,38)), max(SBfinalyKir2(:,38)), max(SBfinalyKir1d1(:,38)) ];

peakCai_BNd1X = [ max(SBfinaly_BNd1XKird2(:,38)), max(SBfinaly_BNd1XKird4(:,38)), max(SBfinaly_BNd1XKird6(:,38)), max(SBfinaly_BNd1XKird8(:,38)), max(SBfinaly_BNd1XKir1(:,38)),...
    max(SBfinaly_BNd1XKir1d2(:,38)), max(SBfinaly_BNd1XKir1d4(:,38)), max(SBfinaly_BNd1XKir1d6(:,38)), max(SBfinaly_BNd1XKir1d8(:,38)), max(SBfinaly_BNd1XKir2(:,38)), max(SBfinaly_BNd1XKir1d1(:,38)) ];

peakCai_BNd2X = [ max(SBfinaly_BNd2XKird2(:,38)), max(SBfinaly_BNd2XKird4(:,38)), max(SBfinaly_BNd2XKird6(:,38)), max(SBfinaly_BNd2XKird8(:,38)), max(SBfinaly_BNd2XKir1(:,38)),...
    max(SBfinaly_BNd2XKir1d2(:,38)), max(SBfinaly_BNd2XKir1d4(:,38)), max(SBfinaly_BNd2XKir1d6(:,38)), max(SBfinaly_BNd2XKir1d8(:,38)), max(SBfinaly_BNd2XKir2(:,38)), max(SBfinaly_BNd2XKir1d1(:,38)) ];

peakCai_BNd4X = [ max(SBfinaly_BNd4XKird2(:,38)), max(SBfinaly_BNd4XKird4(:,38)), max(SBfinaly_BNd4XKird6(:,38)), max(SBfinaly_BNd4XKird8(:,38)), max(SBfinaly_BNd4XKir1(:,38)),...
    max(SBfinaly_BNd4XKir1d2(:,38)), max(SBfinaly_BNd4XKir1d4(:,38)), max(SBfinaly_BNd4XKir1d6(:,38)), max(SBfinaly_BNd4XKir1d8(:,38)), max(SBfinaly_BNd4XKir2(:,38)), max(SBfinaly_BNd4XKir1d1(:,38)) ];

peakCai_BNd6X = [ max(SBfinaly_BNd6XKird2(:,38)), max(SBfinaly_BNd6XKird4(:,38)), max(SBfinaly_BNd6XKird6(:,38)), max(SBfinaly_BNd6XKird8(:,38)), max(SBfinaly_BNd6XKir1(:,38)),...
    max(SBfinaly_BNd6XKir1d2(:,38)), max(SBfinaly_BNd6XKir1d4(:,38)), max(SBfinaly_BNd6XKir1d6(:,38)), max(SBfinaly_BNd6XKir1d8(:,38)), max(SBfinaly_BNd6XKir2(:,38)), max(SBfinaly_BNd6XKir1d1(:,38)) ];

peakCai_BNd8X = [ max(SBfinaly_BNd8XKird2(:,38)), max(SBfinaly_BNd8XKird4(:,38)), max(SBfinaly_BNd8XKird6(:,38)), max(SBfinaly_BNd8XKird8(:,38)), max(SBfinaly_BNd8XKir1(:,38)),...
    max(SBfinaly_BNd8XKir1d2(:,38)), max(SBfinaly_BNd8XKir1d4(:,38)), max(SBfinaly_BNd8XKir1d6(:,38)), max(SBfinaly_BNd8XKir1d8(:,38)), max(SBfinaly_BNd8XKir2(:,38)), max(SBfinaly_BNd8XKir1d1(:,38)) ];

peakCai_BN1X = [ max(SBfinaly_BN1XKird2(:,38)), max(SBfinaly_BN1XKird4(:,38)), max(SBfinaly_BN1XKird6(:,38)), max(SBfinaly_BN1XKird8(:,38)), max(SBfinaly_BN1XKir1(:,38)),...
    max(SBfinaly_BN1XKir1d2(:,38)), max(SBfinaly_BN1XKir1d4(:,38)), max(SBfinaly_BN1XKir1d6(:,38)), max(SBfinaly_BN1XKir1d8(:,38)), max(SBfinaly_BN1XKir2(:,38)), max(SBfinaly_BN1XKir1d1(:,38)) ];



Upstroke_BN0X = [ max(diff(SBfinalyKird2(:,39))./diff(SBfinaltKird2)), max(diff(SBfinalyKird4(:,39))./diff(SBfinaltKird4)), max(diff(SBfinalyKird6(:,39))./diff(SBfinaltKird6)), max(diff(SBfinalyKird8(:,39))./diff(SBfinaltKird8)), max(diff(SBfinalyKir1(:,39))./diff(SBfinaltKir1)),...
    max(diff(SBfinalyKir1d2(:,39))./diff(SBfinaltKir1d2)), max(diff(SBfinalyKir1d4(:,39))./diff(SBfinaltKir1d4)), max(diff(SBfinalyKir1d6(:,39))./diff(SBfinaltKir1d6)), max(diff(SBfinalyKir1d8(:,39))./diff(SBfinaltKir1d8)), max(diff(SBfinalyKir2(:,39))./diff(SBfinaltKir2)), max(diff(SBfinalyKir1d1(:,39))./diff(SBfinaltKir1d1)) ];

Upstroke_BNd1X = [ max(diff(SBfinaly_BNd1XKird2(:,39))./diff(SBfinalt_BNd1XKird2)), max(diff(SBfinaly_BNd1XKird4(:,39))./diff(SBfinalt_BNd1XKird4)), max(diff(SBfinaly_BNd1XKird6(:,39))./diff(SBfinalt_BNd1XKird6)), max(diff(SBfinaly_BNd1XKird8(:,39))./diff(SBfinalt_BNd1XKird8)), max(diff(SBfinaly_BNd1XKir1(:,39))./diff(SBfinalt_BNd1XKir1)),...
    max(diff(SBfinaly_BNd1XKir1d2(:,39))./diff(SBfinalt_BNd1XKir1d2)), max(diff(SBfinaly_BNd1XKir1d4(:,39))./diff(SBfinalt_BNd1XKir1d4)), max(diff(SBfinaly_BNd1XKir1d6(:,39))./diff(SBfinalt_BNd1XKir1d6)), max(diff(SBfinaly_BNd1XKir1d8(:,39))./diff(SBfinalt_BNd1XKir1d8)), max(diff(SBfinaly_BNd1XKir2(:,39))./diff(SBfinalt_BNd1XKir2)), max(diff(SBfinaly_BNd1XKir1d1(:,39))./diff(SBfinalt_BNd1XKir1d1)) ];

Upstroke_BNd2X = [ max(diff(SBfinaly_BNd2XKird2(:,39))./diff(SBfinalt_BNd2XKird2)), max(diff(SBfinaly_BNd2XKird4(:,39))./diff(SBfinalt_BNd2XKird4)), max(diff(SBfinaly_BNd2XKird6(:,39))./diff(SBfinalt_BNd2XKird6)), max(diff(SBfinaly_BNd2XKird8(:,39))./diff(SBfinalt_BNd2XKird8)), max(diff(SBfinaly_BNd2XKir1(:,39))./diff(SBfinalt_BNd2XKir1)),...
    max(diff(SBfinaly_BNd2XKir1d2(:,39))./diff(SBfinalt_BNd2XKir1d2)), max(diff(SBfinaly_BNd2XKir1d4(:,39))./diff(SBfinalt_BNd2XKir1d4)), max(diff(SBfinaly_BNd2XKir1d6(:,39))./diff(SBfinalt_BNd2XKir1d6)), max(diff(SBfinaly_BNd2XKir1d8(:,39))./diff(SBfinalt_BNd2XKir1d8)), max(diff(SBfinaly_BNd2XKir2(:,39))./diff(SBfinalt_BNd2XKir2)), max(diff(SBfinaly_BNd2XKir1d1(:,39))./diff(SBfinalt_BNd2XKir1d1)) ];

Upstroke_BNd4X = [ max(diff(SBfinaly_BNd4XKird2(:,39))./diff(SBfinalt_BNd4XKird2)), max(diff(SBfinaly_BNd4XKird4(:,39))./diff(SBfinalt_BNd4XKird4)), max(diff(SBfinaly_BNd4XKird6(:,39))./diff(SBfinalt_BNd4XKird6)), max(diff(SBfinaly_BNd4XKird8(:,39))./diff(SBfinalt_BNd4XKird8)), max(diff(SBfinaly_BNd4XKir1(:,39))./diff(SBfinalt_BNd4XKir1)),...
    max(diff(SBfinaly_BNd4XKir1d2(:,39))./diff(SBfinalt_BNd4XKir1d2)), max(diff(SBfinaly_BNd4XKir1d4(:,39))./diff(SBfinalt_BNd4XKir1d4)), max(diff(SBfinaly_BNd4XKir1d6(:,39))./diff(SBfinalt_BNd4XKir1d6)), max(diff(SBfinaly_BNd4XKir1d8(:,39))./diff(SBfinalt_BNd4XKir1d8)), max(diff(SBfinaly_BNd4XKir2(:,39))./diff(SBfinalt_BNd4XKir2)), max(diff(SBfinaly_BNd4XKir1d1(:,39))./diff(SBfinalt_BNd4XKir1d1)) ];

Upstroke_BNd6X = [ max(diff(SBfinaly_BNd6XKird2(:,39))./diff(SBfinalt_BNd6XKird2)), max(diff(SBfinaly_BNd6XKird4(:,39))./diff(SBfinalt_BNd6XKird4)), max(diff(SBfinaly_BNd6XKird6(:,39))./diff(SBfinalt_BNd6XKird6)), max(diff(SBfinaly_BNd6XKird8(:,39))./diff(SBfinalt_BNd6XKird8)), max(diff(SBfinaly_BNd6XKir1(:,39))./diff(SBfinalt_BNd6XKir1)),...
    max(diff(SBfinaly_BNd6XKir1d2(:,39))./diff(SBfinalt_BNd6XKir1d2)), max(diff(SBfinaly_BNd6XKir1d4(:,39))./diff(SBfinalt_BNd6XKir1d4)), max(diff(SBfinaly_BNd6XKir1d6(:,39))./diff(SBfinalt_BNd6XKir1d6)), max(diff(SBfinaly_BNd6XKir1d8(:,39))./diff(SBfinalt_BNd6XKir1d8)), max(diff(SBfinaly_BNd6XKir2(:,39))./diff(SBfinalt_BNd6XKir2)), max(diff(SBfinaly_BNd6XKir1d1(:,39))./diff(SBfinalt_BNd6XKir1d1)) ];

Upstroke_BNd8X = [ max(diff(SBfinaly_BNd8XKird2(:,39))./diff(SBfinalt_BNd8XKird2)), max(diff(SBfinaly_BNd8XKird4(:,39))./diff(SBfinalt_BNd8XKird4)), max(diff(SBfinaly_BNd8XKird6(:,39))./diff(SBfinalt_BNd8XKird6)), max(diff(SBfinaly_BNd8XKird8(:,39))./diff(SBfinalt_BNd8XKird8)), max(diff(SBfinaly_BNd8XKir1(:,39))./diff(SBfinalt_BNd8XKir1)),...
    max(diff(SBfinaly_BNd8XKir1d2(:,39))./diff(SBfinalt_BNd8XKir1d2)), max(diff(SBfinaly_BNd8XKir1d4(:,39))./diff(SBfinalt_BNd8XKir1d4)), max(diff(SBfinaly_BNd8XKir1d6(:,39))./diff(SBfinalt_BNd8XKir1d6)), max(diff(SBfinaly_BNd8XKir1d8(:,39))./diff(SBfinalt_BNd8XKir1d8)), max(diff(SBfinaly_BNd8XKir2(:,39))./diff(SBfinalt_BNd8XKir2)), max(diff(SBfinaly_BNd8XKir1d1(:,39))./diff(SBfinalt_BNd8XKir1d1)) ];

Upstroke_BN1X = [ max(diff(SBfinaly_BN1XKird2(:,39))./diff(SBfinalt_BN1XKird2)), max(diff(SBfinaly_BN1XKird4(:,39))./diff(SBfinalt_BN1XKird4)), max(diff(SBfinaly_BN1XKird6(:,39))./diff(SBfinalt_BN1XKird6)), max(diff(SBfinaly_BN1XKird8(:,39))./diff(SBfinalt_BN1XKird8)), max(diff(SBfinaly_BN1XKir1(:,39))./diff(SBfinalt_BN1XKir1)),...
    max(diff(SBfinaly_BN1XKir1d2(:,39))./diff(SBfinalt_BN1XKir1d2)), max(diff(SBfinaly_BN1XKir1d4(:,39))./diff(SBfinalt_BN1XKir1d4)), max(diff(SBfinaly_BN1XKir1d6(:,39))./diff(SBfinalt_BN1XKir1d6)), max(diff(SBfinaly_BN1XKir1d8(:,39))./diff(SBfinalt_BN1XKir1d8)), max(diff(SBfinaly_BN1XKir2(:,39))./diff(SBfinalt_BN1XKir2)), max(diff(SBfinaly_BN1XKir1d1(:,39))./diff(SBfinalt_BN1XKir1d1)) ];




APA_BN0X = [ max(SBfinalyKird2(:,39)) - min(SBfinalyKird2(:,39)), max(SBfinalyKird4(:,39)) - min(SBfinalyKird4(:,39)), max(SBfinalyKird6(:,39)) - min(SBfinalyKird6(:,39)), max(SBfinalyKird8(:,39)) - min(SBfinalyKird8(:,39)), max(SBfinalyKir1(:,39)) - min(SBfinalyKir1(:,39)),...
    max(SBfinalyKir1d2(:,39)) - min(SBfinalyKir1d2(:,39)), max(SBfinalyKir1d4(:,39)) - min(SBfinalyKir1d4(:,39)), max(SBfinalyKir1d6(:,39)) - min(SBfinalyKir1d6(:,39)), max(SBfinalyKir1d8(:,39)) - min(SBfinalyKir1d8(:,39)), max(SBfinalyKir2(:,39)) - min(SBfinalyKir2(:,39)), max(SBfinalyKir1d1(:,39)) - min(SBfinalyKir1d1(:,39)) ];

APA_BNd1X = [ max(SBfinaly_BNd1XKird2(:,39)) - min(SBfinaly_BNd1XKird2(:,39)), max(SBfinaly_BNd1XKird4(:,39)) - min(SBfinaly_BNd1XKird4(:,39)), max(SBfinaly_BNd1XKird6(:,39)) - min(SBfinaly_BNd1XKird6(:,39)), max(SBfinaly_BNd1XKird8(:,39)) - min(SBfinaly_BNd1XKird8(:,39)), max(SBfinaly_BNd1XKir1(:,39)) - min(SBfinaly_BNd1XKir1(:,39)),...
    max(SBfinaly_BNd1XKir1d2(:,39)) - min(SBfinaly_BNd1XKir1d2(:,39)), max(SBfinaly_BNd1XKir1d4(:,39)) - min(SBfinaly_BNd1XKir1d4(:,39)), max(SBfinaly_BNd1XKir1d6(:,39)) - min(SBfinaly_BNd1XKir1d6(:,39)), max(SBfinaly_BNd1XKir1d8(:,39)) - min(SBfinaly_BNd1XKir1d8(:,39)), max(SBfinaly_BNd1XKir2(:,39)) - min(SBfinaly_BNd1XKir2(:,39)), max(SBfinaly_BNd1XKir1d1(:,39)) - min(SBfinaly_BNd1XKir1d1(:,39)) ];

APA_BNd2X = [ max(SBfinaly_BNd2XKird2(:,39)) - min(SBfinaly_BNd2XKird2(:,39)), max(SBfinaly_BNd2XKird4(:,39)) - min(SBfinaly_BNd2XKird4(:,39)), max(SBfinaly_BNd2XKird6(:,39)) - min(SBfinaly_BNd2XKird6(:,39)), max(SBfinaly_BNd2XKird8(:,39)) - min(SBfinaly_BNd2XKird8(:,39)), max(SBfinaly_BNd2XKir1(:,39)) - min(SBfinaly_BNd2XKir1(:,39)),...
    max(SBfinaly_BNd2XKir1d2(:,39)) - min(SBfinaly_BNd2XKir1d2(:,39)), max(SBfinaly_BNd2XKir1d4(:,39)) - min(SBfinaly_BNd2XKir1d4(:,39)), max(SBfinaly_BNd2XKir1d6(:,39)) - min(SBfinaly_BNd2XKir1d6(:,39)), max(SBfinaly_BNd2XKir1d8(:,39)) - min(SBfinaly_BNd2XKir1d8(:,39)), max(SBfinaly_BNd2XKir2(:,39)) - min(SBfinaly_BNd2XKir2(:,39)), max(SBfinaly_BNd2XKir1d1(:,39)) - min(SBfinaly_BNd2XKir1d1(:,39)) ];

APA_BNd4X = [ max(SBfinaly_BNd4XKird2(:,39)) - min(SBfinaly_BNd4XKird2(:,39)), max(SBfinaly_BNd4XKird4(:,39)) - min(SBfinaly_BNd4XKird4(:,39)), max(SBfinaly_BNd4XKird6(:,39)) - min(SBfinaly_BNd4XKird6(:,39)), max(SBfinaly_BNd4XKird8(:,39)) - min(SBfinaly_BNd4XKird8(:,39)), max(SBfinaly_BNd4XKir1(:,39)) - min(SBfinaly_BNd4XKir1(:,39)),...
    max(SBfinaly_BNd4XKir1d2(:,39)) - min(SBfinaly_BNd4XKir1d2(:,39)), max(SBfinaly_BNd4XKir1d4(:,39)) - min(SBfinaly_BNd4XKir1d4(:,39)), max(SBfinaly_BNd4XKir1d6(:,39)) - min(SBfinaly_BNd4XKir1d6(:,39)), max(SBfinaly_BNd4XKir1d8(:,39)) - min(SBfinaly_BNd4XKir1d8(:,39)), max(SBfinaly_BNd4XKir2(:,39)) - min(SBfinaly_BNd4XKir2(:,39)), max(SBfinaly_BNd4XKir1d1(:,39)) - min(SBfinaly_BNd4XKir1d1(:,39)) ];

APA_BNd6X = [ max(SBfinaly_BNd6XKird2(:,39)) - min(SBfinaly_BNd6XKird2(:,39)), max(SBfinaly_BNd6XKird4(:,39)) - min(SBfinaly_BNd6XKird4(:,39)), max(SBfinaly_BNd6XKird6(:,39)) - min(SBfinaly_BNd6XKird6(:,39)), max(SBfinaly_BNd6XKird8(:,39)) - min(SBfinaly_BNd6XKird8(:,39)), max(SBfinaly_BNd6XKir1(:,39)) - min(SBfinaly_BNd6XKir1(:,39)),...
    max(SBfinaly_BNd6XKir1d2(:,39)) - min(SBfinaly_BNd6XKir1d2(:,39)), max(SBfinaly_BNd6XKir1d4(:,39)) - min(SBfinaly_BNd6XKir1d4(:,39)), max(SBfinaly_BNd6XKir1d6(:,39)) - min(SBfinaly_BNd6XKir1d6(:,39)), max(SBfinaly_BNd6XKir1d8(:,39)) - min(SBfinaly_BNd6XKir1d8(:,39)), max(SBfinaly_BNd6XKir2(:,39)) - min(SBfinaly_BNd6XKir2(:,39)), max(SBfinaly_BNd6XKir1d1(:,39)) - min(SBfinaly_BNd6XKir1d1(:,39)) ];

APA_BNd8X = [ max(SBfinaly_BNd8XKird2(:,39)) - min(SBfinaly_BNd8XKird2(:,39)), max(SBfinaly_BNd8XKird4(:,39)) - min(SBfinaly_BNd8XKird4(:,39)), max(SBfinaly_BNd8XKird6(:,39)) - min(SBfinaly_BNd8XKird6(:,39)), max(SBfinaly_BNd8XKird8(:,39)) - min(SBfinaly_BNd8XKird8(:,39)), max(SBfinaly_BNd8XKir1(:,39)) - min(SBfinaly_BNd8XKir1(:,39)),...
    max(SBfinaly_BNd8XKir1d2(:,39)) - min(SBfinaly_BNd8XKir1d2(:,39)), max(SBfinaly_BNd8XKir1d4(:,39)) - min(SBfinaly_BNd8XKir1d4(:,39)), max(SBfinaly_BNd8XKir1d6(:,39)) - min(SBfinaly_BNd8XKir1d6(:,39)), max(SBfinaly_BNd8XKir1d8(:,39)) - min(SBfinaly_BNd8XKir1d8(:,39)), max(SBfinaly_BNd8XKir2(:,39)) - min(SBfinaly_BNd8XKir2(:,39)), max(SBfinaly_BNd8XKir1d1(:,39)) - min(SBfinaly_BNd8XKir1d1(:,39)) ];

APA_BN1X = [ max(SBfinaly_BN1XKird2(:,39)) - min(SBfinaly_BN1XKird2(:,39)), max(SBfinaly_BN1XKird4(:,39)) - min(SBfinaly_BN1XKird4(:,39)), max(SBfinaly_BN1XKird6(:,39)) - min(SBfinaly_BN1XKird6(:,39)), max(SBfinaly_BN1XKird8(:,39)) - min(SBfinaly_BN1XKird8(:,39)), max(SBfinaly_BN1XKir1(:,39)) - min(SBfinaly_BN1XKir1(:,39)),...
    max(SBfinaly_BN1XKir1d2(:,39)) - min(SBfinaly_BN1XKir1d2(:,39)), max(SBfinaly_BN1XKir1d4(:,39)) - min(SBfinaly_BN1XKir1d4(:,39)), max(SBfinaly_BN1XKir1d6(:,39)) - min(SBfinaly_BN1XKir1d6(:,39)), max(SBfinaly_BN1XKir1d8(:,39)) - min(SBfinaly_BN1XKir1d8(:,39)), max(SBfinaly_BN1XKir2(:,39)) - min(SBfinaly_BN1XKir2(:,39)), max(SBfinaly_BN1XKir1d1(:,39)) - min(SBfinaly_BN1XKir1d1(:,39)) ];




MatrixToTable_peakCai  = [ multiplierFactors', peakCai_BN0X', peakCai_BNd1X', peakCai_BNd2X', peakCai_BNd4X', peakCai_BNd6X', peakCai_BNd8X', peakCai_BN1X' ]; 
MatrixToTable_Upstroke = [ multiplierFactors', Upstroke_BN0X', Upstroke_BNd1X', Upstroke_BNd2X', Upstroke_BNd4X', Upstroke_BNd6X', Upstroke_BNd8X', Upstroke_BN1X' ];
MatrixToTable_APA = [ multiplierFactors', APA_BN0X', APA_BNd1X', APA_BNd2X', APA_BNd4X', APA_BNd6X', APA_BNd8X', APA_BN1X' ];

%%

%FULL
figure
plot( [.25, .5, .75, 1., 1.25, 1.5], [ max(SBfinalyKird25(:,38)), max(SBfinalyKird5(:,38)), max(SBfinalyKird75(:,38)), max(SBfinalyKir1(:,38)), max(SBfinalyKir1d25(:,38)), max(SBfinalyKir1d5(:,38)) ], 'DisplayName', '0X BacNav' )
hold on
plot( [.25, .5, .75, 1., 1.25, 1.5], [ max(SBfinaly_BNd1XKird25(:,38)), max(SBfinaly_BNd1XKird5(:,38)), max(SBfinaly_BNd1XKird75(:,38)), max(SBfinaly_BNd1XKir1(:,38)), max(SBfinaly_BNd1XKir1d25(:,38)), max(SBfinaly_BNd1XKir1d5(:,38)) ], 'DisplayName', '.1X BacNav' )
plot( [.25, .5, .75, 1., 1.25, 1.5], [ max(SBfinaly_BNd2XKird25(:,38)), max(SBfinaly_BNd2XKird5(:,38)), max(SBfinaly_BNd2XKird75(:,38)), max(SBfinaly_BNd2XKir1(:,38)), max(SBfinaly_BNd2XKir1d25(:,38)), max(SBfinaly_BNd2XKir1d5(:,38)) ], 'DisplayName', '.2X BacNav' )
plot( [.25, .5, .75, 1., 1.25, 1.5], [ max(SBfinaly_BNd4XKird25(:,38)), max(SBfinaly_BNd4XKird5(:,38)), max(SBfinaly_BNd4XKird75(:,38)), max(SBfinaly_BNd4XKir1(:,38)), max(SBfinaly_BNd4XKir1d25(:,38)), max(SBfinaly_BNd4XKir1d5(:,38)) ], 'DisplayName', '.4X BacNav' )
plot( [.25, .5, .75, 1., 1.25, 1.5], [ max(SBfinaly_BNd6XKird25(:,38)), max(SBfinaly_BNd6XKird5(:,38)), max(SBfinaly_BNd6XKird75(:,38)), max(SBfinaly_BNd6XKir1(:,38)), max(SBfinaly_BNd6XKir1d25(:,38)), max(SBfinaly_BNd6XKir1d5(:,38)) ], 'DisplayName', '.6X BacNav' )
plot( [.25, .5, .75, 1., 1.25, 1.5], [ max(SBfinaly_BNd8XKird25(:,38)), max(SBfinaly_BNd8XKird5(:,38)), max(SBfinaly_BNd8XKird75(:,38)), max(SBfinaly_BNd8XKir1(:,38)), max(SBfinaly_BNd8XKir1d25(:,38)), max(SBfinaly_BNd8XKir1d5(:,38)) ], 'DisplayName', '.8X BacNav' )
plot( [.25, .5, .75, 1., 1.25, 1.5], [ max(SBfinaly_BN1XKird25(:,38)), max(SBfinaly_BN1XKird5(:,38)), max(SBfinaly_BN1XKird75(:,38)), max(SBfinaly_BN1XKir1(:,38)), max(SBfinaly_BN1XKir1d25(:,38)), max(SBfinaly_BN1XKir1d5(:,38)) ], 'DisplayName', '1X BacNav' )
legend()
xlabel('Kir2.1 conductance multiplier')
ylabel('Peak Ca_i (uMol)')
title('Peak Ca_i vs. conductance factor for Kir2.1 channel')




figure
plot( [.25, .5, .75, 1., 1.25, 1.5], [ max(diff(SBfinalyKird25(:,39))./diff(SBfinaltKird25)), max(diff(SBfinalyKird5(:,39))./diff(SBfinaltKird5)), max(diff(SBfinalyKird75(:,39))./diff(SBfinaltKird75)), max(diff(SBfinalyKir1(:,39))./diff(SBfinaltKir1)), max(diff(SBfinalyKir1d25(:,39))./diff(SBfinaltKir1d25)), max(diff(SBfinalyKir1d5(:,39))./diff(SBfinaltKir1d5)) ], 'DisplayName', '0X BacNav' )
hold on
plot( [.25, .5, .75, 1., 1.25, 1.5], [ max(diff(SBfinaly_BNd1XKird25(:,39))./diff(SBfinalt_BNd1XKird25)), max(diff(SBfinaly_BNd1XKird5(:,39))./diff(SBfinalt_BNd1XKird5)), max(diff(SBfinaly_BNd1XKird75(:,39))./diff(SBfinalt_BNd1XKird75)), max(diff(SBfinaly_BNd1XKir1(:,39))./diff(SBfinalt_BNd1XKir1)), max(diff(SBfinaly_BNd1XKir1d25(:,39))./diff(SBfinalt_BNd1XKir1d25)), max(diff(SBfinaly_BNd1XKir1d5(:,39))./diff(SBfinalt_BNd1XKir1d5)) ], 'DisplayName', '.1X BacNav' )
plot( [.25, .5, .75, 1., 1.25, 1.5], [ max(diff(SBfinaly_BNd2XKird25(:,39))./diff(SBfinalt_BNd2XKird25)), max(diff(SBfinaly_BNd2XKird5(:,39))./diff(SBfinalt_BNd2XKird5)), max(diff(SBfinaly_BNd2XKird75(:,39))./diff(SBfinalt_BNd2XKird75)), max(diff(SBfinaly_BNd2XKir1(:,39))./diff(SBfinalt_BNd2XKir1)), max(diff(SBfinaly_BNd2XKir1d25(:,39))./diff(SBfinalt_BNd2XKir1d25)), max(diff(SBfinaly_BNd2XKir1d5(:,39))./diff(SBfinalt_BNd2XKir1d5)) ], 'DisplayName', '.2X BacNav' )
plot( [.25, .5, .75, 1., 1.25, 1.5], [ max(diff(SBfinaly_BNd4XKird25(:,39))./diff(SBfinalt_BNd4XKird25)), max(diff(SBfinaly_BNd4XKird5(:,39))./diff(SBfinalt_BNd4XKird5)), max(diff(SBfinaly_BNd4XKird75(:,39))./diff(SBfinalt_BNd4XKird75)), max(diff(SBfinaly_BNd4XKir1(:,39))./diff(SBfinalt_BNd4XKir1)), max(diff(SBfinaly_BNd4XKir1d25(:,39))./diff(SBfinalt_BNd4XKir1d25)), max(diff(SBfinaly_BNd4XKir1d5(:,39))./diff(SBfinalt_BNd4XKir1d5)) ], 'DisplayName', '.4X BacNav' )
plot( [.25, .5, .75, 1., 1.25, 1.5], [ max(diff(SBfinaly_BNd6XKird25(:,39))./diff(SBfinalt_BNd6XKird25)), max(diff(SBfinaly_BNd6XKird5(:,39))./diff(SBfinalt_BNd6XKird5)), max(diff(SBfinaly_BNd6XKird75(:,39))./diff(SBfinalt_BNd6XKird75)), max(diff(SBfinaly_BNd6XKir1(:,39))./diff(SBfinalt_BNd6XKir1)), max(diff(SBfinaly_BNd6XKir1d25(:,39))./diff(SBfinalt_BNd6XKir1d25)), max(diff(SBfinaly_BNd6XKir1d5(:,39))./diff(SBfinalt_BNd6XKir1d5)) ], 'DisplayName', '.6X BacNav' )
plot( [.25, .5, .75, 1., 1.25, 1.5], [ max(diff(SBfinaly_BNd8XKird25(:,39))./diff(SBfinalt_BNd8XKird25)), max(diff(SBfinaly_BNd8XKird5(:,39))./diff(SBfinalt_BNd8XKird5)), max(diff(SBfinaly_BNd8XKird75(:,39))./diff(SBfinalt_BNd8XKird75)), max(diff(SBfinaly_BNd8XKir1(:,39))./diff(SBfinalt_BNd8XKir1)), max(diff(SBfinaly_BNd8XKir1d25(:,39))./diff(SBfinalt_BNd8XKir1d25)), max(diff(SBfinaly_BNd8XKir1d5(:,39))./diff(SBfinalt_BNd8XKir1d5)) ], 'DisplayName', '.8X BacNav' )
plot( [.25, .5, .75, 1., 1.25, 1.5], [ max(diff(SBfinaly_BN1XKird25(:,39))./diff(SBfinalt_BN1XKird25)), max(diff(SBfinaly_BN1XKird5(:,39))./diff(SBfinalt_BN1XKird5)), max(diff(SBfinaly_BN1XKird75(:,39))./diff(SBfinalt_BN1XKird75)), max(diff(SBfinaly_BN1XKir1(:,39))./diff(SBfinalt_BN1XKir1)), max(diff(SBfinaly_BN1XKir1d25(:,39))./diff(SBfinalt_BN1XKir1d25)), max(diff(SBfinaly_BN1XKir1d5(:,39))./diff(SBfinalt_BN1XKir1d5)) ], 'DisplayName', '1X BacNav' )
legend()
xlabel('Kir2.1 conductance multiplier')
ylabel('AP Upstroke (mV/ms)')
title('AP Upstroke vs. conductance factor for Kir2.1 channel')



%CROPPED

figure
plot( [.25, .5, .75, 1.], [ max(SBfinalyKird25(:,38)), max(SBfinalyKird5(:,38)), max(SBfinalyKird75(:,38)), max(SBfinalyKir1(:,38)) ], 'DisplayName', '0X BacNav' )
hold on
plot( [.25, .5, .75, 1.], [ max(SBfinaly_BNd1XKird25(:,38)), max(SBfinaly_BNd1XKird5(:,38)), max(SBfinaly_BNd1XKird75(:,38)), max(SBfinaly_BNd1XKir1(:,38)) ], 'DisplayName', '.1X BacNav' )
plot( [.25, .5, .75, 1.], [ max(SBfinaly_BNd2XKird25(:,38)), max(SBfinaly_BNd2XKird5(:,38)), max(SBfinaly_BNd2XKird75(:,38)), max(SBfinaly_BNd2XKir1(:,38)) ], 'DisplayName', '.2X BacNav' )
plot( [.25, .5, .75, 1.], [ max(SBfinaly_BNd4XKird25(:,38)), max(SBfinaly_BNd4XKird5(:,38)), max(SBfinaly_BNd4XKird75(:,38)), max(SBfinaly_BNd4XKir1(:,38)) ], 'DisplayName', '.4X BacNav' )
plot( [.25, .5, .75, 1.], [ max(SBfinaly_BNd6XKird25(:,38)), max(SBfinaly_BNd6XKird5(:,38)), max(SBfinaly_BNd6XKird75(:,38)), max(SBfinaly_BNd6XKir1(:,38)) ], 'DisplayName', '.6X BacNav' )
plot( [.25, .5, .75, 1.], [ max(SBfinaly_BNd8XKird25(:,38)), max(SBfinaly_BNd8XKird5(:,38)), max(SBfinaly_BNd8XKird75(:,38)), max(SBfinaly_BNd8XKir1(:,38)) ], 'DisplayName', '.8X BacNav' )
plot( [.25, .5, .75, 1.], [ max(SBfinaly_BN1XKird25(:,38)), max(SBfinaly_BN1XKird5(:,38)), max(SBfinaly_BN1XKird75(:,38)), max(SBfinaly_BN1XKir1(:,38)) ], 'DisplayName', '1X BacNav' )
legend()
xlabel('Kir2.1 conductance multiplier')
ylabel('Peak Ca_i (uMol)')
title('Peak Ca_i vs. conductance factor for Kir2.1 channel')




figure
plot( [.25, .5, .75, 1.], [ max(diff(SBfinalyKird25(:,39))./diff(SBfinaltKird25)), max(diff(SBfinalyKird5(:,39))./diff(SBfinaltKird5)), max(diff(SBfinalyKird75(:,39))./diff(SBfinaltKird75)), max(diff(SBfinalyKir1(:,39))./diff(SBfinaltKir1)) ], 'DisplayName', '0X BacNav' )
hold on
plot( [.25, .5, .75, 1.], [ max(diff(SBfinaly_BNd1XKird25(:,39))./diff(SBfinalt_BNd1XKird25)), max(diff(SBfinaly_BNd1XKird5(:,39))./diff(SBfinalt_BNd1XKird5)), max(diff(SBfinaly_BNd1XKird75(:,39))./diff(SBfinalt_BNd1XKird75)), max(diff(SBfinaly_BNd1XKir1(:,39))./diff(SBfinalt_BNd1XKir1)) ], 'DisplayName', '.1X BacNav' )
plot( [.25, .5, .75, 1.], [ max(diff(SBfinaly_BNd2XKird25(:,39))./diff(SBfinalt_BNd2XKird25)), max(diff(SBfinaly_BNd2XKird5(:,39))./diff(SBfinalt_BNd2XKird5)), max(diff(SBfinaly_BNd2XKird75(:,39))./diff(SBfinalt_BNd2XKird75)), max(diff(SBfinaly_BNd2XKir1(:,39))./diff(SBfinalt_BNd2XKir1)) ], 'DisplayName', '.2X BacNav' )
plot( [.25, .5, .75, 1.], [ max(diff(SBfinaly_BNd4XKird25(:,39))./diff(SBfinalt_BNd4XKird25)), max(diff(SBfinaly_BNd4XKird5(:,39))./diff(SBfinalt_BNd4XKird5)), max(diff(SBfinaly_BNd4XKird75(:,39))./diff(SBfinalt_BNd4XKird75)), max(diff(SBfinaly_BNd4XKir1(:,39))./diff(SBfinalt_BNd4XKir1)) ], 'DisplayName', '.4X BacNav' )
plot( [.25, .5, .75, 1.], [ max(diff(SBfinaly_BNd6XKird25(:,39))./diff(SBfinalt_BNd6XKird25)), max(diff(SBfinaly_BNd6XKird5(:,39))./diff(SBfinalt_BNd6XKird5)), max(diff(SBfinaly_BNd6XKird75(:,39))./diff(SBfinalt_BNd6XKird75)), max(diff(SBfinaly_BNd6XKir1(:,39))./diff(SBfinalt_BNd6XKir1)) ], 'DisplayName', '.6X BacNav' )
plot( [.25, .5, .75, 1.], [ max(diff(SBfinaly_BNd8XKird25(:,39))./diff(SBfinalt_BNd8XKird25)), max(diff(SBfinaly_BNd8XKird5(:,39))./diff(SBfinalt_BNd8XKird5)), max(diff(SBfinaly_BNd8XKird75(:,39))./diff(SBfinalt_BNd8XKird75)), max(diff(SBfinaly_BNd8XKir1(:,39))./diff(SBfinalt_BNd8XKir1)) ], 'DisplayName', '.8X BacNav' )
plot( [.25, .5, .75, 1.], [ max(diff(SBfinaly_BN1XKird25(:,39))./diff(SBfinalt_BN1XKird25)), max(diff(SBfinaly_BN1XKird5(:,39))./diff(SBfinalt_BN1XKird5)), max(diff(SBfinaly_BN1XKird75(:,39))./diff(SBfinalt_BN1XKird75)), max(diff(SBfinaly_BN1XKir1(:,39))./diff(SBfinalt_BN1XKir1)) ], 'DisplayName', '1X BacNav' )
legend()
xlabel('Kir2.1 conductance multiplier')
ylabel('AP Upstroke (mV/ms)')
title('AP Upstroke vs. conductance factor for Kir2.1 channel')









%%


tic



%CAFFEINE EFFECT
param.caffeineBool = false;

param.protocol = 'pace1';%'vclamp';
param.ratePace = 1.e-3;%.001 is 1Hz since 1/.001=1000 so every 1000ms

param.GlycosidesBool = false;
numBeats = 1001;
CycleLength = 1000;
diffRateBool = false;
%Multiplier = 1.
param.GKir = 1.;
param.DADbool = false;



%Multiplier = .05
param.NCXMult = .05;

%Healthy
[SBfinalyNCXMult_d05, SBfinaltNCXMult_d05, SBfinalcurrentsNCXMult_d05, timeImportNCXMult_d05, yImportNCXMult_d05, currentsImportNCXMult_d05] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);
%BacNav .25X and .5X and 1X
[SBfinaly_BNd1XNCXMult_d05, SBfinalt_BNd1XNCXMult_d05, SBfinalcurrents_BNd1XNCXMult_d05, timeImport_BNd1XNCXMult_d05, yImport_BNd1XNCXMult_d05, currentsImport_BNd1XNCXMult_d05] = ShannonBers_rabbit(0.1, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd2XNCXMult_d05, SBfinalt_BNd2XNCXMult_d05, SBfinalcurrents_BNd2XNCXMult_d05, timeImport_BNd2XNCXMult_d05, yImport_BNd2XNCXMult_d05, currentsImport_BNd2XNCXMult_d05] = ShannonBers_rabbit(0.2, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd4XNCXMult_d05, SBfinalt_BNd4XNCXMult_d05, SBfinalcurrents_BNd4XNCXMult_d05, timeImport_BNd4XNCXMult_d05, yImport_BNd4XNCXMult_d05, currentsImport_BNd4XNCXMult_d05] = ShannonBers_rabbit(0.4, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd6XNCXMult_d05, SBfinalt_BNd6XNCXMult_d05, SBfinalcurrents_BNd6XNCXMult_d05, timeImport_BNd6XNCXMult_d05, yImport_BNd6XNCXMult_d05, currentsImport_BNd6XNCXMult_d05] = ShannonBers_rabbit(0.6, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd8XNCXMult_d05, SBfinalt_BNd8XNCXMult_d05, SBfinalcurrents_BNd8XNCXMult_d05, timeImport_BNd8XNCXMult_d05, yImport_BNd8XNCXMult_d05, currentsImport_BNd8XNCXMult_d05] = ShannonBers_rabbit(0.8, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BN1XNCXMult_d05, SBfinalt_BN1XNCXMult_d05, SBfinalcurrents_BN1XNCXMult_d05, timeImport_BN1XNCXMult_d05, yImport_BN1XNCXMult_d05, currentsImport_BN1XNCXMult_d05] = ShannonBers_rabbit(1., 0, numBeats, CycleLength, param, diffRateBool);





%Multiplier = .1
param.NCXMult = .1;

%Healthy
[SBfinalyNCXMult_d1, SBfinaltNCXMult_d1, SBfinalcurrentsNCXMult_d1, timeImportNCXMult_d1, yImportNCXMult_d1, currentsImportNCXMult_d1] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);
%BacNav .25X and .5X and 1X
[SBfinaly_BNd1XNCXMult_d1, SBfinalt_BNd1XNCXMult_d1, SBfinalcurrents_BNd1XNCXMult_d1, timeImport_BNd1XNCXMult_d1, yImport_BNd1XNCXMult_d1, currentsImport_BNd1XNCXMult_d1] = ShannonBers_rabbit(0.1, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd2XNCXMult_d1, SBfinalt_BNd2XNCXMult_d1, SBfinalcurrents_BNd2XNCXMult_d1, timeImport_BNd2XNCXMult_d1, yImport_BNd2XNCXMult_d1, currentsImport_BNd2XNCXMult_d1] = ShannonBers_rabbit(0.2, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd4XNCXMult_d1, SBfinalt_BNd4XNCXMult_d1, SBfinalcurrents_BNd4XNCXMult_d1, timeImport_BNd4XNCXMult_d1, yImport_BNd4XNCXMult_d1, currentsImport_BNd4XNCXMult_d1] = ShannonBers_rabbit(0.4, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd6XNCXMult_d1, SBfinalt_BNd6XNCXMult_d1, SBfinalcurrents_BNd6XNCXMult_d1, timeImport_BNd6XNCXMult_d1, yImport_BNd6XNCXMult_d1, currentsImport_BNd6XNCXMult_d1] = ShannonBers_rabbit(0.6, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd8XNCXMult_d1, SBfinalt_BNd8XNCXMult_d1, SBfinalcurrents_BNd8XNCXMult_d1, timeImport_BNd8XNCXMult_d1, yImport_BNd8XNCXMult_d1, currentsImport_BNd8XNCXMult_d1] = ShannonBers_rabbit(0.8, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BN1XNCXMult_d1, SBfinalt_BN1XNCXMult_d1, SBfinalcurrents_BN1XNCXMult_d1, timeImport_BN1XNCXMult_d1, yImport_BN1XNCXMult_d1, currentsImport_BN1XNCXMult_d1] = ShannonBers_rabbit(1., 0, numBeats, CycleLength, param, diffRateBool);




%Multiplier = .2
param.NCXMult = .2;

%Healthy
[SBfinalyNCXMult_d2, SBfinaltNCXMult_d2, SBfinalcurrentsNCXMult_d2, timeImportNCXMult_d2, yImportNCXMult_d2, currentsImportNCXMult_d2] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);
%BacNav .25X and .5X and 1X
[SBfinaly_BNd1XNCXMult_d2, SBfinalt_BNd1XNCXMult_d2, SBfinalcurrents_BNd1XNCXMult_d2, timeImport_BNd1XNCXMult_d2, yImport_BNd1XNCXMult_d2, currentsImport_BNd1XNCXMult_d2] = ShannonBers_rabbit(0.1, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd2XNCXMult_d2, SBfinalt_BNd2XNCXMult_d2, SBfinalcurrents_BNd2XNCXMult_d2, timeImport_BNd2XNCXMult_d2, yImport_BNd2XNCXMult_d2, currentsImport_BNd2XNCXMult_d2] = ShannonBers_rabbit(0.2, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd4XNCXMult_d2, SBfinalt_BNd4XNCXMult_d2, SBfinalcurrents_BNd4XNCXMult_d2, timeImport_BNd4XNCXMult_d2, yImport_BNd4XNCXMult_d2, currentsImport_BNd4XNCXMult_d2] = ShannonBers_rabbit(0.4, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd6XNCXMult_d2, SBfinalt_BNd6XNCXMult_d2, SBfinalcurrents_BNd6XNCXMult_d2, timeImport_BNd6XNCXMult_d2, yImport_BNd6XNCXMult_d2, currentsImport_BNd6XNCXMult_d2] = ShannonBers_rabbit(0.6, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd8XNCXMult_d2, SBfinalt_BNd8XNCXMult_d2, SBfinalcurrents_BNd8XNCXMult_d2, timeImport_BNd8XNCXMult_d2, yImport_BNd8XNCXMult_d2, currentsImport_BNd8XNCXMult_d2] = ShannonBers_rabbit(0.8, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BN1XNCXMult_d2, SBfinalt_BN1XNCXMult_d2, SBfinalcurrents_BN1XNCXMult_d2, timeImport_BN1XNCXMult_d2, yImport_BN1XNCXMult_d2, currentsImport_BN1XNCXMult_d2] = ShannonBers_rabbit(1., 0, numBeats, CycleLength, param, diffRateBool);



%Multiplier = .4
param.NCXMult = .4;

%Healthy
[SBfinalyNCXMult_d4, SBfinaltNCXMult_d4, SBfinalcurrentsNCXMult_d4, timeImportNCXMult_d4, yImportNCXMult_d4, currentsImportNCXMult_d4] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);
%BacNav .25X and .5X and 1X
[SBfinaly_BNd1XNCXMult_d4, SBfinalt_BNd1XNCXMult_d4, SBfinalcurrents_BNd1XNCXMult_d4, timeImport_BNd1XNCXMult_d4, yImport_BNd1XNCXMult_d4, currentsImport_BNd1XNCXMult_d4] = ShannonBers_rabbit(0.1, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd2XNCXMult_d4, SBfinalt_BNd2XNCXMult_d4, SBfinalcurrents_BNd2XNCXMult_d4, timeImport_BNd2XNCXMult_d4, yImport_BNd2XNCXMult_d4, currentsImport_BNd2XNCXMult_d4] = ShannonBers_rabbit(0.2, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd4XNCXMult_d4, SBfinalt_BNd4XNCXMult_d4, SBfinalcurrents_BNd4XNCXMult_d4, timeImport_BNd4XNCXMult_d4, yImport_BNd4XNCXMult_d4, currentsImport_BNd4XNCXMult_d4] = ShannonBers_rabbit(0.4, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd6XNCXMult_d4, SBfinalt_BNd6XNCXMult_d4, SBfinalcurrents_BNd6XNCXMult_d4, timeImport_BNd6XNCXMult_d4, yImport_BNd6XNCXMult_d4, currentsImport_BNd6XNCXMult_d4] = ShannonBers_rabbit(0.6, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd8XNCXMult_d4, SBfinalt_BNd8XNCXMult_d4, SBfinalcurrents_BNd8XNCXMult_d4, timeImport_BNd8XNCXMult_d4, yImport_BNd8XNCXMult_d4, currentsImport_BNd8XNCXMult_d4] = ShannonBers_rabbit(0.8, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BN1XNCXMult_d4, SBfinalt_BN1XNCXMult_d4, SBfinalcurrents_BN1XNCXMult_d4, timeImport_BN1XNCXMult_d4, yImport_BN1XNCXMult_d4, currentsImport_BN1XNCXMult_d4] = ShannonBers_rabbit(1., 0, numBeats, CycleLength, param, diffRateBool);





%Multiplier = .5
param.NCXMult = .5;

%Healthy
[SBfinalyNCXMult_d5, SBfinaltNCXMult_d5, SBfinalcurrentsNCXMult_d5, timeImportNCXMult_d5, yImportNCXMult_d5, currentsImportNCXMult_d5] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);
%BacNav .25X and .5X and 1X
[SBfinaly_BNd1XNCXMult_d5, SBfinalt_BNd1XNCXMult_d5, SBfinalcurrents_BNd1XNCXMult_d5, timeImport_BNd1XNCXMult_d5, yImport_BNd1XNCXMult_d5, currentsImport_BNd1XNCXMult_d5] = ShannonBers_rabbit(0.1, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd2XNCXMult_d5, SBfinalt_BNd2XNCXMult_d5, SBfinalcurrents_BNd2XNCXMult_d5, timeImport_BNd2XNCXMult_d5, yImport_BNd2XNCXMult_d5, currentsImport_BNd2XNCXMult_d5] = ShannonBers_rabbit(0.2, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd4XNCXMult_d5, SBfinalt_BNd4XNCXMult_d5, SBfinalcurrents_BNd4XNCXMult_d5, timeImport_BNd4XNCXMult_d5, yImport_BNd4XNCXMult_d5, currentsImport_BNd4XNCXMult_d5] = ShannonBers_rabbit(0.4, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd6XNCXMult_d5, SBfinalt_BNd6XNCXMult_d5, SBfinalcurrents_BNd6XNCXMult_d5, timeImport_BNd6XNCXMult_d5, yImport_BNd6XNCXMult_d5, currentsImport_BNd6XNCXMult_d5] = ShannonBers_rabbit(0.6, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd8XNCXMult_d5, SBfinalt_BNd8XNCXMult_d5, SBfinalcurrents_BNd8XNCXMult_d5, timeImport_BNd8XNCXMult_d5, yImport_BNd8XNCXMult_d5, currentsImport_BNd8XNCXMult_d5] = ShannonBers_rabbit(0.8, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BN1XNCXMult_d5, SBfinalt_BN1XNCXMult_d5, SBfinalcurrents_BN1XNCXMult_d5, timeImport_BN1XNCXMult_d5, yImport_BN1XNCXMult_d5, currentsImport_BN1XNCXMult_d5] = ShannonBers_rabbit(1., 0, numBeats, CycleLength, param, diffRateBool);










%Multiplier = .01
param.NCXMult = .01;

%Healthy
[SBfinalyNCXMult_d01, SBfinaltNCXMult_d01, SBfinalcurrentsNCXMult_d01, timeImportNCXMult_d01, yImportNCXMult_d01, currentsImportNCXMult_d01] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);
%BacNav .25X and .5X and 1X
[SBfinaly_BNd1XNCXMult_d01, SBfinalt_BNd1XNCXMult_d01, SBfinalcurrents_BNd1XNCXMult_d01, timeImport_BNd1XNCXMult_d01, yImport_BNd1XNCXMult_d01, currentsImport_BNd1XNCXMult_d01] = ShannonBers_rabbit(0.1, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd2XNCXMult_d01, SBfinalt_BNd2XNCXMult_d01, SBfinalcurrents_BNd2XNCXMult_d01, timeImport_BNd2XNCXMult_d01, yImport_BNd2XNCXMult_d01, currentsImport_BNd2XNCXMult_d01] = ShannonBers_rabbit(0.2, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd4XNCXMult_d01, SBfinalt_BNd4XNCXMult_d01, SBfinalcurrents_BNd4XNCXMult_d01, timeImport_BNd4XNCXMult_d01, yImport_BNd4XNCXMult_d01, currentsImport_BNd4XNCXMult_d01] = ShannonBers_rabbit(0.4, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd6XNCXMult_d01, SBfinalt_BNd6XNCXMult_d01, SBfinalcurrents_BNd6XNCXMult_d01, timeImport_BNd6XNCXMult_d01, yImport_BNd6XNCXMult_d01, currentsImport_BNd6XNCXMult_d01] = ShannonBers_rabbit(0.6, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd8XNCXMult_d01, SBfinalt_BNd8XNCXMult_d01, SBfinalcurrents_BNd8XNCXMult_d01, timeImport_BNd8XNCXMult_d01, yImport_BNd8XNCXMult_d01, currentsImport_BNd8XNCXMult_d01] = ShannonBers_rabbit(0.8, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BN1XNCXMult_d01, SBfinalt_BN1XNCXMult_d01, SBfinalcurrents_BN1XNCXMult_d01, timeImport_BN1XNCXMult_d01, yImport_BN1XNCXMult_d01, currentsImport_BN1XNCXMult_d01] = ShannonBers_rabbit(1., 0, numBeats, CycleLength, param, diffRateBool);




%Multiplier = 0.
param.NCXMult = 0.;

%Healthy
[SBfinalyNCXMult_0, SBfinaltNCXMult_0, SBfinalcurrentsNCXMult_0, timeImportNCXMult_0, yImportNCXMult_0, currentsImportNCXMult_0] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);
%BacNav .25X and .5X and 1X
[SBfinaly_BNd1XNCXMult_0, SBfinalt_BNd1XNCXMult_0, SBfinalcurrents_BNd1XNCXMult_0, timeImport_BNd1XNCXMult_0, yImport_BNd1XNCXMult_0, currentsImport_BNd1XNCXMult_0] = ShannonBers_rabbit(0.1, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd2XNCXMult_0, SBfinalt_BNd2XNCXMult_0, SBfinalcurrents_BNd2XNCXMult_0, timeImport_BNd2XNCXMult_0, yImport_BNd2XNCXMult_0, currentsImport_BNd2XNCXMult_0] = ShannonBers_rabbit(0.2, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd4XNCXMult_0, SBfinalt_BNd4XNCXMult_0, SBfinalcurrents_BNd4XNCXMult_0, timeImport_BNd4XNCXMult_0, yImport_BNd4XNCXMult_0, currentsImport_BNd4XNCXMult_0] = ShannonBers_rabbit(0.4, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd6XNCXMult_0, SBfinalt_BNd6XNCXMult_0, SBfinalcurrents_BNd6XNCXMult_0, timeImport_BNd6XNCXMult_0, yImport_BNd6XNCXMult_0, currentsImport_BNd6XNCXMult_0] = ShannonBers_rabbit(0.6, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd8XNCXMult_0, SBfinalt_BNd8XNCXMult_0, SBfinalcurrents_BNd8XNCXMult_0, timeImport_BNd8XNCXMult_0, yImport_BNd8XNCXMult_0, currentsImport_BNd8XNCXMult_0] = ShannonBers_rabbit(0.8, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BN1XNCXMult_0, SBfinalt_BN1XNCXMult_0, SBfinalcurrents_BN1XNCXMult_0, timeImport_BN1XNCXMult_0, yImport_BN1XNCXMult_0, currentsImport_BN1XNCXMult_0] = ShannonBers_rabbit(1., 0, numBeats, CycleLength, param, diffRateBool);










%Multiplier = .6
param.NCXMult = .6;

%Healthy
[SBfinalyNCXMult_d6, SBfinaltNCXMult_d6, SBfinalcurrentsNCXMult_d6, timeImportNCXMult_d6, yImportNCXMult_d6, currentsImportNCXMult_d6] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);
%BacNav .25X and .5X and 1X
[SBfinaly_BNd1XNCXMult_d6, SBfinalt_BNd1XNCXMult_d6, SBfinalcurrents_BNd1XNCXMult_d6, timeImport_BNd1XNCXMult_d6, yImport_BNd1XNCXMult_d6, currentsImport_BNd1XNCXMult_d6] = ShannonBers_rabbit(0.1, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd2XNCXMult_d6, SBfinalt_BNd2XNCXMult_d6, SBfinalcurrents_BNd2XNCXMult_d6, timeImport_BNd2XNCXMult_d6, yImport_BNd2XNCXMult_d6, currentsImport_BNd2XNCXMult_d6] = ShannonBers_rabbit(0.2, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd4XNCXMult_d6, SBfinalt_BNd4XNCXMult_d6, SBfinalcurrents_BNd4XNCXMult_d6, timeImport_BNd4XNCXMult_d6, yImport_BNd4XNCXMult_d6, currentsImport_BNd4XNCXMult_d6] = ShannonBers_rabbit(0.4, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd6XNCXMult_d6, SBfinalt_BNd6XNCXMult_d6, SBfinalcurrents_BNd6XNCXMult_d6, timeImport_BNd6XNCXMult_d6, yImport_BNd6XNCXMult_d6, currentsImport_BNd6XNCXMult_d6] = ShannonBers_rabbit(0.6, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd8XNCXMult_d6, SBfinalt_BNd8XNCXMult_d6, SBfinalcurrents_BNd8XNCXMult_d6, timeImport_BNd8XNCXMult_d6, yImport_BNd8XNCXMult_d6, currentsImport_BNd8XNCXMult_d6] = ShannonBers_rabbit(0.8, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BN1XNCXMult_d6, SBfinalt_BN1XNCXMult_d6, SBfinalcurrents_BN1XNCXMult_d6, timeImport_BN1XNCXMult_d6, yImport_BN1XNCXMult_d6, currentsImport_BN1XNCXMult_d6] = ShannonBers_rabbit(1., 0, numBeats, CycleLength, param, diffRateBool);




%Multiplier = .8
param.NCXMult = .8;

%Healthy
[SBfinalyNCXMult_d8, SBfinaltNCXMult_d8, SBfinalcurrentsNCXMult_d8, timeImportNCXMult_d8, yImportNCXMult_d8, currentsImportNCXMult_d8] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);
%BacNav .25X and .5X and 1X
[SBfinaly_BNd1XNCXMult_d8, SBfinalt_BNd1XNCXMult_d8, SBfinalcurrents_BNd1XNCXMult_d8, timeImport_BNd1XNCXMult_d8, yImport_BNd1XNCXMult_d8, currentsImport_BNd1XNCXMult_d8] = ShannonBers_rabbit(0.1, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd2XNCXMult_d8, SBfinalt_BNd2XNCXMult_d8, SBfinalcurrents_BNd2XNCXMult_d8, timeImport_BNd2XNCXMult_d8, yImport_BNd2XNCXMult_d8, currentsImport_BNd2XNCXMult_d8] = ShannonBers_rabbit(0.2, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd4XNCXMult_d8, SBfinalt_BNd4XNCXMult_d8, SBfinalcurrents_BNd4XNCXMult_d8, timeImport_BNd4XNCXMult_d8, yImport_BNd4XNCXMult_d8, currentsImport_BNd4XNCXMult_d8] = ShannonBers_rabbit(0.4, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd6XNCXMult_d8, SBfinalt_BNd6XNCXMult_d8, SBfinalcurrents_BNd6XNCXMult_d8, timeImport_BNd6XNCXMult_d8, yImport_BNd6XNCXMult_d8, currentsImport_BNd6XNCXMult_d8] = ShannonBers_rabbit(0.6, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd8XNCXMult_d8, SBfinalt_BNd8XNCXMult_d8, SBfinalcurrents_BNd8XNCXMult_d8, timeImport_BNd8XNCXMult_d8, yImport_BNd8XNCXMult_d8, currentsImport_BNd8XNCXMult_d8] = ShannonBers_rabbit(0.8, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BN1XNCXMult_d8, SBfinalt_BN1XNCXMult_d8, SBfinalcurrents_BN1XNCXMult_d8, timeImport_BN1XNCXMult_d8, yImport_BN1XNCXMult_d8, currentsImport_BN1XNCXMult_d8] = ShannonBers_rabbit(1., 0, numBeats, CycleLength, param, diffRateBool);





%Multiplier = 1.
param.NCXMult = 1.;

%Healthy
[SBfinalyNCXMult_1, SBfinaltNCXMult_1, SBfinalcurrentsNCXMult_1, timeImportNCXMult_1, yImportNCXMult_1, currentsImportNCXMult_1] = ShannonBers_rabbit(0., 0, numBeats, CycleLength, param, diffRateBool);
%BacNav .25X and .5X and 1X
[SBfinaly_BNd1XNCXMult_1, SBfinalt_BNd1XNCXMult_1, SBfinalcurrents_BNd1XNCXMult_1, timeImport_BNd1XNCXMult_1, yImport_BNd1XNCXMult_1, currentsImport_BNd1XNCXMult_1] = ShannonBers_rabbit(0.1, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd2XNCXMult_1, SBfinalt_BNd2XNCXMult_1, SBfinalcurrents_BNd2XNCXMult_1, timeImport_BNd2XNCXMult_1, yImport_BNd2XNCXMult_1, currentsImport_BNd2XNCXMult_1] = ShannonBers_rabbit(0.2, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd4XNCXMult_1, SBfinalt_BNd4XNCXMult_1, SBfinalcurrents_BNd4XNCXMult_1, timeImport_BNd4XNCXMult_1, yImport_BNd4XNCXMult_1, currentsImport_BNd4XNCXMult_1] = ShannonBers_rabbit(0.4, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd6XNCXMult_1, SBfinalt_BNd6XNCXMult_1, SBfinalcurrents_BNd6XNCXMult_1, timeImport_BNd6XNCXMult_1, yImport_BNd6XNCXMult_1, currentsImport_BNd6XNCXMult_1] = ShannonBers_rabbit(0.6, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BNd8XNCXMult_1, SBfinalt_BNd8XNCXMult_1, SBfinalcurrents_BNd8XNCXMult_1, timeImport_BNd8XNCXMult_1, yImport_BNd8XNCXMult_1, currentsImport_BNd8XNCXMult_1] = ShannonBers_rabbit(0.8, 0, numBeats, CycleLength, param, diffRateBool);
[SBfinaly_BN1XNCXMult_1, SBfinalt_BN1XNCXMult_1, SBfinalcurrents_BN1XNCXMult_1, timeImport_BN1XNCXMult_1, yImport_BN1XNCXMult_1, currentsImport_BN1XNCXMult_1] = ShannonBers_rabbit(1., 0, numBeats, CycleLength, param, diffRateBool);






toc


%%

tic
%[APDduration, indStart, indEnd] = calculateAPD(t,v,threshPerc)



%APD, peak Cai, peak negative INCX, peak positive INCX, positive area INCX,
%negative area INCX


tau_m = [.1, .2, .5, 1, 2, 5, 10];
tau_h = [.1, .2, .5, 1, 2, 5, 10];

threshAPD = 80;

APD_Matrix = zeros(length(tau_m), length(tau_h));
peakCai_Matrix = zeros(length(tau_m), length(tau_h));
peakN_INCX_Matrix = zeros(length(tau_m), length(tau_h));
peakP_INCX_Matrix = zeros(length(tau_m), length(tau_h));
PArea_INCX_Matrix = zeros(length(tau_m), length(tau_h));
NArea_INCX_Matrix = zeros(length(tau_m), length(tau_h));


for i=1:length(tau_m)
    for j=1:length(tau_h)

        param.taum_scale = tau_m(i);
        param.tauh_scale = tau_h(j);

        [SBfinaly_BN1X_tauScaled, SBfinalt_BN1X_tauScaled, SBfinalcurrents_BN1X_tauScaled, timeImport_BN1X_tauScaled, yImport_BN1X_tauScaled, currentsImport_BN1X_tauScaled] = ShannonBers_rabbit(1., 0, numBeats, CycleLength, param, diffRateBool);

        [APDduration, indStart, indEnd] = calculateAPD(SBfinalt_BN1X_tauScaled,SBfinaly_BN1X_tauScaled(:,39),threshAPD);
        [positive_area, negative_area] = calculateAreaUnderComplexCurve(SBfinalt_BN1X_tauScaled, SBfinalcurrents_BN1X_tauScaled(:,6));

        APD_Matrix(i,j) = APDduration;
        peakCai_Matrix(i,j) = max(SBfinaly_BN1X_tauScaled(:,38));
        peakN_INCX_Matrix(i,j) = min(SBfinalcurrents_BN1X_tauScaled(:,6));
        peakP_INCX_Matrix(i,j) = max(SBfinalcurrents_BN1X_tauScaled(:,6));
        PArea_INCX_Matrix(i,j) = positive_area;
        NArea_INCX_Matrix(i,j) = negative_area;
         
    end
end

toc





