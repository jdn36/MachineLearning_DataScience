%NOW GET THE CV BASED ON MAX DVM/DT
%propagation stuff
rhoex = 1.5;%Ohms-cm interstitial resistivity 300 for anisotropy
rhoey = 1.5;%Ohms-cm interstitial resistivity 1000 for anisotropy
a = .002;%cm %radius of myocyte
dX = .005;%cm
dY = .005;%cm
Cm = 1;
% L = 2;%cm
rex = rhoex*dX/(pi*(a)^2);%mOhms
rey = rhoey*dY/(pi*(a)^2);%mOhms
Rxx = rex.*2.*pi.*a.*dX;
Ryy = rey.*2.*pi.*a.*dY;

%BACNAV STUFF/BRUGADA
boolBrugada = false;
percBacNav = 0.;


%Simulation Specific stuff
nodes = 6;
numNodesSt = 75;
dimen = 1;
iStMag = -18;%pA/pF
dt = .007; %ms
tFinal = 1000;%ms
t = 0:dt:tFinal;%ms
startOfPulse = 0.;%ms
I_begin = [];
I_begin(1) = startOfPulse;
numPaces = 1000;
paceInterval = 1000;%ms basically if 1000ms it is 1Hz
for i = 2:numPaces+2
    I_begin(i) = I_begin(i-1) + paceInterval;
end
I_stim = 0; %pA/pF
I_dur = 2; %ms
Icounter = 1;
inAP = 0;

%Initializing
Vm = []; h = []; d = []; x_r = []; Na_i = []; K_i = []; Ca_rel = []; o_i = []; u_i = []; 
Cmdn_Ca_i = []; Csqn_Ca_i = []; v = []; m = []; j = []; f = []; x_s = []; Ca_i = []; Ca_up = [];
o_a = []; u_a = []; f_Ca = []; Trpn_Ca_i = []; u = []; w = []; dVm_dtVec = []; dVm_dtMax = []; LATvec = [];
%BacNav stuff
hBN = []; mBN = [];

%1D
    Vm = zeros(nodes,length(t)); 
    h = zeros(1,nodes);
    d = zeros(1,nodes);
    x_r = zeros(1,nodes);
    Na_i = zeros(1,nodes);
    K_i = zeros(1,nodes);
    Ca_rel = zeros(1,nodes);
    o_i = zeros(1,nodes);
    u_i = zeros(1,nodes);
    Cmdn_Ca_i = zeros(1,nodes);
    Csqn_Ca_i = zeros(1,nodes);
    v = zeros(1,nodes);
    m = zeros(1,nodes);
    j = zeros(1,nodes);
    f = zeros(1,nodes);
    x_s = zeros(1,nodes);
    Ca_i = zeros(1,nodes);
    Ca_up = zeros(1,nodes);
    o_a = zeros(1,nodes);
    u_a = zeros(1,nodes);
    f_Ca = zeros(1,nodes);
    Trpn_Ca_i = zeros(1,nodes);
    u = zeros(1,nodes);
    w = zeros(1,nodes);
    dVm_dtVec = zeros(1,nodes);
    dVm_dtMax = zeros(1,nodes);
    LATvec = zeros(1,nodes);
    
    %BacNav stuff
    hBN = zeros(1,nodes);
    mBN = zeros(1,nodes);
    
    pointsToLoop = nodes;

%IC
Vm(:,1) = -81.2; %mV
% Vmtrial(1) = -81.2; %mV
h(:) = .965;
d(:) = .000137;
x_r(:) = .0000329;
Na_i(:) = 11.2;
K_i(:) = 139;
Ca_rel(:) = 1.49;
o_i(:) = .999;
u_i(:) = .999;
Cmdn_Ca_i(:) = .00205;
Csqn_Ca_i(:) = 6.51;
v(:) = 1;
m(:) = .00291;
j(:) = .978;
f(:) = .999;
x_s(:) = .0187;
Ca_i(:) = .000102;
Ca_up(:) = 1.49;
o_a(:) = .0304;
u_a(:) = .00496;
f_Ca(:) = .775;
Trpn_Ca_i(:) = .0118;
u(:) = 0.0;
w(:) = .999;
hBN(:) = 0.8231; %IC for BacNav
mBN(:) = 0.0000094; %IC for BacNav

g_Na_BacNav = 21.6399; %nS/pF max I_Na_BacNav conductance   (21.6399 is original)  instead of 10 to affect peak of notch in AP (the higher it is, the higher the peak is)
tauhParam = 84.6609;%parameter within calculation of tau_h_BN to affect AP shape, the higher it is the shorter the APD is... original is 84.6609 instead of 90.

I_BacNavV = zeros(pointsToLoop,length(t));

% i2=1;
% i3=1;

%node 1 is healthy, node 2 is BacNav with healthy
tic
%double for loop, outer is time inner is nodes
    for i1 = 1:length(t)-1
        %defining where the stimulus is happening
        if t(i1) >= I_begin(Icounter) && t(i1) <= (I_begin(Icounter)+I_dur)
            I_stim = iStMag;
            inAP = 1;
    %         display("here")
        else
            I_stim = 0;
            if inAP == 1
                Icounter = Icounter+1;
                inAP = 0;
            end
            
        end
        
        for i2 = 1:1:pointsToLoop

            %getting the states vector for current iteration
            statesCMpr = [hBN(i2), mBN(i2), Vm(i2,i1), h(i2), d(i2), x_r(i2), Na_i(i2), K_i(i2), Ca_rel(i2),...
                o_i(i2), u_i(i2), Cmdn_Ca_i(i2), Csqn_Ca_i(i2), v(i2),m(i2), j(i2), f(i2), x_s(i2),...
                Ca_i(i2), Ca_up(i2), o_a(i2), u_a(i2), f_Ca(i2), Trpn_Ca_i(i2), u(i2), w(i2)];
        
            %calling the specific cell model for one cell one step at a time
            if i2 == 1
                percBacNav = 0;
            elseif i2 == 2
                percBacNav = 0.2;
            elseif i2 == 3
                percBacNav = 0.5;
            elseif i2 == 4
                percBacNav = 1.;
            elseif i2 == 5
                percBacNav = 1.5;
            elseif i2 == 6
                percBacNav = 2.;
            end
    
            [dVm_dt, statesCMpo, I_bacnav] = CourtemancheFunction(statesCMpr, I_stim, dt, boolBrugada, percBacNav, g_Na_BacNav, tauhParam);
            Vm(i2,i1+1) = Vm(i2,i1) + dt.*(dVm_dt);
            I_BacNavV(i2,i1+1) = I_bacnav;
        
            %getting the states variables for next iteration
            hBN(i2) = statesCMpo(1); mBN(i2) = statesCMpo(2); h(i2) = statesCMpo(4); d(i2) = statesCMpo(5);
            x_r(i2) = statesCMpo(6); Na_i(i2) = statesCMpo(7); K_i(i2) = statesCMpo(8); Ca_rel(i2) = statesCMpo(9); o_i(i2) = statesCMpo(10);
            u_i(i2) = statesCMpo(11); Cmdn_Ca_i(i2) = statesCMpo(12); Csqn_Ca_i(i2) = statesCMpo(13); v(i2) = statesCMpo(14); m(i2) = statesCMpo(15);
            j(i2) = statesCMpo(16); f(i2) = statesCMpo(17); x_s(i2) = statesCMpo(18); Ca_i(i2) = statesCMpo(19); Ca_up(i2) = statesCMpo(20);
            o_a(i2)=statesCMpo(21);u_a(i2)=statesCMpo(22);f_Ca(i2)=statesCMpo(23);Trpn_Ca_i(i2)=statesCMpo(24);u(i2)=statesCMpo(25);w(i2)=statesCMpo(26);
            
        end
    end


toc

%%
%plotting animation
    
        %node 1
        plot(t((length(t) - (ceil((4*paceInterval)/dt)):1:length(t)))./1000, Vm(1,(length(t) - (ceil((4*paceInterval)/dt)):1:length(t))))
        hold on
        %node 2
        plot(t((length(t) - (ceil((4*paceInterval)/dt)):1:length(t)))./1000, Vm(2,(length(t) - (ceil((4*paceInterval)/dt)):1:length(t))))
        %node 3
        plot(t((length(t) - (ceil((4*paceInterval)/dt)):1:length(t)))./1000, Vm(3,(length(t) - (ceil((4*paceInterval)/dt)):1:length(t))))
        %node 4
        plot(t((length(t) - (ceil((4*paceInterval)/dt)):1:length(t)))./1000, Vm(4,(length(t) - (ceil((4*paceInterval)/dt)):1:length(t))))
        %node 5
        plot(t((length(t) - (ceil((4*paceInterval)/dt)):1:length(t)))./1000, Vm(5,(length(t) - (ceil((4*paceInterval)/dt)):1:length(t))))
        %node 6
        plot(t((length(t) - (ceil((4*paceInterval)/dt)):1:length(t)))./1000, Vm(6,(length(t) - (ceil((4*paceInterval)/dt)):1:length(t))))
        


        %node 1
        plot(t./1000, Vm(1,:))
        hold on
        %node 2
        plot(t./1000, Vm(2,:))
        %node 3
        plot(t./1000, Vm(3,:))
        %node 4
        plot(t./1000, Vm(4,:))
        %node 5
        plot(t./1000, Vm(5,:))
        %node 6
        plot(t./1000, Vm(6,:))


        figure
        %node 1
        plot(t./1000, I_BacNavV(1,:))
        hold on
        %node 2
        plot(t./1000, I_BacNavV(2,:))
        %node 3
        plot(t./1000, I_BacNavV(3,:))
        %node 4
        plot(t./1000, I_BacNavV(4,:))
        %node 5
        plot(t./1000, I_BacNavV(5,:))
        %node 6
        plot(t./1000, I_BacNavV(6,:))


        title('Comparison of BacNav influence at steady-state')
        xlabel('time (s)')
        ylabel('Vm (mV)')
        legend('Control','.2X','.5X','1X','1.5X','2X')
    
