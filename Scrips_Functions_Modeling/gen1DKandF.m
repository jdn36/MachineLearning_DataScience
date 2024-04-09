%DON'T USE HUGE STIMULUS
%ALWAYS CHOOSE AN ODD NUMBER FOR NODES --- node number in terms of
%electrical points NOT material points

%Four parts: Parameters, Calculations (electrical and mechanical), Creation of Final Vectors,
%Displaying simulation


%Part 1: Parameters

%simulation time and equation parameters
tic;
dT = .01;
endOfTime = 670;
durationOfPulse = 10;
startOfSecPulse = 360;
t = 0:dT:endOfTime;
Cm = 1;
rhoex = 40;%Ohms-cm interstitial resistivity 300 for anisotropy
rhoey = 20;%Ohms-cm interstitial resistivity 1000 for anisotropy
a = .0010;%cm %radius of myocyte
dX = .0035;%cm
dY = .0035;%cm
nodes = 41;
% L = 2;%cm
rex = rhoex*dX/(pi*(a)^2);%mOhms
rey = rhoey*dY/(pi*(a)^2);%mOhms
percentageChangeI = 1;

%mechanical stuff
kT = 1.5;
kStiffRat = .5;
aThresh = .05;
r0 = dX*2;
dT2 = .001;
c_mech = .00001; %.0001
d_mech = c_mech*30; %.001; %.0035;
scalingTa = .00045;%.0047; %.0047

%adapting the threshold value
accThreshTot = round(((1.05e-4)*(((nodes+1)/2)*((nodes+1)/2))),5);%*2;
threshPerc = 1;

minVal = -84;
maxVal = 29;


Vm = zeros(nodes*nodes,2); %1st: x, 2nd: y, 3rd: time
epsi = (1)*ones((((nodes+1)/2)*((nodes+1)/2)),1);
Vmsecondary  = zeros(nodes*nodes,101); %1st: x, 2nd: y, 3rd: time
Vmnondim  = zeros(nodes*nodes,2); %1st: x, 2nd: y, 3rd: time
Ta = zeros((((nodes+1)/2)*((nodes+1)/2)),2);
Tasecondary = zeros((((nodes+1)/2)*((nodes+1)/2)),101);
massX = zeros((((nodes+1)/2)*((nodes+1)/2)),2);
massY = zeros((((nodes+1)/2)*((nodes+1)/2)),2);
massXsecondary = zeros((((nodes+1)/2)*((nodes+1)/2)),101);
massYsecondary = zeros((((nodes+1)/2)*((nodes+1)/2)),101);
massXdd = zeros((((nodes+1)/2)*((nodes+1)/2)),1);
massYdd = zeros((((nodes+1)/2)*((nodes+1)/2)),1);
newTime = zeros(1,101);

Ta(:,1) = 0;
Vm(:,1) = -84;
Vmsecondary(:,1) = -84;
Vmnondim(:,1) = 0;


%Creation of mass lattice
leftValuex = -floor(((nodes+1)/2)/2)*r0;
bottomValuey = -floor(((nodes+1)/2)/2)*r0;
rightValuex = floor(((nodes+1)/2)/2)*r0;
upperValuey = floor(((nodes+1)/2)/2)*r0;

massX( 1:((nodes+1)/2):(((nodes+1)/2)^2)-(nodes+1)/2+1, 1 ) = leftValuex;
massY( (((nodes+1)/2)^2)-(nodes+1)/2+1:1:(((nodes+1)/2)^2), 1 ) = bottomValuey;

for coli = 1:1:((nodes+1)/2)-1%ceil((((nodes+1)/2)^2)/2):1:(((nodes+1)/2)^2)-1
  
    massX( coli+1:((nodes+1)/2):end, 1 ) = massX( coli:((nodes+1)/2):end, 1 ) + r0;
    massY( (((nodes+1)/2)^2)-(coli+1)*(nodes+1)/2+1:1:(((nodes+1)/2)^2)-(coli)*((nodes+1)/2), 1 ) = massY( (((nodes+1)/2)^2)-(coli)*(nodes+1)/2+1:1:(((nodes+1)/2)^2)-(coli-1)*((nodes+1)/2), 1 ) + r0;
    
end

massXsecondary(:,1) = massX(:,1);
massYsecondary(:,1) = massY(:,1);

% firstMassX = massX(1,1);
% lastMassX = massX(end,1);

VNa = 50;
gbarNa = 4;
gNaC = .003;
gbars = .09;

%to create infinite resistances to simulate blockages in the tissue
%blockageValue = 1 no blocks
blockageValue = 1;

%initializing vectors
% alphaX1 = zeros(nodes*nodes,1);
% alphaM = zeros(nodes*nodes,1);
% alphaH = zeros(nodes*nodes,1);
% alphaJ = zeros(nodes*nodes,1);
% alphaD = zeros(nodes*nodes,1);
% alphaF = zeros(nodes*nodes,1);
% betaX1 = zeros(nodes*nodes,1);
% betaM = zeros(nodes*nodes,1);
% betaH = zeros(nodes*nodes,1);
% n = zeros(nodes*nodes,1);
% m = zeros(nodes*nodes,1);
% h = zeros(nodes*nodes,1);
CaI = zeros(nodes*nodes,1);
CaI(:,1)=3.*10.^(-7);%mole/L

alphaX1 = gatingPF(.0005,.083,50,0,0,.057,1,Vm(:,1));
alphaM = gatingPF(0,0,47,-1,47,-.1,-1,Vm(:,1));
alphaH = gatingPF(.126,-.25,77,0,0,0,0,Vm(:,1));
alphaJ = gatingPF(.055,-.25,78,0,0,-.2,1,Vm(:,1));
alphaD = gatingPF(.095,-.01,-5,0,0,-.072,1,Vm(:,1));
alphaF = gatingPF(.012,-.008,28,0,0,.15,1,Vm(:,1));
betaX1 = gatingPF(.0013,-.06,20,0,0,-.04,1,Vm(:,1));
betaM = gatingPF(40,-.056,72,0,0,0,0,Vm(:,1));
betaH = gatingPF(1.7,0,22.5,0,0,-.082,1,Vm(:,1));
betaJ = gatingPF(.3,0,32,0,0,-.1,1,Vm(:,1));
betaD = gatingPF(.07,-.017,44,0,0,.05,1,Vm(:,1));
betaF = gatingPF(.0065,-.02,30,0,0,-.2,1,Vm(:,1));
%n_inf,m_inf,h_inf
x1 = (alphaX1) ./ (alphaX1 + betaX1);
m = (alphaM) ./ (alphaM + betaM);
h = (alphaH) ./ (alphaH + betaH);
jota = (alphaJ) ./ (alphaJ + betaJ);
d = (alphaD) ./ (alphaD + betaD);
f = (alphaF) ./ (alphaF + betaF);

% iK = [];
% iNa = [];
% iAn = [];
% iK(1) = (1.2*(exp((-Vm(1) - 90) / 50)) + .015*exp((Vm(1) + 90) / 60))*(Vm(1) - VK) + 1.2*((n(1))^4)*(Vm(1) - VK);
% iNa(1) = (400*(m(1)^3)*h(1) + .140)*(Vm(1) - VNa);
% iAn(1) = 0;

%After Cl- current introduction:
% iAn(1) = .075*(Vm(1) - VCl);

iSt = 30;
indexSecondary = 2;
indexPrimary = 1;

A = zeros(nodes*nodes,nodes*nodes);

    
    
    g = zeros(nodes*nodes,1);
    
    Ainternal = A; 
    
    iSt2 = 0;
    

%Matrix stuff for the mass and Ta operations

Axyddn = zeros(((nodes+1)/2).^2,((nodes+1)/2).^2);
Axyddnprev = zeros(((nodes+1)/2).^2,((nodes+1)/2).^2);

    for j = ((nodes+1)/2+2) : 1 : (((nodes+1)/2)^2 - ((nodes+1)/2) - 1)
        %kStiffRat
        Axyddnprev(j,j-((nodes+1)/2)-1) = kStiffRat.*d_mech/dT2;
        Axyddnprev(j,j-((nodes+1)/2)) = d_mech/dT2;
        Axyddnprev(j,j-((nodes+1)/2)+1) = kStiffRat.*d_mech/dT2;
        
        Axyddnprev(j,j-1) = d_mech/dT2;
        Axyddnprev(j,j) = -4.*d_mech/dT2 - 4.*kStiffRat.*d_mech/dT2;
        Axyddnprev(j,j+1) = d_mech/dT2;
        
        Axyddnprev(j,j+((nodes+1)/2)-1) = kStiffRat.*d_mech/dT2; 
        Axyddnprev(j,j+((nodes+1)/2)) = d_mech/dT2;
        Axyddnprev(j,j+((nodes+1)/2)+1) = kStiffRat.*d_mech/dT2; 
    end
    

%%

%Part 2: Calculations (electrical and mechanical)

%Eulers Method solving
for i=1:length(t)
   
   
   %accounting for both type of K current: gK1(Vm - VK) + gK2(Vm - VK)
   %gK1: g based on an instantaneous f(Vm)->decreases when Vm increases
   %gK1 NOT DEPENDING ON ACTIVATION OR INACTIVATION (n^4)
   %gK2: g f(Vm)->increases when Vm increases. n^4 dependent
   %Na current: gNa(Vm - VNa) with gNa being dependent on m^3*h
   
   iK1 = .35.*((4.*((exp(.04.*(Vm(:,indexPrimary)+85))) - (1))./((exp(.08.*(Vm(:,indexPrimary)+53))) + (exp(.04.*(Vm(:,indexPrimary)+53))))) + (.2.*((Vm(:,indexPrimary)+23)./(1-(exp(-.04.*(Vm(:,indexPrimary)+23)))))));
   iX1 = x1.*(.8.*((exp(.04.*(Vm(:,indexPrimary)+77)) - 1)./(exp(.04.*(Vm(:,indexPrimary)+35)))));
   iNa = (gbarNa.*m.^3.*h.*jota +gNaC).*(Vm(:,indexPrimary)-VNa);
   iS = (percentageChangeI).*(gbars.*d.*f.*(Vm(:,indexPrimary) - (-82.3 - 13.0287.*(log(CaI)))));
   
   alphaX1 = gatingPF(.0005,.083,50,0,0,.057,1,Vm(:,indexPrimary));
   alphaM = gatingPF(0,0,47,-1,47,-.1,-1,Vm(:,indexPrimary));
   alphaH = gatingPF(.126,-.25,77,0,0,0,0,Vm(:,indexPrimary));
   alphaJ = gatingPF(.055,-.25,78,0,0,-.2,1,Vm(:,indexPrimary));
   alphaD = gatingPF(.095,-.01,-5,0,0,-.072,1,Vm(:,indexPrimary));
   alphaF = gatingPF(.012,-.008,28,0,0,.15,1,Vm(:,indexPrimary));
   betaX1 = gatingPF(.0013,-.06,20,0,0,-.04,1,Vm(:,indexPrimary));
   betaM = gatingPF(40,-.056,72,0,0,0,0,Vm(:,indexPrimary));
   betaH = gatingPF(1.7,0,22.5,0,0,-.082,1,Vm(:,indexPrimary));
   betaJ = gatingPF(.3,0,32,0,0,-.1,1,Vm(:,indexPrimary));
   betaD = gatingPF(.07,-.017,44,0,0,.05,1,Vm(:,indexPrimary));
   betaF = gatingPF(.0065,-.02,30,0,0,-.2,1,Vm(:,indexPrimary));
   
   x1 = x1 + dT.*((alphaX1).*(1 - x1) - (betaX1).*x1);
   m = m + dT.*(alphaM.*(1 - m) - betaM.*m);
   h = h + dT.*(alphaH.*(1 - h) - betaH.*h);
   jota = jota + dT.*((alphaJ).*(1 - jota) - (betaJ).*jota);
   d = d + dT.*(alphaD.*(1 - d) - betaD.*d);
   f = f + dT.*(alphaF.*(1 - f) - betaF.*f);
   CaI = CaI+dT.*((-10.^(-7)).*iS + .07.*((10.^(-7)) - CaI));
   
%    Ainternal(1:nodes*nodes+1:end) = A(1:nodes*nodes+1:end)';
%     b = zeros(nodes*nodes,1);
    
%     b(1:end) = -gK(1:end).*VK - gNa(1:end).*VNa - .075.*VCl - Vm(1:end,indexPrimary).*Cm./dT;
%     
%     %applying stimulus only on this node
%     %applying stimulus towards the middle right
% %     
% %     b(ceil(nodes*nodes/2 + 5*nodes + 4*nodes/5):ceil(nodes*nodes/2 + 5*nodes + 4*nodes/5)+4) = -gK(ceil(nodes*nodes/2 + 5*nodes + 4*nodes/5):ceil(nodes*nodes/2 + 5*nodes + 4*nodes/5)+4).*VK - gNa(ceil(nodes*nodes/2 + 5*nodes + 4*nodes/5):ceil(nodes*nodes/2 + 5*nodes + 4*nodes/5)+4).*VNa - .075.*VCl - iSt - Vm(ceil(nodes*nodes/2 + 5*nodes + 4*nodes/5):ceil(nodes*nodes/2 + 5*nodes + 4*nodes/5)+4,indexPrimary).*Cm./dT;%-gK(1).*VK - gNa(1).*VNa - .075.*VCl - Vm(1,i).*Cm./dT;
% %     b(ceil(nodes*nodes/2 + 6*nodes + 4*nodes/5):ceil(nodes*nodes/2 + 6*nodes + 4*nodes/5)+4) = -gK(ceil(nodes*nodes/2 + 6*nodes + 4*nodes/5):ceil(nodes*nodes/2 + 6*nodes + 4*nodes/5)+4).*VK - gNa(ceil(nodes*nodes/2 + 6*nodes + 4*nodes/5):ceil(nodes*nodes/2 + 6*nodes + 4*nodes/5)+4).*VNa - .075.*VCl - iSt - Vm(ceil(nodes*nodes/2 + 6*nodes + 4*nodes/5):ceil(nodes*nodes/2 + 6*nodes + 4*nodes/5)+4,indexPrimary).*Cm./dT;%-gK(1).*VK - gNa(1).*VNa - .075.*VCl - Vm(1,i).*Cm./dT;
% %     b(ceil(nodes*nodes/2 + 7*nodes + 4*nodes/5):ceil(nodes*nodes/2 + 7*nodes + 4*nodes/5)+4) = -gK(ceil(nodes*nodes/2 + 7*nodes + 4*nodes/5):ceil(nodes*nodes/2 + 7*nodes + 4*nodes/5)+4).*VK - gNa(ceil(nodes*nodes/2 + 7*nodes + 4*nodes/5):ceil(nodes*nodes/2 + 7*nodes + 4*nodes/5)+4).*VNa - .075.*VCl - iSt - Vm(ceil(nodes*nodes/2 + 7*nodes + 4*nodes/5):ceil(nodes*nodes/2 + 7*nodes + 4*nodes/5)+4,indexPrimary).*Cm./dT;%-gK(1).*VK - gNa(1).*VNa - .075.*VCl - Vm(1,i).*Cm./dT;
% %     b(ceil(nodes*nodes/2 + 8*nodes + 4*nodes/5):ceil(nodes*nodes/2 + 8*nodes + 4*nodes/5)+4) = -gK(ceil(nodes*nodes/2 + 8*nodes + 4*nodes/5):ceil(nodes*nodes/2 + 8*nodes + 4*nodes/5)+4).*VK - gNa(ceil(nodes*nodes/2 + 8*nodes + 4*nodes/5):ceil(nodes*nodes/2 + 8*nodes + 4*nodes/5)+4).*VNa - .075.*VCl - iSt - Vm(ceil(nodes*nodes/2 + 8*nodes + 4*nodes/5):ceil(nodes*nodes/2 + 8*nodes + 4*nodes/5)+4,indexPrimary).*Cm./dT;%-gK(1).*VK - gNa(1).*VNa - .075.*VCl - Vm(1,i).*Cm./dT;
% %     
%     
%     %applying stimulus at the corner
%     b(1:4) = -gK(1:4).*VK - gNa(1:4).*VNa - .075.*VCl - iSt - Vm(1:4,indexPrimary).*Cm./dT;%-gK(1).*VK - gNa(1).*VNa - .075.*VCl - Vm(1,i).*Cm./dT;
%     b(nodes+1:nodes+4) = -gK(nodes+1:nodes+4).*VK - gNa(nodes+1:nodes+4).*VNa - .075.*VCl - iSt - Vm(nodes+1:nodes+4,indexPrimary).*Cm./dT;%-gK(1).*VK - gNa(1).*VNa - .075.*VCl - Vm(1,i).*Cm./dT;
%     b(2*nodes+1:2*nodes+4) = -gK(2*nodes+1:2*nodes+4).*VK - gNa(2*nodes+1:2*nodes+4).*VNa - .075.*VCl - iSt - Vm(2*nodes+1:2*nodes+4,indexPrimary).*Cm./dT;%-gK(1).*VK - gNa(1).*VNa - .075.*VCl - Vm(1,i).*Cm./dT;
%     b(3*nodes+1:3*nodes+4) = -gK(3*nodes+1:3*nodes+4).*VK - gNa(3*nodes+1:3*nodes+4).*VNa - .075.*VCl - iSt - Vm(3*nodes+1:3*nodes+4,indexPrimary).*Cm./dT;%-gK(1).*VK - gNa(1).*VNa - .075.*VCl - Vm(1,i).*Cm./dT;
%     
     %applying stimulus only on this node
%     b(floor(nodes/2 -3):floor(nodes/2 +3)) = -gK(floor(nodes/2 -3):floor(nodes/2 +3)).*VK - gNa(floor(nodes/2 -3):floor(nodes/2 +3)).*VNa - .075.*VCl - iSt - Vm(floor(nodes/2 -3):floor(nodes/2 +3),i).*Cm./dT;%-gK(1).*VK - gNa(1).*VNa - .075.*VCl - Vm(1,i).*Cm./dT;
    
    
    %b(3:end-1) = -gK(3:end-1).*VK - gNa(3:end-1).*VNa - .075.*VCl - Vm(3:end-1,i).*Cm./dT;
    %b(end) = -gK(end).*VK - gNa(end).*VNa - .075.*VCl +iSt - Vm(end,i).*Cm./dT;
    
    %Ax + g = b
    %getting g
    
    

    
    %stimulus applied here: Plane wave
    %from above
    if (t(i) >= 0 && t(i) <= durationOfPulse)
       iSt = 30;%250;
       g(1:end) = (dT./Cm).*(-iK1-iX1-iNa-iS);
       %plane wave above
%        g(1:nodes) = (dT./Cm).*(-iK1(1:nodes)-iX1(1:nodes)-iNa(1:nodes)-iS(1:nodes) + iSt);
%        g(nodes+1:2*nodes) = (dT./Cm).*(-iK1(nodes+1:2*nodes)-iX1(nodes+1:2*nodes)-iNa(nodes+1:2*nodes)-iS(nodes+1:2*nodes) + iSt);
       %plane wave below
%        g(nodes*nodes-nodes+1:nodes*nodes) = (dT./Cm).*(-iK1(nodes*nodes-nodes+1:nodes*nodes)-iX1(nodes*nodes-nodes+1:nodes*nodes)-iNa(nodes*nodes-nodes+1:nodes*nodes)-iS(nodes*nodes-nodes+1:nodes*nodes) + iSt);
%        g(nodes*nodes-2*nodes+1:nodes*nodes-nodes) = (dT./Cm).*(-iK1(nodes*nodes-2*nodes+1:nodes*nodes-nodes)-iX1(nodes*nodes-2*nodes+1:nodes*nodes-nodes)-iNa(nodes*nodes-2*nodes+1:nodes*nodes-nodes)-iS(nodes*nodes-2*nodes+1:nodes*nodes-nodes) + iSt);
       %half plane wave below
%        g(nodes*nodes-nodes+1:nodes*nodes-ceil(nodes/2)) = (dT./Cm).*(-iK1(nodes*nodes-nodes+1:nodes*nodes-ceil(nodes/2))-iX1(nodes*nodes-nodes+1:nodes*nodes-ceil(nodes/2))-iNa(nodes*nodes-nodes+1:nodes*nodes-ceil(nodes/2))-iS(nodes*nodes-nodes+1:nodes*nodes-ceil(nodes/2)) + iSt);%-gK(1).*VK - gNa(1).*VNa - .075.*VCl - Vm(1,i).*Cm./dT;
%        g(nodes*nodes-2*nodes+1:nodes*nodes-nodes-ceil(nodes/2)) = (dT./Cm).*(-iK1(nodes*nodes-2*nodes+1:nodes*nodes-nodes-ceil(nodes/2))-iX1(nodes*nodes-2*nodes+1:nodes*nodes-nodes-ceil(nodes/2))-iNa(nodes*nodes-2*nodes+1:nodes*nodes-nodes-ceil(nodes/2))-iS(nodes*nodes-2*nodes+1:nodes*nodes-nodes-ceil(nodes/2)) + iSt);%-gK(1).*VK - gNa(1).*VNa - .075.*VCl - Vm(1,i).*Cm./dT; 
       %from somewhere in the middle
       g(nodes*nodes-10*nodes+7:nodes*nodes-9*nodes-ceil(nodes/2)) = (dT./Cm).*(-iK1(nodes*nodes-10*nodes+7:nodes*nodes-9*nodes-ceil(nodes/2))-iX1(nodes*nodes-10*nodes+7:nodes*nodes-9*nodes-ceil(nodes/2))-iNa(nodes*nodes-10*nodes+7:nodes*nodes-9*nodes-ceil(nodes/2))-iS(nodes*nodes-10*nodes+7:nodes*nodes-9*nodes-ceil(nodes/2)) + iSt);
       g(nodes*nodes-9*nodes+7:nodes*nodes-8*nodes-ceil(nodes/2)) = (dT./Cm).*(-iK1(nodes*nodes-9*nodes+7:nodes*nodes-8*nodes-ceil(nodes/2))-iX1(nodes*nodes-9*nodes+7:nodes*nodes-8*nodes-ceil(nodes/2))-iNa(nodes*nodes-9*nodes+7:nodes*nodes-8*nodes-ceil(nodes/2))-iS(nodes*nodes-9*nodes+7:nodes*nodes-8*nodes-ceil(nodes/2)) + iSt);
  
%    elseif (t(i) > durationOfPulse)
%        g(1:end) = (dT./Cm).*(-iK1-iX1-iNa-iS);
       
   %For a second pulse:
   elseif (t(i) > durationOfPulse && t(i) <= startOfSecPulse)
       g(1:end) = (dT./Cm).*(-iK1-iX1-iNa-iS);
       
   elseif (t(i) > startOfSecPulse && t(i) <= startOfSecPulse + durationOfPulse)
       %half plane wave below
       g(nodes*nodes-nodes+1:nodes*nodes-ceil(nodes/2)) = (dT./Cm).*(-iK1(nodes*nodes-nodes+1:nodes*nodes-ceil(nodes/2))-iX1(nodes*nodes-nodes+1:nodes*nodes-ceil(nodes/2))-iNa(nodes*nodes-nodes+1:nodes*nodes-ceil(nodes/2))-iS(nodes*nodes-nodes+1:nodes*nodes-ceil(nodes/2)) + iSt);%-gK(1).*VK - gNa(1).*VNa - .075.*VCl - Vm(1,i).*Cm./dT;
       g(nodes*nodes-2*nodes+1:nodes*nodes-nodes-ceil(nodes/2)) = (dT./Cm).*(-iK1(nodes*nodes-2*nodes+1:nodes*nodes-nodes-ceil(nodes/2))-iX1(nodes*nodes-2*nodes+1:nodes*nodes-nodes-ceil(nodes/2))-iNa(nodes*nodes-2*nodes+1:nodes*nodes-nodes-ceil(nodes/2))-iS(nodes*nodes-2*nodes+1:nodes*nodes-nodes-ceil(nodes/2)) + iSt);%-gK(1).*VK - gNa(1).*VNa - .075.*VCl - Vm(1,i).*Cm./dT; 
   
    elseif (t(i) > startOfSecPulse + durationOfPulse)
       g(1:end) = (dT./Cm).*(-iK1-iX1-iNa-iS);
       
   end
    
   
    
    

% %second stimulus to cause reentry
%     g(nodes*5+1:nodes:nodes*25+1) = (dT./Cm).*(gK(nodes*5+1:nodes:nodes*25+1).*VK + gNa(nodes*5+1:nodes:nodes*25+1).*VNa + .075.*VCl + iSt2);%-gK(1).*VK - gNa(1).*VNa - .075.*VCl - Vm(1,i).*Cm./dT;
%     g(nodes*5+2:nodes:nodes*25+2) = (dT./Cm).*(gK(nodes*5+2:nodes:nodes*25+2).*VK + gNa(nodes*5+2:nodes:nodes*25+2).*VNa + .075.*VCl + iSt2);%-gK(1).*VK - gNa(1).*VNa - .075.*VCl - Vm(1,i).*Cm./dT;
%     g(nodes*nodes-nodes+10:-nodes:nodes*nodes-ceil(nodes/2)*nodes+10) = (dT./Cm).*(gK(nodes*nodes-nodes+10:-nodes:nodes*nodes-ceil(nodes/2)*nodes+10).*VK + gNa(nodes*nodes-nodes+10:-nodes:nodes*nodes-ceil(nodes/2)*nodes+10).*VNa + .075.*VCl + iSt2);%-gK(1).*VK - gNa(1).*VNa - .075.*VCl - Vm(1,i).*Cm./dT;
%     g(nodes*nodes-nodes+9:-nodes:nodes*nodes-ceil(nodes/2)*nodes+9) = (dT./Cm).*(gK(nodes*nodes-nodes+9:-nodes:nodes*nodes-ceil(nodes/2)*nodes+9).*VK + gNa(nodes*nodes-nodes+9:-nodes:nodes*nodes-ceil(nodes/2)*nodes+9).*VNa + .075.*VCl + iSt2);%-gK(1).*VK - gNa(1).*VNa - .075.*VCl - Vm(1,i).*Cm./dT;
%     g(nodes*nodes-nodes+8:-nodes:nodes*nodes-ceil(nodes/2)*nodes+8) = (dT./Cm).*(gK(nodes*nodes-nodes+8:-nodes:nodes*nodes-ceil(nodes/2)*nodes+8).*VK + gNa(nodes*nodes-nodes+8:-nodes:nodes*nodes-ceil(nodes/2)*nodes+8).*VNa + .075.*VCl + iSt2);%-gK(1).*VK - gNa(1).*VNa - .075.*VCl - Vm(1,i).*Cm./dT;
%     g(nodes*nodes-nodes+7:-nodes:nodes*nodes-ceil(nodes/2)*nodes+7) = (dT./Cm).*(gK(nodes*nodes-nodes+7:-nodes:nodes*nodes-ceil(nodes/2)*nodes+7).*VK + gNa(nodes*nodes-nodes+7:-nodes:nodes*nodes-ceil(nodes/2)*nodes+7).*VNa + .075.*VCl + iSt2);%-gK(1).*VK - gNa(1).*VNa - .075.*VCl - Vm(1,i).*Cm./dT;
%     g(nodes*nodes-nodes+6:-nodes:nodes*nodes-ceil(nodes/2)*nodes+6) = (dT./Cm).*(gK(nodes*nodes-nodes+6:-nodes:nodes*nodes-ceil(nodes/2)*nodes+6).*VK + gNa(nodes*nodes-nodes+6:-nodes:nodes*nodes-ceil(nodes/2)*nodes+6).*VNa + .075.*VCl + iSt2);%-gK(1).*VK - gNa(1).*VNa - .075.*VCl - Vm(1,i).*Cm./dT;
%     g(nodes*nodes-nodes+5:-nodes:nodes*nodes-ceil(nodes/2)*nodes+5) = (dT./Cm).*(gK(nodes*nodes-nodes+5:-nodes:nodes*nodes-ceil(nodes/2)*nodes+5).*VK + gNa(nodes*nodes-nodes+5:-nodes:nodes*nodes-ceil(nodes/2)*nodes+5).*VNa + .075.*VCl + iSt2);%-gK(1).*VK - gNa(1).*VNa - .075.*VCl - Vm(1,i).*Cm./dT;
%     g(nodes*nodes-nodes+4:-nodes:nodes*nodes-ceil(nodes/2)*nodes+4) = (dT./Cm).*(gK(nodes*nodes-nodes+4:-nodes:nodes*nodes-ceil(nodes/2)*nodes+4).*VK + gNa(nodes*nodes-nodes+4:-nodes:nodes*nodes-ceil(nodes/2)*nodes+4).*VNa + .075.*VCl + iSt2);%-gK(1).*VK - gNa(1).*VNa - .075.*VCl - Vm(1,i).*Cm./dT;
%     g(nodes*nodes-nodes+3:-nodes:nodes*nodes-ceil(nodes/2)*nodes+3) = (dT./Cm).*(gK(nodes*nodes-nodes+3:-nodes:nodes*nodes-ceil(nodes/2)*nodes+3).*VK + gNa(nodes*nodes-nodes+3:-nodes:nodes*nodes-ceil(nodes/2)*nodes+3).*VNa + .075.*VCl + iSt2);%-gK(1).*VK - gNa(1).*VNa - .075.*VCl - Vm(1,i).*Cm./dT;
%     g(nodes*nodes-nodes+2:-nodes:nodes*nodes-ceil(nodes/2)*nodes+2) = (dT./Cm).*(gK(nodes*nodes-nodes+2:-nodes:nodes*nodes-ceil(nodes/2)*nodes+2).*VK + gNa(nodes*nodes-nodes+2:-nodes:nodes*nodes-ceil(nodes/2)*nodes+2).*VNa + .075.*VCl + iSt2);%-gK(1).*VK - gNa(1).*VNa - .075.*VCl - Vm(1,i).*Cm./dT;
%     g(nodes*nodes-nodes+1:-nodes:nodes*nodes-ceil(nodes/2)*nodes+1) = (dT./Cm).*(gK(nodes*nodes-nodes+1:-nodes:nodes*nodes-ceil(nodes/2)*nodes+1).*VK + gNa(nodes*nodes-nodes+1:-nodes:nodes*nodes-ceil(nodes/2)*nodes+1).*VNa + .075.*VCl + iSt2);%-gK(1).*VK - gNa(1).*VNa - .075.*VCl - Vm(1,i).*Cm./dT;
%     
    



    %start of matrix creation for ELECTRICAL PROPAGATION
    %dead tissue only present in some inner nodes, not corners or edges:
    %resistance values of "infinity"
   
    %UPDATING A DUE TO CHANGING DX AND DY FOR INNER AND EDGE NODES (NO
    %CORNERS)
    
        for j = 2:1:nodes-1
            for k = 2:1:nodes-1
                Epoint = j*nodes - (nodes-k);   
                if rem(j,2) == 0 %this is for rows WITH NO mass points on them
                    if rem(k,2) == 0%this is for columns BETWEEN mass points
                        Mpoint = ((Epoint-nodes-1) - nodes*((j-2)/2)+ ((j-2)/2 + 1))/2;
                        dXloc1 = (abs( ((massX(Mpoint,1) + massX(Mpoint+((nodes+1)./2),1))./2) - ((massX(Mpoint+1,1) + massX(Mpoint+1+((nodes+1)./2),1))./2) )./2);
                        dXloc2 = (abs( ((massX(Mpoint,1) + massX(Mpoint+((nodes+1)./2),1))./2) - ((massX(Mpoint+1,1) + massX(Mpoint+1+((nodes+1)./2),1))./2) )./2);
                        dYloc1 = (abs( ((massY(Mpoint,1) + massY(Mpoint+1,1))./2) - ((massY(Mpoint+((nodes+1)/2),1) + massY(Mpoint+((nodes+1)/2)+1,1))./2) )./2);
                        dYloc2 = (abs( ((massY(Mpoint,1) + massY(Mpoint+1,1))./2) - ((massY(Mpoint+((nodes+1)/2),1) + massY(Mpoint+((nodes+1)/2)+1,1))./2) )./2); 
                    else%this is for columns WITH mass points
                        Mpoint = ((Epoint-nodes) - nodes*((j-2)/2)+ ((j-2)/2 + 1))/2;
                        dXloc1 = (abs( ((massX(Mpoint,1) + massX(Mpoint+((nodes+1)./2),1))./2) - ((massX(Mpoint-1,1) + massX(Mpoint-1+((nodes+1)./2),1))./2) )./2);
                        dXloc2 = (abs( ((massX(Mpoint,1) + massX(Mpoint+((nodes+1)./2),1))./2) - ((massX(Mpoint+1,1) + massX(Mpoint+1+((nodes+1)./2),1))./2) )./2);
                        dYloc1 = (abs(massY(Mpoint,1) - massY(Mpoint+((nodes+1)/2),1))./2);
                        dYloc2 = (abs(massY(Mpoint,1) - massY(Mpoint+((nodes+1)/2),1))./2);
                    end
                else %this is for rows WITH mass points on them 
                    if rem(k,2) == 0%this is for columns BETWEEN mass points
                        Mpoint = ((Epoint-1) - nodes*((j-1)/2)+ ((j-1)/2 + 1))/2;
                        dXloc1 = (abs(massX(Mpoint,1) - massX(Mpoint+1,1))./2);
                        dXloc2 = (abs(massX(Mpoint,1) - massX(Mpoint+1,1))./2);
                        dYloc1 = (abs(massY(Mpoint,1) - ((massY(Mpoint-((nodes+1)/2),1) + massY(Mpoint-((nodes+1)/2) + 1,1))./2) )./2);
                        dYloc2 = (abs(massY(Mpoint,1) - ((massY(Mpoint+((nodes+1)/2),1) + massY(Mpoint+((nodes+1)/2) + 1,1))./2) )./2); 
                    else%this is for columns WITH mass points
                        Mpoint = ((Epoint) - nodes*((j-1)/2)+ ((j-1)/2 + 1))/2;
                        dXloc1 = (abs(massX(Mpoint,1) - massX(Mpoint-1,1))./2);
                        dXloc2 = (abs(massX(Mpoint,1) - massX(Mpoint+1,1))./2);
                        dYloc1 = (abs(massY(Mpoint,1) - massY(Mpoint-((nodes+1)/2),1))./2);
                        dYloc2 = (abs(massY(Mpoint,1) - massY(Mpoint+((nodes+1)/2),1))./2);
                    end
                end
                
                %if statements creating the holes in the middle of the conductive
        %tissue
%         if (j > ceil(nodes*nodes/2 - 3*nodes + nodes/5)) && (j < ceil(nodes*nodes/2 - 3*nodes + 4*nodes/5))
%             blockageValue = 10000000000;
%         elseif (j > ceil(nodes*nodes/2 - 2*nodes + nodes/5)) && (j < ceil(nodes*nodes/2 - 2*nodes + 4*nodes/5))
%             blockageValue = 10000000000;
%         elseif (j > ceil(nodes*nodes/2 - 1*nodes + nodes/5)) && (j < ceil(nodes*nodes/2 - 1*nodes + 4*nodes/5))
%             blockageValue = 10000000000;
%         elseif (j > ceil(nodes*nodes/2 + nodes/5)) && (j < ceil(nodes*nodes/2 + 4*nodes/5))
%             blockageValue = 10000000000;
%         elseif (j > ceil(nodes*nodes/2 + 1*nodes + nodes/5)) && (j < ceil(nodes*nodes/2 + 1*nodes + 4*nodes/5))
%             blockageValue = 10000000000;
%         elseif (j > ceil(nodes*nodes/2 + 2*nodes + nodes/5)) && (j < ceil(nodes*nodes/2 + 2*nodes + 4*nodes/5))
%             blockageValue = 10000000000;
%         elseif (j > ceil(nodes*nodes/2 + 3*nodes + nodes/5)) && (j < ceil(nodes*nodes/2 + 3*nodes + 4*nodes/5))
%             blockageValue = 10000000000;
%         else
%             blockageValue = 1;
%         end
        blockageValue = 1; %no blockage
        A(Epoint,Epoint-nodes) = dT./(blockageValue * rey.*2.*pi.*a.*dYloc1.*Cm); %in y
        A(Epoint,Epoint-1) = dT./(blockageValue * rex.*2.*pi.*a.*dXloc1.*Cm); %in x
        A(Epoint,Epoint) = (dT./Cm).*( -1./(blockageValue * rex.*2.*pi.*a.*dXloc1) -1./(blockageValue * rex.*2.*pi.*a.*dXloc2) -1./(blockageValue * rey.*2.*pi.*a.*dYloc1) -1./(blockageValue * rey.*2.*pi.*a.*dYloc2) + Cm./dT);
        A(Epoint,Epoint+1) = dT./(blockageValue * rex.*2.*pi.*a.*dXloc2.*Cm); %in x
        A(Epoint,Epoint+nodes) = dT./(blockageValue * rey.*2.*pi.*a.*dYloc2.*Cm); %in y
            
            end
        end
    
    %rows for edges (without corners): for loops; make the whole row zero first and then add
    %the points

    %top edge:
    Mpoint = 1;
    for k=2:1:nodes-1
        
        if rem(k,2) == 0
            dXloc1 = (abs(massX(Mpoint,1) - massX(Mpoint+1,1))./2);
            dXloc2 = (abs(massX(Mpoint,1) - massX(Mpoint+1,1))./2);
            dYloc2 = (abs( ((massY(Mpoint,1) + massY(Mpoint+1,1))./2) - ((massY(Mpoint+((nodes+1)/2),1) + massY(Mpoint+((nodes+1)/2)+1,1))./2) )./2);
            Mpoint = Mpoint+1;
        else
            dXloc1 = (abs(massX(Mpoint,1) - massX(Mpoint-1,1))./2);
            dXloc2 = (abs(massX(Mpoint,1) - massX(Mpoint+1,1))./2);
            dYloc2 = (abs(massY(Mpoint,1) - massY(Mpoint+((nodes+1)/2),1))./2);
            Mpoint = Mpoint+1;
        end
        
        A(k,:) = 0;
        A(k,k-1) = (dT./Cm).*(1./(rex.*2.*pi.*a.*dXloc1));
        A(k,k) = (dT./Cm).*( -1./(rex.*2.*pi.*a.*dXloc1) -1./(rex.*2.*pi.*a.*dXloc2) -1./(rey.*2.*pi.*a.*dYloc2) + Cm./dT);
        A(k,k+1) = (dT./Cm).*(1./(rex.*2.*pi.*a.*dXloc2));
        A(k,k+nodes) = (dT./Cm).*(1./(rey.*2.*pi.*a.*dYloc2));
        
    end
    
    %bottom edge:
    Mpoint = ((nodes+1)./2).^2 - (nodes+1)./2 + 1;
    for k=nodes*nodes-nodes+2:1:nodes*nodes-1
        
        if rem(k,2) == 0
            dXloc1 = (abs(massX(Mpoint,1) - massX(Mpoint+1,1))./2);
            dXloc2 = (abs(massX(Mpoint,1) - massX(Mpoint+1,1))./2);
            dYloc1 = (abs( ((massY(Mpoint,1) + massY(Mpoint+1,1))./2) - ((massY(Mpoint-((nodes+1)/2),1) + massY(Mpoint-((nodes+1)/2)+1,1))./2) )./2);
            Mpoint = Mpoint+1;
        else
            dXloc1 = (abs(massX(Mpoint,1) - massX(Mpoint-1,1))./2);
            dXloc2 = (abs(massX(Mpoint,1) - massX(Mpoint+1,1))./2);
            dYloc1 = (abs(massY(Mpoint,1) - massY(Mpoint-((nodes+1)/2),1))./2);
            Mpoint = Mpoint+1;
        end
        
        A(k,:) = 0;
        A(k,k-1) = (dT./Cm).*(1./(rex.*2.*pi.*a.*dXloc1));
        A(k,k) = (dT./Cm).*( -1./(rex.*2.*pi.*a.*dXloc1) -1./(rex.*2.*pi.*a.*dXloc2) -1./(rey.*2.*pi.*a.*dYloc1) + Cm./dT);
        A(k,k+1) = (dT./Cm).*(1./(rex.*2.*pi.*a.*dXloc2));
        A(k,k-nodes) = (dT./Cm).*(1./(rey.*2.*pi.*a.*dYloc1));
        
    end
    
    %right edge:
    Mpoint = ((nodes+1)./2);
    for k= 2*nodes :nodes:nodes*nodes-nodes
        
        if rem(k,2) == 0
            dXloc1 = (abs( ((massX(Mpoint,1) + massX(Mpoint+((nodes+1)./2),1))./2) - ((massX(Mpoint-1,1) + massX(Mpoint-1+((nodes+1)./2),1))./2) )./2);
            dYloc1 = (abs(massY(Mpoint,1) - massY(Mpoint+((nodes+1)/2),1))./2);
            dYloc2 = (abs(massY(Mpoint,1) - massY(Mpoint+((nodes+1)/2),1))./2);
            Mpoint = Mpoint+((nodes+1)./2);
        else
            dXloc1 = (abs(massX(Mpoint,1) - massX(Mpoint-1,1))./2);
            dYloc1 = (abs(massY(Mpoint,1) - massY(Mpoint-((nodes+1)/2),1))./2);
            dYloc2 = (abs(massY(Mpoint,1) - massY(Mpoint+((nodes+1)/2),1))./2);
            Mpoint = Mpoint+((nodes+1)./2);
        end
        
        A(k,:) = 0;
        A(k,k-nodes) = (dT./Cm).*(1./(rey.*2.*pi.*a.*dYloc1));
        A(k,k) = (dT./Cm).*( -1./(rey.*2.*pi.*a.*dYloc1) -1./(rey.*2.*pi.*a.*dYloc2) -1./(rex.*2.*pi.*a.*dXloc1) + Cm./dT);
        A(k,k+nodes) = (dT./Cm).*(1./(rey.*2.*pi.*a.*dYloc2));
        A(k,k-1) = (dT./Cm).*(1./(rex.*2.*pi.*a.*dXloc1));
        
    end
    
    %left edge:
    Mpoint = 1;
    for k= nodes+1 :nodes:nodes*nodes-2*nodes+1
        
        if rem(k,2) == 0
            dXloc2 = (abs( ((massX(Mpoint,1) + massX(Mpoint+((nodes+1)./2),1))./2) - ((massX(Mpoint+1,1) + massX(Mpoint+1+((nodes+1)./2),1))./2) )./2);
            dYloc1 = (abs(massY(Mpoint,1) - massY(Mpoint+((nodes+1)/2),1))./2);
            dYloc2 = (abs(massY(Mpoint,1) - massY(Mpoint+((nodes+1)/2),1))./2);
            Mpoint = Mpoint+((nodes+1)./2);
        else
            dXloc2 = (abs(massX(Mpoint,1) - massX(Mpoint+1,1))./2);
            dYloc1 = (abs(massY(Mpoint,1) - massY(Mpoint-((nodes+1)/2),1))./2);
            dYloc2 = (abs(massY(Mpoint,1) - massY(Mpoint+((nodes+1)/2),1))./2);
            Mpoint = Mpoint+((nodes+1)./2);
        end
        
        A(k,:) = 0;
        A(k,k-nodes) = (dT./Cm).*(1./(rey.*2.*pi.*a.*dYloc1));
        A(k,k) = (dT./Cm).*( -1./(rey.*2.*pi.*a.*dYloc1) -1./(rey.*2.*pi.*a.*dYloc2) -1./(rex.*2.*pi.*a.*dXloc2) + Cm./dT);
        A(k,k+nodes) = (dT./Cm).*(1./(rey.*2.*pi.*a.*dYloc2));
        A(k,k+1) = (dT./Cm).*(1./(rex.*2.*pi.*a.*dXloc2));
    
    end
       
    
    %rows for corner points
    %These ones don't deform because of fixed boundary conditions therefore
    %no update on dX and dY needed.
    
    %first corner point 1 left upper
    A(1,1) = (dT./Cm).*(-1./(rex.*2.*pi.*a.*dX) -1./(rey.*2.*pi.*a.*dY) + Cm./dT);%-2./(ri.*2.*pi.*a.*dX) - gK(1) -gNa(1) - Cm./dT;
    A(1,2) = (dT./Cm).*(1./(rex.*2.*pi.*a.*dX));
    A(1,1+nodes) = (dT./Cm).*(1./(rey.*2.*pi.*a.*dY));
    %second corner point nodes; right upper
    A(nodes,nodes) = (dT./Cm).*(-1./(rex.*2.*pi.*a.*dX) -1./(rey.*2.*pi.*a.*dY) + Cm./dT);%-2./(ri.*2.*pi.*a.*dX) - gK(1) -gNa(1) - Cm./dT;
    A(nodes,nodes-1) = (dT./Cm).*(1./(rex.*2.*pi.*a.*dX));
    A(nodes,nodes+nodes) = (dT./Cm).*(1./(rey.*2.*pi.*a.*dY));
    %third corner point nodes*nodes-nodes+1 left lower
    A(end - nodes + 1,end - nodes + 1) = (dT./Cm).*(-1./(rex.*2.*pi.*a.*dX) -1./(rey.*2.*pi.*a.*dY) + Cm./dT);%-2./(ri.*2.*pi.*a.*dX) - gK(1) -gNa(1) - Cm./dT;
    A(end - nodes + 1,end - nodes + 2) = (dT./Cm).*(1./(rex.*2.*pi.*a.*dX));
    A(end - nodes + 1,end - 2*nodes + 1) = (dT./Cm).*(1./(rey.*2.*pi.*a.*dY));
    %4th corner point nodes*nodes right lower
    A(end,end) = (dT./Cm).*(-1./(rex.*2.*pi.*a.*dX) -1./(rey.*2.*pi.*a.*dY) + Cm./dT);%-2./(ri.*2.*pi.*a.*dX) - gK(end) -gNa(end) - Cm./dT;
    A(end,end - 1) = (dT./Cm).*(1./(rex.*2.*pi.*a.*dX));
    A(end,end - nodes) = (dT./Cm).*(1./(rey.*2.*pi.*a.*dY));

    Ainternal = A;

    b = Ainternal*Vm(:,indexPrimary) + g;
    Vm(:,indexPrimary+1) = b;
    Vmnondim(:,indexPrimary+1) = (Vm(:,indexPrimary+1)-minVal)./(maxVal-minVal);
    
    
    
    %end of matrix operations
    
    
    % Start of mechanical part
    %
    %
    
    %Obtaining values for Ta from the Vm array
    for indrowm=1:((nodes+1)/2)
        for indnodes=1:2:nodes
            
            indcole = (indrowm*2-1)*nodes - (nodes-indnodes);
            if (Vmnondim(indcole,indexPrimary) < aThresh)
                epsi(((indcole+indrowm-nodes*(indrowm-1))/2),1) = 1;
            elseif (Vmnondim(indcole,indexPrimary) >= aThresh)
                epsi(((indcole+indrowm-nodes*(indrowm-1))/2),1) = .1;
            end
            
            
            
            Ta(((indcole+indrowm-nodes*(indrowm-1))/2),indexPrimary+1) = Ta(((indcole+indrowm-nodes*(indrowm-1))/2),indexPrimary) + (dT.*(epsi(((indcole+indrowm-nodes*(indrowm-1))/2),1).*(kT.*Vmnondim(indcole,indexPrimary) - Ta(((indcole+indrowm-nodes*(indrowm-1))/2),indexPrimary))));
            %disp(Ta(((indcole+indrowm-nodes*(indrowm-1))/2),i+1))
        end
    end
    
    temporaryMassX = zeros(((nodes+1)/2).^2,3);
    temporaryMassX(:,1) = massX(:,1);
    temporaryMassY = zeros(((nodes+1)/2).^2,3);
    temporaryMassY(:,1) = massY(:,1);
    
    %same matrices for x and y
    for j = (((nodes+1)/2)+2) : 1 : (((nodes+1)/2)^2 - ((nodes+1)/2) - 1)
        
        Axyddn(j,j-((nodes+1)/2)-1) = kStiffRat.*c_mech/(sqrt(2).*r0) - kStiffRat.*d_mech/dT2 + kStiffRat.*((-c_mech) ./ (((temporaryMassX(j,1) - temporaryMassX(j-((nodes+1)/2)-1,1)).^2 + (temporaryMassY(j,1) - temporaryMassY(j-((nodes+1)/2)-1,1)).^2).^(1./2)));
        Axyddn(j,j-((nodes+1)/2)) = c_mech/(r0) - d_mech/dT2 + ((-c_mech + .5.*(scalingTa.*Ta(j,2)) + .5.*(scalingTa.*Ta(j-((nodes+1)/2),2))) ./ (((temporaryMassX(j,1) - temporaryMassX(j-((nodes+1)/2),1)).^2 + (temporaryMassY(j,1) - temporaryMassY(j-((nodes+1)/2),1)).^2).^(1./2)));
        Axyddn(j,j-((nodes+1)/2)+1) = kStiffRat.*c_mech/(sqrt(2).*r0) - kStiffRat.*d_mech/dT2 + kStiffRat.*((-c_mech) ./ (((temporaryMassX(j,1) - temporaryMassX(j-((nodes+1)/2)+1,1)).^2 + (temporaryMassY(j,1) - temporaryMassY(j-((nodes+1)/2)+1,1)).^2).^(1./2)));
        
        Axyddn(j,j-1)             = c_mech/r0 - d_mech/dT2 + ((-c_mech + .5.*(scalingTa.*Ta(j,2)) + .5.*(scalingTa.*Ta(j-1,2))) ./ (((temporaryMassX(j,1) - temporaryMassX(j-1,1)).^2 + (temporaryMassY(j,1) - temporaryMassY(j-1,1)).^2).^(1./2)));
        Axyddn(j,j)               = -4*c_mech/r0 - kStiffRat*4*c_mech/(sqrt(2)*r0) + 4*d_mech/dT2 + kStiffRat*4*d_mech/dT2 - ((-c_mech + .5.*(scalingTa.*Ta(j,2)) + .5.*(scalingTa.*Ta(j+((nodes+1)/2),2))) ./ (((temporaryMassX(j,1) - temporaryMassX(j+((nodes+1)/2),1)).^2 + (temporaryMassY(j,1) - temporaryMassY(j+((nodes+1)/2),1)).^2).^(1./2))) - ((-c_mech + .5.*(scalingTa.*Ta(j,2)) + .5.*(scalingTa.*Ta(j+1,2))) ./ (((temporaryMassX(j,1) - temporaryMassX(j+1,1)).^2 + (temporaryMassY(j,1) - temporaryMassY(j+1,1)).^2).^(1./2))) - ((-c_mech + .5.*(scalingTa.*Ta(j,2)) + .5.*(scalingTa.*Ta(j-1,2))) ./ (((temporaryMassX(j,1) - temporaryMassX(j-1,1)).^2 + (temporaryMassY(j,1) - temporaryMassY(j-1,1)).^2).^(1./2))) - ((-c_mech + .5.*(scalingTa.*Ta(j,2)) + .5.*(scalingTa.*Ta(j-((nodes+1)/2),2))) ./ (((temporaryMassX(j,1) - temporaryMassX(j-((nodes+1)/2),1)).^2 + (temporaryMassY(j,1) - temporaryMassY(j-((nodes+1)/2),1)).^2).^(1./2))) - (( -c_mech * kStiffRat ) ./ (((temporaryMassX(j,1) - temporaryMassX(j-((nodes+1)/2) - 1,1)).^2 + (temporaryMassY(j,1) - temporaryMassY(j-((nodes+1)/2)-1 ,1)).^2).^(1./2))) - (( -c_mech * kStiffRat ) ./ (((temporaryMassX(j,1) - temporaryMassX(j-((nodes+1)/2) + 1,1)).^2 + (temporaryMassY(j,1) - temporaryMassY(j-((nodes+1)/2) + 1 ,1)).^2).^(1./2))) - (( -c_mech * kStiffRat ) ./ (((temporaryMassX(j,1) - temporaryMassX(j+((nodes+1)/2) - 1,1)).^2 + (temporaryMassY(j,1) - temporaryMassY(j+((nodes+1)/2)-1 ,1)).^2).^(1./2))) - (( -c_mech * kStiffRat ) ./ (((temporaryMassX(j,1) - temporaryMassX(j+((nodes+1)/2) + 1,1)).^2 + (temporaryMassY(j,1) - temporaryMassY(j+((nodes+1)/2)+1 ,1)).^2).^(1./2)));
        Axyddn(j,j+1)             = c_mech/r0 - d_mech/dT2 + ((-c_mech + .5.*(scalingTa.*Ta(j,2)) + .5.*(scalingTa.*Ta(j+1,2))) ./ (((temporaryMassX(j,1) - temporaryMassX(j+1,1)).^2 + (temporaryMassY(j,1) - temporaryMassY(j+1,1)).^2).^(1./2)));
        
        Axyddn(j,j+((nodes+1)/2)-1) = kStiffRat.*c_mech/(sqrt(2).*r0) - kStiffRat.*d_mech/dT2 + kStiffRat.*((-c_mech) ./ (((temporaryMassX(j,1) - temporaryMassX(j+((nodes+1)/2)-1,1)).^2 + (temporaryMassY(j,1) - temporaryMassY(j+((nodes+1)/2)-1,1)).^2).^(1./2)));
        Axyddn(j,j+((nodes+1)/2)) = c_mech/r0 - d_mech/dT2 + ((-c_mech + .5.*(scalingTa.*Ta(j,2)) + .5.*(scalingTa.*Ta(j+((nodes+1)/2),2))) ./ (((temporaryMassX(j,1) - temporaryMassX(j+((nodes+1)/2),1)).^2 + (temporaryMassY(j,1) - temporaryMassY(j+((nodes+1)/2),1)).^2).^(1./2)));
        Axyddn(j,j+((nodes+1)/2)+1) = kStiffRat.*c_mech/(sqrt(2).*r0) - kStiffRat.*d_mech/dT2 + kStiffRat.*((-c_mech) ./ (((temporaryMassX(j,1) - temporaryMassX(j+((nodes+1)/2)+1,1)).^2 + (temporaryMassY(j,1) - temporaryMassY(j+((nodes+1)/2)+1,1)).^2).^(1./2)));
        
    end
    
    massXdd(:,1) = Axyddn*temporaryMassX(:,1);
    massYdd(:,1) = Axyddn*temporaryMassY(:,1);
    
    sumat = sum(sqrt((massXdd(:,1)).^2 + (massYdd(:,1)).^2));%sum(massXdd(:,1)) + sum(massYdd(:,1));
    
    temporaryMassX(:,2) = temporaryMassX(:,1) + .5.*(dT2.^2).*(massXdd(:,1));
    temporaryMassY(:,2) = temporaryMassY(:,1) + .5.*(dT2.^2).*(massYdd(:,1));
    
    %setting back the boundaries to be fixed
    temporaryMassX(1:(nodes+1)/2:end,2) = leftValuex;
    temporaryMassX((nodes+1)/2:(nodes+1)/2:end,2) = rightValuex;
    temporaryMassX(1:1:(nodes+1)/2,2) = massX(1:1:(nodes+1)/2,1);
    temporaryMassX(((nodes+1)/2)^2-((nodes+1)/2)+1:1:((nodes+1)/2)^2,2) = massX(((nodes+1)/2)^2-((nodes+1)/2)+1:1:((nodes+1)/2)^2,1);
    temporaryMassY(1:(nodes+1)/2:end,2) = massY(1:(nodes+1)/2:end,1);
    temporaryMassY((nodes+1)/2:(nodes+1)/2:end,2) = massY((nodes+1)/2:(nodes+1)/2:end,1);
    temporaryMassY(1:1:(nodes+1)/2,2) = upperValuey;
    temporaryMassY(((nodes+1)/2)^2-((nodes+1)/2)+1:1:((nodes+1)/2)^2,2) = bottomValuey;
    
    %iterative process for the convervence of the new (x,y) setup
    
    VerT = 2;
    markerWhile = 0;
    threshPerc = 1;
    
       
       while round(sumat,5) > accThreshTot*threshPerc %|| sumat <= -accThreshTot
              
           for j = ((nodes+1)/2+2) : 1 : (((nodes+1)/2)^2 - ((nodes+1)/2) - 1)
       
                Axyddn(j,j-((nodes+1)/2)-1) = kStiffRat.*c_mech/(sqrt(2).*r0) - kStiffRat.*d_mech/dT2 + kStiffRat.*((-c_mech) ./ (((temporaryMassX(j,VerT) - temporaryMassX(j-((nodes+1)/2)-1,VerT)).^2 + (temporaryMassY(j,VerT) - temporaryMassY(j-((nodes+1)/2)-1,VerT)).^2).^(1./2)));
                Axyddn(j,j-((nodes+1)/2)) =   c_mech/r0 - d_mech/dT2 + ((-c_mech + .5.*(scalingTa.*Ta(j,2)) + .5.*(scalingTa.*Ta(j-((nodes+1)/2),2))) ./ (((temporaryMassX(j,VerT) - temporaryMassX(j-((nodes+1)/2),VerT)).^2 + (temporaryMassY(j,VerT) - temporaryMassY(j-((nodes+1)/2),VerT)).^2).^(1./2)));
                Axyddn(j,j-((nodes+1)/2)+1) = kStiffRat.*c_mech/(sqrt(2).*r0) - kStiffRat.*d_mech/dT2 + kStiffRat.*((-c_mech) ./ (((temporaryMassX(j,VerT) - temporaryMassX(j-((nodes+1)/2)+1,VerT)).^2 + (temporaryMassY(j,VerT) - temporaryMassY(j-((nodes+1)/2)+1,VerT)).^2).^(1./2)));
                
                Axyddn(j,j-1)             = c_mech/r0 - d_mech/dT2 + ((-c_mech + .5.*(scalingTa.*Ta(j,2)) + .5.*(scalingTa.*Ta(j-1,2))) ./ (((temporaryMassX(j,VerT) - temporaryMassX(j-1,VerT)).^2 + (temporaryMassY(j,VerT) - temporaryMassY(j-1,VerT)).^2).^(1./2)));
                Axyddn(j,j)               = -4*c_mech/r0 - kStiffRat*4*c_mech/(sqrt(2)*r0) + 4*d_mech/dT2 + kStiffRat*4*d_mech/dT2 - ((-c_mech + .5.*(scalingTa.*Ta(j,2)) + .5.*(scalingTa.*Ta(j+((nodes+1)/2),2))) ./ (((temporaryMassX(j,VerT) - temporaryMassX(j+((nodes+1)/2),VerT)).^2 + (temporaryMassY(j,VerT) - temporaryMassY(j+((nodes+1)/2),VerT)).^2).^(1./2))) - ((-c_mech + .5.*(scalingTa.*Ta(j,2)) + .5.*(scalingTa.*Ta(j+1,2))) ./ (((temporaryMassX(j,VerT) - temporaryMassX(j+1,VerT)).^2 + (temporaryMassY(j,VerT) - temporaryMassY(j+1,VerT)).^2).^(1./2))) - ((-c_mech + .5.*(scalingTa.*Ta(j,2)) + .5.*(scalingTa.*Ta(j-1,2))) ./ (((temporaryMassX(j,VerT) - temporaryMassX(j-1,VerT)).^2 + (temporaryMassY(j,VerT) - temporaryMassY(j-1,VerT)).^2).^(1./2))) - ((-c_mech + .5.*(scalingTa.*Ta(j,2)) + .5.*(scalingTa.*Ta(j-((nodes+1)/2),2))) ./ (((temporaryMassX(j,VerT) - temporaryMassX(j-((nodes+1)/2),VerT)).^2 + (temporaryMassY(j,VerT) - temporaryMassY(j-((nodes+1)/2),VerT)).^2).^(1./2))) - (( -c_mech * kStiffRat ) ./ (((temporaryMassX(j,VerT) - temporaryMassX(j-((nodes+1)/2) - 1,VerT)).^2 + (temporaryMassY(j,VerT) - temporaryMassY(j-((nodes+1)/2)-1 ,VerT)).^2).^(1./2))) - (( -c_mech * kStiffRat ) ./ (((temporaryMassX(j,VerT) - temporaryMassX(j-((nodes+1)/2) + 1,VerT)).^2 + (temporaryMassY(j,VerT) - temporaryMassY(j-((nodes+1)/2) + 1 ,VerT)).^2).^(1./2))) - (( -c_mech * kStiffRat ) ./ (((temporaryMassX(j,VerT) - temporaryMassX(j+((nodes+1)/2) - 1,VerT)).^2 + (temporaryMassY(j,VerT) - temporaryMassY(j+((nodes+1)/2)-1 ,VerT)).^2).^(1./2))) - (( -c_mech * kStiffRat ) ./ (((temporaryMassX(j,VerT) - temporaryMassX(j+((nodes+1)/2) + 1,VerT)).^2 + (temporaryMassY(j,VerT) - temporaryMassY(j+((nodes+1)/2)+1 ,VerT)).^2).^(1./2)));
                Axyddn(j,j+1)             = c_mech/r0 - d_mech/dT2 + ((-c_mech + .5.*(scalingTa.*Ta(j,2)) + .5.*(scalingTa.*Ta(j+1,2))) ./ (((temporaryMassX(j,VerT) - temporaryMassX(j+1,VerT)).^2 + (temporaryMassY(j,VerT) - temporaryMassY(j+1,VerT)).^2).^(1./2)));
                
                Axyddn(j,j+((nodes+1)/2)-1) = kStiffRat.*c_mech/(sqrt(2).*r0) - kStiffRat.*d_mech/dT2 + kStiffRat.*((-c_mech) ./ (((temporaryMassX(j,VerT) - temporaryMassX(j+((nodes+1)/2)-1,VerT)).^2 + (temporaryMassY(j,VerT) - temporaryMassY(j+((nodes+1)/2)-1,VerT)).^2).^(1./2)));
                Axyddn(j,j+((nodes+1)/2)) = c_mech/r0 - d_mech/dT2 + ((-c_mech + .5.*(scalingTa.*Ta(j,2)) + .5.*(scalingTa.*Ta(j+((nodes+1)/2),2))) ./ (((temporaryMassX(j,VerT) - temporaryMassX(j+((nodes+1)/2),VerT)).^2 + (temporaryMassY(j,VerT) - temporaryMassY(j+((nodes+1)/2),VerT)).^2).^(1./2)));
                Axyddn(j,j+((nodes+1)/2)+1) = kStiffRat.*c_mech/(sqrt(2).*r0) - kStiffRat.*d_mech/dT2 + kStiffRat.*((-c_mech) ./ (((temporaryMassX(j,VerT) - temporaryMassX(j+((nodes+1)/2)+1,VerT)).^2 + (temporaryMassY(j,VerT) - temporaryMassY(j+((nodes+1)/2)+1,VerT)).^2).^(1./2)));
                
           end
           
           massXdd(:,1) = Axyddn*temporaryMassX(:,VerT) + Axyddnprev*temporaryMassX(:,VerT-1);
           temporaryMassX(:,VerT+1) = 2.*temporaryMassX(:,VerT) - temporaryMassX(:,VerT-1) + (dT2.^2).*(massXdd(:,1));
           massYdd(:,1) = Axyddn*temporaryMassY(:,VerT) + Axyddnprev*temporaryMassY(:,VerT-1);
           temporaryMassY(:,VerT+1) = 2.*temporaryMassY(:,VerT) - temporaryMassY(:,VerT-1) + (dT2.^2).*(massYdd(:,1));
           
           %setting back the boundaries to be fixed
            temporaryMassX(1:(nodes+1)/2:end,VerT+1) = leftValuex;
            temporaryMassX((nodes+1)/2:(nodes+1)/2:end,VerT+1) = rightValuex;
            temporaryMassX(1:1:(nodes+1)/2,VerT+1) = massX(1:1:(nodes+1)/2,1);
            temporaryMassX(((nodes+1)/2)^2-((nodes+1)/2)+1:1:((nodes+1)/2)^2,VerT+1) = massX(((nodes+1)/2)^2-((nodes+1)/2)+1:1:((nodes+1)/2)^2,1);
            temporaryMassY(1:(nodes+1)/2:end,VerT+1) = massY(1:(nodes+1)/2:end,1);
            temporaryMassY((nodes+1)/2:(nodes+1)/2:end,VerT+1) = massY((nodes+1)/2:(nodes+1)/2:end,1);
            temporaryMassY(1:1:(nodes+1)/2,VerT+1) = upperValuey;
            temporaryMassY(((nodes+1)/2)^2-((nodes+1)/2)+1:1:((nodes+1)/2)^2,VerT+1) = bottomValuey;
   
 
           
           sumat = sum(sqrt((massXdd(:,1)).^2 + (massYdd(:,1)).^2));%sum(massXdd(:,1)) + sum(massYdd(:,1));
           
           temporaryMassX(:,VerT-1) = temporaryMassX(:,VerT);
           temporaryMassX(:,VerT) = temporaryMassX(:,VerT+1);
           temporaryMassY(:,VerT-1) = temporaryMassY(:,VerT);
           temporaryMassY(:,VerT) = temporaryMassY(:,VerT+1);
           
           markerWhile = markerWhile + 1;
           
           %adjusting the threshold: max number of iterations nodes*nodes
           if markerWhile > 1.5e3 && markerWhile <= 2.5e3
               threshPerc = 1.1;
%                disp(['index point: ',num2str(indexSecondary),'; convergence iteration: ',num2str(markerWhile)])
           elseif markerWhile > 2.5e3 && markerWhile <= 3e3
               threshPerc = 1.15;
%                disp(['index point: ',num2str(indexSecondary),'; convergence iteration: ',num2str(markerWhile)])
           elseif markerWhile > 3e3 && markerWhile <= 4e3
               threshPerc = 1.2;
%                disp(['index point: ',num2str(indexSecondary),'; convergence iteration: ',num2str(markerWhile)])
           elseif markerWhile > 4e3 && markerWhile <= 5e3
               threshPerc = 1.3;
%                disp(['index point: ',num2str(indexSecondary),'; convergence iteration: ',num2str(markerWhile)])
           elseif markerWhile > 5e3 && markerWhile <= 6e3
               threshPerc = 1.5;
%                disp(['index point: ',num2str(indexSecondary),'; convergence iteration: ',num2str(markerWhile)])
           elseif markerWhile > 6e3
               disp(['No convergence: iteration marker value of ', num2str(markerWhile)])
           else
               threshPerc = 1;
           end
           
       end
       
       massX(:,indexPrimary+1)= temporaryMassX(:,end-1);
       massY(:,indexPrimary+1)= temporaryMassY(:,end-1);
    

   
   if rem(i, floor(.008*length(t))) == 0
      %Vmsecondary(:,indexSecondary) = Vm(:,indexPrimary+1);
      Vmsecondary(:,indexSecondary) = Vm(:,indexPrimary+1);
      Tasecondary(:, indexSecondary) = Ta(:,indexPrimary+1);
      massXsecondary(:, indexSecondary) = massX(:,indexPrimary+1);
      massYsecondary(:, indexSecondary) = massY(:,indexPrimary+1);
      
      newTime(indexSecondary) = t(i);
      indexSecondary = indexSecondary + 1;
      disp(['start of time point: ',num2str(newTime(indexSecondary-1)),' ms out of ',num2str(endOfTime),' ms'])
   end
   
   Vm(:,indexPrimary) = b;
   Vmnondim(:,indexPrimary) = Vmnondim(:,indexPrimary+1);
   Ta(:, indexPrimary) = Ta(:, indexPrimary+1);
   massX(:, indexPrimary) = massX(:, indexPrimary+1);
   massY(:, indexPrimary) = massY(:, indexPrimary+1);
   
   
   
   
end


%%

%Part 3: Creation of Final Vectors


% Vm(:,end) = [];

Vmfinal = [];
rown = 1;
coln = 1;
indi = 1;
[rowinSec,columninSec] = size(Vmsecondary);

%putting Membrane Voltage into a 3D array
for i = 1:1:columninSec
    for j = 1:nodes*nodes
        Vmfinal(rown,coln,indi) = Vmsecondary(j,i);
        coln = coln + 1;
        if coln > nodes
            rown = rown + 1;
            coln = 1;
        end
    end
    rown = 1;
    coln = 1;
    indi = indi + 1;
end

Tafinal = zeros(((nodes+1)/2),((nodes+1)/2),101);
rown = 1;
coln = 1;
indi = 1;
[rowinSec,columninSec] = size(Tasecondary);

%putting Mass and Ta into a 3D array
for i = 1:1:columninSec
    for j = 1:(((nodes+1)/2)*((nodes+1)/2))
        Tafinal(rown,coln,indi) = Tasecondary(j,i);
        
        
        coln = coln + 1;
        
        if coln > (nodes+1)/2
            rown = rown + 1;
            coln = 1;
        end
    end
    rown = 1;
    coln = 1;
    indi = indi + 1;
end


clear Vm;
toc;

%%
%Part 4: Displaying Simulations


figure
pause(5);
for i = 1:1:length(Vmfinal(1,1,:))
    subplot(3,1,1)
    imagesc((Vmfinal(:,:,i)))
    xlim([1 nodes])
    ylim([1 nodes])
    title(['Voltage as f(x,y) at time = ', num2str(newTime(i)),' ms; i=',num2str(i)])
    xlabel('node in x')
    ylabel('node in y')
    colorbar('southoutside')
    caxis([minVal maxVal])
    
    subplot(3,1,2)
    plot((massXsecondary(:,i)), (massYsecondary(:,i)), 'rs')
    hold on
    plot(massXsecondary(1:(nodes+1)/2:end,i), massYsecondary(1:(nodes+1)/2:end,i), 'bs')
    plot(massXsecondary((nodes+1)/2:(nodes+1)/2:end,i), massYsecondary((nodes+1)/2:(nodes+1)/2:end,i), 'bs')
    plot(massXsecondary(1:1:(nodes+1)/2,i), massYsecondary(1:1:(nodes+1)/2,i), 'bs')
    plot(massXsecondary(((nodes+1)/2)^2-((nodes+1)/2)+1:1:((nodes+1)/2)^2,i), massYsecondary(((nodes+1)/2)^2-((nodes+1)/2)+1:1:((nodes+1)/2)^2,i), 'bs')
    ylim([1*bottomValuey 1*upperValuey]); xlim([1*leftValuex 1*rightValuex])
    title(['Mass lattice arrangement (x,y) at time = ', num2str(newTime(i))])
    xlabel('x value'); ylabel('y value')
    hold off

    subplot(3,1,3)
    imagesc(Tafinal(:,:,i))
    xlim([1 (nodes+1)/2])
    ylim([1 (nodes+1)/2])
    title(['Active Tension as f(x,y) at time = ', num2str(newTime(i)),' ms; i=',num2str(i)])
    xlabel('node in x')
    ylabel('node in y')
    colorbar('southoutside')
    caxis([0 2])
    
    pause(.05)
end


%For conduction velocity
% 
% NOI1 = 19;
% NOI2 = 149;
% peakNOI1 = max(Vmsecondary((nodes*nodes)-((NOI1+1)*nodes)+(nodes/2),:));
% peakNOI2 = max(Vmsecondary((nodes*nodes)-((NOI2+1)*nodes)+(nodes/2),:));
% iforpeakNOI1 = find(Vmsecondary((nodes*nodes)-((NOI1+1)*nodes)+(nodes/2),:) == peakNOI1);
% iforpeakNOI2 = find(Vmsecondary((nodes*nodes)-((NOI2+1)*nodes)+(nodes/2),:) == peakNOI2);
% tforpeakNOI1 = newTime(iforpeakNOI1);
% tforpeakNOI2 = newTime(iforpeakNOI2);
% 
% CondVel = (((dY + a).*(abs(NOI1-NOI2)))./(abs(tforpeakNOI1 - tforpeakNOI2)))*10;%m/s
% 
% figure
% hold on
% plot(newTime,Vmsecondary(nodes*nodes-NOI1*nodes+nodes/2,:),'r')
% plot(newTime,Vmsecondary(nodes*nodes-NOI2*nodes+nodes/2,:),'b')
% xlabel('time (ms)')
% ylabel('Vm (mV)')
% legend('1st Node of Interest','2nd Node of Interest')
% title('Propagation in Time for Two Nodes')