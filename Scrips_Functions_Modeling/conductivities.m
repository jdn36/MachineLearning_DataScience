%Calculating conductivities
display('Calculating conductivities to Use...')
% CV = K*sqrt( ( sigmaI * sigmaE )/( sigmaI + sigmaE ) )

% (based on sim) K = 39.55759 (or 29.0965, not entirely sure)
% %these original ones work great and give a CV of 33.962 cm/s
%Clerc 1976, Henriquez 1993
% sigmaFI_og = (1000)*(1/574);%mS/cm
% sigmaFE_og = (1000)*(1/160);%mS/cm
% sigmaSI_og = (1000)*(1/5171);%mS/cm
% sigmaSE_og = (1000)*(1/423);%mS/cm

%If we want 2X the CV, the speed up factor for the conductivities has
%to be 4X. It is a squared relationship.
%This new one with updated speed up factor works and gives a CV of around
%57 cm/s, Use these ones.

%get K experimentally through simulations
%AS OF NOW, K IS 48.9358 FOR CV CALCULATIONS USING BI SIGMAS AND CV WAS 63 CM/S WITH sigmaFI_Bi = 2.3519 AND
%sigmaFE_Bi = 5.6250 AND FI = .6 AND FE = .4 AND CHI = 1500 SO SPEED UP FACTOR IS .961

K_exp = 46.5905;%39.55759;48.9358;29.0965

speedUpFactor = ( 3.8 )^2;
sigmaFI_og = speedUpFactor*(1000)*(1/574);%mS/cm
sigmaFE_og = speedUpFactor*(1000)*(1/160);%mS/cm
sigmaSI_og = speedUpFactor*(1000)*(1/5171);%mS/cm
sigmaSE_og = speedUpFactor*(1000)*(1/423);%mS/cm

%fraction of intracellular and extracellular space per cell element
fi = .6;
fe = .4;

sigmaFI_Bi = sigmaFI_og*fi
sigmaFE_Bi = sigmaFE_og*fe
sigmaSI_Bi = sigmaSI_og*fi
sigmaSE_Bi = sigmaSE_og*fe

CV_expected = K_exp*sqrt( ( sigmaFI_Bi * sigmaFE_Bi )/( sigmaFI_Bi + sigmaFE_Bi ) )

% sigmaF_Mono = (sigmaFI_og * sigmaFE_og)/(sigmaFI_og + sigmaFE_og)
% sigmaS_Mono = (sigmaSI_og * sigmaSE_og)/(sigmaSI_og + sigmaSE_og)

sigmaF_Mono = (sigmaFI_Bi * sigmaFE_Bi)/(sigmaFI_Bi + sigmaFE_Bi)
sigmaS_Mono = (sigmaSI_Bi * sigmaSE_Bi)/(sigmaSI_Bi + sigmaSE_Bi)


%chi = area/volume = 2pi*a*l / (pi*a^2*l/fi) = 2*fi/a or 2/a if fi assumed
%to be 0
%Courtemanche cell is 16 um in diameters
chi_og = 2500; %cm^-1
chiBi = chi_og*fi






