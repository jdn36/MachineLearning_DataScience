

%matrixTissuePointsExport = makeFullyTransmural_2DGrid(meshX, meshY, meshZ, fibX, fibY, fibZ, nameOfFullTransmuralCompact, dx)

nameOfCoordsFile = '../SimResults/IdealMediumGeometryAllPoints.csv';

coordPoints = csvread(nameOfCoordsFile);
meshX = coordPoints(:,3); meshY = coordPoints(:,4); meshZ = coordPoints(:,5);

nameOfFibrosisCoordsFile = '../SimResults/BMES_CompactFibrosis_FinerCoordsIdealMediumGeometry_35perc_V1.csv';
coordFibrosisPoints = csvread(nameOfFibrosisCoordsFile);
fibX = coordFibrosisPoints(:,1); fibY = coordFibrosisPoints(:,2); fibZ = coordFibrosisPoints(:,3);

figure
plot3(meshX, meshY, meshZ,'ro')
hold on
plot3(fibX, fibY, fibZ,'ko')

%% Compact
%V1
nameOfFibrosisCoordsFile = '../SimResults/BMES_CompactFibrosis_FinerCoordsIdealMediumGeometry_20perc_V1.csv';
coordFibrosisPoints = csvread(nameOfFibrosisCoordsFile);
fibX = coordFibrosisPoints(:,1); fibY = coordFibrosisPoints(:,2); fibZ = coordFibrosisPoints(:,3);
nameOfFullTransmuralCompact = '../SimResults/BMES_CompactFibrosisFT_FinerCoordsIdealMediumGeometry_20perc_V1.csv';
matrixTissuePointsExport = makeFullyTransmural_2DGrid(meshX, meshY, meshZ, fibX, fibY, fibZ, nameOfFullTransmuralCompact, .025);
%V2
nameOfFibrosisCoordsFile = '../SimResults/BMES_CompactFibrosis_FinerCoordsIdealMediumGeometry_20perc_V2.csv';
coordFibrosisPoints = csvread(nameOfFibrosisCoordsFile);
fibX = coordFibrosisPoints(:,1); fibY = coordFibrosisPoints(:,2); fibZ = coordFibrosisPoints(:,3);
nameOfFullTransmuralCompact = '../SimResults/BMES_CompactFibrosisFT_FinerCoordsIdealMediumGeometry_20perc_V2.csv';
matrixTissuePointsExport = makeFullyTransmural_2DGrid(meshX, meshY, meshZ, fibX, fibY, fibZ, nameOfFullTransmuralCompact, .025);
%V3
nameOfFibrosisCoordsFile = '../SimResults/BMES_CompactFibrosis_FinerCoordsIdealMediumGeometry_20perc_V3.csv';
coordFibrosisPoints = csvread(nameOfFibrosisCoordsFile);
fibX = coordFibrosisPoints(:,1); fibY = coordFibrosisPoints(:,2); fibZ = coordFibrosisPoints(:,3);
nameOfFullTransmuralCompact = '../SimResults/BMES_CompactFibrosisFT_FinerCoordsIdealMediumGeometry_20perc_V3.csv';
matrixTissuePointsExport = makeFullyTransmural_2DGrid(meshX, meshY, meshZ, fibX, fibY, fibZ, nameOfFullTransmuralCompact, .025);
%V4
nameOfFibrosisCoordsFile = '../SimResults/BMES_CompactFibrosis_FinerCoordsIdealMediumGeometry_20perc_V4.csv';
coordFibrosisPoints = csvread(nameOfFibrosisCoordsFile);
fibX = coordFibrosisPoints(:,1); fibY = coordFibrosisPoints(:,2); fibZ = coordFibrosisPoints(:,3);
nameOfFullTransmuralCompact = '../SimResults/BMES_CompactFibrosisFT_FinerCoordsIdealMediumGeometry_20perc_V4.csv';
matrixTissuePointsExport = makeFullyTransmural_2DGrid(meshX, meshY, meshZ, fibX, fibY, fibZ, nameOfFullTransmuralCompact, .025);
%V5
nameOfFibrosisCoordsFile = '../SimResults/BMES_CompactFibrosis_FinerCoordsIdealMediumGeometry_20perc_V5.csv';
coordFibrosisPoints = csvread(nameOfFibrosisCoordsFile);
fibX = coordFibrosisPoints(:,1); fibY = coordFibrosisPoints(:,2); fibZ = coordFibrosisPoints(:,3);
nameOfFullTransmuralCompact = '../SimResults/BMES_CompactFibrosisFT_FinerCoordsIdealMediumGeometry_20perc_V5.csv';
matrixTissuePointsExport = makeFullyTransmural_2DGrid(meshX, meshY, meshZ, fibX, fibY, fibZ, nameOfFullTransmuralCompact, .025);

%% Diffuse
%V1
nameOfFibrosisCoordsFile = '../SimResults/BMES_DiffuseFibrosis_FinerCoordsIdealMediumGeometry_20perc_V1.csv';
coordFibrosisPoints = csvread(nameOfFibrosisCoordsFile);
fibX = coordFibrosisPoints(:,1); fibY = coordFibrosisPoints(:,2); fibZ = coordFibrosisPoints(:,3);
nameOfFullTransmuralDiffuse = '../SimResults/BMES_DiffuseFibrosisFT_FinerCoordsIdealMediumGeometry_20perc_V1.csv';
matrixTissuePointsExport = makeFullyTransmural_2DGrid(meshX, meshY, meshZ, fibX, fibY, fibZ, nameOfFullTransmuralDiffuse, .025);
%V2
nameOfFibrosisCoordsFile = '../SimResults/BMES_DiffuseFibrosis_FinerCoordsIdealMediumGeometry_20perc_V2.csv';
coordFibrosisPoints = csvread(nameOfFibrosisCoordsFile);
fibX = coordFibrosisPoints(:,1); fibY = coordFibrosisPoints(:,2); fibZ = coordFibrosisPoints(:,3);
nameOfFullTransmuralDiffuse = '../SimResults/BMES_DiffuseFibrosisFT_FinerCoordsIdealMediumGeometry_20perc_V2.csv';
matrixTissuePointsExport = makeFullyTransmural_2DGrid(meshX, meshY, meshZ, fibX, fibY, fibZ, nameOfFullTransmuralDiffuse, .025);
%V3
nameOfFibrosisCoordsFile = '../SimResults/BMES_DiffuseFibrosis_FinerCoordsIdealMediumGeometry_20perc_V3.csv';
coordFibrosisPoints = csvread(nameOfFibrosisCoordsFile);
fibX = coordFibrosisPoints(:,1); fibY = coordFibrosisPoints(:,2); fibZ = coordFibrosisPoints(:,3);
nameOfFullTransmuralDiffuse = '../SimResults/BMES_DiffuseFibrosisFT_FinerCoordsIdealMediumGeometry_20perc_V3.csv';
matrixTissuePointsExport = makeFullyTransmural_2DGrid(meshX, meshY, meshZ, fibX, fibY, fibZ, nameOfFullTransmuralDiffuse, .025);
%V4
nameOfFibrosisCoordsFile = '../SimResults/BMES_DiffuseFibrosis_FinerCoordsIdealMediumGeometry_20perc_V4.csv';
coordFibrosisPoints = csvread(nameOfFibrosisCoordsFile);
fibX = coordFibrosisPoints(:,1); fibY = coordFibrosisPoints(:,2); fibZ = coordFibrosisPoints(:,3);
nameOfFullTransmuralDiffuse = '../SimResults/BMES_DiffuseFibrosisFT_FinerCoordsIdealMediumGeometry_20perc_V4.csv';
matrixTissuePointsExport = makeFullyTransmural_2DGrid(meshX, meshY, meshZ, fibX, fibY, fibZ, nameOfFullTransmuralDiffuse, .025);
%V5
nameOfFibrosisCoordsFile = '../SimResults/BMES_DiffuseFibrosis_FinerCoordsIdealMediumGeometry_20perc_V5.csv';
coordFibrosisPoints = csvread(nameOfFibrosisCoordsFile);
fibX = coordFibrosisPoints(:,1); fibY = coordFibrosisPoints(:,2); fibZ = coordFibrosisPoints(:,3);
nameOfFullTransmuralDiffuse = '../SimResults/BMES_DiffuseFibrosisFT_FinerCoordsIdealMediumGeometry_20perc_V5.csv';
matrixTissuePointsExport = makeFullyTransmural_2DGrid(meshX, meshY, meshZ, fibX, fibY, fibZ, nameOfFullTransmuralDiffuse, .025);

%% Interstitial
%V1
nameOfFibrosisCoordsFile = '../SimResults/BMES_InterstitialFibrosis_FinerCoordsIdealMediumGeometry_20perc_V1.csv';
coordFibrosisPoints = csvread(nameOfFibrosisCoordsFile);
fibX = coordFibrosisPoints(:,1); fibY = coordFibrosisPoints(:,2); fibZ = coordFibrosisPoints(:,3);
nameOfFullTransmuralInterstitial = '../SimResults/BMES_InterstitialFibrosisFT_FinerCoordsIdealMediumGeometry_20perc_V1.csv';
matrixTissuePointsExport = makeFullyTransmural_2DGrid(meshX, meshY, meshZ, fibX, fibY, fibZ, nameOfFullTransmuralInterstitial, .025);
%V2
nameOfFibrosisCoordsFile = '../SimResults/BMES_InterstitialFibrosis_FinerCoordsIdealMediumGeometry_20perc_V2.csv';
coordFibrosisPoints = csvread(nameOfFibrosisCoordsFile);
fibX = coordFibrosisPoints(:,1); fibY = coordFibrosisPoints(:,2); fibZ = coordFibrosisPoints(:,3);
nameOfFullTransmuralInterstitial = '../SimResults/BMES_InterstitialFibrosisFT_FinerCoordsIdealMediumGeometry_20perc_V2.csv';
matrixTissuePointsExport = makeFullyTransmural_2DGrid(meshX, meshY, meshZ, fibX, fibY, fibZ, nameOfFullTransmuralInterstitial, .025);
%V3
nameOfFibrosisCoordsFile = '../SimResults/BMES_InterstitialFibrosis_FinerCoordsIdealMediumGeometry_20perc_V3.csv';
coordFibrosisPoints = csvread(nameOfFibrosisCoordsFile);
fibX = coordFibrosisPoints(:,1); fibY = coordFibrosisPoints(:,2); fibZ = coordFibrosisPoints(:,3);
nameOfFullTransmuralInterstitial = '../SimResults/BMES_InterstitialFibrosisFT_FinerCoordsIdealMediumGeometry_20perc_V3.csv';
matrixTissuePointsExport = makeFullyTransmural_2DGrid(meshX, meshY, meshZ, fibX, fibY, fibZ, nameOfFullTransmuralInterstitial, .025);
%V4
nameOfFibrosisCoordsFile = '../SimResults/BMES_InterstitialFibrosis_FinerCoordsIdealMediumGeometry_20perc_V4.csv';
coordFibrosisPoints = csvread(nameOfFibrosisCoordsFile);
fibX = coordFibrosisPoints(:,1); fibY = coordFibrosisPoints(:,2); fibZ = coordFibrosisPoints(:,3);
nameOfFullTransmuralInterstitial = '../SimResults/BMES_InterstitialFibrosisFT_FinerCoordsIdealMediumGeometry_20perc_V4.csv';
matrixTissuePointsExport = makeFullyTransmural_2DGrid(meshX, meshY, meshZ, fibX, fibY, fibZ, nameOfFullTransmuralInterstitial, .025);
%V5
nameOfFibrosisCoordsFile = '../SimResults/BMES_InterstitialFibrosis_FinerCoordsIdealMediumGeometry_20perc_V5.csv';
coordFibrosisPoints = csvread(nameOfFibrosisCoordsFile);
fibX = coordFibrosisPoints(:,1); fibY = coordFibrosisPoints(:,2); fibZ = coordFibrosisPoints(:,3);
nameOfFullTransmuralInterstitial = '../SimResults/BMES_InterstitialFibrosisFT_FinerCoordsIdealMediumGeometry_20perc_V5.csv';
matrixTissuePointsExport = makeFullyTransmural_2DGrid(meshX, meshY, meshZ, fibX, fibY, fibZ, nameOfFullTransmuralInterstitial, .025);


%% Patchy

%V1
nameOfFibrosisCoordsFile = '../SimResults/BMES_PatchyFibrosis_FinerCoordsIdealMediumGeometry_20perc_V1.csv';
coordFibrosisPoints = csvread(nameOfFibrosisCoordsFile);
fibX = coordFibrosisPoints(:,1); fibY = coordFibrosisPoints(:,2); fibZ = coordFibrosisPoints(:,3);
nameOfFullTransmuralPatchy = '../SimResults/BMES_PatchyFibrosisFT_FinerCoordsIdealMediumGeometry_20perc_V1.csv';
matrixTissuePointsExport = makeFullyTransmural_2DGrid(meshX, meshY, meshZ, fibX, fibY, fibZ, nameOfFullTransmuralPatchy, .025);
%V2
nameOfFibrosisCoordsFile = '../SimResults/BMES_PatchyFibrosis_FinerCoordsIdealMediumGeometry_20perc_V2.csv';
coordFibrosisPoints = csvread(nameOfFibrosisCoordsFile);
fibX = coordFibrosisPoints(:,1); fibY = coordFibrosisPoints(:,2); fibZ = coordFibrosisPoints(:,3);
nameOfFullTransmuralPatchy = '../SimResults/BMES_PatchyFibrosisFT_FinerCoordsIdealMediumGeometry_20perc_V2.csv';
matrixTissuePointsExport = makeFullyTransmural_2DGrid(meshX, meshY, meshZ, fibX, fibY, fibZ, nameOfFullTransmuralPatchy, .025);
%V3
nameOfFibrosisCoordsFile = '../SimResults/BMES_PatchyFibrosis_FinerCoordsIdealMediumGeometry_20perc_V3.csv';
coordFibrosisPoints = csvread(nameOfFibrosisCoordsFile);
fibX = coordFibrosisPoints(:,1); fibY = coordFibrosisPoints(:,2); fibZ = coordFibrosisPoints(:,3);
nameOfFullTransmuralPatchy = '../SimResults/BMES_PatchyFibrosisFT_FinerCoordsIdealMediumGeometry_20perc_V3.csv';
matrixTissuePointsExport = makeFullyTransmural_2DGrid(meshX, meshY, meshZ, fibX, fibY, fibZ, nameOfFullTransmuralPatchy, .025);
%V4
nameOfFibrosisCoordsFile = '../SimResults/BMES_PatchyFibrosis_FinerCoordsIdealMediumGeometry_20perc_V4.csv';
coordFibrosisPoints = csvread(nameOfFibrosisCoordsFile);
fibX = coordFibrosisPoints(:,1); fibY = coordFibrosisPoints(:,2); fibZ = coordFibrosisPoints(:,3);
nameOfFullTransmuralPatchy = '../SimResults/BMES_PatchyFibrosisFT_FinerCoordsIdealMediumGeometry_20perc_V4.csv';
matrixTissuePointsExport = makeFullyTransmural_2DGrid(meshX, meshY, meshZ, fibX, fibY, fibZ, nameOfFullTransmuralPatchy, .025);
%V5
nameOfFibrosisCoordsFile = '../SimResults/BMES_PatchyFibrosis_FinerCoordsIdealMediumGeometry_20perc_V5.csv';
coordFibrosisPoints = csvread(nameOfFibrosisCoordsFile);
fibX = coordFibrosisPoints(:,1); fibY = coordFibrosisPoints(:,2); fibZ = coordFibrosisPoints(:,3);
nameOfFullTransmuralPatchy = '../SimResults/BMES_PatchyFibrosisFT_FinerCoordsIdealMediumGeometry_20perc_V5.csv';
matrixTissuePointsExport = makeFullyTransmural_2DGrid(meshX, meshY, meshZ, fibX, fibY, fibZ, nameOfFullTransmuralPatchy, .025);





