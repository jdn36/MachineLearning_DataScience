



%Grouping in 12 different "Categories":
% - Compact 10%,
% - Compact 35%,
% - Compact 60%,
% - Diffuse 10%,
% - Diffuse 35%,
% - Diffuse 60%,
% - Interstitial 10%,
% - Interstitial 35%,
% - Interstitial 60%,
% - Patchy 10%,
% - Patchy 35%,
% - Patchy 60%,

%Out of these 12 "categories":
%- Characterize signals (means, stds, etc) using deflection count, EGM
%duration, P-P Amplitude, etc.
%- Compare the variability of characteristics accross different patient's geometries
%- Compare the variability of characteristics and behavior within each type
%but across different degrees of fibrosis
%- Compare the variability of characteristics and behavior within each
%degree of fibrosis but different types
%See how simulated signals compare to real patient signals?
%See how LAT's of simulated signals compare to real patient signals LAT's?
%See how lumped characteristics of simulated signals compare to real patient signals lumped characteristics?

tic

% load('../SimResults/csPacingDatabase.mat')
%SIMULATED DATA
normalizeBool = false;
phie_dt = .1;
numDataPoints = 100;
levelDec = 4;
vecLevel = 1:levelDec;

%Variable: output_folder, exportVe_Name, fibrosis_coords_matrix, 
typesFibrosis = {'Compact', 'Diffuse', 'Interstitial', 'Patchy'};

each_PPASimulated = [];
each_downstrokeSimulated = [];
each_upstrokeRiseSimulated = [];
each_upstrokeBackSimulated = [];
each_PositiveDDSimulated = [];
each_NegativeDDSimulated = [];
each_baselineSimulated = [];
each_SignalWSimulated = [];
each_EGMDurationSimulated = [];
each_deflectionCountSimulated = [];
each_AlldistanceToProjSimulated = [];
each_baselineMeanSDSWSimulated = [];
patient_lumped = [];
cdAll_eachSimulated = [];
eachSimulated_dataToClassify = [];

%preparing arrays
global_PPASimulated = [];
global_downstrokeSimulated = [];
global_upstrokeRiseSimulated = [];
global_upstrokeBackSimulated = [];
global_PositiveDDSimulated = [];
global_NegativeDDSimulated = [];
global_baselineSimulated = [];
global_SignalWSimulated = [];
global_EGMDurationSimulated = [];
global_deflectionCountSimulated = [];
global_AlldistanceToProjSimulated = [];
global_baselineMeanSDSWSimulated = [];
global_simDataCut = [];
global_simDataCutOG = [];


global_PCA1vect_VAR_Persistent = [];
global_PCA2vect_VAR_Persistent = [];
global_PCA3vect_VAR_Persistent = [];
global_PCA4vect_VAR_Persistent = [];
global_PCA5vect_VAR_Persistent = [];
global_idxPCAVAR = [];
global_CPCAVAR = [];
global_DPCAVAR = [];
global_dists_Patient_VAR = [];
global_dists_Simulated_common_VAR = [];
global_dists_Simulated_different_VAR = [];

global_PCA1vect_WD_Persistent = [];
global_PCA2vect_WD_Persistent = [];
global_PCA3vect_WD_Persistent = [];
global_PCA4vect_WD_Persistent = [];
global_PCA5vect_WD_Persistent = [];
global_idxWD = [];
global_CWD = [];
global_DWD = [];
global_dists_Patient_WD = [];
global_dists_Simulated_common_WD = [];
global_dists_Simulated_different_WD = [];

each_LATs_Persistent = [];
global_LATs_Persistent = [];


global_loc_similarity_score = [];
global_loc_xcorr = [];
global_loc_LAT_error = [];

global_PCA1vect_sim_VAR_Persistent = [];
global_PCA2vect_sim_VAR_Persistent = [];
global_PCA3vect_sim_VAR_Persistent = [];
global_PCA4vect_sim_VAR_Persistent = [];
global_PCA5vect_sim_VAR_Persistent = [];

global_PCA1vect_sim_WD_Persistent = [];
global_PCA2vect_sim_WD_Persistent = [];
global_PCA3vect_sim_WD_Persistent = [];
global_PCA4vect_sim_WD_Persistent = [];
global_PCA5vect_sim_WD_Persistent = [];


fType = [];
patientNumber = [];
fPercentage = [];
patternNumber = [];

fType_lumped = [];
patientNumber_lumped = [];
fPercentage_lumped = [];
patternNumber_lumped = [];
similarity_score_lumped = [];
LAT_error_lumped = [];
xcorr_lumped = [];

fTypeNAN = [];
patientNumberNAN = [];
fPercentageNAN = [];
patternNumberNAN = [];


geometryPoints_Paroxysmal2 = csvread('../SimResults/P1simsStuff/Paroxysmal2_AllPoints_Refined.csv');
meshX_Paroxysmal2 = geometryPoints_Paroxysmal2(:,3); meshY_Paroxysmal2 = geometryPoints_Paroxysmal2(:,4); meshZ_Paroxysmal2 = geometryPoints_Paroxysmal2(:,5);
electrodePoints_Paroxysmal2 = csvread('../SimResults/P1simsStuff/P1_Paroxysmal2_RealGeometryElecPoints_OG.csv');
elecX_Paroxysmal2 = electrodePoints_Paroxysmal2(:,1); elecY_Paroxysmal2 = electrodePoints_Paroxysmal2(:,2); elecZ_Paroxysmal2 = electrodePoints_Paroxysmal2(:,3);

geometryPoints_Paroxysmal3 = csvread('../SimResults/P1simsStuff/Paroxysmal3_AllPoints_Refined.csv');
meshX_Paroxysmal3 = geometryPoints_Paroxysmal3(:,3); meshY_Paroxysmal3 = geometryPoints_Paroxysmal3(:,4); meshZ_Paroxysmal3 = geometryPoints_Paroxysmal3(:,5);
electrodePoints_Paroxysmal3 = csvread('../SimResults/P1simsStuff/P1_Paroxysmal3_RealGeometryElecPoints_OG.csv');
elecX_Paroxysmal3 = electrodePoints_Paroxysmal3(:,1); elecY_Paroxysmal3 = electrodePoints_Paroxysmal3(:,2); elecZ_Paroxysmal3 = electrodePoints_Paroxysmal3(:,3);




geometryPoints_Persistent6 = csvread('../SimResults/P1simsStuff/Persistent6_AllPoints_Refined.csv');
meshX_Persistent6 = geometryPoints_Persistent6(:,3); meshY_Persistent6 = geometryPoints_Persistent6(:,4); meshZ_Persistent6 = geometryPoints_Persistent6(:,5);
electrodePoints_Persistent6 = csvread('../SimResults/P1simsStuff/P1_Persistent6_RealGeometryElecPoints_OG.csv');
elecX_Persistent6 = electrodePoints_Persistent6(:,1); elecY_Persistent6 = electrodePoints_Persistent6(:,2); elecZ_Persistent6 = electrodePoints_Persistent6(:,3);


geometryPoints_Persisitent7 = csvread('../SimResults/P1simsStuff/Persistent7_AllPoints_Refined.csv');
meshX_Persisitent7 = geometryPoints_Persisitent7(:,3); meshY_Persisitent7 = geometryPoints_Persisitent7(:,4); meshZ_Persisitent7 = geometryPoints_Persisitent7(:,5);
electrodePoints_Persistent7 = csvread('../SimResults/P1simsStuff/P1_Persistent7_RealGeometryElecPoints_OG.csv');
elecX_Persistent7 = electrodePoints_Persistent7(:,1); elecY_Persistent7 = electrodePoints_Persistent7(:,2); elecZ_Persistent7 = electrodePoints_Persistent7(:,3);


geometryPoints_Persisitent8 = csvread('../SimResults/P1simsStuff/Persistent8_AllPoints_Refined.csv');
meshX_Persisitent8 = geometryPoints_Persisitent8(:,3); meshY_Persisitent8 = geometryPoints_Persisitent8(:,4); meshZ_Persisitent8 = geometryPoints_Persisitent8(:,5);
electrodePoints_Persistent8 = csvread('../SimResults/P1simsStuff/P1_Persistent8_RealGeometryElecPoints_OG.csv');
elecX_Persistent8 = electrodePoints_Persistent8(:,1); elecY_Persistent8 = electrodePoints_Persistent8(:,2); elecZ_Persistent8 = electrodePoints_Persistent8(:,3);


geometryPoints_Persisitent9 = csvread('../SimResults/P1simsStuff/Persistent9_AllPoints_Refined.csv');
meshX_Persisitent9 = geometryPoints_Persisitent9(:,3); meshY_Persisitent9 = geometryPoints_Persisitent9(:,4); meshZ_Persisitent9 = geometryPoints_Persisitent9(:,5);
electrodePoints_Persistent9 = csvread('../SimResults/P1simsStuff/P1_Persistent9_RealGeometryElecPoints_OG.csv');
elecX_Persistent9 = electrodePoints_Persistent9(:,1); elecY_Persistent9 = electrodePoints_Persistent9(:,2); elecZ_Persistent9 = electrodePoints_Persistent9(:,3);



for pN1 = 6:9
    patientNum = pN1;
    loc_geometryPoints = csvread(['../SimResults/P1simsStuff/Persistent',num2str(patientNum),'_AllPoints_Refined.csv']);
    loc_meshX = loc_geometryPoints(:,3); loc_meshY = loc_geometryPoints(:,4); loc_meshZ = loc_geometryPoints(:,5);
    loc_electrodePoints = csvread(['../SimResults/P1simsStuff/P1_Persistent',num2str(patientNum),'_RealGeometryElecPoints_OG.csv']);
    loc_elecX = loc_electrodePoints(:,1); loc_elecY = loc_electrodePoints(:,2); loc_elecZ = loc_electrodePoints(:,3);
    ptCloudSurface = pointCloud([loc_meshX, loc_meshY, loc_meshZ]);


    actTimesCurrent = expDB.dataSet{patientNum}.map.continuousMapping.revisedLat(:,1);

    elec1=[];
    locSignalNorm_Matrix = [];
    %fullEgms -> point in space, which electrode, time points
    for i=1:length(loc_elecX)
        elec_1=[squeeze(expDB.dataSet{patientNum}.map.continuousMapping.fullEgms(i,1,1:140))];
        elec1 = [elec1,elec_1];

        [EGMduration, SD_curve, aboveValues, DeflectionCount, baselineMeanSDSW, bumpinessFactor, xtoCheck, ytoCheck] = SD_SlidingWindowF([1:1:length( elec1(:,i) )], elec1(:,i), 1, 40, 15, 70, 3, numDataPoints);
        locSignalNorm = ( elec1(:,i) - min( elec1(:,i) ) ) ./ ( max( elec1(:,i) ) - min( elec1(:,i) ) );
        [dVeMono_dt, downstrokeMono, midSignalIndMono, signalDiffMonoPointsTime, upstrokeRiseMono, upstrokeBackMono, riseSignalIndMono, backSignalIndMono, PositiveDeflectionDuration_Mono, NegativeDeflectionDuration_Mono, baselineMono, signalWidth, timePeak, timeTrough] = EGM_Characterizing_Function(locSignalNorm, 1:140, 1, 1, EGMduration );
    
        localPPA = max(ytoCheck) - min(ytoCheck);
        
        each_PPASimulated = [each_PPASimulated; localPPA];
        each_downstrokeSimulated = [each_downstrokeSimulated; downstrokeMono];
        each_upstrokeRiseSimulated = [each_upstrokeRiseSimulated; upstrokeRiseMono];
        each_upstrokeBackSimulated = [each_upstrokeBackSimulated; upstrokeBackMono];
        each_PositiveDDSimulated = [each_PositiveDDSimulated; PositiveDeflectionDuration_Mono];
        each_NegativeDDSimulated = [each_NegativeDDSimulated; NegativeDeflectionDuration_Mono];
        each_baselineSimulated = [each_baselineSimulated; baselineMono];
        each_SignalWSimulated = [each_SignalWSimulated; signalWidth];
        each_EGMDurationSimulated = [each_EGMDurationSimulated; EGMduration];
    
        loc_normytoCheck = ( ytoCheck - min(ytoCheck) )./( max(ytoCheck) - min(ytoCheck) );
        locSignalNorm_Matrix = [locSignalNorm_Matrix; reshape(loc_normytoCheck, [1, length(loc_normytoCheck)])];

        [c,l] = wavedec(loc_normytoCheck,levelDec,"rbio4.4");
        approx = appcoef(c,l,"rbio4.4");
        [cd1,cd2,cd3,cd4] = detcoef(c,l,vecLevel);
        cdAll_eachSimulated = [ cdAll_eachSimulated, cd4 ];

        % loc_simDataCut = [loc_simDataCut, loc_normytoCheck];
        % loc_simDataCutOG = [loc_simDataCutOG, ytoCheck];
    
        each_deflectionCountSimulated = [each_deflectionCountSimulated; DeflectionCount];

        [indLoc,distsLoc] = findNearestNeighbors(ptCloudSurface,[loc_elecX(i), loc_elecY(i), loc_elecZ(i)],1);
        each_AlldistanceToProjSimulated = [each_AlldistanceToProjSimulated;  distsLoc ];
    
        each_baselineMeanSDSWSimulated = [each_baselineMeanSDSWSimulated; baselineMeanSDSW];
        patient_lumped = [patient_lumped; patientNum];

    end

    eachSimulated_dataToClassify = [ eachSimulated_dataToClassify; each_downstrokeSimulated, each_upstrokeRiseSimulated, each_upstrokeBackSimulated, each_SignalWSimulated, each_baselineSimulated, each_EGMDurationSimulated, each_AlldistanceToProjSimulated, each_PPASimulated, each_deflectionCountSimulated];
    loc_dataToClassify = eachSimulated_dataToClassify( find(patient_lumped == patientNum) , : );
    loc_dataToClassify_Scaled = loc_dataToClassify;

    maxScalingData = [];
    minScalingData = [];
    % % % % SCALING TO MAKE ALL SCALES MATCH
    for i = 1:length( loc_dataToClassify(1,:) )
        if max(loc_dataToClassify(:,i)) - min(loc_dataToClassify(:,i)) > 0
            loc_dataToClassify_Scaled(:,i) = ( loc_dataToClassify(:,i) - min(loc_dataToClassify(:,i)) )./( max(loc_dataToClassify(:,i)) - min(loc_dataToClassify(:,i)) );
            
            maxScalingData(i,2) = max(loc_dataToClassify(:,i));
            minScalingData(i,2) = min(loc_dataToClassify(:,i));
        end
    end
    
    
    maxScalingData = max(maxScalingData');
    minScalingData = min(minScalingData');
    
    [coeffVAR,scoreVAR,latentVAR,tsquaredVAR,explainedVAR,muVAR] = pca(loc_dataToClassify_Scaled);
    
    PCA1vectVAR = [];
    PCA2vectVAR = [];
    PCA3vectVAR = [];
    PCA4vectVAR = [];
    PCA5vectVAR = [];
    
    for i = 1: length(loc_dataToClassify_Scaled(:,1))
        PCA1vectVAR(i) = coeffVAR(:,1)'*loc_dataToClassify_Scaled(i,:)';
        PCA2vectVAR(i) = coeffVAR(:,2)'*loc_dataToClassify_Scaled(i,:)';
        PCA3vectVAR(i) = coeffVAR(:,3)'*loc_dataToClassify_Scaled(i,:)';
        PCA4vectVAR(i) = coeffVAR(:,4)'*loc_dataToClassify_Scaled(i,:)';
        PCA5vectVAR(i) = coeffVAR(:,5)'*loc_dataToClassify_Scaled(i,:)';
    end
    
    opts = statset('Display','final');
    [idxPCAVAR,CPCAVAR,sumdPCAVAR,DPCAVAR] = kmeans([PCA1vectVAR', PCA2vectVAR', PCA3vectVAR', PCA4vectVAR', PCA5vectVAR'],2,'Distance','cosine','EmptyAction','drop','MaxIter',1000,'Replicates',10,'Options',opts);
    
    indexCluster1_Var = find(idxPCAVAR == 1);
    indexCluster2_Var = find(idxPCAVAR == 2);

    global_PCA1vect_VAR_Persistent = [global_PCA1vect_VAR_Persistent; PCA1vectVAR'];
    global_PCA2vect_VAR_Persistent = [global_PCA2vect_VAR_Persistent; PCA2vectVAR'];
    global_PCA3vect_VAR_Persistent = [global_PCA3vect_VAR_Persistent; PCA3vectVAR'];
    global_PCA4vect_VAR_Persistent = [global_PCA4vect_VAR_Persistent; PCA4vectVAR'];
    global_PCA5vect_VAR_Persistent = [global_PCA5vect_VAR_Persistent; PCA5vectVAR'];
    global_idxPCAVAR = [global_idxPCAVAR; idxPCAVAR];
    global_CPCAVAR = [global_CPCAVAR; CPCAVAR];
    global_DPCAVAR = [global_DPCAVAR; DPCAVAR];






    [coeffWD,scoreWD,latentWD,tsquaredWD,explainedWD,muWD] = pca(cdAll_eachSimulated(:,find(patient_lumped == patientNum))');

    PCA1vectWD = [];
    PCA2vectWD = [];
    PCA3vectWD = [];
    PCA4vectWD = [];
    PCA5vectWD = [];
    
    for i = 1:length(cdAll_eachSimulated(1,find(patient_lumped == patientNum)))
        PCA1vectWD(i) = coeffWD(:,1)'*cdAll_eachSimulated(:,i);
        PCA2vectWD(i) = coeffWD(:,2)'*cdAll_eachSimulated(:,i);
        PCA3vectWD(i) = coeffWD(:,3)'*cdAll_eachSimulated(:,i);
        PCA4vectWD(i) = coeffWD(:,4)'*cdAll_eachSimulated(:,i);
        PCA5vectWD(i) = coeffWD(:,5)'*cdAll_eachSimulated(:,i);
    end
    
    opts = statset('Display','final');
    [idxWD,CWD,sumdWD,DWD] = kmeans([PCA1vectWD', PCA2vectWD', PCA3vectWD', PCA4vectWD', PCA5vectWD'],3,'Distance','cosine','EmptyAction','drop','MaxIter',1000,'Replicates',10,'Options',opts);
    
    indexCluster1_WD = find(idxWD == 1);
    indexCluster2_WD = find(idxWD == 2);
    indexCluster3_WD = find(idxWD == 3);

    global_PCA1vect_WD_Persistent = [global_PCA1vect_WD_Persistent; PCA1vectWD'];
    global_PCA2vect_WD_Persistent = [global_PCA2vect_WD_Persistent; PCA2vectWD'];
    global_PCA3vect_WD_Persistent = [global_PCA3vect_WD_Persistent; PCA3vectWD'];
    global_PCA4vect_WD_Persistent = [global_PCA4vect_WD_Persistent; PCA4vectWD'];
    global_PCA5vect_WD_Persistent = [global_PCA5vect_WD_Persistent; PCA5vectWD'];
    global_idxWD = [global_idxWD; idxWD];
    global_CWD = [global_CWD; CWD];
    global_DWD = [global_DWD; DWD];



    for tF = 1:length(typesFibrosis)
        for fP = [10,35,60]
            fibrosisPercentage = fP;
            fibrosisType = typesFibrosis{tF};


            for pN2 = 1:10
                patternNum = pN2;
                locSignals_unfiltered = csvread([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNum), '_', fibrosisType,'FibrosisNT', num2str(fibrosisPercentage), '_pattern', num2str(patternNum), '_2000points.csv' ]);
                locSignals_unfiltered = locSignals_unfiltered(1:1400,:);

                if ~isnan(mean(mean(locSignals_unfiltered')))
                    locSignals = UpsampleAndFilterFunction(locSignals_unfiltered', phie_dt, 5, 200, 1000, 10,normalizeBool)';
                    loc_actTimesSim = [];
                    loc_similarity_score = [];
                    loc_LAT_error = [];
                    loc_xcorr = [];
    
                    for i = 1:length( locSignals(1,:) )
                        [EGMduration, SD_curve, aboveValues, DeflectionCount, baselineMeanSDSW, bumpinessFactor, xtoCheck, ytoCheck] = SD_SlidingWindowF([1:1:length( locSignals(:,i) )], locSignals(:,i), 1, 40, 15, 70, 3, numDataPoints);
                        locSignalNorm =  ( locSignals(:,i) - min( locSignals(:,i) ) ) ./ ( max( locSignals(:,i) ) - min( locSignals(:,i) ) );
                        [dVeMono_dt, downstrokeMono, midSignalIndMono, signalDiffMonoPointsTime, upstrokeRiseMono, upstrokeBackMono, riseSignalIndMono, backSignalIndMono, PositiveDeflectionDuration_Mono, NegativeDeflectionDuration_Mono, baselineMono, signalWidth, timePeak, timeTrough] = EGM_Characterizing_Function(locSignalNorm, 1:140, 1, 1, EGMduration );
                    
                        localPPA = max(ytoCheck) - min(ytoCheck);
                        
                        global_PPASimulated = [global_PPASimulated; localPPA];
                        global_downstrokeSimulated = [global_downstrokeSimulated; downstrokeMono];
                        global_upstrokeRiseSimulated = [global_upstrokeRiseSimulated; upstrokeRiseMono];
                        global_upstrokeBackSimulated = [global_upstrokeBackSimulated; upstrokeBackMono];
                        global_PositiveDDSimulated = [global_PositiveDDSimulated; PositiveDeflectionDuration_Mono];
                        global_NegativeDDSimulated = [global_NegativeDDSimulated; NegativeDeflectionDuration_Mono];
                        global_baselineSimulated = [global_baselineSimulated; baselineMono];
                        global_SignalWSimulated = [global_SignalWSimulated; signalWidth];
                        global_EGMDurationSimulated = [global_EGMDurationSimulated; EGMduration];
                    
                        loc_normytoCheck = ( ytoCheck - min(ytoCheck) )./( max(ytoCheck) - min(ytoCheck) );
                        % loc_simDataCut = [loc_simDataCut, loc_normytoCheck];
                        % loc_simDataCutOG = [loc_simDataCutOG, ytoCheck];
                    
                        global_deflectionCountSimulated = [global_deflectionCountSimulated; DeflectionCount];
    
                        [indLoc,distsLoc] = findNearestNeighbors(ptCloudSurface,[loc_elecX(i), loc_elecY(i), loc_elecZ(i)],1);
                        global_AlldistanceToProjSimulated = [global_AlldistanceToProjSimulated;  distsLoc ];
                    
                        global_baselineMeanSDSWSimulated = [global_baselineMeanSDSWSimulated; baselineMeanSDSW];


                        % globalSimulated_dataToClassify = [ globalSimulated_dataToClassify; downstrokeMono, upstrokeRiseMono, upstrokeBackMono, signalWidth, baselineMono, EGMduration, distsLoc, localPPA, DeflectionCount];
                        loc_dataToClassify = [downstrokeMono, upstrokeRiseMono, upstrokeBackMono, signalWidth, baselineMono, EGMduration, distsLoc, localPPA, DeflectionCount];
                        loc_dataToClassify_Scaled = loc_dataToClassify;
                        for iteri = 1:length( loc_dataToClassify(1,:) )
                                loc_dataToClassify_Scaled(iteri) = ( loc_dataToClassify(:,iteri) - minScalingData(iteri) )./( maxScalingData(iteri) - minScalingData(iteri) );
                        end

                        [c,l] = wavedec(loc_normytoCheck,levelDec,"rbio4.4");
                        approx = appcoef(c,l,"rbio4.4");
                        [cd1,cd2,cd3,cd4] = detcoef(c,l,vecLevel);
                    
                        global_PCA1vect_sim_WD_Persistent = [global_PCA1vect_sim_WD_Persistent; coeffWD(:,1)'*cd4];
                        global_PCA2vect_sim_WD_Persistent = [global_PCA2vect_sim_WD_Persistent; coeffWD(:,2)'*cd4];
                        global_PCA3vect_sim_WD_Persistent = [global_PCA3vect_sim_WD_Persistent; coeffWD(:,3)'*cd4];
                        global_PCA4vect_sim_WD_Persistent = [global_PCA4vect_sim_WD_Persistent; coeffWD(:,4)'*cd4];
                        global_PCA5vect_sim_WD_Persistent = [global_PCA5vect_sim_WD_Persistent; coeffWD(:,5)'*cd4];

                        global_PCA1vect_sim_VAR_Persistent = [global_PCA1vect_sim_VAR_Persistent; coeffVAR(:,1)'*loc_dataToClassify_Scaled'];
                        global_PCA2vect_sim_VAR_Persistent = [global_PCA2vect_sim_VAR_Persistent; coeffVAR(:,2)'*loc_dataToClassify_Scaled'];
                        global_PCA3vect_sim_VAR_Persistent = [global_PCA3vect_sim_VAR_Persistent; coeffVAR(:,3)'*loc_dataToClassify_Scaled'];
                        global_PCA4vect_sim_VAR_Persistent = [global_PCA4vect_sim_VAR_Persistent; coeffVAR(:,4)'*loc_dataToClassify_Scaled'];
                        global_PCA5vect_sim_VAR_Persistent = [global_PCA5vect_sim_VAR_Persistent; coeffVAR(:,5)'*loc_dataToClassify_Scaled'];



                        %PCA-BASED
                        %WAVELET DECOMPOSITION
                        currIndex = idxWD(i);
                        dists_Patient_WD = DWD(i,currIndex);
                        locDists = [];
                        for j = 1:length(CWD(:,1))
                            locDists(j) = sqrt( (coeffWD(:,1)'*cd4 - CWD(j,1)).^2 + (coeffWD(:,2)'*cd4 - CWD(j,2)).^2 + (coeffWD(:,3)'*cd4 - CWD(j,3)).^2 + (coeffWD(:,4)'*cd4 - CWD(j,4)).^2 + (coeffWD(:,5)'*cd4 - CWD(j,5)).^2  );
                        end
                        dists_Simulated_common_WD = locDists(currIndex);
                        dists_Simulated_different_WD = min(locDists);

                        global_dists_Patient_WD = [global_dists_Patient_WD; dists_Patient_WD];
                        global_dists_Simulated_common_WD = [global_dists_Simulated_common_WD; dists_Simulated_common_WD];
                        global_dists_Simulated_different_WD = [global_dists_Simulated_different_WD; dists_Simulated_different_WD];

                        
                        %CHARACTERISTIC PARAMETERS
                        currIndex = idxPCAVAR(i);
                        dists_Patient_VAR = DPCAVAR(i,currIndex);
                        locDists = [];
                        for j = 1:length(CPCAVAR(:,1))
                            locDists(j) = sqrt( (coeffVAR(:,1)'*loc_dataToClassify_Scaled' - CPCAVAR(j,1)).^2 + (coeffVAR(:,2)'*loc_dataToClassify_Scaled' - CPCAVAR(j,2)).^2 + (coeffVAR(:,3)'*loc_dataToClassify_Scaled' - CPCAVAR(j,3)).^2 + (coeffVAR(:,4)'*loc_dataToClassify_Scaled' - CPCAVAR(j,4)).^2 + (coeffVAR(:,5)'*loc_dataToClassify_Scaled' - CPCAVAR(j,5)).^2  );
                        end
                        dists_Simulated_common_VAR = locDists(currIndex);
                        dists_Simulated_different_VAR = min(locDists);

                        global_dists_Patient_VAR = [global_dists_Patient_VAR; dists_Patient_VAR];
                        global_dists_Simulated_common_VAR = [global_dists_Simulated_common_VAR; dists_Simulated_common_VAR];
                        global_dists_Simulated_different_VAR = [global_dists_Simulated_different_VAR; dists_Simulated_different_VAR];
                    
                        

                        %SIGNAL-SIGNAL BASED
                        dVvec = diff(locSignals(:,i));
                        minVal =  min(dVvec);
                        indFind = find(dVvec == minVal);
                        indFind = indFind(1);
                        loc_actTimesSim(i) = indFind; 
                        correlation_coefficient = corrcoef( locSignalNorm_Matrix(i,:), reshape(loc_normytoCheck, [1,length(loc_normytoCheck)]) );
                        [xc,lags] = xcorr( locSignalNorm_Matrix(i,:), reshape(loc_normytoCheck, [1,length(loc_normytoCheck)]), 10,'normalized' );

                        loc_similarity_score(i) = correlation_coefficient(1, 2);
                        loc_xcorr(i) = max(xc);
                        loc_LAT_error(i) = (abs( actTimesCurrent(i) - loc_actTimesSim(i) )./actTimesCurrent(i))*100;

                        global_loc_similarity_score = [global_loc_similarity_score; correlation_coefficient(1, 2)];
                        global_loc_xcorr = [global_loc_xcorr; max(xc)];
                        global_loc_LAT_error = [global_loc_LAT_error; (abs( actTimesCurrent(i) - loc_actTimesSim(i) )./actTimesCurrent(i))*100];

                        each_LATs_Persistent = [ each_LATs_Persistent; actTimesCurrent(i) ];
                        global_LATs_Persistent = [ global_LATs_Persistent; loc_actTimesSim(i) ];

                        fType = [fType; tF];
                        patientNumber = [patientNumber; patientNum];
                        fPercentage = [fPercentage; fibrosisPercentage];
                        patternNumber = [patternNumber; patternNum];
                    
                    end

                    similarity_score_lumped = [similarity_score_lumped; mean(loc_similarity_score)];
                    LAT_error_lumped = [LAT_error_lumped; mean(loc_LAT_error)];
                    xcorr_lumped = [xcorr_lumped; mean(loc_xcorr)];
                    fType_lumped = [fType_lumped; tF];
                    patientNumber_lumped = [patientNumber_lumped; patientNum];
                    fPercentage_lumped = [fPercentage_lumped; fibrosisPercentage];
                    patternNumber_lumped = [patternNumber_lumped; patternNum];
                    
                else
                    fTypeNAN = [fTypeNAN; tF];
                    patientNumberNAN = [patientNumberNAN; patientNum];
                    fPercentageNAN = [fPercentageNAN; fibrosisPercentage];
                    patternNumberNAN = [patternNumberNAN; patternNum];
                end

            end

        end
    end
end


Compact10_Bin = [];
Diffuse10_Bin = [];
Interstitial10_Bin = [];
Patchy10_Bin = [];
Compact35_Bin = [];
Diffuse35_Bin = [];
Interstitial35_Bin = [];
Patchy35_Bin = [];
Compact60_Bin = [];
Diffuse60_Bin = [];
Interstitial60_Bin = [];
Patchy60_Bin = [];

Compact10_Bin_patNum = [];
Diffuse10_Bin_patNum = [];
Interstitial10_Bin_patNum = [];
Patchy10_Bin_patNum = [];
Compact35_Bin_patNum = [];
Diffuse35_Bin_patNum = [];
Interstitial35_Bin_patNum = [];
Patchy35_Bin_patNum = [];
Compact60_Bin_patNum = [];
Diffuse60_Bin_patNum = [];
Interstitial60_Bin_patNum = [];
Patchy60_Bin_patNum = [];

Compact10_Bin_pattNum = [];
Diffuse10_Bin_pattNum = [];
Interstitial10_Bin_pattNum = [];
Patchy10_Bin_pattNum = [];
Compact35_Bin_pattNum = [];
Diffuse35_Bin_pattNum = [];
Interstitial35_Bin_pattNum = [];
Patchy35_Bin_pattNum = [];
Compact60_Bin_pattNum = [];
Diffuse60_Bin_pattNum = [];
Interstitial60_Bin_pattNum = [];
Patchy60_Bin_pattNum = [];


Compact10_Bin_PC_WD_Persistent = [];
Diffuse10_Bin_PC_WD_Persistent = [];
Interstitial10_Bin_PC_WD_Persistent = [];
Patchy10_Bin_PC_WD_Persistent = [];
Compact35_Bin_PC_WD_Persistent = [];
Diffuse35_Bin_PC_WD_Persistent = [];
Interstitial35_Bin_PC_WD_Persistent = [];
Patchy35_Bin_PC_WD_Persistent = [];
Compact60_Bin_PC_WD_Persistent = [];
Diffuse60_Bin_PC_WD_Persistent = [];
Interstitial60_Bin_PC_WD_Persistent = [];
Patchy60_Bin_PC_WD_Persistent = [];


Compact10_Bin_PC_VAR_Persistent = [];
Diffuse10_Bin_PC_VAR_Persistent = [];
Interstitial10_Bin_PC_VAR_Persistent = [];
Patchy10_Bin_PC_VAR_Persistent = [];
Compact35_Bin_PC_VAR_Persistent = [];
Diffuse35_Bin_PC_VAR_Persistent = [];
Interstitial35_Bin_PC_VAR_Persistent = [];
Patchy35_Bin_PC_VAR_Persistent = [];
Compact60_Bin_PC_VAR_Persistent = [];
Diffuse60_Bin_PC_VAR_Persistent = [];
Interstitial60_Bin_PC_VAR_Persistent = [];
Patchy60_Bin_PC_VAR_Persistent = [];

Compact10_Bin_comparison = [];
Diffuse10_Bin_comparison = [];
Interstitial10_Bin_comparison = [];
Patchy10_Bin_comparison = [];
Compact35_Bin_comparison = [];
Diffuse35_Bin_comparison = [];
Interstitial35_Bin_comparison = [];
Patchy35_Bin_comparison = [];
Compact60_Bin_comparison = [];
Diffuse60_Bin_comparison = [];
Interstitial60_Bin_comparison = [];
Patchy60_Bin_comparison = [];


for i = 1:length(patternNumber)
    if fType(i) == 1 && fPercentage(i) == 10
        Compact10_Bin_patNum = [Compact10_Bin_patNum; patientNumber(i)];
        Compact10_Bin_pattNum = [Compact10_Bin_pattNum; patternNumber(i)];
        Compact10_Bin = [Compact10_Bin; global_downstrokeSimulated(i), global_upstrokeRiseSimulated(i), global_upstrokeBackSimulated(i), global_SignalWSimulated(i), global_baselineSimulated(i), global_EGMDurationSimulated(i), global_AlldistanceToProjSimulated(i), global_PPASimulated(i), global_deflectionCountSimulated(i) ];
        Compact10_Bin_PC_WD_Persistent = [ Compact10_Bin_PC_WD_Persistent; global_PCA1vect_sim_WD_Persistent(i), global_PCA2vect_sim_WD_Persistent(i), global_PCA3vect_sim_WD_Persistent(i), global_PCA4vect_sim_WD_Persistent(i), global_PCA5vect_sim_WD_Persistent(i) ];
        Compact10_Bin_PC_VAR_Persistent = [ Compact10_Bin_PC_VAR_Persistent; global_PCA1vect_sim_VAR_Persistent(i), global_PCA2vect_sim_VAR_Persistent(i), global_PCA3vect_sim_VAR_Persistent(i), global_PCA4vect_sim_VAR_Persistent(i), global_PCA5vect_sim_VAR_Persistent(i) ];
        Compact10_Bin_comparison = [Compact10_Bin_comparison; global_dists_Patient_WD(i), global_dists_Simulated_common_WD(i), global_dists_Simulated_different_WD(i), global_dists_Patient_VAR(i), global_dists_Simulated_common_VAR(i), global_dists_Simulated_different_VAR(i), global_loc_similarity_score(i), global_loc_xcorr(i), global_loc_LAT_error(i)];
    elseif fType(i) == 2 && fPercentage(i) == 10
        Diffuse10_Bin_patNum = [Diffuse10_Bin_patNum; patientNumber(i)];
        Diffuse10_Bin_pattNum = [Diffuse10_Bin_pattNum; patternNumber(i)];
        Diffuse10_Bin = [Diffuse10_Bin; global_downstrokeSimulated(i), global_upstrokeRiseSimulated(i), global_upstrokeBackSimulated(i), global_SignalWSimulated(i), global_baselineSimulated(i), global_EGMDurationSimulated(i), global_AlldistanceToProjSimulated(i), global_PPASimulated(i), global_deflectionCountSimulated(i) ];
        Diffuse10_Bin_PC_WD_Persistent = [ Diffuse10_Bin_PC_WD_Persistent; global_PCA1vect_sim_WD_Persistent(i), global_PCA2vect_sim_WD_Persistent(i), global_PCA3vect_sim_WD_Persistent(i), global_PCA4vect_sim_WD_Persistent(i), global_PCA5vect_sim_WD_Persistent(i) ];
        Diffuse10_Bin_PC_VAR_Persistent = [ Diffuse10_Bin_PC_VAR_Persistent; global_PCA1vect_sim_VAR_Persistent(i), global_PCA2vect_sim_VAR_Persistent(i), global_PCA3vect_sim_VAR_Persistent(i), global_PCA4vect_sim_VAR_Persistent(i), global_PCA5vect_sim_VAR_Persistent(i) ];
        Diffuse10_Bin_comparison = [Diffuse10_Bin_comparison; global_dists_Patient_WD(i), global_dists_Simulated_common_WD(i), global_dists_Simulated_different_WD(i), global_dists_Patient_VAR(i), global_dists_Simulated_common_VAR(i), global_dists_Simulated_different_VAR(i), global_loc_similarity_score(i), global_loc_xcorr(i), global_loc_LAT_error(i)];
    elseif fType(i) == 3 && fPercentage(i) == 10
        Interstitial10_Bin_patNum = [Interstitial10_Bin_patNum; patientNumber(i)];
        Interstitial10_Bin_pattNum = [Interstitial10_Bin_pattNum; patternNumber(i)];
        Interstitial10_Bin = [Interstitial10_Bin; global_downstrokeSimulated(i), global_upstrokeRiseSimulated(i), global_upstrokeBackSimulated(i), global_SignalWSimulated(i), global_baselineSimulated(i), global_EGMDurationSimulated(i), global_AlldistanceToProjSimulated(i), global_PPASimulated(i), global_deflectionCountSimulated(i) ];
        Interstitial10_Bin_PC_WD_Persistent = [ Interstitial10_Bin_PC_WD_Persistent; global_PCA1vect_sim_WD_Persistent(i), global_PCA2vect_sim_WD_Persistent(i), global_PCA3vect_sim_WD_Persistent(i), global_PCA4vect_sim_WD_Persistent(i), global_PCA5vect_sim_WD_Persistent(i) ];
        Interstitial10_Bin_PC_VAR_Persistent = [ Interstitial10_Bin_PC_VAR_Persistent; global_PCA1vect_sim_VAR_Persistent(i), global_PCA2vect_sim_VAR_Persistent(i), global_PCA3vect_sim_VAR_Persistent(i), global_PCA4vect_sim_VAR_Persistent(i), global_PCA5vect_sim_VAR_Persistent(i) ];
        Interstitial10_Bin_comparison = [Interstitial10_Bin_comparison; global_dists_Patient_WD(i), global_dists_Simulated_common_WD(i), global_dists_Simulated_different_WD(i), global_dists_Patient_VAR(i), global_dists_Simulated_common_VAR(i), global_dists_Simulated_different_VAR(i), global_loc_similarity_score(i), global_loc_xcorr(i), global_loc_LAT_error(i)];
    elseif fType(i) == 4 && fPercentage(i) == 10
        Patchy10_Bin_patNum = [Patchy10_Bin_patNum; patientNumber(i)];
        Patchy10_Bin_pattNum = [Patchy10_Bin_pattNum; patternNumber(i)];
        Patchy10_Bin = [Patchy10_Bin; global_downstrokeSimulated(i), global_upstrokeRiseSimulated(i), global_upstrokeBackSimulated(i), global_SignalWSimulated(i), global_baselineSimulated(i), global_EGMDurationSimulated(i), global_AlldistanceToProjSimulated(i), global_PPASimulated(i), global_deflectionCountSimulated(i) ];
        Patchy10_Bin_PC_WD_Persistent = [ Patchy10_Bin_PC_WD_Persistent; global_PCA1vect_sim_WD_Persistent(i), global_PCA2vect_sim_WD_Persistent(i), global_PCA3vect_sim_WD_Persistent(i), global_PCA4vect_sim_WD_Persistent(i), global_PCA5vect_sim_WD_Persistent(i) ];
        Patchy10_Bin_PC_VAR_Persistent = [ Patchy10_Bin_PC_VAR_Persistent; global_PCA1vect_sim_VAR_Persistent(i), global_PCA2vect_sim_VAR_Persistent(i), global_PCA3vect_sim_VAR_Persistent(i), global_PCA4vect_sim_VAR_Persistent(i), global_PCA5vect_sim_VAR_Persistent(i) ];
        Patchy10_Bin_comparison = [Patchy10_Bin_comparison; global_dists_Patient_WD(i), global_dists_Simulated_common_WD(i), global_dists_Simulated_different_WD(i), global_dists_Patient_VAR(i), global_dists_Simulated_common_VAR(i), global_dists_Simulated_different_VAR(i), global_loc_similarity_score(i), global_loc_xcorr(i), global_loc_LAT_error(i)];

    elseif fType(i) == 1 && fPercentage(i) == 35
        Compact35_Bin_patNum = [Compact35_Bin_patNum; patientNumber(i)];
        Compact35_Bin_pattNum = [Compact35_Bin_pattNum; patternNumber(i)];
        Compact35_Bin = [Compact35_Bin; global_downstrokeSimulated(i), global_upstrokeRiseSimulated(i), global_upstrokeBackSimulated(i), global_SignalWSimulated(i), global_baselineSimulated(i), global_EGMDurationSimulated(i), global_AlldistanceToProjSimulated(i), global_PPASimulated(i), global_deflectionCountSimulated(i) ];
        Compact35_Bin_PC_WD_Persistent = [ Compact35_Bin_PC_WD_Persistent; global_PCA1vect_sim_WD_Persistent(i), global_PCA2vect_sim_WD_Persistent(i), global_PCA3vect_sim_WD_Persistent(i), global_PCA4vect_sim_WD_Persistent(i), global_PCA5vect_sim_WD_Persistent(i) ];
        Compact35_Bin_PC_VAR_Persistent = [ Compact35_Bin_PC_VAR_Persistent; global_PCA1vect_sim_VAR_Persistent(i), global_PCA2vect_sim_VAR_Persistent(i), global_PCA3vect_sim_VAR_Persistent(i), global_PCA4vect_sim_VAR_Persistent(i), global_PCA5vect_sim_VAR_Persistent(i) ];
        Compact35_Bin_comparison = [Compact35_Bin_comparison; global_dists_Patient_WD(i), global_dists_Simulated_common_WD(i), global_dists_Simulated_different_WD(i), global_dists_Patient_VAR(i), global_dists_Simulated_common_VAR(i), global_dists_Simulated_different_VAR(i), global_loc_similarity_score(i), global_loc_xcorr(i), global_loc_LAT_error(i)];
    elseif fType(i) == 2 && fPercentage(i) == 35
        Diffuse35_Bin_patNum = [Diffuse35_Bin_patNum; patientNumber(i)];
        Diffuse35_Bin_pattNum = [Diffuse35_Bin_pattNum; patternNumber(i)];
        Diffuse35_Bin = [Diffuse35_Bin; global_downstrokeSimulated(i), global_upstrokeRiseSimulated(i), global_upstrokeBackSimulated(i), global_SignalWSimulated(i), global_baselineSimulated(i), global_EGMDurationSimulated(i), global_AlldistanceToProjSimulated(i), global_PPASimulated(i), global_deflectionCountSimulated(i) ];
        Diffuse35_Bin_PC_WD_Persistent = [ Diffuse35_Bin_PC_WD_Persistent; global_PCA1vect_sim_WD_Persistent(i), global_PCA2vect_sim_WD_Persistent(i), global_PCA3vect_sim_WD_Persistent(i), global_PCA4vect_sim_WD_Persistent(i), global_PCA5vect_sim_WD_Persistent(i) ];
        Diffuse35_Bin_PC_VAR_Persistent = [ Diffuse35_Bin_PC_VAR_Persistent; global_PCA1vect_sim_VAR_Persistent(i), global_PCA2vect_sim_VAR_Persistent(i), global_PCA3vect_sim_VAR_Persistent(i), global_PCA4vect_sim_VAR_Persistent(i), global_PCA5vect_sim_VAR_Persistent(i) ];
        Diffuse35_Bin_comparison = [Diffuse35_Bin_comparison; global_dists_Patient_WD(i), global_dists_Simulated_common_WD(i), global_dists_Simulated_different_WD(i), global_dists_Patient_VAR(i), global_dists_Simulated_common_VAR(i), global_dists_Simulated_different_VAR(i), global_loc_similarity_score(i), global_loc_xcorr(i), global_loc_LAT_error(i)];
    elseif fType(i) == 3 && fPercentage(i) == 35
        Interstitial35_Bin_patNum = [Interstitial35_Bin_patNum; patientNumber(i)];
        Interstitial35_Bin_pattNum = [Interstitial35_Bin_pattNum; patternNumber(i)];
        Interstitial35_Bin = [Interstitial35_Bin; global_downstrokeSimulated(i), global_upstrokeRiseSimulated(i), global_upstrokeBackSimulated(i), global_SignalWSimulated(i), global_baselineSimulated(i), global_EGMDurationSimulated(i), global_AlldistanceToProjSimulated(i), global_PPASimulated(i), global_deflectionCountSimulated(i) ];
        Interstitial35_Bin_PC_WD_Persistent = [ Interstitial35_Bin_PC_WD_Persistent; global_PCA1vect_sim_WD_Persistent(i), global_PCA2vect_sim_WD_Persistent(i), global_PCA3vect_sim_WD_Persistent(i), global_PCA4vect_sim_WD_Persistent(i), global_PCA5vect_sim_WD_Persistent(i) ];
        Interstitial35_Bin_PC_VAR_Persistent = [ Interstitial35_Bin_PC_VAR_Persistent; global_PCA1vect_sim_VAR_Persistent(i), global_PCA2vect_sim_VAR_Persistent(i), global_PCA3vect_sim_VAR_Persistent(i), global_PCA4vect_sim_VAR_Persistent(i), global_PCA5vect_sim_VAR_Persistent(i) ];
        Interstitial35_Bin_comparison = [Interstitial35_Bin_comparison; global_dists_Patient_WD(i), global_dists_Simulated_common_WD(i), global_dists_Simulated_different_WD(i), global_dists_Patient_VAR(i), global_dists_Simulated_common_VAR(i), global_dists_Simulated_different_VAR(i), global_loc_similarity_score(i), global_loc_xcorr(i), global_loc_LAT_error(i)];
    elseif fType(i) == 4 && fPercentage(i) == 35
        Patchy35_Bin_patNum = [Patchy35_Bin_patNum; patientNumber(i)];
        Patchy35_Bin_pattNum = [Patchy35_Bin_pattNum; patternNumber(i)];
        Patchy35_Bin = [Patchy35_Bin; global_downstrokeSimulated(i), global_upstrokeRiseSimulated(i), global_upstrokeBackSimulated(i), global_SignalWSimulated(i), global_baselineSimulated(i), global_EGMDurationSimulated(i), global_AlldistanceToProjSimulated(i), global_PPASimulated(i), global_deflectionCountSimulated(i) ];
        Patchy35_Bin_PC_WD_Persistent = [ Patchy35_Bin_PC_WD_Persistent; global_PCA1vect_sim_WD_Persistent(i), global_PCA2vect_sim_WD_Persistent(i), global_PCA3vect_sim_WD_Persistent(i), global_PCA4vect_sim_WD_Persistent(i), global_PCA5vect_sim_WD_Persistent(i) ];
        Patchy35_Bin_PC_VAR_Persistent = [ Patchy35_Bin_PC_VAR_Persistent; global_PCA1vect_sim_VAR_Persistent(i), global_PCA2vect_sim_VAR_Persistent(i), global_PCA3vect_sim_VAR_Persistent(i), global_PCA4vect_sim_VAR_Persistent(i), global_PCA5vect_sim_VAR_Persistent(i) ];
        Patchy35_Bin_comparison = [Patchy35_Bin_comparison; global_dists_Patient_WD(i), global_dists_Simulated_common_WD(i), global_dists_Simulated_different_WD(i), global_dists_Patient_VAR(i), global_dists_Simulated_common_VAR(i), global_dists_Simulated_different_VAR(i), global_loc_similarity_score(i), global_loc_xcorr(i), global_loc_LAT_error(i)];

    elseif fType(i) == 1 && fPercentage(i) == 60
        Compact60_Bin_patNum = [Compact60_Bin_patNum; patientNumber(i)];
        Compact60_Bin_pattNum = [Compact60_Bin_pattNum; patternNumber(i)];
        Compact60_Bin = [Compact60_Bin; global_downstrokeSimulated(i), global_upstrokeRiseSimulated(i), global_upstrokeBackSimulated(i), global_SignalWSimulated(i), global_baselineSimulated(i), global_EGMDurationSimulated(i), global_AlldistanceToProjSimulated(i), global_PPASimulated(i), global_deflectionCountSimulated(i) ];
        Compact60_Bin_PC_WD_Persistent = [ Compact60_Bin_PC_WD_Persistent; global_PCA1vect_sim_WD_Persistent(i), global_PCA2vect_sim_WD_Persistent(i), global_PCA3vect_sim_WD_Persistent(i), global_PCA4vect_sim_WD_Persistent(i), global_PCA5vect_sim_WD_Persistent(i) ];
        Compact60_Bin_PC_VAR_Persistent = [ Compact60_Bin_PC_VAR_Persistent; global_PCA1vect_sim_VAR_Persistent(i), global_PCA2vect_sim_VAR_Persistent(i), global_PCA3vect_sim_VAR_Persistent(i), global_PCA4vect_sim_VAR_Persistent(i), global_PCA5vect_sim_VAR_Persistent(i) ];
        Compact60_Bin_comparison = [Compact60_Bin_comparison; global_dists_Patient_WD(i), global_dists_Simulated_common_WD(i), global_dists_Simulated_different_WD(i), global_dists_Patient_VAR(i), global_dists_Simulated_common_VAR(i), global_dists_Simulated_different_VAR(i), global_loc_similarity_score(i), global_loc_xcorr(i), global_loc_LAT_error(i)];
    elseif fType(i) == 2 && fPercentage(i) == 60
        Diffuse60_Bin_patNum = [Diffuse60_Bin_patNum; patientNumber(i)];
        Diffuse60_Bin_pattNum = [Diffuse60_Bin_pattNum; patternNumber(i)];
        Diffuse60_Bin = [Diffuse60_Bin; global_downstrokeSimulated(i), global_upstrokeRiseSimulated(i), global_upstrokeBackSimulated(i), global_SignalWSimulated(i), global_baselineSimulated(i), global_EGMDurationSimulated(i), global_AlldistanceToProjSimulated(i), global_PPASimulated(i), global_deflectionCountSimulated(i) ];
        Diffuse60_Bin_PC_WD_Persistent = [ Diffuse60_Bin_PC_WD_Persistent; global_PCA1vect_sim_WD_Persistent(i), global_PCA2vect_sim_WD_Persistent(i), global_PCA3vect_sim_WD_Persistent(i), global_PCA4vect_sim_WD_Persistent(i), global_PCA5vect_sim_WD_Persistent(i) ];
        Diffuse60_Bin_PC_VAR_Persistent = [ Diffuse60_Bin_PC_VAR_Persistent; global_PCA1vect_sim_VAR_Persistent(i), global_PCA2vect_sim_VAR_Persistent(i), global_PCA3vect_sim_VAR_Persistent(i), global_PCA4vect_sim_VAR_Persistent(i), global_PCA5vect_sim_VAR_Persistent(i) ];
        Diffuse60_Bin_comparison = [Diffuse60_Bin_comparison; global_dists_Patient_WD(i), global_dists_Simulated_common_WD(i), global_dists_Simulated_different_WD(i), global_dists_Patient_VAR(i), global_dists_Simulated_common_VAR(i), global_dists_Simulated_different_VAR(i), global_loc_similarity_score(i), global_loc_xcorr(i), global_loc_LAT_error(i)];
    elseif fType(i) == 3 && fPercentage(i) == 60
        Interstitial60_Bin_patNum = [Interstitial60_Bin_patNum; patientNumber(i)];
        Interstitial60_Bin_pattNum = [Interstitial60_Bin_pattNum; patternNumber(i)];
        Interstitial60_Bin = [Interstitial60_Bin; global_downstrokeSimulated(i), global_upstrokeRiseSimulated(i), global_upstrokeBackSimulated(i), global_SignalWSimulated(i), global_baselineSimulated(i), global_EGMDurationSimulated(i), global_AlldistanceToProjSimulated(i), global_PPASimulated(i), global_deflectionCountSimulated(i) ];
        Interstitial60_Bin_PC_WD_Persistent = [ Interstitial60_Bin_PC_WD_Persistent; global_PCA1vect_sim_WD_Persistent(i), global_PCA2vect_sim_WD_Persistent(i), global_PCA3vect_sim_WD_Persistent(i), global_PCA4vect_sim_WD_Persistent(i), global_PCA5vect_sim_WD_Persistent(i) ];
        Interstitial60_Bin_PC_VAR_Persistent = [ Interstitial60_Bin_PC_VAR_Persistent; global_PCA1vect_sim_VAR_Persistent(i), global_PCA2vect_sim_VAR_Persistent(i), global_PCA3vect_sim_VAR_Persistent(i), global_PCA4vect_sim_VAR_Persistent(i), global_PCA5vect_sim_VAR_Persistent(i) ];
        Interstitial60_Bin_comparison = [Interstitial60_Bin_comparison; global_dists_Patient_WD(i), global_dists_Simulated_common_WD(i), global_dists_Simulated_different_WD(i), global_dists_Patient_VAR(i), global_dists_Simulated_common_VAR(i), global_dists_Simulated_different_VAR(i), global_loc_similarity_score(i), global_loc_xcorr(i), global_loc_LAT_error(i)];
    elseif fType(i) == 4 && fPercentage(i) == 60
        Patchy60_Bin_patNum = [Patchy60_Bin_patNum; patientNumber(i)];
        Patchy60_Bin_pattNum = [Patchy60_Bin_pattNum; patternNumber(i)];
        Patchy60_Bin = [Patchy60_Bin; global_downstrokeSimulated(i), global_upstrokeRiseSimulated(i), global_upstrokeBackSimulated(i), global_SignalWSimulated(i), global_baselineSimulated(i), global_EGMDurationSimulated(i), global_AlldistanceToProjSimulated(i), global_PPASimulated(i), global_deflectionCountSimulated(i) ];
        Patchy60_Bin_PC_WD_Persistent = [ Patchy60_Bin_PC_WD_Persistent; global_PCA1vect_sim_WD_Persistent(i), global_PCA2vect_sim_WD_Persistent(i), global_PCA3vect_sim_WD_Persistent(i), global_PCA4vect_sim_WD_Persistent(i), global_PCA5vect_sim_WD_Persistent(i) ];
        Patchy60_Bin_PC_VAR_Persistent = [ Patchy60_Bin_PC_VAR_Persistent; global_PCA1vect_sim_VAR_Persistent(i), global_PCA2vect_sim_VAR_Persistent(i), global_PCA3vect_sim_VAR_Persistent(i), global_PCA4vect_sim_VAR_Persistent(i), global_PCA5vect_sim_VAR_Persistent(i) ];
        Patchy60_Bin_comparison = [Patchy60_Bin_comparison; global_dists_Patient_WD(i), global_dists_Simulated_common_WD(i), global_dists_Simulated_different_WD(i), global_dists_Patient_VAR(i), global_dists_Simulated_common_VAR(i), global_dists_Simulated_different_VAR(i), global_loc_similarity_score(i), global_loc_xcorr(i), global_loc_LAT_error(i)];
    end
end





toc




nonworking = [patientNumberNAN, fTypeNAN, fPercentageNAN, patternNumberNAN];



%Idea: map in space PC's of averages maybe of characteristics of 12
%categories to see where they lie in space in comparison to the other
%categories, like to see which categories are more similar to each other,
%i.e. is degree more relevant or type in terms of similarities.


%Average characteristics for each category
table_all_categories_characteristics = {};
for i = 1:length(Compact10_Bin(1,:))
    table_all_categories_characteristics{1,i} =  createMeanStdForCell(Compact10_Bin(:,i), 1:10, 1:10, 2 );
end
for i = 1:length(Compact35_Bin(1,:))
    table_all_categories_characteristics{2,i} =  createMeanStdForCell(Compact35_Bin(:,i), 1:10, 1:10, 2 );
end
for i = 1:length(Compact60_Bin(1,:))
    table_all_categories_characteristics{3,i} =  createMeanStdForCell(Compact60_Bin(:,i), 1:10, 1:10, 2 );
end

for i = 1:length(Diffuse10_Bin(1,:))
    table_all_categories_characteristics{4,i} =  createMeanStdForCell(Diffuse10_Bin(:,i), 1:10, 1:10, 2 );
end
for i = 1:length(Diffuse35_Bin(1,:))
    table_all_categories_characteristics{5,i} =  createMeanStdForCell(Diffuse35_Bin(:,i), 1:10, 1:10, 2 );
end
for i = 1:length(Diffuse60_Bin(1,:))
    table_all_categories_characteristics{6,i} =  createMeanStdForCell(Diffuse60_Bin(:,i), 1:10, 1:10, 2 );
end

for i = 1:length(Interstitial10_Bin(1,:))
    table_all_categories_characteristics{7,i} =  createMeanStdForCell(Interstitial10_Bin(:,i), 1:10, 1:10, 2 );
end
for i = 1:length(Interstitial35_Bin(1,:))
    table_all_categories_characteristics{8,i} =  createMeanStdForCell(Interstitial35_Bin(:,i), 1:10, 1:10, 2 );
end
for i = 1:length(Interstitial60_Bin(1,:))
    table_all_categories_characteristics{9,i} =  createMeanStdForCell(Interstitial60_Bin(:,i), 1:10, 1:10, 2 );
end

for i = 1:length(Patchy10_Bin(1,:))
    table_all_categories_characteristics{10,i} =  createMeanStdForCell(Patchy10_Bin(:,i), 1:10, 1:10, 2 );
end
for i = 1:length(Patchy35_Bin(1,:))
    table_all_categories_characteristics{11,i} =  createMeanStdForCell(Patchy35_Bin(:,i), 1:10, 1:10, 2 );
end
for i = 1:length(Patchy60_Bin(1,:))
    table_all_categories_characteristics{12,i} =  createMeanStdForCell(Patchy60_Bin(:,i), 1:10, 1:10, 2 );
end

table_all_categories_characteristics = cell2table(table_all_categories_characteristics);





%Variability between patients
table_variability_between_patients = {};
dummyLoc = 1;
for i = [6,8,9]
    table_variability_between_patients{1,dummyLoc} =  createMeanStdForCell(Compact10_Bin(:,i), Compact10_Bin_patNum, 6, 2 );
    table_variability_between_patients{1,dummyLoc+3} =  createMeanStdForCell(Compact10_Bin(:,i), Compact10_Bin_patNum, 7, 2 );
    table_variability_between_patients{1,dummyLoc+6} =  createMeanStdForCell(Compact10_Bin(:,i), Compact10_Bin_patNum, 8, 2 );
    table_variability_between_patients{1,dummyLoc+9} =  createMeanStdForCell(Compact10_Bin(:,i), Compact10_Bin_patNum, 9, 2 );
    dummyLoc = dummyLoc + 1;
end
dummyLoc = 1;
for i = [6,8,9]
    table_variability_between_patients{2,dummyLoc} =  createMeanStdForCell(Compact35_Bin(:,i), Compact35_Bin_patNum, 6, 2 );
    table_variability_between_patients{2,dummyLoc+3} =  createMeanStdForCell(Compact35_Bin(:,i), Compact35_Bin_patNum, 7, 2 );
    table_variability_between_patients{2,dummyLoc+6} =  createMeanStdForCell(Compact35_Bin(:,i), Compact35_Bin_patNum, 8, 2 );
    table_variability_between_patients{2,dummyLoc+9} =  createMeanStdForCell(Compact35_Bin(:,i), Compact35_Bin_patNum, 9, 2 );
    dummyLoc = dummyLoc + 1;
end
dummyLoc = 1;
for i = [6,8,9]
    table_variability_between_patients{3,dummyLoc} =  createMeanStdForCell(Compact60_Bin(:,i), Compact60_Bin_patNum, 6, 2 );
    table_variability_between_patients{3,dummyLoc+3} =  createMeanStdForCell(Compact60_Bin(:,i), Compact60_Bin_patNum, 7, 2 );
    table_variability_between_patients{3,dummyLoc+6} =  createMeanStdForCell(Compact60_Bin(:,i), Compact60_Bin_patNum, 8, 2 );
    table_variability_between_patients{3,dummyLoc+9} =  createMeanStdForCell(Compact60_Bin(:,i), Compact60_Bin_patNum, 9, 2 );
    dummyLoc = dummyLoc + 1;
end

dummyLoc = 1;
for i = [6,8,9]
    table_variability_between_patients{4,dummyLoc} =  createMeanStdForCell(Diffuse10_Bin(:,i), Diffuse10_Bin_patNum, 6, 2 );
    table_variability_between_patients{4,dummyLoc+3} =  createMeanStdForCell(Diffuse10_Bin(:,i), Diffuse10_Bin_patNum, 7, 2 );
    table_variability_between_patients{4,dummyLoc+6} =  createMeanStdForCell(Diffuse10_Bin(:,i), Diffuse10_Bin_patNum, 8, 2 );
    table_variability_between_patients{4,dummyLoc+9} =  createMeanStdForCell(Diffuse10_Bin(:,i), Diffuse10_Bin_patNum, 9, 2 );
    dummyLoc = dummyLoc + 1;
end
dummyLoc = 1;
for i = [6,8,9]
    table_variability_between_patients{5,dummyLoc} =  createMeanStdForCell(Diffuse35_Bin(:,i), Diffuse35_Bin_patNum, 6, 2 );
    table_variability_between_patients{5,dummyLoc+3} =  createMeanStdForCell(Diffuse35_Bin(:,i), Diffuse35_Bin_patNum, 7, 2 );
    table_variability_between_patients{5,dummyLoc+6} =  createMeanStdForCell(Diffuse35_Bin(:,i), Diffuse35_Bin_patNum, 8, 2 );
    table_variability_between_patients{5,dummyLoc+9} =  createMeanStdForCell(Diffuse35_Bin(:,i), Diffuse35_Bin_patNum, 9, 2 );
    dummyLoc = dummyLoc + 1;
end
dummyLoc = 1;
for i = [6,8,9]
    table_variability_between_patients{6,dummyLoc} =  createMeanStdForCell(Diffuse60_Bin(:,i), Diffuse60_Bin_patNum, 6, 2 );
    table_variability_between_patients{6,dummyLoc+3} =  createMeanStdForCell(Diffuse60_Bin(:,i), Diffuse60_Bin_patNum, 7, 2 );
    table_variability_between_patients{6,dummyLoc+6} =  createMeanStdForCell(Diffuse60_Bin(:,i), Diffuse60_Bin_patNum, 8, 2 );
    table_variability_between_patients{6,dummyLoc+9} =  createMeanStdForCell(Diffuse60_Bin(:,i), Diffuse60_Bin_patNum, 9, 2 );
    dummyLoc = dummyLoc + 1;
end

dummyLoc = 1;
for i = [6,8,9]
    table_variability_between_patients{7,dummyLoc} =  createMeanStdForCell(Interstitial10_Bin(:,i), Interstitial10_Bin_patNum, 6, 2 );
    table_variability_between_patients{7,dummyLoc+3} =  createMeanStdForCell(Interstitial10_Bin(:,i), Interstitial10_Bin_patNum, 7, 2 );
    table_variability_between_patients{7,dummyLoc+6} =  createMeanStdForCell(Interstitial10_Bin(:,i), Interstitial10_Bin_patNum, 8, 2 );
    table_variability_between_patients{7,dummyLoc+9} =  createMeanStdForCell(Interstitial10_Bin(:,i), Interstitial10_Bin_patNum, 9, 2 );
    dummyLoc = dummyLoc + 1;
end
dummyLoc = 1;
for i = [6,8,9]
    table_variability_between_patients{8,dummyLoc} =  createMeanStdForCell(Interstitial35_Bin(:,i), Interstitial35_Bin_patNum, 6, 2 );
    table_variability_between_patients{8,dummyLoc+3} =  createMeanStdForCell(Interstitial35_Bin(:,i), Interstitial35_Bin_patNum, 7, 2 );
    table_variability_between_patients{8,dummyLoc+6} =  createMeanStdForCell(Interstitial35_Bin(:,i), Interstitial35_Bin_patNum, 8, 2 );
    table_variability_between_patients{8,dummyLoc+9} =  createMeanStdForCell(Interstitial35_Bin(:,i), Interstitial35_Bin_patNum, 9, 2 );
    dummyLoc = dummyLoc + 1;
end
dummyLoc = 1;
for i = [6,8,9]
    table_variability_between_patients{9,dummyLoc} =  createMeanStdForCell(Interstitial60_Bin(:,i), Interstitial60_Bin_patNum, 6, 2 );
    table_variability_between_patients{9,dummyLoc+3} =  createMeanStdForCell(Interstitial60_Bin(:,i), Interstitial60_Bin_patNum, 7, 2 );
    table_variability_between_patients{9,dummyLoc+6} =  createMeanStdForCell(Interstitial60_Bin(:,i), Interstitial60_Bin_patNum, 8, 2 );
    table_variability_between_patients{9,dummyLoc+9} =  createMeanStdForCell(Interstitial60_Bin(:,i), Interstitial60_Bin_patNum, 9, 2 );
    dummyLoc = dummyLoc + 1;
end

dummyLoc = 1;
for i = [6,8,9]
    table_variability_between_patients{10,dummyLoc} =  createMeanStdForCell(Patchy10_Bin(:,i), Patchy10_Bin_patNum, 6, 2 );
    table_variability_between_patients{10,dummyLoc+3} =  createMeanStdForCell(Patchy10_Bin(:,i), Patchy10_Bin_patNum, 7, 2 );
    table_variability_between_patients{10,dummyLoc+6} =  createMeanStdForCell(Patchy10_Bin(:,i), Patchy10_Bin_patNum, 8, 2 );
    table_variability_between_patients{10,dummyLoc+9} =  createMeanStdForCell(Patchy10_Bin(:,i), Patchy10_Bin_patNum, 9, 2 );
    dummyLoc = dummyLoc + 1;
end
dummyLoc = 1;
for i = [6,8,9]
    table_variability_between_patients{11,dummyLoc} =  createMeanStdForCell(Patchy35_Bin(:,i), Patchy35_Bin_patNum, 6, 2 );
    table_variability_between_patients{11,dummyLoc+3} =  createMeanStdForCell(Patchy35_Bin(:,i), Patchy35_Bin_patNum, 7, 2 );
    table_variability_between_patients{11,dummyLoc+6} =  createMeanStdForCell(Patchy35_Bin(:,i), Patchy35_Bin_patNum, 8, 2 );
    table_variability_between_patients{11,dummyLoc+9} =  createMeanStdForCell(Patchy35_Bin(:,i), Patchy35_Bin_patNum, 9, 2 );
    dummyLoc = dummyLoc + 1;
end
dummyLoc = 1;
for i = [6,8,9]
    table_variability_between_patients{12,dummyLoc} =  createMeanStdForCell(Patchy60_Bin(:,i), Patchy60_Bin_patNum, 6, 2 );
    table_variability_between_patients{12,dummyLoc+3} =  createMeanStdForCell(Patchy60_Bin(:,i), Patchy60_Bin_patNum, 7, 2 );
    table_variability_between_patients{12,dummyLoc+6} =  createMeanStdForCell(Patchy60_Bin(:,i), Patchy60_Bin_patNum, 8, 2 );
    table_variability_between_patients{12,dummyLoc+9} =  createMeanStdForCell(Patchy60_Bin(:,i), Patchy60_Bin_patNum, 9, 2 );
    dummyLoc = dummyLoc + 1;
end

table_variability_between_patients = cell2table(table_variability_between_patients);






%Variability between different densities but the same type of fibrosis
table_variability_between_densitites_and_type = {};

dummyLoc = 1;
for i = [6,8,9]
    table_variability_between_densitites_and_type{1,dummyLoc} =  createMeanStdForCell(Compact10_Bin(:,i), Compact10_Bin_patNum, Compact10_Bin_patNum, 2 );
    table_variability_between_densitites_and_type{1,dummyLoc+3} =  createMeanStdForCell(Compact35_Bin(:,i), Compact35_Bin_patNum, Compact35_Bin_patNum, 2 );
    table_variability_between_densitites_and_type{1,dummyLoc+6} =  createMeanStdForCell(Compact60_Bin(:,i), Compact60_Bin_patNum, Compact60_Bin_patNum, 2 );
    dummyLoc = dummyLoc + 1;
end

dummyLoc = 1;
for i = [6,8,9]
    table_variability_between_densitites_and_type{2,dummyLoc} =  createMeanStdForCell(Diffuse10_Bin(:,i), Diffuse10_Bin_patNum, Diffuse10_Bin_patNum, 2 );
    table_variability_between_densitites_and_type{2,dummyLoc+3} =  createMeanStdForCell(Diffuse35_Bin(:,i), Diffuse35_Bin_patNum, Diffuse35_Bin_patNum, 2 );
    table_variability_between_densitites_and_type{2,dummyLoc+6} =  createMeanStdForCell(Diffuse60_Bin(:,i), Diffuse60_Bin_patNum, Diffuse60_Bin_patNum, 2 );
    dummyLoc = dummyLoc + 1;
end

dummyLoc = 1;
for i = [6,8,9]
    table_variability_between_densitites_and_type{3,dummyLoc} =  createMeanStdForCell(Interstitial10_Bin(:,i), Interstitial10_Bin_patNum, Interstitial10_Bin_patNum, 2 );
    table_variability_between_densitites_and_type{3,dummyLoc+3} =  createMeanStdForCell(Interstitial35_Bin(:,i), Interstitial35_Bin_patNum, Interstitial35_Bin_patNum, 2 );
    table_variability_between_densitites_and_type{3,dummyLoc+6} =  createMeanStdForCell(Interstitial60_Bin(:,i), Interstitial60_Bin_patNum, Interstitial60_Bin_patNum, 2 );
    dummyLoc = dummyLoc + 1;
end

dummyLoc = 1;
for i = [6,8,9]
    table_variability_between_densitites_and_type{4,dummyLoc} =  createMeanStdForCell(Patchy10_Bin(:,i), Patchy10_Bin_patNum, Patchy10_Bin_patNum, 2 );
    table_variability_between_densitites_and_type{4,dummyLoc+3} =  createMeanStdForCell(Patchy35_Bin(:,i), Patchy35_Bin_patNum, Patchy35_Bin_patNum, 2 );
    table_variability_between_densitites_and_type{4,dummyLoc+6} =  createMeanStdForCell(Patchy60_Bin(:,i), Patchy60_Bin_patNum, Patchy60_Bin_patNum, 2 );
    dummyLoc = dummyLoc + 1;
end

table_variability_between_densitites_and_type = cell2table(table_variability_between_densitites_and_type);






%patient characteristics
finalAllPatient_Data = [ each_downstrokeSimulated, each_upstrokeRiseSimulated, each_upstrokeBackSimulated, each_SignalWSimulated, each_baselineSimulated, each_EGMDurationSimulated, each_AlldistanceToProjSimulated, each_PPASimulated, each_deflectionCountSimulated];

table_PersistentPat_Data = {};
for i = 1:length(finalAllPatient_Data(1,:))
    table_PersistentPat_Data{i,1} = createMeanStdForCell(finalAllPatient_Data(:,i), patient_lumped, 6, 2 );
    table_PersistentPat_Data{i,2} = createMeanStdForCell(finalAllPatient_Data(:,i), patient_lumped, 7, 2 );
    table_PersistentPat_Data{i,3} = createMeanStdForCell(finalAllPatient_Data(:,i), patient_lumped, 8, 2 );
    table_PersistentPat_Data{i,4} = createMeanStdForCell(finalAllPatient_Data(:,i), patient_lumped, 9, 2 );
end
table_PersistentPat_Data = cell2table(table_PersistentPat_Data);


%% PLOTTING AND DATA GATHERING
loc_actTimes_2  = expDB.dataSet{2}.map.continuousMapping.revisedLat(:,1);
loc_actTimes_3  = expDB.dataSet{3}.map.continuousMapping.revisedLat(:,1);

loc_actTimes_6  = expDB.dataSet{6}.map.continuousMapping.revisedLat(:,1);
loc_actTimes_7  = expDB.dataSet{7}.map.continuousMapping.revisedLat(:,1);
loc_actTimes_8  = expDB.dataSet{8}.map.continuousMapping.revisedLat(:,1);
loc_actTimes_9  = expDB.dataSet{9}.map.continuousMapping.revisedLat(:,1);



meshX_Paroxysmal2 = expDB.dataSet{2}.geometry.verticesTrimmed(:,1)./10; meshY_Paroxysmal2 = expDB.dataSet{2}.geometry.verticesTrimmed(:,2)./10; meshZ_Paroxysmal2 = expDB.dataSet{2}.geometry.verticesTrimmed(:,3)./10;
meshX_Paroxysmal3 = expDB.dataSet{3}.geometry.verticesTrimmed(:,1)./10; meshY_Paroxysmal3 = expDB.dataSet{3}.geometry.verticesTrimmed(:,2)./10; meshZ_Paroxysmal3 = expDB.dataSet{3}.geometry.verticesTrimmed(:,3)./10;

meshX_Persistent6 = expDB.dataSet{6}.geometry.verticesTrimmed(:,1)./10; meshY_Persistent6 = expDB.dataSet{6}.geometry.verticesTrimmed(:,2)./10; meshZ_Persistent6 = expDB.dataSet{6}.geometry.verticesTrimmed(:,3)./10;
meshX_Persistent7 = expDB.dataSet{7}.geometry.verticesTrimmed(:,1)./10; meshY_Persistent7 = expDB.dataSet{7}.geometry.verticesTrimmed(:,2)./10; meshZ_Persistent7 = expDB.dataSet{7}.geometry.verticesTrimmed(:,3)./10;
meshX_Persistent8 = expDB.dataSet{8}.geometry.verticesTrimmed(:,1)./10; meshY_Persistent8 = expDB.dataSet{8}.geometry.verticesTrimmed(:,2)./10; meshZ_Persistent8 = expDB.dataSet{8}.geometry.verticesTrimmed(:,3)./10;
meshX_Persistent9 = expDB.dataSet{9}.geometry.verticesTrimmed(:,1)./10; meshY_Persistent9 = expDB.dataSet{9}.geometry.verticesTrimmed(:,2)./10; meshZ_Persistent9 = expDB.dataSet{9}.geometry.verticesTrimmed(:,3)./10;



% colorVec = { sscanf('0000ff','%2x%2x%2x',[1 3])/255, sscanf('3aa0ff','%2x%2x%2x',[1 3])/255, sscanf('5df2ff','%2x%2x%2x',[1 3])/255,...
%              sscanf('36aa74','%2x%2x%2x',[1 3])/255, sscanf('00ff00','%2x%2x%2x',[1 3])/255, sscanf('fff121','%2x%2x%2x',[1 3])/255,...
%              sscanf('ffaa00','%2x%2x%2x',[1 3])/255, sscanf('ff5500','%2x%2x%2x',[1 3])/255, sscanf('ff0000','%2x%2x%2x',[1 3])/255 };


colorVec2 = { sscanf('142c5f','%2x%2x%2x',[1 3])/255, sscanf('124e63','%2x%2x%2x',[1 3])/255,...
             sscanf('2e6a57','%2x%2x%2x',[1 3])/255, sscanf('a36f63','%2x%2x%2x',[1 3])/255,...
             sscanf('f6a895','%2x%2x%2x',[1 3])/255, sscanf('af8094','%2x%2x%2x',[1 3])/255, sscanf('ffbad7','%2x%2x%2x',[1 3])/255 };

colorVec = { sscanf('142c5f','%2x%2x%2x',[1 3])/255, sscanf('114d81','%2x%2x%2x',[1 3])/255, sscanf('124e63','%2x%2x%2x',[1 3])/255,...
             sscanf('1b4034','%2x%2x%2x',[1 3])/255, sscanf('2e6a57','%2x%2x%2x',[1 3])/255, sscanf('a36f63','%2x%2x%2x',[1 3])/255,...
             sscanf('f6a895','%2x%2x%2x',[1 3])/255, sscanf('af8094','%2x%2x%2x',[1 3])/255, sscanf('ffbad7','%2x%2x%2x',[1 3])/255 };





%Patient 2
projLAT_patient2 = projLATs_to_Surface(elecX_Paroxysmal2, elecY_Paroxysmal2, elecZ_Paroxysmal2, meshX_Paroxysmal2, meshY_Paroxysmal2, meshZ_Paroxysmal2, loc_actTimes_2);
numBins = 9;
[binEdges2, binIndices_patient2, remainingIndices2] = partitionToBins_Isochrones(projLAT_patient2, numBins, [] );

figure
plot3(meshX_Paroxysmal2(binIndices_patient2{1}), meshY_Paroxysmal2(binIndices_patient2{1}), meshZ_Paroxysmal2(binIndices_patient2{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Paroxysmal2(binIndices_patient2{2}), meshY_Paroxysmal2(binIndices_patient2{2}), meshZ_Paroxysmal2(binIndices_patient2{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Paroxysmal2(binIndices_patient2{3}), meshY_Paroxysmal2(binIndices_patient2{3}), meshZ_Paroxysmal2(binIndices_patient2{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Paroxysmal2(binIndices_patient2{4}), meshY_Paroxysmal2(binIndices_patient2{4}), meshZ_Paroxysmal2(binIndices_patient2{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Paroxysmal2(binIndices_patient2{5}), meshY_Paroxysmal2(binIndices_patient2{5}), meshZ_Paroxysmal2(binIndices_patient2{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Paroxysmal2(binIndices_patient2{6}), meshY_Paroxysmal2(binIndices_patient2{6}), meshZ_Paroxysmal2(binIndices_patient2{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Paroxysmal2(binIndices_patient2{7}), meshY_Paroxysmal2(binIndices_patient2{7}), meshZ_Paroxysmal2(binIndices_patient2{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Paroxysmal2(binIndices_patient2{8}), meshY_Paroxysmal2(binIndices_patient2{8}), meshZ_Paroxysmal2(binIndices_patient2{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Paroxysmal2(binIndices_patient2{9}), meshY_Paroxysmal2(binIndices_patient2{9}), meshZ_Paroxysmal2(binIndices_patient2{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);


%Patient 3
projLAT_patient3 = projLATs_to_Surface(elecX_Paroxysmal3, elecY_Paroxysmal3, elecZ_Paroxysmal3, meshX_Paroxysmal3, meshY_Paroxysmal3, meshZ_Paroxysmal3, loc_actTimes_3);
numBins = 7;
[binEdges3, binIndices_patient3, remainingIndices3] = partitionToBins_Isochrones(projLAT_patient3, numBins, [] );

figure
plot3(meshX_Paroxysmal3(binIndices_patient3{1}), meshY_Paroxysmal3(binIndices_patient3{1}), meshZ_Paroxysmal3(binIndices_patient3{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec2{1}, 'MarkerEdgeColor', colorVec2{1} );
hold on
plot3(meshX_Paroxysmal3(binIndices_patient3{2}), meshY_Paroxysmal3(binIndices_patient3{2}), meshZ_Paroxysmal3(binIndices_patient3{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec2{2}, 'MarkerEdgeColor', colorVec2{2} );
plot3(meshX_Paroxysmal3(binIndices_patient3{3}), meshY_Paroxysmal3(binIndices_patient3{3}), meshZ_Paroxysmal3(binIndices_patient3{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec2{3}, 'MarkerEdgeColor', colorVec2{3} );
plot3(meshX_Paroxysmal3(binIndices_patient3{4}), meshY_Paroxysmal3(binIndices_patient3{4}), meshZ_Paroxysmal3(binIndices_patient3{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec2{4}, 'MarkerEdgeColor', colorVec2{4} );
plot3(meshX_Paroxysmal3(binIndices_patient3{5}), meshY_Paroxysmal3(binIndices_patient3{5}), meshZ_Paroxysmal3(binIndices_patient3{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec2{5}, 'MarkerEdgeColor', colorVec2{5} );
plot3(meshX_Paroxysmal3(binIndices_patient3{6}), meshY_Paroxysmal3(binIndices_patient3{6}), meshZ_Paroxysmal3(binIndices_patient3{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec2{6}, 'MarkerEdgeColor', colorVec2{6} );
plot3(meshX_Paroxysmal3(binIndices_patient3{7}), meshY_Paroxysmal3(binIndices_patient3{7}), meshZ_Paroxysmal3(binIndices_patient3{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec2{7}, 'MarkerEdgeColor', colorVec2{7} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);




%Patient 6
figure
plot3(meshX_Persistent6, meshY_Persistent6, meshZ_Persistent6, 'ro', 'DisplayName', 'tissue', 'MarkerFaceColor','r', 'MarkerSize',20,'LineWidth',3);
hold on
scatter3(elecX_Persistent6, elecY_Persistent6, elecZ_Persistent6, 230, loc_actTimes_6, 'filled', 'MarkerEdgeColor','k','LineWidth',3);
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (msec)';
% colormap(jet);
caxis([min(loc_actTimes_6)-5, max(loc_actTimes_6)+5]);
% title('H1A Patient');
set(gca, 'FontSize', 35);


projLAT_patient6 = projLATs_to_Surface(elecX_Persistent6, elecY_Persistent6, elecZ_Persistent6, meshX_Persistent6, meshY_Persistent6, meshZ_Persistent6, loc_actTimes_6);
figure
scatter3(meshX_Persistent6, meshY_Persistent6, meshZ_Persistent6, 190, projLAT_patient6, 'filled', 'MarkerEdgeColor','k');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (msec)';
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

numBins = 9;
[binEdges6, binIndices_patient6, remainingIndices6] = partitionToBins_Isochrones(projLAT_patient6, numBins, [] );

figure
plot3(meshX_Persistent6(binIndices_patient6{1}), meshY_Persistent6(binIndices_patient6{1}), meshZ_Persistent6(binIndices_patient6{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Persistent6(binIndices_patient6{2}), meshY_Persistent6(binIndices_patient6{2}), meshZ_Persistent6(binIndices_patient6{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Persistent6(binIndices_patient6{3}), meshY_Persistent6(binIndices_patient6{3}), meshZ_Persistent6(binIndices_patient6{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Persistent6(binIndices_patient6{4}), meshY_Persistent6(binIndices_patient6{4}), meshZ_Persistent6(binIndices_patient6{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Persistent6(binIndices_patient6{5}), meshY_Persistent6(binIndices_patient6{5}), meshZ_Persistent6(binIndices_patient6{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Persistent6(binIndices_patient6{6}), meshY_Persistent6(binIndices_patient6{6}), meshZ_Persistent6(binIndices_patient6{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Persistent6(binIndices_patient6{7}), meshY_Persistent6(binIndices_patient6{7}), meshZ_Persistent6(binIndices_patient6{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Persistent6(binIndices_patient6{8}), meshY_Persistent6(binIndices_patient6{8}), meshZ_Persistent6(binIndices_patient6{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Persistent6(binIndices_patient6{9}), meshY_Persistent6(binIndices_patient6{9}), meshZ_Persistent6(binIndices_patient6{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);




%Patient 7
figure
plot3(meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, 'ro', 'DisplayName', 'tissue', 'MarkerFaceColor','r', 'MarkerSize',20,'LineWidth',3);
hold on
scatter3(elecX_Persistent7, elecY_Persistent7, elecZ_Persistent7, 230, loc_actTimes_7, 'filled', 'MarkerEdgeColor','k','LineWidth',3);
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (msec)';
% colormap(jet);
caxis([min(loc_actTimes_7)-5, max(loc_actTimes_7)+5]);
% title('H1A Patient');
set(gca, 'FontSize', 35);




projLAT_patient7 = projLATs_to_Surface(elecX_Persistent7, elecY_Persistent7, elecZ_Persistent7, meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, loc_actTimes_7);
figure
scatter3(meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, 190, projLAT_patient7, 'filled', 'MarkerEdgeColor','k');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (msec)';
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

numBins = 9;
[binEdges7, binIndices_patient7, remainingIndices7] = partitionToBins_Isochrones(projLAT_patient7, numBins, [] );

figure
plot3(meshX_Persistent7(binIndices_patient7{1}), meshY_Persistent7(binIndices_patient7{1}), meshZ_Persistent7(binIndices_patient7{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Persistent7(binIndices_patient7{2}), meshY_Persistent7(binIndices_patient7{2}), meshZ_Persistent7(binIndices_patient7{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Persistent7(binIndices_patient7{3}), meshY_Persistent7(binIndices_patient7{3}), meshZ_Persistent7(binIndices_patient7{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Persistent7(binIndices_patient7{4}), meshY_Persistent7(binIndices_patient7{4}), meshZ_Persistent7(binIndices_patient7{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Persistent7(binIndices_patient7{5}), meshY_Persistent7(binIndices_patient7{5}), meshZ_Persistent7(binIndices_patient7{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Persistent7(binIndices_patient7{6}), meshY_Persistent7(binIndices_patient7{6}), meshZ_Persistent7(binIndices_patient7{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Persistent7(binIndices_patient7{7}), meshY_Persistent7(binIndices_patient7{7}), meshZ_Persistent7(binIndices_patient7{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Persistent7(binIndices_patient7{8}), meshY_Persistent7(binIndices_patient7{8}), meshZ_Persistent7(binIndices_patient7{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Persistent7(binIndices_patient7{9}), meshY_Persistent7(binIndices_patient7{9}), meshZ_Persistent7(binIndices_patient7{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);



%Patient 8
figure
plot3(meshX_Persistent8, meshY_Persistent8, meshZ_Persistent8, 'ro', 'DisplayName', 'tissue', 'MarkerFaceColor','r', 'MarkerSize',20,'LineWidth',3);
hold on
scatter3(elecX_Persistent8, elecY_Persistent8, elecZ_Persistent8, 230, loc_actTimes_8, 'filled', 'MarkerEdgeColor','k','LineWidth',3);
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (msec)';
% colormap(jet);
caxis([30, max(loc_actTimes_8)+5]);
% title('H1A Patient');
set(gca, 'FontSize', 35);


projLAT_patient8 = projLATs_to_Surface(elecX_Persistent8, elecY_Persistent8, elecZ_Persistent8, meshX_Persistent8, meshY_Persistent8, meshZ_Persistent8, loc_actTimes_8);
figure
scatter3(meshX_Persistent8, meshY_Persistent8, meshZ_Persistent8, 190, projLAT_patient8, 'filled', 'MarkerEdgeColor','k');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (msec)';
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

numBins = 9;
[binEdges8, binIndices_patient8, remainingIndices8] = partitionToBins_Isochrones(projLAT_patient8, numBins, [] );

figure
plot3(meshX_Persistent8(binIndices_patient8{1}), meshY_Persistent8(binIndices_patient8{1}), meshZ_Persistent8(binIndices_patient8{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Persistent8(binIndices_patient8{2}), meshY_Persistent8(binIndices_patient8{2}), meshZ_Persistent8(binIndices_patient8{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Persistent8(binIndices_patient8{3}), meshY_Persistent8(binIndices_patient8{3}), meshZ_Persistent8(binIndices_patient8{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Persistent8(binIndices_patient8{4}), meshY_Persistent8(binIndices_patient8{4}), meshZ_Persistent8(binIndices_patient8{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Persistent8(binIndices_patient8{5}), meshY_Persistent8(binIndices_patient8{5}), meshZ_Persistent8(binIndices_patient8{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Persistent8(binIndices_patient8{6}), meshY_Persistent8(binIndices_patient8{6}), meshZ_Persistent8(binIndices_patient8{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Persistent8(binIndices_patient8{7}), meshY_Persistent8(binIndices_patient8{7}), meshZ_Persistent8(binIndices_patient8{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Persistent8(binIndices_patient8{8}), meshY_Persistent8(binIndices_patient8{8}), meshZ_Persistent8(binIndices_patient8{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Persistent8(binIndices_patient8{9}), meshY_Persistent8(binIndices_patient8{9}), meshZ_Persistent8(binIndices_patient8{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);



%Patient 9
figure
plot3(meshX_Persistent9, meshY_Persistent9, meshZ_Persistent9, 'ro', 'DisplayName', 'tissue', 'MarkerFaceColor','r', 'MarkerSize',20,'LineWidth',3);
hold on
scatter3(elecX_Persistent9, elecY_Persistent9, elecZ_Persistent9, 230, loc_actTimes_9, 'filled', 'MarkerEdgeColor','k','LineWidth',3);
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (msec)';
% colormap(jet);
caxis([min(loc_actTimes_9)-5, max(loc_actTimes_9)+5]);
% title('H1A Patient');
set(gca, 'FontSize', 35);


projLAT_patient9 = projLATs_to_Surface(elecX_Persistent9, elecY_Persistent9, elecZ_Persistent9, meshX_Persistent9, meshY_Persistent9, meshZ_Persistent9, loc_actTimes_9);
figure
scatter3(meshX_Persistent9, meshY_Persistent9, meshZ_Persistent9, 190, projLAT_patient9, 'filled', 'MarkerEdgeColor','k');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (msec)';
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

numBins = 9;
[binEdges9, binIndices_patient9, remainingIndices9] = partitionToBins_Isochrones(projLAT_patient9, numBins, [] );

figure
plot3(meshX_Persistent9(binIndices_patient9{1}), meshY_Persistent9(binIndices_patient9{1}), meshZ_Persistent9(binIndices_patient9{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Persistent9(binIndices_patient9{2}), meshY_Persistent9(binIndices_patient9{2}), meshZ_Persistent9(binIndices_patient9{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Persistent9(binIndices_patient9{3}), meshY_Persistent9(binIndices_patient9{3}), meshZ_Persistent9(binIndices_patient9{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Persistent9(binIndices_patient9{4}), meshY_Persistent9(binIndices_patient9{4}), meshZ_Persistent9(binIndices_patient9{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Persistent9(binIndices_patient9{5}), meshY_Persistent9(binIndices_patient9{5}), meshZ_Persistent9(binIndices_patient9{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Persistent9(binIndices_patient9{6}), meshY_Persistent9(binIndices_patient9{6}), meshZ_Persistent9(binIndices_patient9{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Persistent9(binIndices_patient9{7}), meshY_Persistent9(binIndices_patient9{7}), meshZ_Persistent9(binIndices_patient9{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Persistent9(binIndices_patient9{8}), meshY_Persistent9(binIndices_patient9{8}), meshZ_Persistent9(binIndices_patient9{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Persistent9(binIndices_patient9{9}), meshY_Persistent9(binIndices_patient9{9}), meshZ_Persistent9(binIndices_patient9{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);










%% SIGNALS PLOTS


figure
for i = 1:2:180
indToUse = i;
subplot(2,2,1)
plot(signals_Paroxysmal2(indexCluster1Paroxysmal2_WD(67),:)', 'r-', 'DisplayName','1-Patient', 'LineWidth',3)
hold on
plot(signals_Paroxysmal2_Sim(indexCluster1Paroxysmal2_WD(67),:)', 'r--', 'DisplayName','1-Simulated', 'LineWidth',3)
hold off
% xlabel('Time (ms)')
% ylabel('Vm_n')
% legend()
% set(gca,'FontSize',35)
subplot(2,2,2)
plot(signals_Paroxysmal2(indexCluster2Paroxysmal2_WD(17),:)', 'b-', 'DisplayName','2-Patient', 'LineWidth',3)
hold on
plot(signals_Paroxysmal2_Sim(indexCluster2Paroxysmal2_WD(17),:)', 'b--', 'DisplayName','2-Simulated', 'LineWidth',3)
hold off
% xlabel('Time (ms)')
% ylabel('Vm_n')
% legend()
% set(gca,'FontSize',35)
subplot(2,2,3)
plot(signals_Paroxysmal2(indexCluster3Paroxysmal2_WD(51),:)', 'm-', 'DisplayName','3-Patient', 'LineWidth',3)
hold on
plot(signals_Paroxysmal2_Sim(indexCluster3Paroxysmal2_WD(51),:)', 'm--', 'DisplayName','3-Simulated', 'LineWidth',3)
% title(i)
hold off
% xlabel('Time (ms)')
% ylabel('Vm_n')
% legend()
% set(gca,'FontSize',35)
subplot(2,2,4)
plot(signals_Paroxysmal2(indexCluster3Paroxysmal2_WD(51),:)', 'm-', 'DisplayName','3-Patient', 'LineWidth',3)
hold on
plot(signals_Paroxysmal2_Sim(indexCluster3Paroxysmal2_WD(51),:)', 'm--', 'DisplayName','3-Simulated', 'LineWidth',3)
% title(i)
hold off
% xlabel('Time (ms)')
% ylabel('Vm_n')
% legend()
% set(gca,'FontSize',35)
pause(3)
end







%% SIGNAL-SIGNAL COMPARISON TO CLOSEST PATIENT SIGNAL THAT APPLIES


%Persistent 6
patientNumberToTest = 6;

metricsLumped_6 = [];
metricsIndividual_6 = [];
dummyInd = 1;
dummyInd2 = 1;

for degreeNumberToTest = [10,35,60]
    for fibTypeNumberToTest = 1:4

        loc_each_indOfInterest = find( patient_lumped == patientNumberToTest );

        for patternNumberToTest = 1:10
            
            loc_indicesOfInterest_individual = find( patternNumber == patternNumberToTest & fType == fibTypeNumberToTest & patientNumber == patientNumberToTest & fPercentage == degreeNumberToTest );
            
            if ~isempty(loc_indicesOfInterest_individual)
                ppA_error_loc = abs((100.*((each_PPASimulated(loc_each_indOfInterest) - global_PPASimulated(loc_indicesOfInterest_individual))./(each_PPASimulated(loc_each_indOfInterest)))));
                EGMD_error_loc = abs((100.*((each_EGMDurationSimulated(loc_each_indOfInterest) - global_EGMDurationSimulated(loc_indicesOfInterest_individual))./(each_EGMDurationSimulated(loc_each_indOfInterest)))));
                dC_error_loc = abs((100.*((each_deflectionCountSimulated(loc_each_indOfInterest) - global_deflectionCountSimulated(loc_indicesOfInterest_individual))./(each_deflectionCountSimulated(loc_each_indOfInterest)))));
    
                metricsIndividual_6(dummyInd2, :) = [ fibTypeNumberToTest, degreeNumberToTest, patternNumberToTest, mean(ppA_error_loc), mean(EGMD_error_loc), mean(dC_error_loc),...
                                                      mean(global_loc_LAT_error(loc_indicesOfInterest_individual)), mean(global_loc_similarity_score(loc_indicesOfInterest_individual)), mean(global_loc_xcorr(loc_indicesOfInterest_individual)) ];
                dummyInd2 = dummyInd2 + 1;
            end

        end
    end
end

dummyInd = 1;
for degreeNumberToTest = [10,35,60]
    for fibTypeNumberToTest = 1:4
        lumpedIndices = find( metricsIndividual_6(:,1)==fibTypeNumberToTest & metricsIndividual_6(:,2)==degreeNumberToTest );
        metricsLumped_6(dummyInd,:) = [ fibTypeNumberToTest, degreeNumberToTest, mean(metricsIndividual_6(lumpedIndices, 4:9), 1) ];
        dummyInd = dummyInd + 1;
    end
end






%Persistent 7
patientNumberToTest = 7;

metricsLumped_7 = [];
metricsIndividual_7 = [];
dummyInd = 1;
dummyInd2 = 1;

for degreeNumberToTest = [10,35,60]
    for fibTypeNumberToTest = 1:4

        loc_each_indOfInterest = find( patient_lumped == patientNumberToTest );

        for patternNumberToTest = 1:10
            
            loc_indicesOfInterest_individual = find( patternNumber == patternNumberToTest & fType == fibTypeNumberToTest & patientNumber == patientNumberToTest & fPercentage == degreeNumberToTest );
            
            if ~isempty(loc_indicesOfInterest_individual)
                ppA_error_loc = abs((100.*((each_PPASimulated(loc_each_indOfInterest) - global_PPASimulated(loc_indicesOfInterest_individual))./(each_PPASimulated(loc_each_indOfInterest)))));
                EGMD_error_loc = abs((100.*((each_EGMDurationSimulated(loc_each_indOfInterest) - global_EGMDurationSimulated(loc_indicesOfInterest_individual))./(each_EGMDurationSimulated(loc_each_indOfInterest)))));
                dC_error_loc = abs((100.*((each_deflectionCountSimulated(loc_each_indOfInterest) - global_deflectionCountSimulated(loc_indicesOfInterest_individual))./(each_deflectionCountSimulated(loc_each_indOfInterest)))));
    
                metricsIndividual_7(dummyInd2, :) = [ fibTypeNumberToTest, degreeNumberToTest, patternNumberToTest, mean(ppA_error_loc), mean(EGMD_error_loc), mean(dC_error_loc),...
                                                      mean(global_loc_LAT_error(loc_indicesOfInterest_individual)), mean(global_loc_similarity_score(loc_indicesOfInterest_individual)), mean(global_loc_xcorr(loc_indicesOfInterest_individual)) ];
                dummyInd2 = dummyInd2 + 1;
            end

        end
    end
end

dummyInd = 1;
for degreeNumberToTest = [10,35,60]
    for fibTypeNumberToTest = 1:4
        lumpedIndices = find( metricsIndividual_7(:,1)==fibTypeNumberToTest & metricsIndividual_7(:,2)==degreeNumberToTest );
        metricsLumped_7(dummyInd,:) = [ fibTypeNumberToTest, degreeNumberToTest, mean(metricsIndividual_7(lumpedIndices, 4:9), 1) ];
        dummyInd = dummyInd + 1;
    end
end





%Persistent 8
patientNumberToTest = 8;

metricsLumped_8 = [];
metricsIndividual_8 = [];
dummyInd = 1;
dummyInd2 = 1;

for degreeNumberToTest = [10,35,60]
    for fibTypeNumberToTest = 1:4

        loc_each_indOfInterest = find( patient_lumped == patientNumberToTest );

        for patternNumberToTest = 1:10
            
            loc_indicesOfInterest_individual = find( patternNumber == patternNumberToTest & fType == fibTypeNumberToTest & patientNumber == patientNumberToTest & fPercentage == degreeNumberToTest );
            
            if ~isempty(loc_indicesOfInterest_individual)
                ppA_error_loc = abs((100.*((each_PPASimulated(loc_each_indOfInterest) - global_PPASimulated(loc_indicesOfInterest_individual))./(each_PPASimulated(loc_each_indOfInterest)))));
                EGMD_error_loc = abs((100.*((each_EGMDurationSimulated(loc_each_indOfInterest) - global_EGMDurationSimulated(loc_indicesOfInterest_individual))./(each_EGMDurationSimulated(loc_each_indOfInterest)))));
                dC_error_loc = abs((100.*((each_deflectionCountSimulated(loc_each_indOfInterest) - global_deflectionCountSimulated(loc_indicesOfInterest_individual))./(each_deflectionCountSimulated(loc_each_indOfInterest)))));
    
                metricsIndividual_8(dummyInd2, :) = [ fibTypeNumberToTest, degreeNumberToTest, patternNumberToTest, mean(ppA_error_loc), mean(EGMD_error_loc), mean(dC_error_loc),...
                                                      mean(global_loc_LAT_error(loc_indicesOfInterest_individual)), mean(global_loc_similarity_score(loc_indicesOfInterest_individual)), mean(global_loc_xcorr(loc_indicesOfInterest_individual)) ];
                dummyInd2 = dummyInd2 + 1;
            end

        end
    end
end

dummyInd = 1;
for degreeNumberToTest = [10,35,60]
    for fibTypeNumberToTest = 1:4
        lumpedIndices = find( metricsIndividual_8(:,1)==fibTypeNumberToTest & metricsIndividual_8(:,2)==degreeNumberToTest );
        metricsLumped_8(dummyInd,:) = [ fibTypeNumberToTest, degreeNumberToTest, mean(metricsIndividual_8(lumpedIndices, 4:9), 1) ];
        dummyInd = dummyInd + 1;
    end
end





%Persistent 9
patientNumberToTest = 9;

metricsLumped_9 = [];
metricsIndividual_9 = [];
dummyInd = 1;
dummyInd2 = 1;

for degreeNumberToTest = [10,35,60]
    for fibTypeNumberToTest = 1:4

        loc_each_indOfInterest = find( patient_lumped == patientNumberToTest );

        for patternNumberToTest = 1:10
            
            loc_indicesOfInterest_individual = find( patternNumber == patternNumberToTest & fType == fibTypeNumberToTest & patientNumber == patientNumberToTest & fPercentage == degreeNumberToTest );
            
            if ~isempty(loc_indicesOfInterest_individual)
                ppA_error_loc = abs((100.*((each_PPASimulated(loc_each_indOfInterest) - global_PPASimulated(loc_indicesOfInterest_individual))./(each_PPASimulated(loc_each_indOfInterest)))));
                EGMD_error_loc = abs((100.*((each_EGMDurationSimulated(loc_each_indOfInterest) - global_EGMDurationSimulated(loc_indicesOfInterest_individual))./(each_EGMDurationSimulated(loc_each_indOfInterest)))));
                dC_error_loc = abs((100.*((each_deflectionCountSimulated(loc_each_indOfInterest) - global_deflectionCountSimulated(loc_indicesOfInterest_individual))./(each_deflectionCountSimulated(loc_each_indOfInterest)))));
    
                metricsIndividual_9(dummyInd2, :) = [ fibTypeNumberToTest, degreeNumberToTest, patternNumberToTest, mean(ppA_error_loc), mean(EGMD_error_loc), mean(dC_error_loc),...
                                                      mean(global_loc_LAT_error(loc_indicesOfInterest_individual)), mean(global_loc_similarity_score(loc_indicesOfInterest_individual)), mean(global_loc_xcorr(loc_indicesOfInterest_individual)) ];
                dummyInd2 = dummyInd2 + 1;
            end

        end
    end
end

dummyInd = 1;
for degreeNumberToTest = [10,35,60]
    for fibTypeNumberToTest = 1:4
        lumpedIndices = find( metricsIndividual_9(:,1)==fibTypeNumberToTest & metricsIndividual_9(:,2)==degreeNumberToTest );
        metricsLumped_9(dummyInd,:) = [ fibTypeNumberToTest, degreeNumberToTest, mean(metricsIndividual_9(lumpedIndices, 4:9), 1) ];
        dummyInd = dummyInd + 1;
    end
end






%Checking
indErrorBased6 = find( mean(metricsIndividual_6(:,5:7),2) == min(mean(metricsIndividual_6(:,5:7),2))  );
indCorrBased6 = find( mean(metricsIndividual_6(:,8:9),2) == max(mean(metricsIndividual_6(:,8:9),2))  );
minimum_ErrorBased_individual6 = metricsIndividual_6( indErrorBased6, : );
maximum_CorrBased_individual6 = metricsIndividual_6( indCorrBased6, : );
indErrorBased6_L = find( mean(metricsLumped_6(:,4:6),2) == min(mean(metricsLumped_6(:,4:6),2))  );
indCorrBased6_L = find( mean(metricsLumped_6(:,7:8),2) == max(mean(metricsLumped_6(:,7:8),2))  );
minimum_ErrorBased_Lumped6 = metricsLumped_6( indErrorBased6_L, : );
maximum_CorrBased_Lumped6 = metricsLumped_6( indCorrBased6_L, : );


indErrorBased7 = find( mean(metricsIndividual_7(:,5:7),2) == min(mean(metricsIndividual_7(:,5:7),2))  );
indCorrBased7 = find( mean(metricsIndividual_7(:,8:9),2) == max(mean(metricsIndividual_7(:,8:9),2))  );
minimum_ErrorBased_individual7 = metricsIndividual_7( indErrorBased7, : );
maximum_CorrBased_individual7 = metricsIndividual_7( indCorrBased7, : );
indErrorBased7_L = find( mean(metricsLumped_7(:,4:6),2) == min(mean(metricsLumped_7(:,4:6),2))  );
indCorrBased7_L = find( mean(metricsLumped_7(:,7:8),2) == max(mean(metricsLumped_7(:,7:8),2))  );
minimum_ErrorBased_Lumped7 = metricsLumped_7( indErrorBased7_L, : );
maximum_CorrBased_Lumped7 = metricsLumped_7( indCorrBased7_L, : );


indErrorBased8 = find( mean(metricsIndividual_8(:,5:7),2) == min(mean(metricsIndividual_8(:,5:7),2))  );
indCorrBased8 = find( mean(metricsIndividual_8(:,8:9),2) == max(mean(metricsIndividual_8(:,8:9),2))  );
minimum_ErrorBased_individual8 = metricsIndividual_8( indErrorBased8, : );
maximum_CorrBased_individual8 = metricsIndividual_8( indCorrBased8, : );
indErrorBased8_L = find( mean(metricsLumped_8(:,4:6),2) == min(mean(metricsLumped_8(:,4:6),2))  );
indCorrBased8_L = find( mean(metricsLumped_8(:,7:8),2) == max(mean(metricsLumped_8(:,7:8),2))  );
minimum_ErrorBased_Lumped8 = metricsLumped_8( indErrorBased8_L, : );
maximum_CorrBased_Lumped8 = metricsLumped_8( indCorrBased8_L, : );


indErrorBased9 = find( mean(metricsIndividual_9(:,5:7),2) == min(mean(metricsIndividual_9(:,5:7),2))  );
indCorrBased9 = find( mean(metricsIndividual_9(:,8:9),2) == max(mean(metricsIndividual_9(:,8:9),2))  );
minimum_ErrorBased_individual9 = metricsIndividual_9( indErrorBased9, : );
maximum_CorrBased_individual9 = metricsIndividual_9( indCorrBased9, : );
indErrorBased9_L = find( mean(metricsLumped_9(:,4:6),2) == min(mean(metricsLumped_9(:,4:6),2))  );
indCorrBased9_L = find( mean(metricsLumped_9(:,7:8),2) == max(mean(metricsLumped_9(:,7:8),2))  );
minimum_ErrorBased_Lumped9 = metricsLumped_9( indErrorBased9_L, : );
maximum_CorrBased_Lumped9 = metricsLumped_9( indCorrBased9_L, : );
 



table_smallest_error_lumped = [ minimum_ErrorBased_Lumped6; maximum_CorrBased_Lumped6; minimum_ErrorBased_Lumped7; maximum_CorrBased_Lumped7;
                                minimum_ErrorBased_Lumped8; maximum_CorrBased_Lumped8; minimum_ErrorBased_Lumped9; maximum_CorrBased_Lumped9];

table_smallest_error_individual = [ minimum_ErrorBased_individual6; maximum_CorrBased_individual6; minimum_ErrorBased_individual7; maximum_CorrBased_individual7;
                                minimum_ErrorBased_individual8; maximum_CorrBased_individual8; minimum_ErrorBased_individual9; maximum_CorrBased_individual9];



%Persistent 6
 %Patient 6 Patchy 10% pattern 10
% Tissue with electrodes
patientNumToUse = 6;
patternNumToUse = 10;
degreeNumToUse = 10;
fibTypeNumToUse = 4;
loc_actTimesSim = calculateLATs_from_matrix([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNumToUse), '_', typesFibrosis{fibTypeNumToUse},'FibrosisNT', num2str(degreeNumToUse), '_pattern', num2str(patternNumToUse), '_2000points.csv' ], phie_dt, normalizeBool);

projLAT_patient_loc = projLATs_to_Surface(elecX_Persistent6, elecY_Persistent6, elecZ_Persistent6, meshX_Persistent6, meshY_Persistent6, meshZ_Persistent6, loc_actTimesSim);
numBins = 9;
[binEdges_loc_sm6, binIndices_patient_loc, remainingIndicesLoc6] = partitionToBins_Isochrones(projLAT_patient_loc, numBins, binEdges6 );

figure
plot3(meshX_Persistent6(binIndices_patient_loc{1}), meshY_Persistent6(binIndices_patient_loc{1}), meshZ_Persistent6(binIndices_patient_loc{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Persistent6(binIndices_patient_loc{2}), meshY_Persistent6(binIndices_patient_loc{2}), meshZ_Persistent6(binIndices_patient_loc{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Persistent6(binIndices_patient_loc{3}), meshY_Persistent6(binIndices_patient_loc{3}), meshZ_Persistent6(binIndices_patient_loc{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Persistent6(binIndices_patient_loc{4}), meshY_Persistent6(binIndices_patient_loc{4}), meshZ_Persistent6(binIndices_patient_loc{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Persistent6(binIndices_patient_loc{5}), meshY_Persistent6(binIndices_patient_loc{5}), meshZ_Persistent6(binIndices_patient_loc{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Persistent6(binIndices_patient_loc{6}), meshY_Persistent6(binIndices_patient_loc{6}), meshZ_Persistent6(binIndices_patient_loc{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Persistent6(binIndices_patient_loc{7}), meshY_Persistent6(binIndices_patient_loc{7}), meshZ_Persistent6(binIndices_patient_loc{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Persistent6(binIndices_patient_loc{8}), meshY_Persistent6(binIndices_patient_loc{8}), meshZ_Persistent6(binIndices_patient_loc{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Persistent6(binIndices_patient_loc{9}), meshY_Persistent6(binIndices_patient_loc{9}), meshZ_Persistent6(binIndices_patient_loc{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );

plot3(meshX_Persistent6(remainingIndicesLoc6), meshY_Persistent6(remainingIndicesLoc6), meshZ_Persistent6(remainingIndicesLoc6), 'o', 'MarkerSize', 25, 'MarkerFaceColor', sscanf('8b0000','%2x%2x%2x',[1 3])/255, 'MarkerEdgeColor', sscanf('8b0000','%2x%2x%2x',[1 3])/255 );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);


loc_Signals_6chosen = csvread([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNumToUse), '_', typesFibrosis{fibTypeNumToUse},'FibrosisNT', num2str(degreeNumToUse), '_pattern', num2str(patternNumToUse), '_2000points.csv' ]);
loc_Signals_6chosen = loc_Signals_6chosen(1:1400,:);
loc_Signals_6chosen = UpsampleAndFilterFunction(loc_Signals_6chosen', phie_dt, 5, 200, 1000, 10,normalizeBool)';
% loc_Signals_6chosen = (loc_Signals_6chosen - min(loc_Signals_6chosen))./(max(loc_Signals_6chosen) - min(loc_Signals_6chosen));


%Persistent 7
 %Patient 7 Patchy 10% pattern 10
% Tissue with electrodes
patientNumToUse = 7;
patternNumToUse = 10;
degreeNumToUse = 10;
fibTypeNumToUse = 4;
loc_actTimesSim = calculateLATs_from_matrix([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNumToUse), '_', typesFibrosis{fibTypeNumToUse},'FibrosisNT', num2str(degreeNumToUse), '_pattern', num2str(patternNumToUse), '_2000points.csv' ], phie_dt, normalizeBool);

projLAT_patient_loc = projLATs_to_Surface(elecX_Persistent7, elecY_Persistent7, elecZ_Persistent7, meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, loc_actTimesSim);
numBins = 9;
[binEdges_loc_sm7, binIndices_patient_loc, remainingIndicesLoc7] = partitionToBins_Isochrones(projLAT_patient_loc, numBins, binEdges7 );

figure
plot3(meshX_Persistent7(binIndices_patient_loc{1}), meshY_Persistent7(binIndices_patient_loc{1}), meshZ_Persistent7(binIndices_patient_loc{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Persistent7(binIndices_patient_loc{2}), meshY_Persistent7(binIndices_patient_loc{2}), meshZ_Persistent7(binIndices_patient_loc{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Persistent7(binIndices_patient_loc{3}), meshY_Persistent7(binIndices_patient_loc{3}), meshZ_Persistent7(binIndices_patient_loc{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Persistent7(binIndices_patient_loc{4}), meshY_Persistent7(binIndices_patient_loc{4}), meshZ_Persistent7(binIndices_patient_loc{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Persistent7(binIndices_patient_loc{5}), meshY_Persistent7(binIndices_patient_loc{5}), meshZ_Persistent7(binIndices_patient_loc{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Persistent7(binIndices_patient_loc{6}), meshY_Persistent7(binIndices_patient_loc{6}), meshZ_Persistent7(binIndices_patient_loc{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Persistent7(binIndices_patient_loc{7}), meshY_Persistent7(binIndices_patient_loc{7}), meshZ_Persistent7(binIndices_patient_loc{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Persistent7(binIndices_patient_loc{8}), meshY_Persistent7(binIndices_patient_loc{8}), meshZ_Persistent7(binIndices_patient_loc{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Persistent7(binIndices_patient_loc{9}), meshY_Persistent7(binIndices_patient_loc{9}), meshZ_Persistent7(binIndices_patient_loc{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );

plot3(meshX_Persistent7(remainingIndicesLoc7), meshY_Persistent7(remainingIndicesLoc7), meshZ_Persistent7(remainingIndicesLoc7), 'o', 'MarkerSize', 25, 'MarkerFaceColor', sscanf('8b0000','%2x%2x%2x',[1 3])/255, 'MarkerEdgeColor', sscanf('8b0000','%2x%2x%2x',[1 3])/255 );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);

loc_Signals_7chosen = csvread([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNumToUse), '_', typesFibrosis{fibTypeNumToUse},'FibrosisNT', num2str(degreeNumToUse), '_pattern', num2str(patternNumToUse), '_2000points.csv' ]);
loc_Signals_7chosen = loc_Signals_7chosen(1:1400,:);
loc_Signals_7chosen = UpsampleAndFilterFunction(loc_Signals_7chosen', phie_dt, 5, 200, 1000, 10,normalizeBool)';
% loc_Signals_7chosen = (loc_Signals_7chosen - min(loc_Signals_7chosen))./(max(loc_Signals_7chosen) - min(loc_Signals_7chosen));


%Persistent 8
 %Patient 8 Patchy 10% pattern 4
% Tissue with electrodes
patientNumToUse = 8;
patternNumToUse = 4;
degreeNumToUse = 10;
fibTypeNumToUse = 4;
loc_actTimesSim = calculateLATs_from_matrix([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNumToUse), '_', typesFibrosis{fibTypeNumToUse},'FibrosisNT', num2str(degreeNumToUse), '_pattern', num2str(patternNumToUse), '_2000points.csv' ], phie_dt, normalizeBool);

projLAT_patient_loc = projLATs_to_Surface(elecX_Persistent8, elecY_Persistent8, elecZ_Persistent8, meshX_Persistent8, meshY_Persistent8, meshZ_Persistent8, loc_actTimesSim);
numBins = 9;
[binEdges_loc_sm8, binIndices_patient_loc, remainingIndicesLoc8] = partitionToBins_Isochrones(projLAT_patient_loc, numBins, binEdges8 );

figure
plot3(meshX_Persistent8(binIndices_patient_loc{1}), meshY_Persistent8(binIndices_patient_loc{1}), meshZ_Persistent8(binIndices_patient_loc{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Persistent8(binIndices_patient_loc{2}), meshY_Persistent8(binIndices_patient_loc{2}), meshZ_Persistent8(binIndices_patient_loc{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Persistent8(binIndices_patient_loc{3}), meshY_Persistent8(binIndices_patient_loc{3}), meshZ_Persistent8(binIndices_patient_loc{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Persistent8(binIndices_patient_loc{4}), meshY_Persistent8(binIndices_patient_loc{4}), meshZ_Persistent8(binIndices_patient_loc{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Persistent8(binIndices_patient_loc{5}), meshY_Persistent8(binIndices_patient_loc{5}), meshZ_Persistent8(binIndices_patient_loc{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Persistent8(binIndices_patient_loc{6}), meshY_Persistent8(binIndices_patient_loc{6}), meshZ_Persistent8(binIndices_patient_loc{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Persistent8(binIndices_patient_loc{7}), meshY_Persistent8(binIndices_patient_loc{7}), meshZ_Persistent8(binIndices_patient_loc{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Persistent8(binIndices_patient_loc{8}), meshY_Persistent8(binIndices_patient_loc{8}), meshZ_Persistent8(binIndices_patient_loc{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Persistent8(binIndices_patient_loc{9}), meshY_Persistent8(binIndices_patient_loc{9}), meshZ_Persistent8(binIndices_patient_loc{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );

plot3(meshX_Persistent8(remainingIndicesLoc8), meshY_Persistent8(remainingIndicesLoc8), meshZ_Persistent8(remainingIndicesLoc8), 'o', 'MarkerSize', 25, 'MarkerFaceColor', sscanf('8b0000','%2x%2x%2x',[1 3])/255, 'MarkerEdgeColor', sscanf('8b0000','%2x%2x%2x',[1 3])/255 );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);

loc_Signals_8chosen = csvread([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNumToUse), '_', typesFibrosis{fibTypeNumToUse},'FibrosisNT', num2str(degreeNumToUse), '_pattern', num2str(patternNumToUse), '_2000points.csv' ]);
loc_Signals_8chosen = loc_Signals_8chosen(1:1400,:);
loc_Signals_8chosen = UpsampleAndFilterFunction(loc_Signals_8chosen', phie_dt, 5, 200, 1000, 10,normalizeBool)';
% loc_Signals_8chosen = (loc_Signals_8chosen - min(loc_Signals_8chosen))./(max(loc_Signals_8chosen) - min(loc_Signals_8chosen));


%Persistent 9
 %Patient 9 Compact 35% pattern 5
% Tissue with electrodes
patientNumToUse = 9;
patternNumToUse = 5;
degreeNumToUse = 35;
fibTypeNumToUse = 1;
loc_actTimesSim = calculateLATs_from_matrix([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNumToUse), '_', typesFibrosis{fibTypeNumToUse},'FibrosisNT', num2str(degreeNumToUse), '_pattern', num2str(patternNumToUse), '_2000points.csv' ], phie_dt, normalizeBool);

projLAT_patient_loc = projLATs_to_Surface(elecX_Persistent9, elecY_Persistent9, elecZ_Persistent9, meshX_Persistent9, meshY_Persistent9, meshZ_Persistent9, loc_actTimesSim);
numBins = 9;
[binEdges_loc_sm9, binIndices_patient_loc, remainingIndicesLoc9] = partitionToBins_Isochrones(projLAT_patient_loc, numBins, binEdges9 );

figure
plot3(meshX_Persistent9(binIndices_patient_loc{1}), meshY_Persistent9(binIndices_patient_loc{1}), meshZ_Persistent9(binIndices_patient_loc{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Persistent9(binIndices_patient_loc{2}), meshY_Persistent9(binIndices_patient_loc{2}), meshZ_Persistent9(binIndices_patient_loc{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Persistent9(binIndices_patient_loc{3}), meshY_Persistent9(binIndices_patient_loc{3}), meshZ_Persistent9(binIndices_patient_loc{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Persistent9(binIndices_patient_loc{4}), meshY_Persistent9(binIndices_patient_loc{4}), meshZ_Persistent9(binIndices_patient_loc{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Persistent9(binIndices_patient_loc{5}), meshY_Persistent9(binIndices_patient_loc{5}), meshZ_Persistent9(binIndices_patient_loc{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Persistent9(binIndices_patient_loc{6}), meshY_Persistent9(binIndices_patient_loc{6}), meshZ_Persistent9(binIndices_patient_loc{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Persistent9(binIndices_patient_loc{7}), meshY_Persistent9(binIndices_patient_loc{7}), meshZ_Persistent9(binIndices_patient_loc{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Persistent9(binIndices_patient_loc{8}), meshY_Persistent9(binIndices_patient_loc{8}), meshZ_Persistent9(binIndices_patient_loc{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Persistent9(binIndices_patient_loc{9}), meshY_Persistent9(binIndices_patient_loc{9}), meshZ_Persistent9(binIndices_patient_loc{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );

plot3(meshX_Persistent9(remainingIndicesLoc9), meshY_Persistent9(remainingIndicesLoc9), meshZ_Persistent9(remainingIndicesLoc9), 'o', 'MarkerSize', 25, 'MarkerFaceColor', sscanf('8b0000','%2x%2x%2x',[1 3])/255, 'MarkerEdgeColor', sscanf('8b0000','%2x%2x%2x',[1 3])/255 );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);


loc_Signals_9chosen = csvread([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNumToUse), '_', typesFibrosis{fibTypeNumToUse},'FibrosisNT', num2str(degreeNumToUse), '_pattern', num2str(patternNumToUse), '_2000points.csv' ]);
loc_Signals_9chosen = loc_Signals_9chosen(1:1400,:);
loc_Signals_9chosen = UpsampleAndFilterFunction(loc_Signals_9chosen', phie_dt, 5, 200, 1000, 10,normalizeBool)';
% loc_Signals_9chosen = (loc_Signals_9chosen - min(loc_Signals_9chosen))./(max(loc_Signals_9chosen) - min(loc_Signals_9chosen));




%Showing sample signals for all categories using patient 7
loc_Signals_chosen_show = [];
dummyInd = 1;
for i = 1:4
    for j = [10,35,60]
        loc_Signals_ch = csvread([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(7), '_', typesFibrosis{i},'FibrosisNT', num2str(j), '_pattern', num2str(4), '_2000points.csv' ]);
        loc_Signals_ch = loc_Signals_ch(1:1400,:);
        loc_Signals_ch = UpsampleAndFilterFunction(loc_Signals_ch', phie_dt, 5, 200, 1000, 10,normalizeBool)';
        % loc_Signals_chosen_show(dummyInd, :, :) = (loc_Signals_ch - min(loc_Signals_ch))./(max(loc_Signals_ch) - min(loc_Signals_ch));
        loc_Signals_chosen_show(dummyInd, :, :) = loc_Signals_ch;
        dummyInd = dummyInd + 1;
    end
end

elec1 = [];
for i=1:length(loc_Signals_chosen_show(1,1,:))
    elec_1=[squeeze(expDB.dataSet{7}.map.continuousMapping.fullEgms(i,1,1:140))];
    elec1 = [elec1, ((elec_1 - min(elec_1))./(max(elec_1) - min(elec_1))) ];
end





elec6 = [];
for i=1:length(expDB.dataSet{6}.map.continuousMapping.fullEgms(:,1,1))
    elec_1=[squeeze(expDB.dataSet{6}.map.continuousMapping.fullEgms(i,1,1:140))];
    % elec6 = [elec6, ((elec_1 - min(elec_1))./(max(elec_1) - min(elec_1))) ];
    elec6 = [elec6, elec_1 ];
end

elec7 = [];
for i=1:length(expDB.dataSet{7}.map.continuousMapping.fullEgms(:,1,1))
    elec_1=[squeeze(expDB.dataSet{7}.map.continuousMapping.fullEgms(i,1,1:140))];
    % elec7 = [elec7, ((elec_1 - min(elec_1))./(max(elec_1) - min(elec_1))) ];
    elec7 = [elec7, elec_1 ];
end

elec8 = [];
for i=1:length(expDB.dataSet{8}.map.continuousMapping.fullEgms(:,1,1))
    elec_1=[squeeze(expDB.dataSet{8}.map.continuousMapping.fullEgms(i,1,1:140))];
    % elec8 = [elec8, ((elec_1 - min(elec_1))./(max(elec_1) - min(elec_1))) ];
    elec8 = [elec8, elec_1 ];
end

elec9 = [];
for i=1:length(expDB.dataSet{9}.map.continuousMapping.fullEgms(:,1,1))
    elec_1=[squeeze(expDB.dataSet{9}.map.continuousMapping.fullEgms(i,1,1:140))];
    % elec9 = [elec9, ((elec_1 - min(elec_1))./(max(elec_1) - min(elec_1))) ];
    elec9 = [elec9, elec_1 ];
end






indToUse = 123;
figure
% for i = 1:2:400
% indToUse = i;
subplot(2,2,1)
plot(loc_Signals_chosen_show(1,:,indToUse+100)', 'r-', 'DisplayName','Simulated', 'LineWidth',3)
hold on
plot(loc_Signals_chosen_show(1,:,indToUse+350)', 'k-', 'DisplayName','Simulated', 'LineWidth',3)
plot(loc_Signals_chosen_show(1,:,indToUse)', 'b-', 'DisplayName','Simulated', 'LineWidth',3)
hold off
% title(i)
xlabel('Time(msec)')
ylabel('Ve (mV)')
% legend()
set(gca,'FontSize',30)
subplot(2,2,2)
plot(loc_Signals_chosen_show(4,:,indToUse+100)', 'r-', 'DisplayName','Simulated', 'LineWidth',3)
hold on
plot(loc_Signals_chosen_show(4,:,indToUse+350)', 'k-', 'DisplayName','Simulated', 'LineWidth',3)
plot(loc_Signals_chosen_show(4,:,indToUse)', 'b-', 'DisplayName','Simulated', 'LineWidth',3)
hold off
% title(i)
xlabel('Time(msec)')
ylabel('Ve (mV)')
% legend()
set(gca,'FontSize',30)
subplot(2,2,3)
% plot(elec1(:, indToUse), 'r-', 'DisplayName','Patient', 'LineWidth',4)
plot(loc_Signals_chosen_show(7,:,indToUse+100)', 'r-', 'DisplayName','Simulated', 'LineWidth',3)
hold on
plot(loc_Signals_chosen_show(7,:,indToUse+350)', 'k-', 'DisplayName','Simulated', 'LineWidth',3)
plot(loc_Signals_chosen_show(7,:,indToUse)', 'b-', 'DisplayName','Simulated', 'LineWidth',3)
hold off
% title(i)
xlabel('Time(msec)')
ylabel('Ve (mV)')
% legend()
set(gca,'FontSize',30)
subplot(2,2,4)
% plot(elec1(:, indToUse), 'r-', 'DisplayName','Patient', 'LineWidth',4)
plot(loc_Signals_chosen_show(10,:,indToUse+100)', 'r-', 'DisplayName','Simulated', 'LineWidth',3)
hold on
plot(loc_Signals_chosen_show(10,:,indToUse+350)', 'k-', 'DisplayName','Simulated', 'LineWidth',3)
plot(loc_Signals_chosen_show(10,:,indToUse)', 'b-', 'DisplayName','Simulated', 'LineWidth',3)
hold off
% title(i)
xlabel('Time(msec)')
ylabel('Ve (mV)')
% legend()
set(gca,'FontSize',30)
% pause(3)
% end









figure
% pause(5)
% for i = 1:2:400
subplot(2,2,1)
indToUse = 133;
yyaxis right
plot(loc_Signals_6chosen(:,indToUse)', 'r-', 'DisplayName','Simulated', 'LineWidth',3)
ylabel('Ve (mV)')
ylim([min(loc_Signals_6chosen(:,indToUse)') max(loc_Signals_6chosen(:,indToUse)')])
hold on
yyaxis left
plot(elec6(:,indToUse)', 'b-', 'DisplayName','Patient', 'LineWidth',3)
hold off
% title(i)
ylim([min(elec6(:,indToUse)') max(elec6(:,indToUse)')])
xlabel('Time(msec)')
ylabel('Ve (mV)')
legend()
set(gca,'FontSize',30)

subplot(2,2,2)
indToUse = 149;
yyaxis right
plot(loc_Signals_7chosen(:,indToUse)', 'r-', 'DisplayName','Simulated', 'LineWidth',3)
ylabel('Ve (mV)')
ylim([min(loc_Signals_7chosen(:,indToUse)') max(loc_Signals_7chosen(:,indToUse)')])
hold on
yyaxis left
plot(elec7(:,indToUse)', 'b-', 'DisplayName','Patient', 'LineWidth',3)
hold off
% title(i)
ylim([min(elec7(:,indToUse)') max(elec7(:,indToUse)')])
xlabel('Time(msec)')
ylabel('Ve (mV)')
legend()
set(gca,'FontSize',30)

subplot(2,2,3)
indToUse = 57;
yyaxis right
plot(loc_Signals_8chosen(:,indToUse)', 'r-', 'DisplayName','Simulated', 'LineWidth',3)
ylabel('Ve (mV)')
ylim([min(loc_Signals_8chosen(:,indToUse)') max(loc_Signals_8chosen(:,indToUse)')])
hold on
yyaxis left
plot(elec8(:,indToUse)', 'b-', 'DisplayName','Patient', 'LineWidth',3)
hold off
% title(i)
ylim([min(elec8(:,indToUse)') max(elec8(:,indToUse)')])
xlabel('Time(msec)')
ylabel('Ve (mV)')
legend()
set(gca,'FontSize',30)

subplot(2,2,4)
indToUse = 53;
yyaxis right
plot(loc_Signals_9chosen(:,indToUse)', 'r-', 'DisplayName','Simulated', 'LineWidth',3)
ylabel('Ve (mV)')
ylim([min(loc_Signals_9chosen(:,indToUse)') max(loc_Signals_9chosen(:,indToUse)')])
hold on
yyaxis left
plot(elec9(:,indToUse)', 'b-', 'DisplayName','Patient', 'LineWidth',3)
hold off
% title(i)
ylim([min(elec9(:,indToUse)') max(elec9(:,indToUse)')])
xlabel('Time(msec)')
ylabel('Ve (mV)')
legend()
set(gca,'FontSize',30)

% pause(5)
% end



%%



%COMPACT
 %Patient 7 Compact 10% pattern 6
% Tissue with electrodes
patientNumToUse = 7;
patternNumToUse = 6;
degreeNumToUse = 10;
fibTypeNumToUse = 1;
loc_actTimesSim = calculateLATs_from_matrix([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNumToUse), '_', typesFibrosis{fibTypeNumToUse},'FibrosisNT', num2str(degreeNumToUse), '_pattern', num2str(patternNumToUse), '_2000points.csv' ], phie_dt, normalizeBool);

projLAT_patient_loc = projLATs_to_Surface(elecX_Persistent7, elecY_Persistent7, elecZ_Persistent7, meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, loc_actTimesSim);
figure
scatter3(meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, 190, projLAT_patient_loc, 'filled', 'MarkerEdgeColor','k');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (msec)';
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

numBins = 9;
[binEdges_loc, binIndices_patient_loc] = partitionToBins_Isochrones(projLAT_patient_loc, numBins, binEdges7 );

figure
plot3(meshX_Persistent7(binIndices_patient_loc{1}), meshY_Persistent7(binIndices_patient_loc{1}), meshZ_Persistent7(binIndices_patient_loc{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Persistent7(binIndices_patient_loc{2}), meshY_Persistent7(binIndices_patient_loc{2}), meshZ_Persistent7(binIndices_patient_loc{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Persistent7(binIndices_patient_loc{3}), meshY_Persistent7(binIndices_patient_loc{3}), meshZ_Persistent7(binIndices_patient_loc{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Persistent7(binIndices_patient_loc{4}), meshY_Persistent7(binIndices_patient_loc{4}), meshZ_Persistent7(binIndices_patient_loc{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Persistent7(binIndices_patient_loc{5}), meshY_Persistent7(binIndices_patient_loc{5}), meshZ_Persistent7(binIndices_patient_loc{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Persistent7(binIndices_patient_loc{6}), meshY_Persistent7(binIndices_patient_loc{6}), meshZ_Persistent7(binIndices_patient_loc{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Persistent7(binIndices_patient_loc{7}), meshY_Persistent7(binIndices_patient_loc{7}), meshZ_Persistent7(binIndices_patient_loc{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Persistent7(binIndices_patient_loc{8}), meshY_Persistent7(binIndices_patient_loc{8}), meshZ_Persistent7(binIndices_patient_loc{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Persistent7(binIndices_patient_loc{9}), meshY_Persistent7(binIndices_patient_loc{9}), meshZ_Persistent7(binIndices_patient_loc{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);




 %Patient 7 Compact 35% pattern 6
% Tissue with electrodes
patientNumToUse = 7;
patternNumToUse = 6;
degreeNumToUse = 35;
fibTypeNumToUse = 1;
loc_actTimesSim = calculateLATs_from_matrix([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNumToUse), '_', typesFibrosis{fibTypeNumToUse},'FibrosisNT', num2str(degreeNumToUse), '_pattern', num2str(patternNumToUse), '_2000points.csv' ], phie_dt, normalizeBool);

projLAT_patient_loc = projLATs_to_Surface(elecX_Persistent7, elecY_Persistent7, elecZ_Persistent7, meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, loc_actTimesSim);
figure
scatter3(meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, 190, projLAT_patient_loc, 'filled', 'MarkerEdgeColor','k');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (msec)';
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

numBins = 9;
[binEdges_loc, binIndices_patient_loc] = partitionToBins_Isochrones(projLAT_patient_loc, numBins, binEdges7 );

figure
plot3(meshX_Persistent7(binIndices_patient_loc{1}), meshY_Persistent7(binIndices_patient_loc{1}), meshZ_Persistent7(binIndices_patient_loc{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Persistent7(binIndices_patient_loc{2}), meshY_Persistent7(binIndices_patient_loc{2}), meshZ_Persistent7(binIndices_patient_loc{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Persistent7(binIndices_patient_loc{3}), meshY_Persistent7(binIndices_patient_loc{3}), meshZ_Persistent7(binIndices_patient_loc{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Persistent7(binIndices_patient_loc{4}), meshY_Persistent7(binIndices_patient_loc{4}), meshZ_Persistent7(binIndices_patient_loc{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Persistent7(binIndices_patient_loc{5}), meshY_Persistent7(binIndices_patient_loc{5}), meshZ_Persistent7(binIndices_patient_loc{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Persistent7(binIndices_patient_loc{6}), meshY_Persistent7(binIndices_patient_loc{6}), meshZ_Persistent7(binIndices_patient_loc{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Persistent7(binIndices_patient_loc{7}), meshY_Persistent7(binIndices_patient_loc{7}), meshZ_Persistent7(binIndices_patient_loc{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Persistent7(binIndices_patient_loc{8}), meshY_Persistent7(binIndices_patient_loc{8}), meshZ_Persistent7(binIndices_patient_loc{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Persistent7(binIndices_patient_loc{9}), meshY_Persistent7(binIndices_patient_loc{9}), meshZ_Persistent7(binIndices_patient_loc{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);




 %Patient 7 Compact 60% pattern 6
% Tissue with electrodes
patientNumToUse = 7;
patternNumToUse = 6;
degreeNumToUse = 60;
fibTypeNumToUse = 1;
loc_actTimesSim = calculateLATs_from_matrix([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNumToUse), '_', typesFibrosis{fibTypeNumToUse},'FibrosisNT', num2str(degreeNumToUse), '_pattern', num2str(patternNumToUse), '_2000points.csv' ], phie_dt, normalizeBool);

projLAT_patient_loc = projLATs_to_Surface(elecX_Persistent7, elecY_Persistent7, elecZ_Persistent7, meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, loc_actTimesSim);
figure
scatter3(meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, 190, projLAT_patient_loc, 'filled', 'MarkerEdgeColor','k');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (msec)';
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

numBins = 9;
[binEdges_loc, binIndices_patient_loc] = partitionToBins_Isochrones(projLAT_patient_loc, numBins, binEdges7 );

figure
plot3(meshX_Persistent7(binIndices_patient_loc{1}), meshY_Persistent7(binIndices_patient_loc{1}), meshZ_Persistent7(binIndices_patient_loc{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Persistent7(binIndices_patient_loc{2}), meshY_Persistent7(binIndices_patient_loc{2}), meshZ_Persistent7(binIndices_patient_loc{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Persistent7(binIndices_patient_loc{3}), meshY_Persistent7(binIndices_patient_loc{3}), meshZ_Persistent7(binIndices_patient_loc{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Persistent7(binIndices_patient_loc{4}), meshY_Persistent7(binIndices_patient_loc{4}), meshZ_Persistent7(binIndices_patient_loc{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Persistent7(binIndices_patient_loc{5}), meshY_Persistent7(binIndices_patient_loc{5}), meshZ_Persistent7(binIndices_patient_loc{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Persistent7(binIndices_patient_loc{6}), meshY_Persistent7(binIndices_patient_loc{6}), meshZ_Persistent7(binIndices_patient_loc{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Persistent7(binIndices_patient_loc{7}), meshY_Persistent7(binIndices_patient_loc{7}), meshZ_Persistent7(binIndices_patient_loc{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Persistent7(binIndices_patient_loc{8}), meshY_Persistent7(binIndices_patient_loc{8}), meshZ_Persistent7(binIndices_patient_loc{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Persistent7(binIndices_patient_loc{9}), meshY_Persistent7(binIndices_patient_loc{9}), meshZ_Persistent7(binIndices_patient_loc{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);




 %Patient 7 Compact 60% pattern 3
% Tissue with electrodes
patientNumToUse = 7;
patternNumToUse = 3;
degreeNumToUse = 60;
fibTypeNumToUse = 1;
loc_actTimesSim = calculateLATs_from_matrix([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNumToUse), '_', typesFibrosis{fibTypeNumToUse},'FibrosisNT', num2str(degreeNumToUse), '_pattern', num2str(patternNumToUse), '_2000points.csv' ], phie_dt, normalizeBool);

projLAT_patient_loc = projLATs_to_Surface(elecX_Persistent7, elecY_Persistent7, elecZ_Persistent7, meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, loc_actTimesSim);
figure
scatter3(meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, 190, projLAT_patient_loc, 'filled', 'MarkerEdgeColor','k');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (msec)';
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

numBins = 9;
[binEdges_loc, binIndices_patient_loc] = partitionToBins_Isochrones(projLAT_patient_loc, numBins, binEdges7 );

figure
plot3(meshX_Persistent7(binIndices_patient_loc{1}), meshY_Persistent7(binIndices_patient_loc{1}), meshZ_Persistent7(binIndices_patient_loc{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Persistent7(binIndices_patient_loc{2}), meshY_Persistent7(binIndices_patient_loc{2}), meshZ_Persistent7(binIndices_patient_loc{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Persistent7(binIndices_patient_loc{3}), meshY_Persistent7(binIndices_patient_loc{3}), meshZ_Persistent7(binIndices_patient_loc{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Persistent7(binIndices_patient_loc{4}), meshY_Persistent7(binIndices_patient_loc{4}), meshZ_Persistent7(binIndices_patient_loc{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Persistent7(binIndices_patient_loc{5}), meshY_Persistent7(binIndices_patient_loc{5}), meshZ_Persistent7(binIndices_patient_loc{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Persistent7(binIndices_patient_loc{6}), meshY_Persistent7(binIndices_patient_loc{6}), meshZ_Persistent7(binIndices_patient_loc{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Persistent7(binIndices_patient_loc{7}), meshY_Persistent7(binIndices_patient_loc{7}), meshZ_Persistent7(binIndices_patient_loc{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Persistent7(binIndices_patient_loc{8}), meshY_Persistent7(binIndices_patient_loc{8}), meshZ_Persistent7(binIndices_patient_loc{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Persistent7(binIndices_patient_loc{9}), meshY_Persistent7(binIndices_patient_loc{9}), meshZ_Persistent7(binIndices_patient_loc{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);




 %Patient 7 Compact 60% pattern 1
% Tissue with electrodes
patientNumToUse = 7;
patternNumToUse = 1;
degreeNumToUse = 60;
fibTypeNumToUse = 1;
loc_actTimesSim = calculateLATs_from_matrix([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNumToUse), '_', typesFibrosis{fibTypeNumToUse},'FibrosisNT', num2str(degreeNumToUse), '_pattern', num2str(patternNumToUse), '_2000points.csv' ], phie_dt, normalizeBool);

projLAT_patient_loc = projLATs_to_Surface(elecX_Persistent7, elecY_Persistent7, elecZ_Persistent7, meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, loc_actTimesSim);
figure
scatter3(meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, 190, projLAT_patient_loc, 'filled', 'MarkerEdgeColor','k');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (msec)';
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

numBins = 9;
[binEdges_loc, binIndices_patient_loc] = partitionToBins_Isochrones(projLAT_patient_loc, numBins, binEdges7 );

figure
plot3(meshX_Persistent7(binIndices_patient_loc{1}), meshY_Persistent7(binIndices_patient_loc{1}), meshZ_Persistent7(binIndices_patient_loc{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Persistent7(binIndices_patient_loc{2}), meshY_Persistent7(binIndices_patient_loc{2}), meshZ_Persistent7(binIndices_patient_loc{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Persistent7(binIndices_patient_loc{3}), meshY_Persistent7(binIndices_patient_loc{3}), meshZ_Persistent7(binIndices_patient_loc{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Persistent7(binIndices_patient_loc{4}), meshY_Persistent7(binIndices_patient_loc{4}), meshZ_Persistent7(binIndices_patient_loc{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Persistent7(binIndices_patient_loc{5}), meshY_Persistent7(binIndices_patient_loc{5}), meshZ_Persistent7(binIndices_patient_loc{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Persistent7(binIndices_patient_loc{6}), meshY_Persistent7(binIndices_patient_loc{6}), meshZ_Persistent7(binIndices_patient_loc{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Persistent7(binIndices_patient_loc{7}), meshY_Persistent7(binIndices_patient_loc{7}), meshZ_Persistent7(binIndices_patient_loc{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Persistent7(binIndices_patient_loc{8}), meshY_Persistent7(binIndices_patient_loc{8}), meshZ_Persistent7(binIndices_patient_loc{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Persistent7(binIndices_patient_loc{9}), meshY_Persistent7(binIndices_patient_loc{9}), meshZ_Persistent7(binIndices_patient_loc{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);






%DIFFUSE
 %Patient 7 Diffuse 10% pattern 6
% Tissue with electrodes
patientNumToUse = 7;
patternNumToUse = 6;
degreeNumToUse = 10;
fibTypeNumToUse = 2;
loc_actTimesSim = calculateLATs_from_matrix([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNumToUse), '_', typesFibrosis{fibTypeNumToUse},'FibrosisNT', num2str(degreeNumToUse), '_pattern', num2str(patternNumToUse), '_2000points.csv' ], phie_dt, normalizeBool);

projLAT_patient_loc = projLATs_to_Surface(elecX_Persistent7, elecY_Persistent7, elecZ_Persistent7, meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, loc_actTimesSim);
figure
scatter3(meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, 190, projLAT_patient_loc, 'filled', 'MarkerEdgeColor','k');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (msec)';
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

numBins = 9;
[binEdges_loc, binIndices_patient_loc] = partitionToBins_Isochrones(projLAT_patient_loc, numBins, binEdges7 );

figure
plot3(meshX_Persistent7(binIndices_patient_loc{1}), meshY_Persistent7(binIndices_patient_loc{1}), meshZ_Persistent7(binIndices_patient_loc{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Persistent7(binIndices_patient_loc{2}), meshY_Persistent7(binIndices_patient_loc{2}), meshZ_Persistent7(binIndices_patient_loc{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Persistent7(binIndices_patient_loc{3}), meshY_Persistent7(binIndices_patient_loc{3}), meshZ_Persistent7(binIndices_patient_loc{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Persistent7(binIndices_patient_loc{4}), meshY_Persistent7(binIndices_patient_loc{4}), meshZ_Persistent7(binIndices_patient_loc{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Persistent7(binIndices_patient_loc{5}), meshY_Persistent7(binIndices_patient_loc{5}), meshZ_Persistent7(binIndices_patient_loc{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Persistent7(binIndices_patient_loc{6}), meshY_Persistent7(binIndices_patient_loc{6}), meshZ_Persistent7(binIndices_patient_loc{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Persistent7(binIndices_patient_loc{7}), meshY_Persistent7(binIndices_patient_loc{7}), meshZ_Persistent7(binIndices_patient_loc{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Persistent7(binIndices_patient_loc{8}), meshY_Persistent7(binIndices_patient_loc{8}), meshZ_Persistent7(binIndices_patient_loc{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Persistent7(binIndices_patient_loc{9}), meshY_Persistent7(binIndices_patient_loc{9}), meshZ_Persistent7(binIndices_patient_loc{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);




 %Patient 7 Diffuse 35% pattern 6
% Tissue with electrodes
patientNumToUse = 7;
patternNumToUse = 6;
degreeNumToUse = 35;
fibTypeNumToUse = 2;
loc_actTimesSim = calculateLATs_from_matrix([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNumToUse), '_', typesFibrosis{fibTypeNumToUse},'FibrosisNT', num2str(degreeNumToUse), '_pattern', num2str(patternNumToUse), '_2000points.csv' ], phie_dt, normalizeBool);

projLAT_patient_loc = projLATs_to_Surface(elecX_Persistent7, elecY_Persistent7, elecZ_Persistent7, meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, loc_actTimesSim);
figure
scatter3(meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, 190, projLAT_patient_loc, 'filled', 'MarkerEdgeColor','k');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (msec)';
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

numBins = 9;
[binEdges_loc, binIndices_patient_loc] = partitionToBins_Isochrones(projLAT_patient_loc, numBins, binEdges7 );

figure
plot3(meshX_Persistent7(binIndices_patient_loc{1}), meshY_Persistent7(binIndices_patient_loc{1}), meshZ_Persistent7(binIndices_patient_loc{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Persistent7(binIndices_patient_loc{2}), meshY_Persistent7(binIndices_patient_loc{2}), meshZ_Persistent7(binIndices_patient_loc{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Persistent7(binIndices_patient_loc{3}), meshY_Persistent7(binIndices_patient_loc{3}), meshZ_Persistent7(binIndices_patient_loc{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Persistent7(binIndices_patient_loc{4}), meshY_Persistent7(binIndices_patient_loc{4}), meshZ_Persistent7(binIndices_patient_loc{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Persistent7(binIndices_patient_loc{5}), meshY_Persistent7(binIndices_patient_loc{5}), meshZ_Persistent7(binIndices_patient_loc{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Persistent7(binIndices_patient_loc{6}), meshY_Persistent7(binIndices_patient_loc{6}), meshZ_Persistent7(binIndices_patient_loc{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Persistent7(binIndices_patient_loc{7}), meshY_Persistent7(binIndices_patient_loc{7}), meshZ_Persistent7(binIndices_patient_loc{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Persistent7(binIndices_patient_loc{8}), meshY_Persistent7(binIndices_patient_loc{8}), meshZ_Persistent7(binIndices_patient_loc{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Persistent7(binIndices_patient_loc{9}), meshY_Persistent7(binIndices_patient_loc{9}), meshZ_Persistent7(binIndices_patient_loc{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);




 %Patient 7 Diffuse 60% pattern 6
% Tissue with electrodes
patientNumToUse = 7;
patternNumToUse = 6;
degreeNumToUse = 60;
fibTypeNumToUse = 2;
loc_actTimesSim = calculateLATs_from_matrix([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNumToUse), '_', typesFibrosis{fibTypeNumToUse},'FibrosisNT', num2str(degreeNumToUse), '_pattern', num2str(patternNumToUse), '_2000points.csv' ], phie_dt, normalizeBool);

projLAT_patient_loc = projLATs_to_Surface(elecX_Persistent7, elecY_Persistent7, elecZ_Persistent7, meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, loc_actTimesSim);
figure
scatter3(meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, 190, projLAT_patient_loc, 'filled', 'MarkerEdgeColor','k');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (msec)';
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

numBins = 9;
[binEdges_loc, binIndices_patient_loc] = partitionToBins_Isochrones(projLAT_patient_loc, numBins, binEdges7 );

figure
plot3(meshX_Persistent7(binIndices_patient_loc{1}), meshY_Persistent7(binIndices_patient_loc{1}), meshZ_Persistent7(binIndices_patient_loc{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Persistent7(binIndices_patient_loc{2}), meshY_Persistent7(binIndices_patient_loc{2}), meshZ_Persistent7(binIndices_patient_loc{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Persistent7(binIndices_patient_loc{3}), meshY_Persistent7(binIndices_patient_loc{3}), meshZ_Persistent7(binIndices_patient_loc{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Persistent7(binIndices_patient_loc{4}), meshY_Persistent7(binIndices_patient_loc{4}), meshZ_Persistent7(binIndices_patient_loc{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Persistent7(binIndices_patient_loc{5}), meshY_Persistent7(binIndices_patient_loc{5}), meshZ_Persistent7(binIndices_patient_loc{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Persistent7(binIndices_patient_loc{6}), meshY_Persistent7(binIndices_patient_loc{6}), meshZ_Persistent7(binIndices_patient_loc{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Persistent7(binIndices_patient_loc{7}), meshY_Persistent7(binIndices_patient_loc{7}), meshZ_Persistent7(binIndices_patient_loc{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Persistent7(binIndices_patient_loc{8}), meshY_Persistent7(binIndices_patient_loc{8}), meshZ_Persistent7(binIndices_patient_loc{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Persistent7(binIndices_patient_loc{9}), meshY_Persistent7(binIndices_patient_loc{9}), meshZ_Persistent7(binIndices_patient_loc{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);









%INTERSTITIAL
 %Patient 7 Interstitial 10% pattern 6
% Tissue with electrodes
patientNumToUse = 7;
patternNumToUse = 6;
degreeNumToUse = 10;
fibTypeNumToUse = 3;
loc_actTimesSim = calculateLATs_from_matrix([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNumToUse), '_', typesFibrosis{fibTypeNumToUse},'FibrosisNT', num2str(degreeNumToUse), '_pattern', num2str(patternNumToUse), '_2000points.csv' ], phie_dt, normalizeBool);

projLAT_patient_loc = projLATs_to_Surface(elecX_Persistent7, elecY_Persistent7, elecZ_Persistent7, meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, loc_actTimesSim);
figure
scatter3(meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, 190, projLAT_patient_loc, 'filled', 'MarkerEdgeColor','k');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (msec)';
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

numBins = 9;
[binEdges_loc, binIndices_patient_loc] = partitionToBins_Isochrones(projLAT_patient_loc, numBins, binEdges7 );

figure
plot3(meshX_Persistent7(binIndices_patient_loc{1}), meshY_Persistent7(binIndices_patient_loc{1}), meshZ_Persistent7(binIndices_patient_loc{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Persistent7(binIndices_patient_loc{2}), meshY_Persistent7(binIndices_patient_loc{2}), meshZ_Persistent7(binIndices_patient_loc{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Persistent7(binIndices_patient_loc{3}), meshY_Persistent7(binIndices_patient_loc{3}), meshZ_Persistent7(binIndices_patient_loc{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Persistent7(binIndices_patient_loc{4}), meshY_Persistent7(binIndices_patient_loc{4}), meshZ_Persistent7(binIndices_patient_loc{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Persistent7(binIndices_patient_loc{5}), meshY_Persistent7(binIndices_patient_loc{5}), meshZ_Persistent7(binIndices_patient_loc{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Persistent7(binIndices_patient_loc{6}), meshY_Persistent7(binIndices_patient_loc{6}), meshZ_Persistent7(binIndices_patient_loc{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Persistent7(binIndices_patient_loc{7}), meshY_Persistent7(binIndices_patient_loc{7}), meshZ_Persistent7(binIndices_patient_loc{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Persistent7(binIndices_patient_loc{8}), meshY_Persistent7(binIndices_patient_loc{8}), meshZ_Persistent7(binIndices_patient_loc{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Persistent7(binIndices_patient_loc{9}), meshY_Persistent7(binIndices_patient_loc{9}), meshZ_Persistent7(binIndices_patient_loc{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);




 %Patient 7 Interstitial 35% pattern 6
% Tissue with electrodes
patientNumToUse = 7;
patternNumToUse = 6;
degreeNumToUse = 35;
fibTypeNumToUse = 3;
loc_actTimesSim = calculateLATs_from_matrix([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNumToUse), '_', typesFibrosis{fibTypeNumToUse},'FibrosisNT', num2str(degreeNumToUse), '_pattern', num2str(patternNumToUse), '_2000points.csv' ], phie_dt, normalizeBool);

projLAT_patient_loc = projLATs_to_Surface(elecX_Persistent7, elecY_Persistent7, elecZ_Persistent7, meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, loc_actTimesSim);
figure
scatter3(meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, 190, projLAT_patient_loc, 'filled', 'MarkerEdgeColor','k');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (msec)';
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

numBins = 9;
[binEdges_loc, binIndices_patient_loc] = partitionToBins_Isochrones(projLAT_patient_loc, numBins, binEdges7 );

figure
plot3(meshX_Persistent7(binIndices_patient_loc{1}), meshY_Persistent7(binIndices_patient_loc{1}), meshZ_Persistent7(binIndices_patient_loc{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Persistent7(binIndices_patient_loc{2}), meshY_Persistent7(binIndices_patient_loc{2}), meshZ_Persistent7(binIndices_patient_loc{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Persistent7(binIndices_patient_loc{3}), meshY_Persistent7(binIndices_patient_loc{3}), meshZ_Persistent7(binIndices_patient_loc{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Persistent7(binIndices_patient_loc{4}), meshY_Persistent7(binIndices_patient_loc{4}), meshZ_Persistent7(binIndices_patient_loc{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Persistent7(binIndices_patient_loc{5}), meshY_Persistent7(binIndices_patient_loc{5}), meshZ_Persistent7(binIndices_patient_loc{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Persistent7(binIndices_patient_loc{6}), meshY_Persistent7(binIndices_patient_loc{6}), meshZ_Persistent7(binIndices_patient_loc{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Persistent7(binIndices_patient_loc{7}), meshY_Persistent7(binIndices_patient_loc{7}), meshZ_Persistent7(binIndices_patient_loc{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Persistent7(binIndices_patient_loc{8}), meshY_Persistent7(binIndices_patient_loc{8}), meshZ_Persistent7(binIndices_patient_loc{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Persistent7(binIndices_patient_loc{9}), meshY_Persistent7(binIndices_patient_loc{9}), meshZ_Persistent7(binIndices_patient_loc{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);




 %Patient 7 Interstitial 60% pattern 6
% Tissue with electrodes
patientNumToUse = 7;
patternNumToUse = 6;
degreeNumToUse = 60;
fibTypeNumToUse = 3;
loc_actTimesSim = calculateLATs_from_matrix([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNumToUse), '_', typesFibrosis{fibTypeNumToUse},'FibrosisNT', num2str(degreeNumToUse), '_pattern', num2str(patternNumToUse), '_2000points.csv' ], phie_dt, normalizeBool);

projLAT_patient_loc = projLATs_to_Surface(elecX_Persistent7, elecY_Persistent7, elecZ_Persistent7, meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, loc_actTimesSim);
figure
scatter3(meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, 190, projLAT_patient_loc, 'filled', 'MarkerEdgeColor','k');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (msec)';
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

numBins = 9;
[binEdges_loc, binIndices_patient_loc] = partitionToBins_Isochrones(projLAT_patient_loc, numBins, binEdges7 );

figure
plot3(meshX_Persistent7(binIndices_patient_loc{1}), meshY_Persistent7(binIndices_patient_loc{1}), meshZ_Persistent7(binIndices_patient_loc{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Persistent7(binIndices_patient_loc{2}), meshY_Persistent7(binIndices_patient_loc{2}), meshZ_Persistent7(binIndices_patient_loc{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Persistent7(binIndices_patient_loc{3}), meshY_Persistent7(binIndices_patient_loc{3}), meshZ_Persistent7(binIndices_patient_loc{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Persistent7(binIndices_patient_loc{4}), meshY_Persistent7(binIndices_patient_loc{4}), meshZ_Persistent7(binIndices_patient_loc{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Persistent7(binIndices_patient_loc{5}), meshY_Persistent7(binIndices_patient_loc{5}), meshZ_Persistent7(binIndices_patient_loc{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Persistent7(binIndices_patient_loc{6}), meshY_Persistent7(binIndices_patient_loc{6}), meshZ_Persistent7(binIndices_patient_loc{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Persistent7(binIndices_patient_loc{7}), meshY_Persistent7(binIndices_patient_loc{7}), meshZ_Persistent7(binIndices_patient_loc{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Persistent7(binIndices_patient_loc{8}), meshY_Persistent7(binIndices_patient_loc{8}), meshZ_Persistent7(binIndices_patient_loc{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Persistent7(binIndices_patient_loc{9}), meshY_Persistent7(binIndices_patient_loc{9}), meshZ_Persistent7(binIndices_patient_loc{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);







%PATCHY
 %Patient 7 Patchy 10% pattern 6
% Tissue with electrodes
patientNumToUse = 7;
patternNumToUse = 6;
degreeNumToUse = 10;
fibTypeNumToUse = 4;
loc_actTimesSim = calculateLATs_from_matrix([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNumToUse), '_', typesFibrosis{fibTypeNumToUse},'FibrosisNT', num2str(degreeNumToUse), '_pattern', num2str(patternNumToUse), '_2000points.csv' ], phie_dt, normalizeBool);

projLAT_patient_loc = projLATs_to_Surface(elecX_Persistent7, elecY_Persistent7, elecZ_Persistent7, meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, loc_actTimesSim);
figure
scatter3(meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, 190, projLAT_patient_loc, 'filled', 'MarkerEdgeColor','k');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (msec)';
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

numBins = 9;
[binEdges_loc, binIndices_patient_loc] = partitionToBins_Isochrones(projLAT_patient_loc, numBins, binEdges7 );

figure
plot3(meshX_Persistent7(binIndices_patient_loc{1}), meshY_Persistent7(binIndices_patient_loc{1}), meshZ_Persistent7(binIndices_patient_loc{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Persistent7(binIndices_patient_loc{2}), meshY_Persistent7(binIndices_patient_loc{2}), meshZ_Persistent7(binIndices_patient_loc{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Persistent7(binIndices_patient_loc{3}), meshY_Persistent7(binIndices_patient_loc{3}), meshZ_Persistent7(binIndices_patient_loc{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Persistent7(binIndices_patient_loc{4}), meshY_Persistent7(binIndices_patient_loc{4}), meshZ_Persistent7(binIndices_patient_loc{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Persistent7(binIndices_patient_loc{5}), meshY_Persistent7(binIndices_patient_loc{5}), meshZ_Persistent7(binIndices_patient_loc{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Persistent7(binIndices_patient_loc{6}), meshY_Persistent7(binIndices_patient_loc{6}), meshZ_Persistent7(binIndices_patient_loc{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Persistent7(binIndices_patient_loc{7}), meshY_Persistent7(binIndices_patient_loc{7}), meshZ_Persistent7(binIndices_patient_loc{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Persistent7(binIndices_patient_loc{8}), meshY_Persistent7(binIndices_patient_loc{8}), meshZ_Persistent7(binIndices_patient_loc{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Persistent7(binIndices_patient_loc{9}), meshY_Persistent7(binIndices_patient_loc{9}), meshZ_Persistent7(binIndices_patient_loc{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);




 %Patient 7 Patchy 35% pattern 6
% Tissue with electrodes
patientNumToUse = 7;
patternNumToUse = 6;
degreeNumToUse = 35;
fibTypeNumToUse = 4;
loc_actTimesSim = calculateLATs_from_matrix([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNumToUse), '_', typesFibrosis{fibTypeNumToUse},'FibrosisNT', num2str(degreeNumToUse), '_pattern', num2str(patternNumToUse), '_2000points.csv' ], phie_dt, normalizeBool);

projLAT_patient_loc = projLATs_to_Surface(elecX_Persistent7, elecY_Persistent7, elecZ_Persistent7, meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, loc_actTimesSim);
figure
scatter3(meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, 190, projLAT_patient_loc, 'filled', 'MarkerEdgeColor','k');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (msec)';
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

numBins = 9;
[binEdges_loc, binIndices_patient_loc] = partitionToBins_Isochrones(projLAT_patient_loc, numBins, binEdges7 );

figure
plot3(meshX_Persistent7(binIndices_patient_loc{1}), meshY_Persistent7(binIndices_patient_loc{1}), meshZ_Persistent7(binIndices_patient_loc{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Persistent7(binIndices_patient_loc{2}), meshY_Persistent7(binIndices_patient_loc{2}), meshZ_Persistent7(binIndices_patient_loc{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Persistent7(binIndices_patient_loc{3}), meshY_Persistent7(binIndices_patient_loc{3}), meshZ_Persistent7(binIndices_patient_loc{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Persistent7(binIndices_patient_loc{4}), meshY_Persistent7(binIndices_patient_loc{4}), meshZ_Persistent7(binIndices_patient_loc{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Persistent7(binIndices_patient_loc{5}), meshY_Persistent7(binIndices_patient_loc{5}), meshZ_Persistent7(binIndices_patient_loc{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Persistent7(binIndices_patient_loc{6}), meshY_Persistent7(binIndices_patient_loc{6}), meshZ_Persistent7(binIndices_patient_loc{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Persistent7(binIndices_patient_loc{7}), meshY_Persistent7(binIndices_patient_loc{7}), meshZ_Persistent7(binIndices_patient_loc{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Persistent7(binIndices_patient_loc{8}), meshY_Persistent7(binIndices_patient_loc{8}), meshZ_Persistent7(binIndices_patient_loc{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Persistent7(binIndices_patient_loc{9}), meshY_Persistent7(binIndices_patient_loc{9}), meshZ_Persistent7(binIndices_patient_loc{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);




 %Patient 7 Patchy 60% pattern 6
% Tissue with electrodes
patientNumToUse = 7;
patternNumToUse = 6;
degreeNumToUse = 60;
fibTypeNumToUse = 4;
loc_actTimesSim = calculateLATs_from_matrix([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNumToUse), '_', typesFibrosis{fibTypeNumToUse},'FibrosisNT', num2str(degreeNumToUse), '_pattern', num2str(patternNumToUse), '_2000points.csv' ], phie_dt, normalizeBool);

projLAT_patient_loc = projLATs_to_Surface(elecX_Persistent7, elecY_Persistent7, elecZ_Persistent7, meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, loc_actTimesSim);
figure
scatter3(meshX_Persistent7, meshY_Persistent7, meshZ_Persistent7, 190, projLAT_patient_loc, 'filled', 'MarkerEdgeColor','k');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (msec)';
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

numBins = 9;
[binEdges_loc, binIndices_patient_loc] = partitionToBins_Isochrones(projLAT_patient_loc, numBins, binEdges7 );

figure
plot3(meshX_Persistent7(binIndices_patient_loc{1}), meshY_Persistent7(binIndices_patient_loc{1}), meshZ_Persistent7(binIndices_patient_loc{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Persistent7(binIndices_patient_loc{2}), meshY_Persistent7(binIndices_patient_loc{2}), meshZ_Persistent7(binIndices_patient_loc{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Persistent7(binIndices_patient_loc{3}), meshY_Persistent7(binIndices_patient_loc{3}), meshZ_Persistent7(binIndices_patient_loc{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Persistent7(binIndices_patient_loc{4}), meshY_Persistent7(binIndices_patient_loc{4}), meshZ_Persistent7(binIndices_patient_loc{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Persistent7(binIndices_patient_loc{5}), meshY_Persistent7(binIndices_patient_loc{5}), meshZ_Persistent7(binIndices_patient_loc{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Persistent7(binIndices_patient_loc{6}), meshY_Persistent7(binIndices_patient_loc{6}), meshZ_Persistent7(binIndices_patient_loc{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Persistent7(binIndices_patient_loc{7}), meshY_Persistent7(binIndices_patient_loc{7}), meshZ_Persistent7(binIndices_patient_loc{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Persistent7(binIndices_patient_loc{8}), meshY_Persistent7(binIndices_patient_loc{8}), meshZ_Persistent7(binIndices_patient_loc{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Persistent7(binIndices_patient_loc{9}), meshY_Persistent7(binIndices_patient_loc{9}), meshZ_Persistent7(binIndices_patient_loc{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);











 %Patient 6 Patchy 35% pattern 6
% Tissue with electrodes
patientNumToUse = 6;
patternNumToUse = 6;
degreeNumToUse = 35;
fibTypeNumToUse = 4;
loc_actTimesSim = calculateLATs_from_matrix([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNumToUse), '_', typesFibrosis{fibTypeNumToUse},'FibrosisNT', num2str(degreeNumToUse), '_pattern', num2str(patternNumToUse), '_2000points.csv' ], phie_dt, normalizeBool);

projLAT_patient_loc = projLATs_to_Surface(elecX_Persistent6, elecY_Persistent6, elecZ_Persistent6, meshX_Persistent6, meshY_Persistent6, meshZ_Persistent6, loc_actTimesSim);
figure
scatter3(meshX_Persistent6, meshY_Persistent6, meshZ_Persistent6, 190, projLAT_patient_loc, 'filled', 'MarkerEdgeColor','k');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (msec)';
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

numBins = 9;
[binEdges_loc, binIndices_patient_loc] = partitionToBins_Isochrones(projLAT_patient_loc, numBins, binEdges6 );

figure
plot3(meshX_Persistent6(binIndices_patient_loc{1}), meshY_Persistent6(binIndices_patient_loc{1}), meshZ_Persistent6(binIndices_patient_loc{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Persistent6(binIndices_patient_loc{2}), meshY_Persistent6(binIndices_patient_loc{2}), meshZ_Persistent6(binIndices_patient_loc{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Persistent6(binIndices_patient_loc{3}), meshY_Persistent6(binIndices_patient_loc{3}), meshZ_Persistent6(binIndices_patient_loc{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Persistent6(binIndices_patient_loc{4}), meshY_Persistent6(binIndices_patient_loc{4}), meshZ_Persistent6(binIndices_patient_loc{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Persistent6(binIndices_patient_loc{5}), meshY_Persistent6(binIndices_patient_loc{5}), meshZ_Persistent6(binIndices_patient_loc{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Persistent6(binIndices_patient_loc{6}), meshY_Persistent6(binIndices_patient_loc{6}), meshZ_Persistent6(binIndices_patient_loc{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Persistent6(binIndices_patient_loc{7}), meshY_Persistent6(binIndices_patient_loc{7}), meshZ_Persistent6(binIndices_patient_loc{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Persistent6(binIndices_patient_loc{8}), meshY_Persistent6(binIndices_patient_loc{8}), meshZ_Persistent6(binIndices_patient_loc{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Persistent6(binIndices_patient_loc{9}), meshY_Persistent6(binIndices_patient_loc{9}), meshZ_Persistent6(binIndices_patient_loc{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);





 %Patient 8 Patchy 35% pattern 6
% Tissue with electrodes
patientNumToUse = 8;
patternNumToUse = 6;
degreeNumToUse = 35;
fibTypeNumToUse = 4;
loc_actTimesSim = calculateLATs_from_matrix([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNumToUse), '_', typesFibrosis{fibTypeNumToUse},'FibrosisNT', num2str(degreeNumToUse), '_pattern', num2str(patternNumToUse), '_2000points.csv' ], phie_dt, normalizeBool);

projLAT_patient_loc = projLATs_to_Surface(elecX_Persistent8, elecY_Persistent8, elecZ_Persistent8, meshX_Persistent8, meshY_Persistent8, meshZ_Persistent8, loc_actTimesSim);
figure
scatter3(meshX_Persistent8, meshY_Persistent8, meshZ_Persistent8, 190, projLAT_patient_loc, 'filled', 'MarkerEdgeColor','k');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (msec)';
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

numBins = 9;
[binEdges_loc, binIndices_patient_loc] = partitionToBins_Isochrones(projLAT_patient_loc, numBins, binEdges8 );

figure
plot3(meshX_Persistent8(binIndices_patient_loc{1}), meshY_Persistent8(binIndices_patient_loc{1}), meshZ_Persistent8(binIndices_patient_loc{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Persistent8(binIndices_patient_loc{2}), meshY_Persistent8(binIndices_patient_loc{2}), meshZ_Persistent8(binIndices_patient_loc{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Persistent8(binIndices_patient_loc{3}), meshY_Persistent8(binIndices_patient_loc{3}), meshZ_Persistent8(binIndices_patient_loc{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Persistent8(binIndices_patient_loc{4}), meshY_Persistent8(binIndices_patient_loc{4}), meshZ_Persistent8(binIndices_patient_loc{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Persistent8(binIndices_patient_loc{5}), meshY_Persistent8(binIndices_patient_loc{5}), meshZ_Persistent8(binIndices_patient_loc{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Persistent8(binIndices_patient_loc{6}), meshY_Persistent8(binIndices_patient_loc{6}), meshZ_Persistent8(binIndices_patient_loc{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Persistent8(binIndices_patient_loc{7}), meshY_Persistent8(binIndices_patient_loc{7}), meshZ_Persistent8(binIndices_patient_loc{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Persistent8(binIndices_patient_loc{8}), meshY_Persistent8(binIndices_patient_loc{8}), meshZ_Persistent8(binIndices_patient_loc{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Persistent8(binIndices_patient_loc{9}), meshY_Persistent8(binIndices_patient_loc{9}), meshZ_Persistent8(binIndices_patient_loc{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);



 %Patient 9 Patchy 35% pattern 6
% Tissue with electrodes
patientNumToUse = 9;
patternNumToUse = 6;
degreeNumToUse = 35;
fibTypeNumToUse = 4;
loc_actTimesSim = calculateLATs_from_matrix([ '../SimResults/P1simsStuff/tempFolder/FibrosisSimulationsVes/Ve_Catheter_MonoAp_EXPLICIT_INTRACELLULAR_Courtemanche_monodomain_EI_Persistent', num2str(patientNumToUse), '_', typesFibrosis{fibTypeNumToUse},'FibrosisNT', num2str(degreeNumToUse), '_pattern', num2str(patternNumToUse), '_2000points.csv' ], phie_dt, normalizeBool);

projLAT_patient_loc = projLATs_to_Surface(elecX_Persistent9, elecY_Persistent9, elecZ_Persistent9, meshX_Persistent9, meshY_Persistent9, meshZ_Persistent9, loc_actTimesSim);
figure
scatter3(meshX_Persistent9, meshY_Persistent9, meshZ_Persistent9, 190, projLAT_patient_loc, 'filled', 'MarkerEdgeColor','k');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (msec)';
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

numBins = 9;
[binEdges_loc, binIndices_patient_loc] = partitionToBins_Isochrones(projLAT_patient_loc, numBins, binEdges9 );

figure
plot3(meshX_Persistent9(binIndices_patient_loc{1}), meshY_Persistent9(binIndices_patient_loc{1}), meshZ_Persistent9(binIndices_patient_loc{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(meshX_Persistent9(binIndices_patient_loc{2}), meshY_Persistent9(binIndices_patient_loc{2}), meshZ_Persistent9(binIndices_patient_loc{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(meshX_Persistent9(binIndices_patient_loc{3}), meshY_Persistent9(binIndices_patient_loc{3}), meshZ_Persistent9(binIndices_patient_loc{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(meshX_Persistent9(binIndices_patient_loc{4}), meshY_Persistent9(binIndices_patient_loc{4}), meshZ_Persistent9(binIndices_patient_loc{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(meshX_Persistent9(binIndices_patient_loc{5}), meshY_Persistent9(binIndices_patient_loc{5}), meshZ_Persistent9(binIndices_patient_loc{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(meshX_Persistent9(binIndices_patient_loc{6}), meshY_Persistent9(binIndices_patient_loc{6}), meshZ_Persistent9(binIndices_patient_loc{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(meshX_Persistent9(binIndices_patient_loc{7}), meshY_Persistent9(binIndices_patient_loc{7}), meshZ_Persistent9(binIndices_patient_loc{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(meshX_Persistent9(binIndices_patient_loc{8}), meshY_Persistent9(binIndices_patient_loc{8}), meshZ_Persistent9(binIndices_patient_loc{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(meshX_Persistent9(binIndices_patient_loc{9}), meshY_Persistent9(binIndices_patient_loc{9}), meshZ_Persistent9(binIndices_patient_loc{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);








%%


% 
% 
% fibPoints_Persisitent9 = csvread('../SimResults/P1simsStuff/FibrosisPatterns/P1_Persistent9_CompactFibrosis_P1_10perc.csv');
% fibX_Persisitent9 = fibPoints_Persisitent9(:,3); fibY_Persisitent9 = fibPoints_Persisitent9(:,4); fibZ_Persisitent9 = fibPoints_Persisitent9(:,5);
% 
% plot3(meshX_Persisitent9, meshY_Persisitent9, meshZ_Persisitent9,'yo')
% hold on
% plot3(fibX_Persisitent9, fibY_Persisitent9, fibZ_Persisitent9,'ro')



%CHANGE THE DISTANCES OF THE ELECTRODE GRID TO GET SIMILAR DISTANCES TO
%PATIENT DATA AND VARY IT TO GET DIFFERENT DISTANCES. ALSO, GET PEAK TO
%PEAK AMPLITUDE AS ANOTHER VARIABLE TO INCLUDE IN MATRIX OF VARIABLES. KEEP
%DOWNSTROKE, THE TWO UPSTROKES, SIGNAL WIDTH, BASELINE, EGM DURATION, AND
%DISTANCE AS VARIABLES FOR THE MATRIX. USE FASTER CV NEW SIMS

%FOR NOW, TALK ABOUT WAVELET DECOMPOSITION ONE AND TABLES. YOU COULD TALK
%ABOUT VARIABILITY OF SIGNALS AS WELL

