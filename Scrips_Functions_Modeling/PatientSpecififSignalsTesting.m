

% % CHECKING SIGNALS IN CLUSTER
intervalCheck = 4;
figure
subplot(2,2,1)
plot(AllSignals_Persistent(:, indexCluster1Persistent_Var(1:intervalCheck:end))) 
hold on
% plot(simDataBi(1:100, 1:2:end), 'k--','LineWidth', 2)
title('1')

subplot(2,2,2)
plot(AllSignals_Persistent(:, indexCluster2Persistent_Var(1:intervalCheck:end)))
hold on
plot(AllSignals_Persistent(:, indexCluster4Persistent_Var(1:intervalCheck:end)))
% plot(simDataBi(1:100, 1:2:end), 'k--','LineWidth', 2)
title('2')

subplot(2,2,3)
plot(AllSignals_Persistent(:, indexCluster3Persistent_Var(1:intervalCheck:end)))
hold on
plot(AllSignals_Persistent(:, indexCluster6Persistent_Var(1:intervalCheck:end)))
plot(AllSignals_Persistent(:, indexCluster5Persistent_Var(1:intervalCheck:end)))
% plot(simDataBi(1:100, 1:2:end), 'k--','LineWidth', 2)
title('3')

subplot(2,2,4)
plot(AllSignals_Persistent(:, indexCluster7Persistent_Var(1:intervalCheck:end)))
hold on
plot(AllSignals_Persistent(:, indexCluster8Persistent_Var(1:intervalCheck:end)))
% plot(simDataBiCut(:, find(idxPCAnew_Persistent_Bi == 8)), 'k-','LineWidth', 2)
% plot(simDataMonoCut(:, find(idxPCAnew_Persistent_Mono == 8)), 'g-','LineWidth', 2)
title('4')

sgti = sgtitle('Persistent AFib sample signals in clusters') ;
sgti.FontSize = 30;










figure
subplot(2,2,1)
plot(AllSignals_Paroxysmal(:, indexCluster1Paroxysmal_Var(1:intervalCheck:end)))
hold on
plot(AllSignals_Paroxysmal(:, indexCluster3Paroxysmal_Var(1:intervalCheck:end)))
plot(AllSignals_Paroxysmal(:, indexCluster6Paroxysmal_Var(1:intervalCheck:end)))
% plot(simDataBi(1:100, 1:2:end), 'k--','LineWidth', 2)
title('1')

subplot(2,2,2)
plot(AllSignals_Paroxysmal(:, indexCluster2Paroxysmal_Var(1:intervalCheck:end)))
hold on
plot(AllSignals_Paroxysmal(:, indexCluster7Paroxysmal_Var(1:intervalCheck:end)))
% plot(simDataBiCut(:, find(idxPCAnew_Paroxysmal_Bi == 2)), 'k-','LineWidth', 2)
% plot(simDataMonoCut(:, find(idxPCAnew_Paroxysmal_Mono == 2)), 'g-','LineWidth', 2)
title('2')

subplot(2,2,3)
plot(AllSignals_Paroxysmal(:, indexCluster4Paroxysmal_Var(1:intervalCheck:end)))
hold on
plot(AllSignals_Paroxysmal(:, indexCluster8Paroxysmal_Var(1:intervalCheck:end)))
% plot(simDataBi(1:100, 1:2:end), 'k--','LineWidth', 2)
title('3')

subplot(2,2,4)
plot(AllSignals_Paroxysmal(:, indexCluster5Paroxysmal_Var(1:intervalCheck:end)))
hold on
% plot(simDataBi(1:100, 1:2:end), 'k--','LineWidth', 2)
title('4')

sgti = sgtitle('Paroxysmal AFib sample signals in clusters') ;
sgti.FontSize = 30;






























% % CHECKING SIGNALS IN CLUSTER
intervalCheck = 10;
figure
subplot(2,2,1)
plot(AllSignals_Persistent(:, indexCluster1Persistent_Var(1:intervalCheck:end)),'LineWidth',2) 
hold on
% plot(simDataBi(1:100, 1:2:end), 'k--','LineWidth', 2)
title('1')
set(gca,'FontSize',20)
subplot(2,2,2)
plot(AllSignals_Persistent(:, indexCluster2Persistent_Var(1:intervalCheck:end)),'LineWidth',2)
% plot(simDataBi(1:100, 1:2:end), 'k--','LineWidth', 2)
title('2')
set(gca,'FontSize',20)
subplot(2,2,3)
plot(AllSignals_Persistent(:, indexCluster3Persistent_Var(1:intervalCheck:end)),'LineWidth',2)
hold on
% plot(simDataBi(1:100, 1:2:end), 'k--','LineWidth', 2)
title('3')
set(gca,'FontSize',20)
subplot(2,2,4)
plot(AllSignals_Persistent(:, indexCluster4Persistent_Var(1:intervalCheck:end)),'LineWidth',2)
hold on
% plot(simDataBiCut(:, find(idxPCAnew_Persistent_Bi == 8)), 'k-','LineWidth', 2)
% plot(simDataMonoCut(:, find(idxPCAnew_Persistent_Mono == 8)), 'g-','LineWidth', 2)
title('4')
set(gca,'FontSize',20)
sgti = sgtitle('Persistent AFib sample signals in clusters') ;
sgti.FontSize = 30;




figure
subplot(2,2,1)
plot(AllSignals_Paroxysmal(:, indexCluster1Paroxysmal_Var(1:intervalCheck:end)),'LineWidth',2)
hold on
% plot(simDataBi(1:100, 1:2:end), 'k--','LineWidth', 2)
title('1')
set(gca,'FontSize',20)
subplot(2,2,2)
plot(AllSignals_Paroxysmal(:, indexCluster2Paroxysmal_Var(1:intervalCheck:end)),'LineWidth',2)
hold on
% plot(simDataBiCut(:, find(idxPCAnew_Paroxysmal_Bi == 2)), 'k-','LineWidth', 2)
% plot(simDataMonoCut(:, find(idxPCAnew_Paroxysmal_Mono == 2)), 'g-','LineWidth', 2)
title('2')
set(gca,'FontSize',20)
subplot(2,2,3)
plot(AllSignals_Paroxysmal(:, indexCluster3Paroxysmal_Var(1:intervalCheck:end)),'LineWidth',2)
hold on
% plot(simDataBi(1:100, 1:2:end), 'k--','LineWidth', 2)
title('3')
set(gca,'FontSize',20)
subplot(2,2,4)
plot(AllSignals_Paroxysmal(:, indexCluster4Paroxysmal_Var(1:intervalCheck:end)),'LineWidth',2)
hold on
% plot(simDataBi(1:100, 1:2:end), 'k--','LineWidth', 2)
title('4')
set(gca,'FontSize',20)
sgti = sgtitle('Paroxysmal AFib sample signals in clusters') ;
sgti.FontSize = 30;














inttttt = 10;
figure
subplot(2,1,1)
plot(mean(AllSignals_Persistent(:, indexCluster1Persistent_Var(1:inttttt:end)),2), 'LineWidth', 4 )
hold on
plot(mean(AllSignals_Persistent(:, indexCluster2Persistent_Var(1:inttttt:end)),2), 'LineWidth', 4 )
plot(mean(AllSignals_Persistent(:, indexCluster3Persistent_Var(1:inttttt:end)),2), 'LineWidth', 4 )
plot(mean(AllSignals_Persistent(:, indexCluster4Persistent_Var(1:inttttt:end)),2), 'LineWidth', 4 )
legend('1','2','3','4')
title('Persistent')
subplot(2,1,2)
plot(mean(AllSignals_Paroxysmal(:, indexCluster1Paroxysmal_Var(1:inttttt:end)),2), 'LineWidth', 4 )
hold on
plot(mean(AllSignals_Paroxysmal(:, indexCluster2Paroxysmal_Var(1:inttttt:end)),2), 'LineWidth', 4 )
plot(mean(AllSignals_Paroxysmal(:, indexCluster3Paroxysmal_Var(1:inttttt:end)),2), 'LineWidth', 4 )
plot(mean(AllSignals_Paroxysmal(:, indexCluster4Paroxysmal_Var(1:inttttt:end)),2), 'LineWidth', 4 )
legend('1','2','3','4')
title('Paroxysmal')






figure
subplot(2,1,1)
plot(AllSignals_Persistent(:, indexCluster1Persistent_Var(10)),'*--','DisplayName','1','LineWidth',3)
hold on
plot(AllSignals_Persistent(:, indexCluster2Persistent_Var(20)),'v--','DisplayName','2','LineWidth',3)
plot(AllSignals_Persistent(:, indexCluster3Persistent_Var(25)),'+--','DisplayName','3','LineWidth',3)
plot(AllSignals_Persistent(:, indexCluster4Persistent_Var(10)),'o--','DisplayName','4','LineWidth',3)
xlabel('time (ms)')
ylabel('normalized potential')
title('Persistent sample ojective signals: option 1')
set(gca,'FontSize',20)
legend()
subplot(2,1,2)
plot(AllSignals_Paroxysmal(:, indexCluster1Paroxysmal_Var(15)),'*--','DisplayName','1','LineWidth',3)
hold on
plot(AllSignals_Paroxysmal(:, indexCluster2Paroxysmal_Var(15)),'v--','DisplayName','2','LineWidth',3)
plot(AllSignals_Paroxysmal(:, indexCluster3Paroxysmal_Var(20)),'+--','DisplayName','3','LineWidth',3)
plot(AllSignals_Paroxysmal(:, indexCluster4Paroxysmal_Var(15)),'o--','DisplayName','4','LineWidth',3)
xlabel('time (ms)')
ylabel('normalized potential')
title('Paroxysmal sample ojective signals: option 1')
set(gca,'FontSize',20)
legend()









%OPTION 2
intervalCheck = 7;
figure;
subplot(3,1,1)
plot(signalsAllNarrow_Persistent(:,1:intervalCheck:end),'LineWidth',2); 
hold on
% plot(simDataBi(1:100, 1:2:end), 'k--','LineWidth', 2)
title('Narrow Upstroke')
set(gca,'FontSize',20)
subplot(3,1,2)
plot(signalsAllInBetween_Persistent(:,1:intervalCheck:end),'LineWidth',2); 
hold on
% plot(simDataBi(1:100, 1:2:end), 'k--','LineWidth', 2)
title('In-between Upstroke')
set(gca,'FontSize',20)
subplot(3,1,3)
plot(signalsAllWide_Persistent(:,1:intervalCheck:end),'LineWidth',2); 
hold on
% plot(simDataBi(1:100, 1:2:end), 'k--','LineWidth', 2)
title('Wide Upstroke')
set(gca,'FontSize',20)
sgti = sgtitle('Persistent AFib sample signals in clusters: option 2') ;
sgti.FontSize = 30;



figure;
subplot(3,1,1)
plot(signalsAllNarrow_Paroxysmal(:,1:intervalCheck:end),'LineWidth',2); 
hold on
% plot(simDataBi(1:100, 1:2:end), 'k--','LineWidth', 2)
title('Narrow Upstroke')
set(gca,'FontSize',20)
subplot(3,1,2)
plot(signalsAllInBetween_Paroxysmal(:,1:intervalCheck:end),'LineWidth',2); 
hold on
% plot(simDataBi(1:100, 1:2:end), 'k--','LineWidth', 2)
title('In-between Upstroke')
set(gca,'FontSize',20)
subplot(3,1,3)
plot(signalsAllWide_Paroxysmal(:,1:intervalCheck:end),'LineWidth',2); 
hold on
% plot(simDataBi(1:100, 1:2:end), 'k--','LineWidth', 2)
title('Wide Upstroke')
set(gca,'FontSize',20)
sgti = sgtitle('Paroxysmal AFib sample signals in clusters: option 2') ;
sgti.FontSize = 30;











figure
subplot(2,1,1)
plot(signalsAllNarrow_Persistent(:, 10),'*--','DisplayName','Narrow','LineWidth',3)
hold on
plot(signalsAllInBetween_Persistent(:, 20),'v--','DisplayName','In-between','LineWidth',3)
plot(signalsAllWide_Persistent(:, 25),'o--','DisplayName','Wide','LineWidth',3)
xlabel('time (ms)')
ylabel('normalized potential')
title('Persistent sample ojective signals: option 2')
set(gca,'FontSize',20)
legend()
subplot(2,1,2)
plot(signalsAllNarrow_Paroxysmal(:, 10),'*--','DisplayName','Narrow','LineWidth',3)
hold on
plot(signalsAllInBetween_Paroxysmal(:, 15),'v--','DisplayName','In-between','LineWidth',3)
plot(signalsAllWide_Paroxysmal(:, 14),'o--','DisplayName','Wide','LineWidth',3)
xlabel('time (ms)')
ylabel('normalized potential')
title('Paroxysmal sample ojective signals: option 2')
set(gca,'FontSize',20)
legend()



tableOption2 = [ size(signalsAllNarrow_Paroxysmal,2)./size(AllSignals_Paroxysmal,2), size(signalsAllInBetween_Paroxysmal,2)./size(AllSignals_Paroxysmal,2), size(signalsAllWide_Paroxysmal,2)./size(AllSignals_Paroxysmal,2);...
                 size(signalsAllNarrow_Persistent,2)./size(AllSignals_Persistent,2), size(signalsAllInBetween_Persistent,2)./size(AllSignals_Persistent,2), size(signalsAllWide_Persistent,2)./size(AllSignals_Persistent,2)].*100



tableOption1 = [ size(indexCluster1Paroxysmal_Var,1)./size(AllSignals_Paroxysmal,2), size(indexCluster2Paroxysmal_Var,1)./size(AllSignals_Paroxysmal,2), size(indexCluster3Paroxysmal_Var,1)./size(AllSignals_Paroxysmal,2), size(indexCluster4Paroxysmal_Var,1)./size(AllSignals_Paroxysmal,2);...
                 size(indexCluster1Persistent_Var,1)./size(AllSignals_Persistent,2), size(indexCluster2Persistent_Var,1)./size(AllSignals_Persistent,2), size(indexCluster3Persistent_Var,1)./size(AllSignals_Persistent,2), size(indexCluster4Persistent_Var,1)./size(AllSignals_Persistent,2)].*100


%%

%Hotelling's T-Test square for statistically significant difference between
%Persistent and Paroxysmal AFib
alphaInterest = .05;

X_bothPP = [];
dataToClassifyAllButOthers_Paroxysmal_less = dataToClassifyAllButOthers_Paroxysmal;
dataToClassifyAllButOthers_Persistent_less = dataToClassifyAllButOthers_Persistent;


% dataToClassifyAllButOthers_Paroxysmal_less = dataToClassifyAllButOthers_Paroxysmal_less(1:size(dataToClassifyAllButOthers_Persistent_less,1));

% num = 0;
% for i = 1:size(dataToClassifyAllButOthers_Persistent_less,1)
%     num = num + ( dataToClassifyAllButOthers_Paroxysmal_less(i,:) - mean(dataToClassifyAllButOthers_Paroxysmal_less,1) ) .* ( dataToClassifyAllButOthers_Persistent_less(i,:) - mean(dataToClassifyAllButOthers_Persistent_less,1) );
% end


desiredNumberSamples = 300;
idxParoxysmal_fewer = randperm(numel(dataToClassifyAllButOthers_Paroxysmal_less(:,1)));
dataToClassifyAllButOthers_Paroxysmal_less = dataToClassifyAllButOthers_Paroxysmal_less(idxParoxysmal_fewer(1:desiredNumberSamples),:);
idxPersistent_fewer = randperm(numel(dataToClassifyAllButOthers_Persistent_less(:,1)));
dataToClassifyAllButOthers_Persistent_less = dataToClassifyAllButOthers_Persistent_less(idxPersistent_fewer(1:desiredNumberSamples),:);

[T_squared, F_stats, df1, df2, p_value_determined] = HotellingCustom(dataToClassifyAllButOthers_Paroxysmal_less, dataToClassifyAllButOthers_Persistent_less);

for i = 1:1:size(dataToClassifyAllButOthers_Paroxysmal_less,1)
    X_bothPP = [ X_bothPP; 1, dataToClassifyAllButOthers_Paroxysmal_less(i,:) ];
    % dataToClassifyAllButOthers_Paroxysmal_less = [ dataToClassifyAllButOthers_Paroxysmal_less; dataToClassifyAllButOthers_Paroxysmal(i,:)];
end
for i = 1:1:size(dataToClassifyAllButOthers_Persistent_less,1)
    X_bothPP = [ X_bothPP; 2, dataToClassifyAllButOthers_Persistent_less(i,:) ];
    % dataToClassifyAllButOthers_Persistent_less = [ dataToClassifyAllButOthers_Persistent_less; dataToClassifyAllButOthers_Persistent(i,:)];
end

% X_bothPP = X_bothPP(7:end,:);
[MBox] = MBoxtest(X_bothPP,alphaInterest);

[T2Hot2iheOUT] = T2Hot2ihe(X_bothPP,alphaInterest);


if p_value_determined < alphaInterest
    disp("we CAN reject the null hypothesis that the mean vector for sample 1 equals the mean vector for sample 2 given the evidence as usual:")
    varparam =  "T^2: "  + num2str(T_squared) + "; F: " + num2str(F_stats) + "; d.f.: " + num2str(df1) + "," + num2str(df2) + "; p < alpha : " + num2str(p_value_determined) + " < " + num2str(alphaInterest);
    disp(varparam)

else
    disp("we CANNOT reject the null hypothesis that the mean vector for sample 1 equals the mean vector for sample 2 given the evidence as usual:")
    varparam =  "T^2: "  + num2str(T_squared) + "; F: " + num2str(F_stats) + "; d.f.: " + num2str(df1) + "," + num2str(df2) + "; p > alpha : " + num2str(p_value_determined) + " > " + num2str(alphaInterest);
    disp(varparam)
end




if T2Hot2iheOUT(end) < alphaInterest %&& p_value_determined < alphaInterest
    disp("we CAN reject the null hypothesis that the mean vector for sample 1 equals the mean vector for sample 2 given the evidence as usual:")
    varparam =  "T^2: "  + T2Hot2iheOUT(4) + "; F: " + T2Hot2iheOUT(5) + "; d.f.: " + num2str(df1) + "," + num2str(df2) + "; p < alpha : " + T2Hot2iheOUT(end) + " < " + num2str(alphaInterest);
    disp(varparam)

else
    disp("we CANNOT reject the null hypothesis that the mean vector for sample 1 equals the mean vector for sample 2 given the evidence as usual:")
    varparam =  "T^2: "  + T2Hot2iheOUT(4) + "; F: " + T2Hot2iheOUT(5) + "; d.f.: " + num2str(df1) + "," + num2str(df2) + "; p > alpha : " + T2Hot2iheOUT(end) + " > " + num2str(alphaInterest);
    disp(varparam)
end






%%

figure
testvectorsignals = [];
numElecParox = 0;
numElecPersist = 0;

for m = 1:10
patie



end

%% TRYING TO GUESSTIMATE CV

Point1rand = [-4.41, -4.14, 8.19];
Point2rand = [-4.74, -1.91, 9.49];
LAT1rand = 68;
LAT2rand = 97;

CVguesstimate = (sqrt( (Point1rand(1) - Point2rand(1)).^2 + (Point1rand(2) - Point2rand(2)).^2 + (Point1rand(3) - Point2rand(3)).^2  )/abs(LAT1rand - LAT2rand))*1000 %CM/S


%% PAROXYSMAL PATIENTS CHECKING


%PATIENTS 7 AND 9 FOR PERSISTENT MESHES AND 2 AND 3 FOR PAROXYSMAL MESHES



%CHECKING THINGS


%Simulated signals
simulatedParoxysmal2 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_Paroxysmal2_NoFibrosis_2000points.csv");
simulatedParoxysmal2_downsampled = UpsampleAndFilterFunction(simulatedParoxysmal2', .1, 5, 200, 1000, 10,false)';
actTimesSim2 = [];

simulatedParoxysmal3 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_Paroxysmal3_NoFibrosis_2000points.csv");
simulatedParoxysmal3_downsampled = UpsampleAndFilterFunction(simulatedParoxysmal3', .1, 5, 200, 1000, 10,false)';
actTimesSim3 = [];

simulatedParoxysmal2_downsampled_norm = [];
simulatedParoxysmal3_downsampled_norm = [];

% simulatedParoxysmal2_NF = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Paroxysmal2_faster_NoFibers_1000points.csv");
% simulatedParoxysmal2_NF_downsampled = UpsampleAndFilterFunction(simulatedParoxysmal2_NF', .1, 5, 200, 1000, 10,false)';
% actTimesSim_NF = [];

for i = 1:length(simulatedParoxysmal2(1,:))
    dVvec = diff(simulatedParoxysmal2_downsampled(:,i));
    minVal =  min(dVvec);
    indFind = find(dVvec == minVal);
    indFind = indFind(1);
    actTimesSim2(i) = indFind; 

    simulatedParoxysmal2_downsampled_norm(:,i) = ( simulatedParoxysmal2_downsampled(:,i) - min(simulatedParoxysmal2_downsampled(:,i)) )./( max(simulatedParoxysmal2_downsampled(:,i)) - min(simulatedParoxysmal2_downsampled(:,i)) );
    % dVvec2 = diff(simulatedParoxysmal2_NF_downsampled(:,i));
    % minVal2 =  min(dVvec2);
    % indFind2 = find(dVvec2 == minVal2);
    % if length(indFind2) > 0
    %     indFind2 = indFind2(1);
    % else
    %     indFind2 = 0;
    % end
    % actTimesSim_NF(i) = indFind2; 
end
for i = 1:length(simulatedParoxysmal3(1,:))
    dVvec2 = diff(simulatedParoxysmal3_downsampled(:,i));
    minVal2 =  min(dVvec2);
    indFind2 = find(dVvec2 == minVal2);
    if length(indFind2) > 0
        indFind2 = indFind2(1);
    else
        indFind2 = 0;
    end
    actTimesSim3(i) = indFind2; 

    simulatedParoxysmal3_downsampled_norm(:,i) = ( simulatedParoxysmal3_downsampled(:,i) - min(simulatedParoxysmal3_downsampled(:,i)) )./( max(simulatedParoxysmal3_downsampled(:,i)) - min(simulatedParoxysmal3_downsampled(:,i)) );
end


actTimes2 = expDB.dataSet{2}.map.continuousMapping.revisedLat(:,1);
actTimes3 = expDB.dataSet{3}.map.continuousMapping.revisedLat(:,1);


clear colorbar;


figure
subplot(2,1,1)

plot3(patient2xSurfaceCm, patient2ySurfaceCm, patient2zSurfaceCm, 'ro', 'DisplayName', 'tissue', 'MarkerFaceColor','r', 'MarkerSize',20,'LineWidth',3);
hold on
scatter3(patient2xCm, patient2yCm, patient2zCm, 230, actTimes2, 'filled', 'MarkerEdgeColor','k','LineWidth',3);
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (ms)';
% colormap(jet);
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

subplot(2,1,2)

plot3(patient2xSurfaceCm, patient2ySurfaceCm, patient2zSurfaceCm, 'ro', 'DisplayName', 'tissue', 'MarkerFaceColor','r', 'MarkerSize',20,'LineWidth',3);
hold on
scatter3(patient2xCm, patient2yCm, patient2zCm, 230, actTimesSim2, 'filled', 'MarkerEdgeColor','k','LineWidth',3);
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (ms)';
% colormap(jet);
caxis([40, 130]);
% title('H1A Simulation');
set(gca, 'FontSize', 35);




figure
subplot(2,1,1)

plot3(patient3xSurfaceCm, patient3ySurfaceCm, patient3zSurfaceCm, 'ro', 'DisplayName', 'tissue', 'MarkerFaceColor','r', 'MarkerSize',20,'LineWidth',3);
hold on
scatter3(patient3xCm, patient3yCm, patient3zCm, 230, actTimes3, 'filled', 'MarkerEdgeColor','k','LineWidth',3);
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (ms)';
% colormap(jet);
caxis([20, 110]);
% title('H1B Patient');
set(gca, 'FontSize', 35);

subplot(2,1,2)

plot3(patient3xSurfaceCm, patient3ySurfaceCm, patient3zSurfaceCm, 'ro', 'DisplayName', 'tissue', 'MarkerFaceColor','r', 'MarkerSize',20,'LineWidth',3);
hold on
scatter3(patient3xCm, patient3yCm, patient3zCm, 230, actTimesSim3, 'filled', 'MarkerEdgeColor','k','LineWidth',3);
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (ms)';
% colormap(jet);
caxis([20, 110]);
% title('H1B Simulation');
set(gca, 'FontSize', 35);






sameLATindex2 = [];
sameLATCorr2 = [];
sameLATindex3 = [];
sameLATCorr3 = [];

threshLATallowed = 10;

for i = 1:length(patient2xCm)
    if abs(actTimes2(i) - actTimesSim2(i)) < threshLATallowed
        sameLATindex2 = [ sameLATindex2, i ];

        actSignal2 = squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(i,1,1:200));
        correlation_coefficient = corrcoef(actSignal2, simulatedParoxysmal2_downsampled(:,i));
        similarity_score = correlation_coefficient(1, 2);
        sameLATCorr2 = [sameLATCorr2, similarity_score];
    end
end
for i = 1:length(patient3xCm)
    if abs(actTimes3(i) - actTimesSim3(i)) < threshLATallowed
        sameLATindex3 = [ sameLATindex3, i ];

        actSignal3 = squeeze(expDB.dataSet{3}.map.continuousMapping.fullEgms(i,1,1:200));
        correlation_coefficient = corrcoef(actSignal3, simulatedParoxysmal3_downsampled(:,i));
        similarity_score = correlation_coefficient(1, 2);
        sameLATCorr3 = [sameLATCorr3, similarity_score];
    end
end

percentageSameLATs2 = (length(sameLATCorr2)/length(patient2xCm))*100
percentageSameLATs3 = (length(sameLATCorr3)/length(patient3xCm))*100

meanCorr2 = mean(sameLATCorr2)
meanCorr3 = mean(sameLATCorr3)



figure;
subplot(2,2,1)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(174,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedParoxysmal2_downsampled(:,174), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 1")
subplot(2,2,2)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(175,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedParoxysmal2_downsampled(:,175), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 2")
subplot(2,2,3)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(139,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedParoxysmal2_downsampled(:,139), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 3")
subplot(2,2,4)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(177,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedParoxysmal2_downsampled(:,177), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 4")
sgtitle('Patient 2: Other signals')






figure;
subplot(2,2,1)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{3}.map.continuousMapping.fullEgms(93,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedParoxysmal3_downsampled(:,93), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 1")
subplot(2,2,2)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{3}.map.continuousMapping.fullEgms(165,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedParoxysmal3_downsampled(:,165), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 2")
subplot(2,2,3)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{3}.map.continuousMapping.fullEgms(232,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedParoxysmal3_downsampled(:,232), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 3")
subplot(2,2,4)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{3}.map.continuousMapping.fullEgms(630,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedParoxysmal3_downsampled(:,630), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 4")
sgtitle('Patient 3: Other signals')





figure;
subplot(2,2,1)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{3}.map.continuousMapping.fullEgms(175,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedParoxysmal3_downsampled(:,175), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 5")
subplot(2,2,2)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{3}.map.continuousMapping.fullEgms(343,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedParoxysmal3_downsampled(:,343), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 6")
subplot(2,2,3)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{3}.map.continuousMapping.fullEgms(832,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedParoxysmal3_downsampled(:,832), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 7")
subplot(2,2,4)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{3}.map.continuousMapping.fullEgms(469,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedParoxysmal3_downsampled(:,469), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 8")
sgtitle('Patient 3: Other signals')


%% MAPPING LAT TO TISSUE AND ISOCHRONES


projLAT2_patient = [];
projLAT2_simulated = [];
ptCloudSur = pointCloud([patient2xCm, patient2yCm, patient2zCm]);

for i = 1:length(patient2xSurfaceCm)

    [indiT,distsT] = findNearestNeighbors(ptCloudSur,[ patient2xSurfaceCm(i), patient2ySurfaceCm(i), patient2zSurfaceCm(i) ], 1);
    projLAT2_patient(i) = actTimes2(indiT);
    projLAT2_simulated(i) = actTimesSim2(indiT);

end

figure
scatter3(patient2xSurfaceCm, patient2ySurfaceCm, patient2zSurfaceCm, 190, projLAT2_patient, 'filled', 'MarkerEdgeColor','k');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (ms)';
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

figure
scatter3(patient2xSurfaceCm, patient2ySurfaceCm, patient2zSurfaceCm, 190, projLAT2_simulated, 'filled', 'MarkerEdgeColor','k');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (ms)';
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);


numBins = 9;
binEdges2 = linspace( min(min(projLAT2_simulated),min(projLAT2_patient)), max(max(projLAT2_simulated),max(projLAT2_patient)), numBins+1 );
binEdges2(end) = binEdges2(end) + 1;
binIndices_patient2 = {};
binIndices_simulated2 = {};

for i = 1:length(binEdges2)-1
    binIndices_patient2{i} = find( projLAT2_patient >= binEdges2(i) & projLAT2_patient < binEdges2(i+1)  );
    binIndices_simulated2{i} = find( projLAT2_simulated >= binEdges2(i) & projLAT2_simulated < binEdges2(i+1)  );
end


% 
% colorVec = { sscanf('0000ff','%2x%2x%2x',[1 3])/255, sscanf('3aa0ff','%2x%2x%2x',[1 3])/255, sscanf('5df2ff','%2x%2x%2x',[1 3])/255,...
%              sscanf('36aa74','%2x%2x%2x',[1 3])/255, sscanf('00ff00','%2x%2x%2x',[1 3])/255, sscanf('fff121','%2x%2x%2x',[1 3])/255,...
%              sscanf('ffaa00','%2x%2x%2x',[1 3])/255, sscanf('ff5500','%2x%2x%2x',[1 3])/255, sscanf('ff0000','%2x%2x%2x',[1 3])/255 };
% 

figure
plot3(patient2xSurfaceCm(binIndices_patient2{1}), patient2ySurfaceCm(binIndices_patient2{1}), patient2zSurfaceCm(binIndices_patient2{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(patient2xSurfaceCm(binIndices_patient2{2}), patient2ySurfaceCm(binIndices_patient2{2}), patient2zSurfaceCm(binIndices_patient2{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(patient2xSurfaceCm(binIndices_patient2{3}), patient2ySurfaceCm(binIndices_patient2{3}), patient2zSurfaceCm(binIndices_patient2{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(patient2xSurfaceCm(binIndices_patient2{4}), patient2ySurfaceCm(binIndices_patient2{4}), patient2zSurfaceCm(binIndices_patient2{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(patient2xSurfaceCm(binIndices_patient2{5}), patient2ySurfaceCm(binIndices_patient2{5}), patient2zSurfaceCm(binIndices_patient2{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(patient2xSurfaceCm(binIndices_patient2{6}), patient2ySurfaceCm(binIndices_patient2{6}), patient2zSurfaceCm(binIndices_patient2{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(patient2xSurfaceCm(binIndices_patient2{7}), patient2ySurfaceCm(binIndices_patient2{7}), patient2zSurfaceCm(binIndices_patient2{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(patient2xSurfaceCm(binIndices_patient2{8}), patient2ySurfaceCm(binIndices_patient2{8}), patient2zSurfaceCm(binIndices_patient2{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(patient2xSurfaceCm(binIndices_patient2{9}), patient2ySurfaceCm(binIndices_patient2{9}), patient2zSurfaceCm(binIndices_patient2{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);


figure
plot3(patient2xSurfaceCm(binIndices_simulated2{1}), patient2ySurfaceCm(binIndices_simulated2{1}), patient2zSurfaceCm(binIndices_simulated2{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{1}, 'MarkerEdgeColor', colorVec{1} );
hold on
plot3(patient2xSurfaceCm(binIndices_simulated2{2}), patient2ySurfaceCm(binIndices_simulated2{2}), patient2zSurfaceCm(binIndices_simulated2{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{2}, 'MarkerEdgeColor', colorVec{2} );
plot3(patient2xSurfaceCm(binIndices_simulated2{3}), patient2ySurfaceCm(binIndices_simulated2{3}), patient2zSurfaceCm(binIndices_simulated2{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{3}, 'MarkerEdgeColor', colorVec{3} );
plot3(patient2xSurfaceCm(binIndices_simulated2{4}), patient2ySurfaceCm(binIndices_simulated2{4}), patient2zSurfaceCm(binIndices_simulated2{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{4}, 'MarkerEdgeColor', colorVec{4} );
plot3(patient2xSurfaceCm(binIndices_simulated2{5}), patient2ySurfaceCm(binIndices_simulated2{5}), patient2zSurfaceCm(binIndices_simulated2{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{5}, 'MarkerEdgeColor', colorVec{5} );
plot3(patient2xSurfaceCm(binIndices_simulated2{6}), patient2ySurfaceCm(binIndices_simulated2{6}), patient2zSurfaceCm(binIndices_simulated2{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{6}, 'MarkerEdgeColor', colorVec{6} );
plot3(patient2xSurfaceCm(binIndices_simulated2{7}), patient2ySurfaceCm(binIndices_simulated2{7}), patient2zSurfaceCm(binIndices_simulated2{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{7}, 'MarkerEdgeColor', colorVec{7} );
plot3(patient2xSurfaceCm(binIndices_simulated2{8}), patient2ySurfaceCm(binIndices_simulated2{8}), patient2zSurfaceCm(binIndices_simulated2{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
plot3(patient2xSurfaceCm(binIndices_simulated2{9}), patient2ySurfaceCm(binIndices_simulated2{9}), patient2zSurfaceCm(binIndices_simulated2{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);












projLAT3_patient = [];
projLAT3_simulated = [];
ptCloudSur = pointCloud([patient3xCm, patient3yCm, patient3zCm]);

for i = 1:length(patient3xSurfaceCm)

    [indiT,distsT] = findNearestNeighbors(ptCloudSur,[ patient3xSurfaceCm(i), patient3ySurfaceCm(i), patient3zSurfaceCm(i) ], 1);
    projLAT3_patient(i) = actTimes3(indiT);
    projLAT3_simulated(i) = actTimesSim3(indiT);

end

figure
scatter3(patient3xSurfaceCm, patient3ySurfaceCm, patient3zSurfaceCm, 190, projLAT3_patient, 'filled', 'MarkerEdgeColor','k');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (ms)';
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

figure
scatter3(patient3xSurfaceCm, patient3ySurfaceCm, patient3zSurfaceCm, 190, projLAT3_simulated, 'filled', 'MarkerEdgeColor','k');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LATs (ms)';
caxis([40, 130]);
% title('H1A Patient');
set(gca, 'FontSize', 35);



numBins = 7;
binEdges3 = linspace( min(min(projLAT3_simulated),min(projLAT3_patient)), max(max(projLAT3_simulated),max(projLAT3_patient)), numBins+1 );
binEdges3(end) = binEdges3(end) + 1;
binIndices_patient3 = {};
binIndices_simulated3 = {};

for i = 1:length(binEdges3)-1
    binIndices_patient3{i} = find( projLAT3_patient >= binEdges3(i) & projLAT3_patient < binEdges3(i+1)  );
    binIndices_simulated3{i} = find( projLAT3_simulated >= binEdges3(i) & projLAT3_simulated < binEdges3(i+1)  );
end



% colorVec = { sscanf('0000ff','%2x%2x%2x',[1 3])/255, sscanf('3aa0ff','%2x%2x%2x',[1 3])/255,...
%              sscanf('00ff00','%2x%2x%2x',[1 3])/255, sscanf('fff121','%2x%2x%2x',[1 3])/255,...
%              sscanf('ffaa00','%2x%2x%2x',[1 3])/255, sscanf('ff5500','%2x%2x%2x',[1 3])/255, sscanf('ff0000','%2x%2x%2x',[1 3])/255 };


figure
plot3(patient3xSurfaceCm(binIndices_patient3{1}), patient3ySurfaceCm(binIndices_patient3{1}), patient3zSurfaceCm(binIndices_patient3{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec2{1}, 'MarkerEdgeColor', colorVec2{1} );
hold on
plot3(patient3xSurfaceCm(binIndices_patient3{2}), patient3ySurfaceCm(binIndices_patient3{2}), patient3zSurfaceCm(binIndices_patient3{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec2{2}, 'MarkerEdgeColor', colorVec2{2} );
plot3(patient3xSurfaceCm(binIndices_patient3{3}), patient3ySurfaceCm(binIndices_patient3{3}), patient3zSurfaceCm(binIndices_patient3{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec2{3}, 'MarkerEdgeColor', colorVec2{3} );
plot3(patient3xSurfaceCm(binIndices_patient3{4}), patient3ySurfaceCm(binIndices_patient3{4}), patient3zSurfaceCm(binIndices_patient3{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec2{4}, 'MarkerEdgeColor', colorVec2{4} );
plot3(patient3xSurfaceCm(binIndices_patient3{5}), patient3ySurfaceCm(binIndices_patient3{5}), patient3zSurfaceCm(binIndices_patient3{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec2{5}, 'MarkerEdgeColor', colorVec2{5} );
plot3(patient3xSurfaceCm(binIndices_patient3{6}), patient3ySurfaceCm(binIndices_patient3{6}), patient3zSurfaceCm(binIndices_patient3{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec2{6}, 'MarkerEdgeColor', colorVec2{6} );
plot3(patient3xSurfaceCm(binIndices_patient3{7}), patient3ySurfaceCm(binIndices_patient3{7}), patient3zSurfaceCm(binIndices_patient3{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec2{7}, 'MarkerEdgeColor', colorVec2{7} );
% plot3(patient3xSurfaceCm(binIndices_patient3{8}), patient3ySurfaceCm(binIndices_patient3{8}), patient3zSurfaceCm(binIndices_patient3{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
% plot3(patient3xSurfaceCm(binIndices_patient3{9}), patient3ySurfaceCm(binIndices_patient3{9}), patient3zSurfaceCm(binIndices_patient3{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);


figure
plot3(patient3xSurfaceCm(binIndices_simulated3{1}), patient3ySurfaceCm(binIndices_simulated3{1}), patient3zSurfaceCm(binIndices_simulated3{1}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec2{1}, 'MarkerEdgeColor', colorVec2{1} );
hold on
plot3(patient3xSurfaceCm(binIndices_simulated3{2}), patient3ySurfaceCm(binIndices_simulated3{2}), patient3zSurfaceCm(binIndices_simulated3{2}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec2{2}, 'MarkerEdgeColor', colorVec2{2} );
plot3(patient3xSurfaceCm(binIndices_simulated3{3}), patient3ySurfaceCm(binIndices_simulated3{3}), patient3zSurfaceCm(binIndices_simulated3{3}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec2{3}, 'MarkerEdgeColor', colorVec2{3} );
plot3(patient3xSurfaceCm(binIndices_simulated3{4}), patient3ySurfaceCm(binIndices_simulated3{4}), patient3zSurfaceCm(binIndices_simulated3{4}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec2{4}, 'MarkerEdgeColor', colorVec2{4} );
plot3(patient3xSurfaceCm(binIndices_simulated3{5}), patient3ySurfaceCm(binIndices_simulated3{5}), patient3zSurfaceCm(binIndices_simulated3{5}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec2{5}, 'MarkerEdgeColor', colorVec2{5} );
plot3(patient3xSurfaceCm(binIndices_simulated3{6}), patient3ySurfaceCm(binIndices_simulated3{6}), patient3zSurfaceCm(binIndices_simulated3{6}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec2{6}, 'MarkerEdgeColor', colorVec2{6} );
plot3(patient3xSurfaceCm(binIndices_simulated3{7}), patient3ySurfaceCm(binIndices_simulated3{7}), patient3zSurfaceCm(binIndices_simulated3{7}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec2{7}, 'MarkerEdgeColor', colorVec2{7} );
% plot3(patient3xSurfaceCm(binIndices_simulated3{8}), patient3ySurfaceCm(binIndices_simulated3{8}), patient3zSurfaceCm(binIndices_simulated3{8}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{8}, 'MarkerEdgeColor', colorVec{8} );
% plot3(patient3xSurfaceCm(binIndices_simulated3{9}), patient3ySurfaceCm(binIndices_simulated3{9}), patient3zSurfaceCm(binIndices_simulated3{9}), 'o', 'MarkerSize', 25, 'MarkerFaceColor', colorVec{9}, 'MarkerEdgeColor', colorVec{9} );
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% title('H1A Patient');
set(gca, 'FontSize', 35);










% 
% indP2 = [find( abs(patient2xCm - electrodesP2(2,1)) < .01 & abs(patient2yCm - electrodesP2(2,2)) < .01 & abs(patient2zCm - electrodesP2(2,3)) < .01 ),...
%          find( abs(patient2xCm - electrodesP2(3,1)) < .01 & abs(patient2yCm - electrodesP2(3,2)) < .01 & abs(patient2zCm - electrodesP2(3,3)) < .01 ),...
%          find( abs(patient2xCm - electrodesP2(4,1)) < .01 & abs(patient2yCm - electrodesP2(4,2)) < .01 & abs(patient2zCm - electrodesP2(4,3)) < .01 ),...
%          find( abs(patient2xCm - electrodesP2(5,1)) < .01 & abs(patient2yCm - electrodesP2(5,2)) < .01 & abs(patient2zCm - electrodesP2(5,3)) < .01 )];
% 
% indP3 = [find( abs(patient3xCm - electrodesP3(2,1)) < .01 & abs(patient3yCm - electrodesP3(2,2)) < .01 & abs(patient3zCm - electrodesP3(2,3)) < .01 ),...
%          find( abs(patient3xCm - electrodesP3(3,1)) < .01 & abs(patient3yCm - electrodesP3(3,2)) < .01 & abs(patient3zCm - electrodesP3(3,3)) < .01 ),...
%          find( abs(patient3xCm - electrodesP3(4,1)) < .01 & abs(patient3yCm - electrodesP3(4,2)) < .01 & abs(patient3zCm - electrodesP3(4,3)) < .01 ),...
%          find( abs(patient3xCm - electrodesP3(5,1)) < .01 & abs(patient3yCm - electrodesP3(5,2)) < .01 & abs(patient3zCm - electrodesP3(5,3)) < .01 )];
% 


indP2 = [177, 470, 191, 504];
indP3 = [175, 674, 193, 698, 694, 706];




figure
plot( signalsP2(:,indP2(1)), 'DisplayName', 'Patient', 'LineWidth', 4 )
hold on
plot( simulatedParoxysmal2_downsampled(:,indP2(1)), 'DisplayName', 'Simulation', 'LineWidth', 4 )
legend()
xlabel('time (ms)');
ylabel('Ve (mV)');
% title('H1A Patient');
set(gca, 'FontSize', 35);

figure
plot( signalsP2(:,indP2(2)), 'DisplayName', 'Patient', 'LineWidth', 4 )
hold on
plot( simulatedParoxysmal2_downsampled(:,indP2(2)), 'DisplayName', 'Simulation', 'LineWidth', 4 )
legend()
xlabel('time (ms)');
ylabel('Ve (mV)');
% title('H1A Patient');
set(gca, 'FontSize', 35);

figure
plot( signalsP2(:,indP2(3)), 'DisplayName', 'Patient', 'LineWidth', 4 )
hold on
plot( simulatedParoxysmal2_downsampled(:,indP2(3)), 'DisplayName', 'Simulation', 'LineWidth', 4 )
legend()
xlabel('time (ms)');
ylabel('Ve (mV)');
% title('H1A Patient');
set(gca, 'FontSize', 35);

figure
plot( signalsP2(:,indP2(4)), 'DisplayName', 'Patient', 'LineWidth', 4 )
hold on
plot( simulatedParoxysmal2_downsampled(:,indP2(4)), 'DisplayName', 'Simulation', 'LineWidth', 4 )
legend()
xlabel('time (ms)');
ylabel('Ve (mV)');
% title('H1A Patient');
set(gca, 'FontSize', 35);




figure
plot( signalsP3(:,indP3(1)), 'DisplayName', 'Patient', 'LineWidth', 4 )
hold on
plot( simulatedParoxysmal3_downsampled(:,indP3(1)), 'DisplayName', 'Simulation', 'LineWidth', 4 )
legend()
xlabel('time (ms)');
ylabel('Ve (mV)');
% title('H1A Patient');
set(gca, 'FontSize', 35);

figure
plot( signalsP3(:,indP3(2)), 'DisplayName', 'Patient', 'LineWidth', 4 )
hold on
plot( simulatedParoxysmal3_downsampled(:,indP3(2)), 'DisplayName', 'Simulation', 'LineWidth', 4 )
legend()
xlabel('time (ms)');
ylabel('Ve (mV)');
% title('H1A Patient');
set(gca, 'FontSize', 35);

figure
plot( signalsP3(:,indP3(3)), 'DisplayName', 'Patient', 'LineWidth', 4 )
hold on
plot( simulatedParoxysmal3_downsampled(:,indP3(3)), 'DisplayName', 'Simulation', 'LineWidth', 4 )
legend()
xlabel('time (ms)');
ylabel('Ve (mV)');
% title('H1A Patient');
set(gca, 'FontSize', 35);

figure
plot( signalsP3(:,indP3(4)), 'DisplayName', 'Patient', 'LineWidth', 4 )
hold on
plot( simulatedParoxysmal3_downsampled(:,indP3(4)), 'DisplayName', 'Simulation', 'LineWidth', 4 )
legend()
xlabel('time (ms)');
ylabel('Ve (mV)');
% title('H1A Patient');
set(gca, 'FontSize', 35);







figure
plot( signalsP3(:,indP3(5)), 'DisplayName', 'Patient', 'LineWidth', 4 )
hold on
plot( simulatedParoxysmal3_downsampled(:,indP3(5)), 'DisplayName', 'Simulation', 'LineWidth', 4 )
legend()
xlabel('time (ms)');
ylabel('Ve (mV)');
% title('H1A Patient');
set(gca, 'FontSize', 35);

figure
plot( signalsP3(:,indP3(6)), 'DisplayName', 'Patient', 'LineWidth', 4 )
hold on
plot( simulatedParoxysmal3_downsampled(:,indP3(6)), 'DisplayName', 'Simulation', 'LineWidth', 4 )
legend()
xlabel('time (ms)');
ylabel('Ve (mV)');
% title('H1A Patient');
set(gca, 'FontSize', 35);












figure
pause(4)
for i = 650:2:length(signalsP3(1,:))
    subplot(1,2,1)
    plot( signalsP3(:,i), 'DisplayName', 'Patient', 'LineWidth', 4 )
    hold on
    plot( simulatedParoxysmal3_downsampled(:,i), 'DisplayName', 'Simulation', 'LineWidth', 4 )
    hold off
    title(i)
    subplot(1,2,2)
    plot3(patient3xSurfaceCm, patient3ySurfaceCm, patient3zSurfaceCm, 'ro')
    hold on
    scatter3(patient3xCm(i), patient3yCm(i), patient3zCm(i), 100, 'bo', 'filled')
    hold off
    pause(3)
end





%% STUFF FOR TABLES: WAVELET DECOMPOSITION OR VAR + PCA

%simulatedParoxysmal3_downsampled
timePointsToTake = 140;
levelDec = 4;
vecLevel = 1:levelDec;


patientNum = 2;
egmsParoxysmal2=[];
%fullEgms -> point in space, which electrode, time points
for i=1:length(patient2xCm)
    elec_1=[squeeze(expDB.dataSet{patientNum}.map.continuousMapping.fullEgms(i,1,1:timePointsToTake))];
    egmsParoxysmal2 = [egmsParoxysmal2,elec_1];
end


patientNum = 3;
egmsParoxysmal3=[];
%fullEgms -> point in space, which electrode, time points
for i=1:length(patient3xCm)
    elec_1=[squeeze(expDB.dataSet{patientNum}.map.continuousMapping.fullEgms(i,1,1:timePointsToTake))];
    egmsParoxysmal3 = [egmsParoxysmal3,elec_1];
end


cdAll_Paroxysmal2 = [];
caAll_Paroxysmal2 = [];
coefAll_Paroxysmal2 = [];
cdAll_Paroxysmal2_Sim = [];
caAll_Paroxysmal2_Sim = [];
coefAll_Paroxysmal2_Sim = [];

cdAll_Paroxysmal3 = [];
caAll_Paroxysmal3 = [];
coefAll_Paroxysmal3 = [];
cdAll_Paroxysmal3_Sim = [];
caAll_Paroxysmal3_Sim = [];
coefAll_Paroxysmal3_Sim = [];






EGMDuration_Paroxysmal2 = [];
deflectionCount_Paroxysmal2 = [];
signals_Paroxysmal2 = [];
downstroke_Paroxysmal2 = [];
upstrokeRise_Paroxysmal2 = [];
upstrokeBack_Paroxysmal2 = [];
signalWidth_Paroxysmal2 = [];
baseline_Paroxysmal2 = [];
PPA_Paroxysmal2 = [];
AlldistanceToProj_Paroxysmal2 = [];

EGMDuration_Paroxysmal2_Sim = [];
deflectionCount_Paroxysmal2_Sim = [];
signals_Paroxysmal2_Sim = [];
downstroke_Paroxysmal2_Sim = [];
upstrokeRise_Paroxysmal2_Sim = [];
upstrokeBack_Paroxysmal2_Sim = [];
signalWidth_Paroxysmal2_Sim = [];
baseline_Paroxysmal2_Sim = [];
PPA_Paroxysmal2_Sim = [];
AlldistanceToProj_Paroxysmal2_Sim = [];

ptCloud_Paroxysmal2 = pointCloud([patient2xSurfaceCm,patient2ySurfaceCm,patient2zSurfaceCm]);
for i = 1:length(egmsParoxysmal2)
    
    [EGMduration, SD_curve, aboveValues, DeflectionCount, baselineMeanSDSW, bumpinessFactor, xtoCheck, ytoCheck] = SD_SlidingWindowF([1:1:length(egmsParoxysmal2(:,i))], egmsParoxysmal2(:,i), 1, 40, 15, 70, 3, 100);
    locSignal =  ( egmsParoxysmal2(:,i) - min(egmsParoxysmal2(:,i)) ) ./ ( max(egmsParoxysmal2(:,i)) - min(egmsParoxysmal2(:,i)) );
    locSignalCut =  ( ytoCheck - min(ytoCheck) ) ./ ( max(ytoCheck) - min(ytoCheck) );
    locSignalCutOG =  ytoCheck;
    [dVeMono_dt, downstrokeMono, midSignalIndMono, signalDiffMonoPointsTime, upstrokeRiseMono, upstrokeBackMono, riseSignalIndMono, backSignalIndMono, PositiveDeflectionDuration_Mono, NegativeDeflectionDuration_Mono, baselineMono, signalWidth, timePeak, timeTrough] = EGM_Characterizing_Function(locSignal, 1:timePointsToTake, 1, 1, EGMduration);
    localPPA = max(ytoCheck) - min(ytoCheck);

    [indices,distVal] = findNearestNeighbors(ptCloud_Paroxysmal2,[patient2xCm(i), patient2yCm(i), patient2zCm(i)],1);

    EGMDuration_Paroxysmal2 = [EGMDuration_Paroxysmal2; EGMduration];
    deflectionCount_Paroxysmal2 = [deflectionCount_Paroxysmal2; DeflectionCount];
    signals_Paroxysmal2 = [ signals_Paroxysmal2; locSignalCutOG' ];
    downstroke_Paroxysmal2 = [ downstroke_Paroxysmal2; downstrokeMono ];
    upstrokeRise_Paroxysmal2 = [ upstrokeRise_Paroxysmal2; upstrokeRiseMono ];
    upstrokeBack_Paroxysmal2 = [ upstrokeBack_Paroxysmal2; upstrokeBackMono ];
    signalWidth_Paroxysmal2 = [ signalWidth_Paroxysmal2; signalWidth ];
    baseline_Paroxysmal2 = [ baseline_Paroxysmal2; baselineMono ];
    PPA_Paroxysmal2 = [ PPA_Paroxysmal2; localPPA ];
    AlldistanceToProj_Paroxysmal2 = [ AlldistanceToProj_Paroxysmal2; distVal ];

    [c,l] = wavedec(locSignalCut,levelDec,"rbio4.4");
    approx = appcoef(c,l,"rbio4.4");
    [cd1,cd2,cd3,cd4] = detcoef(c,l,vecLevel);
    cdAll_Paroxysmal2 = [ cdAll_Paroxysmal2, cd4 ];
    caAll_Paroxysmal2 = [ caAll_Paroxysmal2, approx ];
    coefAll_Paroxysmal2 = [ coefAll_Paroxysmal2, [cd4; approx] ];



    [EGMduration, SD_curve, aboveValues, DeflectionCount, baselineMeanSDSW, bumpinessFactor, xtoCheck, ytoCheck] = SD_SlidingWindowF([1:1:length(simulatedParoxysmal2_downsampled(1:timePointsToTake,i))], simulatedParoxysmal2_downsampled(1:timePointsToTake,i), 1, 40, 15, 70, 3, 100);
    locSignal =  ( simulatedParoxysmal2_downsampled(1:140,i) - min(simulatedParoxysmal2_downsampled(1:140,i)) ) ./ ( max(simulatedParoxysmal2_downsampled(1:140,i)) - min(simulatedParoxysmal2_downsampled(1:140,i)) );
    locSignalCut =  ( ytoCheck - min(ytoCheck) ) ./ ( max(ytoCheck) - min(ytoCheck) );
    locSignalCutOG =  ytoCheck;
    [dVeMono_dt, downstrokeMono, midSignalIndMono, signalDiffMonoPointsTime, upstrokeRiseMono, upstrokeBackMono, riseSignalIndMono, backSignalIndMono, PositiveDeflectionDuration_Mono, NegativeDeflectionDuration_Mono, baselineMono, signalWidth, timePeak, timeTrough] = EGM_Characterizing_Function(locSignal, 1:timePointsToTake, 1, 1, EGMduration);
    localPPA = max(ytoCheck) - min(ytoCheck);

    EGMDuration_Paroxysmal2_Sim = [EGMDuration_Paroxysmal2_Sim; EGMduration];
    deflectionCount_Paroxysmal2_Sim = [deflectionCount_Paroxysmal2_Sim; DeflectionCount];
    signals_Paroxysmal2_Sim = [ signals_Paroxysmal2_Sim; locSignalCutOG' ];
    downstroke_Paroxysmal2_Sim = [ downstroke_Paroxysmal2_Sim; downstrokeMono ];
    upstrokeRise_Paroxysmal2_Sim = [ upstrokeRise_Paroxysmal2_Sim; upstrokeRiseMono ];
    upstrokeBack_Paroxysmal2_Sim = [ upstrokeBack_Paroxysmal2_Sim; upstrokeBackMono ];
    signalWidth_Paroxysmal2_Sim = [ signalWidth_Paroxysmal2_Sim; signalWidth ];
    baseline_Paroxysmal2_Sim = [ baseline_Paroxysmal2_Sim; baselineMono ];
    PPA_Paroxysmal2_Sim = [ PPA_Paroxysmal2_Sim; localPPA ];
    AlldistanceToProj_Paroxysmal2_Sim = [ AlldistanceToProj_Paroxysmal2_Sim; distVal ];

    [c,l] = wavedec(locSignalCut,levelDec,"rbio4.4");
    approx = appcoef(c,l,"rbio4.4");
    [cd1,cd2,cd3,cd4] = detcoef(c,l,vecLevel);
    cdAll_Paroxysmal2_Sim = [ cdAll_Paroxysmal2_Sim, cd4 ];
    caAll_Paroxysmal2_Sim = [ caAll_Paroxysmal2_Sim, approx ];
    coefAll_Paroxysmal2_Sim = [ coefAll_Paroxysmal2_Sim, [cd4; approx] ];

end









EGMDuration_Paroxysmal3 = [];
deflectionCount_Paroxysmal3 = [];
signals_Paroxysmal3 = [];
downstroke_Paroxysmal3 = [];
upstrokeRise_Paroxysmal3 = [];
upstrokeBack_Paroxysmal3 = [];
signalWidth_Paroxysmal3 = [];
baseline_Paroxysmal3 = [];
PPA_Paroxysmal3 = [];
AlldistanceToProj_Paroxysmal3 = [];

EGMDuration_Paroxysmal3_Sim = [];
deflectionCount_Paroxysmal3_Sim = [];
signals_Paroxysmal3_Sim = [];
downstroke_Paroxysmal3_Sim = [];
upstrokeRise_Paroxysmal3_Sim = [];
upstrokeBack_Paroxysmal3_Sim = [];
signalWidth_Paroxysmal3_Sim = [];
baseline_Paroxysmal3_Sim = [];
PPA_Paroxysmal3_Sim = [];
AlldistanceToProj_Paroxysmal3_Sim = [];

ptCloud_Paroxysmal3 = pointCloud([patient3xSurfaceCm,patient3ySurfaceCm,patient3zSurfaceCm]);
for i = 1:length(egmsParoxysmal3)
    
    [EGMduration, SD_curve, aboveValues, DeflectionCount, baselineMeanSDSW, bumpinessFactor, xtoCheck, ytoCheck] = SD_SlidingWindowF([1:1:length(egmsParoxysmal3(:,i))], egmsParoxysmal3(:,i), 1, 40, 15, 70, 3, 100);
    locSignal =  ( egmsParoxysmal3(:,i) - min(egmsParoxysmal3(:,i)) ) ./ ( max(egmsParoxysmal3(:,i)) - min(egmsParoxysmal3(:,i)) );
    locSignalCut =  ( ytoCheck - min(ytoCheck) ) ./ ( max(ytoCheck) - min(ytoCheck) );
    locSignalCutOG =  ytoCheck;
    [dVeMono_dt, downstrokeMono, midSignalIndMono, signalDiffMonoPointsTime, upstrokeRiseMono, upstrokeBackMono, riseSignalIndMono, backSignalIndMono, PositiveDeflectionDuration_Mono, NegativeDeflectionDuration_Mono, baselineMono, signalWidth, timePeak, timeTrough] = EGM_Characterizing_Function(locSignal, 1:timePointsToTake, 1, 1, EGMduration);
    localPPA = max(ytoCheck) - min(ytoCheck);

    [indices,distVal] = findNearestNeighbors(ptCloud_Paroxysmal3,[patient3xCm(i), patient3yCm(i), patient3zCm(i)],1);

    EGMDuration_Paroxysmal3 = [EGMDuration_Paroxysmal3; EGMduration];
    deflectionCount_Paroxysmal3 = [deflectionCount_Paroxysmal3; DeflectionCount];
    signals_Paroxysmal3 = [ signals_Paroxysmal3; locSignalCutOG' ];
    downstroke_Paroxysmal3 = [ downstroke_Paroxysmal3; downstrokeMono ];
    upstrokeRise_Paroxysmal3 = [ upstrokeRise_Paroxysmal3; upstrokeRiseMono ];
    upstrokeBack_Paroxysmal3 = [ upstrokeBack_Paroxysmal3; upstrokeBackMono ];
    signalWidth_Paroxysmal3 = [ signalWidth_Paroxysmal3; signalWidth ];
    baseline_Paroxysmal3 = [ baseline_Paroxysmal3; baselineMono ];
    PPA_Paroxysmal3 = [ PPA_Paroxysmal3; localPPA ];
    AlldistanceToProj_Paroxysmal3 = [ AlldistanceToProj_Paroxysmal3; distVal ];

    [c,l] = wavedec(locSignalCut,levelDec,"rbio4.4");
    approx = appcoef(c,l,"rbio4.4");
    [cd1,cd2,cd3,cd4] = detcoef(c,l,vecLevel);
    cdAll_Paroxysmal3 = [ cdAll_Paroxysmal3, cd4 ];
    caAll_Paroxysmal3 = [ caAll_Paroxysmal3, approx ];
    coefAll_Paroxysmal3 = [ coefAll_Paroxysmal3, [cd4; approx] ];





    [EGMduration, SD_curve, aboveValues, DeflectionCount, baselineMeanSDSW, bumpinessFactor, xtoCheck, ytoCheck] = SD_SlidingWindowF([1:1:length(simulatedParoxysmal3_downsampled(1:timePointsToTake,i))], simulatedParoxysmal3_downsampled(1:timePointsToTake,i), 1, 40, 15, 70, 3, 100);
    locSignal =  ( simulatedParoxysmal3_downsampled(1:140,i) - min(simulatedParoxysmal3_downsampled(1:140,i)) ) ./ ( max(simulatedParoxysmal3_downsampled(1:140,i)) - min(simulatedParoxysmal3_downsampled(1:140,i)) );
    locSignalCut =  ( ytoCheck - min(ytoCheck) ) ./ ( max(ytoCheck) - min(ytoCheck) );
    locSignalCutOG =  ytoCheck;
    [dVeMono_dt, downstrokeMono, midSignalIndMono, signalDiffMonoPointsTime, upstrokeRiseMono, upstrokeBackMono, riseSignalIndMono, backSignalIndMono, PositiveDeflectionDuration_Mono, NegativeDeflectionDuration_Mono, baselineMono, signalWidth, timePeak, timeTrough] = EGM_Characterizing_Function(locSignal, 1:timePointsToTake, 1, 1, EGMduration);
    localPPA = max(ytoCheck) - min(ytoCheck);

    EGMDuration_Paroxysmal3_Sim = [EGMDuration_Paroxysmal3_Sim; EGMduration];
    deflectionCount_Paroxysmal3_Sim = [deflectionCount_Paroxysmal3_Sim; DeflectionCount];
    signals_Paroxysmal3_Sim = [ signals_Paroxysmal3_Sim; locSignalCutOG' ];
    downstroke_Paroxysmal3_Sim = [ downstroke_Paroxysmal3_Sim; downstrokeMono ];
    upstrokeRise_Paroxysmal3_Sim = [ upstrokeRise_Paroxysmal3_Sim; upstrokeRiseMono ];
    upstrokeBack_Paroxysmal3_Sim = [ upstrokeBack_Paroxysmal3_Sim; upstrokeBackMono ];
    signalWidth_Paroxysmal3_Sim = [ signalWidth_Paroxysmal3_Sim; signalWidth ];
    baseline_Paroxysmal3_Sim = [ baseline_Paroxysmal3_Sim; baselineMono ];
    PPA_Paroxysmal3_Sim = [ PPA_Paroxysmal3_Sim; localPPA ];
    AlldistanceToProj_Paroxysmal3_Sim = [ AlldistanceToProj_Paroxysmal3_Sim; distVal ];

    [c,l] = wavedec(locSignalCut,levelDec,"rbio4.4");
    approx = appcoef(c,l,"rbio4.4");
    [cd1,cd2,cd3,cd4] = detcoef(c,l,vecLevel);
    cdAll_Paroxysmal3_Sim = [ cdAll_Paroxysmal3_Sim, cd4 ];
    caAll_Paroxysmal3_Sim = [ caAll_Paroxysmal3_Sim, approx ];
    coefAll_Paroxysmal3_Sim = [ coefAll_Paroxysmal3_Sim, [cd4; approx] ];


end






[coeffWD_Paroxysmal2,scoreWD_Paroxysmal2,latentWD_Paroxysmal2,tsquaredWD_Paroxysmal2,explainedWD_Paroxysmal2,muWD_Paroxysmal2] = pca(cdAll_Paroxysmal2');
[coeffWD_Paroxysmal3,scoreWD_Paroxysmal3,latentWD_Paroxysmal3,tsquaredWD_Paroxysmal3,explainedWD_Paroxysmal3,muWD_Paroxysmal3] = pca(cdAll_Paroxysmal3');

% coeffWD
% explainedWD
sum(explainedWD_Paroxysmal2(1:5))
sum(explainedWD_Paroxysmal3(1:5))

evaluationWD_2 = evalclusters(scoreWD_Paroxysmal2(:,1:5),"kmeans","CalinskiHarabasz","KList",1:20)
figure; plot(evaluationWD_2)
evaluationWD_3 = evalclusters(scoreWD_Paroxysmal3(:,1:5),"kmeans","CalinskiHarabasz","KList",1:20)
figure; plot(evaluationWD_3)


PCA1vectWD_Paroxysmal2 = [];
PCA2vectWD_Paroxysmal2 = [];
PCA3vectWD_Paroxysmal2 = [];
PCA4vectWD_Paroxysmal2 = [];
PCA5vectWD_Paroxysmal2 = [];

PCA1vectWD_Paroxysmal2_Sim = [];
PCA2vectWD_Paroxysmal2_Sim = [];
PCA3vectWD_Paroxysmal2_Sim = [];
PCA4vectWD_Paroxysmal2_Sim = [];
PCA5vectWD_Paroxysmal2_Sim = [];

for i = 1:length(cdAll_Paroxysmal2(1,:))
    PCA1vectWD_Paroxysmal2(i) = coeffWD_Paroxysmal2(:,1)'*cdAll_Paroxysmal2(:,i);
    PCA2vectWD_Paroxysmal2(i) = coeffWD_Paroxysmal2(:,2)'*cdAll_Paroxysmal2(:,i);
    PCA3vectWD_Paroxysmal2(i) = coeffWD_Paroxysmal2(:,3)'*cdAll_Paroxysmal2(:,i);
    PCA4vectWD_Paroxysmal2(i) = coeffWD_Paroxysmal2(:,4)'*cdAll_Paroxysmal2(:,i);
    PCA5vectWD_Paroxysmal2(i) = coeffWD_Paroxysmal2(:,5)'*cdAll_Paroxysmal2(:,i);

    PCA1vectWD_Paroxysmal2_Sim(i) = coeffWD_Paroxysmal2(:,1)'*cdAll_Paroxysmal2_Sim(:,i);
    PCA2vectWD_Paroxysmal2_Sim(i) = coeffWD_Paroxysmal2(:,2)'*cdAll_Paroxysmal2_Sim(:,i);
    PCA3vectWD_Paroxysmal2_Sim(i) = coeffWD_Paroxysmal2(:,3)'*cdAll_Paroxysmal2_Sim(:,i);
    PCA4vectWD_Paroxysmal2_Sim(i) = coeffWD_Paroxysmal2(:,4)'*cdAll_Paroxysmal2_Sim(:,i);
    PCA5vectWD_Paroxysmal2_Sim(i) = coeffWD_Paroxysmal2(:,5)'*cdAll_Paroxysmal2_Sim(:,i);
end

opts = statset('Display','final');
[idxWD_Paroxysmal2,CWD_Paroxysmal2,sumdWD_Paroxysmal2,DWD_Paroxysmal2] = kmeans([PCA1vectWD_Paroxysmal2', PCA2vectWD_Paroxysmal2', PCA3vectWD_Paroxysmal2', PCA4vectWD_Paroxysmal2', PCA5vectWD_Paroxysmal2'],3,'Distance','cosine','EmptyAction','drop','MaxIter',1000,'Replicates',10,'Options',opts);

indexCluster1Paroxysmal2_WD = find(idxWD_Paroxysmal2 == 1);
indexCluster2Paroxysmal2_WD = find(idxWD_Paroxysmal2 == 2);
indexCluster3Paroxysmal2_WD = find(idxWD_Paroxysmal2 == 3);
indexCluster4Paroxysmal2_WD = find(idxWD_Paroxysmal2 == 4);
indexCluster5Paroxysmal2_WD = find(idxWD_Paroxysmal2 == 5);
indexCluster6Paroxysmal2_WD = find(idxWD_Paroxysmal2 == 6);
indexCluster7Paroxysmal2_WD = find(idxWD_Paroxysmal2 == 7);
% indexCluster8Paroxysmal_WD = find(idxWD_Paroxysmal == 8);




PCA1vectWD_Paroxysmal3 = [];
PCA2vectWD_Paroxysmal3 = [];
PCA3vectWD_Paroxysmal3 = [];
PCA4vectWD_Paroxysmal3 = [];
PCA5vectWD_Paroxysmal3 = [];

PCA1vectWD_Paroxysmal3_Sim = [];
PCA2vectWD_Paroxysmal3_Sim = [];
PCA3vectWD_Paroxysmal3_Sim = [];
PCA4vectWD_Paroxysmal3_Sim = [];
PCA5vectWD_Paroxysmal3_Sim = [];

for i = 1:length(cdAll_Paroxysmal3(1,:))
    PCA1vectWD_Paroxysmal3(i) = coeffWD_Paroxysmal3(:,1)'*cdAll_Paroxysmal3(:,i);
    PCA2vectWD_Paroxysmal3(i) = coeffWD_Paroxysmal3(:,2)'*cdAll_Paroxysmal3(:,i);
    PCA3vectWD_Paroxysmal3(i) = coeffWD_Paroxysmal3(:,3)'*cdAll_Paroxysmal3(:,i);
    PCA4vectWD_Paroxysmal3(i) = coeffWD_Paroxysmal3(:,4)'*cdAll_Paroxysmal3(:,i);
    PCA5vectWD_Paroxysmal3(i) = coeffWD_Paroxysmal3(:,5)'*cdAll_Paroxysmal3(:,i);

    PCA1vectWD_Paroxysmal3_Sim(i) = coeffWD_Paroxysmal3(:,1)'*cdAll_Paroxysmal3_Sim(:,i);
    PCA2vectWD_Paroxysmal3_Sim(i) = coeffWD_Paroxysmal3(:,2)'*cdAll_Paroxysmal3_Sim(:,i);
    PCA3vectWD_Paroxysmal3_Sim(i) = coeffWD_Paroxysmal3(:,3)'*cdAll_Paroxysmal3_Sim(:,i);
    PCA4vectWD_Paroxysmal3_Sim(i) = coeffWD_Paroxysmal3(:,4)'*cdAll_Paroxysmal3_Sim(:,i);
    PCA5vectWD_Paroxysmal3_Sim(i) = coeffWD_Paroxysmal3(:,5)'*cdAll_Paroxysmal3_Sim(:,i);
end

opts = statset('Display','final');
[idxWD_Paroxysmal3,CWD_Paroxysmal3,sumdWD_Paroxysmal3,DWD_Paroxysmal3] = kmeans([PCA1vectWD_Paroxysmal3', PCA2vectWD_Paroxysmal3', PCA3vectWD_Paroxysmal3', PCA4vectWD_Paroxysmal3', PCA5vectWD_Paroxysmal3'],3,'Distance','cosine','EmptyAction','drop','MaxIter',1000,'Replicates',10,'Options',opts);

indexCluster1Paroxysmal3_WD = find(idxWD_Paroxysmal3 == 1);
indexCluster2Paroxysmal3_WD = find(idxWD_Paroxysmal3 == 2);
indexCluster3Paroxysmal3_WD = find(idxWD_Paroxysmal3 == 3);
indexCluster4Paroxysmal3_WD = find(idxWD_Paroxysmal3 == 4);
indexCluster5Paroxysmal3_WD = find(idxWD_Paroxysmal3 == 5);
indexCluster6Paroxysmal3_WD = find(idxWD_Paroxysmal3 == 6);
indexCluster7Paroxysmal3_WD = find(idxWD_Paroxysmal3 == 7);
% indexCluster8Paroxysmal_WD = find(idxWD_Paroxysmal == 8);










dataToClassify_Paroxysmal2 = [ downstroke_Paroxysmal2, upstrokeRise_Paroxysmal2, upstrokeBack_Paroxysmal2, signalWidth_Paroxysmal2, baseline_Paroxysmal2, EGMDuration_Paroxysmal2, AlldistanceToProj_Paroxysmal2, PPA_Paroxysmal2, deflectionCount_Paroxysmal2];
% dataToClassify_Paroxysmal2 = [ downstroke_Paroxysmal2, upstrokeRise_Paroxysmal2, upstrokeBack_Paroxysmal2, signalWidth_Paroxysmal2, baseline_Paroxysmal2, EGMDuration_Paroxysmal2, AlldistanceToProj_Paroxysmal2, PPA_Paroxysmal2];
dataToClassify_Paroxysmal2_Scaled = dataToClassify_Paroxysmal2;

maxScalingData = [];
minScalingData = [];
% % % % SCALING TO MAKE ALL SCALES MATCH
for i = 1:length( dataToClassify_Paroxysmal2(1,:) )
    if max(dataToClassify_Paroxysmal2(:,i)) - min(dataToClassify_Paroxysmal2(:,i)) > 0
        dataToClassify_Paroxysmal2_Scaled(:,i) = ( dataToClassify_Paroxysmal2(:,i) - min(dataToClassify_Paroxysmal2(:,i)) )./( max(dataToClassify_Paroxysmal2(:,i)) - min(dataToClassify_Paroxysmal2(:,i)) );
        
        maxScalingData(i,2) = max(dataToClassify_Paroxysmal2(:,i));
        minScalingData(i,2) = min(dataToClassify_Paroxysmal2(:,i));
    end
end


maxScalingData = max(maxScalingData');
minScalingData = min(minScalingData');

[coeffVAR_Paroxysmal2,scoreVAR_Paroxysmal2,latentVAR_Paroxysmal2,tsquaredVAR_Paroxysmal2,explainedVAR_Paroxysmal2,muVAR_Paroxysmal2] = pca(dataToClassify_Paroxysmal2_Scaled);

sum( explainedVAR_Paroxysmal2(1:5) )

%More important variables are 8, 4, 6, 7

evaluationPCAVAR_Paroxysmal2 = evalclusters(scoreVAR_Paroxysmal2(:,1:5),"kmeans","CalinskiHarabasz","KList",1:20)
figure; plot(evaluationPCAVAR_Paroxysmal2)

dataToClassify_Paroxysmal2_Sim = [ downstroke_Paroxysmal2_Sim, upstrokeRise_Paroxysmal2_Sim, upstrokeBack_Paroxysmal2_Sim, signalWidth_Paroxysmal2_Sim, baseline_Paroxysmal2_Sim, EGMDuration_Paroxysmal2_Sim, AlldistanceToProj_Paroxysmal2_Sim, PPA_Paroxysmal2_Sim, deflectionCount_Paroxysmal2_Sim]; 
dataToClassify_Paroxysmal2_Scaled_Sim = dataToClassify_Paroxysmal2_Sim;
for i = 1:length( dataToClassify_Paroxysmal2_Sim(1,:) )
        dataToClassify_Paroxysmal2_Scaled_Sim(:,i) = ( dataToClassify_Paroxysmal2_Sim(:,i) - minScalingData(i) )./( maxScalingData(i) - minScalingData(i) );
end


PCA1vectVAR_Paroxysmal2 = [];
PCA2vectVAR_Paroxysmal2 = [];
PCA3vectVAR_Paroxysmal2 = [];
PCA4vectVAR_Paroxysmal2 = [];
PCA5vectVAR_Paroxysmal2 = [];

PCA1vectVAR_Paroxysmal2_Sim = [];
PCA2vectVAR_Paroxysmal2_Sim = [];
PCA3vectVAR_Paroxysmal2_Sim = [];
PCA4vectVAR_Paroxysmal2_Sim = [];
PCA5vectVAR_Paroxysmal2_Sim = [];

for i = 1: length(dataToClassify_Paroxysmal2_Scaled(:,1))
    PCA1vectVAR_Paroxysmal2(i) = coeffVAR_Paroxysmal2(:,1)'*dataToClassify_Paroxysmal2_Scaled(i,:)';
    PCA2vectVAR_Paroxysmal2(i) = coeffVAR_Paroxysmal2(:,2)'*dataToClassify_Paroxysmal2_Scaled(i,:)';
    PCA3vectVAR_Paroxysmal2(i) = coeffVAR_Paroxysmal2(:,3)'*dataToClassify_Paroxysmal2_Scaled(i,:)';
    PCA4vectVAR_Paroxysmal2(i) = coeffVAR_Paroxysmal2(:,4)'*dataToClassify_Paroxysmal2_Scaled(i,:)';
    PCA5vectVAR_Paroxysmal2(i) = coeffVAR_Paroxysmal2(:,5)'*dataToClassify_Paroxysmal2_Scaled(i,:)';

    PCA1vectVAR_Paroxysmal2_Sim(i) = coeffVAR_Paroxysmal2(:,1)'*dataToClassify_Paroxysmal2_Scaled_Sim(i,:)';
    PCA2vectVAR_Paroxysmal2_Sim(i) = coeffVAR_Paroxysmal2(:,2)'*dataToClassify_Paroxysmal2_Scaled_Sim(i,:)';
    PCA3vectVAR_Paroxysmal2_Sim(i) = coeffVAR_Paroxysmal2(:,3)'*dataToClassify_Paroxysmal2_Scaled_Sim(i,:)';
    PCA4vectVAR_Paroxysmal2_Sim(i) = coeffVAR_Paroxysmal2(:,4)'*dataToClassify_Paroxysmal2_Scaled_Sim(i,:)';
    PCA5vectVAR_Paroxysmal2_Sim(i) = coeffVAR_Paroxysmal2(:,5)'*dataToClassify_Paroxysmal2_Scaled_Sim(i,:)';
end


opts = statset('Display','final');
[idxPCAVAR_Paroxysmal2,CPCAVAR_Paroxysmal2,sumdPCAVAR_Paroxysmal2,DPCAVAR_Paroxysmal2] = kmeans([PCA1vectVAR_Paroxysmal2', PCA2vectVAR_Paroxysmal2', PCA3vectVAR_Paroxysmal2', PCA4vectVAR_Paroxysmal2', PCA5vectVAR_Paroxysmal2'],2,'Distance','cosine','EmptyAction','drop','MaxIter',1000,'Replicates',10,'Options',opts);

indexCluster1Paroxysmal2_Var = find(idxPCAVAR_Paroxysmal2 == 1);
indexCluster2Paroxysmal2_Var = find(idxPCAVAR_Paroxysmal2 == 2);
indexCluster3Paroxysmal2_Var = find(idxPCAVAR_Paroxysmal2 == 3);
indexCluster4Paroxysmal2_Var = find(idxPCAVAR_Paroxysmal2 == 4);
indexCluster5Paroxysmal2_Var = find(idxPCAVAR_Paroxysmal2 == 5);
indexCluster6Paroxysmal2_Var = find(idxPCAVAR_Paroxysmal2 == 6);
indexCluster7Paroxysmal2_Var = find(idxPCAVAR_Paroxysmal2 == 7);
indexCluster8Paroxysmal2_Var = find(idxPCAVAR_Paroxysmal2 == 8);






dataToClassify_Paroxysmal3 = [ downstroke_Paroxysmal3, upstrokeRise_Paroxysmal3, upstrokeBack_Paroxysmal3, signalWidth_Paroxysmal3, baseline_Paroxysmal3, EGMDuration_Paroxysmal3, AlldistanceToProj_Paroxysmal3, PPA_Paroxysmal3, deflectionCount_Paroxysmal3];
dataToClassify_Paroxysmal3_Scaled = dataToClassify_Paroxysmal3;

maxScalingData = [];
minScalingData = [];
% % % % SCALING TO MAKE ALL SCALES MATCH
for i = 1:length( dataToClassify_Paroxysmal3(1,:) )
    if max(dataToClassify_Paroxysmal3(:,i)) - min(dataToClassify_Paroxysmal3(:,i)) > 0
        dataToClassify_Paroxysmal3_Scaled(:,i) = ( dataToClassify_Paroxysmal3(:,i) - min(dataToClassify_Paroxysmal3(:,i)) )./( max(dataToClassify_Paroxysmal3(:,i)) - min(dataToClassify_Paroxysmal3(:,i)) );
        
        maxScalingData(i,2) = max(dataToClassify_Paroxysmal3(:,i));
        minScalingData(i,2) = min(dataToClassify_Paroxysmal3(:,i));
    end
end


maxScalingData = max(maxScalingData');
minScalingData = min(minScalingData');

[coeffVAR_Paroxysmal3,scoreVAR_Paroxysmal3,latentVAR_Paroxysmal3,tsquaredVAR_Paroxysmal3,explainedVAR_Paroxysmal3,muVAR_Paroxysmal3] = pca(dataToClassify_Paroxysmal3_Scaled);

sum( explainedVAR_Paroxysmal3(1:5) )

%More important variables are 8, 4, 6, 7

evaluationPCAVAR_Paroxysmal3 = evalclusters(scoreVAR_Paroxysmal3(:,1:5),"kmeans","CalinskiHarabasz","KList",1:20)
figure; plot(evaluationPCAVAR_Paroxysmal3)

dataToClassify_Paroxysmal3_Sim = [ downstroke_Paroxysmal3_Sim, upstrokeRise_Paroxysmal3_Sim, upstrokeBack_Paroxysmal3_Sim, signalWidth_Paroxysmal3_Sim, baseline_Paroxysmal3_Sim, EGMDuration_Paroxysmal3_Sim, AlldistanceToProj_Paroxysmal3_Sim, PPA_Paroxysmal3_Sim, deflectionCount_Paroxysmal3_Sim]; 
dataToClassify_Paroxysmal3_Scaled_Sim = dataToClassify_Paroxysmal3_Sim;
for i = 1:length( dataToClassify_Paroxysmal3_Sim(1,:) )
        dataToClassify_Paroxysmal3_Scaled_Sim(:,i) = ( dataToClassify_Paroxysmal3_Sim(:,i) - minScalingData(i) )./( maxScalingData(i) - minScalingData(i) );
end


PCA1vectVAR_Paroxysmal3 = [];
PCA2vectVAR_Paroxysmal3 = [];
PCA3vectVAR_Paroxysmal3 = [];
PCA4vectVAR_Paroxysmal3 = [];
PCA5vectVAR_Paroxysmal3 = [];

PCA1vectVAR_Paroxysmal3_Sim = [];
PCA2vectVAR_Paroxysmal3_Sim = [];
PCA3vectVAR_Paroxysmal3_Sim = [];
PCA4vectVAR_Paroxysmal3_Sim = [];
PCA5vectVAR_Paroxysmal3_Sim = [];

for i = 1: length(dataToClassify_Paroxysmal3_Scaled(:,1))
    PCA1vectVAR_Paroxysmal3(i) = coeffVAR_Paroxysmal3(:,1)'*dataToClassify_Paroxysmal3_Scaled(i,:)';
    PCA2vectVAR_Paroxysmal3(i) = coeffVAR_Paroxysmal3(:,2)'*dataToClassify_Paroxysmal3_Scaled(i,:)';
    PCA3vectVAR_Paroxysmal3(i) = coeffVAR_Paroxysmal3(:,3)'*dataToClassify_Paroxysmal3_Scaled(i,:)';
    PCA4vectVAR_Paroxysmal3(i) = coeffVAR_Paroxysmal3(:,4)'*dataToClassify_Paroxysmal3_Scaled(i,:)';
    PCA5vectVAR_Paroxysmal3(i) = coeffVAR_Paroxysmal3(:,5)'*dataToClassify_Paroxysmal3_Scaled(i,:)';

    PCA1vectVAR_Paroxysmal3_Sim(i) = coeffVAR_Paroxysmal3(:,1)'*dataToClassify_Paroxysmal3_Scaled_Sim(i,:)';
    PCA2vectVAR_Paroxysmal3_Sim(i) = coeffVAR_Paroxysmal3(:,2)'*dataToClassify_Paroxysmal3_Scaled_Sim(i,:)';
    PCA3vectVAR_Paroxysmal3_Sim(i) = coeffVAR_Paroxysmal3(:,3)'*dataToClassify_Paroxysmal3_Scaled_Sim(i,:)';
    PCA4vectVAR_Paroxysmal3_Sim(i) = coeffVAR_Paroxysmal3(:,4)'*dataToClassify_Paroxysmal3_Scaled_Sim(i,:)';
    PCA5vectVAR_Paroxysmal3_Sim(i) = coeffVAR_Paroxysmal3(:,5)'*dataToClassify_Paroxysmal3_Scaled_Sim(i,:)';
end


opts = statset('Display','final');
[idxPCAVAR_Paroxysmal3,CPCAVAR_Paroxysmal3,sumdPCAVAR_Paroxysmal3,DPCAVAR_Paroxysmal3] = kmeans([PCA1vectVAR_Paroxysmal3', PCA2vectVAR_Paroxysmal3', PCA3vectVAR_Paroxysmal3', PCA4vectVAR_Paroxysmal3', PCA5vectVAR_Paroxysmal3'],2,'Distance','cosine','EmptyAction','drop','MaxIter',1000,'Replicates',10,'Options',opts);

indexCluster1Paroxysmal3_Var = find(idxPCAVAR_Paroxysmal3 == 1);
indexCluster2Paroxysmal3_Var = find(idxPCAVAR_Paroxysmal3 == 2);
indexCluster3Paroxysmal3_Var = find(idxPCAVAR_Paroxysmal3 == 3);
indexCluster4Paroxysmal3_Var = find(idxPCAVAR_Paroxysmal3 == 4);
indexCluster5Paroxysmal3_Var = find(idxPCAVAR_Paroxysmal3 == 5);
indexCluster6Paroxysmal3_Var = find(idxPCAVAR_Paroxysmal3 == 6);
indexCluster7Paroxysmal3_Var = find(idxPCAVAR_Paroxysmal3 == 7);
indexCluster8Paroxysmal3_Var = find(idxPCAVAR_Paroxysmal3 == 8);





figure
scatter3(PCA1vectWD_Paroxysmal2(indexCluster1Paroxysmal2_WD), PCA2vectWD_Paroxysmal2(indexCluster1Paroxysmal2_WD), PCA3vectWD_Paroxysmal2(indexCluster1Paroxysmal2_WD), 1800.*abs(PCA4vectWD_Paroxysmal2(indexCluster1Paroxysmal2_WD)), 'r', 'filled', 'MarkerEdgeColor', 'k','DisplayName','1')
hold on
scatter3(PCA1vectWD_Paroxysmal2(indexCluster2Paroxysmal2_WD), PCA2vectWD_Paroxysmal2(indexCluster2Paroxysmal2_WD), PCA3vectWD_Paroxysmal2(indexCluster2Paroxysmal2_WD), 1800.*abs(PCA4vectWD_Paroxysmal2(indexCluster2Paroxysmal2_WD)), 'b', 'filled', 'MarkerEdgeColor', 'k','DisplayName','2')
scatter3(PCA1vectWD_Paroxysmal2(indexCluster3Paroxysmal2_WD), PCA2vectWD_Paroxysmal2(indexCluster3Paroxysmal2_WD), PCA3vectWD_Paroxysmal2(indexCluster3Paroxysmal2_WD), 1800.*abs(PCA4vectWD_Paroxysmal2(indexCluster3Paroxysmal2_WD)), 'm', 'filled', 'MarkerEdgeColor', 'k','DisplayName','3')
scatter3(PCA1vectWD_Paroxysmal2_Sim, PCA2vectWD_Paroxysmal2_Sim, PCA3vectWD_Paroxysmal2_Sim, 1800.*abs(PCA4vectWD_Paroxysmal2_Sim), 'c', 'filled', 'MarkerEdgeColor', 'k','DisplayName','Simulation')
legend()
xlabel('PC 1')
ylabel('PC 2')
zlabel('PC 3')
set(gca,'FontSize',35)


figure
scatter3(PCA1vectWD_Paroxysmal3(indexCluster1Paroxysmal3_WD), PCA2vectWD_Paroxysmal3(indexCluster1Paroxysmal3_WD), PCA3vectWD_Paroxysmal3(indexCluster1Paroxysmal3_WD), 1800.*abs(PCA4vectWD_Paroxysmal3(indexCluster1Paroxysmal3_WD)), 'r', 'filled', 'MarkerEdgeColor', 'k','DisplayName','1')
hold on
scatter3(PCA1vectWD_Paroxysmal3(indexCluster2Paroxysmal3_WD), PCA2vectWD_Paroxysmal3(indexCluster2Paroxysmal3_WD), PCA3vectWD_Paroxysmal3(indexCluster2Paroxysmal3_WD), 1800.*abs(PCA4vectWD_Paroxysmal3(indexCluster2Paroxysmal3_WD)), 'b', 'filled', 'MarkerEdgeColor', 'k','DisplayName','2')
scatter3(PCA1vectWD_Paroxysmal3(indexCluster3Paroxysmal3_WD), PCA2vectWD_Paroxysmal3(indexCluster3Paroxysmal3_WD), PCA3vectWD_Paroxysmal3(indexCluster3Paroxysmal3_WD), 1800.*abs(PCA4vectWD_Paroxysmal3(indexCluster3Paroxysmal3_WD)), 'm', 'filled', 'MarkerEdgeColor', 'k','DisplayName','3')
scatter3(PCA1vectWD_Paroxysmal3_Sim, PCA2vectWD_Paroxysmal3_Sim, PCA3vectWD_Paroxysmal3_Sim, 1800.*abs(PCA4vectWD_Paroxysmal3_Sim), 'c', 'filled', 'MarkerEdgeColor', 'k','DisplayName','Simulation')
legend()
xlabel('PC 1')
ylabel('PC 2')
zlabel('PC 3')
set(gca,'FontSize',35)





figure
scatter3(PCA1vectVAR_Paroxysmal2(indexCluster1Paroxysmal2_Var), PCA2vectVAR_Paroxysmal2(indexCluster1Paroxysmal2_Var), PCA3vectVAR_Paroxysmal2(indexCluster1Paroxysmal2_Var), 1800.*abs(PCA4vectVAR_Paroxysmal2(indexCluster1Paroxysmal2_Var)), 'r', 'filled', 'MarkerEdgeColor', 'k','DisplayName','1')
hold on
scatter3(PCA1vectVAR_Paroxysmal2(indexCluster2Paroxysmal2_Var), PCA2vectVAR_Paroxysmal2(indexCluster2Paroxysmal2_Var), PCA3vectVAR_Paroxysmal2(indexCluster2Paroxysmal2_Var), 1800.*abs(PCA4vectVAR_Paroxysmal2(indexCluster2Paroxysmal2_Var)), 'b', 'filled', 'MarkerEdgeColor', 'k','DisplayName','2')
scatter3(PCA1vectVAR_Paroxysmal2_Sim, PCA2vectVAR_Paroxysmal2_Sim, PCA3vectVAR_Paroxysmal2_Sim, 1800.*abs(PCA4vectVAR_Paroxysmal2_Sim), 'c', 'filled', 'MarkerEdgeColor', 'k','DisplayName','Simulation')
legend()
xlabel('PC 1')
ylabel('PC 2')
zlabel('PC 3')
set(gca,'FontSize',35)


figure
scatter3(PCA1vectVAR_Paroxysmal3(indexCluster1Paroxysmal3_Var), PCA2vectVAR_Paroxysmal3(indexCluster1Paroxysmal3_Var), PCA3vectVAR_Paroxysmal3(indexCluster1Paroxysmal3_Var), 1800.*abs(PCA4vectVAR_Paroxysmal3(indexCluster1Paroxysmal3_Var)), 'r', 'filled', 'MarkerEdgeColor', 'k','DisplayName','1')
hold on
scatter3(PCA1vectVAR_Paroxysmal3(indexCluster2Paroxysmal3_Var), PCA2vectVAR_Paroxysmal3(indexCluster2Paroxysmal3_Var), PCA3vectVAR_Paroxysmal3(indexCluster2Paroxysmal3_Var), 1800.*abs(PCA4vectVAR_Paroxysmal3(indexCluster2Paroxysmal3_Var)), 'b', 'filled', 'MarkerEdgeColor', 'k','DisplayName','2')
scatter3(PCA1vectVAR_Paroxysmal3_Sim, PCA2vectVAR_Paroxysmal3_Sim, PCA3vectVAR_Paroxysmal3_Sim, 1800.*abs(PCA4vectVAR_Paroxysmal3_Sim), 'c', 'filled', 'MarkerEdgeColor', 'k','DisplayName','Simulation')
legend()
xlabel('PC 1')
ylabel('PC 2')
zlabel('PC 3')
set(gca,'FontSize',35)


figure
% for i = 1:2:180
indToUse = i;
subplot(3,1,1)
yyaxis left
plot(signals_Paroxysmal2(indexCluster1Paroxysmal2_WD(67),:)', 'r-', 'DisplayName','1-Patient', 'LineWidth',3)
ylabel('Ve (mV)')
ylim([min(signals_Paroxysmal2(indexCluster1Paroxysmal2_WD(67),:)) max(signals_Paroxysmal2(indexCluster1Paroxysmal2_WD(67),:))])
hold on
yyaxis right
plot(signals_Paroxysmal2_Sim(indexCluster1Paroxysmal2_WD(67),:)', 'r--', 'DisplayName','1-Simulated', 'LineWidth',3)
hold off
xlabel('Time (msec)')
ylabel('Ve (mV)')
ylim([min(signals_Paroxysmal2_Sim(indexCluster1Paroxysmal2_WD(67),:)) max(signals_Paroxysmal2_Sim(indexCluster1Paroxysmal2_WD(67),:))])
legend()
set(gca,'FontSize',35)
subplot(3,1,2)
yyaxis left
plot(signals_Paroxysmal2(indexCluster2Paroxysmal2_WD(17),:)', 'b-', 'DisplayName','1-Patient', 'LineWidth',3)
ylabel('Ve (mV)')
ylim([min(signals_Paroxysmal2(indexCluster2Paroxysmal2_WD(17),:)) max(signals_Paroxysmal2(indexCluster2Paroxysmal2_WD(17),:))])
hold on
yyaxis right
plot(signals_Paroxysmal2_Sim(indexCluster2Paroxysmal2_WD(17),:)', 'b--', 'DisplayName','1-Simulated', 'LineWidth',3)
hold off
xlabel('Time (msec)')
ylabel('Ve (mV)')
ylim([min(signals_Paroxysmal2_Sim(indexCluster2Paroxysmal2_WD(17),:)) max(signals_Paroxysmal2_Sim(indexCluster2Paroxysmal2_WD(17),:))])
legend()
set(gca,'FontSize',35)
subplot(3,1,3)
yyaxis left
plot(signals_Paroxysmal2(indexCluster3Paroxysmal2_WD(51),:)', 'm-', 'DisplayName','1-Patient', 'LineWidth',3)
ylabel('Ve (mV)')
ylim([min(signals_Paroxysmal2(indexCluster3Paroxysmal2_WD(51),:)) max(signals_Paroxysmal2(indexCluster3Paroxysmal2_WD(51),:))])
hold on
yyaxis right
plot(signals_Paroxysmal2_Sim(indexCluster3Paroxysmal2_WD(51),:)', 'm--', 'DisplayName','1-Simulated', 'LineWidth',3)
hold off
xlabel('Time (msec)')
ylabel('Ve (mV)')
ylim([min(signals_Paroxysmal2_Sim(indexCluster3Paroxysmal2_WD(51),:)) max(signals_Paroxysmal2_Sim(indexCluster3Paroxysmal2_WD(51),:))])
legend()
set(gca,'FontSize',35)
% pause(3)
% end







figure
% for i = 1:2:180
indToUse = i;
subplot(3,1,1)
yyaxis left
plot(signals_Paroxysmal3(indexCluster1Paroxysmal3_WD(55),:)', 'r-', 'DisplayName','1-Patient', 'LineWidth',3)
ylabel('Ve (mV)')
ylim([min(signals_Paroxysmal3(indexCluster1Paroxysmal3_WD(55),:)) max(signals_Paroxysmal3(indexCluster1Paroxysmal3_WD(55),:))])
hold on
yyaxis right
plot(signals_Paroxysmal3_Sim(indexCluster1Paroxysmal3_WD(55),:)', 'r--', 'DisplayName','1-Simulated', 'LineWidth',3)
hold off
xlabel('Time (msec)')
ylabel('Ve (mV)')
ylim([min(signals_Paroxysmal3_Sim(indexCluster1Paroxysmal3_WD(55),:)) max(signals_Paroxysmal3_Sim(indexCluster1Paroxysmal3_WD(55),:))])
legend()
set(gca,'FontSize',35)
subplot(3,1,2)
yyaxis left
plot(signals_Paroxysmal3(indexCluster2Paroxysmal3_WD(15),:)', 'b-', 'DisplayName','1-Patient', 'LineWidth',3)
ylabel('Ve (mV)')
ylim([min(signals_Paroxysmal3(indexCluster2Paroxysmal3_WD(15),:)) max(signals_Paroxysmal3(indexCluster2Paroxysmal3_WD(15),:))])
hold on
yyaxis right
plot(signals_Paroxysmal3_Sim(indexCluster2Paroxysmal3_WD(15),:)', 'b--', 'DisplayName','1-Simulated', 'LineWidth',3)
hold off
xlabel('Time (msec)')
ylabel('Ve (mV)')
ylim([min(signals_Paroxysmal3_Sim(indexCluster2Paroxysmal3_WD(15),:)) max(signals_Paroxysmal3_Sim(indexCluster2Paroxysmal3_WD(15),:))])
legend()
set(gca,'FontSize',35)
subplot(3,1,3)
yyaxis left
plot(signals_Paroxysmal3(indexCluster3Paroxysmal3_WD(35),:)', 'm-', 'DisplayName','1-Patient', 'LineWidth',3)
ylabel('Ve (mV)')
ylim([min(signals_Paroxysmal3(indexCluster3Paroxysmal3_WD(35),:)) max(signals_Paroxysmal3(indexCluster3Paroxysmal3_WD(35),:))])
hold on
yyaxis right
plot(signals_Paroxysmal3_Sim(indexCluster3Paroxysmal3_WD(35),:)', 'm--', 'DisplayName','1-Simulated', 'LineWidth',3)
hold off
xlabel('Time (msec)')
ylabel('Ve (mV)')
ylim([min(signals_Paroxysmal3_Sim(indexCluster3Paroxysmal3_WD(35),:)) max(signals_Paroxysmal3_Sim(indexCluster3Paroxysmal3_WD(35),:))])
legend()
set(gca,'FontSize',35)
% pause(3)
% end






figure
% for i = 1:2:180
indToUse = i;
subplot(2,1,1)
yyaxis left
plot(signals_Paroxysmal2(indexCluster1Paroxysmal2_Var(27),:)', 'r-', 'DisplayName','1-Patient', 'LineWidth',3)
ylabel('Ve (mV)')
ylim([min(signals_Paroxysmal2(indexCluster1Paroxysmal2_Var(27),:)) max(signals_Paroxysmal2(indexCluster1Paroxysmal2_Var(27),:))])
hold on
yyaxis right
plot(signals_Paroxysmal2_Sim(indexCluster1Paroxysmal2_Var(27),:)', 'r--', 'DisplayName','1-Simulated', 'LineWidth',3)
hold off
xlabel('Time (msec)')
ylabel('Ve (mV)')
ylim([min(signals_Paroxysmal2_Sim(indexCluster1Paroxysmal2_Var(27),:)) max(signals_Paroxysmal2_Sim(indexCluster1Paroxysmal2_Var(27),:))])
legend()
set(gca,'FontSize',35)
subplot(2,1,2)
yyaxis left
plot(signals_Paroxysmal2(indexCluster2Paroxysmal2_Var(59),:)', 'b-', 'DisplayName','1-Patient', 'LineWidth',3)
ylabel('Ve (mV)')
ylim([min(signals_Paroxysmal2(indexCluster2Paroxysmal2_Var(59),:)) max(signals_Paroxysmal2(indexCluster2Paroxysmal2_Var(59),:))])
hold on
yyaxis right
plot(signals_Paroxysmal2_Sim(indexCluster2Paroxysmal2_Var(59),:)', 'b--', 'DisplayName','1-Simulated', 'LineWidth',3)
hold off
xlabel('Time (msec)')
ylabel('Ve (mV)')
ylim([min(signals_Paroxysmal2_Sim(indexCluster2Paroxysmal2_Var(59),:)) max(signals_Paroxysmal2_Sim(indexCluster2Paroxysmal2_Var(59),:))])
legend()
set(gca,'FontSize',35)
% pause(3)
% end







figure
% for i = 1:2:180
indToUse = i;
subplot(2,1,1)
yyaxis left
plot(signals_Paroxysmal3(indexCluster1Paroxysmal3_Var(97),:)', 'r-', 'DisplayName','1-Patient', 'LineWidth',3)
ylabel('Ve (mV)')
ylim([min(signals_Paroxysmal3(indexCluster1Paroxysmal3_Var(97),:)) max(signals_Paroxysmal3(indexCluster1Paroxysmal3_Var(97),:))])
hold on
yyaxis right
plot(signals_Paroxysmal3_Sim(indexCluster1Paroxysmal3_Var(97),:)', 'r--', 'DisplayName','1-Simulated', 'LineWidth',3)
hold off
xlabel('Time (msec)')
ylabel('Ve (mV)')
ylim([min(signals_Paroxysmal3_Sim(indexCluster1Paroxysmal3_Var(97),:)) max(signals_Paroxysmal3_Sim(indexCluster1Paroxysmal3_Var(97),:))])
legend()
set(gca,'FontSize',35)
subplot(2,1,2)
yyaxis left
plot(signals_Paroxysmal3(indexCluster2Paroxysmal3_Var(13),:)', 'b-', 'DisplayName','1-Patient', 'LineWidth',3)
ylabel('Ve (mV)')
ylim([min(signals_Paroxysmal3(indexCluster2Paroxysmal3_Var(13),:)) max(signals_Paroxysmal3(indexCluster2Paroxysmal3_Var(13),:))])
hold on
yyaxis right
plot(signals_Paroxysmal3_Sim(indexCluster2Paroxysmal3_Var(13),:)', 'b--', 'DisplayName','1-Simulated', 'LineWidth',3)
hold off
xlabel('Time (msec)')
ylabel('Ve (mV)')
ylim([min(signals_Paroxysmal3_Sim(indexCluster2Paroxysmal3_Var(13),:)) max(signals_Paroxysmal3_Sim(indexCluster2Paroxysmal3_Var(13),:))])
legend()
set(gca,'FontSize',35)
% pause(3)
% end







%% SIGNAL TO SIGNAL COMPARISON


loc_similarity_score_Paroxysmal2 = [];
loc_xcorr_Paroxysmal2 = [];
loc_similarity_score_Paroxysmal2_OG = [];
loc_xcorr_Paroxysmal2_OG = [];
loc_LAT_error_Paroxysmal2 = [];

for i = 1:length(signals_Paroxysmal2(:,1))
    correlation_coefficient = corrcoef(signals_Paroxysmal2(i,:), signals_Paroxysmal2_Sim(i,:));
    [xc,lags] = xcorr(signals_Paroxysmal2(i,:), signals_Paroxysmal2_Sim(i,:),10,'normalized');
    loc_similarity_score_Paroxysmal2(i) = correlation_coefficient(1, 2);
    loc_xcorr_Paroxysmal2(i) = max(xc);

    correlation_coefficient = corrcoef(egmsParoxysmal2(:,i), simulatedParoxysmal2_downsampled(1:timePointsToTake,i));
    [xc,lags] = xcorr(egmsParoxysmal2(:,i), simulatedParoxysmal2_downsampled(1:timePointsToTake,i),10,'normalized');
    loc_similarity_score_Paroxysmal2_OG(i) = correlation_coefficient(1, 2);
    loc_xcorr_Paroxysmal2_OG(i) = max(xc);

    loc_LAT_error_Paroxysmal2(i) = (abs( actTimes2(i) - actTimesSim2(i) )./actTimes2(i))*100;
end

loc_similarity_score_Paroxysmal3 = [];
loc_xcorr_Paroxysmal3 = [];
loc_similarity_score_Paroxysmal3_OG = [];
loc_xcorr_Paroxysmal3_OG = [];
loc_LAT_error_Paroxysmal3 = [];

for i = 1:length(signals_Paroxysmal3(:,1))
    correlation_coefficient = corrcoef(signals_Paroxysmal3(i,:), signals_Paroxysmal3_Sim(i,:));
    [xc,lags] = xcorr(signals_Paroxysmal3(i,:), signals_Paroxysmal3_Sim(i,:),10,'normalized');
    loc_similarity_score_Paroxysmal3(i) = correlation_coefficient(1, 2);
    loc_xcorr_Paroxysmal3(i) = max(xc);

    correlation_coefficient = corrcoef(egmsParoxysmal3(:,i), simulatedParoxysmal3_downsampled(1:timePointsToTake,i));
    [xc,lags] = xcorr(egmsParoxysmal3(:,i), simulatedParoxysmal3_downsampled(1:timePointsToTake,i),10,'normalized');
    loc_similarity_score_Paroxysmal3_OG(i) = correlation_coefficient(1, 2);
    loc_xcorr_Paroxysmal3_OG(i) = max(xc);

    loc_LAT_error_Paroxysmal3(i) = (abs( actTimes3(i) - actTimesSim3(i) )./actTimes3(i))*100;
end



ppA_error_Paroxysmal2 = abs((100.*((dataToClassify_Paroxysmal2(:,8) - dataToClassify_Paroxysmal2_Sim(:,8))./(dataToClassify_Paroxysmal2(:,8)))));
dC_error_Paroxysmal2 = abs((100.*((dataToClassify_Paroxysmal2(:,9) - dataToClassify_Paroxysmal2_Sim(:,9))./(dataToClassify_Paroxysmal2(:,9)))));
EGMD_error_Paroxysmal2 = abs((100.*((dataToClassify_Paroxysmal2(:,6) - dataToClassify_Paroxysmal2_Sim(:,6))./(dataToClassify_Paroxysmal2(:,6)))));
similarityScore_mean_Paroxysmal2 = (loc_similarity_score_Paroxysmal2);
xcorr_mean_Paroxysmal2 = (loc_xcorr_Paroxysmal2);
similarityScore_mean_Paroxysmal2_OG = (loc_similarity_score_Paroxysmal2_OG);
xcorr_mean_Paroxysmal2_OG = (loc_xcorr_Paroxysmal2_OG);

ppA_error_Paroxysmal3 = abs((100.*((dataToClassify_Paroxysmal3(:,8) - dataToClassify_Paroxysmal3_Sim(:,8))./(dataToClassify_Paroxysmal3(:,8)))));
dC_error_Paroxysmal3 = abs((100.*((dataToClassify_Paroxysmal3(:,9) - dataToClassify_Paroxysmal3_Sim(:,9))./(dataToClassify_Paroxysmal3(:,9)))));
EGMD_error_Paroxysmal3 = abs((100.*((dataToClassify_Paroxysmal3(:,6) - dataToClassify_Paroxysmal3_Sim(:,6))./(dataToClassify_Paroxysmal3(:,6)))));
similarityScore_mean_Paroxysmal3 = (loc_similarity_score_Paroxysmal3);
xcorr_mean_Paroxysmal3 = (loc_xcorr_Paroxysmal3);
similarityScore_mean_Paroxysmal3_OG = (loc_similarity_score_Paroxysmal3_OG);
xcorr_mean_Paroxysmal3_OG = (loc_xcorr_Paroxysmal3_OG);

table_Paroxysmal2_SSC_MeanCh = {};
for i = 1:9
    table_Paroxysmal2_SSC_MeanCh{i,1} =  createMeanStdForCell(dataToClassify_Paroxysmal2(:,i), 1:10, 1:10, 2 );
    table_Paroxysmal2_SSC_MeanCh{i,2} =  createMeanStdForCell(dataToClassify_Paroxysmal2_Sim(:,i), 1:10, 1:10, 2 );
    table_Paroxysmal2_SSC_MeanCh{i,3} =  createMeanStdForCell(dataToClassify_Paroxysmal3(:,i), 1:10, 1:10, 2 );
    table_Paroxysmal2_SSC_MeanCh{i,4} =  createMeanStdForCell(dataToClassify_Paroxysmal3_Sim(:,i), 1:10, 1:10, 2 );
end
table_Paroxysmal2_SSC_MeanCh = cell2table(table_Paroxysmal2_SSC_MeanCh);


table_Paroxysmal2_SSC_ = {};
table_Paroxysmal2_SSC_{1,1} =  createMeanStdForCell(ppA_error_Paroxysmal2, 1:10, 1:10, 2 );
table_Paroxysmal2_SSC_{2,1} =  createMeanStdForCell(dC_error_Paroxysmal2, 1:10, 1:10, 2 );
table_Paroxysmal2_SSC_{3,1} =  createMeanStdForCell(EGMD_error_Paroxysmal2, 1:10, 1:10, 2 );
table_Paroxysmal2_SSC_{4,1} =  createMeanStdForCell(similarityScore_mean_Paroxysmal2, 1:10, 1:10, 2 );
table_Paroxysmal2_SSC_{5,1} =  createMeanStdForCell(similarityScore_mean_Paroxysmal2_OG, 1:10, 1:10, 2 );
table_Paroxysmal2_SSC_{6,1} =  createMeanStdForCell(xcorr_mean_Paroxysmal2, 1:10, 1:10, 2 );
table_Paroxysmal2_SSC_{7,1} =  createMeanStdForCell(xcorr_mean_Paroxysmal2_OG, 1:10, 1:10, 2 );
table_Paroxysmal2_SSC_{8,1} =  createMeanStdForCell(loc_LAT_error_Paroxysmal2, 1:10, 1:10, 2 );
table_Paroxysmal2_SSC_{1,2} =  createMeanStdForCell(ppA_error_Paroxysmal3, 1:10, 1:10, 2 );
table_Paroxysmal2_SSC_{2,2} =  createMeanStdForCell(dC_error_Paroxysmal3, 1:10, 1:10, 2 );
table_Paroxysmal2_SSC_{3,2} =  createMeanStdForCell(EGMD_error_Paroxysmal3, 1:10, 1:10, 2 );
table_Paroxysmal2_SSC_{4,2} =  createMeanStdForCell(similarityScore_mean_Paroxysmal3, 1:10, 1:10, 2 );
table_Paroxysmal2_SSC_{5,2} =  createMeanStdForCell(similarityScore_mean_Paroxysmal3_OG, 1:10, 1:10, 2 );
table_Paroxysmal2_SSC_{6,2} =  createMeanStdForCell(xcorr_mean_Paroxysmal3, 1:10, 1:10, 2 );
table_Paroxysmal2_SSC_{7,2} =  createMeanStdForCell(xcorr_mean_Paroxysmal3_OG, 1:10, 1:10, 2 );
table_Paroxysmal2_SSC_{8,2} =  createMeanStdForCell(loc_LAT_error_Paroxysmal3, 1:10, 1:10, 2 );
table_Paroxysmal2_SSC_ = cell2table(table_Paroxysmal2_SSC_);


%correlation coefficient NORMALIZED
figure
plot3(patient2xSurfaceCm, patient2ySurfaceCm, patient2zSurfaceCm, 'ro', 'DisplayName', 'tissue', 'MarkerFaceColor','r', 'MarkerSize',20,'LineWidth',3);
hold on
scatter3(patient2xCm, patient2yCm, patient2zCm, 230, loc_similarity_score_Paroxysmal2, 'filled', 'MarkerEdgeColor','k','LineWidth',3);
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'Correlation Coeff.';
% colormap(jet);
caxis([0, 1]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

%correlation coefficient OG
figure
plot3(patient2xSurfaceCm, patient2ySurfaceCm, patient2zSurfaceCm, 'ro', 'DisplayName', 'tissue', 'MarkerFaceColor','r', 'MarkerSize',20,'LineWidth',3);
hold on
scatter3(patient2xCm, patient2yCm, patient2zCm, 230, loc_similarity_score_Paroxysmal2_OG, 'filled', 'MarkerEdgeColor','k','LineWidth',3);
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'Correlation Coeff.';
% colormap(jet);
caxis([0, 1]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

%cross correlation peak NORMALIZED
figure
plot3(patient2xSurfaceCm, patient2ySurfaceCm, patient2zSurfaceCm, 'ro', 'DisplayName', 'tissue', 'MarkerFaceColor','r', 'MarkerSize',20,'LineWidth',3);
hold on
scatter3(patient2xCm, patient2yCm, patient2zCm, 230, loc_xcorr_Paroxysmal2, 'filled', 'MarkerEdgeColor','k','LineWidth',3);
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'X Correlation';
% colormap(jet);
caxis([0, 1]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

%cross correlation peak OG
figure
plot3(patient2xSurfaceCm, patient2ySurfaceCm, patient2zSurfaceCm, 'ro', 'DisplayName', 'tissue', 'MarkerFaceColor','r', 'MarkerSize',20,'LineWidth',3);
hold on
scatter3(patient2xCm, patient2yCm, patient2zCm, 230, loc_xcorr_Paroxysmal2_OG, 'filled', 'MarkerEdgeColor','k','LineWidth',3);
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'X Correlation';
% colormap(jet);
caxis([0, 1]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

%LAT error
figure
plot3(patient2xSurfaceCm, patient2ySurfaceCm, patient2zSurfaceCm, 'ro', 'DisplayName', 'tissue', 'MarkerFaceColor','r', 'MarkerSize',20,'LineWidth',3);
hold on
scatter3(patient2xCm, patient2yCm, patient2zCm, 230, loc_LAT_error_Paroxysmal2, 'filled', 'MarkerEdgeColor','k','LineWidth',3);
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LAT error %';
% colormap(jet);
caxis([0, 50]);
% title('H1A Patient');
set(gca, 'FontSize', 35);










%correlation coefficient NORMALIZED
figure
plot3(patient3xSurfaceCm, patient3ySurfaceCm, patient3zSurfaceCm, 'ro', 'DisplayName', 'tissue', 'MarkerFaceColor','r', 'MarkerSize',20,'LineWidth',3);
hold on
scatter3(patient3xCm, patient3yCm, patient3zCm, 230, loc_similarity_score_Paroxysmal3, 'filled', 'MarkerEdgeColor','k','LineWidth',3);
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'Correlation Coeff.';
% colormap(jet);
caxis([0, 1]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

%correlation coefficient OG
figure
plot3(patient3xSurfaceCm, patient3ySurfaceCm, patient3zSurfaceCm, 'ro', 'DisplayName', 'tissue', 'MarkerFaceColor','r', 'MarkerSize',20,'LineWidth',3);
hold on
scatter3(patient3xCm, patient3yCm, patient3zCm, 230, loc_similarity_score_Paroxysmal3_OG, 'filled', 'MarkerEdgeColor','k','LineWidth',3);
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'Correlation Coeff.';
% colormap(jet);
caxis([0, 1]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

%cross correlation peak NORMALIZED
figure
plot3(patient3xSurfaceCm, patient3ySurfaceCm, patient3zSurfaceCm, 'ro', 'DisplayName', 'tissue', 'MarkerFaceColor','r', 'MarkerSize',20,'LineWidth',3);
hold on
scatter3(patient3xCm, patient3yCm, patient3zCm, 230, loc_xcorr_Paroxysmal3, 'filled', 'MarkerEdgeColor','k','LineWidth',3);
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'X Correlation';
% colormap(jet);
caxis([0, 1]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

%cross correlation peak OG
figure
plot3(patient3xSurfaceCm, patient3ySurfaceCm, patient3zSurfaceCm, 'ro', 'DisplayName', 'tissue', 'MarkerFaceColor','r', 'MarkerSize',20,'LineWidth',3);
hold on
scatter3(patient3xCm, patient3yCm, patient3zCm, 230, loc_xcorr_Paroxysmal3_OG, 'filled', 'MarkerEdgeColor','k','LineWidth',3);
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'X Correlation';
% colormap(jet);
caxis([0, 1]);
% title('H1A Patient');
set(gca, 'FontSize', 35);

%LAT error
figure
plot3(patient3xSurfaceCm, patient3ySurfaceCm, patient3zSurfaceCm, 'ro', 'DisplayName', 'tissue', 'MarkerFaceColor','r', 'MarkerSize',20,'LineWidth',3);
hold on
scatter3(patient3xCm, patient3yCm, patient3zCm, 230, loc_LAT_error_Paroxysmal3, 'filled', 'MarkerEdgeColor','k','LineWidth',3);
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'LAT error %';
% colormap(jet);
caxis([0, 50]);
% title('H1A Patient');
set(gca, 'FontSize', 35);



%% PCA-BASED COMPARISON


%idxWD_Paroxysmal3,CWD_Paroxysmal3,sumdWD_Paroxysmal3,DWD_Paroxysmal3
distsVect_Patient_Paroxysmal3_WD = [];
distsVect_Simulated_common_Paroxysmal3_WD = [];
distsVect_Simulated_different_Paroxysmal3_WD = [];
percentage_Simulated_different_Paroxysmal3_WD = 0;

distsVect_Patient_Paroxysmal3_VAR = [];
distsVect_Simulated_common_Paroxysmal3_VAR = [];
distsVect_Simulated_different_Paroxysmal3_VAR = [];
percentage_Simulated_different_Paroxysmal3_VAR = 0;

distsVect_Patient_Paroxysmal2_WD = [];
distsVect_Simulated_common_Paroxysmal2_WD = [];
distsVect_Simulated_different_Paroxysmal2_WD = [];
percentage_Simulated_different_Paroxysmal2_WD = 0;

distsVect_Patient_Paroxysmal2_VAR = [];
distsVect_Simulated_common_Paroxysmal2_VAR = [];
distsVect_Simulated_different_Paroxysmal2_VAR = [];
percentage_Simulated_different_Paroxysmal2_VAR = 0;

for i = 1:length(DWD_Paroxysmal3(:,1))
    %WAVELET DECOMPOSITION
    currIndex = idxWD_Paroxysmal3(i);
    distsVect_Patient_Paroxysmal3_WD(i) = DWD_Paroxysmal3(i,currIndex);
    locDists = [];
    for j = 1:length(CWD_Paroxysmal3(:,1))
        locDists(j) = sqrt( (PCA1vectWD_Paroxysmal3_Sim(i) - CWD_Paroxysmal3(j,1)).^2 + (PCA2vectWD_Paroxysmal3_Sim(i) - CWD_Paroxysmal3(j,2)).^2 + (PCA3vectWD_Paroxysmal3_Sim(i) - CWD_Paroxysmal3(j,3)).^2 + (PCA4vectWD_Paroxysmal3_Sim(i) - CWD_Paroxysmal3(j,4)).^2 + (PCA5vectWD_Paroxysmal3_Sim(i) - CWD_Paroxysmal3(j,5)).^2  );
    end
    distsVect_Simulated_common_Paroxysmal3_WD(i) = locDists(currIndex);
    distsVect_Simulated_different_Paroxysmal3_WD(i) = min(locDists);
    if distsVect_Simulated_common_Paroxysmal3_WD(i) ~= distsVect_Simulated_different_Paroxysmal3_WD(i)
        percentage_Simulated_different_Paroxysmal3_WD = percentage_Simulated_different_Paroxysmal3_WD + 1;
    end


    %CHARACTERISTIC PARAMETERS
    currIndex = idxPCAVAR_Paroxysmal3(i);
    distsVect_Patient_Paroxysmal3_VAR(i) = DPCAVAR_Paroxysmal3(i,currIndex);
    locDists = [];
    for j = 1:length(CPCAVAR_Paroxysmal3(:,1))
        locDists(j) = sqrt( (PCA1vectVAR_Paroxysmal3_Sim(i) - CPCAVAR_Paroxysmal3(j,1)).^2 + (PCA2vectVAR_Paroxysmal3_Sim(i) - CPCAVAR_Paroxysmal3(j,2)).^2 + (PCA3vectVAR_Paroxysmal3_Sim(i) - CPCAVAR_Paroxysmal3(j,3)).^2 + (PCA4vectVAR_Paroxysmal3_Sim(i) - CPCAVAR_Paroxysmal3(j,4)).^2 + (PCA5vectVAR_Paroxysmal3_Sim(i) - CPCAVAR_Paroxysmal3(j,5)).^2  );
    end
    distsVect_Simulated_common_Paroxysmal3_VAR(i) = locDists(currIndex);
    distsVect_Simulated_different_Paroxysmal3_VAR(i) = min(locDists);
    if distsVect_Simulated_common_Paroxysmal3_VAR(i) ~= distsVect_Simulated_different_Paroxysmal3_VAR(i)
        percentage_Simulated_different_Paroxysmal3_VAR = percentage_Simulated_different_Paroxysmal3_VAR + 1;
    end

end

for i = 1:length(DWD_Paroxysmal2(:,1))
    %WAVELET DECOMPOSITION
    currIndex = idxWD_Paroxysmal2(i);
    distsVect_Patient_Paroxysmal2_WD(i) = DWD_Paroxysmal2(i,currIndex);
    locDists = [];
    for j = 1:length(CWD_Paroxysmal2(:,1))
        locDists(j) = sqrt( (PCA1vectWD_Paroxysmal2_Sim(i) - CWD_Paroxysmal2(j,1)).^2 + (PCA2vectWD_Paroxysmal2_Sim(i) - CWD_Paroxysmal2(j,2)).^2 + (PCA3vectWD_Paroxysmal2_Sim(i) - CWD_Paroxysmal2(j,3)).^2 + (PCA4vectWD_Paroxysmal2_Sim(i) - CWD_Paroxysmal2(j,4)).^2 + (PCA5vectWD_Paroxysmal2_Sim(i) - CWD_Paroxysmal2(j,5)).^2  );
    end
    distsVect_Simulated_common_Paroxysmal2_WD(i) = locDists(currIndex);
    distsVect_Simulated_different_Paroxysmal2_WD(i) = min(locDists);
    if distsVect_Simulated_common_Paroxysmal2_WD(i) ~= distsVect_Simulated_different_Paroxysmal2_WD(i)
        percentage_Simulated_different_Paroxysmal2_WD = percentage_Simulated_different_Paroxysmal2_WD + 1;
    end


    %CHARACTERISTIC PARAMETERS
    currIndex = idxPCAVAR_Paroxysmal2(i);
    distsVect_Patient_Paroxysmal2_VAR(i) = DPCAVAR_Paroxysmal2(i,currIndex);
    locDists = [];
    for j = 1:length(CPCAVAR_Paroxysmal2(:,1))
        locDists(j) = sqrt( (PCA1vectVAR_Paroxysmal2_Sim(i) - CPCAVAR_Paroxysmal2(j,1)).^2 + (PCA2vectVAR_Paroxysmal2_Sim(i) - CPCAVAR_Paroxysmal2(j,2)).^2 + (PCA3vectVAR_Paroxysmal2_Sim(i) - CPCAVAR_Paroxysmal2(j,3)).^2 + (PCA4vectVAR_Paroxysmal2_Sim(i) - CPCAVAR_Paroxysmal2(j,4)).^2 + (PCA5vectVAR_Paroxysmal2_Sim(i) - CPCAVAR_Paroxysmal2(j,5)).^2  );
    end
    distsVect_Simulated_common_Paroxysmal2_VAR(i) = locDists(currIndex);
    distsVect_Simulated_different_Paroxysmal2_VAR(i) = min(locDists);
    if distsVect_Simulated_common_Paroxysmal2_VAR(i) ~= distsVect_Simulated_different_Paroxysmal2_VAR(i)
        percentage_Simulated_different_Paroxysmal2_VAR = percentage_Simulated_different_Paroxysmal2_VAR + 1;
    end

end



percentage_Simulated_different_Paroxysmal3_WD = (percentage_Simulated_different_Paroxysmal3_WD/(length(idxWD_Paroxysmal3)))*100;
percentage_Simulated_different_Paroxysmal3_VAR = (percentage_Simulated_different_Paroxysmal3_VAR/(length(idxPCAVAR_Paroxysmal3)))*100;
percentage_Simulated_different_Paroxysmal2_WD = (percentage_Simulated_different_Paroxysmal2_WD/(length(idxWD_Paroxysmal2)))*100;
percentage_Simulated_different_Paroxysmal2_VAR = (percentage_Simulated_different_Paroxysmal2_VAR/(length(idxPCAVAR_Paroxysmal2)))*100;


table_Paroxysmal2_PCA_approach = {};
table_Paroxysmal2_PCA_approach{1,1} =  createMeanStdForCell(distsVect_Patient_Paroxysmal2_WD, 1:10, 1:10, 2 );
table_Paroxysmal2_PCA_approach{1,2} =  createMeanStdForCell(distsVect_Patient_Paroxysmal2_VAR, 1:10, 1:10, 2 );
table_Paroxysmal2_PCA_approach{2,1} =  createMeanStdForCell(distsVect_Simulated_common_Paroxysmal2_WD, 1:10, 1:10, 2 );
table_Paroxysmal2_PCA_approach{2,2} =  createMeanStdForCell(distsVect_Simulated_common_Paroxysmal2_VAR, 1:10, 1:10, 2 );
table_Paroxysmal2_PCA_approach{3,1} =  createMeanStdForCell(distsVect_Simulated_different_Paroxysmal2_WD, 1:10, 1:10, 2 );
table_Paroxysmal2_PCA_approach{3,2} =  createMeanStdForCell(distsVect_Simulated_different_Paroxysmal2_VAR, 1:10, 1:10, 2 );
table_Paroxysmal2_PCA_approach{4,1} = num2str(percentage_Simulated_different_Paroxysmal2_WD);
table_Paroxysmal2_PCA_approach{4,2} = num2str(percentage_Simulated_different_Paroxysmal2_VAR);
table_Paroxysmal2_PCA_approach = cell2table(table_Paroxysmal2_PCA_approach);

table_Paroxysmal3_PCA_approach = {};
table_Paroxysmal3_PCA_approach{1,1} =  createMeanStdForCell(distsVect_Patient_Paroxysmal3_WD, 1:10, 1:10, 2 );
table_Paroxysmal3_PCA_approach{1,2} =  createMeanStdForCell(distsVect_Patient_Paroxysmal3_VAR, 1:10, 1:10, 2 );
table_Paroxysmal3_PCA_approach{2,1} =  createMeanStdForCell(distsVect_Simulated_common_Paroxysmal3_WD, 1:10, 1:10, 2 );
table_Paroxysmal3_PCA_approach{2,2} =  createMeanStdForCell(distsVect_Simulated_common_Paroxysmal3_VAR, 1:10, 1:10, 2 );
table_Paroxysmal3_PCA_approach{3,1} =  createMeanStdForCell(distsVect_Simulated_different_Paroxysmal3_WD, 1:10, 1:10, 2 );
table_Paroxysmal3_PCA_approach{3,2} =  createMeanStdForCell(distsVect_Simulated_different_Paroxysmal3_VAR, 1:10, 1:10, 2 );
table_Paroxysmal3_PCA_approach{4,1} = num2str(percentage_Simulated_different_Paroxysmal3_WD);
table_Paroxysmal3_PCA_approach{4,2} = num2str(percentage_Simulated_different_Paroxysmal3_VAR);
table_Paroxysmal3_PCA_approach = cell2table(table_Paroxysmal3_PCA_approach);





%% WORKING SELECTING NODES TO STIMULATE
%FIRST WE SELECT TWO POINTS IN THE GENERAL REGION OF WHERE ACTIVATION
%HAPPENS FIRST OR ON EDGES FOR THE CODE TO SELECT ALL THE POINTS IN THAT
%GENERAL DIRECTION UNTIL IT REACHES THE CONNECTION BETWEEN THE ORIGINAL TWO
%POINTS. THEN, IT WILL CHECK WHERE ACTIVATION HAPPENS FIRST BASED ON LAT'S
%AND THEN SELECT A SET OF POINTS IN A RADIUS THAT ARE CLOSER TO THAT FIRST
%ACTIVTION POINTS FROM THE EGMS.

actTimesInterest = actTimes3;

electrodeXInterest = patient3xCm;
electrodeYInterest = patient3yCm;
electrodeZInterest = patient3zCm;

surfaceXInterest = patient3xSurfaceCm;
surfaceYInterest = patient3ySurfaceCm;
surfaceZInterest = patient3zSurfaceCm;

ptCloudSurface = pointCloud( [reshape(surfaceXInterest,[length(surfaceXInterest),1]), reshape(surfaceYInterest,[length(surfaceYInterest),1]), reshape(surfaceZInterest,[length(surfaceZInterest),1])] );

%CORNER POINTS WHERE STIMULUS WILL TAKE PLACE
% %PATIENT 3 Paroxysmal
point1ind = find( abs( surfaceXInterest - (-1.4724) ) < .001 & abs( surfaceYInterest - (-5.6861) ) < .001 & abs( surfaceZInterest - (11.2253) ) < .001 );
point2ind = find( abs( surfaceXInterest - (.8552) ) < .001 & abs( surfaceYInterest - (-4.1922) ) < .001 & abs( surfaceZInterest - (11.2943) ) < .001 );

% %PATIENT 2 Paroxysmal
% point1ind = find( abs( surfaceXInterest - (-2.0308) ) < .001 & abs( surfaceYInterest - (-4.9813) ) < .001 & abs( surfaceZInterest - (9.8328) ) < .001 );
% point2ind = find( abs( surfaceXInterest - (.2207) ) < .001 & abs( surfaceYInterest - (-3.3766) ) < .001 & abs( surfaceZInterest - (9.728) ) < .001 );


% stimulusRadius = .2;
allowedThreshAngle = 90; %degrees, to allow more points to be selected
allowedNumPoints = 15; %too small it won't select properly in the desired direction, too big then it will be too many points to stimulate

roitDirection = [ surfaceXInterest(point1ind) - surfaceXInterest(point2ind), surfaceYInterest(point1ind) - surfaceYInterest(point2ind), surfaceZInterest(point1ind) - surfaceZInterest(point2ind)  ];
% roitDirection = roitDirection./norm(roitDirection);

timeStim = min( actTimesInterest );
indexTimeStim = find( actTimesInterest <= timeStim );

initVar = true;
currPoint = point1ind;
newInd = [];
distVect = [];

while initVar
    [indices,dists] = findNearestNeighbors(ptCloudSurface,[ surfaceXInterest(currPoint), surfaceYInterest(currPoint), surfaceZInterest(currPoint) ],allowedNumPoints);
    workAngle = [];
    dummind = 1;
    for i = 2:length(indices)
        currDirect = [ surfaceXInterest(currPoint) - surfaceXInterest(indices(i)), surfaceYInterest(currPoint) - surfaceYInterest(indices(i)), surfaceZInterest(currPoint) - surfaceZInterest(indices(i))  ];
        % currDirect = currDirect./(norm(currDirect));
        currAngle = rad2deg(atan2(norm(cross(currDirect,roitDirection)), dot(currDirect,roitDirection)));
            if currAngle <= allowedThreshAngle
                workAngle(dummind,1) = currAngle;
                workAngle(dummind,2) = indices(i);
                distVect(dummind) = dists(i);
                dummind = dummind + 1;
            end
    end

    if length(workAngle) > 0
        disp("Moving through...");
        newInd = [ newInd, currPoint, workAngle(:,2)'];
        % minAngleInd = find( workAngle(:,1) == min(workAngle(:,1)) );
        currPoint = workAngle(find( workAngle(:,1) == min(workAngle(:,1)) ),2);
    else
        disp("No good points to follow");
        break
    end
    if ismember(point2ind,workAngle)
        initVar = false;
        disp("Reached the end... connected all points")
    end
end

stimulusRadius = mean(distVect)*1.2
stimSurfaceIndices = [];

for i = 1:length(indexTimeStim)
    [indices,dists] = findNearestNeighbors(ptCloudSurface,[ electrodeXInterest(indexTimeStim(i)), electrodeYInterest(indexTimeStim(i)), electrodeZInterest(indexTimeStim(i)) ],1);
    [moreInd, moreDists] =  findNeighborsInRadius( ptCloudSurface,[ surfaceXInterest(indices), surfaceYInterest(indices), surfaceZInterest(indices) ], stimulusRadius);
    stimSurfaceIndices = [ stimSurfaceIndices; reshape(moreInd, [length(moreInd),1]) ];
end

newInd = [newInd, stimSurfaceIndices'];
newInd = unique(newInd);


figure
plot3(surfaceXInterest, surfaceYInterest, surfaceZInterest, 'ro','DisplayName','tissue')
hold on
plot3(surfaceXInterest(newInd), surfaceYInterest(newInd), surfaceZInterest(newInd), 'yo','DisplayName','Stimulus')
% scatter3(electrodeXInterest,electrodeYInterest,electrodeZInterest,80,actTimesInterest,'filled')
colorbar
% caxis([40, 130]);
title('Patient Data')

size(newInd)

dlmwrite('../SimResults/P1simsStuff/P1_Paroxysmal3_RealGeometry_StimulusPoints.csv', round([ surfaceXInterest(newInd), surfaceYInterest(newInd), surfaceZInterest(newInd) ],4), 'delimiter', ',', 'precision', '%.6f');



%END


%%


indexInterest = 200;
figure
plot((1+min(actTimes2)):(100+min(actTimes2)),simulatedParoxysmal2_downsampled(:,indexInterest),'DisplayName','Simulated','LineWidth',3)
hold on
plot(elec1(:,indexInterest),'DisplayName','Patient','LineWidth',3)
legend()
title("Attempt to match EGMs")
xlabel('time (ms)')
ylabel('potential (mV)')
set(gca,'FontSize',20)



figure
pause(5)
for i = 1:1:length(simulatedParoxysmal2_downsampled(:,1))
    %scatter3(gridX(extrind),gridY(extrind),gridZ(extrind),20,VeMono_multi(i,extrind))
    subplot(2,1,1)
    scatter3(x,y,z,80,simulatedParoxysmal2_downsampled(i,:),'filled')
%     scatter(gridX,gridY,20,VeMono_multi(i,:))
    colorbar
    caxis([min(min(simulatedParoxysmal2_downsampled)), max(max(simulatedParoxysmal2_downsampled))]);
    title([ 'Simulated: ' num2str(i*1) ' msec'])
    subplot(2,1,2)
    scatter3(x,y,z,80,elec1(i,:),'filled')
%     scatter(gridX,gridY,20,VeBi_multi(i,:))
    title([ 'Patient: ' num2str(i*1) ' msec'])
    colorbar
    caxis([min(min(elec1)), max(max(elec1))]);
    pause(.75)
end











%% CHOOSING WAVEFORM

indextry = 5;%5,106
figure;
subplot(2,1,1);
plot( signalsAllInBetween_Paroxysmal(:,indextry), "LineWidth", 4 );
title('Patient 2: Paroxysmal signal 1')
set(gca,'FontSize',20)
subplot(2,1,2);
plot( signalsAllNarrow_Paroxysmal(:,indextry), "LineWidth", 4 );
title('Patient 2: Paroxysmal signal 2')
set(gca,'FontSize',20)



patient2Location = [ AllxNarrow_Paroxysmal(5), AllyNarrow_Paroxysmal(5), AllzNarrow_Paroxysmal(5); AllxInBetween_Paroxysmal(5), AllyInBetween_Paroxysmal(5), AllzInBetween_Paroxysmal(5) ];
patient3Location = [ AllxNarrow_Paroxysmal(106), AllyNarrow_Paroxysmal(106), AllzNarrow_Paroxysmal(106); AllxInBetween_Paroxysmal(106), AllyInBetween_Paroxysmal(106), AllzInBetween_Paroxysmal(106) ];



threshDistToUse = .4;
ind_P2_diffPoint1 = find( AllPatientVecTotal_Paroxysmal == 2 & abs( AllxTotal_Paroxysmal - P2_diffpoint1(1,1) ) < threshDistToUse & abs( AllyTotal_Paroxysmal - P2_diffpoint1(1,2) ) < threshDistToUse & abs( AllzTotal_Paroxysmal - P2_diffpoint1(1,3) ) < threshDistToUse )
threshDistToUse = .6;
ind_P2_diffPoint2 = find( AllPatientVecTotal_Paroxysmal == 2 & abs( AllxTotal_Paroxysmal - P2_diffpoint2(1,1) ) < threshDistToUse & abs( AllyTotal_Paroxysmal - P2_diffpoint2(1,2) ) < threshDistToUse & abs( AllzTotal_Paroxysmal - P2_diffpoint2(1,3) ) < threshDistToUse )
threshDistToUse = .001;
ind_P2_diffPoint3 = find( AllPatientVecTotal_Paroxysmal == 2 & abs( AllxTotal_Paroxysmal - P2_diffpoint3(1,1) ) < threshDistToUse & abs( AllyTotal_Paroxysmal - P2_diffpoint3(1,2) ) < threshDistToUse & abs( AllzTotal_Paroxysmal - P2_diffpoint3(1,3) ) < threshDistToUse )
threshDistToUse = .001;
ind_P2_diffPoint4 = find( AllPatientVecTotal_Paroxysmal == 2 & abs( AllxTotal_Paroxysmal - P2_diffpoint4(1,1) ) < threshDistToUse & abs( AllyTotal_Paroxysmal - P2_diffpoint4(1,2) ) < threshDistToUse & abs( AllzTotal_Paroxysmal - P2_diffpoint4(1,3) ) < threshDistToUse )




threshDistToUse = .001;
ind_P3_diffPoint1 = find( AllPatientVecTotal_Paroxysmal == 3 & abs( AllxTotal_Paroxysmal - P3_diffpoint1(1,1) ) < threshDistToUse & abs( AllyTotal_Paroxysmal - P3_diffpoint1(1,2) ) < threshDistToUse & abs( AllzTotal_Paroxysmal - P3_diffpoint1(1,3) ) < threshDistToUse )
threshDistToUse = .001;
ind_P3_diffPoint2 = find( AllPatientVecTotal_Paroxysmal == 3 & abs( AllxTotal_Paroxysmal - P3_diffpoint2(1,1) ) < threshDistToUse & abs( AllyTotal_Paroxysmal - P3_diffpoint2(1,2) ) < threshDistToUse & abs( AllzTotal_Paroxysmal - P3_diffpoint2(1,3) ) < threshDistToUse )
threshDistToUse = .001;
ind_P3_diffPoint3 = find( AllPatientVecTotal_Paroxysmal == 3 & abs( AllxTotal_Paroxysmal - P3_diffpoint3(1,1) ) < threshDistToUse & abs( AllyTotal_Paroxysmal - P3_diffpoint3(1,2) ) < threshDistToUse & abs( AllzTotal_Paroxysmal - P3_diffpoint3(1,3) ) < threshDistToUse )
threshDistToUse = .001;
ind_P3_diffPoint4 = find( AllPatientVecTotal_Paroxysmal == 3 & abs( AllxTotal_Paroxysmal - P3_diffpoint4(1,1) ) < threshDistToUse & abs( AllyTotal_Paroxysmal - P3_diffpoint4(1,2) ) < threshDistToUse & abs( AllzTotal_Paroxysmal - P3_diffpoint4(1,3) ) < threshDistToUse )




%NEW SIGNALS AND LOCATIONS
figure
subplot(2,2,1)
plot3(patient2xSurfaceCm, patient2ySurfaceCm, patient2zSurfaceCm, 'ro','DisplayName','tissue')
hold on
plot3(patient2xCm, patient2yCm, patient2zCm, 'bo','DisplayName','electrodes')
plot3(AllxTotal_Paroxysmal(5), AllyTotal_Paroxysmal(5), AllzTotal_Paroxysmal(5),'ko','MarkerFaceColor','k','MarkerSize',20,'DisplayName','1')
plot3(AllxTotal_Paroxysmal(352), AllyTotal_Paroxysmal(352), AllzTotal_Paroxysmal(352),'square','MarkerFaceColor','k','MarkerSize',20,'DisplayName','2')
plot3(AllxTotal_Paroxysmal(87), AllyTotal_Paroxysmal(87), AllzTotal_Paroxysmal(87),'diamond','MarkerFaceColor','k','MarkerSize',20,'DisplayName','3')
plot3(AllxTotal_Paroxysmal(509), AllyTotal_Paroxysmal(509), AllzTotal_Paroxysmal(509),'kv','MarkerFaceColor','k','MarkerSize',20,'DisplayName','4')
title('Patient 2')
legend()

subplot(2,2,2)
plot3(patient3xSurfaceCm, patient3ySurfaceCm, patient3zSurfaceCm, 'ro','DisplayName','tissue')
hold on
plot3(patient3xCm, patient3yCm, patient3zCm, 'bo','DisplayName','electrodes')
plot3(AllxTotal_Paroxysmal(460), AllyTotal_Paroxysmal(460), AllzTotal_Paroxysmal(460),'ko','MarkerFaceColor','k','MarkerSize',20,'DisplayName','1')
plot3(AllxTotal_Paroxysmal(145), AllyTotal_Paroxysmal(145), AllzTotal_Paroxysmal(145),'square','MarkerFaceColor','k','MarkerSize',20,'DisplayName','2')
plot3(AllxTotal_Paroxysmal(106), AllyTotal_Paroxysmal(106), AllzTotal_Paroxysmal(106),'diamond','MarkerFaceColor','k','MarkerSize',20,'DisplayName','3')
plot3(AllxTotal_Paroxysmal(453), AllyTotal_Paroxysmal(453), AllzTotal_Paroxysmal(453),'kv','MarkerFaceColor','k','MarkerSize',20,'DisplayName','4')
title('Patient 3')
legend()

subplot(2,2,3)
plot(AllSignals_Paroxysmal(:,5),'DisplayName','1','LineWidth',3);
hold on;
plot(AllSignals_Paroxysmal(:,352),'DisplayName','2','LineWidth',3);
plot(AllSignals_Paroxysmal(:,87),'DisplayName','3','LineWidth',3);
plot(AllSignals_Paroxysmal(:,509),'DisplayName','4','LineWidth',3);
legend()
title("Patient 2")
xlabel('time (ms)')
ylabel('normalized voltage')

subplot(2,2,4)
plot(AllSignals_Paroxysmal(:,460),'DisplayName','1','LineWidth',3);
hold on;
plot(AllSignals_Paroxysmal(:,145),'DisplayName','2','LineWidth',3);
plot(AllSignals_Paroxysmal(:,106),'DisplayName','3','LineWidth',3);
plot(AllSignals_Paroxysmal(:,453),'DisplayName','4','LineWidth',3);
legend()
title("Patient 3")
xlabel('time (ms)')
ylabel('normalized voltage')





%SIMULATED VS PATIENT SIGNALS HERE
innerIndexH1A_1 = find( abs(patient2xCm - AllxTotal_Paroxysmal(5)) < .001 & abs(patient2yCm - AllyTotal_Paroxysmal(5)) < .001 & abs(patient2zCm - AllzTotal_Paroxysmal(5)) < .001 );
innerIndexH1A_2 = find( abs(patient2xCm - AllxTotal_Paroxysmal(352)) < .001 & abs(patient2yCm - AllyTotal_Paroxysmal(352)) < .001 & abs(patient2zCm - AllzTotal_Paroxysmal(352)) < .001 );
innerIndexH1A_3 = find( abs(patient2xCm - AllxTotal_Paroxysmal(87)) < .001 & abs(patient2yCm - AllyTotal_Paroxysmal(87)) < .001 & abs(patient2zCm - AllzTotal_Paroxysmal(87)) < .001 );
innerIndexH1A_4 = find( abs(patient2xCm - AllxTotal_Paroxysmal(509)) < .001 & abs(patient2yCm - AllyTotal_Paroxysmal(509)) < .001 & abs(patient2zCm - AllzTotal_Paroxysmal(509)) < .001 );

innerIndexH1B_1 = find( abs(patient3xCm - AllxTotal_Paroxysmal(460)) < .001 & abs(patient3yCm - AllyTotal_Paroxysmal(460)) < .001 & abs(patient3zCm - AllzTotal_Paroxysmal(460)) < .001 );
innerIndexH1B_2 = find( abs(patient3xCm - AllxTotal_Paroxysmal(145)) < .001 & abs(patient3yCm - AllyTotal_Paroxysmal(145)) < .001 & abs(patient3zCm - AllzTotal_Paroxysmal(145)) < .001 );
innerIndexH1B_3 = find( abs(patient3xCm - AllxTotal_Paroxysmal(106)) < .001 & abs(patient3yCm - AllyTotal_Paroxysmal(106)) < .001 & abs(patient3zCm - AllzTotal_Paroxysmal(106)) < .001 );
innerIndexH1B_4 = find( abs(patient3xCm - AllxTotal_Paroxysmal(453)) < .001 & abs(patient3yCm - AllyTotal_Paroxysmal(453)) < .001 & abs(patient3zCm - AllzTotal_Paroxysmal(453)) < .001 );






figure
% subplot(2,2,1)
plot(Full_AllSignalsOG_Paroxysmal(:,5), '-', 'DisplayName','Patient','LineWidth',5);
hold on
plot(simulatedParoxysmal2_downsampled(:,innerIndexH1A_1), ':', 'DisplayName','Simulated','LineWidth',5)
title("H1A sample signal 4")
xlabel('Time (ms)')
ylabel('Potential (mV)')
legend()
set(gca,'FontSize',35)

figure
% subplot(2,2,2)
plot(Full_AllSignalsOG_Paroxysmal(:,352), '-', 'DisplayName','Patient','LineWidth',5);
hold on
plot(simulatedParoxysmal2_downsampled(:,innerIndexH1A_2), ':','DisplayName','Simulated','LineWidth',5)
title("H1A sample signal 1")
xlabel('Time (ms)')
ylabel('Potential (mV)')
legend()
set(gca,'FontSize',35)

figure
% subplot(2,2,3)
plot(Full_AllSignalsOG_Paroxysmal(:,87),'DisplayName','3 patient','LineWidth',3);
hold on
plot(simulatedParoxysmal2_downsampled(:,innerIndexH1A_3),'DisplayName','3 simulated','LineWidth',3)
title("Patient 2 (H1A): signal 3")
xlabel('time (ms)')
ylabel('potential (mV)')
legend()

figure
% subplot(2,2,4)
plot(Full_AllSignalsOG_Paroxysmal(:,509), '-', 'DisplayName','Patient','LineWidth',5);
hold on
plot(simulatedParoxysmal2_downsampled(:,innerIndexH1A_4), ':','DisplayName','Simulated','LineWidth',5)
title("H1A sample signal 2")
xlabel('Time (ms)')
ylabel('Potential (mV)')
legend()
set(gca,'FontSize',35)




figure
% subplot(2,2,1)
plot(Full_AllSignalsOG_Paroxysmal(:,460),'DisplayName','1 patient','LineWidth',3);
hold on
plot(simulatedParoxysmal3_downsampled(:,innerIndexH1B_1),'DisplayName','1 simulated','LineWidth',3)
title("Patient 3 (H1B): signal 1")
xlabel('time (ms)')
ylabel('potential (mV)')
legend()

figure
% subplot(2,2,2)
plot(Full_AllSignalsOG_Paroxysmal(:,145),'DisplayName','2 patient','LineWidth',3);
hold on
plot(simulatedParoxysmal3_downsampled(:,innerIndexH1B_2),'DisplayName','2 simulated','LineWidth',3)
title("Patient 3 (H1B): signal 2")
xlabel('time (ms)')
ylabel('potential (mV)')
legend()

figure
% subplot(2,2,3)
plot(Full_AllSignalsOG_Paroxysmal(:,106), '-', 'DisplayName','Patient','LineWidth',5);
hold on
plot(simulatedParoxysmal3_downsampled(:,innerIndexH1B_3), ':','DisplayName','Simulated','LineWidth',5)
title("H1B sample signal 1")
xlabel('Time (ms)')
ylabel('Potential (mV)')
legend()
set(gca,'FontSize',35)

figure
% subplot(2,2,4)
plot(Full_AllSignalsOG_Paroxysmal(:,453), '-', 'DisplayName','Patient','LineWidth',5);
hold on
plot(simulatedParoxysmal3_downsampled(:,innerIndexH1B_4), ':','DisplayName','Simulated','LineWidth',5)
title("H1B sample signal 2")
xlabel('Time (ms)')
ylabel('Potential (mV)')
legend()
set(gca,'FontSize',35)












H1A_SignalCutOG_1 = AllSignalsOG_Paroxysmal(:,5);
H1A_SignalCutOG_2 = AllSignalsOG_Paroxysmal(:,352);
H1A_SignalCutOG_3 = AllSignalsOG_Paroxysmal(:,87);
H1A_SignalCutOG_4 = AllSignalsOG_Paroxysmal(:,509);

H1B_SignalCutOG_1 = AllSignalsOG_Paroxysmal(:,460);
H1B_SignalCutOG_2 = AllSignalsOG_Paroxysmal(:,145);
H1B_SignalCutOG_3 = AllSignalsOG_Paroxysmal(:,106);
H1B_SignalCutOG_4 = AllSignalsOG_Paroxysmal(:,453);

H1A_SignalOG_1 = Full_AllSignalsOG_Paroxysmal(:,5);
H1A_SignalOG_2 = Full_AllSignalsOG_Paroxysmal(:,352);
H1A_SignalOG_3 = Full_AllSignalsOG_Paroxysmal(:,87);
H1A_SignalOG_4 = Full_AllSignalsOG_Paroxysmal(:,509);

H1B_SignalOG_1 = Full_AllSignalsOG_Paroxysmal(:,460);
H1B_SignalOG_2 = Full_AllSignalsOG_Paroxysmal(:,145);
H1B_SignalOG_3 = Full_AllSignalsOG_Paroxysmal(:,106);
H1B_SignalOG_4 = Full_AllSignalsOG_Paroxysmal(:,453);

H1A_locations = [ AllxTotal_Paroxysmal(5), AllyTotal_Paroxysmal(5), AllzTotal_Paroxysmal(5);...
                  AllxTotal_Paroxysmal(352), AllyTotal_Paroxysmal(352), AllzTotal_Paroxysmal(352);...
                  AllxTotal_Paroxysmal(87), AllyTotal_Paroxysmal(87), AllzTotal_Paroxysmal(87);...
                  AllxTotal_Paroxysmal(509), AllyTotal_Paroxysmal(509), AllzTotal_Paroxysmal(509);...
                  AllxTotal_Paroxysmal(5), AllyTotal_Paroxysmal(5), AllzTotal_Paroxysmal(5)];

H1B_locations = [ AllxTotal_Paroxysmal(460), AllyTotal_Paroxysmal(460), AllzTotal_Paroxysmal(460);...
                  AllxTotal_Paroxysmal(145), AllyTotal_Paroxysmal(145), AllzTotal_Paroxysmal(145);...
                  AllxTotal_Paroxysmal(106), AllyTotal_Paroxysmal(106), AllzTotal_Paroxysmal(106);...
                  AllxTotal_Paroxysmal(453), AllyTotal_Paroxysmal(453), AllzTotal_Paroxysmal(453);...
                  AllxTotal_Paroxysmal(460), AllyTotal_Paroxysmal(460), AllzTotal_Paroxysmal(460)];




% 
% save("signalsNormal.mat",'H1A_SignalCutOG_1','H1A_SignalCutOG_2','H1A_SignalCutOG_3','H1A_SignalCutOG_4',...
%     'H1A_SignalOG_1','H1A_SignalOG_2','H1A_SignalOG_3','H1A_SignalOG_4',...
%     'H1B_SignalCutOG_1','H1B_SignalCutOG_2','H1B_SignalCutOG_3','H1B_SignalCutOG_4',...
%     'H1B_SignalOG_1','H1B_SignalOG_2','H1B_SignalOG_3','H1B_SignalOG_4','H1A_locations','H1B_locations')
% 

% 
% dlmwrite('../SimResults/P1simsStuff/P1_Paroxysmal2_RealGeometryElecPoints_Reduced.csv', round(H1A_locations,4), 'delimiter', ',', 'precision', '%.6f');
% 


%%


pers7ids = find( patientFractionated_AllPatients_Persistent == 7 & deflectionCountFractionated_AllPatients_Persistent > 10 );
pers9ids = find( patientFractionated_AllPatients_Persistent == 9 & deflectionCountFractionated_AllPatients_Persistent > 10 );
pers6ids = find( patientFractionated_AllPatients_Persistent == 6 & deflectionCountFractionated_AllPatients_Persistent > 8 );
pers8ids = find( patientFractionated_AllPatients_Persistent == 8 & deflectionCountFractionated_AllPatients_Persistent > 10 );


interint = 38;
threshDistToUse = .001;
ind_P6_diffPoint1 = find( abs( patient6xCm - xFractionated_AllPatients_Persistent(pers6ids(interint)) ) < threshDistToUse & abs( patient6yCm - yFractionated_AllPatients_Persistent(pers6ids(interint)) ) < threshDistToUse & abs( patient6zCm - zFractionated_AllPatients_Persistent(pers6ids(interint)) ) < threshDistToUse )
ind_P6_diffPoint1_All = pers6ids(interint);
threshDistToUse = .001;
ind_P6_diffPoint2 = find( abs( patient6xCm - xFractionated_AllPatients_Persistent(pers6ids(interint+10)) ) < threshDistToUse & abs( patient6yCm - yFractionated_AllPatients_Persistent(pers6ids(interint+10)) ) < threshDistToUse & abs( patient6zCm - zFractionated_AllPatients_Persistent(pers6ids(interint+10)) ) < threshDistToUse )
ind_P6_diffPoint2_All = pers6ids(interint+10);
threshDistToUse = .001;
ind_P6_diffPoint3 = find( abs( patient6xCm - xFractionated_AllPatients_Persistent(pers6ids(interint-17)) ) < threshDistToUse & abs( patient6yCm - yFractionated_AllPatients_Persistent(pers6ids(interint-17)) ) < threshDistToUse & abs( patient6zCm - zFractionated_AllPatients_Persistent(pers6ids(interint-17)) ) < threshDistToUse )
ind_P6_diffPoint3_All = pers6ids(interint-17);
threshDistToUse = .001;
ind_P6_diffPoint4 = find( abs( patient6xCm - xFractionated_AllPatients_Persistent(pers6ids(interint-10)) ) < threshDistToUse & abs( patient6yCm - yFractionated_AllPatients_Persistent(pers6ids(interint-10)) ) < threshDistToUse & abs( patient6zCm - zFractionated_AllPatients_Persistent(pers6ids(interint-10)) ) < threshDistToUse )
ind_P6_diffPoint4_All = pers6ids(interint-10);%15

%one is for big vector and the other one for small individual vectors
interint = 22;
threshDistToUse = .001;
ind_P7_diffPoint1 = find( abs( patient7xCm - xFractionated_AllPatients_Persistent(pers7ids(interint)) ) < threshDistToUse & abs( patient7yCm - yFractionated_AllPatients_Persistent(pers7ids(interint)) ) < threshDistToUse & abs( patient7zCm - zFractionated_AllPatients_Persistent(pers7ids(interint)) ) < threshDistToUse )
ind_P7_diffPoint1_All = pers7ids(interint);
threshDistToUse = .001;
ind_P7_diffPoint2 = find( abs( patient7xCm - xFractionated_AllPatients_Persistent(pers7ids(interint+8)) ) < threshDistToUse & abs( patient7yCm - yFractionated_AllPatients_Persistent(pers7ids(interint+8)) ) < threshDistToUse & abs( patient7zCm - zFractionated_AllPatients_Persistent(pers7ids(interint+8)) ) < threshDistToUse )
ind_P7_diffPoint2_All = pers7ids(interint+8);
threshDistToUse = .001;
ind_P7_diffPoint3 = find( abs( patient7xCm - xFractionated_AllPatients_Persistent(pers7ids(interint-20)) ) < threshDistToUse & abs( patient7yCm - yFractionated_AllPatients_Persistent(pers7ids(interint-20)) ) < threshDistToUse & abs( patient7zCm - zFractionated_AllPatients_Persistent(pers7ids(interint-20)) ) < threshDistToUse )
ind_P7_diffPoint3_All = pers7ids(interint-20);
threshDistToUse = .001;
ind_P7_diffPoint4 = find( abs( patient7xCm - xFractionated_AllPatients_Persistent(pers7ids(interint-10)) ) < threshDistToUse & abs( patient7yCm - yFractionated_AllPatients_Persistent(pers7ids(interint-10)) ) < threshDistToUse & abs( patient7zCm - zFractionated_AllPatients_Persistent(pers7ids(interint-10)) ) < threshDistToUse )
ind_P7_diffPoint4_All = pers7ids(interint-10);


interint = 4;
threshDistToUse = .001;
ind_P8_diffPoint1 = find( abs( patient8xCm - xFractionated_AllPatients_Persistent(pers8ids(interint)) ) < threshDistToUse & abs( patient8yCm - yFractionated_AllPatients_Persistent(pers8ids(interint)) ) < threshDistToUse & abs( patient8zCm - zFractionated_AllPatients_Persistent(pers8ids(interint)) ) < threshDistToUse )
ind_P8_diffPoint1_All = pers8ids(interint);
threshDistToUse = .001;
ind_P8_diffPoint2 = find( abs( patient8xCm - xFractionated_AllPatients_Persistent(pers8ids(interint+1)) ) < threshDistToUse & abs( patient8yCm - yFractionated_AllPatients_Persistent(pers8ids(interint+1)) ) < threshDistToUse & abs( patient8zCm - zFractionated_AllPatients_Persistent(pers8ids(interint+1)) ) < threshDistToUse )
ind_P8_diffPoint2_All = pers8ids(interint+1);
threshDistToUse = .001;
ind_P8_diffPoint3 = find( abs( patient8xCm - xFractionated_AllPatients_Persistent(pers8ids(interint-2)) ) < threshDistToUse & abs( patient8yCm - yFractionated_AllPatients_Persistent(pers8ids(interint-2)) ) < threshDistToUse & abs( patient8zCm - zFractionated_AllPatients_Persistent(pers8ids(interint-2)) ) < threshDistToUse )
ind_P8_diffPoint3_All = pers8ids(interint-2);
threshDistToUse = .001;
ind_P8_diffPoint4 = find( abs( patient8xCm - xFractionated_AllPatients_Persistent(pers8ids(interint-1)) ) < threshDistToUse & abs( patient8yCm - yFractionated_AllPatients_Persistent(pers8ids(interint-1)) ) < threshDistToUse & abs( patient8zCm - zFractionated_AllPatients_Persistent(pers8ids(interint-1)) ) < threshDistToUse )
ind_P8_diffPoint4_All = pers8ids(interint-1);


interint = 26;
threshDistToUse = .001;
ind_P9_diffPoint1 = find( abs( patient9xCm - xFractionated_AllPatients_Persistent(pers9ids(interint)) ) < threshDistToUse & abs( patient9yCm - yFractionated_AllPatients_Persistent(pers9ids(interint)) ) < threshDistToUse & abs( patient9zCm - zFractionated_AllPatients_Persistent(pers9ids(interint)) ) < threshDistToUse )
ind_P9_diffPoint1_All = pers9ids(interint);
threshDistToUse = .001;
ind_P9_diffPoint2 = find( abs( patient9xCm - xFractionated_AllPatients_Persistent(pers9ids(interint+10)) ) < threshDistToUse & abs( patient9yCm - yFractionated_AllPatients_Persistent(pers9ids(interint+10)) ) < threshDistToUse & abs( patient9zCm - zFractionated_AllPatients_Persistent(pers9ids(interint+10)) ) < threshDistToUse )
ind_P9_diffPoint2_All = pers9ids(interint+10);
threshDistToUse = .001;
ind_P9_diffPoint3 = find( abs( patient9xCm - xFractionated_AllPatients_Persistent(pers9ids(interint-20)) ) < threshDistToUse & abs( patient9yCm - yFractionated_AllPatients_Persistent(pers9ids(interint-20)) ) < threshDistToUse & abs( patient9zCm - zFractionated_AllPatients_Persistent(pers9ids(interint-20)) ) < threshDistToUse )
ind_P9_diffPoint3_All = pers9ids(interint-20);
threshDistToUse = .001;
ind_P9_diffPoint4 = find( abs( patient9xCm - xFractionated_AllPatients_Persistent(pers9ids(interint-10)) ) < threshDistToUse & abs( patient9yCm - yFractionated_AllPatients_Persistent(pers9ids(interint-10)) ) < threshDistToUse & abs( patient9zCm - zFractionated_AllPatients_Persistent(pers9ids(interint-10)) ) < threshDistToUse )
ind_P9_diffPoint4_All = pers9ids(interint-10);


%PERSISTENT 6
figure
subplot(2,2,1)
plot( signalsFractionatedOG_AllPatients_Persistent((ind_P6_diffPoint1_All),:)', "LineWidth", 5 );
title('Target 1') 
xlabel('time (ms)')
ylabel('Potential (mV)')
set(gca,'FontSize',30)
subplot(2,2,2)
plot( signalsFractionatedOG_AllPatients_Persistent((ind_P6_diffPoint2_All),:)', "LineWidth", 5 );
title('Target 2') 
xlabel('time (ms)')
ylabel('Potential (mV)')
set(gca,'FontSize',30)
subplot(2,2,3)
plot( signalsFractionatedOG_AllPatients_Persistent((ind_P6_diffPoint3_All),:)', "LineWidth", 5 );
title('Target 3') 
xlabel('time (ms)')
ylabel('Potential (mV)')
set(gca,'FontSize',30)
subplot(2,2,4)
plot( signalsFractionatedOG_AllPatients_Persistent((ind_P6_diffPoint4_All),:)', "LineWidth", 5 );
title('Target 4') 
xlabel('time (ms)')
ylabel('Potential (mV)')
set(gca,'FontSize',30)
sgti = sgtitle('H2B Targets') ;
sgti.FontSize = 35;






%PERSISTENT 7
figure
subplot(2,2,1)
plot( signalsFractionatedOG_AllPatients_Persistent((ind_P7_diffPoint1_All),:)', "LineWidth", 5 );
title('Target 1') 
xlabel('time (ms)')
ylabel('Potential (mV)')
set(gca,'FontSize',30)
subplot(2,2,2)
plot( signalsFractionatedOG_AllPatients_Persistent((ind_P7_diffPoint2_All),:)', "LineWidth", 5 );
title('Target 2') 
xlabel('time (ms)')
ylabel('Potential (mV)')
set(gca,'FontSize',30)
subplot(2,2,3)
plot( signalsFractionatedOG_AllPatients_Persistent((ind_P7_diffPoint3_All),:)', "LineWidth", 5 );
title('Target 3') 
xlabel('time (ms)')
ylabel('Potential (mV)')
set(gca,'FontSize',30)
subplot(2,2,4)
plot( signalsFractionatedOG_AllPatients_Persistent((ind_P7_diffPoint4_All),:)', "LineWidth", 5 );
title('Target 4') 
xlabel('time (ms)')
ylabel('Potential (mV)')
set(gca,'FontSize',30)
sgti = sgtitle('H2C Targets') ;
sgti.FontSize = 35;





%PERSISTENT 8
figure
subplot(2,2,1)
plot( signalsFractionatedOG_AllPatients_Persistent((ind_P8_diffPoint1_All),:)', "LineWidth", 5 );
title('Target 1') 
xlabel('time (ms)')
ylabel('Potential (mV)')
set(gca,'FontSize',30)
subplot(2,2,2)
plot( signalsFractionatedOG_AllPatients_Persistent((ind_P8_diffPoint2_All),:)', '-', "LineWidth", 5 );
title('Target 2') 
xlabel('time (ms)')
ylabel('Potential (mV)')
set(gca,'FontSize',30)
subplot(2,2,3)
plot( signalsFractionatedOG_AllPatients_Persistent((ind_P8_diffPoint3_All),:)', "LineWidth", 5 );
title('Target 3') 
xlabel('time (ms)')
ylabel('Potential (mV)')
set(gca,'FontSize',30)
subplot(2,2,4)
plot( signalsFractionatedOG_AllPatients_Persistent((ind_P8_diffPoint4_All),:)', "LineWidth", 5 );
title('Target 4') 
xlabel('time (ms)')
ylabel('Potential (mV)')
set(gca,'FontSize',30)
sgti = sgtitle('H2D Targets') ;
sgti.FontSize = 35;




%PERSISTENT 9
figure
subplot(2,2,1)
plot( signalsFractionatedOG_AllPatients_Persistent((ind_P9_diffPoint1_All),:)', "LineWidth", 5 );
title('Target 1') 
xlabel('time (ms)')
ylabel('Potential (mV)')
set(gca,'FontSize',30)
subplot(2,2,2)
plot( signalsFractionatedOG_AllPatients_Persistent((ind_P9_diffPoint2_All),:)', "LineWidth", 5 );
title('Target 2') 
xlabel('time (ms)')
ylabel('Potential (mV)')
set(gca,'FontSize',30)
subplot(2,2,3)
plot( signalsFractionatedOG_AllPatients_Persistent((ind_P9_diffPoint3_All),:)', "LineWidth", 5 );
title('Target 3') 
xlabel('time (ms)')
ylabel('Potential (mV)')
set(gca,'FontSize',30)
subplot(2,2,4)
plot( signalsFractionatedOG_AllPatients_Persistent((ind_P9_diffPoint4_All),:)', "LineWidth", 5 );
title('Target 4') 
xlabel('time (ms)')
ylabel('Potential (mV)')
set(gca,'FontSize',30)
sgti = sgtitle('H2E Targets') ;
sgti.FontSize = 35;








interint = 48;
figure;
% subplot(2,2,1);
plot( signalsFractionated_AllPatients_Persistent(pers6ids(48),:)', "LineWidth", 4 );
title('Target signal 1: Persistent AFib') 
set(gca,'FontSize',20)
figure
% subplot(2,2,2);
plot( signalsFractionated_AllPatients_Persistent(pers7ids(2),:)', "LineWidth", 4 );
title('Target signal 2: Persistent AFib')
set(gca,'FontSize',20)
figure
% subplot(2,2,3);
plot( signalsFractionated_AllPatients_Persistent(pers8ids(2),:)', "LineWidth", 4 );
title('Target signal 3: Persistent AFib') 
set(gca,'FontSize',20)
figure
% subplot(2,2,4);
plot( signalsFractionated_AllPatients_Persistent(pers9ids(6),:)', "LineWidth", 4 );
title('Target signal 4: Persistent AFib')
set(gca,'FontSize',20)














%% PERSISTENT PATIENTS TO CHECK: COMPACT FIBROSIS

%CHECKING THINGS


%Simulated signals
simulatedPersistent_Compact15_6 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent6_CompactFibrosisNT15_2000points.csv");
simulatedPersistent_Compact15_6_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Compact15_6', .1, 5, 200, 1000, 10,false)';
actTimesSim_Compact15_6 = [];

simulatedPersistent_Compact15_7 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent7_CompactFibrosisNT15_2000points.csv");
simulatedPersistent_Compact15_7_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Compact15_7', .1, 5, 200, 1000, 10,false)';
actTimesSim_Compact15_7 = [];

simulatedPersistent_Compact15_6_downsampled_norm = [];
simulatedPersistent_Compact15_7_downsampled_norm = [];

simulatedPersistent_Compact15_8 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent8_CompactFibrosisNT15_2000points.csv");
simulatedPersistent_Compact15_8_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Compact15_8', .1, 5, 200, 1000, 10,false)';
actTimesSim_Compact15_8 = [];

simulatedPersistent_Compact15_9 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent9_CompactFibrosisNT15_2000points.csv");
simulatedPersistent_Compact15_9_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Compact15_9', .1, 5, 200, 1000, 10,false)';
actTimesSim_Compact15_9 = [];

simulatedPersistent_Compact15_8_downsampled_norm = [];
simulatedPersistent_Compact15_9_downsampled_norm = [];

% simulatedParoxysmal2_NF = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Paroxysmal2_faster_NoFibers_1000points.csv");
% simulatedParoxysmal2_NF_downsampled = UpsampleAndFilterFunction(simulatedParoxysmal2_NF', .1, 5, 200, 1000, 10,false)';
% actTimesSim_NF = [];

for i = 1:length(patient6xCm)
    dVvec = diff(simulatedPersistent_Compact15_6_downsampled(:,i));
    minVal =  min(dVvec);
    indFind = find(dVvec == minVal);
    indFind = indFind(1);
    actTimesSim_Compact15_6(i) = indFind; 

    simulatedPersistent_Compact15_6_downsampled_norm(:,i) = ( simulatedPersistent_Compact15_6_downsampled(:,i) - min(simulatedPersistent_Compact15_6_downsampled(:,i)) )./( max(simulatedPersistent_Compact15_6_downsampled(:,i)) - min(simulatedPersistent_Compact15_6_downsampled(:,i)) );
end
for i = 1:length(patient7xCm)
    dVvec2 = diff(simulatedPersistent_Compact15_7_downsampled(:,i));
    minVal2 =  min(dVvec2);
    indFind2 = find(dVvec2 == minVal2);
    if length(indFind2) > 0
        indFind2 = indFind2(1);
    else
        indFind2 = 0;
    end
    actTimesSim_Compact15_7(i) = indFind2; 

    simulatedPersistent_Compact15_7_downsampled_norm(:,i) = ( simulatedPersistent_Compact15_7_downsampled(:,i) - min(simulatedPersistent_Compact15_7_downsampled(:,i)) )./( max(simulatedPersistent_Compact15_7_downsampled(:,i)) - min(simulatedPersistent_Compact15_7_downsampled(:,i)) );
end
for i = 1:length(patient8xCm)
    dVvec2 = diff(simulatedPersistent_Compact15_8_downsampled(:,i));
    minVal2 =  min(dVvec2);
    indFind2 = find(dVvec2 == minVal2);
    if length(indFind2) > 0
        indFind2 = indFind2(1);
    else
        indFind2 = 0;
    end
    actTimesSim_Compact15_8(i) = indFind2; 

    simulatedPersistent_Compact15_8_downsampled_norm(:,i) = ( simulatedPersistent_Compact15_8_downsampled(:,i) - min(simulatedPersistent_Compact15_8_downsampled(:,i)) )./( max(simulatedPersistent_Compact15_8_downsampled(:,i)) - min(simulatedPersistent_Compact15_8_downsampled(:,i)) );
end
for i = 1:length(patient9xCm)
    dVvec2 = diff(simulatedPersistent_Compact15_9_downsampled(:,i));
    minVal2 =  min(dVvec2);
    indFind2 = find(dVvec2 == minVal2);
    if length(indFind2) > 0
        indFind2 = indFind2(1);
    else
        indFind2 = 0;
    end
    actTimesSim_Compact15_9(i) = indFind2; 

    simulatedPersistent_Compact15_9_downsampled_norm(:,i) = ( simulatedPersistent_Compact15_9_downsampled(:,i) - min(simulatedPersistent_Compact15_9_downsampled(:,i)) )./( max(simulatedPersistent_Compact15_9_downsampled(:,i)) - min(simulatedPersistent_Compact15_9_downsampled(:,i)) );
end





%Simulated signals
simulatedPersistent_Compact25_6 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent6_CompactFibrosisNT25_2000points.csv");
simulatedPersistent_Compact25_6_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Compact25_6', .1, 5, 200, 1000, 10,false)';
actTimesSim_Compact25_6 = [];

simulatedPersistent_Compact25_7 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent7_CompactFibrosisNT25_2000points.csv");
simulatedPersistent_Compact25_7_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Compact25_7', .1, 5, 200, 1000, 10,false)';
actTimesSim_Compact25_7 = [];

simulatedPersistent_Compact25_6_downsampled_norm = [];
simulatedPersistent_Compact25_7_downsampled_norm = [];

simulatedPersistent_Compact25_8 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent8_CompactFibrosisNT25_2000points.csv");
simulatedPersistent_Compact25_8_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Compact25_8', .1, 5, 200, 1000, 10,false)';
actTimesSim_Compact25_8 = [];

simulatedPersistent_Compact25_9 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent9_CompactFibrosisNT25_2000points.csv");
simulatedPersistent_Compact25_9_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Compact25_9', .1, 5, 200, 1000, 10,false)';
actTimesSim_Compact25_9 = [];

simulatedPersistent_Compact25_8_downsampled_norm = [];
simulatedPersistent_Compact25_9_downsampled_norm = [];

% simulatedParoxysmal2_NF = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Paroxysmal2_faster_NoFibers_1000points.csv");
% simulatedParoxysmal2_NF_downsampled = UpsampleAndFilterFunction(simulatedParoxysmal2_NF', .1, 5, 200, 1000, 10,false)';
% actTimesSim_NF = [];

for i = 1:length(patient6xCm)
    dVvec = diff(simulatedPersistent_Compact25_6_downsampled(:,i));
    minVal =  min(dVvec);
    indFind = find(dVvec == minVal);
    indFind = indFind(1);
    actTimesSim_Compact25_6(i) = indFind; 

    simulatedPersistent_Compact25_6_downsampled_norm(:,i) = ( simulatedPersistent_Compact25_6_downsampled(:,i) - min(simulatedPersistent_Compact25_6_downsampled(:,i)) )./( max(simulatedPersistent_Compact25_6_downsampled(:,i)) - min(simulatedPersistent_Compact25_6_downsampled(:,i)) );
end
for i = 1:length(patient7xCm)
    dVvec2 = diff(simulatedPersistent_Compact25_7_downsampled(:,i));
    minVal2 =  min(dVvec2);
    indFind2 = find(dVvec2 == minVal2);
    if length(indFind2) > 0
        indFind2 = indFind2(1);
    else
        indFind2 = 0;
    end
    actTimesSim_Compact25_7(i) = indFind2; 

    simulatedPersistent_Compact25_7_downsampled_norm(:,i) = ( simulatedPersistent_Compact25_7_downsampled(:,i) - min(simulatedPersistent_Compact25_7_downsampled(:,i)) )./( max(simulatedPersistent_Compact25_7_downsampled(:,i)) - min(simulatedPersistent_Compact25_7_downsampled(:,i)) );
end
for i = 1:length(patient8xCm)
    dVvec2 = diff(simulatedPersistent_Compact25_8_downsampled(:,i));
    minVal2 =  min(dVvec2);
    indFind2 = find(dVvec2 == minVal2);
    if length(indFind2) > 0
        indFind2 = indFind2(1);
    else
        indFind2 = 0;
    end
    actTimesSim_Compact25_8(i) = indFind2; 

    simulatedPersistent_Compact25_8_downsampled_norm(:,i) = ( simulatedPersistent_Compact25_8_downsampled(:,i) - min(simulatedPersistent_Compact25_8_downsampled(:,i)) )./( max(simulatedPersistent_Compact25_8_downsampled(:,i)) - min(simulatedPersistent_Compact25_8_downsampled(:,i)) );
end
for i = 1:length(patient9xCm)
    dVvec2 = diff(simulatedPersistent_Compact25_9_downsampled(:,i));
    minVal2 =  min(dVvec2);
    indFind2 = find(dVvec2 == minVal2);
    if length(indFind2) > 0
        indFind2 = indFind2(1);
    else
        indFind2 = 0;
    end
    actTimesSim_Compact25_9(i) = indFind2; 

    simulatedPersistent_Compact25_9_downsampled_norm(:,i) = ( simulatedPersistent_Compact25_9_downsampled(:,i) - min(simulatedPersistent_Compact25_9_downsampled(:,i)) )./( max(simulatedPersistent_Compact25_9_downsampled(:,i)) - min(simulatedPersistent_Compact25_9_downsampled(:,i)) );
end


actTimes6 = expDB.dataSet{6}.map.continuousMapping.revisedLat(:,1);
actTimes7 = expDB.dataSet{7}.map.continuousMapping.revisedLat(:,1);
actTimes8 = expDB.dataSet{8}.map.continuousMapping.revisedLat(:,1);
actTimes9 = expDB.dataSet{9}.map.continuousMapping.revisedLat(:,1);

figure
subplot(2,2,1)
plot3(patient6xSurfaceCm, patient6ySurfaceCm, patient6zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient6xCm,patient6yCm,patient6zCm,80,actTimes6,'filled')
colorbar
caxis([40, 140]);
title('H2B Patient')
set(gca,'FontSize',25)
subplot(2,2,2)
plot3(patient6xSurfaceCm, patient6ySurfaceCm, patient6zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient6xCm,patient6yCm,patient6zCm,80,actTimesSim_Compact15_6,'filled')
colorbar
caxis([40, 140]);
title('H2B Simulation')
set(gca,'FontSize',25)
subplot(2,2,3)
plot3(patient7xSurfaceCm, patient7ySurfaceCm, patient7zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient7xCm,patient7yCm,patient7zCm,80,actTimes7,'filled')
colorbar
caxis([50, 160]);
title('H2C Patient')
set(gca,'FontSize',25)
subplot(2,2,4)
plot3(patient7xSurfaceCm, patient7ySurfaceCm, patient7zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient7xCm,patient7yCm,patient7zCm,80,actTimesSim_Compact15_7,'filled')
colorbar
caxis([50, 160]);
title('H2C Simulation')
set(gca,'FontSize',25)





figure
subplot(2,2,1)
plot3(patient8xSurfaceCm, patient8ySurfaceCm, patient8zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient8xCm,patient8yCm,patient8zCm,80,actTimes8,'filled')
colorbar
caxis([30, 110]);
title('H2D Patient')
set(gca,'FontSize',25)
subplot(2,2,2)
plot3(patient8xSurfaceCm, patient8ySurfaceCm, patient8zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient8xCm,patient8yCm,patient8zCm,80,actTimesSim_Compact15_8,'filled')
colorbar
caxis([30, 110]);
title('H2D Simulation')
set(gca,'FontSize',25)
subplot(2,2,3)
plot3(patient9xSurfaceCm, patient9ySurfaceCm, patient9zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient9xCm,patient9yCm,patient9zCm,80,actTimes9,'filled')
colorbar
caxis([30, 130]);
title('H2E Patient')
set(gca,'FontSize',25)
subplot(2,2,4)
plot3(patient9xSurfaceCm, patient9ySurfaceCm, patient9zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient9xCm,patient9yCm,patient9zCm,80,actTimesSim_Compact15_9,'filled')
colorbar
caxis([30, 130]);
title('H2E Simulation')
set(gca,'FontSize',25)




%COMPACT FIBROSIS
figure;
subplot(2,2,1)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{6}.map.continuousMapping.fullEgms(ind_P6_diffPoint1,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Compact15_6_downsampled(:,ind_P6_diffPoint1), 'b', "LineWidth", 3 )
legend('Patient','Simulated')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 1")
set(gca,'FontSize',30)
subplot(2,2,2)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );

plot( squeeze(expDB.dataSet{6}.map.continuousMapping.fullEgms(ind_P6_diffPoint2,1,1:200)), '-', "LineWidth", 5 );
hold on;
plot( simulatedPersistent_Compact15_6_downsampled(:,ind_P6_diffPoint2), ':', "LineWidth", 5 )
legend('Patient','Simulated')
xlabel('Time (ms)')
ylabel('Potential (mV)')
title("H2B sample signal 1")
set(gca,'FontSize',35)

subplot(2,2,3)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{6}.map.continuousMapping.fullEgms(ind_P6_diffPoint3,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Compact15_6_downsampled(:,ind_P6_diffPoint3), 'b', "LineWidth", 3 )
legend('Patient','Simulated')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 3")
set(gca,'FontSize',30)
subplot(2,2,4)

% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{6}.map.continuousMapping.fullEgms(ind_P6_diffPoint4,1,1:200)), '-', "LineWidth", 5 );
hold on;
plot( simulatedPersistent_Compact15_6_downsampled(:,ind_P6_diffPoint4), ':', "LineWidth", 5 )
legend('Patient','Simulated')
xlabel('Time (ms)')
ylabel('Potential (mV)')
title("H2B sample signal 2")
set(gca,'FontSize',35)

sgti = sgtitle('H2B signals') ;
sgti.FontSize = 35;




%COMPACT FIBROSIS
figure;
subplot(2,2,1)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{7}.map.continuousMapping.fullEgms(ind_P7_diffPoint1,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Compact15_7_downsampled(:,ind_P7_diffPoint1), 'b', "LineWidth", 3 )
legend('Patient','Simulated')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 1")
set(gca,'FontSize',30)
subplot(2,2,2)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{7}.map.continuousMapping.fullEgms(ind_P7_diffPoint2,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Compact15_7_downsampled(:,ind_P7_diffPoint2), 'b', "LineWidth", 3 )
legend('Patient','Simulated')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 2")
set(gca,'FontSize',30)
subplot(2,2,3)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{7}.map.continuousMapping.fullEgms(ind_P7_diffPoint3,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Compact15_7_downsampled(:,ind_P7_diffPoint3), 'b', "LineWidth", 3 )
legend('Patient','Simulated')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 3")
set(gca,'FontSize',30)
subplot(2,2,4)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{7}.map.continuousMapping.fullEgms(ind_P7_diffPoint4,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Compact15_7_downsampled(:,ind_P7_diffPoint4), 'b', "LineWidth", 3 )
legend('Patient','Simulated')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 4")
set(gca,'FontSize',30)
sgti = sgtitle('H2C signals') ;
sgti.FontSize = 35;








%COMPACT FIBROSIS
figure;
subplot(2,2,1)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{8}.map.continuousMapping.fullEgms(ind_P8_diffPoint1,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Compact15_8_downsampled(:,ind_P8_diffPoint1), 'b', "LineWidth", 3 )
legend('Patient','Simulated')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 1")
set(gca,'FontSize',30)
subplot(2,2,2)

% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{8}.map.continuousMapping.fullEgms(ind_P8_diffPoint2,1,1:200)), '-', "LineWidth", 5 );
hold on;
plot( simulatedPersistent_Compact15_8_downsampled(:,ind_P8_diffPoint2), ':', "LineWidth", 5 )
legend('Patient','Simulated')
xlabel('Time (ms)')
ylabel('Potential (mV)')
title("H2D sample signal 1")
set(gca,'FontSize',35)


subplot(2,2,3)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{8}.map.continuousMapping.fullEgms(ind_P8_diffPoint3,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Compact15_8_downsampled(:,ind_P8_diffPoint3), 'b', "LineWidth", 3 )
legend('Patient','Simulated')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 3")
set(gca,'FontSize',30)
subplot(2,2,4)

% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{8}.map.continuousMapping.fullEgms(ind_P8_diffPoint4,1,1:200)), '-', "LineWidth", 5 );
hold on;
plot( simulatedPersistent_Compact15_8_downsampled(:,ind_P8_diffPoint4), ':', "LineWidth", 5 )
legend('Patient','Simulated')
xlabel('Time (ms)')
ylabel('Potential (mV)')
title("H2D sample signal 2")
set(gca,'FontSize',35)

sgti = sgtitle('H2D signals') ;
sgti.FontSize = 35;






%COMPACT FIBROSIS
figure;
subplot(2,2,1)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{9}.map.continuousMapping.fullEgms(ind_P9_diffPoint1,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Compact15_9_downsampled(:,ind_P9_diffPoint1), 'b', "LineWidth", 3 )
legend('Patient','Simulated')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 1")
set(gca,'FontSize',30)
subplot(2,2,2)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{9}.map.continuousMapping.fullEgms(ind_P9_diffPoint2,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Compact15_9_downsampled(:,ind_P9_diffPoint2), 'b', "LineWidth", 3 )
legend('Patient','Simulated')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 2")
set(gca,'FontSize',30)
subplot(2,2,3)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{9}.map.continuousMapping.fullEgms(ind_P9_diffPoint3,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Compact15_9_downsampled(:,ind_P9_diffPoint3), 'b', "LineWidth", 3 )
legend('Patient','Simulated')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 3")
set(gca,'FontSize',30)
subplot(2,2,4)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{9}.map.continuousMapping.fullEgms(ind_P9_diffPoint4,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Compact15_9_downsampled(:,ind_P9_diffPoint4), 'b', "LineWidth", 3 )
legend('Patient','Simulated')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 4")
set(gca,'FontSize',30)
sgti = sgtitle('H2E signals') ;
sgti.FontSize = 35;




%% PERSISTENT PATIENTS TO CHECK: DIFFUSE FIBROSIS

%CHECKING THINGS


%Simulated signals
simulatedPersistent_Diffuse15_6 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent6_DiffuseFibrosisNT15_2000points.csv");
simulatedPersistent_Diffuse15_6_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Diffuse15_6', .1, 5, 200, 1000, 10,false)';
actTimesSim_Diffuse15_6 = [];

simulatedPersistent_Diffuse15_7 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent7_DiffuseFibrosisNT15_2000points.csv");
simulatedPersistent_Diffuse15_7_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Diffuse15_7', .1, 5, 200, 1000, 10,false)';
actTimesSim_Diffuse15_7 = [];

simulatedPersistent_Diffuse15_6_downsampled_norm = [];
simulatedPersistent_Diffuse15_7_downsampled_norm = [];

simulatedPersistent_Diffuse15_8 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent8_DiffuseFibrosisNT15_2000points.csv");
simulatedPersistent_Diffuse15_8_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Diffuse15_8', .1, 5, 200, 1000, 10,false)';
actTimesSim_Diffuse15_8 = [];

simulatedPersistent_Diffuse15_9 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent9_DiffuseFibrosisNT15_2000points.csv");
simulatedPersistent_Diffuse15_9_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Diffuse15_9', .1, 5, 200, 1000, 10,false)';
actTimesSim_Diffuse15_9 = [];

simulatedPersistent_Diffuse15_8_downsampled_norm = [];
simulatedPersistent_Diffuse15_9_downsampled_norm = [];


for i = 1:length(patient6xCm)
    dVvec = diff(simulatedPersistent_Diffuse15_6_downsampled(:,i));
    minVal =  min(dVvec);
    indFind = find(dVvec == minVal);
    indFind = indFind(1);
    actTimesSim_Diffuse15_6(i) = indFind; 

    simulatedPersistent_Diffuse15_6_downsampled_norm(:,i) = ( simulatedPersistent_Diffuse15_6_downsampled(:,i) - min(simulatedPersistent_Diffuse15_6_downsampled(:,i)) )./( max(simulatedPersistent_Diffuse15_6_downsampled(:,i)) - min(simulatedPersistent_Diffuse15_6_downsampled(:,i)) );
end
for i = 1:length(patient7xCm)
    dVvec2 = diff(simulatedPersistent_Diffuse15_7_downsampled(:,i));
    minVal2 =  min(dVvec2);
    indFind2 = find(dVvec2 == minVal2);
    if length(indFind2) > 0
        indFind2 = indFind2(1);
    else
        indFind2 = 0;
    end
    actTimesSim_Diffuse15_7(i) = indFind2; 

    simulatedPersistent_Diffuse15_7_downsampled_norm(:,i) = ( simulatedPersistent_Diffuse15_7_downsampled(:,i) - min(simulatedPersistent_Diffuse15_7_downsampled(:,i)) )./( max(simulatedPersistent_Diffuse15_7_downsampled(:,i)) - min(simulatedPersistent_Diffuse15_7_downsampled(:,i)) );
end
for i = 1:length(patient8xCm)
    dVvec2 = diff(simulatedPersistent_Diffuse15_8_downsampled(:,i));
    minVal2 =  min(dVvec2);
    indFind2 = find(dVvec2 == minVal2);
    if length(indFind2) > 0
        indFind2 = indFind2(1);
    else
        indFind2 = 0;
    end
    actTimesSim_Diffuse15_8(i) = indFind2; 

    simulatedPersistent_Diffuse15_8_downsampled_norm(:,i) = ( simulatedPersistent_Diffuse15_8_downsampled(:,i) - min(simulatedPersistent_Diffuse15_8_downsampled(:,i)) )./( max(simulatedPersistent_Diffuse15_8_downsampled(:,i)) - min(simulatedPersistent_Diffuse15_8_downsampled(:,i)) );
end
for i = 1:length(patient9xCm)
    dVvec2 = diff(simulatedPersistent_Diffuse15_9_downsampled(:,i));
    minVal2 =  min(dVvec2);
    indFind2 = find(dVvec2 == minVal2);
    if length(indFind2) > 0
        indFind2 = indFind2(1);
    else
        indFind2 = 0;
    end
    actTimesSim_Diffuse15_9(i) = indFind2; 

    simulatedPersistent_Diffuse15_9_downsampled_norm(:,i) = ( simulatedPersistent_Diffuse15_9_downsampled(:,i) - min(simulatedPersistent_Diffuse15_9_downsampled(:,i)) )./( max(simulatedPersistent_Diffuse15_9_downsampled(:,i)) - min(simulatedPersistent_Diffuse15_9_downsampled(:,i)) );
end




%Simulated signals
simulatedPersistent_Diffuse25_6 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent6_DiffuseFibrosisNT25_2000points.csv");
simulatedPersistent_Diffuse25_6_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Diffuse25_6', .1, 5, 200, 1000, 10,false)';
actTimesSim_Diffuse25_6 = [];

simulatedPersistent_Diffuse25_7 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent7_DiffuseFibrosisNT25_2000points.csv");
simulatedPersistent_Diffuse25_7_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Diffuse25_7', .1, 5, 200, 1000, 10,false)';
actTimesSim_Diffuse25_7 = [];

simulatedPersistent_Diffuse25_6_downsampled_norm = [];
simulatedPersistent_Diffuse25_7_downsampled_norm = [];

simulatedPersistent_Diffuse25_8 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent8_DiffuseFibrosisNT25_2000points.csv");
simulatedPersistent_Diffuse25_8_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Diffuse25_8', .1, 5, 200, 1000, 10,false)';
actTimesSim_Diffuse25_8 = [];

simulatedPersistent_Diffuse25_9 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent9_DiffuseFibrosisNT25_2000points.csv");
simulatedPersistent_Diffuse25_9_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Diffuse25_9', .1, 5, 200, 1000, 10,false)';
actTimesSim_Diffuse25_9 = [];

simulatedPersistent_Diffuse25_8_downsampled_norm = [];
simulatedPersistent_Diffuse25_9_downsampled_norm = [];


for i = 1:length(patient6xCm)
    dVvec = diff(simulatedPersistent_Diffuse25_6_downsampled(:,i));
    minVal =  min(dVvec);
    indFind = find(dVvec == minVal);
    indFind = indFind(1);
    actTimesSim_Diffuse25_6(i) = indFind; 

    simulatedPersistent_Diffuse25_6_downsampled_norm(:,i) = ( simulatedPersistent_Diffuse25_6_downsampled(:,i) - min(simulatedPersistent_Diffuse25_6_downsampled(:,i)) )./( max(simulatedPersistent_Diffuse25_6_downsampled(:,i)) - min(simulatedPersistent_Diffuse25_6_downsampled(:,i)) );
end
for i = 1:length(patient7xCm)
    dVvec2 = diff(simulatedPersistent_Diffuse25_7_downsampled(:,i));
    minVal2 =  min(dVvec2);
    indFind2 = find(dVvec2 == minVal2);
    if length(indFind2) > 0
        indFind2 = indFind2(1);
    else
        indFind2 = 0;
    end
    actTimesSim_Diffuse25_7(i) = indFind2; 

    simulatedPersistent_Diffuse25_7_downsampled_norm(:,i) = ( simulatedPersistent_Diffuse25_7_downsampled(:,i) - min(simulatedPersistent_Diffuse25_7_downsampled(:,i)) )./( max(simulatedPersistent_Diffuse25_7_downsampled(:,i)) - min(simulatedPersistent_Diffuse25_7_downsampled(:,i)) );
end
for i = 1:length(patient8xCm)
    dVvec2 = diff(simulatedPersistent_Diffuse25_8_downsampled(:,i));
    minVal2 =  min(dVvec2);
    indFind2 = find(dVvec2 == minVal2);
    if length(indFind2) > 0
        indFind2 = indFind2(1);
    else
        indFind2 = 0;
    end
    actTimesSim_Diffuse25_8(i) = indFind2; 

    simulatedPersistent_Diffuse25_8_downsampled_norm(:,i) = ( simulatedPersistent_Diffuse25_8_downsampled(:,i) - min(simulatedPersistent_Diffuse25_8_downsampled(:,i)) )./( max(simulatedPersistent_Diffuse25_8_downsampled(:,i)) - min(simulatedPersistent_Diffuse25_8_downsampled(:,i)) );
end
for i = 1:length(patient9xCm)
    dVvec2 = diff(simulatedPersistent_Diffuse25_9_downsampled(:,i));
    minVal2 =  min(dVvec2);
    indFind2 = find(dVvec2 == minVal2);
    if length(indFind2) > 0
        indFind2 = indFind2(1);
    else
        indFind2 = 0;
    end
    actTimesSim_Diffuse25_9(i) = indFind2; 

    simulatedPersistent_Diffuse25_9_downsampled_norm(:,i) = ( simulatedPersistent_Diffuse25_9_downsampled(:,i) - min(simulatedPersistent_Diffuse25_9_downsampled(:,i)) )./( max(simulatedPersistent_Diffuse25_9_downsampled(:,i)) - min(simulatedPersistent_Diffuse25_9_downsampled(:,i)) );
end




figure
subplot(2,2,1)
plot3(patient6xSurfaceCm, patient6ySurfaceCm, patient6zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient6xCm,patient6yCm,patient6zCm,80,actTimes6,'filled')
colorbar
caxis([40, 140]);
title('Patient 6 Data')
subplot(2,2,2)
plot3(patient6xSurfaceCm, patient6ySurfaceCm, patient6zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient6xCm,patient6yCm,patient6zCm,80,actTimesSim_Diffuse15_6,'filled')
colorbar
caxis([40, 140]);
title('Simulation 6 Data')
subplot(2,2,3)
plot3(patient7xSurfaceCm, patient7ySurfaceCm, patient7zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient7xCm,patient7yCm,patient7zCm,80,actTimes7,'filled')
colorbar
caxis([50, 160]);
title('Patient 7 Data')
subplot(2,2,4)
plot3(patient7xSurfaceCm, patient7ySurfaceCm, patient7zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient7xCm,patient7yCm,patient7zCm,80,actTimesSim_Diffuse15_7,'filled')
colorbar
caxis([50, 160]);
title('Simulation 7 Data')






figure
subplot(2,2,1)
plot3(patient8xSurfaceCm, patient8ySurfaceCm, patient8zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient8xCm,patient8yCm,patient8zCm,80,actTimes8,'filled')
colorbar
caxis([30, 110]);
title('Patient 8 Data')
subplot(2,2,2)
plot3(patient8xSurfaceCm, patient8ySurfaceCm, patient8zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient8xCm,patient8yCm,patient8zCm,80,actTimesSim_Diffuse15_8,'filled')
colorbar
caxis([30, 110]);
title('Simulation 8 Data')
subplot(2,2,3)
plot3(patient9xSurfaceCm, patient9ySurfaceCm, patient9zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient9xCm,patient9yCm,patient9zCm,80,actTimes9,'filled')
colorbar
caxis([30, 130]);
title('Patient 9 Data')
subplot(2,2,4)
plot3(patient9xSurfaceCm, patient9ySurfaceCm, patient9zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient9xCm,patient9yCm,patient9zCm,80,actTimesSim_Diffuse15_9,'filled')
colorbar
caxis([30, 130]);
title('Simulation 9 Data')





%DIFFUSE FIBROSIS
figure;
subplot(2,2,1)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{6}.map.continuousMapping.fullEgms(ind_P6_diffPoint1,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Diffuse15_6_downsampled(:,ind_P6_diffPoint1), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 1")
subplot(2,2,2)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{6}.map.continuousMapping.fullEgms(ind_P6_diffPoint2,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Diffuse15_6_downsampled(:,ind_P6_diffPoint2), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 2")
subplot(2,2,3)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{6}.map.continuousMapping.fullEgms(ind_P6_diffPoint3,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Diffuse15_6_downsampled(:,ind_P6_diffPoint3), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 3")
subplot(2,2,4)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{6}.map.continuousMapping.fullEgms(ind_P6_diffPoint4,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Diffuse15_6_downsampled(:,ind_P6_diffPoint4), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 4")
sgtitle('H2B signals')





figure;
subplot(2,2,1)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{7}.map.continuousMapping.fullEgms(ind_P7_diffPoint1,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Diffuse15_7_downsampled(:,ind_P7_diffPoint1), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 1")
subplot(2,2,2)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{7}.map.continuousMapping.fullEgms(ind_P7_diffPoint2,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Diffuse15_7_downsampled(:,ind_P7_diffPoint2), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 2")
subplot(2,2,3)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{7}.map.continuousMapping.fullEgms(ind_P7_diffPoint3,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Diffuse15_7_downsampled(:,ind_P7_diffPoint3), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 3")
subplot(2,2,4)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{7}.map.continuousMapping.fullEgms(ind_P7_diffPoint4,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Diffuse15_7_downsampled(:,ind_P7_diffPoint4), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 4")
sgtitle('H2C signals')





figure;
subplot(2,2,1)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{8}.map.continuousMapping.fullEgms(ind_P8_diffPoint1,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Diffuse15_8_downsampled(:,ind_P8_diffPoint1), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 1")
subplot(2,2,2)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{8}.map.continuousMapping.fullEgms(ind_P8_diffPoint2,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Diffuse15_8_downsampled(:,ind_P8_diffPoint2), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 2")
subplot(2,2,3)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{8}.map.continuousMapping.fullEgms(ind_P8_diffPoint3,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Diffuse15_8_downsampled(:,ind_P8_diffPoint3), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 3")
subplot(2,2,4)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{8}.map.continuousMapping.fullEgms(ind_P8_diffPoint4,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Diffuse15_8_downsampled(:,ind_P8_diffPoint4), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 4")
sgtitle('H2D signals')






figure;
subplot(2,2,1)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{9}.map.continuousMapping.fullEgms(ind_P9_diffPoint1,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Diffuse15_9_downsampled(:,ind_P9_diffPoint1), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 1")
subplot(2,2,2)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{9}.map.continuousMapping.fullEgms(ind_P9_diffPoint2,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Diffuse15_9_downsampled(:,ind_P9_diffPoint2), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 2")
subplot(2,2,3)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{9}.map.continuousMapping.fullEgms(ind_P9_diffPoint3,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Diffuse15_9_downsampled(:,ind_P9_diffPoint3), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 3")
subplot(2,2,4)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{9}.map.continuousMapping.fullEgms(ind_P9_diffPoint4,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Diffuse15_9_downsampled(:,ind_P9_diffPoint4), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 4")
sgtitle('H2E signals')




%% PERSISTENT PATIENTS TO CHECK: INTERSTITIAL FIBROSIS



%CHECKING THINGS


%Simulated signals
simulatedPersistent_Interstitial15_6 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent6_InterstitialFibrosisNT15_2000points.csv");
simulatedPersistent_Interstitial15_6_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Interstitial15_6', .1, 5, 200, 1000, 10,false)';
actTimesSim_Interstitial15_6 = [];

simulatedPersistent_Interstitial15_7 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent7_InterstitialFibrosisNT15_2000points.csv");
simulatedPersistent_Interstitial15_7_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Interstitial15_7', .1, 5, 200, 1000, 10,false)';
actTimesSim_Interstitial15_7 = [];

simulatedPersistent_Interstitial15_6_downsampled_norm = [];
simulatedPersistent_Interstitial15_7_downsampled_norm = [];

simulatedPersistent_Interstitial15_8 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent8_InterstitialFibrosisNT15_2000points.csv");
simulatedPersistent_Interstitial15_8_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Interstitial15_8', .1, 5, 200, 1000, 10,false)';
actTimesSim_Interstitial15_8 = [];

simulatedPersistent_Interstitial15_9 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent9_InterstitialFibrosisNT15_2000points.csv");
simulatedPersistent_Interstitial15_9_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Interstitial15_9', .1, 5, 200, 1000, 10,false)';
actTimesSim_Interstitial15_9 = [];

simulatedPersistent_Interstitial15_8_downsampled_norm = [];
simulatedPersistent_Interstitial15_9_downsampled_norm = [];


for i = 1:length(patient6xCm)
    dVvec = diff(simulatedPersistent_Interstitial15_6_downsampled(:,i));
    minVal =  min(dVvec);
    indFind = find(dVvec == minVal);
    indFind = indFind(1);
    actTimesSim_Interstitial15_6(i) = indFind; 

    simulatedPersistent_Interstitial15_6_downsampled_norm(:,i) = ( simulatedPersistent_Interstitial15_6_downsampled(:,i) - min(simulatedPersistent_Interstitial15_6_downsampled(:,i)) )./( max(simulatedPersistent_Interstitial15_6_downsampled(:,i)) - min(simulatedPersistent_Interstitial15_6_downsampled(:,i)) );
end
for i = 1:length(patient7xCm)
    dVvec2 = diff(simulatedPersistent_Interstitial15_7_downsampled(:,i));
    minVal2 =  min(dVvec2);
    indFind2 = find(dVvec2 == minVal2);
    if length(indFind2) > 0
        indFind2 = indFind2(1);
    else
        indFind2 = 0;
    end
    actTimesSim_Interstitial15_7(i) = indFind2; 

    simulatedPersistent_Interstitial15_7_downsampled_norm(:,i) = ( simulatedPersistent_Interstitial15_7_downsampled(:,i) - min(simulatedPersistent_Interstitial15_7_downsampled(:,i)) )./( max(simulatedPersistent_Interstitial15_7_downsampled(:,i)) - min(simulatedPersistent_Interstitial15_7_downsampled(:,i)) );
end
for i = 1:length(patient8xCm)
    dVvec2 = diff(simulatedPersistent_Interstitial15_8_downsampled(:,i));
    minVal2 =  min(dVvec2);
    indFind2 = find(dVvec2 == minVal2);
    if length(indFind2) > 0
        indFind2 = indFind2(1);
    else
        indFind2 = 0;
    end
    actTimesSim_Interstitial15_8(i) = indFind2; 

    simulatedPersistent_Interstitial15_8_downsampled_norm(:,i) = ( simulatedPersistent_Interstitial15_8_downsampled(:,i) - min(simulatedPersistent_Interstitial15_8_downsampled(:,i)) )./( max(simulatedPersistent_Interstitial15_8_downsampled(:,i)) - min(simulatedPersistent_Interstitial15_8_downsampled(:,i)) );
end
for i = 1:length(patient9xCm)
    dVvec2 = diff(simulatedPersistent_Interstitial15_9_downsampled(:,i));
    minVal2 =  min(dVvec2);
    indFind2 = find(dVvec2 == minVal2);
    if length(indFind2) > 0
        indFind2 = indFind2(1);
    else
        indFind2 = 0;
    end
    actTimesSim_Interstitial15_9(i) = indFind2; 

    simulatedPersistent_Interstitial15_9_downsampled_norm(:,i) = ( simulatedPersistent_Interstitial15_9_downsampled(:,i) - min(simulatedPersistent_Interstitial15_9_downsampled(:,i)) )./( max(simulatedPersistent_Interstitial15_9_downsampled(:,i)) - min(simulatedPersistent_Interstitial15_9_downsampled(:,i)) );
end




%Simulated signals
simulatedPersistent_Interstitial25_6 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent6_InterstitialFibrosisNT25_2000points.csv");
simulatedPersistent_Interstitial25_6_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Interstitial25_6', .1, 5, 200, 1000, 10,false)';
actTimesSim_Interstitial25_6 = [];

simulatedPersistent_Interstitial25_7 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent7_InterstitialFibrosisNT25_2000points.csv");
simulatedPersistent_Interstitial25_7_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Interstitial25_7', .1, 5, 200, 1000, 10,false)';
actTimesSim_Interstitial25_7 = [];

simulatedPersistent_Interstitial25_6_downsampled_norm = [];
simulatedPersistent_Interstitial25_7_downsampled_norm = [];

simulatedPersistent_Interstitial25_8 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent8_InterstitialFibrosisNT25_2000points.csv");
simulatedPersistent_Interstitial25_8_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Interstitial25_8', .1, 5, 200, 1000, 10,false)';
actTimesSim_Interstitial25_8 = [];

simulatedPersistent_Interstitial25_9 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent9_InterstitialFibrosisNT25_2000points.csv");
simulatedPersistent_Interstitial25_9_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Interstitial25_9', .1, 5, 200, 1000, 10,false)';
actTimesSim_Interstitial25_9 = [];

simulatedPersistent_Interstitial25_8_downsampled_norm = [];
simulatedPersistent_Interstitial25_9_downsampled_norm = [];


for i = 1:length(patient6xCm)
    dVvec = diff(simulatedPersistent_Interstitial25_6_downsampled(:,i));
    minVal =  min(dVvec);
    indFind = find(dVvec == minVal);
    indFind = indFind(1);
    actTimesSim_Interstitial25_6(i) = indFind; 

    simulatedPersistent_Interstitial25_6_downsampled_norm(:,i) = ( simulatedPersistent_Interstitial25_6_downsampled(:,i) - min(simulatedPersistent_Interstitial25_6_downsampled(:,i)) )./( max(simulatedPersistent_Interstitial25_6_downsampled(:,i)) - min(simulatedPersistent_Interstitial25_6_downsampled(:,i)) );
end
for i = 1:length(patient7xCm)
    dVvec2 = diff(simulatedPersistent_Interstitial25_7_downsampled(:,i));
    minVal2 =  min(dVvec2);
    indFind2 = find(dVvec2 == minVal2);
    if length(indFind2) > 0
        indFind2 = indFind2(1);
    else
        indFind2 = 0;
    end
    actTimesSim_Interstitial25_7(i) = indFind2; 

    simulatedPersistent_Interstitial25_7_downsampled_norm(:,i) = ( simulatedPersistent_Interstitial25_7_downsampled(:,i) - min(simulatedPersistent_Interstitial25_7_downsampled(:,i)) )./( max(simulatedPersistent_Interstitial25_7_downsampled(:,i)) - min(simulatedPersistent_Interstitial25_7_downsampled(:,i)) );
end
for i = 1:length(patient8xCm)
    dVvec2 = diff(simulatedPersistent_Interstitial25_8_downsampled(:,i));
    minVal2 =  min(dVvec2);
    indFind2 = find(dVvec2 == minVal2);
    if length(indFind2) > 0
        indFind2 = indFind2(1);
    else
        indFind2 = 0;
    end
    actTimesSim_Interstitial25_8(i) = indFind2; 

    simulatedPersistent_Interstitial25_8_downsampled_norm(:,i) = ( simulatedPersistent_Interstitial25_8_downsampled(:,i) - min(simulatedPersistent_Interstitial25_8_downsampled(:,i)) )./( max(simulatedPersistent_Interstitial25_8_downsampled(:,i)) - min(simulatedPersistent_Interstitial25_8_downsampled(:,i)) );
end
for i = 1:length(patient9xCm)
    dVvec2 = diff(simulatedPersistent_Interstitial25_9_downsampled(:,i));
    minVal2 =  min(dVvec2);
    indFind2 = find(dVvec2 == minVal2);
    if length(indFind2) > 0
        indFind2 = indFind2(1);
    else
        indFind2 = 0;
    end
    actTimesSim_Interstitial25_9(i) = indFind2; 

    simulatedPersistent_Interstitial25_9_downsampled_norm(:,i) = ( simulatedPersistent_Interstitial25_9_downsampled(:,i) - min(simulatedPersistent_Interstitial25_9_downsampled(:,i)) )./( max(simulatedPersistent_Interstitial25_9_downsampled(:,i)) - min(simulatedPersistent_Interstitial25_9_downsampled(:,i)) );
end





figure
subplot(2,2,1)
plot3(patient6xSurfaceCm, patient6ySurfaceCm, patient6zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient6xCm,patient6yCm,patient6zCm,80,actTimes6,'filled')
colorbar
caxis([40, 140]);
title('Patient 6 Data')
subplot(2,2,2)
plot3(patient6xSurfaceCm, patient6ySurfaceCm, patient6zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient6xCm,patient6yCm,patient6zCm,80,actTimesSim_Interstitial15_6,'filled')
colorbar
caxis([40, 140]);
title('Simulation 6 Data')
subplot(2,2,3)
plot3(patient7xSurfaceCm, patient7ySurfaceCm, patient7zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient7xCm,patient7yCm,patient7zCm,80,actTimes7,'filled')
colorbar
caxis([50, 160]);
title('Patient 7 Data')
subplot(2,2,4)
plot3(patient7xSurfaceCm, patient7ySurfaceCm, patient7zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient7xCm,patient7yCm,patient7zCm,80,actTimesSim_Interstitial15_7,'filled')
colorbar
caxis([50, 160]);
title('Simulation 7 Data')






figure
subplot(2,2,1)
plot3(patient8xSurfaceCm, patient8ySurfaceCm, patient8zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient8xCm,patient8yCm,patient8zCm,80,actTimes8,'filled')
colorbar
caxis([30, 110]);
title('Patient 8 Data')
subplot(2,2,2)
plot3(patient8xSurfaceCm, patient8ySurfaceCm, patient8zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient8xCm,patient8yCm,patient8zCm,80,actTimesSim_Interstitial15_8,'filled')
colorbar
caxis([30, 110]);
title('Simulation 8 Data')
subplot(2,2,3)
plot3(patient9xSurfaceCm, patient9ySurfaceCm, patient9zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient9xCm,patient9yCm,patient9zCm,80,actTimes9,'filled')
colorbar
caxis([30, 130]);
title('Patient 9 Data')
subplot(2,2,4)
plot3(patient9xSurfaceCm, patient9ySurfaceCm, patient9zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient9xCm,patient9yCm,patient9zCm,80,actTimesSim_Interstitial15_9,'filled')
colorbar
caxis([30, 130]);
title('Simulation 9 Data')





%INTERSTITIAL FIBROSIS
figure;
subplot(2,2,1)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{6}.map.continuousMapping.fullEgms(ind_P6_diffPoint1,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Interstitial15_6_downsampled(:,ind_P6_diffPoint1), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 1")
subplot(2,2,2)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{6}.map.continuousMapping.fullEgms(ind_P6_diffPoint2,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Interstitial15_6_downsampled(:,ind_P6_diffPoint2), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 2")
subplot(2,2,3)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{6}.map.continuousMapping.fullEgms(ind_P6_diffPoint3,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Interstitial15_6_downsampled(:,ind_P6_diffPoint3), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 3")
subplot(2,2,4)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{6}.map.continuousMapping.fullEgms(ind_P6_diffPoint4,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Interstitial15_6_downsampled(:,ind_P6_diffPoint4), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 4")
sgtitle('H2B signals')





figure;
subplot(2,2,1)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{7}.map.continuousMapping.fullEgms(ind_P7_diffPoint1,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Interstitial15_7_downsampled(:,ind_P7_diffPoint1), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 1")
subplot(2,2,2)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{7}.map.continuousMapping.fullEgms(ind_P7_diffPoint2,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Interstitial15_7_downsampled(:,ind_P7_diffPoint2), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 2")
subplot(2,2,3)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{7}.map.continuousMapping.fullEgms(ind_P7_diffPoint3,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Interstitial15_7_downsampled(:,ind_P7_diffPoint3), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 3")
subplot(2,2,4)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{7}.map.continuousMapping.fullEgms(ind_P7_diffPoint4,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Interstitial15_7_downsampled(:,ind_P7_diffPoint4), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 4")
sgtitle('H2C signals')





figure;
subplot(2,2,1)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{8}.map.continuousMapping.fullEgms(ind_P8_diffPoint1,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Interstitial15_8_downsampled(:,ind_P8_diffPoint1), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 1")
subplot(2,2,2)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{8}.map.continuousMapping.fullEgms(ind_P8_diffPoint2,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Interstitial15_8_downsampled(:,ind_P8_diffPoint2), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 2")
subplot(2,2,3)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{8}.map.continuousMapping.fullEgms(ind_P8_diffPoint3,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Interstitial15_8_downsampled(:,ind_P8_diffPoint3), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 3")
subplot(2,2,4)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{8}.map.continuousMapping.fullEgms(ind_P8_diffPoint4,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Interstitial15_8_downsampled(:,ind_P8_diffPoint4), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 4")
sgtitle('H2D signals')






figure;
subplot(2,2,1)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{9}.map.continuousMapping.fullEgms(ind_P9_diffPoint1,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Interstitial15_9_downsampled(:,ind_P9_diffPoint1), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 1")
subplot(2,2,2)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{9}.map.continuousMapping.fullEgms(ind_P9_diffPoint2,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Interstitial15_9_downsampled(:,ind_P9_diffPoint2), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 2")
subplot(2,2,3)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{9}.map.continuousMapping.fullEgms(ind_P9_diffPoint3,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Interstitial15_9_downsampled(:,ind_P9_diffPoint3), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 3")
subplot(2,2,4)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{9}.map.continuousMapping.fullEgms(ind_P9_diffPoint4,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Interstitial15_9_downsampled(:,ind_P9_diffPoint4), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Target 4")
sgtitle('H2E signals')



%% PERSISTENT PATIENTS TO CHECK: PATCHY FIBROSIS



%CHECKING THINGS


%Simulated signals
simulatedPersistent_Patchy15_6 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent6_PatchyFibrosisNT15_2000points.csv");
simulatedPersistent_Patchy15_6_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Patchy15_6', .1, 5, 200, 1000, 10,false)';
actTimesSim_Patchy15_6 = [];

simulatedPersistent_Patchy15_7 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent7_PatchyFibrosisNT15_2000points.csv");
simulatedPersistent_Patchy15_7_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Patchy15_7', .1, 5, 200, 1000, 10,false)';
actTimesSim_Patchy15_7 = [];

simulatedPersistent_Patchy15_6_downsampled_norm = [];
simulatedPersistent_Patchy15_7_downsampled_norm = [];

simulatedPersistent_Patchy15_8 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent8_PatchyFibrosisNT15_2000points.csv");
simulatedPersistent_Patchy15_8_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Patchy15_8', .1, 5, 200, 1000, 10,false)';
actTimesSim_Patchy15_8 = [];

simulatedPersistent_Patchy15_9 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent9_PatchyFibrosisNT15_2000points.csv");
simulatedPersistent_Patchy15_9_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Patchy15_9', .1, 5, 200, 1000, 10,false)';
actTimesSim_Patchy15_9 = [];

simulatedPersistent_Patchy15_8_downsampled_norm = [];
simulatedPersistent_Patchy15_9_downsampled_norm = [];


for i = 1:length(patient6xCm)
    dVvec = diff(simulatedPersistent_Patchy15_6_downsampled(:,i));
    minVal =  min(dVvec);
    indFind = find(dVvec == minVal);
    indFind = indFind(1);
    actTimesSim_Patchy15_6(i) = indFind; 

    simulatedPersistent_Patchy15_6_downsampled_norm(:,i) = ( simulatedPersistent_Patchy15_6_downsampled(:,i) - min(simulatedPersistent_Patchy15_6_downsampled(:,i)) )./( max(simulatedPersistent_Patchy15_6_downsampled(:,i)) - min(simulatedPersistent_Patchy15_6_downsampled(:,i)) );
end
for i = 1:length(patient7xCm)
    dVvec2 = diff(simulatedPersistent_Patchy15_7_downsampled(:,i));
    minVal2 =  min(dVvec2);
    indFind2 = find(dVvec2 == minVal2);
    if length(indFind2) > 0
        indFind2 = indFind2(1);
    else
        indFind2 = 0;
    end
    actTimesSim_Patchy15_7(i) = indFind2; 

    simulatedPersistent_Patchy15_7_downsampled_norm(:,i) = ( simulatedPersistent_Patchy15_7_downsampled(:,i) - min(simulatedPersistent_Patchy15_7_downsampled(:,i)) )./( max(simulatedPersistent_Patchy15_7_downsampled(:,i)) - min(simulatedPersistent_Patchy15_7_downsampled(:,i)) );
end
for i = 1:length(patient8xCm)
    dVvec2 = diff(simulatedPersistent_Patchy15_8_downsampled(:,i));
    minVal2 =  min(dVvec2);
    indFind2 = find(dVvec2 == minVal2);
    if length(indFind2) > 0
        indFind2 = indFind2(1);
    else
        indFind2 = 0;
    end
    actTimesSim_Patchy15_8(i) = indFind2; 

    simulatedPersistent_Patchy15_8_downsampled_norm(:,i) = ( simulatedPersistent_Patchy15_8_downsampled(:,i) - min(simulatedPersistent_Patchy15_8_downsampled(:,i)) )./( max(simulatedPersistent_Patchy15_8_downsampled(:,i)) - min(simulatedPersistent_Patchy15_8_downsampled(:,i)) );
end
for i = 1:length(patient9xCm)
    dVvec2 = diff(simulatedPersistent_Patchy15_9_downsampled(:,i));
    minVal2 =  min(dVvec2);
    indFind2 = find(dVvec2 == minVal2);
    if length(indFind2) > 0
        indFind2 = indFind2(1);
    else
        indFind2 = 0;
    end
    actTimesSim_Patchy15_9(i) = indFind2; 

    simulatedPersistent_Patchy15_9_downsampled_norm(:,i) = ( simulatedPersistent_Patchy15_9_downsampled(:,i) - min(simulatedPersistent_Patchy15_9_downsampled(:,i)) )./( max(simulatedPersistent_Patchy15_9_downsampled(:,i)) - min(simulatedPersistent_Patchy15_9_downsampled(:,i)) );
end




%Simulated signals
simulatedPersistent_Patchy25_6 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent6_PatchyFibrosisNT25_2000points.csv");
simulatedPersistent_Patchy25_6_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Patchy25_6', .1, 5, 200, 1000, 10,false)';
actTimesSim_Patchy25_6 = [];

simulatedPersistent_Patchy25_7 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent7_PatchyFibrosisNT25_2000points.csv");
simulatedPersistent_Patchy25_7_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Patchy25_7', .1, 5, 200, 1000, 10,false)';
actTimesSim_Patchy25_7 = [];

simulatedPersistent_Patchy25_6_downsampled_norm = [];
simulatedPersistent_Patchy25_7_downsampled_norm = [];

simulatedPersistent_Patchy25_8 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent8_PatchyFibrosisNT25_2000points.csv");
simulatedPersistent_Patchy25_8_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Patchy25_8', .1, 5, 200, 1000, 10,false)';
actTimesSim_Patchy25_8 = [];

simulatedPersistent_Patchy25_9 = csvread("../SimResults/P1simsStuff/Ve_Catheter_MonoAp_SBDF1_Courtemanche_monodomain_NoFibrosis_Persistent9_PatchyFibrosisNT25_2000points.csv");
simulatedPersistent_Patchy25_9_downsampled = UpsampleAndFilterFunction(simulatedPersistent_Patchy25_9', .1, 5, 200, 1000, 10,false)';
actTimesSim_Patchy25_9 = [];

simulatedPersistent_Patchy25_8_downsampled_norm = [];
simulatedPersistent_Patchy25_9_downsampled_norm = [];


for i = 1:length(patient6xCm)
    dVvec = diff(simulatedPersistent_Patchy25_6_downsampled(:,i));
    minVal =  min(dVvec);
    indFind = find(dVvec == minVal);
    indFind = indFind(1);
    actTimesSim_Patchy25_6(i) = indFind; 

    simulatedPersistent_Patchy25_6_downsampled_norm(:,i) = ( simulatedPersistent_Patchy25_6_downsampled(:,i) - min(simulatedPersistent_Patchy25_6_downsampled(:,i)) )./( max(simulatedPersistent_Patchy25_6_downsampled(:,i)) - min(simulatedPersistent_Patchy25_6_downsampled(:,i)) );
end
for i = 1:length(patient7xCm)
    dVvec2 = diff(simulatedPersistent_Patchy25_7_downsampled(:,i));
    minVal2 =  min(dVvec2);
    indFind2 = find(dVvec2 == minVal2);
    if length(indFind2) > 0
        indFind2 = indFind2(1);
    else
        indFind2 = 0;
    end
    actTimesSim_Patchy25_7(i) = indFind2; 

    simulatedPersistent_Patchy25_7_downsampled_norm(:,i) = ( simulatedPersistent_Patchy25_7_downsampled(:,i) - min(simulatedPersistent_Patchy25_7_downsampled(:,i)) )./( max(simulatedPersistent_Patchy25_7_downsampled(:,i)) - min(simulatedPersistent_Patchy25_7_downsampled(:,i)) );
end
for i = 1:length(patient8xCm)
    dVvec2 = diff(simulatedPersistent_Patchy25_8_downsampled(:,i));
    minVal2 =  min(dVvec2);
    indFind2 = find(dVvec2 == minVal2);
    if length(indFind2) > 0
        indFind2 = indFind2(1);
    else
        indFind2 = 0;
    end
    actTimesSim_Patchy25_8(i) = indFind2; 

    simulatedPersistent_Patchy25_8_downsampled_norm(:,i) = ( simulatedPersistent_Patchy25_8_downsampled(:,i) - min(simulatedPersistent_Patchy25_8_downsampled(:,i)) )./( max(simulatedPersistent_Patchy25_8_downsampled(:,i)) - min(simulatedPersistent_Patchy25_8_downsampled(:,i)) );
end
for i = 1:length(patient9xCm)
    dVvec2 = diff(simulatedPersistent_Patchy25_9_downsampled(:,i));
    minVal2 =  min(dVvec2);
    indFind2 = find(dVvec2 == minVal2);
    if length(indFind2) > 0
        indFind2 = indFind2(1);
    else
        indFind2 = 0;
    end
    actTimesSim_Patchy25_9(i) = indFind2; 

    simulatedPersistent_Patchy25_9_downsampled_norm(:,i) = ( simulatedPersistent_Patchy25_9_downsampled(:,i) - min(simulatedPersistent_Patchy25_9_downsampled(:,i)) )./( max(simulatedPersistent_Patchy25_9_downsampled(:,i)) - min(simulatedPersistent_Patchy25_9_downsampled(:,i)) );
end




figure
subplot(2,2,1)
plot3(patient6xSurfaceCm, patient6ySurfaceCm, patient6zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient6xCm,patient6yCm,patient6zCm,80,actTimes6,'filled')
colorbar
caxis([40, 140]);
title('Patient 6 Data')
subplot(2,2,2)
plot3(patient6xSurfaceCm, patient6ySurfaceCm, patient6zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient6xCm,patient6yCm,patient6zCm,80,actTimesSim_Patchy15_6,'filled')
colorbar
caxis([40, 140]);
title('Simulation 6 Data')
subplot(2,2,3)
plot3(patient7xSurfaceCm, patient7ySurfaceCm, patient7zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient7xCm,patient7yCm,patient7zCm,80,actTimes7,'filled')
colorbar
caxis([50, 160]);
title('Patient 7 Data')
subplot(2,2,4)
plot3(patient7xSurfaceCm, patient7ySurfaceCm, patient7zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient7xCm,patient7yCm,patient7zCm,80,actTimesSim_Patchy15_7,'filled')
colorbar
caxis([50, 160]);
title('Simulation 7 Data')






figure
subplot(2,2,1)
plot3(patient8xSurfaceCm, patient8ySurfaceCm, patient8zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient8xCm,patient8yCm,patient8zCm,80,actTimes8,'filled')
colorbar
caxis([30, 110]);
title('Patient 8 Data')
subplot(2,2,2)
plot3(patient8xSurfaceCm, patient8ySurfaceCm, patient8zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient8xCm,patient8yCm,patient8zCm,80,actTimesSim_Patchy15_8,'filled')
colorbar
caxis([30, 110]);
title('Simulation 8 Data')
subplot(2,2,3)
plot3(patient9xSurfaceCm, patient9ySurfaceCm, patient9zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient9xCm,patient9yCm,patient9zCm,80,actTimes9,'filled')
colorbar
caxis([30, 130]);
title('Patient 9 Data')
subplot(2,2,4)
plot3(patient9xSurfaceCm, patient9ySurfaceCm, patient9zSurfaceCm, 'ro','DisplayName','tissue')
hold on
scatter3(patient9xCm,patient9yCm,patient9zCm,80,actTimesSim_Patchy15_9,'filled')
colorbar
caxis([30, 130]);
title('Simulation 9 Data')





%PATCHY FIBROSIS
figure;
subplot(2,2,1)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{6}.map.continuousMapping.fullEgms(ind_P6_diffPoint1,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Patchy15_6_downsampled(:,ind_P6_diffPoint1), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 1")
subplot(2,2,2)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{6}.map.continuousMapping.fullEgms(ind_P6_diffPoint2,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Patchy15_6_downsampled(:,ind_P6_diffPoint2), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 2")
subplot(2,2,3)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{6}.map.continuousMapping.fullEgms(ind_P6_diffPoint3,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Patchy15_6_downsampled(:,ind_P6_diffPoint3), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 3")
subplot(2,2,4)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{6}.map.continuousMapping.fullEgms(ind_P6_diffPoint4,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Patchy15_6_downsampled(:,ind_P6_diffPoint4), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 4")
sgtitle('Patient 6: signals')





figure;
subplot(2,2,1)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{7}.map.continuousMapping.fullEgms(ind_P7_diffPoint1,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Patchy15_7_downsampled(:,ind_P7_diffPoint1), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 1")
subplot(2,2,2)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{7}.map.continuousMapping.fullEgms(ind_P7_diffPoint2,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Patchy15_7_downsampled(:,ind_P7_diffPoint2), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 2")
subplot(2,2,3)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{7}.map.continuousMapping.fullEgms(ind_P7_diffPoint3,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Patchy15_7_downsampled(:,ind_P7_diffPoint3), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 3")
subplot(2,2,4)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{7}.map.continuousMapping.fullEgms(ind_P7_diffPoint4,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Patchy15_7_downsampled(:,ind_P7_diffPoint4), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 4")
sgtitle('Patient 7: signals')





figure;
subplot(2,2,1)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{8}.map.continuousMapping.fullEgms(ind_P8_diffPoint1,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Patchy15_8_downsampled(:,ind_P8_diffPoint1), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 1")
subplot(2,2,2)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{8}.map.continuousMapping.fullEgms(ind_P8_diffPoint2,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Patchy15_8_downsampled(:,ind_P8_diffPoint2), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 2")
subplot(2,2,3)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{8}.map.continuousMapping.fullEgms(ind_P8_diffPoint3,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Patchy15_8_downsampled(:,ind_P8_diffPoint3), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 3")
subplot(2,2,4)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{8}.map.continuousMapping.fullEgms(ind_P8_diffPoint4,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Patchy15_8_downsampled(:,ind_P8_diffPoint4), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 4")
sgtitle('Patient 8: signals')






figure;
subplot(2,2,1)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{9}.map.continuousMapping.fullEgms(ind_P9_diffPoint1,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Patchy15_9_downsampled(:,ind_P9_diffPoint1), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 1")
subplot(2,2,2)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{9}.map.continuousMapping.fullEgms(ind_P9_diffPoint2,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Patchy15_9_downsampled(:,ind_P9_diffPoint2), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 2")
subplot(2,2,3)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{9}.map.continuousMapping.fullEgms(ind_P9_diffPoint3,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Patchy15_9_downsampled(:,ind_P9_diffPoint3), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 3")
subplot(2,2,4)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{9}.map.continuousMapping.fullEgms(ind_P9_diffPoint4,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedPersistent_Patchy15_9_downsampled(:,ind_P9_diffPoint4), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 4")
sgtitle('Patient 9: signals')


%%

threshLATallowed = 10;%ms

sameLATindex_Compact15_6 = [];
sameLATCorr_Compact15_6 = [];

for i = 1:length(patient6xCm)
    if abs(actTimes6(i) - actTimesSim_Compact15_6(i)) < threshLATallowed
        sameLATindex_Compact15_6 = [ sameLATindex_Compact15_6, i ];

        actSignal = squeeze(expDB.dataSet{6}.map.continuousMapping.fullEgms(i,1,1:200));
        correlation_coefficient = corrcoef(actSignal, simulatedPersistent_Compact15_6_downsampled(:,i));
        similarity_score = correlation_coefficient(1, 2);
        sameLATCorr2 = [sameLATCorr2, similarity_score];
    end
end

percentageSameLATs2 = (length(sameLATCorr2)/length(patient2xCm))*100
meanCorr2 = mean(sameLATCorr2)









%%
sameLATindex2 = [];
sameLATCorr2 = [];
sameLATindex3 = [];
sameLATCorr3 = [];

threshLATallowed = 10;

for i = 1:length(patient2xCm)
    if abs(actTimes2(i) - actTimesSim2(i)) < threshLATallowed
        sameLATindex2 = [ sameLATindex2, i ];

        actSignal2 = squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(i,1,1:200));
        correlation_coefficient = corrcoef(actSignal2, simulatedParoxysmal2_downsampled(:,i));
        similarity_score = correlation_coefficient(1, 2);
        sameLATCorr2 = [sameLATCorr2, similarity_score];
    end
end
for i = 1:length(patient3xCm)
    if abs(actTimes3(i) - actTimesSim3(i)) < threshLATallowed
        sameLATindex3 = [ sameLATindex3, i ];

        actSignal3 = squeeze(expDB.dataSet{3}.map.continuousMapping.fullEgms(i,1,1:200));
        correlation_coefficient = corrcoef(actSignal3, simulatedParoxysmal3_downsampled(:,i));
        similarity_score = correlation_coefficient(1, 2);
        sameLATCorr3 = [sameLATCorr3, similarity_score];
    end
end

percentageSameLATs2 = (length(sameLATCorr2)/length(patient2xCm))*100
percentageSameLATs3 = (length(sameLATCorr3)/length(patient3xCm))*100

meanCorr2 = mean(sameLATCorr2)
meanCorr3 = mean(sameLATCorr3)




figure;
subplot(2,2,1)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(174,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedParoxysmal2_downsampled(:,174), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 1")
subplot(2,2,2)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(175,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedParoxysmal2_downsampled(:,175), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 2")
subplot(2,2,3)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(139,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedParoxysmal2_downsampled(:,139), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 3")
subplot(2,2,4)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(177,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedParoxysmal2_downsampled(:,177), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 4")
sgtitle('Patient 2: Other signals')






figure;
subplot(2,2,1)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{3}.map.continuousMapping.fullEgms(93,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedParoxysmal3_downsampled(:,93), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 1")
subplot(2,2,2)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{3}.map.continuousMapping.fullEgms(165,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedParoxysmal3_downsampled(:,165), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 2")
subplot(2,2,3)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{3}.map.continuousMapping.fullEgms(232,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedParoxysmal3_downsampled(:,232), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 3")
subplot(2,2,4)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{3}.map.continuousMapping.fullEgms(630,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedParoxysmal3_downsampled(:,630), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 4")
sgtitle('Patient 3: Other signals')





figure;
subplot(2,2,1)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{3}.map.continuousMapping.fullEgms(175,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedParoxysmal3_downsampled(:,175), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 5")
subplot(2,2,2)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{3}.map.continuousMapping.fullEgms(343,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedParoxysmal3_downsampled(:,343), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 6")
subplot(2,2,3)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{3}.map.continuousMapping.fullEgms(832,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedParoxysmal3_downsampled(:,832), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 7")
subplot(2,2,4)
% plot( squeeze(expDB.dataSet{2}.map.continuousMapping.fullEgms(sameLATindex2(find(sameLATCorr2 > .65)),1,1:200)) );
plot( squeeze(expDB.dataSet{3}.map.continuousMapping.fullEgms(469,1,1:200)), 'r', "LineWidth", 3 );
hold on;
plot( simulatedParoxysmal3_downsampled(:,469), 'b', "LineWidth", 3 )
legend('Patient signal','Simulated Signal')
xlabel('time (ms)')
ylabel('Potential (mV)')
title("Random signal 8")
sgtitle('Patient 3: Other signals')




%%

indicestrial = find( deflectionCountFractionated_AllPatients_Paroxysmal <= 3 & patientFractionated_AllPatients_Paroxysmal == 2 );


figure;
for i = 1:2:length(indicestrial)
    plot(signalsFractionated_AllPatients_Paroxysmal(indicestrial(i),:)')
    title(deflectionCountFractionated_AllPatients_Paroxysmal(indicestrial(i)))
    pause(2)
end




indicestrial2 = find( deflectionCountAll_Paroxysmal <= 3 & AllPatientVecTotal_Paroxysmal == 2 );



figure;
for i = 1:2:length(indicestrial2)
    plot(AllSignalsOG_Paroxysmal(:,indicestrial2(i)))
    title(deflectionCountAll_Paroxysmal(indicestrial2(i)))
    pause(2)
end




indicestrial3 = find( patientComplexFractionated_AllPatients_Paroxysmal == 2 );



figure;
for i = 1:2:length(indicestrial3)
    plot(signalsComplexFractionatedOG_AllPatients_Paroxysmal(indicestrial3(i),:)')
    % title(deflectionCountAll_Paroxysmal(indicestrial3(i)))
    pause(2)
end





%signalsAllOther_Paroxysmal, signalsAllOther_Persistent
%AllPatientVecOther_Paroxysmal, AllPatientVecOther_Persistent


matrixCount_perPatient = [  length(find(AllPatientVecTotal_Persistent == 1)), length(find(patientFractionated_AllPatients_Persistent == 1)), (length(find(patientComplexFractionated_AllPatients_Persistent == 1))+length(find(AllPatientVecOther_Persistent == 1))); ...
                 length(find(AllPatientVecTotal_Paroxysmal == 2)), length(find(patientFractionated_AllPatients_Paroxysmal == 2)), (length(find(patientComplexFractionated_AllPatients_Paroxysmal == 2))+length(find(AllPatientVecOther_Paroxysmal == 2))); ...
                 length(find(AllPatientVecTotal_Paroxysmal == 3)), length(find(patientFractionated_AllPatients_Paroxysmal == 3)), (length(find(patientComplexFractionated_AllPatients_Paroxysmal == 3))+length(find(AllPatientVecOther_Paroxysmal == 3))); ...
                 length(find(AllPatientVecTotal_Paroxysmal == 4)), length(find(patientFractionated_AllPatients_Paroxysmal == 4)), (length(find(patientComplexFractionated_AllPatients_Paroxysmal == 4))+length(find(AllPatientVecOther_Paroxysmal == 4))); ...
                 length(find(AllPatientVecTotal_Paroxysmal == 5)), length(find(patientFractionated_AllPatients_Paroxysmal == 5)), (length(find(patientComplexFractionated_AllPatients_Paroxysmal == 5))+length(find(AllPatientVecOther_Paroxysmal == 5))); ...
                 length(find(AllPatientVecTotal_Persistent == 6)), length(find(patientFractionated_AllPatients_Persistent == 6)), (length(find(patientComplexFractionated_AllPatients_Persistent == 6))+length(find(AllPatientVecOther_Persistent == 6))); ...
                 length(find(AllPatientVecTotal_Persistent == 7)), length(find(patientFractionated_AllPatients_Persistent == 7)), (length(find(patientComplexFractionated_AllPatients_Persistent == 7))+length(find(AllPatientVecOther_Persistent == 7))); ...
                 length(find(AllPatientVecTotal_Persistent == 8)), length(find(patientFractionated_AllPatients_Persistent == 8)), (length(find(patientComplexFractionated_AllPatients_Persistent == 8))+length(find(AllPatientVecOther_Persistent == 8))); ...
                 length(find(AllPatientVecTotal_Persistent == 9)), length(find(patientFractionated_AllPatients_Persistent == 9)), (length(find(patientComplexFractionated_AllPatients_Persistent == 9))+length(find(AllPatientVecOther_Persistent == 9))); ...
                 length(find(AllPatientVecTotal_Paroxysmal == 10)), length(find(patientFractionated_AllPatients_Paroxysmal == 10)), (length(find(patientComplexFractionated_AllPatients_Paroxysmal == 10))+length(find(AllPatientVecOther_Paroxysmal == 10)))]


matrixCount_All = [ length(AllPatientVecTotal_Paroxysmal), length(patientFractionated_AllPatients_Paroxysmal), (length(patientComplexFractionated_AllPatients_Paroxysmal)+length(AllPatientVecOther_Paroxysmal));...
                    length(AllPatientVecTotal_Persistent), length(patientFractionated_AllPatients_Persistent), (length(patientComplexFractionated_AllPatients_Persistent)+length(AllPatientVecOther_Persistent));
                    (length(AllPatientVecTotal_Paroxysmal) + length(AllPatientVecTotal_Persistent)), (length(patientFractionated_AllPatients_Paroxysmal)+length(patientFractionated_AllPatients_Persistent)), ((length(patientComplexFractionated_AllPatients_Paroxysmal)+length(AllPatientVecOther_Paroxysmal))+(length(patientComplexFractionated_AllPatients_Persistent)+length(AllPatientVecOther_Persistent))) ]


characteristics_perPatient_mean = {};
for i = [2,3,4,5,10]
    for j = 1:7
        switch j
            case 1
                characteristics_perPatient_mean{i,j} =  createMeanStdForCell(downstrokeAll_Paroxysmal, AllPatientVecTotal_Paroxysmal, i, 2 );
            case 2
                characteristics_perPatient_mean{i,j} =  createMeanStdForCell(upstrokeRiseAll_Paroxysmal, AllPatientVecTotal_Paroxysmal, i, 2 );
            case 3
                characteristics_perPatient_mean{i,j} =  createMeanStdForCell(upstrokeBackAll_Paroxysmal, AllPatientVecTotal_Paroxysmal, i, 2 );
            case 4
                characteristics_perPatient_mean{i,j} =  createMeanStdForCell(signalWidthAll_Paroxysmal, AllPatientVecTotal_Paroxysmal, i, 2 );
            case 5
                characteristics_perPatient_mean{i,j} =  createMeanStdForCell(baselineAll_Paroxysmal, AllPatientVecTotal_Paroxysmal, i, 2 );
            case 6
                characteristics_perPatient_mean{i,j} =  createMeanStdForCell(EGMDurationAll_Paroxysmal, AllPatientVecTotal_Paroxysmal, i, 2 );
            case 7
                characteristics_perPatient_mean{i,j} =  createMeanStdForCell(PPAAll_Paroxysmal, AllPatientVecTotal_Paroxysmal, i, 2 );
        end
    end
end
for i = [1,6,7,8,9]
    for j = 1:7
        switch j
            case 1
                characteristics_perPatient_mean{i,j} =  createMeanStdForCell(downstrokeAll_Persistent, AllPatientVecTotal_Persistent, i, 2 );
            case 2
                characteristics_perPatient_mean{i,j} =  createMeanStdForCell(upstrokeRiseAll_Persistent, AllPatientVecTotal_Persistent, i, 2 );
            case 3
                characteristics_perPatient_mean{i,j} =  createMeanStdForCell(upstrokeBackAll_Persistent, AllPatientVecTotal_Persistent, i, 2 );
            case 4
                characteristics_perPatient_mean{i,j} =  createMeanStdForCell(signalWidthAll_Persistent, AllPatientVecTotal_Persistent, i, 2 );
            case 5
                characteristics_perPatient_mean{i,j} =  createMeanStdForCell(baselineAll_Persistent, AllPatientVecTotal_Persistent, i, 2 );
            case 6
                characteristics_perPatient_mean{i,j} =  createMeanStdForCell(EGMDurationAll_Persistent, AllPatientVecTotal_Persistent, i, 2 );
            case 7
                characteristics_perPatient_mean{i,j} =  createMeanStdForCell(PPAAll_Persistent, AllPatientVecTotal_Persistent, i, 2 );
        end
    end
end

characteristics_perPatient_mean = cell2table(characteristics_perPatient_mean);








characteristicsFract_perPatient_mean = {};
for i = [2,3,4,5,10]
    for j = 1:8
        switch j
            case 1
                characteristicsFract_perPatient_mean{i,j} =  createMeanStdForCell(downstrokeFractionated_AllPatients_Paroxysmal, patientFractionated_AllPatients_Paroxysmal, i, 2 );
            case 2
                characteristicsFract_perPatient_mean{i,j} =  createMeanStdForCell(upstrokeRiseFractionated_AllPatients_Paroxysmal, patientFractionated_AllPatients_Paroxysmal, i, 2 );
            case 3
                characteristicsFract_perPatient_mean{i,j} =  createMeanStdForCell(upstrokeBackFractionated_AllPatients_Paroxysmal, patientFractionated_AllPatients_Paroxysmal, i, 2 );
            case 4
                characteristicsFract_perPatient_mean{i,j} =  createMeanStdForCell(signalWidthFractionated_AllPatients_Paroxysmal, patientFractionated_AllPatients_Paroxysmal, i, 2 );
            case 5
                characteristicsFract_perPatient_mean{i,j} =  createMeanStdForCell(baselineFractionated_AllPatients_Paroxysmal, patientFractionated_AllPatients_Paroxysmal, i, 2 );
            case 6
                characteristicsFract_perPatient_mean{i,j} =  createMeanStdForCell(EGMDurationFractionated_AllPatients_Paroxysmal, patientFractionated_AllPatients_Paroxysmal, i, 2 );
            case 7
                characteristicsFract_perPatient_mean{i,j} =  createMeanStdForCell(PPAFractionated_AllPatients_Paroxysmal, patientFractionated_AllPatients_Paroxysmal, i, 2 );
            case 8
                characteristicsFract_perPatient_mean{i,j} =  createMeanStdForCell(deflectionCountFractionated_AllPatients_Paroxysmal, patientFractionated_AllPatients_Paroxysmal, i, 2 );
        end
    end
end
for i = [1,6,7,8,9]
    for j = 1:8
        switch j
            case 1
                characteristicsFract_perPatient_mean{i,j} =  createMeanStdForCell(downstrokeFractionated_AllPatients_Persistent, patientFractionated_AllPatients_Persistent, i, 2 );
            case 2
                characteristicsFract_perPatient_mean{i,j} =  createMeanStdForCell(upstrokeRiseFractionated_AllPatients_Persistent, patientFractionated_AllPatients_Persistent, i, 2 );
            case 3
                characteristicsFract_perPatient_mean{i,j} =  createMeanStdForCell(upstrokeBackFractionated_AllPatients_Persistent, patientFractionated_AllPatients_Persistent, i, 2 );
            case 4
                characteristicsFract_perPatient_mean{i,j} =  createMeanStdForCell(signalWidthFractionated_AllPatients_Persistent, patientFractionated_AllPatients_Persistent, i, 2 );
            case 5
                characteristicsFract_perPatient_mean{i,j} =  createMeanStdForCell(baselineFractionated_AllPatients_Persistent, patientFractionated_AllPatients_Persistent, i, 2 );
            case 6
                characteristicsFract_perPatient_mean{i,j} =  createMeanStdForCell(EGMDurationFractionated_AllPatients_Persistent, patientFractionated_AllPatients_Persistent, i, 2 );
            case 7
                characteristicsFract_perPatient_mean{i,j} =  createMeanStdForCell(PPAFractionated_AllPatients_Persistent, patientFractionated_AllPatients_Persistent, i, 2 );
            case 8
                characteristicsFract_perPatient_mean{i,j} =  createMeanStdForCell(deflectionCountFractionated_AllPatients_Persistent, patientFractionated_AllPatients_Persistent, i, 2 );
        end
    end
end

characteristicsFract_perPatient_mean = cell2table(characteristicsFract_perPatient_mean);








characteristics_AllPatients_mean = {};
i = 1;
for j = 1:7
    switch j
        case 1
            characteristics_AllPatients_mean{i,j} =  createMeanStdForCell(downstrokeAll_Paroxysmal, AllPatientVecTotal_Paroxysmal, 1:10, 2 );
        case 2
            characteristics_AllPatients_mean{i,j} =  createMeanStdForCell(upstrokeRiseAll_Paroxysmal, AllPatientVecTotal_Paroxysmal, 1:10, 2 );
        case 3
            characteristics_AllPatients_mean{i,j} =  createMeanStdForCell(upstrokeBackAll_Paroxysmal, AllPatientVecTotal_Paroxysmal, 1:10, 2 );
        case 4
            characteristics_AllPatients_mean{i,j} =  createMeanStdForCell(signalWidthAll_Paroxysmal, AllPatientVecTotal_Paroxysmal, 1:10, 2 );
        case 5
            characteristics_AllPatients_mean{i,j} =  createMeanStdForCell(baselineAll_Paroxysmal, AllPatientVecTotal_Paroxysmal, 1:10, 2 );
        case 6
            characteristics_AllPatients_mean{i,j} =  createMeanStdForCell(EGMDurationAll_Paroxysmal, AllPatientVecTotal_Paroxysmal, 1:10, 2 );
        case 7
            characteristics_AllPatients_mean{i,j} =  createMeanStdForCell(PPAAll_Paroxysmal, AllPatientVecTotal_Paroxysmal, 1:10, 2 );
    end
end
i = 2;
for j = 1:7
    switch j
        case 1
            characteristics_AllPatients_mean{i,j} =  createMeanStdForCell(downstrokeAll_Persistent, AllPatientVecTotal_Persistent, 1:10, 2 );
        case 2
            characteristics_AllPatients_mean{i,j} =  createMeanStdForCell(upstrokeRiseAll_Persistent, AllPatientVecTotal_Persistent, 1:10, 2 );
        case 3
            characteristics_AllPatients_mean{i,j} =  createMeanStdForCell(upstrokeBackAll_Persistent, AllPatientVecTotal_Persistent, 1:10, 2 );
        case 4
            characteristics_AllPatients_mean{i,j} =  createMeanStdForCell(signalWidthAll_Persistent, AllPatientVecTotal_Persistent, 1:10, 2 );
        case 5
            characteristics_AllPatients_mean{i,j} =  createMeanStdForCell(baselineAll_Persistent, AllPatientVecTotal_Persistent, 1:10, 2 );
        case 6
            characteristics_AllPatients_mean{i,j} =  createMeanStdForCell(EGMDurationAll_Persistent, AllPatientVecTotal_Persistent, 1:10, 2 );
        case 7
            characteristics_AllPatients_mean{i,j} =  createMeanStdForCell(PPAAll_Persistent, AllPatientVecTotal_Persistent, 1:10, 2 );
    end
end

characteristics_AllPatients_mean = cell2table(characteristics_AllPatients_mean);





characteristicsFract_AllPatients_mean = {};
i = 1;
for j = 1:8
    switch j
        case 1
            characteristicsFract_AllPatients_mean{i,j} =  createMeanStdForCell(downstrokeFractionated_AllPatients_Paroxysmal, patientFractionated_AllPatients_Paroxysmal, 1:10, 2 );
        case 2
            characteristicsFract_AllPatients_mean{i,j} =  createMeanStdForCell(upstrokeRiseFractionated_AllPatients_Paroxysmal, patientFractionated_AllPatients_Paroxysmal, 1:10, 2 );
        case 3
            characteristicsFract_AllPatients_mean{i,j} =  createMeanStdForCell(upstrokeBackFractionated_AllPatients_Paroxysmal, patientFractionated_AllPatients_Paroxysmal, 1:10, 2 );
        case 4
            characteristicsFract_AllPatients_mean{i,j} =  createMeanStdForCell(signalWidthFractionated_AllPatients_Paroxysmal, patientFractionated_AllPatients_Paroxysmal, 1:10, 2 );
        case 5
            characteristicsFract_AllPatients_mean{i,j} =  createMeanStdForCell(baselineFractionated_AllPatients_Paroxysmal, patientFractionated_AllPatients_Paroxysmal, 1:10, 2 );
        case 6
            characteristicsFract_AllPatients_mean{i,j} =  createMeanStdForCell(EGMDurationFractionated_AllPatients_Paroxysmal, patientFractionated_AllPatients_Paroxysmal, 1:10, 2 );
        case 7
            characteristicsFract_AllPatients_mean{i,j} =  createMeanStdForCell(PPAFractionated_AllPatients_Paroxysmal, patientFractionated_AllPatients_Paroxysmal, 1:10, 2 );
        case 8
            characteristicsFract_AllPatients_mean{i,j} =  createMeanStdForCell(deflectionCountFractionated_AllPatients_Paroxysmal, patientFractionated_AllPatients_Paroxysmal, 1:10, 2 );
    end
end
i = 2;
for j = 1:8
    switch j
        case 1
            characteristicsFract_AllPatients_mean{i,j} =  createMeanStdForCell(downstrokeFractionated_AllPatients_Persistent, patientFractionated_AllPatients_Persistent, 1:10, 2 );
        case 2
            characteristicsFract_AllPatients_mean{i,j} =  createMeanStdForCell(upstrokeRiseFractionated_AllPatients_Persistent, patientFractionated_AllPatients_Persistent, 1:10, 2 );
        case 3
            characteristicsFract_AllPatients_mean{i,j} =  createMeanStdForCell(upstrokeBackFractionated_AllPatients_Persistent, patientFractionated_AllPatients_Persistent, 1:10, 2 );
        case 4
            characteristicsFract_AllPatients_mean{i,j} =  createMeanStdForCell(signalWidthFractionated_AllPatients_Persistent, patientFractionated_AllPatients_Persistent, 1:10, 2 );
        case 5
            characteristicsFract_AllPatients_mean{i,j} =  createMeanStdForCell(baselineFractionated_AllPatients_Persistent, patientFractionated_AllPatients_Persistent, 1:10, 2 );
        case 6
            characteristicsFract_AllPatients_mean{i,j} =  createMeanStdForCell(EGMDurationFractionated_AllPatients_Persistent, patientFractionated_AllPatients_Persistent, 1:10, 2 );
        case 7
            characteristicsFract_AllPatients_mean{i,j} =  createMeanStdForCell(PPAFractionated_AllPatients_Persistent, patientFractionated_AllPatients_Persistent, 1:10, 2 );
        case 8
            characteristicsFract_AllPatients_mean{i,j} =  createMeanStdForCell(deflectionCountFractionated_AllPatients_Persistent, patientFractionated_AllPatients_Persistent, 1:10, 2 );
    end
end

characteristicsFract_AllPatients_mean = cell2table(characteristicsFract_AllPatients_mean);






%% CHECKING FIBROSIS PATTERNS


Pers7_compactFibrosisPoints = csvread('./P1simsStuff/updatedFibrosisPatterns/P1_Persistent7_CompactFibrosis_P6_60perc.csv');
Pers7_diffuseFibrosisPoints = csvread('./P1simsStuff/updatedFibrosisPatterns/P1_Persistent7_DiffuseFibrosis_P6_60perc.csv');
Pers7_interstitialFibrosisPoints = csvread('./P1simsStuff/updatedFibrosisPatterns/P1_Persistent7_InterstitialFibrosis_P6_60perc.csv');
Pers7_patchyFibrosisPoints = csvread('./P1simsStuff/updatedFibrosisPatterns/P1_Persistent7_PatchyFibrosis_P6_60perc.csv');

%geometryPoints_Persistent7


figure
plot3(geometryPoints_Persisitent7(:,3), geometryPoints_Persisitent7(:,4), geometryPoints_Persisitent7(:,5), 'o', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'y')
hold on
plot3(Pers7_compactFibrosisPoints(:,1), Pers7_compactFibrosisPoints(:,2), Pers7_compactFibrosisPoints(:,3), 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r')


figure
plot3(geometryPoints_Persisitent7(:,3), geometryPoints_Persisitent7(:,4), geometryPoints_Persisitent7(:,5), 'o', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'y')
hold on
plot3(Pers7_diffuseFibrosisPoints(:,1), Pers7_diffuseFibrosisPoints(:,2), Pers7_diffuseFibrosisPoints(:,3), 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r')


figure
plot3(geometryPoints_Persisitent7(:,3), geometryPoints_Persisitent7(:,4), geometryPoints_Persisitent7(:,5), 'o', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'y')
hold on
plot3(Pers7_interstitialFibrosisPoints(:,1), Pers7_interstitialFibrosisPoints(:,2), Pers7_interstitialFibrosisPoints(:,3), 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r')


figure
plot3(geometryPoints_Persisitent7(:,3), geometryPoints_Persisitent7(:,4), geometryPoints_Persisitent7(:,5), 'o', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'y')
hold on
plot3(Pers7_patchyFibrosisPoints(:,1), Pers7_patchyFibrosisPoints(:,2), Pers7_patchyFibrosisPoints(:,3), 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r')







