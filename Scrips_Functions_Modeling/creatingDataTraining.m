
%70% of the signals for training, 20% for teting, and 10% for final
%testing: this is seen in patterns

matrixTraining = [];
vectorTraining = [];
matrixTest = [];
vectorTest = [];
matrixOther = [];
vectorOther = [];

% output vector will be [1, 0, 0] for 10%, [0, 1, 0] for 35%, and [0, 0, 1] for 60%


%TRAINING SET
for i = [6,7,8,9]
    for j = 1:7
        
        %10% DENSITY
        indexOfInterestC = find(Compact10_Bin_patNum == i & Compact10_Bin_pattNum == j);
        indexOfInterestD = find(Diffuse10_Bin_patNum == i & Diffuse10_Bin_pattNum == j);
        indexOfInterestI = find(Interstitial10_Bin_patNum == i & Interstitial10_Bin_pattNum == j);
        indexOfInterestP = find(Patchy10_Bin_patNum == i & Patchy10_Bin_pattNum == j);

        matrixTraining = [matrixTraining; mean((Compact10_Bin(indexOfInterestC,:)),1); mean((Diffuse10_Bin(indexOfInterestD,:)),1); mean((Interstitial10_Bin(indexOfInterestI,:)),1); mean((Patchy10_Bin(indexOfInterestP,:)),1); ];
        vectorTraining = [vectorTraining; [1, 0, 0]; [1, 0, 0]; [1, 0, 0]; [1, 0, 0]];

        %35% DENSITY
        indexOfInterestC = find(Compact35_Bin_patNum == i & Compact35_Bin_pattNum == j);
        indexOfInterestD = find(Diffuse35_Bin_patNum == i & Diffuse35_Bin_pattNum == j);
        indexOfInterestI = find(Interstitial35_Bin_patNum == i & Interstitial35_Bin_pattNum == j);
        indexOfInterestP = find(Patchy35_Bin_patNum == i & Patchy35_Bin_pattNum == j);

        matrixTraining = [matrixTraining; mean((Compact35_Bin(indexOfInterestC,:)),1); mean((Diffuse35_Bin(indexOfInterestD,:)),1); mean((Interstitial35_Bin(indexOfInterestI,:)),1); mean((Patchy35_Bin(indexOfInterestP,:)),1); ];
        vectorTraining = [vectorTraining; [0, 1, 0]; [0, 1, 0]; [0, 1, 0]; [0, 1, 0]];

        %60% DENSITY
        indexOfInterestC = find(Compact60_Bin_patNum == i & Compact60_Bin_pattNum == j);
        indexOfInterestD = find(Diffuse60_Bin_patNum == i & Diffuse60_Bin_pattNum == j);
        indexOfInterestI = find(Interstitial60_Bin_patNum == i & Interstitial60_Bin_pattNum == j);
        indexOfInterestP = find(Patchy60_Bin_patNum == i & Patchy60_Bin_pattNum == j);

        matrixTraining = [matrixTraining; mean((Compact60_Bin(indexOfInterestC,:)),1); mean((Diffuse60_Bin(indexOfInterestD,:)),1); mean((Interstitial60_Bin(indexOfInterestI,:)),1); mean((Patchy60_Bin(indexOfInterestP,:)),1); ];
        vectorTraining = [vectorTraining; [0, 0, 1]; [0, 0, 1]; [0, 0, 1]; [0, 0, 1]];
    
    end
end





%TEST SET
for i = [6,7,8,9]
    for j = 8:9
        
        %10% DENSITY
        indexOfInterestC = find(Compact10_Bin_patNum == i & Compact10_Bin_pattNum == j);
        indexOfInterestD = find(Diffuse10_Bin_patNum == i & Diffuse10_Bin_pattNum == j);
        indexOfInterestI = find(Interstitial10_Bin_patNum == i & Interstitial10_Bin_pattNum == j);
        indexOfInterestP = find(Patchy10_Bin_patNum == i & Patchy10_Bin_pattNum == j);

        matrixTest = [matrixTest; mean((Compact10_Bin(indexOfInterestC,:)),1); mean((Diffuse10_Bin(indexOfInterestD,:)),1); mean((Interstitial10_Bin(indexOfInterestI,:)),1); mean((Patchy10_Bin(indexOfInterestP,:)),1); ];
        vectorTest = [vectorTest; [1, 0, 0]; [1, 0, 0]; [1, 0, 0]; [1, 0, 0]];

        %35% DENSITY
        indexOfInterestC = find(Compact35_Bin_patNum == i & Compact35_Bin_pattNum == j);
        indexOfInterestD = find(Diffuse35_Bin_patNum == i & Diffuse35_Bin_pattNum == j);
        indexOfInterestI = find(Interstitial35_Bin_patNum == i & Interstitial35_Bin_pattNum == j);
        indexOfInterestP = find(Patchy35_Bin_patNum == i & Patchy35_Bin_pattNum == j);

        matrixTest = [matrixTest; mean((Compact35_Bin(indexOfInterestC,:)),1); mean((Diffuse35_Bin(indexOfInterestD,:)),1); mean((Interstitial35_Bin(indexOfInterestI,:)),1); mean((Patchy35_Bin(indexOfInterestP,:)),1); ];
        vectorTest = [vectorTest; [0, 1, 0]; [0, 1, 0]; [0, 1, 0]; [0, 1, 0]];

        %60% DENSITY
        indexOfInterestC = find(Compact60_Bin_patNum == i & Compact60_Bin_pattNum == j);
        indexOfInterestD = find(Diffuse60_Bin_patNum == i & Diffuse60_Bin_pattNum == j);
        indexOfInterestI = find(Interstitial60_Bin_patNum == i & Interstitial60_Bin_pattNum == j);
        indexOfInterestP = find(Patchy60_Bin_patNum == i & Patchy60_Bin_pattNum == j);

        matrixTest = [matrixTest; mean((Compact60_Bin(indexOfInterestC,:)),1); mean((Diffuse60_Bin(indexOfInterestD,:)),1); mean((Interstitial60_Bin(indexOfInterestI,:)),1); mean((Patchy60_Bin(indexOfInterestP,:)),1); ];
        vectorTest = [vectorTest; [0, 0, 1]; [0, 0, 1]; [0, 0, 1]; [0, 0, 1]];
    
    end
end



%OTHER SET
for i = [6,7,8,9]
        j = 10;
        
        %10% DENSITY
        indexOfInterestC = find(Compact10_Bin_patNum == i & Compact10_Bin_pattNum == j);
        indexOfInterestD = find(Diffuse10_Bin_patNum == i & Diffuse10_Bin_pattNum == j);
        indexOfInterestI = find(Interstitial10_Bin_patNum == i & Interstitial10_Bin_pattNum == j);
        indexOfInterestP = find(Patchy10_Bin_patNum == i & Patchy10_Bin_pattNum == j);

        matrixOther = [matrixOther; mean((Compact10_Bin(indexOfInterestC,:)),1); mean((Diffuse10_Bin(indexOfInterestD,:)),1); mean((Interstitial10_Bin(indexOfInterestI,:)),1); mean((Patchy10_Bin(indexOfInterestP,:)),1); ];
        vectorOther = [vectorOther; [1, 0, 0]; [1, 0, 0]; [1, 0, 0]; [1, 0, 0]];

        %35% DENSITY
        indexOfInterestC = find(Compact35_Bin_patNum == i & Compact35_Bin_pattNum == j);
        indexOfInterestD = find(Diffuse35_Bin_patNum == i & Diffuse35_Bin_pattNum == j);
        indexOfInterestI = find(Interstitial35_Bin_patNum == i & Interstitial35_Bin_pattNum == j);
        indexOfInterestP = find(Patchy35_Bin_patNum == i & Patchy35_Bin_pattNum == j);

        matrixOther = [matrixOther; mean((Compact35_Bin(indexOfInterestC,:)),1); mean((Diffuse35_Bin(indexOfInterestD,:)),1); mean((Interstitial35_Bin(indexOfInterestI,:)),1); mean((Patchy35_Bin(indexOfInterestP,:)),1); ];
        vectorOther = [vectorOther; [0, 1, 0]; [0, 1, 0]; [0, 1, 0]; [0, 1, 0]];

        %60% DENSITY
        indexOfInterestC = find(Compact60_Bin_patNum == i & Compact60_Bin_pattNum == j);
        indexOfInterestD = find(Diffuse60_Bin_patNum == i & Diffuse60_Bin_pattNum == j);
        indexOfInterestI = find(Interstitial60_Bin_patNum == i & Interstitial60_Bin_pattNum == j);
        indexOfInterestP = find(Patchy60_Bin_patNum == i & Patchy60_Bin_pattNum == j);

        matrixOther = [matrixOther; mean((Compact60_Bin(indexOfInterestC,:)),1); mean((Diffuse60_Bin(indexOfInterestD,:)),1); mean((Interstitial60_Bin(indexOfInterestI,:)),1); mean((Patchy60_Bin(indexOfInterestP,:)),1); ];
        vectorOther = [vectorOther; [0, 0, 1]; [0, 0, 1]; [0, 0, 1]; [0, 0, 1]];
    
end
            


matrixTraining = [matrixTraining(:,1:6), matrixTraining(:,8:9)];
matrixTest = [matrixTest(:,1:6), matrixTest(:,8:9)];
matrixOther = [matrixOther(:,1:6), matrixOther(:,8:9)];


dlmwrite('./fibrosisDegreTraining.csv', round(matrixTraining,4), 'delimiter', ',', 'precision', '%.4f');
dlmwrite('./fibrosisDegreTest.csv', round(matrixTest,4), 'delimiter', ',', 'precision', '%.4f');
dlmwrite('./fibrosisDegreOther.csv', round(matrixOther,4), 'delimiter', ',', 'precision', '%.4f');

dlmwrite('./fibrosisDegreTrainingAns.csv', round(vectorTraining,4), 'delimiter', ',', 'precision', '%.4f');
dlmwrite('./fibrosisDegreTestAns.csv', round(vectorTest,4), 'delimiter', ',', 'precision', '%.4f');
dlmwrite('./fibrosisDegreOtherAns.csv', round(vectorOther,4), 'delimiter', ',', 'precision', '%.4f');





