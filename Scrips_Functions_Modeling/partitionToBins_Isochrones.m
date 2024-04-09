function [binEdges2, binIndices_patient2, remainingIndices] = partitionToBins_Isochrones(projLAT2_patient, numBins, alreadyPartitions )


    if length(alreadyPartitions) > 0
        binEdges2 = alreadyPartitions;
        disp('Prescribed edges...')
    else
        binEdges2 = linspace( min(projLAT2_patient), max(projLAT2_patient), numBins+1 );
        binEdges2(1) = binEdges2(1) - 1;
        binEdges2(end) = binEdges2(end) + 1;
        disp('Obtaining edges...')
    end

    binIndices_patient2 = {};
    allIndicesUsed = [];
    
    for i = 1:length(binEdges2)-1
        binIndices_patient2{i} = find( projLAT2_patient >= binEdges2(i) & projLAT2_patient < binEdges2(i+1)  );
        allIndicesUsed = [allIndicesUsed, find( projLAT2_patient >= binEdges2(i) & projLAT2_patient < binEdges2(i+1)  )];
    end

    allIndicesUsed = unique(allIndicesUsed);
    fullIndices = 1:length(projLAT2_patient);

    remainingIndices = setdiff(1:numel(fullIndices), allIndicesUsed);

    

end