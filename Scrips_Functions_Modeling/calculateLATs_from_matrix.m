function loc_actTimesSim = calculateLATs_from_matrix(nameOfFile, phie_dt, normalizeBool)
    locSignals_unfiltered = csvread(nameOfFile);
    locSignals_unfiltered = locSignals_unfiltered(1:1400,:);
    locSignals = UpsampleAndFilterFunction(locSignals_unfiltered', phie_dt, 5, 200, 1000, 10,normalizeBool)';
    loc_actTimesSim = [];
    
    
    for i = 1:length(locSignals(1,:))
        dVvec = diff(locSignals(:,i));
        minVal =  min(dVvec);
        indFind = find(dVvec == minVal);
        indFind = indFind(1);
        loc_actTimesSim(i) = indFind; 
    end
end


