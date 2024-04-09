function [ tauVec, GlobalAT ] = DELTAmethodTau( egmcourt, dtresolution )

    numTriConnections = length(egmcourt(:,1,1));
    tauVec = [];
    GlobalAT = [];

    for i = 1:numTriConnections

%         i
        U = squeeze(max(egmcourt(i,1,:)) - min(egmcourt(i,1,:)));
        if U ~= 1.
            egmcourt(i,1,:) = (egmcourt(i,1,:) - min(egmcourt(i,1,:)))./(max(egmcourt(i,1,:)) - min(egmcourt(i,1,:)));
            egmcourt(i,2,:) = (egmcourt(i,2,:) - min(egmcourt(i,2,:)))./(max(egmcourt(i,2,:)) - min(egmcourt(i,2,:)));
            egmcourt(i,3,:)=egmcourt(i,1,:) - egmcourt(i,2,:);
            U = squeeze(max(egmcourt(i,1,:)) - min(egmcourt(i,1,:)));
        end

%         U

        dVmdt1 = squeeze(diff(egmcourt(i,1,:))./(dtresolution));
        dVmdt12 = squeeze(diff(egmcourt(i,2,:))./(dtresolution));

%         size((dVmdt1))
%         size((dVmdt12))
% 
%         find( ((dVmdt1) - .97*min((dVmdt1))) < .00001 )

        indicesElec1 = find( ((dVmdt1) - .97*min((dVmdt1))) < .00001 );
        indicesElec2 = find( ((dVmdt12) - .97*min((dVmdt12))) < .00001 );

        if indicesElec1(1) <= indicesElec2(1)
            tau = (U/(min((dVmdt1))))*asin(min(egmcourt(i,3,:))/U);
            AT = find( egmcourt(i,3,:) == min(egmcourt(i,3,:)) )*dtresolution;
        else
            tau = (U/(min((dVmdt1))))*asin(max(egmcourt(i,3,:))/U);
            AT = find( egmcourt(i,3,:) == max(egmcourt(i,3,:)) )*dtresolution;
        end
        tauVec(i) = tau;
        GlobalAT(i) = AT(1);

    end

end