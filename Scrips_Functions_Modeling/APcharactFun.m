
function [UpstrokeVec, APAVec, APDVec, RMPVec] = APcharactFun(VmVec, tVec)
%VmVec has to be a cell of vectors with info FOR THE VERY LAST AP ONLY.
%Same thing for tVec. 
%MAKE SURE THAT DIMENSIONS MATCH FOR VMVEC AND TVEC (individual samples
%should either  be both row vectors or both column vectors)

    UpstrokeVec = [];
    APAVec = [];
    APDVec = [];
    RMPVec = [];

    for i = 1:length(VmVec(1,:))
        
        UpstrokeVec(i) = max( diff( VmVec{i} ) ./ diff( tVec{i} ) );
        APAVec(i) = max(VmVec{i}) - min(VmVec{i});
        Normali1 = find( VmVec{i} <= (max(VmVec{i}) - .8*APAVec(i)) ); Normali2 = find( abs(diff(Normali1)) ==  max(abs(diff(Normali1))));
        APDVec(i) = tVec{i}( Normali1(Normali2+2) ) - tVec{i}( find( diff( VmVec{i} ) ./ diff( tVec{i} ) ==  UpstrokeVec(i) ) );
        RMPVec(i) = min( VmVec{i} );

    end
    
    
end



