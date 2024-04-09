%% TRYING TO TRIM LAPW FROM WHOLE GEOMETRY

%The steps for extracting the surface points describing the LAPW are as
%follows:
% 1 - The surface points for the whole LA are obtained and stored in
% variables xsur, ysur, zsur which are the x,y,z coordinate vectors for the
% whole LA geometry.
% 2 - A point cloud object is created from these coordinates and downsampled
% to better appreciate the whole geometry. With this downsampled point
% cloud, using knowledge of where the LAPW is found (the wall where the
% four Pulmonary Veins come and join), the coordinates for 8 different points that determine
% the boundaries of where the LAPW is found are obtained. These 8 different
% points correspond to four points describing the corners of the LAPW
% (boundaries) which create the edges of the LAPW and four points in
% between these edges.
% 3 - Once the coordinates are obtained by finding the closest points from
% the original point cloud to these points, their inidices are stored in
% ind1 through ind8. It is IMPORTANT to note that the order on which these
% indices representing the points are stored in variables ind1 through
% ind8. The variable ind1 corresponds to the first corner and then the rest
% of the points are registered in CLOCKWISE direction from the surface.
% This means that ind1, ind2, and ind3 correspond to the first edge; ind3,
% ind4, and ind5 correspond to the second edge next to edge 1 in the
% clockwise direction; ind5, ind6, and ind7 correspond to the third edge;
% and ind7, ind8, and ind1 correspond to the fourth and final edge
% representing the boundaries of the LAPW.
% 4 - Once the indices representing the boundaries are obtained, the
% function Extraction_LAPW is called. This functions takes the eight
% indices of importance, the original and downsampled pointCloud objects in
% that order, and a value representing the width of elements around the
% boundaries of the LAPW to include (a number to high includes too many
% points, but a number too low makes the middle of the LAPW to look
% sparse). The output of the function are the indices of the LAPW only from the
% original matrix containing the whole LA.


factorToDownSample = 20;

xsur = expDB.dataSet{7}.geometry.vertices(:,1);
ysur = expDB.dataSet{7}.geometry.vertices(:,2);
zsur = expDB.dataSet{7}.geometry.vertices(:,3);

actxsur = expDB.dataSet{7}.geometry.verticesTrimmed(:,1);
actysur = expDB.dataSet{7}.geometry.verticesTrimmed(:,2);
actzsur = expDB.dataSet{7}.geometry.verticesTrimmed(:,3);

ptCloudSur = pointCloud([xsur, ysur, zsur]);%create point cloud again once done
ptCloudSur2 = pcdownsample(ptCloudSur,'nonuniformGridSample',factorToDownSample);


newXsur = ptCloudSur2.Location(:,1);
newYsur = ptCloudSur2.Location(:,2);
newZsur = ptCloudSur2.Location(:,3);

figure; plot3(newXsur, newYsur, newZsur, 'ko', 'MarkerFaceColor','r', 'MarkerSize', 15)
xlabel('X'); ylabel('Y'); zlabel('Z')
set(gca,'FontSize',35)
% figure; pcshow(ptCloudSur, MarkerSize=20)





% [meshSurf,radiiSurf] = pc2surfacemesh(ptCloudSur,"ball-pivot");
% % figure
% surfaceMeshShow(meshSurf, BackgroundColor="White")

%%

threshToUse = .75;
ind1 = find( abs(newXsur - 60.402) < threshToUse & abs(newYsur - -36.352) < threshToUse & abs(newZsur - 132.851) < threshToUse  )
ind3 = find( abs(newXsur - 34.92) < threshToUse & abs(newYsur - -48.031) < threshToUse & abs(newZsur - 135.717) < threshToUse  )
ind5 = find( abs(newXsur - 46.11) < threshToUse & abs(newYsur - -71.024) < threshToUse & abs(newZsur - 103.082) < threshToUse  )
ind7 = find( abs(newXsur - 64.647) < threshToUse & abs(newYsur - -60.218) < threshToUse & abs(newZsur - 100.547) < threshToUse  )

ind2 = find( abs(newXsur - 46.107) < threshToUse & abs(newYsur - -47.251) < threshToUse & abs(newZsur - 135.023) < threshToUse  )
ind4 = find( abs(newXsur - 36.537) < threshToUse & abs(newYsur - -47.479) < threshToUse & abs(newZsur - 120.373) < threshToUse  )
ind6 = find( abs(newXsur - 57.365) < threshToUse & abs(newYsur - -64.969) < threshToUse & abs(newZsur - 100.305) < threshToUse  )
ind8 = find( abs(newXsur - 61.419) < threshToUse & abs(newYsur - -39.779) < threshToUse & abs(newZsur - 110.009) < threshToUse  )


ind1 = ind1(1);
ind2 = ind2(1);
ind3 = ind3(1);
ind4 = ind4(1);
ind5 = ind5(1);
ind6 = ind6(1);
ind7 = ind7(1);
ind8 = ind8(1);


figure;
% plot3(newXsur, newYsur, newZsur, 'ko', 'MarkerFaceColor','r', 'MarkerSize', 10)
hold on
% plot3(actxsur, actysur, actzsur, 'mo', 'MarkerFaceColor','m', 'MarkerSize', 10)
plot3(xsur, ysur, zsur, 'ko', 'MarkerFaceColor','r', 'MarkerSize', 8)
plot3(newXsur(ind1), newYsur(ind1), newZsur(ind1), 'ko', 'MarkerFaceColor','b', 'MarkerSize', 20)
plot3(newXsur(ind2), newYsur(ind2), newZsur(ind2), 'ko', 'MarkerFaceColor','b', 'MarkerSize', 20)
plot3(newXsur(ind3), newYsur(ind3), newZsur(ind3), 'ko', 'MarkerFaceColor','b', 'MarkerSize', 20)
plot3(newXsur(ind4), newYsur(ind4), newZsur(ind4), 'ko', 'MarkerFaceColor','b', 'MarkerSize', 20)
plot3(newXsur(ind5), newYsur(ind5), newZsur(ind5), 'ko', 'MarkerFaceColor','b', 'MarkerSize', 20)
plot3(newXsur(ind6), newYsur(ind6), newZsur(ind6), 'ko', 'MarkerFaceColor','b', 'MarkerSize', 20)
plot3(newXsur(ind7), newYsur(ind7), newZsur(ind7), 'ko', 'MarkerFaceColor','b', 'MarkerSize', 20)
plot3(newXsur(ind8), newYsur(ind8), newZsur(ind8), 'ko', 'MarkerFaceColor','b', 'MarkerSize', 20)
xlabel('X'); ylabel('Y'); zlabel('Z')
set(gca,'FontSize',35)
%%
%Number of elements around points of interest that can be added to final
%indices
numTestElem = 10;

allLAPWInd = Extraction_LAPW( ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8, ptCloudSur, ptCloudSur2, numTestElem  );





