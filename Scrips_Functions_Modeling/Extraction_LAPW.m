function allLAPWInd = Extraction_LAPW( ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8, ptCloudSur, ptCloudSur2, numTestElem  )

%This function delineates the edges of the LAPW using the provided 8
%indices creating the boundaries for the LAPW. For this, basic anatomical
%knowledge is recquired to identify where the PVs are and delimit the
%points where the LAPW could be found. Using the edges, all the points
%within those edges are found by "sweeping" in a "straight" line from side
%to side within the point cloud data.



%INPUTS:
% ind1 through ind8. It is IMPORTANT to note that the order on which these
% indices representing the points are stored in variables ind1 through
% ind8. The variable ind1 corresponds to the first corner and then the rest
% of the points are registered in CLOCKWISE direction from the surface.
% This means that ind1, ind2, and ind3 correspond to the first edge; ind3,
% ind4, and ind5 correspond to the second edge next to edge 1 in the
% clockwise direction; ind5, ind6, and ind7 correspond to the third edge;
% and ind7, ind8, and ind1 correspond to the fourth and final edge
% representing the boundaries of the LAPW.
%ptCloudSur - original point cloud containing surface points of the whole LA geometry
%ptCloudSur2 - downsampled point cloud containing surface points of the whole LA geometry
%numTestElem - Number of elements around points of interest that can be added to final indices




%OUTPUT: allLAPWInd - indices of the extracted points of interest of the downsampled point cloud




%preparing vectors and matrices
xsur = ptCloudSur.Location(:,1);
ysur = ptCloudSur.Location(:,2);
zsur = ptCloudSur.Location(:,3);

newXsur = ptCloudSur2.Location(:,1);
newYsur = ptCloudSur2.Location(:,2);
newZsur = ptCloudSur2.Location(:,3);


%first edge: first half
firstCorner = ind1;
secondCorner = ind2;
firstEdgeInd = [];
[indFirst,dists] = findNearestNeighbors(ptCloudSur,[ newXsur(firstCorner), newYsur(firstCorner), newZsur(firstCorner) ],1);
firstEdgeInd = [ firstEdgeInd; indFirst(1) ];
[indLast,dists] = findNearestNeighbors(ptCloudSur,[ newXsur(secondCorner), newYsur(secondCorner), newZsur(secondCorner) ],1);



dummInd = 1;
lastDistance = 10000;
currentDistance = 100;
while currentDistance < lastDistance
    disp("Moving through..."+num2str(dummInd));
    [indices,dists] = findNearestNeighbors(ptCloudSur,[ xsur(firstEdgeInd(dummInd)), ysur(firstEdgeInd(dummInd)), zsur(firstEdgeInd(dummInd)) ], 10);
    locVecD = [];
    for i=2:length(indices)
        locVecD(i-1) = sqrt( ( xsur( indLast ) - xsur( indices(i) ) )^2 + ( ysur( indLast ) - ysur( indices(i) ) )^2 + ( zsur( indLast ) - zsur( indices(i) ) )^2 );
    end
    nextInd = find(locVecD == min(locVecD)) + 1;
    firstEdgeInd = [ firstEdgeInd; indices(nextInd(1)) ];
    dummInd = dummInd + 1;
    lastDistance = sqrt( ( xsur( indLast ) - xsur( firstEdgeInd(dummInd-1) ) )^2 + ( ysur( indLast ) - ysur( firstEdgeInd(dummInd-1) ) )^2 + ( zsur( indLast ) - zsur( firstEdgeInd(dummInd-1) ) )^2 );
    currentDistance = sqrt( ( xsur( indLast ) - xsur( firstEdgeInd(dummInd) ) )^2 + ( ysur( indLast ) - ysur( firstEdgeInd(dummInd) ) )^2 + ( zsur( indLast ) - zsur( firstEdgeInd(dummInd) ) )^2 );
end


%first edge: second half
firstCorner = ind2;
secondCorner = ind3;
secondEdgeInd = [];
[indFirst,dists] = findNearestNeighbors(ptCloudSur,[ newXsur(firstCorner), newYsur(firstCorner), newZsur(firstCorner) ],1);
secondEdgeInd = [ secondEdgeInd; indFirst(1) ];
[indLast,dists] = findNearestNeighbors(ptCloudSur,[ newXsur(secondCorner), newYsur(secondCorner), newZsur(secondCorner) ],1);



dummInd = 1;
lastDistance = 10000;
currentDistance = 100;
while currentDistance < lastDistance
    disp("Moving through..."+num2str(dummInd));
    [indices,dists] = findNearestNeighbors(ptCloudSur,[ xsur(secondEdgeInd(dummInd)), ysur(secondEdgeInd(dummInd)), zsur(secondEdgeInd(dummInd)) ], 10);
    locVecD = [];
    for i=2:length(indices)
        locVecD(i-1) = sqrt( ( xsur( indLast ) - xsur( indices(i) ) )^2 + ( ysur( indLast ) - ysur( indices(i) ) )^2 + ( zsur( indLast ) - zsur( indices(i) ) )^2 );
    end
    nextInd = find(locVecD == min(locVecD)) + 1;
    secondEdgeInd = [ secondEdgeInd; indices(nextInd(1)) ];
    dummInd = dummInd + 1;
    lastDistance = sqrt( ( xsur( indLast ) - xsur( secondEdgeInd(dummInd-1) ) )^2 + ( ysur( indLast ) - ysur( secondEdgeInd(dummInd-1) ) )^2 + ( zsur( indLast ) - zsur( secondEdgeInd(dummInd-1) ) )^2 );
    currentDistance = sqrt( ( xsur( indLast ) - xsur( secondEdgeInd(dummInd) ) )^2 + ( ysur( indLast ) - ysur( secondEdgeInd(dummInd) ) )^2 + ( zsur( indLast ) - zsur( secondEdgeInd(dummInd) ) )^2 );
end


%second edge: first half
firstCorner = ind3;
secondCorner = ind4;
thirdEdgeInd = [];
[indFirst,dists] = findNearestNeighbors(ptCloudSur,[ newXsur(firstCorner), newYsur(firstCorner), newZsur(firstCorner) ],1);
thirdEdgeInd = [ thirdEdgeInd; indFirst(1) ];
[indLast,dists] = findNearestNeighbors(ptCloudSur,[ newXsur(secondCorner), newYsur(secondCorner), newZsur(secondCorner) ],1);

dummInd = 1;
lastDistance = 10000;
currentDistance = 100;
while currentDistance < lastDistance
    disp("Moving through..."+num2str(dummInd));
    [indices,dists] = findNearestNeighbors(ptCloudSur,[ xsur(thirdEdgeInd(dummInd)), ysur(thirdEdgeInd(dummInd)), zsur(thirdEdgeInd(dummInd)) ], 10);
    locVecD = [];
    for i=2:length(indices)
        locVecD(i-1) = sqrt( ( xsur( indLast ) - xsur( indices(i) ) )^2 + ( ysur( indLast ) - ysur( indices(i) ) )^2 + ( zsur( indLast ) - zsur( indices(i) ) )^2 );
    end
    nextInd = find(locVecD == min(locVecD)) + 1;
    thirdEdgeInd = [ thirdEdgeInd; indices(nextInd(1)) ];
    dummInd = dummInd + 1;
    lastDistance = sqrt( ( xsur( indLast ) - xsur( thirdEdgeInd(dummInd-1) ) )^2 + ( ysur( indLast ) - ysur( thirdEdgeInd(dummInd-1) ) )^2 + ( zsur( indLast ) - zsur( thirdEdgeInd(dummInd-1) ) )^2 );
    currentDistance = sqrt( ( xsur( indLast ) - xsur( thirdEdgeInd(dummInd) ) )^2 + ( ysur( indLast ) - ysur( thirdEdgeInd(dummInd) ) )^2 + ( zsur( indLast ) - zsur( thirdEdgeInd(dummInd) ) )^2 );
end



%second edge: second half
firstCorner = ind4;
secondCorner = ind5;
fourthEdgeInd = [];
[indFirst,dists] = findNearestNeighbors(ptCloudSur,[ newXsur(firstCorner), newYsur(firstCorner), newZsur(firstCorner) ],1);
fourthEdgeInd = [ fourthEdgeInd; indFirst(1) ];
[indLast,dists] = findNearestNeighbors(ptCloudSur,[ newXsur(secondCorner), newYsur(secondCorner), newZsur(secondCorner) ],1);

dummInd = 1;
lastDistance = 10000;
currentDistance = 100;
while currentDistance < lastDistance
    disp("Moving through..."+num2str(dummInd));
    [indices,dists] = findNearestNeighbors(ptCloudSur,[ xsur(fourthEdgeInd(dummInd)), ysur(fourthEdgeInd(dummInd)), zsur(fourthEdgeInd(dummInd)) ], 10);
    locVecD = [];
    for i=2:length(indices)
        locVecD(i-1) = sqrt( ( xsur( indLast ) - xsur( indices(i) ) )^2 + ( ysur( indLast ) - ysur( indices(i) ) )^2 + ( zsur( indLast ) - zsur( indices(i) ) )^2 );
    end
    nextInd = find(locVecD == min(locVecD)) + 1;
    fourthEdgeInd = [ fourthEdgeInd; indices(nextInd(1)) ];
    dummInd = dummInd + 1;
    lastDistance = sqrt( ( xsur( indLast ) - xsur( fourthEdgeInd(dummInd-1) ) )^2 + ( ysur( indLast ) - ysur( fourthEdgeInd(dummInd-1) ) )^2 + ( zsur( indLast ) - zsur( fourthEdgeInd(dummInd-1) ) )^2 );
    currentDistance = sqrt( ( xsur( indLast ) - xsur( fourthEdgeInd(dummInd) ) )^2 + ( ysur( indLast ) - ysur( fourthEdgeInd(dummInd) ) )^2 + ( zsur( indLast ) - zsur( fourthEdgeInd(dummInd) ) )^2 );
end



%third edge: first half
firstCorner = ind5;
secondCorner = ind6;
fifthEdgeInd = [];
[indFirst,dists] = findNearestNeighbors(ptCloudSur,[ newXsur(firstCorner), newYsur(firstCorner), newZsur(firstCorner) ],1);
fifthEdgeInd = [ fifthEdgeInd; indFirst(1) ];
[indLast,dists] = findNearestNeighbors(ptCloudSur,[ newXsur(secondCorner), newYsur(secondCorner), newZsur(secondCorner) ],1);

dummInd = 1;
lastDistance = 10000;
currentDistance = 100;
while currentDistance < lastDistance
    disp("Moving through..."+num2str(dummInd));
    [indices,dists] = findNearestNeighbors(ptCloudSur,[ xsur(fifthEdgeInd(dummInd)), ysur(fifthEdgeInd(dummInd)), zsur(fifthEdgeInd(dummInd)) ], 10);
    locVecD = [];
    for i=2:length(indices)
        locVecD(i-1) = sqrt( ( xsur( indLast ) - xsur( indices(i) ) )^2 + ( ysur( indLast ) - ysur( indices(i) ) )^2 + ( zsur( indLast ) - zsur( indices(i) ) )^2 );
    end
    nextInd = find(locVecD == min(locVecD)) + 1;
    fifthEdgeInd = [ fifthEdgeInd; indices(nextInd(1)) ];
    dummInd = dummInd + 1;
    lastDistance = sqrt( ( xsur( indLast ) - xsur( fifthEdgeInd(dummInd-1) ) )^2 + ( ysur( indLast ) - ysur( fifthEdgeInd(dummInd-1) ) )^2 + ( zsur( indLast ) - zsur( fifthEdgeInd(dummInd-1) ) )^2 );
    currentDistance = sqrt( ( xsur( indLast ) - xsur( fifthEdgeInd(dummInd) ) )^2 + ( ysur( indLast ) - ysur( fifthEdgeInd(dummInd) ) )^2 + ( zsur( indLast ) - zsur( fifthEdgeInd(dummInd) ) )^2 );
end



%third edge: second half
firstCorner = ind6;
secondCorner = ind7;
sixthEdgeInd = [];
[indFirst,dists] = findNearestNeighbors(ptCloudSur,[ newXsur(firstCorner), newYsur(firstCorner), newZsur(firstCorner) ],1);
sixthEdgeInd = [ sixthEdgeInd; indFirst(1) ];
[indLast,dists] = findNearestNeighbors(ptCloudSur,[ newXsur(secondCorner), newYsur(secondCorner), newZsur(secondCorner) ],1);

dummInd = 1;
lastDistance = 10000;
currentDistance = 100;
while currentDistance < lastDistance
    disp("Moving through..."+num2str(dummInd));
    [indices,dists] = findNearestNeighbors(ptCloudSur,[ xsur(sixthEdgeInd(dummInd)), ysur(sixthEdgeInd(dummInd)), zsur(sixthEdgeInd(dummInd)) ], 10);
    locVecD = [];
    for i=2:length(indices)
        locVecD(i-1) = sqrt( ( xsur( indLast ) - xsur( indices(i) ) )^2 + ( ysur( indLast ) - ysur( indices(i) ) )^2 + ( zsur( indLast ) - zsur( indices(i) ) )^2 );
    end
    nextInd = find(locVecD == min(locVecD)) + 1;
    sixthEdgeInd = [ sixthEdgeInd; indices(nextInd(1)) ];
    dummInd = dummInd + 1;
    lastDistance = sqrt( ( xsur( indLast ) - xsur( sixthEdgeInd(dummInd-1) ) )^2 + ( ysur( indLast ) - ysur( sixthEdgeInd(dummInd-1) ) )^2 + ( zsur( indLast ) - zsur( sixthEdgeInd(dummInd-1) ) )^2 );
    currentDistance = sqrt( ( xsur( indLast ) - xsur( sixthEdgeInd(dummInd) ) )^2 + ( ysur( indLast ) - ysur( sixthEdgeInd(dummInd) ) )^2 + ( zsur( indLast ) - zsur( sixthEdgeInd(dummInd) ) )^2 );
end



%fourth edge: first half
firstCorner = ind7;
secondCorner = ind8;
seventhEdgeInd = [];
[indFirst,dists] = findNearestNeighbors(ptCloudSur,[ newXsur(firstCorner), newYsur(firstCorner), newZsur(firstCorner) ],1);
seventhEdgeInd = [ seventhEdgeInd; indFirst(1) ];
[indLast,dists] = findNearestNeighbors(ptCloudSur,[ newXsur(secondCorner), newYsur(secondCorner), newZsur(secondCorner) ],1);

dummInd = 1;
lastDistance = 10000;
currentDistance = 100;
while currentDistance < lastDistance
    disp("Moving through..."+num2str(dummInd));
    [indices,dists] = findNearestNeighbors(ptCloudSur,[ xsur(seventhEdgeInd(dummInd)), ysur(seventhEdgeInd(dummInd)), zsur(seventhEdgeInd(dummInd)) ], 10);
    locVecD = [];
    for i=2:length(indices)
        locVecD(i-1) = sqrt( ( xsur( indLast ) - xsur( indices(i) ) )^2 + ( ysur( indLast ) - ysur( indices(i) ) )^2 + ( zsur( indLast ) - zsur( indices(i) ) )^2 );
    end
    nextInd = find(locVecD == min(locVecD)) + 1;
    seventhEdgeInd = [ seventhEdgeInd; indices(nextInd(1)) ];
    dummInd = dummInd + 1;
    lastDistance = sqrt( ( xsur( indLast ) - xsur( seventhEdgeInd(dummInd-1) ) )^2 + ( ysur( indLast ) - ysur( seventhEdgeInd(dummInd-1) ) )^2 + ( zsur( indLast ) - zsur( seventhEdgeInd(dummInd-1) ) )^2 );
    currentDistance = sqrt( ( xsur( indLast ) - xsur( seventhEdgeInd(dummInd) ) )^2 + ( ysur( indLast ) - ysur( seventhEdgeInd(dummInd) ) )^2 + ( zsur( indLast ) - zsur( seventhEdgeInd(dummInd) ) )^2 );
end



%fourth edge: second half
firstCorner = ind8;
secondCorner = ind1;
eigthEdgeInd = [];
[indFirst,dists] = findNearestNeighbors(ptCloudSur,[ newXsur(firstCorner), newYsur(firstCorner), newZsur(firstCorner) ],1);
eigthEdgeInd = [ eigthEdgeInd; indFirst(1) ];
[indLast,dists] = findNearestNeighbors(ptCloudSur,[ newXsur(secondCorner), newYsur(secondCorner), newZsur(secondCorner) ],1);

dummInd = 1;
lastDistance = 10000;
currentDistance = 100;
while currentDistance < lastDistance
    disp("Moving through..."+num2str(dummInd));
    [indices,dists] = findNearestNeighbors(ptCloudSur,[ xsur(eigthEdgeInd(dummInd)), ysur(eigthEdgeInd(dummInd)), zsur(eigthEdgeInd(dummInd)) ], 10);
    locVecD = [];
    for i=2:length(indices)
        locVecD(i-1) = sqrt( ( xsur( indLast ) - xsur( indices(i) ) )^2 + ( ysur( indLast ) - ysur( indices(i) ) )^2 + ( zsur( indLast ) - zsur( indices(i) ) )^2 );
    end
    nextInd = find(locVecD == min(locVecD)) + 1;
    eigthEdgeInd = [ eigthEdgeInd; indices(nextInd(1)) ];
    dummInd = dummInd + 1;
    lastDistance = sqrt( ( xsur( indLast ) - xsur( eigthEdgeInd(dummInd-1) ) )^2 + ( ysur( indLast ) - ysur( eigthEdgeInd(dummInd-1) ) )^2 + ( zsur( indLast ) - zsur( eigthEdgeInd(dummInd-1) ) )^2 );
    currentDistance = sqrt( ( xsur( indLast ) - xsur( eigthEdgeInd(dummInd) ) )^2 + ( ysur( indLast ) - ysur( eigthEdgeInd(dummInd) ) )^2 + ( zsur( indLast ) - zsur( eigthEdgeInd(dummInd) ) )^2 );
end


% 
% figure;
% % plot3(newXsur, newYsur, newZsur, 'ko', 'MarkerFaceColor','r', 'MarkerSize', 10)
% plot3(actxsur, actysur, actzsur, 'mo', 'MarkerFaceColor','m', 'MarkerSize', 10)
% hold on
% plot3(xsur, ysur, zsur, 'ko', 'MarkerFaceColor','r', 'MarkerSize', 8)
% plot3(newXsur(ind1), newYsur(ind1), newZsur(ind1), 'ko', 'MarkerFaceColor','b', 'MarkerSize', 15)
% plot3(newXsur(ind2), newYsur(ind2), newZsur(ind2), 'ko', 'MarkerFaceColor','b', 'MarkerSize', 15)
% plot3(newXsur(ind3), newYsur(ind3), newZsur(ind3), 'ko', 'MarkerFaceColor','b', 'MarkerSize', 15)
% plot3(newXsur(ind4), newYsur(ind4), newZsur(ind4), 'ko', 'MarkerFaceColor','b', 'MarkerSize', 15)
% plot3(xsur(firstEdgeInd), ysur(firstEdgeInd), zsur(firstEdgeInd), 'co', 'MarkerFaceColor','c', 'MarkerSize', 15)
% plot3(xsur(secondEdgeInd), ysur(secondEdgeInd), zsur(secondEdgeInd), 'co', 'MarkerFaceColor','c', 'MarkerSize', 15)
% plot3(xsur(thirdEdgeInd), ysur(thirdEdgeInd), zsur(thirdEdgeInd), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 15)
% plot3(xsur(fourthEdgeInd), ysur(fourthEdgeInd), zsur(fourthEdgeInd), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 15)
% plot3(xsur(fifthEdgeInd), ysur(fifthEdgeInd), zsur(fifthEdgeInd), 'yo', 'MarkerFaceColor','y', 'MarkerSize', 15)
% plot3(xsur(sixthEdgeInd), ysur(sixthEdgeInd), zsur(sixthEdgeInd), 'yo', 'MarkerFaceColor','y', 'MarkerSize', 15)
% plot3(xsur(seventhEdgeInd), ysur(seventhEdgeInd), zsur(seventhEdgeInd), 'go', 'MarkerFaceColor','g', 'MarkerSize', 15)
% plot3(xsur(eigthEdgeInd), ysur(eigthEdgeInd), zsur(eigthEdgeInd), 'go', 'MarkerFaceColor','g', 'MarkerSize', 15)
% xlabel('x')
% ylabel('y')
% zlabel('z')


%now that edges are obtained, all the points within those edges will be
%found by "sweeping" from side to side in two different directions.


finEdge1 = [firstEdgeInd; secondEdgeInd];
finEdge2 = [thirdEdgeInd; fourthEdgeInd];
finEdge3 = [fifthEdgeInd; sixthEdgeInd];
finEdge4 = [seventhEdgeInd; eigthEdgeInd];

insideEdgeInd = [];
% numTestElem = 9;


%Sweep from left to right
for i = 1:length(finEdge4)

    firstCorner = finEdge4(i);
    if i <= length(finEdge2)
        secondCorner = finEdge2(i);
    else
        secondCorner = finEdge2(end);
    end
    locInsideInd = [];
    [indFirst,dists] = findNearestNeighbors(ptCloudSur,[ xsur(firstCorner), ysur(firstCorner), zsur(firstCorner) ],1);
    locInsideInd = [ locInsideInd; indFirst(1) ];
    [indLast,dists] = findNearestNeighbors(ptCloudSur,[ xsur(secondCorner), ysur(secondCorner), zsur(secondCorner) ],1);
    locInsideIndTot = [];

    dummInd = 1;
    lastDistance = 10000;
    currentDistance = 100;
    while currentDistance < lastDistance
        disp("Moving through..."+num2str(dummInd));
        [indices,dists] = findNearestNeighbors(ptCloudSur,[ xsur(locInsideInd(dummInd)), ysur(locInsideInd(dummInd)), zsur(locInsideInd(dummInd)) ], numTestElem);
        locVecD = [];
        for i=2:length(indices)
            locVecD(i-1) = sqrt( ( xsur( indLast ) - xsur( indices(i) ) )^2 + ( ysur( indLast ) - ysur( indices(i) ) )^2 + ( zsur( indLast ) - zsur( indices(i) ) )^2 );
        end
        nextInd = find(locVecD == min(locVecD)) + 1;
        locInsideInd = [ locInsideInd; indices(nextInd(1)) ];
        locInsideIndTot = [locInsideIndTot; indices];
        dummInd = dummInd + 1;
        lastDistance = sqrt( ( xsur( indLast ) - xsur( locInsideInd(dummInd-1) ) )^2 + ( ysur( indLast ) - ysur( locInsideInd(dummInd-1) ) )^2 + ( zsur( indLast ) - zsur( locInsideInd(dummInd-1) ) )^2 );
        currentDistance = sqrt( ( xsur( indLast ) - xsur( locInsideInd(dummInd) ) )^2 + ( ysur( indLast ) - ysur( locInsideInd(dummInd) ) )^2 + ( zsur( indLast ) - zsur( locInsideInd(dummInd) ) )^2 );
    end

    insideEdgeInd = [insideEdgeInd; locInsideIndTot];

end




%Sweep from bottom to top
for i = 1:length(finEdge3)

    firstCorner = finEdge3(i);
    if i <= length(finEdge1)
        secondCorner = finEdge1(i);
    else
        secondCorner = finEdge1(end);
    end
    locInsideInd = [];
    [indFirst,dists] = findNearestNeighbors(ptCloudSur,[ xsur(firstCorner), ysur(firstCorner), zsur(firstCorner) ],1);
    locInsideInd = [ locInsideInd; indFirst(1) ];
    [indLast,dists] = findNearestNeighbors(ptCloudSur,[ xsur(secondCorner), ysur(secondCorner), zsur(secondCorner) ],1);
    locInsideIndTot = [];

    dummInd = 1;
    lastDistance = 10000;
    currentDistance = 100;
    while currentDistance < lastDistance
        disp("Moving through..."+num2str(dummInd));
        [indices,dists] = findNearestNeighbors(ptCloudSur,[ xsur(locInsideInd(dummInd)), ysur(locInsideInd(dummInd)), zsur(locInsideInd(dummInd)) ], numTestElem);
        locVecD = [];
        for i=2:length(indices)
            locVecD(i-1) = sqrt( ( xsur( indLast ) - xsur( indices(i) ) )^2 + ( ysur( indLast ) - ysur( indices(i) ) )^2 + ( zsur( indLast ) - zsur( indices(i) ) )^2 );
        end
        nextInd = find(locVecD == min(locVecD)) + 1;
        locInsideInd = [ locInsideInd; indices(nextInd(1)) ];
        locInsideIndTot = [locInsideIndTot; indices];
        dummInd = dummInd + 1;
        lastDistance = sqrt( ( xsur( indLast ) - xsur( locInsideInd(dummInd-1) ) )^2 + ( ysur( indLast ) - ysur( locInsideInd(dummInd-1) ) )^2 + ( zsur( indLast ) - zsur( locInsideInd(dummInd-1) ) )^2 );
        currentDistance = sqrt( ( xsur( indLast ) - xsur( locInsideInd(dummInd) ) )^2 + ( ysur( indLast ) - ysur( locInsideInd(dummInd) ) )^2 + ( zsur( indLast ) - zsur( locInsideInd(dummInd) ) )^2 );
    end

    insideEdgeInd = [insideEdgeInd; locInsideIndTot];

end


%putting it all together
insideEdgeInd = unique(insideEdgeInd);

figure;
% plot3(newXsur, newYsur, newZsur, 'ko', 'MarkerFaceColor','r', 'MarkerSize', 10)
% plot3(actxsur, actysur, actzsur, 'mo', 'MarkerFaceColor','m', 'MarkerSize', 10)
hold on
plot3(xsur, ysur, zsur, 'ko', 'MarkerFaceColor','r', 'MarkerSize', 8)
plot3(xsur(insideEdgeInd), ysur(insideEdgeInd), zsur(insideEdgeInd), 'mo', 'MarkerFaceColor','m', 'MarkerSize', 15)
plot3(xsur(finEdge1), ysur(finEdge1), zsur(finEdge1), 'co', 'MarkerFaceColor','c', 'MarkerSize', 15)
plot3(xsur(finEdge2), ysur(finEdge2), zsur(finEdge2), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 15)
plot3(xsur(finEdge3), ysur(finEdge3), zsur(finEdge3), 'yo', 'MarkerFaceColor','y', 'MarkerSize', 15)
plot3(xsur(finEdge4), ysur(finEdge4), zsur(finEdge4), 'go', 'MarkerFaceColor','g', 'MarkerSize', 15)
plot3(newXsur(ind1), newYsur(ind1), newZsur(ind1), 'bo', 'MarkerFaceColor','b', 'MarkerSize', 20)
plot3(newXsur(ind3), newYsur(ind3), newZsur(ind3), 'bo', 'MarkerFaceColor','b', 'MarkerSize', 20)
plot3(newXsur(ind5), newYsur(ind5), newZsur(ind5), 'bo', 'MarkerFaceColor','b', 'MarkerSize', 20)
plot3(newXsur(ind7), newYsur(ind7), newZsur(ind7), 'bo', 'MarkerFaceColor','b', 'MarkerSize', 20)
xlabel('X'); ylabel('Y'); zlabel('Z')
set(gca,'FontSize',35)



allLAPWInd = [insideEdgeInd; finEdge1; finEdge2; finEdge3; finEdge4];
allLAPWInd = unique(allLAPWInd);



end