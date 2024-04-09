%Map data from original mesh to generated new meshes and create subdomains
%Three sections: one section to do tissue-only mapping, another one to do
%tissue-bath mapping and subdomain determining, and another one to do
%multiple regions and fiber mapping
%Final extra section to save the CSV's of all of these

%NOTE: RIGHT NOW TO GENERATE THE DIFFERENT CASES FOR THE MESHES, WE ARE
%BUILDING ON TOP OF THE PREVIOUS CASE: FROM TISSUE TO TISSUE WITH BATH,
%FROM TISSUE WITH BATH TO TISSUE WITH BATH AND FAT... AND IT COULD KEEP
%GOING ON

%Section for tissue-only mapping
%original points from paraview and original mesh

lapwOG = csvread('../SimResults/P1simsStuff/DataAll_M75_transformed_Persistent8.csv'); %This is the CSV data of the clipped region from ParaView containing EVERYTHING and all the points NOT JUST SURFACE

xog = lapwOG(:,12);
yog = lapwOG(:,13);
zog = lapwOG(:,14);

fibx = lapwOG(:,1);
fiby = lapwOG(:,2);
fibz = lapwOG(:,3);

sheetx = lapwOG(:,4);
sheety = lapwOG(:,5);
sheetz = lapwOG(:,6);

crossfx = lapwOG(:,7);
crossfy = lapwOG(:,8);
crossfz = lapwOG(:,9);

% scatter3(xog,yog,zog,14,crossfz,'filled');
% colorbar

%refined points to map
lapwNM = csvread('../SimResults/P1simsStuff/Persistent8_AllPoints_Refined_second.csv'); %this is the CSV data from the REFINED libMesh mesh containing ALL the points that are to be mapped with fiber stuff.
%MAKE SURE TO USE THE REFINED MESH PRODUCT OF LIBMESH REFINEMENT TO GET
%EXACT NODES OF THE NEW MESH THAT WILL BE USED TO RUN THE ACTUAL SIMS
%THIS APPLIES FOR ALL THE MESHES: TISSUE ONLY, TISSUE-BATH, ETC.
xnm = lapwNM(:,12);
ynm = lapwNM(:,13);
znm = lapwNM(:,14);










plot3(xnm,ynm,znm,'bo','MarkerSize',7)
hold on
plot3(xog,yog,zog,'ro','MarkerSize',7)
%plotting to check how the points look like

%actual mapping for tissue
nmfibx = [];
nmfiby = [];
nmfibz = [];
nmsheetx = [];
nmsheety = [];
nmsheetz = [];
nmcrossfx = [];
nmcrossfy = [];
nmcrossfz = [];

maxDist = 0.;
dataSetCloud = pointCloud([xog, yog, zog]);
farPointIndexNM = 0;
farPointIndexOG = 0;

%mapping fiber data from original clipped data from patient meshes to new
%generated mesh based on closest point to the new mesh points
for i = 1:length(xnm)
    [indices,dists] = findNearestNeighbors(dataSetCloud,[xnm(i), ynm(i), znm(i)],1);
    nmfibx(i) = fibx(indices);
    nmfiby(i) = fiby(indices);
    nmfibz(i) = fibz(indices);

    nmsheetx(i) = sheetx(indices);
    nmsheety(i) = sheety(indices);
    nmsheetz(i) = sheetz(indices);

    nmcrossfx(i) = crossfx(indices);
    nmcrossfy(i) = crossfy(indices);
    nmcrossfz(i) = crossfz(indices);

    locDistPoints = sqrt( (xnm(i) - xog(indices)).^2 + (ynm(i) - yog(indices)).^2 + (znm(i) - zog(indices)).^2 );
    %this part is done to know where the max distance from original to new data is in the new mesh that will have bath
    %and other layers as well. If a point in the cloud of points of the
    %tissue with bath points is farther than this maxDist, the it is not
    %tissue any more but bath or anything else. This is done later, but
    %that maxDist threshold is calculated here with the OG clipped points
    %and the new tissueOnly mesh. 
    if locDistPoints > maxDist
        maxDist = locDistPoints;
        farPointIndexNM = i;
        farPointIndexOG = indices;
    end

end


figure
subplot(2,1,1);
scatter3(xnm,ynm,znm,14,nmfibz,'filled');
subplot(2,1,2);
scatter3(xog,yog,zog,14,fibz,'filled');




 

%%
%Section for Tissue-Bath mapping


lapwNMBath = csvread('M71data/M71_LAPW_tissueBath_refinedPoints.csv');
%this is the CSV data from the REFINED libMesh mesh containing ALL the points
% %for BOTH tissue AND something else (i.e. bath) that are to be mapped with
% %fiber stuff and subdomain ID.

xnmb = lapwNMBath(:,2);
ynmb = lapwNMBath(:,3);
znmb = lapwNMBath(:,4);

nmbfibx = [];
nmbfiby = [];
nmbfibz = [];
nmbsheetx = [];
nmbsheety = [];
nmbsheetz = [];
nmbcrossfx = [];
nmbcrossfy = [];
nmbcrossfz = [];

plot3(xnm,ynm,znm,'ro','MarkerSize',7)
hold on
plot3(xnmb,ynmb,znmb,'bo','MarkerSize',7)

%like libMesh: 1 tissue, 2 bath, fat = 3, fibrosis = 0
SubdomainIDInd = [];
dataSetCloudBath = pointCloud([xnm, ynm, znm]);
tissueCounter = 1;
tissueIndices = [];
bathCounter = 1;
bathIndices = [];
maxDistTB = 0.;

%Mapping of fiber stuff to new mesh as well as subdomain checking
for i = 1:length(xnmb)
    [indices,dists] = findNearestNeighbors(dataSetCloudBath,[xnmb(i), ynmb(i), znmb(i)],1);
    
    locDistPoints = sqrt( (xnmb(i) - xnm(indices)).^2 + (ynmb(i) - ynm(indices)).^2 + (znmb(i) - znm(indices)).^2 );
    
    if locDistPoints > maxDist*.5 %We are in the bath region or fat or anything else but tissue
        SubdomainIDInd(i) = 2;
        bathIndices(bathCounter) = i; 
        bathCounter = bathCounter + 1;
    else %We are in the tissue region
        SubdomainIDInd(i) = 1;
        tissueIndices(tissueCounter) = i; 
        [indices2,dists2] = findNearestNeighbors(dataSetCloud,[xnmb(i), ynmb(i), znmb(i)],1);
        nmbfibx(tissueCounter) = fibx(indices2);
        nmbfiby(tissueCounter) = fiby(indices2);
        nmbfibz(tissueCounter) = fibz(indices2);
    
        nmbsheetx(tissueCounter) = sheetx(indices2);
        nmbsheety(tissueCounter) = sheety(indices2);
        nmbsheetz(tissueCounter) = sheetz(indices2);
    
        nmbcrossfx(tissueCounter) = crossfx(indices2);
        nmbcrossfy(tissueCounter) = crossfy(indices2);
        nmbcrossfz(tissueCounter) = crossfz(indices2);
        tissueCounter = tissueCounter + 1;

        if locDistPoints > maxDistTB
            maxDistTB = locDistPoints;
        end

    end

end
% 
% scatter3(xnmb,ynmb,znmb,14,SubdomainIDInd,'filled');
% hold on
% plot3(xnm,ynm,znm, 'ro')

xnmactualBath = xnmb(bathIndices);
ynmactualBath = ynmb(bathIndices);
znmactualBath = znmb(bathIndices);
%this points will be used to actually get to identify the bath in the cases
%where there is more than just tissue and bath (like fat). Creat pointCloud
%of this and find nearest neighbor in radius to differentiate between fat
%and bath.

plot3(xnmb(tissueIndices),ynmb(tissueIndices),znmb(tissueIndices), 'ro')
hold on
plot3(xnmb(bathIndices),ynmb(bathIndices),znmb(bathIndices), 'bo')
% plot3(xnm,ynm,znm, 'ro')
plot3(xnmb,ynmb,znmb, 'ko')

subplot(2,1,1);
scatter3(xog,yog,zog,14,fibx,'filled');
subplot(2,1,2);
scatter3(xnmb(tissueIndices),ynmb(tissueIndices),znmb(tissueIndices),14,nmbfibx,'filled');
% 






%%
%Section for tissue and other regions not only bath mapping





lapwNMBathFat = csvread('M75data/tissueBathFat/M75RefinedPointsMatlabTissueBathFat.csv');
%this is the CSV data from the REFINED libMesh mesh containing ALL the points
% %for BOTH tissue AND something else (i.e. bath) that are to be mapped with
% %fiber stuff and subdomain ID.

xnmbf = lapwNMBathFat(:,2);
ynmbf = lapwNMBathFat(:,3);
znmbf = lapwNMBathFat(:,4);

nmbffibx = [];
nmbffiby = [];
nmbffibz = [];
nmbfsheetx = [];
nmbfsheety = [];
nmbfsheetz = [];
nmbfcrossfx = [];
nmbfcrossfy = [];
nmbfcrossfz = [];

% plot3(xnm,ynm,znm,'ro','MarkerSize',7)
% hold on
% plot3(xnmbf,ynmbf,znmbf,'bo','MarkerSize',7)

%like libMesh: 1 tissue, 2 bath, fat = 3, fibrosis = 0
SubdomainIDInd = [];
dataSetCloudBath = pointCloud([xnm, ynm, znm]);
tissueCounter = 1;
tissueIndices = [];
bathCounter = 1;
bathIndices = [];
fatCounter = 1;
fatIndices = [];
radiusOfBathPoints = .18;
% maxDistTB = 0.;

dataSetCloudActualBath = pointCloud([xnmactualBath, ynmactualBath, znmactualBath]);

%Mapping of fiber stuff to new mesh as well as subdomain checking
for i = 1:length(xnmbf)
    [indices,dists] = findNearestNeighbors(dataSetCloud,[xnmbf(i), ynmbf(i), znmbf(i)],1);
    
%     locDistPoints = sqrt( (xnmbf(i) - xnm(indices)).^2 + (ynmbf(i) - ynm(indices)).^2 + (znmbf(i) - znm(indices)).^2 );
    
    if dists > maxDistTB %We are in the bath region or fat or anything else but tissue

        %Use the previous mesh of only bath to check for points close to
        %this in a radius of the bath of the previous mesh. If it is then,
        %the subdomain is bath, if not then it's bath.

        [idx3,dists3] = findNeighborsInRadius(dataSetCloudActualBath,[xnmbf(i), ynmbf(i), znmbf(i)],radiusOfBathPoints);

        if(isempty(idx3))
            SubdomainIDInd(i) = 3;
            fatIndices(fatCounter) = i; 
            fatCounter = fatCounter + 1;
        else
            SubdomainIDInd(i) = 2;
            bathIndices(bathCounter) = i; 
            bathCounter = bathCounter + 1;
        end
        
    else %We are in the tissue region
        SubdomainIDInd(i) = 1;
        tissueIndices(tissueCounter) = i; 
        [indices2,dists2] = findNearestNeighbors(dataSetCloud,[xnmbf(i), ynmbf(i), znmbf(i)],1);
        nmbffibx(tissueCounter) = fibx(indices2);
        nmbffiby(tissueCounter) = fiby(indices2);
        nmbffibz(tissueCounter) = fibz(indices2);
    
        nmbfsheetx(tissueCounter) = sheetx(indices2);
        nmbfsheety(tissueCounter) = sheety(indices2);
        nmbfsheetz(tissueCounter) = sheetz(indices2);
    
        nmbfcrossfx(tissueCounter) = crossfx(indices2);
        nmbfcrossfy(tissueCounter) = crossfy(indices2);
        nmbfcrossfz(tissueCounter) = crossfz(indices2);
        tissueCounter = tissueCounter + 1;

%         if locDistPoints > maxDistTB
%             maxDistTB = locDistPoints;
%         end

    end

end

% subplot(2,1,1);
% scatter3(xnmbf(tissueIndices),ynmbf(tissueIndices),znmbf(tissueIndices),14,nmbfcrossfz,'filled');
% subplot(2,1,2);
% scatter3(xog,yog,zog,14,crossfz,'filled');

% plot3(xnmbf(tissueIndices),ynmbf(tissueIndices),znmbf(tissueIndices),'ro');
scatter3(xnmbf(tissueIndices),ynmbf(tissueIndices),znmbf(tissueIndices),14,nmbfcrossfz,'filled');
hold on;
plot3(xnmbf(bathIndices),ynmbf(bathIndices),znmbf(bathIndices),'bo');
plot3(xnmbf(fatIndices),ynmbf(fatIndices),znmbf(fatIndices),'ko');







%%
%extra section for saving the files into csv's that will be loaded in
%libMesh

%Writing files to be read in libMesh
numberofDigits = 8.;
multiplier_round = 10.^numberofDigits;
% identifierCol = 1:length(xnmEx);



%OnlyTissue
xnmEx = (floor(xnm.*multiplier_round))./multiplier_round;
ynmEx = (floor(ynm.*multiplier_round))./multiplier_round;
znmEx = (floor(znm.*multiplier_round))./multiplier_round;
matrixTissuePointsExport = [xnmEx, ynmEx, znmEx];
% matrixTissuePointsExportId = [identifierCol',xnmEx, ynmEx, znmEx];
dlmwrite('../SimResults/P1simsStuff/Persistent8_RealGeometry_tissueOnlyCoordPoints_second.csv', round(matrixTissuePointsExport,4), 'delimiter', ',', 'precision', '%.6f');
% dlmwrite('M75tissueOnlyCoordPointsId.csv', matrixTissuePointsExportId, 'delimiter', ',', 'precision', '%.6f');
matrixTissueFiberDirectionsExport = [nmfibx', nmfiby', nmfibz',nmsheetx', nmsheety', nmsheetz',nmcrossfx', nmcrossfy', nmcrossfz'];
dlmwrite('../SimResults/P1simsStuff/Persistent8_RealGeometry_tissueOnlyFiberDirections_second.csv', matrixTissueFiberDirectionsExport, 'delimiter', ',', 'precision', '%.6f');






%Tissue-Bath
xnmbEx = (floor(xnmb.*multiplier_round))./multiplier_round;
ynmbEx = (floor(ynmb.*multiplier_round))./multiplier_round;
znmbEx = (floor(znmb.*multiplier_round))./multiplier_round;
matrixTissueBathTissuePointsExport = [xnmbEx(tissueIndices), ynmbEx(tissueIndices), znmbEx(tissueIndices)];
dlmwrite('M75tissueBathCoordPoints.csv', round(matrixTissueBathTissuePointsExport,4), 'delimiter', ',', 'precision', '%.4f');
matrixTissueBathFiberDirectionsExport = [nmbfibx', nmbfiby', nmbfibz', nmbsheetx', nmbsheety', nmbsheetz',nmbcrossfx', nmbcrossfy', nmbcrossfz'];
dlmwrite('M75tissueBathFiberDirections.csv', matrixTissueBathFiberDirectionsExport, 'delimiter', ',', 'precision', '%.6f');






%Tissue-Bath_Fat
xnmbfEx = (floor(xnmbf.*multiplier_round))./multiplier_round;
ynmbfEx = (floor(ynmbf.*multiplier_round))./multiplier_round;
znmbfEx = (floor(znmbf.*multiplier_round))./multiplier_round;
matrixTissueBathFatTissuePointsExport = [xnmbfEx(tissueIndices), ynmbfEx(tissueIndices), znmbfEx(tissueIndices)];
dlmwrite('M75tissueBathFatCoordPoints.csv', round(matrixTissueBathFatTissuePointsExport,4), 'delimiter', ',', 'precision', '%.4f');
matrixTissueBathFatFiberDirectionsExport = [nmbffibx', nmbffiby', nmbffibz', nmbfsheetx', nmbfsheety', nmbfsheetz',nmbfcrossfx', nmbfcrossfy', nmbfcrossfz'];
dlmwrite('M75tissueBathFatFiberDirections.csv', matrixTissueBathFatFiberDirectionsExport, 'delimiter', ',', 'precision', '%.6f');
matrixTissueBathFatBathPointsExport = [xnmbfEx(bathIndices), ynmbfEx(bathIndices), znmbfEx(bathIndices)];
dlmwrite('M75tissueBathFatBathPoints.csv', round(matrixTissueBathFatBathPointsExport,4), 'delimiter', ',', 'precision', '%.4f');
matrixTissueBathFatFatPointsExport = [xnmbfEx(fatIndices), ynmbfEx(fatIndices), znmbfEx(fatIndices)];
dlmwrite('M75tissueBathFatFatPoints.csv', round(matrixTissueBathFatFatPointsExport,4), 'delimiter', ',', 'precision', '%.4f');



%%
%Whole LA stuff

LAWhole = csvread('M75data/M75WholeLAPoints.csv');
matrixWholeLAFiber = LAWhole(:,1:9);
matrixWholeLAPoints = LAWhole(:,12:14);
matrixWholeLAPoints = (floor(matrixWholeLAPoints.*multiplier_round))./multiplier_round;
dlmwrite('M75WholeLACoordPoints.csv', round(matrixWholeLAPoints,4), 'delimiter', ',', 'precision', '%.4f');
dlmwrite('M75WholeLAFibers.csv', matrixWholeLAFiber, 'delimiter', ',', 'precision', '%.6f');

%%
%grounded Nodes
groundcoords = csvread('../SimResults/P1simsStuff/SurfacePoints_ForEGMs_ParaView.csv');
groundcoordsX = groundcoords(:,4);
groundcoordsY = groundcoords(:,5);
groundcoordsZ = groundcoords(:,6);

groundcoordsXEx = (floor(groundcoordsX.*multiplier_round))./multiplier_round;
groundcoordsYEx = (floor(groundcoordsY.*multiplier_round))./multiplier_round;
groundcoordsZEx = (floor(groundcoordsZ.*multiplier_round))./multiplier_round;
matrixGroundedPointsExport = [groundcoordsXEx, groundcoordsYEx, groundcoordsZEx];
dlmwrite('../SimResults/P1simsStuff/SurfacePoints_ForEGMsV2.csv', matrixGroundedPointsExport, 'delimiter', ',', 'precision', '%.6f');


