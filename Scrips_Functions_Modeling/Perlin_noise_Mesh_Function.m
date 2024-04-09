%Function to create fibrosis patterns

function [meshX, meshY, meshZ,tissueInd,fibroticInd] = Perlin_noise_Mesh_Function(fibrosisDensity_desired, lb, gamma, R, phi_fiber_orientation, psi_fiberZ_orientation, d, ld, f, L, meshX, meshY, meshZ, TrialVal, fiberBool, nameOfFiberFile, nameOfFibrosisFile)
%function Perlin_noiseF(WhiteNoise, DimToUse, fibrosisDensity_desired, lb, gamma, R, phi_fiber_orientation, psi_fiberZ_orientation, d, ld, f, L, maxX, minX, maxY, minY, maxZ, minZ, numElX, numElY, numElZ,caseName)
%Jakes et. al. 2019 Perlin Noise
%CHOOSE DIMENSION TO USE: 2D OR 3D



%INPUTS:
%fibrosisDensity_desired - desired fibrosis density in percentage
%lb - Base size of fibrotic obstacles (mm)
%gamma - Roughness of base noise field
%R - Anisotropy of base obstacles MAYBE DOUBLE THEM SINCE NOW IT'S 3D AND NOW IT'S NOT SQRT(R) BUT R
%phi_fiber_orientation - first angle for fiber orientation when not using prescribed fibers
%psi_fiberZ_orientation - second angle for fiber orientation when not using prescribed fibers
%d - effect of density variation field in final pattern
%ld - Size of density variation features (mm)
%f - effect of fiber orientation field in final pattern
%L - distance between fibers (mm)
%meshX - x coordinate points
%meshY - y coordinate points
%meshZ - z coordinate points
%TrialVal - factor to make the noise field size bigger (3 used typically)
%fiberBool - boolean to determine if prescribing fiber or not
%nameOfFiberFile - string with the name of the file containing fiber orientation information to create fiber field
%nameOfFibrosisFile - to save coordinates of fibrosis points give a name of fibrosis file


%meshX - x coordinate points
%meshY - y coordinate points
%meshZ - z coordinate points
%tissueInd - indices of coordinate points that are tissue domain
%fibroticInd - indices of coordinate points that are fibrosis domain



%GIVE EVERYTHING IN MM

%PARAMETERS TO USE
%phi_fiber can change based on whatever desired for fiber orientation
%                     lb    gamma    R    phi_fiber    d     ld     f     L
%   Interstitial      .24    .96    1.89     34       .32   4.67   .30   .31
%   Compact           .96    .59    2.47     -9       .44   2.03    0     -
%   Diffuse           .07    .44    2.17     11       .49   1.22    0     -
%   Patchy            .32    .78    2.50     68       .42   2.10   .38   .31

%Density (according to table, but just suggested, also can vary based on
%desired density):
%Interstitial: 9.6%
%Compact: 47.2%
%Diffuse: 22.0%
%Patchy: 26.9%

%Base Noise field
%lb = .96;             %[.01, 2]   Base size of fibrotic obstacles (mm)
lb = lb/10;           %now in cm
%R = 2.47;                %[.5, 50]   Anisotropy of base obstacles MAYBE DOUBLE THEM SINCE NOW IT'S 3D AND NOW IT'S NOT SQRT(R) BUT R
%gamma = .59;           %[0, .99]   Roughness of base noise field

%Density Variation Noise field
%ld = 2.03;               %[1, 8]   Size of density variation features (mm)
ld = ld/10;           %now in cm

%Fiber selection Noise field
%L = .31;               %[.3, 2]   Base separation distance of fibres (mm)
L = L/10;    %now in cm
%phi_fiber_orientation = 11;                             %[-pi/2, pi/2]   Fibre orientation in degrees
phi_fiber_orientation = phi_fiber_orientation*pi/180;          %now in radians

%Others
%f = 0.;                     %[0, .4]   Extent of pattern coincidence with fibres
%d = .44;                    %[0, .5]   Magnitude of local density variation

Seed = 1;

%Fixed: Fiber selection Noise field
s = 15;   %Sharpening factor
phi_phase_modulation = 10*pi;   %Strength of phase modulation
Beta_dissimilarity = 1;%.25;   %Dissimilarity between fibres
Beta_wiggliness = 1;%.25;   %Fibre ?wiggliness?

% Beta_dissimilarity = Beta_dissimilarity*1;%*5
% Beta_wiggliness = Beta_wiggliness*1;%*60

Nf = 4; %number of octaves of Perlin noise

%Geometry parameters
%maxX = 1.; %cm
%minX = -1.; %cm
%maxY = 1.; %cm
%minY = -1.; %cm

desiredH = .025;%cm

%numElX = 80;
%numElY = 80;


numElTot = length(meshX);
numElPDim = ceil((numElTot^(1/3))*3.25);
maxX = max(meshX);
maxY = max(meshY);
maxZ = max(meshZ);
minX = min(meshX);
minY = min(meshY);
minZ = min(meshZ);


width = numElPDim;
height = numElPDim;

fiberInfoMat = [];

if fiberBool
    fiberInfoMat = csvread(nameOfFiberFile);
    disp('Using prescribed fibers...reading file...')
end

%numEl = ceil((maxX-minX)/(desiredH));

%if 3D:
%maxZ = .1; %cm
%minZ = -.1; %cm
%numElZ = 8;
base = numElPDim;

numPlanes = 6; %number of planes to show in images
%Extra angle parameter for 3D
%psi_fiberZ_orientation = 30;                            %[-pi/2, pi/2]   Fibre orientation in degrees in z direction
psi_fiberZ_orientation = psi_fiberZ_orientation*pi/180;        %now in radians

%White noise parameter
% TrialVal = 10; %factor to use to make the white noise grid bigger so that when obtaining noise values the grid is large enough

%%
        %PARAMETERS TO USE 3D
        IMG_W = width;    
        IMG_H = height;
        IMG_B = base;

        A_R = [R.^(-2./6), 0,            0;...
               0,          R.^(1./6),    0;...
               0,          0,         R.^(1./6)];
        

        R_phi = [sin(psi_fiberZ_orientation)*cos(phi_fiber_orientation), sin(psi_fiberZ_orientation)*sin(phi_fiber_orientation), cos(psi_fiberZ_orientation);...
               cos(psi_fiberZ_orientation)*cos(phi_fiber_orientation), cos(psi_fiberZ_orientation)*sin(phi_fiber_orientation), -sin(psi_fiberZ_orientation);...
               -sin(phi_fiber_orientation), cos(phi_fiber_orientation), 0];

       
        %CREATING GEOMETRY
        dx = (maxX - minX)/numElPDim;
        dy = (maxY - minY)/numElPDim;
        dz = (maxZ - minZ)/numElPDim;

        if dx > .1 & dy > .1 & dz > .1
            disp('The mesh interval is not small enough...')
        end

%         numPlanesX = numElX/numPlanes;
%         numPlanesY = numElY/numPlanes;
%         numPlanesZ = numElZ/numPlanes;
        
%         FunGrid = ones(numElY,numElX,numElZ);
        
        
        limB = 1.1*(max([minX maxX minY maxY minZ maxZ]));
        
        
%         slice(X,Y,Z,FunGrid, minX+dx/2:dx*numPlanesX:maxX-dx/2, minY+dy/2:dy*numPlanesY:maxY-dy/2, minZ+dz/2:dz*numPlanesZ:maxZ-dz/2)
%         FunGridAlt = reshape(FunGrid,1,[]);



        
        %PLOTTING HERE
        % figure
        % plot3(meshX,meshY,meshZ,'ro')





%         set(gca,'clim',[min(FunGridAlt) max(FunGridAlt)])
%         axis([-limB limB -limB limB -limB limB])
        % title("Grid to use")


        %%


        %IMAGING WHITE NOISE 

        CMAP = [0 0 0 ; hsv(255)];%HSV colormap
        img = ones(IMG_H, IMG_W, IMG_B);
        % initiate figure
        

        % initiate image
        % hImg = image('XData',1:IMG_W, 'YData',1:IMG_H, 'ZData',1:IMG_B, ...
        %     'CData',img, 'CDataMapping','direct');

        
        WhiteNoise = zeros(IMG_W*TrialVal,IMG_H*TrialVal, IMG_B*TrialVal);

        % draw the whitenoise
        for x = 1:IMG_W*TrialVal
          for y = 1:IMG_H*TrialVal
              for z = 1:IMG_B*TrialVal
                clr = (1-rand()*2);
                WhiteNoise(x,y,z) = clr;
              end
          end
        end
        
        % randomVectorGrid = zeros(IMG_W*TrialVal,IMG_H*TrialVal, 2);
        %    
        % draw the whitenoise
        % for x = 1:IMG_W*TrialVal
        %   for y = 1:IMG_H*TrialVal
        %         v= 2*rand(2,1)-1;
        %         n=v/sqrt(v(1)^2+v(2)^2);
        %         randomVectorGrid(x,y,:) = n;
        %   end
        % end

        for x = 1:IMG_H
          for y = 1:IMG_W
              for z = 1:IMG_B
                img(x,y,z) = (WhiteNoise(x,y,z)+1)*125;
              end
          end
        end
        % set(hImg, 'CData',img)
        % drawnow  


%         slice(X,Y,Z,img, minX+dx/2:dx*2:maxX-dx/2, minY+dy/2:dy*2:maxY-dy/2, minZ+dz/2:dz*2:maxZ-dz/2)
        imgAlt = reshape(img,1,[]);
%         scatter3(xAlt,yAlt,zAlt,20,imgAlt)





        %PLOTTING HERE
%         set(gca,'clim',[min(imgAlt) max(imgAlt)])
% %         axis([-limB limB -limB limB -limB limB])
%         title("White Noise 3D Field")
%         %hold on
%         plot3(meshX,meshY,meshZ,'ro')


        %%



        %CREATING PERLIN NOISE FIELDS FOR THE THREE DIFFERENT NOISE FIELDS
        %JUST DO MAPPING FROM GRID TO NOISE FIELD
        % draw the perlin noise

        %three noisefields:
        % - Nb_field: Fibrotic Obstacles
        % - Nd_field: Density Variation
        % - F_field: Fiber Selection



        % Nb_perlinMatrix = zeros(IMG_H, IMG_W, IMG_B);
        % Nd_perlinMatrix = zeros(IMG_H, IMG_W, IMG_B);
        % F_perlinMatrix = zeros(IMG_H, IMG_W, IMG_B);

        Nb_field_img = ones(numElTot,1);
        Nd_field_img = ones(numElTot,1);
        F_field_img = ones(numElTot,1);

        % I_perlinMatrix = zeros(IMG_H, IMG_W);
        I_final = ones(numElTot,1);

        % [gX,gY] = gradient(WhiteNoise);
        % gX = rescale(gX); 
        % gY = rescale(gY); 



        for x = 1:numElTot
            %for y = 1:IMG_W
                %for z = 1:IMG_B
                    %N_b BASE PERLIN NOISEFIELD
                    % offsetV1 = -.5 + (.5+.5) .* rand(1,1);
                    %perlin = perlinNoise2d(x, y, 0.1, 4, IMG_W, IMG_H, WhiteNoise);
                    xV_Nb = (1./lb).*([meshX(x), meshY(x), meshZ(x)]*A_R);
            %         xV_Nb = (1./lb).*([x, y]*A_R);
                    perlinNb = perlinNoise3d(xV_Nb(1), xV_Nb(2), xV_Nb(3), gamma, Nf, IMG_H, IMG_W, IMG_B, WhiteNoise);
                    %perlin = perlinNoise2d(50*X(x,y)+2, 10*Y(x,y)+2, 0.1, 4, IMG_W, IMG_H, WhiteNoise);
            %         Nb_perlinMatrix(x,y) = perlinNb;
                    cRGB_Nb = floor(1+ (perlinNb - 0) / (1 - 0) * (256-1)); 
                    %cRGB_Nb = floor(5.3857 + (perlinNb - (-4.3230)) / (5.3857 - (-4.3230)) * (256-1)); 
                    Nb_field_img(x) = cRGB_Nb;


                    %N_d DENSITY VARIATION PERLIN NOISEFIELD
                    xV_Nd = (1./ld).*([meshX(x), meshY(x), meshZ(x)]);
                    perlinNd = perlinNoise3d(xV_Nd(1), xV_Nd(2), xV_Nd(3), .5, Nf, IMG_H, IMG_W, IMG_B, WhiteNoise);
            %         Nd_perlinMatrix(x,y) = perlinNd;
                    cRGB_Nd = floor(1+ (perlinNd - 0) / (1 - 0) * (256-1)); 
                    %cRGB_Nd = floor(1.7932 + (perlinNd - (-1.3009)) / (1.7932 - -1.3009) * (256-1));
                    Nd_field_img(x) = cRGB_Nd;


                    %F FIBER SELECTION NOISEFIELD

                    %With rotation included
                    if fiberBool
                        fx = fiberInfoMat(x,1);
                        fy = fiberInfoMat(x,2);
                        fz = fiberInfoMat(x,3);

                        sx = fiberInfoMat(x,4);
                        sy = fiberInfoMat(x,5);
                        sz = fiberInfoMat(x,6);

                        cx = fiberInfoMat(x,7);
                        cy = fiberInfoMat(x,8);
                        cz = fiberInfoMat(x,9);

                        % phi_fiber_orientation = atan(fz/fx);
                        % psi_fiberZ_orientation = atan(fy/fx);
                        % if fy > fz
                        %     if fy > fx
                        %         phi_fiber_orientation = atan(fy/fz);
                        %         psi_fiberZ_orientation = atan(fx/fz);
                        %     end
                        % end
                        % if fx > fz
                        %     if fx > fy
                        %         phi_fiber_orientation = atan(fx/fy);
                        %         psi_fiberZ_orientation = atan(fz/fy);
                        %     end
                        % end

                        %TAKING ADVANTAGE OF PRESCRIBED FIBERS
                        R_phi = [fx, sx, cx; fy, sy, cy; fz, sz, cz];
                        
                        % R_phi = [sin(psi_fiberZ_orientation)*cos(phi_fiber_orientation), sin(psi_fiberZ_orientation)*sin(phi_fiber_orientation), cos(psi_fiberZ_orientation);...
                        % cos(psi_fiberZ_orientation)*cos(phi_fiber_orientation), cos(psi_fiberZ_orientation)*sin(phi_fiber_orientation), -sin(psi_fiberZ_orientation);...
                        % -sin(phi_fiber_orientation), cos(phi_fiber_orientation), 0];
                    end
                    % xV_F = (( [ meshX(x), meshY(x), meshZ(x) ] * (R_phi)) * [ Beta_wiggliness/L, 0, 0; 0, (Beta_dissimilarity/(L)), 0; 0, 0, (Beta_wiggliness/L) ] );


                    % Z_field = (phi_phase_modulation).*perlinNoise3d(xV_F(1), xV_F(2), xV_F(3), .5, Nf, IMG_H, IMG_W, IMG_B, WhiteNoise);

                    % perlinF = (.5 + .5*cos( (2*pi/L) * ( xV_F(2) ) + Z_field ) ).^2;
                    % perlinF = ((.5 + .5*cos( (2*pi/L) * ( xV_F(2) ) ) ).^2 + (.5 + .5*cos( (2*pi/L) * ( xV_F(1) ) ) ).^2 + (.5 + .5*cos( (2*pi/L) * ( xV_F(3) ) ) ).^2)./3;
                    % perlinF = (.5 + .5*cos( (2*pi/L) * ( xV_F(3) ) ) ).^2;

                    dotProduct = dot( [ meshX(x), meshY(x), meshZ(x) ], [-fy,fx,fz]./norm([-fy,fx,fz]) );
                    perlinF  = mod(dotProduct, L + L ) <= L;

                    % display('different')
            %       F_perlinMatrix(x,y) = perlinF;
                    %cRGB_F = floor(1+ (perlinF - 0) / (1 - 0) * (256-1));
                    cRGB_F = floor(1+ (perlinF - 0) / (1 - 0) * (256-1));
                    F_field_img(x) = cRGB_F;


                    %FINAL IMAGE (DETERMINE IF ROTATION OR NO ROTATION)
                    %if rotation, rotate first all the points in X,Y
                    I_perlin = ((perlinNb)*((1 - f) + (f*perlinF))) + d*perlinNd;
                    cRGB_I = floor(1+ (I_perlin - 0) / (1 - 0) * (256-1)); 
                    I_final(x) = cRGB_I;

                %end
            %end
        end
        % 
        % Nb_perlinMatrix = rescale(Nb_perlinMatrix); 
        % Nd_perlinMatrix = rescale(Nd_perlinMatrix); 
        % F_perlinMatrix = rescale(F_perlinMatrix); 
        % 
        % for x = 1:IMG_W
        %     for y = 1:IMG_H
        %         cRGB_Nb = floor(1+ (Nb_perlinMatrix(x,y) - 0) / (1 - 0) * (256-1)); 
        %         %cRGB_Nb = floor(5.3857 + (perlinNb - (-4.3230)) / (5.3857 - (-4.3230)) * (256-1)); 
        %         Nb_field_img(x,y) = cRGB_Nb;
        %         
        %         
        %         
        %         cRGB_Nd = floor(1+ (Nd_perlinMatrix(x,y) - 0) / (1 - 0) * (256-1)); 
        %         %cRGB_Nd = floor(1.7932 + (perlinNd - (-1.3009)) / (1.7932 - -1.3009) * (256-1));
        %         Nd_field_img(x,y) = cRGB_Nd;
        %         
        %         
        %         
        %         %cRGB_F = floor(1+ (perlinF - 0) / (1 - 0) * (256-1));
        %         cRGB_F = floor(1+ (F_perlinMatrix(x,y) - 0) / (1 - 0) * (256-1));
        %         F_field_img(x,y) = cRGB_F;
        %         
        %         
        %         %FINAL IMAGE (DETERMINE IF ROTATION OR NO ROTATION)
        %         %if rotation, rotate first all the points in X,Y
        %         I_perlinMatrix(x,y) = ((Nb_perlinMatrix(x,y))*((1 - f) + (f*F_perlinMatrix(x,y)))) + d*Nd_perlinMatrix(x,y);
        %         cRGB_I = floor(1+ (I_perlinMatrix(x,y) - 0) / (1 - 0) * (256-1)); 
        %         I_final(x,y) = cRGB_I;
        %         
        %         
        %     end
        % end




        %DETERMINING THRESHOLD FOR FIBROSIS IMAGING
        % meanV = mean(mean(mean(I_final)));
        % stdV = std(std(I_final));
        % maxV = max(max(I_final));
        % minV = min(min(I_final));

        colorFibrosis = 46;% YELLOW COLOR
        colorTissue = 18;% BLUE COLOR

        %fibroticThreshold = 100;
        fibroticThreshold = prctile(I_final,(100-fibrosisDensity_desired));
        fibCount = 0;
        totCount = 0;

        %histogram(I_final, 20)

        for x = 1:numElTot
            %for y = 1:IMG_W
                %for z = 1:IMG_B
                    if I_final(x) >= fibroticThreshold
                        I_final(x) = colorFibrosis;
                        fibCount = fibCount + 1;
                        totCount = totCount + 1;
                    else
                        I_final(x) = colorTissue;
                        totCount = totCount + 1;
                    end
                %end
            %end
        end

        fibrosis_density_Produced = fibCount/totCount


        %%


        fibroticInd = find(I_final == colorFibrosis);
        tissueInd = find(I_final == colorTissue);




% 
%         % FIGURES PLOTTING
%         %PLOTTING HERE
%         % initiate figure for Nb
%         hFig2 = figure('Pointer','crosshair', 'DoubleBuffer','on');
%         % initiate axis and colormap
%         % initiate image
%         Nb_field_imgAlt = Nb_field_img;%reshape(Nb_field_img,1,[]);
%         scatter3(meshX,meshY,meshZ,20,Nb_field_imgAlt)
%         set(gca,'clim',[min(Nb_field_imgAlt) max(Nb_field_imgAlt)])
% %         axis([-limB limB -limB limB -limB limB])
%         title("Base Noise field in 3D with points")
% 
%         % initiate figure for Nd
%         hFig3 = figure('Pointer','crosshair', 'DoubleBuffer','on');
%         % initiate axis
%         % initiate image
%         Nd_field_imgAlt = Nd_field_img;%reshape(Nd_field_img,1,[]);
%         scatter3(meshX,meshY,meshZ,20,Nd_field_imgAlt)
%         set(gca,'clim',[min(Nd_field_imgAlt) max(Nd_field_imgAlt)])
% %         axis([-limB limB -limB limB -limB limB])
%         title("Density Variation field in 3D with points")
% 
%         % initiate figure for F
%         hFig4 = figure('Pointer','crosshair', 'DoubleBuffer','on');
%         % initiate axis and colormap
%         % initiate image
%         F_field_imgAlt = F_field_img;%reshape(F_field_img,1,[]);
%         scatter3(meshX,meshY,meshZ,30,F_field_imgAlt)
% %         set(gca,'clim',[min(F_field_imgAlt) max(F_field_imgAlt)])
% %         axis([-limB limB -limB limB -limB limB])
%         title("Fiber Direction field in 3D with points")
% 
% 
% 
%         % % initiate figure for final Image
%         % hFig5 = figure('Pointer','crosshair', 'DoubleBuffer','on');
%         % % initiate axis and colormap
%         % % initiate image
%         %plotting in grid itself


%         %PLOTTING HERE
%         figure
%         plot3(meshX(tissueInd),meshY(tissueInd),meshZ(tissueInd),"oy")
%         hold on
%         plot3(meshX(fibroticInd),meshY(fibroticInd),meshZ(fibroticInd),"or")
%         grid on
% %         axis([-limB limB -limB limB -limB limB])
%         title("Image of Fibrosis (red) and Tissue (yellow) in points")
%         saveas(hFig5,'FibrosisImage.png');


        %%


        %EXPORTING GRID POINTS
        
        %EXPORTING CSV TO SEND TO LIBMESH CODE AND READ EACH POINT OF INTEREST
        %(WHICH IS EACH LINE OF EACH CSV FILE)
        %coordsTissue = [round(X(tissueInd),2), round(Y(tissueInd),2)];
        %coordsFibrosis = [round(X(fibroticInd),2), round(Y(fibroticInd),2)];
        numberofDigits = 8.;
        multiplier_round = 10.^numberofDigits;
        % coordsTissue = [ (floor(X(tissueInd).*multiplier_round))./multiplier_round, (floor(Y(tissueInd).*multiplier_round))./multiplier_round ];
        
        
        coordsFibrosis = [ (floor(meshX(fibroticInd).*multiplier_round))./multiplier_round, (floor(meshY(fibroticInd).*multiplier_round))./multiplier_round, (floor(meshZ(fibroticInd).*multiplier_round))./multiplier_round ];
        
        
        %dlmwrite('coords_Tissue.csv', coordsTissue, 'delimiter', ',', 'precision', '%.6f');
        dlmwrite(nameOfFibrosisFile, coordsFibrosis, 'delimiter', ',', 'precision', '%.6f');

        % these indices are the ones that will be used to set the specific
        %conductivities differently.


        
%%







%FUNCTIONS TO BE USED

%SHARED FUNCTION
 function r = cosineInterpolate(a, b, x)
    ft = x*pi;
    f1 = (1-cos(ft))*0.5;
    r = (a*(1-f1)+b*f1);
 end

%2D FUNCTIONS
function re = smoothedNoise2d(x,y, IMG_W, IMG_H, WhiteNoise)
    
    centerX = floor(IMG_W / 2);
    centerY = floor(IMG_H / 2);
    
    x = centerX + x;
    y = centerY + y;
    
    if x < 0
        x = -x;
%     elseif x >= IMG_W
%         x = IMG_W-1;
    end
    if y < 0
        y = -y;
%     elseif y >= IMG_H
%         y = IMG_H-1;
    end
    
    if x == 1 || x == 0
        x = 2;
    end
    if y == 1 || y == 0
        y = 2;
    end
    
    
%     x
%     y


    
    corners = (WhiteNoise(x-1,y-1) + WhiteNoise(x+1,y-1) + WhiteNoise(x-1,y+1) + WhiteNoise(x+1,y+1))/16;
    sides = (WhiteNoise(x-1,y) + WhiteNoise(x+1,y) + WhiteNoise(x,y+1) + WhiteNoise(x,y-1))/8;
    center = WhiteNoise(x,y)/4;
    re = (corners + sides + center);
    
end

function rb = interpolatedNoise2d(x,y, IMG_W, IMG_H, WhiteNoise)
    integerX = floor(x);
    integerY = floor(y);

    fractionalX = x-integerX;
    fractionalY = y-integerY;

    v1 = smoothedNoise2d(integerX,integerY, IMG_W, IMG_H, WhiteNoise);
    v2 = smoothedNoise2d(integerX+1,integerY, IMG_W, IMG_H, WhiteNoise);
    v3 = smoothedNoise2d(integerX, integerY+1, IMG_W, IMG_H, WhiteNoise);
    v4 = smoothedNoise2d(integerX+1,integerY+1, IMG_W, IMG_H, WhiteNoise);

    i1 = cosineInterpolate(v1,v2,fractionalX);
    i2 = cosineInterpolate(v3,v4, fractionalX);

    rb = cosineInterpolate(i1,i2,fractionalY);
end

function res = perlinNoise2d(x,y, gamma, Nf, IMG_W, IMG_H, WhiteNoise, g)
    P_total = 0;
    
    for i = 1:Nf
        frequency = 2^(i-1); 
        amplitude = gamma^(i-1);
        offsetV = -.5 + (.5+.5).*rand(1,1);
%         offsetVector(i) = offsetV;
%         offsetV = 0;
        P_total = P_total + interpolatedNoise2d(x*frequency + offsetV,y*frequency + offsetV, IMG_W, IMG_H, WhiteNoise)*amplitude;
%         P_total = P_total + Perlin_Noise_Jakes(x*frequency + offsetV, y*frequency + offsetV, g, IMG_W, IMG_H)*amplitude;
        
    end

    res = 1/2 + ((sqrt(2)*(1 - gamma))/(2*(1 - gamma^Nf))) * P_total;
end






%3D FUNCTIONS

function re = smoothedNoise3d(x,y,z, IMG_W, IMG_H, IMG_B, WhiteNoise)
    
    [widthloc, heightloc, baseloc] = size(WhiteNoise);

    centerX = floor(IMG_W / 2);
    centerY = floor(IMG_H / 2);
    centerZ = floor(IMG_B / 2);
    
    x = centerX + x;
    y = centerY + y;
    z = centerZ + z;
    
    if x < 0
        x = -x;
%     elseif x >= IMG_W
%         x = IMG_W-1;
    end
    if y < 0
        y = -y;
%     elseif y >= IMG_H
%         y = IMG_H-1;
    end
    if z < 0
        z = -z;
%     elseif y >= IMG_H
%         y = IMG_H-1;
    end
    
    
    if x == 1 || x == 0
        x = 2;
    end
    if y == 1 || y == 0
        y = 2;
    end
    if z == 1 || z == 0
        z = 2;
    end



    if x + 1 > widthloc
        x = widthloc-1;
    end
    if y + 1 > heightloc
        y = heightloc-1;
    end
    if z + 1 > baseloc
        z = baseloc-1;
    end


%     x
%     y
%     z


    
    corners = ( WhiteNoise(x-1,y-1,z-1) + WhiteNoise(x+1,y-1,z-1) + WhiteNoise(x-1,y+1,z-1) + WhiteNoise(x+1,y+1,z-1) + WhiteNoise(x-1,y-1,z+1) + WhiteNoise(x+1,y-1,z+1) + WhiteNoise(x-1,y+1,z+1) + WhiteNoise(x+1,y+1,z+1) )/32;
    sides = ( WhiteNoise(x-1,y,z) + WhiteNoise(x+1,y,z) + WhiteNoise(x,y+1,z) + WhiteNoise(x,y-1,z) + WhiteNoise(x,y,z+1) + WhiteNoise(x,y,z-1) )/12;
    center = WhiteNoise(x,y,z)/4;
    re = (corners + sides + center);
    
end

function rb = interpolatedNoise3d(x,y,z, IMG_W, IMG_H, IMG_B, WhiteNoise)
    integerX = floor(x);
    integerY = floor(y);
    integerZ = floor(z);

    fractionalX = x-integerX;
    fractionalY = y-integerY;
    fractionalZ = z-integerZ;

    %in lower face z
    v1 = smoothedNoise3d(integerX, integerY, integerZ, IMG_W, IMG_H, IMG_B, WhiteNoise);
    v2 = smoothedNoise3d(integerX+1, integerY, integerZ, IMG_W, IMG_H, IMG_B, WhiteNoise);
    v3 = smoothedNoise3d(integerX, integerY+1, integerZ, IMG_W, IMG_H, IMG_B, WhiteNoise);
    v4 = smoothedNoise3d(integerX+1, integerY+1, integerZ, IMG_W, IMG_H, IMG_B, WhiteNoise);

    %in upper face z
    v5 = smoothedNoise3d(integerX, integerY, integerZ+1, IMG_W, IMG_H, IMG_B, WhiteNoise);
    v6 = smoothedNoise3d(integerX+1, integerY, integerZ+1, IMG_W, IMG_H, IMG_B, WhiteNoise);
    v7 = smoothedNoise3d(integerX, integerY+1, integerZ+1, IMG_W, IMG_H, IMG_B, WhiteNoise);
    v8 = smoothedNoise3d(integerX+1, integerY+1, integerZ+1, IMG_W, IMG_H, IMG_B, WhiteNoise);

    %in lower face z
    i1 = cosineInterpolate(v1,v2,fractionalX);
    i2 = cosineInterpolate(v3,v4, fractionalX);
    rb1 = cosineInterpolate(i1,i2,fractionalY);
    
    %in upper face z
    i3 = cosineInterpolate(v5,v6,fractionalX);
    i4 = cosineInterpolate(v7,v8, fractionalX);
    rb2 = cosineInterpolate(i3,i4,fractionalY);
    
    %Final interpolation
    rb = cosineInterpolate(rb1,rb2,fractionalZ);
    
end

function res = perlinNoise3d(x,y,z, gamma, Nf, IMG_W, IMG_H, IMG_B, WhiteNoise)
    P_total = 0;
    
    for i = 1:Nf
        frequency = 2^(i-1); 
        amplitude = gamma^(i-1);
        offsetV = -.5 + (.5+.5).*rand(1,1);
%         offsetVector(i) = offsetV;
%         offsetV = 0;
        P_total = P_total + interpolatedNoise3d(x*frequency + offsetV,y*frequency + offsetV, z*frequency + offsetV, IMG_W, IMG_H, IMG_B, WhiteNoise)*amplitude;
%         P_total = P_total + Perlin_Noise_Jakes(x*frequency + offsetV, y*frequency + offsetV, g, IMG_W, IMG_H)*amplitude;
        
    end

    res = 1/2 + ((sqrt(2)*(1 - gamma))/(2*(1 - gamma^Nf))) * P_total;
end

end