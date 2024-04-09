%Function to get EGM duration and deflection count based on sliding standard deviation window. Cuts the signal 

function [EGMduration, DeflectionCount, xtoCheckSend, ytoCheckSend] = SD_SlidingWindowFunction(x, y, dt, windowWidth, thresholdPercent, factorWindow, numDataPoints)

%EGM duration based on: A novel algorithm for 3-D visualization of electrogram duration for substrate-mapping in patients with ischemic heart disease and ventricular tachycardia
%masjedi et. al. 2021


%INPUT:
%x - time vector
%y - signal vector
%dt - time step (msec)
%windowWidth - parameter corresponding to the sliding window width to use based on Masjedi's paper (40ms used)
%thresholdPercent - percentage used based on Masjedi's papaer to detect duration based on sliding standard deviation window (15% used)
%factorWindow - factor to make window smaller to detect derivatives to calculate number of deflections (3 used)
%numDataPoints - number of data points to export



%OUTPUT:
%EGMduration - EGM duration in msec
%DeflectionCount - deflectoin count
%xtoCheckSend - time vector of exported points
%ytoCheckSend - signal vector of exported points


windowWidth = (windowWidth/dt);

locVec = [];
SD_curve = [];

for i = 1:length(x)
    
    if i-(windowWidth/2) < 1
        locVec = y( 1 : i+(windowWidth/2) );
    elseif i+ceil(windowWidth/2) > length(x)
        locVec = y( i-(windowWidth/2) : length(x) );
    else
        locVec = y( i-(windowWidth/2) : i+(windowWidth/2) );
    end
    
    locSTD = std(locVec);
    SD_curve(i) = locSTD;
    
end

threshValue = (thresholdPercent/100)*max(SD_curve);
aboveValues = find(SD_curve >= threshValue);
EGMduration = x(aboveValues(end)) - x(aboveValues(1));

EGMduration

belowValues = find(SD_curve < threshValue);
ynorm = ( y - min(y) )./( max(y) - min(y) );
baselineMeanSDSW = mean( ynorm(belowValues) );



xtoCheck = x(aboveValues);
ytoCheck = y(aboveValues);
SDtoCheck = SD_curve(aboveValues);

avgSlope = [];
avgSlope2 = [];

for i = 1:length(xtoCheck)
    if i-ceil(windowWidth/factorWindow) < 1 && i+ceil(windowWidth/factorWindow) > length(xtoCheck)
        locVec = ytoCheck( 1 : length(xtoCheck) );
    elseif i-ceil(windowWidth/factorWindow) < 1
        locVec = ytoCheck( 1 : i+ceil(windowWidth/factorWindow) );
    elseif i+ceil(windowWidth/factorWindow) > length(xtoCheck)
        locVec = ytoCheck( i-ceil(windowWidth/factorWindow) : length(xtoCheck) );
    else
        locVec = ytoCheck( i-ceil(windowWidth/factorWindow) : i+ceil(windowWidth/factorWindow) );
    end
    avgSlope(i) = mean(sign(diff(locVec))); 
    avgSlope2(i) = mean((diff(locVec)));
end

derSlope = diff(avgSlope);
totCount = 0;


derSlope = (derSlope - min(derSlope))./(max(derSlope) - min(derSlope));
maxValSlope = max(derSlope);
minValSlope = min(derSlope);
meanValSlope = mean(derSlope);
% for i=2:length(derSlope)
%     if abs( derSlope(i) - derSlope(i-1) ) > (maxValSlope - minValSlope)*.4
%         totCount = totCount + 1;
%     end
% end
% DeflectionCount = ceil(totCount/2)
for i=2:length(derSlope)
    % if sign( (derSlope(i)-meanValSlope) - (derSlope(i-1)-meanValSlope) ) ~= sign( (derSlope(i)-meanValSlope) - (derSlope(i+1)-meanValSlope) )
    if sign( ( round( derSlope(i)-meanValSlope, 1 ) ) ) ~= sign( ( round( derSlope(i-1)-meanValSlope, 1 ) ) )
        % if abs((derSlope(i)-meanValSlope) - (derSlope(i-1)-meanValSlope)) > 0.005 && abs((derSlope(i)-meanValSlope) - (derSlope(i+1)-meanValSlope)) > 0.005 
            totCount = totCount + 1;
        % end
    end
end



normavgSlope = (avgSlope2 - min(avgSlope2))./(max(avgSlope2) - min(avgSlope2));
meannormavgSlope = mean(normavgSlope);
shiftmeanNAS = normavgSlope - meannormavgSlope;
minNAS = min(normavgSlope - meannormavgSlope);
maxNAS = max(normavgSlope - meannormavgSlope);

trackingChange = 0;
currentSlope = 0;
previousSlope = 0;
newDeflectionCount = 0;
allowedPercentage = .1;

timePointsDeflections = [];
yPointsDeflections = [];


for i = 1:length(shiftmeanNAS)

    if i == 1
        currentSlope = sign(diff( shiftmeanNAS(i:i+1) ));
        previousSlope = sign(diff( shiftmeanNAS(i:i+1) ));

        if sign(shiftmeanNAS(i)) ~= sign(shiftmeanNAS(i+1)) 
            trackingChange = 1;
            % disp('crossed base')
        end
    else
        previousSlope = currentSlope;
        currentSlope = sign(diff( shiftmeanNAS(i-1:i) ));

        if sign(shiftmeanNAS(i)) ~= sign(shiftmeanNAS(i-1)) || abs(round(shiftmeanNAS(i),1)) < allowedPercentage
            trackingChange = 1;
            % disp('crossed base')
        end
    end

    if trackingChange == 1
        if previousSlope ~= currentSlope
            if shiftmeanNAS(i) >= allowedPercentage*maxNAS || shiftmeanNAS(i) <= allowedPercentage*minNAS
                newDeflectionCount = newDeflectionCount + 1;
                trackingChange = 0;
                timePointsDeflections = [timePointsDeflections, xtoCheck(i)];
                yPointsDeflections = [yPointsDeflections, ytoCheck(i)];
            end
        end
        if i == length(shiftmeanNAS) && trackingChange == 1
            if shiftmeanNAS(i) >= allowedPercentage*maxNAS || shiftmeanNAS(i) <= allowedPercentage*minNAS
                newDeflectionCount = newDeflectionCount + 1;
                trackingChange = 0;
                timePointsDeflections = [timePointsDeflections, xtoCheck(i)];
                yPointsDeflections = [yPointsDeflections, ytoCheck(i)];
            end
        end
    end
    
    
end

DeflectionCount = newDeflectionCount;
DeflectionCount2 = floor(totCount/2)


extraPointsSide = floor( (numDataPoints - length(ytoCheck)) / 2 );
initPoint = 1;
finPoint = 1;

if aboveValues(1) - extraPointsSide <= 0
    initPoint = 1;
    finPoint = numDataPoints;
elseif aboveValues(end) + extraPointsSide > length(y)
    initPoint = length(y) - numDataPoints;
    finPoint = length(y);
else
    initPoint = aboveValues(1) - extraPointsSide;
    finPoint = aboveValues(end) + extraPointsSide;
end

% initPoint
% finPoint

if length(x) == numDataPoints
    xtoCheckSend = x;
    ytoCheckSend = y;
else
    xtoCheckSend = x(initPoint:finPoint);
    ytoCheckSend = y(initPoint:finPoint);
    
    if length(ytoCheckSend) < numDataPoints
        xtoCheckSend(end+1) = xtoCheckSend(end);
        ytoCheckSend(end+1) = ytoCheckSend(end);
    end
    
    if length(ytoCheckSend) < numDataPoints
        xtoCheckSend(end+1) = xtoCheckSend(end);
        ytoCheckSend(end+1) = ytoCheckSend(end);
    end



    if length(ytoCheckSend) > numDataPoints
        xtoCheckSend = xtoCheckSend(1:numDataPoints);
        ytoCheckSend = ytoCheckSend(1:numDataPoints);
    end


end
% 
% vq = interp1(1:1:length(ytocheck),ytocheck,length(ytocheck)/numDataPoints:length(ytocheck)/numDataPoints:length(ytocheck) )






%%FINAL METHOD TO CALCULATE DEFLECTION COUNT LEULLOCHE ET. AL. 2007 AND JACQUEMET 2009

trackingChange1 = 0;
newDeflectionCount = 0;
allowedPercentage = .1;

timePointsDeflections2 = [];
yPointsDeflections2 = [];

maxytoCheckSend = max(ytoCheck);
minytoCheckSend = min(ytoCheck);

peak1 = 0;
peak2 = 0;

currLocSlope = 0;
prevLocSlope = 0;

% size(ytoCheckSend)
% size(ytoCheck)

for i = 1:length(ytoCheck)-1
    if i == 1
        currLocSlope = sign(diff(ytoCheck(i:i+1))./dt);
        prevLocSlope = sign(diff(ytoCheck(i:i+1))./dt);
        peak1 = ytoCheck(i);
        peak2 = ytoCheck(i);
    else
        prevLocSlope = currLocSlope;
        currLocSlope = sign(diff(ytoCheck(i:i+1))./dt);
    end

    if abs(round(currLocSlope,1)) <= .1 || currLocSlope ~= prevLocSlope
        trackingChange1 = 1;
        peak2 = ytoCheck(i);
    end

    if trackingChange1 == 1
        if abs(peak1 - peak2) > abs((maxytoCheckSend - minytoCheckSend)*allowedPercentage) && (currLocSlope ~= prevLocSlope || abs(round(currLocSlope,1)) <= .1)
            newDeflectionCount = newDeflectionCount + 1;
            peak1 = ytoCheck(i);
            peak2 = ytoCheck(i);
            trackingChange1 = 0;

            timePointsDeflections2 = [timePointsDeflections2, i];
            yPointsDeflections2 = [yPointsDeflections2, ytoCheck(i)];
        end
    end
end

DeflectionCount = newDeflectionCount

bumpinessFactor = 0;

if DeflectionCount2 > 0
    bumpinessFactor = abs(DeflectionCount - DeflectionCount2)/DeflectionCount2
end
% 
% figure
% subplot(4,1,1)
% plot(x,y)
% hold on
% plot(x(aboveValues(1)), y(aboveValues(1)),'*r')
% plot(x(aboveValues(end)), y(aboveValues(end)),'*r')
% plot(timePointsDeflections, yPointsDeflections, 'gv')
% title1 = "EGM duration is "+ num2str(EGMduration) + " ms with " + DeflectionCount + " deflections and 'bumpiness' factor of " + bumpinessFactor;
% title(title1)
% subplot(4,1,2)
% plot(x,SD_curve, 'b')
% hold on
% plot(x(aboveValues(1)), SD_curve(aboveValues(1)),'*r')
% plot(x(aboveValues(end)), SD_curve(aboveValues(end)),'*r')
% subplot(4,1,3)
% plot(xtoCheck(2:end), derSlope - meanValSlope)
% hold on
% plot(xtoCheck(2:end), 0.*ones(1,length(xtoCheck(2:end))), 'r')
% axis([min(x) max(x) min(derSlope - meanValSlope) max(derSlope - meanValSlope)])
% subplot(4,1,4)
% plot(xtoCheck, shiftmeanNAS)
% hold on
% plot(xtoCheck, 0.*ones(1,length(xtoCheck)), 'r')
% % plot(timePointsDeflections, 0.*ones(1,length(timePointsDeflections)), 'gv')
% axis([min(x) max(x) minNAS maxNAS])
% 
% 
% 
% figure
% plot(ytoCheck)
% hold on
% plot(timePointsDeflections2, yPointsDeflections2, 'gv')
% title(["number of data points is ", num2str(length(xtoCheck))])
% 
% 

