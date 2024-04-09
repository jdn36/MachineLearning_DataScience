
%Function to characterize the signal and get characteristic parameters
function [dVe_dt, downstroke, midSignalInd, signalDiffPointsTime, upstrokeRise,...
    upstrokeBack, riseSignalInd, backSignalInd, PositiveDeflectionDuration,...
    NegativeDeflectionDuration, baseline, signalWidth, timePeak, timeTrough ] = EGM_Characterizing_Function(sampleShiftNorm, tvectShiftNorm, dt, dnsizfact, EGMduration)


    %sampleShiftNorm - original signal vector
    %tvectShiftNorm - time vector based on original signal
    %dt - dt or time step of the original signal
    %dnsizfact - the factor needed to multiply dt to get to 1kHz resolution so new time step of 1ms
    %EGMduration - EGM duration in ms

    %OUTPUT:
    %dVe_dt - EGM derivative
    %downstroke - downstroke of normalized signal
    %midSignalInd - index of data point where the downstroke happened
    %signalDiffPointsTime - time vector for the points containing EGM
    %upstrokeRise - rising upstroke of normalized signal (positive deflection)
    %upstrokeBack - recovery upstroke of normalized signal (negative deflection)
    %riseSignalInd - index of data point where the rising upstroke happened
    %backSignalInd - index of data point where the recovery upstroke happened
    %PositiveDeflectionDuration - the duration of the positive deflection (time length between rising upstroke and downstroke)
    %NegativeDeflectionDuration - the duration of the negative deflection (time length between recovery upstroke and downstroke)
    %baseline - baseline of normalized signal
    %signalWidth - time width between peak and trough of the signal
    %timePeak - time that the peak happened
    %timeTrough - time that the trough happened



    %characterize normalized egm's more
    %THIS SHOULD WORK WELL WITH HEALTHY EGM'S

    %first derivative
    dVe_dt = diff(sampleShiftNorm)./(dt.*dnsizfact);

    %downstroke calculation
    downstroke = min(dVe_dt);

    %finding the activation point of the EGM to split the signal into
    %"positive" and "negative" deflection and to also select the part of
    %the signal that contains the + and - deflection "only"
    midSignalInd = find(dVe_dt==downstroke);
    midSignalInd = midSignalInd(1);

    if midSignalInd-ceil(EGMduration/2) < 1
        firstValueM = 1;
    else
        firstValueM = midSignalInd-ceil(EGMduration/2);
    end
    if midSignalInd+ceil(EGMduration/2) > length(dVe_dt)
        lastValueM = length(dVe_dt);
    else
        lastValueM = midSignalInd+ceil(EGMduration/2);
    end


    %cropping the data
    if firstValueM(1) < 1
        firstValueM = firstValueM(2);
        if lastValueM(2) > length(dVe_dt)
            lastValueM = lastValueM(1);
        else
            lastValueM = lastValueM(2);
        end
    else
        firstValueM = firstValueM(1);
        lastValueM = lastValueM(1);
    end
    % if lastValueM(1) > length(dVe_dt)
    %     lastValueM = length(dVe_dt);
    % else
    %     lastValueM = lastValueM(1);
    % end
    
    

    % if firstValueM < 1
    %     firstValueM = 1;
    % end
    % if lastValueM > length(sampleShiftNorm)*(dt*dnsizfact)-1
    %     lastValueM = length(sampleShiftNorm)*(dt*dnsizfact)-1;
    % end

    % firstValueM
    % lastValueM

    signalDiffPoints = dVe_dt( firstValueM:lastValueM );
    signalDiffPointsTime = tvectShiftNorm( firstValueM:lastValueM );
    downstroke = min(signalDiffPoints);

    %finding the upstrokes and where they are happening
    upstrokeRise = max( signalDiffPoints(1:find(signalDiffPoints==downstroke))  );
    upstrokeBack = max( signalDiffPoints(find(signalDiffPoints==downstroke):end)  );
    riseSignalInd = find(signalDiffPoints == upstrokeRise);
    backSignalInd = find(signalDiffPoints == upstrokeBack);

    %duration of "individual" deflections:
    %Positive deflection duration is defined as the time period between the
    %point where the first upstroke happens and the activation time point
    %Negative deflection duration is defined as the time period between the
    %activation time point and the point where the second upstroke happens
    %baseline is a measure of symmetry in x axis. If perfectly symmetrical
    %(positive and negative deflections have the same relative amplitude),
    %the baseline value is .5; i.e. if the negative deflection is stronger
    %(larger relative amplitude) the baseline is shifted upwards > .5
    PositiveDeflectionDuration = tvectShiftNorm(midSignalInd(1)) - signalDiffPointsTime( riseSignalInd(1) );
    NegativeDeflectionDuration = signalDiffPointsTime( backSignalInd(1) ) - tvectShiftNorm(midSignalInd(1));
    baseline = mean( [ sampleShiftNorm( 1:firstValueM )' sampleShiftNorm(lastValueM:end)' ] );

    midSignalInd = midSignalInd(1);
    riseSignalInd = riseSignalInd(1);
    backSignalInd = backSignalInd(1);

    %calculating width between max and min values of EGM (peak - trough duration)
    peakPoint = max(sampleShiftNorm);
    troughPoint = min(sampleShiftNorm);
    timePeak = tvectShiftNorm(find( sampleShiftNorm == peakPoint ));
    timeTrough = tvectShiftNorm(find( sampleShiftNorm == troughPoint ));
    signalWidth = abs( timeTrough(1) - timePeak(1) );



    % 
    % 
    % figure
    % plot(tvectShiftNorm, sampleShiftNorm,'DisplayName','Signal','LineWidth',4)
    % hold on
    % plot(signalDiffPointsTime(riseSignalInd), sampleShiftNorm( find(tvectShiftNorm == signalDiffPointsTime(riseSignalInd)) ),'m*','MarkerSize',12,'DisplayName','Rising Upstroke')
    % plot(tvectShiftNorm(midSignalInd), sampleShiftNorm(midSignalInd),'mv','MarkerSize',12,'DisplayName','Activation (Downstroke)')
    % plot(signalDiffPointsTime(backSignalInd), sampleShiftNorm( find(tvectShiftNorm == signalDiffPointsTime(backSignalInd)) ),'mo','MarkerSize',12,'DisplayName','Back upstroke')
    % plot(timePeak, 1, 'm+','MarkerSize',12, 'DisplayName','Peak')
    % plot(timeTrough, 0, 'm+','MarkerSize',12, 'DisplayName','Trough')
    % % title([pointsDescription, 'Normalized and Shifted Signal'])
    % xlabel('time (ms)')
    % ylabel('Potential (mV)')
    % legend()
    % set(gca,'FontSize',17)
    % 



end


