%function to classify signals based on characteristics. Classifying into
%smooth biphasic, complex multiphasic, and signals that could not be
%classified.

function [patient_SmoothBiphasic, patient_NotClassified, patient_ComplexMultiphasic, signalCharacteristics_SmoothBiphasic, signalCharacteristics_NotClassified, signalCharacteristics_ComplexMultiphasic] =  SignalClassified(elec1, surfaceCoordinates, electrodeCoordinates, timePointsToTake, numDataPoints, dt, DeflectionCountAllowed, EGMdurationAllowed, signalWidthAllowed, thresholdAmplitude)

%INPUTS:
%elec1 - matrix with all the signals where rows are time points and columns observations
%surfaceCoordinates - matrix with surface points where rows are point indices and columns are x,y,z respectively
%electrodeCoordinates - matrix with electrode geometry points where rows are point indices and columns are x,y,z respectively
%timePointsToTake - ending of time points (total signal length)
%numDataPoints - the length of the cut signal that contains the EGM to minimize points describing the signal
%dt - dt or time step (optimized for 1kHz resolution so 1ms)
%DeflectionCountAllowed - threshold for smooth biphasic classification (less than 4 deflections is typically considered smooth biphasic)
%EGMdurationAllowed - threshold for EGM duration for smooth biphasic; (less than 80ms was considered smooth biphasic)
%signalWidthAllowed - threshold for signal width for smooth biphasic; (less than 25ms was considered smooth biphasic)
%thresholdAmplitude - threshold for amplitude for complex signals; (less than .5mV was considered not classified (typically that is the threshold for fractionated EGMs as well))


%OUTPUTS:
%patient_SmoothBiphasic - indices of the signals that are smooth biphasic
%patient_ComplexMultiphasic - indices of the signals that are complex multiphasic
%patient_NotClassified - indices of the signals that were not classified
%signalCharacteristics_SmoothBiphasic - matrix of characteristics of smooth biphasic signals
%signalCharacteristics_ComplexMultiphasic - matrix of characteristics of complex multiphasic signals
%signalCharacteristics_NotClassified - matrix of characteristics of signals not classified

%The characteristic matrix has rows corresponding to the indices found in
%their respective vector coming from EGM observation and the columns are:
%
%    1               2                                             3                                4           5            6               7                        8                  9
%[downstroke, rising upstroke (positive deflection), recovery upstroke (negative deflection), signal width, baseline, EGM duration, distance to electrode, peak-to-peak amplitude, deflection count ];



    EGMDuration_SmoothBiphasic = [];
    deflectionCount_SmoothBiphasic = [];
    signals_SmoothBiphasic = [];
    signals_SmoothBiphasicOG = [];
    downstroke_SmoothBiphasic = [];
    upstrokeRise_SmoothBiphasic = [];
    upstrokeBack_SmoothBiphasic = [];
    signalWidth_SmoothBiphasic = [];
    baseline_SmoothBiphasic = [];
    PPA_SmoothBiphasic = [];
    AlldistanceToProj_SmoothBiphasic = [];
    patient_SmoothBiphasic = [];

    EGMDuration_NotClassified = [];
    deflectionCount_NotClassified = [];
    signals_NotClassified = [];
    signals_NotClassifiedOG = [];
    downstroke_NotClassified = [];
    baseline_NotClassified = [];
    PPA_NotClassified = [];
    AlldistanceToProj_NotClassified = [];
    patient_NotClassified = [];

    EGMDuration_ComplexMultiphasic = [];
    deflectionCount_ComplexMultiphasic = [];
    signals_ComplexMultiphasic = [];
    signals_ComplexMultiphasicOG = [];
    downstroke_ComplexMultiphasic = [];
    upstrokeRise_ComplexMultiphasic = [];
    upstrokeBack_ComplexMultiphasic = [];
    signalWidth_ComplexMultiphasic = [];
    baseline_ComplexMultiphasic = [];
    PPA_ComplexMultiphasic = [];
    AlldistanceToProj_ComplexMultiphasic = [];
    patient_ComplexMultiphasic = [];



    ptCloud = pointCloud([surfaceCoordinates(:,1), surfaceCoordinates(:,2), surfaceCoordinates(:,3)]);%create point cloud
    tvectShiftNorm = 1:dt:timePointsToTake;

    for i=1:length(elec1(1,:))
        [EGMduration,DeflectionCount, xtoCheck, ytoCheck] = SD_SlidingWindowFunction([dt:dt:length(elec1(:,i))], elec1(:,i), dt, 40, 15, 3, numDataPoints);
        locSignal =  ( elec1(:,i) - min(elec1(:,i)) ) ./ ( max(elec1(:,i)) - min(elec1(:,i)) );
        [dVeMono_dt, downstrokeMono, midSignalIndMono, signalDiffMonoPointsTime, upstrokeRiseMono, upstrokeBackMono, riseSignalIndMono, backSignalIndMono, PositiveDeflectionDuration_Mono, NegativeDeflectionDuration_Mono, baselineMono, signalWidth, timePeak, timeTrough] = EGM_Characterizing_Function(locSignal, tvectShiftNorm, dt, 1/dt, EGMduration);
        [indices,dists] = findNearestNeighbors(ptCloud,[electrodeCoordinates(i,1), electrodeCoordinates(i,2), electrodeCoordinates(i,3)],1);
        
        if (DeflectionCount < DeflectionCountAllowed  && EGMduration < EGMdurationAllowed) && signalWidth <= signalWidthAllowed

            locSignalCut =  ( ytoCheck - min(ytoCheck) ) ./ ( max(ytoCheck) - min(ytoCheck) );
            locSignalCutOG =  ytoCheck;
            localPPA = max(ytoCheck) - min(ytoCheck);

            EGMDuration_SmoothBiphasic = [EGMDuration_SmoothBiphasic; EGMduration];
            deflectionCount_SmoothBiphasic = [deflectionCount_SmoothBiphasic; DeflectionCount];
            signals_SmoothBiphasic = [ signals_SmoothBiphasic; locSignalCut' ];
            signals_SmoothBiphasicOG = [ signals_SmoothBiphasicOG; locSignalCutOG' ];
            downstroke_SmoothBiphasic = [ downstroke_SmoothBiphasic; downstrokeMono ];
            upstrokeRise_SmoothBiphasic = [ upstrokeRise_SmoothBiphasic; upstrokeRiseMono ];
            upstrokeBack_SmoothBiphasic = [ upstrokeBack_SmoothBiphasic; upstrokeBackMono ];
            signalWidth_SmoothBiphasic = [ signalWidth_SmoothBiphasic; signalWidth ];
            baseline_SmoothBiphasic = [ baseline_SmoothBiphasic; baselineMono ];
            PPA_SmoothBiphasic = [ PPA_SmoothBiphasic; localPPA ];
            AlldistanceToProj_SmoothBiphasic = [ AlldistanceToProj_SmoothBiphasic; dists ];
            patient_SmoothBiphasic = [ patient_SmoothBiphasic; i ];
            
        else

            [EGMduration, DeflectionCount, xtoCheck, ytoCheck] = SD_SlidingWindowFunction([dt:dt:length(elec1(:,i))], elec1(:,i), dt, 40, 15, 3, length(elec1(:,1)));
            locSignal =  ( elec1(:,i) - min(elec1(:,i)) ) ./ ( max(elec1(:,i)) - min(elec1(:,i)) );
            locSignalCut =  ( ytoCheck - min(ytoCheck) ) ./ ( max(ytoCheck) - min(ytoCheck) );
            locSignalCutOG =  ytoCheck;
            [dVeMono_dt, downstrokeMono, midSignalIndMono, signalDiffMonoPointsTime, upstrokeRiseMono, upstrokeBackMono, riseSignalIndMono, backSignalIndMono, PositiveDeflectionDuration_Mono, NegativeDeflectionDuration_Mono, baselineMono, signalWidth, timePeak, timeTrough] = EGM_Characterizing_Function(locSignal, tvectShiftNorm, dt, 1/dt, EGMduration);
            localPPA = max(ytoCheck) - min(ytoCheck);

            %CHECK FOR FIRST AND LAST VALUE OF THE NORMALIZED SIGNALS TO
            %SEE HOW FAR APART THE AMPLITUDE IS (IF THE SIGNAL IS KIND OF
            %NORMAL AND THERE ARE BUMPS THE BASELINE SHOULD NOT CHANGE
            %AFTER THE BIPHASIC PART). MAYBE TRY THRESHOLD OF .55
            %DIFFERENCE. ALSO CHECK THE PEAK-PEAK AMPLITUDE, SET A
            %THRESHOLD FOR ANYTHING BELOW THAT TO BE CONSIDERED NOT AN EGM
            %BUT MAYBE JUST NOISE. LAST THING IS, WHERE IS THE TROUGH OR
            %PEAK IN COMPARISON TO UPSTROKE RISE AND UPSTROKE BACK: IF THE
            %TIME DIFFERENCE BETWEEN WHERE THESE POINTS ARE OCURRING IS TOO
            %LARGE OR THE UPSTROKE BACK FOR INSTANCE HAPPENS BEFORE THE
            %TROUGH, THEN THE SIGNAL IS JUST NOT A BIPHASIC SIGNAL AT ALL.
            %IF ALL OF THESE THINGS ARE TRUE, THEN THE SIGNAL IS NEITHER
            %BIPHASIC NOR MULTIPHASIC. SAVE IN A DIFFERENT BIN.

            

            if abs(locSignalCut(1) - locSignalCut(end)) > thresholdAmplitude & ( (abs(ytoCheck(timePeak(1)) - ytoCheck(timeTrough(1))) ~= localPPA) | localPPA < thresholdAmplitude | signalDiffMonoPointsTime(backSignalIndMono(1)) < timeTrough(1))
                %CFAE or things that cannot be identified
                EGMDuration_NotClassified = [EGMDuration_NotClassified; EGMduration];
                deflectionCount_NotClassified = [deflectionCount_NotClassified; DeflectionCount];
                signals_NotClassified = [ signals_NotClassified; locSignalCut' ];
                signals_NotClassifiedOG = [ signals_NotClassifiedOG; locSignalCutOG' ];
                downstroke_NotClassified = [ downstroke_NotClassified; downstrokeMono ];
                baseline_NotClassified = [ baseline_NotClassified; baselineMono ];
                PPA_NotClassified = [ PPA_NotClassified; localPPA ];
                AlldistanceToProj_NotClassified = [ AlldistanceToProj_NotClassified; dists ];
                patient_NotClassified = [ patient_NotClassified; i ];
            else
                %multiphasic complex signals
                EGMDuration_ComplexMultiphasic = [EGMDuration_ComplexMultiphasic; EGMduration];
                deflectionCount_ComplexMultiphasic = [deflectionCount_ComplexMultiphasic; DeflectionCount];
                signals_ComplexMultiphasic = [ signals_ComplexMultiphasic; locSignalCut' ];
                signals_ComplexMultiphasicOG = [ signals_ComplexMultiphasicOG; locSignalCutOG' ];
                downstroke_ComplexMultiphasic = [ downstroke_ComplexMultiphasic; downstrokeMono ];
                upstrokeRise_ComplexMultiphasic = [ upstrokeRise_ComplexMultiphasic; upstrokeRiseMono ];
                upstrokeBack_ComplexMultiphasic = [ upstrokeBack_ComplexMultiphasic; upstrokeBackMono ];
                signalWidth_ComplexMultiphasic = [ signalWidth_ComplexMultiphasic; signalWidth ];
                baseline_ComplexMultiphasic = [ baseline_ComplexMultiphasic; baselineMono ];
                PPA_ComplexMultiphasic = [ PPA_ComplexMultiphasic; localPPA ];
                AlldistanceToProj_ComplexMultiphasic = [ AlldistanceToProj_ComplexMultiphasic; dists ];
                patient_ComplexMultiphasic = [ patient_ComplexMultiphasic; i ];
            end
            
        end
    end


    signalCharacteristics_SmoothBiphasic = [downstroke_SmoothBiphasic, upstrokeRise_SmoothBiphasic, upstrokeBack_SmoothBiphasic, signalWidth_SmoothBiphasic, baseline_SmoothBiphasic, EGMDuration_SmoothBiphasic, AlldistanceToProj_SmoothBiphasic, PPA_SmoothBiphasic, deflectionCount_SmoothBiphasic ];
    signalCharacteristics_NotClassified = [downstroke_NotClassified, baseline_NotClassified, EGMDuration_NotClassified, AlldistanceToProj_NotClassified, PPA_NotClassified, deflectionCount_NotClassified ];
    signalCharacteristics_ComplexMultiphasic = [downstroke_ComplexMultiphasic, upstrokeRise_ComplexMultiphasic, upstrokeBack_ComplexMultiphasic, signalWidth_ComplexMultiphasic, baseline_ComplexMultiphasic, EGMDuration_ComplexMultiphasic, AlldistanceToProj_ComplexMultiphasic, PPA_ComplexMultiphasic, deflectionCount_ComplexMultiphasic ];
   

end





