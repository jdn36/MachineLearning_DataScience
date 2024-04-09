function [dominantFreq, freqRang, FFTvalAbs, S1] = FreqCalcFourier(y, FsOrig, paddingBool, minimumdF)

    LOrig = length(y);  %number of data points
%     FsOrig = FsOrig.*1000;
%     TOrig = 1/(FsOrig);             % Sampling period      s     
                % Length of signal    data points we have (technically ms)
%     tOrig = (0:LOrig-1)*TOrig;        % Time vector    s

    % FsInt = 1000;            % Sampling frequency     Hz               
    % TInt = 1/FsInt;             % Sampling period      s     
    % LInt = 2000;             % Length of signal    msx10^-1
    % tInt = (0:LInt-1)*TInt;        % Time vector    s

    S1 = y;
    % S1Int = interp1(tOrig, S1 , tInt ,'spline');

    if paddingBool == 1
        while FsOrig/LOrig > minimumdF
            S1 = [S1; zeros(200,1); flip(S1)];
%             S1 = [S1; flip(S1); zeros(300,1)];
            LOrig = length(S1);
        end
    else

    end

    [FFTval, freqRang] = posFFT(S1, FsOrig);

%     freqRang = freqRang(2:end-1);
%     FFTval = FFTval(2:end-1);

    % subplot(2,1,1)
    % plot(S1)
    % subplot(2,1,2)
    % stem(freqRang, abs(FFTval))
    % title("Single-Sided Amplitude Spectrum of X(t)")
    % xlabel("f (Hz)")
    % ylabel("|Power(f)|")

    dominantFreq = freqRang(find(abs(FFTval) == max(abs(FFTval))));
    FFTvalAbs = FFTval;


end