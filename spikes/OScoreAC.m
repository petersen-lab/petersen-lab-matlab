%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:  Function that takes as input an auto-correlation histogram of spikes and 
%               computes the oscillation score and the peak oscillation frequency in a given band.
%
% Author:       Ovidiu F. Jurjut, 13.10.2011
%
% Disclaimer:   This code is freely usable for non-profit scientific purposes.
%               I do not warrant that the code is bug free. Use it at your own risk!
%
% Article:      The Oscillation Score: An Efficient Method for Estimating Oscillation 
%               Strength in Neuronal Activity
%               Muresan et al. 2008, Journal of Neurophysiology 99: 1333-1353               
%               http://jn.physiology.org/content/99/3/1333
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:        AutoCorrelogram - array containing the auto-correlation histogram
%               LowBoundFrequency - low boundary of the frequency band of interest in Hz
%               HighBoundFrequency - high boundary of the frequency band of interest in Hz
%               SamplingFrequency - sampling frequency of the signal in Hz
%
% Output:       OscScore - oscillation score for the specified frequency band. See paper
%                   for ranges of the score specific to different frequency bands.
%               OscFreq - peak oscillating frequency in the specified frequency band in Hz
%               AutoCorrelogramWithoutPeak - array containing the smoothed autocorrelogram
%                   with no central peak
%               Spectrum - array containing the frequency spectrum of the smoothed peakless
%                   autocorrelogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [OscScore, OscFreq, AutoCorrelogramWithoutPeak, Spectrum] = ...
    OScoreAC(AutoCorrelogram, LowBoundFrequency, HighBoundFrequency, SamplingFrequency)

    CorrelogramSize = size(AutoCorrelogram,2);
    HalfCorrelogramSize = floor(CorrelogramSize/2)+1;
    
    % Checking size of input autocorrelogram
    if(CorrelogramSize <= 0)
        OscScore = -1; 
        return;
    end
    
    % Checking frequency band input parameters
    if((LowBoundFrequency <= 0)||(HighBoundFrequency > SamplingFrequency/2)||(LowBoundFrequency > HighBoundFrequency))
        OscScore = -1; 
        return;
    end
    
    MirroredAutoCorrelogram = flipdim(AutoCorrelogram,2);
    
    % Computing "high" gaussian kernel with sigma as a function of the HighBoundFrequency
    KernelSTD = min([2 134.0/(1.5*HighBoundFrequency)]) * (SamplingFrequency/1000.0);
    GaussianKernel = -round(3*KernelSTD):round(3*KernelSTD);
    GaussianKernel = (1.0/(sqrt(2.0*pi) * KernelSTD) * exp(-(GaussianKernel.^2)/(2.0*KernelSTD*KernelSTD)));
    HalfKernelSize = floor(size(GaussianKernel,2)/2);

    % Smoothing the autocorrelogram with the "high" gaussian kernel (mirror signal at edges with HalfKernelSize+1 elements to reduce border effects)
    PaddedAutoCorrelogram =  [MirroredAutoCorrelogram(end-HalfKernelSize:end), AutoCorrelogram, MirroredAutoCorrelogram(1:HalfKernelSize+1)];
    HighSmoothedAutoCorrelogram = conv(PaddedAutoCorrelogram, GaussianKernel);
    HighSmoothedAutoCorrelogram = HighSmoothedAutoCorrelogram((size(GaussianKernel,2)+1):(end-size(GaussianKernel,2)));
    
    % Computing "low" gaussian kernel with sigma as a function of the LowBoundFrequency
    KernelSTD = 2 * 134.0/(1.5*LowBoundFrequency) * (SamplingFrequency/1000.0);
    GaussianKernel = -round(3*KernelSTD):round(3*KernelSTD);
    GaussianKernel = (1.0/(sqrt(2.0*pi) * KernelSTD) * exp(-(GaussianKernel.^2)/(2.0*KernelSTD*KernelSTD)));
    HalfKernelSize = floor(size(GaussianKernel,2)/2);

    % Smoothing the autocorrelogram with the "low" gaussian kernel (mirror signal at edges with HalfKernelSize+1 elements to reduce border effects)
    PaddedAutoCorrelogram =  [MirroredAutoCorrelogram(end-HalfKernelSize:end), AutoCorrelogram, MirroredAutoCorrelogram(1:HalfKernelSize+1)];
    LowSmoothedAutoCorrelogram = conv(PaddedAutoCorrelogram, GaussianKernel);
    LowSmoothedAutoCorrelogram = LowSmoothedAutoCorrelogram((size(GaussianKernel,2)+1):(end-size(GaussianKernel,2)));

    % Cutting central peak based on the low-smoothed autocorrelogram
        % Compute the scaling factor of the Y axis in the autocorrelogram such as to match the units of the time axis
        Min = min(LowSmoothedAutoCorrelogram);
        Min = min([Min 0]); % in case the shift predictor is not extracted minimum is 0
        Max = max(LowSmoothedAutoCorrelogram);
        if (Max - Min ~= 0)
            Ratio = CorrelogramSize/(Max - Min);
        else
            if (Max == 0)
                OscScore = -1; 
                return;
            else
                Ratio = 1;
            end
        end
    
        % Finding cutting points
        LimitSlopeLeft = tan(10.0 * pi / 180.0);
        Derivative = (LowSmoothedAutoCorrelogram(2:end)-LowSmoothedAutoCorrelogram(1:end-1))*Ratio;
        LeftCutIndex = find(Derivative(1:HalfCorrelogramSize-1) < LimitSlopeLeft,1,'last');
        RightCutIndex = HalfCorrelogramSize + (HalfCorrelogramSize-LeftCutIndex);
        
        % Cutting the peak and filling in with the mean values of the left and right cut points
        AutoCorrelogramWithoutPeak = HighSmoothedAutoCorrelogram;
        AutoCorrelogramWithoutPeak(LeftCutIndex:RightCutIndex) = mean([HighSmoothedAutoCorrelogram(LeftCutIndex),HighSmoothedAutoCorrelogram(RightCutIndex)]);
    
    % Applying Fast Fourier Transform on smoothed correlogram without peak (manual multiplication with Blackman window)
    BlackManWin = blackman(size(AutoCorrelogramWithoutPeak,2))';
    Spectrum = abs(fft(AutoCorrelogramWithoutPeak.*BlackManWin, 2^nextpow2(HalfCorrelogramSize)));
    SpectrumSize = size(Spectrum,2)/2;
    Spectrum = Spectrum(1:SpectrumSize);
    SpectrumIntegral = mean(Spectrum);
    
    % Search in the frequency spectrum for the peak frequency in the band
    LowFrequencyIndex = round(LowBoundFrequency * SpectrumSize/(SamplingFrequency/2));
    HighFrequencyIndex = round(HighBoundFrequency * SpectrumSize/(SamplingFrequency/2));
    [PeakFreqVal, PeakFreqIndex] = max(Spectrum(LowFrequencyIndex:HighFrequencyIndex));
    
    % Calculate oscillation score and peak frequency in the band
    if(SpectrumIntegral>0)
        OscScore = PeakFreqVal / SpectrumIntegral;
        OscFreq = ((PeakFreqIndex-1+LowFrequencyIndex-1) * (SamplingFrequency/2))/SpectrumSize;
    else
        OscScore = -1; 
        return;
    end
end