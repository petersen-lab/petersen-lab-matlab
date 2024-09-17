%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:  Function that takes as input a neuron's spike times organized by trials
%               and computes its oscillation score, the confidence of the score's estimation 
%               and the peak oscillation frequency in a given band.
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
% Input:        SpikeTimesPerTrials - cell array of size (1 x Trial_Count) where each cell 
%                   contains an array of spike times corresponding to one trial. Spike times 
%                   are aligned to the start of the trial.
%               TrialLength - duration of trial in sample units
%               LowBoundFrequency - low boundary of the frequency band of interest in Hz
%               HighBoundFrequency - high boundary of the frequency band of interest in Hz
%               SamplingFrequency - sampling frequency of the signal in Hz
%
% Output:       OscScore - oscillation score for the specified frequency band. See paper
%                   for ranges of the score specific to different frequency bands.
%               CnfScore - confidence score in estimating the oscillation score, range [0..1]
%               OscFreq - peak oscillating frequency in the specified frequency band in Hz
%               AutoCorrelogram - array containing the autocorrelogram computed on all trials
%               AutoCorrelogramWithoutPeak - array containing the autocorrelogram computed
%                   on all trials, smoothed and with no central peak
%               Spectrum - array containing the frequency spectrum of the smoothed peakless
%                   autocorrelogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [OscScore, CnfScore, OscFreq, AutoCorrelogram, AutoCorrelogramWithoutPeak, Spectrum] = ...
    OScoreSpikes(SpikeTimesPerTrials, TrialLength, LowBoundFrequency, HighBoundFrequency, SamplingFrequency)

    TrialCount = size(SpikeTimesPerTrials,2);
    if(TrialCount<1 || TrialLength<=0)
        OscScore = -1; 
        return;
    end
    
    CorrelationWindow = floor(2^ceil(max(log2(3*SamplingFrequency/LowBoundFrequency),log2(SamplingFrequency/4.0))));
    AutoCorrelogramSize = 2 * CorrelationWindow + 1;
    
    % compute autocorrelogram and oscillation score for each individual trial
    TrialAutoCorrelograms = zeros(TrialCount,AutoCorrelogramSize);
    TrialOScores = zeros(TrialCount,1);
    for i=1:TrialCount
        Trial = zeros(1,TrialLength);
        Trial(SpikeTimesPerTrials{i}) = 1;
        TrialAutoCorrelograms(i,:) = xcorr(Trial,Trial,CorrelationWindow);
        [TrialOScores(i),Aux] = OScoreAC(TrialAutoCorrelograms(i,:), LowBoundFrequency, HighBoundFrequency, SamplingFrequency);
    end
    
    % compute oscillation score on the autocorrelogram for all trials
    AutoCorrelogram = sum(TrialAutoCorrelograms,1);
    [OscScore, OscFreq, AutoCorrelogramWithoutPeak, Spectrum] = OScoreAC(AutoCorrelogram, LowBoundFrequency, HighBoundFrequency, SamplingFrequency);

    % compute confidence score
    CV = std(TrialOScores(TrialOScores>=0))/mean(TrialOScores(TrialOScores>=0));
    if(TrialCount>1)
        CnfScore = 1.0 /(1.0 + CV);
    else
        CnfScore = 0;
    end;
    
end