function [generalisedPhase, params] = generalisedPhaseForPointProcess(times, options)
% [generalisedPhase, params] = generalisedPhaseForPointProcess(times, <options>)
%
% Calculates generalised wide-band phase and instantaneous frequency for a
% given point process signal using generalized_phase_vector function of
% generalized-phase toolbox (https://github.com/mullerlab/generalized-phase).
%
% Args:
%   times (numeric, required, positional): a shape-(M, N) numeric array
%     where either M = 1 or N = 1. This input vector contains times when
%     discrete events occur (i.e., spike times).
%   freqRange (numeric, optional, keyword): a shape-(1, 2) vector array
%     with frequency cutoff values for evaluating the wide-band phase of
%     the signal (default=[2 12]).
%   passbandFilter (logical, optional, keyword): a shape-(1, 1) logical
%     scalar controlling whether to bandpass filter the signal within the
%     freqRange (default=true).
%   stepsize (numeric, optional, keyword): a shape-(1, 1) scalar
%     representing the step size (in seconds) used for binning spikes to
%     produce a continuous spiking rate trace (default=0.002).
%   convolutionPoints (numeric, optional, keyword): a shape-(1, 1) scalar
%     representing points of Gaussian convolution (gausswin) used to smooth
%     binned spiking rate (default=25 sample points).
%   showPhase (logical, optional, keyword): a shape-(1, 1) logical scalar
%     for ploting the wide-band phase of the oscillation/fluctuation
%     (default=false). As part of displaying the pahse, the convolved
%     spiking rate and the filtered spiking rate are also displayed.
%   smoothFreq (logical, optional, keyword): a shape-(1, 1) logical scalar
%     for smoothing instantaneous frequency estimates (default=false).
%
% Returns:
%   generalisedPhase (struct): a structure with the following fields:
%     phase (numeric): a shape-(1, L) numeric array of wide-band
%       oscillation/fluctuation phases estimated using
%       generalized_phase_vector function. The number of samples is
%       determined by the stepsize input parameter and correspond to the
%       length of the convolved spiking rate vector.
%     phaseUnwrapped (numeric): a shape-(1, L) numeric array of
%       unwrapped wide-band oscillation/fluctuation phases estimated using
%       generalized_phase_vector function.
%     instantFreq (numeric): a shape-(1, L) numeric array of wide-band
%       oscillation/fluctuation instantenous frequencies estimated using
%       generalized_phase_vector function.
%     amplitude (numeric): a shape-(1, L) numeric array of amplitude
%       (envelope) of the filtered and Hilbert-transformed spiking signal.
%     spikingRate (numeric): a shape-(1, L) numeric array of convolved
%       population spiking rate (in spikes/s; continuous).
%     filteredSpikes (numeric): a shape-(1, L) numeric array of convolved
%       and filtered population spiking rate (continuous).
%     timestamps (numeric): a shape-(1, L) numeric array of timestamps
%       corresponding to wide-band oscillation/fluctuation phase estimate
%       array.
%     instantFreqPerCycle (numeric): a shape-(1, K) numeric array of
%       wide-band instantenous frequencies estimated per individual
%       oscillation/fluctuation cycle based on the phase info. The number
%       of samples corresponds to the number of detected cycles.
%     cycleEndpoints (numeric): a shape-(1, K) numeric array of wide-band
%       oscillation/fluctuation cycle enpoint times.
%   params (struct): a structure with its fields corresponding to the
%     function input parameters including the default values.
%
% Dependencies:
%   CellExplorer (https://cellexplorer.org/).
%   petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).
%   mullerlab/generalized-phase
%     (https://github.com/mullerlab/generalized-phase).
%   mullerlab/wave-matlab (https://github.com/mullerlab/wave-matlab).
%   Circular Statistics Toolbox
%     (https://github.com/circstat/circstat-matlab).
%   smoothn function on Matlab File Exchange
%     (https://se.mathworks.com/matlabcentral/fileexchange/25634-smoothn).
%
% Author:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  times {mustBeVector,mustBeNumeric}
  options.freqRange (1,2) {mustBeVector,mustBePositive} = [2 12]
  options.passbandFilter (1,1) {islogical} = true
  options.stepsize (1,1) {mustBePositive} = 0.002
  options.convolutionPoints (1,1) {mustBeNumeric,mustBePositive} = 25
  options.showPhase (1,1) {islogical} = false
  options.smoothFreq (1,1) {islogical} = false
end

% Get convolved population rate
populationRate.numcells = 1;
populationRate.times = {times};
[spikesPresentation, spikeTimeBins] = ...
  spikes_convolution(populationRate, options.stepsize, options.convolutionPoints);
spikingRate = spikesPresentation./options.stepsize; % Turn into the spiking rate

% Filter the population rate
convolvedSR = round(1/(spikeTimeBins(2)-spikeTimeBins(1)));
if options.passbandFilter
  Wn_WB = [options.freqRange(1)/(convolvedSR/2) ...
    options.freqRange(2)/(convolvedSR/2)]; % Band-pass frequencies normalised by the Nyquist frequency
  [bWB,aWB] = butter(3,Wn_WB); % Apply butterworth filter
  spikingRateFiltered = filtfilt(bWB,aWB,spikingRate)';
else
  spikingRateFiltered = spikingRate(:);
end
spikingRateFiltered = spikingRateFiltered - mean(spikingRateFiltered);

% Obtain wide-band oscillation/fluctuation phase
[phase, instantFreq] = generalized_phase_vector(spikingRateFiltered, ...
  round(1/options.stepsize), options.freqRange(1));
phase = angle(phase);
phaseUnwrapped = unwrap(phase);

% Calculate instantaneous frequency
cycleDurations = diff(find(diff(ceil((phaseUnwrapped+pi)/(2*pi))))); % Cycle endpoints at pi
instantFreqPerCycle = convolvedSR./cycleDurations; % Find points where phase switches from negative to positive and use the distance between these points to calculate the instantaneous frequency
cycleEndpoints = cumsum(cycleDurations)/convolvedSR; % Get corresponding times
instantFreqPerCycle(instantFreqPerCycle > options.freqRange(2)) = nan; % Remove values outside the theta frequency range
if options.smoothFreq
  instantFreqPerCycle = nanconv(instantFreqPerCycle,ce_gausswin(7)/sum(ce_gausswin(7)),'edge'); % Smooth it
end

% Obtain the amplitude of the Hilbert-transformed signal
amplitude = abs(hilbert(spikingRateFiltered));

% Display phase and spiking rate
if options.showPhase
  % Phase
  figure
  yyaxis left
  plot(spikeTimeBins, phase, '-', 'color','#4DBEEE')
  ylabel('Signal phase (rad)')
  
  % Spiking rate
  hold on
  yyaxis right
  plot(spikeTimeBins, spikingRate, 'c-')
  if options.passbandFilter
    plot(spikeTimeBins, spikingRateFiltered - mean(spikingRateFiltered) + mean(spikingRate), 'b-')
  end
  ylabel('(Mean-subtracted) spiking rate (spikes/s)')
  
  title('Spiking Wide-band Phase Graph')
  if options.passbandFilter
    legend('Phase', 'Spiking Rate', 'Filtered spiking rate')
  else
    legend('Phase', 'Spiking Rate')
  end
end

% Assign output
generalisedPhase.phase = phase;
generalisedPhase.phaseUnwrapped = phaseUnwrapped;
generalisedPhase.instantFreq = instantFreq;
generalisedPhase.amplitude = amplitude;
generalisedPhase.spikingRate = spikingRate;
generalisedPhase.filteredSpikes = spikingRateFiltered;
generalisedPhase.timestamps = spikeTimeBins;
generalisedPhase.instantFreqPerCycle = instantFreqPerCycle;
generalisedPhase.cycleEndpoints = cycleEndpoints;
params = rmfield(options,'showPhase');