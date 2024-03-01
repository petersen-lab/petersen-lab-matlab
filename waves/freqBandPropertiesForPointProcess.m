function [freqBandProperties, params] = freqBandPropertiesForPointProcess(times, freqRange, options)
% [freqBandProperties, params] = freqBandPropertiesForPointProcess(times, freqRange, <options>)
%
% Calculates various properties for a given frequency band including a
% spectrogram, instantaneous phase, instantaneous frequency, oscillation
% amplitude, and power for a given point process signal using multi-taper
% estimation and the Chronux toolbox (http://chronux.org/).
%
% Args:
%   times (numeric, required, positional): a shape-(M, N) numeric array
%     where either M = 1 or N = 1. This input vector contains times when
%     discrete events occur (i.e., spike times).
%   freqRange (numeric, required, positional): a shape-(1, 2) numeric array
%     defining the frequency band limits in Hz.
%   sr (numeric, optional, keyword): a shape-(1, 1) scalar corresponding to
%     the data sampling rate (in Hz) (default = 500Hz).
%   window (numeric, optional, keyword): a shape-(1, 1) scalar representing
%     duration (in s) of the spectrogram time window for frequency band power
%     calculation (default = 1).
%   overlap (numeric, optional, keyword): a shape-(1, 1) scalar
%     representing overlap (in s) between successive spectrogram windows
%     (default = 0.5). Ideally, should be equal to window/2.
%   step (numeric, optional, keyword): a shape-(1, 1) scalar representing
%     step (in s) between successive spectrogram windows (default = 0.5).
%     Ideally, should be equal to window/2.
%   tapers (numeric, optional, keyword): a shape-(1, 2) vector array. The
%     values give relative resolution and order of the tapers [NW K] used
%     in spectrogram calculatios (default = [3 5]).
%   pad (numeric, optional, keyword): a shape-(1, 1) scalar representing
%     FFT padding (can take values -1,0,1,2...). -1 corresponds to no
%     padding, 0 corresponds to padding to the next highest power of 2 etc.
%     e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad
%     the FFT to 512 points, if pad=1, we pad to 1024 points etc. Setting
%     pad to larger values, gives denser spectrogram frequency grids.
%     Defaults to 0.
%   showSpectrogram (logical, optional, keyword): a shape-(1, 1) logical
%     scalar for ploting the spectrogram (default = false).
%   parent (figure handle, optional, keyword): a shape-(1, 1) parent figure
%     or uipanel handle (default = figure, which creates a new figure).
%   cutoffs (numeric, optional, keyword): a shape-(1, 2) vector array with
%     cutoff values for the spectrogram color plot (default = [0 13]).
%   filterSpectrogram (logical, optional, keyword): a shape-(1, 1) logical
%     scalar for carrying out spectrogram filtering (default = false).
%   showPower (logical, optional, keyword): a shape-(1, 1) logical scalar
%     for ploting the frequency band power (default = false). If
%     'showSpectrogram' input parameter is set to true, power is overlayed
%     the spectrogram and scaled to the figure window size having arbitrary
%     units.
%   stepsize (numeric, optional, keyword): a shape-(1, 1) scalar
%     representing the step size (in seconds) used for binning spikes prior
%     to producing a continuous spiking rate trace (default = 0.002).
%   convolutionPoints (numeric, optional, keyword): a shape-(1, 1) scalar
%     representing points of gaussian convolution (gausswin) used to smooth
%     binned spiking rate (default = 25 sample points).
%   showPhase (logical, optional, keyword): a shape-(1, 1) logical scalar
%     for ploting the frequency band oscillation phase (default = false). As part of
%     displaying the pahse, the frequency band power, the convolved
%     spiking rate and the filtered spiking rate are also displayed.
%   smoothFrequencies (logical, optional, keyword): a shape-(1, 1) logical
%     scalar for smoothing instantaneous frequency estimates
%     (default = false).
%
% Returns:
%   freqBandProperties (struct): a structure with the following fields:
%     spectrogram (numeric): a shape-(M, N) numeric array containing the
%       spectrogram. M corresponds to frequencies and N corresponds to
%       time.
%     spectrogramTimestamps (numeric): a shape-(1, N) numeric array with
%       spectrogram timestamps.
%     spectrogramFrequencies (numeric): a shape-(M, 1) numeric array with
%       spectrogram frequencies.
%     spectrogramPower (numeric): a shape-(1, N) numeric array with power
%       values of maximal frequencies at spectrogram timestamps.
%     moderatePowerPeriods (logical): a shape-(1, N) logical array with
%       true values representing points when frequency band power is
%       moderate (above the mean power but below mean + norminv(0.95)*SD).
%     highPowerPeriods (logical): a shape-(1, N) logical array with true
%       values representing points when frequency band power is high (above
%       mean + norminv(0.95)*SD).
%     spectrogramMaxFrequency (numeric): a shape-(1, N) numeric array with
%       maximal frequencies at spectrogram timestamps.
%     instFrequency (numeric): a shape-(1, K) numeric array of
%       instantaneous frequencies.
%     instTimestamps (numeric): a shape-(1, K) numeric array of timestamps
%       corresponding to instantaneous frequency estimates.
%     instPhase (numeric): a shape-(1, L) numeric array of frequency band
%       oscillation phases estimated using Hilbert transform. The number of
%       samples is determined by the stepsize input parameter and
%       correspond to the length of the convolved spiking rate vector.
%     instPhaseUnwrapped (numeric): a shape-(1, L) numeric array of
%       unwrapped frequency band oscillation phases estimated using Hilbert
%       transform.
%     amplitude (numeric): a shape-(1, L) numeric array of amplitude
%       (envelope) of the filtered and Hilbert-transformed spiking signal.
%     instPhaseTimestamps (numeric): a shape-(1, L) numeric array of
%       timestamps corresponding to frequency band oscillation phase
%       estimate array.
%     spikingRate (numeric): a shape-(1, L) numeric array of convolved
%       population spiking rate (in spikes/s; continuous).
%     convolvedFilteredSpikes (numeric): a shape-(1, L) numeric array of
%       convolved and filtered population spiking rate (continuous). The
%       band-pass filtering frequency is between 4 and 11 Hz.
%   params (struct): a structure with its fields corresponding to the
%     function input parameters including the default values.
%
% Dependencies:
%   CellExplorer (https://cellexplorer.org/).
%   petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).
%   FMA Toolbox (https://github.com/michael-zugaro/FMAToolbox).
%   Chronux Toolbox (http://chronux.org/).
%   The latter two toolboxes are contained within buzcode/externalPackages
%     (https://github.com/buzsakilab/buzcode).
%
% Comments:
%   The time displacement between successive short time spectra can be
%     supplied either as a 'step' (explicit time difference) or as an
%     'overlap' (between successive time windows).
%   The function was adapted from FMA Toolbox MTPointSpectrogram function.
%     Type 'help MTPointSpectrogram' for more info.
%
% Examples:
%   freqBandProperties = freqBandPropertiesForPointProcess( ...
%     populationRate.times{1}, freqRange=[1 4], pad=3);
%   freqBandProperties = freqBandPropertiesForPointProcess( ...
%     populationRate.times{1}, freqRange=[1 4], pad=3, ...
%     showSpectrogram=true, showPower=true, showPhase=true);
%
% Author:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  times {mustBeVector}
  freqRange {mustBeVector,mustBePositive}
  options.sr (1,1) {mustBeNumeric} = 500
  options.window (1,1) {mustBeNumeric} = 1
  options.overlap (1,1) {mustBeNumeric} = 0.5
  options.step (1,1) {mustBeNumeric} = 0.5
  options.tapers (1,2) {mustBeVector} = [3 5]
  options.pad (1,1) {mustBeNumeric} = 0
  options.showSpectrogram (1,1) {islogical} = false
  options.parent (1,1) {mustBeUnderlyingType(options.parent,'matlab.ui.Figure')} = figure
  options.cutoffs (1,2) {mustBeVector} = [0 13]
  options.filterSpectrogram (1,1) {islogical} = false
  options.stepsize (1,1) {mustBeNumeric} = 0.002
  options.convolutionPoints (1,1) {mustBeNumeric} = 25
  options.showPower (1,1) {islogical} = false
  options.showPhase (1,1) {islogical} = false
  options.smoothFrequencies (1,1) {islogical} = false
end

% Parse input
if options.showSpectrogram
  options.show = 'on';
else
  close(gcf);
  options.parent = [];
  options.show = 'off';
end

% Spectrogram
[spectrogram, t, f] = MTPointSpectrogram( ...
  times, 'frequency',options.sr, 'range',freqRange, ...
  'window',options.window, 'overlap',options.overlap, 'step',options.step, ...
  'tapers',options.tapers, 'pad',options.pad, 'show',options.show, ...
  'parent',options.parent, 'cutoffs',options.cutoffs);

% Filter noise in the spectrogram
if options.filterSpectrogram
  spectrogram = medfilt2(abs(spectrogram),[2,4]); % Median filter with the window size of 2-by-4
  %h = 1/10*ones(4,1);
  %H = h*h';
  %spectrogram = filter2(H,spectrogram); % Filter the spectrogram again: No longer applied
end

% Calculate power and frequency
[power, maxIndex] = max(spectrogram);
thetaFrequency = f(maxIndex);

% Get convolved population rate
populationRate.numcells = 1;
populationRate.times = {times};
[spikesPresentation, spikeTimeBins] = spikes_convolution(populationRate, options.stepsize, options.convolutionPoints);
spikingRate = spikesPresentation./options.stepsize; % Turn into the spiking rate

% Filter the population rate
convolvedSR = round(1/(spikeTimeBins(2)-spikeTimeBins(1)));
Wn_theta = [freqRange(1)/(convolvedSR/2) freqRange(2)/(convolvedSR/2)]; % Band-pass frequencies normalised by the Nyquist frequency
[btheta,atheta] = butter(3,Wn_theta); % Apply butterworth filter
spikingRateFiltered = filtfilt(btheta,atheta,spikingRate)';

% Obtain frequency band oscillation phase
hilbert1 = hilbert(spikingRateFiltered); % Apply Hilbert transform
populationRatePhase = atan2(imag(hilbert1), real(hilbert1));
populationRatePhaseUnwrapped = unwrap(populationRatePhase);

% Obtain the amplitude of the Hilbert-transformed signal
amplitude = abs(hilbert1);

% Calculate instantaneous frequency
%cycleDurations = diff(find(diff(populationRatePhase>0)==1)); % Cycles endpoints at 0 rad crossings
cycleDurations = diff(find(diff(ceil((populationRatePhaseUnwrapped+pi)/(2*pi))))); % Cycle endpoints at pi (more appropriate)
instFreq = (convolvedSR)./cycleDurations; % Find points where phase switches from negative to positive and use the distance between these points to calculate the instantaneous frequency
instTime = cumsum(cycleDurations)/convolvedSR; % Get corresponding times
instFreq(instFreq>freqRange(2)) = nan; % Remove values outside the frequency band limits
if options.smoothFrequencies
  instFreq = nanconv(instFreq,ce_gausswin(7)/sum(ce_gausswin(7)),'edge'); % Not sure why this is needed, makes it less instantaneous :)
end

% Find periods with moderate and high frequency band power
meanPower = mean(power);
stdPower = std(power);
highPowerThr = meanPower + norminv(0.95)*stdPower;
highPowerPeriods = power > highPowerThr;
moderatePowerPeriods = power > meanPower & ~highPowerPeriods;

% Display power
if options.showPower
  if options.showSpectrogram
    gcf
  else
    figure
  end
  hold on
  yyaxis right
  plot(t, power, 'k-')
  ylabel('Power (db/Hz)')
  hold off
end

% Display phase and spiking rate
if options.showPhase
  % Phase
  figure
  yyaxis left
  plot(spikeTimeBins, populationRatePhase, '-', 'color','#4DBEEE')
  ylabel('Oscillation phase (rad)')
  
  % Spiking rate
  hold on
  yyaxis right
  plot(spikeTimeBins, spikingRate, 'c-')
  plot(spikeTimeBins, spikingRateFiltered - mean(spikingRateFiltered) + mean(spikingRate), 'b-')
  yLim = ylim;
  scaledPower = yLim(1) + yLim(2)*(power./max(power));
  plot(t, scaledPower, 'g-')
  plot(t(find(moderatePowerPeriods,1)), scaledPower(find(moderatePowerPeriods,1)), 'y-')
  plot(t(find(highPowerPeriods,1)), scaledPower(find(highPowerPeriods,1)), 'r-')
  for i = 1:numel(t)-1
    if moderatePowerPeriods(i) && moderatePowerPeriods(i+1)
      plot([mean(t(max([1 i-1]):i)) t(i:i+1)'], [mean(scaledPower(max([1 i-1]):i)) scaledPower(i:i+1)], 'y-')
    elseif moderatePowerPeriods(i)
      plot([mean(t(i-1:i)) t(i) mean(t(i:i+1))], [mean(scaledPower(i-1:i)) scaledPower(i) mean(scaledPower(i:i+1))], 'y-')
    end
    if highPowerPeriods(i) && highPowerPeriods(i+1)
      plot([mean(t(max([1 i-1]):i)) t(i:i+1)'], [mean(scaledPower(max([1 i-1]):i)) scaledPower(i:i+1)], 'r-')
    elseif highPowerPeriods(i)
      plot([mean(t(i-1:i)) t(i) mean(t(i:i+1))], [mean(scaledPower(i-1:i)) scaledPower(i) mean(scaledPower(i:i+1))], 'r-')
    end
  end
  ylabel('Convolved spiking rate (spikes/s)')
  
  title('Oscillation Power and Phase Graph')
  legend('Phase', 'Spiking Rate', 'Filtered spiking rate', 'Low power periods', 'Moderate power periods', 'High power periods')
end


% Assign output
freqBandProperties.spectrogram = spectrogram;
freqBandProperties.spectrogramTimestamps = t';
freqBandProperties.spectrogramFrequencies = f;
freqBandProperties.spectrogramPower = power;
freqBandProperties.moderatePowerPeriods = logical(moderatePowerPeriods);
freqBandProperties.highPowerPeriods = logical(highPowerPeriods);
freqBandProperties.spectrogramMaxFrequency = thetaFrequency';
freqBandProperties.instFrequency = instFreq';
freqBandProperties.instTimestamps = instTime';
freqBandProperties.instPhase = populationRatePhase';
freqBandProperties.instPhaseUnwrapped = populationRatePhaseUnwrapped';
freqBandProperties.amplitude = amplitude';
freqBandProperties.instPhaseTimestamps = spikeTimeBins;
freqBandProperties.spikingRate = spikingRate;
freqBandProperties.spikingRateFiltered = spikingRateFiltered';
params = rmfield(options,'show');
params.freqRange = freqRange;