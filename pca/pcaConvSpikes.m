function [pcaOut, pcaOutFilt, scoreCorr] = pcaConvSpikes(spikeTimes, options)
% [pcaOut, pcaOutFilt, scoreCorr] = pcaConvSpikes(spikeTimes, <options>)
%
% Function performs the principal component analysis (PCA) on unit spiking
% data. PCA is performed on convolved unfiltered and filtered data.
%
% Args:
%   spikeTimes (cell, required, positional): a shape-(M, 1) cell array of
%     numeric spike time vectors.
%   intervals (numeric, optional, keyword): a shape-(N, 2) numeric array
%     of time intervals of interest with rows corresponding to individual
%     intervals, while the first and the second columns corresponding to
%     start and end times, respectively. By default, intervals are not used
%     and the entire duration of the signal is selected (all spike times).
%   includeUnits (logical | numeric, optional, keyword): a shape-(M, 1)
%     logical array indicating which units should be included (true) or
%     excluded (false) from the analysis. Alternatively, a shape-(L, 1)
%     numeric linear index array can be suppplied. By default, all units
%     are included.
%   normalise (logical, optional, keyword): a shape-(1, 1) logical scalar
%     determinig whether individual unit activity should be normalised
%     (default).
%   samplingInterval (numeric, optional, keyword): a shape-(1, 1) numeric
%     scalar setting the sampling interval of the convolved spiking data
%     (default=0.002).
%   convPoints (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     setting the size of the Gaussian convolution window in terms of
%     sample points (default=25). Combined with the default sampling
%     interval value, the convolved-only signal is effectively low-pass
%     filtered at ~20 Hz.
%   freqRange (numeric, optional, keyword): a shape-(1, 2) numeric array of
%     band-pass filter frequency range. PCA is also carried out on a
%     band-pass filtered signal. This parameter allows changing the default
%     frequency pass-band of [4 12] Hz. Setting the passband frequency to
%     an empty array allows to turn off band-pass filtering.
%
% Returns:
%   pcaOut (struct): a shape-(1, 1) scalar structure containing the output
%     of the PCA performed on the convolved signal. The fields are exactly
%     the same as those of the native Matlab's pca function with an
%     additional timestamps field.
%   pcaOutFilt (struct): a shape-(1, 1) scalar structure containing the
%     output of the PCA performed on the convolved and band-pass filtered
%     signal. The fields are exactly the same as those of the native
%     Matlab's pca function with an additional timestamps field.
%   scoreCorr (struct): a shape-(1, 1) scalar structure containing the
%     correlation coefficients and corresponding p-values for a match
%     between band-pass filtered and non-filtered PC scores. Individual
%     entries of fields rho and p correspond to individual components.
%     Spearman correlation analysis is used.
%
% Comments:
%   This is the wrapper of the Matlab's native pca function. To learn more
%   about this function, type 'help pca'.
%
% Dependencies:
%   petersen-lab/petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  spikeTimes (:,1) {mustBeA(spikeTimes,'cell'),mustBeNonempty}
  options.intervals (:,2) {mustBeNumeric} = []
  options.includeUnits (:,1) {mustBeNumericOrListedType(options.includeUnits,'logical')} = (1:numel(spikeTimes))'
  options.normalise (1,1) {mustBeA(options.normalise,'logical')} = true
  options.samplingInterval (1,1) {mustBeNumeric,mustBePositive} = 0.002
  options.convPoints (1,1) {mustBeNumeric,mustBePositive} = 25
  options.freqRange (1,2) {mustBeNumeric} = [4 12]
end

% Convolve data
emptyChSpikeTimes = cellfun('isempty', spikeTimes);
includeUnits = options.includeUnits & ~emptyChSpikeTimes;
lastSpikeTime = max(cellfun(@max, spikeTimes));
[convChSpikeTimes, convTimeBins, convParams] = convolveSpikes( ...
  spikeTimes, stepSize=options.samplingInterval, ...
  convolutionPoints=options.convPoints, endTime=lastSpikeTime);

% Filter data
if ~isempty(options.freqRange)
  filtConvChSpikeTimes = bandpassFilterTimeSeries(convChSpikeTimes, ...
    sampleRate=round(1/convParams.stepSize), frequencyRange=options.freqRange);
end

% Select time intervals
if ~isempty(options.intervals)
  [~, inds] = selectArrayValues(convTimeBins, options.intervals);
else
  inds = round(convTimeBins./options.samplingInterval);
end

% PCA
signal = zeros(size(convChSpikeTimes(:,inds)'));
if options.normalise
  for unit = 1:size(convChSpikeTimes(:,inds)',2)
    signal(:,unit) = (convChSpikeTimes(unit,inds) - ...
      mean(convChSpikeTimes(unit,inds)))./max([realmin std(convChSpikeTimes(unit,inds))]);
  end
else
  for unit = 1:size(convChSpikeTimes(:,inds)',2)
    signal(:,unit) = (convChSpikeTimes(unit,inds) - ...
      mean(convChSpikeTimes(unit,inds)));
  end
end
[pcaOut.coeff, pcaOut.score, pcaOut.latent, pcaOut.tsquared, ...
  pcaOut.explained, pcaOut.mu] = pca(signal(:,includeUnits));
pcaOut.timestamps = convTimeBins(inds)';

if ~isempty(options.freqRange)
  signal = zeros(size(filtConvChSpikeTimes(:,inds)'));
  if options.normalise
    for unit = 1:size(convChSpikeTimes(:,inds)',2)
      signal(:,unit) = (filtConvChSpikeTimes(unit,inds) - ...
        mean(filtConvChSpikeTimes(unit,inds)))./max([realmin std(filtConvChSpikeTimes(unit,inds))]);
    end
  else
    for unit = 1:size(convChSpikeTimes(:,inds)',2)
      signal(:,unit) = (filtConvChSpikeTimes(unit,inds) - ...
        mean(filtConvChSpikeTimes(unit,inds)));
    end
  end
  [pcaOutFilt.coeff, pcaOutFilt.score, pcaOutFilt.latent, pcaOutFilt.tsquared, ...
    pcaOutFilt.explained, pcaOutFilt.mu] = pca(signal(:,includeUnits));
  pcaOutFilt.timestamps = convTimeBins(inds)';
else
  pcaOutFilt.coeff = [];
  pcaOutFilt.score = [];
  pcaOutFilt.latent = [];
  pcaOutFilt.tsquared = [];
  pcaOutFilt.explained = [];
  pcaOutFilt.mu = [];
end

% Quantify match between non-filtered and filtered PC scores
if ~isempty(options.freqRange)
  nComponents = size(pcaOut.coeff,2);
  scoreCorr.rho = zeros(1,nComponents);
  scoreCorr.p = zeros(1,nComponents);
  for comp = 1:nComponents
    [scoreCorr.rho(comp), scoreCorr.p(comp)] = corr( ...
      pcaOut.coeff(:,comp), pcaOutFilt.coeff(:,comp), 'type','Spearman');
  end
else
  scoreCorr.rho = [];
  scoreCorr.p = [];
end