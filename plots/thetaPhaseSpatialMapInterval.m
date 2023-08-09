function fullCoherence = thetaPhaseSpatialMapInterval(intervals, ...
  unitSpikeTimes, populationSpikeTimes, maxChan, figTitle, ...
  resamplingInterval, frequencyRange, options)
% fullCoherence = thetaPhaseSpatialMapInterval(intervals, ...
%   unitSpikeTimes, populationSpikeTimes, maxChan, figTitle, ...
%   resamplingInterval, frequencyRange, <options>)
%
% A wrapper function for thetaPhaseSpatialMap function. Function calculates
% phase-electrode correlations and related statistics and plots the
% correlations for phase data within specified time intervals.
%
% Args:
%   intervals (numeric, required, positional): a shape-(N, 2) numeric array
%     of time intervals. Each row corresponds to individual time intervals
%     of interest with the first element being the start time and the
%     second element being the end time.
%   unitSpikeTimes (cell, required, positional): a shape-(M, 1) cell array
%     of shape-(1, L) numeric arrays of spike times where M corresponds to
%     individual units and L corresponds to spike times.
%   populationSpikeTimes (numeric, required, positional): a shape-(1, L)
%     numeric array of population (units + MUAs) spike times where L
%     corresponds to spike times.
%   maxChan (numeric, required, positional): a shape-(M, 1) numeric array
%     of recording channels with the largest amplitude waveforms for
%     corresponding units.
%   figTitle (char, required, positional): a shape-(1, L) character array
%     with the title for individual figures. The title will be appended by
%     the frequency value of a figure.
%   resamplingInterval (numeric, required, positional): a shape-(1, 1)
%     numeric scalar representing bin size in seconds for generating spike
%     count vectors.
%   frequencyRange (numeric, required, positional): a shape-(1, 2) numeric
%     array with the frequency range for estimating coherence and phase
%     values.
%   coherenceRange (char, optional, keyword): a shape-(1, K) character
%     array controlling coherence range for including units to correlation
%     analyses. Can take one of the two values:
%       'full' - include all units irrespictive of their coherence values
%                (default).
%       'aboveMean' - include only the top 50% most coherent units.
%       'high' - include only the top 25% most coherent units.
%   parallelise (logical, optional, keyword): a shape-(1, 1) logical scalar
%     for performing coherence analysis using the parfor loop
%     (default = true).
%   include (logical, optional, keyword): a shape-(M, 1) logical array
%     matching the shape of the unitSpikeTimes array and marking units to
%     be included in the analysis correlation analyses. By default, all
%     values are included.
%   figPath (char, optional, keyword): a shape-(1, J) character array
%     specifying the full folder path for saving figures. If left empty,
%     figures are not drawn and not saved (default = '').
%
% Returns:
%   fullCoherence (struct): a shape-(1, 1) scalar structure with the
%     following fields:
%     coherence (numeric): a shape-(G, H) numeric array containing
%       coherence values for the signal with respect to the reference
%       (range = [0 1]).
%     coherenceConf (numeric): a shape-(G, H) numeric array containing
%       coherence 95% confidence interval. Add/subtract this interval to
%       actual coherence values to get upper and lower intervals
%       (coherence +/- coherenceConf).
%     phase (numeric): a shape-(G, H) numeric array containing phase radian
%       values for the signal with respect to the reference. Negative phase
%       indicates lag, whereas positive phase indicates lead.
%     phaseConf (cell): a shape-(G, 1) cell array containing phase upper
%       and lower 95% confidence intervals (rad). Each cell contains a
%       shape-(2, L) numeric array.
%     frequency (numeric): a shape-(G, H) numeric array containing
%       frequency values corresponding to coherence and phase estimates.
%     rateAdjustedCoherence (numeric): a shape-(G, H) numeric array
%       containing firing rate adjusted coherence.
%     rateAdjustedCoherenceConf (numeric): a shape-(G, H) numeric array
%       containing firing rate adjusted coherence 95% confidence interval.
%       Add/subtract this interval to actual coherence values to get upper
%       and lower intervals
%       (rateAdjustedCoherence +/- rateAdjustedCoherenceConf).
%     kappaSignal (numeric): a shape-(G, H) numeric array with coherence
%       firing rate adjustment factor kappa1 corresponding to the primary
%       signal.
%     kappaReference (numeric): a shape-(G, H) numeric array with
%       coherence firing rate adjustment factor kappa2 corresponding to the
%       secondary (reference) signal.
%
% Comments:
%   The top left corner shows the ratio of units with significant phase and
%   the total numer of units with non-zero spike counts. The bottom
%   right corner shows the correlation coeffiecient and its significance
%   p-value.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com)

arguments
  intervals (:,2) {mustBeNumeric,mustBeNonnegative}
  unitSpikeTimes (:,1) {mustBeA(unitSpikeTimes,'cell')}
  populationSpikeTimes {mustBeVector,mustBeNonnegative,mustBeNonempty}
  maxChan (:,1) {mustBeVector,mustBeNonnegative,mustBeNonempty}
  figTitle (1,:) {mustBeNonzeroLengthText}
  resamplingInterval (1,1) {mustBeNumeric,mustBePositive}
  frequencyRange (1,2) {mustBeNumeric,mustBeNonnegative}
  options.coherenceRange {mustBeMember(options.coherenceRange,{'full','aboveMean','high'})} = 'full'
  options.parallelise (1,1) {mustBeA(options.parallelise,'logical')} = true
  options.include (:,:) {mustBeA(options.include,'logical')} = []
  options.figPath (1,:) {mustBeText} = '';
end

%% Parse input
if strcmpi(options.coherenceRange,'aboveMean')
  n = 51;
elseif strcmpi(options.coherenceRange,'high')
  n = 76;
end

if isempty(options.include)
  options.include = true(size(unitSpikeTimes));
end

%% Work out intervals of interest
if ~isempty(intervals)
  unitSpikeTimes = selectCellValues(unitSpikeTimes, intervals);
  populationSpikeTimes = selectArrayValues(populationSpikeTimes, intervals);
else
  unitSpikeTimes = [];
  populationSpikeTimes = [];
end

%% Carry out coherence and correlation analyses
if isempty(unitSpikeTimes) || isempty(populationSpikeTimes)
  warning('Empty spike count vectors supplied for coherence analysis. Quiting...');
else
  % Calculate coherence
  fullCoherence = coherence(unitSpikeTimes, populationSpikeTimes, ...
    intervals=intervals, stepsize=resamplingInterval, ...
    range=frequencyRange, parallelise=options.parallelise);

  % Label coherent units
  if strcmpi(options.coherenceRange,'full')
    includeUnits = true(size(fullCoherence.coherence));
  else
    quantiles = quantile(fullCoherence.coherence,100);
    coherenceThr = quantiles(n,:);
    includeUnits = false(size(fullCoherence.coherence));
    for f = 1:numel(coherenceThr)
      includeUnits(fullCoherence.coherence(:,f) >= coherenceThr(f),f) = true;
    end
  end
  options.include = repmat(options.include, 1, size(fullCoherence.coherence,2));
  includeUnits = logical(options.include & includeUnits);
  if ~any(includeUnits)
    return
  end

  % plot phase spatial map
  if ~isempty(options.figPath)
    thetaPhaseSpatialMap(fullCoherence.phase, maxChan, ...
      fullCoherence.frequency(1,:), unitSpikeTimes, figTitle, ...
      include=includeUnits, figPath=options.figPath);
  end
end