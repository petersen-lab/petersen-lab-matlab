function thetaPhaseSpatialMapInterval(intervals, unitSpikeTimes, ...
  populationSpikeTimes, maxChan, figTitle, resamplingInterval, ...
  frequencyRange, options)
% thetaPhaseSpatialMapInterval(intervals, unitSpikeTimes, ...
%   populationSpikeTimes, maxChan, figTitle, resamplingInterval, ...
%   frequencyRange, <options>)
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
%   figPath (char, optional, keyword): a shape-(1, J) character array
%     specifying the full folder path for saving figures. If left empty,
%     figures are not saved (default = []).
%
% Returns:
%   None.
%
% Comments:
%   The top left corner shows the ratio of units with significant phase and
%   the total numer of units with non-zero spike counts. The bottom
%   rightcorner shows the correlation coeffiecient and its significance
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
  options.figPath (1,:) {mustBeNonzeroLengthText} = [];
end

%% Parse input
if strcmpi(options.coherenceRange,'aboveMean')
  n = 51;
elseif strcmpi(options.coherenceRange,'high')
  n = 76;
end

%% Work out intervals of interest
if ~isempty(intervals)
  unitSpikeTimes = selectCellValues(unitSpikeTimes, intervals);
  populationSpikeTimes = selectArrayValues(populationSpikeTimes, intervals);
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

  % plot phase spatial map
  thetaPhaseSpatialMap(fullCoherence.phase, maxChan, ...
    fullCoherence.frequency(1,:), unitSpikeTimes, figTitle, ...
    include=includeUnits, figPath=options.figPath);
end