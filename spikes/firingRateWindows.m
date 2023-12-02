function [maxFiringRate, firingRateEvolution] = firingRateWindows(spikeTimes, options)
% [maxFiringRate, firingRateEvolution] = firingRateWindows(spikeTimes, samplingRate, <options>)
%
% Function calculates the maximum firing rate over a given sliding time
% window. It also reports firing rates for all windows across the entirety
% of the recording.
%
% Args:
%   spikeTimes (numeric, required, positional): a shape-(1, N) numeric
%     array of spike times in seconds.
%   samplingRate (numeric, required, positional): a shape-(1, 1) numeric
%     scalar corresponding to the data sampling rate in Hz.
%   windowSize (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     indicating the size of time window in seconds over which the firing
%     rate is calculated (default = 20 minutes or 1200 seconds).
%   stepSize (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     indicating a time step size between two consequtive windows. By
%     default the step size is not fixed but corresponds to consecutive
%     inter-spike intervals.
%   startTime (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     indicating the start time of the first time window (default = 0).
%
% Returns:
%   maxFiringRate (numeric): a shape-(1, 1) numeric scalar indicating the
%     highest firing rate across all time windows (spikes per second).
%   firingRateEvolution (numeric): a shape-(1, M) numeric array containing
%     firing rates for all consequtive time windows across the entirety
%     of the recording.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  spikeTimes (1,:) {mustBeNumeric,mustBeVector,mustBeNonempty}
  options.windowSize (1,1) {mustBeNumeric,mustBePositive} = 1200
  options.stepSize (1,1) {mustBeNumeric} = 0
  options.startTime (1,1) {mustBeNumeric} = 0
end

% Parse input
spikeTimes = spikeTimes(:)';

% Set start and end times for all windows
if options.stepSize
  startTimes = options.startTime:options.stepSize:spikeTimes(end);
  endTimes = options.startTime+options.windowSize:options.stepSize:spikeTimes(end);
  startTimes = startTimes(1:numel(endTimes)+1);
  endTimes = [endTimes spikeTimes(end)];
else
  startTimes = spikeTimes;
  endTimes = startTimes + options.windowSize;
  startTimes = startTimes(endTimes <= spikeTimes(end));
  endTimes = endTimes(endTimes <= spikeTimes(end));
end

% Calculate firing rate for all time windows
firingRateEvolution = zeros(size(startTimes));
nWindows = numel(startTimes);
for window = 1:nWindows
  firingRateEvolution(window) = sum(spikeTimes >= startTimes(window) & spikeTimes <= endTimes(window))/options.windowSize;
end

% Calculate highest firing rate
maxFiringRate = max(firingRateEvolution);