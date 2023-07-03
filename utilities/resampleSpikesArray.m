function [resampledSpikes, spikeTimeBins, resampledSpikeTimes] = resampleSpikesArray(spikeTimes, options)
% [resampledSpikes, spikeTimeBins, resampledSpikeTimes] = resampleSpikesArray(spikeTimes, <stepsize>)
%
% Function resamples spike times series.
%
% Args:
%   spikeTimes (cell, required, positional): a shape-(K, 1) cell
%     array of shape-(1, N) numeric arrays of spike times where N
%     corresponds to spike times.
%   stepsize (numeric, optional, keyword): a shape-(1, 1) nunmeric scalar
%     corresponding to the new sample interval (default = 0.002).
% Returns:
%   resampledSpikes (numeric): a shape-(K, 1) numeric array of shape-(1, M)
%     resampled spike counts with the time bin of the same size as the
%     input stepsize variable.
%   spikeTimeBins (numeric): a shape-(1, M) numeric array with time bins
%     corresponding to the length of individual elements of resampledSpikes
%     output variable.
%   resampledSpikeTimes (numeric): a shape-(K, 1) numeric array of
%     shape-(1, L) numeric arrays of resampled spike times.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  spikeTimes (:,1) {mustBeA(spikeTimes,'cell')}
  options.stepsize (1,1) {mustBeNumeric,mustBePositive} = 0.002
end


%% Generating continuous resampled spike container
spikeTimeBins = 0:options.stepsize:max(cellfun(@(x) max(x),spikeTimes))+options.stepsize;
resampledSpikes = zeros(numel(spikeTimes),numel(spikeTimeBins));
resampledSpikeTimes = cell(numel(spikeTimes),1);


%% Filling in the container variables
for unit = 1:numel(spikeTimes)
  [resampledSpikesUnit, spikeTimeBins, resampledSpikeTimes{unit}] = resampleSpikes(spikeTimes{unit}, stepsize=options.stepsize);
  resampledSpikes(unit,1:numel(resampledSpikesUnit)) = resampledSpikesUnit;
end