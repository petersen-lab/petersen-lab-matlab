function [resampledSpikes, spikeTimeBins, resampledSpikeTimes] = resampleSpikesArray(spikeTimes, options)
% [resampledSpikes, spikeTimeBins, resampledSpikeTimes] = resampleSpikesArray(spikeTimes, <stepsize>)
%
% Function resamples spike times series.
%
% Args:
%   spikeTimes (cell, required, positional): a shape-(K, 1) cell
%     array of shape-(1, N) numeric arrays of spike times where N
%     corresponds to spike times. N might be different in different cells.
%   stepsize (numeric, optional, keyword): a shape-(1, 1) nunmeric scalar
%     corresponding to the new sample interval (default = 0.002).
%   startTime (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     corresponding to the start time bin (default = stepsize).
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
% Dependencies:
%   petersen-lab/petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  spikeTimes (:,1) {mustBeA(spikeTimes,'cell'),mustBeNonempty}
  options.stepsize (1,1) {mustBeNumeric,mustBePositive} = 0.002
  options.startTime (1,1) {mustBeNumeric} = 0
end


%% Generating continuous resampled spike container
spikeTimeBins = max([1 floor(options.startTime/options.stepsize)])*options.stepsize:options.stepsize:max(cellfun(@(x) getMaxSpikeTime(x),spikeTimes))+options.stepsize;
resampledSpikes = zeros(numel(spikeTimes),numel(spikeTimeBins));
resampledSpikeTimes = cell(numel(spikeTimes),1);


%% Filling in the container variables
for unit = 1:numel(spikeTimes)
  if ~isempty(spikeTimes{unit})
    [resampledSpikesUnit, ~, resampledSpikeTimes{unit}] = resampleSpikes( ...
      spikeTimes{unit}, stepsize=options.stepsize, startTime=options.startTime);
    resampledSpikes(unit,1:numel(resampledSpikesUnit)) = resampledSpikesUnit;
  end
end