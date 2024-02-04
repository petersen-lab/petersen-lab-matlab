function [resampledSpikeCounts, timeBins] = resampleSpikeCounts(spikeCounts, options)
% [resampledSpikeCounts, timeBins] = resampleSpikeCounts(spikeCounts, <options>)
%
% Function resamples a spike count vectors in a row matrix.
%
% Args:
%   spikeCounts (numeric, required, positional): a shape-(M, N) numeric
%     array containing spike counts per sample point. Individual spike
%     count vectors correspond to individual rows in a matrix.
%   stepsize (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     corresponding to the sampling interval used in the spikeCounts vector
%     (default = 0.002).
%   newStepsize (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     corresponding to the new sampling interval (default = 0.05).
%
% Returns:
%   resampledSpikeCounts (numeric) a shape-(M, L) numeric array containing
%     resampled spike count vectors (rows in a matrix).
%   timeBins (numeric) a shape-(1, L) numeric array time bins corresponding
%     to the columns of resampledSpikeCounts matrix.
%
% Dependencies:
%   None.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  spikeCounts (1,:) {mustBeNumeric,mustBeNonempty,mustBeNonnegative}
  options.stepsize (1,1) {mustBeNumeric,mustBePositive} = 0.002
  options.newStepsize (1,1) {mustBeNumeric,mustBePositive} = 0.05
end

% parse input
sr = 1/options.stepsize;
newsr = 1/options.newStepsize;

% Re-bin the spikes
biggestIndex = size(spikeCounts,2);
sampleDuration = ceil(sr/newsr);
newBiggestIndex = ceil(biggestIndex*(newsr/sr));
resampledSpikeCounts = zeros(size(spikeCounts,1), newBiggestIndex);
for s = 1:newBiggestIndex
  iEnd = (s-1)*sampleDuration+sampleDuration;
  if iEnd == biggestIndex+sampleDuration
      iSum = 0;
  elseif iEnd > biggestIndex
    iSum = sum(spikeCounts(:,(s-1)*sampleDuration+1:biggestIndex),2);
    iSum = iSum*(sampleDuration/(biggestIndex-((s-1)*sampleDuration)));
  else
    iSum = sum(spikeCounts(:,(s-1)*sampleDuration+(1:sampleDuration)),2);
  end
  resampledSpikeCounts(:,s) = iSum;
end

% Calculate time bins
timeBins = (1:size(resampledSpikeCounts,1)).*options.newStepsize;