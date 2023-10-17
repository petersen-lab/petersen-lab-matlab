function [binCounts, binLocs, totalCounts, phaseMeans] = phaseHistrogram(phase, options)
% [valueCounts, binLocs] = phaseHistrogram(phase, <options>)
%
% Function produces a phase histogram (bin counts of phase values, not an
% actual graph).
%
% Args:
%   phase (numeric, required, positional): a shape-(M, N) numeric array of
%     phase values in radians. The histogram bins values column-wise.
%   centre (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     representing the centre of the new phase range for the phase
%     histogram (default = 0).
%   nBins (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     indicating the numer of histogram bins (default = 10).
%
% Returns:
%   binCounts (numeric): a shape-(nBins, N) numeric array with phase value
%     counts in each bin.
%   binLocs (numeric): a shape-(1, nBins) numeric array with phase
%     histogram bin locations in radians.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  phase (:,:) {mustBeNumeric}
  options.centre (1,1) {mustBeNumeric,mustBeNonNan,mustBeReal} = 0
  options.nBins (1,1) {mustBeNumeric,mustBePositive} = 10
end

% Recentre phase values
[phase, phaseRange] = recentrePhase(phase, options.centre);

% Work out bin locations
binSize = (phaseRange(2) - phaseRange(1))/options.nBins;
binLocs = phaseRange(1)+binSize/2:binSize:phaseRange(2);

% Bin phase values
if size(phase,1) == 1
  binCounts = hist([phase; phase], binLocs)./2; % deal with the singular case
else
  binCounts = hist(phase, binLocs); %#ok<*HIST>
end

% Calculate phase means and all values per column (excludes NaNs)
[phaseMeans, ~, totalCounts] = datamean(phase, 'circularNP');