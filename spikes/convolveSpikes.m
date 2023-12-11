function [spikesPresentation, timeBins, parameters] = convolveSpikes(spikeTimes, options)
% [spikes_presentation, timeBins, parameters] = convolveSpikes(spikeTimes, <options>)
%
% Gaussian convolution of the spikes raster into continuous rates.
%
% Args:
%   spikeTimes (cell | numeric, required, positional): a shape-(N, 1) cell
%     array of numeric spike time arrays corresponding to individual units
%     or probe recording channels or a single shape-(1, M) numeric array of
%     spike times.
%   stepSize (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     with the sampling interval for convolved continuous spike traces
%     (default=0.002).
%   convolutionPoints (numeric, optional, keyword): a shape-(1, 1) numeric
%     scalar with Gaussian convolution sample points (gausswin;
%     default=25).
%   startTime (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     with the start time bin (default = stepSize).
%   endTime (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     with the end time bin (default = last spike time).
%
% Returns:
%   spikesPresentation (numeric): a shape-(N, M) numeric array of convolved
%     continuous spike time arrays corresponding to individual units or
%     probe recording channels. If the input variable spikeTimes contains
%     empty cells, the corresponding row vector of spikesPresentation would
%     contain NaNs.
%   timeBins (numeric): a shape-(1, M) numeric time array corresponding to
%     columns of spikesPresentation.
%   parameters (struct): a shape-(1, 1) scalar structure with the following
%     fields: stepSize, convolutionPoints, and startTime. These fields
%     contain user supplied input parameters or their default values.
%
% Dependencies:
%   petersen-lab/petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab).
%   CellExplorer (https://github.com/petersenpeter/CellExplorer).
%
% Comments:
%   This is a more generic version of
%   petersen-lab-matlab/spikes/spikes_convolution.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  spikeTimes (:,:) {mustBeNumericOrListedType(spikeTimes,'cell')}
  options.stepSize (1,1) {mustBeNumeric,mustBePositive} = 0.002
  options.convolutionPoints (1,1) {mustBeNumeric,mustBePositive} = 25
  options.startTime (:,:) {mustBeNumeric,mustBeScalarOrEmpty} = []
  options.endTime (:,:) {mustBeNumeric,mustBeScalarOrEmpty} = []
end

% Parse input
if ~iscell(spikeTimes)
  spikeTimes = {spikeTimes};
end
if isempty(options.startTime) || ~isnumeric(options.startTime)
  options.startTime = options.stepSize;
end

% Find data limits
nDataRows = numel(spikeTimes);
if isempty(options.endTime)
  options.endTime = ceil(max(cellfun(@(x) getMaxSpikeTime(x), spikeTimes)));
else
  assert(options.endTime >= max(cellfun(@(x) getMaxSpikeTime(x), spikeTimes)));
end

% Generate continuous representation of raster actvity
timeBins = max([1 floor(options.startTime/options.stepSize)])*options.stepSize:options.stepSize:options.endTime+options.stepSize;
spikesPresentation = zeros(numel(spikeTimes), numel(timeBins));

for row = 1:nDataRows
  if isempty(spikeTimes{row})
    spikesPresentation(row,:) = nan(size(timeBins));
  else
    idx = ceil(spikeTimes{row}/options.stepSize) - timeBins(1)/options.stepSize + 1; % Assign bin indices to spikes
    idx(idx == 0) = 1;
    [spkCounts,idx] = groupcounts(idx(:)); % Count spikes within bin indices
    spikesPresentation(row,idx) = spkCounts; % Mark spike presentations

    % Convolving the spike times with an n-point gaussian function
    spikesPresentation(row,:) = nanconv(spikesPresentation(row,:), ...
      gausswin(options.convolutionPoints)'/sum(gausswin(options.convolutionPoints)),'edge');
  end
end
assert(numel(timeBins) == size(spikesPresentation,2));

% List processing paprameters
parameters.stepSize = options.stepSize;
parameters.convolutionPoints = options.convolutionPoints;
parameters.startTime = options.startTime;
parameters.endTime = options.endTime;