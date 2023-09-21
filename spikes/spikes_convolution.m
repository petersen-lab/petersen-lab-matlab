function [spikes_presentation,time_bins,parameters] = spikes_convolution(spikes, stepsize, convolution_points, startTime)
% Gaussian convolution of the spikes raster into continuous rates 
% Inputs
%   spikes              % Spikes struct with numcells and times fields
%   stepsize            % step size of continuous traces
%   convolution_points  % points of gaussian convolution (gausswin)
%   startTime           % start time bin (default = stepsize)
%
% Output
%   spikes_presentation
%   time_bins
%   parameters (stepsize, convolution_points, startTime)

if nargin < 2
    % Setting default convolution points (steps)
    stepsize = 0.002;
end
if nargin < 3
    % Setting default convolution points (steps)
    convolution_points = 50;
end
if nargin < 4
  startTime = stepsize;
end
if isfield(spikes,'times')
    spikes.spindices = generateSpinDices(spikes.times);
else
    error('Spikes structure is missing spike times.');
end

% Generating continuous representation of the raster actvity
time_bins = max([1 floor(startTime/stepsize)])*stepsize:stepsize:ceil(max(spikes.spindices(:,1)));
spikes_presentation = zeros(spikes.numcells,numel(time_bins));

for i = 1:spikes.numcells
    idx = round(spikes.times{i}/stepsize) - time_bins(1)/stepsize + 1; % Assign bin indices to spikes
    idx(idx == 0) = 1;
    [spkCounts,idx] = groupcounts(idx); % Count spikes within bin indices
    spikes_presentation(i,idx) = spkCounts; % Mark spike presentations
    
    % Convoluting the spike times with a n points gaussian convolution
    spikes_presentation(i,:) = nanconv(spikes_presentation(i,:),gausswin(convolution_points)'/sum(gausswin(convolution_points)),'edge');
end

% List processing paprameters
parameters.stepsize = stepsize;
parameters.convolution_points = convolution_points;
parameters.startTime = startTime;