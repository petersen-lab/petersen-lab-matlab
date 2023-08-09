function contaminationPercentage = refractoryContamination(spikeTimes)
% contaminationPercentage = refractoryContamination(spikeTimes)
%
% Function calculates the contamination percentage of the spiking
% autocorrelagram (ACG) 1.5 ms refractory period around zero lag relative
% to the ACG shoulder (1000-1500 ms away from the zero lag).
%
% Args:
%   spikeTimes (numeric, required, positional): a shape-(1, N) numeric
%     array of spike times in seconds.
%
% Returns:
%   contaminationPercentage (numeric): a shape-(1, 1) numeric scalar with
%     the refractory period contamination percentage.
%
% Dependencies:
%   CellExplorer (https://cellexplorer.org/).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  spikeTimes (1,:) {mustBeNumeric,mustBeVector}
end

% Parse input
spikeTimes = spikeTimes(:)';

% Calculate the autocorrelogram (ACG)
binSize = 0.0005; % in seconds
duration = 3; % in seconds
[acg,t] = CCG(spikeTimes, ones(size(spikeTimes)), 'binSize',binSize, 'duration',duration); % Calls CellExplorer CCG function (make sure no conflict with other CCG instances on the Matlab path)

% Calculate the refractory period contamination percentage relative to the ACG shoulder.
contaminationPercentage = mean(acg(t>=-0.0015 & t<=0.0015))/mean(acg(t>=1 & t<=1.5));