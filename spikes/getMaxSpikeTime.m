function maxTime = getMaxSpikeTime(spikeTimes)
% maxTime = getMaxSpikeTime(spikeTimes)
%
% Function picks maximum spike time from a given spike times vector. An
% time pick is replaced by 0.
%
% Args:
%   spikeTimes (numeric, required, positional): a shape-(1, N) numeric
%     array of spike times where N corresponds to individual spike times.
%
% Returns:
%   maxTime (numeric): a shape-(1, 1) numeric scalar corresponding to the
%     maximum spike time. An empty pick is replaced by 0.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  spikeTimes
end

if numel(size(spikeTimes)) > 2 || (size(spikeTimes,1) > 1 && size(spikeTimes,2) > 1)
  error('Value must be a 1-by-n vector or an n-by-1 vector.');
end

maxTime = max(spikeTimes);
if isempty(maxTime)
  maxTime = 0;
end