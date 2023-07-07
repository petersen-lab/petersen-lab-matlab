function intervals = logical2intervals(logicalVec)
% intervals = logical2intervals(logicalVec)
%
% Function finds onset and offset indices of true values in a logical
% vector.
%
% Args:
%   logicalVec (logical, required, positional): a shape-(1, N) logical
%     array.
%
% Returns:
%   intervals (numeric): a shape-(M, 2) numeric array of indices. Each row
%     corresponds to individual intervals of true values with the first
%     element being the onset index and the second element being the offset
%     index.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  logicalVec {mustBeVector,mustBeA(logicalVec,'logical')}
end

% Parse input
logicalVec = logicalVec(:)';

% Find onsets
changes = [0 diff(logicalVec)];
onsets = changes;
onsets(onsets < 0) = 0;
onsets = find(onsets);

% Find offsets
offsets = -changes;
offsets(offsets < 0) = 0;
offsets = find(offsets) - 1;

% Check for missing onsets
if offsets(1) < onsets(1)
  onsets = [1 onsets];
end

% Check for missing offsets
if numel(offsets) < numel(onsets)
  offsets = [offsets numel(logicalVec)];
end

% Assign output
intervals = [onsets(:) offsets(:)];