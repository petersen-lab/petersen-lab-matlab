function overlappingIntervals = intervalOverlap(intervals1, intervals2)
% overlappingIntervals = intervalOverlap(intervals1, interval2)
%
% Functtion selects overlapping intervals.
%
% Args:
%   intervals1 (numeric, required, positional): a shape-(N, 2) numeric
%     array of time intervals. Each row corresponds to individual time
%     intervals with the first element being the start time and the second
%     element being the end time.
%   intervals2 (numeric, required, positional): a shape-(M, 2) numeric
%     array of time intervals. Each row corresponds to individual time
%     intervals with the first element being the start time and the second
%     element being the end time.
%
% Returns:
%   overlappingIntervals (numeric): a shape-(L, 2) numeric array of time
%     intervals selected from intervals1 so that they all overlap with any
%     of the intervals in intervals2.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  intervals1 (:,2) {mustBeNumeric,mustBeNonnegative}
  intervals2 (:,2) {mustBeNumeric,mustBeNonnegative}
end

% Select overlapping intervals
startTimesAll = [];
stopTimesAll = [];
for interv = 1:size(intervals2,1)
  startTimes = selectArrayValues(intervals1(:,1), intervals2(interv,:));
  stopTimes = selectArrayValues(intervals1(:,2), intervals2(interv,:));
  if ~isempty(startTimes) || ~isempty(stopTimes)
    if numel(startTimes) > numel(stopTimes)
      stopTimes = [stopTimes; intervals2(interv,2)]; %#ok<*AGROW>
    elseif numel(startTimes) < numel(stopTimes)
      startTimes = [intervals2(interv,1); startTimes];
    elseif numel(startTimes) == numel(stopTimes) && startTimes(1) > stopTimes(1)
      startTimes = [intervals2(interv,1); startTimes];
      stopTimes = [stopTimes; intervals2(interv,2)];
    end
    assert(startTimes(1) < stopTimes(1) && startTimes(end) < stopTimes(end));
    startTimesAll = [startTimesAll; startTimes];
    stopTimesAll = [stopTimesAll; stopTimes];
  end
end

% Assign output
overlappingIntervals = [startTimes(:) stopTimes(:)];