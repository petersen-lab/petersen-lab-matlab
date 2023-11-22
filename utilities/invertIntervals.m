function invertedIntervals = invertIntervals(intervals, endTime, options)
% invertedIntervals = invertIntervals(intervals, endTime, <startTime>)
%
% Function inverts time intervals.
%
% Args:
%   intervals (numeric, required, positional): a shape-(N, 2) numeric
%     array of time intervals. Each row corresponds to individual time
%     intervals with the first element being the start time and the second
%     element being the end time.
%   endTime (numeric, required, positional): a shape-(1, 1) numeric scalar
%     defining the last time sample. No inverse intervals will be produced
%     that exceed this time limit.
%   startTime(numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     defining the first time sample. Any time prior to this point will be
%     excluded (default=0).
%
% Returns:
%   invertedIntervals (numeric): a shape-(M, 2) numeric array of time
%     intervals refering to times outside the boundaries of original
%     intervals.
%
% Dependencies:
%   petersen-lab/petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  intervals (:,2) {mustBeNumeric,mustBeNonnegative}
  endTime (1,1) {mustBeNumeric,mustBePositive}
  options.startTime (1,1) {mustBeNumeric,mustBeNonnegative} = 0
end

% Invert intervals
invertedIntervals = [intervals(1:end-1,2) intervals(2:end,1)];
if intervals(1,1) > 0
  invertedIntervals = [0 intervals(1,1); invertedIntervals];
end
if intervals(end,end) < endTime
  invertedIntervals = [invertedIntervals; intervals(end,end) endTime];
end

% Crop intervals
if ~isempty(invertedIntervals)
  invertedIntervals = invertedIntervals( ...
    logical(sum(invertedIntervals >= options.startTime,2) & ...
    sum(invertedIntervals <= endTime,2)),:);
  if ~isempty(invertedIntervals)
    invertedIntervals(1,1) = max([invertedIntervals(1,1) options.startTime]);
    invertedIntervals(end,end) = min([invertedIntervals(end,end) endTime]);
  end
end