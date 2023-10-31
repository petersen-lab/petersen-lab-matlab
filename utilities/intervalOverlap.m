function overlappingIntervals = intervalOverlap(intervals1, intervals2)
% overlappingIntervals = intervalOverlap(intervals1, interval2)
%
% Function selects overlapping intervals.
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
% Dependencies:
%   petersen-lab/petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).
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
if ~isempty(intervals1) && ~isempty(intervals2)
  for interv = 1:size(intervals2,1)

    % Cut out non-overlapping outer times
    intervals1_i2 = intervals1;
    intervals1_i2(:,1) = max([intervals1_i2(:,1) ...
      repmat(min(min(intervals2(interv,:))),size(intervals1_i2(:,1)))],[],2);
    intervals1_i2(:,2) = min([intervals1_i2(:,2) ...
      repmat(max(max(intervals2(interv,:))),size(intervals1_i2(:,2)))],[],2);
    intervals1_i2 = intervals1_i2(intervals1_i2(:,1) < intervals1_i2(:,2),:);

    % Cut out non-overlapping inner times
    if ~isempty(intervals1_i2)
      startTimes = selectArrayValues(intervals1_i2(:,1), intervals2(interv,:));
      stopTimes = selectArrayValues(intervals1_i2(:,2), intervals2(interv,:));
      if ~isempty(startTimes) || ~isempty(stopTimes)
        if numel(startTimes) > numel(stopTimes)
          stopTimes = [stopTimes; intervals2(interv,2)]; %#ok<*AGROW>
        elseif numel(startTimes) < numel(stopTimes)
          startTimes = [intervals2(interv,1); startTimes];
        elseif numel(startTimes) == numel(stopTimes) && startTimes(1) > stopTimes(1)
          startTimes = [intervals2(interv,1); startTimes];
          stopTimes = [stopTimes; intervals2(interv,2)];
        end
        assert(numel(startTimes)==numel(stopTimes) && ...
          startTimes(1)<stopTimes(1) && startTimes(end)<stopTimes(end));
        startTimesAll = [startTimesAll; startTimes];
        stopTimesAll = [stopTimesAll; stopTimes];
      end
    end
  end
end

% Assign output
overlappingIntervals = [startTimesAll(:) stopTimesAll(:)];