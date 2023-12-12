function sessionIntervals = getTimeIntervals(dataFile, epochsOfInterest, options)
% sessionIntervals = getTimeIntervals(dataFile, epochsOfInterest, <options>)
%
% Function extracts time intervals of interest for a single recording
% session (located in the ms-preliminary repository). The function works
% with data saved in the CellExplorer format (https://cellexplorer.org/).
%
% Args:
%   dataFile (char, required, positional): a shape-(1, M) character array
%     with the generic data filename in the following format:
%     <data-folder>/<animalID>/<session-basename>/<session-basename>.*.mat.
%   epochsOfInterest (char | cell, required, positional): a shape-(1, N)
%     character array with the name of an epoch of interest (e.g.,
%     'ThetaMaze_AlternativeRunning'). Alternatively, a shape-(L, 1) cell
%     array of character arrays of epochs.
%   thetaPower (char, optional, keyword): a shape-(1, K) character array
%     describing the theta power range of the data interval. The following
%     options are available:
%     'moderate' - intervals with above-average theta power only.
%     'high' - intervals with the theta power above 75th percentile.
%     '' - no restrictions on theta frequency range power (default).
%   minIntervalLength (numeric, optional, keyword): a shape-(1, 1) numeric
%     scalar defining the minimal interval duration (default=0). Smaller
%     intervals will be excluded.
%   onlyTrials (logical, optional, keyword): a shape-(1, 1) logical scalar
%     controlling whether intervals should only be limited to trial periods
%     (default=false).
%   onlyHighSpeed (logical, optional, keyword): a shape-(1, 1) logical
%     scalar controlling whether intervals should only be limited to
%     periods with the running speed of 10 cm/s (default=false).
%   sleepState (char, optional, keyword): a shape-(1, 1) character array
%     determining whether intervals should be limited to particular brain
%     states of interest. The following options are available:
%     'wake' - periods when the animal is awake.
%     'ma' - micro-arousal periods.
%     'nrem' - non-REM sleep periods.
%     'rem' - REM sleep periods.
%     '' - no brain state restrictions (default).
%   excludeNoise (logical, optional, keyword): a shape-(1, 1) logical
%     scalar controlling whether intervals containing LFP noise
%     (saturations) should be excluded (default=false).
%
% Returns:
%   sessionIntervals (numeric): a shape-(J, 2) numeric array of selected
%     intervals of interest with the first column corresponding to the
%     onsets while the second one to the offsets.
%
% Dependencies:
%   petersen-lab/petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  dataFile (1,:) {mustBeA(dataFile,'char')}
  epochsOfInterest {mustBeCharOrListedType(epochsOfInterest,'cell')}
  options.thetaPower (1,:) {mustBeA(options.thetaPower,'char')} = '';
  options.minIntervalLength (1,1) {mustBeNumeric,mustBeNonnegative} = 0;
  options.onlyTrials (1,1) {mustBeA(options.onlyTrials,'logical')} = false
  options.onlyHighSpeed (1,1) {mustBeA(options.onlyHighSpeed,'logical')} = false
  options.sleepState {mustBeMember(options.sleepState, {'wake','ma','nrem','rem',''})} = ''
  options.excludeNoise (1,1) {mustBeA(options.excludeNoise,'logical')} = false;
end

% Parse input
if ischar(epochsOfInterest)
  epochsOfInterest = {epochsOfInterest};
end

% Load data
sessionFile = strrep(dataFile, '*', 'session');
if exist(sessionFile,'file')
  load(sessionFile); %#ok<*LOAD>
else
  sessionIntervals = [];
  warning(['Session file ' sessionFile ' does not exist. Terminating function call.'])
  return
end
if options.excludeNoise
  lfpNoiseFile = strrep(dataFile, '*', 'lfpNoise.events');
  if exist(lfpNoiseFile,'file')
    load(lfpNoiseFile);
  else
    options.excludeNoise = false;
    warning(['LFP noise file ' lfpNoiseFile ' does not exist.'])
  end
end
if strcmpi(options.thetaPower,'moderate') || strcmpi(options.thetaPower,'high')
  thetaPowerFile = strrep(dataFile, '*', 'thetaPower.timeseriesCollection');
  if exist(thetaPowerFile,'file')
    load(thetaPowerFile);
  else
    options.thetaPower = '';
    warning(['Theta power file ' thetaPowerFile ' does not exist.'])
  end
end
if options.onlyTrials || options.onlyHighSpeed
  circularTrackFile = strrep(dataFile, '*', 'circular_track.behavior');
  if exist(circularTrackFile,'file')
    load(circularTrackFile);
  else
    options.onlyTrials = false;
    options.onlyHighSpeed = false;
    warning(['Circular track behaviour file ' circularTrackFile ' does not exist. Behaviour data is missing.'])
  end
end
if ~isempty(options.sleepState)
  eegStatesFile = strrep(dataFile, '*', 'eegstates');
  if exist(eegStatesFile,'file')
    load(eegStatesFile);
  else
    options.sleepState = '';
    warning(['EEG states file ' eegStatesFile ' does not exist. EEG states data is missing.'])
  end
end

% Proceed with interval selection
sessionIntervals = [];
for epoch = 1:numel(session.epochs)
  if isfield(session.epochs{epoch},'behavioralParadigm') && ...
      ismember(session.epochs{epoch}.behavioralParadigm, epochsOfInterest)

    % Select time intervals corresponding to the sessions of interest
    interval = [session.epochs{epoch}.startTime session.epochs{epoch}.stopTime];

    % Exclude intervals with noisy (saturated) LFPs
    if options.excludeNoise && exist('lfpNoise','var')
      cleanIntervals = invertIntervals(lfpNoise.timestamps, ...
        session.epochs{epoch}.stopTime, startTime=session.epochs{epoch}.startTime);
      if ~isempty(cleanIntervals)
        interval = intervalOverlap(interval, cleanIntervals);
      end
    end

    % Select time intervals corresponding to increased theta power
    if (strcmpi(options.thetaPower,'moderate') || strcmpi(options.thetaPower,'high')) && ~isempty(interval)
      if strcmpi(options.thetaPower,'moderate')
        increasedThetaPowerIntervals = thetaPower.data(:,3)' | thetaPower.data(:,4)';
      elseif strcmpi(options.thetaPower,'high')
        increasedThetaPowerIntervals = logical(thetaPower.data(:,4)');
      end
      increasedThetaPowerIntervals = logical2intervals(increasedThetaPowerIntervals);
      increasedThetaPowerIntervals(:,2) = increasedThetaPowerIntervals(:,2) + 1;
      increasedThetaPowerIntervals(end,2) = min([increasedThetaPowerIntervals(end,2) numel(thetaPower.timestamps)]);
      increasedThetaPowerIntervals = thetaPower.timestamps(increasedThetaPowerIntervals);
      interval = intervalOverlap(interval, increasedThetaPowerIntervals);
    end

    % Select time intervals corresponding to trials only
    if options.onlyTrials && ~isempty(interval)
      if exist('circular_track','var') && isfield(circular_track,'trials')
        trialIntervals = [circular_track.trials.alternation.start ...
          circular_track.trials.alternation.end];
      else
        trialIntervals = [];
      end
      interval = intervalOverlap(interval, trialIntervals);
    end

    % Select time intervals corresponding to high speeds
    if options.onlyHighSpeed && ~isempty(interval)
      if exist('circular_track','var') && isfield(circular_track,'speed')
        highSpeedIntervals = circular_track.speed >= circular_track.speed_th;
        highSpeedIntervals = logical2intervals(highSpeedIntervals);
        highSpeedIntervals(:,2) = highSpeedIntervals(:,2) + 1;
        highSpeedIntervals(end,2) = min([highSpeedIntervals(end,2) numel(circular_track.timestamps)]);
        highSpeedIntervals = circular_track.timestamps(highSpeedIntervals);
      else
        highSpeedIntervals = [];
      end
      interval = intervalOverlap(interval, highSpeedIntervals);
    end

    % Select intervals corresponding to a brain state of interest
    if ~isempty(options.sleepState)
      if exist('SleepState','var') && isfield(SleepState,'ints')
        switch options.sleepState
          case 'wake'
            interval = intervalOverlap(interval, SleepState.ints.WAKEstate);
          case 'ma'
            interval = intervalOverlap(interval, SleepState.ints.MAstate);
          case 'nrem'
            interval = intervalOverlap(interval, SleepState.ints.NREMstate);
          case 'rem'
            interval = intervalOverlap(interval, SleepState.ints.REMstate);
        end
      end
    end

    % Concatenate all relevant intervals for the session
    if ~isempty(interval) && sum(interval(:,2) - interval(:,1) >= options.minIntervalLength)
      sessionIntervals = [sessionIntervals; ...
        interval(interval(:,2) - interval(:,1) >= options.minIntervalLength,:)]; %#ok<*AGROW>
    end
  end
end