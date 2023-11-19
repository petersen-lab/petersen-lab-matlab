function saveLFPNoisePeriods(lfpNoisePeriods, options)
% saveLFPNoisePeriods(lfpNoisePeriods, <outputFilename>)
%
% Save LFP recording noise periods in CellExplorer format
% (https://cellexplorer.org/).
%
% Args:
%   lfpNoisePeriods (numeric, required, positional): a shape-(M, 2) numeric
%     array of LFP recording noise periods.
%   outputFilename (char, optional, keyword): a shape-(1, N) character
%     array with a full path filename where the LFP recording noise period
%     data should be saved. The filename should end with
%     '.lfpNoise.events.mat'. If ending is different,
%     '<basefolder-name>.lfpNoise.events.mat' will be appended to the
%     filename. If left empty, the output data will be saved in a file
%     named 'lfpNoise.events.mat' in the current working directory.
%
% Returns:
%   None.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  lfpNoisePeriods (:,2) {mustBeNumeric,mustBeNonnegative}
  options.outputFilename (1,:) {mustBeA(options.outputFilename,'char')} = ''
end

% Parse input
for iPeriod = 1:size(lfpNoisePeriods)
  assert(lfpNoisePeriods(iPeriod,2) > lfpNoisePeriods(iPeriod,1));
  if iPeriod > 1
    assert(lfpNoisePeriods(iPeriod,1) > lfpNoisePeriods(iPeriod-1,1));
  end
end

% Organise data within the events container
nPeriods = size(lfpNoisePeriods, 1);
lfpNoise.timestamps = lfpNoisePeriods;
lfpNoise.eventID = ones(nPeriods,1);
lfpNoise.eventIDlabels = repmat({'noise'}, [nPeriods 1]);
lfpNoise.center = mean(lfpNoisePeriods')'; %#ok<UDIM>
lfpNoise.duration = lfpNoisePeriods(:,2) - lfpNoisePeriods(:,1);
lfpNoise.detectorinfo.detectorname = 'Detected using visualiseSaturations + manual curation; saved using saveLFPNoiseData.';
lfpNoise.detectorinfo.detectiondate = datetime;

% Save the data
if isempty(options.outputFilename)
  options.outputFilename = 'lfpNoise.events.mat';
elseif ~contains(options.outputFilename, 'lfpNoise.events.mat', 'IgnoreCase',true)
  [~, basename] = fileparts(options.outputFilename);
  if contains(options.outputFilename, '.mat')
    options.outputFilename = strrep(options.outputFilename, '.mat', ...
      [basename '.lfpNoise.events.mat']);
  else
    options.outputFilename = fullfile(options.outputFilename, ...
      [basename '.lfpNoise.events.mat']);
  end
end
save(options.outputFilename, 'lfpNoise', '-v7.3');