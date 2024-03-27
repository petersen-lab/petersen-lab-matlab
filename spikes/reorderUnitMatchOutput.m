function reorderUnitMatchOutput(outputDir, channelOrder, channelOrderAllSessions, options)
% reorderUnitMatchOutput(outputDir, channelOrder, channelOrderAllSessions, <options>)
%
% Function reorders channels in UnitMatch output data.
%
% Args:
%   outputDir (char, required, positional): a shape-(1, N) character array
%     with a full path to the UnitMatch output folder where
%     PreparedData.mat file is stored and where there is also a folder with
%     extracted waveforms called RawWaveforms.
%   channelOrder (numeric, required, positional): a shape-(1, M) numeric
%     array with correct channel order (where the current channels should
%     be moved): e.g., channels1to384(channelOrder).
%   channelOrderAllSessions (numeric, required, positional): a shape-(L, M)
%     numeric array with correct channel orders for all recording sessions.
%     Rows correspond to consequtive recording sessions and columns
%     correspond to the recording channels on the probe.
%   reorderClusinfo (logical, optional, keyword): a shape-(1, 1) logical
%     scalar with true value setting the channel info in clusinfo variable
%     to be reordered (default).
%   reorderWaveforms (logical, optional, keyword): a shape-(1, 1) logical
%     scalar with true value setting channels of extracted waveforms in
%     RawWaveforms folder to be reordered (default=false).
%   verbose (logical, optional, keyword): a shape-(1, 1) logical scalar to
%     set the printing of function execution stages (default=false).
%
% Returns:
%   None
%
% Dependencies:
%   kwikteam/npy-matlab (https://github.com/kwikteam/npy-matlab).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  outputDir (1,:) {mustBeA(outputDir,'char'),mustBeVector}
  channelOrder (1,:) {mustBePositive,mustBeVector}
  channelOrderAllSessions (:,:) {mustBePositive}
  options.reorderClusinfo (1,1) {mustBeA(options.reorderClusinfo,'logical')} = true
  options.reorderWaveforms (1,1) {mustBeA(options.reorderWaveforms,'logical')} = false
  options.verbose (1,1) {mustBeA(options.verbose,'logical')} = false
end

% Load PreparedData.mat file
filename = fullfile(outputDir, 'PreparedData.mat');
if options.verbose
  disp(['Reordering channel info in ' filename]);
end
load(filename); %#ok<*LOAD> 

% Reorder clusinfo
if options.reorderClusinfo
  for u = 1:numel(clusinfo.ch)
    clusinfo.depth(u) = Params.AllChannelPos{clusinfo.RecSesID(u)}(clusinfo.ch(u)+1,2);
    clusinfo.ch(u) = channelOrder(clusinfo.ch(u)+1)-1;
  end
end

% Reorder Params
for session = 1:size(channelOrderAllSessions,1)
  Params.AllChannelPos{session} = ...
    Params.AllChannelPos{session}(channelOrderAllSessions(session,:),:);
  if isfield(Params, 'Coordinates')
    Params.Coordinates{session} = ...
      Params.Coordinates{session}(channelOrderAllSessions(session,:),:);
  end
end

% Save PreparedData.mat file
save(filename, 'clusinfo','Params','sp', '-v7.3'); %#ok<*USENS,*NODEF> 

% Reorder waveforms
if options.reorderWaveforms
  waveformFolder = fullfile(outputDir, 'RawWaveforms');
  if options.verbose
    disp(['Reordering waveforms in ' waveformFolder]);
  end
  waveforms = dir(waveformFolder);
  channels = sort(channelOrder);
  reorderedChannels = channels(channelOrder);
  for iWave = 1:numel(waveforms)
    if endsWith(waveforms(iWave).name, '.npy')
      waveformFilename = fullfile(waveforms(iWave).folder, waveforms(iWave).name);
      data = readNPY(waveformFilename);
      data = data(:,reorderedChannels,:);
      writeNPY(data, waveformFilename);
    end
  end
end