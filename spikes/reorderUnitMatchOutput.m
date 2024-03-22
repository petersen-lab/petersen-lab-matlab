function reorderUnitMatchOutput(outputDir, channelOrder, channelOrderAllSessions)
% reorderUnitMatchOutput(outputDir, channelOrder, channelOrderAllSessions)
%
% Function reorders channels in UnitMatch output data
%
% Args:
%   outputDir
%   channelOrder
%   channelOrderAllSessions
%
% Returns:
%   None
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  outputDir
  channelOrder
  channelOrderAllSessions
end

% Load PreparedData.mat file
filename = fullfile(outputDir, 'PreparedData.mat');
load(filename); %#ok<*LOAD> 

% Reorder clusinfo
for u = 1:numel(clusinfo.ch)
  clusinfo.depth(u) = Params.AllChannelPos{clusinfo.RecSesID(u)}(clusinfo.ch(u)+1,2);
  clusinfo.ch(u) = channelOrder(clusinfo.ch(u)+1)-1;
end

% Reorder Params
for session = 1:size(channelOrderAllSessions)
  channels = sort(channelOrderAllSessions(session,:));
  reorderedChannels = channels(channelOrderAllSessions(session,:));
  Params.Coordinates{session} = Params.Coordinates{session}(reorderedChannels,:);
end

% Save PreparedData.mat file
save(filename, 'clusinfo','Params','sp', '-v7.3'); %#ok<*USENS,*NODEF> 

% Reorder waveforms
waveformFolder = fullfile(outputDir, 'RawWaveforms');
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