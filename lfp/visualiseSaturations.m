function saturationDetectionOutput = visualiseSaturations(lfpFile, options)
% saturationDetectionOutput = visualiseSaturations(lfpFile, <options>)
% 
% Function takes in an binary file with LFP data and detects and displays
% LFP amplitude saturations in LFP recording channels of interest. Use the
% output of the function and the figures it produces to decide which
% periods of LFP and spiking data you would like to exclude from your data
% analysis.
%
% Args:
%   lfpFile (char, required, positional): a shape-(1, P) character array
%     with the full path and a basename of the LFP file in the form:
%     <path-to-data-folder>\<basename>.lfp.
%   chanMapFile (char, optional, keyword): a shape-(1, Q) character array
%     with the full path to the channel map matlab file (kilosort format).
%     Supply this file path if chanMap.mat file does not exist in the same
%     folder as the LFP file or is named differently.
%   lfpSamplingInterval (numeric, optional, keyword): a shape-(1, 1) LFP
%     sampling interval. If left empty, the function will try to discern it
%     by loading sessionInfo data from the same folder location as the LFP
%     file. If the sessionInfo file does not exist, the sampling interval
%     is assumed to be 0.0004 seconds by default.
%   channels (numeric, optional, keyword): a shape-(R, 1) numeric array of
%     LFP channels of interest for which signal saturations should be
%     detected. By default, select 10 LFP channels across the entire
%     electrode array.
%
% Returns:
%   saturationDetectionOutput (struct): a shape-(1, 1) structure scalar
%     accumulating the output of detectLFPsaturations function for each LFP
%     channel of interest. For more explanation about individual strucuture
%     fields, type help detectLFPsaturations.
%
% Comments:
%   visualiseSaturations is a wrapper of detectLFPsaturations function
%   making it compatible with the CellExplorer format
%   (https://cellexplorer.org/datastructure/data-structure-and-format/).
%
% Dependencies:
%   dervinism/lfpProcessingFunctions
%     (https://github.com/dervinism/lfpProcessingFunctions).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  lfpFile (1,:) {mustBeA(lfpFile,'char')}
  options.chanMapFile (1,:) {mustBeA(options.chanMapFile,'char')} = ''
  options.lfpSamplingInterval (:,:) {mustBeNumeric} = []
  options.channels (:,:) {mustBeNumeric} = []
end

% Load LFP data
fR = fopen(lfpFile, 'r');
lfp = fread(fR, 'int16=>int16');

% Load the probe channel map file
if isempty(options.chanMapFile)
  options.chanMapFile = fullfile(fileparts(lfpFile), 'chanMap.mat');
end  
load(options.chanMapFile, 'connected');
nCh = sum(connected);

% Reshape LFP data
try
  lfp = reshape(lfp,nCh,[]);
catch
  lfp = reshape(lfp,nCh+1,[]); % In case there is an extra analog channel
end

% Work out the sampling interval
sessionInfoFile = strrep(lfpFile, 'lfp', 'sessionInfo.mat');
if isempty(options.lfpSamplingInterval)
  if exist(sessionInfoFile,'file') && contains(sessionInfoFile, 'sessionInfo.mat')
    load(sessionInfoFile); %#ok<LOAD> 
    options.lfpSamplingInterval = 1/sessionInfo.lfpSampleRate;
  else
    options.lfpSamplingInterval = 0.0004;
  end
end

% Work out which LFP channels need to be evaluated
if isempty(options.channels)
  options.channels = floor(nCh/10):floor(nCh/10):nCh;
end

% Detect LFP saturations
nCh = numel(options.channels);
nSamples = size(lfp,2);
saturationDetectionOutput.LFPsaturations = zeros(nCh,nSamples);
saturationDetectionOutput.time = zeros(nCh,nSamples);
saturationDetectionOutput.nSaturations = zeros(nCh,1);
saturationDetectionOutput.fSaturations = zeros(nCh,1);
saturationDetectionOutput.meanSatDuration = zeros(nCh,1);
saturationDetectionOutput.f = zeros(nCh,2);
for iCh = 1:nCh
  [saturationDetectionOutput.LFPsaturations(iCh,:), ...
    saturationDetectionOutput.time(iCh,:), ...
    saturationDetectionOutput.nSaturations(iCh), ...
    saturationDetectionOutput.fSaturations(iCh), ...
    saturationDetectionOutput.meanSatDuration(iCh), ...
    saturationDetectionOutput.f(iCh,:)] = detectLFPsaturations( ...
    double(lfp(options.channels(iCh),:)), options.lfpSamplingInterval, ...
    'hist2', true);
  for iFig = 1:numel(saturationDetectionOutput.f(iCh,:))
    fH = figure(saturationDetectionOutput.f(iCh,iFig));
    if iFig == 1
      title(['LFP histogram: ch#' num2str(options.channels(iCh))]);
      set(fH, 'Name',['LFP_histogram_ch' num2str(options.channels(iCh))]);
    elseif iFig == 2
      title(['LFP saturations: ch#' num2str(options.channels(iCh))]);
      set(fH, 'Name',['LFP_saturations_ch' num2str(options.channels(iCh))]);
    end
  end
end