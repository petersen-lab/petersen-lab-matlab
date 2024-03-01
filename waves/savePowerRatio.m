function filename = savePowerRatio(freqBandFilename1, freqBandFilename2, options)
% filename = savePowerRatio(freqBandFilename1, freqBandFilename2, <outputFilename>)
%
% Function loads spectrograms of two frequency bands and saves the ratio of
% the two over the time as in powerRatioOverTime =
% power_in_freqBandFilename1/power_in_freqBandFilename2.
%
% Args:
%   freqBandFilename1 (char, required, positional): a shape-(1, N)
%     character array containing the MAT filename where the spectrogram of
%     the first frequency band is stored. The structure of this file should
%     have been generated using either instThetaForPointProcess or
%     freqBandPropertiesForPointProcess functions.
%   freqBandFilename2 (char, required, positional): a shape-(1, M)
%     character array containing the MAT filename where the spectrogram of
%     the second frequency band is stored. The structure of this file
%     should have been generated using either instThetaForPointProcess or
%     freqBandPropertiesForPointProcess functions.
%   outputFilename (char, optional, positional): a shape-(1, K) character
%     array containing the output filename (full path) where the power
%     ratio over time should be saved.
%
% Returns:
%   filename (char): a shape-(1, L or K) character array containing the
%     output filename (full path) where the power ratio over time was
%     saved. If the output file has been supplied, then it corresponds to
%     outputFilename. The inner structure of this file follows the
%     timeseries container structure of the CellExplorer format.
%
% Dependencies:
%   petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).
%   CellExplorer (https://cellexplorer.org/).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  freqBandFilename1 (1,:) {mustBeVector,mustBeA(freqBandFilename1,'char')}
  freqBandFilename2 (1,:) {mustBeVector,mustBeA(freqBandFilename2,'char')}
  options.outputFilename (1,:) {mustBeVector,mustBeA(options.outputFilename,'char')} = 'powerRatio.timeseries.mat'
end

% Parse input
if ~endsWith(freqBandFilename1, 'timeseries.mat')
  error('freqBandFilename1 should be a timeseries container.');
end
if ~endsWith(freqBandFilename2, 'timeseries.mat')
  error('freqBandFilename2 should be a timeseries container.');
end
if ~endsWith(options.outputFilename, 'timeseries.mat')
  error('outputFilename is a timeseries container and therefore should reflect that in its name (should end with timeseries.mat).');
end

% Load data
data1 = load(freqBandFilename1);
fieldName = fieldnames(data1);
data1 = data1.(fieldName{1});
data2 = load(freqBandFilename2);
fieldName = fieldnames(data2);
data2 = data2.(fieldName{1});

% Calculate the power ratio
powerSignal1 = sum(data1.data'); %#ok<*UDIM> 
powerSignal2 = sum(data2.data');
assert(numel(powerSignal1) == numel(powerSignal2), ...
  'The power signals in the two supplied files are of different lengths.');
powerRatio = powerSignal1./powerSignal2;

% Save the power ratio
filename = options.outputFilename;
powerRatioTimeseries.data = powerRatio';
powerRatioTimeseries.timestamps = data1.timestamps;
powerRatioTimeseries.precision = class(powerRatio);
powerRatioTimeseries.units = 'a.u.';
powerRatioTimeseries.nChannels = 1;
powerRatioTimeseries.channelNames = 'powerRatio';
powerRatioTimeseries.sr = round(1/(data1.timestamps(2)-data1.timestamps(1)));
powerRatioTimeseries.nSamples = numel(powerRatio);
powerRatioTimeseries.description = ['The power ratio of two signals. Description of signal 1 is as follows: ' ...
  data1.description ' Description of signal 2 is as follows: ' data2.description];
powerRatioTimeseries.processingInfo.params = [];
powerRatioTimeseries.processingInfo.function = 'petersen-lab-matlab/waves/savePowerRatio';
powerRatioTimeseries.processingInfo.date = datetime;
powerRatioTimeseries.processingInfo.username = getenv('username');
powerRatioTimeseries.processingInfo.hostname = getenv('computername');
save(filename, 'powerRatioTimeseries', '-v7.3');