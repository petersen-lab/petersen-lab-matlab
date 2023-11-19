function saveOscScore(spikesFile, options)
% saveOscScore(spikesFile, <outputFilename>)
%
% Detect and save unit oscillation scores and frequencies in the
% CellExplorer format (https://cellexplorer.org/).
%
% Args:
%   spikesFile (char, required, positional): a shape-(1, M) character
%     array with the filename storing unit spike times.
%   freqRange (numeric, optional, keyword): a shape-(1, 2) numeric array
%     defining frequency range over which to calculate the oscillation
%     score (default=[4 11]);
%   sr (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     frequency for sampling the signal in Hz (default=500).
%   outputFilename (char, optional, keyword): a shape-(1, N) character
%     array with a full path filename where the unit oscillation scores
%     should be saved. The filename should end with
%     '.oscillationScore.cellinfo.mat'. If ending is different,
%     '<basefolder-name>.oscillationScore.cellinfo.mat' will
%     be appended to the filename. If left empty, the output data will be
%     saved in one of the two different ways:
%       (1) If spikesFile follows CellExplorer convention and ends with
%           'spikes.cellinfo.mat', the output file will be saved in the
%           same folder in a file ending with
%           'oscillationScore.cellinfo.mat'.
%       (2) Otherwise, the output is saved in a current working directory
%           in the file named 'oscillationScore.cellinfo.mat'.
%
% Dependencies:
%   The Oscillation Score (https://www.raulmuresan.ro/sources/oscore/).
%   petersen-lab/petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).
%
% Returns:
%   None.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  spikesFile (1,:) {mustBeA(spikesFile,'char')}
  options.freqRange (1,2) {mustBePositive,mustBeVector} = [4 11]
  options.sr (1,1) {mustBePositive} = 500
  options.outputFilename (1,:) {mustBeA(options.outputFilename,'char')} = ''
end

% Parse input
if ~exist(spikesFile, 'file')
  warning('The supplied spikesFile does not exist. Terminating function call.')
  return
end

% Load the unit spiking data
load(spikesFile); %#ok<*LOAD>

% Detect oscillation scores and frequencies of individual units
[oscScore, oscFreq] = multiOscScore(spikes.times, freqRange=options.freqRange, ...
  sr=options.sr);

% Organise oscillation score data within the cellinfo container
oscillationScore.data = [oscScore, oscFreq];
oscillationScore.definition = ['Unit oscillation scores (column 1) and ' ...
  'frequencies (column 2) calculated as in Muresan et al. (2008) ' ...
  'The Oscillation Score: An Efficient Method for Estimating Oscillation ' ...
  'Strength in Neuronal Activity , Journal of Neurophysiology 99: 1333-1353 ' ...
  '(https://www.raulmuresan.ro/sources/oscore/).'];
oscillationScore.processingInfo.params = options;
oscillationScore.processingInfo.function = 'petersen-lab-matlab/spikes/multiOscScore';
oscillationScore.processingInfo.date = datetime;
oscillationScore.processingInfo.username = getenv('username');
oscillationScore.processingInfo.hostname = getenv('computername');

% Save unit oscillation scores
if isempty(options.outputFilename)
  if contains(spikesFile, 'spikes.cellinfo.mat')
    options.outputFilename = strrep(spikesFile, 'spikes', 'oscillationScore');
  else
    options.outputFilename = 'oscillationScore.cellinfo.mat';
  end
elseif ~contains(options.outputFilename, 'oscillationScore.cellinfo.mat', ...
    'IgnoreCase',true)
  [~, basename] = fileparts(options.outputFilename);
  if contains(options.outputFilename, '.mat')
    options.outputFilename = strrep(options.outputFilename, '.mat', ...
      [basename '.oscillationScore.cellinfo.mat']);
  else
    options.outputFilename = fullfile(options.outputFilename, ...
      [basename '.oscillationScore.cellinfo.mat']);
  end
end
save(options.outputFilename, 'oscillationScore', '-v7.3');