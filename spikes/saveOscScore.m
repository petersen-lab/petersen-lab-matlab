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
%     score (default=[4 12]);
%   sr (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     frequency for sampling the signal in Hz (default=500).
%   shuffle (logical | struct, optional, keyword): a shape-(1, 1) logical
%     scalar to indicate whether the 95% significance level cut-off should
%     be calculated based on randomly shuffling the data (default=false).
%     Alternatively, one can supply a shape-(1, 1) structure scalar
%     specifying the parameters of the random number generator to be used
%     for data shuffling (guarantees reproducibility). The structure should
%     have fields 'Type' and 'Seed' which specify the type of the random
%     number generator and the numeric seed. For more information, check
%     the documentation of the Matlab's rng function.
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
% Comments:
%   If shuffle parameter is enabled, the shuffled data will be saved in a
%   separate file with the 'oscillationScore' part replaced by
%   'oscillationScoreShuffled'.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  spikesFile (1,:) {mustBeA(spikesFile,'char')}
  options.freqRange (1,2) {mustBePositive,mustBeVector} = [4 12]
  options.sr (1,1) {mustBePositive} = 500
  options.shuffle (1,1) {mustBeLogicalOrListedType(options.shuffle,'struct')} = false;
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
[oscScore, oscFreq, ~, shuffledOscScore] = multiOscScore(spikes.times, freqRange=options.freqRange, ...
  sr=options.sr, shuffle=options.shuffle);

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

oscillationScoreShuffled = oscillationScore;
oscillationScoreShuffled.data = shuffledOscScore;
oscillationScoreShuffled.definition = ['Shuffled unit oscillation scores ' ...
  'calculated as in Muresan et al. (2008) ' ...
  'The Oscillation Score: An Efficient Method for Estimating Oscillation ' ...
  'Strength in Neuronal Activity , Journal of Neurophysiology 99: 1333-1353 ' ...
  '(https://www.raulmuresan.ro/sources/oscore/).'];

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

options.outputFilename = strrep(options.outputFilename, ...
  'oscillationScore', 'oscillationScoreShuffled');
save(options.outputFilename, 'oscillationScoreShuffled', '-v7.3');