function [oscScore, oscFreq, oscScoreHist, shuffledOscScore, ...
  shuffledOscScoreHist, significanceCutoff, options] = multiOscScore( ...
  spikeTimes, options)
% [oscScore, oscFreq, oscScoreHist, shuffledOscScore, ...
%   shuffledOscScoreHist, significanceCutoff, options] = multiOscScore( ...
%   spikeTimes, options)
%
% Function calculates multiple oscillation scores and related statistics.
%
% Args:
%   spikeTimes (cell | numeric, required, positional): a shape-(M, 1) cell
%     array of numeric spike time vectors or a shape-(1, N) single numeric
%     spike time vector.
%   freqRange (numeric, optional, keyword): a shape-(1, 2) numeric array
%     defining frequency range over which to calculate the oscillation
%     score (default=[4 12]);
%   sampleRate(numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     frequency of sampling rate from the raw data (default=30000).
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
%
% Returns:
%   oscScore (numeric): a shape-(M, 1) numeric array of ocillation scores
%     for each spikeTimes entry (cell).
%   oscFreq (numeric): a shape-(M, 1) numeric array of ocillation
%     frequencies for each spikeTimes entry (cell).
%   oscScoreHist (numeric): a shape-(2, K) numeric array of oscillation
%     score bin values (the first row) and corresponding entry counts (the
%     second row) comprising oscillation score histogram for the input
%     data. Empty entries are ommited from the histogram calculation.
%   shuffledOscScore (numeric): a shape-(M, 1) numeric array of shuffled
%     ocillation scores for each spikeTimes entry (cell).
%   shuffledOscScoreHist (numeric): a shape-(2, L) numeric array of
%     shuffled oscillation score bin values (the first row) and
%     corresponding entry counts (the second row) comprising oscillation
%     score histogram for the shuffled input data. Empty entries are
%     ommited from the histogram calculation.
%   significanceCutoff (numeric): a shape-(1, 1) numeric scalar marking 95%
%     oscillation score confidence interval. Spike times series with
%     oscillation scores above this value are considered to be oscillatory
%     within the defined freqRange.
%   options (struct): a shape-(1, 1) scalar structure containing parameter
%     values used in the oscillation score calculations. The structure
%     fields correspond to the optional input variables: 'freqRange', 'sr',
%     and 'shuffle'. If data shuffling was used to estimate the
%     significance threshold, 'shuffle' field will hold information
%     desribing the type and the seed of the random number generator (type
%     'doc rng' for more info) used in the data shuffling procedure. You
%     have to use the same random number generator parameters if you want
%     to reproduce the same results.
%
% Comments:
%   The shuffling procedure is done by pooling spike times from all entries
%   together and then randomly drawing new spike times from this pool to
%   replace spike times of individual entries. Shuffled entries have the
%   same number of spikes as the original ones. The oscillation scores for
%   the shuffled entries are then calculated and the 95th percentile of
%   the distribution of these shuffle oscillation scores is taken to
%   represent the significance cut-off level.
%
% Dependencies:
%   The Oscillation Score (https://www.raulmuresan.ro/sources/oscore/).
%   petersen-lab/petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).


arguments
  spikeTimes (:,:) {mustBeNumericOrListedType(spikeTimes,'cell')}
  options.freqRange (1,2) {mustBePositive,mustBeVector} = [4 12]
  options.sampleRate (1,1) {mustBeInteger,mustBePositive} = 30000 % Sample rate of the data in Hz
  options.sr (1,1) {mustBePositive} = 500 % Sample rate of the autocorrelogram
  options.shuffle (1,1) {mustBeLogicalOrListedType(options.shuffle,'struct')} = false;
end

% Parse input
if ~iscell(spikeTimes)
  spikeTimes = {spikeTimes};
end

% Calculate oscillation score
nEntries = numel(spikeTimes);
oscScore = zeros(nEntries,1);
oscFreq = zeros(nEntries,1);
for entry = 1:nEntries
  if isempty(spikeTimes{entry})
    oscScore(entry) = 0;
  else
    [oscScore(entry), ~, oscFreq(entry)] = OScoreSpikes( ...
      {round(spikeTimes{entry}*options.sampleRate)}, ...
      round(spikeTimes{entry}(end)*options.sampleRate), ...
      options.freqRange(1)+1, options.freqRange(2), options.sr);
  end
end

% Bin oscillation scores into a histogram
if nEntries > 1
  [oscScoreHist, scoreBins] = hist(oscScore(logical(oscScore)), round(max(oscScore))+1); %#ok<*HIST>
  oscScoreHist = [scoreBins; oscScoreHist];
else
  oscScoreHist = [];
end

% Establish significance cutoff
if (isstruct(options.shuffle) || options.shuffle) && nEntries > 1
  % Initialise the random number generator
  if isstruct(options.shuffle)
    rng(options.shuffle.Seed, options.shuffle.Type);
    options.rngParams = options.shuffle;
  else
      rng('shuffle');
      options.rngParams = rng;
      options.rngParams = rmfield(options.rngParams, 'State');
  end
  nShuffles = 10; % Number of shuffles to perform
  allShuffledScores = [];

  for i = 1:nShuffles
      % Shuffle and calculate scores
      randSpikeTimes = shuffleSpikeTimes(spikeTimes, customRNG=options.rngParams);
      shuffledOscScore = zeros(nEntries,1);
      for entry = 1:nEntries
          if ~isempty(randSpikeTimes{entry})
              shuffledOscScore(entry) = OScoreSpikes( ...
                  {round(randSpikeTimes{entry}*options.sampleRate)}, ...
                  round(randSpikeTimes{entry}(end)*options.sampleRate), ...
                  options.freqRange(1), options.freqRange(2), options.sr);
          end
      end
      allShuffledScores = [allShuffledScores; shuffledOscScore];
  end

  % Calculate significance cutoff from all shuffled scores
  significanceCutoff = prctile(allShuffledScores(allShuffledScores > 0), 95);

  % You might want to keep the last shuffled scores for the histogram
  shuffledOscScore = allShuffledScores((end-nEntries+1):end);
  [shuffledOscScoreHist, scoreBins] = hist(shuffledOscScore(shuffledOscScore > 0), round(max(shuffledOscScore))+1);
  shuffledOscScoreHist = [scoreBins; shuffledOscScoreHist];
end