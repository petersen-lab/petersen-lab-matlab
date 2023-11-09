function [randSpikeTimes, rngParams] = shuffleSpikeTimes(spikeTimes, options)
% [randSpikeTimes, rngParams] = shuffleSpikeTimes(spikeTimes, <customRNG>)
%
% Randomly shuffles spike times between spike times series but preserves
% spike counts of individual spike times series.
%
% Args:
%   spikeTimes (cell, required, positional): a shape-(M, 1) cell array of
%     numeric spike time vectors.
%   customRNG (struct, optional, keyword): a shape-(1, 1) structure scalar
%     specifying the parameters of the random number generator to be used
%     for data shuffling (guarantees reproducibility). The structure should
%     have fields 'Type' and 'Seed' which specify the type of the random
%     number generator and the numeric seed. For more information, check
%     the documentation of the Matlab's rng function.
%
% Returns:
%   randSpikeTimes (cell): a shape-(M, 1) cell array of randomly shuffled
%     numeric spike time vectors.
%   rngParams (struct): a shape-(1, 1) structure scalar specifying the
%     parameters of the random number generator used for data shuffling.
%     The structure contains fields 'Type' and 'Seed' which specify the
%     type of the random number generator used and its numeric seed. For
%     more information, check the documentation of the Matlab's rng
%     function.
%
% Dependencies:
%   petersen-lab/petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  spikeTimes (:,:) {mustBeA(spikeTimes,'cell')}
  options.customRNG (:,:) = []
end

% Parse input
if ~isempty(options.customRNG) && ~isstruct(options.customRNG)
  error('customRNG parameter can only be a structure or an empty array.');
end

% Initialise the random number generator
if isempty(options.customRNG)
  rng('shuffle');
  options.customRNG = rng;
  options.customRNG = rmfield(options.customRNG, 'State');
else
  rng(options.customRNG.Seed, options.customRNG.Type);
end

% Get population rate
for entry = 1:numel(spikeTimes)
  if ~isempty(spikeTimes{entry})
    if size(spikeTimes{entry},2) == 1
      populationRate = sort(concatenateCells(spikeTimes));
    elseif size(spikeTimes{entry},1) == 1
      populationRate = sort(concatenateCells(spikeTimes,2));
    else
      error('Individual cells of spikeTimes contain non-vector arrays.');
    end
    break
  elseif entry == numel(spikeTimes)
    warning('Empty spikeTimes cell array.');
    randSpikeTimes = spikeTimes;
    rngParams = options.customRNG;
    return
  end
end

% Shuffle the data
randSpikeTimes = cell(size(spikeTimes));
for entry = 1:numel(spikeTimes)
  nTotalSpikes = numel(populationRate);
  nSpikes = numel(spikeTimes{entry});
  randSpikeInds = randperm(nTotalSpikes, nSpikes);
  randSpikeTimes{entry} = populationRate(randSpikeInds);
  populationRate = populationRate(~ismember(1:nTotalSpikes, randSpikeInds));
end
rngParams = options.customRNG;