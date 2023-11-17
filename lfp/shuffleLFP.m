function [randLFP, rngParams] = shuffleLFP(lfp, options)
% [randLFP, rngParams] = shuffleLFP(lfp, <options>)
%
% Randomly shuffles LFP recording traces. First, the algorithm shuffles
% channel positions. Then, every channel trace is cut into a number of
% equal length temporal segments and the segments temporal positions are
% shuffled. The latter shuffling procedure is carried out separately for
% each individual recording channel. The number of temporal segments is the
% same as the number of recording channels.
%
% Args:
%   lfp (numeric, required, positional): a shape-(M, N) numeric array of
%     LFP timeseries for individual recording channels on the silicon
%     probe. Rows correspond to individual timeseries.
%   randLevel (char, optional, keyword): a shape-(1, 1) character scalar
%     setting the level of LFP randomisation. Two levels are available:
%     '1' - shuffling recording channel positions only;
%     '2' - shuffling both channel positions and temporal segments within
%           individual channels (default).
%   customRNG (struct, optional, keyword): a shape-(1, 1) structure scalar
%     specifying the parameters of the random number generator to be used
%     for data shuffling (guarantees reproducibility). The structure should
%     have fields 'Type' and 'Seed' which specify the type of the random
%     number generator and the numeric seed. For more information, check
%     the documentation of the Matlab's rng function.
%
% Returns:
%   randLFP (numeric): a shape-(M, N) numeric array of randomly shuffled
%     LFP traces.
%   rngParams (struct): a shape-(1, 1) structure scalar specifying the
%     parameters of the random number generator used for data shuffling.
%     The structure contains fields 'Type' and 'Seed' which specify the
%     type of the random number generator used and its numeric seed. For
%     more information, check the documentation of the Matlab's rng
%     function.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  lfp (:,:) {mustBeNumeric}
  options.randLevel (1,1) {mustBeMember(options.randLevel,{'1','2'})} = '2'
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

% Shuffle recording channel positions
nCh = size(lfp,1);
shuffledCh = randperm(nCh);
lfp = lfp(shuffledCh,:);

% Shuffle individual channels
if strcmpi(options.randLevel,'2')
  randLFP = zeros(size(lfp));
  nSegments = nCh;
  nSamples = size(lfp,2);
  segmentLength = floor(nSamples/nSegments);
  for iCh = 1:nCh
    shuffledSegments = randperm(nSegments);
    for iSegment = 1:nSegments
      randLFP(iCh, max([1 (iSegment-1)*segmentLength+1]):iSegment*segmentLength) = ...
        lfp(iCh, max([1 (shuffledSegments(iSegment)-1)*segmentLength+1]): ...
        shuffledSegments(iSegment)*segmentLength);
    end
  end
else
  randLFP = lfp;
end
rngParams = options.customRNG;