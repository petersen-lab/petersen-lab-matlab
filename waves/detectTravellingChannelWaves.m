function travellingThetaWave = detectTravellingChannelWaves(chSpikeTimes, xcoords, ycoords, options)
% travellingWaveData = detectTravellingChannelWaves(chSpikeTimes, xcoords, ycoords, <options>)
%
% Function detects travelling recording channel waves and their direction
% across spiking channels. To detect travelling waves across LFP channels,
% use petersen-lab-matlab/waves/detectTravellingLFPWaves function.
%
% Args:
%   chSpikeTimes (cell, required, positional): a shape-(M, 1) cell array of
%     numeric spike time vectors for individual recording channels on the
%     silicon probe.
%   xcoords (numeric, required, positional): a shape-(1, M) numeric array
%     x-axis coordinates of recording channels.
%   ycoords (numeric, required, positional): a shape-(1, M) numeric array
%     y-axis coordinates of recording channels.
%   samplingInterval (numeric, optional, keyword): a shape-(1, 1) sampling
%     interval for producing convolved spike times series and oscillation
%     score calculations (default=0.002).
%   channelOrder (numeric, optional, keyword): a shape-(1, M) numeric array
%     for re-ordering recording channels. By default, no re-ordering is
%     carried out.
%   channelsOI (numeric, optional, keyword): a shape-(1, N) numeric array
%     for selecting a smaller number of recording channels. By default, all
%     channels are included.
%   freqRange (numeric, optional, keyword): a shape-(1, 2) numeric array
%     defining frequency range over which to calculate the oscillation
%     score (default=[4 11]);
%   axis (char, optional, keyword): a shape-(1, K) character array
%     indicating axes along which to calculate the travelling wave
%     direction. The following three options are available:
%       'vertical' - calculate the direction along the y-axis,
%       'horizontal' - calculate the direction along the x-axis,
%       'all' - calculate direction vertically and horizontally (default).
%   omitnans (logical, optional, keyword): a shape-(1, 1) logical scalar
%     flag to indicate whether empty electrode rows and columns containing
%     only NaNs should be ommited from travelling wave calculations
%     (default=true).
%   firingRateTh (numeric, optional, keyword): a shape-(1, 1) numeric
%     scalar firing rate threshold for excluding recording channels with
%     the firing rate below this value in Hz (default=0)
%   oscillationTh (numeric, optional, keyword): a shape-(1, 1) numeric
%     scalar representing oscillation score threshold for excluding
%     channels with the score that is equal or below this value
%     (default=0). If left empty, shuffling procedure is used to work out a
%     significant oscillation score threshold.
%   electrodePattern (char, optional, keyword): a shape-(1, L) character
%     array representing probe electrode pattern. Currently only two types
%     of patterns are defined:
%       'linear' - describing electrodes with overall increasing only or
%                  decreasing only x and y coordinates (default).
%       'checkerboard' - Neuropixels 1.0 electrode pattern where recording
%                        channels are arranged in a checkerboard pattern.
%
% Returns:
%   travellingThetaWave (struct): a shape(1, 1) scalar structure with the
%     following fields (output of the WaveMonk/GradientMethod function):
%     phaseGradDir (numeric): a shape-(1, J) numeric array with phase
%       gradient directionality (PGD) index values. PGD varries between 0
%       and 1 with high values being indicative of travelling waves. The
%       number of samples in this vector is determined by the
%       samplingInterval input variable.
%     centreWeightedWaveDir (numeric): a shape-(1, J) numeric array with
%       centre-weighted travelling wave direction in radians. Not sure
%       about its meaning (see WaveMonk/GradientMethod).
%     waveSpeed (numeric): a shape-(1, J) numeric array indicating the
%       instantaneous speed of the travelling wave in meters per second.
%     wavelength (numeric): a shape-(1, J) numeric array with instantaneous
%       travelling wavelengths in meters.
%     netWaveDir (numeric): a shape-(1, J) numeric array containing
%       travelling wave direction in radians.
%     timestamps (numeric): a shape-(1, J) numeric array of timestamps
%       corresponding to data samples in the above output vectors.
%
% Dependencies
%   CellExplorer (https://cellexplorer.org/)
%   petersen-lab/petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).
%   The Oscillation Score (https://www.raulmuresan.ro/sources/oscore/).
%   erfanzabeh/WaveMonk (https://github.com/erfanzabeh/WaveMonk).
%   Circular Statistics Toolbox
%     (https://uk.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  chSpikeTimes (:,1) {mustBeA(chSpikeTimes,'cell'),mustBeVector}
  xcoords (1,:) {mustBeNumeric,mustBeVector}
  ycoords (1,:) {mustBeNumeric,mustBeVector}
  options.samplingInterval (1,1) {mustBePositive} = 0.002
  options.channelOrder (1,:) {mustBePositive,mustBeVector} = 1:numel(xcoords)
  options.channelsOI (1,:) {mustBePositive,mustBeVector} = 1:numel(xcoords)
  options.freqRange (1,2) {mustBePositive,mustBeVector} = [4 11]
  options.axis (1,:) {mustBeMember(options.axis,{'all','vertical','horizontal'})} = 'all'
  options.omitnans (1,1) {mustBeLogical} = true
  options.firingRateTh (1,1) {mustBeNumeric,mustBeNonnegative} = 0
  options.oscillationTh (1,1) {mustBeNumeric} = []
  options.electrodePattern (1,:) {mustBeMember(options.electrodePattern,{'linear','checkerboard'})} = 'linear'
end

% Re-order channels
chSpikeTimes = chSpikeTimes(options.channelOrder);

% Threshold channels based on their firing rates
emptyChSpikeTimes = cellfun('isempty', chSpikeTimes);
lastSpikeTime = max(cellfun(@max, chSpikeTimes(~emptyChSpikeTimes)));
chFiringRates = zeros(size(emptyChSpikeTimes));
chFiringRates(~emptyChSpikeTimes) = cellfun(@(x) numel(x)/lastSpikeTime, chSpikeTimes(~emptyChSpikeTimes));
highFiringChannels = find(chFiringRates >= options.firingRateTh);
channelsOI = options.channelsOI(ismember(options.channelsOI, highFiringChannels));
if numel(channelsOI) < 2
  warning('Too few high firing channels to perform travelling wave detection. Terminating function call.');
  travellingThetaWave.phaseGradDir = []; travellingThetaWave.centreWeightedWaveDir = [];
  travellingThetaWave.waveSpeed = []; travellingThetaWave.wavelength = [];
  travellingThetaWave.netWaveDir = []; travellingThetaWave.timestamps = [];
  return
end

% Threshold channels based on their oscillation scores
if isempty(options.oscillationTh)
  [oscScore, ~, ~, ~, ~, options.oscillationTh] = multiOscScore(chSpikeTimes, ...
    freqRange=options.freqRange, sr=round(1/options.samplingInterval), shuffle=true);
else
  oscScore = multiOscScore(chSpikeTimes, freqRange=options.freqRange, ...
    sr=round(1/options.samplingInterval));
end
oscillatingChannels = find(oscScore > options.oscillationTh);
channelsOI = channelsOI(ismember(channelsOI, oscillatingChannels));
if numel(channelsOI) < 2
  warning('Too few oscillating channels to perform travelling wave detection. Terminating function call.');
  travellingThetaWave.phaseGradDir = []; travellingThetaWave.centreWeightedWaveDir = [];
  travellingThetaWave.waveSpeed = []; travellingThetaWave.wavelength = [];
  travellingThetaWave.netWaveDir = []; travellingThetaWave.timestamps = [];
  return
end

% Remove unwanted channels
chSpikeTimes = chSpikeTimes(channelsOI);
nCh = numel(chSpikeTimes);
old2newChMapping = options.channelOrder(channelsOI);

% Get electrode dimensions
% y-axis: yChanInd, yIndInit, yInd, yDim. Each variable has a different meaning and should not be confused.
[~, yChanInd] = unique(ycoords(old2newChMapping),'stable'); % Unique ycoords of relevant channels. Ranges across all channels: nHorzDim*nVertDim. For Npx: [1 384].
[~, ~, yIndInit] = unique(ycoords,'stable'); % All vertical locations. Ranges across all nVertDims. For Npx: [1 192].
[~, ~, yOrd] = unique(yIndInit(old2newChMapping)); % Vertical channel order
uniqueVertLocs = logical([1 diff(yOrd)']); % Index for unique (incremental) vertical locations
yInd = round(cumsum(uniqueVertLocs)); % New vertical channel indices. Ranges across all nVertDims. For Npx: [1 192].
yDim = numel(unique(yInd)); % Number of relevant unique vertical locations. Up to nVertDims. For Npx: Up to 192.
% x-axis: xChanInd, xIndInit, xInd, xDim. Each variable has a different meaning and should not be confused (similar to y-axis).
[~, xChanInd] = unique(xcoords(old2newChMapping),'stable');
uniquexcoords = unique(xcoords);
[~, xIndInit] = ismember(xcoords, uniquexcoords);
xInd = xIndInit(old2newChMapping);
xDim = numel(unique(xIndInit));

% Collapse channel spike times along specific dimensions or NaN values
reshapedChSpikeTimes = cell(yDim,xDim);
for ch = 1:nCh
  reshapedChSpikeTimes(yInd(ch),xInd(ch)) = chSpikeTimes(ch);
end
if strcmpi(options.axis, 'all') && options.omitnans
  % In case of omiting rows and columns containing only NaNs but not
  % collapsing over specified dimensions
  [~, nonEmptyVertInds] = collapseCell(reshapedChSpikeTimes, dim=1);
  [~, nonEmptyHorzInds] = collapseCell(reshapedChSpikeTimes, dim=2);
  nonEmptyInds = ismember(yInd, nonEmptyVertInds) & ...
    ismember(xInd, nonEmptyHorzInds);
  chSpikeTimes = chSpikeTimes(nonEmptyInds);
  [~, yChanInd] = unique(ycoords(nonEmptyInds),'stable');
  yInd = yInd(nonEmptyInds);
  yDim = numel(nonEmptyVertInds);
  xChanInd = xChanInd(nonEmptyHorzInds);
  xDim = numel(nonEmptyHorzInds);
  xInd = xInd(nonEmptyInds);
else
  % In case of omiting rows and columns containing only NaNs and/or
  % collapsing over specified dimensions
  if strcmpi(options.axis, 'vertical')
    [chSpikeTimes, nonEmptyVertInds] = collapseCell( ...
      reshapedChSpikeTimes, dim=1, sortElements='ascend');
    nonEmptyHorzInds = 1;
    if options.omitnans
      yInd = yInd(ismember(yInd, nonEmptyVertInds));
    end
    xInd = ones(size(yInd));
    xChanInd = xInd;
  elseif strcmpi(options.axis, 'horizontal')
    [chSpikeTimes, nonEmptyHorzInds] = collapseCell( ...
      reshapedChSpikeTimes, dim=2, sortElements='ascend');
    nonEmptyVertInds = 1;
    if options.omitnans
      xInd = xInd(ismember(xInd, nonEmptyHorzInds));
    end
    yInd = ones(size(xInd));
    yChanInd = yInd;
  end
  if options.omitnans
    % In case of omiting rows and columns containing only NaNs
    chSpikeTimes = chSpikeTimes(nonEmptyVertInds,nonEmptyHorzInds);
  end
  if strcmpi(options.axis, 'vertical') || ...
      strcmpi(options.axis, 'horizontal') || options.omitnans
    [yDim,xDim] = size(chSpikeTimes);
  end
end

% Calculate hypothetical (because some channels were potentially removed) electrode pitch size
vertDistances = abs(diff(ycoords(yChanInd)));
vertSpacing = median(vertDistances(logical(vertDistances)));
horzDistances = abs(diff(xcoords(xChanInd)));
horzSpacing = median(horzDistances(logical(horzDistances)));
if strcmpi(options.axis, 'vertical')
  pitchSize = vertSpacing/1E6;
elseif strcmpi(options.axis, 'horizontal')
  pitchSize = horzSpacing/1E6;
else
  pitchSize = (horzSpacing + vertSpacing)/2E6; % meters
end

% Convolve spike times with a Gaussian function
[convChSpikeTimes, convTimeBins, convParams] = convolveSpikes( ...
  chSpikeTimes, stepSize=options.samplingInterval, convolutionPoints=50);

% Band-pass filter the convolved population rate at theta frequency range
filtConvChSpikeTimes = bandpassFilterTimeSeries(convChSpikeTimes, ...
  sampleRate=round(1/convParams.stepSize), frequencyRange=[6 11]);

% Calculate theta phase
hilbert1 = hilbert(filtConvChSpikeTimes'); % Apply Hilbert transform
chSpikeTimesPhase = atan2(imag(hilbert1), real(hilbert1))';

% Reshape recording channel spike times matrix
nChunks = ceil(0.05/options.samplingInterval); % So we don't run out of memory
nSamples = numel(convTimeBins);
phaseGradDir = zeros(size(convTimeBins));
centreWeightedWaveDir = zeros(size(convTimeBins));
waveSpeed = zeros(size(convTimeBins));
wavelength = zeros(size(convTimeBins));
netWaveDir = zeros(size(convTimeBins));
for chunk = 1:nChunks
  if chunk == 1
    disp('   Analysing data...')
  end
  disp(['      chunk #' num2str(chunk)]);
  chunkSize = ceil(nSamples/nChunks) + 1;
  inds = max([1 chunkSize*(chunk-1)-1]):min([chunkSize*chunk nSamples]);
  if yDim == 1 % To prevent triggering error in GradientMethod
    reshapedChSpikeTimesPhase = nan(2,xDim,numel(inds));
    for ch = 1:size(chSpikeTimesPhase,1)
      reshapedChSpikeTimesPhase(:, ch,:) = [chSpikeTimesPhase(ch,inds); chSpikeTimesPhase(ch,inds)];
    end
  elseif xDim == 1 % To prevent triggering error in GradientMethod
    reshapedChSpikeTimesPhase = nan(yDim,2,numel(inds));
    for ch = 1:size(chSpikeTimesPhase,1)
      reshapedChSpikeTimesPhase(ch,1,:) = chSpikeTimesPhase(ch,inds);
      reshapedChSpikeTimesPhase(ch,2,:) = chSpikeTimesPhase(ch,inds);
    end
  else
    reshapedChSpikeTimesPhase = nan(yDim,xDim,numel(inds));
    for ch = 1:size(chSpikeTimesPhase,1)
      reshapedChSpikeTimesPhase(yInd(ch), xInd(ch),:) = chSpikeTimesPhase(ch,inds);
    end
  end

  % Detect travelling waves
  [phaseGradDir_chunk, centreWeightedWaveDir_chunk, waveSpeed_chunk, ...
    wavelength_chunk, netWaveDir_chunk] = GradientMethod( ...
    SmoothContinuesPhase(reshapedChSpikeTimesPhase), ...
    reshapedChSpikeTimesPhase, pitchSize, round(1/convParams.stepSize), ...
    options.electrodePattern);
  if chunk == 1
    leftInds = inds;
    rightInds = 1:numel(inds);
  else
    leftInds = inds(2:end);
    rightInds = 2:numel(inds);
  end
  phaseGradDir(leftInds) = phaseGradDir_chunk(rightInds);
  centreWeightedWaveDir(leftInds) = centreWeightedWaveDir_chunk(rightInds);
  waveSpeed(leftInds) = waveSpeed_chunk(rightInds);
  wavelength(leftInds) = wavelength_chunk(rightInds);
  netWaveDir(leftInds) = netWaveDir_chunk(rightInds);

  % Detect travelling waves in shuffled data
  % Write the shuffling function and re-run GradientMethod function with
  % shuffled data to determine PGD cut-off value.
end

% Assign output
travellingThetaWave.phaseGradDir = phaseGradDir;
travellingThetaWave.centreWeightedWaveDir = centreWeightedWaveDir;
travellingThetaWave.waveSpeed = waveSpeed;
travellingThetaWave.wavelength = wavelength;
travellingThetaWave.netWaveDir = netWaveDir;
travellingThetaWave.timestamps = convTimeBins;