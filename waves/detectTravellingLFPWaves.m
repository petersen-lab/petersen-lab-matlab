function [travellingWave, options] = detectTravellingLFPWaves(lfp, xcoords, ycoords, options)
% [travellingWave, options] = detectTravellingLFPWaves(lfp, xcoords, ycoords, <options>)
%
% Function detects travelling recording channel LFP waves and their
% direction across spiking channels. To detect travelling waves across LFP
% channels, use petersen-lab-matlab/waves/detectTravellingLFPWaves function.
%
% Args:
%   lfp (numeric, required, positional): a shape-(M, N) numeric array of
%     LFP timeseries for individual recording channels on the silicon
%     probe. Rows correspond to individual timeseries.
%   xcoords (numeric, required, positional): a shape-(1, M) numeric array
%     x-axis coordinates of recording channels.
%   ycoords (numeric, required, positional): a shape-(1, M) numeric array
%     y-axis coordinates of recording channels.
%   lfpSamplingInterval (numeric, optional, keyword): a shape-(1, 1) LFP
%     sampling interval (default=0.0004).
%   samplingInterval (numeric, optional, keyword): a shape-(1, 1) sampling
%     interval for downsampling the LFP signal (default=0.002).
%   channelOrder (numeric, optional, keyword): a shape-(1, M) numeric array
%     for re-ordering recording channels. By default, no re-ordering is
%     carried out.
%   channelsOI (numeric, optional, keyword): a shape-(1, N) numeric array
%     for selecting a smaller number of recording channels. By default, all
%     channels are included.
%   freqRange (numeric, optional, keyword): a shape-(1, 2) numeric array
%     defining frequency range over which to calculate the oscillation
%     score and filter data (default=[4 12]);
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
%   oscillationTh (logical, optional, keyword): a shape-(1, 1) logical
%     scalar setting whether to exclude non-oscillatory channels (true) or
%     include all channels (false; default).
%   electrodePattern (char, optional, keyword): a shape-(1, L) character
%     array representing probe electrode pattern. Currently only two types
%     of patterns are defined:
%       'linear' - describing electrodes with overall increasing only or
%                  decreasing only x and y coordinates (default).
%       'checkerboard' - Neuropixels 1.0 electrode pattern where recording
%                        channels are arranged in a checkerboard pattern.
%   pgdTh (numeric | struct, optional, keyword): a shape-(1, 1) numeric
%     scalar setting the phase gradient directionality index significance
%     threshold (default=0). If left empty, the threshold is estimated
%     based on random data shuffling. Alternatively, one can supply a
%     shape-(1, 1) structure scalar specifying the parameters of the random
%     number generator to be used for data shuffling (guarantees
%     reproducibility). The structure should have fields 'Type' and 'Seed'
%     which specify the type of the random number generator and the numeric
%     seed. For more information, check the documentation of the Matlab's
%     rng function.
%   envTh (numeric | struct, optional, keyword): a shape-(1, 1) numeric
%     scalar setting the oscillation envelope significance threshold
%     (default=0). The threshold allows determining oscillatory periods in
%     the LFP data. If left empty, the threshold is estimated based on
%     random data shuffling. Alternatively, one can supply a shape-(1, 1)
%     structure scalar specifying the parameters of the random number
%     generator to be used for data shuffling (guarantees reproducibility).
%     The structure should have fields 'Type' and 'Seed' which specify the
%     type of the random number generator and the numeric seed. For more
%     information, check the documentation of the Matlab's rng function.
%   cycleDemarcationMethod (char, optional, keyword): a shape-(1, F)
%     character array determining the method of oscillation cycle
%     demarcation. The following two methods are available:
%     'medianCh' - a median of oscillatory cycles across all LFP channels
%                  (default).
%     'maxPowerCh' - oscillatory cycles of the most oscillatory LFP
%                    recording channel only.
%   verbose (logical, optional, keyword): a shape-(1, 1) logical scalar
%     controlling whether the function should produce processing reports
%     (default=false);
%
% Returns:
%   travellingWave (struct): a shape(1, 1) scalar structure with the
%     following fields (output of the WaveMonk/GradientMethod function):
%     phaseGradDir (numeric): a shape-(1, J) numeric array with phase
%       gradient directionality (PGD) index values. PGD varries between 0
%       and 1 with high values being indicative of travelling waves. The
%       number of samples in this vector is determined by the
%       samplingInterval input variable.
%     shuffledPhaseGradDir (numeric): a shape-(1, J) numeric array with PGD
%       index values based on shuffled LFP data.
%     centreWeightedWaveDir (numeric): a shape-(1, J) numeric array with
%       centre-weighted travelling wave direction in radians (direction at
%       the centre of the electrode; see WaveMonk/GradientMethod). Positive
%       phase indicates a spread from lower order to higher order channels.
%       For example, if the lower order channels are deeper, then the wave
%       travels in the ventrodorsal direction. Negative phase in this case
%       would indicate a spread in the dorsoventral direction.
%     waveSpeed (numeric): a shape-(1, J) numeric array indicating the
%       instantaneous speed of the travelling wave in meters per second.
%     wavelength (numeric): a shape-(1, J) numeric array with instantaneous
%       travelling wavelengths in meters.
%     netWaveDir (numeric): a shape-(1, J) numeric array containing
%       travelling wave direction in radians. Positive phase indicates a
%       spread from higher order to lower order channels. For example, if
%       the lower order channels are deeper, then the wave travels in the
%       dorsoventral direction. Negative phase in this case would indicate
%       a spread in the ventrodorsal direction.
%     shuffledNetWaveDir (numeric): a shape-(1, J) numeric array containing
%       travelling wave direction in radians calculated based on the
%       shuffled LFP data.
%     lfpOscEnvelope (numeric): a shape-(1, J) numeric array of the average
%       Hilbert transform envelope of the band-pass (freqRange) filtered
%       LFP traces.
%     shuffledLFPOscEnvelope (numeric): a shape-(1, J) numeric array of the
%       average Hilbert transform envelope of the band-pass (freqRange)
%       filtered shuffled LFP traces.
%     travellingWaveLocs (logical): a shape-(1, J) logical array marking
%       locations where the phase gradient directionality index exceeds
%       pgdTh threshold and also the LFP signal is oscillatory as indicated
%       by lfpOscEnvelope exceeding the envTh parameter.
%     travellingWaveCycleLocs (logical): a shape-(1, J) logical array
%       marking full oscillation cycles containing travelling waves as
%       marked by travellingWaveLocs.
%     cycleNumbers (numeric): a shape-(1, J) numeric array of data sample
%       labels in terms of successive oscillation cycles (cycle order
%       number). Cycles are determined based either on the most oscillatory
%       LFP channel or either on the median cycle of all LFP channels
%       (default).
%     timestamps (numeric): a shape-(1, J) numeric array of timestamps
%       corresponding to data samples in the above output vectors.
%     shuffledPGDHist (numeric): a shape-(2, I) numeric array with PGD bins
%       (1st row) and PGD counts (2nd row) based on the shuffled LFP data.
%     shuffledEnvelopeHist (numeric): a shape-(2, H) numeric array with
%       Hilbert transform envelope bins (1st row) and the corresponding
%       counts (2nd row) for shuffled LFP data.
%     netWaveDirHist (numeric): a shape-(2, G) numeric array with estimated
%       electrode channel net spiking phase gradient direction bins (1st
%       row) and the corresponding counts (2nd row). This is the histogram
%       of the netWaveDir vector.
%     travellingWaveDirHist (numeric): a shape-(2, G) numeric array with
%       estimated travelling wave-only direction bins (1st row) and the
%       corresponding counts (2nd row).
%     nontravellingWaveDirHist (numeric): a shape-(2, G) numeric array with
%       bins over estimated electrode channel net LFP phase gradient
%       direction outside detected travelling waves (1st row) and the
%       corresponding counts (2nd row). the sum of travellingWaveDirHist
%       and nontravellingWaveDirHist counts should equal netWaveDirHist
%       count.
%     shuffledNetWaveDirHist (numeric): same as netWaveDirHist but for the
%       shuffled LFP data.
%     shuffledTravellingWaveDirHist (numeric): same as
%       travellingWaveDirHist but for the shuffled LFP data.
%     shuffledNontravellingWaveDirHist (numeric): same as
%       nontravellingWaveDirHist but for the shuffled LFP data.
%     travellingWaveCount (numeric): a shape-(1, 1) numeric scalar of the
%       total count of travelling waves within freqRange.
%     proportionPerOscCycle (numeric): a shape-(1, 1) numeric scalar with
%       the propotion of all oscillation cycles that are travelling waves
%       within freqRange.
%     incidence (numeric): a shape-(1, 1) numeric scalar with the incidence
%       rate of travelling waves (per second).
%   options (struct): a shape-(1, 1) structure scalar containing all
%     optional input parameters. There is also one extra field 'rngs'
%     containing random number generator parameters ('Type' and 'Seed')
%     used to shuffle LFP data to estimate the significance threshold of
%     the PGD index ('pgdTh') as well as to estimate the Hilbert transform
%     envelope significance threshold ('envTh') for the shuffled LFP data.
%     If both thresholds were estimated based on shuffled data, their
%     random number generator parameters should be the same. This info is
%     useful if you want to replicate the travelling wave detection results.
%
% Dependencies
%   CellExplorer (https://cellexplorer.org/)
%   petersen-lab/petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).
%   erfanzabeh/WaveMonk (https://github.com/erfanzabeh/WaveMonk).
%   Circular Statistics Toolbox
%     (https://uk.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  lfp (:,:) {mustBeNumeric}
  xcoords (1,:) {mustBeNumeric,mustBeVector}
  ycoords (1,:) {mustBeNumeric,mustBeVector}
  options.lfpSamplingInterval (1,1) {mustBePositive} = 0.0004
  options.samplingInterval (1,1) {mustBePositive} = 0.002
  options.channelOrder (1,:) {mustBePositive,mustBeVector} = 1:numel(xcoords)
  options.channelsOI (1,:) {mustBePositive,mustBeVector} = 1:numel(xcoords)
  options.freqRange (1,2) {mustBePositive,mustBeVector} = [4 12]
  options.axis (1,:) {mustBeMember(options.axis,{'all','vertical','horizontal'})} = 'all'
  options.omitnans (1,1) {islogical} = true
  options.oscillationTh (1,1) {islogical} = false
  options.electrodePattern (1,:) {mustBeMember(options.electrodePattern,{'linear','checkerboard'})} = 'linear'
  options.pgdTh (:,:) {mustBeNumericOrListedType(options.pgdTh,'struct')} = 0
  options.envTh (:,:) {mustBeNumericOrListedType(options.envTh,'struct')} = 0
  options.cycleDemarcationMethod (1,:) {mustBeMember(options.cycleDemarcationMethod,{'medianCh','maxPowerCh'})} = 'medianCh'
  options.verbose (1,1) {islogical} = false
end

% Parse input
if isstruct(options.pgdTh) && isstruct(options.envTh)
  assert(all(isfield(options.pgdTh,{'Type','Seed'})))
  assert(all(isfield(options.envTh,{'Type','Seed'})))
  assert(strcmpi(options.pgdTh.Type,options.envTh.Type))
  assert(options.pgdTh.Seed == options.envTh.Seed)
end

% Downsample LFP traces
originalTimestamps = (1:size(lfp,2))*options.lfpSamplingInterval;
downsampledTimestamps = options.samplingInterval:options.samplingInterval:originalTimestamps(end);
lfp = interp1(originalTimestamps, single(lfp'), downsampledTimestamps)';

% Re-order channels
lfp = lfp(options.channelOrder,:);

% Threshold channels based on their oscillation scores
if options.oscillationTh
  oscillatingChannels = detectLFPOsc(lfp, options.freqRange, ...
    samplingFreq=1/options.samplingInterval);
  oscillatingChannels = find(oscillatingChannels);
  options.channelsOI = options.channelsOI(ismember(options.channelsOI, ...
    oscillatingChannels));
end
if numel(options.channelsOI) < 2
  warning('Too few oscillating channels to perform travelling wave detection. Terminating function call.');
  travellingWave.phaseGradDir = []; travellingWave.shuffledPhaseGradDir = [];
  travellingWave.centreWeightedWaveDir = []; travellingWave.waveSpeed = [];
  travellingWave.wavelength = []; travellingWave.netWaveDir = [];
  travellingWave.shuffledNetWaveDir = []; travellingWave.lfpOscEnvelope = [];
  travellingWave.shuffledLFPOscEnvelope = []; travellingWave.travellingWaveLocs = [];
  travellingWave.travellingWaveCycleLocs = []; travellingWave.cycleNumbers = [];
  travellingWave.timestamps = []; travellingWave.shuffledPGDHist = [];
  travellingWave.shuffledEnvelopeHist = []; travellingWave.netWaveDirHist = [];
  travellingWave.travellingWaveDirHist = []; travellingWave.nontravellingWaveDirHist = [];
  travellingWave.shuffledNetWaveDirHist = []; travellingWave.shuffledTravellingWaveDirHist = [];
  travellingWave.shuffledNontravellingWaveDirHist = []; travellingWave.travellingWaveCount = [];
  travellingWave.proportionPerOscCycle = []; travellingWave.incidence = [];
  options.rngs.pgdTh = []; options.rngs.envTh = [];
  return
end

% Remove unwanted channels
lfp = lfp(options.channelsOI,:);
nCh = size(lfp,1);
nSamples = size(lfp,2);
old2newChMapping = options.channelOrder(options.channelsOI);

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
reshapedLFP = cell(yDim,xDim);
for ch = 1:nCh
  reshapedLFP(yInd(ch),xInd(ch)) = {lfp(ch,:)};
end
if strcmpi(options.axis, 'all') && options.omitnans
  % In case of omiting rows and columns containing only NaNs but not
  % collapsing over specified dimensions
  [~, nonEmptyVertInds] = collapseCell(reshapedLFP, dim=1);
  [~, nonEmptyHorzInds] = collapseCell(reshapedLFP, dim=2);
  nonEmptyInds = ismember(yInd, nonEmptyVertInds) & ...
    ismember(xInd, nonEmptyHorzInds);
  lfp = lfp(nonEmptyInds,:);
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
    [reshapedLFP, nonEmptyVertInds] = collapseCell( ...
      reshapedLFP, dim=1, sortElements='none');
    lfp = zeros(numel(reshapedLFP),nSamples);
    for iCh = 1:numel(reshapedLFP)
      lfp(iCh,:) = mean(reshape(reshapedLFP{iCh}', ...
        nSamples, numel(reshapedLFP{iCh})/nSamples)'); %#ok<*UDIM> 
    end
    nonEmptyHorzInds = 1; %#ok<*NASGU> 
    if options.omitnans
      lfp = lfp(nonEmptyVertInds,:);
      yInd = yInd(ismember(yInd, nonEmptyVertInds));
    end
    xInd = ones(size(yInd));
    xChanInd = xInd;
    yDim = size(lfp,1);
    xDim = 1;
  elseif strcmpi(options.axis, 'horizontal')
    [reshapedLFP, nonEmptyHorzInds] = collapseCell( ...
      reshapedLFP, dim=2, sortElements='none');
    lfp = zeros(numel(reshapedLFP),nSamples);
    for iCh = 1:numel(reshapedLFP)
      lfp(iCh,:) = mean(reshape(reshapedLFP{iCh}', ...
        nSamples, numel(reshapedLFP{iCh})/nSamples)'); %#ok<*UDIM> 
    end
    nonEmptyVertInds = 1;
    if options.omitnans
      lfp = lfp(nonEmptyHorzInds,:);
      xInd = xInd(ismember(xInd, nonEmptyHorzInds));
    end
    yInd = ones(size(xInd));
    yChanInd = yInd;
    yDim = 1;
    xDim = size(lfp,1);
  end
end

% Generate randomly shuffled LFP recording
randLevel = '2';
if isstruct(options.pgdTh)
  [randLFP, options.rngs.pgdTh] = shuffleLFP(lfp, randLevel=randLevel, ...
    customRNG=options.pgdTh);
  if isempty(options.envTh)
    options.rngs.envTh = options.rngs.pgdTh;
  else
    options.rngs.envTh = [];
  end
elseif isstruct(options.envTh)
  [randLFP, options.rngs.envTh] = shuffleLFP(lfp, randLevel=randLevel, ...
    customRNG=options.envTh);
  if isempty(options.pgdTh)
    options.rngs.pgdTh = options.rngs.envTh;
  else
    options.rngs.pgdTh = [];
  end
elseif isempty(options.pgdTh) || isempty(options.envTh)
  [randLFP, rngParams] = shuffleLFP(lfp, randLevel=randLevel);
  if isempty(options.pgdTh)
    options.rngs.pgdTh = rngParams;
  else
    options.rngs.pgdTh = [];
  end
  if isempty(options.envTh)
    options.rngs.envTh = rngParams;
  else
    options.rngs.envTh = [];
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

% Band-pass filter the LFP recording traces at the frequency range of interest
filtLFP = bandpassFilterTimeSeries(lfp, ...
  sampleRate=round(1/options.samplingInterval), frequencyRange=options.freqRange);
if exist('randLFP','var')
  filtRandLFP = bandpassFilterTimeSeries(randLFP, ...
    sampleRate=round(1/options.samplingInterval), frequencyRange=options.freqRange);
end

% Calculate oscillation phase and envelope
hilbert1 = hilbert(filtLFP'); % Apply Hilbert transform
lfpPhase = atan2(imag(hilbert1), real(hilbert1))';
lfpOscEnvelope = abs(hilbert1)';
if exist('randLFP','var')
  hilbert1 = hilbert(filtRandLFP');
  randLFPPhase = atan2(imag(hilbert1), real(hilbert1))';
  randLFPOscEnvelope = abs(hilbert1)';
  randLFPOscEnvMean = mean(randLFPOscEnvelope);
end

% Estimate Significance threshold for the oscillation envelope
if isempty(options.envTh) || isstruct(options.envTh)
  flatRandLFPOscEnvelope = reshape(randLFPOscEnvelope,numel(randLFPOscEnvelope),1);
  [randEnvelopeHist, envBins] = hist(flatRandLFPOscEnvelope, 99);
  randEnvelopeHist = [envBins; randEnvelopeHist];
  options.envTh = prctile(flatRandLFPOscEnvelope,95);
else
  randEnvelopeHist = [];
end

% Reshape recording channel spike times matrix
nChunks = ceil(0.05/options.samplingInterval); % So we don't run out of memory
phaseGradDir = zeros(1,nSamples);
randPhaseGradDir = zeros(1,nSamples);
centreWeightedWaveDir = zeros(1,nSamples);
waveSpeed = zeros(1,nSamples);
wavelength = zeros(1,nSamples);
netWaveDir = zeros(1,nSamples);
randNetWaveDir = zeros(1,nSamples);
for chunk = 1:nChunks
  if chunk == 1 && options.verbose
    disp('   Analysing data...')
  end
  if options.verbose
    disp(['      chunk #' num2str(chunk)]);
  end
  chunkSize = ceil(nSamples/nChunks) + 1;
  inds = max([1 chunkSize*(chunk-1)-1]):min([chunkSize*chunk nSamples]);
  reshapedLFPPhase = reshape2ElectrodeDimensions(lfpPhase, xDim, yDim, ...
    xInd, yInd, inds);

  % Calculate the phase gradient
  [phaseGradDir_chunk, centreWeightedWaveDir_chunk, waveSpeed_chunk, ...
    wavelength_chunk, netWaveDir_chunk] = GradientMethod( ...
    SmoothContinuesPhase(reshapedLFPPhase), reshapedLFPPhase, pitchSize, ...
    round(1/options.samplingInterval), options.electrodePattern);
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

  % Calculate the shuffled phase gradient
  if isempty(options.pgdTh) || isstruct(options.pgdTh)
    reshapedLFPPhase = reshape2ElectrodeDimensions(randLFPPhase, ...
      xDim, yDim, xInd, yInd, inds);
    [phaseGradDir_chunk, ~, ~, ~, netWaveDir_chunk] = GradientMethod( ...
    SmoothContinuesPhase(reshapedLFPPhase), reshapedLFPPhase, pitchSize, ...
    round(1/options.samplingInterval), options.electrodePattern);
    randPhaseGradDir(leftInds) = phaseGradDir_chunk(rightInds);
    randNetWaveDir(leftInds) = netWaveDir_chunk(rightInds);
  end
end

% Estimate PGD index significance cutoff
if isempty(options.pgdTh) || isstruct(options.pgdTh)
  % Bin shuffled PGD index into a histogram
  [shuffledPGDHist, pgdBins] = hist(randPhaseGradDir, 99);
  shuffledPGDHist = [pgdBins; shuffledPGDHist];
  % Calculate significance cutoff
  pgdSignificanceCutoff = prctile(randPhaseGradDir,95);
else
  shuffledPGDHist = [];
  pgdSignificanceCutoff = options.pgdTh;
end

% Number oscillatory cycles
if strcmpi(options.cycleDemarcationMethod,'medianCh')
  cycleNumbers = round(median(ceil((unwrap(lfpPhase')'+pi)./(2*pi))));
elseif strcmpi(options.cycleDemarcationMethod,'maxPowerCh')
  [oscillatingChannels, ~, oscMaxPower] = detectLFPOsc(lfp, options.freqRange, ...
    samplingFreq=1/options.samplingInterval);
  oscillatingChannels = find(oscillatingChannels);
  [~, maxPowerCh] = max(oscMaxPower(oscillatingChannels));
  maxPowerCh = oscillatingChannels(maxPowerCh);
  cycleNumbers = median(ceil((unwrap(lfpPhase(maxPowerCh,:))'+pi)./(2*pi)));
end
cycleNumbers(cycleNumbers == 0) = 1;

% Detect travelling waves
lfpOscEnvMean = mean(lfpOscEnvelope);
oscCycleLocs = lfpOscEnvMean>options.envTh;
oscCycleCount = numel(unique(cycleNumbers(oscCycleLocs)));
travellingWaveLocs = phaseGradDir>pgdSignificanceCutoff & lfpOscEnvMean>options.envTh;
travellingWaveOscCycles = unique(cycleNumbers(travellingWaveLocs));
travellingWaveCycleLocs = ismember(cycleNumbers, travellingWaveOscCycles);
travellingWaveCount = numel(travellingWaveOscCycles);
proportionPerOscCycle = travellingWaveCount/oscCycleCount;
incidence = travellingWaveCount/downsampledTimestamps(end);

% Calculate travelling/nontravelling wave direction distributions
if strcmpi(options.axis, 'vertical') || strcmpi(options.axis, 'horizontal')
  [netWaveDirHist, netWaveDirBins] = hist(netWaveDir(~isnan(netWaveDir)), ...
    sort(unique((~isnan(netWaveDir))))); %#ok<*HIST> 
  selectNetWaveDir = netWaveDir(travellingWaveLocs);
  [netWaveDirHist_travel, netWaveDirBins_travel] = hist( ...
    selectNetWaveDir(~isnan(selectNetWaveDir)), ...
    sort(unique(selectNetWaveDir(~isnan(selectNetWaveDir)))));
  selectNetWaveDir = netWaveDir(~travellingWaveLocs);
  [netWaveDirHist_nontravel, netWaveDirBins_nontravel] = hist( ...
    selectNetWaveDir(~isnan(selectNetWaveDir)), ...
    sort(unique(selectNetWaveDir(~isnan(selectNetWaveDir)))));
else
  [netWaveDirHist, netWaveDirBins] = hist(netWaveDir, 99);
  [netWaveDirHist_travel, netWaveDirBins_travel] = hist( ...
    netWaveDir(travellingWaveLocs), 99);
  [netWaveDirHist_nontravel, netWaveDirBins_nontravel] = hist( ...
    netWaveDir(~travellingWaveLocs), 99);
end
netWaveDirHist = [netWaveDirBins; netWaveDirHist];
netWaveDirHist_travel = [netWaveDirBins_travel; netWaveDirHist_travel];
netWaveDirHist_nontravel = [netWaveDirBins_nontravel; netWaveDirHist_nontravel];
if isempty(options.pgdTh) || isstruct(options.pgdTh)
  if strcmpi(options.axis, 'vertical') || strcmpi(options.axis, 'horizontal')
    [randNetWaveDirHist, randNetWaveDirBins] = hist( ...
      randNetWaveDir(~isnan(randNetWaveDir)), ...
      sort(unique(randNetWaveDir(~isnan(randNetWaveDir)))));
    selectRandNetWaveDir = randNetWaveDir(travellingWaveLocs);
    [randNetWaveDirHist_travel, randNetWaveDirBins_travel] = hist( ...
      selectRandNetWaveDir(~isnan(selectRandNetWaveDir)), ...
      sort(unique(selectRandNetWaveDir(~isnan(selectRandNetWaveDir)))));
    selectRandNetWaveDir = randNetWaveDir(~travellingWaveLocs);
    [randNetWaveDirHist_nontravel, randNetWaveDirBins_nontravel] = hist( ...
      selectRandNetWaveDir(~isnan(selectRandNetWaveDir)), ...
      sort(unique(selectRandNetWaveDir(~isnan(selectRandNetWaveDir)))));
  else
    [randNetWaveDirHist, randNetWaveDirBins] = hist(randNetWaveDir, 99);
    [randNetWaveDirHist_travel, randNetWaveDirBins_travel] = hist( ...
      randNetWaveDir(travellingWaveLocs), 99);
    [randNetWaveDirHist_nontravel, randNetWaveDirBins_nontravel] = hist( ...
      randNetWaveDir(~travellingWaveLocs), 99);
  end
  randNetWaveDirHist = [randNetWaveDirBins; randNetWaveDirHist];
  randNetWaveDirHist_travel = [randNetWaveDirBins_travel; randNetWaveDirHist_travel];
  randNetWaveDirHist_nontravel = [randNetWaveDirBins_nontravel; ...
    randNetWaveDirHist_nontravel];
  options.pgdTh = pgdSignificanceCutoff;
else
  randNetWaveDirHist = [];
  randNetWaveDirHist_travel = [];
  randNetWaveDirHist_nontravel = [];
end

% Assign output
% Timeseries
travellingWave.phaseGradDir = phaseGradDir;
travellingWave.shuffledPhaseGradDir = randPhaseGradDir;
travellingWave.centreWeightedWaveDir = centreWeightedWaveDir;
travellingWave.waveSpeed = waveSpeed;
travellingWave.wavelength = wavelength;
travellingWave.netWaveDir = netWaveDir;
travellingWave.shuffledNetWaveDir = randNetWaveDir;
travellingWave.lfpOscEnvelope = lfpOscEnvMean;
travellingWave.shuffledLFPOscEnvelope = randLFPOscEnvMean;
travellingWave.travellingWaveLocs = travellingWaveLocs;
travellingWave.travellingWaveCycleLocs = travellingWaveCycleLocs;
travellingWave.cycleNumbers = cycleNumbers;
travellingWave.timestamps = downsampledTimestamps;
% Histograms
travellingWave.shuffledPGDHist = shuffledPGDHist;
travellingWave.shuffledEnvelopeHist = randEnvelopeHist;
travellingWave.netWaveDirHist = netWaveDirHist;
travellingWave.travellingWaveDirHist = netWaveDirHist_travel;
travellingWave.nontravellingWaveDirHist = netWaveDirHist_nontravel;
travellingWave.shuffledNetWaveDirHist = randNetWaveDirHist;
travellingWave.shuffledTravellingWaveDirHist = randNetWaveDirHist_travel;
travellingWave.shuffledNontravellingWaveDirHist = randNetWaveDirHist_nontravel;
% Scalar measures
travellingWave.travellingWaveCount = travellingWaveCount;
travellingWave.proportionPerOscCycle = proportionPerOscCycle;
travellingWave.incidence = incidence;
end



%% Local functions
function electrodeArrangedTimeSeries = reshape2ElectrodeDimensions(verticallyArrangedTimeSeries, xDim, yDim, xInd, yInd, inds)
% electrodeArrangedTimeSeries = reshape2ElectrodeDimensions(verticallyArrangedTimeSeries, xDim, yDim, xInd, yInd, inds)
%
% A helper function to detectTravellingLFPWaves

if yDim == 1 % To prevent triggering error in GradientMethod
  electrodeArrangedTimeSeries = nan(2,xDim,numel(inds));
  for ch = 1:size(verticallyArrangedTimeSeries,1)
    electrodeArrangedTimeSeries(:, ch,:) = [verticallyArrangedTimeSeries(ch,inds); verticallyArrangedTimeSeries(ch,inds)];
  end
elseif xDim == 1 % To prevent triggering error in GradientMethod
  electrodeArrangedTimeSeries = nan(yDim,2,numel(inds));
  for ch = 1:size(verticallyArrangedTimeSeries,1)
    electrodeArrangedTimeSeries(ch,1,:) = verticallyArrangedTimeSeries(ch,inds);
    electrodeArrangedTimeSeries(ch,2,:) = verticallyArrangedTimeSeries(ch,inds);
  end
else
  electrodeArrangedTimeSeries = nan(yDim,xDim,numel(inds));
  for ch = 1:size(verticallyArrangedTimeSeries,1)
    electrodeArrangedTimeSeries(yInd(ch), xInd(ch),:) = verticallyArrangedTimeSeries(ch,inds);
  end
end
end