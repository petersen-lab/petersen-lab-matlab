function [travellingSpikingWaveCollectionFile, travellingSpikingWaveTimeseriesFiles] = saveTravellingThetaChannelWaves(dataFile, options)
% [travellingSpikingWaveCollectionFile, travellingSpikingWaveTimeseriesFiles] = saveTravellingThetaChannelWaves(dataFile, options)
%
% Function detects travelling recording channel spiking waves, their
% direction, and saves this data in accordance with the CellExplorer
% format. It's a wrapper for detectTravellingChannelWaves function that
% enables it to work with the CellExplorer format (see
% detectTravellingChannelWaves for more info).
%
% Args:
%   dataFile (char, required, positional): a shape-(P, 1) character array
%     with the full path and a basename of the data file in the form:
%     <path-to-data-folder>\<basename>.*.mat.
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
%   firingRateTh (numeric, optional, keyword): a shape-(1, 1) numeric
%     scalar firing rate threshold for excluding recording channels with
%     the firing rate below this value in Hz (default=0)
%   oscillationTh (numeric | struct, optional, keyword): a shape-(1, 1)
%     numeric scalar representing oscillation score threshold for excluding
%     channels with the score that is equal or below this value
%     (default=0). If left empty, shuffling procedure is used to work out a
%     significant oscillation score threshold. Alternatively, one can
%     supply a shape-(1, 1) structure scalar specifying the parameters of
%     the random number generator to be used for data shuffling (guarantees
%     reproducibility). The structure should have fields 'Type' and 'Seed'
%     which specify the type of the random number generator and the numeric
%     seed. For more information, check the documentation of the Matlab's
%     rng function.
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
%     the population firing rate. If left empty, the threshold is estimated
%     based on randomly generated population firing rate. Alternatively,
%     one can supply a shape-(1, 1) structure scalar specifying the
%     parameters of the random number generator to be used for population
%     firing rate generation (guarantees reproducibility). The structure
%     should have fields 'Type' and 'Seed' which specify the type of the
%     random number generator and the numeric seed. For more information,
%     check the documentation of the Matlab's rng function.
%   verbose (logical, optional, keyword): a shape-(1, 1) logical scalar
%     controlling whether the function should produce processing reports
%     (default=false);
%
% Returns:
%   travellingSpikingWaveCollectionFile (char): a shape-(1, L) character array
%     with the name of the timeseries collection file where all the
%     timeseries output variables of the detectTravellingChannelWaves
%     function are saved.
%   travellingSpikingWaveTimeseriesFiles (cell): a shape-(J, 1) cell array of MAT
%     filenames containing individual timeseries vectors from a breakup
%     collection. Individual files follow CellExplorer timeseries format
%     (https://cellexplorer.org/datastructure/data-structure-and-format/#time-series).
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
%   Martynas Dervinis (martynas,dervinis@gmail.com).

arguments
  dataFile (1,:) {mustBeA(dataFile,'char')}
  options.samplingInterval (1,1) {mustBePositive} = 0.002
  options.channelOrder (:,:) {mustBeNumeric} = []
  options.channelsOI (:,:) {mustBeNumeric} = []
  options.freqRange (1,2) {mustBePositive,mustBeVector} = [4 12]
  options.axis (1,:) {mustBeMember(options.axis,{'all','vertical','horizontal'})} = 'all'
  options.omitnans (1,1) {islogical} = true
  options.firingRateTh (1,1) {mustBeNumeric,mustBeNonnegative} = 0
  options.oscillationTh (:,:) {mustBeNumericOrListedType(options.oscillationTh,'struct')} = 0
  options.pgdTh (:,:) {mustBeNumericOrListedType(options.pgdTh,'struct')} = 0
  options.envTh (:,:) {mustBeNumericOrListedType(options.envTh,'struct')} = 0
  options.verbose (1,1) {islogical} = false
end

% Load spiking data
spikesFile = strrep(dataFile, '*', 'spikes.cellinfo');
if ~exist(spikesFile, 'file')
  if options.verbose
    disp('   No data.')
    disp('Done.')
  end
  return
end
if options.verbose
  disp('   Loading data...')
end
load(spikesFile); %#ok<*LOAD>
spikesFile = strrep(dataFile, '*', 'muaSpikes.cellinfo');
load(spikesFile);

% Load the probe channel map file
chanMapFile = fullfile(fileparts(dataFile), 'chanMap.mat');
load(chanMapFile, 'connected','xcoords','ycoords');
nCh = sum(connected);
if all(diff(xcoords) >= 0) % Test if the electrode pattern is linear
  electrodePattern = 'linear';
else
  electrodePattern = 'checkerboard'; % this has to be re-written in the future
end

% Work out channel order and determine channels of interest
if isempty(options.channelOrder)
  options.channelOrder = 1:numel(xcoords);
end
if isempty(options.channelsOI)
  options.channelsOI = 1:numel(xcoords);
end

% Get channel spike times
if options.verbose
  disp('   Preprocessing data...')
end
chSpikeTimes = getChannelSpikeTimes([spikes.times'; muaSpikes.times'], ...
  [spikes.maxWaveformCh1'; muaSpikes.maxWaveformCh1'], nCh);

% Detect travelling theta waves
[travellingSpikingThetaWave, params] = detectTravellingChannelWaves( ...
  chSpikeTimes, xcoords, ycoords, samplingInterval=options.samplingInterval, ...
  channelOrder=options.channelOrder, channelsOI=options.channelsOI, ...
  freqRange=options.freqRange, axis=options.axis, ...
  omitnans=options.omitnans, firingRateTh=options.firingRateTh, ...
  oscillationTh=options.oscillationTh, ...
  electrodePattern=electrodePattern, pgdTh=options.pgdTh, ...
  envTh=options.envTh, verbose=options.verbose);
if isempty(travellingSpikingThetaWave.phaseGradDir)
  if options.verbose
    disp('   Too few oscillating/strongly firing channels to perform phase gradient detection. Terminating function call.')
    disp('Done.')
  end
  return
end

% Save calculations
% Full timeseries collection
if options.verbose
  disp('   Saving data...')
end
travellingSpikingThetaWave2save.data = [travellingSpikingThetaWave.phaseGradDir; ...
  travellingSpikingThetaWave.shuffledPhaseGradDir; ...
  travellingSpikingThetaWave.centreWeightedWaveDir; ...
  travellingSpikingThetaWave.waveSpeed; travellingSpikingThetaWave.wavelength; ...
  travellingSpikingThetaWave.netWaveDir; travellingSpikingThetaWave.shuffledNetWaveDir; ...
  travellingSpikingThetaWave.populationRateOscEnvelope; ...
  travellingSpikingThetaWave.randPopulationRateOscEnvelope; ...
  travellingSpikingThetaWave.travellingWaveLocs; ...
  travellingSpikingThetaWave.travellingWaveCycleLocs; ...
  travellingSpikingThetaWave.cycleNumbers]';
travellingSpikingThetaWave2save.timestamps = travellingSpikingThetaWave.timestamps';
travellingSpikingThetaWave2save.precision = class(travellingSpikingThetaWave.phaseGradDir);
travellingSpikingThetaWave2save.units = {'a.u.', 'a.u.', 'rad', 'm/s', 'm', ...
  'rad', 'rad', 'a.u.', 'a.u.', 'bool', 'bool', '#'};
travellingSpikingThetaWave2save.nChannels = size(travellingSpikingThetaWave2save.data,2);
travellingSpikingThetaWave2save.channelNames = {'pgdIndex', 'pgdIndexShuffled', ...
  'centreWaveDir', 'waveSpeed', 'wavelength', 'netWaveDir', 'netWaveDirShuffled', ...
  'oscillationEnvelope', 'oscillationEnvelopeRandom', 'travellingWaveLocs', ...
  'travellingWaveCycleLocs', 'cycleNumbers'};
travellingSpikingThetaWave2save.sr = round(1/params.samplingInterval);
travellingSpikingThetaWave2save.nSamples = numel(travellingSpikingThetaWave.phaseGradDir);
travellingSpikingThetaWave2save.description = [ ...
  'The 1st data column contains spiking travelling theta (4-11 Hz) wave phase gradient directionality index based on Zabeh et al., Nat Commun. ' ...
  'Rows correspond to timestamps which are downsampled compared to the original data. ' ...
  'The 2nd data column is the same as the 1st one but for shuffled data. ' ...
  'The 3rd data column contains wave direction at the centre of the electrode (Time_Dire output variable). ' ...
  'The 4th data column contains the wave travelling speed (Speed). ' ...
  'The 5th data column corresponds to the length of the theta wave (Wavelength_Grad). ' ...
  'The 6th data column corresponds to overal (net) direction of the traveling wave (NetGradDir). ' ...
  'The 7th data column is the same as the 6th but for shuffled data. ' ...
  'The 8th data column is the population firing rate Hilbert transform oscillation envelope. ' ...
  'The 9th data column is the same as the 8th one but for a random population firing rate. ' ...
  'The 10th data column marks detected phase gradient locations. ' ...
  'The 11th data column marks theta travelling waves in full cycles. ' ...
  'The 12th data column marks theta cycle order number as detected by the Hilbert transform.'];
travellingSpikingThetaWave2save.processingInfo.params = params;
travellingSpikingThetaWave2save.processingInfo.function = [ ...
  'petersen-lab/petersen-lab-matlab/waves/detectTravellingChannelWaves (https://github.com/petersen-lab/petersen-lab-matlab); ' ...
  'erfanzabeh/WaveMonk/GradientMethod (https://github.com/erfanzabeh/WaveMonk)'];
travellingSpikingThetaWave2save.processingInfo.date = datetime;
travellingSpikingThetaWave2save.processingInfo.username = getenv('username');
travellingSpikingThetaWave2save.processingInfo.hostname = getenv('computername');
travellingSpikingThetaWave = travellingSpikingThetaWave2save;
travellingSpikingWaveCollectionFile = strrep(dataFile, '*', 'travellingSpikingThetaWave.timeseriesCollection');
save(travellingSpikingWaveCollectionFile, 'travellingSpikingThetaWave', '-v7.3');
% Individual timeseries
travellingSpikingWaveTimeseriesFiles = breakupTimeseriesCollection( ...
  travellingSpikingThetaWave, basename=strrep(travellingSpikingWaveCollectionFile, ...
  '.travellingSpikingThetaWave.timeseriesCollection.mat', ''));
if options.verbose
  disp('Done.')
end