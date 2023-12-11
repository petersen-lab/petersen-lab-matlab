function [travellingLFPWaveCollectionFile, travellingLFPWaveTimeseriesFiles] = saveTravellingThetaLFPWaves(lfpFile, options)
% [travellingLFPWaveCollectionFile, travellingLFPWaveTimeseriesFiles] = saveTravellingThetaLFPWaves(lfpFile, <options>)
%
% Function detects travelling recording channel LFP waves, their direction,
% and saves this data in accordance with the CellExplorer format. It's a
% wrapper for detectTravellingLFPWaves function that enables it to work
% with the CellExplorer format (see detectTravellingLFPWaves for more info).
%
% Args:
%   lfpFile (char, required, positional): a shape-(1, P) character array
%     with the full path and a basename of the LFP file in the form:
%     <path-to-data-folder>\<basename>.lfp.
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
%     score (default=[4 12]);
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
%   travellingLFPWaveCollectionFile (char): a shape-(1, L) character array
%     with the name of the timeseries collection file where all the
%     timeseries output variables of the detectTravellingChannelWaves
%     function are saved.
%   travellingLFPWaveTimeseriesFiles (cell): a shape-(J, 1) cell array of MAT
%     filenames containing individual timeseries vectors from a breakup
%     collection. Individual files follow CellExplorer timeseries format
%     (https://cellexplorer.org/datastructure/data-structure-and-format/#time-series).
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
%   Martynas Dervinis (martynas,dervinis@gmail.com).

arguments
  lfpFile (1,:) {mustBeA(lfpFile,'char')}
  options.lfpSamplingInterval (1,1) {mustBePositive} = 0.0004
  options.samplingInterval (1,1) {mustBePositive} = 0.002
  options.channelOrder (:,:) {mustBeNumeric} = []
  options.channelsOI (:,:) {mustBeNumeric} = []
  options.freqRange (1,2) {mustBePositive,mustBeVector} = [4 12]
  options.axis (1,:) {mustBeMember(options.axis,{'all','vertical','horizontal'})} = 'all'
  options.omitnans (1,1) {islogical} = true
  options.oscillationTh (1,1) {islogical} = false
  options.pgdTh (:,:) {mustBeNumericOrListedType(options.pgdTh,'struct')} = 0
  options.envTh (:,:) {mustBeNumericOrListedType(options.envTh,'struct')} = 0
  options.cycleDemarcationMethod (1,:) {mustBeMember(options.cycleDemarcationMethod,{'medianCh','maxPowerCh'})} = 'medianCh'
  options.verbose (1,1) {islogical} = false
end

% Load LFP data
if ~exist(lfpFile, 'file')
  if options.verbose
    disp('   No data.')
    disp('Done.')
  end
  return
end
if options.verbose
  disp('   Loading data...')
end
fR = fopen(lfpFile, 'r');
lfp = fread(fR, 'int16=>int16');

% Load the probe channel map file
chanMapFile = fullfile(fileparts(lfpFile), 'chanMap.mat');
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

% Reshape LFP data
lfp = reshape(lfp,nCh,[]);

% Detect travelling theta waves
[travellingLFPThetaWave, params] = detectTravellingLFPWaves(lfp, ...
  xcoords, ycoords, lfpSamplingInterval=options.lfpSamplingInterval, ...
  samplingInterval=options.samplingInterval, ...
  channelOrder=options.channelOrder, channelsOI=options.channelsOI, ...
  freqRange=options.freqRange, axis=options.axis, ...
  omitnans=options.omitnans, oscillationTh=options.oscillationTh, ...
  electrodePattern=electrodePattern, pgdTh=options.pgdTh, ...
  envTh=options.envTh, cycleDemarcationMethod=options.cycleDemarcationMethod, ...
  verbose=options.verbose);
if isempty(travellingLFPThetaWave.phaseGradDir)
  if options.verbose
    disp('   No data.')
    disp('Done.')
  end
  return
end

% Save calculations
% Full timeseries collection
if options.verbose
  disp('   Saving data...')
end
travellingLFPThetaWave2save.data = [travellingLFPThetaWave.phaseGradDir; ...
  travellingLFPThetaWave.shuffledPhaseGradDir; ...
  travellingLFPThetaWave.centreWeightedWaveDir; ...
  travellingLFPThetaWave.waveSpeed; travellingLFPThetaWave.wavelength; ...
  travellingLFPThetaWave.netWaveDir; travellingLFPThetaWave.shuffledNetWaveDir; ...
  travellingLFPThetaWave.lfpOscEnvelope; ...
  travellingLFPThetaWave.shuffledLFPOscEnvelope; ...
  travellingLFPThetaWave.travellingWaveLocs; ...
  travellingLFPThetaWave.travellingWaveCycleLocs; ...
  travellingLFPThetaWave.cycleNumbers]';
travellingLFPThetaWave2save.timestamps = travellingLFPThetaWave.timestamps';
travellingLFPThetaWave2save.precision = class(travellingLFPThetaWave.phaseGradDir);
travellingLFPThetaWave2save.units = {'a.u.', 'a.u.', 'rad', 'm/s', 'm', ...
  'rad', 'rad', 'a.u.', 'a.u.', 'bool', 'bool', '#'};
travellingLFPThetaWave2save.nChannels = size(travellingLFPThetaWave2save.data,2);
travellingLFPThetaWave2save.channelNames = {'pgdIndex', 'pgdIndexShuffled', ...
  'centreWaveDir', 'waveSpeed', 'wavelength', 'netWaveDir', 'netWaveDirShuffled', ...
  'oscillationEnvelope', 'oscillationEnvelopeShuffled', 'travellingWaveLocs', ...
  'travellingWaveCycleLocs', 'cycleNumbers'};
travellingLFPThetaWave2save.sr = round(1/params.samplingInterval);
travellingLFPThetaWave2save.nSamples = numel(travellingLFPThetaWave.phaseGradDir);
travellingLFPThetaWave2save.description = [ ...
  'The 1st data column contains LFP travelling theta (4-11 Hz) wave phase gradient directionality index based on Zabeh et al., Nat Commun. ' ...
  'Rows correspond to timestamps which are downsampled compared to the original data. ' ...
  'The 2nd data column is the same as the 1st one but for shuffled data. ' ...
  'The 3rd data column contains wave direction at the centre of the electrode (Time_Dire output variable). ' ...
  'The 4th data column contains the wave travelling speed (Speed). ' ...
  'The 5th data column corresponds to the length of the theta wave (Wavelength_Grad). ' ...
  'The 6th data column corresponds to overal (net) direction of the traveling wave (NetGradDir). ' ...
  'The 7th data column is the same as the 6th but for shuffled data. ' ...
  'The 8th data column is the average Hilbert transform oscillation envelope of all LFP recording channels. ' ...
  'The 9th data column is the same as the 8th one but for the shuffled LFP data. ' ...
  'The 10th data column marks detected phase gradient locations. ' ...
  'The 11th data column marks theta travelling waves in full cycles. ' ...
  'The 12th data column marks theta cycle order number as detected by the Hilbert transform.'];
travellingLFPThetaWave2save.processingInfo.params = params;
travellingLFPThetaWave2save.processingInfo.function = [ ...
  'petersen-lab/petersen-lab-matlab/waves/detectTravellingLFPWaves (https://github.com/petersen-lab/petersen-lab-matlab); ' ...
  'erfanzabeh/WaveMonk/GradientMethod (https://github.com/erfanzabeh/WaveMonk)'];
travellingLFPThetaWave2save.processingInfo.date = datetime;
travellingLFPThetaWave2save.processingInfo.username = getenv('username');
travellingLFPThetaWave2save.processingInfo.hostname = getenv('computername');
travellingLFPThetaWave = travellingLFPThetaWave2save;
travellingLFPWaveCollectionFile = strrep(lfpFile, 'lfp', 'travellingLFPThetaWave.timeseriesCollection.mat');
save(travellingLFPWaveCollectionFile, 'travellingLFPThetaWave', '-v7.3');
% Individual timeseries
travellingLFPWaveTimeseriesFiles = breakupTimeseriesCollection( ...
  travellingLFPThetaWave, basename=strrep(travellingLFPWaveCollectionFile, ...
  '.travellingLFPThetaWave.timeseriesCollection.mat', ''));
if options.verbose
  disp('Done.')
end