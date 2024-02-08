function savedFiles = saveGeneralisedPhase(populationRateFolderOrFile, options)
% savedFiles = saveGeneralisedPhase(populationRateFolderOrFile, <options>)
%
% Function estimates and saves wide-band frequency oscillation/fluctuation
% phase given the population firing rate file or folder.
%
% Args:
%   populationRateFolderOrFile (char, required, positional): a shape-(1, N)
%     character array with the population firing rate filename (must end
%     with populationRate.cellinfo.mat) or a folder name where the file is.
%     In either case full path has to be supplied.
%   freqRange (numeric, optional, keyword): a shape-(1, 2) vector array
%     with frequency cutoff values for evaluating the wide-band phase of
%     the signal (default=[2 12]).
%   passbandFilter (logical, optional, keyword): a shape-(1, 1) logical
%     scalar controlling whether to bandpass filter the signal within the
%     freqRange (default=true).
%   stepsize (numeric, optional, keyword): a shape-(1, 1) scalar
%     representing the step size (in seconds) used for binning spikes to
%     produce a continuous spiking rate trace (default=0.002).
%   convolutionPoints (numeric, optional, keyword): a shape-(1, 1) scalar
%     representing points of Gaussian convolution (gausswin) used to smooth
%     binned spiking rate (default=25 sample points).
%   showPhase (logical, optional, keyword): a shape-(1, 1) logical scalar
%     for ploting the wide-band phase of the oscillation/fluctuation
%     (default=false). As part of displaying the pahse, the convolved
%     spiking rate and the filtered spiking rate are also displayed.
%   smoothFreq (logical, optional, keyword): a shape-(1, 1) logical scalar
%     for smoothing instantaneous frequency estimates (default = false).
%
% Returns:
%   savedFiles (cell): a shape-(3, 1) cell array of filename character
%     arrays. The first file contains instantenous wideband frequency
%     info, the second one contains instantaneous wideband phase info,
%     while the third one stores wideband oscillation/fluctuation amplitude.
%
% Dependencies:
%   CellExplorer (https://cellexplorer.org/).
%   petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).
%   mullerlab/generalized-phase
%     (https://github.com/mullerlab/generalized-phase).
%   mullerlab/wave-matlab (https://github.com/mullerlab/wave-matlab).
%   Circular Statistics Toolbox
%     (https://github.com/circstat/circstat-matlab).
%   smoothn function on Matlab File Exchange
%     (https://se.mathworks.com/matlabcentral/fileexchange/25634-smoothn).
%
% Author:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  populationRateFolderOrFile (1,:) {mustBeA(populationRateFolderOrFile,'char')}
  options.freqRange (1,2) {mustBeVector,mustBePositive} = [2 12]
  options.passbandFilter (1,1) {islogical} = true
  options.stepsize (1,1) {mustBePositive} = 0.002
  options.convolutionPoints (1,1) {mustBeNumeric,mustBePositive} = 25
  options.showPhase (1,1) {islogical} = false
  options.smoothFreq (1,1) {islogical} = false
end

% Load spiking file
if endsWith(populationRateFolderOrFile, 'populationRate.cellinfo.mat')
  load(populationRateFolderOrFile); %#ok<*LOAD>
  [file2find.folder, file2find.name] = fileparts(populationRateFolderOrFile);
  file2find.name = [file2find.name '.mat'];
elseif exist(populationRateFolderOrFile, 'dir')
  file2find = dir(fullfile(populationRateFolderOrFile, '*populationRate.cellinfo.mat'));
  load(fullfile(file2find.folder, file2find.name));
else
  error(['The supplied filename is not recognised as a legitimate population firing rate file.' ...
         'The filename must end with populationRate.cellinfo.mat.'])
end
basename = strrep(file2find.name, '.populationRate.cellinfo.mat', '');

% Estimate wide-band frequency oscillation/fluctuation phase and instantaneous frequency
[generalisedPhase, parameters] = generalisedPhaseForPointProcess(populationRate.times{1}, ...
  freqRange=options.freqRange, passbandFilter=options.passbandFilter, ...
  stepsize=options.stepsize, convolutionPoints=options.convolutionPoints, ...
  showPhase=options.showPhase, smoothFreq=options.smoothFreq);

% Instantaneous wide-band frequency
widebandInstantFrequency.data = generalisedPhase.instantFreqPerCycle';
widebandInstantFrequency.timestamps = generalisedPhase.cycleEndpoints';
widebandInstantFrequency.precision = class(generalisedPhase.instantFreqPerCycle);
widebandInstantFrequency.units = 'Hz';
widebandInstantFrequency.nChannels = 1;
widebandInstantFrequency.sr = 1/(generalisedPhase.cycleEndpoints(end)/numel(generalisedPhase.cycleEndpoints));
widebandInstantFrequency.nSamples = numel(generalisedPhase.cycleEndpoints);
widebandInstantFrequency.description = ['A vector of instantaneous wide-band oscillation/fluctuation frequencies. The number of samples ' ...
                                        'corresponds to the number of oscillation/fluctuation cycles in the population firing rate. ' ...
                                        'Wideband frequency range: ' num2str(options.freqRange(1)) ' to ' num2str(options.freqRange(2)) ' Hz.'];
widebandInstantFrequency.processingInfo.params = parameters;
widebandInstantFrequency.processingInfo.function = 'petersen-lab-matlab/waves/generalisedPhaseForPointProcess';
widebandInstantFrequency.processingInfo.date = datetime;
widebandInstantFrequency.processingInfo.username = getenv('username');
widebandInstantFrequency.processingInfo.hostname = getenv('computername');
filenameFreq = fullfile(file2find.folder, [basename '.widebandInstantFrequency.timeseries.mat']);
save(filenameFreq, 'widebandInstantFrequency', '-v7.3');

% Instantaneous wide-band phase
widebandInstantPhase.data = generalisedPhase.phase';
widebandInstantPhase.timestamps = generalisedPhase.timestamps';
widebandInstantPhase.precision = class(generalisedPhase.phase);
widebandInstantPhase.units = 'rad';
widebandInstantPhase.nChannels = 1;
widebandInstantPhase.sr = round(1/(generalisedPhase.timestamps(2)-generalisedPhase.timestamps(1)));
widebandInstantPhase.nSamples = numel(generalisedPhase.timestamps);
widebandInstantPhase.description = ['A vector of instantaneous wide-band oscillation/fluctuation phases. The number of samples ' ...
                                    'corresponds to the number of samples in the convolved population firing rate. ' ...
                                    'Wideband frequency range: ' num2str(options.freqRange(1)) ' to ' num2str(options.freqRange(2)) ' Hz.'];
widebandInstantPhase.processingInfo.params = parameters;
widebandInstantPhase.processingInfo.function = 'petersen-lab-matlab/waves/generalisedPhaseForPointProcess';
widebandInstantPhase.processingInfo.date = datetime;
widebandInstantPhase.processingInfo.username = getenv('username');
widebandInstantPhase.processingInfo.hostname = getenv('computername');
filenamePhase = fullfile(file2find.folder, [basename '.widebandInstantPhase.timeseries.mat']);
save(filenamePhase, 'widebandInstantPhase', '-v7.3');

% Wide-band oscillation/fluctuation amplitude
widebandAmplitude.data = generalisedPhase.amplitude';
widebandAmplitude.timestamps = generalisedPhase.timestamps';
widebandAmplitude.precision = class(generalisedPhase.amplitude);
widebandAmplitude.units = 'a.u.';
widebandAmplitude.nChannels = 1;
widebandAmplitude.sr = round(1/(generalisedPhase.timestamps(2)-generalisedPhase.timestamps(1)));
widebandAmplitude.nSamples = numel(generalisedPhase.timestamps);
widebandAmplitude.description = ['A vector of wide-band oscillation/fluctuation amplitude. The number of samples ' ...
                                    'corresponds to the number of samples in the convolved population firing rate. ' ...
                                    'Wideband frequency range: ' num2str(options.freqRange(1)) ' to ' num2str(options.freqRange(2)) ' Hz.'];
widebandAmplitude.processingInfo.params = parameters;
widebandAmplitude.processingInfo.function = 'petersen-lab-matlab/waves/generalisedPhaseForPointProcess';
widebandAmplitude.processingInfo.date = datetime;
widebandAmplitude.processingInfo.username = getenv('username');
widebandAmplitude.processingInfo.hostname = getenv('computername');
filenameAmplitude = fullfile(file2find.folder, [basename '.widebandAmplitude.timeseries.mat']);
save(filenameAmplitude, 'widebandAmplitude', '-v7.3');

% Assign output
savedFiles = {filenameFreq; filenamePhase; filenameAmplitude};