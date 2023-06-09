function [fullCoherence, half1Coherence, half2Coherence] = coherence(timesSignal, timesReference, options)
% [coh, phi] = coherence(timesSignal, timesReference, <options>)
%
% Function calculates coherence and phase of (a) signal(s) with respect to
% a reference signal.
%
% Args:
%   spikeTimes (cell or numeric, required, positional): a shape-(1, K) cell
%     array of shape-(1, N) numeric arrays of spike times where N
%     corresponds to spike times. Alternatively, one can supply a
%     shape-(1, N) numeric array of spike times (in the case of a single
%     signal vector).
%   timesReference (numeric, required, positional): a shape-(1, N) numeric
%     array of reference spike times where N corresponds to spike times.
%   intervals (numeric, required, positional): a shape-(N, 2) numeric array
%     of time intervals. Each row corresponds to individual time intervals
%     of interest with the first element being the start time and the
%     second element being the end time (default = []).
%   stepsize (numeric, optional, keyword): a shape-(1, 1) nunmeric scalar
%     corresponding to the new sampling interval in seconds after
%     resampling prior to coherence analysis (default = 0.002).
%   startTime (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     corresponding to the start time bin used to generate spike count
%     vectors (or continuous resampled signal) for coherence analysis
%     (default = 0).
%   range (numeric, optional, keyword): a shape-(1, 2) numeric array with
%     the frequency range for estimating coherence and phase values.
%     Default values are [0 0.5/options.stepsize].
%   typespk1 (char, optional, keyword): a shape-(1, M) character array
%     desribing the type of the signal. Can take one of the two values:
%       'pb' - point process (default);
%       'c' - continuous.
%   typespk2 (char, optional, keyword): a shape-(1, M) character array
%     desribing the type of the reference. Can take one of the two values:
%       'pb' - point process (default);
%       'c' - continuous.
%   winfactor (numeric, optional, keyword): a shape-(1, 1) numeric scalar.
%     Each phase/coherence estimation window is at least this many times
%     than 1/(highest frequency). Default is 5.
%   freqfactor (numeric, optional, keyword): a shape-(1, 1) numeric scalar.
%     Each phase/coherence estimation window is at least
%     opt.winfactor/opt.freqfacor times than 1/(lowest frequency). It has
%     to be > 1; default = 2.
%   tapers (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     corresponding to the number of tapers used in phase/coherence
%     calculations (default = 3).
%   decimate (logical, optional, keyword): a shape-(1, 1) logical scalar
%     controlling whether to decimate signals for low frequencies to reduce
%     runtime (default = false).
%   monotoneFreq (logical, optional, keyword): a shape-(1, 1) logical
%     scalar used to remove repeating frequencies caused by transitions
%     across different phase/coherence estimation window sizes
%     (default = true).
%   jack (logical, optional, keyword): a shape-(1, 1) logical scalar for
%     using jackknife error estimates (default = false).
%   pad (numeric, optional, keyword): a shape-(1, 1) scalar representing
%     FFT padding (can take values -1,0,1,2...). -1 corresponds to no
%     padding, 0 corresponds to padding to the next highest power of 2 etc.
%     e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad
%     the FFT to 512 points, if pad=1, we pad to 1024 points etc. Setting
%     pad to larger values, gives denser frequency grids. Defaults to 0.
%   rateAdjust (logical, optional, keyword): a shape-(1, 1) logical scalar
%     to adjusting coherence for firing rate. This only applies to point
%     process signals (default=true).
%   fullCoherence (logical, optional, keyword): a shape-(1, 1) logical
%     scalar for performing coherence analysis on the full signal duration
%     (default = true).
%   halfCoherence (logical, optional, keyword): a shape-(1, 1) logical
%     scalar for performing coherence analysis on signal halves
%     (default = false).
%   parallelise (logical, optional, keyword): a shape-(1, 1) logical
%     scalar for performing coherence analysis using the parfor loop. This
%     might be faster for multiple signals with high sampling frequencies.
%     Otherwise this option does not provide a substantial runtime
%     reduction (default = false).
%
% Returns:
%   fullCoherence (struct): a structure with the following fields:
%     coherence (numeric): a shape-(K, L) numeric array containing
%       coherence values for the signal with respect to the reference
%       (range = [0 1]).
%     coherenceConf (numeric): a shape-(K, L) numeric array containing
%       coherence 95% confidence interval. Add/subtract this interval to
%       actual coherence values to get upper and lower intervals
%       (coherence +/- coherenceConf).
%     phase (numeric): a shape-(K, L) numeric array containing phase radian
%       values for the signal with respect to the reference. Negative phase
%       indicates lag, whereas positive phase indicates lead.
%     phaseConf (cell): a shape-(K, 1) cell array containing phase upper
%       and lower 95% confidence intervals (rad). Each cell contains a
%       shape-(2, L) numeric array.
%     frequency (numeric): a shape-(K, L) numeric array containing
%       frequency values corresponding to coherence and phase estimates.
%     rateAdjustedCoherence (numeric): a shape-(K, L) numeric array
%       containing firing rate adjusted coherence.
%     rateAdjustedCoherenceConf (numeric): a shape-(K, L) numeric array
%       containing firing rate adjusted coherence 95% confidence interval.
%       Add/subtract this interval to actual coherence values to get upper
%       and lower intervals
%       (rateAdjustedCoherence +/- rateAdjustedCoherenceConf).
%     kappaSignal (numeric): a shape-(K, L) numeric array with coherence
%       firing rate adjustment factor kappa1 corresponding to the primary
%       signal.
%     kappaReference (numeric): a shape-(K, L) numeric array with
%       coherence firing rate adjustment factor kappa2 corresponding to the
%       secondary (reference) signal.
%   half1Coherence (struct): a structure with the following fields:
%     coherence (numeric): a shape-(K, J) numeric array containing
%       coherence values for the first half of the signal with respect to
%       the first half of the reference (range = [0 1]).
%     coherenceConf (numeric): a shape-(K, J) numeric array containing
%       coherence 95% confidence interval corresponding to the first half
%       of the signal. Add/subtract this interval to actual coherence
%       values to get upper and lower intervals
%       (coherence +/- coherenceConf).
%     phase (numeric): a shape-(K, J) numeric array containing phase radian
%       values for the first half of the signal with respect to the first
%       half of the reference. Negative phase indicates lag, whereas
%       positive phase indicates lead.
%     phaseConf (cell): a shape-(K, 1) cell array containing phase upper
%       and lower 95% confidence intervals (rad) for the first half
%       of the signal. Each cell contains a shape-(2, J) numeric array.
%     frequency (numeric): a shape-(K, J) numeric array containing
%       frequency values corresponding to the first half coherence and
%       phase estimates.
%     rateAdjustedCoherence (numeric): a shape-(K, J) numeric array
%       containing firing rate adjusted coherence.
%     rateAdjustedCoherenceConf (numeric): a shape-(K, J) numeric array
%       containing firing rate adjusted coherence 95% confidence interval.
%       Add/subtract this interval to actual coherence values to get upper
%       and lower intervals
%       (rateAdjustedCoherence +/- rateAdjustedCoherenceConf).
%     kappaSignal (numeric): a shape-(K, J) numeric array with coherence
%       firing rate adjustment factor kappa1 corresponding to the primary
%       signal.
%     kappaReference (numeric): a shape-(K, J) numeric array with
%       coherence firing rate adjustment factor kappa2 corresponding to the
%       secondary (reference) signal.
%   half2Coherence (struct): a structure with the same fields as
%     half1Coherence but for respective 2nd halves of the signal and the
%     reference comparisons.
%
% Dependencies:
%   Chronux Toolbox (http://chronux.org/).
%   Circular Statistics Toolbox (https://github.com/circstat/circstat-matlab).
%   petersen-lab-matlab repository (https://github.com/petersen-lab/petersen-lab-matlab).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  timesSignal (1,:) {mustBeNumericOrListedType(timesSignal,'cell')}
  timesReference {mustBeVector,mustBeNonnegative,mustBeNonempty}
  options.intervals (:,2) {mustBeNumeric,mustBeNonnegative} = [];
  options.stepsize (1,1) {mustBeNumeric,mustBePositive} = 0.002
  options.startTime (1,1) {mustBeNumeric,mustBeNonnegative} = 0
  options.range (1,2) {mustBeNumeric,mustBeNonnegative} = [0 0]
  options.typespk1 {mustBeMember(options.typespk1,{'pb','c'})} = 'pb'
  options.typespk2 {mustBeMember(options.typespk2,{'pb','c'})} = 'pb'
  options.winfactor (1,1) {mustBeNumeric,mustBePositive} = 5
  options.freqfactor (1,1) {mustBeNumeric,mustBeGreaterThan(options.freqfactor,1)} = 2
  options.tapers (1,1) {mustBeNumeric,mustBePositive} = 3
  options.decimate (1,1) {mustBeA(options.decimate,'logical')} = false
  options.monotoneFreq (1,1) {mustBeA(options.monotoneFreq,'logical')} = true
  options.jack (1,1) {mustBeA(options.jack,'logical')} = false
  options.pad (1,1) {mustBeNumeric} = 0
  options.rateAdjust (1,1) {mustBeA(options.rateAdjust,'logical')} = true
  options.fullCoherence (1,1) {mustBeA(options.fullCoherence,'logical')} = true
  options.halfCoherence (1,1) {mustBeA(options.halfCoherence,'logical')} = false
  options.parallelise (1,1) {mustBeA(options.parallelise,'logical')} = false
end

% Parse input
if ~iscell(timesSignal)
  timesSignal = {timesSignal};
end
timesReference = timesReference(:)';
if options.range(2) == 0
  options.range(2) = 0.5/options.stepsize;
end

% Resample signals (unit rates) and the reference (population rate)
downsampledSignal = resampleSpikesArray(timesSignal, stepsize=options.stepsize, startTime=options.startTime);
[downsampledReference, spikeTimeBins] = resampleSpikes(timesReference, stepsize=options.stepsize, startTime=options.startTime);

% Find indices of times falling within intervals of interest
[~, includeIdx] = selectArrayValues(spikeTimeBins, options.intervals);

% Calculate phase and coherence for all signals
nUnits = numel(timesSignal);
fullCoherence_temp = cell(nUnits,1); % Initialise temporary containers
half1Coherence_temp = cell(nUnits,1);
half2Coherence_temp = cell(nUnits,1);
if options.parallelise
  parfor unit = 1:nUnits
    [fullCoherence_temp{unit}, half1Coherence_temp{unit}, half2Coherence_temp{unit}] = coherenceCalc( ...
      downsampledSignal(unit,includeIdx), downsampledReference(includeIdx), ...
      range=options.range, samplingInterval=options.stepsize, ...
      typespk1=options.typespk1, typespk2=options.typespk2, ...
      winfactor=options.winfactor, freqfactor=options.freqfactor, ...
      tapers=options.tapers, decimate=options.decimate, ...
      monotoneFreq = options.monotoneFreq, jack=options.jack, ...
      pad=options.pad, fullCoherence=options.fullCoherence, ...
      halfCoherence=options.halfCoherence); %#ok<*PFBNS>
  end
else
  for unit = 1:nUnits
    [fullCoherence_temp{unit}, half1Coherence_temp{unit}, half2Coherence_temp{unit}] = coherenceCalc( ...
      downsampledSignal(unit,includeIdx), downsampledReference(includeIdx), ...
      range=options.range, samplingInterval=options.stepsize, ...
      typespk1=options.typespk1, typespk2=options.typespk2, ...
      winfactor=options.winfactor, freqfactor=options.freqfactor, ...
      tapers=options.tapers, decimate=options.decimate, ...
      monotoneFreq = options.monotoneFreq, jack=options.jack, ...
      pad=options.pad, fullCoherence=options.fullCoherence, ...
      halfCoherence=options.halfCoherence);
  end
end

% Repackage cell arrays into matrices
if options.fullCoherence
  fullCoherence_temp = cell2mat(fullCoherence_temp); % Convert into a structure of comma-separated lists
  fullCoherence.coherence = vertcat(fullCoherence_temp.coherence); % Convert lists into matrices
  fullCoherence.coherenceConf = vertcat(fullCoherence_temp.coherenceConf);
  fullCoherence.phase = vertcat(fullCoherence_temp.phase);
  fullCoherence.phaseConf = vertcat({fullCoherence_temp.phaseConf})';
  fullCoherence.frequency = vertcat(fullCoherence_temp.frequency);
else
  fullCoherence.coherence = [];
  fullCoherence.coherenceConf = [];
  fullCoherence.phase = [];
  fullCoherence.phaseConf = [];
  fullCoherence.frequency = [];
end

if options.halfCoherence
  half1Coherence_temp = cell2mat(half1Coherence_temp); % Convert into a structure of comma-separated lists
  half1Coherence.coherence = vertcat(half1Coherence_temp.coherence); % Convert lists into matrices
  half1Coherence.coherenceConf = vertcat(half1Coherence_temp.coherenceConf);
  half1Coherence.phase = vertcat(half1Coherence_temp.phase);
  half1Coherence.phaseConf = vertcat({half1Coherence_temp.phaseConf})';
  half1Coherence.frequency = vertcat(half1Coherence_temp.frequency);
  half2Coherence_temp = cell2mat(half2Coherence_temp); % Convert into a structure of comma-separated lists
  half2Coherence.coherence = vertcat(half2Coherence_temp.coherence); % Convert lists into matrices
  half2Coherence.coherenceConf = vertcat(half2Coherence_temp.coherenceConf);
  half2Coherence.phase = vertcat(half2Coherence_temp.phase);
  half2Coherence.phaseConf = vertcat({half2Coherence_temp.phaseConf})';
  half2Coherence.frequency = vertcat(half2Coherence_temp.frequency);
else
  half1Coherence.coherence = [];
  half1Coherence.coherenceConf = [];
  half1Coherence.phase = [];
  half1Coherence.phaseConf = [];
  half1Coherence.frequency = [];
  half2Coherence.coherence = [];
  half2Coherence.coherenceConf = [];
  half2Coherence.phase = [];
  half2Coherence.phaseConf = [];
  half2Coherence.frequency = [];
end

% Rate-adjust coherence
if options.rateAdjust && (strcmpi(options.typespk1,'pb') || strcmpi(options.typespk2,'pb')) ...
    && (options.fullCoherence || options.halfCoherence)
  % Reference (population rate)
  if strcmpi(options.typespk2,'pb') % Mean firing rates
    [mfrFullReference, mfrHalvesReference] = rateCalc( ...
      downsampledReference(includeIdx), samplingInterval=options.stepsize);
  else
    mfrFullReference = 1;
    mfrHalvesReference = [1 1];
  end
  [fullPSDReference, half1PSDReference, half2PSDReference] = psdCalc( ...
    downsampledReference(includeIdx), range=options.range, ...
    samplingInterval=options.stepsize, typespk1=options.typespk2, ...
    winfactor=options.winfactor, freqfactor=options.freqfactor, ...
    tapers=options.tapers, decimate=options.decimate, ...
    monotoneFreq = options.monotoneFreq, jack=options.jack, ...
    pad=options.pad, fullPSD=options.fullCoherence, ...
    halfPSD=options.halfCoherence); % PSDs

  % Signal (unit rates)
  mfrFullSignal = zeros(nUnits,1); % Initialise containers
  mfrHalvesSignal = zeros(nUnits,2);
  fullPSDSignal = cell(nUnits,1);
  half1PSDSignal = cell(nUnits,1);
  half2PSDSignal = cell(nUnits,1);
  if options.fullCoherence
    fullCoherence.rateAdjustedCoherence = zeros(nUnits,numel(fullPSDReference.psd));
    fullCoherence.rateAdjustedCoherenceConf = zeros(nUnits,numel(fullPSDReference.psd));
    fullCoherence.kappaSignal = zeros(nUnits,numel(fullPSDReference.psd));
    fullCoherence.kappaReference = zeros(nUnits,numel(fullPSDReference.psd));
  else
    fullCoherence.rateAdjustedCoherence = [];
    fullCoherence.rateAdjustedCoherenceConf = [];
    fullCoherence.kappaSignal = [];
    fullCoherence.kappaReference = [];
  end
  if options.halfCoherence
    half1Coherence.rateAdjustedfCoherence = zeros(nUnits,numel(half1PSDReference.psd));
    half1Coherence.rateAdjustedfCoherenceConf = zeros(nUnits,numel(half1PSDReference.psd));
    half1Coherence.kappaSignal = zeros(nUnits,numel(half1PSDReference.psd));
    half1Coherence.kappaReference = zeros(nUnits,numel(half1PSDReference.psd));
    half2Coherence.rateAdjustedfCoherence = zeros(nUnits,numel(half2PSDReference.psd));
    half2Coherence.rateAdjustedfCoherenceConf = zeros(nUnits,numel(half2PSDReference.psd));
    half2Coherence.kappaSignal = zeros(nUnits,numel(half2PSDReference.psd));
    half2Coherence.kappaReference = zeros(nUnits,numel(half2PSDReference.psd));
  else
    half1Coherence.rateAdjustedfCoherence = [];
    half1Coherence.rateAdjustedfCoherenceConf = [];
    half1Coherence.kappaSignal = [];
    half1Coherence.kappaReference = [];
    half2Coherence.rateAdjustedfCoherence = [];
    half2Coherence.rateAdjustedfCoherenceConf = [];
    half2Coherence.kappaSignal = [];
    half2Coherence.kappaReference = [];
  end
  if options.parallelise
    fullCoherence_rateAdjustedCoherence = zeros(nUnits,numel(fullPSDReference.psd)); % These are temporary containers to enable running the parfor loop
    fullCoherence_kappaSignal = zeros(nUnits,numel(fullPSDReference.psd));
    fullCoherence_kappaReference = zeros(nUnits,numel(fullPSDReference.psd));
    fullCoherence_rateAdjustedCoherenceConf = zeros(nUnits,numel(fullPSDReference.psd));
    half1Coherence_rateAdjustedCoherence = zeros(nUnits,numel(half1PSDReference.psd));
    half1Coherence_kappaSignal = zeros(nUnits,numel(half1PSDReference.psd));
    half1Coherence_kappaReference = zeros(nUnits,numel(half1PSDReference.psd));
    half1Coherence_rateAdjustedCoherenceConf = zeros(nUnits,numel(half1PSDReference.psd));
    half2Coherence_rateAdjustedCoherence = zeros(nUnits,numel(half2PSDReference.psd));
    half2Coherence_kappaSignal = zeros(nUnits,numel(half2PSDReference.psd));
    half2Coherence_kappaReference = zeros(nUnits,numel(half2PSDReference.psd));
    half2Coherence_rateAdjustedCoherenceConf = zeros(nUnits,numel(half2PSDReference.psd));
    parfor unit = 1:nUnits
      if strcmpi(options.typespk1,'pb') % Mean firing rates
        [mfrFullSignal(unit), mfrHalvesSignal(unit,:)] = rateCalc( ...
          downsampledSignal(unit,includeIdx), samplingInterval=options.stepsize); %#ok<*PFOUS> 
      else
        mfrFullSignal(unit) = 1;
        mfrHalvesSignal(unit,:) = [1 1];
      end
      [fullPSDSignal{unit}, half1PSDSignal{unit}, half2PSDSignal{unit}] = psdCalc( ...
        downsampledSignal(unit,includeIdx), range=options.range, ...
        samplingInterval=options.stepsize, typespk1=options.typespk1, ...
        winfactor=options.winfactor, freqfactor=options.freqfactor, ...
        tapers=options.tapers, decimate=options.decimate, ...
        monotoneFreq=options.monotoneFreq, jack=options.jack, ...
        pad=options.pad, fullPSD=options.fullCoherence, ...
        halfPSD=options.halfCoherence); % PSDs
      % Full signal
      if options.fullCoherence
        [fullCoherence_rateAdjustedCoherence(unit,:), fullCoherence_kappaSignal(unit,:), ...
          fullCoherence_kappaReference(unit,:)] = coherenceRateAdjustment( ...
          mfrFullSignal(unit), mfrFullReference, fullPSDSignal{unit}.psd, fullPSDReference.psd, ...
          fullCoherence.coherence(unit,:), samplingInterval=options.stepsize); % Full coherence
        fullCoherence_rateAdjustedCoherenceConf(unit,:) = ... % Full coherence 95% confidence intervals
          fullCoherence_kappaSignal(unit,:).*fullCoherence_kappaReference(unit,:).*fullCoherence.coherenceConf(unit,:);
      end
      if options.halfCoherence
        % The first half of the signal
        mfrHalvesSignal_unit = mfrHalvesSignal(unit,:);
        [half1Coherence_rateAdjustedCoherence(unit,:), half1Coherence_kappaSignal(unit,:), ...
          half1Coherence_kappaReference(unit,:)] = coherenceRateAdjustment( ...
          mfrHalvesSignal_unit(1), mfrHalvesReference(1), half1PSDSignal{unit}.psd, half1PSDReference.psd, ...
          half1Coherence.coherence(unit,:), samplingInterval=options.stepsize); % Full coherence
        half1Coherence_rateAdjustedCoherenceConf(unit,:) = ... % Full coherence 95% confidence intervals
          half1Coherence_kappaSignal(unit,:).*half1Coherence_kappaReference(unit,:).*half1Coherence.coherenceConf(unit,:);
        % The second half of the signal
        [half2Coherence_rateAdjustedCoherence(unit,:), half2Coherence_kappaSignal(unit,:), ...
          half2Coherence_kappaReference(unit,:)] = coherenceRateAdjustment( ...
          mfrHalvesSignal_unit(2), mfrHalvesReference(2), half2PSDSignal{unit}.psd, half2PSDReference.psd, ...
          half2Coherence.coherence(unit,:), samplingInterval=options.stepsize); % Full coherence
        half2Coherence_rateAdjustedCoherenceConf(unit,:) = ... % Full coherence 95% confidence intervals
          half2Coherence_kappaSignal(unit,:).*half2Coherence_kappaReference(unit,:).*half2Coherence.coherenceConf(unit,:);
      end
    end
    if options.fullCoherence
      fullCoherence.rateAdjustedCoherence = fullCoherence_rateAdjustedCoherence; % Converting from parfor to the regular format
      fullCoherence.kappaSignal = fullCoherence_kappaSignal;
      fullCoherence.kappaReference = fullCoherence_kappaReference;
      fullCoherence.rateAdjustedCoherenceConf = fullCoherence_rateAdjustedCoherenceConf;
    end
    if options.halfCoherence
      half1Coherence.rateAdjustedCoherence = half1Coherence_rateAdjustedCoherence;
      half1Coherence.kappaSignal = half1Coherence_kappaSignal;
      half1Coherence.kappaReference = half1Coherence_kappaReference;
      half1Coherence.rateAdjustedCoherenceConf = half1Coherence_rateAdjustedCoherenceConf;
      half2Coherence.rateAdjustedCoherence = half2Coherence_rateAdjustedCoherence;
      half2Coherence.kappaSignal = half2Coherence_kappaSignal;
      half2Coherence.kappaReference = half2Coherence_kappaReference;
      half2Coherence.rateAdjustedCoherenceConf = half2Coherence_rateAdjustedCoherenceConf;
    end
  else
    for unit = 1:nUnits
      if strcmpi(options.typespk1,'pb') % Mean firing rates
        [mfrFullSignal(unit), mfrHalvesSignal(unit,:)] = rateCalc(downsampledSignal(unit,:), ...
          samplingInterval=options.stepsize);
      else
        mfrFullSignal(unit) = 1;
        mfrHalvesSignal(unit,:) = [1 1];
      end
      [fullPSDSignal{unit}, half1PSDSignal{unit}, half2PSDSignal{unit}] = psdCalc( ...
        downsampledSignal(unit,includeIdx), range=options.range, ...
        samplingInterval=options.stepsize, typespk1=options.typespk1, ...
        winfactor=options.winfactor, freqfactor=options.freqfactor, ...
        tapers=options.tapers, decimate=options.decimate, ...
        monotoneFreq=options.monotoneFreq, jack=options.jack, ...
        pad=options.pad, fullPSD=options.fullCoherence, ...
        halfPSD=options.halfCoherence); % PSDs
      % Full signal
      if options.fullCoherence
        [fullCoherence.rateAdjustedCoherence(unit,:), fullCoherence.kappaSignal(unit,:), ...
          fullCoherence.kappaReference(unit,:)] = coherenceRateAdjustment( ...
          mfrFullSignal(unit), mfrFullReference, fullPSDSignal{unit}.psd, fullPSDReference.psd, ...
          fullCoherence.coherence(unit,:), samplingInterval=options.stepsize); % Full coherence
        fullCoherence.rateAdjustedCoherenceConf(unit,:) = ... % Full coherence 95% confidence intervals
          fullCoherence.kappaSignal(unit,:).*fullCoherence.kappaReference(unit,:).*fullCoherence.coherenceConf(unit,:);
      end
      if options.halfCoherence
        % The first half of the signal
        [half1Coherence.rateAdjustedCoherence(unit,:), half1Coherence.kappaSignal(unit,:), ...
          half1Coherence.kappaReference(unit,:)] = coherenceRateAdjustment( ...
          mfrHalvesSignal(unit,1), mfrHalvesReference(1), half1PSDSignal{unit}.psd, half1PSDReference.psd, ...
          half1Coherence.coherence(unit,:), samplingInterval=options.stepsize); % Full coherence
        half1Coherence.rateAdjustedCoherenceConf(unit,:) = ... % Full coherence 95% confidence intervals
          half1Coherence.kappaSignal(unit,:).*half1Coherence.kappaReference(unit,:).*half1Coherence.coherenceConf(unit,:);
        % The second half of the signal
        [half2Coherence.rateAdjustedCoherence(unit,:), half2Coherence.kappaSignal(unit,:), ...
          half2Coherence.kappaReference(unit,:)] = coherenceRateAdjustment( ...
          mfrHalvesSignal(unit,2), mfrHalvesReference(2), half2PSDSignal{unit}.psd, half2PSDReference.psd, ...
          half2Coherence.coherence(unit,:), samplingInterval=options.stepsize); % Full coherence
        half2Coherence.rateAdjustedCoherenceConf(unit,:) = ... % Full coherence 95% confidence intervals
          half2Coherence.kappaSignal(unit,:).*half2Coherence.kappaReference(unit,:).*half2Coherence.coherenceConf(unit,:);
      end
    end
  end
else
  fullCoherence.rateAdjustedCoherence = [];
  fullCoherence.rateAdjustedCoherenceConf = [];
  fullCoherence.kappaSignal = [];
  fullCoherence.kappaReference = [];
  half1Coherence.rateAdjustedfCoherence = [];
  half1Coherence.rateAdjustedfCoherenceConf = [];
  half1Coherence.kappaSignal = [];
  half1Coherence.kappaReference = [];
  half2Coherence.rateAdjustedfCoherence = [];
  half2Coherence.rateAdjustedfCoherenceConf = [];
  half2Coherence.kappaSignal = [];
  half2Coherence.kappaReference = [];
end
end



%% Local functions
function [fullCoherence, half1Coherence, half2Coherence] = coherenceCalc(signal, reference, options)
% [fullCoherence, half1Coherence, half2Coherence] = coherenceCalc(signal, reference, <options>)
%
% Function calculates full and half interval phase and coherence of a
% signal with respect to the reference.
%
% Args:
%   signal (numeric, required, positional): a shape-(1, N) numeric array of
%     signal spike counts or a continuous signal.
%   reference (numeric, required, positional): a shape-(1, N) numeric array
%     of reference spike counts or a continuous reference signal.
%   range (numeric, optional, keyword): a shape-(1, 2) numeric array with
%     the frequency range for estimating coherence and phase values.
%     Default values are [0 0.5/options.samplingInterval].
%   samplingInterval (numeric, optional, keyword): a shape-(1, 1) numeric
%     scalar corresponding to the sampling interval in seconds
%     (default = 0.002).
%   typespk1 (char, optional, keyword): a shape-(1, M) character array
%     desribing the type of the signal. Can take one of the two values:
%       'pb' - point process (default);
%       'c' - continuous.
%   typespk2 (char, optional, keyword): a shape-(1, M) character array
%     desribing the type of the reference. Can take one of the two values:
%       'pb' - point process (default);
%       'c' - continuous.
%   winfactor (numeric, optional, keyword): a shape-(1, 1) numeric scalar.
%     Each phase/coherence estimation window is at least this many times
%     than 1/(highest frequency). Default is 5.
%   freqfactor (numeric, optional, keyword): a shape-(1, 1) numeric scalar.
%     Each phase/coherence estimation window is at least
%     opt.winfactor/opt.freqfacor times than 1/(lowest frequency). It has
%     to be > 1; default = 2.
%   tapers (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     corresponding to the number of tapers used in phase/coherence
%     calculations (default = 3).
%   decimate (logical, optional, keyword): a shape-(1, 1) logical scalar
%     controlling whether to decimate signals for low frequencies to reduce
%     runtime (default = false).
%   monotoneFreq (logical, optional, keyword): a shape-(1, 1) logical
%     scalar used to remove repeating frequencies caused by transitions
%     across different phase/coherence estimation window sizes
%     (default = true).
%   jack (logical, optional, keyword): a shape-(1, 1) logical scalar for
%     using jackknife error estimates (default = false).
%   pad (numeric, optional, keyword): a shape-(1, 1) scalar representing
%     FFT padding (can take values -1,0,1,2...). -1 corresponds to no
%     padding, 0 corresponds to padding to the next highest power of 2 etc.
%     e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad
%     the FFT to 512 points, if pad=1, we pad to 1024 points etc. Setting
%     pad to larger values, gives denser frequency grids. Defaults to 0.
%   fullCoherence (logical, optional, keyword): a shape-(1, 1) logical
%     scalar for performing coherence analysis on the full signal duration
%     (default = true).
%   halfCoherence (logical, optional, keyword): a shape-(1, 1) logical
%     scalar for performing coherence analysis on signal halves
%     (default = false).
%
% Returns:
%   fullCoherence (struct): a structure with the following fields:
%     coherence (numeric): a shape-(1, L) numeric array containing
%       coherence values for the signal with respect to the reference
%       (range = [0 1]).
%     coherenceConf (numeric): a shape-(1, L) numeric array containing
%       coherence 95% confidence intervals. Add/subtract this interval to
%       actual coherence values to get upper and lower intervals
%       (coherence +/- coherenceConf).
%     phase (numeric): a shape-(1, L) numeric array containing phase radian
%       values for the signal with respect to the reference. Negative phase
%       indicates lag, whereas positive phase indicates lead.
%     phaseConf (numeric): a shape-(2, L) numeric array containing phase
%       upper and lower 95% confidence intervals (rad).
%     frequency (numeric): a shape-(1, L) numeric array containing
%       frequency values corresponding to coherence and phase estimates.
%   half1Coherence (struct): a structure with the following fields:
%     coherence (numeric): a shape-(1, K) numeric array containing
%       coherence values for the first half of the signal with respect to
%       the first half of the reference (range = [0 1]).
%     coherenceConf (numeric): a shape-(1, K) numeric array containing
%       coherence 95% confidence interval corresponding to the first half
%       of the signal. Add/subtract this interval to actual coherence
%       values to get upper and lower intervals
%       (coherence +/- coherenceConf).
%     phase (numeric): a shape-(1, K) numeric array containing phase radian
%       values for the first half of the signal with respect to the first
%       half of the reference. Negative phase indicates lag, whereas
%       positive phase indicates lead.
%     phaseConf (numeric): a shape-(2, K) numeric array containing phase
%       upper and lower 95% confidence intervals (rad) for the first half
%       of the signal.
%     frequency (numeric): a shape-(1, K) numeric array containing
%       frequency values corresponding to the first half coherence and
%       phase estimates.
%   half2Coherence (struct): a structure with the same fields as
%     half1Coherence but for respective 2nd halves of the signal and the
%     reference comparisons.
%
% Comments:
%   If supplied input signals are not of the same lengths, one of the
%   signals is end-padded by zeros.
%
% Dependencies:
%   Chronux Toolbox (http://chronux.org/).
%   Circular Statistics Toolbox (https://github.com/circstat/circstat-matlab).
%   petersen-lab-matlab repository (https://github.com/petersen-lab/petersen-lab-matlab).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  signal {mustBeVector,mustBeNonnegative}
  reference {mustBeVector,mustBeNonnegative}
  options.range (1,2) {mustBeNumeric} = [0 0]
  options.samplingInterval (1,1) {mustBeNumeric,mustBePositive} = 0.002
  options.typespk1 {mustBeMember(options.typespk1,{'pb','c'})} = 'pb'
  options.typespk2 {mustBeMember(options.typespk2,{'pb','c'})} = 'pb'
  options.winfactor (1,1) {mustBeNumeric,mustBePositive} = 5
  options.freqfactor (1,1) {mustBeNumeric,mustBeGreaterThan(options.freqfactor,1)} = 2
  options.tapers (1,1) {mustBeNumeric,mustBePositive} = 3
  options.decimate (1,1) {mustBeA(options.decimate,'logical')} = false
  options.monotoneFreq (1,1) {mustBeA(options.monotoneFreq,'logical')} = true
  options.jack (1,1) {mustBeA(options.jack,'logical')} = false
  options.pad (1,1) {mustBeNumeric} = 0
  options.fullCoherence (1,1) {mustBeA(options.fullCoherence,'logical')} = true
  options.halfCoherence (1,1) {mustBeA(options.halfCoherence,'logical')} = false
end

% Parse input
if ~options.fullCoherence && ~options.halfCoherence
  warning('Both full and half signal duration coherence calculations disabled upon calling coherenceCalc function.');
end

signal = signal(:)';
reference = reference(:)';
if options.range(2) == 0
  options.range(2) = 0.5/options.samplingInterval;
end
options.minFreq = options.range(1); % A parameter used by freqDependentWindowCoherence
options.maxFreq = options.range(2); % A parameter used by freqDependentWindowCoherence

[signal, reference] = padSignals(signal, reference);
signal = double(signal);
reference = double(reference);
signal(isnan(signal)) = 0;
reference(isnan(reference)) = 0;

% Preprocess input
tStart = 1;
tEnd = numel(signal);
tMid = round(tEnd/2);
if isfield(options,'typespk1') && strcmpi(options.typespk1,'c')
  signal_1sthalf = signal(tStart:tMid) - mean(signal(tStart:tMid), 'omitnan');
  signal_2ndhalf = signal(tMid+1:tEnd) - mean(signal(tMid+1:tEnd), 'omitnan');
  signal = signal - mean(signal, 'omitnan');
else
  signal_1sthalf = signal(tStart:tMid);
  signal_2ndhalf = signal(tMid+1:tEnd);
end
if isfield(options,'typespk2') && strcmpi(options.typespk2,'c')
  reference_1sthalf = reference(tStart:tMid) - mean(reference(tStart:tMid), 'omitnan');
  reference_2ndhalf = reference(tMid+1:tEnd) - mean(reference(tMid+1:tEnd), 'omitnan');
  reference = reference - mean(reference, 'omitnan');
else
  reference_1sthalf = reference(tStart:tMid);
  reference_2ndhalf = reference(tMid+1:tEnd);
end

% Calculate half phase and coherence
if options.halfCoherence
  if sum(signal_1sthalf) && sum(reference_1sthalf) % 1st halves of both signals are not empty
    [freq_1sthalf, coh_1sthalf, phi_1sthalf, coh_1sthalf_conf, phi_1sthalf_confU, phi_1sthalf_confL] = ...
      freqDependentWindowCoherence(reference_1sthalf', signal_1sthalf', options.samplingInterval, [], options);
  end
  if sum(signal_2ndhalf) && sum(reference_2ndhalf) % 2nd halves of both signals are not empty
    [freq_2ndhalf, coh_2ndhalf, phi_2ndhalf, coh_2ndhalf_conf, phi_2ndhalf_confU, phi_2ndhalf_confL] = ...
      freqDependentWindowCoherence(reference_2ndhalf', signal_2ndhalf', options.samplingInterval, [], options);
    if ~sum(signal_1sthalf) || ~sum(reference_1sthalf) || ~exist('freq_1sthalf','var') % The case when one of the first halves of two signals was empty
      freq_1sthalf = NaN(size(freq_2ndhalf)); coh_1sthalf = NaN(size(coh_2ndhalf)); phi_1sthalf = NaN(size(phi_2ndhalf));
    end
  else % 2nd halves of any one of the two signals are empty
    if sum(signal_1sthalf)
      if ~exist('freq_1sthalf','var') % The case when one of the first halves of two signals was empty
        freq_1sthalf = freqDependentWindowCoherence(reference_1sthalf', [], options.samplingInterval, [], options);
        coh_1sthalf = NaN(size(freq_1sthalf)); phi_1sthalf = NaN(size(freq_1sthalf));
      end
      freq_2ndhalf = freq_1sthalf; coh_2ndhalf = NaN(size(coh_1sthalf)); phi_2ndhalf = NaN(size(phi_1sthalf));
    else % Both signals are totally empty
      freq_1sthalf = freqDependentWindowCoherence(reference_1sthalf', [], options.samplingInterval, [], options);
      coh_1sthalf = NaN(size(freq_1sthalf)); phi_1sthalf = NaN(size(freq_1sthalf));
      freq_2ndhalf = freq_1sthalf; coh_2ndhalf = NaN(size(coh_1sthalf)); phi_2ndhalf = NaN(size(phi_1sthalf));
    end
  end

  % Eliminate phase values with no defined confidence intervals
  if exist('phi_1sthalf','var') && exist('phi_1sthalf_confU','var') && exist('phi_1sthalf_confL','var')
    phi_1sthalf(isnan(phi_1sthalf_confU)) = NaN;
    phi_1sthalf(isnan(phi_1sthalf_confL)) = NaN;
  else
    phi_1sthalf = NaN(size(freq_1sthalf));
    phi_1sthalf_confU = NaN(size(freq_1sthalf)); phi_1sthalf_confL = NaN(size(freq_1sthalf));
    coh_1sthalf = NaN(size(freq_1sthalf)); coh_1sthalf_conf = NaN(size(freq_1sthalf));
  end
  if exist('phi_2ndhalf','var') && exist('phi_2ndhalf_confU','var') && exist('phi_2ndhalf_confL','var')
    phi_2ndhalf(isnan(phi_2ndhalf_confU)) = NaN;
    phi_2ndhalf(isnan(phi_2ndhalf_confL)) = NaN;
  else
    phi_2ndhalf = NaN(size(freq_2ndhalf));
    phi_2ndhalf_confU = NaN(size(freq_2ndhalf)); phi_2ndhalf_confL = NaN(size(freq_2ndhalf));
    coh_2ndhalf = NaN(size(freq_2ndhalf)); coh_2ndhalf_conf = NaN(size(freq_2ndhalf));
  end
end

% Calculate full phase and coherence
if options.fullCoherence
  if sum(signal) && sum(reference)
    [freq, coh, phi, coh_conf, phi_confU, phi_confL] = freqDependentWindowCoherence(reference', signal', options.samplingInterval, [], options);
    % Eliminate phase values with no defined confidence intervals
    phi(isnan(phi_confU)) = NaN;
    phi(isnan(phi_confL)) = NaN;
  else
    freq = freqDependentWindowCoherence(reference', [], options.samplingInterval, [], options);
    coh = NaN(size(freq)); phi = NaN(size(freq)); coh_conf = NaN(size(freq)); phi_confU = NaN(size(freq)); phi_confL = NaN(size(freq));
  end
else
  freq = []; coh = []; phi = []; coh_conf = []; phi_confU = []; phi_confL = [];
end

% Assign output structures
if options.halfCoherence
  half1Coherence.coherence = coh_1sthalf(:)';
  half1Coherence.coherenceConf = coh_1sthalf_conf(:)';
  half1Coherence.phase = phi_1sthalf(:)';
  half1Coherence.phaseConf = [phi_1sthalf_confU(:)'; phi_1sthalf_confL(:)'];
  if exist('freq_1sthalf','var')
    half1Coherence.frequency = freq_1sthalf(:)';
  else
    half1Coherence.frequency = freq_2ndhalf(:)';
  end
  half2Coherence.coherence = coh_2ndhalf(:)';
  half2Coherence.coherenceConf = coh_2ndhalf_conf(:)';
  half2Coherence.phase = phi_2ndhalf(:)';
  half2Coherence.phaseConf = [phi_2ndhalf_confU(:)'; phi_2ndhalf_confL(:)'];
  if exist('freq_2ndhalf','var')
    half2Coherence.frequency = freq_2ndhalf(:)';
  else
    half2Coherence.frequency = freq_1sthalf(:)';
  end
else
  half1Coherence.coherence = [];
  half1Coherence.coherenceConf = [];
  half1Coherence.phase = [];
  half1Coherence.phaseConf = [];
  half1Coherence.frequency = [];
  half2Coherence.coherence = [];
  half2Coherence.coherenceConf = [];
  half2Coherence.phase = [];
  half2Coherence.phaseConf = [];
  half2Coherence.frequency = [];
end

if options.fullCoherence
  fullCoherence.coherence = coh(:)';
  fullCoherence.coherenceConf = coh_conf(:)';
  fullCoherence.phase = phi(:)';
  fullCoherence.phaseConf = [phi_confU(:)'; phi_confL(:)'];
  fullCoherence.frequency = freq(:)';
else
  fullCoherence.coherence = [];
  fullCoherence.coherenceConf = [];
  fullCoherence.phase = [];
  fullCoherence.phaseConf = [];
  fullCoherence.frequency = [];
end
end


function [fullPSD, half1PSD, half2PSD] = psdCalc(signal, options)
% [psdFull, psdHalf] = psdCalc(signal, options)
%
% Function calculates a power spectral density (PSD) of a signal.
%
% Args:
%   signal (numeric, required, positional): a shape-(1, N) numeric array of
%     signal spike counts or a continuous signal.
%   range (numeric, optional, keyword): a shape-(1, 2) numeric array with
%     the frequency range for estimating PSD values. Default values are
%     [0 0.5/options.samplingInterval].
%   samplingInterval (numeric, optional, keyword): a shape-(1, 1) numeric
%     scalar corresponding to the sampling interval in seconds
%     (default = 0.002).
%   typespk1 (char, optional, keyword): a shape-(1, M) character array
%     desribing the type of the signal. Can take one of the two values:
%       'pb' - point process (default);
%       'c' - continuous.
%   winfactor (numeric, optional, keyword): a shape-(1, 1) numeric scalar.
%     Each PSD estimation window is at least this many times than
%     1/(highest frequency). Default is 5.
%   freqfactor (numeric, optional, keyword): a shape-(1, 1) numeric scalar.
%     Each PSD estimation window is at least opt.winfactor/opt.freqfacor
%     times than 1/(lowest frequency). It has to be > 1; default = 2.
%   tapers (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     corresponding to the number of tapers used in PSD calculations
%     (default = 3).
%   decimate (logical, optional, keyword): a shape-(1, 1) logical scalar
%     controlling whether to decimate signals for low frequencies to reduce
%     runtime (default = false).
%   monotoneFreq (logical, optional, keyword): a shape-(1, 1) logical
%     scalar used to remove repeating frequencies caused by transitions
%     across different PSD estimation window sizes (default = true).
%   jack (logical, optional, keyword): a shape-(1, 1) logical scalar for
%     using jackknife error estimates (default = false).
%   pad (numeric, optional, keyword): a shape-(1, 1) scalar representing
%     FFT padding (can take values -1,0,1,2...). -1 corresponds to no
%     padding, 0 corresponds to padding to the next highest power of 2 etc.
%     e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad
%     the FFT to 512 points, if pad=1, we pad to 1024 points etc. Setting
%     pad to larger values, gives denser frequency grids. Defaults to 0.
%   fullPSD (logical, optional, keyword): a shape-(1, 1) logical scalar for
%     performing PSD analysis on the full signal duration (default = true).
%   halfPSD (logical, optional, keyword): a shape-(1, 1) logical scalar for
%     performing PSD analysis on signal halves (default = false).
%
% Returns:
%   fullPSD (struct): a structure with the following fields:
%     psd (numeric): a shape-(1, L) numeric array containing PSD values for
%       the full duration signal.
%     psdConf (numeric): a shape-(2, L) numeric array containing PSD upper
%       and lower 95% confidence intervals.
%     frequency (numeric): a shape-(1, L) numeric array containing
%       frequency values corresponding to the PSD estimate.
%   half1PSD (struct): a structure with the following fields:
%     psd (numeric): a shape-(1, K) numeric array containing PSD values
%       for the first half of the signal.
%     psdConf (numeric): a shape-(2, K) numeric array containing PSD upper
%       and lower 95% confidence intervals for the first half of the signal.
%     frequency (numeric): a shape-(1, K) numeric array containing
%       frequency values corresponding to the first half PSD estimate.
%   half2PSD (struct): a structure with the same fields as half1PSD but for
%     respective 2nd halves of the signal and the reference comparisons.
%
% Dependencies:
%   Chronux Toolbox (http://chronux.org/).
%   Circular Statistics Toolbox (https://github.com/circstat/circstat-matlab).
%   petersen-lab-matlab repository (https://github.com/petersen-lab/petersen-lab-matlab).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  signal {mustBeVector,mustBeNonnegative}
  options.range (1,2) {mustBeNumeric} = [0 0]
  options.samplingInterval (1,1) {mustBeNumeric,mustBePositive} = 0.002
  options.typespk1 {mustBeMember(options.typespk1,{'pb','c'})} = 'pb'
  options.winfactor (1,1) {mustBeNumeric,mustBePositive} = 5
  options.freqfactor (1,1) {mustBeNumeric,mustBeGreaterThan(options.freqfactor,1)} = 2
  options.tapers (1,1) {mustBeNumeric,mustBePositive} = 3
  options.decimate (1,1) {mustBeA(options.decimate,'logical')} = false
  options.monotoneFreq (1,1) {mustBeA(options.monotoneFreq,'logical')} = true
  options.jack (1,1) {mustBeA(options.jack,'logical')} = false
  options.pad (1,1) {mustBeNumeric} = 0
  options.fullPSD (1,1) {mustBeA(options.fullPSD,'logical')} = true
  options.halfPSD (1,1) {mustBeA(options.halfPSD,'logical')} = false
end

% Parse input
if ~options.fullPSD && ~options.halfPSD
  warning('Both full and half signal duration PSD calculations disabled upon calling psdCalc function.');
end
options.minFreq = options.range(1); % A parameter used by freqDependentWindowCoherence
options.maxFreq = options.range(2); % A parameter used by freqDependentWindowCoherence

% Calculate half PSD
if options.halfPSD
  tStart = 1;
  tEnd = numel(signal);
  tMid = round((tEnd+tStart)/2);
  [freq_1sthalf, psd_1sthalf, ~, psd_1sthalf_conf] = freqDependentWindowCoherence(signal(tStart:tMid), [], options.samplingInterval, [], options);
  [freq_2ndhalf, psd_2ndhalf, ~, psd_2ndhalf_conf] = freqDependentWindowCoherence(signal(tMid+1:tEnd), [], options.samplingInterval, [], options);
  assert(numel(freq_1sthalf) == numel(freq_2ndhalf) && max(abs(freq_1sthalf - freq_2ndhalf)) < 1e-9) % Assert that frequency bins match between the two half PSD estimates
end

% Calculate full PSD
if options.fullPSD
  [freq, psd, ~, psd_conf] = freqDependentWindowCoherence(signal, [], options.samplingInterval, [], options);
end

% Assign output structures
if options.halfPSD
  half1PSD.psd = psd_1sthalf(:)';
  half1PSD.psdConf = psd_1sthalf_conf;
  half1PSD.frequency = freq_1sthalf(:)';
  half2PSD.psd = psd_2ndhalf(:)';
  half2PSD.psdConf = psd_2ndhalf_conf;
  half2PSD.frequency = freq_2ndhalf(:)';
else
  half1PSD.psd = [];
  half1PSD.psdConf = [];
  half1PSD.frequency = [];
  half2PSD.psd = [];
  half2PSD.psdConf = [];
  half2PSD.frequency = [];
end

if options.fullPSD
  fullPSD.psd = psd(:)';
  fullPSD.psdConf = psd_conf;
  fullPSD.frequency = freq(:)';
else
  fullPSD.psd = [];
  fullPSD.psdConf = [];
  fullPSD.frequency = [];
end
end


function [mfr, mfrHalves] = rateCalc(signal, options)
% [mfr, mfrHalves] = rateCalc(signal, <samplingInterval>)
%
% Function calculates mean firing rate for the full duration and half
% duration signals.
%
% Args:
%   signal (numeric, required, positional): a shape-(1, N) numeric array of
%     signal spike counts or a continuous signal.
%   samplingInterval (numeric, optional, keyword): a shape-(1, 1) numeric
%     scalar corresponding to the sampling interval in seconds
%     (default = 0.002).
%
% Returns:
%   mfr (numeric): a shape-(1, 1) numeric scalar corresponding to the mean
%     firing rate.
%   mfrHalves (numeric): a shape-(2, 1) numeric array corresponding to the
%     mean firing rates of the two halves of the signal.

arguments
  signal {mustBeVector,mustBeNonnegative}
  options.samplingInterval (1,1) {mustBeNumeric,mustBePositive} = 0.002
end

% Calculate half mean firing rates
tStart = 1;
tEnd = numel(signal);
tMid = round((tEnd+tStart)/2);
mfrHalves = [mean(signal(tStart:tMid)); mean(signal(tMid+1:tEnd))]./options.samplingInterval;

% Calculate the mean firing rate for the full duration signal
mfr = mean(signal)/options.samplingInterval;
end


function [adjustedCoherence, kappa1, kappa2] = coherenceRateAdjustment(mfr1, mfr2, psd1, psd2, coh, options)
% [adjustedCoherence, kappa1, kappa2] = coherenceRateAdjustment(mfr1, mfr2, psd1, psd2, coh, <samplingInterval>)
%
% Computes the rate adjustment factors for coherence and adjusts it
% accordingly.
%
% Args:
%   mfr1 (numeric, required, positional): a shape-(1, 1) numeric scalar
%     with the mean firing rate for the corrected signal.
%   mfr2 (numeric, required, positional): a shape-(1, 1) numeric scalar
%     with the mean firing rate for the secondary signal (mfr2 = 1 for
%     LFP).
%   psd1 (numeric, required, positional): a shape-(1, N) numeric array
%     containing PSD values for the corrected signal.
%   psd2 (numeric, required, positional): a shape-(1, N) numeric array
%     containing PSD values for the secondary (reference) signal.
%   coh (numeric, required, positional): a shape-(1, N) numeric array
%     containing coherence values between the two signals.
%   samplingInterval (numeric, optional, keyword): a shape-(1, 1) numeric
%     scalar corresponding to the sampling interval in seconds
%     (default = 0.002).
%
% Returns:
%   adjustedCoherence (numeric): a shape-(1, N) numeric array containing
%     rate-adjusted coherence values between the two signals.
%   kappa1 (numeric): a shape-(1, N) numeric scalar with coherence firing
%     rate adjustment factor kappa1 corresponding to the primary signal.
%   kappa2 (numeric): a shape-(1, N) numeric scalar with coherence firing
%     rate adjustment factor kappa2 corresponding to the secondary
%     (reference) signal.
%
% Comments:
%   Use: C_n*y (f) = kappa*C_ny (f) in case of LFP comparison and
%        C_n1*n2* (f) = kappa1*kappa2*C_n1n2 (f) in case of two spiking
%        rates.
%
% References:
%   Aoi, MC, Lepage, KQ, Kramer, MA, Eden, UT (2015) Rate-adjusted
%     spike–LFP coherence comparisons from spike-train statistics,
%     240:141-153, Journal of Neuroscience Methods.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  mfr1 (1,1) {mustBeNumeric,mustBeNonnegative}
  mfr2 (1,1) {mustBeNumeric,mustBeNonnegative}
  psd1 {mustBeVector,mustBeNonnegative}
  psd2 {mustBeVector,mustBeNonnegative}
  coh {mustBeVector}
  options.samplingInterval (1,1) {mustBeNumeric,mustBePositive} = 0.002
end

% Calculate rate adjustment factors
kappa1 = kappaCalc(mfr1, mfr2, psd1, samplingInterval=options.samplingInterval);
kappa2 = kappaCalc(mfr2, mfr1, psd2, samplingInterval=options.samplingInterval);

% Adjust coherence
adjustedCoherence = kappa1.*kappa2.*coh;
end


function kappa = kappaCalc(mfr1, mfr2, psd1, options)
% kappa = kappaCalc(mfr1, mfr2, spectrum1, <options>)
%
% Function calculates rate adjustment coefficient kappa used to correct for
% the spike-field or spike-spike coherence when the conditions have
% different firing rates (see Aoi et al., 2015).
%
% Args:
%   mfr1 (numeric, required, positional): a shape-(1, 1) numeric scalar
%     with the mean firing rate for the corrected signal.
%   mfr2 (numeric, required, positional): a shape-(1, 1) numeric scalar
%     with the mean firing rate for the secondary signal (mfr2 = 1 for
%     LFP).
%   psd1 (numeric, required, positional): a shape-(1, N) numeric array
%     containing PSD values for the corrected signal.
%   beta (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     corresponding to homogenous Poisson noise rate (default = 0).
%   samplingInterval (numeric, optional, keyword): a shape-(1, 1) numeric
%     scalar corresponding to the sampling interval in seconds
%     (default = 0.002).
%
% Returns:
%   kappa (numeric): a shape-(1, N) numeric scalar with coherence
%     firing rate adjustment factor kappa.
%
% Comments:
%   Use: C_n*y (f) = kappa*C_ny (f) in case of LFP comparison and
%        C_n1*n2* (f) = kappa1*kappa2*C_n1n2 (f) in case of two spiking
%        rates.
%
% References:
%   Aoi, MC, Lepage, KQ, Kramer, MA, Eden, UT (2015) Rate-adjusted
%     spike–LFP coherence comparisons from spike-train statistics,
%     240:141-153, Journal of Neuroscience Methods.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  mfr1 (1,1) {mustBeNumeric,mustBeNonnegative}
  mfr2 (1,1) {mustBeNumeric,mustBeNonnegative}
  psd1 {mustBeVector,mustBeNonnegative}
  options.beta (1,1) {mustBeNumeric} = 0
  options.samplingInterval (1,1) {mustBeNumeric,mustBePositive} = 0.002
end

alpha = mfr2/mfr1;
kappa_temp = 1 + ((options.samplingInterval^2)*((1/alpha - 1)*mfr1 + options.beta/(alpha^2)))./psd1;
kappa_temp(kappa_temp < 0) = NaN;
kappa = 1./sqrt(kappa_temp);
end