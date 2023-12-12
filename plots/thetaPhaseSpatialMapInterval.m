function [fullCoherence, thetaPhaseTopography, fullInterpCoherence, ...
  interpThetaPhaseTopography] = thetaPhaseSpatialMapInterval(intervals, ...
  unitSpikeTimes, populationSpikeTimes, maxChan, figTitle, ...
  resamplingInterval, frequencyRange, instantThetaFrequency, ...
  instantThetaFreqTimes, options)
% fullCoherence = thetaPhaseSpatialMapInterval(intervals, ...
%   unitSpikeTimes, populationSpikeTimes, maxChan, figTitle, ...
%   resamplingInterval, frequencyRange, <options>)
%
% A wrapper function for thetaPhaseSpatialMap function. Function calculates
% phase-electrode correlations and related statistics and plots the
% correlations for phase data within specified time intervals.
%
% Args:
%   intervals (numeric, required, positional): a shape-(N, 2) numeric array
%     of time intervals. Each row corresponds to individual time intervals
%     of interest with the first element being the start time and the
%     second element being the end time.
%   unitSpikeTimes (cell, required, positional): a shape-(M, 1) cell array
%     of shape-(1, L) numeric arrays of spike times where M corresponds to
%     individual units and L corresponds to spike times.
%   populationSpikeTimes (numeric, required, positional): a shape-(1, L)
%     numeric array of population (units + MUAs) spike times where L
%     corresponds to spike times.
%   maxChan (numeric, required, positional): a shape-(M, 1) numeric array
%     of recording channels with the largest amplitude waveforms for
%     corresponding units.
%   figTitle (char, required, positional): a shape-(1, L) character array
%     with the title for individual figures. The title will be appended by
%     the frequency value of a figure.
%   resamplingInterval (numeric, required, positional): a shape-(1, 1)
%     numeric scalar representing bin size in seconds for generating spike
%     count vectors.
%   frequencyRange (numeric, required, positional): a shape-(1, 2) numeric
%     array with the frequency range for estimating coherence and phase
%     values.
%   instantThetaFrequency (numeric, required, positional): a shape-(1, I)
%     numeric array of instantaneous theta frequencies.
%   instantThetaFreqTimes (numeric, required, positional): a shape-(1, I)
%     numeric array of timestamps corresponding to instantThetaFrequency.
%   coherenceRange (char, optional, keyword): a shape-(1, K) character
%     array controlling coherence range for including units to correlation
%     analyses. Can take one of the two values:
%       'full' - include all units irrespictive of their coherence values
%                (default).
%       'aboveMean' - include only the top 50% most coherent units.
%       'high' - include only the top 25% most coherent units.
%   parallelise (logical, optional, keyword): a shape-(1, 1) logical scalar
%     for performing coherence analysis using the parfor loop
%     (default = true).
%   include (logical, optional, keyword): a shape-(M, 1) logical array
%     matching the shape of the unitSpikeTimes array and marking units to
%     be included in the analysis correlation analyses. By default, all
%     values are included.
%   maxThetaFrequency (numeric, optional, keyword): a shape-(2, E) numeric
%     array with the first row corresponding to timestamps and the second
%     row corresponding to the theta frequency with the largest power.
%   figPath (char, optional, keyword): a shape-(1, J) character array
%     specifying the full folder path for saving figures. If left empty,
%     figures are not drawn and not saved (default = '').
%
% Returns:
%   fullCoherence (struct): a shape-(1, 1) scalar structure with the
%     following fields:
%     coherence (numeric): a shape-(G, H) numeric array containing
%       coherence values for the signal with respect to the reference
%       (range = [0 1]).
%     coherenceConf (numeric): a shape-(G, H) numeric array containing
%       coherence 95% confidence interval. Add/subtract this interval to
%       actual coherence values to get upper and lower intervals
%       (coherence +/- coherenceConf).
%     phase (numeric): a shape-(G, H) numeric array containing phase radian
%       values for the signal with respect to the reference. Negative phase
%       indicates lag, whereas positive phase indicates lead.
%     phaseConf (cell): a shape-(G, 1) cell array containing phase upper
%       and lower 95% confidence intervals (rad). Each cell contains a
%       shape-(2, L) numeric array.
%     frequency (numeric): a shape-(G, H) numeric array containing
%       frequency values corresponding to coherence and phase estimates.
%     rateAdjustedCoherence (numeric): a shape-(G, H) numeric array
%       containing firing rate adjusted coherence.
%     rateAdjustedCoherenceConf (numeric): a shape-(G, H) numeric array
%       containing firing rate adjusted coherence 95% confidence interval.
%       Add/subtract this interval to actual coherence values to get upper
%       and lower intervals
%       (rateAdjustedCoherence +/- rateAdjustedCoherenceConf).
%     kappaSignal (numeric): a shape-(G, H) numeric array with coherence
%       firing rate adjustment factor kappa1 corresponding to the primary
%       signal.
%     kappaReference (numeric): a shape-(G, H) numeric array with
%       coherence firing rate adjustment factor kappa2 corresponding to the
%       secondary (reference) signal.
%   thetaPhaseTopography (struct): a shape-(1, 1) scalar structure with the
%     following fields:
%     frequency (numeric): a shape-(1, H) numeric array containing
%       frequency values corresponding to columns of phase estimates.
%     phase (numeric): a shape-(G, H) numeric array containing phase radian
%       values for the signal with respect to the reference. Negative phase
%       indicates lag, whereas positive phase indicates lead.
%     maxChan (numeric): a shape-(M, 1) numeric array of recording channels
%       with the largest amplitude waveforms for corresponding units.
%     includeUnits (logical): a shape-(G, H) logical array corresponding to
%       the shape of the phase array and indicating which phase values were
%       included in correlation analyses.
%     r (numeric): a shape-(1, H) numeric array with correlation
%       coefficients for a circular-linear correlation analysis between
%       individual phase array columns and the maxChan array.
%     pval (numeric): a shape-(1, H) numeric array with correlation
%       coefficient significance p-values for a circular-linear correlation
%       analysis between individual phase array columns and the maxChan
%       array.
%     coefficients (numeric): a shape-(H, 2) numeric array with linear
%       regression equation coeffiecients corresponding to correlation
%       coefficients (slopes in column 1 and y-intercepts in column 2). The
%       first dimension of the array corresponds to frequencies.
%     fIndMostCoh (numeric): a shape-(1, 1) numeric scalar indexing
%       frequency with most coherent phase values on average across all
%       units.
%     fIndMostSignificant (numeric): a shape-(1, 1) numeric scalar indexing
%       frequency that is the most significant (the lowest p-value).
%   fullInterpCoherence (struct): a shape-(1, 1) scalar structure with the
%     same fields as fullCoherence except for mostCoherentFrequency field.
%     The difference here is that all values are interpolated with 0.01 Hz
%     intersample interval as follows:
%     frequencyRange(1):0.01:frequencyRange(2).
%   interpThetaPhaseTopography (struct): a shape-(1, 1) scalar structure
%     with the following fields:
%     frequency (numeric): a shape-(1, F) numeric array containing
%       frequency values corresponding to columns of interpolated phase
%       estimates.
%     phase (numeric): a shape-(G, F) numeric array containing interpolated
%       phase radian values for the signal with respect to the reference.
%       Negative phase indicates lag, whereas positive phase indicates lead.
%     maxChan (numeric): a shape-(M, 1) numeric array of recording channels
%       with the largest amplitude waveforms for corresponding units.
%     includeUnits (logical): a shape-(G, F) logical array corresponding to
%       the shape of the interpolated phase array and indicating which
%       phase values were included in correlation analyses.
%     r (numeric): a shape-(1, 2) numeric array with correlation
%       coefficients for a circular-linear correlation analysis between
%       certain phase array columns and the maxChan array. The first value
%       corresponds to the frequency with the largest spectrogram power.
%       The second value corresponds to the average instantaneous frequency
%       over the period of interest estimated by the Hilbert Transform.
%     pval (numeric): a shape-(1, 2) numeric array with correlation
%       coefficient significance p-values for a circular-linear correlation
%       analysis between individual certain phase array columns and the
%       maxChan array.  The first value corresponds to the frequency with
%       the largest spectrogram power. The second value corresponds to the
%       average instantaneous frequency over the period of interest
%       estimated by the Hilbert Transform.
%     coefficients (numeric): a shape-(2, 2) numeric array with linear
%       regression equation coeffiecients corresponding to correlation
%       coefficients (slopes in column 1 and y-intercepts in column 2). The
%       first dimension of the array corresponds to frequencies. The first
%       frequency has the largest histogram power, whereas the second one
%       is the average instantaneous frequency estimated by Hilbert
%       Transform.
%     fIndMostPower (numeric): a shape-(1, 1) numeric scalar indexing
%       frequency with the largest power within the theta frequency band
%       for the population rate.
%     fIndMeanInstant (numeric): a shape-(1, 1) numeric scalar indexing
%       average instantenous frequency.
%
% Comments:
%   The top left corner shows the ratio of units with significant phase and
%   the total numer of units with non-zero spike counts. The top right
%   corner shows the fitted line equation. The bottom right corner shows
%   the correlation coeffiecient and its significance p-value.
%
% Dependencies:
%   CellExplorer (https://cellexplorer.org/).
%   Circular Statistics Toolbox
%     (https://uk.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics).
%   dervinism/circStatNP (https://github.com/dervinism/circStatNP).
%   petersen-lab/petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).
%   FMA Toolbox (https://github.com/michael-zugaro/FMAToolbox).
%   Chronux Toolbox (http://chronux.org/).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com)

arguments
  intervals (:,2) {mustBeNumeric,mustBeNonnegative}
  unitSpikeTimes (:,1) {mustBeA(unitSpikeTimes,'cell')}
  populationSpikeTimes {mustBeVector,mustBeNonnegative,mustBeNonempty}
  maxChan (:,1) {mustBeVector,mustBeNonnegative,mustBeNonempty}
  figTitle (1,:) {mustBeNonzeroLengthText}
  resamplingInterval (1,1) {mustBeNumeric,mustBePositive}
  frequencyRange (1,2) {mustBeNumeric,mustBeNonnegative}
  instantThetaFrequency (1,:) {mustBeVector,mustBeNonempty}
  instantThetaFreqTimes (1,:) {mustBeVector,mustBeNonnegative,mustBeNonempty}
  options.coherenceRange {mustBeMember(options.coherenceRange,{'full','aboveMean','high'})} = 'full'
  options.parallelise (1,1) {mustBeA(options.parallelise,'logical')} = true
  options.include (:,:) {mustBeA(options.include,'logical')} = []
  options.maxThetaFrequency (2,:) {mustBeNumeric,mustBeNonnegative} = []
  options.figPath (1,:) {mustBeText} = ''
end

%% Parse input
if strcmpi(options.coherenceRange,'aboveMean')
  n = 51;
elseif strcmpi(options.coherenceRange,'high')
  n = 76;
end

if isempty(options.include)
  options.include = true(size(unitSpikeTimes));
end

%% Work out intervals of interest
if ~isempty(intervals)
  unitSpikeTimes = selectCellValues(unitSpikeTimes, intervals);
  populationSpikeTimes = selectArrayValues(populationSpikeTimes, intervals);
else
  unitSpikeTimes = [];
  populationSpikeTimes = [];
end

%% Carry out coherence and correlation analyses
if isempty(unitSpikeTimes) || isempty(populationSpikeTimes)
  warning('Empty spike count vectors supplied for coherence analysis. Quiting...');
  fullCoherence = [];
  thetaPhaseTopography = [];
  fullInterpCoherence = [];
  interpThetaPhaseTopography = [];
else
  % Calculate coherence
  freqGrid = frequencyRange(1):0.01:frequencyRange(2);
  [fullCoherence, ~, ~, fullInterpCoherence] = coherence(unitSpikeTimes, ...
    populationSpikeTimes, intervals=intervals, stepsize=resamplingInterval, ...
    freqRange=frequencyRange, freqGrid=freqGrid, parallelise=options.parallelise);

  % Label coherent units
  if strcmpi(options.coherenceRange,'full')
    includeUnits = true(size(fullCoherence.coherence));
    includeUnitsInterp = true(size(fullInterpCoherence.coherence));
  else
    quantiles = quantile(fullCoherence.rateAdjustedCoherence,100);
    coherenceThr = quantiles(n,:);
    includeUnits = false(size(fullCoherence.rateAdjustedCoherence));
    for f = 1:numel(coherenceThr)
      includeUnits(fullCoherence.rateAdjustedCoherence(:,f) >= coherenceThr(f),f) = true;
    end
    quantiles = quantile(fullInterpCoherence.rateAdjustedCoherence,100);
    coherenceThr = quantiles(n,:);
    includeUnitsInterp = false(size(fullInterpCoherence.rateAdjustedCoherence));
    for f = 1:numel(coherenceThr)
      includeUnitsInterp(fullInterpCoherence.rateAdjustedCoherence(:,f) >= coherenceThr(f),f) = true;
    end
  end
  include = repmat(options.include, 1, size(fullCoherence.coherence,2));
  includeUnits = logical(include & includeUnits);
  if ~any(includeUnits)
    return
  end
  include = repmat(options.include, 1, size(fullInterpCoherence.coherence,2));
  includeUnitsInterp = logical(include & includeUnitsInterp);

  % plot phase spatial map
  if ~isempty(options.figPath)
    % All frequencies
    [r, pval, coefficients] = thetaPhaseSpatialMap(fullCoherence.phase, ...
      maxChan, fullCoherence.frequency(1,:), unitSpikeTimes, figTitle, ...
      include=includeUnits, figPath=options.figPath);

    % Most coherent frequency (coherence analysis output)
    [~, maxIdx] = max(fullCoherence.rateAdjustedCoherence(includeUnits),[],2,'omitnan','linear');
    frequency = fullCoherence.frequency(includeUnits);
    mostCoherentFrequency = mode(frequency(maxIdx));
    fIndMostCoh = find(frequency(1,:) == mostCoherentFrequency);
    figTitleMostCoh = [figTitle ' (most coherent freq)'];
    thetaPhaseSpatialMap(fullCoherence.phase(:,fIndMostCoh), maxChan, ...
      fullCoherence.frequency(1,fIndMostCoh), unitSpikeTimes, figTitleMostCoh, ...
      include=includeUnits(:,fIndMostCoh), figPath=options.figPath);

    % Most significant frequency (coherence analysis output)
    [~, fIndMostSignificant] = min(pval);
    figTitleMostSignificant = [figTitle ' (most significant freq)'];
    thetaPhaseSpatialMap(fullCoherence.phase(:,fIndMostSignificant), ...
      maxChan, fullCoherence.frequency(1,fIndMostSignificant), ...
      unitSpikeTimes, figTitleMostSignificant, ...
      include=includeUnits(:,fIndMostSignificant), figPath=options.figPath);
    
    % Most powerful frequency (instantaneous theta frequency analysis output)
    if ~isempty(options.maxThetaFrequency)
      maxThetaFrequency = [];
      for interval = 1:size(intervals,1)
        maxThetaFrequency = [maxThetaFrequency options.maxThetaFrequency(2, ...
          options.maxThetaFrequency(1,:) >= intervals(interval,1) & ...
          options.maxThetaFrequency(1,:) <= intervals(interval,2))]; %#ok<*AGROW>
      end
      maxThetaFrequency = mean(maxThetaFrequency);
      [~,fIndMostPower] = min(abs(fullInterpCoherence.frequency(1,:) - maxThetaFrequency));
      figTitleMostPower = [figTitle ' (most powerful freq)'];
      if any(~isnan(fullInterpCoherence.phase(:,fIndMostPower)))
        [rMostPower, pvalMostPower, coefficientsMostPower] = ...
          thetaPhaseSpatialMap(fullInterpCoherence.phase(:,fIndMostPower), ...
          maxChan, fullInterpCoherence.frequency(1,fIndMostPower), ...
          unitSpikeTimes, figTitleMostPower, ...
          include=includeUnitsInterp(:,fIndMostPower), figPath=options.figPath);
      else
        rMostPower = NaN;
        pvalMostPower = NaN;
        coefficientsMostPower = [NaN NaN];
      end
    end

    % Average frequency estimated using Hilbert Transform (instantaneous theta frequency analysis output)
    [~, instantThetaFrequencyInds] = selectArrayValues(instantThetaFreqTimes, intervals);
    instantThetaFrequency = instantThetaFrequency(instantThetaFrequencyInds);
    meanInstantThetaFrequency = mean(instantThetaFrequency, 'omitnan');
    [~,fIndMeanInstant] = min(abs(fullInterpCoherence.frequency(1,:) - meanInstantThetaFrequency));
    figTitleMeanInstant = [figTitle ' (mean instant freq)'];
    if any(~isnan(fullInterpCoherence.phase(:,fIndMeanInstant)))
      [rInst, pvalInst, coefficientsInst] = thetaPhaseSpatialMap( ...
        fullInterpCoherence.phase(:,fIndMeanInstant), maxChan, ...
        fullInterpCoherence.frequency(1,fIndMeanInstant), unitSpikeTimes, ...
        figTitleMeanInstant, include=includeUnitsInterp(:,fIndMeanInstant), ...
        figPath=options.figPath);
    else
        rInst = NaN;
        pvalInst = NaN;
        coefficientsInst = [NaN NaN];
    end
  end

  % Assign output
  thetaPhaseTopography.frequency = fullCoherence.frequency(1,:);
  thetaPhaseTopography.phase = fullCoherence.phase;
  thetaPhaseTopography.includeUnits = includeUnits;
  thetaPhaseTopography.maxChan = maxChan;
  thetaPhaseTopography.r = r;
  thetaPhaseTopography.pval = pval;
  thetaPhaseTopography.coefficients = coefficients;
  thetaPhaseTopography.fIndMostCoh = fIndMostCoh;
  thetaPhaseTopography.fIndMostSignificant = fIndMostSignificant;

  interpThetaPhaseTopography.frequency = freqGrid;
  interpThetaPhaseTopography.phase = fullInterpCoherence.phase;
  interpThetaPhaseTopography.includeUnits = includeUnitsInterp;
  interpThetaPhaseTopography.maxChan = maxChan;
  interpThetaPhaseTopography.r = [rInst rMostPower];
  interpThetaPhaseTopography.pval = [pvalInst pvalMostPower];
  interpThetaPhaseTopography.coefficients = [coefficientsInst; coefficientsMostPower];
  interpThetaPhaseTopography.fIndMostPower = fIndMostPower;
  interpThetaPhaseTopography.fIndMeanInstant = fIndMeanInstant;
end