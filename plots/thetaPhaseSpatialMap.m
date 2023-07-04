function thetaPhaseSpatialMap(phase, maxChan, frequencies, unitTimes, figTitle)
% thetaPhaseSpatialMap(phase, maxChan, frequencies, unitTimes, figTitle)
%
% Function calculates phase-electrode correlations and related statistics
% and plots the correlations.
%
% Args:
%   phase (numeric, required, positional): a shape-(M, N) numeric array of
%     phase values in radians. The first dimension corresponds to
%     individual units while the second dimension corresponds to different
%     frequencies.
%   maxChan (numeric, required, positional): a shape-(M, 1) numeric array
%     of recording channels with the largest amplitude waveforms for
%     corresponding units.
%   frequencies (numeric, required, positional): a shape-(1, N) numeric
%     array of frequencies.
%   unitTimes (cell, required, positional): a shape-(M, 1) cell array of
%     shape-(1, K) numeric arrays of spike times corresponding to
%     individual units. K may have different values for different units.
%   figTitle (char, required, positional): a shape-(1, L) character array
%     with the title for individual figures. The title will be appended by
%     the frequency value of a figure.
%
% Returns:
%   None.
%
% Comments:
%   The top left corner shows the ratio of units with significant phase and
%   the total numer of units with non-zero spike counts. The bottom
%   rightcorner shows the correlation coeffiecient and its significance
%   p-value.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com)

arguments
  phase {mustBeNumeric}
  maxChan (:,1) {mustBeVector}
  frequencies (1,:) {mustBeVector}
  unitTimes (:,1) {mustBeA(unitTimes,'cell')}
  figTitle (1,:) {mustBeNonzeroLengthText}
end

for f = 1:numel(frequencies)
  % Calculate correlations
  [r, pval] = corrLinearCircular(phase(:,f), maxChan, type='circlinearnp');

  % Plot
  figure; plot(maxChan, phase(:,f), '.', 'MarkerSize',10)

  % Label axes and figures
  xlabel('Electrode #')
  ylabel('Phase (rad)')
  title([figTitle '. Frequency: ' num2str(frequencies(f))])

  % Print stats
  str = ['r=' num2str(r), ', p=' num2str(pval)];
  xLim = xlim;
  xAxisLength = xLim(2) - xLim(1);
  yLim = ylim;
  yAxisLength = yLim(2) - yLim(1);
  text(xLim(1)+0.6*xAxisLength, yLim(1)+0.05*yAxisLength, str)

  % Print significant unit counts
  nonEmptyUnits = ~cellfun(@isempty,unitTimes);
  unitCount = sum(nonEmptyUnits);
  significantUnitCount = sum(~isnan(phase(nonEmptyUnits,f)));
  str = ['n=' num2str(significantUnitCount), '/' num2str(unitCount)];
  text(xLim(1)+0.05*xAxisLength, yLim(1)+0.95*yAxisLength, str)
end