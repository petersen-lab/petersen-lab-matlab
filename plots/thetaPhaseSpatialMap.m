function thetaPhaseSpatialMap(phase, maxChan, frequencies, unitTimes, figTitle, options)
% thetaPhaseSpatialMap(phase, maxChan, frequencies, unitTimes, figTitle, <include>)
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
%   include (logical, optional, keyword): a shape-(M, N) logical array
%     matching the shape of the phase array and marking values to be
%     included in the analysis corresponding to most coherent units. By
%     default, all values are included.
%   figPath (char, optional, keyword): a shape-(1, J) character array
%     specifying the full folder path for saving figures. If left empty,
%     figures are not saved (default = '').
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
  frequencies (1,:) {mustBeVector,mustBePositive}
  unitTimes (:,1) {mustBeA(unitTimes,'cell')}
  figTitle (1,:) {mustBeNonzeroLengthText}
  options.include (:,:) {mustBeA(options.include,'logical')} = []
  options.figPath (1,:) {mustBeText} = '';
end

% Parse input
if isempty(options.include)
  options.include = true(size(phase));
end

for f = 1:numel(frequencies)
  % Calculate correlations
  [r, pval] = corrLinearCircular(phase(options.include(:,f),f), ...
    maxChan(options.include(:,f)), type='circlinearnp');

  % Plot data
  fH = figure; plot(maxChan(options.include(:,f)), ...
    phase(options.include(:,f),f), '.', 'MarkerSize',10)

  % Draw the fitted line
  [~, slope, coefficients] = fitLine(maxChan(options.include(:,f)), ...
    phase(options.include(:,f),f), type='linear-circular');
  xLim = xlim;
  xAxisLength = xLim(2) - xLim(1);
  xAxisStep = xAxisLength/10000;
  x = xLim(1):xAxisStep:xLim(2);
  yFit = x.*slope + coefficients(2);
  hold on; plot(x, yFit, 'k--'); hold off;

  % Label axes and figures
  xlabel('Electrode #')
  ylabel('Phase (rad)')
  tH = title([figTitle ', Frequency: ' num2str(frequencies(f))]);

  % Print stats
  str = ['r=' num2str(r), ', p=' num2str(pval)];
  yLim = ylim;
  yAxisLength = yLim(2) - yLim(1);
  text(xLim(1)+0.6*xAxisLength, yLim(1)+0.05*yAxisLength, str)

  % Print significant unit counts
  nonEmptyUnits = ~cellfun(@isempty,unitTimes);
  unitCount = sum(nonEmptyUnits & options.include(:,f));
  significantUnitCount = sum(~isnan(phase(nonEmptyUnits,f)) & ...
    options.include(nonEmptyUnits,f));
  str = ['n=' num2str(significantUnitCount), '/' num2str(unitCount)];
  text(xLim(1)+0.05*xAxisLength, yLim(1)+0.95*yAxisLength, str)

  % Save figures
  if ~isempty(options.figPath)
    if ~exist(options.figPath,'dir')
      mkdir(options.figPath);
    end
    filename = strrep(tH.String,' ','_');
    filename = strrep(filename, ',','_');
    filename = strrep(filename, '-','_');
    filename = strrep(filename, ':','' );
    filename = strrep(filename, '^','' );
    filename = strrep(filename, '{','' );
    filename = strrep(filename, '}','' );
    filename = strrep(filename, '.','p');
    filename = [options.figPath filesep filename '.fig']; %#ok<*AGROW> 
    savefig(fH,filename,'compact');
  end
end