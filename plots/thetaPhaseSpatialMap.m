function [r, pval, coefficients] = thetaPhaseSpatialMap(phase, maxChan, frequencies, unitTimes, figTitle, options)
% thetaPhaseSpatialMap(phase, maxChan, frequencies, unitTimes, figTitle, <options>)
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
%   fitFunc (char, optional, keyword): a shape-(1, J) character array
%     specifying which function should be used to fit the linear-circular
%     regression line. The available options are:
%       'fma' - CircularRegression function from FMA Toolbox (default)
%               (https://github.com/michael-zugaro/FMAToolbox).
%       'pp' - CircularLinearRegression used by Peter Petersen
%              (https://github.com/petersenpeter/Code_Petersen_Buzsaki_Neuron_2020/blob/master/CircularLinearRegression.m).
%       'fma&pp' - both functions.
%   figPath (char, optional, keyword): a shape-(1, K) character array
%     specifying the full folder path for saving figures. If left empty,
%     figures are not saved (default = '').
%
% Returns:
%   r (numeric): a shape-(1, N) numeric array of correlation coefficient
%     for all phase frequencies.
%   pval (numeric): a shape-(1, N) numeric array of p-value of
%     corresponding correlation coefficients for all phase frequencies.
%   coefficients (numeric): a shape-(N, 2) numeric array of the fitted line
%     equation coeffients (slope and y-intercept).
%
% Comments:
%   The top left corner shows the ratio of units with significant phase and
%   the total numer of units with non-zero spike counts. The bottom
%   rightcorner shows the correlation coeffiecient and its significance
%   p-value.
%
% Dependencies:
%   CellExplorer (https://cellexplorer.org/).
%   Circular Statistics Toolbox
%     (https://uk.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics).
%   dervinism/circStatNP (https://github.com/dervinism/circStatNP).
%   petersen-lab/petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).
%   FMA Toolbox (https://github.com/michael-zugaro/FMAToolbox).
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
  options.fitFunc (1,:) {mustBeMember(options.fitFunc, {'fma','pp','fma&pp'})} = 'fma'
  options.figPath (1,:) {mustBeText} = ''
end

% Parse input
if isempty(options.include)
  options.include = true(size(phase));
end

nFreq = numel(frequencies);
r = zeros(1,nFreq);
pval = zeros(1,nFreq);
coefficients = zeros(nFreq,2);
for f = 1:nFreq
  if any(~isnan(phase(options.include(:,f),f))) && ~isempty(phase(options.include(:,f),f))
    % Calculate correlations
    [r(f), pval(f)] = corrLinearCircular(phase(options.include(:,f),f), ...
      maxChan(options.include(:,f)), type='circlinearnp');

    % Plot data
    fH = figure; plot(maxChan(options.include(:,f)), ...
      phase(options.include(:,f),f), '.', 'MarkerSize',10)

    % Draw the fitted line
    xLim = xlim;
    xAxisLength = xLim(2) - xLim(1);
    xAxisStep = xAxisLength/10000;
    x = xLim(1):xAxisStep:xLim(2);
    if strcmpi(options.fitFunc, 'fma') || strcmpi(options.fitFunc, 'fma&pp')
      [~, slope, coefficients(f,:)] = fitLine(maxChan(options.include(:,f)), ...
        phase(options.include(:,f),f), type='linear-circular-fma');
      yFit = x.*slope + coefficients(f,2);
      hold on; p1 = plot(x, yFit, 'k--'); hold off;
    end
    if strcmpi(options.fitFunc, 'pp') || strcmpi(options.fitFunc, 'fma&pp')
      [~, slope, coefficients(f,:)] = fitLine(maxChan(options.include(:,f)), ...
        phase(options.include(:,f),f), type='linear-circular-pp');
      yFit = x.*slope + coefficients(f,2);
      hold on; p2 = plot(x, yFit, 'b--'); hold off;
    end
    if strcmpi(options.fitFunc, 'fma&pp')
      legend([p1,p2], {'fma','pp'})
    end

    % Adjust the correlation direction based on the slope of the regression line
    if slope < 0
      r(f) = -abs(r(f));
    end

    % Label axes and figures
    xlabel('Electrode #')
    ylabel('Phase (rad)')
    tH = title([figTitle ', Freq: ' num2str(round(frequencies(f),1))], 'Interpreter','none');

    % Print stats
    str = ['r=' num2str(r(f)), ', p=' num2str(pval(f))];
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

    % Print the line equation
    if coefficients(f,2) >= 0
      lineEqStr = [num2str(coefficients(f,1)) 'x+' num2str(coefficients(f,2))];
    else
      lineEqStr = [num2str(coefficients(f,1)) 'x' num2str(coefficients(f,2))];
    end
    text(xLim(1)+0.6*xAxisLength, yLim(1)+0.95*yAxisLength, lineEqStr)

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
      title('');
      saveas(fH,filename(1:end-4),'png');
      close(fH);
    end
  end
end