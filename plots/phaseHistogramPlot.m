function fH = phaseHistogramPlot(binCounts, options)
% phaseHistogramPlot(binCounts, binLocs, options)
%
% Function creates a phase histogram plot.
%
% Args:
%   binCounts (numeric, required, positional): a shape-(1, N) numeric array
%     of phase value counts.
%   binLocs (numeric, optional, keyword): a shape-(1, N) numeric array
%     of phase bin values corresponding to phase value counts. By default,
%     the function assumes 9 locations centred around 0 origin of the
%     x-axis.
%   centre (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     representing the central histogram value on the x-axis in radians
%     (default = 0).
%   dataMean (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     representing the mean data value in radians. If supplied, the mean
%     data value will be indicated in the figure as a vertical dotted red
%     line. By default, mean value is not displayed.
%   figText (cell, optional, keyword): a shape-(1, M) cell array of strings
%     to be displayed in the figure corners. Can take up to 4 values
%     corresponding to each figure corner as in
%     {'top-left','top-right','bottom-right','bottom-left'}. By default,
%     no text is displayed.
%   figTitle (char, optional, keyword): a shape-(1, L) character array
%     containing figure title (default = 'Phase Histogram').
%   figPath (char, optional, keyword): a shape-(1, K) character array
%     containing the full folder path for saving the figure. If left empty,
%     the figure will not be saved (default).
%
% Returns:
%   fH (graphics): a shape-(1, 1) graphics object containing the figure
%     handle.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  binCounts (1,:) {mustBeVector,mustBeNumeric,mustBeNonnegative}
  options.binLocs (1,:) {mustBeNumeric} = []
  options.centre (1,1) {mustBeNumeric,mustBeNonNan,mustBeReal} = 0
  options.dataMean (1,:) {mustBeNumeric} = []
  options.figText (1,:) {mustBeUnderlyingType(options.figText,'cell')} = {}
  options.figTitle (1,:) {mustBeText} = 'Phase Histogram'
  options.figPath (1,:) {mustBeText} = ''
end

% Parse input
if isempty(options.binLocs)
  phaseRange = [-pi pi];
  binSize = (phaseRange(2) - phaseRange(1))/options.nBins;
  options.binLocs = phaseRange(1)+binSize/2:binSize:phaseRange(2);
end

% Drawing parameters
fontSize = 14;

% Determine x-axis tick labels
switch options.centre
  case 0
    xTickLocs = [-pi, -pi/2, 0, pi/2, pi];
    xTickLabels = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
  otherwise
    xTickLocs = 'auto';
    xTickLabels = 'auto';
end

% Plot the histogram
fH = figure; bar(options.binLocs, binCounts, ...
  'FaceColor',[.7 .7 .7], 'EdgeColor',[.7 .7 .7]);

% Get axes' dimensions
xLim = xlim;
xAxisLength = xLim(2) - xLim(1);
yLim = ylim;
yAxisLength = yLim(2) - yLim(1);

% Draw the data mean
if ~isempty(options.dataMean)
  hold on; plot([options.dataMean options.dataMean],yLim, 'r:'); hold off
end

% Label the x-axis
xticks(xTickLocs);
xticklabels(xTickLabels);

% Add the corner text messages
if ~isempty(options.figText)
  % Corner coordinates
  xPos = [xLim(1)+0.05*xAxisLength xLim(1)+0.6*xAxisLength];
  yPos = [yLim(1)+0.05*yAxisLength yLim(1)+0.95*yAxisLength];
  cornerCoords = [xPos(1) yPos(2); xPos(2) yPos(2); xPos(2) yPos(1); xPos(1) yPos(1);];
  % Place the text messages
  for iTxt = 1:numel(options.figText)
    if ~isempty(options.figText{iTxt})
      text(cornerCoords(iTxt,1), cornerCoords(iTxt,2), options.figText{iTxt});
    end
  end
end

% Label axes
xlabel('Phase (rad)', 'FontSize',fontSize, 'FontWeight','bold');
ylabel('Spike count', 'FontSize',fontSize, 'FontWeight','bold');

% Add the figure title
if ~isempty(options.figTitle)
  tH = title(options.figTitle, 'Interpreter','none');
else
  tH.String = 'Phase Histogram';
end

% Save the figure
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
end