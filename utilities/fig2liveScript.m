function [figH, axes] = fig2liveScript(fig, options)
% [figH, axes] = fig2liveEditor(fig, <options>)
%
% Function converts a regular matlab figure into a Live Script figure.
%
% Args:
%   fig (graphics | char, required, positional): a shape-(1,1) graphics
%     object corresponding to a figure handle or a shape-(1,N) character
%     array containing the full figure file path.
%   xlim (numeric, optional, keyword): a shape-(1,2) numeric array with
%     x-axis limits (default = []).
%   ylim (numeric, optional, keyword): a shape-(1,2) numeric array with
%     y-axis limits (default = []).
%   legendLocation (char, optional, keyword): a shape-(1,M) character array
%     defining legend position. It could take one of the following values:
%     'NorthEast', 'NorthWest', 'SouthEast', 'SouthWest', or 'None' (no
%     legend; default).
%   figSize (numeric, optional, keyword): a shape-(1,1) or shape-(1,2)
%     numeric array defining the wdth and the height of the figure. If a
%     scalar is supplied, it acts both as height and width (default = []).
%     Increase the size if figure axes labels don't fit in a live script
%     document.
%   tight (logical, optional, keyword): a shape-(1,1) logical scalar for
%     tightly fitting the figures (default = false). It can prevent the
%     figure and its axes' labels being cut off in live script documents.
%     However, this property messes up legends if multiple figure are
%     displayed side by side.
%
% Returns:
%   figH (graphics): a shape-(1,1) graphics object handle containing the
%     Live Script figure.
%   axes (graphics): a shape-(1,1) graphics object containing the axes
%     handle from the Live Script figure.
%
% Dependencies:
%   petersen-lab/petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  fig {mustBeFigureOrListedType(fig,'char')}
  options.xlim {mustBeNumeric,mustBeNonnegative} = []
  options.ylim {mustBeNumeric,mustBeNonnegative} = []
  options.legendLocation {mustBeMember(options.legendLocation,{'NorthEast','NorthWest','SouthEast','SouthWest','None'})} = 'None'
  options.figSize {mustBeNumeric,mustBeNonnegative} = []
  options.tight {mustBeA(options.tight,'logical')} = false
end

% Open the figure
if ischar(fig)
  fig = openfig(fig);
end

% Copy old figure to the new figure
figH = figure;
for child = 1:numel(fig.Children)
  if strcmp(fig.Children(child).Type,'axes')
    set(fig.Children(child), 'Parent',figH);
    axes = gca;
    remainingAxes = numel(fig.Children);
    if remainingAxes <= child
      break
    end
  end
end

% Adjust the new figure
if ~isempty(options.xlim) % Adjust x axis limits
  xlim(options.xlim);
end

if ~isempty(options.ylim) % Adjust y axis limits
  ylim(options.ylim);
end

if ~strcmpi(options.legendLocation,'None') % Add the legend
  for child = 1:numel(figH.Children)
    if strcmp(figH.Children(child).Type,'legend')
      set(figH.Children(child), 'Location',options.legendLocation);
      break
    end
  end
end

if ~isempty(options.figSize) % Adjust figure size
  if numel(options.figSize) == 1
    options.figSize = [options.figSize options.figSize];
  end
  resizeFig(figH, gca, options.figSize(1), options.figSize(2), [0 0], [0 0], 0);
end

if options.tight % Uncover any occluding axes
  figH = tightfig(figH);
end

% Close the old figure
if ischar(fig)
  close(fig);
end