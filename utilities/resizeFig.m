function paperSize = resizeFig(f, ax, width, height, label, margin, gap)
% paperSize = resizeFig(f, ax, width, height, label, margin, gap)
%
% Function resizes and crops a figure.
%
% Args:
%   f (graphics, required, positional): a shape-(1,1) graphics object
%     corresponding to a figure handle.
%   ax (graphics, required, positional): a shape-(1,M) graphics object
%     containing axes objects corresponding to the figure handle;
%   width (numeric, required, positional): a shape-(1,1) numeric scalar
%     corresponding to the figure width in centimeters.
%   height (numeric, required, positional): a shape-(1,1) numeric scalar
%     corresponding to the figure height in centimeters.
%   label (numeric, required, positional): a shape-(1,2) numeric array
%     corresponding to the size of y and x label areas in centimeters.
%   margin (numeric, required, positional): a shape-(1,2) numeric array
%     corresponding to the size of right and top figure margins in
%     centimeters.
%   gap (numeric, optional, positional): a shape-(1,1) numeric scalar
%     corresponding to the distance between graphs, if there are more than
%     one axes object to display in the figure (default = 0)
%
% Returns:
%   paperSize (numeric): a shape-(1,2) numeric array with figure
%     dimensions [width height].
%
% Dependencies:
%   petersen-lab/petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com)

arguments
  f {mustBeFigureOrListedType(f)}
  ax {mustBeFigureOrListedType(ax)}
  width {mustBeNumeric,mustBePositive}
  height {mustBeNumeric,mustBePositive}
  label {mustBeVector,mustBeNumeric,mustBeNonnegative}
  margin {mustBeVector,mustBeNumeric,mustBeNonnegative}
  gap {mustBeNumeric,mustBeNonnegative} = 0
end

% Figure dimensions
paperSize = [label(1)+width+margin(1) label(2)+height+margin(2)];

% Convert to pixels
width = width*40;
height = height*40;
label = label*40;
margin = margin*40;
gap = gap*40;
screenSize = get(0, 'ScreenSize');
screenSize = [screenSize(3) screenSize(4) screenSize(3) screenSize(4)];

fPos(1) = (screenSize(1)-width)/2 - label(1);
fPos(2) = (screenSize(2)-height)/2 - label(2);
fPos(3) = label(1)+width+margin(1);
fPos(4) = label(2)+height+margin(2);
set(f, 'Units', 'pixels', 'Position', fPos);
set(f, 'Units', 'normalized');

% Position axes
for i = 1:length(ax)
    axPos(1) = label(1);
    axPos(2) = label(2) + (length(ax)-i)*((height-gap*(length(ax)-1))/length(ax)+gap);
    axPos(3) = width;
    axPos(4) = (height-gap*(length(ax)-1))/length(ax);
    set(ax(i), 'Units', 'pixels', 'Position', axPos);
    set(ax(i), 'Units', 'normalized');
end