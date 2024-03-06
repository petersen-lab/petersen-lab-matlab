function [figH, ax1, ax2, varargout] = adjacentFigs(figPath1, figPath2, options)
% [figH, ax1, ax2, varargout] = adjacentFigs(figPath1, figPath2, <options>)
%
% Function takes two or three separate figures and positions them
% horizontally next to each other as subplots of the same figure. The
% function is useful when creating a live script document.
%
% Args:
%   figPath1 (char, required, positional): a shape-(1,N) character array
%     containing the full figure file path (supports wildcards).
%   figPath2 (char, required, positional): a shape-(1,M) character array
%     containing the full figure file path (supports wildcards).
%   xlim1 (numeric, optional, keyword): a shape-(1,2) numeric array with
%     figure 1 x-axis limits (default = []).
%   xlim2 (numeric, optional, keyword): a shape-(1,2) numeric array with
%     figure 2 x-axis limits (default = []).
%   ylim1 (numeric, optional, keyword): a shape-(1,2) numeric array with
%     figure 1 y-axis limits (default = []).
%   ylim2 (numeric, optional, keyword): a shape-(1,2) numeric array with
%     figure 2 y-axis limits (default = []).
%   legendLocation1 (char, optional, keyword): a shape-(1,L) character
%     array defining legend position for the figure 1. It could take one of
%     the following values: 'NorthEast', 'NorthWest', 'SouthEast',
%     'SouthWest', or 'None' (no legend; default).
%   legendLocation2 (char, optional, keyword): a shape-(1,K) character
%     array defining legend position for the figure 2. It could take one of
%     the following values: 'NorthEast', 'NorthWest', 'SouthEast',
%     'SouthWest', or 'None' (no legend; default).
%   colormap1 (char, optional, keyword): a shape-(1,H) character array
%     controlling the choice of the first figure's colormap. Only two
%     options are currently available: 'hsv' (default) and 'parula'.
%   colormap2 (char, optional, keyword): a shape-(1,F) character array
%     controlling the choice of the second figure's colormap. Only two
%     options are currently available: 'hsv' (default) and 'parula'.
%   figPath3 (char, optional, keyword): a shape-(1,J) character array
%     containing the full figure file path (supports wildcards;
%     default = []).
%   xlim3 (numeric, optional, keyword): a shape-(1,2) numeric array with
%     figure 3 x-axis limits (default = [])
%   ylim3 (numeric, optional, keyword): a shape-(1,2) numeric array with
%     figure 3 y-axis limits (default = []).
%   legendLocation3 (char, optional, keyword): a shape-(1,I) character
%     array defining legend position for the figure 3. It could take one of
%     the following values: 'NorthEast', 'NorthWest', 'SouthEast',
%     'SouthWest', or 'None' (no legend; default).
%   colormap3 (char, optional, keyword): a shape-(1,G) character array
%     controlling the choice of the third figure's colormap. Only two
%     options are currently available: 'hsv' (default) and 'parula'.
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
%   ax1 (graphics): a shape-(1,1) graphics object containing the axes
%     handle 1 from the Live Script figure.
%   ax2 (graphics): a shape-(1,1) graphics object containing the axes
%     handle 2 from the Live Script figure.
%   ax3 (graphics): a shape-(1,1) graphics object containing the axes
%     handle 3 from the Live Script figure (in the case of three figures).
%
% Dependencies:
%   petersen-lab/petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  figPath1 (1,:) {mustBeVector,mustBeText}
  figPath2 (1,:) {mustBeVector,mustBeText}
  options.xlim1 (1,:) {mustBeNumeric,mustBeNonnegative} = []
  options.xlim2 (1,:) {mustBeNumeric,mustBeNonnegative} = []
  options.ylim1 (1,:) {mustBeNumeric,mustBeNonnegative} = []
  options.ylim2 (1,:) {mustBeNumeric,mustBeNonnegative} = []
  options.legendLocation1 (1,:) {mustBeMember(options.legendLocation1,{'NorthEast','NorthWest','SouthEast','SouthWest','None'})} = 'None'
  options.legendLocation2 (1,:) {mustBeMember(options.legendLocation2,{'NorthEast','NorthWest','SouthEast','SouthWest','None'})} = 'None'
  options.colormap1 (1,:) {mustBeVector,mustBeText} = 'hsv'
  options.colormap2 (1,:) {mustBeVector,mustBeText} = 'hsv'
  options.figPath3 (1,:) {mustBeText} = ''
  options.xlim3 (1,:) {mustBeNumeric,mustBeNonnegative} = []
  options.ylim3 (1,:) {mustBeNumeric,mustBeNonnegative} = []
  options.legendLocation3 (1,:) {mustBeMember(options.legendLocation3,{'NorthEast','NorthWest','SouthEast','SouthWest','None'})} = 'None'
  options.colormap3 (1,:) {mustBeVector,mustBeText} = 'hsv'
  options.figSize (1,:) {mustBeNumeric,mustBeNonnegative} = []
  options.figPosition (1,:) {mustBeNumeric,mustBeNonnegative} = []
  options.tight (1,1) {mustBeA(options.tight,'logical')} = false
end

% Parse input
if isempty(options.figPosition)
  if isempty(options.figPath3)
    options.figPosition = [0 0 1200 500];
  else
    options.figPosition = [0 0 1800 500];
  end
else
  options.tight = false;
end

% Draw figure 1
figPath1 = dir(figPath1);
figPath1 = fullfile(figPath1(1).folder, figPath1(1).name);
[figH, ax1] = fig2liveScript(figPath1, xlim=options.xlim1, ...
  ylim=options.ylim1, legendLocation=options.legendLocation1, ...
  figSize=options.figSize, tight=options.tight);
if strcmpi(options.colormap1,'hsv')
  set(ax1, 'Colormap',hsv);
elseif strcmpi(options.colormap1,'parula')
  set(ax1, 'Colormap',parula);
end

% Draw figure 2
figPath2 = dir(figPath2);
figPath2 = fullfile(figPath2(1).folder, figPath2(1).name);
[~, ax2] = fig2liveScript(figPath2, xlim=options.xlim2, ...
  ylim=options.ylim2, legendLocation=options.legendLocation2, ...
  figSize=options.figSize, tight=options.tight);
if strcmpi(options.colormap2,'hsv')
  set(ax2, 'Colormap',hsv);
elseif strcmpi(options.colormap2,'parula')
  set(ax2, 'Colormap',parula);
end
set(ax2, 'Parent',figH);

% Draw figure 3
if ~isempty(options.figPath3)
  figPath3 = dir(options.figPath3);
  figPath3 = fullfile(figPath3(1).folder, figPath3(1).name);
  [~, ax3] = fig2liveScript(figPath3, xlim=options.xlim3, ...
  ylim=options.ylim3, legendLocation=options.legendLocation3, ...
  figSize=options.figSize, tight=options.tight);
  if strcmpi(options.colormap3,'hsv')
    set(ax3, 'Colormap',hsv);
  elseif strcmpi(options.colormap3,'parula')
    set(ax3, 'Colormap',parula);
  end
  set(ax3, 'Parent',figH);
  varargout{1} = ax3;
end

% Position all figures
if exist('ax3','var')
  subplot(1,3,1,ax1);
  subplot(1,3,2,ax2);
  subplot(1,3,3,ax3);
else
  subplot(1,2,1,ax1);
  subplot(1,2,2,ax2);
end
set(figH, 'Position', options.figPosition);