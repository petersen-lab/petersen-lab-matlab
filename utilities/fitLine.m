function [yFit, slope, coefficients] = fitLine(x, y, options)
% [yFit, slope, coefficients] = fitLine(x, y, <type>)
%
% Function fits a line to the data y given x.
% 
% Args:
%   x (numeric, required, positional): a shape-(1, N) numeric array of data
%     values (x-axis; could be independent variable).
%   y (numeric, required, positional): a shape-(1, N) numeric array of data
%     values (y-axis; could be dependent variable).
%   type (char, optional, keyword): a shape-(1, M) character array
%     desribing data types. Could be one of the following values:
%       'linear-linear' - meaning both data types are linear (default).
%       'linear-circular' - meaning that x is linear and y is circular.
% 
% Returns:
%   yFit (numeric): a shape-(1, L) numeric array of fitted y values
%     corresponding to input variable x.
%   slope (numeric): a shape-(1, 1) numeric scalar corresponding to the
%     slope of the fitted line.
%   coefficients (numeric): a shape-(1, 2) numeric array containing
%     coefficients a (slope) and b (y-intercept) of the fitted line
%     (y = a*x + b).
%
% Comments:
%   If x contains non-unique values, they are adjusted adding random small
%   values of noise to x values. The input x and y vectors are further
%   sorted to make x monotonicly increasing. If y vector contains NaNs,
%   they are omitted.
%
%   You may need to replace the following line in CircularRegression
%   function:
%     TSE = norm(phi-CircularMean(phi))^2;
%   with the following line:
%     TSE = norm(phi-circ_mean(phi))^2;
%
% Dependencies:
%   FMA Toolbox (https://github.com/michael-zugaro/FMAToolbox).
%   Circular Statistics Toolbox (https://github.com/circstat/circstat-matlab).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  x (1,:) {mustBeVector}
  y (1,:) {mustBeVector}
  options.type {mustBeMember(options.type,{'linear-linear','linear-circular'})} = 'linear-linear'
end

% Adjust for non-unique x-values
if numel(x) - numel(unique(x)) > 0
  xAdjusted = x + 1e-9 * rand(1, length(x));
else
  xAdjusted = x;
end

% Make sure x values are monotonicly increasing
[xAdjusted, monotonicIdx] = sort(xAdjusted);
yAdjusted = y(monotonicIdx);

% Exclude NaN values
points2exclude = isnan(yAdjusted);
xAdjusted = xAdjusted(~points2exclude);
yAdjusted = yAdjusted(~points2exclude);

% Fit the line
if strcmpi(options.type, 'linear-linear')
  coefficients = polyfit(xAdjusted, yAdjusted, 1);
  yFit = polyval(coefficients , xAdjusted);
  slope = coefficients(1);
elseif strcmpi(options.type, 'linear-circular')
  coefficients = CircularRegression(xAdjusted, yAdjusted);
  yFit = coefficients(1).*xAdjusted + coefficients(2);
  slope = coefficients(1);
end