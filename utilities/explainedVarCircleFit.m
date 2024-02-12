function explainedVariance = explainedVarCircleFit(x, y, xc, yc, R)
% explainedVariance = explainedVarCircleFit(x, y, xc, yc, R)
%
% Function calculates variance explained by a fitted circle to the data.
%
% Args:
%   x, y (numeric, required, positional): a shape-(1, N) numeric arrays
%     with data coordinates on the x and y axes.
%   xc, yc (numeric, required, positional): a shape-(1, 1) numeric scalars
%     with centre of the fitted circle in x and y coordinates.
%   R (numeric, required, positional): a shape-(1, 1) numeric scalar with
%     the radius of the fitted circle.
%
% Returns:
%   explainedVariance (numeric): a shape-(1, 1) numeric scalar with
%     explained variance proportion by the circular fit.

arguments
  x (1,:) {mustBeVector,mustBeNumeric}
  y (1,:) {mustBeVector,mustBeNumeric}
  xc (1,1) {mustBeNumeric}
  yc (1,1) {mustBeNumeric}
  R (1,1) {mustBePositive}
end

% Calculate the explained variance
sse = sum((vecnorm([x-xc y-yc]') - R).^2);
sst = sum(vecnorm([x-mean(x) y-mean(y)]').^2);
explainedVariance = max([0 1-sse/sst]);
